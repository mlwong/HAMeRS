#include "apps/Navier-Stokes/NavierStokesSpecialBoundaryConditions.hpp"

/*
 * Set the data on the patch physical boundary to some values, depending on the flow problems
 * and flow models.
 */
void
NavierStokesSpecialBoundaryConditions::setSpecialBoundaryConditions(
    hier::Patch& patch,
    const std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& conservative_variables,
    const double fill_time,
    const hier::IntVector& ghost_width_to_fill)
{
    if ((d_project_name != "2D smooth isopycnic Rayleigh-Taylor instability") &&  // single-mode 2 species
        (d_project_name != "2D smooth isopycnic Rayleigh-Taylor instability 3 species") && // single-mode 3 species
        (d_project_name != "2D smooth multi-mode isopycnic Rayleigh-Taylor instability")) // multi-mode 3 species
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Can only initialize data for 'project_name' = '2D smooth isopycnic Rayleigh-Taylor instability' or"
            << "'2D smooth isopycnic Rayleigh-Taylor instability 3 species' or '2D smooth multi-mode isopycnic Rayleigh-Taylor instability'!\n"
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
    
    if (d_flow_model_type != FLOW_MODEL::FOUR_EQN_CONSERVATIVE)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Flow model should be conservative four-equation models!"
            << std::endl);
    }
    
    if (d_flow_model->getNumberOfSpecies() != 2 && d_flow_model->getNumberOfSpecies() != 3)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Number of species should be 2 or 3!"
            << std::endl);
    }
    
    // Get the number of ghost cells of gradient.
    hier::IntVector num_ghosts = conservative_variables[0]->getGhostCellWidth();
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of patch plus
    // ghost cells.
    const hier::Box ghost_box = conservative_variables[0]->getGhostBox();
    const hier::IntVector ghostcell_dims = ghost_box.numberCells();
    
    const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
        HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(patch_geom);
#endif
    
    TBOX_ASSERT(d_special_boundary_conditions_db != nullptr);
    TBOX_ASSERT(d_special_boundary_conditions_db->keyExists("gravity"));
    TBOX_ASSERT(d_special_boundary_conditions_db->keyExists("species_mass"));
    
    std::vector<double> gravity_vector = d_special_boundary_conditions_db->getDoubleVector("gravity");
    std::vector<double> W_vector = d_special_boundary_conditions_db->getDoubleVector("species_mass"); // molecular mass of mixing fluids
    
    if( d_project_name == "2D smooth isopycnic Rayleigh-Taylor instability")
    {
        TBOX_ASSERT(static_cast<int>(W_vector.size()) == 2);
        
        if (patch_geom->getTouchesRegularBoundary(0, 0) ||
            patch_geom->getTouchesRegularBoundary(0, 1))
        {
            const double* const dx = patch_geom->getDx();
            const double* const patch_xlo = patch_geom->getXLower();
            
            HAMERS_SHARED_PTR<pdat::CellData<double> > partial_density = conservative_variables[0];
            HAMERS_SHARED_PTR<pdat::CellData<double> > momentum        = conservative_variables[1];
            HAMERS_SHARED_PTR<pdat::CellData<double> > total_energy    = conservative_variables[2];
            
            double* rho_Y_0 = partial_density->getPointer(0);
            double* rho_Y_1 = partial_density->getPointer(1);
            double* rho_u   = momentum->getPointer(0);
            double* rho_v   = momentum->getPointer(1);
            double* E       = total_energy->getPointer(0);
            
            const double gamma = double(7)/double(5);
            const double g = gravity_vector[0]; // gravity
            
            const double p_i = 100000.0; // interface pressure
            const double T_0 = 300.0; 
            
            const double W_1 = W_vector[0]; // molecular mass of heavier gas
            const double W_2 = W_vector[1]; // molecular mass of lighter gas
            
            const double R_u = 8.31446261815324; // universal gas constant
            const double R_1 = R_u/W_1;          // gas constant of heavier gas
            const double R_2 = R_u/W_2;
            
            double lambda = 701.53278340668; // wavelength of single-mode perturbation
            double eta_0  = 0.02*lambda;      // 1% perturbation // DEBUGGING
            
            const double delta = 0.04*lambda; // characteristic length of interface.
            const double shift = 0.0;
            const double rho_1 = p_i/(R_1*T_0);
            const double rho_2 = p_i/(R_2*T_0); 
            
            // Assume it is left boundary first.
            int i_lo = -ghost_width_to_fill[0];
            int i_hi = 0;
            
            for (int bi = 0; bi < 2; bi++) // loop over left and righ boundaries
            {
                if (patch_geom->getTouchesRegularBoundary(0, bi))
                {
                    if (bi == 1) // for right boundary
                    {
                        i_lo = interior_dims[0];
                        i_hi = interior_dims[0] + ghost_width_to_fill[0];
                    }
                    
                    for (int j = -num_ghosts[1]; j < interior_dims[1] + num_ghosts[1]; j++)
                    {
                        for (int i = i_lo; i < i_hi; i++)
                        {
                            const int idx_cell = (i + num_ghosts[0]) +
                                (j + num_ghosts[1])*ghostcell_dims[0];
                            
                            // Compute the coordinates.
                            double x[2];
                            x[0] = patch_xlo[0] + (i + double(1)/double(2))*dx[0];
                            x[1] = patch_xlo[1] + (j + double(1)/double(2))*dx[1];
                            
                            const double x_shifted = x[0] - shift;
                            
                            const double eta = eta_0*cos(2.0*M_PI/lambda*x[1]);
                            
                            const double Z_2_H = 0.5*(1.0 + erf((x_shifted - eta)/delta)); // volume fraction of second species (Z_2)
                            
                            
                            const double rho = rho_1*(1 - Z_2_H) + rho_2*Z_2_H;
                            
                            const double p_H = p_i + 0.5*(rho_1+rho_2)*g*(x_shifted) +
                                0.5*(rho_1-rho_2)*g*(delta*(exp(-pow(x_shifted/delta,2.0))-1.0)/sqrt(M_PI) + x_shifted*erf(x_shifted/delta));
                            
                            const double p = p_H;
                            
                            rho_Y_0[idx_cell] = rho_1*(1.0 - Z_2_H);
                            rho_Y_1[idx_cell] = rho_2*Z_2_H;
                            
                            const double u = 0.0;
                            const double v = 0.0;
                            
                            rho_u[idx_cell] = rho*u;
                            rho_v[idx_cell] = rho*v;
                            E[idx_cell]     = p/(gamma - double(1)) + double(1)/double(2)*rho*(u*u + v*v);
                        }
                    }
                }
            }
        }
    }
    else if(d_project_name == "2D smooth isopycnic Rayleigh-Taylor instability 3 species")
    {
        TBOX_ASSERT(static_cast<int>(W_vector.size()) == 3);
        
        if (patch_geom->getTouchesRegularBoundary(0, 0) ||
            patch_geom->getTouchesRegularBoundary(0, 1))
        {
            const double* const dx = patch_geom->getDx();
            const double* const patch_xlo = patch_geom->getXLower();
            
            HAMERS_SHARED_PTR<pdat::CellData<double> > partial_density = conservative_variables[0];
            HAMERS_SHARED_PTR<pdat::CellData<double> > momentum        = conservative_variables[1];
            HAMERS_SHARED_PTR<pdat::CellData<double> > total_energy    = conservative_variables[2];
            
            double* rho_Y_0 = partial_density->getPointer(0);
            double* rho_Y_1 = partial_density->getPointer(1);
            double* rho_Y_2 = partial_density->getPointer(2);
            double* rho_u   = momentum->getPointer(0);
            double* rho_v   = momentum->getPointer(1);
            double* E       = total_energy->getPointer(0);
            
            const double gamma = double(7)/double(5);
            const double g     = gravity_vector[0]; // gravity
            
            const double p_i = 100000.0; // interface pressure
            const double T_0 = 300.0; 
            
            const double W_1 = W_vector[0]; // molecular mass of first gas
            const double W_2 = W_vector[1]; // molecular mass of second gas
            const double W_3 = W_vector[2]; // molecular mass of third gas
            
            const double R_u = 8.31446261815324; // universal gas constant
            const double R_1 = R_u/W_1;          // gas constant of heavier gas
            const double R_2 = R_u/W_2;
            const double R_3 = R_u/W_3;
            
            
            double lambda = 701.53278340668; // wavelength of single-mode perturbation
            double eta_0  = 0.02*lambda;      // 1% perturbation // DEBUGGING
            
            const double delta = 0.01*lambda; // characteristic length of interface.
            const double shift = lambda/4.0;
            const double rho_1 = p_i/(R_1*T_0);
            const double rho_2 = p_i/(R_2*T_0);
            const double rho_3 = p_i/(R_3*T_0);
            
            // Assume it is left boundary first.
            int i_lo = -ghost_width_to_fill[0];
            int i_hi = 0;
            
            for (int bi = 0; bi < 2; bi++) // loop over left and righ boundaries
            {
                if (patch_geom->getTouchesRegularBoundary(0, bi))
                {
                    if (bi == 1) // for right boundary
                    {
                        i_lo = interior_dims[0];
                        i_hi = interior_dims[0] + ghost_width_to_fill[0];
                    }
                    
                    for (int j = -num_ghosts[1]; j < interior_dims[1] + num_ghosts[1]; j++)
                    {
                        for (int i = i_lo; i < i_hi; i++)
                        {
                            const int idx_cell = (i + num_ghosts[0]) +
                                (j + num_ghosts[1])*ghostcell_dims[0];
                            
                            // Compute the coordinates.
                            double x[2];
                            x[0] = patch_xlo[0] + (i + double(1)/double(2))*dx[0];
                            x[1] = patch_xlo[1] + (j + double(1)/double(2))*dx[1];
                            
                            const double eta   = eta_0*cos(2.0*M_PI/lambda*x[1])*0.0;
                            
                            const double Z_2_H = 0.5*(1.0 + erf(((x[0] - eta + shift)/delta))) - 0.5*(1.0 + erf(((x[0] - eta - shift)/delta)));
                            const double Z_3_H = 0.5*(1.0 + erf(((x[0] - eta - shift)/delta)));
                            const double Z_1_H = 1.0 - Z_2_H - Z_3_H;
                            
                            const double rho = rho_1*Z_1_H + rho_2*Z_2_H + rho_3*Z_3_H;
                            
                            const double ksi_1 = (x[0] + shift)/delta;
                            const double ksi_2 = (x[0] - shift)/delta;
                            
                            const double p = p_i + 0.5*g*(rho_1 + rho_3)*x[0] +
                                0.5*g*delta*(rho_2 - rho_1)*( -(shift/delta)*erf(shift/delta) + ksi_1*erf(ksi_1) + (exp(-ksi_1*ksi_1) - exp(-shift*shift/(delta*delta)))/sqrt(M_PI) ) +
                                0.5*g*delta*(rho_3 - rho_2)*( -(shift/delta)*erf(shift/delta) + ksi_2*erf(ksi_2) + (exp(-ksi_2*ksi_2) - exp(-shift*shift/(delta*delta)))/sqrt(M_PI) );
                            
                            rho_Y_0[idx_cell] = rho_1*Z_1_H;
                            rho_Y_1[idx_cell] = rho_2*Z_2_H;
                            rho_Y_2[idx_cell] = rho_3*Z_3_H;
                            
                            const double u = 0.0;
                            const double v = 0.0;
                            
                            rho_u[idx_cell] = rho*u;
                            rho_v[idx_cell] = rho*v;
                            E[idx_cell]     = p/(gamma - double(1)) + double(1)/double(2)*rho*(u*u + v*v);
                        }
                    }
                }
            }
        }
    }
    else if(d_project_name == "2D smooth multi-mode isopycnic Rayleigh-Taylor instability")
    {
        TBOX_ASSERT(static_cast<int>(W_vector.size()) == 3);
        
        if (patch_geom->getTouchesRegularBoundary(0, 0) ||
            patch_geom->getTouchesRegularBoundary(0, 1))
        {
            const double* const dx = patch_geom->getDx();
            const double* const patch_xlo = patch_geom->getXLower();
            
            HAMERS_SHARED_PTR<pdat::CellData<double> > partial_density = conservative_variables[0];
            HAMERS_SHARED_PTR<pdat::CellData<double> > momentum        = conservative_variables[1];
            HAMERS_SHARED_PTR<pdat::CellData<double> > total_energy    = conservative_variables[2];
            
            double* rho_Y_0 = partial_density->getPointer(0);
            double* rho_Y_1 = partial_density->getPointer(1);
            double* rho_Y_2 = partial_density->getPointer(2);
            double* rho_u   = momentum->getPointer(0);
            double* rho_v   = momentum->getPointer(1);
            double* E       = total_energy->getPointer(0);
            
            const double gamma = double(7)/double(5);
            const double g     = gravity_vector[0]; // gravity
            
            const double p_i = 100000.0; // interface pressure
            const double T_0 = 300.0; 
            
            const double W_1 = W_vector[0]; // molecular mass of first gas
            const double W_2 = W_vector[1]; // molecular mass of second gas
            const double W_3 = W_vector[2]; // molecular mass of third gas
            
            const double R_u = 8.31446261815324; // universal gas constant
            const double R_1 = R_u/W_1;          // gas constant of heavier gas
            const double R_2 = R_u/W_2;
            const double R_3 = R_u/W_3;
            
            const double rho_1 = p_i/(R_1*T_0);
            const double rho_2 = p_i/(R_2*T_0);
            const double rho_3 = p_i/(R_3*T_0);

            const double lambda  = 701.53278340668/4.0;
            const double eta_0   = lambda*0.04;
            const double delta   = 0.04*lambda; // characteristic length of interface
            const int    waven   = 16;          // dominant wave number
            const double width   = 16.0*lambda; // domain size in y direction
            const double shift = lambda/4.0;
            
            double rmod[9]; // random seed
            rmod[0] = 6.031966614958411e+000;
            rmod[1] = 1.273017034173460e+000;
            rmod[2] = 5.934447177754063e+000;
            rmod[3] = 3.101658133166612e+000;
            rmod[4] = 2.294026034817427e+000;
            rmod[5] = 4.916046917518752e+000;
            rmod[6] = 0.571212135466553e+000;
            rmod[7] = 4.966766749458944e+000;
            rmod[8] = 5.027899324302027e+000;
            
            double rmod_2[9]; // random seed
            rmod_2[0] = 2.620226532717789200e+000;
            rmod_2[1] = 4.525932273597345700e+000;
            rmod_2[2] = 7.186381718527406600e-004;
            rmod_2[3] = 1.899611578242180700e+000;
            rmod_2[4] = 9.220944569241362700e-001;
            rmod_2[5] = 5.801805019369201700e-001;
            rmod_2[6] = 1.170307423440345900e+000;
            rmod_2[7] = 2.171222082895173200e+000;
            rmod_2[8] = 2.492963564452900500e+000;
            
            // Assume it is left boundary first.
            int i_lo = -ghost_width_to_fill[0];
            int i_hi = 0;
            
            for (int bi = 0; bi < 2; bi++) // loop over left and righ boundaries
            {
                if (patch_geom->getTouchesRegularBoundary(0, bi))
                {
                    if (bi == 1) // for right boundary
                    {
                        i_lo = interior_dims[0];
                        i_hi = interior_dims[0] + ghost_width_to_fill[0];
                    }
            
                    for (int j = -num_ghosts[1]; j < interior_dims[1] + num_ghosts[1]; j++)
                    {
                        for (int i = i_lo; i < i_hi; i++)
                        {
                            const int idx_cell = (i + num_ghosts[0]) +
                                (j + num_ghosts[1])*ghostcell_dims[0];
                            
                            // Compute the coordinates.
                            double x[2];
                            x[0] = patch_xlo[0] + (i + double(1)/double(2))*dx[0];
                            x[1] = patch_xlo[1] + (j + double(1)/double(2))*dx[1];
                            
                            double eta = 0.0;
                            for (int m = waven - 4; m <= waven + 4; m++)
                            {
                                eta += eta_0/3.0*cos(2.0*M_PI*m/width*x[1] + rmod[m-waven+4]);
                            }
                            double eta_2 = 0.0;
                            for (int m = waven - 4; m <= waven + 4; m++)
                            {
                                eta_2 += eta_0/3.0*cos(2.0*M_PI*m/width*x[1] + rmod_2[m-waven+4]);
                            }
                            
                            const double Z_2_H = 0.5*(1.0 + erf(((x[0] - eta + shift)/delta))) - 0.5*(1.0 + erf(((x[0] - eta_2 - shift)/delta)));
                            const double Z_3_H = 0.5*(1.0 + erf(((x[0] - eta_2 - shift)/delta)));
                            const double Z_1_H = 1.0 - Z_2_H - Z_3_H;
                            
                            const double rho = rho_1*Z_1_H + rho_2*Z_2_H + rho_3*Z_3_H;
                            
                            const double ksi_1 = (x[0] + shift)/delta;
                            const double ksi_2 = (x[0] - shift)/delta;
                            
                            const double p = p_i + 0.5*g*(rho_1 + rho_3)*x[0] +
                                0.5*g*delta*(rho_2 - rho_1)*( -(shift/delta)*erf(shift/delta) + ksi_1*erf(ksi_1) + (exp(-ksi_1*ksi_1) - exp(-shift*shift/(delta*delta)))/sqrt(M_PI) ) +
                                0.5*g*delta*(rho_3 - rho_2)*( -(shift/delta)*erf(shift/delta) + ksi_2*erf(ksi_2) + (exp(-ksi_2*ksi_2) - exp(-shift*shift/(delta*delta)))/sqrt(M_PI) );
                            
                            rho_Y_0[idx_cell] = rho*Z_1_H;
                            rho_Y_1[idx_cell] = rho*Z_2_H;
                            rho_Y_2[idx_cell] = rho*Z_3_H;
                            
                            const double u = 0.0;
                            const double v = 0.0;
                            
                            rho_u[idx_cell] = rho*u;
                            rho_v[idx_cell] = rho*v;
                            E[idx_cell]     = p/(gamma - double(1)) + double(1)/double(2)*rho*(u*u + v*v);
                        }
                    }
                }
            }
        }
    }
}
