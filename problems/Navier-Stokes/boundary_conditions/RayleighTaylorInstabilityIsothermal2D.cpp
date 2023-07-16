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
    if ((d_project_name != "2D discontinuous Rayleigh-Taylor instability") && (d_project_name != "2D smooth Rayleigh-Taylor instability"))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Can only initialize data for 'project_name' = '2D discontinuous Rayleigh-Taylor instability' or"
            << "'2D smooth Rayleigh-Taylor instability'!\n"
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
    
    if (d_flow_model->getNumberOfSpecies() != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Number of species should be 2!"
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
    
    if( d_project_name == "2D discontinuous Rayleigh-Taylor instability")
    {
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
            const double g     = 90.0; // HARD-CODED for now!
            
            const double p_i = 100000.0; // interface pressure
            const double T_0 = 300.0; 
            
            const double W_1 = 0.0600; //0.033280; //0.0400; // molecular mass of heavier gas
            const double W_2 = 0.0200; //0.030720; //0.0240; // molecular mass of lighter gas
            
            const double R_u = 8.31446261815324; // universal gas constant
            const double R_1 = R_u/W_1;          // gas constant of heavier gas
            const double R_2 = R_u/W_2;
            
            /*
             * Update the left boundary conditions.
             */
            
            if (patch_geom->getTouchesRegularBoundary(0, 0))
            {
                for (int j = -num_ghosts[1]; j < interior_dims[1] + num_ghosts[1]; j++)
                {
                    for (int i = -ghost_width_to_fill[0];
                         i < 0;
                         i++)
                    {
                        const int idx_cell = (i + num_ghosts[0]) +
                            (j + num_ghosts[1])*ghostcell_dims[0];
                        
                        // Compute the coordinates.
                        double x[2];
                        x[0] = patch_xlo[0] + (i + double(1)/double(2))*dx[0];
                        x[1] = patch_xlo[1] + (j + double(1)/double(2))*dx[1];
                        
                        const double rho = p_i/(R_1*T_0)*exp((g*x[0])/(R_1*T_0));
                        rho_Y_0[idx_cell] = rho;
                        rho_Y_1[idx_cell] = 0.0;
                        
                        const double p = p_i*exp((g*x[0])/(R_1*T_0));
                            
                        const double u = 0.0;
                        const double v = 0.0;
                            
                        rho_u[idx_cell] = rho*u;
                        rho_v[idx_cell] = rho*v;
                        E[idx_cell]     = p/(gamma - double(1)) + double(1)/double(2)*rho*(u*u + v*v);
                    }
                }
            }
            
            /*
             * Update the right boundary conditions.
             */
            
            if (patch_geom->getTouchesRegularBoundary(0, 1))
            {
                for (int j = -num_ghosts[1]; j < interior_dims[1] + num_ghosts[1]; j++)
                {
                    for (int i = interior_dims[0];
                         i < interior_dims[0] + ghost_width_to_fill[0];
                         i++)
                    {
                        const int idx_cell = (i + num_ghosts[0]) +
                            (j + num_ghosts[1])*ghostcell_dims[0];
                        
                        // Compute the coordinates.
                        double x[2];
                        x[0] = patch_xlo[0] + (i + double(1)/double(2))*dx[0];
                        x[1] = patch_xlo[1] + (j + double(1)/double(2))*dx[1];
                        
                        const double rho = p_i/(R_2*T_0)*exp((g*x[0])/(R_2*T_0));
                        rho_Y_0[idx_cell] = 0.0;
                        rho_Y_1[idx_cell] = rho;
                            
                        const double p = p_i*exp((g*x[0])/(R_2*T_0));
                            
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
    else if(d_project_name == "2D smooth Rayleigh-Taylor instability")
    {
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
            const double g     = 90.0; // HARD-CODED for now!
            
            const double p_i = 100000.0; // interface pressure
            const double T_0 = 300.0; 
            
            const double W_1 = 0.060; // molecular mass of heavier gas
            const double W_2 = 0.020; // molecular mass of lighter gas
            
            const double R_u = 8.31446261815324; // universal gas constant
            const double R_1 = R_u/W_1;          // gas constant of heavier gas
            const double R_2 = R_u/W_2;
            
            const double lambda = 701.53278340668; // wavelength of single-mode perturbation
            const double eta_0  = 0.04*lambda;
            const double delta  = 0.04*lambda; // characteristic length of interface.
            const double shift  = 0.0; // location of interface.
            
            
            // Assume it is left boundary first.
            int i_lo = -ghost_width_to_fill[0];
            int i_hi = 0;
            
            for (int bi = 0; bi < 2; bi++) // loop over left and righ boundaries
            {
                if (patch_geom->getTouchesRegularBoundary(0, bi))
                {
                    // Discretize the domain in x-direction for the approximated integral.
                    const int integral_N_x = 10000;
                    
                    std::vector<double> integral_vector(integral_N_x + 3);
                    
                    std::ifstream f_in;
                    std::string integral_filename = "integral.dat";
                    f_in.open(integral_filename, std::ios::in | std::ios::binary);
                    f_in.read((char*)&integral_vector[0], sizeof(double)*integral_vector.size());
                    f_in.close();
                    
                    const double x_domain_lo = integral_vector[integral_N_x + 0];
                    const double x_domain_hi = integral_vector[integral_N_x + 1];
                    const double dx_uniform  = integral_vector[integral_N_x + 2];
                    
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
                            
                            const double eta = eta_0*cos(2.0*M_PI/lambda*x[1]);
                            
                            double X_2_H = 0.5*(1.0 + erf((x[0] - eta - shift)/delta)); // mass fraction of second species (Y_2)
                            
                            const double R_H   = R_1*(1.0 - X_2_H) + X_2_H*R_2;
                            
                            // const int N_int = 100000; // number of numerical quadrature points
                            // const double dx_p = (x[0] - shift)/(N_int - 1.0);
                            
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
                    
                            
                            // for (int ii = 0; ii < N_int; ii++)
                            //{
                            //    const double x_p = shift + ii*dx_p;  //Bug fixed 3.22.2023
                            //    integral += 1.0/(0.5*(R_2 - R_1)*erf((x_p - shift)/(delta)) + 0.5*(R_1 + R_2))*dx_p;
                            //}
                            p_H = p_i*exp(g/T_0*integral);
                            rho_H = p_H/(R_H*T_0);
                            
                            double rho, p;
                            
                            rho = rho_H;
                            p   = p_H;
                            
                            rho_Y_0[idx_cell] = rho*(1.0 - X_2_H);
                            rho_Y_1[idx_cell] = rho*X_2_H;
                            
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
