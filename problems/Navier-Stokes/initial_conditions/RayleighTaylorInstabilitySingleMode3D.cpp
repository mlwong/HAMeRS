#include "apps/Navier-Stokes/NavierStokesInitialConditions.hpp"

/*
 * Set the data on the patch interior to some initial values.
 */
void
NavierStokesInitialConditions::initializeDataOnPatch(
    hier::Patch& patch,
    const std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& conservative_variables,
    const double data_time,
    const bool initial_time)
{
    // Follow Reckinger, Scott J., Daniel Livescu, and Oleg V. Vasilyev.
    // "Comprehensive numerical methodology for direct numerical simulations of compressible Rayleigh–Taylor instability."
    // Journal of Computational Physics 313 (2016): 181-208.
    // Note that the sign of gravity in the paper is flipped.
    NULL_USE(data_time);
    
    if ((d_project_name != "3D smooth Rayleigh-Taylor instability")
       ) 
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Can only initialize data for 'project_name' = '3D smooth Rayleigh-Taylor instability' !\n "
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
    
    if (initial_time)
    {
        const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
            HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                patch.getPatchGeometry()));
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(patch_geom);
#endif
        
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
	double* rho_w   = momentum->getPointer(2);
        double* E       = total_energy->getPointer(0);
        
        const double gamma = double(7)/double(5); // assume both gases have the same ratio of specific heat ratios
        // const double gamma_0 = double(7)/double(5);
        // const double gamma_1 = double(7)/double(5);
        
              double lambda = 701.53278340668; // wavelength of single-mode perturbation
              double eta_0  = 0.01*lambda;      // 1% perturbation
        // const double eta_0  = 0.0*lambda;      // no perturbation
        
        const double W_1 = 0.04000; // molecular weight of heavier gas
        const double W_2 = 0.02400; // molecular weight of lighter gas
        
        const double p_i = 100000.0; // interface pressure
        const double T_0 = 300.0;    // background temperature
        
        TBOX_ASSERT(d_initial_conditions_db != nullptr);
        TBOX_ASSERT(d_initial_conditions_db->keyExists("gravity"));
        
        std::vector<double> gravity_vector = d_initial_conditions_db->getDoubleVector("gravity");
        const double g = gravity_vector[0]; // gravity
        
        const double R_u = 8.31446261815324; // universal gas constant
        const double R_1 = R_u/W_1;          // gas constant of heavier gas
        const double R_2 = R_u/W_2;          // gas constant of lighter gas
        
        // const double rho_i = p_i/(R_u*T_0)*(W_1 + W_2)/2.0;
        
        if (d_project_name == "3D smooth Rayleigh-Taylor instability")
        {
            const double delta = 0.01*lambda; // characteristic length of interface.
            
            for (int k = 0; k < patch_dims[2]; k++)
            {
                for (int j = 0; j < patch_dims[1]; j++)
                {
		    for (int i = 0; i < patch_dims[0]; i++)
		    {
                        // Compute index into linear data array.
                        int idx_cell = i + j*patch_dims[0] +  k*patch_dims[0]*patch_dims[1];
                        
                        // Compute the coordinates.
                        double x[3];
                        x[0] = patch_xlo[0] + (double(i) + double(1)/double(2))*dx[0];
                        x[1] = patch_xlo[1] + (double(j) + double(1)/double(2))*dx[1];
		        x[2] = patch_xlo[2] + (double(k) + double(1)/double(2))*dx[2];
                        
                        const double eta = eta_0*cos(2.0*M_PI/lambda*x[1])*cos(2.0*M_PI/lambda*x[2]);
                        
                        const double X_2_H = 0.5*(1.0 + erf((x[0] - eta)/delta)); // mass fraction of second species (Y_2)
                        const double R_H   = R_1*(1.0 - X_2_H) + X_2_H*R_2;
                        
                        const int N_int = 10000; // number of numerical quadrature points
                        const double dx_p = x[0]/(N_int - 1.0);
                        
                        double integral = 0.0;
                        for (int ii = 0; ii < N_int; ii++)
                        {
                            const double x_p = x[0] + ii*dx_p;
                            integral += 1.0/(0.5*(R_2 - R_1)*erf((x_p - eta)/delta) + 0.5*(R_1 + R_2))*dx_p;
                        }
                        
                        const double p_H = p_i*exp(g/T_0*integral);
                        const double rho_H = p_H/(R_H*T_0);
                        
                        // Scott's implementation
                        // const double dX_2_H_dx = 1.0/(delta*sqrt(M_PI))*exp(-(x[0]/delta)*(x[0]/delta));
                        // const double dlnR_H_dx = (R_2 - R_1)*dX_2_H_dx;
                        // const double p_H = p_i*exp(g/(R_H*T_0)*(x[0] - 0.5*delta*delta*dlnR_H_dx));
                        // const double rho_H = p_H/(R_H*T_0);
                        
                        // const double X_2 = 0.5*(1.0 + erf((x[0] - eta)/delta)); // mass fraction of second species (Y_2)
                        
                        double rho, p;
                        
                        rho = rho_H;
                        p   = p_H;
                        
                        // if (x[0] < eta)
                        // {
                        //     const double p_1_H   = p_i*exp((g*x[0])/(R_1*T_0));
                        //     const double rho_1_H = p_i/(R_1*T_0)*exp((g*x[0])/(R_1*T_0));
                        //     
                        //     const double p_1   = p_i*exp((g*x[0])/(R_1*T_0));
                        //     const double rho_1 = p_i/(R_1*T_0)*exp((g*x[0])/(R_1*T_0));
                        //     
                        //     p   = p_1*p_H/p_1_H;
                        //     rho = rho_1*rho_H/rho_1_H;
                        //     
                        //     
                        // }
                        // else
                        // {
                        //     const double p_2_H   = p_i*exp((g*x[0])/(R_2*T_0));
                        //     const double rho_2_H = p_i/(R_2*T_0)*exp((g*x[0])/(R_2*T_0));
                        //     
                        //     const double p_2   = p_i*exp((g*x[0])/(R_2*T_0));
                        //     const double rho_2 = p_i/(R_2*T_0)*exp((g*x[0])/(R_2*T_0));
                        //     
                        //     p   = p_2*p_H/p_2_H;
                        //     rho = rho_2*rho_H/rho_2_H;
                        // }
                        
                        rho_Y_0[idx_cell] = rho*(1.0 - X_2_H);
                        rho_Y_1[idx_cell] = rho*X_2_H;
                        
                        const double u = 0.0;
                        const double v = 0.0;
		        const double w = 0.0;
                        
                        rho_u[idx_cell] = rho*u;
                        rho_v[idx_cell] = rho*v;
		        rho_w[idx_cell] = rho*w;
                        E[idx_cell]     = p/(gamma - double(1)) + double(1)/double(2)*rho*(u*u + v*v + w*w);
 
     		    }
		}
            }
        }
    }
}
