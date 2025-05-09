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
    
    if (d_project_name != "2D smooth Isentropic Rayleigh-Taylor instability") 
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Can only initialize data for 'project_name' = "
            << "'2D smooth Isentropic Rayleigh-Taylor instability' \n"
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
        double* E       = total_energy->getPointer(0);
        
        const double gamma = double(7)/double(5); // assume both gases have the same ratio of specific heat ratios
        // const double gamma_0 = double(7)/double(5);
        // const double gamma_1 = double(7)/double(5);
        
        double lambda = 701.53278340668; // wavelength of single-mode perturbation
        //double eta_0  = 0.02*lambda;      // 1% perturbation // DEBUGGING

        TBOX_ASSERT(d_initial_conditions_db != nullptr);
        TBOX_ASSERT(d_initial_conditions_db->keyExists("gravity"));
        TBOX_ASSERT(d_initial_conditions_db->keyExists("species_mass"));

        TBOX_ASSERT(d_initial_conditions_db->keyExists("eta_const"));
        const double eta_0 = d_initial_conditions_db->getDouble("eta_const");   // height of perturbation
        
        TBOX_ASSERT(d_initial_conditions_db->keyExists("delta_const"));
        const double delta = d_initial_conditions_db->getDouble("delta_const"); // characteristic length of interface.

        const double p_i = 100000.0; // interface pressure
        const double T_0 = 300.0;    // interface temperature, left as T_0

        std::vector<double> gravity_vector = d_initial_conditions_db->getDoubleVector("gravity");
        const double g   = gravity_vector[0]; // gravity
        std::vector<double> W_vector = d_initial_conditions_db->getDoubleVector("species_mass"); // molecular mass of mixing fluids                
        const double W_1 = W_vector[0]; // molecular mass of heavier gas
        const double W_2 = W_vector[1]; // molecular mass of lighter gas
        
        const double R_u = 8.31446261815324; // universal gas constant
        const double R_1 = R_u/W_1;          // gas constant of heavier gas
        const double R_2 = R_u/W_2;          // gas constant of lighter gas
        
        
        // const double rho_i = p_i/(R_u*T_0)*(W_1 + W_2)/2.0;
        
        if (d_project_name == "2D smooth Isentropic Rayleigh-Taylor instability")
        {
            const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());

            TBOX_ASSERT(d_initial_conditions_db->keyExists("x_lo"));
            TBOX_ASSERT(d_initial_conditions_db->keyExists("x_up"));
        
            std::vector<double> x_lo_vector = d_initial_conditions_db->getDoubleVector("x_lo");
            double x_domain_lo = x_lo_vector[0]; // Lower end of computational domain
            std::vector<double> x_hi_vector = d_initial_conditions_db->getDoubleVector("x_up");
            double x_domain_hi = x_hi_vector[0]; // Upper end of computational domain
            
            const double shift = 0.0; // location of interface.
            
            std::string integral_filename = "integral.dat";
            const int integral_N_x = 10000;
            const int integral_N_int = 1000000; // number of numerical quadrature points
            
            // Discretize the domain in x-direction for the approximated integral.
            
        
            //double x_domain_lo = -4.0*lambda; // Hard coded but read from input
            //double x_domain_hi =  4.0*lambda; // Hard coded but read from input
            
            x_domain_lo -= 0.1*lambda; // enlarge the domain for domain ghost cells
            x_domain_hi += 0.1*lambda; // enlarge the domain for domain ghost cells
            const double dx_uniform  = (x_domain_hi - x_domain_lo)/double(integral_N_x);
            
            std::vector<double> integral_vector(integral_N_x + 3);
            integral_vector[integral_N_x + 0] = x_domain_lo;
            integral_vector[integral_N_x + 1] = x_domain_hi;
            integral_vector[integral_N_x + 2] = dx_uniform;
            
            std::ifstream f_in;
            f_in.open(integral_filename, std::ios::in | std::ios::binary);
            
            if (!f_in.is_open())
            {
                for (int i = 0; i < integral_N_x; i++)
                {
                    integral_vector[i] = 0.0;
                    const double x_pos = i*dx_uniform + 0.5*dx_uniform + x_domain_lo;
                    const double dx_p = (x_pos - shift)/(double(integral_N_int) - 1.0);
                    for (int ii = 0; ii < integral_N_int; ii++)
                    {
                        const double x_p = shift + ii*dx_p;  //Bug fixed 3.22.2023
                        integral_vector[i] += 1.0/(0.5*(R_2 - R_1)*erf((x_p - shift)/(delta)) + 0.5*(R_1 + R_2))*dx_p;
                    }
                }
                
                mpi.Barrier();
                if (mpi.getRank() == 0)
                {
                    std::ofstream f_out;
                    f_out.open(integral_filename, std::ios::out | std::ios::binary);
                    if (!f_out.is_open())
                    {
                        TBOX_ERROR(d_object_name
                            << ": "
                            << "Failed to open file to output integral!"
                            << std::endl);
                    }
                    
                    f_out.write((char*)&integral_vector[0], sizeof(double)*integral_vector.size());
                    f_out.close();
                }
                mpi.Barrier();
            }
            else
            {
                f_in.read((char*)&integral_vector[0], sizeof(double)*integral_vector.size());
                f_in.close();
            }
            
            for (int j = 0; j < patch_dims[1]; j++)
            {
                for (int i = 0; i < patch_dims[0]; i++)
                {
                    // Compute index into linear data array.
                    int idx_cell = i + j*patch_dims[0];
                    
                    // Compute the coordinates.
                    double x[2];
                    x[0] = patch_xlo[0] + (double(i) + double(1)/double(2))*dx[0];
                    x[1] = patch_xlo[1] + (double(j) + double(1)/double(2))*dx[1];
                    
                    const double eta = eta_0*cos(2.0*M_PI/lambda*x[1]);
                    
                    double X_2_H = 0.5*(1.0 + erf((x[0] - eta - shift)/delta)); // mass fraction of second species (Y_2)
                    // if (std::abs(x[0] - shift)/delta > 5.0*delta)
                    // {
                    //     X_2_H = 0.5*(1.0 + erf((x[0] - shift)/delta));
                    // }
                    
                    const double R_H   = R_1*(1.0 - X_2_H) + X_2_H*R_2;
                    
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
                    
                    // for (int ii = 0; ii < N_int; ii++)
                    // {
                    //     // const double x_p = x[0] + ii*dx_p;  //Bug fixed 3.22.2023 OLD
                    //     const double x_p = shift + ii*dx_p;  //Bug fixed 3.22.2023
                    //     integral += 1.0/(0.5*(R_2 - R_1)*erf((x_p - shift)/(delta)) + 0.5*(R_1 + R_2))*dx_p;
                    // }
                    
                    p_H   = p_i*std::pow(((gamma-1.0)*g*integral/gamma/T_0)+1.0,gamma/(gamma-1.0));
                    rho_H = p_H*std::pow(p_H/p_i,(1.0-gamma)/gamma)/(R_H*T_0);
                    
                    // Scott's implementation
                    // const double dX_2_H_dx = 1.0/(delta*sqrt(M_PI))*exp(-(x[0]/delta)*(x[0]/delta));
                    // const double dlnR_H_dx = (R_2 - R_1)*dX_2_H_dx;
                    // const double p_H = p_i*exp(g/(R_H*T_0)*(x[0] - 0.5*delta*delta*dlnR_H_dx));
                    // const double rho_H = p_H/(R_H*T_0);
                    
                    // const double X_2 = 0.5*(1.0 + erf((x[0] - eta)/delta)); // mass fraction of second species (Y_2)
                    
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
