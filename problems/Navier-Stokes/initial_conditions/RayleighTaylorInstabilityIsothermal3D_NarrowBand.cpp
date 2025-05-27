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
    
    if ((d_project_name != "3D discontinuous Rayleigh-Taylor instability") &&
        (d_project_name != "3D smooth Rayleigh-Taylor instability") &&
        (d_project_name != "3D smooth multi-mode Rayleigh-Taylor instability")
       ) 
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Can only initialize data for 'project_name' = '3D discontinuous Rayleigh-Taylor instability' or "
            << "'3D smooth Rayleigh-Taylor instability' or "
            << "'3D smooth multi-mode Rayleigh-Taylor instability'!\n"
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
        
        // const double gamma = double(7)/double(5); // assume both gases have the same ratio of specific heat ratios, hard coded 

        
        TBOX_ASSERT(d_initial_conditions_db != nullptr);
        TBOX_ASSERT(d_initial_conditions_db->keyExists("gravity"));
        TBOX_ASSERT(d_initial_conditions_db->keyExists("species_mass"));
        TBOX_ASSERT(d_initial_conditions_db->keyExists("species_gamma"));
        TBOX_ASSERT(d_initial_conditions_db->keyExists("width"));
        TBOX_ASSERT(d_initial_conditions_db->keyExists("delta"));
        TBOX_ASSERT(d_initial_conditions_db->keyExists("eta_0"));
        TBOX_ASSERT(d_initial_conditions_db->keyExists("m_min"));
        TBOX_ASSERT(d_initial_conditions_db->keyExists("m_max"));
        
        std::vector<double> gravity_vector = d_initial_conditions_db->getDoubleVector("gravity"); //gravity 
        const double g = gravity_vector[0]; // gravity
        std::vector<double> W_vector = d_initial_conditions_db->getDoubleVector("species_mass"); // molecular mass of mixing fluids                
        const double W_1 = W_vector[0]; // molecular mass of heavier gas
        const double W_2 = W_vector[1]; // molecular mass of lighter gas
        std::vector<double> gamma_vector = d_initial_conditions_db->getDoubleVector("species_gamma"); // specific heat ratio 
        const double gamma = gamma_vector[0]; // specific heat ratio for heavier gas (assume both gases have same specific heat ratios)

        // Read domain size in y direction from input
        const double width = d_initial_conditions_db->getDouble("width"); // width in y direction (assumes Ly = Lz)
            
        // Read characteristic thickness of interface from input
        const double delta = d_initial_conditions_db->getDouble("delta");
        
        // Read characteristic length of interface from input
        const double eta_0 = d_initial_conditions_db->getDouble("eta_0");

        // Read minimum wave number from input
        const int m_min = d_initial_conditions_db->getInteger("m_min");
            
        // Read maximum wave number from input
        const int m_max = d_initial_conditions_db->getInteger("m_max");

        const double domain_width = width;  // wavelength of single-mode perturbation
        const double dominant_m   = (m_min + m_max)/2.0;            // Perturbed between wave numbers (k) 4 to 12
        const double lambda_0     = domain_width/dominant_m;
        

        // std::cout << "domain_width: " << domain_width << ", dominant_m: " << dominant_m << ", lambda_0: " << lambda_0 << ", eta_0: " << eta_0 << std::endl;
        
        const double p_i = 100000.0; // interface pressure
        const double T_0 = 300.0;    // background temperature
        
        
        const double R_u = 8.31446261815324; // universal gas constant
        const double R_1 = R_u/W_1;          // gas constant of heavier gas
        const double R_2 = R_u/W_2;          // gas constant of lighter gas
        
        // const double rho_i = p_i/(R_u*T_0)*(W_1 + W_2)/2.0;
        
        if (d_project_name == "3D discontinuous Rayleigh-Taylor instability")
        {

        }
        else if (d_project_name == "3D smooth Rayleigh-Taylor instability")
        {
            
        }
        else if (d_project_name == "3D smooth multi-mode Rayleigh-Taylor instability")
        {
            const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());

            TBOX_ASSERT(d_initial_conditions_db->keyExists("x_lo"));
            TBOX_ASSERT(d_initial_conditions_db->keyExists("x_up"));
        
            std::vector<double> x_lo_vector = d_initial_conditions_db->getDoubleVector("x_lo");
            double x_domain_lo = x_lo_vector[0]; // Lower end of computational domain
            std::vector<double> x_hi_vector = d_initial_conditions_db->getDoubleVector("x_up");
            double x_domain_hi = x_hi_vector[0]; // Upper end of computational domain
            

            const double shift   = 0.0; // location of interface.


            const double k_0     = 2.0*M_PI/width;
            const double delta_ky   = 1*k_0;
            const double delta_kz   = 1*k_0;
            const double epsilon = 1.0e-15;
            const double k_max   = m_max*k_0;
            const double k_min   = m_min*k_0;
            
            std::vector<double> a_mn((2*m_max+1) * (2*m_max+1)); // Create vector for random numbers a
            // Read random number list a
            std::ifstream f_rand_a;
            std::string random_a_filename = "random_a.dat";
            f_rand_a.open(random_a_filename, std::ios::in | std::ios::binary);
            f_rand_a.read((char*)&a_mn[0], sizeof(double)*a_mn.size());
            f_rand_a.close();

            std::vector<double> b_mn((2*m_max+1) * (2*m_max+1)); // Create vector for random numbers b
            // Read random number list b
            std::ifstream f_rand_b;
            std::string random_b_filename = "random_b.dat";
            f_rand_b.open(random_b_filename, std::ios::in | std::ios::binary);
            f_rand_b.read((char*)&b_mn[0], sizeof(double)*b_mn.size());
            f_rand_b.close();

            std::vector<double> c_mn((2*m_max+1) * (2*m_max+1)); // Create vector for random numbers c
            // Read random number list c
            std::ifstream f_rand_c;
            std::string random_c_filename = "random_c.dat";
            f_rand_c.open(random_c_filename, std::ios::in | std::ios::binary);
            f_rand_c.read((char*)&c_mn[0], sizeof(double)*c_mn.size());
            f_rand_c.close();
        
            std::vector<double> d_mn((2*m_max+1) * (2*m_max+1)); // Create vector for random numbers d
            // Read random number list d
            std::ifstream f_rand_d;
            std::string random_d_filename = "random_d.dat";
            f_rand_d.open(random_d_filename, std::ios::in | std::ios::binary);
            f_rand_d.read((char*)&d_mn[0], sizeof(double)*d_mn.size());
            f_rand_d.close();
        
           
            std::string integral_filename = "integral.dat";
            const int integral_N_x   = 10000;
            const int integral_N_int = 1000000; // number of numerical quadrature points
            
            // Discretize the domain in x-direction for the approximated integral.
            
            x_domain_lo -= 1.1*width; // enlarge the domain for domain ghost cells
            x_domain_hi += 1.1*width; // enlarge the domain for domain ghost cells
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
                    const double dx_p  = (x_pos - shift)/(double(integral_N_int) - 1.0);
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
         
            // Initial perturbation follows:  B. THORNBER1‡, D. DRIKAKIS1, D. L. YOUNGS2 AND R. J. R. WILLIAMS 2010
            // The influence of initial conditions on turbulent mixing due to Richtmyer–Meshkov instability
            // Appendix A

            for (int k = 0; k < patch_dims[2]; k++)
            {
                for (int j = 0; j < patch_dims[1]; j++)
                {
                    for (int i = 0; i < patch_dims[0]; i++)
                    {
                        // Compute index into linear data array.
                        int idx_cell = i + j*patch_dims[0] + k*patch_dims[0]*patch_dims[1];
                        
                        // Compute the coordinates.
                        double x[3];
                        x[0] = patch_xlo[0] + (double(i) + double(1)/double(2))*dx[0];
                        x[1] = patch_xlo[1] + (double(j) + double(1)/double(2))*dx[1];
                        x[2] = patch_xlo[2] + (double(k) + double(1)/double(2))*dx[2];

                        // Initialize the parameters

                        double eta      = 0.0;

                        for (int m = -m_max; m <= m_max; m++)
                        {
                            for (int n = -m_max; n <= m_max; n++)
                            {
                                if (m*m + n*n > m_max*m_max || m*m + n*n < m_min*m_min)
                                {                                 
                                }
                                else
                                {   
                                    int idx_random = (m + m_max) * (2*m_max+1) + (n + m_max) ;
                                    double ratio   = sqrt((1.0 / 4.0 * ( a_mn[idx_random]*a_mn[idx_random] +
                                                           b_mn[idx_random]*b_mn[idx_random] +
                                                           c_mn[idx_random]*c_mn[idx_random] +
                                                           d_mn[idx_random]*d_mn[idx_random] ) *
                                               2.0 * M_PI * sqrt((m * m + n * n)*(k_0*k_0))) / (eta_0 * eta_0  * delta_ky * delta_kz / k_0)); 
                                               // eta_0 is P in the eq. A6 of the ref. paper
                                
                                    if (ratio <= epsilon) {
                                        TBOX_ERROR(d_object_name
                                            << ": "
                                            << "Can not initialize when a = b = c = d = 0."
                                            << std::endl);
                                    }
                                    else {
                                        ratio = 1.0/ratio;
                                    }

                                    // std::cout << "ratio: " << ratio << std::endl;

                                    eta += ratio * (a_mn[idx_random]*cos(k_0 * m * x[1])*cos(k_0 * n * x[2]) +
                                                      b_mn[idx_random]*cos(k_0 * m * x[1])*sin(k_0 * n * x[2]) +
                                                      c_mn[idx_random]*sin(k_0 * m * x[1])*cos(k_0 * n * x[2]) +
                                                      d_mn[idx_random]*sin(k_0 * m * x[1])*sin(k_0 * n * x[2]));
                                    // std::cout << "m: " << m << ", n: " << n << ", idx_random: " << idx_random << ", ratio: " << ratio << ", eta: " << eta << ", amn: " << a_mn[idx_random] << ", bmn: " << b_mn[idx_random] << ", cmn: " << c_mn[idx_random] << ", dmn: " << d_mn[idx_random] << std::endl;
                                }
                            }
                        }
                        
                        // std::cout << "delta: " << delta << std::endl;
                        
                        const double X_2_H  = 0.5*(1.0 + erf((x[0] - eta)/delta)); // mass fraction of second species (Y_2)
                        const double R_H    = R_1*(1.0 - X_2_H) + X_2_H*R_2;
                        
                        const int N_int     = 10000; // number of numerical quadrature points
                        const double dx_p   = x[0]/(N_int - 1.0);
                        

                        // Compute the integral with linear interpolation.
                        const int idx_integral_vector_lo = int(std::floor((x[0] - x_domain_lo)/dx_uniform));
                        const int idx_integral_vector_hi = idx_integral_vector_lo + 1;
                        
                        const double x_integral_vector_lo = double(idx_integral_vector_lo)*dx_uniform + 0.5*dx_uniform + x_domain_lo;
                        const double x_integral_vector_hi = x_integral_vector_lo + dx_uniform;
                        
                        const double weight_lo = ( x_integral_vector_hi - x[0])/dx_uniform;
                        const double weight_hi = (-x_integral_vector_lo + x[0])/dx_uniform;

                        double integral    = weight_lo*integral_vector[idx_integral_vector_lo] + weight_hi*integral_vector[idx_integral_vector_hi];
                        
                        double p_H         = p_i*exp(g/T_0*integral);
                        double rho_H       = p_H/(R_H*T_0);

                        double rho, p;
                        rho         = rho_H;
                        p           = p_H;
                
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