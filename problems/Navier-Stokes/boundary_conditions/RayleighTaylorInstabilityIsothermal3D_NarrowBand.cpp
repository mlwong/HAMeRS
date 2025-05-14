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
    if ((d_project_name != "3D smooth Rayleigh-Taylor instability") &&
       (d_project_name != "3D smooth multi-mode Rayleigh-Taylor instability"))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Can only initialize data for 'project_name' = "
            << "'3D smooth Rayleigh-Taylor instability' or"
            << "'3D smooth multi-mode Rayleigh Taylor instability'!\n"
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
    // Added 9/7/23
    TBOX_ASSERT(d_special_boundary_conditions_db != nullptr);
    TBOX_ASSERT(d_special_boundary_conditions_db->keyExists("gravity"));
    TBOX_ASSERT(d_special_boundary_conditions_db->keyExists("species_mass"));
    
    std::vector<double> gravity_vector = d_special_boundary_conditions_db->getDoubleVector("gravity");
    std::vector<double> W_vector = d_special_boundary_conditions_db->getDoubleVector("species_mass"); // molecular mass of mixing fluids
    
    /*
     * Determine the ghost cell width to fill.
     */
    hier::IntVector gcw_to_fill(tbox::Dimension(3));
    
    // If the ghost fill width is not used, it is set to the ghost cell width of the data.
    if (ghost_width_to_fill == -hier::IntVector::getOne(tbox::Dimension(3)))
    {
        gcw_to_fill = num_ghosts;
    }
    else
    {
        gcw_to_fill = hier::IntVector::min(
            num_ghosts,
            ghost_width_to_fill);
    }
    
    if (d_project_name == "3D smooth Rayleigh-Taylor instability")
    {
        // New code
        const double* const dx = patch_geom->getDx();
        const double* const patch_xlo = patch_geom->getXLower();
        
        HAMERS_SHARED_PTR<pdat::CellData<double> > partial_density = conservative_variables[0];
        HAMERS_SHARED_PTR<pdat::CellData<double> > momentum        = conservative_variables[1];
        HAMERS_SHARED_PTR<pdat::CellData<double> > total_energy    = conservative_variables[2];
        
        double* rho_Y_0 = partial_density->getPointer(0);
        double* rho_Y_1 = partial_density->getPointer(1);
        double* rho_u   = momentum->getPointer(0);
        double* rho_v   = momentum->getPointer(1);
        double* rho_w   = momentum->getPointer(2);
        double* E       = total_energy->getPointer(0);
        
        const double gamma = double(7)/double(5);
        const double g     = gravity_vector[0]; // gravity
        //const double g     = 90.0; // HARD-CODED for now!
        
        const double p_i = 100000.0; // interface pressure
        const double T_0 = 300.0; 
        
        const double W_1 = W_vector[0]; // molecular mass of heavier gas
        const double W_2 = W_vector[1]; // molecular mass of lighter gas
        
        const double R_u = 8.31446261815324; // universal gas constant
        const double R_1 = R_u/W_1;          // gas constant of heavier gas
        const double R_2 = R_u/W_2;
        
        const double lambda = 701.53278340668; // wavelength of single-mode perturbation
        const double eta_0  = 0.02*lambda;
        const double delta  = 0.02*lambda; // characteristic length of interface.
        const double shift  = 0.0; // location of interface.
            
        for (int codim = 1; codim <= d_dim.getValue(); codim++)
        {
            const std::vector<hier::BoundaryBox>& boundary_boxes = patch_geom->getCodimensionBoundaries(codim);
            
            if (!boundary_boxes.empty())
            {
                // Discretize the domain in x-direction for the approximated integral.
                const int integral_N_x = 10000;
                
                std::vector<double>& integral_vector = d_special_boundary_conditions_vector; // Reference to the class member vector that stores the numerical integral.
                
                if (integral_vector.empty())
                {
                    integral_vector.resize(integral_N_x + 3);
                    std::ifstream f_in;
                    std::string integral_filename = "integral.dat";
                    f_in.open(integral_filename, std::ios::in | std::ios::binary);
                    f_in.read((char*)&integral_vector[0], sizeof(double)*integral_vector.size());
                    f_in.close();
                }
                
                const double x_domain_lo = integral_vector[integral_N_x + 0];
                const double x_domain_hi = integral_vector[integral_N_x + 1];
                const double dx_uniform  = integral_vector[integral_N_x + 2];
                
                // Loop over the boundary boxes.
                for (int bi = 0; bi < static_cast<int>(boundary_boxes.size()); bi++)
                {
                    hier::Box fill_box(patch_geom->getBoundaryFillBox(
                        boundary_boxes[bi],
                        interior_box,
                        gcw_to_fill));
                    
                    hier::Index fill_box_lo_idx(fill_box.lower());
                    hier::Index fill_box_hi_idx(fill_box.upper());
                    
                    /*
                     * Offset the indices.
                     */
                    fill_box_lo_idx = fill_box_lo_idx - interior_box.lower();
                    fill_box_hi_idx = fill_box_hi_idx - interior_box.lower();
                    
                    for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                    {
                        for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                        {
                            for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                            {
                                const int idx_cell = (i + num_ghosts[0]) +
                                    (j + num_ghosts[1])*ghostcell_dims[0] +
                                    (k + num_ghosts[2]*ghostcell_dims[0]*ghostcell_dims[1]);
                                
                                // Compute the coordinates.
                                double x[3];
                                x[0] = patch_xlo[0] + (i + double(1)/double(2))*dx[0];
                                x[1] = patch_xlo[1] + (j + double(1)/double(2))*dx[1];
                                x[2] = patch_xlo[2] + (k + double(1)/double(2))*dx[2];
                                
                                const double eta = eta_0*cos(2.0*M_PI/lambda*x[1])*cos(2.0*M_PI/lambda*x[2]);
                                
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
    }
    if (d_project_name == "3D smooth multi-mode Rayleigh-Taylor instability")
    {
        // New code
        const double* const dx = patch_geom->getDx();
        const double* const patch_xlo = patch_geom->getXLower();
        
        HAMERS_SHARED_PTR<pdat::CellData<double> > partial_density = conservative_variables[0];
        HAMERS_SHARED_PTR<pdat::CellData<double> > momentum        = conservative_variables[1];
        HAMERS_SHARED_PTR<pdat::CellData<double> > total_energy    = conservative_variables[2];
        
        double* rho_Y_0 = partial_density->getPointer(0);
        double* rho_Y_1 = partial_density->getPointer(1);
        double* rho_u   = momentum->getPointer(0);
        double* rho_v   = momentum->getPointer(1);
        double* rho_w   = momentum->getPointer(2);
        double* E       = total_energy->getPointer(0);
        
        const double gamma = double(7)/double(5);
        const double g     = gravity_vector[0]; // gravity
        
        const double p_i = 100000.0; // interface pressure
        const double T_0 = 300.0; 
        
        const double W_1 = W_vector[0]; // molecular mass of heavier gas
        const double W_2 = W_vector[1]; // molecular mass of lighter gas

        
        const double R_u = 8.31446261815324; // universal gas constant
        const double R_1 = R_u/W_1;          // gas constant of heavier gas
        const double R_2 = R_u/W_2;
        
        const double lambda = 701.53278340668; // characteristic wave-length
        const double width   = 4.0*lambda;  // domain size in y direction
        const int k_max = 24;
        const int k_min = 8;
        const int k_base = 16;
        // const int dominant_waven    = (k_max + k_min) / 2.0 ;
        const double lambda_0       = width / (k_base);
        const double delta   = lambda_0 * 0.08; // Characteristic length of interface
        const double eta_0   = lambda_0 * 0.04;
        const double shift  = 0.0; // location of interface.

        // const double std_desired = (width/12)*0.04;


        const int delta_ky = 1;
        const int delta_kz = 1;
        const double k_0 = 2.0*M_PI/width;
        const double epsilon = 1.0e-15;

        std::vector<double> a_mn((k_max+1) * (k_max+1)); // Create vector for random number a
        // Read random number list a
        std::ifstream f_rand_a;
        std::string random_a_filename = "random_a.dat";
        f_rand_a.open(random_a_filename, std::ios::in | std::ios::binary);
        f_rand_a.read((char*)&a_mn[0], sizeof(double)*a_mn.size());
        f_rand_a.close();

        std::vector<double> b_mn((k_max+1) * (k_max+1)); // Create vector for random number b
        // Read random number list b
        std::ifstream f_rand_b;
        std::string random_b_filename = "random_b.dat";
        f_rand_b.open(random_b_filename, std::ios::in | std::ios::binary);
        f_rand_b.read((char*)&b_mn[0], sizeof(double)*b_mn.size());
        f_rand_b.close();

        std::vector<double> c_mn((k_max+1) * (k_max+1)); // Create vector for random number c
        // Read random number list c
        std::ifstream f_rand_c;
        std::string random_c_filename = "random_c.dat";
        f_rand_c.open(random_c_filename, std::ios::in | std::ios::binary);
        f_rand_c.read((char*)&c_mn[0], sizeof(double)*c_mn.size());
        f_rand_c.close();
    
        std::vector<double> d_mn((k_max+1) * (k_max+1)); // Create vector for random number d
        // Read random number list d
        std::ifstream f_rand_d;
        std::string random_d_filename = "random_d.dat";
        f_rand_d.open(random_d_filename, std::ios::in | std::ios::binary);
        f_rand_d.read((char*)&d_mn[0], sizeof(double)*d_mn.size());
        f_rand_d.close();

       
        for (int codim = 1; codim <= d_dim.getValue(); codim++)
        {
            const std::vector<hier::BoundaryBox>& boundary_boxes = patch_geom->getCodimensionBoundaries(codim);
            
            if (!boundary_boxes.empty())
            {
                // Discretize the domain in x-direction for the approximated integral.
                const int integral_N_x = 10000;
                
                std::vector<double>& integral_vector = d_special_boundary_conditions_vector; // Reference to the class member vector that stores the numerical integral.
                
                if (integral_vector.empty())
                {
                    integral_vector.resize(integral_N_x + 3);
                    std::ifstream f_in;
                    std::string integral_filename = "integral.dat";
                    f_in.open(integral_filename, std::ios::in | std::ios::binary);
                    f_in.read((char*)&integral_vector[0], sizeof(double)*integral_vector.size());
                    f_in.close();
                }
                
                const double x_domain_lo = integral_vector[integral_N_x + 0];
                const double x_domain_hi = integral_vector[integral_N_x + 1];
                const double dx_uniform  = integral_vector[integral_N_x + 2];
                
                // Loop over the boundary boxes.
                for (int bi = 0; bi < static_cast<int>(boundary_boxes.size()); bi++)
                {
                    hier::Box fill_box(patch_geom->getBoundaryFillBox(
                        boundary_boxes[bi],
                        interior_box,
                        gcw_to_fill));
                    
                    hier::Index fill_box_lo_idx(fill_box.lower());
                    hier::Index fill_box_hi_idx(fill_box.upper());
                    
                    /*
                     * Offset the indices.
                     */
                    fill_box_lo_idx = fill_box_lo_idx - interior_box.lower();
                    fill_box_hi_idx = fill_box_hi_idx - interior_box.lower();
                    
                    for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                    {
                        for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                        {
                            for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                            {
                                const int idx_cell = (i + num_ghosts[0]) +
                                    (j + num_ghosts[1])*ghostcell_dims[0] +
                                    (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                
                                // Compute the coordinates.
                                double x[3];
                                x[0] = patch_xlo[0] + (i + double(1)/double(2))*dx[0];
                                if ( x[0] < x_domain_lo || x[0] > x_domain_hi )
                                {
                                    std::cout<<"SHORT DOMAIN INTEGRAL" <<x[0] <<"; " << x_domain_hi <<"; " << x_domain_lo <<"; " << std::endl;
                                }
                                x[1] = patch_xlo[1] + (j + double(1)/double(2))*dx[1];
                                x[2] = patch_xlo[2] + (k + double(1)/double(2))*dx[2];
                                
                                double eta = 0.0; 
                                int idx_random = 0;       
                                double k_amp = 0.0;
                                double ratio = 0.0;


                                for (int m = 0; m <= k_max; m++)
                                {
                                    for (int n = 0; n <= k_max; n++)
                                    {
                                        if (m*m + n*n > k_max*k_max || m*m + n*n < k_min*k_min)
                                        {
                                            k_amp = 0.0;
                                            ratio = 1.0;                                  
                                        }
                                        else
                                        {
                                            k_amp = 1.0;
                                        }

                                        int idx_random = m + (k_max+1)*n;

                                        ratio = sqrt(1.0 / 4.0 * (a_mn[idx_random]*a_mn[idx_random] + b_mn[idx_random]*b_mn[idx_random] + \
                                        c_mn[idx_random]*c_mn[idx_random] + d_mn[idx_random] * d_mn[idx_random]) * \
                                        2.0 * M_PI * sqrt(m * m + n * n) / (eta_0 * delta_ky * delta_kz));
                                        
                                        if (ratio <= epsilon) {
                                            ratio = 1.0;
                                        }

                                        eta += (k_amp / ratio)*(a_mn[idx_random]*cos(k_0 * m * x[1])*cos(k_0 * n * x[2]) + \
                                                                b_mn[idx_random]*cos(k_0 * m * x[1])*sin(k_0 * n * x[2]) + \
                                                                c_mn[idx_random]*sin(k_0 * m * x[1])*cos(k_0 * n * x[2]) + \
                                                                d_mn[idx_random]*sin(k_0 * m * x[1])*sin(k_0 * n * x[2])); 
                                    }
                                }




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
                                //
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
    }
    if (d_project_name == "TEST 3D smooth multi-mode Rayleigh-Taylor instability")
    {
        // New code
        const double* const dx = patch_geom->getDx();
        const double* const patch_xlo = patch_geom->getXLower();
        
        HAMERS_SHARED_PTR<pdat::CellData<double> > partial_density = conservative_variables[0];
        HAMERS_SHARED_PTR<pdat::CellData<double> > momentum        = conservative_variables[1];
        HAMERS_SHARED_PTR<pdat::CellData<double> > total_energy    = conservative_variables[2];
        
        double* rho_Y_0 = partial_density->getPointer(0);
        double* rho_Y_1 = partial_density->getPointer(1);
        double* rho_u   = momentum->getPointer(0);
        double* rho_v   = momentum->getPointer(1);
        double* rho_w   = momentum->getPointer(2);
        double* E       = total_energy->getPointer(0);
        
        const double gamma = double(7)/double(5);
        const double g     = gravity_vector[0]; // gravity
        
        const double p_i = 100000.0; // interface pressure
        const double T_0 = 300.0; 
        
        const double W_1 = W_vector[0]; // molecular mass of heavier gas
        const double W_2 = W_vector[1]; // molecular mass of lighter gas

        
        const double R_u = 8.31446261815324; // universal gas constant
        const double R_1 = R_u/W_1;          // gas constant of heavier gas
        const double R_2 = R_u/W_2;
        
        const double lambda = 701.53278340668 / 4.0; // characteristic wave-length
        // const double eta_0  = 0.08*lambda;
        const double shift  = 0.0; // location of interface.
        const double width   = 12.0*lambda;  // domain size in y direction

        for (int codim = 1; codim <= d_dim.getValue(); codim++)
        {
            const std::vector<hier::BoundaryBox>& boundary_boxes = patch_geom->getCodimensionBoundaries(codim);
            
            if (!boundary_boxes.empty())
            {
                for (int bi = 0; bi < static_cast<int>(boundary_boxes.size()); bi++)
                {
                    hier::Box fill_box(patch_geom->getBoundaryFillBox(
                        boundary_boxes[bi],
                        interior_box,
                        gcw_to_fill));
                    
                    hier::Index fill_box_lo_idx(fill_box.lower());
                    hier::Index fill_box_hi_idx(fill_box.upper());
                    
                    /*
                     * Offset the indices.
                     */
                    fill_box_lo_idx = fill_box_lo_idx - interior_box.lower();
                    fill_box_hi_idx = fill_box_hi_idx - interior_box.lower();
                    
                    for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                    {
                        for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                        {
                            for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                            {
                                const int idx_cell = (i + num_ghosts[0]) +
                                    (j + num_ghosts[1])*ghostcell_dims[0] +
                                    (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                
                                // Compute the coordinates.
                                double x[3];
                                x[0] = patch_xlo[0] + (i + double(1)/double(2))*dx[0];
                                x[1] = patch_xlo[1] + (j + double(1)/double(2))*dx[1];
                                x[2] = patch_xlo[2] + (k + double(1)/double(2))*dx[2];
                                
                                
                                double rho, p, X_2_H;
                                
                                rho = 1.0;
                                p   = p_i;
                                X_2_H = 0.5;

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
    }    
}
