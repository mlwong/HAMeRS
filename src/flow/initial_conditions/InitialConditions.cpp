#include "flow/initial_conditions/InitialConditions.hpp"

/*
 * Set the data on the patch interior to some initial values,
 * depending on the flow problems and flow models.
 */
void
InitialConditions::initializeDataOnPatch(
    hier::Patch& patch,
    const double data_time,
    const bool initial_time,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
    NULL_USE(data_time);
    
    if (initial_time)
    {
        const double* const domain_xlo = d_grid_geometry->getXLower();
        const double* const domain_xhi = d_grid_geometry->getXUpper();
        
        const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
            BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                patch.getPatchGeometry()));
        
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(patch_geom);
#endif
        
        const double* const dx = patch_geom->getDx();
        const double* const patch_xlo = patch_geom->getXLower();
        
        // Get the dimensions of box that covers the interior of Patch.
        hier::Box patch_box = patch.getBox();
        const hier::IntVector patch_dims = patch_box.numberCells();
        
        // Get the ratio to level zero index space.
        const hier::IntVector ratio_to_level_zero = patch_geom->getRatio();
        
        switch (d_flow_model_label)
        {
            case SINGLE_SPECIES:
            {
                boost::shared_ptr<pdat::CellData<double> > density(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_density, data_context)));
                
                boost::shared_ptr<pdat::CellData<double> > momentum(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_momentum, data_context)));
                
                boost::shared_ptr<pdat::CellData<double> > total_energy(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_total_energy, data_context)));
                
#ifdef DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(density);
                TBOX_ASSERT(momentum);
                TBOX_ASSERT(total_energy);
#endif
                
                if (d_dim == tbox::Dimension(1))
                {
                    // NOT YET IMPLEMENTED.
                    
                    TBOX_ERROR(d_object_name
                        << ": "
                        << "Cannot initialize data for unknown 1D problem"
                        << " for single-species flow with name = '"
                        << d_project_name
                        << "'."
                        << std::endl);
                }
                else if (d_dim == tbox::Dimension(2))
                {
                    if (d_project_name == "2D advection of density wave")
                    {
                        // Initialize data for a 2D density wave advection problem.
                        double* rho   = density->getPointer(0);
                        double* rho_u = momentum->getPointer(0);
                        double* rho_v = momentum->getPointer(1);
                        double* E     = total_energy->getPointer(0);
                        
                        const double x_a = 1.0/3*(domain_xlo[0] + domain_xhi[0]);
                        const double x_b = 2.0/3*(domain_xlo[0] + domain_xhi[0]);
                        
                        const double y_a = 1.0/3*(domain_xlo[1] + domain_xhi[1]);
                        const double y_b = 2.0/3*(domain_xlo[1] + domain_xhi[1]);
                        
                        double gamma = d_equation_of_state->getSpeciesThermodynamicProperty(
                            "gamma",
                            0);
                        
                        // Initial conditions inside the square.
                        double rho_i = 10.0;
                        double u_i   = 1.0;
                        double v_i   = 1.0;
                        double p_i   = 1.0;
                        
                        // Initial conditions outside the square.
                        double rho_o = 1.0;
                        double u_o   = 1.0;
                        double v_o   = 1.0;
                        double p_o   = 1.0;
                        
                        for (int j = 0; j < patch_dims[1]; j++)
                        {
                            for (int i = 0; i < patch_dims[0]; i++)
                            {
                                // Compute index into linear data array.
                                int idx_cell = i + j*patch_dims[0];
                                
                                // Compute the coordinates.
                                double x[2];
                                x[0] = patch_xlo[0] + (i + 0.5)*dx[0];
                                x[1] = patch_xlo[1] + (j + 0.5)*dx[1];
                                
                                if ((x[0] >= x_a) &&
                                    (x[0] <= x_b) &&
                                    (x[1] >= y_a) &&
                                    (x[1] <= y_b))
                                {
                                    rho[idx_cell]   = rho_i;
                                    rho_u[idx_cell] = rho_i*u_i;
                                    rho_v[idx_cell] = rho_i*v_i;
                                    E[idx_cell]     = p_i/(gamma - 1.0) + 0.5*rho_i*(u_i*u_i +
                                        v_i*v_i);
                                }
                                else
                                {
                                    rho[idx_cell]   = rho_o;
                                    rho_u[idx_cell] = rho_o*u_o;
                                    rho_v[idx_cell] = rho_o*v_o;
                                    E[idx_cell]     = p_o/(gamma - 1.0) + 0.5*rho_o*(u_o*u_o +
                                        v_o*v_o);
                                }
                            }
                        }
                    }
                    else if (d_project_name == "2D steady wedge flow")
                    {
                        // Initialize data for a 2D steady wedge flow problem.
                        double* rho   = density->getPointer(0);
                        double* rho_u = momentum->getPointer(0);
                        double* rho_v = momentum->getPointer(1);
                        double* E     = total_energy->getPointer(0);
                        
                        const double gamma = 1.4;
                        const double R     = 287.058;
                        
                        const double p_inf = 1e5;
                        const double T_inf = 300;
                        const double M_inf = 2;
                        const double theta = 10.0/180*M_PI;
                        
                        const double U_inf   = M_inf * sqrt(gamma*R*T_inf);
                        const double rho_inf = p_inf/R/T_inf;
                        
                        const double u_inf = U_inf*cos(theta);
                        const double v_inf = -U_inf*sin(theta);
                        
                        for (int j = 0; j < patch_dims[1]; j++)
                        {
                            for (int i = 0; i < patch_dims[0]; i++)
                            {
                                // Compute index into linear data array.
                                int idx_cell = i + j*patch_dims[0];
                                
                                rho[idx_cell]   = rho_inf;
                                rho_u[idx_cell] = rho_inf*u_inf;
                                rho_v[idx_cell] = rho_inf*v_inf;
                                E[idx_cell]     = p_inf/(gamma - 1.0) + 0.5*rho_inf*(
                                    u_inf*u_inf + v_inf*v_inf);
                            }
                        }
                    }
                    else if (d_project_name == "2D double-Mach reflection")
                    {
                        // Initialize data for a 2D double-Mach reflection problem.
                        double* rho   = density->getPointer(0);
                        double* rho_u = momentum->getPointer(0);
                        double* rho_v = momentum->getPointer(1);
                        double* E     = total_energy->getPointer(0);
                        
                        const double x_0 = 1.0/6.0;
                        
                        const double gamma = 1.4;
                        
                        const double rho_post_shock = 8.0;
                        const double u_post_shock   = 8.25*cos(M_PI/6.0);
                        const double v_post_shock   = -8.25*sin(M_PI/6.0);
                        const double p_post_shock   = 116.5;
                        
                        const double rho_pre_shock  = 1.4;
                        const double u_pre_shock    = 0.0;
                        const double v_pre_shock    = 0.0;
                        const double p_pre_shock    = 1.0;
                        
                        const double rho_u_post_shock = rho_post_shock*u_post_shock;
                        const double rho_v_post_shock = rho_post_shock*v_post_shock;
                        
                        const double rho_u_pre_shock  = rho_pre_shock*u_pre_shock;
                        const double rho_v_pre_shock  = rho_pre_shock*v_pre_shock;
                        
                        const double E_pre_shock = p_pre_shock/(gamma - 1.0) +
                            0.5*rho_pre_shock*(u_pre_shock*u_pre_shock + v_pre_shock*v_pre_shock);
                        
                        const double E_post_shock = p_post_shock/(gamma - 1.0) +
                            0.5*rho_post_shock*(u_post_shock*u_post_shock + v_post_shock*v_post_shock);
                            
                        for (int j = 0; j < patch_dims[1]; j++)
                        {
                            for (int i = 0; i < patch_dims[0]; i++)
                            {
                                // Compute index into linear data array.
                                int idx_cell = i + j*patch_dims[0];
                                
                                // Compute the coordinates.
                                double x[2];
                                x[0] = patch_xlo[0] + (i + 0.5)*dx[0];
                                x[1] = patch_xlo[1] + (j + 0.5)*dx[1];
                                
                                if (x[0] < x_0 + x[1]*sqrt(1.0/3.0))
                                {
                                    rho[idx_cell] = rho_post_shock;
                                    rho_u[idx_cell] = rho_u_post_shock;
                                    rho_v[idx_cell] = rho_v_post_shock;
                                    E[idx_cell] = E_post_shock;
                                }
                                else
                                {
                                    rho[idx_cell] = rho_pre_shock;
                                    rho_u[idx_cell] = rho_u_pre_shock;
                                    rho_v[idx_cell] = rho_v_pre_shock;
                                    E[idx_cell] = E_pre_shock;
                                }
                            }
                        }
                    }
                    else if (d_project_name == "2D shock-vortex interaction")
                    {
                        // Initialize data for a 2D shock-vortex interaction problem.
                        double* rho   = density->getPointer(0);
                        double* rho_u = momentum->getPointer(0);
                        double* rho_v = momentum->getPointer(1);
                        double* E     = total_energy->getPointer(0);
                        
                        double gamma = d_equation_of_state->getSpeciesThermodynamicProperty(
                            "gamma",
                            0);
                        
                        // Vortex strength.
                        const double M_v = 1.0;
                        
                        // Vortex radius.
                        const double R = 1.0;
                        
                        // Vortex center.
                        double x_v[2];
                        x_v[0] = 4.0;
                        x_v[1] = 0.0;
                        
                        // Post-shock condition.
                        const double rho_post = 1.34161490;
                        const double p_post = 1.51333333/gamma;
                        const double u_post = -0.89444445;
                        const double v_post = 0.0;
                        
                        // Pre-shock condition.
                        const double rho_pre = 1.0;
                        const double p_pre = 1.0/gamma;
                        const double u_pre = -1.2;
                        const double v_pre = 0.0;
                        
                        for (int j = 0; j < patch_dims[1]; j++)
                        {
                            for (int i = 0; i < patch_dims[0]; i++)
                            {
                                // Compute index into linear data array.
                                int idx_cell = i + j*patch_dims[0];
                                
                                // Compute the coordinates.
                                double x[2];
                                x[0] = patch_xlo[0] + (i + 0.5)*dx[0];
                                x[1] = patch_xlo[1] + (j + 0.5)*dx[1];
                                
                                if (x[0] < 0)
                                {
                                    rho[idx_cell]     = rho_post;
                                    rho_u[idx_cell]   = rho_post*u_post;
                                    rho_v[idx_cell]   = rho_post*v_post;
                                    E[idx_cell]       = p_post/(gamma - 1.0) +
                                        0.5*rho_post*(u_post*u_post + v_post*v_post);
                                }
                                else
                                {
                                    double r = sqrt(pow(x[0] - x_v[0], 2) + pow(x[1] - x_v[1], 2));
                                    
                                    if (r > 4)
                                    {
                                        rho[idx_cell]   = rho_pre;
                                        rho_u[idx_cell] = rho_pre*u_pre;
                                        rho_v[idx_cell] = rho_pre*v_pre;
                                        E[idx_cell]     = p_pre/(gamma - 1.0) +
                                            0.5*rho_pre*(u_pre*u_pre + v_pre*v_pre);
                                    }
                                    else
                                    {
                                        double p = 1.0/(gamma)*pow((1.0 - 0.5*(gamma - 1.0)*M_v*M_v*exp(1 - pow(r/R, 2))),
                                                                  gamma/(gamma - 1.0));
                                        double u = u_pre - M_v*exp(0.5*(1 - pow(r/R, 2)))*(x[1] - x_v[1]);
                                        double v = v_pre + M_v*exp(0.5*(1 - pow(r/R, 2)))*(x[0] - x_v[0]);
                                        
                                        
                                        rho[idx_cell] = pow((1.0 - 0.5*(gamma - 1.0)*M_v*M_v*exp(1 - pow(r/R, 2))),
                                                             1.0/(gamma - 1.0));
                                        rho_u[idx_cell] = rho[idx_cell]*u;
                                        rho_v[idx_cell] = rho[idx_cell]*v;
                                        E[idx_cell] = p/(gamma - 1.0) +
                                            0.5*rho[idx_cell]*(u*u + v*v);
                                    }
                                }
                            }
                        }
                    }
                    /*
                    else if (d_project_name == "2D Kelvin-Helmholtz instability")
                    {
                        // Radice, David, and Luciano Rezzolla.
                        // "THC: a new high-order finite-difference high-resolution shock-capturing code
                        // for special-relativistic hydrodynamics." Astronomy & Astrophysics 547 (2012): A26.
                        double* rho   = density->getPointer(0);
                        double* rho_u = momentum->getPointer(0);
                        double* rho_v = momentum->getPointer(1);
                        double* E     = total_energy->getPointer(0);
                        
                        double gamma = d_equation_of_state->getSpeciesThermodynamicProperty(
                            "gamma",
                            0);
                        
                        // Some parameters for the problem.
                        const double a       = 0.04; // Characteristic size of the shear layer.
                        const double v_shear = 0.5;  // Correspnding to a relative Lorentz factor of W = 2.29
                        const double A_0     = 0.1;  // Perturbation amplitude.
                        const double sigma   = 0.1;  // Characteristic lengthscale.
                        
                        
                        const double rho_0 = 0.505;
                        const double rho_1 = 0.495;
                        
                        for (int j = 0; j < patch_dims[1]; j++)
                        {
                            for (int i = 0; i < patch_dims[0]; i++)
                            {
                                // Compute index into linear data array.
                                int idx_cell = i + j*patch_dims[0];
                                
                                // Compute the coordinates.
                                double x[2];
                                x[0] = patch_xlo[0] + (i + 0.5)*dx[0];
                                x[1] = patch_xlo[1] + (j + 0.5)*dx[1];
                                
                                double u, v;
                                
                                if (x[1] > 0)
                                {
                                    rho[idx_cell] = rho_0 + rho_1*tanh((x[1] - 0.5)/a);
                                    u = v_shear*tanh((x[1] - 0.5)/a);
                                    v = A_0*v_shear*sin(2.0*M_PI*x[0])*exp(-pow(x[1] - 0.5, 2)/sigma);
                                }
                                else
                                {
                                    rho[idx_cell] = rho_0 - rho_1*tanh((x[1] + 0.5)/a);
                                    u = -v_shear*tanh((x[1] + 0.5)/a);
                                    v = -A_0*v_shear*sin(2.0*M_PI*x[0])*exp(-pow(x[1] + 0.5, 2)/sigma);
                                }
                                
                                double p = 1.0;
                                
                                rho_u[idx_cell] = rho[idx_cell]*u;
                                rho_v[idx_cell] = rho[idx_cell]*v;
                                E[idx_cell]     = p/(gamma - 1.0) + 0.5*rho[idx_cell]*(u*u + v*v);
                            }
                        }
                    }
                    */
                    else if (d_project_name == "2D Kelvin-Helmholtz instability")
                    {
                        // Initialize data for a 2D Kelvin-Helmohltz instability problem.
                        // (McNally, Colin P., Wladimir Lyra, and Jean-Claude Passy.
                        // "A well-posed Kelvin-Helmholtz instability test and comparison."
                        // The Astrophysical Journal Supplement Series 201.2 (2012): 18.)
                        double* rho   = density->getPointer(0);
                        double* rho_u = momentum->getPointer(0);
                        double* rho_v = momentum->getPointer(1);
                        double* E     = total_energy->getPointer(0);
                        
                        double gamma = d_equation_of_state->getSpeciesThermodynamicProperty(
                            "gamma",
                            0);
                        
                        // Some parameters for the problem.
                        const double rho_1 = 1.0;
                        const double rho_2 = 2.0;
                        const double rho_m = 0.5*(rho_1 - rho_2);
                        const double u_1   = 0.5;
                        const double u_2   = -0.5;
                        const double u_m   = 0.5*(u_1 - u_2);
                        const double L     = 0.025;
                        
                        for (int j = 0; j < patch_dims[1]; j++)
                        {
                            for (int i = 0; i < patch_dims[0]; i++)
                            {
                                // Compute index into linear data array.
                                int idx_cell = i + j*patch_dims[0];
                                
                                // Compute the coordinates.
                                double x[2];
                                x[0] = patch_xlo[0] + (i + 0.5)*dx[0];
                                x[1] = patch_xlo[1] + (j + 0.5)*dx[1];
                                
                                double u   = 0.0;
                                double v   = 0.01*sin(4*M_PI*x[0]);
                                double p   = 2.5;
                                
                                if (x[1] < 0.25)
                                {
                                    rho[idx_cell] = rho_1 - rho_m*exp((x[1] - 0.25)/L);
                                    u = u_1 - u_m*exp((x[1] - 0.25)/L);
                                    
                                }
                                else if (x[1] < 0.5)
                                {
                                    rho[idx_cell] = rho_2 + rho_m*exp((-x[1] + 0.25)/L);
                                    u   = u_2 + u_m*exp((-x[1] + 0.25)/L);
                                }
                                else if (x[1] < 0.75)
                                {
                                    rho[idx_cell] = rho_2 + rho_m*exp(-(0.75 - x[1])/L);
                                    u   = u_2 + u_m*exp(-(0.75 - x[1])/L);
                                }
                                else
                                {
                                    rho[idx_cell] = rho_1 - rho_m*exp(-(x[1] - 0.75)/L);
                                    u   = u_1 - u_m*exp(-(x[1] - 0.75)/L);
                                }
                                
                                rho_u[idx_cell] = rho[idx_cell]*u;
                                rho_v[idx_cell] = rho[idx_cell]*v;
                                E[idx_cell]     = p/(gamma - 1.0) +
                                    0.5*rho[idx_cell]*(u*u + v*v);
                            }
                        }
                    }
                    else
                    {
                        TBOX_ERROR(d_object_name
                            << ": "
                            << "Cannot initialize data for unknown 2D problem"
                            << " for single-species flow with name = '"
                            << d_project_name
                            << "'."
                            << std::endl);
                    }
                }
                else if (d_dim == tbox::Dimension(3))
                {
                    if (d_project_name == "3D advection of density wave")
                    {
                        // Initialize data for a 3D density wave advection problem.
                        double* rho   = density->getPointer(0);
                        double* rho_u = momentum->getPointer(0);
                        double* rho_v = momentum->getPointer(1);
                        double* rho_w = momentum->getPointer(2);
                        double* E     = total_energy->getPointer(0);
                        
                        const double x_a = 1.0/3*(domain_xlo[0] + domain_xhi[0]);
                        const double x_b = 2.0/3*(domain_xlo[0] + domain_xhi[0]);
                        
                        const double y_a = 1.0/3*(domain_xlo[1] + domain_xhi[1]);
                        const double y_b = 2.0/3*(domain_xlo[1] + domain_xhi[1]);
                        
                        const double z_a = 1.0/3*(domain_xlo[2] + domain_xhi[2]);
                        const double z_b = 2.0/3*(domain_xlo[2] + domain_xhi[2]);
                        
                        double gamma = d_equation_of_state->getSpeciesThermodynamicProperty(
                            "gamma",
                            0);
                        
                        // Initial conditions inside the cube.
                        double rho_i = 10.0;
                        double u_i   = 1.0;
                        double v_i   = 1.0;
                        double w_i   = 1.0;
                        double p_i   = 1.0;
                        
                        // Initial conditions outside the cube.
                        double rho_o = 1.0;
                        double u_o   = 1.0;
                        double v_o   = 1.0;
                        double w_o   = 1.0;
                        double p_o   = 1.0;
                        
                        for (int k = 0; k < patch_dims[2]; k++)
                        {
                            for (int j = 0; j < patch_dims[1]; j++)
                            {
                                for (int i = 0; i < patch_dims[0]; i++)
                                {
                                    // Compute index into linear data array.
                                    int idx_cell = i +
                                        j*patch_dims[0] +
                                        k*patch_dims[0]*patch_dims[1];
                                    
                                    // Compute the coordinates.
                                    double x[3];
                                    x[0] = patch_xlo[0] + (i + 0.5)*dx[0];
                                    x[1] = patch_xlo[1] + (j + 0.5)*dx[1];
                                    x[2] = patch_xlo[2] + (k + 0.5)*dx[2];
                                    
                                    if ((x[0] >= x_a) &&
                                        (x[0] <= x_b) &&
                                        (x[1] >= y_a) &&
                                        (x[1] <= y_b) &&
                                        (x[2] >= z_a) &&
                                        (x[2] <= z_b))
                                    {
                                        rho[idx_cell]   = rho_i;
                                        rho_u[idx_cell] = rho_i*u_i;
                                        rho_v[idx_cell] = rho_i*v_i;
                                        rho_w[idx_cell] = rho_i*w_i;
                                        E[idx_cell]     = p_i/(gamma - 1.0) + 0.5*rho_i*(u_i*u_i +
                                            v_i*v_i + w_i*w_i);
                                    }
                                    else
                                    {
                                        rho[idx_cell]   = rho_o;
                                        rho_u[idx_cell] = rho_o*u_o;
                                        rho_v[idx_cell] = rho_o*v_o;
                                        rho_w[idx_cell] = rho_o*w_o;
                                        E[idx_cell]     = p_o/(gamma - 1.0) + 0.5*rho_o*(u_o*u_o +
                                            v_o*v_o + w_i*w_i);
                                    }
                                }
                            }
                        }
                    }
                    else if (d_project_name == "3D Taylor-Green vortex")
                    {
                        // Initialize data for a 3D inviscid Taylor-Green problem.
                        double* rho   = density->getPointer(0);
                        double* rho_u = momentum->getPointer(0);
                        double* rho_v = momentum->getPointer(1);
                        double* rho_w = momentum->getPointer(2);
                        double* E     = total_energy->getPointer(0);
                        
                        double gamma = d_equation_of_state->getSpeciesThermodynamicProperty(
                            "gamma",
                            0);
                        
                        for (int k = 0; k < patch_dims[2]; k++)
                        {
                            for (int j = 0; j < patch_dims[1]; j++)
                            {
                                for (int i = 0; i < patch_dims[0]; i++)
                                {
                                    // Compute index into linear data array.
                                    int idx_cell = i +
                                        j*patch_dims[0] +
                                        k*patch_dims[0]*patch_dims[1];
                                    
                                    // Compute the coordinates.
                                    double x[3];
                                    x[0] = patch_xlo[0] + (i + 0.5)*dx[0];
                                    x[1] = patch_xlo[1] + (j + 0.5)*dx[1];
                                    x[2] = patch_xlo[2] + (k + 0.5)*dx[2];
                                    
                                    // Compute density, velocities and pressure.
                                    const double rho_i = 1.0;
                                    const double u_i = sin(x[0])*cos(x[1])*cos(x[2]);
                                    const double v_i = -cos(x[0])*sin(x[1])*cos(x[2]);
                                    const double w_i = 0.0;
                                    const double p_i = 100.0 +
                                        ((cos(2*x[2]) + 2.0)*(cos(2*x[0]) + cos(2*x[1])) - 2.0)/
                                            16.0;
                                    
                                    rho[idx_cell]   = rho_i;
                                    rho_u[idx_cell] = rho_i*u_i;
                                    rho_v[idx_cell] = rho_i*v_i;
                                    rho_w[idx_cell] = rho_i*w_i;
                                    E[idx_cell]     = p_i/(gamma - 1.0) + 0.5*rho_i*(u_i*u_i +
                                        v_i*v_i + w_i*w_i);
                                }
                            }
                        }
                    }
                    else
                    {
                        TBOX_ERROR(d_object_name
                            << ": "
                            << "Cannot initialize data for unknown 3D problem"
                            << " for single-species flow with name = '"
                            << d_project_name
                            << "'."
                            << std::endl);
                    }
                }
                
                break;
            }
            case FOUR_EQN_CONSERVATIVE:
            {
                boost::shared_ptr<pdat::CellData<double> > partial_density(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_partial_density, data_context)));
                
                boost::shared_ptr<pdat::CellData<double> > momentum(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_momentum, data_context)));
                
                boost::shared_ptr<pdat::CellData<double> > total_energy(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_total_energy, data_context)));
                
#ifdef DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(partial_density);
                TBOX_ASSERT(momentum);
                TBOX_ASSERT(total_energy);
#endif
                
                if (d_dim == tbox::Dimension(1))
                {
                    // NOT YET IMPLEMENTED.
                    
                    TBOX_ERROR(d_object_name
                        << ": "
                        << "Cannot initialize data for unknown 1D problem"
                        << " for muti-species flow with conservative four-equation model with name = '"
                        << d_project_name
                        << "'."
                        << std::endl);
                }
                else if (d_dim == tbox::Dimension(2))
                {
                    if (d_project_name == "2D shock-bubble interaction with AMR")
                    {
                        if (d_num_species != 2)
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "Please provide only two-species for the problem"
                                << " '2D shock-bubble interaction with AMR'."
                                << std::endl);
                        }
                        
                        // Characteristic length of the problem.
                        const double D = 1.0;
                        
                        // Compute the characteristic length of the initial interface thickness.
                        const double C_epsilon = 6.0;
                        const double epsilon_i = C_epsilon*0.0025; // epsilon_i for smoothing interface
                        
                        double* rho_Y_1   = partial_density->getPointer(0);
                        double* rho_Y_2   = partial_density->getPointer(1);
                        double* rho_u     = momentum->getPointer(0);
                        double* rho_v     = momentum->getPointer(1);
                        double* E         = total_energy->getPointer(0);
                        
                        // species 0: He
                        // species 1: air
                        
                        // He, pre-shock condition.
                        const double rho_He = 0.1819;
                        const double u_He   = 0.0;
                        const double v_He   = 0.0;
                        const double p_He   = 1.0/1.4;
                        
                        // air, pre-shock condition.
                        const double rho_pre = 1.0;
                        const double u_pre   = 0.0;
                        const double v_pre   = 0.0;
                        const double p_pre   = 1.0/1.4;
                        
                        // air, post-shock condition.
                        const double rho_post = 1.3764;
                        const double u_post   = -0.3336;
                        const double v_post   = 0.0;
                        const double p_post   = 1.5698/1.4;
                        
                        for (int j = 0; j < patch_dims[1]; j++)
                        {
                            for (int i = 0; i < patch_dims[0]; i++)
                            {
                                // Compute index into linear data array.
                                int idx_cell = i + j*patch_dims[0];
                                
                                // Compute the coordinates.
                                double x[2];
                                x[0] = patch_xlo[0] + (i + 0.5)*dx[0];
                                x[1] = patch_xlo[1] + (j + 0.5)*dx[1];
                                
                                if (x[0] > 4.5*D)
                                {
                                    const double Y_1 = 0.0;
                                    const double Y_2 = 1.0;
                                    std::vector<const double*> Y_ptr;
                                    Y_ptr.push_back(&Y_1);
                                    Y_ptr.push_back(&Y_2);
                                    
                                    const double gamma = d_equation_of_state->
                                        getMixtureThermodynamicPropertyWithMassFraction(
                                            "gamma",
                                            Y_ptr);
                                    
                                    rho_Y_1[idx_cell] = 0.0;
                                    rho_Y_2[idx_cell] = rho_post;
                                    rho_u[idx_cell]   = rho_post*u_post;
                                    rho_v[idx_cell]   = rho_post*v_post;
                                    E[idx_cell]       = p_post/(gamma - 1.0) +
                                        0.5*rho_post*(u_post*u_post + v_post*v_post);
                                }
                                else
                                {
                                    // Compute the distance from the initial material interface.
                                    const double dR = sqrt(pow(x[0] - 3.5, 2) + x[1]*x[1]) - 0.5*D;
                                    
                                    const double f_sm = 0.5*(1.0 + erf(dR/epsilon_i));
                                    
                                    // Smooth the primitive quantity.
                                    const double rho_Y_1_i = rho_He*(1.0 - f_sm);
                                    const double rho_Y_2_i = rho_pre*f_sm;
                                    const double u_i       = u_He*(1.0 - f_sm) + u_pre*f_sm;
                                    const double v_i       = v_He*(1.0 - f_sm) + v_pre*f_sm;
                                    const double p_i       = p_He*(1.0 - f_sm) + p_pre*f_sm;
                                    
                                    const double rho_i = rho_Y_1_i + rho_Y_2_i;
                                    const double Y_1_i = rho_Y_1_i/rho_i;
                                    const double Y_2_i = 1.0 - Y_1_i;
                                    
                                    std::vector<const double*> Y_ptr;
                                    Y_ptr.push_back(&Y_1_i);
                                    Y_ptr.push_back(&Y_2_i);
                                    
                                    const double gamma = d_equation_of_state->
                                        getMixtureThermodynamicPropertyWithMassFraction(
                                            "gamma",
                                            Y_ptr);
                                    
                                    rho_Y_1[idx_cell] = rho_Y_1_i;
                                    rho_Y_2[idx_cell] = rho_Y_2_i;
                                    rho_u[idx_cell]   = rho_i*u_i;
                                    rho_v[idx_cell]   = rho_i*v_i;
                                    E[idx_cell]       = p_i/(gamma - 1.0) +
                                        0.5*rho_i*(u_i*u_i + v_i*v_i);
                                }
                            }
                        }
                    }
                    else if (d_project_name == "2D Poggi Richtmyer-Meshkov instability problem")
                    {
                        if (d_num_species != 2)
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "Please provide only two-species for the problem"
                                << " '2D Poggi Richtmyer-Meshkov instability problem'."
                                << std::endl);
                        }
                        
                        // Compute the characteristic length of the initial interface thickness.
                        const double epsilon_i = 0.001; // epsilon_i for smoothing interface
                        
                        double* rho_Y_1   = partial_density->getPointer(0);
                        double* rho_Y_2   = partial_density->getPointer(1);
                        double* rho_u     = momentum->getPointer(0);
                        double* rho_v     = momentum->getPointer(1);
                        double* E         = total_energy->getPointer(0);
                        
                        // species 0: SF6
                        // species 1: air
                        const double gamma_0 = d_equation_of_state->
                            getSpeciesThermodynamicProperty(
                                "gamma",
                                0);
                        
                        const double gamma_1 = d_equation_of_state->
                            getSpeciesThermodynamicProperty(
                                "gamma",
                                1);
                        
                        NULL_USE(gamma_1);
                        
                        // Unshocked SF6.
                        const double rho_unshocked = 5.972856088;
                        const double u_unshocked   = 0.0;
                        const double v_unshocked   = 0.0;
                        const double p_unshocked   = 101325;
                        
                        // Shocked SF6.
                        const double rho_shocked = 11.94060777;
                        const double u_shocked   = 98.92358848;
                        const double v_shocked   = 0.0;
                        const double p_shocked   = 218274.2532;
                        
                        // Air.
                        const double rho_air = 1.184554102;
                        const double u_air   = 0.0;
                        const double v_air   = 0.0;
                        const double p_air   = 101325;
                        
                        // Shock hits the interface after 0.4 ms.
                        const double L_x_shock = 0.120827285;
                        const double L_x_interface = 0.2;
                        
                        for (int j = 0; j < patch_dims[1]; j++)
                        {
                            for (int i = 0; i < patch_dims[0]; i++)
                            {
                                // Compute index into linear data array.
                                int idx_cell = i + j*patch_dims[0];
                                
                                // Compute the coordinates.
                                double x[2];
                                x[0] = patch_xlo[0] + (i + 0.5)*dx[0];
                                x[1] = patch_xlo[1] + (j + 0.5)*dx[1];
                                
                                if (x[0] < L_x_shock)
                                {
                                    rho_Y_1[idx_cell] = rho_shocked;
                                    rho_Y_2[idx_cell] = 0.0;
                                    rho_u[idx_cell] = rho_shocked*u_shocked;
                                    rho_v[idx_cell] = rho_shocked*v_shocked;
                                    E[idx_cell]     = p_shocked/(gamma_0 - 1.0) +
                                        0.5*rho_shocked*(u_shocked*u_shocked + v_shocked*v_shocked);
                                }
                                else
                                {
                                    const double f_sm = 0.5*(1.0 + erf((x[0] - L_x_interface)/epsilon_i));
                                    
                                    // Smooth the primitive quantities;
                                    const double rho_Y_1_i = rho_unshocked*(1.0 - f_sm);
                                    const double rho_Y_2_i = rho_air*f_sm;
                                    const double u_i = u_unshocked*(1.0 - f_sm) + u_air*f_sm;
                                    const double v_i = v_unshocked*(1.0 - f_sm) + v_air*f_sm;
                                    const double p_i = p_unshocked*(1.0 - f_sm) + p_air*f_sm;
                                    
                                    const double rho_i = rho_Y_1_i + rho_Y_2_i;
                                    const double Y_1_i = rho_Y_1_i/rho_i;
                                    const double Y_2_i = 1.0 - Y_1_i;
                                    
                                    std::vector<const double*> Y_ptr;
                                    Y_ptr.reserve(2);
                                    Y_ptr.push_back(&Y_1_i);
                                    Y_ptr.push_back(&Y_2_i);
                                    
                                    const double gamma = d_equation_of_state->
                                        getMixtureThermodynamicPropertyWithMassFraction(
                                            "gamma",
                                            Y_ptr);
                                    
                                    rho_Y_1[idx_cell] = rho_Y_1_i;
                                    rho_Y_2[idx_cell] = rho_Y_2_i;
                                    rho_u[idx_cell]   = rho_i*u_i;
                                    rho_v[idx_cell]   = rho_i*v_i;
                                    E[idx_cell]       = p_i/(gamma - 1.0) +
                                        0.5*rho_i*(u_i*u_i + v_i*v_i);
                                }
                            }
                        }
                    }
                    else if (d_project_name == "2D Poggi Richtmyer-Meshkov instability problem with perturbations")
                    {
                        if (d_num_species != 2)
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "Please provide only two-species for the problem"
                                << " '2D Poggi Richtmyer-Meshkov instability problem'."
                                << std::endl);
                        }
                        
                        // Compute the characteristic length of the initial interface thickness.
                        const double epsilon_i = 0.001; // epsilon_i for smoothing interface
                        
                        double* rho_Y_1   = partial_density->getPointer(0);
                        double* rho_Y_2   = partial_density->getPointer(1);
                        double* rho_u     = momentum->getPointer(0);
                        double* rho_v     = momentum->getPointer(1);
                        double* E         = total_energy->getPointer(0);
                        
                        // species 0: SF6
                        // species 1: air
                        const double gamma_0 = d_equation_of_state->
                            getSpeciesThermodynamicProperty(
                                "gamma",
                                0);
                        
                        const double gamma_1 = d_equation_of_state->
                            getSpeciesThermodynamicProperty(
                                "gamma",
                                1);
                        
                        NULL_USE(gamma_1);
                        
                        // Unshocked SF6.
                        const double rho_unshocked = 5.972856088;
                        const double u_unshocked   = 0.0;
                        const double v_unshocked   = 0.0;
                        const double p_unshocked   = 101325;
                        
                        // Shocked SF6.
                        const double rho_shocked = 11.94060777;
                        const double u_shocked   = 98.92358848;
                        const double v_shocked   = 0.0;
                        const double p_shocked   = 218274.2532;
                        
                        // Air.
                        const double rho_air = 1.184554102;
                        const double u_air   = 0.0;
                        const double v_air   = 0.0;
                        const double p_air   = 101325;
                        
                        // Shock hits the interface after 0.4 ms.
                        const double L_x_shock = 0.120827285;
                        const double L_x_interface = 0.2;
                        
                        // Perturbations due to S mode.
                        double h = 0.001; // Length scale of perturbations.
                        double lambda[] = {0.8, 1.0, 1.25, 1.6, 2.0, 2.5}; // Wavelength of perturbations.

                        for (int j = 0; j < patch_dims[1]; j++)
                        {
                            for (int i = 0; i < patch_dims[0]; i++)
                            {
                                // Compute index into linear data array.
                                int idx_cell = i + j*patch_dims[0];
                                
                                // Compute the coordinates.
                                double x[2];
                                x[0] = patch_xlo[0] + (i + 0.5)*dx[0];
                                x[1] = patch_xlo[1] + (j + 0.5)*dx[1];
                                
                                double x_2 = 0.0;
                                
                                double Chi = 0.0;
                                for (int m = 0; m < 6; m++)
                                {
                                    double kappa = 2*M_PI/(lambda[m]*h);
                                    Chi += h/6.0*(cos(kappa*x[1])*cos(kappa*x_2));
                                }
                                
                                double Chi_phi = computePhiModeLocation(x[1],0.0);
                                
                                if (x[0] < L_x_shock)
                                {
                                    rho_Y_1[idx_cell] = rho_shocked;
                                    rho_Y_2[idx_cell] = 0.0;
                                    rho_u[idx_cell] = rho_shocked*u_shocked;
                                    rho_v[idx_cell] = rho_shocked*v_shocked;
                                    E[idx_cell]     = p_shocked/(gamma_0 - 1.0) +
                                        0.5*rho_shocked*(u_shocked*u_shocked + v_shocked*v_shocked);
                                }
                                else
                                {
                                    const double f_sm = 0.5*(1.0 + erf((x[0] - (L_x_interface + Chi + Chi_phi))/epsilon_i));
                                    
                                    // Smooth the primitive quantities;
                                    const double rho_Y_1_i = rho_unshocked*(1.0 - f_sm);
                                    const double rho_Y_2_i = rho_air*f_sm;
                                    const double u_i = u_unshocked*(1.0 - f_sm) + u_air*f_sm;
                                    const double v_i = v_unshocked*(1.0 - f_sm) + v_air*f_sm;
                                    const double p_i = p_unshocked*(1.0 - f_sm) + p_air*f_sm;
                                    
                                    const double rho_i = rho_Y_1_i + rho_Y_2_i;
                                    const double Y_1_i = rho_Y_1_i/rho_i;
                                    const double Y_2_i = 1.0 - Y_1_i;
                                    
                                    std::vector<const double*> Y_ptr;
                                    Y_ptr.reserve(2);
                                    Y_ptr.push_back(&Y_1_i);
                                    Y_ptr.push_back(&Y_2_i);
                                    
                                    const double gamma = d_equation_of_state->
                                        getMixtureThermodynamicPropertyWithMassFraction(
                                            "gamma",
                                            Y_ptr);
                                    
                                    rho_Y_1[idx_cell] = rho_Y_1_i;
                                    rho_Y_2[idx_cell] = rho_Y_2_i;
                                    rho_u[idx_cell]   = rho_i*u_i;
                                    rho_v[idx_cell]   = rho_i*v_i;
                                    E[idx_cell]       = p_i/(gamma - 1.0) +
                                        0.5*rho_i*(u_i*u_i + v_i*v_i);
                                }
                            }
                        }
                    }
                    else
                    {
                        TBOX_ERROR(d_object_name
                            << ": "
                            << "Cannot initialize data for unknown 2D problem"
                            << " for muti-species flow with conservative four-equation model with name = '"
                            << d_project_name
                            << "'."
                            << std::endl);
                    }
                }
                else if (d_dim == tbox::Dimension(3))
                {
                    if (d_project_name == "3D Poggi Richtmyer-Meshkov instability problem")
                    {
                        if (d_num_species != 2)
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "Please provide only two-species for the problem"
                                << " '3D Poggi Richtmyer-Meshkov instability problem'."
                                << std::endl);
                        }
                        
                        // Compute the characteristic length of the initial interface thickness.
                        const double epsilon_i = 0.002; // epsilon_i for smoothing interface
                        
                        double* rho_Y_1   = partial_density->getPointer(0);
                        double* rho_Y_2   = partial_density->getPointer(1);
                        double* rho_u     = momentum->getPointer(0);
                        double* rho_v     = momentum->getPointer(1);
                        double* rho_w     = momentum->getPointer(2);
                        double* E         = total_energy->getPointer(0);
                        
                        // species 0: SF6
                        // species 1: air
                        const double gamma_0 = d_equation_of_state->
                            getSpeciesThermodynamicProperty(
                                "gamma",
                                0);
                        
                        const double gamma_1 = d_equation_of_state->
                            getSpeciesThermodynamicProperty(
                                "gamma",
                                1);
                        
                        NULL_USE(gamma_1);
                        
                        // Unshocked SF6.
                        const double rho_unshocked = 5.972856088;
                        const double u_unshocked   = 0.0;
                        const double v_unshocked   = 0.0;
                        const double w_unshocked   = 0.0;
                        const double p_unshocked   = 101325;
                        
                        // Shocked SF6.
                        const double rho_shocked = 11.94060777;
                        const double u_shocked   = 98.92358848;
                        const double v_shocked   = 0.0;
                        const double w_shocked   = 0.0;
                        const double p_shocked   = 218274.2532;
                        
                        // Air.
                        const double rho_air = 1.184554102;
                        const double u_air   = 0.0;
                        const double v_air   = 0.0;
                        const double w_air   = 0.0;
                        const double p_air   = 101325;
                        
                        // Shock hits the interface after 0.4 ms.
                        const double L_x_shock = 0.120827285;
                        const double L_x_interface = 0.2;
                        
                        for (int k = 0; k < patch_dims[2]; k++)
                        {
                            for (int j = 0; j < patch_dims[1]; j++)
                            {
                                for (int i = 0; i < patch_dims[0]; i++)
                                {
                                    // Compute the linear index of current cell.
                                    int idx_cell = i +
                                        j*patch_dims[0] +
                                        k*patch_dims[0]*patch_dims[1];
                                    
                                    // Compute the coordinates.
                                    double x[3];
                                    x[0] = patch_xlo[0] + (i + 0.5)*dx[0];
                                    x[1] = patch_xlo[1] + (j + 0.5)*dx[1];
                                    x[2] = patch_xlo[2] + (k + 0.5)*dx[2];
                                    
                                    if (x[0] < L_x_shock)
                                    {
                                        rho_Y_1[idx_cell] = rho_shocked;
                                        rho_Y_2[idx_cell] = 0.0;
                                        rho_u[idx_cell] = rho_shocked*u_shocked;
                                        rho_v[idx_cell] = rho_shocked*v_shocked;
                                        rho_w[idx_cell] = rho_shocked*w_shocked;
                                        E[idx_cell]     = p_shocked/(gamma_0 - 1.0) +
                                            0.5*rho_shocked*(u_shocked*u_shocked +
                                                v_shocked*v_shocked +
                                                w_shocked*w_shocked);
                                    }
                                    else
                                    {
                                        const double f_sm = 0.5*(1.0 + erf((x[0] - L_x_interface)/epsilon_i));
                                        
                                        // Smooth the primitive quantities;
                                        const double rho_Y_1_i = rho_unshocked*(1.0 - f_sm);
                                        const double rho_Y_2_i = rho_air*f_sm;
                                        const double u_i = u_unshocked*(1.0 - f_sm) + u_air*f_sm;
                                        const double v_i = v_unshocked*(1.0 - f_sm) + v_air*f_sm;
                                        const double w_i = w_unshocked*(1.0 - f_sm) + w_air*f_sm;
                                        const double p_i = p_unshocked*(1.0 - f_sm) + p_air*f_sm;
                                        
                                        const double rho_i = rho_Y_1_i + rho_Y_2_i;
                                        const double Y_1_i = rho_Y_1_i/rho_i;
                                        const double Y_2_i = 1.0 - Y_1_i;
                                        
                                        std::vector<const double*> Y_ptr;
                                        Y_ptr.reserve(2);
                                        Y_ptr.push_back(&Y_1_i);
                                        Y_ptr.push_back(&Y_2_i);
                                        
                                        const double gamma = d_equation_of_state->
                                            getMixtureThermodynamicPropertyWithMassFraction(
                                                "gamma",
                                                Y_ptr);
                                        
                                        rho_Y_1[idx_cell] = rho_Y_1_i;
                                        rho_Y_2[idx_cell] = rho_Y_2_i;
                                        rho_u[idx_cell]   = rho_i*u_i;
                                        rho_v[idx_cell]   = rho_i*v_i;
                                        rho_w[idx_cell]   = rho_i*w_i;
                                        E[idx_cell]       = p_i/(gamma - 1.0) +
                                            0.5*rho_i*(u_i*u_i + v_i*v_i + w_i*w_i);
                                    }
                                }
                            }
                        }
                    }
                    else if (d_project_name == "3D Poggi Richtmyer-Meshkov instability problem with perturbations")
                    {
                        if (d_num_species != 2)
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "Please provide only two-species for the problem"
                                << " '3D Poggi Richtmyer-Meshkov instability problem'."
                                << std::endl);
                        }
                        
                        // Compute the characteristic length of the initial interface thickness.
                        const double epsilon_i = 0.001; // epsilon_i for smoothing interface
                        
                        double* rho_Y_1   = partial_density->getPointer(0);
                        double* rho_Y_2   = partial_density->getPointer(1);
                        double* rho_u     = momentum->getPointer(0);
                        double* rho_v     = momentum->getPointer(1);
                        double* rho_w     = momentum->getPointer(2);
                        double* E         = total_energy->getPointer(0);
                        
                        // species 0: SF6
                        // species 1: air
                        const double gamma_0 = d_equation_of_state->
                            getSpeciesThermodynamicProperty(
                                "gamma",
                                0);
                        
                        const double gamma_1 = d_equation_of_state->
                            getSpeciesThermodynamicProperty(
                                "gamma",
                                1);
                        
                        NULL_USE(gamma_1);
                        
                        // Unshocked SF6.
                        const double rho_unshocked = 5.972856088;
                        const double u_unshocked   = 0.0;
                        const double v_unshocked   = 0.0;
                        const double w_unshocked   = 0.0;
                        const double p_unshocked   = 101325;
                        
                        // Shocked SF6.
                        const double rho_shocked = 11.94060777;
                        const double u_shocked   = 98.92358848;
                        const double v_shocked   = 0.0;
                        const double w_shocked   = 0.0;
                        const double p_shocked   = 218274.2532;
                        
                        // Air.
                        const double rho_air = 1.184554102;
                        const double u_air   = 0.0;
                        const double v_air   = 0.0;
                        const double w_air   = 0.0;
                        const double p_air   = 101325;
                        
                        // Shock hits the interface after 0.4 ms.
                        const double L_x_shock = 0.120827285;
                        const double L_x_interface = 0.2;
                        
                        // Perturbations due to S mode.
                        double h = 0.001; // Length scale of perturbations.
                        double lambda[]={0.8, 1.0, 1.25, 1.6, 2.0, 2.5}; // Wavelength of perturbations.
                        
                        for (int k = 0; k < patch_dims[2]; k++)
                        {
                            for (int j = 0; j < patch_dims[1]; j++)
                            {
                                for (int i = 0; i < patch_dims[0]; i++)
                                {
                                    // Compute the linear index of current cell.
                                    int idx_cell = i +
                                        j*patch_dims[0] +
                                        k*patch_dims[0]*patch_dims[1];
                                    
                                    // Compute the coordinates.
                                    double x[3];
                                    x[0] = patch_xlo[0] + (i + 0.5)*dx[0];
                                    x[1] = patch_xlo[1] + (j + 0.5)*dx[1];
                                    x[2] = patch_xlo[2] + (k + 0.5)*dx[2];
                                    
                                    double Chi = 0.0;
                                    for (int m = 0; m < 6; m++)
                                    {
                                        double kappa = 2*M_PI/(lambda[m]*h);
                                        Chi += h/6.0*(cos(kappa*x[1])*cos(kappa*x[2]));
                                    }
                                    
                                    double Chi_phi = computePhiModeLocation(x[1],x[2]);
                                    
                                    if (x[0] < L_x_shock)
                                    {
                                        rho_Y_1[idx_cell] = rho_shocked;
                                        rho_Y_2[idx_cell] = 0.0;
                                        rho_u[idx_cell] = rho_shocked*u_shocked;
                                        rho_v[idx_cell] = rho_shocked*v_shocked;
                                        rho_w[idx_cell] = rho_shocked*w_shocked;
                                        E[idx_cell]     = p_shocked/(gamma_0 - 1.0) +
                                            0.5*rho_shocked*(u_shocked*u_shocked +
                                                v_shocked*v_shocked +
                                                w_shocked*w_shocked);
                                    }
                                    else
                                    {
                                        const double f_sm = 0.5*(1.0 + erf((x[0] - (L_x_interface + Chi + Chi_phi))/epsilon_i));
                                        
                                        // Smooth the primitive quantities;
                                        const double rho_Y_1_i = rho_unshocked*(1.0 - f_sm);
                                        const double rho_Y_2_i = rho_air*f_sm;
                                        const double u_i = u_unshocked*(1.0 - f_sm) + u_air*f_sm;
                                        const double v_i = v_unshocked*(1.0 - f_sm) + v_air*f_sm;
                                        const double w_i = w_unshocked*(1.0 - f_sm) + w_air*f_sm;
                                        const double p_i = p_unshocked*(1.0 - f_sm) + p_air*f_sm;
                                        
                                        const double rho_i = rho_Y_1_i + rho_Y_2_i;
                                        const double Y_1_i = rho_Y_1_i/rho_i;
                                        const double Y_2_i = 1.0 - Y_1_i;
                                        
                                        std::vector<const double*> Y_ptr;
                                        Y_ptr.reserve(2);
                                        Y_ptr.push_back(&Y_1_i);
                                        Y_ptr.push_back(&Y_2_i);
                                        
                                        const double gamma = d_equation_of_state->
                                            getMixtureThermodynamicPropertyWithMassFraction(
                                                "gamma",
                                                Y_ptr);
                                        
                                        rho_Y_1[idx_cell] = rho_Y_1_i;
                                        rho_Y_2[idx_cell] = rho_Y_2_i;
                                        rho_u[idx_cell]   = rho_i*u_i;
                                        rho_v[idx_cell]   = rho_i*v_i;
                                        rho_w[idx_cell]   = rho_i*w_i;
                                        E[idx_cell]       = p_i/(gamma - 1.0) +
                                            0.5*rho_i*(u_i*u_i + v_i*v_i + w_i*w_i);
                                    }
                                }
                            }
                        }
                    }
                    else
                    {
                        TBOX_ERROR(d_object_name
                            << ": "
                            << "Cannot initialize data for unknown 3D problem"
                            << " for muti-species flow with conservative four-equation model with name = '"
                            << d_project_name
                            << "'."
                            << std::endl);
                    }
                }
                
                break;
            }
            case FIVE_EQN_ALLAIRE:
            {
                boost::shared_ptr<pdat::CellData<double> > partial_density(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_partial_density, data_context)));
                
                boost::shared_ptr<pdat::CellData<double> > momentum(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_momentum, data_context)));
                
                boost::shared_ptr<pdat::CellData<double> > total_energy(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_total_energy, data_context)));
                
                boost::shared_ptr<pdat::CellData<double> > volume_fraction(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_volume_fraction, data_context)));
                
#ifdef DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(partial_density);
                TBOX_ASSERT(momentum);
                TBOX_ASSERT(total_energy);
                TBOX_ASSERT(volume_fraction);
#endif
                
                if (d_dim == tbox::Dimension(1))
                {
                    // NOT YET IMPLEMENTED.
                    
                    TBOX_ERROR(d_object_name
                        << ": "
                        << "Cannot initialize data for unknown 1D problem"
                        << " for muti-species flow with five-equation model by Allaire with name = '"
                        << d_project_name
                        << "'."
                        << std::endl);
                }
                else if (d_dim == tbox::Dimension(2))
                {
                    if (d_project_name == "2D shock-bubble interaction")
                    {
                        if (d_num_species != 2)
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "Please provide only two-species for the problem"
                                << " '2D shock-bubble interaction'."
                                << std::endl);
                        }
                        
                        // Characteristic length of the problem.
                        const double D = 1.0;
                        
                        // Compute the characteristic length of the initial interface thickness.
                        const double C_epsilon = 3.0;
                        const double epsilon_i = C_epsilon*sqrt(dx[0]*dx[1]);
                        
                        double* Z_rho_1   = partial_density->getPointer(0);
                        double* Z_rho_2   = partial_density->getPointer(1);
                        double* rho_u     = momentum->getPointer(0);
                        double* rho_v     = momentum->getPointer(1);
                        double* E         = total_energy->getPointer(0);
                        double* Z_1       = volume_fraction->getPointer(0);
                        double* Z_2       = volume_fraction->getPointer(1);
                        
                        // species 0: He
                        // species 1: air
                        const double gamma_0 = d_equation_of_state->
                            getSpeciesThermodynamicProperty(
                                "gamma",
                                0);
                        
                        const double gamma_1 = d_equation_of_state->
                            getSpeciesThermodynamicProperty(
                                "gamma",
                                1);
                        
                        NULL_USE(gamma_0);
                        
                        // He, pre-shock condition.
                        const double rho_He = 0.1819;
                        const double u_He   = 0.0;
                        const double v_He   = 0.0;
                        const double p_He   = 1.0/1.4;
                        const double Z_He   = 1.0;
                        
                        // air, pre-shock condition.
                        const double rho_pre = 1.0;
                        const double u_pre   = 0.0;
                        const double v_pre   = 0.0;
                        const double p_pre   = 1.0/1.4;
                        const double Z_pre   = 0.0;
                        
                        // air, post-shock condition.
                        const double rho_post = 1.3764;
                        const double u_post   = -0.3336;
                        const double v_post   = 0.0;
                        const double p_post   = 1.5698/1.4;
                        const double Z_post   = 0.0;
                        
                        for (int j = 0; j < patch_dims[1]; j++)
                        {
                            for (int i = 0; i < patch_dims[0]; i++)
                            {
                                // Compute index into linear data array.
                                int idx_cell = i + j*patch_dims[0];
                                
                                // Compute the coordinates.
                                double x[2];
                                x[0] = patch_xlo[0] + (i + 0.5)*dx[0];
                                x[1] = patch_xlo[1] + (j + 0.5)*dx[1];
                                
                                if (x[0] > 4.5*D)
                                {
                                    Z_rho_1[idx_cell] = 0.0;
                                    Z_rho_2[idx_cell] = rho_post;
                                    rho_u[idx_cell]   = rho_post*u_post;
                                    rho_v[idx_cell]   = rho_post*v_post;
                                    E[idx_cell]       = p_post/(gamma_1 - 1.0) +
                                        0.5*rho_post*(u_post*u_post + v_post*v_post);
                                    Z_1[idx_cell]     = Z_post;
                                    Z_2[idx_cell]     = 1.0 - Z_post;
                                }
                                else
                                {
                                    // Compute the distance from the initial material interface.
                                    const double dR = sqrt(pow(x[0] - 3.5, 2) + x[1]*x[1]) - 0.5*D;
                                    
                                    const double f_sm = 0.5*(1.0 + erf(dR/epsilon_i));
                                    
                                    // Smooth the primitive quantity.
                                    const double Z_rho_1_i = rho_He*(1.0 - f_sm);
                                    const double Z_rho_2_i = rho_pre*f_sm;
                                    const double u_i       = u_He*(1.0 - f_sm) + u_pre*f_sm;
                                    const double v_i       = v_He*(1.0 - f_sm) + v_pre*f_sm;
                                    const double p_i       = p_He*(1.0 - f_sm) + p_pre*f_sm;
                                    const double Z_1_i     = Z_He*(1.0 - f_sm) + Z_pre*f_sm;
                                    
                                    const double rho_i = Z_rho_1_i + Z_rho_2_i;
                                    const double Z_2_i = 1.0 - Z_1_i;
                                    
                                    std::vector<const double*> Z_ptr;
                                    Z_ptr.push_back(&Z_1_i);
                                    Z_ptr.push_back(&Z_2_i);
                                    
                                    const double gamma = d_equation_of_state->
                                        getMixtureThermodynamicPropertyWithVolumeFraction(
                                            "gamma",
                                            Z_ptr);
                                    
                                    Z_rho_1[idx_cell] = Z_rho_1_i;
                                    Z_rho_2[idx_cell] = Z_rho_2_i;
                                    rho_u[idx_cell]   = rho_i*u_i;
                                    rho_v[idx_cell]   = rho_i*v_i;
                                    E[idx_cell]       = p_i/(gamma - 1.0) +
                                        0.5*rho_i*(u_i*u_i + v_i*v_i);
                                    Z_1[idx_cell]     = Z_1_i;
                                    Z_2[idx_cell]     = Z_2_i;
                                }
                            }
                        }
                    }
                    else if (d_project_name == "2D shock-bubble interaction with AMR")
                    {
                        if (d_num_species != 2)
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "Please provide only two-species for the problem"
                                << " '2D shock-bubble interaction with AMR'."
                                << std::endl);
                        }
                        
                        // Characteristic length of the problem.
                        const double D = 1.0;
                        
                        // Compute the characteristic length of the initial interface thickness.
                        const double C_epsilon = 3.0;
                        const double epsilon_i = C_epsilon*0.0025; // epsilon_i for smoothing interface
                        
                        double* Z_rho_1   = partial_density->getPointer(0);
                        double* Z_rho_2   = partial_density->getPointer(1);
                        double* rho_u     = momentum->getPointer(0);
                        double* rho_v     = momentum->getPointer(1);
                        double* E         = total_energy->getPointer(0);
                        double* Z_1       = volume_fraction->getPointer(0);
                        double* Z_2       = volume_fraction->getPointer(1);
                        
                        // species 0: He
                        // species 1: air
                        const double gamma_0 = d_equation_of_state->
                            getSpeciesThermodynamicProperty(
                                "gamma",
                                0);
                        
                        const double gamma_1 = d_equation_of_state->
                            getSpeciesThermodynamicProperty(
                                "gamma",
                                1);
                        
                        NULL_USE(gamma_0);
                        
                        // He, pre-shock condition.
                        const double rho_He = 0.1819;
                        const double u_He   = 0.0;
                        const double v_He   = 0.0;
                        const double p_He   = 1.0/1.4;
                        const double Z_He   = 1.0;
                        
                        // air, pre-shock condition.
                        const double rho_pre = 1.0;
                        const double u_pre   = 0.0;
                        const double v_pre   = 0.0;
                        const double p_pre   = 1.0/1.4;
                        const double Z_pre   = 0.0;
                        
                        // air, post-shock condition.
                        const double rho_post = 1.3764;
                        const double u_post   = -0.3336;
                        const double v_post   = 0.0;
                        const double p_post   = 1.5698/1.4;
                        const double Z_post   = 0.0;
                        
                        for (int j = 0; j < patch_dims[1]; j++)
                        {
                            for (int i = 0; i < patch_dims[0]; i++)
                            {
                                // Compute index into linear data array.
                                int idx_cell = i + j*patch_dims[0];
                                
                                // Compute the coordinates.
                                double x[2];
                                x[0] = patch_xlo[0] + (i + 0.5)*dx[0];
                                x[1] = patch_xlo[1] + (j + 0.5)*dx[1];
                                
                                if (x[0] > 4.5*D)
                                {
                                    Z_rho_1[idx_cell] = 0.0;
                                    Z_rho_2[idx_cell] = rho_post;
                                    rho_u[idx_cell]   = rho_post*u_post;
                                    rho_v[idx_cell]   = rho_post*v_post;
                                    E[idx_cell]       = p_post/(gamma_1 - 1.0) +
                                        0.5*rho_post*(u_post*u_post + v_post*v_post);
                                    Z_1[idx_cell]     = Z_post;
                                    Z_2[idx_cell]     = 1.0 - Z_post;
                                }
                                else
                                {
                                    // Compute the distance from the initial material interface.
                                    const double dR = sqrt(pow(x[0] - 3.5, 2) + x[1]*x[1]) - 0.5*D;
                                    
                                    const double f_sm = 0.5*(1.0 + erf(dR/epsilon_i));
                                    
                                    // Smooth the primitive quantity.
                                    const double Z_rho_1_i = rho_He*(1.0 - f_sm);
                                    const double Z_rho_2_i = rho_pre*f_sm;
                                    const double u_i       = u_He*(1.0 - f_sm) + u_pre*f_sm;
                                    const double v_i       = v_He*(1.0 - f_sm) + v_pre*f_sm;
                                    const double p_i       = p_He*(1.0 - f_sm) + p_pre*f_sm;
                                    const double Z_1_i     = Z_He*(1.0 - f_sm) + Z_pre*f_sm;
                                    
                                    const double rho_i = Z_rho_1_i + Z_rho_2_i;
                                    const double Z_2_i = 1.0 - Z_1_i;
                                    
                                    std::vector<const double*> Z_ptr;
                                    Z_ptr.push_back(&Z_1_i);
                                    Z_ptr.push_back(&Z_2_i);
                                    
                                    const double gamma = d_equation_of_state->
                                        getMixtureThermodynamicPropertyWithVolumeFraction(
                                            "gamma",
                                            Z_ptr);
                                    
                                    Z_rho_1[idx_cell] = Z_rho_1_i;
                                    Z_rho_2[idx_cell] = Z_rho_2_i;
                                    rho_u[idx_cell]   = rho_i*u_i;
                                    rho_v[idx_cell]   = rho_i*v_i;
                                    E[idx_cell]       = p_i/(gamma - 1.0) +
                                        0.5*rho_i*(u_i*u_i + v_i*v_i);
                                    Z_1[idx_cell]     = Z_1_i;
                                    Z_2[idx_cell]     = Z_2_i;
                                }
                            }
                        }
                    }
                    else if (d_project_name == "2D Richtmyer-Meshkov instability")
                    {
                        // Initialize data for a 2D Richtmyer-Meshkov instability.
                        
                        if (d_num_species != 2)
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "Please provide only two-species for the problem"
                                << " '2D Richtmyer-Meshkov instability'."
                                << std::endl);
                        }
                        
                        const double D = 1.0;
                        
                        double* Z_rho_1   = partial_density->getPointer(0);
                        double* Z_rho_2   = partial_density->getPointer(1);
                        double* rho_u     = momentum->getPointer(0);
                        double* rho_v     = momentum->getPointer(1);
                        double* E         = total_energy->getPointer(0);
                        double* Z_1       = volume_fraction->getPointer(0);
                        double* Z_2       = volume_fraction->getPointer(1);
                        
                        // species 0: SF6
                        // species 1: air
                        const double gamma_0 = d_equation_of_state->
                            getSpeciesThermodynamicProperty(
                                "gamma",
                                0);
                        
                        const double gamma_1 = d_equation_of_state->
                            getSpeciesThermodynamicProperty(
                                "gamma",
                                1);
                        
                        NULL_USE(gamma_0);
                        
                        // SF6, pre-shock condition.
                        const double rho_SF6 = 5.04;
                        const double u_SF6   = 1.24;
                        const double v_SF6   = 0.0;
                        const double p_SF6   = 1.0/1.4;
                        const double Z_SF6   = 1.0;
                        
                        // air, pre-shock condition.
                        const double rho_pre = 1.0;
                        const double u_pre   = 1.24;
                        const double v_pre   = 0.0;
                        const double p_pre   = 1.0/1.4;
                        const double Z_pre   = 0.0;
                        
                        // air, post-shock condition.
                        const double rho_post = 1.4112;
                        const double u_post   = 0.8787;
                        const double v_post   = 0.0;
                        const double p_post   = 1.6272/1.4;
                        const double Z_post   = 0.0;
                        
                        // Compute the characteristic length of the initial interface thickness.
                        const double epsilon_i = (6.0)*sqrt(dx[0]*dx[1]);
                        
                        for (int j = 0; j < patch_dims[1]; j++)
                        {
                            for (int i = 0; i < patch_dims[0]; i++)
                            {
                                // Compute index into linear data array.
                                int idx_cell = i + j*patch_dims[0];
                                
                                // Compute the coordinates.
                                double x[2];
                                x[0] = patch_xlo[0] + (i + 0.5)*dx[0];
                                x[1] = patch_xlo[1] + (j + 0.5)*dx[1];
                                
                                if (x[0] > 0.7*D)
                                {
                                    Z_rho_1[idx_cell] = 0.0;
                                    Z_rho_2[idx_cell] = rho_post;
                                    rho_u[idx_cell]   = rho_post*u_post;
                                    rho_v[idx_cell]   = rho_post*v_post;
                                    E[idx_cell]       = p_post/(gamma_1 - 1.0) +
                                        0.5*rho_post*(u_post*u_post + v_post*v_post);
                                    Z_1[idx_cell]     = Z_post;
                                    Z_2[idx_cell]     = 1.0 - Z_post;
                                }
                                else
                                {
                                    // Compute the distance from the initial material interface.
                                    const double dR = x[0] - (0.4 - 0.1*sin(2*M_PI*(x[1]/D + 0.25)))*D;
                                    
                                    const double f_sm = 0.5*(1.0 + erf(dR/epsilon_i));
                                    
                                    // Smooth the primitive quantity.
                                    const double Z_rho_1_i = rho_SF6*(1.0 - f_sm);
                                    const double Z_rho_2_i = rho_pre*f_sm;
                                    const double u_i       = u_SF6*(1.0 - f_sm) + u_pre*f_sm;
                                    const double v_i       = v_SF6*(1.0 - f_sm) + v_pre*f_sm;
                                    const double p_i       = p_SF6*(1.0 - f_sm) + p_pre*f_sm;
                                    const double Z_1_i     = Z_SF6*(1.0 - f_sm) + Z_pre*f_sm;
                                    
                                    const double rho_i = Z_rho_1_i + Z_rho_2_i;
                                    const double Z_2_i = 1.0 - Z_1_i;
                                    
                                    std::vector<const double*> Z_ptr;
                                    Z_ptr.push_back(&Z_1_i);
                                    Z_ptr.push_back(&Z_2_i);
                                    
                                    const double gamma = d_equation_of_state->
                                        getMixtureThermodynamicPropertyWithVolumeFraction(
                                            "gamma",
                                            Z_ptr);
                                    
                                    Z_rho_1[idx_cell] = Z_rho_1_i;
                                    Z_rho_2[idx_cell] = Z_rho_2_i;
                                    rho_u[idx_cell]   = rho_i*u_i;
                                    rho_v[idx_cell]   = rho_i*v_i;
                                    E[idx_cell]       = p_i/(gamma - 1.0) +
                                        0.5*rho_i*(u_i*u_i + v_i*v_i);
                                    Z_1[idx_cell]     = Z_1_i;
                                    Z_2[idx_cell]     = Z_2_i;
                                }
                            }
                        }
                    }
                    else if (d_project_name == "2D advection of material interface")
                    {
                        // Initialize data for a 2D material interface advection problem.
                        
                        if (d_num_species != 2)
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "Please provide only two-species for the problem"
                                << " '2D advection of material interface'."
                                << std::endl);
                        }
                        
                        double* Z_rho_1   = partial_density->getPointer(0);
                        double* Z_rho_2   = partial_density->getPointer(1);
                        double* rho_u     = momentum->getPointer(0);
                        double* rho_v     = momentum->getPointer(1);
                        double* E         = total_energy->getPointer(0);
                        double* Z_1       = volume_fraction->getPointer(0);
                        double* Z_2       = volume_fraction->getPointer(1);
                        
                        const double x_a = 1.0/3*(domain_xlo[0] + domain_xhi[0]);
                        const double x_b = 2.0/3*(domain_xlo[0] + domain_xhi[0]);
                        
                        const double y_a = 1.0/3*(domain_xlo[1] + domain_xhi[1]);
                        const double y_b = 2.0/3*(domain_xlo[1] + domain_xhi[1]);
                        
                        // material initial conditions.
                        double gamma_m  = d_equation_of_state->
                            getSpeciesThermodynamicProperty(
                                "gamma",
                                0);
                        double rho_m    = 10.0;
                        double u_m      = 0.5;
                        double v_m      = 0.5;
                        double p_m      = 1.0/1.4;
                        
                        // ambient initial conditions.
                        double gamma_a  = d_equation_of_state->
                            getSpeciesThermodynamicProperty(
                                "gamma",
                                1);
                        double rho_a    = 1.0;
                        double u_a      = 0.5;
                        double v_a      = 0.5;
                        double p_a      = 1.0/1.4;
                        
                        for (int j = 0; j < patch_dims[1]; j++)
                        {
                            for (int i = 0; i < patch_dims[0]; i++)
                            {
                                // Compute index into linear data array.
                                int idx_cell = i + j*patch_dims[0];
                                
                                // Compute the coordinates.
                                double x[2];
                                x[0] = patch_xlo[0] + (i + 0.5)*dx[0];
                                x[1] = patch_xlo[1] + (j + 0.5)*dx[1];
                                
                                if ((x[0] >= x_a) &&
                                    (x[0] <= x_b) &&
                                    (x[1] >= y_a) &&
                                    (x[1] <= y_b))
                                {
                                    Z_rho_1[idx_cell] = rho_m;
                                    Z_rho_2[idx_cell] = 0.0;
                                    rho_u[idx_cell]   = rho_m*u_m;
                                    rho_v[idx_cell]   = rho_m*v_m;
                                    E[idx_cell]       = p_m/(gamma_m - 1.0) +
                                        0.5*rho_m*(u_m*u_m + v_m*v_m);
                                    Z_1[idx_cell]     = 1.0;
                                    Z_2[idx_cell]     = 0.0;
                                }
                                else
                                {
                                    Z_rho_1[idx_cell] = 0.0;
                                    Z_rho_2[idx_cell] = rho_a;
                                    rho_u[idx_cell]   = rho_a*u_a;
                                    rho_v[idx_cell]   = rho_a*v_a;
                                    E[idx_cell]       = p_a/(gamma_a - 1.0) +
                                        0.5*rho_a*(u_a*u_a + v_a*v_a);
                                    Z_1[idx_cell]     = 0.0;
                                    Z_2[idx_cell]     = 1.0;
                                }
                            }
                        }
                    }
                    else
                    {
                        TBOX_ERROR(d_object_name
                            << ": "
                            << "Cannot initialize data for unknown 2D problem"
                            << " for muti-species flow with five-equation model by Allaire with name = '"
                            << d_project_name
                            << "'."
                            << std::endl);
                    }
                }
                else if (d_dim == tbox::Dimension(3))
                {
                    if (d_project_name == "3D advection of material interface")
                    {
                        // Initialize data for a 3D material interface advection problem.
                        
                        if (d_num_species != 2)
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "Please provide only two-species for the problem"
                                << " '3D advection of material interface'."
                                << std::endl);
                        }
                        
                        double* Z_rho_1   = partial_density->getPointer(0);
                        double* Z_rho_2   = partial_density->getPointer(1);
                        double* rho_u     = momentum->getPointer(0);
                        double* rho_v     = momentum->getPointer(1);
                        double* rho_w     = momentum->getPointer(2);
                        double* E         = total_energy->getPointer(0);
                        double* Z_1       = volume_fraction->getPointer(0);
                        double* Z_2       = volume_fraction->getPointer(1);
                        
                        const double x_a = 1.0/3*(domain_xlo[0] + domain_xhi[0]);
                        const double x_b = 2.0/3*(domain_xlo[0] + domain_xhi[0]);
                        
                        const double y_a = 1.0/3*(domain_xlo[1] + domain_xhi[1]);
                        const double y_b = 2.0/3*(domain_xlo[1] + domain_xhi[1]);
                        
                        const double z_a = 1.0/3*(domain_xlo[2] + domain_xhi[2]);
                        const double z_b = 2.0/3*(domain_xlo[2] + domain_xhi[2]);
                        
                        // material initial conditions.
                        double gamma_m  = d_equation_of_state->
                            getSpeciesThermodynamicProperty(
                                "gamma",
                                0);
                        double rho_m    = 10.0;
                        double u_m      = 0.5;
                        double v_m      = 0.5;
                        double w_m      = 0.5;
                        double p_m      = 1.0/1.4;
                        
                        // ambient initial conditions.
                        double gamma_a  = d_equation_of_state->
                            getSpeciesThermodynamicProperty(
                                "gamma",
                                1);
                        double rho_a    = 1.0;
                        double u_a      = 0.5;
                        double v_a      = 0.5;
                        double w_a      = 0.5;
                        double p_a      = 1.0/1.4;
                        
                        for (int k = 0; k < patch_dims[2]; k++)
                        {
                            for (int j = 0; j < patch_dims[1]; j++)
                            {
                                for (int i = 0; i < patch_dims[0]; i++)
                                {
                                    // Compute the linear index of current cell.
                                    int idx_cell = i +
                                        j*patch_dims[0] +
                                        k*patch_dims[0]*patch_dims[1];
                                    
                                    // Compute the coordinates.
                                    double x[3];
                                    x[0] = patch_xlo[0] + (i + 0.5)*dx[0];
                                    x[1] = patch_xlo[1] + (j + 0.5)*dx[1];
                                    x[2] = patch_xlo[2] + (k + 0.5)*dx[2];
                                    
                                    if ((x[0] >= x_a) &&
                                        (x[0] <= x_b) &&
                                        (x[1] >= y_a) &&
                                        (x[1] <= y_b) &&
                                        (x[2] >= z_a) &&
                                        (x[2] <= z_b))
                                    {
                                        Z_rho_1[idx_cell] = rho_m;
                                        Z_rho_2[idx_cell] = 0.0;
                                        rho_u[idx_cell]   = rho_m*u_m;
                                        rho_v[idx_cell]   = rho_m*v_m;
                                        rho_w[idx_cell]   = rho_m*w_m;
                                        E[idx_cell]       = p_m/(gamma_m - 1.0) +
                                            0.5*rho_m*(u_m*u_m + v_m*v_m + w_m*w_m);
                                        Z_1[idx_cell]     = 1.0;
                                        Z_2[idx_cell]     = 0.0;
                                    }
                                    else
                                    {
                                        Z_rho_1[idx_cell] = 0.0;
                                        Z_rho_2[idx_cell] = rho_a;
                                        rho_u[idx_cell]   = rho_a*u_a;
                                        rho_v[idx_cell]   = rho_a*v_a;
                                        rho_w[idx_cell]   = rho_a*w_a;
                                        E[idx_cell]       = p_a/(gamma_a - 1.0) +
                                            0.5*rho_a*(u_a*u_a + v_a*v_a + w_a*w_a);
                                        Z_1[idx_cell]     = 0.0;
                                        Z_2[idx_cell]     = 1.0;
                                    }
                                }
                            }
                        }
                    }
                    else if (d_project_name == "3D shock-bubble interaction with AMR")
                    {
                        if (d_num_species != 2)
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "Please provide only two-species for the problem"
                                << " '3D shock-bubble interaction with AMR'."
                                << std::endl);
                        }
                        
                        // Compute the characteristic length of the initial interface thickness.
                        // const double C_epsilon = 0.15;
                        const double C_epsilon = 1.5;
                        const double epsilon_i = C_epsilon*396.875e-6; // epsilon_i for smoothing interface.
                        //const double epsilon_i = C_epsilon*pow(dx[0]*dx[1]*dx[2], 1.0/3.0); // epsilon_i for smoothing interface.
                        
                        double* Z_rho_1   = partial_density->getPointer(0);
                        double* Z_rho_2   = partial_density->getPointer(1);
                        double* rho_u     = momentum->getPointer(0);
                        double* rho_v     = momentum->getPointer(1);
                        double* rho_w     = momentum->getPointer(2);
                        double* E         = total_energy->getPointer(0);
                        double* Z_1       = volume_fraction->getPointer(0);
                        double* Z_2       = volume_fraction->getPointer(1);
                        
                        // species 0: Kr
                        // species 1: air
                        const double gamma_0 = d_equation_of_state->
                            getSpeciesThermodynamicProperty(
                                "gamma",
                                0);
                        
                        const double gamma_1 = d_equation_of_state->
                            getSpeciesThermodynamicProperty(
                                "gamma",
                                1);
                        
                        NULL_USE(gamma_0);
                        
                        // Kr, pre-shock condition.
                        const double rho_Kr = 3.485;
                        const double u_Kr   = 0.0;
                        const double v_Kr   = 0.0;
                        const double w_Kr   = 0.0;
                        const double p_Kr   = 1.013e5;
                        const double Z_Kr   = 1.0;
                        
                        // air, pre-shock condition.
                        const double rho_pre = 1.205;
                        const double u_pre   = 0.0;
                        const double v_pre   = 0.0;
                        const double w_pre   = 0.0;
                        const double p_pre   = 1.013e5;
                        const double Z_pre   = 0.0;
                        
                        // air, post-shock condition.
                        const double rho_post = 2.609923442;
                        const double u_post   = 0.0;
                        const double v_post   = 310.1374648;
                        const double w_post   = 0.0;
                        const double p_post   = 3.166131794e5;
                        const double Z_post   = 0.0;
                        
                        for (int k = 0; k < patch_dims[2]; k++)
                        {
                            for (int j = 0; j < patch_dims[1]; j++)
                            {
                                for (int i = 0; i < patch_dims[0]; i++)
                                {
                                    // Compute the linear index of current cell.
                                    int idx_cell = i +
                                        j*patch_dims[0] +
                                        k*patch_dims[0]*patch_dims[1];
                                    
                                    // Compute the coordinates.
                                    double x[3];
                                    x[0] = patch_xlo[0] + (i + 0.5)*dx[0];
                                    x[1] = patch_xlo[1] + (j + 0.5)*dx[1];
                                    x[2] = patch_xlo[2] + (k + 0.5)*dx[2];
                                    
                                    // If y < 2cm,
                                    if (x[1] < 0.02)
                                    {
                                        Z_rho_1[idx_cell] = 0.0;
                                        Z_rho_2[idx_cell] = rho_post;
                                        rho_u[idx_cell]   = rho_post*u_post;
                                        rho_v[idx_cell]   = rho_post*v_post;
                                        rho_w[idx_cell]   = rho_post*w_post;
                                        E[idx_cell]       = p_post/(gamma_1 - 1.0) +
                                            0.5*rho_post*(u_post*u_post + v_post*v_post + w_post*w_post);
                                        Z_1[idx_cell]     = Z_post;
                                        Z_2[idx_cell]     = 1.0 - Z_post;
                                    }
                                    else
                                    {
                                        // Compute the distance from the initial material interface.
                                        const double dR = sqrt(x[0]*x[0] + pow(x[1] - 0.075, 2) + x[2]*x[2]) - 0.0254;
                                        
                                        const double f_sm = 0.5*(1.0 + erf(dR/epsilon_i));
                                        
                                        // Smooth the primitive quantity.
                                        const double Z_rho_1_i = rho_Kr*(1.0 - f_sm);
                                        const double Z_rho_2_i = rho_pre*f_sm;
                                        const double u_i       = u_Kr*(1.0 - f_sm) + u_pre*f_sm;
                                        const double v_i       = v_Kr*(1.0 - f_sm) + v_pre*f_sm;
                                        const double w_i       = w_Kr*(1.0 - f_sm) + w_pre*f_sm;
                                        const double p_i       = p_Kr*(1.0 - f_sm) + p_pre*f_sm;
                                        const double Z_1_i     = Z_Kr*(1.0 - f_sm) + Z_pre*f_sm;
                                        
                                        const double rho_i = Z_rho_1_i + Z_rho_2_i;
                                        const double Z_2_i = 1.0 - Z_1_i;
                                        
                                        std::vector<const double*> Z_ptr;
                                        Z_ptr.push_back(&Z_1_i);
                                        Z_ptr.push_back(&Z_2_i);
                                        
                                        const double gamma = d_equation_of_state->
                                            getMixtureThermodynamicPropertyWithVolumeFraction(
                                                "gamma",
                                                Z_ptr);
                                        
                                        Z_rho_1[idx_cell] = Z_rho_1_i;
                                        Z_rho_2[idx_cell] = Z_rho_2_i;
                                        rho_u[idx_cell]   = rho_i*u_i;
                                        rho_v[idx_cell]   = rho_i*v_i;
                                        rho_w[idx_cell]   = rho_i*w_i;
                                        E[idx_cell]       = p_i/(gamma - 1.0) +
                                            0.5*rho_i*(u_i*u_i + v_i*v_i + w_i*w_i);
                                        Z_1[idx_cell]     = Z_1_i;
                                        Z_2[idx_cell]     = Z_2_i;
                                    }
                                }
                            }
                        }
                    }
                    else
                    {
                        TBOX_ERROR(d_object_name
                            << ": "
                            << "Cannot initialize data for unknown 3D problem"
                            << " for muti-species flow with five-equation model by Allaire with name = '"
                            << d_project_name
                            << "'."
                            << std::endl);
                    }
                }
                
                break;
            }
        }
    }
}

double
InitialConditions::computePhiModeLocation(
    const double& x_1,
    const double& x_2)
{
    double Chi_phi = 0.0;

    // Characteristics of the perturbations.
    std::vector<double> amp {0.0027295181740957884, 0.0017468916314213045, 0.0012131191884870169, 0.00089127124052107388, 0.00068237954352394709, 
    0.00053916408377200752, 0.00043672290785532614, 0.0003609280230209307, 0.00030327979712175423, 0.00025841592180788524, 
    0.00022281781013026847, 0.0001940990701579227, 0.00017059488588098677, 0.00015111519302952461, 0.00013479102094300188, 
    0.00012097587475216789, 0.00010918072696383153, 9.9030137835674859e-05, 9.0232005755232675e-05, 8.2556315284560707e-05, 
    7.5819949280438558e-05, 6.9875665256852187e-05, 6.4603980451971311e-05, 5.9907120419111959e-05, 5.5704452532567118e-05, 
    5.1929002123106551e-05, 4.8524767539480674e-05, 4.5444631410543818e-05, 4.2648721470246693e-05, 4.0103113668992295e-05, 
    3.7778798257381152e-05, 3.5650849620842951e-05, 3.369775523575047e-05, 3.1900869821426308e-05, 3.0243968688041973e-05, 
    2.8712880200876143e-05, 2.7295181740957884e-05, 2.5979946927740995e-05, 2.4757534458918715e-05, 2.361941091699979e-05, 
    2.2558001438808169e-05, 2.1566563350880302e-05, 2.0639078821140177e-05, 1.9770163325275059e-05, 1.8954987320109639e-05, 
    1.8189208990225995e-05, 1.7468916314213047e-05, 1.6790577003280513e-05, 1.6150995112992828e-05, 1.554727333055629e-05, 
    1.497678010477799e-05, 1.4437120920837226e-05, 1.3926113133141779e-05, 1.3441763861351989e-05, 1.2982250530776638e-05, 
    1.2545903701675557e-05, 1.2131191884870169e-05, 1.1736708085335289e-05, 1.1361157852635955e-05, 1.1003348648408319e-05, 
    1.0662180367561673e-05, 0.0025689582815019182, 0.0016797034917512544, 0.0011803321833927735, 0.00087344581571065217, 
    0.00067188139670050181, 0.00053258891201869035, 0.00043239891866863971, 0.00035796959660272629, 0.00030118821231401801, 
    0.0002568958281501918, 0.00022168675525651075, 0.00019324022471474606, 0.00016993109255071054, 0.00015059410615700898, 
    0.00013437627934010036, 0.00012064168725285254, 0.00010890845582427083, 9.880608775007377e-05, 9.0045960382541487e-05, 
    8.2400548651948333e-05, 7.56885455555158e-05, 6.9764042788390742e-05, 6.4508553597537088e-05, 5.9825055870592621e-05, 
    5.5633491446538352e-05, 5.1867328723910461e-05, 4.8470910971734308e-05, 4.539739166895281e-05, 4.2607112961495234e-05, 
    4.0066321821589558e-05, 3.7746145882050649e-05, 3.5621770624414843e-05, 3.3671773928706715e-05, 3.1877584514987308e-05, 
    3.0223038605904926e-05, 2.8694014970783585e-05, 2.7278132907890455e-05, 2.5964501061553279e-05, 2.474350752721394e-05, 
    2.3606643667855463e-05, 2.2546355593976571e-05, 2.1555918452878882e-05, 2.0629329610549177e-05, 1.9761217550014752e-05, 
    1.8946763898278791e-05, 1.8181636463585602e-05, 1.7461931541596408e-05, 1.6784124052856499e-05, 1.6145024319975089e-05, 
    1.5541740493072105e-05, 1.4971645795520265e-05, 1.4432349896078191e-05, 1.3921673823886714e-05, 1.3437627934010035e-05, 
    1.2978392506844757e-05, 1.2542300627665885e-05, 1.2127823045135412e-05, 1.1733554751620799e-05, 1.1358203065158025e-05, 
    1.1000577024063628e-05, 0.0021836145392766303, 0.0015059410615700904, 0.0010918072696383151, 0.00082400548651948339, 
    0.00064223957037547956, 0.00051379165630038371, 0.00041992587293781359, 0.0003493783262842608, 0.00029508304584819337, 
    0.00025244098719961046, 0.00021836145392766304, 0.00019070869338660528, 0.00016797034917512545, 0.00014905218698133997, 
    0.00013314722800467259, 0.00011965011174118526, 0.00010809972966715993, 9.8139979293331706e-05, 8.9492399150681572e-05, 
    8.193675569518313e-05, 7.5297053078504504e-05, 6.9431304905457253e-05, 6.422395703754795e-05, 5.9580205710139991e-05, 
    5.5421688814127688e-05, 5.168318436157705e-05, 4.8310056178686516e-05, 4.5256259881380943e-05, 4.2482773137677635e-05, 
    3.9956350215491862e-05, 3.7648526539252245e-05, 3.5534817563492765e-05, 3.3594069835025089e-05, 3.1807932108909401e-05, 
    3.0160421813213135e-05, 2.8637567728218113e-05, 2.7227113956067708e-05, 2.5918273463224103e-05, 2.4701521937518443e-05, 
    2.3568424600935031e-05, 2.2511490095635372e-05, 2.1524046715393102e-05, 2.0600137162987083e-05, 1.9734428732730505e-05, 
    1.892213638887895e-05, 1.8158956667581126e-05, 1.7441010697097685e-05, 1.6764794927267797e-05, 1.6127138399384272e-05, 
    1.5525165583196803e-05, 1.4956263967648155e-05, 1.4418055723186734e-05, 1.3908372861634588e-05, 1.3425235409017095e-05, 
    1.2966832180977615e-05, 1.2531503812204481e-05, 1.2117727742933577e-05, 1.1724104908867815e-05, 1.1349347917238203e-05, 
    1.0992270522409414e-05, 0.0024262383769740343, 0.0017468916314213045, 0.0012844791407509591, 0.00097049535078961343, 
    0.00075297053078504509, 0.00059825055870592629, 0.00048524767539480672, 0.00040066321821589549, 0.0003359406983502509, 
    0.00028543980905576874, 0.00024534994823332927, 0.00021303556480747613, 0.0001866337213056949, 0.00016480109730389664, 
    0.00014655131136084772, 0.00013114802037697481, 0.00011803321833927734, 0.00010677821707954183, 9.7049535078961349e-05, 
    8.8584768327652359e-05, 8.1175261683146121e-05, 7.4653488522277974e-05, 6.8883739409357429e-05, 6.3755169029974617e-05, 
    5.9176545779854495e-05, 5.5072245631188665e-05, 5.1379165630038368e-05, 4.8044324296515534e-05, 4.5022980191270737e-05, 
    4.2277145000515604e-05, 3.9774399622525144e-05, 3.7486944880285505e-05, 3.539083532052885e-05, 3.3465356923779783e-05, 
    3.1692518712287818e-05, 3.0056635089836621e-05, 2.8543980905576871e-05, 2.7142505149491988e-05, 2.5841592180788528e-05, 
    2.4631861695167862e-05, 2.3505000422783966e-05, 2.2453619941147876e-05, 2.1471136079416234e-05, 2.0551666252015343e-05, 
    1.9689941742800998e-05, 1.8881232505634504e-05, 1.812128248362349e-05, 1.7406253800531131e-05, 1.6732678461889888e-05, 
    1.6097416434033397e-05, 1.5497619157392695e-05, 1.4930697704455594e-05, 1.4394294919424065e-05, 1.3886260981091449e-05, 
    1.3404631916983615e-05, 1.2947610668702227e-05, 1.2513550368347452e-05, 1.2100939536030092e-05, 1.1708388950544937e-05, 
    1.1334619980672883e-05, 1.0978454194452643e-05, 0.0027295181740957884, 0.0025689582815019182, 0.0021836145392766303, 
    0.0017468916314213045, 0.0013647590870478937, 0.0010651778240373807, 0.00083985174587562718, 0.00067188139670050181, 
    0.00054590363481915756, 0.00045022980191270737, 0.0003764852653925226, 0.00031877584514987308, 0.00027295181740957878, 
    0.00023606643667855464, 0.00020600137162987085, 0.00018121282483623492, 0.00016055989259386989, 0.00014318783864109057, 
    0.00012844791407509593, 0.00011584162012077616, 0.0001049814682344534, 9.5562999530705956e-05, 8.73445815710652e-05, 
    8.0132643643179103e-05, 7.3770761462048343e-05, 6.8131498885386278e-05, 6.3110246799902615e-05, 5.8620524544339074e-05, 
    5.459036348191576e-05, 5.0959499166315762e-05, 4.7677173346651321e-05, 4.4700399985192024e-05, 4.1992587293781363e-05, 
    3.9522435100029511e-05, 3.7263046745334992e-05, 3.5191209335642729e-05, 3.3286807001168147e-05, 3.1532339917352067e-05, 
    2.9912527935296314e-05, 2.8413982293775282e-05, 2.7024932416789982e-05, 2.5734997516518925e-05, 2.4534994823332926e-05, 
    2.3416777901089875e-05, 2.2373099787670393e-05, 2.1397496710207062e-05, 2.0484188923795782e-05, 1.962799585866634e-05, 
    1.8824263269626126e-05, 1.8068800490497561e-05, 1.7357826226364313e-05, 1.6687921584078187e-05, 1.6055989259386988e-05, 
    1.5459217977179684e-05, 1.4895051427534998e-05, 1.4361161060681554e-05, 1.3855422203531922e-05, 1.3375893043042148e-05, 
    1.2920796090394263e-05, 1.2488501797407096e-05, 1.2077514044671629e-05, 1.1686457261314587e-05, 1.1314064970345236e-05, 
    1.095916958231684e-05, 0.0017468916314213045, 0.0016797034917512544, 0.0015059410615700904, 0.0012844791407509591, 
    0.0010651778240373807, 0.00087344581571065217, 0.00071593919320545269, 0.00059016609169638674, 0.00049069989646665865, 
    0.00041200274325974164, 0.0003493783262842608, 0.00029912527935296304, 0.00025841592180788524, 0.00022511490095635366, 
    0.00019761217550014754, 0.00017468916314213043, 0.00015541740493072104, 0.00013908372861634591, 0.00012513550368347453, 
    0.00011314064970345238, 0.00010275833126007672, 9.3717362200713755e-05, 8.5800178360574884e-05, 7.8830849793380171e-05, 
    7.2666041240486881e-05, 6.7188139670050178e-05, 6.2299986855253366e-05, 5.7920810060388072e-05, 5.3983054123031653e-05, 
    5.0429896980984539e-05, 4.7213287335710933e-05, 4.429238416382618e-05, 4.1632307707848063e-05, 3.9203133559724069e-05, 
    3.6979077718486553e-05, 3.4937832628426087e-05, 3.3060023304718105e-05, 3.1328759530511194e-05, 2.9729265340730161e-05, 
    2.8248571012634286e-05, 2.6875255868020073e-05, 2.5599232582375509e-05, 2.441156555926921e-05, 2.3304317388224449e-05, 
    2.2270418554580625e-05, 2.1303556480747613e-05, 2.0398080703191319e-05, 1.9548921569173058e-05, 1.8751520302933707e-05, 
    1.8001768666748808e-05, 1.7295956746745589e-05, 1.6630727641101528e-05, 1.6003038030609236e-05, 1.5410123777534442e-05, 
    1.4849469835271206e-05, 1.4318783864109055e-05, 1.3815973041927432e-05, 1.333912363638748e-05, 1.2886482970059787e-05, 
    1.2456443464213523e-05, 1.204752849256072e-05, 1.1658379814610948e-05, 1.1287746390677853e-05, 1.0934474407995146e-05, 
    0.0012131191884870169, 0.0011803321833927735, 0.0010918072696383151, 0.00097049535078961343, 0.00083985174587562718, 
    0.00071593919320545269, 0.00060655959424350857, 0.00051379165630038371, 0.00043672290785532614, 0.0003732674426113898, 
    0.00032111978518773978, 0.00027816745723269182, 0.00024262383769740336, 0.00021303556480747613, 0.00018824263269626127, 
    0.00016732678461889891, 0.00014956263967648157, 0.00013437627934010036, 0.00012131191884870168, 0.0001100057702406363, 
    0.00010016580455397387, 9.1556165168831468e-05, 8.3985174587562726e-05, 7.729608988589842e-05, 7.1359952263942185e-05, 
    6.6070031445586406e-05, 6.1337487058332318e-05, 5.7087961811153742e-05, 5.3258891201869034e-05, 4.9797366916228755e-05, 
    4.6658430326423725e-05, 4.3803701891206235e-05, 4.120027432597416e-05, 3.8819814031584546e-05, 3.6637827840211931e-05, 
    3.4633061685592871e-05, 3.2787005094243701e-05, 3.1083480986144197e-05, 2.9508304584819334e-05, 2.8048998577734488e-05, 
    2.6694554269885457e-05, 2.5435230509919989e-05, 2.4262383769740337e-05, 2.3168324024155227e-05, 2.214619208191309e-05, 
    2.1189854820733919e-05, 2.029381542078653e-05, 1.9453136207364189e-05, 1.8663372130569494e-05, 1.7920513248064266e-05, 
    1.7220934852339357e-05, 1.6561354109037773e-05, 1.5938792257493654e-05, 1.5350541576637124e-05, 1.4794136444963624e-05, 
    1.4267327927321989e-05, 1.3768061407797166e-05, 1.3294456860131694e-05, 1.2844791407509592e-05, 1.2417483874191815e-05, 
    1.2011081074128884e-05, 1.1624245617655738e-05, 1.1255745047817684e-05, 1.0904442143703522e-05, 0.00089127124052107388, 
    0.00087344581571065217, 0.00082400548651948339, 0.00075297053078504509, 0.00067188139670050181, 0.00059016609169638674, 
    0.00051379165630038371, 0.00044563562026053694, 0.00038648044942949213, 0.0003359406983502509, 0.00029310262272169544, 
    0.0002568958281501918, 0.00022628129940690473, 0.00020033160910794777, 0.00017825424810421472, 0.00015938792257493654, 
    0.00014318783864109057, 0.00012920796090394265, 0.00011708388950544936, 0.00010651778240373807, 9.7265681036820957e-05, 
    8.9127124052107375e-05, 8.193675569518313e-05, 7.5557596514762317e-05, 6.9875665256852187e-05, 6.479568365806025e-05, 
    6.0237642462803606e-05, 5.6134049852869691e-05, 5.2427720030651399e-05, 4.9069989646665853e-05, 4.6019273746609698e-05, 
    4.3239891866863962e-05, 4.0701109772164607e-05, 3.8376353941592811e-05, 3.6242564967246973e-05, 3.427966309696437e-05, 
    3.2470104673258443e-05, 3.0798512542688723e-05, 2.9251366902567059e-05, 2.7816745723269173e-05, 2.648410599486514e-05, 
    2.5244098719961047e-05, 2.4088411905974966e-05, 2.3009636873304849e-05, 2.2001154048127263e-05, 2.1057035094278017e-05, 
    2.0171958792393818e-05, 1.9341138523265107e-05, 1.8560259577361926e-05, 1.7825424810421472e-05, 1.7133107408996709e-05, 
    1.6480109730389663e-05, 1.5863527346724523e-05, 1.5280717559668513e-05, 1.4729271765778285e-05, 1.4206991146887641e-05, 
    1.3711865238785746e-05, 1.3242052997432568e-05, 1.279586603736672e-05, 1.237175376360697e-05, 1.1968290157723379e-05, 
    1.1584162012077613e-05, 1.1218158434506195e-05, 1.0869161469769193e-05, 0.00068237954352394709, 0.00067188139670050181, 
    0.00064223957037547956, 0.00059825055870592629, 0.00054590363481915756, 0.00049069989646665865, 0.00043672290785532614, 
    0.00038648044942949213, 0.00034118977176197344, 0.00030118821231401801, 0.00026629445600934517, 0.00023606643667855464, 
    0.00020996293646890679, 0.00018743472440142754, 0.00016797034917512545, 0.00015111519302952461, 0.00013647590870478939, 
    0.00012371753763606972, 0.00011255745047817684, 0.00010275833126007672, 9.412131634813065e-05, 8.647978373372795e-05, 
    7.9693961287468271e-05, 7.3646358828891411e-05, 6.8237954352394695e-05, 6.3385037424575622e-05, 5.9016609169638661e-05, 
    5.5072245631188665e-05, 5.1500342907467712e-05, 4.8256674901141002e-05, 4.5303206209058729e-05, 4.2607112961495234e-05, 
    4.0139973148467472e-05, 3.7877095217287603e-05, 3.5796959660272642e-05, 3.3880753130746793e-05, 3.2111978518773982e-05, 
    3.0476127554454017e-05, 2.896040503019404e-05, 2.7553495763742974e-05, 2.6245367058613349e-05, 2.5027100736694908e-05, 
    2.3890749882676489e-05, 2.2829216301898909e-05, 2.18361453927663e-05, 2.0905835703940928e-05, 2.0033160910794776e-05, 
    1.9213502325355306e-05, 1.8442690365512086e-05, 1.7716953665530471e-05, 1.703287472134657e-05, 1.6387351139036627e-05, 
    1.5777561699975654e-05, 1.5200936576934429e-05, 1.4655131136084768e-05, 1.4138002844134871e-05, 1.364759087047894e-05, 
    1.318209803366514e-05, 1.2739874791578941e-05, 1.2319405017075489e-05, 1.191929333666283e-05, 1.1538253840299238e-05, 
    1.1175099996298006e-05, 1.0828735627456636e-05, 0.00053916408377200752, 0.00053258891201869035, 0.00051379165630038371, 
    0.00048524767539480672, 0.00045022980191270737, 0.00041200274325974164, 0.0003732674426113898, 0.0003359406983502509, 
    0.00030118821231401801, 0.00026958204188600382, 0.00024128337450570505, 0.00021619945933431986, 0.0001940990701579227, 
    0.00017468916314213043, 0.00015766169958676031, 0.00014271990452788437, 0.00012959136731612053, 0.00011803321833927734, 
    0.00010783281675440149, 9.880608775007377e-05, 9.0794783337905634e-05, 8.3663392309449441e-05, 7.729608988589842e-05, 
    7.1593919320545285e-05, 6.6472284300658478e-05, 6.185876881803486e-05, 5.7691269201496192e-05, 5.3916408377200752e-05, 
    5.04881974399221e-05, 4.7366909745696975e-05, 4.4518135357321729e-05, 4.1911987318169487e-05, 3.9522435100029511e-05, 
    3.7326744261138973e-05, 3.5305004677067586e-05, 3.3439732607605368e-05, 3.1715534339529855e-05, 3.0118821231401796e-05, 
    2.8637567728218113e-05, 2.7261105359258809e-05, 2.5979946927740995e-05, 2.4785636087135418e-05, 2.3670618311941793e-05, 
    2.2628129940690472e-05, 2.1652102521334957e-05, 2.0737080145077215e-05, 1.9878147831375788e-05, 1.9070869338660529e-05, 
    1.8311233033766296e-05, 1.7595604667821361e-05, 1.6920686085057192e-05, 1.6283479040094186e-05, 1.5681253423889629e-05, 
    1.511151930295246e-05, 1.4572002264108309e-05, 1.4060621630886222e-05, 1.3575471179836062e-05, 1.311480203769748e-05, 
    1.2677007484915127e-05, 1.2260609428841272e-05, 1.1864246342171315e-05, 1.148666248961931e-05, 1.1126698289307673e-05, 
    1.0783281675440151e-05, 0.00043672290785532614, 0.00043239891866863971, 0.00041992587293781359, 0.00040066321821589549, 
    0.0003764852653925226, 0.0003493783262842608, 0.00032111978518773978, 0.00029310262272169544, 0.00026629445600934517, 
    0.00024128337450570505, 0.00021836145392766304, 0.00019761217550014754, 0.00017898479830136317, 0.00016235052336629221, 
    0.00014754152292409669, 0.00013437627934010036, 0.00012267497411666466, 0.00011226809970573935, 0.00010300068581493541, 
    9.473381949139395e-05, 8.73445815710652e-05, 8.072512159987543e-05, 7.4781319838240759e-05, 6.9431304905457253e-05, 
    6.4603980451971311e-05, 6.0237642462803606e-05, 5.6278725239088414e-05, 5.268068852295852e-05, 4.9403043875036885e-05, 
    4.641051093042786e-05, 4.3672290785532607e-05, 4.1161442776185319e-05, 3.885435123268026e-05, 3.6730269794392443e-05, 
    3.4770932154086478e-05, 3.2960219460779325e-05, 3.1283875920868632e-05, 2.9729265340730161e-05, 2.8285162425863095e-05, 
    2.6941573587620365e-05, 2.5689582815019181e-05, 2.4521218857682551e-05, 2.3429340550178439e-05, 2.2407537601607292e-05, 
    2.1450044590143721e-05, 2.0551666252015343e-05, 1.9707712448345043e-05, 1.8913941440247991e-05, 1.816651031012172e-05, 
    1.7461931541596408e-05, 1.6797034917512545e-05, 1.6168934019079085e-05, 1.5574996713813342e-05, 1.5012819108123966e-05, 
    1.4480202515097018e-05, 1.3975133051370435e-05, 1.3495763530757913e-05, 1.3040397368030042e-05, 1.2607474245246135e-05, 
    1.2195557326314608e-05, 1.1803321833927733e-05, 1.1429544827409738e-05, 1.1073096040956545e-05, 1.0732929659752426e-05, 
    0.0003609280230209307, 0.00035796959660272629, 0.0003493783262842608, 0.0003359406983502509, 0.00031877584514987308, 
    0.00029912527935296304, 0.00027816745723269182, 0.0002568958281501918, 0.00023606643667855464, 0.00021619945933431986, 
    0.00019761217550014754, 0.00018046401151046535, 0.00016480109730389664, 0.00015059410615700898, 0.00013776747881871483, 
    0.00012622049359980523, 0.00011584162012077616, 0.00010651778240373807, 9.8139979293331706e-05, 9.0606412418117458e-05, 
    8.3823974636339002e-05, 7.7708702465360534e-05, 7.2185604604186135e-05, 6.7188139670050178e-05, 6.2657519061022401e-05, 
    5.8541944752724681e-05, 5.4795847911584213e-05, 5.1379165630038368e-05, 4.8256674901141002e-05, 4.539739166895281e-05, 
    4.2774036028925178e-05, 4.0362560799937728e-05, 3.8141738677321052e-05, 3.6092802302093067e-05, 3.4199131390393597e-05, 
    3.2445981267111903e-05, 3.0820247555068885e-05, 2.931026227216954e-05, 2.7905617115356301e-05, 2.6597010222614255e-05, 
    2.5376113181599427e-05, 2.4235455485867154e-05, 2.3168324024155227e-05, 2.2168675525651071e-05, 2.1231060177701801e-05, 
    2.0350554886082297e-05, 1.9522704866129908e-05, 1.8743472440142749e-05, 1.8009192076508293e-05, 1.7316530842796435e-05, 
    1.6662453561820911e-05, 1.6044192059343355e-05, 1.5459217977179684e-05, 1.4905218698133996e-05, 1.4380075991285024e-05, 
    1.3881847039266564e-05, 1.3408747554661532e-05, 1.2959136731612053e-05, 1.2531503812204481e-05, 1.2124456075939094e-05, 
    1.1736708085335289e-05, 1.1367072042043886e-05, 1.1014449126237733e-05, 1.0677821707954185e-05, 0.00030327979712175423, 
    0.00030118821231401801, 0.00029508304584819337, 0.00028543980905576874, 0.00027295181740957878, 0.00025841592180788524, 
    0.00024262383769740336, 0.00022628129940690473, 0.00020996293646890679, 0.0001940990701579227, 0.00017898479830136317, 
    0.00016480109730389664, 0.00015163989856087714, 0.00013952808557678151, 0.00012844791407509593, 0.00011835309155970896, 
    0.00010918072696383153, 0.00010085979396196905, 9.3316860652847451e-05, 8.647978373372795e-05, 8.0279946296934944e-05, 
    7.4653488522277974e-05, 6.9541864308172955e-05, 6.4891962534223806e-05, 6.065595942435084e-05, 5.6791015325790135e-05, 
    5.3258891201869034e-05, 5.0025533545856364e-05, 4.7060658174065318e-05, 4.4337351051302149e-05, 4.1831696154724727e-05, 
    3.9522435100029511e-05, 3.7390659919120393e-05, 3.5419538349985895e-05, 3.3594069835025089e-05, 3.1900869821426308e-05, 
    3.032797971217542e-05, 2.8864699792156389e-05, 2.7501442560159075e-05, 2.622960407539496e-05, 2.5041451138493468e-05, 
    2.393002234823705e-05, 2.2889041292207867e-05, 2.1912840333935083e-05, 2.0996293646890681e-05, 2.0134758315137213e-05, 
    1.9324022471474605e-05, 1.8560259577361926e-05, 1.7839988065985546e-05, 1.7160035672114973e-05, 1.6517507861396601e-05, 
    1.5909759849010057e-05, 1.5334371764583079e-05, 1.4789126578236576e-05, 1.4271990452788436e-05, 1.378109523052465e-05, 
    1.3314722800467258e-05, 1.2871291124530686e-05, 1.2449341729057189e-05, 1.204752849256072e-05, 1.1664607581605931e-05, 
    1.1299428405053716e-05, 1.0950925472801559e-05, 0.00025841592180788524, 0.0002568958281501918, 0.00025244098719961046, 
    0.00024534994823332927, 0.00023606643667855464, 0.00022511490095635366, 0.00021303556480747613, 0.00020033160910794777, 
    0.00018743472440142754, 0.00017468916314213043, 0.00016235052336629221, 0.00015059410615700898, 0.00013952808557678151, 
    0.00012920796090394265, 0.00011965011174118526, 0.00011084337762825535, 0.00010275833126007672, 9.5354346693302656e-05, 
    8.8584768327652359e-05, 8.2400548651948333e-05, 7.6752707883185608e-05, 7.1593919320545285e-05, 6.6879465215210723e-05, 
    6.2567751841737277e-05, 5.8620524544339074e-05, 5.5002885120318151e-05, 5.168318436157705e-05, 4.8632840518410478e-05, 
    4.582611834788312e-05, 4.3239891866863962e-05, 4.0853405786279331e-05, 3.8648044942949217e-05, 3.6607117171443938e-05, 
    3.4715652452728633e-05, 3.2960219460779325e-05, 3.1328759530511194e-05, 2.9810437396267989e-05, 2.8395507662895061e-05, 
    2.7075195775283701e-05, 2.5841592180788528e-05, 2.4687558386394918e-05, 2.3606643667855463e-05, 2.2593011270322097e-05, 
    2.1641373035447287e-05, 2.0746931489564187e-05, 1.9905328525766914e-05, 1.9112599906141191e-05, 1.8365134897196218e-05, 
    1.7659640430866404e-05, 1.6993109255071053e-05, 1.636279160192305e-05, 1.5766169958676034e-05, 1.5200936576934429e-05, 
    1.4664973400111688e-05, 1.4156334128211544e-05, 1.3673228173303888e-05, 1.3214006289117282e-05, 1.277714768447414e-05, 
    1.236124845330671e-05, 1.1965011174118524e-05, 1.1587235549358612e-05, 1.1226809970573934e-05, 1.0882703908679942e-05, 
    0.00022281781013026847, 0.00022168675525651075, 0.00021836145392766304, 0.00021303556480747613, 0.00020600137162987085, 
    0.00019761217550014754, 0.00018824263269626127, 0.00017825424810421472, 0.00016797034917512545, 0.00015766169958676031, 
    0.00014754152292409669, 0.00013776747881871483, 0.00012844791407509593, 0.00011965011174118526, 0.00011140890506513424, 
    0.00010373465744782091, 9.6620112357373032e-05, 9.0045960382541487e-05, 8.3985174587562726e-05, 7.8406267119448137e-05, 
    7.3275655680423861e-05, 6.8559326193928739e-05, 6.422395703754795e-05, 6.0237642462803606e-05, 5.6570324851726182e-05, 
    5.3194020445228524e-05, 5.0082902276986943e-05, 4.7213287335710933e-05, 4.4563562026053681e-05, 4.2114070188556041e-05, 
    3.9846980643734135e-05, 3.7746145882050649e-05, 3.5796959660272642e-05, 3.3986218510142112e-05, 3.2301990225985662e-05, 
    3.0733491052450816e-05, 2.9270972376362341e-05, 2.7905617115356301e-05, 2.6629445600934517e-05, 2.5435230509919989e-05, 
    2.4316420259205239e-05, 2.3267070210725955e-05, 2.2281781013026844e-05, 2.1355643415908366e-05, 2.0484188923795782e-05, 
    1.9663345693621164e-05, 1.8889399128690579e-05, 1.8158956667581126e-05, 1.7468916314213047e-05, 1.6816438500397617e-05, 
    1.6198920914515062e-05, 1.5613975969085667e-05, 1.5059410615700901e-05, 1.4533208248097372e-05, 1.4033512463217423e-05, 
    1.3558612476104504e-05, 1.310693000766285e-05, 1.2677007484915127e-05, 1.2267497411666463e-05, 1.1877152783664023e-05, 
    1.1504818436652424e-05, 1.114942322837187e-05, 1.080997296671599e-05, 0.0001940990701579227, 0.00019324022471474606, 
    0.00019070869338660528, 0.0001866337213056949, 0.00018121282483623492, 0.00017468916314213043, 0.00016732678461889891, 
    0.00015938792257493654, 0.00015111519302952461, 0.00014271990452788437, 0.00013437627934010036, 0.00012622049359980523, 
    0.00011835309155970896, 0.00011084337762825535, 0.00010373465744782091, 9.7049535078961349e-05, 9.0794783337905634e-05, 
    8.496554627535527e-05, 7.9548799245050301e-05, 7.4526093490669983e-05, 6.9875665256852187e-05, 6.5574010188487403e-05, 
    6.1597025085377447e-05, 5.7920810060388072e-05, 5.4522210718517617e-05, 5.1379165630038368e-05, 4.8470910971734308e-05, 
    4.5778082584415727e-05, 4.3282746070894575e-05, 4.0968377847591558e-05, 3.8819814031584546e-05, 3.6823179414445719e-05, 
    3.4965805272644202e-05, 3.3236142150329232e-05, 3.1623671821529775e-05, 3.0118821231401796e-05, 2.8712880200876143e-05, 
    2.7397923955792103e-05, 2.616674103387215e-05, 2.5012766772928185e-05, 2.393002234823705e-05, 2.291305917394156e-05, 
    2.1956908388905286e-05, 2.1057035094278017e-05, 2.0209296985438506e-05, 1.9409907015792277e-05, 1.8655399737519273e-05, 
    1.7942600980087353e-05, 1.7268600547857891e-05, 1.6630727641101528e-05, 1.6026528728635823e-05, 1.5453747624038438e-05, 
    1.4910307540297918e-05, 1.4394294919424065e-05, 1.3903944853719389e-05, 1.3437627934010035e-05, 1.2993838377129607e-05, 
    1.2571183300383597e-05, 1.2168373024667766e-05, 1.1784212300467516e-05, 1.141759236223075e-05, 1.1067483726693515e-05, 
    1.0732929659752426e-05, 0.00017059488588098677, 0.00016993109255071054, 0.00016797034917512545, 0.00016480109730389664, 
    0.00016055989259386989, 0.00015541740493072104, 0.00014956263967648157, 0.00014318783864109057, 0.00013647590870478939, 
    0.00012959136731612053, 0.00012267497411666466, 0.00011584162012077616, 0.00010918072696383153, 0.00010275833126007672, 
    9.6620112357373032e-05, 9.0794783337905634e-05, 8.5297442940493359e-05, 8.0132643643179103e-05, 7.5297053078504504e-05, 
    7.0781670641057701e-05, 6.6573614002336294e-05, 6.2657519061022401e-05, 5.9016609169638661e-05, 5.5633491446538352e-05, 
    5.2490734117226699e-05, 4.9571272174270836e-05, 4.6858681100356884e-05, 4.4337351051302149e-05, 4.1992587293781363e-05, 
    3.9810657051533827e-05, 3.7778798257381152e-05, 3.5885201960174691e-05, 3.4118977176197348e-05, 3.2470104673258443e-05, 
    3.092938440901743e-05, 2.9488380003735729e-05, 2.8139362619544211e-05, 2.6875255868020073e-05, 2.5689582815019181e-05, 
    2.4576415748752174e-05, 2.3530329087032662e-05, 2.2546355593976571e-05, 2.1619945933431988e-05, 2.0746931489564187e-05, 
    1.9923490321867068e-05, 1.9146116083091899e-05, 1.8411589707222853e-05, 1.7716953665530471e-05, 1.7059488588098674e-05, 
    1.6436692053267825e-05, 1.5846259356143906e-05, 1.5286066078240328e-05, 1.4754152292409665e-05, 1.4248708249765942e-05, 
    1.3768061407797166e-05, 1.3310664670994396e-05, 1.2875085726866928e-05, 1.2459997371050673e-05, 1.2064168725285251e-05, 
    1.1686457261314587e-05, 1.1325801552264682e-05, 1.0981214680797741e-05, 0.00015111519302952461, 0.00015059410615700898, 
    0.00014905218698133997, 0.00014655131136084772, 0.00014318783864109057, 0.00013908372861634591, 0.00013437627934010036, 
    0.00012920796090394265, 0.00012371753763606972, 0.00011803321833927734, 0.00011226809970573935, 0.00010651778240373807, 
    0.00010085979396196905, 9.5354346693302656e-05, 9.0045960382541487e-05, 8.496554627535527e-05, 8.0132643643179103e-05, 
    7.5557596514762317e-05, 7.12435412488297e-05, 6.7188139670050178e-05, 6.3385037424575622e-05, 5.9825055870592621e-05, 
    5.6497142025268571e-05, 5.3389108539770914e-05, 5.04881974399221e-05, 4.7781499765352958e-05, 4.5256259881380943e-05, 
    4.2900089180287449e-05, 4.0701109772164607e-05, 3.8648044942949217e-05, 3.6730269794392443e-05, 3.4937832628426087e-05, 
    3.3261455282203057e-05, 3.1692518712287818e-05, 3.0223038605904926e-05, 2.8845634600748099e-05, 2.7553495763742974e-05, 
    2.6340344261479257e-05, 2.5200398606770115e-05, 2.4128337450570501e-05, 2.3119264576777455e-05, 2.2168675525651071e-05, 
    2.1272426101087488e-05, 2.0426702893139672e-05, 1.962799585866634e-05, 1.8873072941025328e-05, 1.8158956667581126e-05, 
    1.7482902636322104e-05, 1.6842379786167605e-05, 1.6235052336629225e-05, 1.5658763279143999e-05, 1.511151930295246e-05, 
    1.4591477041607952e-05, 1.4096930531159654e-05, 1.3626299777077261e-05, 1.3178120333594629e-05, 1.2751033805994922e-05, 
    1.2343779193197459e-05, 1.1955184994670849e-05, 1.1584162012077613e-05, 1.1229696782086041e-05, 1.0890845582427085e-05, 
    0.00013479102094300188, 0.00013437627934010036, 0.00013314722800467259, 0.00013114802037697481, 0.00012844791407509593, 
    0.00012513550368347453, 0.00012131191884870168, 0.00011708388950544936, 0.00011255745047817684, 0.00010783281675440149, 
    0.00010300068581493541, 9.8139979293331706e-05, 9.3316860652847451e-05, 8.8584768327652359e-05, 8.3985174587562726e-05, 
    7.9548799245050301e-05, 7.5297053078504504e-05, 7.12435412488297e-05, 6.7395510471500954e-05, 6.3755169029974617e-05, 
    6.0320843626426263e-05, 5.7087961811153742e-05, 5.4049864833579964e-05, 5.119846516475101e-05, 4.8524767539480674e-05, 
    4.6019273746609698e-05, 4.3672290785532607e-05, 4.1474160290154423e-05, 3.9415424896690079e-05, 3.7486944880285505e-05, 
    3.5679976131971092e-05, 3.3986218510142112e-05, 3.2397841829030132e-05, 3.090749524807687e-05, 2.9508304584819334e-05, 
    2.8193861062319309e-05, 2.6958204188600373e-05, 2.5795800818389018e-05, 2.4701521937518443e-05, 2.3670618311941793e-05, 
    2.2698695834476409e-05, 2.1781691164854174e-05, 2.091584807736236e-05, 2.0097694793158127e-05, 1.9324022471474605e-05, 
    1.8591864957655429e-05, 1.7898479830136321e-05, 1.7241330748335021e-05, 1.661807107516462e-05, 1.6026528728635823e-05, 
    1.5464692204508715e-05, 1.4930697704455594e-05, 1.4422817300374048e-05, 1.3939448064325759e-05, 1.3479102094300188e-05, 
    1.3040397368030042e-05, 1.2622049359980525e-05, 1.2222863360070701e-05, 1.1841727436424244e-05, 1.147760598831343e-05, 
    1.1129533839330432e-05, 1.0796610824606334e-05, 0.00012097587475216789, 0.00012064168725285254, 0.00011965011174118526, 
    0.00011803321833927734, 0.00011584162012077616, 0.00011314064970345238, 0.0001100057702406363, 0.00010651778240373807, 
    0.00010275833126007672, 9.880608775007377e-05, 9.473381949139395e-05, 9.0606412418117458e-05, 8.647978373372795e-05, 
    8.2400548651948333e-05, 7.8406267119448137e-05, 7.4526093490669983e-05, 7.0781670641057701e-05, 6.7188139670050178e-05, 
    6.3755169029974617e-05, 6.0487937376083945e-05, 5.7388029941567157e-05, 5.4454227912135424e-05, 5.168318436157705e-05, 
    4.9069989646665853e-05, 4.6608634776448898e-05, 4.429238416382618e-05, 4.2114070188556041e-05, 4.0066321821589558e-05, 
    3.8141738677321052e-05, 3.6333020620243447e-05, 3.4633061685592871e-05, 3.3035015722793203e-05, 3.1532339917352067e-05, 
    3.0118821231401796e-05, 2.8788589838848126e-05, 2.7536122815594332e-05, 2.6356240667189265e-05, 2.5244098719961047e-05, 
    2.419517495043358e-05, 2.3205255465213927e-05, 2.2270418554580625e-05, 2.1387018014462589e-05, 2.0551666252015343e-05, 
    1.9761217550014752e-05, 1.9012751756870967e-05, 1.8303558585721966e-05, 1.7631122642524264e-05, 1.6993109255071053e-05, 
    1.6387351139036627e-05, 1.5811835910764884e-05, 1.5264694437445863e-05, 1.4744190001867865e-05, 1.4248708249765942e-05, 
    1.3776747881871487e-05, 1.3326912049292831e-05, 1.2897900409194513e-05, 1.2488501797407096e-05, 1.209758747521679e-05, 
    1.1724104908867815e-05, 1.1367072042043886e-05, 1.1025572023613385e-05, 1.0698748355103531e-05, 0.00010918072696383153, 
    0.00010890845582427083, 0.00010809972966715993, 0.00010677821707954183, 0.0001049814682344534, 0.00010275833126007672, 
    0.00010016580455397387, 9.7265681036820957e-05, 9.412131634813065e-05, 9.0794783337905634e-05, 8.73445815710652e-05, 
    8.3823974636339002e-05, 8.0279946296934944e-05, 7.6752707883185608e-05, 7.3275655680423861e-05, 6.9875665256852187e-05, 
    6.6573614002336294e-05, 6.3385037424575622e-05, 6.0320843626426263e-05, 5.7388029941567157e-05, 5.459036348191576e-05, 
    5.1929002123106551e-05, 4.9403043875036885e-05, 4.7010000845567931e-05, 4.4746199575340793e-05, 4.2607112961495234e-05, 
    4.0587630841573054e-05, 3.8682277046530221e-05, 3.6885380731024172e-05, 3.5191209335642729e-05, 3.3594069835025089e-05, 
    3.2088384118686716e-05, 3.0668743529166166e-05, 2.9329946800223377e-05, 2.8067024926434839e-05, 2.6875255868020073e-05, 
    2.5750171453733852e-05, 2.4687558386394918e-05, 2.3683454872848487e-05, 2.2734144084087766e-05, 2.18361453927663e-05, 
    2.098620412567641e-05, 2.0181280399968857e-05, 1.9418537476893117e-05, 1.869532995956019e-05, 1.8009192076508293e-05, 
    1.7357826226364313e-05, 1.6739091907065007e-05, 1.6150995112992828e-05, 1.5591678252600002e-05, 1.5059410615700901e-05, 
    1.4552579402043524e-05, 1.4069681309772104e-05, 1.3609314672961238e-05, 1.317017213073963e-05, 1.2751033805994922e-05, 
    1.2350760968759221e-05, 1.1968290157723379e-05, 1.1602627732606965e-05, 1.1252844830077971e-05, 1.0918072696383152e-05, 
    9.9030137835674859e-05, 9.880608775007377e-05, 9.8139979293331706e-05, 9.7049535078961349e-05, 9.5562999530705956e-05, 
    9.3717362200713755e-05, 9.1556165168831468e-05, 8.9127124052107375e-05, 8.647978373372795e-05, 8.3663392309449441e-05, 
    8.072512159987543e-05, 7.7708702465360534e-05, 7.4653488522277974e-05, 7.1593919320545285e-05, 6.8559326193928739e-05, 
    6.5574010188487403e-05, 6.2657519061022401e-05, 5.9825055870592621e-05, 5.7087961811153742e-05, 5.4454227912135424e-05, 
    5.1929002123106551e-05, 4.9515068917837436e-05, 4.7213287335710933e-05, 4.5022980191270737e-05, 4.2942272158832454e-05, 
    4.0968377847591558e-05, 3.9097843138346116e-05, 3.7326744261138973e-05, 3.5650849620842951e-05, 3.4065749442693139e-05, 
    3.2566958080188372e-05, 3.1149993427626683e-05, 2.9810437396267989e-05, 2.8543980905576871e-05, 2.7346456346607769e-05, 
    2.62138600153257e-05, 2.5142366600767193e-05, 2.4128337450570501e-05, 2.3168324024155227e-05, 2.2259067678660858e-05, 
    2.1397496710207062e-05, 2.0580721388092656e-05, 1.9806027567134973e-05, 1.9070869338660529e-05, 1.837286107931536e-05, 
    1.7709769174992951e-05, 1.7079503631416746e-05, 1.6480109730389663e-05, 1.5909759849010057e-05, 1.5366745526225408e-05, 
    1.4849469835271206e-05, 1.4356440100438072e-05, 1.3886260981091449e-05, 1.3437627934010035e-05, 1.3009321056161041e-05, 
    1.2600199303385058e-05, 1.2209195075631148e-05, 1.1835309155970897e-05, 1.147760598831343e-05, 1.1135209277290312e-05, 
    1.0807297892980108e-05, 9.0232005755232675e-05, 9.0045960382541487e-05, 8.9492399150681572e-05, 8.8584768327652359e-05, 
    8.73445815710652e-05, 8.5800178360574884e-05, 8.3985174587562726e-05, 8.193675569518313e-05, 7.9693961287468271e-05, 
    7.729608988589842e-05, 7.4781319838240759e-05, 7.2185604604186135e-05, 6.9541864308172955e-05, 6.6879465215210723e-05, 
    6.422395703754795e-05, 6.1597025085377447e-05, 5.9016609169638661e-05, 5.6497142025268571e-05, 5.4049864833579964e-05, 
    5.168318436157705e-05, 4.9403043875036885e-05, 4.7213287335710933e-05, 4.5116002877616337e-05, 4.3111836905757763e-05, 
    4.120027432597416e-05, 3.9379883485601989e-05, 3.7648526539252245e-05, 3.6003537333497623e-05, 3.4441869704678708e-05, 
    3.2960219460779325e-05, 3.1555123399951307e-05, 3.0223038605904926e-05, 2.896040503019404e-05, 2.7763694078533134e-05, 
    2.6629445600934517e-05, 2.555429536894828e-05, 2.4534994823332926e-05, 2.3568424600935031e-05, 2.2651603104529365e-05, 
    2.1781691164854174e-05, 2.095599365908475e-05, 2.0171958792393818e-05, 1.9427175616340133e-05, 1.8719370246692075e-05, 
    1.8046401151046534e-05, 1.7406253800531131e-05, 1.6797034917512545e-05, 1.6216966500383444e-05, 1.56643797652556e-05, 
    1.5137709111103159e-05, 1.463548618818117e-05, 1.4156334128211544e-05, 1.3698961977896053e-05, 1.3262159363963743e-05, 
    1.2844791407509592e-05, 1.2445793897273469e-05, 1.2064168725285251e-05, 1.1698979583587625e-05, 1.1349347917238203e-05, 
    1.1014449126237733e-05, 1.0693509007231295e-05, 8.2556315284560707e-05, 8.2400548651948333e-05, 8.193675569518313e-05, 
    8.1175261683146121e-05, 8.0132643643179103e-05, 7.8830849793380171e-05, 7.729608988589842e-05, 7.5557596514762317e-05, 
    7.3646358828891411e-05, 7.1593919320545285e-05, 6.9431304905457253e-05, 6.7188139670050178e-05, 6.4891962534223806e-05, 
    6.2567751841737277e-05, 6.0237642462803606e-05, 5.7920810060388072e-05, 5.5633491446538352e-05, 5.3389108539770914e-05, 
    5.119846516475101e-05, 4.9069989646665853e-05, 4.7010000845567931e-05, 4.5022980191270737e-05, 4.3111836905757763e-05, 
    4.1278157642280354e-05, 3.9522435100029511e-05, 3.7844272777757886e-05, 3.6242564967246973e-05, 3.4715652452728633e-05, 
    3.3261455282203057e-05, 3.1877584514987308e-05, 3.0561435119337026e-05, 2.931026227216954e-05, 2.8121243261772441e-05, 
    2.6991527061515833e-05, 2.5918273463224103e-05, 2.4898683458114377e-05, 2.393002234823705e-05, 2.3009636873304849e-05, 
    2.2134967453387023e-05, 2.1303556480747613e-05, 2.0513053445529644e-05, 1.9761217550014752e-05, 1.9045918353917405e-05, 
    1.8365134897196218e-05, 1.7716953665530471e-05, 1.7099565695196791e-05, 1.6511263056912141e-05, 1.595043491071315e-05, 
    1.541556328469206e-05, 1.4905218698133996e-05, 1.4418055723186734e-05, 1.3952808557678151e-05, 1.3508286664253823e-05, 
    1.308337051693607e-05, 1.2677007484915127e-05, 1.2288207874376087e-05, 1.1916041142027999e-05, 1.1559632288388728e-05, 
    1.1218158434506195e-05, 1.0890845582427085e-05, 7.5819949280438558e-05, 7.56885455555158e-05, 7.5297053078504504e-05, 
    7.4653488522277974e-05, 7.3770761462048343e-05, 7.2666041240486881e-05, 7.1359952263942185e-05, 6.9875665256852187e-05, 
    6.8237954352394695e-05, 6.6472284300658478e-05, 6.4603980451971311e-05, 6.2657519061022401e-05, 6.065595942435084e-05, 
    5.8620524544339074e-05, 5.6570324851726182e-05, 5.4522210718517617e-05, 5.2490734117226699e-05, 5.04881974399221e-05, 
    4.8524767539480674e-05, 4.6608634776448898e-05, 4.4746199575340793e-05, 4.2942272158832454e-05, 4.120027432597416e-05, 
    3.9522435100029511e-05, 3.7909974640219286e-05, 3.6363272927171197e-05, 3.4882021394195378e-05, 3.3465356923779783e-05, 
    3.2111978518773982e-05, 3.0820247555068885e-05, 2.9588272889927241e-05, 2.8413982293775282e-05, 2.7295181740957884e-05, 
    2.622960407539496e-05, 2.5214948490492263e-05, 2.4248912151878185e-05, 2.3329215163211863e-05, 2.2453619941147876e-05, 
    2.1619945933431988e-05, 2.0826080489047502e-05, 2.0069986574233736e-05, 1.9349707924471692e-05, 1.8663372130569494e-05, 
    1.8009192076508293e-05, 1.7385466077043239e-05, 1.6790577003280513e-05, 1.6222990633555952e-05, 1.5681253423889629e-05, 
    1.516398985608771e-05, 1.4669899491277326e-05, 1.4197753831447534e-05, 1.3746393070674412e-05, 1.3314722800467258e-05, 
    1.290171071950742e-05, 1.2506383386464091e-05, 1.2127823045135412e-05, 1.176516454351633e-05, 1.141759236223075e-05, 
    1.1084337762825537e-05, 1.0764676062492632e-05, 6.9875665256852187e-05, 6.9764042788390742e-05, 6.9431304905457253e-05, 
    6.8883739409357429e-05, 6.8131498885386278e-05, 6.7188139670050178e-05, 6.6070031445586406e-05, 6.479568365806025e-05, 
    6.3385037424575622e-05, 6.185876881803486e-05, 6.0237642462803606e-05, 5.8541944752724681e-05, 5.6791015325790135e-05, 
    5.5002885120318151e-05, 5.3194020445228524e-05, 5.1379165630038368e-05, 4.9571272174270836e-05, 4.7781499765352958e-05, 
    4.6019273746609698e-05, 4.429238416382618e-05, 4.2607112961495234e-05, 4.0968377847591558e-05, 3.9379883485601989e-05, 
    3.7844272777757886e-05, 3.6363272927171197e-05, 3.4937832628426087e-05, 3.3568248105713005e-05, 3.2254276798768551e-05, 
    3.099523831478539e-05, 2.9790102855069999e-05, 2.8637567728218113e-05, 2.7536122815594332e-05, 2.648410599486514e-05, 
    2.5479749583157881e-05, 2.4521218857682551e-05, 2.3606643667855463e-05, 2.2734144084087766e-05, 2.1901850945603111e-05, 
    2.1107922080972743e-05, 2.0350554886082297e-05, 1.962799585866634e-05, 1.8938547608643802e-05, 1.8280573790511765e-05, 
    1.7652502338533797e-05, 1.7052827327423904e-05, 1.6480109730389663e-05, 1.5932977302273846e-05, 1.5410123777534442e-05, 
    1.4910307540297918e-05, 1.4432349896078191e-05, 1.3975133051370435e-05, 1.353759788764185e-05, 1.3118741599739446e-05, 
    1.271761525495999e-05, 1.2333321317574871e-05, 1.1965011174118524e-05, 1.1611882686927043e-05, 1.1273177796988285e-05, 
    1.0948180191910906e-05, 6.4603980451971311e-05, 6.4508553597537088e-05, 6.422395703754795e-05, 6.3755169029974617e-05, 
    6.3110246799902615e-05, 6.2299986855253366e-05, 6.1337487058332318e-05, 6.0237642462803606e-05, 5.9016609169638661e-05, 
    5.7691269201496192e-05, 5.6278725239088414e-05, 5.4795847911584213e-05, 5.3258891201869034e-05, 5.168318436157705e-05, 
    5.0082902276986943e-05, 4.8470910971734308e-05, 4.6858681100356884e-05, 4.5256259881380943e-05, 4.3672290785532607e-05, 
    4.2114070188556041e-05, 4.0587630841573054e-05, 3.9097843138346116e-05, 3.7648526539252245e-05, 3.6242564967246973e-05, 
    3.4882021394195378e-05, 3.3568248105713005e-05, 3.2301990225985662e-05, 3.1083480986144197e-05, 2.9912527935296314e-05, 
    2.8788589838848126e-05, 2.7710844407063837e-05, 2.667824727277496e-05, 2.5689582815019181e-05, 2.474350752721394e-05, 
    2.3838586673325664e-05, 2.2973324979238624e-05, 2.214619208191309e-05, 2.1355643415908366e-05, 2.0600137162987083e-05, 
    1.9878147831375788e-05, 1.9188176970796402e-05, 1.8528761470315069e-05, 1.7898479830136321e-05, 1.7295956746745589e-05, 
    1.6719866303802681e-05, 1.6168934019079085e-05, 1.5641937960434319e-05, 1.5137709111103159e-05, 1.4655131136084768e-05, 
    1.4193139676806177e-05, 1.3750721280079538e-05, 1.3326912049292831e-05, 1.2920796090394263e-05, 1.2531503812204481e-05, 
    1.215821012960262e-05, 1.1800132608898301e-05, 1.145652958697078e-05, 1.1126698289307673e-05, 1.080997296671599e-05, 
    5.9907120419111959e-05, 5.9825055870592621e-05, 5.9580205710139991e-05, 5.9176545779854495e-05, 5.8620524544339074e-05, 
    5.7920810060388072e-05, 5.7087961811153742e-05, 5.6134049852869691e-05, 5.5072245631188665e-05, 5.3916408377200752e-05, 
    5.268068852295852e-05, 5.1379165630038368e-05, 5.0025533545856364e-05, 4.8632840518410478e-05, 4.7213287335710933e-05, 
    4.5778082584415727e-05, 4.4337351051302149e-05, 4.2900089180287449e-05, 4.1474160290154423e-05, 4.0066321821589558e-05, 
    3.8682277046530221e-05, 3.7326744261138973e-05, 3.6003537333497623e-05, 3.4715652452728633e-05, 3.3465356923779783e-05, 
    3.2254276798768551e-05, 3.1083480986144197e-05, 2.9953560209555966e-05, 2.8864699792156389e-05, 2.7816745723269173e-05, 
    2.680926383396723e-05, 2.5841592180788528e-05, 2.4912886928427052e-05, 2.4022162148257767e-05, 2.3168324024155227e-05, 
    2.2350199992596015e-05, 2.1566563350880302e-05, 2.0816153853924032e-05, 2.0097694793158127e-05, 1.9409907015792277e-05, 
    1.8751520302933707e-05, 1.812128248362349e-05, 1.7517966620751149e-05, 1.6940376565373393e-05, 1.6387351139036627e-05, 
    1.5857767169764927e-05, 1.5350541576637124e-05, 1.4864632670365084e-05, 1.439904081290228e-05, 1.3952808557678151e-05, 
    1.3525020373345501e-05, 1.311480203769748e-05, 1.2721319774405071e-05, 1.2343779193197459e-05, 1.198142408382239e-05, 
    1.1633535105362978e-05, 1.1299428405053716e-05, 1.0978454194452643e-05, 1.0669995305529591e-05, 5.5704452532567118e-05, 
    5.5633491446538352e-05, 5.5421688814127688e-05, 5.5072245631188665e-05, 5.459036348191576e-05, 5.3983054123031653e-05, 
    5.3258891201869034e-05, 5.2427720030651399e-05, 5.1500342907467712e-05, 5.04881974399221e-05, 4.9403043875036885e-05, 
    4.8256674901141002e-05, 4.7060658174065318e-05, 4.582611834788312e-05, 4.4563562026053681e-05, 4.3282746070894575e-05, 
    4.1992587293781363e-05, 4.0701109772164607e-05, 3.9415424896690079e-05, 3.8141738677321052e-05, 3.6885380731024172e-05, 
    3.5650849620842951e-05, 3.4441869704678708e-05, 3.3261455282203057e-05, 3.2111978518773982e-05, 3.099523831478539e-05, 
    2.9912527935296314e-05, 2.8864699792156389e-05, 2.7852226266283559e-05, 2.6875255868020073e-05, 2.5933664361955227e-05, 
    2.5027100736694908e-05, 2.4155028089343258e-05, 2.3316759629221896e-05, 2.2511490095635372e-05, 2.1738322939538379e-05, 
    2.0996293646890681e-05, 2.0284389589193037e-05, 1.9601566779862034e-05, 1.8946763898278791e-05, 1.8318913920105965e-05, 
    1.7716953665530471e-05, 1.7139831548482185e-05, 1.6586513781060617e-05, 1.6055989259386988e-05, 1.554727333055629e-05, 
    1.5059410615700901e-05, 1.4591477041607952e-05, 1.4142581212931546e-05, 1.3711865238785746e-05, 1.3298505111307131e-05, 
    1.290171071950742e-05, 1.2520725569246736e-05, 1.2154826269282664e-05, 1.1803321833927733e-05, 1.146555284471846e-05, 
    1.114089050651342e-05, 1.0828735627456636e-05, 5.1929002123106551e-05, 5.1867328723910461e-05, 5.168318436157705e-05, 
    5.1379165630038368e-05, 5.0959499166315762e-05, 5.0429896980984539e-05, 4.9797366916228755e-05, 4.9069989646665853e-05, 
    4.8256674901141002e-05, 4.7366909745696975e-05, 4.641051093042786e-05, 4.539739166895281e-05, 4.4337351051302149e-05, 
    4.3239891866863962e-05, 4.2114070188556041e-05, 4.0968377847591558e-05, 3.9810657051533827e-05, 3.8648044942949217e-05, 
    3.7486944880285505e-05, 3.6333020620243447e-05, 3.5191209335642729e-05, 3.4065749442693139e-05, 3.2960219460779325e-05, 
    3.1877584514987308e-05, 3.0820247555068885e-05, 2.9790102855069999e-05, 2.8788589838848126e-05, 2.7816745723269173e-05, 
    2.6875255868020073e-05, 2.5964501061553279e-05, 2.508460125533177e-05, 2.4235455485867154e-05, 2.3416777901089875e-05, 
    2.2628129940690472e-05, 2.1868948815990292e-05, 2.1138572500257799e-05, 2.0436261481297429e-05, 1.9761217550014752e-05, 
    1.9112599906141191e-05, 1.848953885924327e-05, 1.7891147392680304e-05, 1.7316530842796435e-05, 1.6764794927267797e-05, 
    1.6235052336629225e-05, 1.5726428082654883e-05, 1.5238063777227012e-05, 1.4769120996121951e-05, 1.4318783864109055e-05, 
    1.3886260981091449e-05, 1.3470786793810186e-05, 1.3071622503900814e-05, 1.2688056590799713e-05, 1.2319405017075489e-05, 
    1.1965011174118524e-05, 1.1624245617655738e-05, 1.1296505635161048e-05, 1.0981214680797741e-05, 1.0677821707954185e-05, 
    4.8524767539480674e-05, 4.8470910971734308e-05, 4.8310056178686516e-05, 4.8044324296515534e-05, 4.7677173346651321e-05, 
    4.7213287335710933e-05, 4.6658430326423725e-05, 4.6019273746609698e-05, 4.5303206209058729e-05, 4.4518135357321729e-05, 
    4.3672290785532607e-05, 4.2774036028925178e-05, 4.1831696154724727e-05, 4.0853405786279331e-05, 3.9846980643734135e-05, 
    3.8819814031584546e-05, 3.7778798257381152e-05, 3.6730269794392443e-05, 3.5679976131971092e-05, 3.4633061685592871e-05, 
    3.3594069835025089e-05, 3.2566958080188372e-05, 3.1555123399951307e-05, 3.0561435119337026e-05, 2.9588272889927241e-05, 
    2.8637567728218113e-05, 2.7710844407063837e-05, 2.680926383396723e-05, 2.5933664361955227e-05, 2.508460125533177e-05, 
    2.4262383769740337e-05, 2.3467109503241595e-05, 2.2698695834476409e-05, 2.1956908388905286e-05, 2.1241386568838818e-05, 
    2.0551666252015343e-05, 1.9887199811262575e-05, 1.9247373638401327e-05, 1.8631523372667496e-05, 1.8038947040699136e-05, 
    1.7468916314213047e-05, 1.6920686085057192e-05, 1.6393502547121851e-05, 1.5886609961998042e-05, 1.5399256271344362e-05, 
    1.4930697704455594e-05, 1.4480202515097018e-05, 1.4047053967684984e-05, 1.3630552679629404e-05, 1.3230018414278281e-05, 
    1.2844791407509592e-05, 1.2474233300637707e-05, 1.2117727742933577e-05, 1.1774680718666113e-05, 1.1444520646103932e-05, 
    1.1126698289307673e-05, 1.0820686517723644e-05, 4.5444631410543818e-05, 4.539739166895281e-05, 4.5256259881380943e-05, 
    4.5022980191270737e-05, 4.4700399985192024e-05, 4.429238416382618e-05, 4.3803701891206235e-05, 4.3239891866863962e-05, 
    4.2607112961495234e-05, 4.1911987318169487e-05, 4.1161442776185319e-05, 4.0362560799937728e-05, 3.9522435100029511e-05, 
    3.8648044942949217e-05, 3.7746145882050649e-05, 3.6823179414445719e-05, 3.5885201960174691e-05, 3.4937832628426087e-05, 
    3.3986218510142112e-05, 3.3035015722793203e-05, 3.2088384118686716e-05, 3.1149993427626683e-05, 3.0223038605904926e-05, 
    2.931026227216954e-05, 2.8413982293775282e-05, 2.7536122815594332e-05, 2.667824727277496e-05, 2.5841592180788528e-05, 
    2.5027100736694908e-05, 2.4235455485867154e-05, 2.3467109503241595e-05, 2.2722315705271909e-05, 2.2001154048127263e-05, 
    2.1303556480747613e-05, 2.0629329610549177e-05, 1.9978175107745931e-05, 1.9349707924471692e-05, 1.8743472440142749e-05, 
    1.8158956667581126e-05, 1.7595604667821361e-05, 1.7052827327423904e-05, 1.6530011652359049e-05, 1.6026528728635823e-05, 
    1.5541740493072105e-05, 1.5075005448923926e-05, 1.462568345128353e-05, 1.4193139676806177e-05, 1.3776747881871487e-05, 
    1.3375893043042148e-05, 1.2989973463870498e-05, 1.2618402422864087e-05, 1.2260609428841272e-05, 1.1916041142027999e-05, 
    1.1584162012077613e-05, 1.1264454677723139e-05, 1.0956420166967538e-05, 4.2648721470246693e-05, 4.2607112961495234e-05, 
    4.2482773137677635e-05, 4.2277145000515604e-05, 4.1992587293781363e-05, 4.1632307707848063e-05, 4.120027432597416e-05, 
    4.0701109772164607e-05, 4.0139973148467472e-05, 3.9522435100029511e-05, 3.885435123268026e-05, 3.8141738677321052e-05, 
    3.7390659919120393e-05, 3.6607117171443938e-05, 3.5796959660272642e-05, 3.4965805272644202e-05, 3.4118977176197348e-05, 
    3.3261455282203057e-05, 3.2397841829030132e-05, 3.1532339917352067e-05, 3.0668743529166166e-05, 2.9810437396267989e-05, 
    2.896040503019404e-05, 2.8121243261772441e-05, 2.7295181740957884e-05, 2.648410599486514e-05, 2.5689582815019181e-05, 
    2.4912886928427052e-05, 2.4155028089343258e-05, 2.3416777901089875e-05, 2.2698695834476409e-05, 2.2001154048127263e-05, 
    2.132436073512334e-05, 2.0668381819939712e-05, 2.0033160910794776e-05, 1.9418537476893117e-05, 1.8824263269626126e-05, 
    1.8250017043682661e-05, 1.7695417660264425e-05, 1.7160035672114973e-05, 1.6643403500584073e-05, 1.6145024319975089e-05, 
    1.56643797652556e-05, 1.5200936576934429e-05, 1.4754152292409665e-05, 1.4323480087088427e-05, 1.3908372861634588e-05, 
    1.3508286664253823e-05, 1.3122683529306675e-05, 1.2751033805994922e-05, 1.2392818043567709e-05, 1.204752849256072e-05, 
    1.1714670275089221e-05, 1.1393762271206004e-05, 1.1084337762825537e-05, 1.0785944871704769e-05, 4.0103113668992295e-05, 
    4.0066321821589558e-05, 3.9956350215491862e-05, 3.9774399622525144e-05, 3.9522435100029511e-05, 3.9203133559724069e-05, 
    3.8819814031584546e-05, 3.8376353941592811e-05, 3.7877095217287603e-05, 3.7326744261138973e-05, 3.6730269794392443e-05, 
    3.6092802302093067e-05, 3.5419538349985895e-05, 3.4715652452728633e-05, 3.3986218510142112e-05, 3.3236142150329232e-05, 
    3.2470104673258443e-05, 3.1692518712287818e-05, 3.090749524807687e-05, 3.0118821231401796e-05, 2.9329946800223377e-05, 
    2.8543980905576871e-05, 2.7763694078533134e-05, 2.6991527061515833e-05, 2.622960407539496e-05, 2.5479749583157881e-05, 
    2.474350752721394e-05, 2.4022162148257767e-05, 2.3316759629221896e-05, 2.2628129940690472e-05, 2.1956908388905286e-05, 
    2.1303556480747613e-05, 2.0668381819939712e-05, 2.0051556834496147e-05, 1.9453136207364189e-05, 1.8873072941025328e-05, 
    1.8311233033766296e-05, 1.7767408781746382e-05, 1.7241330748335021e-05, 1.6732678461889888e-05, 1.6241089916523838e-05, 
    1.5766169958676034e-05, 1.5307497646523874e-05, 1.4864632670365084e-05, 1.4437120920837226e-05, 1.4024499288867249e-05, 
    1.3626299777077261e-05, 1.3242052997432568e-05, 1.2871291124530686e-05, 1.2513550368347452e-05, 1.2168373024667766e-05, 
    1.1835309155970897e-05, 1.1513917950311788e-05, 1.1203768800803646e-05, 1.0904442143703522e-05, 3.7778798257381152e-05, 
    3.7746145882050649e-05, 3.7648526539252245e-05, 3.7486944880285505e-05, 3.7263046745334992e-05, 3.6979077718486553e-05, 
    3.6637827840211931e-05, 3.6242564967246973e-05, 3.5796959660272642e-05, 3.5305004677067586e-05, 3.4770932154086478e-05, 
    3.4199131390393597e-05, 3.3594069835025089e-05, 3.2960219460779325e-05, 3.2301990225985662e-05, 3.1623671821529775e-05, 
    3.092938440901743e-05, 3.0223038605904926e-05, 2.9508304584819334e-05, 2.8788589838848126e-05, 2.8067024926434839e-05, 
    2.7346456346607769e-05, 2.6629445600934517e-05, 2.5918273463224103e-05, 2.5214948490492263e-05, 2.4521218857682551e-05, 
    2.3838586673325664e-05, 2.3168324024155227e-05, 2.2511490095635372e-05, 2.1868948815990292e-05, 2.1241386568838818e-05, 
    2.0629329610549177e-05, 2.0033160910794776e-05, 1.9453136207364189e-05, 1.8889399128690579e-05, 1.8341995290017899e-05, 
    1.7810885312207425e-05, 1.7295956746745589e-05, 1.6797034917512545e-05, 1.6313892710322231e-05, 1.5846259356143906e-05, 
    1.539382826419902e-05, 1.4956263967648155e-05, 1.4533208248097372e-05, 1.4124285506317143e-05, 1.3729107445939204e-05, 
    1.3347277134942729e-05, 1.2978392506844757e-05, 1.2622049359980525e-05, 1.2277843909342877e-05, 1.1945374941338239e-05, 
    1.1624245617655738e-05, 1.1314064970345236e-05, 1.1014449126237733e-05, 1.0725022295071862e-05, 3.5650849620842951e-05, 
    3.5621770624414843e-05, 3.5534817563492765e-05, 3.539083532052885e-05, 3.5191209335642729e-05, 3.4937832628426087e-05, 
    3.4633061685592871e-05, 3.427966309696437e-05, 3.3880753130746793e-05, 3.3439732607605368e-05, 3.2960219460779325e-05, 
    3.2445981267111903e-05, 3.1900869821426308e-05, 3.1328759530511194e-05, 3.0733491052450816e-05, 3.0118821231401796e-05, 
    2.9488380003735729e-05, 2.8845634600748099e-05, 2.8193861062319309e-05, 2.7536122815594332e-05, 2.6875255868020073e-05, 
    2.62138600153257e-05, 2.555429536894828e-05, 2.4898683458114377e-05, 2.4248912151878185e-05, 2.3606643667855463e-05, 
    2.2973324979238624e-05, 2.2350199992596015e-05, 2.1738322939538379e-05, 2.1138572500257799e-05, 2.0551666252015343e-05, 
    1.9978175107745931e-05, 1.9418537476893117e-05, 1.8873072941025328e-05, 1.8341995290017899e-05, 1.7825424810421472e-05, 
    1.732339975626046e-05, 1.6835886964353361e-05, 1.636279160192305e-05, 1.59039660544547e-05, 1.5459217977179684e-05, 
    1.5028317544918312e-05, 1.4611003942968423e-05, 1.4206991146887641e-05, 1.3815973041927432e-05, 1.3437627934010035e-05, 
    1.3071622503900814e-05, 1.271761525495999e-05, 1.2375259502842903e-05, 1.2044205952987481e-05, 1.1724104908867815e-05, 
    1.1414608150949455e-05, 1.111537052316941e-05, 1.082605126066748e-05, 3.369775523575047e-05, 3.3671773928706715e-05, 
    3.3594069835025089e-05, 3.3465356923779783e-05, 3.3286807001168147e-05, 3.3060023304718105e-05, 3.2787005094243701e-05, 
    3.2470104673258443e-05, 3.2111978518773982e-05, 3.1715534339529855e-05, 3.1283875920868632e-05, 3.0820247555068885e-05, 
    3.032797971217542e-05, 2.9810437396267989e-05, 2.9270972376362341e-05, 2.8712880200876143e-05, 2.8139362619544211e-05, 
    2.7553495763742974e-05, 2.6958204188600373e-05, 2.6356240667189265e-05, 2.5750171453733852e-05, 2.5142366600767193e-05, 
    2.4534994823332926e-05, 2.393002234823705e-05, 2.3329215163211863e-05, 2.2734144084087766e-05, 2.214619208191309e-05, 
    2.1566563350880302e-05, 2.0996293646890681e-05, 2.0436261481297429e-05, 1.9887199811262575e-05, 1.9349707924471692e-05, 
    1.8824263269626126e-05, 1.8311233033766296e-05, 1.7810885312207425e-05, 1.732339975626046e-05, 1.6848877617875238e-05, 
    1.6387351139036627e-05, 1.5938792257493654e-05, 1.5503120619642394e-05, 1.5080210906606566e-05, 1.4669899491277326e-05, 
    1.4271990452788436e-05, 1.3886260981091449e-05, 1.3512466208394991e-05, 1.315034350663433e-05, 1.2799616291187753e-05, 
    1.2459997371050673e-05, 1.2131191884870169e-05, 1.1812899860841926e-05, 1.1504818436652424e-05, 1.1206643773552122e-05, 
    1.0918072696383152e-05, 3.1900869821426308e-05, 3.1877584514987308e-05, 3.1807932108909401e-05, 3.1692518712287818e-05, 
    3.1532339917352067e-05, 3.1328759530511194e-05, 3.1083480986144197e-05, 3.0798512542688723e-05, 3.0476127554454017e-05, 
    3.0118821231401796e-05, 2.9729265340730161e-05, 2.931026227216954e-05, 2.8864699792156389e-05, 2.8395507662895061e-05, 
    2.7905617115356301e-05, 2.7397923955792103e-05, 2.6875255868020073e-05, 2.6340344261479257e-05, 2.5795800818389018e-05, 
    2.5244098719961047e-05, 2.4687558386394918e-05, 2.4128337450570501e-05, 2.3568424600935031e-05, 2.3009636873304849e-05, 
    2.2453619941147876e-05, 2.1901850945603111e-05, 2.1355643415908366e-05, 2.0816153853924032e-05, 2.0284389589193037e-05, 
    1.9761217550014752e-05, 1.9247373638401327e-05, 1.8743472440142749e-05, 1.8250017043682661e-05, 1.7767408781746382e-05, 
    1.7295956746745589e-05, 1.6835886964353361e-05, 1.6387351139036627e-05, 1.595043491071315e-05, 1.5525165583196803e-05, 
    1.511151930295246e-05, 1.4709427681216777e-05, 1.4318783864109055e-05, 1.3939448064325759e-05, 1.3571252574745994e-05, 
    1.3214006289117282e-05, 1.2867498758259462e-05, 1.2531503812204481e-05, 1.2205782779634605e-05, 1.1890087336110159e-05, 
    1.1584162012077613e-05, 1.1287746390677853e-05, 1.1000577024063628e-05, 1.0722389095392245e-05, 3.0243968688041973e-05, 
    3.0223038605904926e-05, 3.0160421813213135e-05, 3.0056635089836621e-05, 2.9912527935296314e-05, 2.9729265340730161e-05, 
    2.9508304584819334e-05, 2.9251366902567059e-05, 2.896040503019404e-05, 2.8637567728218113e-05, 2.8285162425863095e-05, 
    2.7905617115356301e-05, 2.7501442560159075e-05, 2.7075195775283701e-05, 2.6629445600934517e-05, 2.616674103387215e-05, 
    2.5689582815019181e-05, 2.5200398606770115e-05, 2.4701521937518443e-05, 2.419517495043358e-05, 2.3683454872848487e-05, 
    2.3168324024155227e-05, 2.2651603104529365e-05, 2.2134967453387023e-05, 2.1619945933431988e-05, 2.1107922080972743e-05, 
    2.0600137162987083e-05, 2.0097694793158127e-05, 1.9601566779862034e-05, 1.9112599906141191e-05, 1.8631523372667496e-05, 
    1.8158956667581126e-05, 1.7695417660264425e-05, 1.7241330748335021e-05, 1.6797034917512545e-05, 1.636279160192305e-05, 
    1.5938792257493654e-05, 1.5525165583196803e-05, 1.5121984344020986e-05, 1.4729271765778285e-05, 1.4347007485391789e-05, 
    1.3975133051370435e-05, 1.3613556978033856e-05, 1.3262159363963743e-05, 1.2920796090394263e-05, 1.2589302619063883e-05, 
    1.2267497411666463e-05, 1.1955184994670849e-05, 1.1652158694112225e-05, 1.1358203065158025e-05, 1.1073096040956545e-05, 
    1.0796610824606334e-05, 2.8712880200876143e-05, 2.8694014970783585e-05, 2.8637567728218113e-05, 2.8543980905576871e-05, 
    2.8413982293775282e-05, 2.8248571012634286e-05, 2.8048998577734488e-05, 2.7816745723269173e-05, 2.7553495763742974e-05, 
    2.7261105359258809e-05, 2.6941573587620365e-05, 2.6597010222614255e-05, 2.622960407539496e-05, 2.5841592180788528e-05, 
    2.5435230509919989e-05, 2.5012766772928185e-05, 2.4576415748752174e-05, 2.4128337450570501e-05, 2.3670618311941793e-05, 
    2.3205255465213927e-05, 2.2734144084087766e-05, 2.2259067678660858e-05, 2.1781691164854174e-05, 2.1303556480747613e-05, 
    2.0826080489047502e-05, 2.0350554886082297e-05, 1.9878147831375788e-05, 1.9409907015792277e-05, 1.8946763898278791e-05, 
    1.848953885924327e-05, 1.8038947040699136e-05, 1.7595604667821361e-05, 1.7160035672114973e-05, 1.6732678461889888e-05, 
    1.6313892710322231e-05, 1.59039660544547e-05, 1.5503120619642394e-05, 1.511151930295246e-05, 1.4729271765778285e-05, 
    1.4356440100438072e-05, 1.399304414787972e-05, 1.3639066453945224e-05, 1.3294456860131694e-05, 1.2959136731612053e-05, 
    1.2633002830642929e-05, 1.2315930847583931e-05, 1.2007778604765636e-05, 1.1708388950544937e-05, 1.141759236223075e-05, 
    1.1135209277290312e-05, 1.0861052172477647e-05, 2.7295181740957884e-05, 2.7278132907890455e-05, 2.7227113956067708e-05, 
    2.7142505149491988e-05, 2.7024932416789982e-05, 2.6875255868020073e-05, 2.6694554269885457e-05, 2.648410599486514e-05, 
    2.6245367058613349e-05, 2.5979946927740995e-05, 2.5689582815019181e-05, 2.5376113181599427e-05, 2.5041451138493468e-05, 
    2.4687558386394918e-05, 2.4316420259205239e-05, 2.393002234823705e-05, 2.3530329087032662e-05, 2.3119264576777455e-05, 
    2.2698695834476409e-05, 2.2270418554580625e-05, 2.18361453927663e-05, 2.1397496710207062e-05, 2.095599365908475e-05, 
    2.0513053445529644e-05, 2.0069986574233736e-05, 1.962799585866634e-05, 1.9188176970796402e-05, 1.8751520302933707e-05, 
    1.8318913920105965e-05, 1.7891147392680304e-05, 1.7468916314213047e-05, 1.7052827327423904e-05, 1.6643403500584073e-05, 
    1.6241089916523838e-05, 1.5846259356143906e-05, 1.5459217977179684e-05, 1.5080210906606566e-05, 1.4709427681216777e-05, 
    1.4347007485391789e-05, 1.399304414787972e-05, 1.364759087047894e-05, 1.3310664670994396e-05, 1.2982250530776638e-05, 
    1.2662305243703281e-05, 1.2350760968759221e-05, 1.204752849256072e-05, 1.1752500211391983e-05, 1.146555284471846e-05, 
    1.1186549893835198e-05, 1.0915343860418048e-05, 2.5979946927740995e-05, 2.5964501061553279e-05, 2.5918273463224103e-05, 
    2.5841592180788528e-05, 2.5734997516518925e-05, 2.5599232582375509e-05, 2.5435230509919989e-05, 2.5244098719961047e-05, 
    2.5027100736694908e-05, 2.4785636087135418e-05, 2.4521218857682551e-05, 2.4235455485867154e-05, 2.393002234823705e-05, 
    2.3606643667855463e-05, 2.3267070210725955e-05, 2.291305917394156e-05, 2.2546355593976571e-05, 2.2168675525651071e-05, 
    2.1781691164854174e-05, 2.1387018014462589e-05, 2.098620412567641e-05, 2.0580721388092656e-05, 2.0171958792393818e-05, 
    1.9761217550014752e-05, 1.9349707924471692e-05, 1.8938547608643802e-05, 1.8528761470315069e-05, 1.812128248362349e-05, 
    1.7716953665530471e-05, 1.7316530842796435e-05, 1.6920686085057192e-05, 1.6530011652359049e-05, 1.6145024319975089e-05, 
    1.5766169958676034e-05, 1.539382826419902e-05, 1.5028317544918312e-05, 1.4669899491277326e-05, 1.4318783864109055e-05, 
    1.3975133051370435e-05, 1.3639066453945224e-05, 1.3310664670994396e-05, 1.2989973463870498e-05, 1.2677007484915127e-05, 
    1.237175376360697e-05, 1.2074174947617531e-05, 1.1784212300467516e-05, 1.1501788460767085e-05, 1.1226809970573934e-05, 
    1.095916958231684e-05, 1.0698748355103531e-05, 2.4757534458918715e-05, 2.474350752721394e-05, 2.4701521937518443e-05, 
    2.4631861695167862e-05, 2.4534994823332926e-05, 2.441156555926921e-05, 2.4262383769740337e-05, 2.4088411905974966e-05, 
    2.3890749882676489e-05, 2.3670618311941793e-05, 2.3429340550178439e-05, 2.3168324024155227e-05, 2.2889041292207867e-05, 
    2.2593011270322097e-05, 2.2281781013026844e-05, 2.1956908388905286e-05, 2.1619945933431988e-05, 2.1272426101087488e-05, 
    2.091584807736236e-05, 2.0551666252015343e-05, 2.0181280399968857e-05, 1.9806027567134973e-05, 1.9427175616340133e-05, 
    1.9045918353917405e-05, 1.8663372130569494e-05, 1.8280573790511765e-05, 1.7898479830136321e-05, 1.7517966620751149e-05, 
    1.7139831548482185e-05, 1.6764794927267797e-05, 1.6393502547121851e-05, 1.6026528728635823e-05, 1.56643797652556e-05, 
    1.5307497646523874e-05, 1.4956263967648155e-05, 1.4611003942968423e-05, 1.4271990452788436e-05, 1.3939448064325759e-05, 
    1.3613556978033856e-05, 1.3294456860131694e-05, 1.2982250530776638e-05, 1.2677007484915127e-05, 1.2378767229459359e-05, 
    1.208754242610922e-05, 1.1803321833927733e-05, 1.1526073049757881e-05, 1.1255745047817684e-05, 1.0992270522409414e-05, 
    1.0735568039708114e-05, 2.361941091699979e-05, 2.3606643667855463e-05, 2.3568424600935031e-05, 2.3505000422783966e-05, 
    2.3416777901089875e-05, 2.3304317388224449e-05, 2.3168324024155227e-05, 2.3009636873304849e-05, 2.2829216301898909e-05, 
    2.2628129940690472e-05, 2.2407537601607292e-05, 2.2168675525651071e-05, 2.1912840333935083e-05, 2.1641373035447287e-05, 
    2.1355643415908366e-05, 2.1057035094278017e-05, 2.0746931489564187e-05, 2.0426702893139672e-05, 2.0097694793158127e-05, 
    1.9761217550014752e-05, 1.9418537476893117e-05, 1.9070869338660529e-05, 1.8719370246692075e-05, 1.8365134897196218e-05, 
    1.8009192076508293e-05, 1.7652502338533797e-05, 1.7295956746745589e-05, 1.6940376565373393e-05, 1.6586513781060617e-05, 
    1.6235052336629225e-05, 1.5886609961998042e-05, 1.5541740493072105e-05, 1.5200936576934429e-05, 1.4864632670365084e-05, 
    1.4533208248097372e-05, 1.4206991146887641e-05, 1.3886260981091449e-05, 1.3571252574745994e-05, 1.3262159363963743e-05, 
    1.2959136731612053e-05, 1.2662305243703281e-05, 1.237175376360697e-05, 1.208754242610922e-05, 1.1809705458499893e-05, 
    1.1538253840299238e-05, 1.1273177796988285e-05, 1.1014449126237733e-05, 1.0762023357696553e-05, 2.2558001438808169e-05, 
    2.2546355593976571e-05, 2.2511490095635372e-05, 2.2453619941147876e-05, 2.2373099787670393e-05, 2.2270418554580625e-05, 
    2.214619208191309e-05, 2.2001154048127263e-05, 2.18361453927663e-05, 2.1652102521334957e-05, 2.1450044590143721e-05, 
    2.1231060177701801e-05, 2.0996293646890681e-05, 2.0746931489564187e-05, 2.0484188923795782e-05, 2.0209296985438506e-05, 
    1.9923490321867068e-05, 1.962799585866634e-05, 1.9324022471474605e-05, 1.9012751756870967e-05, 1.869532995956019e-05, 
    1.837286107931536e-05, 1.8046401151046534e-05, 1.7716953665530471e-05, 1.7385466077043239e-05, 1.7052827327423904e-05, 
    1.6719866303802681e-05, 1.6387351139036627e-05, 1.6055989259386988e-05, 1.5726428082654883e-05, 1.5399256271344362e-05, 
    1.5075005448923926e-05, 1.4754152292409665e-05, 1.4437120920837226e-05, 1.4124285506317143e-05, 1.3815973041927432e-05, 
    1.3512466208394991e-05, 1.3214006289117282e-05, 1.2920796090394263e-05, 1.2633002830642929e-05, 1.2350760968759221e-05, 
    1.2074174947617531e-05, 1.1803321833927733e-05, 1.1538253840299238e-05, 1.1279000719404084e-05, 1.1025572023613385e-05, 
    1.0777959226439441e-05, 2.1566563350880302e-05, 2.1555918452878882e-05, 2.1524046715393102e-05, 2.1471136079416234e-05, 
    2.1397496710207062e-05, 2.1303556480747613e-05, 2.1189854820733919e-05, 2.1057035094278017e-05, 2.0905835703940928e-05, 
    2.0737080145077215e-05, 2.0551666252015343e-05, 2.0350554886082297e-05, 2.0134758315137213e-05, 1.9905328525766914e-05, 
    1.9663345693621164e-05, 1.9409907015792277e-05, 1.9146116083091899e-05, 1.8873072941025328e-05, 1.8591864957655429e-05, 
    1.8303558585721966e-05, 1.8009192076508293e-05, 1.7709769174992951e-05, 1.7406253800531131e-05, 1.7099565695196791e-05, 
    1.6790577003280513e-05, 1.6480109730389663e-05, 1.6168934019079085e-05, 1.5857767169764927e-05, 1.554727333055629e-05, 
    1.5238063777227012e-05, 1.4930697704455594e-05, 1.462568345128353e-05, 1.4323480087088427e-05, 1.4024499288867249e-05, 
    1.3729107445939204e-05, 1.3437627934010035e-05, 1.315034350663433e-05, 1.2867498758259462e-05, 1.2589302619063883e-05, 
    1.2315930847583931e-05, 1.204752849256072e-05, 1.1784212300467516e-05, 1.1526073049757881e-05, 1.1273177796988285e-05, 
    1.1025572023613385e-05, 1.0783281675440151e-05, 2.0639078821140177e-05, 2.0629329610549177e-05, 2.0600137162987083e-05, 
    2.0551666252015343e-05, 2.0484188923795782e-05, 2.0398080703191319e-05, 2.029381542078653e-05, 2.0171958792393818e-05, 
    2.0033160910794776e-05, 1.9878147831375788e-05, 1.9707712448345043e-05, 1.9522704866129908e-05, 1.9324022471474605e-05, 
    1.9112599906141191e-05, 1.8889399128690579e-05, 1.8655399737519273e-05, 1.8411589707222853e-05, 1.8158956667581126e-05, 
    1.7898479830136321e-05, 1.7631122642524264e-05, 1.7357826226364313e-05, 1.7079503631416746e-05, 1.6797034917512545e-05, 
    1.6511263056912141e-05, 1.6222990633555952e-05, 1.5932977302273846e-05, 1.5641937960434319e-05, 1.5350541576637124e-05, 
    1.5059410615700901e-05, 1.4769120996121951e-05, 1.4480202515097018e-05, 1.4193139676806177e-05, 1.3908372861634588e-05, 
    1.3626299777077261e-05, 1.3347277134942729e-05, 1.3071622503900814e-05, 1.2799616291187753e-05, 1.2531503812204481e-05, 
    1.2267497411666463e-05, 1.2007778604765636e-05, 1.1752500211391983e-05, 1.1501788460767085e-05, 1.1255745047817684e-05, 
    1.1014449126237733e-05, 1.0777959226439441e-05, 1.9770163325275059e-05, 1.9761217550014752e-05, 1.9734428732730505e-05, 
    1.9689941742800998e-05, 1.962799585866634e-05, 1.9548921569173058e-05, 1.9453136207364189e-05, 1.9341138523265107e-05, 
    1.9213502325355306e-05, 1.9070869338660529e-05, 1.8913941440247991e-05, 1.8743472440142749e-05, 1.8560259577361926e-05, 
    1.8365134897196218e-05, 1.8158956667581126e-05, 1.7942600980087353e-05, 1.7716953665530471e-05, 1.7482902636322104e-05, 
    1.7241330748335021e-05, 1.6993109255071053e-05, 1.6739091907065007e-05, 1.6480109730389663e-05, 1.6216966500383444e-05, 
    1.595043491071315e-05, 1.5681253423889629e-05, 1.5410123777534442e-05, 1.5137709111103159e-05, 1.4864632670365084e-05, 
    1.4591477041607952e-05, 1.4318783864109055e-05, 1.4047053967684984e-05, 1.3776747881871487e-05, 1.3508286664253823e-05, 
    1.3242052997432568e-05, 1.2978392506844757e-05, 1.271761525495999e-05, 1.2459997371050673e-05, 1.2205782779634605e-05, 
    1.1955184994670849e-05, 1.1708388950544937e-05, 1.146555284471846e-05, 1.1226809970573934e-05, 1.0992270522409414e-05, 
    1.0762023357696553e-05, 1.8954987320109639e-05, 1.8946763898278791e-05, 1.892213638887895e-05, 1.8881232505634504e-05, 
    1.8824263269626126e-05, 1.8751520302933707e-05, 1.8663372130569494e-05, 1.8560259577361926e-05, 1.8442690365512086e-05, 
    1.8311233033766296e-05, 1.816651031012172e-05, 1.8009192076508293e-05, 1.7839988065985546e-05, 1.7659640430866404e-05, 
    1.7468916314213047e-05, 1.7268600547857891e-05, 1.7059488588098674e-05, 1.6842379786167605e-05, 1.661807107516462e-05, 
    1.6387351139036627e-05, 1.6150995112992828e-05, 1.5909759849010057e-05, 1.56643797652556e-05, 1.541556328469206e-05, 
    1.516398985608771e-05, 1.4910307540297918e-05, 1.4655131136084768e-05, 1.439904081290228e-05, 1.4142581212931546e-05, 
    1.3886260981091449e-05, 1.3630552679629404e-05, 1.3375893043042148e-05, 1.3122683529306675e-05, 1.2871291124530686e-05, 
    1.2622049359980525e-05, 1.2375259502842903e-05, 1.2131191884870169e-05, 1.1890087336110159e-05, 1.1652158694112225e-05, 
    1.141759236223075e-05, 1.1186549893835198e-05, 1.095916958231684e-05, 1.0735568039708114e-05, 1.8189208990225995e-05, 
    1.8181636463585602e-05, 1.8158956667581126e-05, 1.812128248362349e-05, 1.8068800490497561e-05, 1.8001768666748808e-05, 
    1.7920513248064266e-05, 1.7825424810421472e-05, 1.7716953665530471e-05, 1.7595604667821361e-05, 1.7461931541596408e-05, 
    1.7316530842796435e-05, 1.7160035672114973e-05, 1.6993109255071053e-05, 1.6816438500397617e-05, 1.6630727641101528e-05, 
    1.6436692053267825e-05, 1.6235052336629225e-05, 1.6026528728635823e-05, 1.5811835910764884e-05, 1.5591678252600002e-05, 
    1.5366745526225408e-05, 1.5137709111103159e-05, 1.4905218698133996e-05, 1.4669899491277326e-05, 1.4432349896078191e-05, 
    1.4193139676806177e-05, 1.3952808557678151e-05, 1.3711865238785746e-05, 1.3470786793810186e-05, 1.3230018414278281e-05, 
    1.2989973463870498e-05, 1.2751033805994922e-05, 1.2513550368347452e-05, 1.2277843909342877e-05, 1.2044205952987481e-05, 
    1.1812899860841926e-05, 1.1584162012077613e-05, 1.1358203065158025e-05, 1.1135209277290312e-05, 1.0915343860418048e-05, 
    1.0698748355103531e-05, 1.7468916314213047e-05, 1.7461931541596408e-05, 1.7441010697097685e-05, 1.7406253800531131e-05, 
    1.7357826226364313e-05, 1.7295956746745589e-05, 1.7220934852339357e-05, 1.7133107408996709e-05, 1.703287472134657e-05, 
    1.6920686085057192e-05, 1.6797034917512545e-05, 1.6662453561820911e-05, 1.6517507861396601e-05, 1.636279160192305e-05, 
    1.6198920914515062e-05, 1.6026528728635823e-05, 1.5846259356143906e-05, 1.5658763279143999e-05, 1.5464692204508715e-05, 
    1.5264694437445863e-05, 1.5059410615700901e-05, 1.4849469835271206e-05, 1.463548618818117e-05, 1.4418055723186734e-05, 
    1.4197753831447534e-05, 1.3975133051370435e-05, 1.3750721280079538e-05, 1.3525020373345501e-05, 1.3298505111307131e-05, 
    1.3071622503900814e-05, 1.2844791407509592e-05, 1.2618402422864087e-05, 1.2392818043567709e-05, 1.2168373024667766e-05, 
    1.1945374941338239e-05, 1.1724104908867815e-05, 1.1504818436652424e-05, 1.1287746390677853e-05, 1.1073096040956545e-05, 
    1.0861052172477647e-05, 1.6790577003280513e-05, 1.6784124052856499e-05, 1.6764794927267797e-05, 1.6732678461889888e-05, 
    1.6687921584078187e-05, 1.6630727641101528e-05, 1.6561354109037773e-05, 1.6480109730389663e-05, 1.6387351139036627e-05, 
    1.6283479040094186e-05, 1.6168934019079085e-05, 1.6044192059343355e-05, 1.5909759849010057e-05, 1.5766169958676034e-05, 
    1.5613975969085667e-05, 1.5453747624038438e-05, 1.5286066078240328e-05, 1.511151930295246e-05, 1.4930697704455594e-05, 
    1.4744190001867865e-05, 1.4552579402043524e-05, 1.4356440100438072e-05, 1.4156334128211544e-05, 1.3952808557678151e-05, 
    1.3746393070674412e-05, 1.353759788764185e-05, 1.3326912049292831e-05, 1.311480203769748e-05, 1.290171071950742e-05, 
    1.2688056590799713e-05, 1.2474233300637707e-05, 1.2260609428841272e-05, 1.204752849256072e-05, 1.1835309155970897e-05, 
    1.1624245617655738e-05, 1.1414608150949455e-05, 1.1206643773552122e-05, 1.1000577024063628e-05, 1.0796610824606334e-05, 
    1.6150995112992828e-05, 1.6145024319975089e-05, 1.6127138399384272e-05, 1.6097416434033397e-05, 1.6055989259386988e-05, 
    1.6003038030609236e-05, 1.5938792257493654e-05, 1.5863527346724523e-05, 1.5777561699975654e-05, 1.5681253423889629e-05, 
    1.5574996713813342e-05, 1.5459217977179684e-05, 1.5334371764583079e-05, 1.5200936576934429e-05, 1.5059410615700901e-05, 
    1.4910307540297918e-05, 1.4754152292409665e-05, 1.4591477041607952e-05, 1.4422817300374048e-05, 1.4248708249765942e-05, 
    1.4069681309772104e-05, 1.3886260981091449e-05, 1.3698961977896053e-05, 1.3508286664253823e-05, 1.3314722800467258e-05, 
    1.3118741599739446e-05, 1.2920796090394263e-05, 1.2721319774405071e-05, 1.2520725569246736e-05, 1.2319405017075489e-05, 
    1.2117727742933577e-05, 1.1916041142027999e-05, 1.1714670275089221e-05, 1.1513917950311788e-05, 1.1314064970345236e-05, 
    1.111537052316941e-05, 1.0918072696383152e-05, 1.0722389095392245e-05, 1.554727333055629e-05, 1.5541740493072105e-05, 
    1.5525165583196803e-05, 1.5497619157392695e-05, 1.5459217977179684e-05, 1.5410123777534442e-05, 1.5350541576637124e-05, 
    1.5280717559668513e-05, 1.5200936576934429e-05, 1.511151930295246e-05, 1.5012819108123966e-05, 1.4905218698133996e-05, 
    1.4789126578236576e-05, 1.4664973400111688e-05, 1.4533208248097372e-05, 1.4394294919424065e-05, 1.4248708249765942e-05, 
    1.4096930531159654e-05, 1.3939448064325759e-05, 1.3776747881871487e-05, 1.3609314672961238e-05, 1.3437627934010035e-05, 
    1.3262159363963743e-05, 1.308337051693607e-05, 1.290171071950742e-05, 1.271761525495999e-05, 1.2531503812204481e-05, 
    1.2343779193197459e-05, 1.2154826269282664e-05, 1.1965011174118524e-05, 1.1774680718666113e-05, 1.1584162012077613e-05, 
    1.1393762271206004e-05, 1.1203768800803646e-05, 1.1014449126237733e-05, 1.082605126066748e-05, 1.497678010477799e-05, 
    1.4971645795520265e-05, 1.4956263967648155e-05, 1.4930697704455594e-05, 1.4895051427534998e-05, 1.4849469835271206e-05, 
    1.4794136444963624e-05, 1.4729271765778285e-05, 1.4655131136084768e-05, 1.4572002264108309e-05, 1.4480202515097018e-05, 
    1.4380075991285024e-05, 1.4271990452788436e-05, 1.4156334128211544e-05, 1.4033512463217423e-05, 1.3903944853719389e-05, 
    1.3768061407797166e-05, 1.3626299777077261e-05, 1.3479102094300188e-05, 1.3326912049292831e-05, 1.317017213073963e-05, 
    1.3009321056161041e-05, 1.2844791407509592e-05, 1.2677007484915127e-05, 1.2506383386464091e-05, 1.2333321317574871e-05, 
    1.215821012960262e-05, 1.198142408382239e-05, 1.1803321833927733e-05, 1.1624245617655738e-05, 1.1444520646103932e-05, 
    1.1264454677723139e-05, 1.1084337762825537e-05, 1.0904442143703522e-05, 1.0725022295071862e-05, 1.4437120920837226e-05, 
    1.4432349896078191e-05, 1.4418055723186734e-05, 1.4394294919424065e-05, 1.4361161060681554e-05, 1.4318783864109055e-05, 
    1.4267327927321989e-05, 1.4206991146887641e-05, 1.4138002844134871e-05, 1.4060621630886222e-05, 1.3975133051370435e-05, 
    1.3881847039266564e-05, 1.378109523052465e-05, 1.3673228173303888e-05, 1.3558612476104504e-05, 1.3437627934010035e-05, 
    1.3310664670994396e-05, 1.3178120333594629e-05, 1.3040397368030042e-05, 1.2897900409194513e-05, 1.2751033805994922e-05, 
    1.2600199303385058e-05, 1.2445793897273469e-05, 1.2288207874376087e-05, 1.2127823045135412e-05, 1.1965011174118524e-05, 
    1.1800132608898301e-05, 1.1633535105362978e-05, 1.146555284471846e-05, 1.1296505635161048e-05, 1.1126698289307673e-05, 
    1.0956420166967538e-05, 1.0785944871704769e-05, 1.3926113133141779e-05, 1.3921673823886714e-05, 1.3908372861634588e-05, 
    1.3886260981091449e-05, 1.3855422203531922e-05, 1.3815973041927432e-05, 1.3768061407797166e-05, 1.3711865238785746e-05, 
    1.364759087047894e-05, 1.3575471179836062e-05, 1.3495763530757913e-05, 1.3408747554661532e-05, 1.3314722800467258e-05, 
    1.3214006289117282e-05, 1.310693000766285e-05, 1.2993838377129607e-05, 1.2875085726866928e-05, 1.2751033805994922e-05, 
    1.2622049359980525e-05, 1.2488501797407096e-05, 1.2350760968759221e-05, 1.2209195075631148e-05, 1.2064168725285251e-05, 
    1.1916041142027999e-05, 1.176516454351633e-05, 1.1611882686927043e-05, 1.145652958697078e-05, 1.1299428405053716e-05, 
    1.114089050651342e-05, 1.0981214680797741e-05, 1.0820686517723644e-05, 1.3441763861351989e-05, 1.3437627934010035e-05, 
    1.3425235409017095e-05, 1.3404631916983615e-05, 1.3375893043042148e-05, 1.333912363638748e-05, 1.3294456860131694e-05, 
    1.3242052997432568e-05, 1.318209803366514e-05, 1.311480203769748e-05, 1.3040397368030042e-05, 1.2959136731612053e-05, 
    1.2871291124530686e-05, 1.277714768447414e-05, 1.2677007484915127e-05, 1.2571183300383597e-05, 1.2459997371050673e-05, 
    1.2343779193197459e-05, 1.2222863360070701e-05, 1.209758747521679e-05, 1.1968290157723379e-05, 1.1835309155970897e-05, 
    1.1698979583587625e-05, 1.1559632288388728e-05, 1.141759236223075e-05, 1.1273177796988285e-05, 1.1126698289307673e-05, 
    1.0978454194452643e-05, 1.0828735627456636e-05, 1.0677821707954185e-05, 1.2982250530776638e-05, 1.2978392506844757e-05, 
    1.2966832180977615e-05, 1.2947610668702227e-05, 1.2920796090394263e-05, 1.2886482970059787e-05, 1.2844791407509592e-05, 
    1.279586603736672e-05, 1.2739874791578941e-05, 1.2677007484915127e-05, 1.2607474245246135e-05, 1.2531503812204481e-05, 
    1.2449341729057189e-05, 1.236124845330671e-05, 1.2267497411666463e-05, 1.2168373024667766e-05, 1.2064168725285251e-05, 
    1.1955184994670849e-05, 1.1841727436424244e-05, 1.1724104908867815e-05, 1.1602627732606965e-05, 1.147760598831343e-05, 
    1.1349347917238203e-05, 1.1218158434506195e-05, 1.1084337762825537e-05, 1.0948180191910906e-05, 1.080997296671599e-05, 
    1.0669995305529591e-05, 1.2545903701675557e-05, 1.2542300627665885e-05, 1.2531503812204481e-05, 1.2513550368347452e-05, 
    1.2488501797407096e-05, 1.2456443464213523e-05, 1.2417483874191815e-05, 1.237175376360697e-05, 1.2319405017075489e-05, 
    1.2260609428841272e-05, 1.2195557326314608e-05, 1.2124456075939094e-05, 1.204752849256072e-05, 1.1965011174118524e-05, 
    1.1877152783664023e-05, 1.1784212300467516e-05, 1.1686457261314587e-05, 1.1584162012077613e-05, 1.147760598831343e-05, 
    1.1367072042043886e-05, 1.1252844830077971e-05, 1.1135209277290312e-05, 1.1014449126237733e-05, 1.0890845582427085e-05, 
    1.0764676062492632e-05, 1.2131191884870169e-05, 1.2127823045135412e-05, 1.2117727742933577e-05, 1.2100939536030092e-05, 
    1.2077514044671629e-05, 1.204752849256072e-05, 1.2011081074128884e-05, 1.1968290157723379e-05, 1.191929333666283e-05, 
    1.1864246342171315e-05, 1.1803321833927733e-05, 1.1736708085335289e-05, 1.1664607581605931e-05, 1.1587235549358612e-05, 
    1.1504818436652424e-05, 1.141759236223075e-05, 1.1325801552264682e-05, 1.1229696782086041e-05, 1.1129533839330432e-05, 
    1.1025572023613385e-05, 1.0918072696383152e-05, 1.0807297892980108e-05, 1.0693509007231295e-05, 1.1736708085335289e-05, 
    1.1733554751620799e-05, 1.1724104908867815e-05, 1.1708388950544937e-05, 1.1686457261314587e-05, 1.1658379814610948e-05, 
    1.1624245617655738e-05, 1.1584162012077613e-05, 1.1538253840299238e-05, 1.148666248961931e-05, 1.1429544827409738e-05, 
    1.1367072042043886e-05, 1.1299428405053716e-05, 1.1226809970573934e-05, 1.114942322837187e-05, 1.1067483726693515e-05, 
    1.0981214680797741e-05, 1.0890845582427085e-05, 1.0796610824606334e-05, 1.0698748355103531e-05, 1.1361157852635955e-05, 
    1.1358203065158025e-05, 1.1349347917238203e-05, 1.1334619980672883e-05, 1.1314064970345236e-05, 1.1287746390677853e-05, 
    1.1255745047817684e-05, 1.1218158434506195e-05, 1.1175099996298006e-05, 1.1126698289307673e-05, 1.1073096040956545e-05, 
    1.1014449126237733e-05, 1.0950925472801559e-05, 1.0882703908679942e-05, 1.080997296671599e-05, 1.0732929659752426e-05, 
    1.1003348648408319e-05, 1.1000577024063628e-05, 1.0992270522409414e-05, 1.0978454194452643e-05, 1.095916958231684e-05, 
    1.0934474407995146e-05, 1.0904442143703522e-05, 1.0869161469769193e-05, 1.0828735627456636e-05, 1.0783281675440151e-05, 
    1.0732929659752426e-05, 1.0677821707954185e-05, 1.0662180367561673e-05};

    std::vector<double> kx {0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 
        0, 1, 1, 1, 1, 
        1, 1, 1, 1, 1, 
        1, 1, 1, 1, 1, 
        1, 1, 1, 1, 1, 
        1, 1, 1, 1, 1, 
        1, 1, 1, 1, 1, 
        1, 1, 1, 1, 1, 
        1, 1, 1, 1, 1, 
        1, 1, 1, 1, 1, 
        1, 1, 1, 1, 1, 
        1, 1, 1, 1, 1, 
        1, 1, 1, 1, 1, 
        1, 2, 2, 2, 2, 
        2, 2, 2, 2, 2, 
        2, 2, 2, 2, 2, 
        2, 2, 2, 2, 2, 
        2, 2, 2, 2, 2, 
        2, 2, 2, 2, 2, 
        2, 2, 2, 2, 2, 
        2, 2, 2, 2, 2, 
        2, 2, 2, 2, 2, 
        2, 2, 2, 2, 2, 
        2, 2, 2, 2, 2, 
        2, 2, 2, 2, 2, 
        2, 3, 3, 3, 3, 
        3, 3, 3, 3, 3, 
        3, 3, 3, 3, 3, 
        3, 3, 3, 3, 3, 
        3, 3, 3, 3, 3, 
        3, 3, 3, 3, 3, 
        3, 3, 3, 3, 3, 
        3, 3, 3, 3, 3, 
        3, 3, 3, 3, 3, 
        3, 3, 3, 3, 3, 
        3, 3, 3, 3, 3, 
        3, 3, 3, 3, 3, 
        3, 3, 4, 4, 4, 
        4, 4, 4, 4, 4, 
        4, 4, 4, 4, 4, 
        4, 4, 4, 4, 4, 
        4, 4, 4, 4, 4, 
        4, 4, 4, 4, 4, 
        4, 4, 4, 4, 4, 
        4, 4, 4, 4, 4, 
        4, 4, 4, 4, 4, 
        4, 4, 4, 4, 4, 
        4, 4, 4, 4, 4, 
        4, 4, 4, 4, 4, 
        4, 4, 4, 4, 4, 
        4, 5, 5, 5, 5, 
        5, 5, 5, 5, 5, 
        5, 5, 5, 5, 5, 
        5, 5, 5, 5, 5, 
        5, 5, 5, 5, 5, 
        5, 5, 5, 5, 5, 
        5, 5, 5, 5, 5, 
        5, 5, 5, 5, 5, 
        5, 5, 5, 5, 5, 
        5, 5, 5, 5, 5, 
        5, 5, 5, 5, 5, 
        5, 5, 5, 5, 5, 
        5, 5, 5, 5, 5, 
        6, 6, 6, 6, 6, 
        6, 6, 6, 6, 6, 
        6, 6, 6, 6, 6, 
        6, 6, 6, 6, 6, 
        6, 6, 6, 6, 6, 
        6, 6, 6, 6, 6, 
        6, 6, 6, 6, 6, 
        6, 6, 6, 6, 6, 
        6, 6, 6, 6, 6, 
        6, 6, 6, 6, 6, 
        6, 6, 6, 6, 6, 
        6, 6, 6, 6, 6, 
        6, 6, 6, 6, 7, 
        7, 7, 7, 7, 7, 
        7, 7, 7, 7, 7, 
        7, 7, 7, 7, 7, 
        7, 7, 7, 7, 7, 
        7, 7, 7, 7, 7, 
        7, 7, 7, 7, 7, 
        7, 7, 7, 7, 7, 
        7, 7, 7, 7, 7, 
        7, 7, 7, 7, 7, 
        7, 7, 7, 7, 7, 
        7, 7, 7, 7, 7, 
        7, 7, 7, 7, 7, 
        7, 7, 7, 8, 8, 
        8, 8, 8, 8, 8, 
        8, 8, 8, 8, 8, 
        8, 8, 8, 8, 8, 
        8, 8, 8, 8, 8, 
        8, 8, 8, 8, 8, 
        8, 8, 8, 8, 8, 
        8, 8, 8, 8, 8, 
        8, 8, 8, 8, 8, 
        8, 8, 8, 8, 8, 
        8, 8, 8, 8, 8, 
        8, 8, 8, 8, 8, 
        8, 8, 8, 8, 8, 
        8, 8, 9, 9, 9, 
        9, 9, 9, 9, 9, 
        9, 9, 9, 9, 9, 
        9, 9, 9, 9, 9, 
        9, 9, 9, 9, 9, 
        9, 9, 9, 9, 9, 
        9, 9, 9, 9, 9, 
        9, 9, 9, 9, 9, 
        9, 9, 9, 9, 9, 
        9, 9, 9, 9, 9, 
        9, 9, 9, 9, 9, 
        9, 9, 9, 9, 9, 
        9, 9, 9, 9, 9, 
        9, 10, 10, 10, 10, 
        10, 10, 10, 10, 10, 
        10, 10, 10, 10, 10, 
        10, 10, 10, 10, 10, 
        10, 10, 10, 10, 10, 
        10, 10, 10, 10, 10, 
        10, 10, 10, 10, 10, 
        10, 10, 10, 10, 10, 
        10, 10, 10, 10, 10, 
        10, 10, 10, 10, 10, 
        10, 10, 10, 10, 10, 
        10, 10, 10, 10, 10, 
        10, 10, 10, 10, 10, 
        11, 11, 11, 11, 11, 
        11, 11, 11, 11, 11, 
        11, 11, 11, 11, 11, 
        11, 11, 11, 11, 11, 
        11, 11, 11, 11, 11, 
        11, 11, 11, 11, 11, 
        11, 11, 11, 11, 11, 
        11, 11, 11, 11, 11, 
        11, 11, 11, 11, 11, 
        11, 11, 11, 11, 11, 
        11, 11, 11, 11, 11, 
        11, 11, 11, 11, 11, 
        11, 11, 11, 11, 12, 
        12, 12, 12, 12, 12, 
        12, 12, 12, 12, 12, 
        12, 12, 12, 12, 12, 
        12, 12, 12, 12, 12, 
        12, 12, 12, 12, 12, 
        12, 12, 12, 12, 12, 
        12, 12, 12, 12, 12, 
        12, 12, 12, 12, 12, 
        12, 12, 12, 12, 12, 
        12, 12, 12, 12, 12, 
        12, 12, 12, 12, 12, 
        12, 12, 12, 12, 12, 
        12, 12, 13, 13, 13, 
        13, 13, 13, 13, 13, 
        13, 13, 13, 13, 13, 
        13, 13, 13, 13, 13, 
        13, 13, 13, 13, 13, 
        13, 13, 13, 13, 13, 
        13, 13, 13, 13, 13, 
        13, 13, 13, 13, 13, 
        13, 13, 13, 13, 13, 
        13, 13, 13, 13, 13, 
        13, 13, 13, 13, 13, 
        13, 13, 13, 13, 13, 
        13, 13, 13, 13, 13, 
        14, 14, 14, 14, 14, 
        14, 14, 14, 14, 14, 
        14, 14, 14, 14, 14, 
        14, 14, 14, 14, 14, 
        14, 14, 14, 14, 14, 
        14, 14, 14, 14, 14, 
        14, 14, 14, 14, 14, 
        14, 14, 14, 14, 14, 
        14, 14, 14, 14, 14, 
        14, 14, 14, 14, 14, 
        14, 14, 14, 14, 14, 
        14, 14, 14, 14, 14, 
        14, 14, 14, 15, 15, 
        15, 15, 15, 15, 15, 
        15, 15, 15, 15, 15, 
        15, 15, 15, 15, 15, 
        15, 15, 15, 15, 15, 
        15, 15, 15, 15, 15, 
        15, 15, 15, 15, 15, 
        15, 15, 15, 15, 15, 
        15, 15, 15, 15, 15, 
        15, 15, 15, 15, 15, 
        15, 15, 15, 15, 15, 
        15, 15, 15, 15, 15, 
        15, 15, 15, 15, 15, 
        15, 16, 16, 16, 16, 
        16, 16, 16, 16, 16, 
        16, 16, 16, 16, 16, 
        16, 16, 16, 16, 16, 
        16, 16, 16, 16, 16, 
        16, 16, 16, 16, 16, 
        16, 16, 16, 16, 16, 
        16, 16, 16, 16, 16, 
        16, 16, 16, 16, 16, 
        16, 16, 16, 16, 16, 
        16, 16, 16, 16, 16, 
        16, 16, 16, 16, 16, 
        16, 16, 16, 17, 17, 
        17, 17, 17, 17, 17, 
        17, 17, 17, 17, 17, 
        17, 17, 17, 17, 17, 
        17, 17, 17, 17, 17, 
        17, 17, 17, 17, 17, 
        17, 17, 17, 17, 17, 
        17, 17, 17, 17, 17, 
        17, 17, 17, 17, 17, 
        17, 17, 17, 17, 17, 
        17, 17, 17, 17, 17, 
        17, 17, 17, 17, 17, 
        17, 17, 17, 17, 17, 
        18, 18, 18, 18, 18, 
        18, 18, 18, 18, 18, 
        18, 18, 18, 18, 18, 
        18, 18, 18, 18, 18, 
        18, 18, 18, 18, 18, 
        18, 18, 18, 18, 18, 
        18, 18, 18, 18, 18, 
        18, 18, 18, 18, 18, 
        18, 18, 18, 18, 18, 
        18, 18, 18, 18, 18, 
        18, 18, 18, 18, 18, 
        18, 18, 18, 18, 18, 
        18, 18, 19, 19, 19, 
        19, 19, 19, 19, 19, 
        19, 19, 19, 19, 19, 
        19, 19, 19, 19, 19, 
        19, 19, 19, 19, 19, 
        19, 19, 19, 19, 19, 
        19, 19, 19, 19, 19, 
        19, 19, 19, 19, 19, 
        19, 19, 19, 19, 19, 
        19, 19, 19, 19, 19, 
        19, 19, 19, 19, 19, 
        19, 19, 19, 19, 19, 
        19, 19, 19, 19, 20, 
        20, 20, 20, 20, 20, 
        20, 20, 20, 20, 20, 
        20, 20, 20, 20, 20, 
        20, 20, 20, 20, 20, 
        20, 20, 20, 20, 20, 
        20, 20, 20, 20, 20, 
        20, 20, 20, 20, 20, 
        20, 20, 20, 20, 20, 
        20, 20, 20, 20, 20, 
        20, 20, 20, 20, 20, 
        20, 20, 20, 20, 20, 
        20, 20, 20, 20, 20, 
        21, 21, 21, 21, 21, 
        21, 21, 21, 21, 21, 
        21, 21, 21, 21, 21, 
        21, 21, 21, 21, 21, 
        21, 21, 21, 21, 21, 
        21, 21, 21, 21, 21, 
        21, 21, 21, 21, 21, 
        21, 21, 21, 21, 21, 
        21, 21, 21, 21, 21, 
        21, 21, 21, 21, 21, 
        21, 21, 21, 21, 21, 
        21, 21, 21, 21, 21, 
        21, 22, 22, 22, 22, 
        22, 22, 22, 22, 22, 
        22, 22, 22, 22, 22, 
        22, 22, 22, 22, 22, 
        22, 22, 22, 22, 22, 
        22, 22, 22, 22, 22, 
        22, 22, 22, 22, 22, 
        22, 22, 22, 22, 22, 
        22, 22, 22, 22, 22, 
        22, 22, 22, 22, 22, 
        22, 22, 22, 22, 22, 
        22, 22, 22, 22, 22, 
        22, 22, 23, 23, 23, 
        23, 23, 23, 23, 23, 
        23, 23, 23, 23, 23, 
        23, 23, 23, 23, 23, 
        23, 23, 23, 23, 23, 
        23, 23, 23, 23, 23, 
        23, 23, 23, 23, 23, 
        23, 23, 23, 23, 23, 
        23, 23, 23, 23, 23, 
        23, 23, 23, 23, 23, 
        23, 23, 23, 23, 23, 
        23, 23, 23, 23, 23, 
        23, 23, 24, 24, 24, 
        24, 24, 24, 24, 24, 
        24, 24, 24, 24, 24, 
        24, 24, 24, 24, 24, 
        24, 24, 24, 24, 24, 
        24, 24, 24, 24, 24, 
        24, 24, 24, 24, 24, 
        24, 24, 24, 24, 24, 
        24, 24, 24, 24, 24, 
        24, 24, 24, 24, 24, 
        24, 24, 24, 24, 24, 
        24, 24, 24, 24, 24, 
        24, 24, 25, 25, 25, 
        25, 25, 25, 25, 25, 
        25, 25, 25, 25, 25, 
        25, 25, 25, 25, 25, 
        25, 25, 25, 25, 25, 
        25, 25, 25, 25, 25, 
        25, 25, 25, 25, 25, 
        25, 25, 25, 25, 25, 
        25, 25, 25, 25, 25, 
        25, 25, 25, 25, 25, 
        25, 25, 25, 25, 25, 
        25, 25, 25, 25, 25, 
        25, 26, 26, 26, 26, 
        26, 26, 26, 26, 26, 
        26, 26, 26, 26, 26, 
        26, 26, 26, 26, 26, 
        26, 26, 26, 26, 26, 
        26, 26, 26, 26, 26, 
        26, 26, 26, 26, 26, 
        26, 26, 26, 26, 26, 
        26, 26, 26, 26, 26, 
        26, 26, 26, 26, 26, 
        26, 26, 26, 26, 26, 
        26, 26, 26, 26, 26, 
        27, 27, 27, 27, 27, 
        27, 27, 27, 27, 27, 
        27, 27, 27, 27, 27, 
        27, 27, 27, 27, 27, 
        27, 27, 27, 27, 27, 
        27, 27, 27, 27, 27, 
        27, 27, 27, 27, 27, 
        27, 27, 27, 27, 27, 
        27, 27, 27, 27, 27, 
        27, 27, 27, 27, 27, 
        27, 27, 27, 27, 27, 
        27, 27, 27, 27, 28, 
        28, 28, 28, 28, 28, 
        28, 28, 28, 28, 28, 
        28, 28, 28, 28, 28, 
        28, 28, 28, 28, 28, 
        28, 28, 28, 28, 28, 
        28, 28, 28, 28, 28, 
        28, 28, 28, 28, 28, 
        28, 28, 28, 28, 28, 
        28, 28, 28, 28, 28, 
        28, 28, 28, 28, 28, 
        28, 28, 28, 28, 28, 
        28, 28, 29, 29, 29, 
        29, 29, 29, 29, 29, 
        29, 29, 29, 29, 29, 
        29, 29, 29, 29, 29, 
        29, 29, 29, 29, 29, 
        29, 29, 29, 29, 29, 
        29, 29, 29, 29, 29, 
        29, 29, 29, 29, 29, 
        29, 29, 29, 29, 29, 
        29, 29, 29, 29, 29, 
        29, 29, 29, 29, 29, 
        29, 29, 29, 29, 29, 
        30, 30, 30, 30, 30, 
        30, 30, 30, 30, 30, 
        30, 30, 30, 30, 30, 
        30, 30, 30, 30, 30, 
        30, 30, 30, 30, 30, 
        30, 30, 30, 30, 30, 
        30, 30, 30, 30, 30, 
        30, 30, 30, 30, 30, 
        30, 30, 30, 30, 30, 
        30, 30, 30, 30, 30, 
        30, 30, 30, 30, 30, 
        30, 30, 31, 31, 31, 
        31, 31, 31, 31, 31, 
        31, 31, 31, 31, 31, 
        31, 31, 31, 31, 31, 
        31, 31, 31, 31, 31, 
        31, 31, 31, 31, 31, 
        31, 31, 31, 31, 31, 
        31, 31, 31, 31, 31, 
        31, 31, 31, 31, 31, 
        31, 31, 31, 31, 31, 
        31, 31, 31, 31, 31, 
        31, 31, 31, 32, 32, 
        32, 32, 32, 32, 32, 
        32, 32, 32, 32, 32, 
        32, 32, 32, 32, 32, 
        32, 32, 32, 32, 32, 
        32, 32, 32, 32, 32, 
        32, 32, 32, 32, 32, 
        32, 32, 32, 32, 32, 
        32, 32, 32, 32, 32, 
        32, 32, 32, 32, 32, 
        32, 32, 32, 32, 32, 
        32, 32, 32, 32, 33, 
        33, 33, 33, 33, 33, 
        33, 33, 33, 33, 33, 
        33, 33, 33, 33, 33, 
        33, 33, 33, 33, 33, 
        33, 33, 33, 33, 33, 
        33, 33, 33, 33, 33, 
        33, 33, 33, 33, 33, 
        33, 33, 33, 33, 33, 
        33, 33, 33, 33, 33, 
        33, 33, 33, 33, 33, 
        33, 33, 33, 33, 34, 
        34, 34, 34, 34, 34, 
        34, 34, 34, 34, 34, 
        34, 34, 34, 34, 34, 
        34, 34, 34, 34, 34, 
        34, 34, 34, 34, 34, 
        34, 34, 34, 34, 34, 
        34, 34, 34, 34, 34, 
        34, 34, 34, 34, 34, 
        34, 34, 34, 34, 34, 
        34, 34, 34, 34, 34, 
        34, 34, 34, 34, 35, 
        35, 35, 35, 35, 35, 
        35, 35, 35, 35, 35, 
        35, 35, 35, 35, 35, 
        35, 35, 35, 35, 35, 
        35, 35, 35, 35, 35, 
        35, 35, 35, 35, 35, 
        35, 35, 35, 35, 35, 
        35, 35, 35, 35, 35, 
        35, 35, 35, 35, 35, 
        35, 35, 35, 35, 35, 
        35, 35, 35, 36, 36, 
        36, 36, 36, 36, 36, 
        36, 36, 36, 36, 36, 
        36, 36, 36, 36, 36, 
        36, 36, 36, 36, 36, 
        36, 36, 36, 36, 36, 
        36, 36, 36, 36, 36, 
        36, 36, 36, 36, 36, 
        36, 36, 36, 36, 36, 
        36, 36, 36, 36, 36, 
        36, 36, 36, 36, 36, 
        36, 37, 37, 37, 37, 
        37, 37, 37, 37, 37, 
        37, 37, 37, 37, 37, 
        37, 37, 37, 37, 37, 
        37, 37, 37, 37, 37, 
        37, 37, 37, 37, 37, 
        37, 37, 37, 37, 37, 
        37, 37, 37, 37, 37, 
        37, 37, 37, 37, 37, 
        37, 37, 37, 37, 37, 
        37, 37, 37, 37, 38, 
        38, 38, 38, 38, 38, 
        38, 38, 38, 38, 38, 
        38, 38, 38, 38, 38, 
        38, 38, 38, 38, 38, 
        38, 38, 38, 38, 38, 
        38, 38, 38, 38, 38, 
        38, 38, 38, 38, 38, 
        38, 38, 38, 38, 38, 
        38, 38, 38, 38, 38, 
        38, 38, 38, 38, 38, 
        38, 39, 39, 39, 39, 
        39, 39, 39, 39, 39, 
        39, 39, 39, 39, 39, 
        39, 39, 39, 39, 39, 
        39, 39, 39, 39, 39, 
        39, 39, 39, 39, 39, 
        39, 39, 39, 39, 39, 
        39, 39, 39, 39, 39, 
        39, 39, 39, 39, 39, 
        39, 39, 39, 39, 39, 
        39, 39, 40, 40, 40, 
        40, 40, 40, 40, 40, 
        40, 40, 40, 40, 40, 
        40, 40, 40, 40, 40, 
        40, 40, 40, 40, 40, 
        40, 40, 40, 40, 40, 
        40, 40, 40, 40, 40, 
        40, 40, 40, 40, 40, 
        40, 40, 40, 40, 40, 
        40, 40, 40, 40, 40, 
        40, 40, 41, 41, 41, 
        41, 41, 41, 41, 41, 
        41, 41, 41, 41, 41, 
        41, 41, 41, 41, 41, 
        41, 41, 41, 41, 41, 
        41, 41, 41, 41, 41, 
        41, 41, 41, 41, 41, 
        41, 41, 41, 41, 41, 
        41, 41, 41, 41, 41, 
        41, 41, 41, 41, 41, 
        41, 41, 42, 42, 42, 
        42, 42, 42, 42, 42, 
        42, 42, 42, 42, 42, 
        42, 42, 42, 42, 42, 
        42, 42, 42, 42, 42, 
        42, 42, 42, 42, 42, 
        42, 42, 42, 42, 42, 
        42, 42, 42, 42, 42, 
        42, 42, 42, 42, 42, 
        42, 42, 42, 42, 42, 
        42, 43, 43, 43, 43, 
        43, 43, 43, 43, 43, 
        43, 43, 43, 43, 43, 
        43, 43, 43, 43, 43, 
        43, 43, 43, 43, 43, 
        43, 43, 43, 43, 43, 
        43, 43, 43, 43, 43, 
        43, 43, 43, 43, 43, 
        43, 43, 43, 43, 43, 
        43, 43, 43, 43, 44, 
        44, 44, 44, 44, 44, 
        44, 44, 44, 44, 44, 
        44, 44, 44, 44, 44, 
        44, 44, 44, 44, 44, 
        44, 44, 44, 44, 44, 
        44, 44, 44, 44, 44, 
        44, 44, 44, 44, 44, 
        44, 44, 44, 44, 44, 
        44, 44, 44, 44, 44, 
        44, 45, 45, 45, 45, 
        45, 45, 45, 45, 45, 
        45, 45, 45, 45, 45, 
        45, 45, 45, 45, 45, 
        45, 45, 45, 45, 45, 
        45, 45, 45, 45, 45, 
        45, 45, 45, 45, 45, 
        45, 45, 45, 45, 45, 
        45, 45, 45, 45, 45, 
        45, 45, 46, 46, 46, 
        46, 46, 46, 46, 46, 
        46, 46, 46, 46, 46, 
        46, 46, 46, 46, 46, 
        46, 46, 46, 46, 46, 
        46, 46, 46, 46, 46, 
        46, 46, 46, 46, 46, 
        46, 46, 46, 46, 46, 
        46, 46, 46, 46, 46, 
        46, 46, 47, 47, 47, 
        47, 47, 47, 47, 47, 
        47, 47, 47, 47, 47, 
        47, 47, 47, 47, 47, 
        47, 47, 47, 47, 47, 
        47, 47, 47, 47, 47, 
        47, 47, 47, 47, 47, 
        47, 47, 47, 47, 47, 
        47, 47, 47, 47, 47, 
        47, 48, 48, 48, 48, 
        48, 48, 48, 48, 48, 
        48, 48, 48, 48, 48, 
        48, 48, 48, 48, 48, 
        48, 48, 48, 48, 48, 
        48, 48, 48, 48, 48, 
        48, 48, 48, 48, 48, 
        48, 48, 48, 48, 48, 
        48, 48, 48, 48, 49, 
        49, 49, 49, 49, 49, 
        49, 49, 49, 49, 49, 
        49, 49, 49, 49, 49, 
        49, 49, 49, 49, 49, 
        49, 49, 49, 49, 49, 
        49, 49, 49, 49, 49, 
        49, 49, 49, 49, 49, 
        49, 49, 49, 49, 49, 
        49, 50, 50, 50, 50, 
        50, 50, 50, 50, 50, 
        50, 50, 50, 50, 50, 
        50, 50, 50, 50, 50, 
        50, 50, 50, 50, 50, 
        50, 50, 50, 50, 50, 
        50, 50, 50, 50, 50, 
        50, 50, 50, 50, 50, 
        50, 51, 51, 51, 51, 
        51, 51, 51, 51, 51, 
        51, 51, 51, 51, 51, 
        51, 51, 51, 51, 51, 
        51, 51, 51, 51, 51, 
        51, 51, 51, 51, 51, 
        51, 51, 51, 51, 51, 
        51, 51, 51, 51, 51, 
        52, 52, 52, 52, 52, 
        52, 52, 52, 52, 52, 
        52, 52, 52, 52, 52, 
        52, 52, 52, 52, 52, 
        52, 52, 52, 52, 52, 
        52, 52, 52, 52, 52, 
        52, 52, 52, 52, 52, 
        52, 52, 52, 53, 53, 
        53, 53, 53, 53, 53, 
        53, 53, 53, 53, 53, 
        53, 53, 53, 53, 53, 
        53, 53, 53, 53, 53, 
        53, 53, 53, 53, 53, 
        53, 53, 53, 53, 53, 
        53, 53, 53, 53, 54, 
        54, 54, 54, 54, 54, 
        54, 54, 54, 54, 54, 
        54, 54, 54, 54, 54, 
        54, 54, 54, 54, 54, 
        54, 54, 54, 54, 54, 
        54, 54, 54, 54, 54, 
        54, 54, 54, 54, 55, 
        55, 55, 55, 55, 55, 
        55, 55, 55, 55, 55, 
        55, 55, 55, 55, 55, 
        55, 55, 55, 55, 55, 
        55, 55, 55, 55, 55, 
        55, 55, 55, 55, 55, 
        55, 55, 56, 56, 56, 
        56, 56, 56, 56, 56, 
        56, 56, 56, 56, 56, 
        56, 56, 56, 56, 56, 
        56, 56, 56, 56, 56, 
        56, 56, 56, 56, 56, 
        56, 56, 56, 57, 57, 
        57, 57, 57, 57, 57, 
        57, 57, 57, 57, 57, 
        57, 57, 57, 57, 57, 
        57, 57, 57, 57, 57, 
        57, 57, 57, 57, 57, 
        57, 57, 57, 58, 58, 
        58, 58, 58, 58, 58, 
        58, 58, 58, 58, 58, 
        58, 58, 58, 58, 58, 
        58, 58, 58, 58, 58, 
        58, 58, 58, 58, 58, 
        58, 59, 59, 59, 59, 
        59, 59, 59, 59, 59, 
        59, 59, 59, 59, 59, 
        59, 59, 59, 59, 59, 
        59, 59, 59, 59, 59, 
        59, 60, 60, 60, 60, 
        60, 60, 60, 60, 60, 
        60, 60, 60, 60, 60, 
        60, 60, 60, 60, 60, 
        60, 60, 60, 60, 61, 
        61, 61, 61, 61, 61, 
        61, 61, 61, 61, 61, 
        61, 61, 61, 61, 61, 
        61, 61, 61, 61, 62, 
        62, 62, 62, 62, 62, 
        62, 62, 62, 62, 62, 
        62, 62, 62, 62, 62, 
        63, 63, 63, 63, 63, 
        63, 63, 63, 63, 63, 
        63, 63, 64};
    
    std::vector<double> ky {4, 5, 6, 7, 8, 
        9, 10, 11, 12, 13, 
        14, 15, 16, 17, 18, 
        19, 20, 21, 22, 23, 
        24, 25, 26, 27, 28, 
        29, 30, 31, 32, 33, 
        34, 35, 36, 37, 38, 
        39, 40, 41, 42, 43, 
        44, 45, 46, 47, 48, 
        49, 50, 51, 52, 53, 
        54, 55, 56, 57, 58, 
        59, 60, 61, 62, 63, 
        64, 4, 5, 6, 7, 
        8, 9, 10, 11, 12, 
        13, 14, 15, 16, 17, 
        18, 19, 20, 21, 22, 
        23, 24, 25, 26, 27, 
        28, 29, 30, 31, 32, 
        33, 34, 35, 36, 37, 
        38, 39, 40, 41, 42, 
        43, 44, 45, 46, 47, 
        48, 49, 50, 51, 52, 
        53, 54, 55, 56, 57, 
        58, 59, 60, 61, 62, 
        63, 4, 5, 6, 7, 
        8, 9, 10, 11, 12, 
        13, 14, 15, 16, 17, 
        18, 19, 20, 21, 22, 
        23, 24, 25, 26, 27, 
        28, 29, 30, 31, 32, 
        33, 34, 35, 36, 37, 
        38, 39, 40, 41, 42, 
        43, 44, 45, 46, 47, 
        48, 49, 50, 51, 52, 
        53, 54, 55, 56, 57, 
        58, 59, 60, 61, 62, 
        63, 3, 4, 5, 6, 
        7, 8, 9, 10, 11, 
        12, 13, 14, 15, 16, 
        17, 18, 19, 20, 21, 
        22, 23, 24, 25, 26, 
        27, 28, 29, 30, 31, 
        32, 33, 34, 35, 36, 
        37, 38, 39, 40, 41, 
        42, 43, 44, 45, 46, 
        47, 48, 49, 50, 51, 
        52, 53, 54, 55, 56, 
        57, 58, 59, 60, 61, 
        62, 63, 0, 1, 2, 
        3, 4, 5, 6, 7, 
        8, 9, 10, 11, 12, 
        13, 14, 15, 16, 17, 
        18, 19, 20, 21, 22, 
        23, 24, 25, 26, 27, 
        28, 29, 30, 31, 32, 
        33, 34, 35, 36, 37, 
        38, 39, 40, 41, 42, 
        43, 44, 45, 46, 47, 
        48, 49, 50, 51, 52, 
        53, 54, 55, 56, 57, 
        58, 59, 60, 61, 62, 
        63, 0, 1, 2, 3, 
        4, 5, 6, 7, 8, 
        9, 10, 11, 12, 13, 
        14, 15, 16, 17, 18, 
        19, 20, 21, 22, 23, 
        24, 25, 26, 27, 28, 
        29, 30, 31, 32, 33, 
        34, 35, 36, 37, 38, 
        39, 40, 41, 42, 43, 
        44, 45, 46, 47, 48, 
        49, 50, 51, 52, 53, 
        54, 55, 56, 57, 58, 
        59, 60, 61, 62, 63, 
        0, 1, 2, 3, 4, 
        5, 6, 7, 8, 9, 
        10, 11, 12, 13, 14, 
        15, 16, 17, 18, 19, 
        20, 21, 22, 23, 24, 
        25, 26, 27, 28, 29, 
        30, 31, 32, 33, 34, 
        35, 36, 37, 38, 39, 
        40, 41, 42, 43, 44, 
        45, 46, 47, 48, 49, 
        50, 51, 52, 53, 54, 
        55, 56, 57, 58, 59, 
        60, 61, 62, 63, 0, 
        1, 2, 3, 4, 5, 
        6, 7, 8, 9, 10, 
        11, 12, 13, 14, 15, 
        16, 17, 18, 19, 20, 
        21, 22, 23, 24, 25, 
        26, 27, 28, 29, 30, 
        31, 32, 33, 34, 35, 
        36, 37, 38, 39, 40, 
        41, 42, 43, 44, 45, 
        46, 47, 48, 49, 50, 
        51, 52, 53, 54, 55, 
        56, 57, 58, 59, 60, 
        61, 62, 63, 0, 1, 
        2, 3, 4, 5, 6, 
        7, 8, 9, 10, 11, 
        12, 13, 14, 15, 16, 
        17, 18, 19, 20, 21, 
        22, 23, 24, 25, 26, 
        27, 28, 29, 30, 31, 
        32, 33, 34, 35, 36, 
        37, 38, 39, 40, 41, 
        42, 43, 44, 45, 46, 
        47, 48, 49, 50, 51, 
        52, 53, 54, 55, 56, 
        57, 58, 59, 60, 61, 
        62, 63, 0, 1, 2, 
        3, 4, 5, 6, 7, 
        8, 9, 10, 11, 12, 
        13, 14, 15, 16, 17, 
        18, 19, 20, 21, 22, 
        23, 24, 25, 26, 27, 
        28, 29, 30, 31, 32, 
        33, 34, 35, 36, 37, 
        38, 39, 40, 41, 42, 
        43, 44, 45, 46, 47, 
        48, 49, 50, 51, 52, 
        53, 54, 55, 56, 57, 
        58, 59, 60, 61, 62, 
        63, 0, 1, 2, 3, 
        4, 5, 6, 7, 8, 
        9, 10, 11, 12, 13, 
        14, 15, 16, 17, 18, 
        19, 20, 21, 22, 23, 
        24, 25, 26, 27, 28, 
        29, 30, 31, 32, 33, 
        34, 35, 36, 37, 38, 
        39, 40, 41, 42, 43, 
        44, 45, 46, 47, 48, 
        49, 50, 51, 52, 53, 
        54, 55, 56, 57, 58, 
        59, 60, 61, 62, 63, 
        0, 1, 2, 3, 4, 
        5, 6, 7, 8, 9, 
        10, 11, 12, 13, 14, 
        15, 16, 17, 18, 19, 
        20, 21, 22, 23, 24, 
        25, 26, 27, 28, 29, 
        30, 31, 32, 33, 34, 
        35, 36, 37, 38, 39, 
        40, 41, 42, 43, 44, 
        45, 46, 47, 48, 49, 
        50, 51, 52, 53, 54, 
        55, 56, 57, 58, 59, 
        60, 61, 62, 63, 0, 
        1, 2, 3, 4, 5, 
        6, 7, 8, 9, 10, 
        11, 12, 13, 14, 15, 
        16, 17, 18, 19, 20, 
        21, 22, 23, 24, 25, 
        26, 27, 28, 29, 30, 
        31, 32, 33, 34, 35, 
        36, 37, 38, 39, 40, 
        41, 42, 43, 44, 45, 
        46, 47, 48, 49, 50, 
        51, 52, 53, 54, 55, 
        56, 57, 58, 59, 60, 
        61, 62, 0, 1, 2, 
        3, 4, 5, 6, 7, 
        8, 9, 10, 11, 12, 
        13, 14, 15, 16, 17, 
        18, 19, 20, 21, 22, 
        23, 24, 25, 26, 27, 
        28, 29, 30, 31, 32, 
        33, 34, 35, 36, 37, 
        38, 39, 40, 41, 42, 
        43, 44, 45, 46, 47, 
        48, 49, 50, 51, 52, 
        53, 54, 55, 56, 57, 
        58, 59, 60, 61, 62, 
        0, 1, 2, 3, 4, 
        5, 6, 7, 8, 9, 
        10, 11, 12, 13, 14, 
        15, 16, 17, 18, 19, 
        20, 21, 22, 23, 24, 
        25, 26, 27, 28, 29, 
        30, 31, 32, 33, 34, 
        35, 36, 37, 38, 39, 
        40, 41, 42, 43, 44, 
        45, 46, 47, 48, 49, 
        50, 51, 52, 53, 54, 
        55, 56, 57, 58, 59, 
        60, 61, 62, 0, 1, 
        2, 3, 4, 5, 6, 
        7, 8, 9, 10, 11, 
        12, 13, 14, 15, 16, 
        17, 18, 19, 20, 21, 
        22, 23, 24, 25, 26, 
        27, 28, 29, 30, 31, 
        32, 33, 34, 35, 36, 
        37, 38, 39, 40, 41, 
        42, 43, 44, 45, 46, 
        47, 48, 49, 50, 51, 
        52, 53, 54, 55, 56, 
        57, 58, 59, 60, 61, 
        62, 0, 1, 2, 3, 
        4, 5, 6, 7, 8, 
        9, 10, 11, 12, 13, 
        14, 15, 16, 17, 18, 
        19, 20, 21, 22, 23, 
        24, 25, 26, 27, 28, 
        29, 30, 31, 32, 33, 
        34, 35, 36, 37, 38, 
        39, 40, 41, 42, 43, 
        44, 45, 46, 47, 48, 
        49, 50, 51, 52, 53, 
        54, 55, 56, 57, 58, 
        59, 60, 61, 0, 1, 
        2, 3, 4, 5, 6, 
        7, 8, 9, 10, 11, 
        12, 13, 14, 15, 16, 
        17, 18, 19, 20, 21, 
        22, 23, 24, 25, 26, 
        27, 28, 29, 30, 31, 
        32, 33, 34, 35, 36, 
        37, 38, 39, 40, 41, 
        42, 43, 44, 45, 46, 
        47, 48, 49, 50, 51, 
        52, 53, 54, 55, 56, 
        57, 58, 59, 60, 61, 
        0, 1, 2, 3, 4, 
        5, 6, 7, 8, 9, 
        10, 11, 12, 13, 14, 
        15, 16, 17, 18, 19, 
        20, 21, 22, 23, 24, 
        25, 26, 27, 28, 29, 
        30, 31, 32, 33, 34, 
        35, 36, 37, 38, 39, 
        40, 41, 42, 43, 44, 
        45, 46, 47, 48, 49, 
        50, 51, 52, 53, 54, 
        55, 56, 57, 58, 59, 
        60, 61, 0, 1, 2, 
        3, 4, 5, 6, 7, 
        8, 9, 10, 11, 12, 
        13, 14, 15, 16, 17, 
        18, 19, 20, 21, 22, 
        23, 24, 25, 26, 27, 
        28, 29, 30, 31, 32, 
        33, 34, 35, 36, 37, 
        38, 39, 40, 41, 42, 
        43, 44, 45, 46, 47, 
        48, 49, 50, 51, 52, 
        53, 54, 55, 56, 57, 
        58, 59, 60, 61, 0, 
        1, 2, 3, 4, 5, 
        6, 7, 8, 9, 10, 
        11, 12, 13, 14, 15, 
        16, 17, 18, 19, 20, 
        21, 22, 23, 24, 25, 
        26, 27, 28, 29, 30, 
        31, 32, 33, 34, 35, 
        36, 37, 38, 39, 40, 
        41, 42, 43, 44, 45, 
        46, 47, 48, 49, 50, 
        51, 52, 53, 54, 55, 
        56, 57, 58, 59, 60, 
        0, 1, 2, 3, 4, 
        5, 6, 7, 8, 9, 
        10, 11, 12, 13, 14, 
        15, 16, 17, 18, 19, 
        20, 21, 22, 23, 24, 
        25, 26, 27, 28, 29, 
        30, 31, 32, 33, 34, 
        35, 36, 37, 38, 39, 
        40, 41, 42, 43, 44, 
        45, 46, 47, 48, 49, 
        50, 51, 52, 53, 54, 
        55, 56, 57, 58, 59, 
        60, 0, 1, 2, 3, 
        4, 5, 6, 7, 8, 
        9, 10, 11, 12, 13, 
        14, 15, 16, 17, 18, 
        19, 20, 21, 22, 23, 
        24, 25, 26, 27, 28, 
        29, 30, 31, 32, 33, 
        34, 35, 36, 37, 38, 
        39, 40, 41, 42, 43, 
        44, 45, 46, 47, 48, 
        49, 50, 51, 52, 53, 
        54, 55, 56, 57, 58, 
        59, 60, 0, 1, 2, 
        3, 4, 5, 6, 7, 
        8, 9, 10, 11, 12, 
        13, 14, 15, 16, 17, 
        18, 19, 20, 21, 22, 
        23, 24, 25, 26, 27, 
        28, 29, 30, 31, 32, 
        33, 34, 35, 36, 37, 
        38, 39, 40, 41, 42, 
        43, 44, 45, 46, 47, 
        48, 49, 50, 51, 52, 
        53, 54, 55, 56, 57, 
        58, 59, 0, 1, 2, 
        3, 4, 5, 6, 7, 
        8, 9, 10, 11, 12, 
        13, 14, 15, 16, 17, 
        18, 19, 20, 21, 22, 
        23, 24, 25, 26, 27, 
        28, 29, 30, 31, 32, 
        33, 34, 35, 36, 37, 
        38, 39, 40, 41, 42, 
        43, 44, 45, 46, 47, 
        48, 49, 50, 51, 52, 
        53, 54, 55, 56, 57, 
        58, 59, 0, 1, 2, 
        3, 4, 5, 6, 7, 
        8, 9, 10, 11, 12, 
        13, 14, 15, 16, 17, 
        18, 19, 20, 21, 22, 
        23, 24, 25, 26, 27, 
        28, 29, 30, 31, 32, 
        33, 34, 35, 36, 37, 
        38, 39, 40, 41, 42, 
        43, 44, 45, 46, 47, 
        48, 49, 50, 51, 52, 
        53, 54, 55, 56, 57, 
        58, 0, 1, 2, 3, 
        4, 5, 6, 7, 8, 
        9, 10, 11, 12, 13, 
        14, 15, 16, 17, 18, 
        19, 20, 21, 22, 23, 
        24, 25, 26, 27, 28, 
        29, 30, 31, 32, 33, 
        34, 35, 36, 37, 38, 
        39, 40, 41, 42, 43, 
        44, 45, 46, 47, 48, 
        49, 50, 51, 52, 53, 
        54, 55, 56, 57, 58, 
        0, 1, 2, 3, 4, 
        5, 6, 7, 8, 9, 
        10, 11, 12, 13, 14, 
        15, 16, 17, 18, 19, 
        20, 21, 22, 23, 24, 
        25, 26, 27, 28, 29, 
        30, 31, 32, 33, 34, 
        35, 36, 37, 38, 39, 
        40, 41, 42, 43, 44, 
        45, 46, 47, 48, 49, 
        50, 51, 52, 53, 54, 
        55, 56, 57, 58, 0, 
        1, 2, 3, 4, 5, 
        6, 7, 8, 9, 10, 
        11, 12, 13, 14, 15, 
        16, 17, 18, 19, 20, 
        21, 22, 23, 24, 25, 
        26, 27, 28, 29, 30, 
        31, 32, 33, 34, 35, 
        36, 37, 38, 39, 40, 
        41, 42, 43, 44, 45, 
        46, 47, 48, 49, 50, 
        51, 52, 53, 54, 55, 
        56, 57, 0, 1, 2, 
        3, 4, 5, 6, 7, 
        8, 9, 10, 11, 12, 
        13, 14, 15, 16, 17, 
        18, 19, 20, 21, 22, 
        23, 24, 25, 26, 27, 
        28, 29, 30, 31, 32, 
        33, 34, 35, 36, 37, 
        38, 39, 40, 41, 42, 
        43, 44, 45, 46, 47, 
        48, 49, 50, 51, 52, 
        53, 54, 55, 56, 57, 
        0, 1, 2, 3, 4, 
        5, 6, 7, 8, 9, 
        10, 11, 12, 13, 14, 
        15, 16, 17, 18, 19, 
        20, 21, 22, 23, 24, 
        25, 26, 27, 28, 29, 
        30, 31, 32, 33, 34, 
        35, 36, 37, 38, 39, 
        40, 41, 42, 43, 44, 
        45, 46, 47, 48, 49, 
        50, 51, 52, 53, 54, 
        55, 56, 0, 1, 2, 
        3, 4, 5, 6, 7, 
        8, 9, 10, 11, 12, 
        13, 14, 15, 16, 17, 
        18, 19, 20, 21, 22, 
        23, 24, 25, 26, 27, 
        28, 29, 30, 31, 32, 
        33, 34, 35, 36, 37, 
        38, 39, 40, 41, 42, 
        43, 44, 45, 46, 47, 
        48, 49, 50, 51, 52, 
        53, 54, 55, 0, 1, 
        2, 3, 4, 5, 6, 
        7, 8, 9, 10, 11, 
        12, 13, 14, 15, 16, 
        17, 18, 19, 20, 21, 
        22, 23, 24, 25, 26, 
        27, 28, 29, 30, 31, 
        32, 33, 34, 35, 36, 
        37, 38, 39, 40, 41, 
        42, 43, 44, 45, 46, 
        47, 48, 49, 50, 51, 
        52, 53, 54, 55, 0, 
        1, 2, 3, 4, 5, 
        6, 7, 8, 9, 10, 
        11, 12, 13, 14, 15, 
        16, 17, 18, 19, 20, 
        21, 22, 23, 24, 25, 
        26, 27, 28, 29, 30, 
        31, 32, 33, 34, 35, 
        36, 37, 38, 39, 40, 
        41, 42, 43, 44, 45, 
        46, 47, 48, 49, 50, 
        51, 52, 53, 54, 0, 
        1, 2, 3, 4, 5, 
        6, 7, 8, 9, 10, 
        11, 12, 13, 14, 15, 
        16, 17, 18, 19, 20, 
        21, 22, 23, 24, 25, 
        26, 27, 28, 29, 30, 
        31, 32, 33, 34, 35, 
        36, 37, 38, 39, 40, 
        41, 42, 43, 44, 45, 
        46, 47, 48, 49, 50, 
        51, 52, 53, 54, 0, 
        1, 2, 3, 4, 5, 
        6, 7, 8, 9, 10, 
        11, 12, 13, 14, 15, 
        16, 17, 18, 19, 20, 
        21, 22, 23, 24, 25, 
        26, 27, 28, 29, 30, 
        31, 32, 33, 34, 35, 
        36, 37, 38, 39, 40, 
        41, 42, 43, 44, 45, 
        46, 47, 48, 49, 50, 
        51, 52, 53, 0, 1, 
        2, 3, 4, 5, 6, 
        7, 8, 9, 10, 11, 
        12, 13, 14, 15, 16, 
        17, 18, 19, 20, 21, 
        22, 23, 24, 25, 26, 
        27, 28, 29, 30, 31, 
        32, 33, 34, 35, 36, 
        37, 38, 39, 40, 41, 
        42, 43, 44, 45, 46, 
        47, 48, 49, 50, 51, 
        52, 0, 1, 2, 3, 
        4, 5, 6, 7, 8, 
        9, 10, 11, 12, 13, 
        14, 15, 16, 17, 18, 
        19, 20, 21, 22, 23, 
        24, 25, 26, 27, 28, 
        29, 30, 31, 32, 33, 
        34, 35, 36, 37, 38, 
        39, 40, 41, 42, 43, 
        44, 45, 46, 47, 48, 
        49, 50, 51, 52, 0, 
        1, 2, 3, 4, 5, 
        6, 7, 8, 9, 10, 
        11, 12, 13, 14, 15, 
        16, 17, 18, 19, 20, 
        21, 22, 23, 24, 25, 
        26, 27, 28, 29, 30, 
        31, 32, 33, 34, 35, 
        36, 37, 38, 39, 40, 
        41, 42, 43, 44, 45, 
        46, 47, 48, 49, 50, 
        51, 0, 1, 2, 3, 
        4, 5, 6, 7, 8, 
        9, 10, 11, 12, 13, 
        14, 15, 16, 17, 18, 
        19, 20, 21, 22, 23, 
        24, 25, 26, 27, 28, 
        29, 30, 31, 32, 33, 
        34, 35, 36, 37, 38, 
        39, 40, 41, 42, 43, 
        44, 45, 46, 47, 48, 
        49, 50, 0, 1, 2, 
        3, 4, 5, 6, 7, 
        8, 9, 10, 11, 12, 
        13, 14, 15, 16, 17, 
        18, 19, 20, 21, 22, 
        23, 24, 25, 26, 27, 
        28, 29, 30, 31, 32, 
        33, 34, 35, 36, 37, 
        38, 39, 40, 41, 42, 
        43, 44, 45, 46, 47, 
        48, 49, 0, 1, 2, 
        3, 4, 5, 6, 7, 
        8, 9, 10, 11, 12, 
        13, 14, 15, 16, 17, 
        18, 19, 20, 21, 22, 
        23, 24, 25, 26, 27, 
        28, 29, 30, 31, 32, 
        33, 34, 35, 36, 37, 
        38, 39, 40, 41, 42, 
        43, 44, 45, 46, 47, 
        48, 49, 0, 1, 2, 
        3, 4, 5, 6, 7, 
        8, 9, 10, 11, 12, 
        13, 14, 15, 16, 17, 
        18, 19, 20, 21, 22, 
        23, 24, 25, 26, 27, 
        28, 29, 30, 31, 32, 
        33, 34, 35, 36, 37, 
        38, 39, 40, 41, 42, 
        43, 44, 45, 46, 47, 
        48, 0, 1, 2, 3, 
        4, 5, 6, 7, 8, 
        9, 10, 11, 12, 13, 
        14, 15, 16, 17, 18, 
        19, 20, 21, 22, 23, 
        24, 25, 26, 27, 28, 
        29, 30, 31, 32, 33, 
        34, 35, 36, 37, 38, 
        39, 40, 41, 42, 43, 
        44, 45, 46, 47, 0, 
        1, 2, 3, 4, 5, 
        6, 7, 8, 9, 10, 
        11, 12, 13, 14, 15, 
        16, 17, 18, 19, 20, 
        21, 22, 23, 24, 25, 
        26, 27, 28, 29, 30, 
        31, 32, 33, 34, 35, 
        36, 37, 38, 39, 40, 
        41, 42, 43, 44, 45, 
        46, 0, 1, 2, 3, 
        4, 5, 6, 7, 8, 
        9, 10, 11, 12, 13, 
        14, 15, 16, 17, 18, 
        19, 20, 21, 22, 23, 
        24, 25, 26, 27, 28, 
        29, 30, 31, 32, 33, 
        34, 35, 36, 37, 38, 
        39, 40, 41, 42, 43, 
        44, 45, 0, 1, 2, 
        3, 4, 5, 6, 7, 
        8, 9, 10, 11, 12, 
        13, 14, 15, 16, 17, 
        18, 19, 20, 21, 22, 
        23, 24, 25, 26, 27, 
        28, 29, 30, 31, 32, 
        33, 34, 35, 36, 37, 
        38, 39, 40, 41, 42, 
        43, 44, 0, 1, 2, 
        3, 4, 5, 6, 7, 
        8, 9, 10, 11, 12, 
        13, 14, 15, 16, 17, 
        18, 19, 20, 21, 22, 
        23, 24, 25, 26, 27, 
        28, 29, 30, 31, 32, 
        33, 34, 35, 36, 37, 
        38, 39, 40, 41, 42, 
        43, 0, 1, 2, 3, 
        4, 5, 6, 7, 8, 
        9, 10, 11, 12, 13, 
        14, 15, 16, 17, 18, 
        19, 20, 21, 22, 23, 
        24, 25, 26, 27, 28, 
        29, 30, 31, 32, 33, 
        34, 35, 36, 37, 38, 
        39, 40, 41, 42, 0, 
        1, 2, 3, 4, 5, 
        6, 7, 8, 9, 10, 
        11, 12, 13, 14, 15, 
        16, 17, 18, 19, 20, 
        21, 22, 23, 24, 25, 
        26, 27, 28, 29, 30, 
        31, 32, 33, 34, 35, 
        36, 37, 38, 39, 40, 
        41, 0, 1, 2, 3, 
        4, 5, 6, 7, 8, 
        9, 10, 11, 12, 13, 
        14, 15, 16, 17, 18, 
        19, 20, 21, 22, 23, 
        24, 25, 26, 27, 28, 
        29, 30, 31, 32, 33, 
        34, 35, 36, 37, 38, 
        39, 0, 1, 2, 3, 
        4, 5, 6, 7, 8, 
        9, 10, 11, 12, 13, 
        14, 15, 16, 17, 18, 
        19, 20, 21, 22, 23, 
        24, 25, 26, 27, 28, 
        29, 30, 31, 32, 33, 
        34, 35, 36, 37, 38, 
        0, 1, 2, 3, 4, 
        5, 6, 7, 8, 9, 
        10, 11, 12, 13, 14, 
        15, 16, 17, 18, 19, 
        20, 21, 22, 23, 24, 
        25, 26, 27, 28, 29, 
        30, 31, 32, 33, 34, 
        35, 36, 37, 0, 1, 
        2, 3, 4, 5, 6, 
        7, 8, 9, 10, 11, 
        12, 13, 14, 15, 16, 
        17, 18, 19, 20, 21, 
        22, 23, 24, 25, 26, 
        27, 28, 29, 30, 31, 
        32, 33, 34, 35, 0, 
        1, 2, 3, 4, 5, 
        6, 7, 8, 9, 10, 
        11, 12, 13, 14, 15, 
        16, 17, 18, 19, 20, 
        21, 22, 23, 24, 25, 
        26, 27, 28, 29, 30, 
        31, 32, 33, 34, 0, 
        1, 2, 3, 4, 5, 
        6, 7, 8, 9, 10, 
        11, 12, 13, 14, 15, 
        16, 17, 18, 19, 20, 
        21, 22, 23, 24, 25, 
        26, 27, 28, 29, 30, 
        31, 32, 0, 1, 2, 
        3, 4, 5, 6, 7, 
        8, 9, 10, 11, 12, 
        13, 14, 15, 16, 17, 
        18, 19, 20, 21, 22, 
        23, 24, 25, 26, 27, 
        28, 29, 30, 0, 1, 
        2, 3, 4, 5, 6, 
        7, 8, 9, 10, 11, 
        12, 13, 14, 15, 16, 
        17, 18, 19, 20, 21, 
        22, 23, 24, 25, 26, 
        27, 28, 29, 0, 1, 
        2, 3, 4, 5, 6, 
        7, 8, 9, 10, 11, 
        12, 13, 14, 15, 16, 
        17, 18, 19, 20, 21, 
        22, 23, 24, 25, 26, 
        27, 0, 1, 2, 3, 
        4, 5, 6, 7, 8, 
        9, 10, 11, 12, 13, 
        14, 15, 16, 17, 18, 
        19, 20, 21, 22, 23, 
        24, 0, 1, 2, 3, 
        4, 5, 6, 7, 8, 
        9, 10, 11, 12, 13, 
        14, 15, 16, 17, 18, 
        19, 20, 21, 22, 0, 
        1, 2, 3, 4, 5, 
        6, 7, 8, 9, 10, 
        11, 12, 13, 14, 15, 
        16, 17, 18, 19, 0, 
        1, 2, 3, 4, 5, 
        6, 7, 8, 9, 10, 
        11, 12, 13, 14, 15, 
        0, 1, 2, 3, 4, 
        5, 6, 7, 8, 9, 
        10, 11, 0};
    
    std::vector<double> ar_1d {0.62944737278635787, 0.81158387415123845, 0.74602636741298789, 0.82675171227803879, 0.26471849245081902, 
        0.80491919000118095, 0.44300356226590321, 0.093763038409967692, 0.9150136708685952, 0.92977707039855306, 
        0.68477383664490343, 0.94118556352123139, 0.91433389648589114, 0.029248702554317552, 0.60056093777760022, 
        0.71622732274556933, 0.15647743474745002, 0.83147105037813418, 0.58441465911910884, 0.91898485278580599, 
        0.31148139831317367, 0.92857664285162089, 0.69825861173755421, 0.8679864955151011, 0.35747030971554694, 
        0.5154802611566669, 0.48626493624983236, 0.21554596093166367, 0.31095578035511329, 0.65762662437687647, 
        0.41209217603921755, 0.93633430724515865, 0.44615403007822008, 0.90765721873769212, 0.80573643752830493, 
        0.64691565665458528, 0.38965724595163409, 0.36580103987827894, 0.90044409767670985, 0.93110783899418248, 
        0.12251128068720352, 0.23688308581398321, 0.53103357629800474, 0.59039980227412636, 0.62625479089124281, 
        0.02047120842353789, 0.10882759857820101, 0.29262602022252926, 0.41872966171614512, 0.50937336396472177, 
        0.44794984600284327, 0.35940535370734961, 0.31019600794768132, 0.67477652961073886, 0.76200463688324671, 
        0.0032718960357140947, 0.91948791703216215, 0.31922854666773359, 0.17053550195955469, 0.55237612101772604, 
        0.5025341186113057, 0.48980976908146179, 0.011914103330284753, 0.39815344531337193, 0.78180650507159699, 
        0.91858285041088861, 0.094431059927606142, 0.72275111434264172, 0.70141198888188505, 0.48498349175252708, 
        0.68143451196732507, 0.4914356420569379, 0.62856965213763272, 0.5129500625500214, 0.85852724637445554, 
        0.30003246803038253, 0.60680949913758364, 0.49783228404793789, 0.23208935229327832, 0.053422302194541471, 
        0.29668098587400649, 0.66165725579258172, 0.17052818230544853, 0.099447216582279063, 0.83438732765962009, 
        0.42832196235925291, 0.51440045822144254, 0.50745818855699065, 0.23910830604928668, 0.13564328145044224, 
        0.84829142087387277, 0.8920997626667857, 0.061595106017945378, 0.55833446020402233, 0.86802136845836597, 
        0.74018758305253973, 0.13764732174438543, 0.061218717883588347, 0.97619586099751721, 0.32575471120223698, 
        0.6756353836135145, 0.58856908136781394, 0.37756991591039024, 0.057066271012425451, 0.66870254100043813, 
        0.20396388280327304, 0.47405743091971142, 0.3081581969535645, 0.37842900628001552, 0.49630318564741893, 
        0.098916802995004494, 0.83235724400613487, 0.54204606256636234, 0.82667472300333911, 0.69524396206155403, 
        0.65163395497909482, 0.07668487052011419, 0.99226943325377093, 0.84364894249363265, 0.11464346044910734, 
        0.78669445963883122, 0.92379616171010737, 0.99073155173186511, 0.54982092942300476, 0.63460644130686594, 
        0.73738941072701936, 0.83112830897817935, 0.20043470180220702, 0.48025919429869157, 0.60013696044861509, 
        0.13717234507291076, 0.82129518885904584, 0.63630594339429503, 0.47239416695601988, 0.70892203923056596, 
        0.72786288258267251, 0.7385844152801786, 0.15940917473114036, 0.099720403672663993, 0.71009040355254638, 
        0.70606223544378732, 0.24411026297013194, 0.29809523821545825, 0.026499079734106701, 0.19638393249611652, 
        0.84806661661831617, 0.52016769289268394, 0.75336213032966892, 0.63218442343516656, 0.52009494867019446, 
        0.165465861831261, 0.90069113934851575, 0.80543221983056212, 0.88957437944329198, 0.018271815063840124, 
        0.021494723199962262, 0.32456118035724568, 0.8001076928353239, 0.26150643775956994, 0.77759448941242515, 
        0.56050413664227583, 0.22052232607749311, 0.51661742817233458, 0.19217570882377055, 0.80709094966322281, 
        0.73605341478732988, 0.88410118155097028, 0.91226908045960453, 0.15041719015693111, 0.88044091410568837, 
        0.53044017325518733, 0.29368285755585788, 0.6423880803959181, 0.9691931246968899, 0.91395239668438433, 
        0.66201994107459128, 0.29823094991290411, 0.46344477131734063, 0.29549192627261345, 0.098152587138110148, 
        0.094017784572689944, 0.40735838878445363, 0.48938561414831239, 0.62208996993491095, 0.37355086673062998, 
        0.63297768852546055, 0.26303080701932702, 0.25123712145938071, 0.56045487030275365, 0.83774846226842947, 
        0.85877194193745998, 0.55142535721680463, 0.026416735193655283, 0.12828282283816184, 0.10643250114038749, 
        0.38730105596688524, 0.01701731076225399, 0.021543128344219342, 0.63525541664452412, 0.58966283376690609, 
        0.28863626038738333, 0.24278123467946311, 0.62316091656495431, 0.065651177598909705, 0.29854579284623339, 
        0.87800312399977365, 0.75188562298596762, 0.10031268579684438, 0.24495017200245495, 0.17408940906283354, 
        0.58451541453394307, 0.39750733944101868, 0.058153302964818687, 0.53902367957688302, 0.68861758539077811, 
        0.61047142086590145, 0.54815643805520242, 0.65858390570428282, 0.54467140436689299, 0.12860263179220177, 
        0.37779542669917432, 0.84675928420648772, 0.13958521734083207, 0.63036735975172786, 0.80976193735978574, 
        0.95949675671217038, 0.12226005374779358, 0.77776155311880246, 0.48387060817586613, 0.18256030777489585, 
        0.18979214801722866, 0.47557650443830912, 0.20568617876416595, 0.42243156086736589, 0.5565065319655198, 
        0.76516469828838818, 0.40664825356334622, 0.3624433961482354, 0.15166648057238552, 0.015716569322236307, 
        0.82896840581991205, 0.4750355306033347, 0.60202924553947756, 0.94155944487570742, 0.85770827895608925, 
        0.4606617257109058, 0.022782052392841701, 0.15705012204687785, 0.52543284045695704, 0.082302343640137776, 
        0.92617707857382592, 0.093611437477935988, 0.042271661608002997, 0.53681122658295233, 0.022204512159666123, 
        0.24812017634737904, 0.35827108173149536, 0.20896956866281391, 0.2651267029110469, 0.97596400632326552, 
        0.92452226752089572, 0.77033601640495064, 0.82657365527847793, 0.59236774717042406, 0.80257544268885139, 
        0.47625763225856788, 0.32928632007440695, 0.35945590275467598, 0.72689372528926066, 0.44245499716348036, 
        0.78647627678551713, 0.30751469733711922, 0.011652126721459766, 0.55810344646255028, 0.43007415680138816, 
        0.80744112111263266, 0.78184500866157847, 0.33167389452500751, 0.39749166466958896, 0.60438034662814166, 
        0.93891810739072668, 0.48814852073492476, 4.4871180401750976e-05, 0.040155717707879068, 0.80944447613472548, 
        0.21973329684511689, 0.23533277917690931, 0.71888461129242454, 0.61097884905937128, 0.1534430312293702, 
        0.63415506117017206, 0.52013597886256524, 0.77302386615220264, 0.94265169507178781, 0.020197222975552265, 
        0.6641457086354865, 0.9573612992823175, 0.42538894335782818, 0.00094324830968606221, 0.057823250916121394, 
        0.88076226484072162, 0.36394380829812634, 0.91513772499851664, 0.85710907079871523, 0.043299684928567395, 
        0.80653994843826604, 0.63629710771924941, 0.63509418415857266, 0.4448791847336846, 0.7002691150440663, 
        0.31921050581661436, 0.037189885021076341, 0.94594910952772504, 0.29798298542471224, 0.60066115070480297, 
        0.092404582546161018, 0.13521699243307661, 0.65062759080409127, 0.83306037028217195, 0.73365798478567656, 
        0.6532227737619889, 0.21812439535252892, 0.66275948567813914, 0.60672878320488044, 0.87905764166021272, 
        0.2014844587728486, 0.05375166101659179, 0.16640106413842615, 0.31371978194741446, 0.25594671838020844, 
        0.41603184007657013, 0.13669765950255952, 0.9690257487279621, 0.96812744875830758, 0.66566318017068804, 
        0.78756731014267234, 0.25518051988892609, 0.60376319491405073, 0.02062472396795223, 0.32101317321848466, 
        0.9032609295554539, 0.84066407967312751, 0.89464600463841482, 0.4757161910339931, 0.46176114720288819, 
        0.15432876998238432, 0.095741802429689438, 0.88547396855386862, 0.1645117913666756, 0.96610493293971222, 
        0.39709010257586907, 0.40219751180185259, 0.33267770316885126, 0.078252930085713324, 0.39621104036061672, 
        0.33305582680517376, 0.64373509119932448, 0.74397120055965482, 0.99816078952272136, 0.65775786728713581, 
        0.93479835893894392, 0.12239958541932028, 0.76373300090361984, 0.33835060906878756, 0.61913346564009175, 
        0.2621669078722102, 0.078548125479176889, 0.96327590194149937, 0.68719009554687305, 0.71104561169182268, 
        0.28952907374017589, 0.24745557944233676, 0.61815260952739393, 0.14349401404122797, 0.035955877936287406, 
        0.7587767734056754, 0.17901496939011796, 0.54762464049464854, 0.23076175126117837, 0.16597276549534778, 
        0.49638775505537436, 0.41911867144604131, 0.2341817687864467, 0.46943818037994123, 0.64875253337767003, 
        0.96532679944390054, 0.46049758453519529, 0.31224599177003376, 0.16813866655690402, 0.7844619695125139, 
        0.81261630129946583, 0.75930744896380964, 0.6355211187412837, 0.47854400188907076, 0.18871250132866169, 
        0.95497481451953647, 0.14948135957173014, 0.3745622263587689, 0.67703051137650005, 0.64246762649526357, 
        0.15422862179983099, 0.81154132222453068, 0.19704733751348225, 0.058151487283331971, 0.39189862660321584, 
        0.39977569985658312, 0.27706151654367583, 0.93279232786714106, 0.86238780176389751, 0.36080052963900777, 
        0.061728561388253, 0.30889141551413268, 0.18476160591769486, 0.63996244556388127, 0.43671788641176734, 
        0.93729866046218735, 0.062667813131348948, 0.34970863635887994, 0.78874159334195615, 0.22191731749240118, 
        0.55760448364818505, 0.15309416207452342, 0.81835342842512104, 0.46705701844185521, 0.69268656481738677, 
        0.43798939493225819, 0.11982972199655828, 0.054285483521304467, 0.085151268624652054, 0.75074319720836979, 
        0.036104216722208315, 0.88724524909677571, 0.27541819614434804, 0.91538787968316648, 0.51858592903967948, 
        0.35224460772750366, 0.42187085665104651, 0.34361633082843035, 0.39028099910347436, 0.86401446305997887, 
        0.49041968680598935, 0.5519199383515625, 0.33566545402743309, 0.68878431305440913, 0.31107517739791568, 
        0.56103930546271585, 0.35066413149399911, 0.98656937136304501, 0.20434097516359051, 0.22645761095803074, 
        0.83198248826285059, 0.99769788574178553, 0.075101681515342689, 0.15130192036924961, 0.078167267942071961, 
        0.54031945721721852, 0.35505638562644282, 0.56947858952148311, 0.057285692578776803, 0.92847453346176412, 
        0.64825116863293819, 0.4435160667822049, 0.053028014069359308, 0.69455759912353643, 0.31775078590178141, 
        0.21477842753669485, 0.61650948907640446, 0.47685367995388317, 0.5143008033636618, 0.83484868409876478, 
        0.46187682662796337, 0.53100003324287681, 0.62267604641701801, 0.42500365386773753, 0.81777307262693011, 
        0.15241876132601395, 0.36672648658930607, 0.093186229180645563, 0.14854231625762426, 0.28888556286267297, 
        0.29523526034536873, 0.35803350818640389, 0.27157342102816728, 0.89034822621880272, 0.5821301551479543, 
        0.41856340542108961, 0.52753884601240686, 0.76120750440538898, 0.21460788137126929, 0.099724606068207988, 
        0.082549012702264335, 0.32388950381130388, 0.54057102960732029, 0.29956397311778993, 0.32401919671826906, 
        0.16768282006040702, 0.68385830538261794, 0.66583363815043151, 0.48711801554170542, 0.22692147362575055, 
        0.16449832905445416, 0.081478674248819383, 0.7398820647160147, 0.47044194704873998, 0.3638518490378817, 
        0.76157091789161768, 0.87965894068984118, 0.29110374994504706, 0.041073550102224354, 0.27863392208021676, 
        0.089432221053525662, 0.29462296058625537, 0.087771867999278053, 0.44209324115962279, 0.044990611554204252, 
        0.98740924824170406, 0.56264673520073227, 0.78840345349954366, 0.78060507095361165, 0.8728172580497886, 
        0.19084000828474834, 0.10325417586700936, 0.2683676463236575, 0.52700928169762684, 0.25579275922833755, 
        0.54396077110849017, 0.86570714055763909, 0.9454817080060276, 0.61594330114445017, 0.7222515943416894, 
        0.39253267416598958, 0.81235994645026888, 0.050808807718671467, 0.060688436785726774, 0.7222796227866648, 
        0.030293332895796166, 0.21308727756946721, 0.34286227934805202, 0.48251588690841296, 0.040104934780773593, 
        0.30457465744494971, 0.70000549233663389, 0.1721841344629238, 0.47570936454438528, 0.91109181544352302, 
        0.5098665344623583, 0.51442928435807622, 0.11519537399611335, 0.3755921702402143, 0.28154357919627881, 
        0.47268014860240348, 0.21058504944247369, 0.36683173393595681, 0.40809486066853196, 0.11538917323325837, 
        0.96084475289336257, 0.33828423957185794, 0.15138100633372509, 0.45945915313586938, 0.60589240380908849, 
        0.64344236992261972, 0.14015718123346721, 0.77554190851270755, 0.21763400907767338, 0.53822877477659126, 
        0.20641696597276682, 0.61702819177469026, 0.51015419801416706, 0.24520891032979431, 0.56796216807721134, 
        0.58081443593382676, 0.8986078236995938, 0.34486913184958956, 0.34252874090347962, 0.12271003482608878, 
        0.66700119117795054, 0.53770850485922983, 0.66549290901055658, 0.72396095740414435, 0.97974430726300787, 
        0.028846913011408848, 0.76856204625391067, 0.17605211061699499, 0.69049530268791059, 0.6002743542850959, 
        0.18609032572218642, 0.49741143643138286, 0.65116763157231117, 0.57992605988906165, 0.36295150920201569, 
        0.068128254741452787, 0.8200986424588379, 0.7765885116135931, 0.72741490212340265, 0.35730460960037647, 
        0.0096459618206787834, 0.62057918796483968, 0.0099883500195583341, 0.70478355604662268, 0.89005170618762364, 
        0.70142534857801486, 0.12111905470976958, 0.85921773351332642, 0.39333440111045515, 0.16558193035168012, 
        0.63079442295484234, 0.75802780919435553, 0.97782323215917821, 0.99895524928611046, 0.73087718202604912, 
        0.22513293896799746, 0.97990041141766171, 0.055360138676884718, 0.040953229579562489, 0.60269521104390456, 
        0.54431412858791628, 0.0038114176072205908, 0.80170497706400967, 0.14932243826037528, 0.69035637010807327, 
        0.47728058399080364, 0.17197407165295164, 0.50653094802805021, 0.33283243463893597, 0.83303437279475467, 
        0.25191957034316648, 0.32188911589468461, 0.45950371063444218, 0.7815042326506445, 0.96460644576721277, 
        0.5380581706717924, 0.16289297575079553, 0.85662612462837595, 0.16018073151688328, 0.96603412332547745, 
        0.75828085780288368, 0.72542143739933929, 0.031406977575795025, 0.68971134915252641, 0.58118983195813056, 
        0.10458268307755003, 0.25976677012884242, 0.93601796847486618, 0.22942683823428101, 0.27517707545389469, 
        0.90093484191587758, 0.020860021645356497, 0.61497920787585048, 0.7538325049081096, 0.5890116581846403, 
        0.70697017877022028, 0.62185565105477236, 0.91469517817771329, 0.27039583371976406, 0.43626628823914015, 
        0.077193356090680787, 0.390326078888664, 0.0017679730348210132, 0.071602111502226506, 0.10963366940791541, 
        0.75213544480385974, 0.019285413063963697, 0.70599631068163249, 0.74785481172346668, 0.45941133541460477, 
        0.58307728249737245, 0.12995914147640231, 0.28062365032551617, 0.16594209671422799, 0.58804896893551328, 
        0.89586624258633729, 0.83585758580454828, 0.78858114683655711, 0.7159177561920036, 0.66707911824715849, 
        0.24191728787061662, 0.14741952968239636, 0.89584421942826076, 0.86240276921649928, 0.45732336335654233, 
        0.47568330759518029, 0.87319099861436356, 0.72088112607646382, 0.86881023792242673, 0.96879662448194326, 
        0.71787763336773236, 0.57111797853006108, 0.026754837175149815, 0.64479507898826993, 0.20282100652831381, 
        0.73213749802405759, 0.93822090251009693, 0.87828341213909589, 0.39738787082721561, 0.40893233104928717, 
        0.33412743632764941, 0.065863625942295112, 0.29639681293231424, 0.94954363701392741, 0.68441322483866895, 
        0.11806508997739007, 0.70819989854688625, 0.30424161134547778, 0.10794670388979433, 0.89152103111774084, 
        0.6457849324217817, 0.32561612392194883, 0.33834200959339089, 0.79697227566859907, 0.76368960310657874, 
        0.97683585756996116, 0.079964198075858306, 0.4138348386455255, 0.99898324019540929, 0.42430131036972552, 
        0.17095492221378317, 0.070320116749725381, 0.5279141569839132, 0.63640807781534159, 0.79955691960901687, 
        0.64376609222646719, 0.28073017303583936, 0.88659062186341764, 0.043771347322558496, 0.32830205064615092, 
        0.64866194064867844, 0.58210665201372924, 0.81030711800892807, 0.35078235467249375, 0.063063600192005298, 
        0.8242649484792457, 0.791976850441243, 0.49109214740343465, 0.47253491119327706, 0.12372285056327481, 
        0.63161180496894698, 0.19442270067570977, 0.40012601982042217, 0.73175413434263525, 0.57479693328231418, 
        0.78988335088162875, 0.85709437442591097, 0.51502688212656111, 0.89249121558257238, 0.11655588587115218, 
        0.97343359906549498, 0.79438270194714411, 0.60668361726473607, 0.81325896648981399, 0.38526620082415985, 
        0.087884666312517012, 0.79666121275048973, 0.99077945531018341, 0.33581433309500186, 0.40530636822415667, 
        0.87590955736073473, 0.40351205622447228, 0.90729747020363849, 0.010856284915405912, 0.5228517733802267, 
        0.26213999842718705, 0.82021669844781653, 0.83827515373937267, 0.55448107309788641, 0.81026948844714175, 
        0.067543903534000416, 0.78169157591508287, 0.65161771571119664, 0.32380456239565625, 0.41205389394703196, 
        0.49262685540735873, 0.97932676331320723, 0.90310532149355582, 0.33583224314724736, 0.20693596766153899, 
        0.052204931591122383, 0.45941889644645539, 0.41450697063084307, 0.56275410359855438, 0.42404604877165752, 
        0.38506397277303761, 0.11333966992802513, 0.20695841483681487, 0.87681866589207069, 0.56035106298234716, 
        0.3248322718959098, 0.21573181452589152, 0.48250809900443592, 0.79037351605300077, 0.74422324043400967, 
        0.099080214030395908, 0.029541182830081336, 0.78095135836887652, 0.59792055762575869, 0.46868216739193946, 
        0.8973362277752579, 0.85422940180204776, 0.82294508065055916, 0.59670172822790479, 0.88601627914140546, 
        0.36743114481671624, 0.73583408857287447, 0.44544907931353173, 0.77929303871530275, 0.76501429569633417, 
        0.28143584593185111, 0.34237157048639455, 0.30762404519154729, 0.49826292620703816, 0.16637146290975235, 
        0.48006465597555525, 0.53034617050419142, 0.46991508339210442, 0.94119705017322897, 0.73386058350383188, 
        0.82753094027300733, 0.26712676736160224, 0.26160239133996366, 0.37005694532321853, 0.19588327076777823, 
        0.57872788728381086, 0.2646941631242461, 0.58794428098961071, 0.82666690520893549, 0.54386783419821461, 
        0.58865095707048032, 0.22345673790439657, 0.10355706391445496, 0.54209349595379996, 0.28388124079837462, 
        0.031039255203639149, 0.69630894976746682, 0.56386393317600381, 0.79878735527515699, 0.4118673324827431, 
        0.52525396058884155, 0.061744514055856792, 0.81700253732117556, 0.18936916023881745, 0.79030750576848696, 
        0.7754320756879467, 0.56885578148782545, 0.41685936418613823, 0.20706687750177366, 0.92884533443180262, 
        0.135030012059278, 0.38950438923588027, 0.51619855057890818, 0.13471534770579896, 0.31099607960707476, 
        0.78048989855389683, 0.86751969677066487, 0.625078387156627, 0.46764232218472235, 0.59566052042319417, 
        0.024792449132152727, 0.53791652811773827, 0.20798650956424969, 0.4541224117526188, 0.92553073185934398, 
        0.34658982821730633, 0.1408710814962939, 0.096521530190528892, 0.2197143385804321, 0.88119340628344567, 
        0.36837712332226791, 0.54544426172586946, 0.3928659780121897, 0.74933563778164003, 0.7396970992211529, 
        0.81529532256159687, 0.98435941286133022, 0.15378122967166696, 0.31114634987582868, 0.44584504938404868, 
        0.06241858716487747, 0.78236412345390915, 0.26353274705697771, 0.74700026934139419, 0.73139339137285075, 
        0.80281181457800455, 0.71594550313614325, 0.66349740301694449, 0.60750215548608932, 0.36504044970112925, 
        0.36714200170741806, 0.5648733811543587, 0.49791630796852782, 0.78584481057195399, 0.4064464491125821, 
        0.11147588543877318, 0.63113266448469352, 0.57593831493535852, 0.8453063837746464, 0.82760082155913572, 
        0.41343043539386115, 0.11557793350975243, 0.37314202012681741, 0.66759287419569868, 0.24499451855979038, 
        0.97586946990499079, 0.65913595388623336, 0.4844154988559739, 0.20640136273371268, 0.85201046084612408, 
        0.36819213392401795, 0.19522333460767616, 0.96567040278790195, 0.19563202955503067, 0.24134389439915682, 
        0.69126038904145526, 0.2373095911110561, 0.67773205630127831, 0.51622486265483714, 0.74222224383077839, 
        0.29844651022821478, 0.37107141749507444, 0.41170273246430078, 0.061258607713771118, 0.66484677257036773, 
        0.19498038374515869, 0.3293773385895078, 0.40154995333378674, 0.094814916861351906, 0.15470869355907513, 
        0.28078736405552873, 0.11663839973859425, 0.48509073140387815, 0.15133043274861868, 0.14128842284758991, 
        0.7502544825603743, 0.95113196789925203, 0.41962946973854565, 0.36495883420154729, 0.30738026793295026, 
        0.91387184814136813, 0.87146174556976086, 0.084227332291266332, 0.51904320633583079, 0.52779588857295656, 
        0.51865476626219253, 0.48129612995722848, 0.48737668297465198, 0.78815916653446938, 0.36312086094063134, 
        0.073478842812561673, 0.57567358949013125, 0.80296252462378326, 0.64714894785567711, 0.64998052523584082, 
        0.67286018043001361, 0.33197443282222117, 0.78877875070848558, 0.033116416702540841, 0.40540461390095062, 
        0.69281924676119955, 0.90691413977249535, 0.081768162482952933, 0.35946779642093385, 0.92687396390309429, 
        0.61840770258758671, 0.49723774355239425, 0.75962596402583871, 0.050090329525217525, 0.34833274247350166, 
        0.09289887980613698, 0.20223849523360204, 0.16981322677390676, 0.63852447949041125, 0.48922651902389847, 
        0.95892845068363086, 0.84735122524081441, 0.30739977801650586, 0.8652271440971282, 0.6729752629449488, 
        0.84219451178439497, 0.5893157707775063, 0.15478839341329742, 0.11992880847949272, 0.4847725265751246, 
        0.50389278773489998, 0.54266103578899716, 0.8716258252162028, 0.53465902155314882, 0.34240437071307106, 
        0.43042502957168027, 0.28412165686770452, 0.1619034127502339, 0.21847583559165096, 0.63228020575064536, 
        0.36514427268830074, 0.62907954580130254, 0.57814702987791655, 0.70452778068769129, 0.011273235143512306, 
        0.27132277772275382, 0.90178883075627048, 0.11207168996237926, 0.87996236044104803, 0.73349979399863741, 
        0.26237746853802246, 0.28985269624230203, 0.9940065432132954, 0.55165700203374568, 0.30490214593722986, 
        0.2099812838165187, 0.2255091370337301, 0.71562568141899185, 0.94973002857959377, 0.15777549246951739, 
        0.63179942114497756, 0.45155053493890618, 0.25927462696960379, 0.68312017493641219, 0.46845938238662632, 
        0.14205174564875778, 0.64628988474939453, 0.91476804519144661, 0.4693559276141599, 0.84916179047920193, 
        0.55245919060591797, 0.25287238471471007, 0.8249993008468286, 0.28023309649343009, 0.63876622449378329, 
        0.90989778505285157, 0.44634695836619054, 0.30512470883641973, 0.32123364900580764, 0.23226279785605719, 
        0.25469300488693403, 0.95670037073938707, 0.82113997704605701, 0.60111731255762213, 0.49169496868544238, 
        0.62622562722152142, 0.23338736274894156, 0.23455846463289842, 0.15098971940562778, 0.060103409530031238, 
        0.44986048835612991, 0.50274208067605919, 0.096722459096054481, 0.54457434794690474, 0.60889916722614013, 
        0.97220848379194025, 0.94001609946122011, 0.071328381334475388, 0.82584556019821531, 0.60418288111160789, 
        0.97828981940067994, 0.86610748320449971, 0.8787967237690697, 0.96364493272660856, 0.3676772274927107, 
        0.56747296016643722, 0.068275135765456518, 0.77071890186228353, 0.79800979781228021, 0.25187525216099171, 
        0.72426201517488353, 0.56439681257575081, 0.63571784821913258, 0.91636027205409132, 0.78611668289958603, 
        0.23288697017137028, 0.87932202032213347, 0.29108853806534229, 0.17874181988097204, 0.96869883396890311, 
        0.89115837807052523, 0.3532893568670783, 0.97660452462657243, 0.53366277443318855, 0.32660147121704775, 
        0.32476372079896243, 0.51166942641944146, 0.40898549833680642, 0.36035674246100413, 0.055693660837595749, 
        0.17681297318493083, 0.20527643607279433, 0.50104011184747121, 0.16706634852391722, 0.10358502987193496, 
        0.16714123751542886, 0.0236398399169504, 0.8348145460026013, 0.43914026971897946, 0.99231222259373886, 
        0.29093139008615632, 0.94251763036681147, 0.30710247739928076, 0.77308772352061261, 0.09061027001618438, 
        0.17314542195836946, 0.56453586328539918, 0.74869082527474751, 0.38217081286636945, 0.4522088633296637, 
        0.5657441459582464, 0.38757522997379446, 0.98039549547387583, 0.68642667602101914, 0.84466399559255145, 
        0.54190844134784899, 0.91468028812990254, 0.24362772589956228, 0.40867924896673524, 0.45902609100929381, 
        0.55144585867097118, 0.46189053645326994, 0.34606233000823794, 0.04501560454627751, 0.24743282533488498, 
        0.52711013471817925, 0.64575249100336696, 0.65928677139516911, 0.53384333587052968, 0.86895654623653917, 
        0.7842221898312729, 0.63554498879788857, 0.80180943515508307, 0.020472396153106143, 0.61350933906225902, 
        0.79178314758513646, 0.8018207006363709, 0.91166885647467777, 0.11459031160952371, 0.54499013437524813, 
        0.3761198850074623, 0.64203504137132983, 0.32208864350456357, 0.57970872591289591, 0.020305039530500446, 
        0.81272864653042953, 0.25784787730463576, 0.79693222237537564, 0.21829049454729077, 0.89076676955268486, 
        0.0025658264064305492, 0.13655765623150562, 0.99512069902437816, 0.6232051619906458, 0.028696660203964663, 
        0.7888955111347864, 0.72490681046587047, 0.2199901711094292, 0.85471244999624973, 0.83498766483223386, 
        0.42714802318863154, 0.23667476724388004, 0.31342421951730914, 0.87205465337953947, 0.75045191867901484, 
        0.46117072301141415, 0.29295486485162758, 0.6663039713385901, 0.20343554356244908, 0.49964441872127185, 
        0.67044102095626101, 0.3550792052754812, 0.10452323371670991, 0.9582582648677842, 0.098617066036605472, 
        0.33915278078679445, 0.23894311035555971, 0.2787268579955946, 0.51301908700388865, 0.17219850261962066, 
        0.015309791230124548, 0.38948646626522021, 0.94546777015956818, 0.3444900790131864, 0.67560636615715119, 
        0.47814445454705612, 0.90834891275908625, 0.93615474099204321, 0.28626202763491593, 0.32530766857442894, 
        0.43699688170301876, 0.53923386536507256, 0.42225710236065028, 0.24914583398661705, 0.18121730583927187, 
        0.32087593262520375, 0.90489065377226785, 0.30243038297988223, 0.097318839288513592, 0.51819000575977858, 
        0.43009002659235307, 0.71236458401257585, 0.43698460976289333, 0.46210165944748294, 0.72447421496096731, 
        0.67344556349943496, 0.72279656851527996, 0.17641877077898727, 0.26768639909012437, 0.61351908932221133, 
        0.0075615715523118521, 0.020811322553291545, 0.75409744677008783, 0.29371637412208895, 0.10111288685650344, 
        0.92706057368685379, 0.9154044041709144, 0.94591666828126941, 0.62158631375524798, 0.33424060008014989, 
        0.17287922935783584, 0.35022483281023065, 0.27795590161067829, 0.24055685414217098, 0.62230177020057043, 
        0.96148504517171718, 0.83225298343420029, 0.94960333436978073, 0.30269906483070663, 0.53752436767129552, 
        0.19301771375082022, 0.75595896349574088, 0.46312235720543371, 0.48430765977479062, 0.33666952251474136, 
        0.69553197427410707, 0.30398468056777306, 0.7566830913845477, 0.76830611549932026, 0.81144321987170787, 
        0.86008125221497811, 0.201960061930889, 0.90519707594169851, 0.31525299367852044, 0.47193231776032829, 
        0.58936431466761041, 0.08981179648934301, 0.37244692596808182, 0.78726539186720323, 0.89041642016053513, 
        0.3926772405659682, 0.90761688733690438, 0.6090464727324032, 0.44033160503305124, 0.44350654699859349, 
        0.75559814230305444, 0.16486592897482066, 0.85863132922670005, 0.84548914037818435, 0.60074418470221413, 
        0.42810628715927357, 0.087326465682723731, 0.96955247394121957, 0.43135613425324304, 0.67793919466000041, 
        0.13347887798410341, 0.058750568786918533, 0.12142682208911215, 0.46181691293800498, 0.49803693528389226, 
        0.007775546920977261, 0.29361933254556094, 0.38450883593811858, 0.72255072784371088, 0.048854132404874129, 
        0.27508143836250798, 0.57622685679709451, 0.56059164168821929, 0.3370244285070898, 0.73299228067737543, 
        0.95688822559300601, 0.11968141174501934, 0.39836196386102185, 0.87881942774691657, 0.96180727209371808, 
        0.42675922221148133, 0.60164057390306924, 0.79222270286520757, 0.19505315363565923, 0.76803347144763956, 
        0.88746308239158278, 0.09831617483980537, 0.45677364918871333, 0.15351659571601961, 0.94828505783372075, 
        0.10693804343037128, 0.29260391470131109, 0.042405905344836459, 0.25537467844097583, 0.87426933268312457, 
        0.65906564905302978, 0.69817095990891032, 0.25493152020092102, 0.18636915043698066, 0.74510512929511719, 
        0.86700321701421057, 0.33692854872275113, 0.58644708412978996, 0.30770118412517622, 0.85589689755870024, 
        0.18654616982127248, 0.3338630664149973, 0.86745131909185935, 0.62190006447653023, 0.03090345633001923, 
        0.51349842013102442, 0.16590509251446606, 0.94357198597858805, 0.97594940246166395, 0.72829505806240546, 
        0.22223244817514298, 0.090516343921776432, 0.50662560472384155, 0.56884618604835424, 0.76567521116610493, 
        0.82742336258603633, 0.11656984724455, 0.19773620549241078, 0.70224655932490188, 0.79942696907147703, 
        0.099212838695768912, 0.58865532107207308, 0.79930198101816163, 0.52517107874645386, 0.76497261480059597, 
        0.43009956395035753, 0.34645197290451546, 0.32855980883610147, 0.75437001220094824, 0.18536315388637425, 
        0.44942609736861994, 0.43333948127828337, 0.43323123674212738, 0.79239771369899081, 0.65315778462712037, 
        0.21994697437200461, 0.0041941140716268421, 0.38961038534498771, 0.66873800219987545, 0.21925937935146478, 
        0.14947432171281605, 0.34791565815441738, 0.087150798296506071, 0.42759116646626905, 0.7688100905519899, 
        0.44171134163386294, 0.96277445054722777, 0.34955293425657152, 0.12298235152697434, 0.12435964166745994, 
        0.76592636098868438, 0.62936337903716444, 0.35028912273898194, 0.50754377153071717, 0.31457355775344098, 
        0.24861571998082344, 0.093107586165917855, 0.12384030853195349, 0.20835554275199808, 0.20373824093732829, 
        0.030734438423144717, 0.31506108419182199, 0.90183039652991837, 0.44469702865499583, 0.19984050927532482, 
        0.66374267865961989, 0.73132331654258365, 0.87906645612034318, 0.83150589537068087, 0.6722033633405835, 
        0.351560159411902, 0.39654644559070595, 0.9766380177393208, 0.079810187699250523, 0.80925461474435689, 
        0.706970287115537, 0.26228241402991004, 0.71864082285757069, 0.94844326247742661, 0.14167685494436699, 
        0.9937004291498972, 0.10708314699106514, 0.030916909680651639, 0.33863586801946588, 0.13999640321320728, 
        0.016387492792510505, 0.85792583217409457, 0.77547843106543213, 0.87073281789826318, 0.12763008853710334, 
        0.65325901323778357, 0.21093063153393299, 0.22694976341052397, 0.63728147016441827, 0.77247007812521296, 
        0.8622232549744413, 0.61843070398175559, 0.48283549716245533, 0.79573136968059832, 0.18672372077221788, 
        0.0076801636936958406, 0.22561917767111828, 0.63884448446375175, 0.063778338985851502, 0.59584980690313682, 
        0.092213068638956175, 0.14417816035965836, 0.9321054992239981, 0.24011010826694745, 0.39077990383600714, 
        0.44032922403595887, 0.30620961942477276, 0.033980841050620381, 0.1133892640242038, 0.68700956177191097, 
        0.12411209112354959, 0.38960657821691247, 0.14708892857952516, 0.67254084135262815, 0.46277413232386611, 
        0.27993791582845673, 0.091575256293202845, 0.22722020281565536, 0.551109284501071, 0.46854221137187868, 
        0.13944430985168887, 0.38750514826782689, 0.89042697433692064, 0.5684651964634817, 0.41114371626213764, 
        0.78133152095130853, 0.22013868568702999, 0.18180946081067284, 0.081239904067349267, 0.89932002652209175, 
        0.54262483195955835, 0.66837812211379899, 0.96871061460694596, 0.72742173011173139, 0.84386189393513145, 
        0.33808518016311084, 0.00042264856220475266, 0.56401240249969287, 0.14323145080817823, 0.75562169815316049, 
        0.34233246597730593, 0.19917109624562768, 0.88804768522441901, 0.88731396291328646, 0.69499872599505941, 
        0.96075786646883188, 0.12964890860609635, 0.66444295056563529, 0.23478034290816852, 0.040258830621483943, 
        0.7277364458672213, 0.80460416367711374, 0.81610440637353832, 0.78396661172648208, 0.033993516193889795, 
        0.71368795583284772, 0.11874114480600784, 0.99084075210535305, 0.5333639972429749, 0.69741845291656479, 
        0.8336425405074761, 0.97393654956731623, 0.010266203597645074, 0.4571567511649699, 0.79849897615752874, 
        0.015697661659073159, 0.17121825140375679, 0.52577419182148177, 0.8340747017789123, 0.32319238616542845, 
        0.033958029412426205, 0.65790396494910564, 0.87711572866368481, 0.1809663542851434, 0.11873063847831977, 
        0.88383786062256586, 0.31182764051428258, 0.096108581479228405, 0.67939484143933804, 0.065247004948692622, 
        0.10777413158255, 0.36013106016672136, 0.26562018936526521, 0.52141878761290927, 0.15784698491818716, 
        0.73377410934501563, 0.18644647956954885, 0.77476971794990646, 0.11230832654608536, 0.3996311975722, 
        0.19722629237101441, 0.66672712690626867, 0.19274267445278581, 0.21964812373878684, 0.27910221324211082, 
        0.71948928146111135, 0.47973961177511981, 0.82636979826833268, 0.14120532582838941, 0.48543443046028023, 
        0.40488923169776325, 0.15028317659074841, 0.76158548115742586, 0.009866152399081729, 0.41281445507512093, 
        0.51285325463809839, 0.57014016386801214, 0.85182084627814625, 0.2122331460366067, 0.99321175407138518, 
        0.55864621613395204, 0.99739886259880772, 0.62164065801052959, 0.71503189035805703, 0.46384800175324803, 
        0.65021586932843434, 0.72270205681480948, 0.19777122073573006, 0.80211581136545118, 0.87875953176847643, 
        0.55763108812184425, 0.034657249207172391, 0.24797776815049088, 0.047560073931739, 0.47025481290282034, 
        0.8632855590583417, 0.12734584503979463, 0.65229392526999797, 0.94778578369019018, 0.90935654816089895, 
        0.13880696028116546, 0.92311714620732643, 0.52482896800598589, 0.98530267779430614, 0.36007728092027746, 
        0.41190150849776686, 0.29025757321476497, 0.10461968901545982, 0.56378254043344, 0.54473242722608517, 
        0.54394335279909956, 0.25827056205318955, 0.78185764205138986, 0.71275381545652139, 0.19513290037311171, 
        0.36396176500652988, 0.21727086126997919, 0.82039045338896344, 0.8181963756017645, 0.18318881781487573, 
        0.33485718532889774, 0.70612725841983837, 0.11520421391766988, 0.80871095643588742, 0.93364118809574514, 
        0.064852965000265161, 0.43299469300876359, 0.64139631235642125, 0.32693414839787494, 0.62457410284946535, 
        0.35614563363415241, 0.19228657753958656, 0.097132599713927181, 0.90252281443546667, 0.10546426635788508, 
        0.45037719032499868, 0.51699651653177559, 0.51370964118689888, 0.69168110165985119, 0.91283272201604815, 
        0.87132275103320023, 0.63742887441648755, 0.45652369564566819, 0.64837654344894058, 0.27925804591509285, 
        0.62242005471678996, 0.99760320631991317, 0.36716097253491387, 0.39923397279411099, 0.25051036035808094, 
        0.086124350687703188, 0.12192559322466168, 0.42514546503491468, 0.0033182134933804619, 0.52309237133929365, 
        0.52481609757742809, 0.15211180299311411, 0.49532567528343274, 0.29106901196264312, 0.75356096335098766, 
        0.0087957201855333178, 0.30547737455673274, 0.81570463038979124, 0.7043010639349625, 0.6036605978671965, 
        0.34454047491485729, 0.13697763419621833, 0.38880781936910869, 0.48643087345970426, 0.98048270022968098, 
        0.064566142121777359, 0.44121606949407077, 0.89246030703985335, 0.81288653302000435, 0.2146308480517769, 
        0.95028953231025581, 0.34287359314491961, 0.67434127075410899, 0.94299927662252681, 0.88613422912950401, 
        0.099352374533093979, 0.16494060347257422, 0.37327562899260647, 0.4388655051919963, 0.30008150385078158, 
        0.453829101825304, 0.25230466832003673, 0.1631641664257315, 0.76776297444119934, 0.88469127757025556, 
        0.95953044795197306, 0.43035254627876451, 0.18994859761919258, 0.92432206202251965, 0.62844347318284344, 
        0.61392036804594685, 0.31671179070971656, 0.86579579163616094, 0.21866492676480753, 0.4535665840000731, 
        0.69610584030631162, 0.20578231451309659, 0.25055506609751488, 0.73777058591398936, 0.12991856420874637, 
        0.81697366557477835, 0.22925391602589174, 0.97804181541837654, 0.14652076652674473, 0.57945971605152358, 
        0.5292664536982612, 0.10396057307189821, 0.13871636656986408, 0.87719711541830603, 0.0074222287202303505, 
        0.28463046917059698, 0.5574685397457908, 0.67411289106376127, 0.94215046288502702, 0.69274577538634685, 
        0.011998911781096577, 0.44224877770953874, 0.49323444365640268, 0.5261392317742597, 0.91469056323824649, 
        0.24052007215554161, 0.2005242910152556, 0.65479099670466279, 0.81930652371355084, 0.48947559471282243, 
        0.7171410625184198, 0.82213410679558807, 0.39926753445152507, 0.45036471004677492, 0.54022784228500709, 
        0.15210691264270881, 0.62125621001587827, 0.1923132632318687, 0.97687853439949013, 0.82000237002622334, 
        0.35811793470447717, 0.022817877638355855, 0.87878726686351549, 0.45137584709168777, 0.11311149712398394, 
        0.058719804962513278, 0.65996486406639066, 0.71751806814360797, 0.57805784662789805, 0.364333892547543, 
        0.095585092474036459, 0.5044559400998847, 0.78027658849862847, 0.78051526281219274, 0.46023267259119827, 
        0.049274690792621811, 0.94530215395499417, 0.42081737055633894, 0.3762801097049342, 0.41708574470454529, 
        0.70071467474924276, 0.82329484801570674, 0.27855229455212749, 0.48925940411111313, 0.82266831993543388, 
        0.67651117507445147, 0.16943723852664072, 0.89621747079204428, 0.87794214161498152, 0.16928260671022111, 
        0.42978382868271625, 0.65546434689652644, 0.61802711860520443, 0.11494007559423269, 0.21317698726484724, 
        0.65314795808552994, 0.35374218687683778, 0.58479393124003898, 0.36379054769947361, 0.73237802928774798, 
        0.34292577895606136, 0.14198215092481137, 0.660465867947023, 0.70468844569652633, 0.04784056346508736, 
        0.81620483301389979, 0.10435005343166948, 0.93412021450024674, 0.89227414712888775, 0.61012645711780422, 
        0.097250290593103639, 0.23470754088008294, 0.57928740737938145, 0.27142626100041212, 0.064699869997820736, 
        0.42331341196253303, 0.74295303599169338, 0.34262077665554203, 0.30023605079555415, 0.94967229600551506, 
        0.84806527741172877, 0.1740383341655436, 0.17222700445327987, 0.38172714706746591, 0.47233191694640952, 
        0.5175325301604079, 0.99043196225950392, 0.62685711171726166, 0.56229053706952903, 0.60840403794653608, 
        0.98471794635985299, 0.6045231395285775, 0.15154658019338774, 0.45772773610274786, 0.0032928349528966727, 
        0.61798053439808953, 0.28698213303790787, 0.85351313102038517, 0.18198291054969862, 0.82037566145628249, 
        0.61246812766677938, 0.13526441693120472, 0.49831945813601841, 0.92163102670483399, 0.89264997961096526, 
        0.52734664732752234, 0.11764110101912051, 0.6323141110684507, 0.0041023699621061116, 0.035691200468101369, 
        0.98848602128810059, 0.7097033661813501, 0.9248078794239416, 0.35788201795406649, 0.19299722239127859, 
        0.8699581723561034, 0.041030908735036897, 0.53641677667376619, 0.20741950298044909, 0.41015495022557769, 
        0.11711806912110023, 0.51326140178873514, 0.99096211711051563, 0.92486280873326798, 0.070134210512313944, 
        0.92774025994342901, 0.76874824163748157, 0.89710341353958944, 0.39130210872685312, 0.16038366628542478, 
        0.061928904676562979, 0.802416185306287, 0.081100850340488639, 0.13603877828663458, 0.085333974872683127, 
        0.42482961157904331, 0.96665057411953503, 0.60184176411569634, 0.71498135015241893, 0.043051054194273686, 
        0.48632929164073291, 0.26181662235027892, 0.32352982673082864, 0.66078237309137111, 0.44243195912196254, 
        0.60355641097232215, 0.60985693343147829, 0.3463207033000486, 0.76067572075839918, 0.057796269968503067, 
        0.19206125565811782, 0.64153704715310655, 0.93784999229458377, 0.18508852439475953, 0.68897479296689634, 
        0.23065019426169897, 0.24677784393452562, 0.75436349867419428, 0.56970485456604747, 0.070091433370931489, 
        0.62795385364186784, 0.79688827436024079, 0.14152291377034532, 0.33134116075556053, 0.19329420889049254, 
        0.80398161741312424, 0.40413289871935154, 0.24508979585827029, 0.46991185086002885, 0.90820557294242721, 
        0.085626226207102363, 0.080211664718731646, 0.37778072079199609, 0.85753085615083258, 0.636039065661957, 
        0.81402214625864411, 0.073021504475513588, 0.98133497594446117, 0.83005182347157969, 0.2854834782662079, 
        0.99716188359519409, 0.93922945501016453, 0.58305955338936055, 0.0900677099958056, 0.74546793334973316, 
        0.98270462424381932, 0.45415920229705442, 0.29176706788682716, 0.56089189219204938, 0.12668671672941767, 
        0.12689043525546362, 0.90157363928702727, 0.90073619696121954, 0.81779964917389214, 0.18807406288824247, 
        0.51783188966195581, 0.6827382039442873, 0.71442552818121108, 0.92722440190270961, 0.022200427678792156, 
        0.55937979886273403, 0.54758271831742822, 0.073575609025651634, 0.5242194184222948, 0.3048656991024814, 
        0.077536481217001629, 0.27864752439871276, 0.83467208173369034, 0.67685485221933828, 0.43127080833494702, 
        0.15547775238481853, 0.13340204558176216, 0.76848556526218759, 0.21389648924617433, 0.64204969474453577, 
        0.26666716064024309, 0.24800110862468805, 0.34411680605350248, 0.60593063191647389, 0.99895571727178423, 
        0.96195632186429259, 0.74592611561028588, 0.53551970807641469, 0.9527350666835499, 0.21486522080371029, 
        0.77838135742569992, 0.18508098424320218, 0.76815361332392329, 0.096265554953175014, 0.2619938465767655, 
        0.58330803437330636, 0.11811344695135295, 0.91239230435175567, 0.75194816670391362, 0.058473508266086327, 
        0.71379265556438654, 0.91321905686609361, 0.38325029040261116, 0.95797093335007699, 0.4334642030157263, 
        0.73243900001115381, 0.37055936882537455, 0.81890911149878942, 0.22173796476648655, 0.79996559286572166, 
        0.61313249612181608, 0.50884916297650018, 0.30747856667045004, 0.16274917139457412, 0.68856033929021554, 
        0.63800116785368188, 0.2498470556925545, 0.47712084322129544, 0.61022484141390798, 0.86555468487364284, 
        0.90158063301113223, 0.0048459779556786486, 0.51029180940104402, 0.48481014135388367, 0.66225915142022629, 
        0.68699630649959698, 0.085382558897182825, 0.23620099218876178, 0.86436671009410104, 0.6701764711187026, 
        0.79084703268470036, 0.16503704933606023, 0.16549361272358287, 0.70985187793990989, 0.93026859990413846, 
        0.77084019096960654, 0.18453837728803646, 0.92723593910433588, 0.49229588547544068, 0.69034245849778864, 
        0.71218334265026417, 0.21191810664749422, 0.49103831189726299, 0.35169158998835592, 0.19641816430600278, 
        0.18725335719774572, 0.22761772287817528, 0.21960364767590956, 0.66621855466750834, 0.62381585360950331, 
        0.81074226027566909, 0.35362745447747113, 0.53919385854467006, 0.531763995082003, 0.48073093924852461, 
        0.38563682766345853, 0.64815642673276064, 0.65595509827897391, 0.41326413198321599, 0.38126151515450202, 
        0.046059596883986087, 0.34940201501114787, 0.66368502097213611, 0.62058958148768117, 0.11399652255721793, 
        0.4740719325396856, 0.36113240460890039, 0.53269370911523017, 0.087149270773776211, 0.23086621097695215, 
        0.07720287739604359, 0.98340744192390961, 0.51044080335566044, 0.96090969886074862, 0.53043471867430458, 
        0.057118177105330181, 0.89712789286618433, 0.51375094196924453, 0.20395953255102039, 0.71433741824441754, 
        0.97655427851622845, 0.85896894485123498, 0.18097013312236943, 0.99931707556372218, 0.081756645526086524, 
        0.58453875600788119, 0.56143282363513936, 0.34838747453155627, 0.80810116678483102, 0.49506732347660609, 
        0.49701791478658386, 0.086598866170217725, 0.32373536695189364, 0.66466732324105648, 0.10514435605191497, 
        0.91508623392404154, 0.78566663863600317, 0.2869927657459912, 0.092803785983774478, 0.30663649749010191, 
        0.24560558131826071, 0.59324937061679228, 0.4917498030812999, 0.74892753729036055, 0.64478801351804083, 
        0.94969899714299566, 0.17114223815193674, 0.46281493594587442, 0.56274800055192564, 0.26542816973726202, 
        0.48973571248334546, 0.7845343764622148, 0.51479322744066058, 0.74080604832491947, 0.54986423574417342, 
        0.29997216199846766, 0.42583077370764788, 0.85497601435885584, 0.89737248636640898, 0.18533351041277313, 
        0.67420219956453531, 0.67681153934606053, 0.66487813957745834, 0.0044012308189453897, 0.99865895805143867, 
        0.2891856974642375, 0.90584458090233211, 0.57267882975840179, 0.20432175179458945, 0.33266363923941045, 
        0.54079493606360485, 0.87224035500413666, 0.36637756939996535, 0.92422759144572675, 0.12405364473440383, 
        0.88067330884181416, 0.98833135085525825, 0.2206140601397506, 0.60215151393223909, 0.53403693317270373, 
        0.86493735184195719, 0.52652555066748374, 0.65289907940536329, 0.14692717720973225, 0.58516332824045181, 
        0.34191760832150497, 0.553076047748295, 0.37522731596149295, 0.16904695705711115, 0.65982828204304567, 
        0.41907501690546378, 0.19489121587375502, 0.72411459686629231, 0.22947910043279918, 0.98237557665958564, 
        0.59260224962875574, 0.65441816989330981, 0.3517232379055486, 0.5021010599203839, 0.04842874686093479, 
        0.20184954654641629, 0.19887649818282416, 0.60104553160770191, 0.78986245741765848, 0.64288438825278371, 
        0.6821726600908451, 0.29098751119912181, 0.13986107476833309, 0.14447843682536932, 0.4016494615424735, 
        0.48493951397456025, 0.5157683012415708, 0.22174247286154625, 0.14139503000032505, 0.91268931500131734, 
        0.14594273454096851, 0.69944410843734728, 0.44730936740355309, 0.24464723030490076, 0.17672336922707776, 
        0.92693693161398305, 0.82819461006549799, 0.00099778320790111508, 0.043179352996109754, 0.81966799566207293, 
        0.80933295938025029, 0.76877787701972089, 0.12202066638881948, 0.56344522583433254, 0.70306995411318307, 
        0.23963185452174196, 0.47875264121115868, 0.10868757712084953, 0.68799902294987225, 0.60759016253749598, 
        0.39229675444469514, 0.03341086453142883, 0.32437591225795259, 0.59697165904309824, 0.97497515037870586, 
        0.68190489055521097, 0.52624043354789052, 0.40447326533578076, 0.24905667460490521, 0.94740980534996799, 
        0.944611111376914, 0.28739609884280948, 0.72019775370832129, 0.19623320398382882, 0.26386159595491421, 
        0.97047327882966061, 0.11895481197455093, 0.86718383168441804, 0.44068641200149616, 0.031922966375316442, 
        0.27806218536246052, 0.77527373177920822, 0.60252649980215756, 0.20926745208957964, 0.98435060499061056, 
        0.19529676612810376, 0.31771296214586231, 0.80269624606902323, 0.99076355303577257, 0.3063265621893807, 
        0.78312714428897046, 0.92777193840928796, 0.23618247809394743, 0.13428872379599111, 0.92392931561065206, 
        0.49221093883777534, 0.32503219687016349, 0.04662663682401913, 0.48021143703120517, 0.92398753383234555, 
        0.080408075027040171, 0.93945967316046253, 0.39262892682225448, 0.039432315841756216, 0.88193888867103176, 
        0.78007246645642625, 0.33959551497195894, 0.54059760425776093, 0.77210271600962166, 0.37815457309278844, 
        0.54313535572852523, 0.30399433443400636, 0.86767975080334647, 0.44913724027816615, 0.43635952487826635, 
        0.76013251945211247, 0.1113392826015307, 0.51182824181592235, 0.2065927516306969, 0.56653187480807254, 
        0.7721387226431482, 0.95712777031995588, 0.69719335095909796, 0.89870700465891318, 0.067596332418728089, 
        0.34869345010393293, 0.26041024999949935, 0.53940167244465242, 0.15976994656257792, 0.20631264087979306, 
        0.19975819166571407, 0.10314423923216953, 0.92915307863442309, 0.027629670076553037, 0.18453967645472602, 
        0.78390792676888221, 0.08024874211537214, 0.098234492648339167, 0.10228067160883758, 0.61080863401599528, 
        0.40170014557212119, 0.74447100784319575, 0.89561569514219808, 0.56063737718930784, 0.080715979777266966, 
        0.91706719225801914, 0.58009071395887646, 0.096250781914998251, 0.33314362707826706, 0.88180937463983877, 
        0.48181060578717516, 0.013589048111005031, 0.60014917372129073, 0.14561289779742381, 0.66261944257545768, 
        0.50338918170937452, 0.26329831674239901, 0.88363559047976858, 0.96565491817856541, 0.65811181755760129, 
        0.25318194991055187, 0.077493030314438327, 0.30101527076012102, 0.45325906461302834, 0.81102287215884306, 
        0.75514786861531968, 0.97127571210802843, 0.41139474630752559, 0.6401702429915197, 0.85258853639590004, 
        0.86363912634593087, 0.16218645515698293, 0.27430244505044943, 0.30253852952988458, 0.72924402506710595, 
        0.88809494118595822, 0.63371034754770039, 0.057844932656857084, 0.38870112124488321, 0.57519020743633198, 
        0.086559892499702196, 0.4050405443069276, 0.91286911452093777, 0.11091567425305726, 0.82920433837471008, 
        0.88531970491748502, 0.25890097219895858, 0.59235813299691631, 0.38238266571628698, 0.30938427556966897, 
        0.89363333432318903, 0.040380636422954419, 0.90762605039585376, 0.85280872838360011, 0.58593610385385198, 
        0.55005562796393481, 0.82837564248754703, 0.56510129542059917, 0.40893160767024805, 0.69630855537103087, 
        0.69582104462756433, 0.56970918218629918, 0.4583369956973653, 0.54437859036773251, 0.35795356582209781, 
        0.65912360879833498, 0.64436438920195171, 0.14136570096375656, 0.14365927245769505, 0.42796345500123945, 
        0.39826711860241804, 0.59251588656315413, 0.11682188769755442, 0.10756877661755682, 0.068675183269778728, 
        0.44192162545077029, 0.35075063683318941, 0.80732905272707978, 0.81705179769674929, 0.49439388826142605, 
        0.47897698583290782, 0.379275680148218, 0.73633866904897527, 0.75299833536039884, 0.61819429351101762, 
        0.7085358010074172, 0.17008723053743102, 0.85327661983906644, 0.64465244437061231, 0.4458059465369304, 
        0.85171607603497157, 0.014722803076437696, 0.30976579689933703, 0.78024695592337379, 0.077051148938042324, 
        0.43558967088355538, 0.9519150357599282, 0.92714896894903975, 0.34751085387937519, 0.94602724778564085, 
        0.26993474938850515, 0.3817007622081654, 0.7581752308387435, 0.83153140808157922, 0.72904358863699614, 
        0.33576421634331854, 0.79495978479976848, 0.00070244430545329806, 0.23057648443011347, 0.16626593617729646, 
        0.39650785325482785, 0.94133531449127594, 0.055765375477294654, 0.93585430716002915, 0.65428464567262923, 
        0.32002758703566281, 0.69342191776106032, 0.50786086578780321, 0.16298250473626519, 0.87535375814955074, 
        0.90442541522986764, 0.89204468165952289, 0.95876393176001828, 0.36295702839242483, 0.19725705326646925, 
        0.77193993096868829, 0.59249070145169047, 0.23570117427855752, 0.8595729520260833, 0.86144205486399739, 
        0.72798522275982069, 0.57778263931401752, 0.81520308415234077, 0.52426241598538725, 0.51270417476643049, 
        0.7903555746822104, 0.71670552045751568, 0.3963995323591325, 0.46748426652656239, 0.30106124890982522, 
        0.032541131548257995, 0.34722372982223604, 0.3235523935570026, 0.76486903075648738, 0.70436517163678736, 
        0.96047067254236596, 0.92858346059911878, 0.94074580402662433, 0.75227898570636054, 0.065179946451080717, 
        0.31338800660894717, 0.4196288667342658, 0.50907327943258962, 0.11623610271794305, 0.14441480890308789, 
        0.46561188545159982, 0.50747214826612264, 0.79675265001126205, 0.45688806791557202, 0.18633973069102328, 
        0.87663168050739437, 0.489145090124395, 0.066326481217645172, 0.90950974075274038, 0.46450492711980762, 
        0.49983075882082106, 0.8553454404786951, 0.86283530288908872, 0.40119922687827025, 0.18316706831273422, 
        0.59340170928732738, 0.27176639777377476, 0.596740515660517, 0.0034021047284364947, 0.30162427189828378, 
        0.59190997111225374, 0.53325246385686542, 0.20167787686958483, 0.77507524323342714, 0.031531306391828773, 
        0.67568138174046855, 0.84158021599740929, 0.0035441410892882441, 0.44477775547122333, 0.30503992149658821, 
        0.83459760760133506, 0.019678907029470949, 0.9483829676051958, 0.60544211580537222, 0.7776300324221328, 
        0.40529141326275608, 0.20716292887654664, 0.15848863330301577, 0.37704928469909493, 0.38768634186445228, 
        0.81625633345700099, 0.1958227619259072, 0.40963839194007567, 0.38700644194724165, 0.78887773852945031, 
        0.18765520912853639, 0.43454496751048333, 0.68955675812099937, 0.99868266482777535, 0.43280916584323714, 
        0.10162174995919959, 0.74180438300736773, 0.91549314933501447, 0.80944396293265375, 0.73805182726966723, 
        0.66745789624613128, 0.60093674979812484, 0.83576018399727858, 0.72539268310824578, 0.0094646351148972574, 
        0.19008303976540541, 0.65285567693056556, 0.15036722493570687, 0.21243586514556645, 0.57110851916804362, 
        0.039864510112876284, 0.97837161927877303, 0.020169675173210999, 0.38974642237416512, 0.17715634522447021, 
        0.93044651980050674, 0.41433688806293034, 0.60288347555506938, 0.3069960188738452, 0.83336811035225566, 
        0.022212741885719645, 0.26633318429542663, 0.47895900977497363, 0.049480327883989705, 0.60904175491746226, 
        0.63382235809906273, 0.62105703662912215, 0.75261308627127366, 0.64199295798204314, 0.27579524002287181, 
        0.96776039623556964, 0.79190925527946665, 0.030750114239957771, 0.089043876692493473, 0.21288441649705381, 
        0.52087243568035868, 0.71069374375380079, 0.23426310333647549, 0.83070162276554926, 0.46774659756254366, 
        0.33602133475113893, 0.67949906974928975, 0.25655338396688543, 0.65643054006725876, 0.64696274207687954, 
        0.74096027793066654, 0.75976716483575735, 0.91184174306460686, 0.37344003696160555, 0.46754580986718408, 
        0.12565512551068303, 0.24032216618474367, 0.95931487969314877, 0.2020131547767412, 0.11962617146126475, 
        0.68638422491120887, 0.34793143186716158, 0.37187614265307034, 0.78900132780851084, 0.50595164147192051, 
        0.37864196817470153, 0.18226218391565796, 0.41602178812810608, 0.71272480988101505, 0.74264413314604361, 
        0.83368820649308018, 0.07652432326373737, 0.93922195267769726, 0.50640106541770846, 0.40008533905923893, 
        0.57097550767201799, 0.3598093710468373, 0.11458536137985753, 0.70135780518020585, 0.11713098360577123, 
        0.80354732809512597, 0.16096459197362978, 0.28374329442192736, 0.022024965241027372, 0.48807571799561833, 
        0.85833853798479298, 0.066486371579128622, 0.49198429888611028, 0.1375631054438522, 0.40506002250468565, 
        0.19534026402679716, 0.63631984266862363, 0.71250223229336407, 0.16840244907319701, 0.25284243464020473, 
        0.55661090646320432, 0.56201176457633251, 0.044464886050999963, 0.13315457302447697, 0.48260797659710608, 
        0.85910073534317721, 0.69466687279436345, 0.3597604691158145, 0.72669695620284047, 0.71680437715737333, 
        0.60033175066585009, 0.21467984352842406, 0.086089967308776272, 0.67535061593499357, 0.98869386374630897, 
        0.54297029433081834, 0.52957511251067979, 0.15786129527623616, 0.88637363805273517, 0.17149399046845559, 
        0.6516895111319374, 0.45722132113262015, 0.068582099882058323, 0.49387142449717469, 0.83411349657116962, 
        0.51638905943083224, 0.77406146445602619, 0.86240365779256667, 0.63294360894432011, 0.47414530917317754, 
        0.39342970529944332, 0.55398581141205749, 0.0038068383074874035, 0.14900634399474666, 0.22247426154860261, 
        0.71154378272347252, 0.34159421796044631, 0.047184089422321218, 0.40236917978595566, 0.40793815730238925, 
        0.23677824802223535, 0.13536999560297325, 0.77572191062245532, 0.6858985779905773, 0.79759711166437963, 
        0.87800618289311005, 0.63087031187749498, 0.99728434242793784, 0.99381860488506968, 0.82506236736249705, 
        0.47854529759136066, 0.95440265871813312, 0.1518303623882411, 0.31787006218154246, 0.082707866615408854, 
        0.85233816516614569, 0.40300110588152083, 0.32383041229743625, 0.71896021441357916, 0.31904337846913489, 
        0.72375919389982712, 0.015597910705205553, 0.71331295378337978, 0.23137276423477271, 0.39138185365753864, 
        0.25580835558109749, 0.099223691296149052, 0.052764115802231171, 0.89941270700706499, 0.83300487642133825, 
        0.44034215999957493, 0.10598538133059132, 0.17514252799062691, 0.75526827682904352, 0.061798959541961684, 
        0.12516304859655869, 0.4923698799496905, 0.06417906838450782, 0.72165470211574978, 0.066975975237025365, 
        0.0037910257545814563, 0.025138691169109162, 0.54106245050396629, 0.82889553518188563, 0.86523337278218682, 
        0.77678186960915596, 0.53366462866068032, 0.72319151996873443, 0.42347018601583786, 0.74562598608570752, 
        0.87600401383926085, 0.72062144290499908, 0.21219971029799312, 0.96112565943010786, 0.2895880519707521, 
        0.79281955890786771, 0.03553918912640297, 0.9718138496219535, 0.2457606888695909, 0.53780927377773047, 
        0.054867925553268027, 0.44998392271447774, 0.21483158228987076, 0.17673289159050132, 0.13313031993388402, 
        0.51165421298594982, 0.14207929244398509, 0.97964508895644586, 0.21764289808872017, 0.91595037639471588, 
        0.80910695290079881, 0.92881887323970469, 0.77247022265150278, 0.50611726988500694, 0.98216988749332601, 
        0.62984054289607805, 0.71900112690246676, 0.75973288056800525, 0.80924642781895062, 0.29487979750272086, 
        0.18684113557498261, 0.17036327695718501, 0.3353635283690064, 0.2960542301114093, 0.13326186532281614, 
        0.720482481801485, 0.50385901595016747, 0.5164264128255629, 0.30091864495612475, 0.71474812166893331, 
        0.83125905790970211, 0.94417832342568331, 0.93708012888176828, 0.67080910088971724, 0.67142669634172103, 
        0.90028499730881406, 0.091772417781659676, 0.88633967931095747, 0.35705386026002195, 0.61293360758666715, 
        0.20279751050561834, 0.57924093024607681, 0.59837007008863119, 0.90087046884357802, 0.43360273231810886, 
        0.30691364795992904, 0.020689309343057793, 0.94570447402160585, 0.49697981995416396, 0.13568229947761457, 
        0.40207167308560576, 0.48778043787632375, 0.77312758882309396, 0.10639827443439898, 0.63197450659241472, 
        0.80332539651072454, 0.71918690695054477, 0.94474198367764983, 0.79831286851272876, 0.79987100105235109, 
        0.048212019234725956, 0.7596009747057304, 0.6444118163554422, 0.41221518901722876, 0.66271951277409857, 
        0.93033235886454202, 0.51567742077186507, 0.91422435600232821, 0.31425897000630476, 0.27648741884941175, 
        0.31398837605934982, 0.56705720763491829, 0.57240125406108899, 0.44617991909094301, 0.44232219042508558, 
        0.16486285740896012, 0.1579884372634397, 0.81586258649126453, 0.95194534998650293, 0.017708391977269144, 
        0.44346595969624181, 0.32048567169787057, 0.42530078096984814, 0.65819352950104326, 0.20147337143714106, 
        0.39529913691579655, 0.59264712574197009, 0.33265297697901985, 0.11386782118345384, 0.13340946884860738, 
        0.6495211152108098, 0.61359537647921925, 0.2328426632935221, 0.46197878425322636, 0.11935586323181524, 
        0.88956798350861188, 0.42894327149789513, 0.35843931078801394, 0.91876116261417251, 0.550668623917973, 
        0.21545393264421087, 0.89600471866670217, 0.88071675282926831, 0.46257530951441295, 0.97336008796056595, 
        0.54441346418805292, 0.049291400084410464, 0.36179911110790597, 0.16613049013715653, 0.23970159421363757, 
        0.57345922184958242, 0.23412402173171376, 0.94066449697979237, 0.05535786983185309, 0.33325462631242564, 
        0.95168960821997839, 0.1108870772395143, 0.69260727556714485, 0.18387301119617838, 0.075963672452215247, 
        0.65261344169286328, 0.98240639973749033, 0.047895298527217633, 0.85087377926289642, 0.47804372304384968, 
        0.13486011234384176, 0.93755532588581536, 0.6489988396378441, 0.9192185579957679, 0.29269213088804547, 
        0.24085302911942241, 0.046850638239341436, 0.82377962500143753, 0.97029228846425308, 0.6866123055320752, 
        0.056864246940303342, 0.085983389690513112, 0.88061793522838339, 0.31606061125616502, 0.77927060675808524, 
        0.78073478776957628, 0.12445377156377257, 0.43953982078535225, 0.97049262490240595, 0.21751741466880725, 
        0.49250158616310746, 0.73477054630327499, 0.090009782436846208, 0.65560786294184448, 0.6740125420456482, 
        0.66669729406884581, 0.59256913141216749, 0.088835398497546736, 0.74988544308469707, 0.75799982766066543, 
        0.71270153707453199, 0.79955203792465523, 0.56426975212705921, 0.84604183013021617, 0.051569679036255867, 
        0.6700549314693256, 0.061212995916330559, 0.17246012150100953, 0.0054925276999973072, 0.74912008000068964, 
        0.73542934914412261, 0.74095073028217162, 0.20590036260377387, 0.46939654789240159, 0.72960516011386956, 
        0.88378108355204876, 0.084492229795827711, 0.44441928819903542, 0.3220016968078403, 0.19756028614737486, 
        0.053959504389654134, 0.78847142892902866, 0.55672131662662117, 0.86126479466199513, 0.44243094275791006, 
        0.24126033062786933, 0.72934509276607451, 0.16007926117860172, 0.52024505301623369, 0.19530903282735701, 
        0.041191192796015352, 0.79709500367040298, 0.86941847582762044, 0.63577366762223497, 0.41781795780364162, 
        0.48644114080901546, 0.79941956503309886, 0.8695203903206159, 0.32817765730828197, 0.99132985478240343, 
        0.65619120629877203, 0.014870754606352232, 0.26767512770187563, 0.5467205647112594, 0.069663419193794152, 
        0.42103026722086989, 0.86325820510959228, 0.83005850886614718, 0.86332106373199591, 0.18036700169899578, 
        0.75323723538123866, 0.11396681467552749, 0.79787947392859948, 0.29272199642867869, 0.75964376089494401, 
        0.13822205088741657, 0.75006545273098335, 0.30285004579221608, 0.9161579901824608, 0.71532056122321674, 
        0.84681433024716579, 0.48105317860177332, 0.08694974536199962, 0.33649863063105823, 0.39849215120625492, 
        0.14271442478151619, 0.25738668390888342, 0.75552785018010238, 0.32470347925607768, 0.75082722374960298, 
        0.064965962323545012, 0.7173279433266706, 0.86373108553482303, 0.42848481383974368, 0.38402794363402615, 
        0.342327943423125, 0.30484167084452851, 0.062097676690203496, 0.43021406352461034, 0.0096212460461631011, 
        0.024001282989543293, 0.0043233819399977058, 0.87195110664885567, 0.22143663440694916, 0.76570966854484745, 
        0.51915289205752235, 0.36981698925829232, 0.67850723011919212, 0.94028926579606931, 0.56966077439193308, 
        0.52068669283647484, 0.16820569971606547, 0.19409684485751733, 0.02008067551513415, 0.0087115618299036246, 
        0.30273668506342744, 0.4874122892458872, 0.39609306174553605, 0.8207761295957523, 0.6519309649684879, 
        0.22082598215661364, 0.55060970653767893, 0.64126180943204703, 0.78127758056442387, 0.81031667618327408, 
        0.75270200784500108, 0.99958512888950324, 0.72850937211416711, 0.92624495543225471, 0.089364884014494272, 
        0.99523236436714635, 0.022016308492092573, 0.74702435187144345, 0.85955705532692228, 0.97507747365599307, 
        0.84542685909341264, 0.12852183919955573, 0.13707589265131048, 0.32431086948701648, 0.44144685569450282, 
        0.97268886667332977, 0.25187961444835238, 0.84537044078433587, 0.092959296058544805, 0.052226260495749299, 
        0.0069253631196106724, 0.38205783756121514, 0.90167511768181563, 0.96398656174603992, 0.027126207085156562, 
        0.98517963694414434, 0.088312595869498578, 0.14791823821580263, 0.57357044240417787, 0.61350103914144971, 
        0.6655084919406673, 0.45327419920992162, 0.059491387811549945, 0.65816147499865929, 0.023748862676507443, 
        0.10392369250909272, 0.57343021196820776, 0.1756412526454596, 0.71447354930371887, 0.89555294434969901, 
        0.36661583061186986, 0.21711895682385607, 0.56068813043767829, 0.1874401153973102, 0.25981588952496049, 
        0.11061524945090673, 0.74484187592347451, 0.66160333668877258, 0.99800116376649317, 0.16367505760575307, 
        0.023028642247110875, 0.68026231369505941, 0.33367715444751722, 0.96415926189748902, 0.76064928696151624, 
        0.90425061775598548, 0.95170794186438434, 0.93817092266318847, 0.012236312493126622, 0.72544425627765952, 
        0.51424766741378014, 0.66854032743759806, 0.62718238480573651, 0.25794620815126157, 0.99552360287833341, 
        0.24055212054613118, 0.80881139298487037, 0.36070333012936939, 0.2423537385476604, 0.26390942140466822, 
        0.51346095505882006, 0.1428261432077933, 0.96345754755284529, 0.69936403781103351, 0.4331116286218788, 
        0.3649153289004361, 0.28375093803599927, 0.9738535399237005, 0.83202740848762247, 0.49938822727659904, 
        0.62269323063001147, 0.83115379224382924, 0.062507927483714143, 0.60124595829590044, 0.47764648384212705, 
        0.71667403939508212, 0.1242114841899542, 0.2992392551730616, 0.043005005604620949, 0.17480680024818795, 
        0.70838860956095728, 0.81066121448872774, 0.2803884966644441, 0.67412085711866032, 0.13182118200331461, 
        0.8632322739543894, 0.56620376736926215, 0.37137427395797751, 0.067561430034515535, 0.47936380866694073, 
        0.13853634215403821, 0.50245867422384638, 0.36139681695820602, 0.82160446443388113, 0.77044008198882108, 
        0.5891788483317042, 0.85161961433375422, 0.642320098448244, 0.035082296472930263, 0.25401063281861891, 
        0.82636469180180883, 0.32793650767761573, 0.22161435401477392, 0.48001516283230417, 0.63526971941746169, 
        0.20068958635298606, 0.83000587099459433, 0.84471600218010834, 0.89280436936688368, 0.054049905354731864, 
        0.7622934467605047, 0.23971401951799032, 0.62566506603021832, 0.51180816885662006, 0.76884531818714974, 
        0.42529358951829943, 0.24370312383571346, 0.50216077800813363, 0.49429250952145387, 0.53448718905387005, 
        0.90027622990551937, 0.37057713521781821, 0.24055624430200262, 0.49336925993745018, 0.95451129669133872, 
        0.23217297467876929, 0.47958879537255528, 0.75493907202824784, 0.61219196702129319, 0.077757840440601189, 
        0.81807662415155757, 0.12853765054639621, 0.62523422170278642, 0.063379171789465172, 0.2899333564532347, 
        0.37043299048901202, 0.4534829376853875, 0.031545716267076473, 0.58128979582495588, 0.59101483728208715, 
        0.35621215511959115, 0.89502891472395274, 0.60234472472499778, 0.35713727763862768, 0.89201789016812438, 
        0.81688371425033979, 0.81687674115635622, 0.0199060741064081, 0.22980713028011013, 0.36785759822808051, 
        0.84502503696915654, 0.70122813200741341, 0.710946186102976, 0.25902841469197613, 0.244782866162502, 
        0.99510380486274053, 0.034688266491385233, 0.98102243030063074, 0.54693107488832271, 0.20398962715957203, 
        0.39313737902327972, 0.87071848109295469, 0.49532316475672711, 0.15919913753437376, 0.62263484047604667, 
        0.24078942845877727, 0.36186435956643126, 0.97210213345322871, 0.43636181041945443, 0.17363328286020518, 
        0.80273952551417249, 0.46911825033281906, 0.27461258070696526, 0.85231624652786775, 0.75898366719604327, 
        0.96319241783276532, 0.0064011564417014721, 0.95517269111969383, 0.8923369020609786, 0.71825240479631791, 
        0.78694862348037109, 0.068359853767904211, 0.12171347626708462, 0.011087305233202338, 0.86442902302836977, 
        0.79529308474738247, 0.42286938223747805, 0.46190637049296623, 0.18838835086072558, 0.048241951887384005, 
        0.26337792385319436, 0.31122218099299914, 0.87640081160852534, 0.24085037156798883, 0.43431981677134912, 
        0.58963749613402672, 0.12173183681242805, 0.94549956067836693, 0.75236867192069234, 0.22018449254891004, 
        0.59281523942516223, 0.039833649580802488, 0.8923513956997593, 0.72437486223375691, 0.11413066789226112, 
        0.09601831467131694, 0.13372166357848214, 0.36078998242861227, 0.25724284381381612, 0.84354258162862994, 
        0.087298552699942578, 0.90431239967940136, 0.47651398405357925, 0.9239969292080783, 0.90848864871882973, 
        0.48474466892723056, 0.87489913408145892, 0.026727972315308879, 0.51819040184454446, 0.4800694425662646, 
        0.5179483226841497, 0.98668552768301665, 0.28658850015086834, 0.50571568220601937, 0.77990300577466698, 
        0.19409045922880641, 0.13880959935879345, 0.4614360346431805, 0.47764765885749028, 0.81038463683287221, 
        0.098073021417657502, 0.28014902053548307, 0.73592236526236854, 0.09435514035189696, 0.30439919789628744, 
        0.65399666311612892, 0.38384618315957586, 0.19527284198989192, 0.76846263209289645, 0.40115902011028459, 
        0.51625540970482575, 0.51966363143917538, 0.41814852561314075, 0.44512256208901646, 0.98778439394680206, 
        0.25057771770888437, 0.12613443853602457, 0.39140306231664423, 0.41827933878192125, 0.51496738222999894, 
        0.87336722004466893, 0.72038048165407331, 0.20554347912432913, 0.041161791282769089, 0.12999218618971908, 
        0.020761220906909017, 0.46037673465236839, 0.979480284124431, 0.63264860746775198, 0.72331344984804535, 
        0.93473485578708759, 0.33608396082669478, 0.49749399156591423, 0.28873291636328746, 0.66152407469667085, 
        0.90441020144486806, 0.086540470784498336, 0.49717304768836112, 0.15714461970722193, 0.83095319806674506, 
        0.79119305228879044, 0.034986765239921835, 0.11452004447252917, 0.37649009550338342, 0.88937245603666493, 
        0.50758354681742412, 0.73610504047816239, 0.28815576183323266, 0.20825807197059709, 0.77104245497927004, 
        0.9575209027561018, 0.6881614384154926, 0.42385871235553063, 0.49931423492705673, 0.023240729950068983, 
        0.45807223069002623, 0.59476780767240989, 0.56744758912698101, 0.95266902955883981, 0.18647158487310267, 
        0.39119284381306518, 0.93540576639065853, 0.79193966130330162, 0.61994904641544002, 0.9964013522416868, 
        0.42352752663053694, 0.73545357781974685, 0.76338418032400912, 0.92195404654471469, 0.19638905085955538, 
        0.20862701126265648, 0.032864820897775537, 0.98498519181343203, 0.37786672717623504, 0.89204075078951583, 
        0.74707291292843436, 0.77344709495203601, 0.29086131761496037, 0.51612633874659819, 0.12066661968123582, 
        0.22545421509976893, 0.39838181867855016, 0.59628342839031312, 0.59128325345340804, 0.56218453182077366, 
        0.29780313517079993, 0.89140529220004594, 0.41740923395844964, 0.98585688260799231, 0.67504702962022645, 
        0.77286732265590619, 0.82575085177308694, 0.036686039254889558, 0.70361190401670703, 0.61982760620841071, 
        0.62648091146031959, 0.5056002865086735, 0.89162231726546115, 0.21792291290775667, 0.55446485846902083, 
        0.022128178225443218, 0.94449982310345648, 0.98077056313375111, 0.0018798475072103749, 0.33600502607768723, 
        0.65223372011098735, 0.25127278969767852, 0.15026966561265298, 0.501970769050889, 0.69296240044682844, 
        0.28642560679528284, 0.7120983995077097, 0.70121531049049701, 0.32426051938332079, 0.44960701430931538, 
        0.98797805376576764, 0.60382393723499606, 0.0051834305228584121, 0.075685682648242159, 0.74182504924264947, 
        0.44568713733289611, 0.33617186334969085, 0.64234441539258325, 0.10099070974919822, 0.91975015737594323, 
        0.19204334397259859, 0.61714141750477491, 0.96906470972688941, 0.77184889533987544, 0.5723259937737768, 
        0.93073987737506658, 0.097752417193651242, 0.97241088644381057, 0.052578501734181327, 0.90239418145758976, 
        0.5020954425891333, 0.22715439889116595, 0.1371332471478921, 0.66177296949808051, 0.64929363926567163, 
        0.094002775397776128, 0.23887966010464434, 0.8517370856812676, 0.48168725490233344, 0.47526614684565183, 
        0.89383199802405633, 0.020199410864761402, 0.58375836979350137, 0.095652995341882319, 0.69839741509796749, 
        0.21913834107500429, 0.47675302301257338, 0.95287711470370451, 0.046598116430220493, 0.14016939658045069, 
        0.58568383482688491, 0.35319724104591077, 0.77825956327977952, 0.24957989151577631, 0.34019204173192819, 
        0.31578965561076111, 0.63423770261422896, 0.06337062847512942, 0.042244097017762439, 0.54861941797329905, 
        0.7594744584400428, 0.25090030383505013, 0.30670262889890676, 0.33076122832024968, 0.14923312676195044, 
        0.72787528211756491, 0.60287385565869456, 0.34490192739529668, 0.80366274196014409, 0.60168563676209463, 
        0.40343170034003473, 0.0069594594732498294, 0.77980846336935961, 0.0028291797371959859, 0.44600952179205011, 
        0.067916357333699162, 0.14849910757128315, 0.17439119655788238, 0.97047267438991458, 0.40552960685650952, 
        0.01349430297398202, 0.23746780314607974, 0.87013693025228989, 0.28284007867797034, 0.53149751493551478, 
        0.59299189422206022, 0.62755916632843967, 0.2131308128130549, 0.89284965827701246, 0.2498968434550628, 
        0.5499883548353488, 0.66940475825844614, 0.82443779407331896};
     
    std::vector<double> dr_1d {0.36158843638153426, 0.34044066343210289, 0.5915286598097822, 0.5344294727596266, 0.86005370209434973, 
        0.90006789602722881, 0.68358043018622872, 0.42710826719799289, 0.37425822947390053, 0.7177003527244179, 
        0.02417330096686543, 0.44265385458832918, 0.85769027123971031, 0.46420770003252021, 0.4996960414762468, 
        0.18536079378899029, 0.52101660882423384, 0.04170095016607811, 0.5618452797821929, 0.68477580779749392, 
        0.32586581810258841, 0.6324705703329101, 0.5877552002051416, 0.061789592896476231, 0.38095004133085153, 
        0.37515849502748355, 0.97370444677993517, 0.5398677517110837, 0.65916169345038522, 0.41216953090062436, 
        0.19067283729549755, 0.50574824016508413, 0.00655559234952241, 0.73026530802412148, 0.86394461174117532, 
        0.93709188042598557, 0.8024873760579998, 0.093954172910549749, 0.19406018204093201, 0.78591968919963473, 
        0.44833234096124475, 0.22736469592115105, 0.56593638503612564, 0.13324116488328719, 0.62263831741384545, 
        0.15355172048198074, 0.88805838768332857, 0.7429042978833551, 0.01520441588953747, 0.57764633002492172, 
        0.05393871820060081, 0.65760339697688153, 0.35503682387609437, 0.95229307045471123, 0.44357791947371283, 
        0.85433839037514225, 0.50244721174401574, 0.66237724106833729, 0.8446761866622623, 0.34595141325997081, 
        0.60813852514142464, 0.076500666601134082, 0.07341024184138556, 0.64150050170934025, 0.90381301238290757, 
        0.8474542211204561, 0.41734212837565021, 0.53014711232602418, 0.20220888456197628, 0.46375114618678182, 
        0.66502663319107147, 0.99074762069240863, 0.2995019655092801, 0.40790566208853885, 0.86460639593751498, 
        0.37530509613903407, 0.13670783733988423, 0.23830494142094594, 0.26915896962084385, 0.27354275781206994, 
        0.18476111010818941, 0.26260055284841899, 0.063202479219284458, 0.0068276865679237631, 0.82107156479742072, 
        0.58713795585238304, 0.32279264944523423, 0.14825124397671496, 0.026135121140708151, 0.47556202861911778, 
        0.15918692590623795, 0.75665539795086478, 0.87809955523101557, 0.1182470621395959, 0.8314839973250876, 
        0.12647523085145718, 0.078621160583139504, 0.53611580054260166, 0.53381945560212363, 0.17472348497293821, 
        0.082052266355279624, 0.72196304360737318, 0.32167218889648863, 0.29224166778309435, 0.30562740421223533, 
        0.492563795371042, 0.90505629074679406, 0.40359802035212922, 0.68318784444341052, 0.27740588654581555, 
        0.48325818576564883, 0.4118000200588634, 0.40178434574214816, 0.98754842517699704, 0.25130879744831791, 
        0.80299123858531729, 0.36331034384549521, 0.19416584642961765, 0.40440963505536165, 0.74997123742854432, 
        0.22328866268588698, 0.63537548794773202, 0.96235140848665313, 0.72397967057430934, 0.83235823580798796, 
        0.32457601932575142, 0.52774166725571359, 0.364389224135659, 0.96889681094870017, 0.096501782512833367, 
        0.49850251870885343, 0.68370346944817806, 0.66622058602594336, 0.80619521857916721, 0.78975174435556261, 
        0.49018619288144727, 0.45874387567581043, 0.43493920422057042, 0.73313646453592907, 0.10842200096746368, 
        0.017574462429140958, 0.060980581637381848, 0.71943408717097967, 0.35544982143937731, 0.61167651758950425, 
        0.062485738917163136, 0.9117916772263841, 0.86664611080135368, 0.083036054899900957, 0.43667931301072582, 
        0.038199327268969574, 0.36972756269760931, 0.58348322240650452, 0.21632197718813218, 0.34764765524647601, 
        0.7616944047236156, 0.7332107732781501, 0.79518442022532376, 0.91823380626786366, 0.69419600065041553, 
        0.6949238378090985, 0.68889482388674428, 0.82086119414627801, 0.091150877193742152, 0.33779138562057542, 
        0.6626033639011335, 0.580470273827427, 0.42542183558397517, 0.054804989916438007, 0.41717630479149692, 
        0.91611896204052945, 0.01155112459588703, 0.38989376936504283, 0.57962236760393782, 0.52722607282433609, 
        0.5313943334539093, 0.070601202782515138, 0.23877360243190204, 0.23065880033506825, 0.75475219688790096, 
        0.7524129433455764, 0.43108201740836449, 0.47146561365223416, 0.1773837940339249, 0.65796378542542855, 
        0.8702282320994783, 0.20186668367591443, 0.89557723973884906, 0.1423722044769431, 0.4953399513635115, 
        0.3595128837915702, 0.014131362715608331, 0.55669449104810131, 0.87854748531259474, 0.035389421283993361, 
        0.079992173747893247, 0.55788580059273452, 0.80810999635769765, 0.87966916708104304, 0.63901774988453797, 
        0.54295660945207502, 0.60860853796052683, 0.79023630910631693, 0.3686002924482894, 0.31369296249523249, 
        0.980762402975907, 0.93261641375797621, 0.15149345915169454, 0.02003265438749291, 0.16700894093455609, 
        0.83346050483741152, 0.32030925995556703, 0.89538955845921464, 0.11366106826399114, 0.42405049826988828, 
        0.024178749479857542, 0.2352012918119577, 0.57244485047593741, 0.2913173241227125, 0.2387154174307462, 
        0.79257431962622293, 0.2449769624745064, 0.47427838498269415, 0.51742877469399784, 0.24584810074902808, 
        0.045855524315099983, 0.17352450701164779, 0.56441592894826864, 0.71710723424165845, 0.72201650659926786, 
        0.43212197687033815, 0.23078569789024805, 0.55897900368794717, 0.90969462283884983, 0.83920649245252132, 
        0.23036859289295752, 0.67471396718745114, 0.59350283955143457, 0.77236564073173497, 0.68235140264916128, 
        0.28834376749514568, 0.69551712879145167, 0.16556239397797756, 0.17235533048010288, 0.85167597789967164, 
        0.15015653541507756, 0.98004521007958623, 0.61875457163573788, 0.21761671918231618, 0.040172639247550768, 
        0.46311968341786591, 0.48380723200820075, 0.037980928141983172, 0.54531703204574411, 0.90279540351621579, 
        0.66151637091548121, 0.48310738449874013, 0.6041798848876776, 0.21138612836604143, 0.64740729075984071, 
        0.62122992534320831, 0.60447403205986694, 0.41622675281349331, 0.71875570860605609, 0.56216903600276469, 
        0.59240629790531241, 0.98668649986183721, 0.81275094356029975, 0.30123708781545355, 0.56964694776527258, 
        0.51225379571536056, 0.3206243898505976, 0.6042753189891672, 0.013666753271772913, 0.90151596549178015, 
        0.21078556339309773, 0.16892968063621283, 0.2130744919045815, 0.42928573417582472, 0.1969272924213703, 
        0.71738204387595061, 0.84096996432277082, 0.50167777898625077, 0.42881658961021429, 0.59363157577610637, 
        0.71446833974091084, 0.0091003812165424414, 0.22137590296380227, 0.40759347070306018, 0.23332013796030182, 
        0.45736687317598745, 0.77456936005577459, 0.88830335746235356, 0.72356731486565118, 0.72612995745076025, 
        0.15650703490112994, 0.17737165092422824, 0.91828268113729505, 0.50049058594608287, 0.96199422851937122, 
        0.53296838874850749, 0.80754647758341802, 0.23083560156220062, 0.0005451018850557432, 0.14051173774481218, 
        0.95326319348442645, 0.014289637969761237, 0.19823426721600934, 0.98998052364059541, 0.4780039833095151, 
        0.33064942148032017, 0.92851323233903504, 0.34230207275681224, 0.40165028713117201, 0.062253053630644262, 
        0.99707364965515399, 0.76757761880453512, 0.19123228538446346, 0.39758882697305475, 0.90116318281690089, 
        0.07870361358756961, 0.42470497791653994, 0.83074672284723716, 0.16436267033702934, 0.69386260511381481, 
        0.85381138766101894, 0.16112955366169102, 0.42597049520578145, 0.27615983119767051, 0.44964938764031315, 
        0.71662335625167839, 0.30416862666588829, 0.92349121065173345, 0.90713810313343135, 0.58792524633928323, 
        0.53648971814155089, 0.23106375968107362, 0.83784497490463461, 0.20507452241898894, 0.40427645660920364, 
        0.48734991040817532, 0.22979544643187566, 0.49699242165011981, 0.92648052514830193, 0.055754251505268693, 
        0.29014465687078705, 0.44206011919590749, 0.035721819744508609, 0.50866672280032343, 0.40498365115090662, 
        0.30094818648372734, 0.7827866110418551, 0.72221293341874215, 0.58016984647543035, 0.20181217679145291, 
        0.77576604680888672, 0.48694397619233998, 0.93360320873485625, 0.23831470595866189, 0.66930846693811219, 
        0.65239805765343784, 0.31138620241722426, 0.092905168117457748, 0.49734143314213286, 0.91968723088951276, 
        0.53324917485258738, 0.27777412750029096, 0.26691007002473244, 0.97219619244406497, 0.5856867077100838, 
        0.51416788311153994, 0.77265621926089323, 0.055540281006085879, 0.6821729478545504, 0.62184157129741213, 
        0.046984494566103985, 0.76742493550763125, 0.75144122556927972, 0.27037275892600943, 0.80543730181673867, 
        0.81687891247926236, 0.92996795697940327, 0.92050061177440368, 0.97713578245252664, 0.37236319692067466, 
        0.2466216145239355, 0.0086485055977207903, 0.52699135857418855, 0.90225040644097598, 0.4518553131877201, 
        0.40265611719603389, 0.082218348374993644, 0.16458961791049953, 0.32182039743032576, 0.65874917266684552, 
        0.20161319219888085, 0.83955438908876157, 0.54790912537349512, 0.27798660769782679, 0.35087597448833407, 
        0.83283577266401432, 0.025331067210928593, 0.66572964956276959, 0.80922699707429113, 0.44719182453949013, 
        0.23400835139899212, 0.40396588802727074, 0.38342327142963462, 0.76091385365842301, 0.84909661523903779, 
        0.83749444925900685, 0.034653970854451099, 0.74346988853844742, 0.49417768151798747, 0.76792306950523037, 
        0.60744617357318065, 0.75729269052606751, 0.087389986921644569, 0.37075850355715678, 0.2359181698639663, 
        0.58307852570368879, 0.67835895908052013, 0.36047277185353832, 0.16615598539735554, 0.28577961501217097, 
        0.57183709691975393, 0.23454161156948472, 0.35038229169143986, 0.2020469692175344, 0.30738099199236424, 
        0.2711972634768951, 0.65703784633691042, 0.59072313583785352, 0.014666175399778547, 0.29075429090379989, 
        0.55012239994805245, 0.5263903858711223, 0.68966689802839798, 0.63305394163827233, 0.6924567281268148, 
        0.25962629827248263, 0.23354093688084876, 0.72266949612020315, 0.072181614234081204, 0.14109643981106723, 
        0.39061361740356593, 0.92183454635098139, 0.092625958504874095, 0.27315314686345493, 0.14178497846199112, 
        0.85422429409458167, 0.72753113963507698, 0.66032149826721254, 0.642601941665772, 0.51299126475360923, 
        0.50355871012260356, 0.60173124976595682, 0.96588102300796397, 0.4192774049581407, 0.64912770896128191, 
        0.71659478472064642, 0.81882296608685823, 0.92332616038653259, 0.14119940817192078, 0.12575846070089214, 
        0.64667759061159469, 0.027358120751201032, 0.096943831859612084, 0.66944511475965296, 0.012213922935105481, 
        0.070234534952495986, 0.6023856803356511, 0.24633782466817356, 0.94736921269128449, 0.36241807749861521, 
        0.065999068821099405, 0.34645252502405666, 0.20438088122129994, 0.27613454508747415, 0.73015872682865024, 
        0.82762700709471226, 0.28111767458312387, 0.31754828459175832, 0.35066066261972995, 0.48911542841231737, 
        0.68435514637391992, 0.033314335418147323, 0.69626259709193628, 0.23867145142407864, 0.64203880288216597, 
        0.65727124131439951, 0.34004943716433367, 0.93294397381851635, 0.61258519302843384, 0.55562414003558036, 
        0.99954624538951231, 0.87252260507129464, 0.14903376382703759, 0.1913236958206268, 0.19941423223278942, 
        0.77615471218821175, 0.15137845352473289, 0.22709176533224773, 0.97612257151108528, 0.56019844232100868, 
        0.29183784417939473, 0.46751624162374639, 0.41700393206965147, 0.62322091697834581, 0.95428075274876534, 
        0.10119163610303272, 0.51271961321059423, 0.73745308947976262, 0.057221534885241399, 0.82827033628383218, 
        0.94786035411781455, 0.17085191732465921, 0.76204923465264107, 0.85306517786268921, 0.18712154496231581, 
        0.76723060984953628, 0.15104813782469106, 0.21451477176807932, 0.85847283504539673, 0.84954486236639615, 
        0.28415883228340166, 0.79100066420485016, 0.40045014737891282, 0.20839134234333523, 0.83019049844585524, 
        0.57104226489785548, 0.50240220492320575, 0.54669372507718372, 0.40600878203624768, 0.50830517043565493, 
        0.094573267727308563, 0.10696683986850464, 0.26114558340608429, 0.97091390102949804, 0.26855640596523567, 
        0.20091122506164094, 0.81837391416015071, 0.14167326776446254, 0.32916015750587357, 0.91427831662179826, 
        0.12015375543969453, 0.20308518703135192, 0.44052472314104674, 0.35755668334292645, 0.57448197979359383, 
        0.83675324990866984, 0.45104118422880224, 0.73503677293292879, 0.11871413479491166, 0.070744457945359773, 
        0.13939819843919277, 0.54797184004082755, 0.30784084256158328, 0.31545854333305656, 0.67795819083052455, 
        0.13524634611265607, 0.01017216871104365, 0.24933543814008585, 0.039255543828265349, 0.31515719446578339, 
        0.554292880280147, 0.23211227522038613, 0.42311285398246334, 0.038133185216097765, 0.45836044424314348, 
        0.87511847792681308, 0.034508805101846507, 0.80613749725333506, 0.56361471573817412, 0.74643868034580985, 
        0.83461413613829305, 0.069193166858613431, 0.95614035438404255, 0.61656110007109688, 0.64158099841753669, 
        0.66921550629209881, 0.63679306567578897, 0.38287312662241479, 0.57247296154408867, 0.403790280453586, 
        0.53667029411949208, 0.0023001892920386791, 0.8189296130145578, 0.88429419382412333, 0.12649614351420957, 
        0.14450949114233125, 0.13013431955435717, 0.64763313445271731, 0.74780234503483678, 0.3997664449324807, 
        0.99575595033357889, 0.90221471578054513, 0.53259872666399333, 0.50260832437629954, 0.72227032571813066, 
        0.30136091471474891, 0.69731887900504419, 0.0065586434841513608, 0.6173045058243809, 0.26573758688229332, 
        0.37680299538569639, 0.27914043571051361, 0.45864345248713523, 0.71969141617073862, 0.25391075882914693, 
        0.63881725418122337, 0.14661309266515699, 0.67286768723189305, 0.81210347115102044, 0.84531442973335857, 
        0.3229306997166479, 0.16123509479426557, 0.049529685015446301, 0.61064097245106463, 0.06155628194961138, 
        0.54537933671542471, 0.41896995039296225, 0.70273489359030017, 0.316232362446768, 0.26796509578018535, 
        0.54138749706845668, 0.6355439694545959, 0.66729498452709946, 0.70078581676251761, 0.59450520451480493, 
        0.90991780362754504, 0.96818332238103011, 0.9150287003423665, 0.94861542925954412, 0.94222132619141763, 
        0.40480864648305359, 0.050146158438850286, 0.72467842295864915, 0.79280980692982306, 0.6219788891776783, 
        0.32143911489607069, 0.88246187574697066, 0.95141425134126378, 0.78412949396451292, 0.64220154395547424, 
        0.4931028542640723, 0.9010629549361, 0.85743023607293223, 0.021748960857344013, 0.69978764247927305, 
        0.99408288775599951, 0.99121477483564924, 0.085214928118157784, 0.72269635283312939, 0.81827804894012468, 
        0.69070112455522459, 0.75774576927704707, 0.49236447531548788, 0.76502123960218582, 0.018044513672988538, 
        0.6623364720812277, 0.66222300639161302, 0.85602171347747746, 0.66103114231195526, 0.7674746400904191, 
        0.22427421074177056, 0.23486107954879665, 0.45709415206596238, 0.73576628567030844, 0.48300491442178073, 
        0.10425453974470744, 0.41927771566143823, 0.88866256580403413, 0.65176425305673225, 0.51080743207134915, 
        0.2818582487905259, 0.61722501979808864, 0.70674201290626137, 0.20376470615847508, 0.76901238250458182, 
        0.83943809645401135, 0.27906475984088575, 0.65781109928040205, 0.57078058744671711, 0.58208064781334867, 
        0.30937662835275481, 0.94770710387640245, 0.57155257179555052, 0.84512609054597743, 0.015373738048586771, 
        0.66802331023778172, 0.73729262931330086, 0.51956694236503931, 0.85147236134550019, 0.66541514880445241, 
        0.48119767829492854, 0.57395646527582023, 0.044629563367653002, 0.20528579834991456, 0.04178043489281924, 
        0.98780727915709399, 0.2089569155588249, 0.88981729265455289, 0.019115569723019421, 0.12410688045546259, 
        0.54531114477771769, 0.48813313459667262, 0.11419186497800182, 0.89399989497991283, 0.82435631168516754, 
        0.59597191535523675, 0.31116474592238808, 0.93532889310716327, 0.1141336593623834, 0.43960310978438111, 
        0.77918390208608312, 0.56670524743944273, 0.62204043123291841, 0.72267623426146521, 0.76379840146618361, 
        0.84711212338061115, 0.97448887432771669, 0.24568149391401062, 0.66437662204480263, 0.080445591489719526, 
        0.79667513014434088, 0.92146450822050929, 0.86645822354709057, 0.94318376572843454, 0.2781439599048221, 
        0.28841080753349257, 0.86410539581263546, 0.58417600716873563, 0.92079230644852861, 0.061281373552294127, 
        0.69980655467473318, 0.98261384236838056, 0.14587530448060959, 0.91074420594204231, 0.44849406704016781, 
        0.16178342472754759, 0.080515814841560784, 0.41088238312816183, 0.98994223339775789, 0.56503155767248336, 
        0.85371914677101768, 0.98340868521472058, 0.64925668533278835, 0.53467173605575957, 0.99427379069737154, 
        0.54469385588356212, 0.83908441238892961, 0.28399861104356061, 0.78935963404816945, 0.46367817257663413, 
        0.52768752907916827, 0.61102031115108102, 0.79149400479616538, 0.060482447414313079, 0.56187621821765266, 
        0.84541581308544722, 0.3593547359677034, 0.7150876175615628, 0.48030666468975602, 0.75612588541323511, 
        0.62346476634153891, 0.51838767469963676, 0.93662139977558612, 0.28467831416592904, 0.13374182552129121, 
        0.2471809090998669, 0.5749032411202637, 0.58431367906517839, 0.70911303975077766, 0.021715099218688305, 
        0.97430850625247456, 0.62677595037560963, 0.029539086100459855, 0.67645215340059073, 0.71788606280107969, 
        0.4644339748234998, 0.38213383196685347, 0.93101476043073483, 0.022285472783080129, 0.94278082031301769, 
        0.77509724129860391, 0.48642754840163449, 0.27708274267884159, 0.1883658259398866, 0.0027549745274999538, 
        0.13570531857557633, 0.146991497327291, 0.84752111529234009, 0.41882669308414244, 0.12267009464558898, 
        0.26666693664344021, 0.86155171809624664, 0.95553831540446765, 0.81280690798809263, 0.32346988883588668, 
        0.20554579897087755, 0.052364273558146657, 0.28748731500744662, 0.048843116733559011, 0.34204403656669613, 
        0.91929039398759027, 0.82183146898717929, 0.59548491969410278, 0.18155149228177292, 0.82439350850415294, 
        0.79774136222476222, 0.41341071406573726, 0.89682358923213257, 0.0082563611956734118, 0.53675150379735492, 
        0.43402924757178263, 0.54927950841128559, 0.33742071474882152, 0.093497297183514405, 0.47476927636260058, 
        0.019771309406294568, 0.2349716162932054, 0.81096643369445398, 0.93051549793945321, 0.25653470510491005, 
        0.73593765989018678, 0.23660355665730748, 0.23395979187059379, 0.98238749513759771, 0.42634543048157547, 
        0.41238285917872819, 0.070411795010890232, 0.6135779228157987, 0.37887020184834141, 0.89908999366265419, 
        0.63113176860129161, 0.90868335304523429, 0.77008300408250241, 0.67958877159622721, 0.76368950969263594, 
        0.17917029260432615, 0.75954282001870044, 0.14418504452144632, 0.89878026925830401, 0.48722968040624037, 
        0.97973087030263284, 0.30038391791125751, 0.58295525221313982, 0.33165432008960605, 0.94668941344174029, 
        0.24539964307703599, 0.87292421488526273, 0.25298093351738937, 0.6674969291269035, 0.53744384378844323, 
        0.89558273341229566, 0.80351263342702728, 0.58658339980265128, 0.25397135031070928, 0.66410941075308361, 
        0.50766911317954655, 0.24372619452357269, 0.21181363830169864, 0.28144407709957608, 0.82229507821543502, 
        0.31664599272070548, 0.097341188134732803, 0.07890507782242806, 0.29090474687118029, 0.02704168253502881, 
        0.62885293712655588, 0.80563358802600682, 0.072572878708977662, 0.17963497920948623, 0.62565585892032272, 
        0.22266029327590298, 0.89611603726330791, 0.1514538969074426, 0.68468940365754594, 0.00054805275731095726, 
        0.12195022468157224, 0.70188583350862377, 0.94344109820841426, 0.51333957596869539, 0.59221247188568049, 
        0.41288881116511922, 0.76958650796393813, 0.24981673802402748, 0.65778749011354054, 0.68355316390952248, 
        0.33047690311559741, 0.92027995177250843, 0.88623578442140616, 0.77460101523854497, 0.29657489845117135, 
        0.038391896641120971, 0.86695861610175529, 0.79554245106400479, 0.0055397421260174884, 0.54260685771726136, 
        0.87927554115206163, 0.47508507019194246, 0.30213873740495778, 0.73279217910045191, 0.27709147364957221, 
        0.23011381006425213, 0.53139625134940061, 0.3058313643122712, 0.23702186943496706, 0.3999628345191244, 
        0.31972088898646667, 0.83785307342483351, 0.087466373393521701, 0.1150060948037539, 0.091628360768955419, 
        0.89056401837327903, 0.56176209489965045, 0.76480559186346775, 0.96024915435058111, 0.31647023352537951, 
        0.53205489138883011, 0.31439259696121424, 0.2376128083088096, 0.093957675816208619, 0.97967479474417152, 
        0.19816233494097424, 0.20313669671694434, 0.29883490192273232, 0.31455873483555852, 0.013401488850230514, 
        0.40354789050706064, 0.77560520285898016, 0.88988430768763083, 0.80327609890853924, 0.29956557094552361, 
        0.52814178247265464, 0.97591788319568096, 0.74935041716881345, 0.27104564959479971, 0.35246018056450557, 
        0.24848365756501689, 0.72691698293421814, 0.41604649802794103, 0.73305001883140686, 0.34530102885485858, 
        0.59483003825614689, 0.7370303964021836, 0.50231464379921875, 0.16124070001406188, 0.99953786960116653, 
        0.70107253526459323, 0.45233139691146329, 0.74485003792575966, 0.20250228782362423, 0.35762309433747008, 
        0.43141385108737551, 0.12936838843158549, 0.80751810236396171, 0.85021128997401885, 0.010584904411662377, 
        0.25516361854005987, 0.43852785364482094, 0.95217423728068962, 0.14986500399379965, 0.90693107187192346, 
        0.1549372234630344, 0.064531477416013283, 0.95474402648818657, 0.86985191507070736, 0.84791267145028759, 
        0.06828572788854026, 0.26640706200524922, 0.27210764091851125, 0.69725266834709765, 0.70078280273466853, 
        0.29839581672571258, 0.32806780866287744, 0.56805569253743204, 0.026521576747964959, 0.070403974252468249, 
        0.73749489680641922, 0.77278294177659013, 0.34911485125753394, 0.67031927604607211, 0.31297586975086444, 
        0.9678339529929505, 0.95957974378256949, 0.49969397280404237, 0.24914140708438071, 0.45648484630707542, 
        0.0036593694653832554, 0.6996565899185474, 0.61816359275622612, 0.75171242555566375, 0.99441915734177622, 
        0.69409442003631372, 0.068321649286836372, 0.021267210272956527, 0.22956824024046085, 0.37880709256717449, 
        0.99288910509357931, 0.63048489062734236, 0.27683691633271645, 0.10332123859616105, 0.51182845947618505, 
        0.60677017199381034, 0.64794177445614443, 0.70438129879727263, 0.065496489132994773, 0.94139794635224039, 
        0.68248803004457725, 0.8429035872909747, 0.5248022246558266, 0.6351417625374054, 0.18834257038748148, 
        0.067375558724751539, 0.90307243616316968, 0.93001061168711674, 0.53056973355164261, 0.14906718034799438, 
        0.83185005090044295, 0.0091351824192804632, 0.66797524583941414, 0.34800466005665265, 0.40712778798054883, 
        0.11659698306181232, 0.8650463938978199, 0.86204430757828332, 0.66643083537327374, 0.89487652525744354, 
        0.6221768352421635, 0.4209113859609761, 0.94049112335503615, 0.99685357868374225, 0.97490961349808347, 
        0.69982684396221884, 0.916957608906243, 0.060918098468840665, 0.85182696522016776, 0.37635819216854838, 
        0.79034155999121292, 0.6695366231654436, 0.99530248864050241, 0.28041969232981967, 0.60636515017487036, 
        0.50973848963160662, 0.8717555364682148, 0.47370336129244817, 0.79455927747421651, 0.032561950332177902, 
        0.1622310746975606, 0.23742198394532554, 0.77353023002481391, 0.15888497032521531, 0.43230426196569893, 
        0.90363771197859433, 0.56166803489123618, 0.52164797182604006, 0.94148302592554001, 0.40462183671727314, 
        0.98472889045520806, 0.22184136055870862, 0.18382092509765569, 0.50210665037328828, 0.30491944829180184, 
        0.35944582326007679, 0.7926557510083303, 0.071129789025296875, 0.67025951504593295, 0.76687740470744803, 
        0.33291402409618565, 0.69548498382871715, 0.52532212744225326, 0.61402754923970226, 0.26590488731050255, 
        0.42085461572741556, 0.377327078235562, 0.35810427772573861, 0.063296989001548454, 0.74638791521967862, 
        0.89092065218040584, 0.0008015510303622797, 0.13447357094327783, 0.80856498326642412, 0.26036670285335517, 
        0.96607244415391191, 0.17040067904509382, 0.68127303519049787, 0.06237243884067456, 0.090444659387681536, 
        0.64179175811684042, 0.26892339354088124, 0.92591755073676074, 0.068027502903462045, 0.040771495423203685, 
        0.58733216605187222, 0.81457256227141039, 0.76159908572268442, 0.99225292897302553, 0.023058184398851767, 
        0.35692726864159208, 0.13149529368789503, 0.043082706509605995, 0.3589743117039228, 0.20315959691646968, 
        0.82630548474500376, 0.3650356334342193, 0.89346694223075107, 0.80182198272466976, 0.022057848316708872, 
        0.77974331963437926, 0.090523754661113731, 0.37757179794566875, 0.70516602332365896, 0.55512207521603907, 
        0.20189753378643993, 0.79659989635735395, 0.38591225214886049, 0.87789770752856211, 0.56105771110019154, 
        0.83433412593518774, 0.90079067026381132, 0.96725867828662948, 0.77066313363431727, 0.97518010806100652, 
        0.56754534144753688, 0.97714494237183791, 0.28489853116017994, 0.033981879063215592, 0.50891571122540369, 
        0.61250531989325196, 0.81825884370395841, 0.26311666923195465, 0.98441417872505133, 0.20540263843162498, 
        0.042290511897877137, 0.38376462234042408, 0.48888889525525747, 0.67869855513754973, 0.47513615954024147, 
        0.028476122644146118, 0.10645174481479991, 0.31757770886452441, 0.6782836102797718, 0.96498500301195667, 
        0.2529290021466637, 0.63744965878934301, 0.75396596628509327, 0.15993741927993499, 0.34292839314420509, 
        0.46362725466428056, 0.10047820761701431, 0.63896620765419643, 0.35698373882976386, 0.88862301948636024, 
        0.93188563103756916, 0.42696309046397452, 0.84521990593414897, 0.80114458051105575, 0.69321189413842221, 
        0.20860305483574426, 0.66156946937063155, 0.13909537333508037, 0.16756372623438009, 0.45752696651179425, 
        0.18705264217041728, 0.90361491410510086, 0.82397004182966471, 0.90282863302151961, 0.30799759442227082, 
        0.41951158243288722, 0.77340152567914999, 0.57993840658594409, 0.73824690447976415, 0.041031774650466035, 
        0.8109254994616073, 0.19493952880596721, 0.56847749565497385, 0.84252103350673058, 0.86611972353027089, 
        0.20574312194459443, 0.24501565356994126, 0.32986288240001538, 0.58438052888832615, 0.33301574024600167, 
        0.38531772866760061, 0.5923686456689341, 0.91742719177978227, 0.42366491281032825, 0.66618681575124938, 
        0.11444565375558202, 0.26598713672469465, 0.85993434318441708, 0.058661938795776081, 0.25294884789606553, 
        0.36163784545592992, 0.846395046316192, 0.69433171549508388, 0.18855841257094408, 0.37504800231444335, 
        0.38779791673998676, 0.78137574679709565, 0.018658067426799763, 0.61164701883372175, 0.34712218137724538, 
        0.099758290251447779, 0.22243133021885786, 0.79365801452876172, 0.35224016201066188, 0.65679389666158916, 
        0.77982230453365653, 0.44154533088719794, 0.53527233139493458, 0.56788680786718237, 0.93187574365645065, 
        0.1268962045153883, 0.87372742081726673, 0.47582207438243618, 0.13949099968457812, 0.28089391018808563, 
        0.94632274207548828, 0.00083838451544271386, 0.6540298151600179, 0.48203515895740168, 0.90822998672545818, 
        0.5070567918040787, 0.32146428925809611, 0.34118061776607034, 0.31899249513444938, 0.97399175146554429, 
        0.43613571335465129, 0.21776237573047164, 0.93298948067499232, 0.18796919499866283, 0.43261527360235763, 
        0.84267282956370693, 0.96803148650897475, 0.96684174564679259, 0.79261218597212224, 0.73140641922579763, 
        0.60192194443483671, 0.10995770297509799, 0.16225570496687536, 0.74575498681595742, 0.30924148616780234, 
        0.72792982412024876, 0.45080313026203545, 0.68032787477877488, 0.85848951146130448, 0.2424158498760427, 
        0.46366113233606066, 0.69415764243050448, 0.26199915984522004, 0.36725036752156148, 0.91822370505980322, 
        0.0026520215159757665, 0.47721703345084054, 0.97448866221545627, 0.21070637703292894, 0.152902642912496, 
        0.6147576418066707, 0.30993612303716889, 0.75645679027098733, 0.8047458508843226, 0.69553428665484396, 
        0.61483982638694878, 0.5819512972829568, 0.87858979155206907, 0.2203456150930887, 0.40006829147712675, 
        0.46835981315182829, 0.79158149248537812, 0.58515097510830572, 0.56545788550376552, 0.064795386437494384, 
        0.49329531635749735, 0.85809082997669206, 0.25160666234207985, 0.95063705974200219, 0.87591501498094226, 
        0.74077603544336945, 0.09877185542512823, 0.344671687845157, 0.71222166695623135, 0.0031091514908172524, 
        0.90243099192774401, 0.37233528025441598, 0.28326233544103729, 0.57277328722915222, 0.42169901364760287, 
        0.0042639786259024426, 0.6368692319539031, 0.19025658505845033, 0.072849217836399971, 0.33825371876011867, 
        0.1766195474971739, 0.58801248137445428, 0.31358355898129364, 0.074786528065364566, 0.2643526644006704, 
        0.35914123920851138, 0.13555806374137158, 0.30355369226108553, 0.017772503384438032, 0.20308353017638847, 
        0.045031582004028126, 0.86682475409678883, 0.17794302563278763, 0.93818391690895209, 0.56144206178072253, 
        0.45803577736598933, 0.53130123201733293, 0.51316160772564134, 0.68654063409033839, 0.54031369990366729, 
        0.95732185746313814, 0.77728111489228491, 0.2078461175810864, 0.015878698320937135, 0.48381076343033103, 
        0.92606755738264757, 0.94876264606740657, 0.45285195804405931, 0.70407033922847972, 0.70421088525680275, 
        0.40964819821337639, 0.23801295519496368, 0.84717711166819365, 0.17831931393473477, 0.71401889182171097, 
        0.59784024106122735, 0.86049632627608807, 0.99056104424638036, 0.30008116284542696, 0.35706492210498042, 
        0.49275431642142231, 0.68636182574177385, 0.41208245157231649, 0.94627923169478834, 0.81338567420640584, 
        0.59576994458499399, 0.42280525244176625, 0.56681356768451252, 0.24784528804139172, 0.65083092014738608, 
        0.92995482982217514, 0.18905005541346154, 0.50066248943298741, 0.038206316715732713, 0.76167825491947538, 
        0.43863036480946827, 0.19828501861718673, 0.94754772270477994, 0.6896026674627489, 0.66782174143271189, 
        0.6102132655002086, 0.65957607020815034, 0.32385028194107579, 0.34223309906195998, 0.89526736187452238, 
        0.46861813277031139, 0.0010364443590868966, 0.88659987492075043, 0.42045516726487087, 0.24687914601423167, 
        0.77248927669815104, 0.92970141453770938, 0.13497686557516042, 0.83087359801373428, 0.43335282518596152, 
        0.013574759069765019, 0.3438245175025032, 0.5070482821882949, 0.67201164920670764, 0.4925672954461966, 
        0.068850033580639236, 0.12966953407168469, 0.68459170860884133, 0.20096240734367887, 0.87490152006480582, 
        0.78448284244979538, 0.79996286945774964, 0.10092998238093798, 0.14528662635746548, 0.6952383392744641, 
        0.50490799497307126, 0.10526241259272129, 0.065566875393897384, 0.29069773140285471, 0.5462296245447571, 
        0.76336228490081437, 0.46818418695225272, 0.18713558727983015, 0.20835749120413616, 0.28221010332017449, 
        0.74506632796924999, 0.0076160061081371744, 0.37906128898599567, 0.15714892776953948, 0.8872178740776826, 
        0.14612078144728025, 0.93374311765982165, 0.85886775757193146, 0.84996450266583556, 0.28338616089137347, 
        0.48002176735505553, 0.57371929120710208, 0.023164486760669067, 0.125054306779812, 0.36958752828705532, 
        0.81520680930461586, 0.74515793407889008, 0.88587576043940874, 0.80681152679616908, 0.69176501044285454, 
        0.81879122971204721, 0.97731771888047758, 0.047365407490908007, 0.3006905473178878, 0.22970929858813349, 
        0.29860409891477957, 0.52570695950468949, 0.15137196276636344, 0.26383829742446685, 0.44359199222580181, 
        0.67967284559962216, 0.14633008752185184, 0.26324475669932612, 0.6669339630415354, 0.45962907277087961, 
        0.19839852765544586, 0.10850918482730765, 0.11226928633476074, 0.81923175594277797, 0.48876261997414727, 
        0.93476910139419411, 0.14051491313395736, 0.92547544105195745, 0.9515955060722201, 0.044679816026178543, 
        0.81925995719562894, 0.23350415731997076, 0.7689044738095645, 0.4899647813376935, 0.81809968388560295, 
        0.7891208895025108, 0.20296567801884358, 0.2500400036720154, 0.13519457237168253, 0.78902377088999476, 
        0.57166839226421229, 0.99228168227855873, 0.76116292206342062, 0.52975672991897582, 0.51027075530584409, 
        0.28183407027947815, 0.39095103758306426, 0.65124162377762662, 0.76737742189452152, 0.89074591597091723, 
        0.21840932319075912, 0.60264012942332945, 0.68577342404506436, 0.25033079355980759, 0.39797070551930913, 
        0.82826213375913027, 0.062359145576879937, 0.7771330442811335, 0.47266798362074769, 0.53042894939133967, 
        0.67931445544847402, 0.0089199405481250604, 0.69526726717418907, 0.53846368838100722, 0.31590791416939568, 
        0.12589676248015103, 0.41634252470230493, 0.24460996943024771, 0.43181081658475629, 0.4385379318396021, 
        0.17545368290185004, 0.27558857317174201, 0.5627840440050571, 0.72902569739603784, 0.8041390290993391, 
        0.42072782824372057, 0.00088414778039691555, 0.56718077635604636, 0.35412412539521143, 0.70037653996220062, 
        0.39323895923582275, 0.74197391375451027, 0.89189257483142903, 0.77282369848573151, 0.030009225874839229, 
        0.35881749186327139, 0.95358671566726461, 0.74907928671993962, 0.50448610249280734, 0.6541038189123165, 
        0.56286065382048789, 0.61823466648355185, 0.14271731879204475, 0.97108650576741673, 0.34942902256060759, 
        0.73059594274725947, 0.09896507583195735, 0.14454986035715378, 0.58404656396635501, 0.16052678004389254, 
        0.065073088681911084, 0.85140878405603204, 0.79816367548436329, 0.089666796937483806, 0.80224790034073679, 
        0.89634658219198515, 0.61720503736021382, 0.33019084124462728, 0.54263376550265008, 0.64480361529147423, 
        0.30353084201314684, 0.66905880926353345, 0.94373286398348566, 0.91073953055511248, 0.36058142341120081, 
        0.72112422944932808, 0.87818703070486759, 0.36038792889610338, 0.83484749689486804, 0.48661665819768429, 
        0.77123626885803853, 0.84008558023636359, 0.39987311992663477, 0.85321810661241471, 0.53479746992027488, 
        0.830095568269972, 0.45752815918534351, 0.10422071282458001, 0.30249732643386484, 0.6609962587701097, 
        0.062898482719160542, 0.26760171638791963, 0.97180799163461518, 0.059257329537995007, 0.77265194233621171, 
        0.77194312864385828, 0.11491769669964147, 0.31909577632310526, 0.41045405946463065, 0.90073680302979975, 
        0.38857249513111092, 0.58638694595816787, 0.1095246519134776, 0.7585569690238374, 0.11571538143959192, 
        0.50466812030282293, 0.78980158445901605, 0.68368834612978047, 0.73828663470810385, 0.62169213229602271, 
        0.69272033821488721, 0.94219586866980087, 0.98183020301637258, 0.19290348394407464, 0.21809864876761531, 
        0.83784686402110475, 0.46714888674624477, 0.39770609579830318, 0.0088510419391409911, 0.4836739530032379, 
        0.46570872775310534, 0.7664784979169692, 0.49208312651520325, 0.61957821418838277, 0.49046736242332423, 
        0.32571337177324677, 0.16864996087267525, 0.062095790088121383, 0.82546993107986211, 0.65743485232014476, 
        0.37188941687768917, 0.4653498898084909, 0.93896712918490532, 0.63245159466832934, 0.40011815411020413, 
        0.1776302186222749, 0.52702162576915867, 0.60988935463782967, 0.41075855929183036, 0.63890349266831459, 
        0.044665835843709534, 0.40765635342798801, 0.074436223431312287, 0.8504648333585616, 0.56822171540794875, 
        0.99798224298597504, 0.81321248623368336, 0.3600838912565969, 0.029903610693907234, 0.044137762238141809, 
        0.7941628415570825, 0.99376097637675054, 0.28206236125098871, 0.25048519756042165, 0.21326695754941105, 
        0.98467510899919186, 0.090569868064201531, 0.018220748893168048, 0.50644215935720194, 0.90920945067672099, 
        0.68345948952202029, 0.90350016478873596, 0.36735900500570695, 0.5668388608785575, 0.94480011942449127, 
        0.17292816147805601, 0.55608769856323614, 0.4553948873847613, 0.30197359425275616, 0.32923004029983471, 
        0.87755960104052777, 0.070162531435560949, 0.20311976317979452, 0.34091532291189885, 0.1189310564471262, 
        0.73425211044387972, 0.12159250411643563, 0.095285484546552279, 0.20972724830917722, 0.20345601094974097, 
        0.50269730225177289, 0.044700692109150619, 0.019133932697505296, 0.82264213091504645, 0.49829673381437245, 
        0.10488195443995241, 0.2759217921634356, 0.41889058300789106, 0.98523961856283737, 0.86438974018520032, 
        0.81554100722230327, 0.90707994565892225, 0.67440868236499729, 0.94091017454365278, 0.1940133305357441, 
        0.51954576175257405, 0.85940967704449078, 0.39991756978205251, 0.62708964779397447, 0.84658025159836936, 
        0.29105380912190348, 0.73597825118204874, 0.68364073933081326, 0.87570589676526422, 0.4036869577846478, 
        0.82703660584106675, 0.23357447912447848, 0.6524568335507932, 0.30280212947709528, 0.002607188033771024, 
        0.43097862286581456, 0.66111947747060862, 0.63671976496667937, 0.87633880750700777, 0.99934780731976214, 
        0.2807784492675236, 0.98528861444702831, 0.78715778661477231, 0.78641117349009404, 0.26578191229160053, 
        0.52078459882657779, 0.30771959660509673, 0.5007603291995435, 0.22587137171217386, 0.15792325913289296, 
        0.28015457829085189, 0.57510594591321706, 0.46001115275866677, 0.68796472145182053, 0.48093676044374223, 
        0.65220389782797117, 0.63561550525015198, 0.86912763448870178, 0.22070052422588615, 0.40310561855008653, 
        0.77676457449358383, 0.8083525350598042, 0.19566754464155744, 0.62446522071006916, 0.62915634807264054, 
        0.82112606971320257, 0.46255741065149225, 0.80771315244807895, 0.095536425569003214, 0.85862464728276189, 
        0.51744368986499278, 0.46373011964739264, 0.91901676868165816, 0.15095447428680542, 0.080430174085812878, 
        0.90765567578817219, 0.58218851602114263, 0.76733580755119668, 0.29244013578558592, 0.78317782042995221, 
        0.96699384843315794, 0.5033119231079155, 0.21271144026427002, 0.63339054556134156, 0.66010842521376101, 
        0.021921335097626837, 0.52145651223615963, 0.83021587069757152, 0.80194991252871595, 0.57152497985287853, 
        0.09412221195641246, 0.56941896199122044, 0.61111781149553646, 0.49377830399500255, 0.04888355111728182, 
        0.16651716245540826, 0.47890383629631805, 0.83035508073832731, 0.40373429763196667, 0.83425871775375993, 
        0.058963614508009465, 0.46106489353091207, 0.52594095388955919, 0.54434421036593439, 0.95740075178066464, 
        0.75997238499507835, 0.59634036827467529, 0.35166952290075271, 0.33808883327836914, 0.40741227329784091, 
        0.85990408735354062, 0.43607913017102473, 0.66224005731149527, 0.49033320806012415, 0.045731296984582004, 
        0.30688908219897049, 0.93315186610779155, 0.37394502354455783, 0.84712108437584388, 0.58283029392730246, 
        0.2692320414194842, 0.17019877845082498, 0.63332658763327543, 0.84616211550798859, 0.69267446693858981, 
        0.65375198034193205, 0.39808651473338208, 0.23223101998869078, 0.30149789446181541, 0.63474104807897014, 
        0.53253243873994571, 0.25164888483583026, 0.62027504439505954, 0.29300351043676209, 0.99278939098531094, 
        0.434224189300777, 0.27720118950227302, 0.18414915920729014, 0.34942159062849165, 0.977902502443456, 
        0.75352830744425359, 0.47175397273270026, 0.68676595978098831, 0.13070517715060626, 0.66448643620547654, 
        0.28018087709035266, 0.84752390490513374, 0.11385109148158734, 0.45213955155544805, 0.73589064138061921, 
        0.39930969793762228, 0.028190184964336007, 0.63456749837259285, 0.79757092362068738, 0.59683134550950134, 
        0.73050797852401983, 0.35242169200447293, 0.901090096192672, 0.064261121231668605, 0.50462858013062006, 
        0.12544749365422625, 0.33817556633947121, 0.095391421905757445, 0.2181089886334362, 0.7262695205040548, 
        0.2386081003563294, 0.49791123095430212, 0.68660640488420754, 0.88375236421005665, 0.32058607585984999, 
        0.63435283044867585, 0.24490466046366599, 0.94520986435344145, 0.21063718130195519, 0.32352766890189, 
        0.85596744458314911, 0.79684938860256405, 0.70141292142688427, 0.48641616586380643, 0.42900899821699645, 
        0.5598938208578963, 0.4027904564888356, 0.014976447469657872, 0.93533549695846507, 0.047677011338697062, 
        0.98981103430484896, 0.01884275630442267, 0.0069067929859851773, 0.53754575183806197, 0.22383103585822695, 
        0.093481480925684091, 0.73429206369925759, 0.51701265649115125, 0.13047386817170015, 0.29727319822836984, 
        0.59611612624815313, 0.5591035142365901, 0.71578154212751488, 0.80945783827468842, 0.41596022991403947, 
        0.45173951629580422, 0.32110839125739066, 0.45463754629171516, 0.65945082020687651, 0.32803645162286288, 
        0.071703315167209514, 0.65821883643846379, 0.46527583870685985, 0.64768843875440152, 0.13762726718562623, 
        0.048689924403007678, 0.5704418158138993, 0.73869680528400972, 0.89728853463939551, 0.25501295394266998, 
        0.94182489629860378, 0.72761078071013663, 0.38911119737919031, 0.031357914908463735, 0.085164840629544658, 
        0.61696701793946329, 0.58737205937972448, 0.0037341740760725717, 0.44673588426359379, 0.76067011779409066, 
        0.77319552031551875, 0.94057798222424527, 0.88505262575639798, 0.27625626502472223, 0.81881969591726755, 
        0.85058377020237996, 0.63509812759842088, 0.93662467599618449, 0.44987261609174412, 0.71168238010712259, 
        0.27185363321605238, 0.57969511736963053, 0.13250625321797571, 0.24517215222305344, 0.64317223102306587, 
        0.39023834098699162, 0.36125380492633652, 0.56995978276586068, 0.0074728741145511091, 0.47800653613286737, 
        0.46495104831316492, 0.67410006996279126, 0.84215044420444185, 0.55555307466957027, 0.83278256530399308, 
        0.85259689251568105, 0.53910331532409717, 0.63535405338868656, 0.48085107056099985, 0.51649757435592392, 
        0.92238996022107722, 0.067150075800876508, 0.573992976001658, 0.15486383264558823, 0.88745991685918857, 
        0.99737120402351165, 0.96253581430931345, 0.14038836461690507, 0.30705889118306118, 0.11500610621813467, 
        0.40043059499685096, 0.68185435494852542, 0.33051257174518955, 0.36840556654970258, 0.58481872469658236, 
        0.30275112468021192, 0.49985053329386497, 0.30999322941227558, 0.34273481426899655, 0.8549710657842331, 
        0.51214486490740718, 0.42352273315063127, 0.21235975229719606, 0.5320587284743763, 0.69230468854967686, 
        0.80396238375790596, 0.19141429515775576, 0.86292771076786434, 0.56391325011177384, 0.73881698307504218, 
        0.1715642917013489, 0.32241830400263694, 0.56696968064429143, 0.50422798475424413, 0.10885057770611506, 
        0.54082920972883874, 0.98617447205926534, 0.53327879446926851, 0.95633350097118419, 0.21379313078862205, 
        0.49491968848566747, 0.59156615721164774, 0.32455893443321915, 0.82948416005989234, 0.986194904513221, 
        0.49284640072429897, 0.59934158833934559, 0.81559788214201001, 0.94914787077058005, 0.76027301773555456, 
        0.037939526833722592, 0.64394962337637041, 0.27404559748585267, 0.90782075379921379, 0.89385274645370183, 
        0.93321049066397133, 0.86531162305609688, 0.12493041267975036, 0.35832011437896472, 0.7318958831991198, 
        0.73082969954859922, 0.61189532063739538, 0.049528528789790194, 0.88851394116287596, 0.97668015327744739, 
        0.18021164030575654, 0.25762424571504461, 0.54629181510994296, 0.10793808235229507, 0.46755712486568135, 
        0.081807212065081902, 0.13418934376246261, 0.48075313141758369, 0.73257098326766279, 0.16154316281825065, 
        0.013721022485082335, 0.35135088999980035, 0.36937760083940274, 0.11381297794666412, 0.12867384032494322, 
        0.5860470311233521, 0.63111273936822343, 0.50422530741906657, 0.57851188264255948, 0.0025349061688835217, 
        0.11036384539424615, 0.26149947676889296, 0.80402167030737925, 0.50864763612622177, 0.23145028457923789, 
        0.39008008478596579, 0.53393520550255702, 0.46554181526582172, 0.92098739149674236, 0.40687644493402586, 
        0.11272995018738285, 0.93814894217611688, 0.37824942226835878, 0.43576279677518537, 0.11805823913987168, 
        0.066708415607317884, 0.75144800512925936, 0.21380281455099404, 0.083866123998163955, 0.5835095900534224, 
        0.5145463259676375, 0.093384056357003198, 0.28522361406451213, 0.4019955793904737, 0.78155500662269017, 
        0.98677846807159542, 0.19460415960270638, 0.31835213324373912, 0.16001286110375768, 0.81990322024601392, 
        0.27201230422351097, 0.051122502470773545, 0.48075383977900188, 0.89765826528610826, 0.46390025831995918, 
        0.67141721816308397, 0.43920598862091342, 0.48114510404508071, 0.094191919353094011, 0.082536934971938036, 
        0.57622636809798489, 0.7392111106300614, 0.57508106031842487, 0.93885571179149063, 0.6390637806228967, 
        0.86119792563515851, 0.90969624791903048, 0.51870822942412542, 0.98228993980279955, 0.34318706862574655, 
        0.80962646003472449, 0.14483262322331969, 0.6890774193420568, 0.0047348674540561575, 0.13546575654594828, 
        0.62345644852026605, 0.35160707306592087, 0.43207468567762164, 0.10585741714243269, 0.71543652036548533, 
        0.23927288117466561, 0.20686150957791649, 0.15348339138821054, 0.96119699457775987, 0.15516012253888922, 
        0.86439218561709752, 0.7862544395222284, 0.46429951746870435, 0.94104557391221144, 0.21777107739435175, 
        0.43933122435783778, 0.39449534854043633, 0.081955211842961839, 0.90394277363527187, 0.22929578964959041, 
        0.2765672395037222, 0.42483188950137252, 0.63342547178190989, 0.098943810262933418, 0.61326504072648502, 
        0.58034587004866878, 0.43408425362684455, 0.86338152274991331, 0.89013854355520938, 0.27504527867999506, 
        0.15142883756580017, 0.81106269613720428, 0.16535823274630812, 0.69188334496627535, 0.079998527540722009, 
        0.87418273592061713, 0.32191012983029021, 0.21068511687236979, 0.48201891568206578, 0.69582431994108962, 
        0.89011214987981102, 0.24600032543376948, 0.86543963091328813, 0.63683646264820615, 0.15149523043113877, 
        0.62822993680302885, 0.41710328809746722, 0.076671786130667607, 0.30604511797709555, 0.36366856692043315, 
        0.080170840775108498, 0.52821015980695041, 0.94449592482096789, 0.31690691191204579, 0.68239287455322484, 
        0.60532520109059162, 0.18284045633054702, 0.34521624223890068, 0.49202945713596824, 0.49270384069023998, 
        0.65209007177887202, 0.76491039635417679, 0.65192495122978378, 0.25481832797419046, 0.68378051110120541, 
        0.020165709703090506, 0.66847251218121095, 0.42861056485905924, 0.8140745926158226, 0.56292451777082309, 
        0.74192289581592008, 0.57632811276822449, 0.67332988055553966, 0.71860337835283072, 0.046754945948072502, 
        0.045274411470320075, 0.77976626542018268, 0.86984759916527077, 0.018912250314903956, 0.24159589421157524, 
        0.46713487965543221, 0.54000753823173997, 0.95625624860176606, 0.72202222060052024, 0.53901489381830925, 
        0.93953036181630867, 0.22641116794354166, 0.986855579594063, 0.34721306463864021, 0.72567792762360317, 
        0.23049273996933239, 0.12527166503846554, 0.26767085355943077, 0.083212769810237397, 0.37001315112664712, 
        0.68134655319668069, 0.69476753278919201, 0.72609845276550367, 0.41960623271939745, 0.07027106785259396, 
        0.77345763334021256, 0.40175770459681437, 0.64002772503690331, 0.60733616085901376, 0.027912209603009552, 
        0.096876723557826283, 0.58431048109955874, 0.56917774919211972, 0.05297926544096665, 0.14204672711473743, 
        0.15591848252104046, 0.44231261983542547, 0.85372337765756079, 0.18972383653728131, 0.72395472895394852, 
        0.10239787884141993, 0.30515478091402026, 0.39305376281060389, 0.21487990190568373, 0.44219538013309356, 
        0.59912468039231559, 0.59233823586214518, 0.90817380678529758, 0.1113236693353139, 0.086234374693886728, 
        0.19963427658097066, 0.68524366294123795, 0.93759909409645337, 0.62545971676758971, 0.88718879870399725, 
        0.89584427751538542, 0.09403121759370614, 0.62166510231451433, 0.85775840142019133, 0.34543456197705669, 
        0.25533544072034964, 0.18860861803730189, 0.12235215641597152, 0.35729840140629521, 0.069859342512815337, 
        0.90652833625249918, 0.29060403572793336, 0.32195827629736895, 0.79172009343481564, 0.090868516516904307, 
        0.49856743466293918, 0.75024560172426469, 0.093530363756992863, 0.85050112308846715, 0.32668122718667703, 
        0.40730328385613546, 0.837901136616507, 0.3201438168019064, 0.38020983856327084, 0.70744721962314516, 
        0.064196554631402725, 0.083042636161430705, 0.61220837489833957, 0.64953447246255402, 0.61912779362420678, 
        0.94865358336345396, 0.88637007900895215, 0.71412883941416117, 0.65716252489110083, 0.25168956266129805, 
        0.94096958455219148, 0.055335469502435286, 0.35688554073445311, 0.77042206166948279, 0.52787215200737037, 
        0.42183879569753313, 0.6544874114036856, 0.35258775456720048, 0.60222019798358128, 0.40073260156611323, 
        0.55125144629576894, 0.10554166113739205, 0.10938344325435989, 0.46129457328398038, 0.54724964497952699, 
        0.80169420925378754, 0.72366568728933545, 0.58815363277104393, 0.62117325155587522, 0.94205123496913656, 
        0.74522561504650309, 0.73252862741396685, 0.74337808059296084, 0.87051958605511159, 0.45353988968848213, 
        0.88535851197712589, 0.27630100294333015, 0.7450616417625533, 0.2658535236021915, 0.52759619001519065, 
        0.62536352377352289, 0.091298529854120014, 0.48976708989923479, 0.38835278438407572, 0.96890276293233502, 
        0.17497618765805778, 0.92511662181033638, 0.69969534838159064, 0.98411940583042212, 0.26805107828215302, 
        0.28141772965386069, 0.77188076944230688, 0.081635852777233398, 0.16715723351197109, 0.034196994817781823, 
        0.77224575274453344, 0.7012081834512911, 0.13064695467679788, 0.88192238351494789, 0.23793988610914507, 
        0.44473893675691589, 0.80976154212003082, 0.33438079488264072, 0.40719889125496711, 0.19712552807984496, 
        0.69629602039584992, 0.12724578663695785, 0.97463539894983264, 0.54195410641402453, 0.47264069290404453, 
        0.022770464120163592, 0.56983022377629067, 0.30778041820791291, 0.49563418225619893, 0.17276109362944259, 
        0.88849487690136364, 0.21994371727987883, 0.051032118967155427, 0.65064143295434995, 0.39271285772577746, 
        0.64358027863808775, 0.13135952658371819, 0.89123676626730086, 0.47997163846501056, 0.17820710240453175, 
        0.040532769123095003, 0.60269876603103389, 0.52197467550343513, 0.56038695846303344, 0.23461983831820099, 
        0.71174240078477435, 0.43226404034175125, 0.19698714270220474, 0.075229264824387077, 0.41456401301084123, 
        0.1975789078861041, 0.97122919886186931, 0.85071911683354617, 0.18213634989302085, 0.10799695945149446, 
        0.85324120496747957, 0.81021328842647033, 0.24916243270669414, 0.091997728053576866, 0.77664325197046247, 
        0.80892646882421571, 0.26655399720135597, 0.81082328697349504, 0.26110587899311155, 0.97154344625494926, 
        0.36704442646583524, 0.77625852332731982, 0.25892626208345071, 0.87858110420496094, 0.34798984407567612, 
        0.045121826827021927, 0.38890839480204487, 0.032665964533829861, 0.41405649971522029, 0.62722196660673979, 
        0.36838331998250107, 0.37739563754049987, 0.31004625778892048, 0.33260978618899961, 0.72214325971866766, 
        0.5235585085108585, 0.75167515329002232, 0.74247099476069001, 0.6544191938893813, 0.70044796445824553, 
        0.91917675258720433, 0.54041397827651827, 0.75003070127533378, 0.86517304726334943, 0.29358443189194361, 
        0.35180681914328682, 0.28069382496917195, 0.75950206187683333, 0.25272464895356195, 0.53347511518764201, 
        0.6638266826965431, 0.039449367187392115, 0.25491121868584243, 0.42782737159885853, 0.38720778931102218, 
        0.47262578114739418, 0.83200701464940741, 0.23006852750786955, 0.81364987996589444, 0.25539440415673242, 
        0.61594395421784509, 0.55392938012867243, 0.72900101083415514, 0.33283677634222286, 0.7291720013429186, 
        0.53093431666082802, 0.36277253296272449, 0.49524097792040522, 0.59985077710104528, 0.86198029528462916, 
        0.10382147024602273, 0.19237675365596862, 0.50023160163271285, 0.02563299511913808, 0.23042775457355447, 
        0.87718297954974256, 0.57261701835422407, 0.087746728929496198, 0.17871077547983671, 0.80192641210074656, 
        0.8874209260273187, 0.11298150513185545, 0.075625054355092036, 0.7318787993083633, 0.081894448424284283, 
        0.71473638495180913, 0.60396527334011707, 0.68878179741872647, 0.87724446437226455, 0.32214797786449223, 
        0.96279471613733714, 0.41779412687075101, 0.94764765320356337, 0.52928090853955911, 0.51263213958131648, 
        0.36423212209902012, 0.72429097954966748, 0.25962333405109783, 0.71402938531970661, 0.79959666654072015, 
        0.30326319553627124, 0.027379613426073313, 0.359034984669232, 0.40824378240607317, 0.078231886507900583, 
        0.27145540560928416, 0.43946955153427192, 0.84759113292230692, 0.11075005790748627, 0.66858740689775109, 
        0.20250205304563829, 0.84116934182092185, 0.022668771047668512, 0.82828466078969876, 0.81613212231493804, 
        0.9860743486505148, 0.8071327542190263, 0.37370501315628579, 0.57075022207623527, 0.20480051593145476, 
        0.068183293101268694, 0.40373730762198701, 0.73365443096463645, 0.40998209085483861, 0.66674583412969191, 
        0.36579211766256092, 0.78032063971862198, 0.66417592281392013, 0.94318550833485526, 0.56345745512504175, 
        0.412155256786636, 0.92197261516875728, 0.23259024765569891, 0.33872834110881467, 0.92559693449780989, 
        0.99333087651256435, 0.71507129167306327, 0.72481091778158957, 0.44792516778691249, 0.063383073996734263, 
        0.044442604174963707, 0.13523542991037663, 0.33393489272570465, 0.17316689543119801, 0.17130369487473418, 
        0.96784356857184739, 0.88452140440652527, 0.20691422821976091, 0.58264116432672597, 0.18837092556993396, 
        0.38087381491778638, 0.80360581269496745, 0.81374602037643817, 0.36194572750209364, 0.7739232294128966, 
        0.31486334688488293, 0.36903378060316361, 0.052102750006198351, 0.71751168313471991, 0.90184976314753307, 
        0.76525224698757333, 0.12512290781856472, 0.66991360182807669, 0.34971319880588103, 0.26472140032415581, 
        0.5896789320080138, 0.80138256751564452, 0.90356311532149935, 0.9970521541789561, 0.40920830165716948, 
        0.90309089523513264, 0.11451923235720796, 0.57969273693019741, 0.82704123110476724, 0.066508719971205998, 
        0.60815288898100328, 0.12532071322207483, 0.50175290392353755, 0.98153929124139982, 0.046430296284603267, 
        0.49934094694939324, 0.38416348316232618, 0.93389820885109964, 0.58239127881127395, 0.040952505547749096, 
        0.54890774120465924, 0.13439574766160622, 0.99632705214819661, 0.73626932653760302, 0.90935316980588743, 
        0.75224497700439885, 0.62752358848376311, 0.29300542116648476, 0.74366184500196342, 0.83735899428631178, 
        0.31845410844222033, 0.94520229783681686, 0.97035960165777557, 0.078661894253501785, 0.25233016658180629, 
        0.41346977228215187, 0.89482274881378077, 0.23544530231227601, 0.38581710989678974, 0.20414262480224243, 
        0.55051542427646072, 0.18366841135784528, 0.24764038470410488, 0.7012983991915831, 0.54852961843128889, 
        0.59391656865196851, 0.9937620125770894, 0.43737835925848678, 0.4207619500488633, 0.32928413430570802, 
        0.17032893363728885, 0.0034601954514494881, 0.89824605101544686, 0.90636335179970406, 0.46584662337238414, 
        0.23063999075241903, 0.91983454232980155, 0.16586717186193445, 0.12942528342873993, 0.28963476510274733, 
        0.76040164200335925, 0.24906026599028386, 0.24803109267985834, 0.40851717742684768, 0.85064014621449946, 
        0.41258801121525934, 0.53052038290973269, 0.30820523816637291, 0.69698071811132123, 0.67927123626033992, 
        0.68427124743976475, 0.017318539284530754, 0.2066020001186657, 0.67720606300097463, 0.27090864969948814, 
        0.68788115170718123, 0.5645384238924307, 0.4708672912916696, 0.37060207539819157, 0.63360885695726954, 
        0.10505345490921569, 0.34664276677067396, 0.44036742007856366, 0.86352781683379454, 0.20061517407490537, 
        0.24115368631150935, 0.18569289898081331, 0.86298662810544635, 0.58952467579367229, 0.44723655620440206, 
        0.15029522566129772, 0.59953499272606559, 0.68696069055421449, 0.1525399979789499, 0.089501925194469534, 
        0.055934767678013086, 0.62978435841010572, 0.83662057202313211, 0.071812754038174775, 0.93889044859901305, 
        0.13003168490164407, 0.11572992264596316, 0.27756425861305245, 0.93156517341482092, 0.41980089458239234, 
        0.6613521030909324, 0.18676494642070618, 0.21612983865291002, 0.54471917397378911, 0.88745276851521471, 
        0.70945029389688652, 0.2314313489111901, 0.20076207794228473, 0.34916502920349801, 0.11078022097228346, 
        0.40915979356077758, 0.26775999593128752, 0.30190334265677565, 0.26044811107621046, 0.32888416161380807, 
        0.98419267323379178, 0.88872252432413856, 0.29930427626709011, 0.61400411024689339, 0.83915227674529369, 
        0.42261589933516963, 0.1017188550226662, 0.83864830762865439, 0.81990295572287919, 0.4846073955937884, 
        0.1459040845847015, 0.15544658008147438, 0.79904278977813248, 0.56354876601426218, 0.93406376894682852, 
        0.13203724393077954, 0.56964484634900159, 0.05046768924098477, 0.33744629435108986, 0.13681270647575827, 
        0.43580320203875078, 0.83242415623445609, 0.78005597473611665, 0.73058788324621382, 0.76017448920547492, 
        0.78690484216768852, 0.30617108350543254, 0.91944384664488288, 0.0094319271091349854, 0.78890693082644958, 
        0.22854307111251204, 0.41589624428394578, 0.53192686005496292, 0.59810328526841605, 0.2393878734574566, 
        0.18958028449729736, 0.46322805007130707, 0.24488528393514564, 0.60913461976456573, 0.79205748295449951, 
        0.4584806149909697, 0.29718052883052404, 0.050690753290428603, 0.8658242200794144, 0.80711264131569527, 
        0.1982138367308528, 0.53287058156891143, 0.93534332358606065, 0.15971852755055305, 0.6843963046812207, 
        0.11379397624458853, 0.67976597996393484, 0.59009545266398145, 0.24250977362562587, 0.65195822540618598, 
        0.42092036244055442, 0.96299412936199369, 0.40299278030795893, 0.90418274352059913, 0.49801867009116418, 
        0.51345477896587632, 0.084178150889594017, 0.43587425424194626, 0.51025030621966772, 0.42736356435116218, 
        0.92626882046003001, 0.53860930873013646, 0.074642826872167722, 0.5899834955647465, 0.13191741089662212, 
        0.71557605506174848, 0.24883245006804344, 0.58712125357487266, 0.62563130171005454, 0.80758181291833075, 
        0.080751468331655918, 0.63572662446250527, 0.41683329066485597, 0.91358229622671416, 0.70810803916056742, 
        0.53332867452637367, 0.50659783772888844, 0.6594345140598783, 0.52988632990809292, 0.44907784809515294, 
        0.9032132523406331, 0.30665803355024179, 0.40535612438180002, 0.19114758360844308, 0.39552514841002373, 
        0.51462216727100141, 0.28051221187192321, 0.75020085677384074, 0.23436511894453971, 0.28896109257353375, 
        0.27412362676337998, 0.86302342725010806, 0.73432374355672403, 0.08420760396550997, 0.84476922956735989, 
        0.80984239176407713, 0.43660823841120067, 0.22775568513462785, 0.32371588418642872, 0.59995080820640978, 
        0.91997245514489667, 0.33021196870070813, 0.082599678776028806, 0.73793774332999362, 0.11409493411782923, 
        0.95720376607002411, 0.034638896216645731, 0.61597115828304139, 0.47201873231023606, 0.14460879054539832, 
        0.98203029326422864, 0.43656092732235585, 0.10114681654760616, 0.31923353752360906, 0.5064166919011337, 
        0.60947802762553738, 0.94168883564009187, 0.55963453495488191, 0.13470099716817363, 0.84776674629685433, 
        0.49675635138167151, 0.73309094276207931, 0.12894784356539102, 0.081962565639087881, 0.86215710467659123, 
        0.97686173471456383, 0.49780911914330428, 0.3690540169440979, 0.39868610239440949, 0.91592192366925018, 
        0.055802601416877318, 0.48795533246335987, 0.18258001511291022, 0.89502528557720207, 0.83856193996926631, 
        0.75751264524884854, 0.18388966943170759, 0.28069565367185434, 0.43861386495094834, 0.047102183001541009, 
        0.47831329603583539, 0.013841387917107983, 0.71167180140679753, 0.44881256737764885, 0.60177939636268163, 
        0.68542825117271478, 0.25904896044136749, 0.72454400990948797, 0.36952525622027377, 0.26839575857940856, 
        0.71735341762611404, 0.84139550986695233, 0.75229069505815138, 0.15914259751969984, 0.024672785579256562, 
        0.079349011321637075, 0.031354499215009213, 0.45601121961247348, 0.53684032389739311, 0.79906815374407714, 
        0.81739465328882588, 0.20728570446848704, 0.26952703383410759, 0.19718100136913086, 0.33697429252955158, 
        0.78912818183654942, 0.82532895218598923, 0.078020698626848128, 0.14310012798540694, 0.23430698891807089, 
        0.11775175391367454, 0.54829862360052539, 0.79096094009873208, 0.98004468790061927, 0.88169339184933149, 
        0.35469525898456578, 0.55895301369644179, 0.329042827385293, 0.23913722500594026, 0.98577143492824737, 
        0.29601138990480536, 0.07955164063693676, 0.53549049389844061, 0.47966445412741487, 0.7779846903056995, 
        0.71961525182521213, 0.19412157722070167, 0.30950428696038612, 0.83002798751345153, 0.13363383872168688, 
        0.42047848831552082, 0.26376460109185107, 0.40915632102809063, 0.2440575807346812, 0.90493161161709224, 
        0.98921979049818276, 0.58647170173910879, 0.21472791221416765, 0.30473934861431795, 0.43547413584265415, 
        0.94401340965736091, 0.86631647764430486, 0.85412828174178768, 0.8244564787247346, 0.33520287616912703, 
        0.052361790842057498, 0.50671575188337359, 0.085798913173891878, 0.56173570120815319, 0.043764762818563119, 
        0.86389986777201289, 0.70577603935442657, 0.16646859649628842, 0.43941426676177464, 0.19620063428297585, 
        0.92705028659773658, 0.87262699389727461, 0.35423498965888101, 0.80324293923025025, 0.65991074746249168, 
        0.25767240246771128, 0.920476754757523, 0.41848190706713462, 0.28268272042194931, 0.65187827970425616, 
        0.87569611473262943, 0.18668115940444885, 0.073874148181669597, 0.59465140112538895, 0.73909978009828614, 
        0.19588073821427532, 0.95397292581711812, 0.79885287420316775, 0.094107768657010649, 0.88394764700682926, 
        0.78746124669096473, 0.99685058301279406, 0.73265480733433996, 0.23041159318655446, 0.94611103103226868, 
        0.35495472875907663, 0.072447532113448565, 0.80196553415737859, 0.14197628853990607, 0.34824981161347091, 
        0.099015326705326867, 0.15568759023491507, 0.85031138480032897, 0.88531466124903613, 0.39808712998869256, 
        0.043442497564549765, 0.123760470124888, 0.51689045982814208, 0.82544032530556266, 0.65146853939818539, 
        0.11090823216218859, 0.96412512642914838, 0.15653512292433835, 0.53115300721213976, 0.6211808849251379, 
        0.097451670978136784, 0.50006166117944839, 0.9108750635873788, 0.71469930602841236, 0.025126769560798401, 
        0.94385027517387687, 0.29664124253131874, 0.22934200683301875, 0.060699256168083027, 0.15556139529349577, 
        0.82262645433969883, 0.24755943307858508, 0.54247045490073509, 0.15295025755448899, 0.45280755821041407, 
        0.11086831226666627, 0.25503009380284847, 0.0692825070731522, 0.22911568005139116, 0.74690019816138342, 
        0.39930253020590878, 0.19994056876072319, 0.035447028789436619, 0.87635057117532766, 0.5372655588413866, 
        0.76302801687794219, 0.80243992825299371, 0.78055091109395214, 0.93322184522719964, 0.67806998380292205, 
        0.014533449632591289, 0.77256711791088439, 0.01910868607595928, 0.19888344035106775, 0.81956773375962721, 
        0.95644806825533668, 0.30603754515613768, 0.07773574690039986, 0.72751622096710977, 0.4743418557784298, 
        0.64792588569203935, 0.34204938987130751, 0.882601062296831, 0.51178385797226111, 0.91422033035731975, 
        0.02164992936957022, 0.12915971861501307, 0.98743226320192634, 0.54196788002270635, 0.37239615825489536, 
        0.8842547503670517, 0.91185241636173298, 0.62588793934510312, 0.17585344457127827, 0.23167366235433851, 
        0.04625306478077218, 0.78429793378389867, 0.18815871880105317, 0.20877779227520277, 0.81039288708310031, 
        0.3280374444614953, 0.71031632307739256, 0.51023894051367447, 0.24205668302386019, 0.45936077590715141, 
        0.56880140098643794, 0.2674622220872831, 0.59943114029810252, 0.58297194656331097, 0.79251238432694571, 
        0.49141894161116983, 0.072771614076971858, 0.93993038803900886, 0.12976634035479773, 0.56462020211243713, 
        0.7155265736912706, 0.72399863565470679, 0.37339561032714874, 0.34603157700529774, 0.68002680803674109, 
        0.014934427469702527, 0.89619503182361382, 0.55603673789261632, 0.14665522835391243, 0.44007623722652989, 
        0.33298397637373522, 0.2661200959067036, 0.58958482042418492, 0.9226184288649133, 0.45333867099104386, 
        0.74536891860112808, 0.42837043493466509, 0.31369851877092669, 0.53622674980080642, 0.24390788308900513, 
        0.84975246529337745, 0.93351175099946215, 0.21993991324764806, 0.2325696404911739, 0.93886045331716539, 
        0.71513117363018086, 0.20708637076347114, 0.69568081018241656, 0.009232493752688864, 0.98417832485314771, 
        0.83814635492924161, 0.1785904086935981, 0.46473149557146831, 0.69778944227932471, 0.63046655337416135, 
        0.81920506372246904, 0.55588492236634601, 0.51424627285128621, 0.18960495291843427, 0.40380216034248884, 
        0.13092187424275048, 0.16947270384689883, 0.31017080146518561, 0.4084664790012893, 0.67827264195454084, 
        0.99697213326648382, 0.12400170428189061, 0.09296824988802932, 0.10165552364253716, 0.38907352158866026, 
        0.85498818410879807, 0.88842061773607206, 0.80652665865105266, 0.20821257278293714, 0.85241951396703541, 
        0.76610748317634503, 0.3183969842089327, 0.92653911283280221, 0.07673970098439642, 0.47445760245014412, 
        0.63755499469143406, 0.14750667688720775, 0.80366408594223104, 0.39355686664236611, 0.56024522393007814, 
        0.073920833438577249, 0.53814286594898575, 0.27788729125679446, 0.78629649131614321, 0.87861922513322788, 
        0.6484569084627152, 0.16733628195652517, 0.47958368196263201, 0.78590958498089059, 0.94829338682197362, 
        0.72487394324684806, 0.15182513474154358, 0.52913356375801501, 0.048874375303603168, 0.50896271789850189, 
        0.66040311324339673, 0.34545401182269253, 0.23730548549197561, 0.98630782107806159, 0.48045104827528484, 
        0.98347095590414901, 0.74363728552510233, 0.27883447944931783, 0.61785346925488915, 0.42879532162813772, 
        0.64426597746330438, 0.97369749997112631, 0.96600174019679219, 0.92040706086072155, 0.60866106386956731, 
        0.71958457115219865, 0.1325050559951908, 0.5082253619146555, 0.039029532775519238, 0.1718089448997504, 
        0.61853394765828162, 0.006882926956886859, 0.89817073468801256, 0.88780851364589908, 0.32960520016293504, 
        0.26180951350459503, 0.78395563002623803, 0.34679325286811613, 0.37053408099930141, 0.39148958639443499, 
        0.59966286630725252, 0.32121033439903846, 0.039873190703289207, 0.33652717839867141, 0.86834542563935346, 
        0.50704610370100367, 0.022412742284776854, 0.48063577295116633, 0.49166293405796724, 0.69116951716951935, 
        0.076308843979698615, 0.82027226812551857, 0.29858717970705984, 0.87064976543119244, 0.85543891826501661, 
        0.17540072079190416, 0.79647573123893678, 0.26624812055597946, 0.44830893006401484, 0.46780499731406544, 
        0.5469227011454707, 0.9973584398533768, 0.78601474327238385, 0.093004985508627458, 0.15684201100839745, 
        0.36878892229835025, 0.9880337501323786, 0.96778576710600817, 0.92909578216706734, 0.33253047810433833, 
        0.45231635274905613, 0.33191489085604231, 0.045545532705741465, 0.45329126878330328, 0.43675805687113467, 
        0.556035112715205, 0.83787156094190118, 0.55683971918434971, 0.59208101192568252, 0.24819861179130931, 
        0.45041023107418643, 0.66872837539647167, 0.96220923614670006, 0.59576473685556564, 0.061814076263803308, 
        0.24320963698256981, 0.31921253256120252, 0.87229763035390806, 0.52267186814711386, 0.19467582843147668, 
        0.34857943255333246, 0.10215131814962919, 0.89708566686272917, 0.38500812486792024, 0.93088953389487883, 
        0.86295029984633831, 0.2439465469182025, 0.23534038389354994, 0.12506678241911962, 0.6611214192861421, 
        0.91539747710308284, 0.848988693418532, 0.77604879345100275, 0.12489458319860436, 0.61128974908294587, 
        0.55709441825076444, 0.41152609895803671, 0.19352287177169303, 0.17259214622681629, 0.93698829628083846, 
        0.16348685054217915, 0.80039941744509369, 0.66679435644802676, 0.79581241042581841, 0.7075057904706441, 
        0.34312953375928679, 0.27970636362204004, 0.25512909663072603, 0.6742540677350517, 0.22096238339178775, 
        0.59995313283452134, 0.20041260028198749, 0.51023816990448467, 0.40957917667910948, 0.28001182118053158, 
        0.77023430877409216, 0.58071335455681106, 0.039287024177821817, 0.77402488116753321, 0.7352607199399166, 
        0.8719954853076779, 0.84140520907935357, 0.25390832793672402, 0.17733206760924469, 0.2771793719957405, 
        0.71311392549450781, 0.52715006583513913, 0.95232059237845235, 0.5630744901787148, 0.8711664070926477, 
        0.47830433436864794, 0.49276046497104264, 0.43850056931972992, 0.38702523120017429, 0.54663790169616866, 
        0.49443417136791479, 0.36323266634017704, 0.020861914423620176, 0.54613295286668295, 0.14563105208436955, 
        0.90770067016886347, 0.65685539120879066, 0.81447019020369993, 0.50443884524161642, 0.42764953662709471, 
        0.25482068527165103, 0.074344153674740454, 0.74112150942210664, 0.098769325161996768, 0.93901864636334054, 
        0.11571902052057959, 0.18192785965265679, 0.77805950579440686, 0.60029541705555389, 0.67407179383231974, 
        0.92624514540471847, 0.45455067206510202, 0.53973282019943203, 0.28291062245268339, 0.73131460205820997, 
        0.99718919458656496, 0.027041709464680652, 0.22430739052993998, 0.50029448410813826, 0.27049514706717637, 
        0.20188491765571381, 0.85260151929933436, 0.0090541387202227597, 0.21837400949955077, 0.99038307288617533, 
        0.17545752868980546, 0.55492484502821982, 0.31436636758105418, 0.05662938714540422, 0.65154632288432479, 
        0.92486753875478156, 0.37252708657628353, 0.59417011494948713, 0.42964840396437576, 0.97021024893888774, 
        0.81170049121792198, 0.34255395824112367, 0.38886735578155474, 0.96395573388982259, 0.67558518253293287, 
        0.11202107319335242, 0.53354005907423607, 0.36251187030299903, 0.42907345603661029, 0.080438243981869606, 
        0.83853681475777209, 0.9778868393277993, 0.86520223084055492, 0.076953305247298376, 0.80984433910742504, 
        0.2270966926295872, 0.20603310456097867, 0.12063801742071156, 0.69167951122097482, 0.4303900904741047, 
        0.32663699536173807, 0.20451635305373084, 0.31292589342558164, 0.38013182678495228, 0.33671238592439479, 
        0.62360534229972009, 0.79863318559432228, 0.42687209192101561, 0.29034755666804246, 0.071523396151897733, 
        0.9816659560824823, 0.94364383416924258, 0.41906403123625768, 0.81031403518857803, 0.73165876222640591, 
        0.76166135404979718, 0.91055040234045626, 0.11819841811552867, 0.75152043495713938, 0.72990328590567466, 
        0.28969884954374869, 0.26229289492659613, 0.72909850004247589, 0.95801653097992023, 0.84635949003403876, 
        0.24657114541122693, 0.70161267290835538, 0.93180620858563001, 0.56449941500605805, 0.34550984751826341, 
        0.63495265939770174, 0.65276189954977459, 0.35223546526351179, 0.75120096161782191, 0.51431094474025785, 
        0.54074857082835215, 0.28211276142977093, 0.27367181183700184, 0.46511377743512794, 0.32539447543941891, 
        0.82607870219698776, 0.096880695990769761, 0.090671944519704839, 0.94192898055291185, 0.27420374531363767, 
        0.88108242605092446, 0.66160777192357556, 0.36935948481451719, 0.10913462513386829, 0.98792833268666347, 
        0.42382969474728793, 0.24677687048129471, 0.70701393326868356, 0.85093276905783322, 0.078514247545707461, 
        0.25921124122467543, 0.64911503682226535, 0.075071312126054446, 0.63503216589503442, 0.07923834681382691, 
        0.97193237123557741, 0.99003358971350064, 0.66754643224416754, 0.2681223222052902, 0.43642095289942651, 
        0.68115272644445746, 0.44245705855959705, 0.28987689537302641, 0.42559170165839322, 0.3561382146325538, 
        0.68916641004531742, 0.22467371862870622, 0.79146626803991404, 0.77553373586347973, 0.21277334452558461, 
        0.35091548429798114, 0.49557161081035606, 0.89997541697014927, 0.24936434659988871, 0.58610411707425425, 
        0.78068967779620335, 0.13169922076870355, 0.451333207583964, 0.8576694798306328, 0.68310139959906802, 
        0.90095940760899307, 0.36176749598024682, 0.56351162480316708, 0.61473039250843931, 0.46989281361357538, 
        0.79179879425959898, 0.2159721760483122, 0.97119832272227824};

    for (int i = 0; i < 3263; i++)
    {
        Chi_phi += amp[i]/100*(dr_1d[i]/100*sin(2*M_PI*(kx[i]*x_1/0.05 + ky[i]*x_2/0.05)));
    }
    
    return Chi_phi;
}
