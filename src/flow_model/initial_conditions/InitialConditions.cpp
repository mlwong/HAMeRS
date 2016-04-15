#include "flow_model/initial_conditions/InitialConditions.hpp"

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
        
        switch (d_flow_model)
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
            case FOUR_EQN_SHYUE:
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
                
                boost::shared_ptr<pdat::CellData<double> > mass_fraction(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_mass_fraction, data_context)));
                
#ifdef DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(density);
                TBOX_ASSERT(momentum);
                TBOX_ASSERT(total_energy);
                TBOX_ASSERT(mass_fraction);
#endif
                
                if (d_dim == tbox::Dimension(1))
                {
                    // NOT YET IMPLEMENTED.
                    
                    TBOX_ERROR(d_object_name
                        << ": "
                        << "Cannot initialize data for unknown 1D problem"
                        << " for muti-species flow with four-equation model by Shyue with name = '"
                        << d_project_name
                        << "'."
                        << std::endl);
                }
                else if (d_dim == tbox::Dimension(2))
                {
                    if (d_project_name == "2D advection of material interface")
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
                        
                        double* rho       = density->getPointer(0);
                        double* rho_u     = momentum->getPointer(0);
                        double* rho_v     = momentum->getPointer(1);
                        double* E         = total_energy->getPointer(0);
                        double* Y_1       = mass_fraction->getPointer(0);
                        double* Y_2       = mass_fraction->getPointer(1);
                        
                        const double x_a = 1.0/3*(domain_xlo[0] + domain_xhi[0]);
                        const double x_b = 2.0/3*(domain_xlo[0] + domain_xhi[0]);
                        
                        const double y_a = 1.0/3*(domain_xlo[1] + domain_xhi[1]);
                        const double y_b = 2.0/3*(domain_xlo[1] + domain_xhi[1]);
                        
                        // Material initial conditions.
                        double gamma_m  = d_equation_of_state->
                            getSpeciesThermodynamicProperty(
                                "gamma",
                                0);
                        double rho_m    = 10.0;
                        double u_m      = 0.5;
                        double v_m      = 0.5;
                        double p_m      = 1.0/1.4;
                        
                        // Ambient initial conditions.
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
                                    rho[idx_cell]     = rho_m;
                                    rho_u[idx_cell]   = rho_m*u_m;
                                    rho_v[idx_cell]   = rho_m*v_m;
                                    E[idx_cell]       = p_m/(gamma_m - 1.0) +
                                        0.5*rho_m*(u_m*u_m + v_m*v_m);
                                    Y_1[idx_cell]     = 1.0;
                                    Y_2[idx_cell]     = 0.0;
                                }
                                else
                                {
                                    rho[idx_cell]     = rho_a;
                                    rho_u[idx_cell]   = rho_a*u_a;
                                    rho_v[idx_cell]   = rho_a*v_a;
                                    E[idx_cell]       = p_a/(gamma_a - 1.0) +
                                        0.5*rho_a*(u_a*u_a + v_a*v_a);
                                    Y_1[idx_cell]     = 0.0;
                                    Y_2[idx_cell]     = 1.0;
                                }
                            }
                        }
                    }
                    else
                    {
                        TBOX_ERROR(d_object_name
                            << ": "
                            << "Cannot initialize data for unknown 2D problem"
                            << " for muti-species flow with four-equation model by Shyue with name = '"
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
                        
                        double* rho       = density->getPointer(0);
                        double* rho_u     = momentum->getPointer(0);
                        double* rho_v     = momentum->getPointer(1);
                        double* rho_w     = momentum->getPointer(2);
                        double* E         = total_energy->getPointer(0);
                        double* Y_1       = mass_fraction->getPointer(0);
                        double* Y_2       = mass_fraction->getPointer(1);
                        
                        const double x_a = 1.0/3*(domain_xlo[0] + domain_xhi[0]);
                        const double x_b = 2.0/3*(domain_xlo[0] + domain_xhi[0]);
                        
                        const double y_a = 1.0/3*(domain_xlo[1] + domain_xhi[1]);
                        const double y_b = 2.0/3*(domain_xlo[1] + domain_xhi[1]);
                        
                        const double z_a = 1.0/3*(domain_xlo[2] + domain_xhi[2]);
                        const double z_b = 2.0/3*(domain_xlo[2] + domain_xhi[2]);
                        
                        // Material initial conditions.
                        double gamma_m  = d_equation_of_state->
                            getSpeciesThermodynamicProperty(
                                "gamma",
                                0);
                        double rho_m    = 10.0;
                        double u_m      = 0.5;
                        double v_m      = 0.5;
                        double w_m      = 0.5;
                        double p_m      = 1.0/1.4;
                        
                        // Ambient initial conditions.
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
                                        rho[idx_cell]     = rho_m;
                                        rho_u[idx_cell]   = rho_m*u_m;
                                        rho_v[idx_cell]   = rho_m*v_m;
                                        rho_w[idx_cell]   = rho_m*w_m;
                                        E[idx_cell]       = p_m/(gamma_m - 1.0) +
                                            0.5*rho_m*(u_m*u_m + v_m*v_m + w_m*w_m);
                                        Y_1[idx_cell]     = 1.0;
                                        Y_2[idx_cell]     = 0.0;
                                    }
                                    else
                                    {
                                        rho[idx_cell]     = rho_a;
                                        rho_u[idx_cell]   = rho_a*u_a;
                                        rho_v[idx_cell]   = rho_a*v_a;
                                        rho_w[idx_cell]   = rho_a*w_a;
                                        E[idx_cell]       = p_a/(gamma_a - 1.0) +
                                            0.5*rho_a*(u_a*u_a + v_a*v_a + w_a*w_a);
                                        Y_1[idx_cell]     = 0.0;
                                        Y_2[idx_cell]     = 1.0;
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
                            << " for muti-species flow with four-equation model by Shyue with name = '"
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
                        const double C_epsilon = 1.5;
                        // const double epsilon_i = C_epsilon*396.875*10-6; // epsilon_i for smoothing interface.
                        const double epsilon_i = C_epsilon*pow(dx[0]*dx[1]*dx[2], 1.0/3.0); // epsilon_i for smoothing interface.
                        
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
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "d_flow_model '"
                    << d_flow_model
                    << "' not yet implemented."
                    << std::endl);
            }
        }
    }
}
