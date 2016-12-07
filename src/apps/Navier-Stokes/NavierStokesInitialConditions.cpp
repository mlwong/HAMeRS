#include "apps/Navier-Stokes/NavierStokesInitialConditions.hpp"

#include <boost/math/constants/constants.hpp>

/*
 * Set the data on the patch interior to some initial values,
 * depending on the flow problems and flow models.
 */
void
NavierStokesInitialConditions::initializeDataOnPatch(
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
        
        switch (d_flow_model_type)
        {
            case FLOW_MODEL::SINGLE_SPECIES:
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
                        
                        double gamma = 1.4;
                        
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
                    else if (d_project_name == "2D Couette flow in x-direction")
                    {
                        const double gamma = 1.4;
                        const double rho_avg = 1.159750086791891;
                        const double u_w = 69.445;
                        const double p_c = 1.0e5;
                        
                        double* rho   = density->getPointer(0);
                        double* rho_u = momentum->getPointer(0);
                        double* rho_v = momentum->getPointer(1);
                        double* E     = total_energy->getPointer(0);
                        
                        const double H = 1.0;
                        
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
                                
                                double phi = x[1]/H;
                                double u = u_w*phi;
                                
                                rho[idx_cell] = rho_avg;
                                rho_u[idx_cell] = rho[idx_cell]*u;
                                rho_v[idx_cell] = 0.0;
                                E[idx_cell] = p_c/(gamma - 1.0) + 0.5*rho_avg*u*u;
                            }
                        }
                    }
                    else if (d_project_name == "2D Couette flow in y-direction")
                    {
                        const double gamma = 1.4;
                        const double rho_avg = 1.159750086791891;
                        const double v_w = 69.445;
                        const double p_c = 1.0e5;
                        
                        double* rho   = density->getPointer(0);
                        double* rho_u = momentum->getPointer(0);
                        double* rho_v = momentum->getPointer(1);
                        double* E     = total_energy->getPointer(0);
                        
                        const double H = 1.0;
                        
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
                                
                                double phi = x[0]/H;
                                double v = v_w*phi;
                                
                                rho[idx_cell] = rho_avg;
                                rho_u[idx_cell] = 0.0;
                                rho_v[idx_cell] = rho[idx_cell]*v;
                                E[idx_cell] = p_c/(gamma - 1.0) + 0.5*rho_avg*v*v;
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
                        
                        double gamma = 1.4;
                        
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
                        
                        double gamma = 5.0/3.0;
                        
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
                    else if (d_project_name == "3D Couette flow in x-direction")
                    {
                        const double gamma = 1.4;
                        const double rho_avg = 1.159750086791891;
                        const double u_w = 69.445;
                        const double p_c = 1.0e5;
                        
                        const double H = 1.0;
                        
                        double* rho   = density->getPointer(0);
                        double* rho_u = momentum->getPointer(0);
                        double* rho_v = momentum->getPointer(1);
                        double* rho_w = momentum->getPointer(2);
                        double* E     = total_energy->getPointer(0);
                        
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
                                    
                                    double phi = x[1]/H;
                                    double u = u_w*phi;
                                    
                                    rho[idx_cell] = rho_avg;
                                    rho_u[idx_cell] = rho[idx_cell]*u;
                                    rho_v[idx_cell] = 0.0;
                                    rho_w[idx_cell] = 0.0;
                                    E[idx_cell] = p_c/(gamma - 1.0) + 0.5*rho_avg*u*u;
                                }
                            }
                        }
                    }
                    else if (d_project_name == "3D Couette flow in y-direction")
                    {
                        const double gamma = 1.4;
                        const double rho_avg = 1.159750086791891;
                        const double v_w = 69.445;
                        const double p_c = 1.0e5;
                        
                        const double H = 1.0;
                        
                        double* rho   = density->getPointer(0);
                        double* rho_u = momentum->getPointer(0);
                        double* rho_v = momentum->getPointer(1);
                        double* rho_w = momentum->getPointer(2);
                        double* E     = total_energy->getPointer(0);
                        
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
                                    
                                    double phi = x[0]/H;
                                    double v = v_w*phi;
                                    
                                    rho[idx_cell] = rho_avg;
                                    rho_u[idx_cell] = 0.0;
                                    rho_v[idx_cell] = rho[idx_cell]*v;
                                    rho_w[idx_cell] = 0.0;
                                    E[idx_cell] = p_c/(gamma - 1.0) + 0.5*rho_avg*v*v;
                                }
                            }
                        }
                    }
                    else if (d_project_name == "3D Couette flow in z-direction")
                    {
                        const double gamma = 1.4;
                        const double rho_avg = 1.159750086791891;
                        const double w_w = 69.445;
                        const double p_c = 1.0e5;
                        
                        const double H = 1.0;
                        
                        double* rho   = density->getPointer(0);
                        double* rho_u = momentum->getPointer(0);
                        double* rho_v = momentum->getPointer(1);
                        double* rho_w = momentum->getPointer(2);
                        double* E     = total_energy->getPointer(0);
                        
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
                                    
                                    double phi = x[1]/H;
                                    double w = w_w*phi;
                                    
                                    rho[idx_cell] = rho_avg;
                                    rho_u[idx_cell] = 0.0;
                                    rho_v[idx_cell] = 0.0;
                                    rho_w[idx_cell] = rho[idx_cell]*w;
                                    E[idx_cell] = p_c/(gamma - 1.0) + 0.5*rho_avg*w*w;
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
            case FLOW_MODEL::FOUR_EQN_CONSERVATIVE:
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
                    if (d_project_name == "2D binary diffusion in x-direction")
                    {
                        if (d_num_species != 2)
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "Please provide only two-species for the problem"
                                << " '2D binary diffustion'."
                                << std::endl);
                        }
                        
                        // Characteristic length of the problem.
                        const double K = 0.1*0.1;
                        
                        double* rho_Y_1   = partial_density->getPointer(0);
                        double* rho_Y_2   = partial_density->getPointer(1);
                        double* rho_u     = momentum->getPointer(0);
                        double* rho_v     = momentum->getPointer(1);
                        double* E         = total_energy->getPointer(0);
                        
                        // Initial conditions of mixture.
                        const double rho = 1.0;
                        const double u   = 0.0;
                        const double v   = 0.0;
                        const double p   = 1.0;
                        
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
                                
                                const double Y_1 = exp(-x[0]*x[0]/K);
                                const double Y_2 = 1.0 - Y_1;
                                
                                const double gamma = 1.4;
                                
                                rho_Y_1[idx_cell] = rho*Y_1;
                                rho_Y_2[idx_cell] = rho*Y_2;
                                rho_u[idx_cell]   = rho*u;
                                rho_v[idx_cell]   = rho*v;
                                E[idx_cell]       = p/(gamma - 1.0) + 0.5*rho*(u*u + v*v);
                            }
                        }
                    }
                    else if (d_project_name == "2D binary diffusion in y-direction")
                    {
                        if (d_num_species != 2)
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "Please provide only two-species for the problem"
                                << " '2D binary diffustion'."
                                << std::endl);
                        }
                        
                        // Characteristic length of the problem.
                        const double K = 0.1*0.1;
                        
                        double* rho_Y_1   = partial_density->getPointer(0);
                        double* rho_Y_2   = partial_density->getPointer(1);
                        double* rho_u     = momentum->getPointer(0);
                        double* rho_v     = momentum->getPointer(1);
                        double* E         = total_energy->getPointer(0);
                        
                        // Initial conditions of mixture.
                        const double rho = 1.0;
                        const double u   = 0.0;
                        const double v   = 0.0;
                        const double p   = 1.0;
                        
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
                                
                                const double Y_1 = exp(-x[1]*x[1]/K);
                                const double Y_2 = 1.0 - Y_1;
                                
                                const double gamma = 1.4;
                                
                                rho_Y_1[idx_cell] = rho*Y_1;
                                rho_Y_2[idx_cell] = rho*Y_2;
                                rho_u[idx_cell]   = rho*u;
                                rho_v[idx_cell]   = rho*v;
                                E[idx_cell]       = p/(gamma - 1.0) + 0.5*rho*(u*u + v*v);
                            }
                        }
                    }
                    else if (d_project_name == "2D shock-bubble interaction")
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
                                    
                                    const double gamma = 1.4;
                                    
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
                                    
                                    const double c_p_air = 1.4;
                                    const double c_p_He = 1.648;
                                    
                                    const double c_v_air = 1.0;
                                    const double c_v_He = 1.0;
                                    
                                    const double gamma =
                                        (Y_1_i*c_p_He + Y_2_i*c_p_air)/(Y_1_i*c_v_He + Y_2_i*c_v_air);
                                    
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
                    else if (d_project_name == "2D Poggi's RMI 0")
                    {
                        if (d_num_species != 2)
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "Please provide only two-species for the problem"
                                << " '2D Poggi's RMI'."
                                << std::endl);
                        }
                        
                        // Characteristic length of the initial interface thickness.
                        const double epsilon_i = 0.001;
                        
                        double* rho_Y_0   = partial_density->getPointer(0);
                        double* rho_Y_1   = partial_density->getPointer(1);
                        double* rho_u     = momentum->getPointer(0);
                        double* rho_v     = momentum->getPointer(1);
                        double* E         = total_energy->getPointer(0);
                        
                        // species 0: SF6.
                        // species 1: air.
                        const double gamma_0 = 1.09312;
                        const double gamma_1 = 1.39909;
                        
                        const double c_p_SF6 = 668.286;
                        const double c_p_air = 1040.50;
                        
                        const double c_v_SF6 = 611.359;
                        const double c_v_air = 743.697;
                        
                        NULL_USE(gamma_1);
                        
                        // Unshocked SF6.
                        const double rho_unshocked = 5.97286552525647;
                        const double u_unshocked   = 0.0;
                        const double v_unshocked   = 0.0;
                        const double p_unshocked   = 101325.0;
                        
                        // Shocked SF6.
                        const double rho_shocked = 11.9708247309869;
                        const double u_shocked   = 98.9344103891513;
                        const double v_shocked   = 0.0;
                        const double p_shocked   = 218005.430874;
                        
                        // Air.
                        const double rho_air = 1.14560096494103;
                        const double u_air   = 0.0;
                        const double v_air   = 0.0;
                        const double p_air   = 101325.0;
                        
                        // Shock hits the interface after 0.1 ms.
                        const double L_x_shock = 0.180254509478037;
                        const double L_x_interface = 0.2;
                        
                        for (int j = 0; j < patch_dims[1]; j++)
                        {
                            for (int i = 0; i < patch_dims[0]; i++)
                            {
                                // Compute the linear index.
                                int idx_cell = i + j*patch_dims[0];
                                
                                // Compute the coordinates.
                                double x[2];
                                x[0] = patch_xlo[0] + (i + 0.5)*dx[0];
                                x[1] = patch_xlo[1] + (j + 0.5)*dx[1];
                                
                                double S = 0.0;
                                
                                if (x[0] < L_x_shock)
                                {
                                    rho_Y_0[idx_cell] = rho_shocked;
                                    rho_Y_1[idx_cell] = 0.0;
                                    rho_u[idx_cell] = rho_shocked*u_shocked;
                                    rho_v[idx_cell] = rho_shocked*v_shocked;
                                    E[idx_cell]     = p_shocked/(gamma_0 - 1.0) +
                                        0.5*rho_shocked*(u_shocked*u_shocked + v_shocked*v_shocked);
                                }
                                else
                                {
                                    const double f_sm = 0.5*(1.0 + erf((x[0] - (L_x_interface + S))/epsilon_i));
                                    
                                    // Smooth the primitive variables.
                                    const double rho_Y_0_i = rho_unshocked*(1.0 - f_sm);
                                    const double rho_Y_1_i = rho_air*f_sm;
                                    const double u_i = u_unshocked*(1.0 - f_sm) + u_air*f_sm;
                                    const double v_i = v_unshocked*(1.0 - f_sm) + v_air*f_sm;
                                    const double p_i = p_unshocked*(1.0 - f_sm) + p_air*f_sm;
                                    
                                    const double rho_i = rho_Y_0_i + rho_Y_1_i;
                                    const double Y_0_i = rho_Y_0_i/rho_i;
                                    const double Y_1_i = 1.0 - Y_0_i;
                                    
                                    const double gamma_i = (Y_0_i*c_p_SF6 + Y_1_i*c_p_air)/
                                        (Y_0_i*c_v_SF6 + Y_1_i*c_v_air);
                                    
                                    rho_Y_0[idx_cell] = rho_Y_0_i;
                                    rho_Y_1[idx_cell] = rho_Y_1_i;
                                    rho_u[idx_cell]   = rho_i*u_i;
                                    rho_v[idx_cell]   = rho_i*v_i;
                                    E[idx_cell]       = p_i/(gamma_i - 1.0) + 0.5*rho_i*(u_i*u_i + v_i*v_i);
                                }
                            }
                        }
                    }
                    else if (d_project_name == "2D Poggi's RMI 1")
                    {
                        if (d_num_species != 2)
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "Please provide only two-species for the problem"
                                << " '2D Poggi's RMI'."
                                << std::endl);
                        }
                        
                        // Characteristic length of the initial interface thickness.
                        const double epsilon_i = 0.001;
                        
                        double* rho_Y_0   = partial_density->getPointer(0);
                        double* rho_Y_1   = partial_density->getPointer(1);
                        double* rho_u     = momentum->getPointer(0);
                        double* rho_v     = momentum->getPointer(1);
                        double* E         = total_energy->getPointer(0);
                        
                        // species 0: SF6.
                        // species 1: air.
                        const double gamma_0 = 1.09312;
                        const double gamma_1 = 1.39909;
                        
                        const double c_p_SF6 = 668.286;
                        const double c_p_air = 1040.50;
                        
                        const double c_v_SF6 = 611.359;
                        const double c_v_air = 743.697;
                        
                        NULL_USE(gamma_1);
                        
                        // Unshocked SF6.
                        const double rho_unshocked = 5.97286552525647;
                        const double u_unshocked   = 0.0;
                        const double v_unshocked   = 0.0;
                        const double p_unshocked   = 101325.0;
                        
                        // Shocked SF6.
                        const double rho_shocked = 11.9708247309869;
                        const double u_shocked   = 98.9344103891513;
                        const double v_shocked   = 0.0;
                        const double p_shocked   = 218005.430874;
                        
                        // Air.
                        const double rho_air = 1.14560096494103;
                        const double u_air   = 0.0;
                        const double v_air   = 0.0;
                        const double p_air   = 101325.0;
                        
                        // Shock hits the interface after 0.1 ms.
                        const double L_x_shock = 0.180254509478037;
                        const double L_x_interface = 0.2;
                        
                        // Perturbations due to S mode.
                        const double A = 1.0e-3; // Amplitude.
                        
                        for (int j = 0; j < patch_dims[1]; j++)
                        {
                            for (int i = 0; i < patch_dims[0]; i++)
                            {
                                // Compute the linear index.
                                int idx_cell = i + j*patch_dims[0];
                                
                                // Compute the coordinates.
                                double x[2];
                                x[0] = patch_xlo[0] + (i + 0.5)*dx[0];
                                x[1] = patch_xlo[1] + (j + 0.5)*dx[1];
                                
                                double S = 0.0;
                                for (int m = 40; m <= 60; m++)
                                {
                                    S += A*cos(2.0*M_PI*m/0.05*x[1] + tan(7.0*m));
                                }
                                
                                if (x[0] < L_x_shock)
                                {
                                    rho_Y_0[idx_cell] = rho_shocked;
                                    rho_Y_1[idx_cell] = 0.0;
                                    rho_u[idx_cell] = rho_shocked*u_shocked;
                                    rho_v[idx_cell] = rho_shocked*v_shocked;
                                    E[idx_cell]     = p_shocked/(gamma_0 - 1.0) +
                                        0.5*rho_shocked*(u_shocked*u_shocked + v_shocked*v_shocked);
                                }
                                else
                                {
                                    const double f_sm = 0.5*(1.0 + erf((x[0] - (L_x_interface + S))/epsilon_i));
                                    
                                    // Smooth the primitive variables.
                                    const double rho_Y_0_i = rho_unshocked*(1.0 - f_sm);
                                    const double rho_Y_1_i = rho_air*f_sm;
                                    const double u_i = u_unshocked*(1.0 - f_sm) + u_air*f_sm;
                                    const double v_i = v_unshocked*(1.0 - f_sm) + v_air*f_sm;
                                    const double p_i = p_unshocked*(1.0 - f_sm) + p_air*f_sm;
                                    
                                    const double rho_i = rho_Y_0_i + rho_Y_1_i;
                                    const double Y_0_i = rho_Y_0_i/rho_i;
                                    const double Y_1_i = 1.0 - Y_0_i;
                                    
                                    const double gamma_i = (Y_0_i*c_p_SF6 + Y_1_i*c_p_air)/
                                        (Y_0_i*c_v_SF6 + Y_1_i*c_v_air);
                                    
                                    rho_Y_0[idx_cell] = rho_Y_0_i;
                                    rho_Y_1[idx_cell] = rho_Y_1_i;
                                    rho_u[idx_cell]   = rho_i*u_i;
                                    rho_v[idx_cell]   = rho_i*v_i;
                                    E[idx_cell]       = p_i/(gamma_i - 1.0) + 0.5*rho_i*(u_i*u_i + v_i*v_i);
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
                    if (false)
                    {}
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
            case FLOW_MODEL::FIVE_EQN_ALLAIRE:
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
                    if (d_project_name == "2D shock-bubble interaction with smoothed interface")
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
                        const double gamma_0 = 1.648;
                        const double gamma_1 = 1.4;
                        
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
                                    
                                    const double gamma = 1.0/(Z_1_i/(gamma_0 - 1.0) + Z_2_i/(gamma_1 - 1.0)) + 1.0;
                                    
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
                    else if (d_project_name == "2D shock-bubble interaction")
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
                        const double gamma_0 = 1.648;
                        const double gamma_1 = 1.4;
                        
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
                                    
                                    const double gamma = 1.0/(Z_1_i/(gamma_0 - 1.0) + Z_2_i/(gamma_1 - 1.0)) + 1.0;
                                    
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
                    else if (d_project_name == "2D Richtmyer-Meshkov instability with smoothed interface")
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
                        const double gamma_0 = 1.093;
                        const double gamma_1 = 1.4;
                        
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
                                    
                                    const double gamma = 1.0/(Z_1_i/(gamma_0 - 1.0) + Z_2_i/(gamma_1 - 1.0)) + 1.0;
                                    
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
                        double gamma_m  = 1.6;
                        double rho_m    = 10.0;
                        double u_m      = 0.5;
                        double v_m      = 0.5;
                        double p_m      = 1.0/1.4;
                        
                        // ambient initial conditions.
                        double gamma_a  = 1.4;
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
                        double gamma_m  = 1.6;
                        double rho_m    = 10.0;
                        double u_m      = 0.5;
                        double v_m      = 0.5;
                        double w_m      = 0.5;
                        double p_m      = 1.0/1.4;
                        
                        // ambient initial conditions.
                        double gamma_a  = 1.4;
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
