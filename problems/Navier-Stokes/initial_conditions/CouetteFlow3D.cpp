#include "apps/Navier-Stokes/NavierStokesInitialConditions.hpp"

/*
 * Set the data on the patch interior to some initial values.
 */
void
NavierStokesInitialConditions::initializeDataOnPatch(
    hier::Patch& patch,
    const std::vector<boost::shared_ptr<pdat::CellData<double> > >& conservative_variables,
    const double data_time,
    const bool initial_time)
{
    NULL_USE(data_time);
    
    if ((d_project_name != "3D Couette flow in x-direction") &&
        (d_project_name != "3D Couette flow in y-direction") &&
        (d_project_name != "3D Couette flow in z-direction"))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Can only initialize data for 'project_name' = "
            << "'3D Couette flow in x-direction' or "
            << "'3D Couette flow in y-direction' or "
            << "'3D Couette flow in z-direction'"
            << "!\n"
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
    
    if (d_flow_model_type != FLOW_MODEL::SINGLE_SPECIES)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Flow model should be single-species!"
            << std::endl);
    }
    
    if (initial_time)
    {
        const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
            BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
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
         * Initialize data for 3D Couette flow problem.
         */
        
        boost::shared_ptr<pdat::CellData<double> > density      = conservative_variables[0];
        boost::shared_ptr<pdat::CellData<double> > momentum     = conservative_variables[1];
        boost::shared_ptr<pdat::CellData<double> > total_energy = conservative_variables[2];
        
        double* rho   = density->getPointer(0);
        double* rho_u = momentum->getPointer(0);
        double* rho_v = momentum->getPointer(1);
        double* rho_w = momentum->getPointer(2);
        double* E     = total_energy->getPointer(0);
        
        if (d_project_name == "3D Couette flow in x-direction")
        {
            const double gamma   = double(7)/double(5);
            const double rho_avg = double(1.159750086791891);
            const double u_w     = double(69.445);
            const double p_c     = double(1.0e5);
            
            const double H = double(1);
            
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
                        x[0] = patch_xlo[0] + (double(i) + double(1)/double(2))*dx[0];
                        x[1] = patch_xlo[1] + (double(j) + double(1)/double(2))*dx[1];
                        x[2] = patch_xlo[2] + (double(k) + double(1)/double(2))*dx[2];
                        
                        double phi = x[1]/H;
                        double u   = u_w*phi;
                        
                        rho[idx_cell]   = rho_avg;
                        rho_u[idx_cell] = rho[idx_cell]*u;
                        rho_v[idx_cell] = double(0);
                        rho_w[idx_cell] = double(0);
                        E[idx_cell]     = p_c/(gamma - double(1)) + double(1)/double(2)*rho_avg*u*u;
                    }
                }
            }
        }
        else if (d_project_name == "3D Couette flow in y-direction")
        {
            const double gamma   = double(7)/double(5);
            const double rho_avg = double(1.159750086791891);
            const double v_w     = double(69.445);
            const double p_c     = double(1.0e5);
            
            const double H = double(1);
            
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
                        x[0] = patch_xlo[0] + (double(i) + double(1)/double(2))*dx[0];
                        x[1] = patch_xlo[1] + (double(j) + double(1)/double(2))*dx[1];
                        x[2] = patch_xlo[2] + (double(k) + double(1)/double(2))*dx[2];
                        
                        double phi = x[2]/H;
                        double v   = v_w*phi;
                        
                        rho[idx_cell]   = rho_avg;
                        rho_u[idx_cell] = double(0);
                        rho_v[idx_cell] = rho[idx_cell]*v;
                        rho_w[idx_cell] = double(0);
                        E[idx_cell]     = p_c/(gamma - double(1)) + double(1)/double(2)*rho_avg*v*v;
                    }
                }
            }
        }
        else if (d_project_name == "3D Couette flow in z-direction")
        {
            const double gamma   = double(7)/double(5);
            const double rho_avg = double(1.159750086791891);
            const double w_w     = double(69.445);
            const double p_c     = double(1.0e5);
            
            const double H = double(1);
            
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
                        x[0] = patch_xlo[0] + (double(i) + double(1)/double(2))*dx[0];
                        x[1] = patch_xlo[1] + (double(j) + double(1)/double(2))*dx[1];
                        x[2] = patch_xlo[2] + (double(k) + double(1)/double(2))*dx[2];
                        
                        double phi = x[0]/H;
                        double w   = w_w*phi;
                        
                        rho[idx_cell]   = rho_avg;
                        rho_u[idx_cell] = double(0);
                        rho_v[idx_cell] = double(0);
                        rho_w[idx_cell] = rho[idx_cell]*w;
                        E[idx_cell]     = p_c/(gamma - double(1)) + double(1)/double(2)*rho_avg*w*w;
                    }
                }
            }
        }
    }
}
