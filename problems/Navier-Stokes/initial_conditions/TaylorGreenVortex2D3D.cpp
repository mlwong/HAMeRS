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
    NULL_USE(data_time);
    
    if ((d_project_name != "3D Taylor-Green vortex") &&
        (d_project_name != "2D Taylor-Green vortex"))

    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Can only initialize data for 'project_name' = '3D Taylor-Green vortex' or "
            << "Can only initialize data for 'project_name' = '2D Taylor-Green vortex' !\n"
            << "'project_name' = '"
            << d_project_name
            << "' is given."
            << std::endl);
    }
    
    if (d_dim != tbox::Dimension(3) && d_dim != tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Dimension of problem should be 2 or 3!"
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
        const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
            HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                patch.getPatchGeometry()));
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(patch_geom);
#endif
        
        const double* const dx = patch_geom->getDx();
        const double* const patch_xlo = patch_geom->getXLower();
        

        if (d_project_name == "3D Taylor-Green vortex")
        {        
            // Get the dimensions of box that covers the interior of Patch.
            hier::Box patch_box = patch.getBox();
            const hier::IntVector patch_dims = patch_box.numberCells();
        
            /*
             * Initialize data for a 3D inviscid Taylor-Green problem.
             */
        

            HAMERS_SHARED_PTR<pdat::CellData<double> > density      = conservative_variables[0];
            HAMERS_SHARED_PTR<pdat::CellData<double> > momentum     = conservative_variables[1];
            HAMERS_SHARED_PTR<pdat::CellData<double> > total_energy = conservative_variables[2];
        
            double* rho   = density->getPointer(0);
            double* rho_u = momentum->getPointer(0);
            double* rho_v = momentum->getPointer(1);
            double* rho_w = momentum->getPointer(2);
            double* E     = total_energy->getPointer(0);
        
            double gamma = double(5)/double(3);
        
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
                        
                        // Compute density, velocities and pressure.
                        const double rho_i = double(1);
                        const double u_i = sin(x[0])*cos(x[1])*cos(x[2]);
                        const double v_i = -cos(x[0])*sin(x[1])*cos(x[2]);
                        const double w_i = double(0);
                        const double p_i = double(100) + ((cos(double(2)*x[2]) + double(2))*
                            (cos(double(2)*x[0]) + cos(double(2)*x[1])) - double(2))/double(16);
                        
                        rho[idx_cell]   = rho_i;
                        rho_u[idx_cell] = rho_i*u_i;
                        rho_v[idx_cell] = rho_i*v_i;
                        rho_w[idx_cell] = rho_i*w_i;
                        E[idx_cell]     = p_i/(gamma - double(1)) + double(1)/double(2)*rho_i*
                            (u_i*u_i + v_i*v_i + w_i*w_i);
                    }
                }
            }
        }
        else if (d_project_name == "2D Taylor-Green vortex")
        {
            // Get the dimensions of box that covers the interior of Patch.
            hier::Box patch_box = patch.getBox();
            const hier::IntVector patch_dims = patch_box.numberCells();
        
            /*
             * Initialize data for a 2D inviscid Taylor-Green problem.
             */
        

            HAMERS_SHARED_PTR<pdat::CellData<double> > density      = conservative_variables[0];
            HAMERS_SHARED_PTR<pdat::CellData<double> > momentum     = conservative_variables[1];
            HAMERS_SHARED_PTR<pdat::CellData<double> > total_energy = conservative_variables[2];
        
            double* rho   = density->getPointer(0);
            double* rho_u = momentum->getPointer(0);
            double* rho_v = momentum->getPointer(1);
            double* E     = total_energy->getPointer(0);
        
            double gamma = double(5)/double(3);

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
                    
                    // Compute density, velocities and pressure.
                    const double rho_i = double(1);
                    const double u_i = sin(x[0])*cos(x[1]);
                    const double v_i = -cos(x[0])*sin(x[1]);
                    const double p_i = double(100) + ((cos(double(2)*x[0]) + cos(double(2)*x[1])) 
                                        - double(2))/double(16);
                    rho[idx_cell]   = rho_i;
                    rho_u[idx_cell] = rho_i*u_i;
                    rho_v[idx_cell] = rho_i*v_i;
                    E[idx_cell]     = p_i/(gamma - double(1)) + double(1)/double(2)*rho_i*
                        (u_i*u_i + v_i*v_i);
                }
            }
        }
    }
}
