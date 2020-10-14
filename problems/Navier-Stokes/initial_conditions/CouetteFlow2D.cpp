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
    
    if ((d_project_name != "2D Couette flow in x-direction") &&
        (d_project_name != "2D Couette flow in y-direction"))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Can only initialize data for 'project_name' = "
            << "'2D Couette flow in x-direction' or "
            << "'2D Couette flow in y-direction'"
            << "!\n"
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
        
        // Get the dimensions of box that covers the interior of Patch.
        hier::Box patch_box = patch.getBox();
        const hier::IntVector patch_dims = patch_box.numberCells();
        
        /*
         * Initialize data for 2D Couette flow problem.
         */
            
        HAMERS_SHARED_PTR<pdat::CellData<double> > density      = conservative_variables[0];
        HAMERS_SHARED_PTR<pdat::CellData<double> > momentum     = conservative_variables[1];
        HAMERS_SHARED_PTR<pdat::CellData<double> > total_energy = conservative_variables[2];
        
        double* rho   = density->getPointer(0);
        double* rho_u = momentum->getPointer(0);
        double* rho_v = momentum->getPointer(1);
        double* E     = total_energy->getPointer(0);
        
        if (d_project_name == "2D Couette flow in x-direction")
        {
            const double gamma   = double(7)/double(5);
            const double rho_avg = double(1.159750086791891);
            const double u_w     = double(69.445);
            const double p_c     = double(1.0e5);
            
            const double H = double(1);
            
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
                    
                    double phi = x[1]/H;
                    double u   = u_w*phi;
                    
                    rho[idx_cell]   = rho_avg;
                    rho_u[idx_cell] = rho[idx_cell]*u;
                    rho_v[idx_cell] = double(0);
                    E[idx_cell]     = p_c/(gamma - double(1)) + double(1)/double(2)*rho_avg*u*u;
                }
            }
        }
        else if (d_project_name == "2D Couette flow in y-direction")
        {
            const double gamma   = double(7)/double(5);
            const double rho_avg = double(1.159750086791891);
            const double v_w     = double(69.445);
            const double p_c     = double(1.0e5);
            
            const double H = double(1);
            
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
                    
                    double phi = x[0]/H;
                    double v   = v_w*phi;
                    
                    rho[idx_cell]   = rho_avg;
                    rho_u[idx_cell] = double(0);
                    rho_v[idx_cell] = rho[idx_cell]*v;
                    E[idx_cell]     = p_c/(gamma - double(1)) + double(1)/double(2)*rho_avg*v*v;
                }
            }
        }
    }
}
