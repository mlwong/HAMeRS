#include "apps/Euler/EulerInitialConditions.hpp"

/*
 * Set the data on the patch interior to some initial values.
 */
void
EulerInitialConditions::initializeDataOnPatch(
    hier::Patch& patch,
    const std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& conservative_variables,
    const double data_time,
    const bool initial_time)
{
    NULL_USE(data_time);
    
    if ((d_project_name != "2D convergence test single-species") &&
        (d_project_name != "3D convergence test single-species"))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Can only initialize data for 'project_name' = '2D convergence test single-species' or "
            << "'3D convergence test single-species'!\n"
            << "'project_name' = '"
            << d_project_name
            << "' is given."
            << std::endl);
    }
    
    
    if (d_flow_model_type != FLOW_MODEL::SINGLE_SPECIES)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Flow model should be single-species!"
            << std::endl);
    }
    
    if (d_project_name == "2D convergence test single-species")
    {
        if (d_dim != tbox::Dimension(2))
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Dimension of problem should be 2!"
                << std::endl);
        }
        
        if (initial_time)
        {
            const double* const domain_xlo = d_grid_geometry->getXLower();
            const double* const domain_xhi = d_grid_geometry->getXUpper();
            
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
             * Initialize data for a 2D density wave advection problem.
             */
            
            HAMERS_SHARED_PTR<pdat::CellData<double> > density      = conservative_variables[0];
            HAMERS_SHARED_PTR<pdat::CellData<double> > momentum     = conservative_variables[1];
            HAMERS_SHARED_PTR<pdat::CellData<double> > total_energy = conservative_variables[2];
            
            double* rho   = density->getPointer(0);
            double* rho_u = momentum->getPointer(0);
            double* rho_v = momentum->getPointer(1);
            double* E     = total_energy->getPointer(0);
            
            double gamma = double(7)/double(5);
            
            double u = double(1);
            double v = double(1);
            double p = double(1);
            
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
                    
                    rho[idx_cell]   = double(1) + double(1)/double(2)*sin(M_PI*(x[0] + x[1]));
                    rho_u[idx_cell] = rho[idx_cell]*u;
                    rho_v[idx_cell] = rho[idx_cell]*v;
                    E[idx_cell]     = p/(gamma - double(1)) + double(1)/double(2)*rho[idx_cell]*
                        (u*u + v*v);
                }
            }
        }
    }
    else if (d_project_name == "3D convergence test single-species")
    {
        if (d_dim != tbox::Dimension(3))
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Dimension of problem should be 3!"
                << std::endl);
        }
        
        if (initial_time)
        {
            const double* const domain_xlo = d_grid_geometry->getXLower();
            const double* const domain_xhi = d_grid_geometry->getXUpper();
            
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
             * Initialize data for a 2D density wave advection problem.
             */
            
            HAMERS_SHARED_PTR<pdat::CellData<double> > density      = conservative_variables[0];
            HAMERS_SHARED_PTR<pdat::CellData<double> > momentum     = conservative_variables[1];
            HAMERS_SHARED_PTR<pdat::CellData<double> > total_energy = conservative_variables[2];
            
            double* rho   = density->getPointer(0);
            double* rho_u = momentum->getPointer(0);
            double* rho_v = momentum->getPointer(1);
            double* rho_w = momentum->getPointer(2);
            double* E     = total_energy->getPointer(0);
            
            double gamma = double(7)/double(5);
            
            double u = double(1);
            double v = double(1);
            double w = double(1);
            double p = double(1);
            
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
                        
                        rho[idx_cell]   = double(1) + double(1)/double(2)*sin(M_PI*(x[0] + x[1] + x[2]));
                        rho_u[idx_cell] = rho[idx_cell]*u;
                        rho_v[idx_cell] = rho[idx_cell]*v;
                        rho_w[idx_cell] = rho[idx_cell]*w;
                        E[idx_cell]     = p/(gamma - double(1)) + double(1)/double(2)*rho[idx_cell]*
                            (u*u + v*v + w*w);
                    }
                }
            }
        }
    }
}
