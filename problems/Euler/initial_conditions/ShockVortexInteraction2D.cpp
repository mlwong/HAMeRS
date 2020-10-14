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
    
    if (d_project_name != "2D shock-vortex interaction")
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Can only initialize data for 'project_name' = '2D shock-vortex interaction'!\n"
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
         * Initialize data for a 2D shock-vortex interaction problem.
         */
        
        HAMERS_SHARED_PTR<pdat::CellData<double> > density      = conservative_variables[0];
        HAMERS_SHARED_PTR<pdat::CellData<double> > momentum     = conservative_variables[1];
        HAMERS_SHARED_PTR<pdat::CellData<double> > total_energy = conservative_variables[2];
        
        double* rho   = density->getPointer(0);
        double* rho_u = momentum->getPointer(0);
        double* rho_v = momentum->getPointer(1);
        double* E     = total_energy->getPointer(0);
        
        const double gamma = double(7)/double(5);
        
        // Vortex strength.
        const double M_v = double(1);
        
        // Vortex radius.
        const double R = double(1);
        
        // Vortex center.
        double x_v[2];
        x_v[0] = double(4);
        x_v[1] = double(0);
        
        // Post-shock condition.
        const double rho_post = double(1.34161490);
        const double p_post   = double(1.51333333)/gamma;
        const double u_post   = double(-0.89444445);
        const double v_post   = double(0);
        
        // Pre-shock condition.
        const double rho_pre = double(1);
        const double p_pre   = double(1)/gamma;
        const double u_pre   = double(-6)/double(5);
        const double v_pre   = double(0);
        
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
                
                if (x[0] < 0)
                {
                    rho[idx_cell]     = rho_post;
                    rho_u[idx_cell]   = rho_post*u_post;
                    rho_v[idx_cell]   = rho_post*v_post;
                    E[idx_cell]       = p_post/(gamma - double(1)) + double(1)/double(2)*rho_post*
                        (u_post*u_post + v_post*v_post);
                }
                else
                {
                    double r = sqrt(pow(x[0] - x_v[0], 2) + pow(x[1] - x_v[1], 2));
                    
                    if (r > double(4))
                    {
                        rho[idx_cell]   = rho_pre;
                        rho_u[idx_cell] = rho_pre*u_pre;
                        rho_v[idx_cell] = rho_pre*v_pre;
                        E[idx_cell]     = p_pre/(gamma - double(1)) + double(1)/double(2)*rho_pre*
                            (u_pre*u_pre + v_pre*v_pre);
                    }
                    else
                    {
                        double p = double(1)/(gamma)*
                            pow((double(1) - double(1)/double(2)*(gamma - double(1))*M_v*M_v*exp(1 - pow(r/R, 2))),
                                gamma/(gamma - double(1)));
                        double u = u_pre - M_v*exp(double(1)/double(2)*(double(1) - pow(r/R, 2)))*(x[1] - x_v[1]);
                        double v = v_pre + M_v*exp(double(1)/double(2)*(double(1) - pow(r/R, 2)))*(x[0] - x_v[0]);
                        
                        rho[idx_cell] =
                            pow((double(1) - double(1)/double(2)*(gamma - double(1))*M_v*M_v*exp(1 - pow(r/R, 2))),
                                double(1)/(gamma - double(1)));
                        rho_u[idx_cell] = rho[idx_cell]*u;
                        rho_v[idx_cell] = rho[idx_cell]*v;
                        E[idx_cell] = p/(gamma - double(1)) + double(1)/double(2)*rho[idx_cell]*
                            (u*u + v*v);
                    }
                }
            }
        }
    }
}
