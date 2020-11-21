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
    
    if ((d_project_name != "3D vortex leaving domain in x-direction") &&
        (d_project_name != "3D vortex leaving domain in y-direction") &&
        (d_project_name != "3D vortex leaving domain in z-direction"))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Can only initialize data for 'project_name' = "
            << "'3D vortex leaving domain in x-direction', "
            << "'3D vortex leaving domain in y-direction' or "
            << "'3D vortex leaving domain in z-direction'"
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
         * Initialize data for 3D vortex leaving domain problem.
         */
            
        HAMERS_SHARED_PTR<pdat::CellData<double> > density      = conservative_variables[0];
        HAMERS_SHARED_PTR<pdat::CellData<double> > momentum     = conservative_variables[1];
        HAMERS_SHARED_PTR<pdat::CellData<double> > total_energy = conservative_variables[2];
        
        double* rho   = density->getPointer(0);
        double* rho_u = momentum->getPointer(0);
        double* rho_v = momentum->getPointer(1);
        double* rho_w = momentum->getPointer(2);
        double* E     = total_energy->getPointer(0);
        
        if (d_project_name == "3D vortex leaving domain in x-direction")
        {
            const double gamma = double(7)/double(5);
            
            const double L       = double(1.0);
            const double Gamma_v = double(0.024);
            const double R_v     = L/double(10);
            const double Ma      = double(0.283);
            
            const double x_v = double(0);
            const double y_v = double(0);
            
            const double rho_inf = double(1);
            const double p_inf   = double(1.0)/gamma;
            const double u_inf   = Ma*sqrt(gamma*p_inf/rho_inf);
            
            const double c = sqrt(gamma*p_inf/rho_inf);
            const double Gamma_normalized = Gamma_v/(c*R_v);
            
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
                        
                        const double r = sqrt(pow(x[0] - x_v, 2) + pow(x[1] - y_v, 2));
                        
                        const double exp_factor = exp(-pow(r/R_v, 2));
                        const double exp_factor_half  = sqrt(exp_factor);
                        
                        const double p_vortex   = p_inf*exp(-double(1)/double(2)*gamma*Gamma_normalized*Gamma_normalized*exp_factor);
                        const double rho_vortex = rho_inf/p_inf*p_vortex;
                        
                        const double u_vortex = u_inf - exp_factor_half * (x[1] - y_v) * Gamma_v/pow(R_v, 2);
                        const double v_vortex =         exp_factor_half * (x[0] - x_v) * Gamma_v/pow(R_v, 2);
                        rho[idx_cell]   = rho_vortex;
                        rho_u[idx_cell] = rho_vortex*u_vortex;
                        rho_v[idx_cell] = rho_vortex*v_vortex;
                        rho_w[idx_cell] = double(0);
                        E[idx_cell]     = p_vortex/(gamma - double(1)) + double(1)/double(2)*rho_vortex*(u_vortex*u_vortex + v_vortex*v_vortex);
                    }
                }
            }
        }
        else if (d_project_name == "3D vortex leaving domain in y-direction")
        {
            const double gamma = double(7)/double(5);
            
            const double L       = double(1.0);
            const double Gamma_v = double(0.024);
            const double R_v     = L/double(10);
            const double Ma      = -double(0.283);
            
            const double x_v = double(0);
            const double y_v = double(0);
            
            const double rho_inf = double(1);
            const double p_inf   = double(1.0)/gamma;
            const double v_inf   = Ma*sqrt(gamma*p_inf/rho_inf);
            
            const double c = sqrt(gamma*p_inf/rho_inf);
            const double Gamma_normalized = Gamma_v/(c*R_v);
            
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
                        
                        const double r = sqrt(pow(x[0] - x_v, 2) + pow(x[1] - y_v, 2));
                        
                        const double exp_factor = exp(-pow(r/R_v, 2));
                        const double exp_factor_half  = sqrt(exp_factor);
                        
                        const double p_vortex   = p_inf*exp(-double(1)/double(2)*gamma*Gamma_normalized*Gamma_normalized*exp_factor);
                        const double rho_vortex = rho_inf/p_inf*p_vortex;
                        
                        const double v_vortex = v_inf - exp_factor_half * (x[0] - x_v) * Gamma_v/pow(R_v, 2);
                        const double u_vortex =         exp_factor_half * (x[1] - y_v) * Gamma_v/pow(R_v, 2);
                        rho[idx_cell]   = rho_vortex;
                        rho_u[idx_cell] = rho_vortex*u_vortex;
                        rho_v[idx_cell] = rho_vortex*v_vortex;
                        rho_w[idx_cell] = double(0);
                        E[idx_cell]     = p_vortex/(gamma - double(1)) + double(1)/double(2)*rho_vortex*(u_vortex*u_vortex + v_vortex*v_vortex);
                    }
                }
            }
        }
        else if (d_project_name == "3D vortex leaving domain in z-direction")
        {
            const double gamma = double(7)/double(5);
            
            const double L       = double(1.0);
            const double Gamma_v = double(0.024);
            const double R_v     = L/double(10);
            const double Ma      = double(0.283);
            
            const double x_v = double(0);
            const double z_v = double(0);
            
            const double rho_inf = double(1);
            const double p_inf   = double(1.0)/gamma;
            const double w_inf   = Ma*sqrt(gamma*p_inf/rho_inf);
            
            const double c = sqrt(gamma*p_inf/rho_inf);
            const double Gamma_normalized = Gamma_v/(c*R_v);
            
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
                        
                        const double r = sqrt(pow(x[0] - x_v, 2) + pow(x[2] - z_v, 2));
                        
                        const double exp_factor = exp(-pow(r/R_v, 2));
                        const double exp_factor_half  = sqrt(exp_factor);
                        
                        const double p_vortex   = p_inf*exp(-double(1)/double(2)*gamma*Gamma_normalized*Gamma_normalized*exp_factor);
                        const double rho_vortex = rho_inf/p_inf*p_vortex;
                        
                        const double w_vortex = w_inf - exp_factor_half * (x[0] - x_v) * Gamma_v/pow(R_v, 2);
                        const double u_vortex =         exp_factor_half * (x[2] - z_v) * Gamma_v/pow(R_v, 2);
                        rho[idx_cell]   = rho_vortex;
                        rho_u[idx_cell] = rho_vortex*u_vortex;
                        rho_v[idx_cell] = double(0);
                        rho_w[idx_cell] = rho_vortex*w_vortex;
                        E[idx_cell]     = p_vortex/(gamma - double(1)) + double(1)/double(2)*rho_vortex*(u_vortex*u_vortex + w_vortex*w_vortex);
                    }
                }
            }
        }
    }    
}
