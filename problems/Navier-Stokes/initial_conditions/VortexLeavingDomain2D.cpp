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
    
    if ((d_project_name != "2D vortex leaving domain in x-direction") &&
        (d_project_name != "2D vortex leaving domain in y-direction"))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Can only initialize data for 'project_name' = "
            << "'2D vortex leaving domain in x-direction' or "
            << "'2D vortex leaving domain in y-direction'"
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
    
    if (d_flow_model_type == FLOW_MODEL::SINGLE_SPECIES)
    {
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
             * Initialize data for 2D vortex leaving domain problem.
             */
                
            HAMERS_SHARED_PTR<pdat::CellData<double> > density      = conservative_variables[0];
            HAMERS_SHARED_PTR<pdat::CellData<double> > momentum     = conservative_variables[1];
            HAMERS_SHARED_PTR<pdat::CellData<double> > total_energy = conservative_variables[2];
            
            double* rho   = density->getPointer(0);
            double* rho_u = momentum->getPointer(0);
            double* rho_v = momentum->getPointer(1);
            double* E     = total_energy->getPointer(0);
            
            if (d_project_name == "2D vortex leaving domain in x-direction")
            {
                const double gamma = double(7)/double(5);
                
                const double L       = double(1);
                const double Gamma_v = double(0.024);
                const double R_v     = L/double(10);
                const double Ma      = -double(0.283);
                
                const double x_v = double(0);
                const double y_v = double(0);
                
                const double rho_inf = double(1);
                const double p_inf   = double(1)/gamma;
                const double u_inf   = Ma*sqrt(gamma*p_inf/rho_inf);
                
                const double c = sqrt(gamma*p_inf/rho_inf);
                const double Gamma_normalized = Gamma_v/(c*R_v);
                
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
                        E[idx_cell]     = p_vortex/(gamma - double(1)) + double(1)/double(2)*rho_vortex*(u_vortex*u_vortex + v_vortex*v_vortex);
                    }
                }
            }
            else if (d_project_name == "2D vortex leaving domain in y-direction")
            {
                const double gamma = double(7)/double(5);
                
                const double L       = double(1);
                const double Gamma_v = double(0.024);
                const double R_v     = L/double(10);
                const double Ma      = -double(0.283);
                
                const double x_v = double(0);
                const double y_v = double(0);
                
                const double rho_inf = double(1);
                const double p_inf   = double(1)/gamma;
                const double v_inf   = Ma*sqrt(gamma*p_inf/rho_inf);
                
                const double c = sqrt(gamma*p_inf/rho_inf);
                const double Gamma_normalized = Gamma_v/(c*R_v);
                
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
                        E[idx_cell]     = p_vortex/(gamma - double(1)) + double(1)/double(2)*rho_vortex*(u_vortex*u_vortex + v_vortex*v_vortex);
                    }
                }
            }
        }
    }
    else if (d_flow_model_type == FLOW_MODEL::FOUR_EQN_CONSERVATIVE)
    {
        if (d_flow_model->getNumberOfSpecies() != 2)
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Number of species should be 2!"
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
             * Initialize data for 2D vortex leaving domain problem.
             */
                
            HAMERS_SHARED_PTR<pdat::CellData<double> > partial_densities = conservative_variables[0];
            HAMERS_SHARED_PTR<pdat::CellData<double> > momentum          = conservative_variables[1];
            HAMERS_SHARED_PTR<pdat::CellData<double> > total_energy      = conservative_variables[2];
            
            double* rho_Y_0 = partial_densities->getPointer(0);
            double* rho_Y_1 = partial_densities->getPointer(1);
            double* rho_u   = momentum->getPointer(0);
            double* rho_v   = momentum->getPointer(1);
            double* E       = total_energy->getPointer(0);
            
            if (d_project_name == "2D vortex leaving domain in x-direction")
            {
                const double c_p_0 = double(7);
                const double c_p_1 = double(7);
                
                const double c_v_0 = double(5);
                const double c_v_1 = double(5);
                
                const double L       = double(1);
                const double Gamma_v = double(0.024);
                const double R_v     = L/double(10);
                const double Ma      = double(0.283);
                
                const double x_v = double(0);
                const double y_v = double(0);
                
                const double rho_inf = double(1);
                const double Y_0_inf = double(1)/double(2);
                const double Y_1_inf = double(1)/double(2);
                
                const double c_p_mixture = Y_0_inf*c_p_0 + Y_1_inf*c_p_1;
                const double c_v_mixture = Y_0_inf*c_v_0 + Y_1_inf*c_v_1;
                
                const double gamma_mixture = c_p_mixture/c_v_mixture;
                
                const double p_inf = double(1)/gamma_mixture;
                const double u_inf = Ma*sqrt(gamma_mixture*p_inf/rho_inf);
                
                const double c = sqrt(gamma_mixture*p_inf/rho_inf);
                const double Gamma_normalized = Gamma_v/(c*R_v);
                
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
                        
                        const double r = sqrt(pow(x[0] - x_v, 2) + pow(x[1] - y_v, 2));
                        
                        const double exp_factor = exp(-pow(r/R_v, 2));
                        const double exp_factor_half = sqrt(exp_factor);
                        
                        const double p_vortex   = p_inf*exp(-double(1)/double(2)*gamma_mixture*Gamma_normalized*Gamma_normalized*exp_factor);
                        const double rho_vortex = rho_inf/p_inf*p_vortex;
                        
                        const double u_vortex = u_inf - exp_factor_half * (x[1] - y_v) * Gamma_v/pow(R_v, 2);
                        const double v_vortex =         exp_factor_half * (x[0] - x_v) * Gamma_v/pow(R_v, 2);
                        
                        rho_Y_0[idx_cell] = rho_vortex*Y_0_inf;
                        rho_Y_1[idx_cell] = rho_vortex*Y_1_inf;
                        rho_u[idx_cell]   = rho_vortex*u_vortex;
                        rho_v[idx_cell]   = rho_vortex*v_vortex;
                        E[idx_cell]       = p_vortex/(gamma_mixture - double(1)) + double(1)/double(2)*rho_vortex*(u_vortex*u_vortex + v_vortex*v_vortex);
                    }
                }
            }
            else if (d_project_name == "2D vortex leaving domain in y-direction")
            {
                const double c_p_0 = double(7);
                const double c_p_1 = double(7);
                
                const double c_v_0 = double(5);
                const double c_v_1 = double(5);
                
                const double L       = double(1);
                const double Gamma_v = double(0.024);
                
                const double R_v     = L/double(10);
                const double Ma      = -double(0.283);
                
                const double x_v = double(0);
                const double y_v = double(0);
                
                const double rho_inf = double(1);
                const double Y_0_inf = double(1)/double(2);
                const double Y_1_inf = double(1)/double(2);
                const double c_p_mixture = Y_0_inf*c_p_0 + Y_1_inf*c_p_1;
                const double c_v_mixture = Y_0_inf*c_v_0 + Y_1_inf*c_v_1;
                const double gamma_mixture = c_p_mixture/c_v_mixture;
                const double p_inf = double(1)/gamma_mixture;
                const double v_inf = Ma*sqrt(gamma_mixture*p_inf/rho_inf);
                const double c = sqrt(gamma_mixture*p_inf/rho_inf);
                const double Gamma_normalized = Gamma_v/(c*R_v);
                
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
                        
                        const double r = sqrt(pow(x[0] - x_v, 2) + pow(x[1] - y_v, 2));
                        
                        const double exp_factor = exp(-pow(r/R_v, 2));
                        const double exp_factor_half  = sqrt(exp_factor);
                         
                        const double p_vortex   = p_inf*exp(-double(1)/double(2)*gamma_mixture*Gamma_normalized*Gamma_normalized*exp_factor);
                        const double rho_vortex = rho_inf/p_inf*p_vortex;
                        
                        const double v_vortex = v_inf - exp_factor_half * (x[0] - x_v) * Gamma_v/pow(R_v, 2);
                        const double u_vortex =         exp_factor_half * (x[1] - y_v) * Gamma_v/pow(R_v, 2);
                        
                        rho_Y_0[idx_cell] = rho_vortex*Y_0_inf;
                        rho_Y_1[idx_cell] = rho_vortex*Y_1_inf;
                        rho_u[idx_cell]   = rho_vortex*u_vortex;
                        rho_v[idx_cell]   = rho_vortex*v_vortex;
                        E[idx_cell]       = p_vortex/(gamma_mixture - double(1)) + double(1)/double(2)*rho_vortex*(u_vortex*u_vortex + v_vortex*v_vortex);
                    }
                }
            }
        }
    }
}
