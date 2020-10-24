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
    
    if (d_project_name != "2D Rayleigh-Taylor instability")
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Can only initialize data for 'project_name' = '2D Rayleigh-Taylor instability'!\n"
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
    
    if (d_flow_model_type != FLOW_MODEL::FOUR_EQN_CONSERVATIVE)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Flow model should be conservative four-equation models!"
            << std::endl);
    }
    
    if (d_flow_model->getNumberOfSpecies() != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Number of species should be 2!"
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
         * Initialize data for a 2D Rayleigh-Taylor instability problem.
         */
        
        boost::shared_ptr<pdat::CellData<double> > partial_density = conservative_variables[0];
        boost::shared_ptr<pdat::CellData<double> > momentum        = conservative_variables[1];
        boost::shared_ptr<pdat::CellData<double> > total_energy    = conservative_variables[2];
        
        double* rho_Y_0 = partial_density->getPointer(0);
        double* rho_Y_1 = partial_density->getPointer(1);
        double* rho_u   = momentum->getPointer(0);
        double* rho_v   = momentum->getPointer(1);
        double* E       = total_energy->getPointer(0);
        
        const double gamma = double(7)/double(5);
        // const double gamma_0 = double(7)/double(5);
        // const double gamma_1 = double(7)/double(5);
        
        const double lambda = 701.53278340668;
        const double eta_0  = 0.01*lambda;
        
        const double W_1 = 0.03328;
        const double W_2 = 0.03072;
        
        const double p_i = 100000.0;
        const double T_0 = 300.0;
        
        const double g   = 10.0;
        const double R_u = 8.31446261815324;
        
        const double R_1 = R_u/W_1;
        const double R_2 = R_u/W_2;
        
        // const double rho_i = p_i/(R_u*T_0)*(W_1 + W_2)/2.0;
        
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
                
                // const double delta = eta_0*cos(2.0*M_PI/lambda*x[1]);
                const double delta = 0.0;
                
                if (x[0] < delta) // heavier fluid
                {
                    const double rho = p_i/(R_1*T_0)*exp((g*x[0])/(R_1*T_0));
                    rho_Y_0[idx_cell] = rho;
                    rho_Y_1[idx_cell] = 0.0;
                    
                    const double p = p_i*exp((g*x[0])/(R_1*T_0));
                    
                    const double u = 0.0;
                    const double v = 0.0;
                    
                    rho_u[idx_cell] = rho*u;
                    rho_v[idx_cell] = rho*v;
                    E[idx_cell]     = p/(gamma - double(1)) + double(1)/double(2)*rho*(u*u + v*v);
                }
                else // lighter fluid
                {
                    const double rho = p_i/(R_2*T_0)*exp((g*x[0])/(R_2*T_0));
                    rho_Y_0[idx_cell] = 0.0;
                    rho_Y_1[idx_cell] = rho;
                    
                    const double p = p_i*exp((g*x[0])/(R_2*T_0));
                    
                    const double u = 0.0;
                    const double v = 0.0;
                    
                    rho_u[idx_cell] = rho*u;
                    rho_v[idx_cell] = rho*v;
                    E[idx_cell]     = p/(gamma - double(1)) + double(1)/double(2)*rho*(u*u + v*v);
                }
            }
        }
    }
}
