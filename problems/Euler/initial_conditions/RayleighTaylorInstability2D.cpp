#include "apps/Euler/EulerInitialConditions.hpp"

/*
 * Set the data on the patch interior to some initial values.
 */
void
EulerInitialConditions::initializeDataOnPatch(
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
    
    if ((d_flow_model_type != FLOW_MODEL::SINGLE_SPECIES) && (d_flow_model_type != FLOW_MODEL::FOUR_EQN_CONSERVATIVE))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Flow model should be single-species or conservative four-equation models!"
            << std::endl);
    }
    
    if (d_flow_model_type == FLOW_MODEL::FOUR_EQN_CONSERVATIVE)
    {
        if (d_flow_model->getNumberOfSpecies() != 2)
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Number of species should be 2!"
                << std::endl);
        }
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
        
        if (d_flow_model_type == FLOW_MODEL::SINGLE_SPECIES)
        {
            boost::shared_ptr<pdat::CellData<double> > density      = conservative_variables[0];
            boost::shared_ptr<pdat::CellData<double> > momentum     = conservative_variables[1];
            boost::shared_ptr<pdat::CellData<double> > total_energy = conservative_variables[2];
            
            double* rho   = density->getPointer(0);
            double* rho_u = momentum->getPointer(0);
            double* rho_v = momentum->getPointer(1);
            double* E     = total_energy->getPointer(0);
            
            const double gamma = double(5)/double(3);
            
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
                    
                    if (x[1] < double(1)/double(2))
                    {
                        rho[idx_cell] = double(2);
                        
                        const double p = double(2)*x[1] + double(1);
                        const double u = double(0);
                        
                        const double a = sqrt(gamma*p/rho[idx_cell]);
                        const double v = -double(25)/double(1000)*a*cos(double(8)*M_PI*x[0]);
                        
                        rho_u[idx_cell] = rho[idx_cell]*u;
                        rho_v[idx_cell] = rho[idx_cell]*v;
                        E[idx_cell]     = p/(gamma - double(1)) + double(1)/double(2)*rho[idx_cell]*(u*u + v*v);
                    }
                    else
                    {
                        rho[idx_cell] = double(1);
                        
                        const double p = x[1] + double(3)/double(2);
                        const double u = double(0);
                        
                        const double a = sqrt(gamma*p/rho[idx_cell]);
                        const double v = -double(25)/double(1000)*a*cos(double(8)*M_PI*x[0]);
                        
                        rho_u[idx_cell] = rho[idx_cell]*u;
                        rho_v[idx_cell] = rho[idx_cell]*v;
                        E[idx_cell]     = p/(gamma - double(1)) + double(1)/double(2)*rho[idx_cell]*(u*u + v*v);
                    }
                }
            }
        }
        else if (d_flow_model_type == FLOW_MODEL::FOUR_EQN_CONSERVATIVE)
        {
            boost::shared_ptr<pdat::CellData<double> > partial_density = conservative_variables[0];
            boost::shared_ptr<pdat::CellData<double> > momentum     = conservative_variables[1];
            boost::shared_ptr<pdat::CellData<double> > total_energy = conservative_variables[2];
            
            double* rho_Y_0   = partial_density->getPointer(0);
            double* rho_Y_1   = partial_density->getPointer(1);
            double* rho_u = momentum->getPointer(0);
            double* rho_v = momentum->getPointer(1);
            double* E     = total_energy->getPointer(0);
            
            const double gamma_0 = double(5)/double(3);
            const double gamma_1 = double(7)/double(5);
            
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
                    
                    if (x[1] < double(1)/double(2))
                    {
                        rho_Y_0[idx_cell] = double(2);
                        rho_Y_1[idx_cell] = double(0);
                        
                        const double rho = rho_Y_0[idx_cell] + rho_Y_1[idx_cell];
                        const double p   = double(2)*x[1] + double(1);
                        const double u   = double(0);
                        
                        const double a = sqrt(gamma_0*p/rho);
                        const double v = -double(25)/double(1000)*a*cos(double(8)*M_PI*x[0]);
                        
                        rho_u[idx_cell] = rho*u;
                        rho_v[idx_cell] = rho*v;
                        E[idx_cell]     = p/(gamma_0 - double(1)) + double(1)/double(2)*rho*(u*u + v*v);
                    }
                    else
                    {
                        rho_Y_0[idx_cell] = double(0);
                        rho_Y_1[idx_cell] = double(1);
                        
                        const double rho = rho_Y_0[idx_cell] + rho_Y_1[idx_cell];
                        const double p   = x[1] + double(3)/double(2);
                        const double u   = double(0);
                        
                        const double a = sqrt(gamma_1*p/rho);
                        const double v = -double(25)/double(1000)*a*cos(double(8)*M_PI*x[0]);
                        
                        rho_u[idx_cell] = rho*u;
                        rho_v[idx_cell] = rho*v;
                        E[idx_cell]     = p/(gamma_1 - double(1)) + double(1)/double(2)*rho*(u*u + v*v);
                    }
                }
            }
        }
        
    }
}
