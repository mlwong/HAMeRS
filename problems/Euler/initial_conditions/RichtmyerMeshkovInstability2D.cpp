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
    
    if (d_project_name != "2D Richtmyer-Meshkov instability")
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Can only initialize data for 'project_name' = '2D Richtmyer-Meshkov instability'!\n"
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
            << "Flow model should be conservative four-equation model!"
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
         * Initialize data for a 2D Richtmyer-Meshkov instability problem.
         */
        
        boost::shared_ptr<pdat::CellData<double> > partial_density = conservative_variables[0];
        boost::shared_ptr<pdat::CellData<double> > momentum        = conservative_variables[1];
        boost::shared_ptr<pdat::CellData<double> > total_energy    = conservative_variables[2];
        
        double* rho_Y_0   = partial_density->getPointer(0);
        double* rho_Y_1   = partial_density->getPointer(1);
        double* rho_u     = momentum->getPointer(0);
        double* rho_v     = momentum->getPointer(1);
        double* E         = total_energy->getPointer(0);
        
        // Define and compute the characteristic lengths of the problem.
        const double D = double(0.001);
        const double epsilon_i = double(6)/double(128)*D;
        
        // species 0: SF6
        // species 1: air
        const double gamma_0 = double(1.09312);
        const double gamma_1 = double(1.39909);
        
        const double c_p_0 = double(668.286);
        const double c_p_1 = double(1040.50);
        
        const double c_v_0 = double(611.359);
        const double c_v_1 = double(743.697);
        
        NULL_USE(gamma_0);
        
        // SF6, pre-shock condition.
        const double rho_SF6 = double(5.972856);
        const double u_SF6   = double(436.201332);
        const double v_SF6   = double(0);
        const double p_SF6   = double(101325);
        
        // air, pre-shock condition.
        const double rho_pre = double(1.145598);
        const double u_pre   = double(436.201332);
        const double v_pre   = double(0);
        const double p_pre   = double(101325);
        
        // air, post-shock condition.
        const double rho_post = double(1.616874);
        const double u_post   = double(309.060123);
        const double v_post   = double(0);
        const double p_post   = double(164859);
        
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
                
                if (x[0] > double(7)/double(10)*D)
                {
                    rho_Y_0[idx_cell] = double(0);
                    rho_Y_1[idx_cell] = rho_post;
                    rho_u[idx_cell]   = rho_post*u_post;
                    rho_v[idx_cell]   = rho_post*v_post;
                    E[idx_cell]       = p_post/(gamma_1 - double(1)) + double(1)/double(2)*rho_post*
                        (u_post*u_post + v_post*v_post);
                }
                else
                {
                    // Compute the distance from the initial material interface.
                    const double dR = x[0] - (double(2)/double(5) - double(1)/double(10)*
                        sin(2*M_PI*(x[1]/D + double(1)/double(4))))*D;
                    
                    const double f_sm = double(1)/double(2)*(double(1) + erf(dR/epsilon_i));
                    
                    // Smooth the primitive quantity.
                    const double rho_Y_0_i = rho_SF6*(double(1) - f_sm);
                    const double rho_Y_1_i = rho_pre*f_sm;
                    const double u_i       = u_SF6*(double(1) - f_sm) + u_pre*f_sm;
                    const double v_i       = v_SF6*(double(1) - f_sm) + v_pre*f_sm;
                    const double p_i       = p_SF6*(double(1) - f_sm) + p_pre*f_sm;
                    
                    const double rho_i = rho_Y_0_i + rho_Y_1_i;
                    const double Y_0_i = rho_Y_0_i/rho_i;
                    const double Y_1_i = double(1) - Y_0_i;
                    
                    const double gamma = (Y_0_i*c_p_0 + Y_1_i*c_p_1)/(Y_0_i*c_v_0 + Y_1_i*c_v_1);
                    
                    rho_Y_0[idx_cell] = rho_Y_0_i;
                    rho_Y_1[idx_cell] = rho_Y_1_i;
                    rho_u[idx_cell]   = rho_i*u_i;
                    rho_v[idx_cell]   = rho_i*v_i;
                    E[idx_cell]       = p_i/(gamma - double(1)) + double(1)/double(2)*rho_i*
                        (u_i*u_i + v_i*v_i);
                }
            }
        }
    }
}
