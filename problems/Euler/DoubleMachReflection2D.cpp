#include "apps/Euler/EulerInitialConditions.hpp"

/*
 * Set the data on the patch interior to some initial values.
 */
void
EulerInitialConditions::initializeDataOnPatch(
    hier::Patch& patch,
    const std::vector<boost::shared_ptr<pdat::CellVariable<double> > >& conservative_variables,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double data_time,
    const bool initial_time)
{
    NULL_USE(data_time);
    
    if (d_project_name != "2D double-Mach reflection")
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Can only initialize data for 'project_name' = '2D double-Mach reflection'!\n"
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
         * Initialize data for a 2D double-Mach reflection problem.
         */
        
        boost::shared_ptr<pdat::CellVariable<double> > var_density      = conservative_variables[0];
        boost::shared_ptr<pdat::CellVariable<double> > var_momentum     = conservative_variables[1];
        boost::shared_ptr<pdat::CellVariable<double> > var_total_energy = conservative_variables[2];
        
        boost::shared_ptr<pdat::CellData<double> > density(
            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                patch.getPatchData(var_density, data_context)));
        
        boost::shared_ptr<pdat::CellData<double> > momentum(
            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                patch.getPatchData(var_momentum, data_context)));
        
        boost::shared_ptr<pdat::CellData<double> > total_energy(
            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                patch.getPatchData(var_total_energy, data_context)));
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(density);
        TBOX_ASSERT(momentum);
        TBOX_ASSERT(total_energy);
#endif
        
        double* rho   = density->getPointer(0);
        double* rho_u = momentum->getPointer(0);
        double* rho_v = momentum->getPointer(1);
        double* E     = total_energy->getPointer(0);
        
        const double x_0 = double(1)/double(6);
        
        double gamma = double(7)/double(5);
        
        const double rho_post_shock = double(8);
        const double u_post_shock   = double(33)/double(4)*cos(M_PI/double(6));
        const double v_post_shock   = -double(33)/double(4)*sin(M_PI/double(6));
        const double p_post_shock   = double(233)/double(2);
        
        const double rho_pre_shock  = double(7)/double(5);
        const double u_pre_shock    = double(0);
        const double v_pre_shock    = double(0);
        const double p_pre_shock    = double(1);
        
        const double rho_u_post_shock = rho_post_shock*u_post_shock;
        const double rho_v_post_shock = rho_post_shock*v_post_shock;
        
        const double rho_u_pre_shock  = rho_pre_shock*u_pre_shock;
        const double rho_v_pre_shock  = rho_pre_shock*v_pre_shock;
        
        const double E_pre_shock = p_pre_shock/(gamma - double(1)) +
            double(1)/double(2)*rho_pre_shock*(u_pre_shock*u_pre_shock + v_pre_shock*v_pre_shock);
        
        const double E_post_shock = p_post_shock/(gamma - double(1)) +
            double(1)/double(2)*rho_post_shock*(u_post_shock*u_post_shock + v_post_shock*v_post_shock);
            
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
                
                if (x[0] < x_0 + x[1]*sqrt(double(1)/double(3)))
                {
                    rho[idx_cell] = rho_post_shock;
                    rho_u[idx_cell] = rho_u_post_shock;
                    rho_v[idx_cell] = rho_v_post_shock;
                    E[idx_cell] = E_post_shock;
                }
                else
                {
                    rho[idx_cell] = rho_pre_shock;
                    rho_u[idx_cell] = rho_u_pre_shock;
                    rho_v[idx_cell] = rho_v_pre_shock;
                    E[idx_cell] = E_pre_shock;
                }
            }
        }
    }
}
