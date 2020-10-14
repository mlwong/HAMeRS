#include "apps/Euler/EulerSpecialBoundaryConditions.hpp"

/*
 * Set the data on the patch physical boundary to some values, depending on the flow problems
 * and flow models.
 */
void
EulerSpecialBoundaryConditions::setSpecialBoundaryConditions(
    hier::Patch& patch,
    const std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& conservative_variables,
    const double fill_time,
    const hier::IntVector& ghost_width_to_fill)
{
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
    
    // Get the number of ghost cells of gradient.
    hier::IntVector num_ghosts = conservative_variables[0]->getGhostCellWidth();
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of patch plus
    // ghost cells.
    const hier::Box ghost_box = conservative_variables[0]->getGhostBox();
    const hier::IntVector ghostcell_dims = ghost_box.numberCells();
    
    const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
        HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(patch_geom);
#endif
    
    if (patch_geom->getTouchesRegularBoundary(1, 0) ||
        patch_geom->getTouchesRegularBoundary(1, 1))
    {
        const double* const dx = patch_geom->getDx();
        const double* const patch_xlo = patch_geom->getXLower();
        
        HAMERS_SHARED_PTR<pdat::CellData<double> > density = conservative_variables[0];
        HAMERS_SHARED_PTR<pdat::CellData<double> > momentum = conservative_variables[1];
        HAMERS_SHARED_PTR<pdat::CellData<double> > total_energy = conservative_variables[2];
        
        double* rho   = density->getPointer(0);
        double* rho_u = momentum->getPointer(0);
        double* rho_v = momentum->getPointer(1);
        double* E     = total_energy->getPointer(0);
        
        const double x_0 = double(1)/double(6);
        
        double gamma = double(7)/double(5);
        
        const double rho_post_shock = double(8);
        const double u_post_shock = double(33)/double(4)*cos(M_PI/double(6));
        const double v_post_shock = -double(33)/double(4)*sin(M_PI/double(6));
        const double p_post_shock = double(233)/double(2);
        
        const double rho_pre_shock = double(7)/double(5);
        const double u_pre_shock = double(0);
        const double v_pre_shock = double(0);
        const double p_pre_shock = double(1);
        
        const double rho_u_post_shock = rho_post_shock*u_post_shock;
        const double rho_v_post_shock = rho_post_shock*v_post_shock;
        
        const double rho_u_pre_shock = rho_pre_shock*u_pre_shock;
        const double rho_v_pre_shock = rho_pre_shock*v_pre_shock;
        
        const double E_pre_shock = p_pre_shock/(gamma - double(1)) +
            double(1)/double(2)*rho_pre_shock*(u_pre_shock*u_pre_shock + v_pre_shock*v_pre_shock);
        
        const double E_post_shock = p_post_shock/(gamma - double(1)) +
            double(1)/double(2)*rho_post_shock*(u_post_shock*u_post_shock + v_post_shock*v_post_shock);
        
        /*
         * Update the bottom boundary conditions.
         */
        
        if (patch_geom->getTouchesRegularBoundary(1, 0))
        {
            for (int i = 0; i < interior_dims[0]; i++)
            {
                for (int j = -ghost_width_to_fill[1];
                     j < 0;
                     j++)
                {
                    const int idx_cell = (i + num_ghosts[0]) +
                        (j + num_ghosts[1])*ghostcell_dims[0];
                    
                    // Compute the coordinates.
                    double x[2];
                    x[0] = patch_xlo[0] + (i + double(1)/double(2))*dx[0];
                    x[1] = patch_xlo[1] + (j + double(1)/double(2))*dx[1];
                    
                    if (x[0] < x_0)
                    {
                        rho[idx_cell] = rho_post_shock;
                        rho_u[idx_cell] = rho_u_post_shock;
                        rho_v[idx_cell] = rho_v_post_shock;
                        E[idx_cell] = E_post_shock;
                    }
                    else
                    {
                        const int idx_mirror_cell = (i + num_ghosts[0]) +
                            (-j + num_ghosts[1] - 1)*ghostcell_dims[0];
                        
                        rho[idx_cell] = rho[idx_mirror_cell];
                        rho_u[idx_cell] = rho_u[idx_mirror_cell];
                        rho_v[idx_cell] = -rho_v[idx_mirror_cell];
                        E[idx_cell] = E[idx_mirror_cell];
                    }
                }
            }
        }
        
        /*
         * Update the top boundary conditions.
         */
        
        if (patch_geom->getTouchesRegularBoundary(1, 1))
        {
            const double x_s = x_0 + (double(1) + double(20)*fill_time)/sqrt(double(3));
            
            for (int i = 0; i < interior_dims[0]; i++)
            {
                for (int j = interior_dims[1];
                     j < interior_dims[1] + ghost_width_to_fill[1];
                     j++)
                {
                    const int idx_cell = (i + num_ghosts[0]) +
                        (j + num_ghosts[1])*ghostcell_dims[0];
                    
                    // Compute the coordinates.
                    double x[2];
                    x[0] = patch_xlo[0] + (i + double(1)/double(2))*dx[0];
                    x[1] = patch_xlo[1] + (j + double(1)/double(2))*dx[1];
                    
                    if (x[0] >= x_s)
                    {
                        rho[idx_cell] = rho_pre_shock;
                        rho_u[idx_cell] = rho_u_pre_shock;
                        rho_v[idx_cell] = rho_v_pre_shock;
                        E[idx_cell] = E_pre_shock;
                    }
                    else
                    {
                        rho[idx_cell] = rho_post_shock;
                        rho_u[idx_cell] = rho_u_post_shock;
                        rho_v[idx_cell] = rho_v_post_shock;
                        E[idx_cell] = E_post_shock;
                    }
                }
            }
        }
    }
}
