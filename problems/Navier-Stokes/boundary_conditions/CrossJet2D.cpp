#include "apps/Navier-Stokes/NavierStokesSpecialBoundaryConditions.hpp"

/*
 * Set the data on the patch physical boundary to some values, depending on the flow problems
 * and flow models.
 */
void
NavierStokesSpecialBoundaryConditions::setSpecialBoundaryConditions(
    hier::Patch& patch,
    const std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& conservative_variables,
    const double fill_time,
    const hier::IntVector& ghost_width_to_fill)
{
    if (d_project_name != "2D cross jet")
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Can only initialize data for 'project_name' = "
            << "'2D constant uniform flow'!\n"
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
            << "Flow model should be conservative single-species models!"
            << std::endl);
    }
    
    if (d_flow_model->getNumberOfSpecies() != 1)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Number of species should be 1!"
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

    TBOX_ASSERT(d_special_boundary_conditions_db != nullptr);
    TBOX_ASSERT(d_special_boundary_conditions_db->keyExists("p_inf"));
    TBOX_ASSERT(d_special_boundary_conditions_db->keyExists("u_inf"));
    TBOX_ASSERT(d_special_boundary_conditions_db->keyExists("v_inf"));
    TBOX_ASSERT(d_special_boundary_conditions_db->keyExists("gamma"));
    TBOX_ASSERT(d_special_boundary_conditions_db->keyExists("x_a"));
    TBOX_ASSERT(d_special_boundary_conditions_db->keyExists("x_b"));
    TBOX_ASSERT(d_special_boundary_conditions_db->keyExists("x_al"));
    TBOX_ASSERT(d_special_boundary_conditions_db->keyExists("x_br"));
    TBOX_ASSERT(d_special_boundary_conditions_db->keyExists("x_j_c"));
    TBOX_ASSERT(d_special_boundary_conditions_db->keyExists("V_jet"));
    
    const double p_inf = d_special_boundary_conditions_db->getDouble("p_inf");
    const double u_inf = d_special_boundary_conditions_db->getDouble("u_inf"); 
    const double v_inf = d_special_boundary_conditions_db->getDouble("v_inf");
    const double V_jet = d_special_boundary_conditions_db->getDouble("V_jet");
    const double gamma = d_special_boundary_conditions_db->getDouble("gamma");
    const double x_a   = d_special_boundary_conditions_db->getDouble("x_a");
    const double x_b   = d_special_boundary_conditions_db->getDouble("x_b");
    const double x_al  = d_special_boundary_conditions_db->getDouble("x_al");
    const double x_br  = d_special_boundary_conditions_db->getDouble("x_br");
    const double x_j_c = d_special_boundary_conditions_db->getDouble("x_j_c");
    
    /*
     * Determine the ghost cell width to fill.
     */
    hier::IntVector gcw_to_fill(tbox::Dimension(2));
    
    // If the ghost fill width is not used, it is set to the ghost cell width of the data.
    if (ghost_width_to_fill == -hier::IntVector::getOne(tbox::Dimension(2)))
    {
        gcw_to_fill = num_ghosts;
    }
    else
    {
        gcw_to_fill = hier::IntVector::min(
            num_ghosts,
            ghost_width_to_fill);
    }
    
    if (d_project_name == "2D cross jet")
    {
        // New code
        const double* const dx = patch_geom->getDx();
        const double* const patch_xlo = patch_geom->getXLower();
        
        HAMERS_SHARED_PTR<pdat::CellData<Real> > density      = conservative_variables[0];
        HAMERS_SHARED_PTR<pdat::CellData<Real> > momentum     = conservative_variables[1];
        HAMERS_SHARED_PTR<pdat::CellData<Real> > total_energy = conservative_variables[2];
        
        Real* rho         = density->getPointer(0);
        Real* rho_u       = momentum->getPointer(0);
        Real* rho_v       = momentum->getPointer(1);
        Real* E           = total_energy->getPointer(0);
        
        for (int codim = 1; codim <= d_dim.getValue(); codim++)
        {
            const std::vector<hier::BoundaryBox>& boundary_boxes = patch_geom->getCodimensionBoundaries(codim);
            
            if (!boundary_boxes.empty())
            {
                // Loop over the boundary boxes.
                for (int bi = 0; bi < static_cast<int>(boundary_boxes.size()); bi++)
                {
                    hier::Box fill_box(patch_geom->getBoundaryFillBox(
                        boundary_boxes[bi],
                        interior_box,
                        gcw_to_fill));
                    
                    hier::Index fill_box_lo_idx(fill_box.lower());
                    hier::Index fill_box_hi_idx(fill_box.upper());
                    
                    /*
                     * Offset the indices.
                     */
                    fill_box_lo_idx = fill_box_lo_idx - interior_box.lower();
                    fill_box_hi_idx = fill_box_hi_idx - interior_box.lower();
                    
                    for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                    {
                        for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                        {
                            const int idx_cell = (i + num_ghosts[0]) +
                                (j + num_ghosts[1])*ghostcell_dims[0];
                            
                            const int idx_mirror_cell = (i + num_ghosts[0]) +
                            (-j + num_ghosts[1] - 1)*ghostcell_dims[0];
                            
                            // Compute the coordinates.
                            double x[2];
                            x[0] = patch_xlo[0] + (i + double(1)/double(2))*dx[0];
                            x[1] = patch_xlo[1] + (j + double(1)/double(2))*dx[1];
                            double u = 0.0;
                            if(x[1] < 0.0)
                            {
                                if (x[0] < x_a)
                                {
                                    u               = (u_inf*(x_a-x[0])/(x_a - x_al));
                                    rho[idx_cell]   = rho[idx_mirror_cell];
                                    rho_u[idx_cell] = 2*rho[idx_cell]*u-rho_u[idx_mirror_cell];
                                    rho_v[idx_cell] = -rho_v[idx_mirror_cell];
                                    E[idx_cell]     = E[idx_mirror_cell];
                                }
                                else if (x[0] > x_b)
                                {
                                    u               = (u_inf*(x[0]-x_b)/(x_br - x_b));
                                    rho[idx_cell]   = rho[idx_mirror_cell];
                                    rho_u[idx_cell] = 2*rho[idx_cell]*u-rho_u[idx_mirror_cell];
                                    rho_v[idx_cell] = -rho_v[idx_mirror_cell];
                                    E[idx_cell]     = E[idx_mirror_cell];
                                }
                                else
                                {
                                    const double r     = fabs(x[0]-x_j_c);
                                    const double r_0   = 0.5;
                                    const double u_ref = 0.0;
                                    
                                    const double theta_0 = 0.1;
                                    const double v_ref = V_jet*0.5*(1.0-tanh(r_0/(4.0*theta_0)*(r/r_0-r_0/r)));
                                    
                                    rho[idx_cell]   = rho[idx_mirror_cell];
                                    
                                    const double rho_ref   = rho[idx_cell];
                                    const double rho_u_ref = rho_ref * u_ref;
                                    const double rho_v_ref = rho_ref * v_ref;
                                    const double E_ref     = p_inf/(gamma - double(1)) + double(1)/double(2)*rho_ref*(u_ref*u_ref + v_ref*v_ref);
                                    
                                    rho_u[idx_cell] = rho_u_ref;
                                    rho_v[idx_cell] = rho_v_ref;
                                    E[idx_cell]     = E_ref;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
