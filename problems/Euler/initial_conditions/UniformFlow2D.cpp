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
    
    if (d_project_name != "2D uniform flow")
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Can only initialize data for 'project_name' = '2D uniform flow'!\n"
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
        
        // Initial conditions.
        double rho_inf = double(1);
        double u_inf   = double(1);
        double v_inf   = double(1);
        double p_inf   = double(1);
        
        for (int j = 0; j < patch_dims[1]; j++)
        {
            for (int i = 0; i < patch_dims[0]; i++)
            {
                // Compute index into linear data array.
                int idx_cell = i + j*patch_dims[0];
                
                rho[idx_cell]   = rho_inf;
                rho_u[idx_cell] = rho_inf*u_inf;
                rho_v[idx_cell] = rho_inf*v_inf;
                E[idx_cell]     = p_inf/(gamma - double(1)) + double(1)/double(2)*rho_inf*
                    (u_inf*u_inf + v_inf*v_inf);
            }
        }
    }
}
