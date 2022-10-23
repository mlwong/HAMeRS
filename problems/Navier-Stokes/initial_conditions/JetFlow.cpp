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
    // Uniform multi-species jet
    NULL_USE(data_time);
    
    if ((d_project_name != "2D jet") && (d_project_name != "3D jet"))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Can only initialize data for 'project_name' = '2D jet' or '3D jet'!\n"
            << "'project_name' = '"
            << d_project_name
            << "' is given."
            << std::endl);
    }
    
    if (d_dim != tbox::Dimension(2) && d_dim != tbox::Dimension(3)) 
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Dimension of problem should be 2 or 3!"
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
         * Initialize data for a 2D Rayleigh-Taylor instability problem (At = 0.04, M = 0.3).
         */
        
        HAMERS_SHARED_PTR<pdat::CellData<double> > density         = conservative_variables[0];
        HAMERS_SHARED_PTR<pdat::CellData<double> > momentum        = conservative_variables[1];
        HAMERS_SHARED_PTR<pdat::CellData<double> > total_energy    = conservative_variables[2];
        
        double* rho     = density->getPointer(0);
        double* rho_u   = momentum->getPointer(0);
        double* rho_v   = momentum->getPointer(1);
        // if(d_project_name == "3D jet")
        //{
        double* rho_w   = momentum->getPointer(2);
        //}

        double* E       = total_energy->getPointer(0);
        
        const double gamma = double(7)/double(5); // assume both gases have the same ratio of specific heat ratios
        // const double gamma_0 = double(7)/double(5);
        // const double gamma_1 = double(7)/double(5);
        
        //double lambda = 701.53278340668; // wavelength of single-mode perturbation
        //double eta_0  = 0.01*lambda;      // 1% perturbation
        // const double eta_0  = 0.0*lambda;      // no perturbation
        
        const double W_0   = 1.00000; // molecular weight of heavier gas
        
        const double p_ref = 100000.0; // interface pressure
        const double T_ref = 300.0;    // background temperature
        
        TBOX_ASSERT(d_initial_conditions_db != nullptr);
        TBOX_ASSERT(d_initial_conditions_db->keyExists("gravity"));
        
        // std::vector<double> gravity_vector = d_initial_conditions_db->getDoubleVector("gravity");
        // const double g = gravity_vector[0]; // gravity
        
        const double R_u = 8.31446261815324; // universal gas constant
        const double R_0 = R_u/W_0;          // gas constant of heavier gas
        
        // const double rho_i = p_i/(R_u*T_0)*(W_1 + W_2)/2.0;
        
        if (d_project_name == "2D jet")
        {
            const double u_ref = 0.0;
            const double v_ref = 0.0;
            const double Z_ref = 0.0;
                    
            const double rho_ref = p_ref/(R_0*T_ref);

            const double rho_u_ref = rho_ref * u_ref;
            const double rho_v_ref = rho_ref * v_ref;
            const double E_ref     = p_ref/(gamma - double(1)) + double(1)/double(2)*rho_ref*(u_ref*u_ref + v_ref*v_ref);

            for (int j = 0; j < patch_dims[1]; j++)
            {
                for (int i = 0; i < patch_dims[0]; i++)
                {
                    // Compute index into linear data array.
                    int idx_cell = i + j*patch_dims[0];
                    
                    rho[idx_cell]   = rho_ref;
                    rho_u[idx_cell] = rho_u_ref;
                    rho_v[idx_cell] = rho_v_ref;
                    E[idx_cell]     = E_ref;
                }
            }
        }
        else if (d_project_name == "3D jet")
        {
            const double u_ref = 0.0;
            const double v_ref = 0.0;
            const double w_ref = 0.0;
            const double Z_ref = 0.0;

            const double rho_ref = p_ref/(R_0*T_ref);
               
            const double rho_u_ref = rho_ref * u_ref;
            const double rho_v_ref = rho_ref * v_ref;
            const double rho_w_ref = rho_ref * w_ref;

            const double E_ref     = p_ref/(gamma - double(1)) + double(1)/double(2)*rho_ref*(u_ref*u_ref + v_ref*v_ref + w_ref*w_ref);

            for (int k = 0; k < patch_dims[2]; k++)
            {
                for (int j = 0; j < patch_dims[1]; j++)
                {
                    for (int i = 0; i < patch_dims[0]; i++)
                    {
                        // Compute index into linear data array.
                        int idx_cell    =  i + j*patch_dims[0] +  k*patch_dims[0]*patch_dims[1];

                        rho[idx_cell]   = rho_ref;
                        rho_u[idx_cell] = rho_u_ref;
                        rho_v[idx_cell] = rho_v_ref;
                        rho_w[idx_cell] = rho_w_ref;
                        E[idx_cell]     = E_ref;
                    }
                }
            }
        } 
    }
}
