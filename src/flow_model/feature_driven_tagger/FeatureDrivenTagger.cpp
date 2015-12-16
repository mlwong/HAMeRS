#include "flow_model/feature_driven_tagger/FeatureDrivenTagger.hpp"

#define EPSILON 1e-40

FeatureDrivenTagger::FeatureDrivenTagger(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
    const hier::IntVector& num_ghosts,
    const FLOW_MODEL& flow_model,
    const int& num_species,
    const boost::shared_ptr<EquationOfState>& equation_of_state,
    const boost::shared_ptr<tbox::Database>& feature_driven_tagger_db):
        d_object_name(object_name),
        d_dim(dim),
        d_grid_geometry(grid_geometry),
        d_num_ghosts(num_ghosts),
        d_flow_model(flow_model),
        d_num_species(num_species),
        d_equation_of_state(equation_of_state),
        d_density(NULL),
        d_partial_density(NULL),
        d_momentum(NULL),
        d_total_energy(NULL),
        d_mass_fraction(NULL),
        d_volume_fraction(NULL),
        d_variables_set(false),
        d_shock_Jameson_tol(0.01),
        d_shock_Ducros_tol(0.01),
        d_shock_Larsson_tol(0.7),
        d_density_Jameson_tol(0.01),
        d_vorticity_Q_criterion_tol(0.01)
{
    if (feature_driven_tagger_db != nullptr)
    {
        std::vector<std::string> refinement_keys = feature_driven_tagger_db->getAllKeys();
        int num_keys = static_cast<int>(refinement_keys.size());
        
        if (feature_driven_tagger_db->keyExists("refinement_criteria"))
        {
            d_refinement_criteria = feature_driven_tagger_db->getStringVector("refinement_criteria");
        }
        else if (feature_driven_tagger_db->keyExists("d_refinement_criteria"))
        {
            d_refinement_criteria = feature_driven_tagger_db->getStringVector("d_refinement_criteria");
        }
        else
        {
            TBOX_WARNING(
                d_object_name << ": "
                              << "No key 'refinement_criteria'/'d_refinement_criteria' found in data for"
                              << " Refinement_data. No refinement will occur."
                              << std::endl);
        }
        
        std::vector<std::string> ref_keys_defined(num_keys);
        int def_key_cnt = 0;
        boost::shared_ptr<tbox::Database> error_db;
        for (int i = 0; i < num_keys; i++)
        {
            std::string error_key = refinement_keys[i];
            error_db.reset();
            
            if (!((error_key == "refinement_criteria") || (error_key == "d_refinement_criteria")))
            {
                if (!((error_key == "SHOCK_JAMESON") ||
                      (error_key == "SHOCK_DUCROS") ||
                      (error_key == "SHOCK_LARSSON") ||
                      (error_key == "DENSITY_JAMESON") ||
                      (error_key == "VORTICITY_Q_CRITERION")))
                {
                    TBOX_ERROR(d_object_name
                               << ": "
                               << "Unknown refinement criteria: "
                               << error_key
                               << "\nin input."
                               << std::endl);
                }
                else
                {
                    error_db = feature_driven_tagger_db->getDatabase(error_key);
                    ref_keys_defined[def_key_cnt] = error_key;
                    ++def_key_cnt;
                }
                
                if (error_db && error_key == "SHOCK_JAMESON")
                {                    
                    if (error_db->keyExists("shock_Jameson_tol"))
                    {
                        d_shock_Jameson_tol = error_db->getDouble("shock_Jameson_tol");
                    }
                    else if (error_db->keyExists("d_shock_Jameson_tol"))
                    {
                        d_shock_Jameson_tol = error_db->getDouble("d_shock_Jameson_tol");
                    }
                    else
                    {
                        TBOX_ERROR(d_object_name
                                   << ": "
                                   << "No key 'shock_Jameson_tol'/'d_shock_Jameson_tol' found in data for "
                                   << error_key
                                   << "."
                                   << std::endl);
                    }
                }
                else if (error_db && error_key == "SHOCK_DUCROS")
                {                    
                    if (error_db->keyExists("shock_Ducros_tol"))
                    {
                        d_shock_Ducros_tol = error_db->getDouble("shock_Ducros_tol");
                    }
                    else if (error_db->keyExists("d_shock_Ducros_tol"))
                    {
                        d_shock_Ducros_tol = error_db->getDouble("d_shock_Ducros_tol");
                    }
                    else
                    {
                        TBOX_ERROR(d_object_name
                                   << ": "
                                   << "No key 'shock_Ducros_tol'/'d_shock_Ducros_tol' found in data for "
                                   << error_key
                                   << "."
                                   << std::endl);
                    }
                }
                else if (error_db && error_key == "SHOCK_LARSSON")
                {                    
                    if (error_db->keyExists("shock_Larsson_tol"))
                    {
                        d_shock_Larsson_tol = error_db->getDouble("shock_Larsson_tol");
                    }
                    else if (error_db->keyExists("d_shock_Larsson_tol"))
                    {
                        d_shock_Larsson_tol = error_db->getDouble("d_shock_Larsson_tol");
                    }
                    else
                    {
                        TBOX_ERROR(d_object_name
                                   << ": "
                                   << "No key 'shock_Larsson_tol'/'d_shock_Larsson_tol' found in data for "
                                   << error_key
                                   << "."
                                   << std::endl);
                    }
                }
                else if (error_db && error_key == "DENSITY_JAMESON")
                {                    
                    if (error_db->keyExists("density_Jameson_tol"))
                    {
                        d_density_Jameson_tol = error_db->getDouble("density_Jameson_tol");
                    }
                    else if (error_db->keyExists("d_density_Jameson_tol"))
                    {
                        d_density_Jameson_tol = error_db->getDouble("d_density_Jameson_tol");
                    }
                    else
                    {
                        TBOX_ERROR(d_object_name
                                   << ": "
                                   << "No key 'density_Jameson_tol'/'d_density_Jameson_tol' found in data for "
                                   << error_key
                                   << "."
                                   << std::endl);
                    }
                }
                else if (error_db && error_key == "VORTICITY_Q_CRITERION")
                {        
                    if (error_db->keyExists("vorticity_Q_criterion_tol"))
                    {
                        d_vorticity_Q_criterion_tol = error_db->getDouble("vorticity_Q_criterion_tol");
                    }
                    else if (error_db->keyExists("d_vorticity_Q_criterion_tol"))
                    {
                        d_vorticity_Q_criterion_tol = error_db->getDouble("d_vorticity_Q_criterion_tol");
                    }
                    else
                    {
                        TBOX_ERROR(d_object_name
                                   << ": "
                                   << "No key 'vorticity_Q_criterion_tol'/'d_vorticity_Q_criterion_tol' found in data for "
                                   << error_key
                                   << "."
                                   << std::endl);
                    }
                }
            }
        } // loop over refine criteria
        
        /*
         * Check that input is found for each string identifier in key list.
         */
        for (int k0 = 0;
             k0 < static_cast<int>(d_refinement_criteria.size());
             k0++)
        {
            std::string use_key = d_refinement_criteria[k0];
            bool key_found = false;
            for (int k1 = 0; k1 < def_key_cnt; k1++)
            {
                std::string def_key = ref_keys_defined[k1];
                if (def_key == use_key)
                {
                    key_found = true;
                }
            }
            
            if (!key_found)
            {
                TBOX_ERROR(d_object_name
                           << ": "
                           << "No input found for specified refine criteria: "
                           << d_refinement_criteria[k0]
                           << "."
                           << std::endl);
            }
        }
    }
    else
    {
        TBOX_WARNING(d_object_name
            << ": "
            << "Key data 'Feature_driven_tagger' not found in input/restart database."
            << " No refinement will occur."
            << std::endl);
    }
}


/*
 * Print all characteristics of the feature driven tagger class.
 */
void
FeatureDrivenTagger::printClassData(std::ostream& os) const
{
    NULL_USE(os);
}


/*
 * Put the characteristics of the feature driven tagger into the restart
 * database.
 */
void
FeatureDrivenTagger::putToRestart(
    const boost::shared_ptr<tbox::Database>& restart_db) const
{
    if (d_refinement_criteria.size() > 0)
    {
        restart_db->putStringVector("d_refinement_criteria", d_refinement_criteria);
    }
    for (int i = 0; i < static_cast<int>(d_refinement_criteria.size()); i++)
    {
        if (d_refinement_criteria[i] == "SHOCK_JAMESON")
        {
            boost::shared_ptr<tbox::Database> error_db =
                restart_db->putDatabase("SHOCK_JAMESON");
                
            error_db->putDouble("d_shock_Jameson_tol",
                d_shock_Jameson_tol);
        }
        else if (d_refinement_criteria[i] == "SHOCK_DUCROS")
        {
            boost::shared_ptr<tbox::Database> error_db =
                restart_db->putDatabase("SHOCK_DUCROS");
            
            error_db->putDouble("d_shock_Ducros_tol",
                d_shock_Ducros_tol);
        }
        else if (d_refinement_criteria[i] == "SHOCK_LARSSON")
        {
            boost::shared_ptr<tbox::Database> error_db =
                restart_db->putDatabase("SHOCK_LARSSON");
            
            error_db->putDouble("d_shock_Larsson_tol",
                d_shock_Larsson_tol);
        }
        else if (d_refinement_criteria[i] == "DENSITY_JAMESON")
        {
            boost::shared_ptr<tbox::Database> error_db =
                restart_db->putDatabase("DENSITY_JAMESON");
            
            error_db->putDouble("d_density_Jameson_tol",
                d_density_Jameson_tol);
        }
        else if (d_refinement_criteria[i] == "VORTICITY_Q_CRITERION")
        {
            boost::shared_ptr<tbox::Database> error_db =
                restart_db->putDatabase("VORTICITY_Q_CRITERION");
            
            error_db->putDouble("d_vorticity_Q_criterion_tol",
                d_vorticity_Q_criterion_tol);
        }
    }
}


/*
 * Tag cells for refinement using different feature-based detectors.
 */
void
FeatureDrivenTagger::tagCells(
   hier::Patch& patch,
   const double regrid_time,
   const bool initial_error,
   const bool uses_richardson_extrapolation_too,
   boost::shared_ptr<pdat::CellData<int> > tags,
   const boost::shared_ptr<hier::VariableContext>& data_context)
{
    NULL_USE(uses_richardson_extrapolation_too);
    
    const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(patch_geom);
#endif
    
    const double* dx = patch_geom->getDx();
    
    // Get the dimensions of box that covers the interior of Patch.
    hier::Box dummy_box = patch.getBox();
    const hier::Box interior_box = dummy_box;
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of Patch plus
    // ghost cells.
    dummy_box.grow(d_num_ghosts);
    const hier::Box ghost_box = dummy_box;
    const hier::IntVector ghostcell_dims = ghost_box.numberCells();
    
    for (int ncrit = 0;
             ncrit < static_cast<int>(d_refinement_criteria.size());
             ncrit++)
    {
        std::string ref = d_refinement_criteria[ncrit];
        
        // Get the pointer of the tags.
        int* tag_ptr  = tags->getPointer();
        
        if (ref == "SHOCK_JAMESON")
        {
            // Get the cell data of the time-dependent variables.
            boost::shared_ptr<pdat::CellData<double> > density(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_density, data_context)));
            
            boost::shared_ptr<pdat::CellData<double> > momentum(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_momentum, data_context)));
            
            boost::shared_ptr<pdat::CellData<double> > total_energy(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_total_energy, data_context)));
            
            // Allocate temporary patch data.
            boost::shared_ptr<pdat::CellData<double> > pressure(
                new pdat::CellData<double>(interior_box, 1, d_num_ghosts));
            
            if (d_dim == tbox::Dimension(1))
            {
                // NOT YET IMPLEMENTED
            }
            else if (d_dim == tbox::Dimension(2))
            {
                // Get the arrays of conservative variables.
                double* rho   = density->getPointer(0);
                double* rho_u = momentum->getPointer(0);
                double* rho_v = momentum->getPointer(1);
                double* E     = total_energy->getPointer(0);
                
                // Get the array of pressure.
                double* p     = pressure->getPointer(0);
                
                // Compute the field of pressure.
                for (int j = -d_num_ghosts[1]; j < interior_dims[1] + d_num_ghosts[1]; j++)
                {
                    for (int i = -d_num_ghosts[0]; i < interior_dims[0] + d_num_ghosts[0]; i++)
                    {
                        // Compute index into linear data array.
                        const int idx = (i + d_num_ghosts[0]) +
                            (j + d_num_ghosts[1])*ghostcell_dims[0];
                        
                        std::vector<const double*> m_ptr;
                        m_ptr.push_back(&rho_u[idx]);
                        m_ptr.push_back(&rho_v[idx]);
                        
                        p[idx] = d_equation_of_state->getPressure(
                            &rho[idx],
                            m_ptr,
                            &E[idx]);
                    }
                }
                
                // Use the Jameson's switch with threshold to tag cells.
                for (int j = 0; j < interior_dims[1]; j++)
                {
                    for (int i = 0; i < interior_dims[0]; i++)
                    {
                        // Compute indices.
                        const int idx = (i + d_num_ghosts[0]) +
                            (j + d_num_ghosts[1])*ghostcell_dims[0];
                        
                        const int idx_x_L = (i - 1 + d_num_ghosts[0]) +
                            (j + d_num_ghosts[1])*ghostcell_dims[0];
                        
                        const int idx_x_R = (i + 1 + d_num_ghosts[0]) +
                            (j + d_num_ghosts[1])*ghostcell_dims[0];
                        
                        const int idx_y_B = (i + d_num_ghosts[0]) +
                            (j - 1 + d_num_ghosts[1])*ghostcell_dims[0];
                        
                        const int idx_y_T = (i + d_num_ghosts[0]) +
                            (j + 1 + d_num_ghosts[1])*ghostcell_dims[0];
                        
                        const int idx_nghost = i + j*interior_dims[0];
                        
                        const double psi_x = fabs(p[idx_x_R] - 2*p[idx] + p[idx_x_L])/
                            (p[idx_x_R] + 2*p[idx] + p[idx_x_L] + EPSILON);
                        
                        const double psi_y = fabs(p[idx_y_T] - 2*p[idx] + p[idx_y_B])/
                            (p[idx_y_T] + 2*p[idx] + p[idx_y_B] + EPSILON);
                        
                        const double psi = fmax(psi_x, psi_y);
                        
                        if (psi > d_shock_Jameson_tol)
                        {
                            tag_ptr[idx_nghost] |= 1;
                        }
                    }
                }
            }
            else if (d_dim == tbox::Dimension(3))
            {
                // Get the arrays of conservative variables.
                double* rho   = density->getPointer(0);
                double* rho_u = momentum->getPointer(0);
                double* rho_v = momentum->getPointer(1);
                double* rho_w = momentum->getPointer(2);
                double* E     = total_energy->getPointer(0);
                
                // Get the array of pressure.
                double* p     = pressure->getPointer(0);
                
                for (int k = -d_num_ghosts[2]; k < interior_dims[2] + d_num_ghosts[2]; k++)
                {
                    for (int j = -d_num_ghosts[1]; j < interior_dims[1] + d_num_ghosts[1]; j++)
                    {
                        for (int i = -d_num_ghosts[0]; i < interior_dims[0] + d_num_ghosts[0]; i++)
                        {
                            // Compute index into linear data array.
                            const int idx = (i + d_num_ghosts[0]) +
                                (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            
                            std::vector<const double*> m_ptr;
                            m_ptr.push_back(&rho_u[idx]);
                            m_ptr.push_back(&rho_v[idx]);
                            m_ptr.push_back(&rho_w[idx]);
                            
                            p[idx] = d_equation_of_state->getPressure(
                                &rho[idx],
                                m_ptr,
                                &E[idx]);
                        }
                    }
                }
                
                // Use the Jameson's switch with threshold to tag cells.
                for (int k = 0; k < interior_dims[2]; k++)
                {
                    for (int j = 0; j < interior_dims[1]; j++)
                    {
                        for (int i = 0; i < interior_dims[0]; i++)
                        {
                            // Compute indices.
                            const int idx = (i + d_num_ghosts[0]) +
                                (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            
                            const int idx_x_L = (i - 1 + d_num_ghosts[0]) +
                                (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            
                            const int idx_x_R = (i + 1 + d_num_ghosts[0]) +
                                (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            
                            const int idx_y_B = (i + d_num_ghosts[0]) +
                                (j - 1 + d_num_ghosts[1])*ghostcell_dims[0] +
                                (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            
                            const int idx_y_T = (i + d_num_ghosts[0]) +
                                (j + 1 + d_num_ghosts[1])*ghostcell_dims[0] +
                                (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            
                            const int idx_z_B = (i + d_num_ghosts[0]) +
                                (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                (k - 1 + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            
                            const int idx_z_F = (i + d_num_ghosts[0]) +
                                (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                (k + 1 + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            
                            const int idx_nghost = i + j*interior_dims[0] + k*interior_dims[0]*interior_dims[1];
                            
                            const double psi_x = fabs(p[idx_x_R] - 2*p[idx] + p[idx_x_L])/
                                (p[idx_x_R] + 2*p[idx] + p[idx_x_L] + EPSILON);
                            
                            const double psi_y = fabs(p[idx_y_T] - 2*p[idx] + p[idx_y_B])/
                                (p[idx_y_T] + 2*p[idx] + p[idx_y_B] + EPSILON);
                            
                            const double psi_z = fabs(p[idx_z_F] - 2*p[idx] + p[idx_z_B])/
                                (p[idx_z_F] + 2*p[idx] + p[idx_z_B] + EPSILON);
                            
                            const double psi = fmax(fmax(psi_x, psi_y), psi_z);
                            
                            if (psi > d_shock_Jameson_tol)
                            {
                                tag_ptr[idx_nghost] |= 1;
                            }
                        }
                    }
                }
            }
        }
        else if (ref == "SHOCK_DUCROS")
        {
            // Get the cell data of the time-dependent variables.
            boost::shared_ptr<pdat::CellData<double> > density(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_density, data_context)));
            
            boost::shared_ptr<pdat::CellData<double> > momentum(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_momentum, data_context)));
            
            boost::shared_ptr<pdat::CellData<double> > total_energy(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_total_energy, data_context)));
            
            // Allocate temporary patch data.
            boost::shared_ptr<pdat::CellData<double> > pressure(
                new pdat::CellData<double>(interior_box, 1, d_num_ghosts));
            
            boost::shared_ptr<pdat::CellData<double> > velocity(
                new pdat::CellData<double>(interior_box, d_dim.getValue(), d_num_ghosts));
            
            if (d_dim == tbox::Dimension(1))
            {
                // NOT YET IMPLEMENTED
            }
            else if (d_dim == tbox::Dimension(2))
            {
                // Get the arrays of conservative variables.
                double* rho   = density->getPointer(0);
                double* rho_u = momentum->getPointer(0);
                double* rho_v = momentum->getPointer(1);
                double* E     = total_energy->getPointer(0);
                
                // Get the array of pressure.
                double* p     = pressure->getPointer(0);
                double* u     = velocity->getPointer(0);
                double* v     = velocity->getPointer(1);
                
                // Compute the field of pressure.
                for (int j = -d_num_ghosts[1]; j < interior_dims[1] + d_num_ghosts[1]; j++)
                {
                    for (int i = -d_num_ghosts[0]; i < interior_dims[0] + d_num_ghosts[0]; i++)
                    {
                        // Compute index into linear data array.
                        const int idx = (i + d_num_ghosts[0]) +
                            (j + d_num_ghosts[1])*ghostcell_dims[0];
                        
                        std::vector<const double*> m_ptr;
                        m_ptr.push_back(&rho_u[idx]);
                        m_ptr.push_back(&rho_v[idx]);
                        
                        p[idx] = d_equation_of_state->getPressure(
                            &rho[idx],
                            m_ptr,
                            &E[idx]);
                        
                        u[idx] = rho_u[idx]/rho[idx];
                        v[idx] = rho_v[idx]/rho[idx];
                    }
                }
                
                // Use the Ducros switch with threshold to tag cells.
                for (int j = 0; j < interior_dims[1]; j++)
                {
                    for (int i = 0; i < interior_dims[0]; i++)
                    {
                        // Compute indices.
                        const int idx = (i + d_num_ghosts[0]) +
                            (j + d_num_ghosts[1])*ghostcell_dims[0];
                        
                        const int idx_x_L = (i - 1 + d_num_ghosts[0]) +
                            (j + d_num_ghosts[1])*ghostcell_dims[0];
                        
                        const int idx_x_R = (i + 1 + d_num_ghosts[0]) +
                            (j + d_num_ghosts[1])*ghostcell_dims[0];
                        
                        const int idx_y_B = (i + d_num_ghosts[0]) +
                            (j - 1 + d_num_ghosts[1])*ghostcell_dims[0];
                        
                        const int idx_y_T = (i + d_num_ghosts[0]) +
                            (j + 1 + d_num_ghosts[1])*ghostcell_dims[0];
                        
                        const int idx_nghost = i + j*interior_dims[0];
                        
                        // Compute the Jameson's sensor.
                        const double psi_x = fabs(p[idx_x_R] - 2*p[idx] + p[idx_x_L])/
                            (p[idx_x_R] + 2*p[idx] + p[idx_x_L] + EPSILON);
                        
                        const double psi_y = fabs(p[idx_y_T] - 2*p[idx] + p[idx_y_B])/
                            (p[idx_y_T] + 2*p[idx] + p[idx_y_B] + EPSILON);
                        
                        const double psi = fmax(psi_x, psi_y);
                        
                        // Compute another indicator.
                        
                        const double dudx = (u[idx_x_R] - u[idx_x_L])/(2*dx[0]);
                        const double dudy = (u[idx_y_T] - u[idx_y_B])/(2*dx[1]);
                        
                        const double dvdx = (v[idx_x_R] - v[idx_x_L])/(2*dx[0]);
                        const double dvdy = (v[idx_y_T] - v[idx_y_B])/(2*dx[1]);
                        
                        const double theta = dudx + dvdy;
                        const double Omega = fabs(dvdx - dudy);
                        
                        const double phi = theta*theta/(theta*theta + Omega*Omega + 1e-40);
                        
                        if (psi*phi > d_shock_Ducros_tol)
                        {
                            tag_ptr[idx_nghost] |= 1;
                        }
                    }
                }
            }
            else if (d_dim == tbox::Dimension(3))
            {
                // Get the arrays of conservative variables.
                double* rho   = density->getPointer(0);
                double* rho_u = momentum->getPointer(0);
                double* rho_v = momentum->getPointer(1);
                double* rho_w = momentum->getPointer(2);
                double* E     = total_energy->getPointer(0);
                
                // Get the array of pressure.
                double* p     = pressure->getPointer(0);
                double* u     = velocity->getPointer(0);
                double* v     = velocity->getPointer(1);
                double* w     = velocity->getPointer(2);
                
                for (int k = -d_num_ghosts[2]; k < interior_dims[2] + d_num_ghosts[2]; k++)
                {
                    for (int j = -d_num_ghosts[1]; j < interior_dims[1] + d_num_ghosts[1]; j++)
                    {
                        for (int i = -d_num_ghosts[0]; i < interior_dims[0] + d_num_ghosts[0]; i++)
                        {
                            // Compute index into linear data array.
                            const int idx = (i + d_num_ghosts[0]) +
                                (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            
                            std::vector<const double*> m_ptr;
                            m_ptr.push_back(&rho_u[idx]);
                            m_ptr.push_back(&rho_v[idx]);
                            m_ptr.push_back(&rho_w[idx]);
                            
                            p[idx] = d_equation_of_state->getPressure(
                                &rho[idx],
                                m_ptr,
                                &E[idx]);
                            
                            u[idx] = rho_u[idx]/rho[idx];
                            v[idx] = rho_v[idx]/rho[idx];
                            w[idx] = rho_w[idx]/rho[idx];
                        }
                    }
                }
                
                // Use the Ducros switch with threshold to tag cells.
                for (int k = 0; k < interior_dims[2]; k++)
                {
                    for (int j = 0; j < interior_dims[1]; j++)
                    {
                        for (int i = 0; i < interior_dims[0]; i++)
                        {
                            // Compute indices.
                            const int idx = (i + d_num_ghosts[0]) +
                                (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            
                            const int idx_x_L = (i - 1 + d_num_ghosts[0]) +
                                (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            
                            const int idx_x_R = (i + 1 + d_num_ghosts[0]) +
                                (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            
                            const int idx_y_B = (i + d_num_ghosts[0]) +
                                (j - 1 + d_num_ghosts[1])*ghostcell_dims[0] +
                                (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            
                            const int idx_y_T = (i + d_num_ghosts[0]) +
                                (j + 1 + d_num_ghosts[1])*ghostcell_dims[0] +
                                (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            
                            const int idx_z_B = (i + d_num_ghosts[0]) +
                                (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                (k - 1 + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            
                            const int idx_z_F = (i + d_num_ghosts[0]) +
                                (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                (k + 1 + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            
                            const int idx_nghost = i + j*interior_dims[0] + k*interior_dims[0]*interior_dims[1];
                            
                            // Compute Jameson's sensor.
                            
                            const double psi_x = fabs(p[idx_x_R] - 2*p[idx] + p[idx_x_L])/
                                (p[idx_x_R] + 2*p[idx] + p[idx_x_L] + EPSILON);
                            
                            const double psi_y = fabs(p[idx_y_T] - 2*p[idx] + p[idx_y_B])/
                                (p[idx_y_T] + 2*p[idx] + p[idx_y_B] + EPSILON);
                            
                            const double psi_z = fabs(p[idx_z_F] - 2*p[idx] + p[idx_z_B])/
                                (p[idx_z_F] + 2*p[idx] + p[idx_z_B] + EPSILON);
                            
                            const double psi = fmax(fmax(psi_x, psi_y), psi_z);
                            
                            // Compute another indicator.
                            
                            const double dudx = (u[idx_x_R] - u[idx_x_L])/(2*dx[0]);
                            const double dudy = (u[idx_y_T] - u[idx_y_B])/(2*dx[1]);
                            const double dudz = (u[idx_z_F] - u[idx_z_B])/(2*dx[2]);
                            
                            const double dvdx = (v[idx_x_R] - v[idx_x_L])/(2*dx[0]);
                            const double dvdy = (v[idx_y_T] - v[idx_y_B])/(2*dx[1]);
                            const double dvdz = (v[idx_z_F] - v[idx_z_B])/(2*dx[2]);
                            
                            const double dwdx = (w[idx_x_R] - w[idx_x_L])/(2*dx[0]);
                            const double dwdy = (w[idx_y_T] - w[idx_y_B])/(2*dx[1]);
                            const double dwdz = (w[idx_z_F] - w[idx_z_B])/(2*dx[2]);
                            
                            const double theta = dudx + dvdy + dwdz;;
                            
                            const double Omega = sqrt(pow(dwdy - dvdz, 2) +
                                pow(dudz - dwdx, 2) +
                                pow(dvdx - dudy, 2));
                            
                            const double phi = theta*theta/(theta*theta + Omega*Omega + 1e-40);
                            
                            if (psi*phi > d_shock_Ducros_tol)
                            {
                                tag_ptr[idx_nghost] |= 1;
                            }
                        }
                    }
                }
            }
        }
        else if (ref == "SHOCK_LARSSON")
        {
            // Get the cell data of the time-dependent variables.
            boost::shared_ptr<pdat::CellData<double> > density(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_density, data_context)));
            
            boost::shared_ptr<pdat::CellData<double> > momentum(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_momentum, data_context)));
            
            // Allocate temporary patch data.
            boost::shared_ptr<pdat::CellData<double> > velocity(
                new pdat::CellData<double>(interior_box, d_dim.getValue(), d_num_ghosts));
            
            if (d_dim == tbox::Dimension(1))
            {
                // NOT YET IMPLEMENTED
            }
            else if (d_dim == tbox::Dimension(2))
            {
                // Get the array of conservative variables.
                double* rho   = density->getPointer(0);
                double* rho_u = momentum->getPointer(0);
                double* rho_v = momentum->getPointer(1);
                
                // Get the arrays of temporary patch data.
                double* u     = velocity->getPointer(0);
                double* v     = velocity->getPointer(1);
                
                // Compute the field of velocities.
                for (int j = -d_num_ghosts[1]; j < interior_dims[1] + d_num_ghosts[1]; j++)
                {
                    for (int i = -d_num_ghosts[0]; i < interior_dims[0] + d_num_ghosts[0]; i++)
                    {
                        // Compute index into linear data array.
                        const int idx = (i + d_num_ghosts[0]) +
                            (j + d_num_ghosts[1])*ghostcell_dims[0];
                        
                        u[idx] = rho_u[idx]/rho[idx];
                        v[idx] = rho_v[idx]/rho[idx];
                    }
                }
                
                // Use the Larsson's switch with threshold to tag cells.
                for (int j = 0; j < interior_dims[1]; j++)
                {
                    for (int i = 0; i < interior_dims[0]; i++)
                    {
                        // Compute indices.
                        const int idx_x_L = (i - 1 + d_num_ghosts[0]) +
                            (j + d_num_ghosts[1])*ghostcell_dims[0];
                        
                        const int idx_x_R = (i + 1 + d_num_ghosts[0]) +
                            (j + d_num_ghosts[1])*ghostcell_dims[0];
                        
                        const int idx_y_B = (i + d_num_ghosts[0]) +
                            (j - 1 + d_num_ghosts[1])*ghostcell_dims[0];
                        
                        const int idx_y_T = (i + d_num_ghosts[0]) +
                            (j + 1 + d_num_ghosts[1])*ghostcell_dims[0];
                        
                        const int idx_nghost = i + j*interior_dims[0];
                        
                        const double dudx = (u[idx_x_R] - u[idx_x_L])/(2*dx[0]);
                        const double dudy = (u[idx_y_T] - u[idx_y_B])/(2*dx[1]);
                        
                        const double dvdx = (v[idx_x_R] - v[idx_x_L])/(2*dx[0]);
                        const double dvdy = (v[idx_y_T] - v[idx_y_B])/(2*dx[1]);
                        
                        const double theta = dudx + dvdy;
                        const double Omega = fabs(dvdx - dudy);
                        
                        const double s = -theta/(fabs(theta) + Omega + EPSILON);
                        
                        if (s > d_shock_Larsson_tol)
                        {
                            tag_ptr[idx_nghost] |= 1;
                        }
                    }
                }
            }
            else if (d_dim == tbox::Dimension(3))
            {
                // Get the arrays of conservative variables.
                double* rho   = density->getPointer(0);
                double* rho_u = momentum->getPointer(0);
                double* rho_v = momentum->getPointer(1);
                double* rho_w = momentum->getPointer(2);
                
                // Get the arrays of temporary patch data.
                double* u     = velocity->getPointer(0);
                double* v     = velocity->getPointer(1);
                double* w     = velocity->getPointer(2);
                
                // Compute the field of velocities.
                for (int k = -d_num_ghosts[2]; k < interior_dims[2] + d_num_ghosts[2]; k++)
                {
                    for (int j = -d_num_ghosts[1]; j < interior_dims[1] + d_num_ghosts[1]; j++)
                    {
                        for (int i = -d_num_ghosts[0]; i < interior_dims[0] + d_num_ghosts[0]; i++)
                        {
                            // Compute index into linear data array.
                            const int idx = (i + d_num_ghosts[0]) +
                                (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            
                            u[idx] = rho_u[idx]/rho[idx];
                            v[idx] = rho_v[idx]/rho[idx];
                            w[idx] = rho_w[idx]/rho[idx];
                        }
                    }
                }
                
                // Use the Larsson's switch with threshold to tag cells.
                for (int k = 0; k < interior_dims[2]; k++)
                {
                    for (int j = 0; j < interior_dims[1]; j++)
                    {
                        for (int i = 0; i < interior_dims[0]; i++)
                        {
                            // Compute indices.
                            const int idx_x_L = (i - 1 + d_num_ghosts[0]) +
                                (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            
                            const int idx_x_R = (i + 1 + d_num_ghosts[0]) +
                                (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            
                            const int idx_y_B = (i + d_num_ghosts[0]) +
                                (j - 1 + d_num_ghosts[1])*ghostcell_dims[0] +
                                (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            
                            const int idx_y_T = (i + d_num_ghosts[0]) +
                                (j + 1 + d_num_ghosts[1])*ghostcell_dims[0] +
                                (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            
                            const int idx_z_B = (i + d_num_ghosts[0]) +
                                (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                (k - 1 + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            
                            const int idx_z_F = (i + d_num_ghosts[0]) +
                                (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                (k + 1 + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            
                            const int idx_nghost = i + j*interior_dims[0] + k*interior_dims[0]*interior_dims[1];
                            
                            const double dudx = (u[idx_x_R] - u[idx_x_L])/(2*dx[0]);
                            const double dudy = (u[idx_y_T] - u[idx_y_B])/(2*dx[1]);
                            const double dudz = (u[idx_z_F] - u[idx_z_B])/(2*dx[2]);
                            
                            const double dvdx = (v[idx_x_R] - v[idx_x_L])/(2*dx[0]);
                            const double dvdy = (v[idx_y_T] - v[idx_y_B])/(2*dx[1]);
                            const double dvdz = (v[idx_z_F] - v[idx_z_B])/(2*dx[2]);
                            
                            const double dwdx = (w[idx_x_R] - w[idx_x_L])/(2*dx[0]);
                            const double dwdy = (w[idx_y_T] - w[idx_y_B])/(2*dx[1]);
                            const double dwdz = (w[idx_z_F] - w[idx_z_B])/(2*dx[2]);
                            
                            const double theta = dudx + dvdy + dwdz;;
                            
                            const double Omega = sqrt(pow(dwdy - dvdz, 2) +
                                pow(dudz - dwdx, 2) +
                                pow(dvdx - dudy, 2));
                            
                            const double s = -theta/(fabs(theta) + Omega + EPSILON);
                            
                            if (s > d_shock_Larsson_tol)
                            {
                                tag_ptr[idx_nghost] |= 1;
                            }
                        }
                    }
                }
            }
        }
        else if (ref == "DENSITY_JAMESON")
        {
            switch (d_flow_model)
            {
                case SINGLE_SPECIES:
                {
                    // Get the cell data of the time-dependent variables.
                    boost::shared_ptr<pdat::CellData<double> > density(
                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                            patch.getPatchData(d_density, data_context)));
                    
                    if (d_dim == tbox::Dimension(1))
                    {
                        // NOT YET IMPLEMENTED
                    }
                    else if (d_dim == tbox::Dimension(2))
                    {
                        // Get the arrays of density.
                        double* rho   = density->getPointer(0);
                        
                        // Use the Jameson's switch with threshold to tag cells.
                        for (int j = 0; j < interior_dims[1]; j++)
                        {
                            for (int i = 0; i < interior_dims[0]; i++)
                            {
                                // Compute indices.
                                const int idx = (i + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0];
                                
                                const int idx_x_L = (i - 1 + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0];
                                
                                const int idx_x_R = (i + 1 + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0];
                                
                                const int idx_y_B = (i + d_num_ghosts[0]) +
                                    (j - 1 + d_num_ghosts[1])*ghostcell_dims[0];
                                
                                const int idx_y_T = (i + d_num_ghosts[0]) +
                                    (j + 1 + d_num_ghosts[1])*ghostcell_dims[0];
                                
                                const int idx_nghost = i + j*interior_dims[0];
                                
                                const double psi_x = fabs(rho[idx_x_R] - 2*rho[idx] + rho[idx_x_L])/
                                    (rho[idx_x_R] + 2*rho[idx] + rho[idx_x_L] + EPSILON);
                                
                                const double psi_y = fabs(rho[idx_y_T] - 2*rho[idx] + rho[idx_y_B])/
                                    (rho[idx_y_T] + 2*rho[idx] + rho[idx_y_B] + EPSILON);
                                
                                const double psi = fmax(psi_x, psi_y);
                                
                                if (psi > d_density_Jameson_tol)
                                {
                                    tag_ptr[idx_nghost] |= 1;
                                }
                            }
                        }
                    }
                    else if (d_dim == tbox::Dimension(3))
                    {
                        // Get the arrays of density.
                        double* rho   = density->getPointer(0);
                        
                        // Use the Jameson's switch with threshold to tag cells.
                        for (int k = 0; k < interior_dims[2]; k++)
                        {
                            for (int j = 0; j < interior_dims[1]; j++)
                            {
                                for (int i = 0; i < interior_dims[0]; i++)
                                {
                                    // Compute indices.
                                    const int idx = (i + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const int idx_x_L = (i - 1 + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const int idx_x_R = (i + 1 + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const int idx_y_B = (i + d_num_ghosts[0]) +
                                        (j - 1 + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const int idx_y_T = (i + d_num_ghosts[0]) +
                                        (j + 1 + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const int idx_z_B = (i + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k - 1 + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const int idx_z_F = (i + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + 1 + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const int idx_nghost = i + j*interior_dims[0] + k*interior_dims[0]*interior_dims[1];
                                    
                                    const double psi_x = fabs(rho[idx_x_R] - 2*rho[idx] + rho[idx_x_L])/
                                        (rho[idx_x_R] + 2*rho[idx] + rho[idx_x_L] + EPSILON);
                                    
                                    const double psi_y = fabs(rho[idx_y_T] - 2*rho[idx] + rho[idx_y_B])/
                                        (rho[idx_y_T] + 2*rho[idx] + rho[idx_y_B] + EPSILON);
                                    
                                    const double psi_z = fabs(rho[idx_z_F] - 2*rho[idx] + rho[idx_z_B])/
                                        (rho[idx_z_F] + 2*rho[idx] + rho[idx_z_B] + EPSILON);
                                    
                                    const double psi = fmax(fmax(psi_x, psi_y), psi_z);
                                    
                                    if (psi > d_density_Jameson_tol)
                                    {
                                        tag_ptr[idx_nghost] |= 1;
                                    }
                                }
                            }
                        }
                    }
                    
                    break;
                }
                case FOUR_EQN_SHYUE:
                {
                    // Get the cell data of the time-dependent variables.
                    boost::shared_ptr<pdat::CellData<double> > density(
                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                            patch.getPatchData(d_density, data_context)));
                    
                    if (d_dim == tbox::Dimension(1))
                    {
                        // NOT YET IMPLEMENTED
                    }
                    else if (d_dim == tbox::Dimension(2))
                    {
                        // Get the arrays of density.
                        double* rho   = density->getPointer(0);
                        
                        // Use the Jameson's switch with threshold to tag cells.
                        for (int j = 0; j < interior_dims[1]; j++)
                        {
                            for (int i = 0; i < interior_dims[0]; i++)
                            {
                                // Compute indices.
                                const int idx = (i + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0];
                                
                                const int idx_x_L = (i - 1 + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0];
                                
                                const int idx_x_R = (i + 1 + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0];
                                
                                const int idx_y_B = (i + d_num_ghosts[0]) +
                                    (j - 1 + d_num_ghosts[1])*ghostcell_dims[0];
                                
                                const int idx_y_T = (i + d_num_ghosts[0]) +
                                    (j + 1 + d_num_ghosts[1])*ghostcell_dims[0];
                                
                                const int idx_nghost = i + j*interior_dims[0];
                                
                                const double psi_x = fabs(rho[idx_x_R] - 2*rho[idx] + rho[idx_x_L])/
                                    (rho[idx_x_R] + 2*rho[idx] + rho[idx_x_L] + EPSILON);
                                
                                const double psi_y = fabs(rho[idx_y_T] - 2*rho[idx] + rho[idx_y_B])/
                                    (rho[idx_y_T] + 2*rho[idx] + rho[idx_y_B] + EPSILON);
                                
                                const double psi = fmax(psi_x, psi_y);
                                
                                if (psi > d_density_Jameson_tol)
                                {
                                    tag_ptr[idx_nghost] |= 1;
                                }
                            }
                        }
                    }
                    else if (d_dim == tbox::Dimension(3))
                    {
                        // Get the arrays of density.
                        double* rho   = density->getPointer(0);
                        
                        // Use the Jameson's switch with threshold to tag cells.
                        for (int k = 0; k < interior_dims[2]; k++)
                        {
                            for (int j = 0; j < interior_dims[1]; j++)
                            {
                                for (int i = 0; i < interior_dims[0]; i++)
                                {
                                    // Compute indices.
                                    const int idx = (i + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const int idx_x_L = (i - 1 + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const int idx_x_R = (i + 1 + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const int idx_y_B = (i + d_num_ghosts[0]) +
                                        (j - 1 + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const int idx_y_T = (i + d_num_ghosts[0]) +
                                        (j + 1 + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const int idx_z_B = (i + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k - 1 + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const int idx_z_F = (i + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + 1 + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const int idx_nghost = i + j*interior_dims[0] + k*interior_dims[0]*interior_dims[1];
                                    
                                    const double psi_x = fabs(rho[idx_x_R] - 2*rho[idx] + rho[idx_x_L])/
                                        (rho[idx_x_R] + 2*rho[idx] + rho[idx_x_L] + EPSILON);
                                    
                                    const double psi_y = fabs(rho[idx_y_T] - 2*rho[idx] + rho[idx_y_B])/
                                        (rho[idx_y_T] + 2*rho[idx] + rho[idx_y_B] + EPSILON);
                                    
                                    const double psi_z = fabs(rho[idx_z_F] - 2*rho[idx] + rho[idx_z_B])/
                                        (rho[idx_z_F] + 2*rho[idx] + rho[idx_z_B] + EPSILON);
                                    
                                    const double psi = fmax(fmax(psi_x, psi_y), psi_z);
                                    
                                    if (psi > d_density_Jameson_tol)
                                    {
                                        tag_ptr[idx_nghost] |= 1;
                                    }
                                }
                            }
                        }
                    }
                    
                    break;
                }
                case FIVE_EQN_ALLAIRE:
                {
                    // Get the cell data of the time-dependent variables.
                    boost::shared_ptr<pdat::CellData<double> > partial_density(
                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                            patch.getPatchData(d_partial_density, data_context)));
                    
                    // Allocate temporary patch data.
                    boost::shared_ptr<pdat::CellData<double> > density(
                        new pdat::CellData<double>(interior_box, 1, d_num_ghosts));
                    
                    if (d_dim == tbox::Dimension(1))
                    {
                        // NOT YET IMPLEMENTED
                    }
                    else if (d_dim == tbox::Dimension(2))
                    {
                        // Get the arrays of conservative variables.
                        std::vector<double*> Z_rho;
                        for (int si = 0; si < d_num_species; si++)
                        {
                            Z_rho.push_back(partial_density->getPointer(si));
                        }
                        
                        // Get the array of pressure.
                        double* rho     = density->getPointer(0);
                        
                        // Compute the field of density.
                        for (int j = -d_num_ghosts[1]; j < interior_dims[1] + d_num_ghosts[1]; j++)
                        {
                            for (int i = -d_num_ghosts[0]; i < interior_dims[0] + d_num_ghosts[0]; i++)
                            {
                                // Compute index into linear data array.
                                const int idx = (i + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0];
                                
                                std::vector<const double*> Z_rho_ptr;
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Z_rho_ptr.push_back(&Z_rho[si][idx]);
                                }
                                
                                rho[idx] = d_equation_of_state->
                                    getTotalDensity(
                                        Z_rho_ptr);
                            }
                        }
                        
                        // Use the Jameson's switch with threshold to tag cells.
                        for (int j = 0; j < interior_dims[1]; j++)
                        {
                            for (int i = 0; i < interior_dims[0]; i++)
                            {
                                // Compute indices.
                                const int idx = (i + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0];
                                
                                const int idx_x_L = (i - 1 + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0];
                                
                                const int idx_x_R = (i + 1 + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0];
                                
                                const int idx_y_B = (i + d_num_ghosts[0]) +
                                    (j - 1 + d_num_ghosts[1])*ghostcell_dims[0];
                                
                                const int idx_y_T = (i + d_num_ghosts[0]) +
                                    (j + 1 + d_num_ghosts[1])*ghostcell_dims[0];
                                
                                const int idx_nghost = i + j*interior_dims[0];
                                
                                const double psi_x = fabs(rho[idx_x_R] - 2*rho[idx] + rho[idx_x_L])/
                                    (rho[idx_x_R] + 2*rho[idx] + rho[idx_x_L] + EPSILON);
                                
                                const double psi_y = fabs(rho[idx_y_T] - 2*rho[idx] + rho[idx_y_B])/
                                    (rho[idx_y_T] + 2*rho[idx] + rho[idx_y_B] + EPSILON);
                                
                                const double psi = fmax(psi_x, psi_y);
                                
                                if (psi > d_density_Jameson_tol)
                                {
                                    tag_ptr[idx_nghost] |= 1;
                                }
                            }
                        }
                    }
                    else if (d_dim == tbox::Dimension(3))
                    {
                        // Get the arrays of conservative variables.
                        std::vector<double*> Z_rho;
                        for (int si = 0; si < d_num_species; si++)
                        {
                            Z_rho.push_back(partial_density->getPointer(si));
                        }
                        
                        // Get the array of pressure.
                        double* rho     = density->getPointer(0);
                        
                        // Compute the field of density.
                        for (int k = -d_num_ghosts[2]; k < interior_dims[2] + d_num_ghosts[2]; k++)
                        {
                            for (int j = -d_num_ghosts[1]; j < interior_dims[1] + d_num_ghosts[1]; j++)
                            {
                                for (int i = -d_num_ghosts[0]; i < interior_dims[0] + d_num_ghosts[0]; i++)
                                {
                                    // Compute index into linear data array.
                                    const int idx = (i + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    std::vector<const double*> Z_rho_ptr;
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        Z_rho_ptr.push_back(&Z_rho[si][idx]);
                                    }
                                    
                                    rho[idx] = d_equation_of_state->
                                        getTotalDensity(
                                            Z_rho_ptr);
                                }
                            }
                        }
                        
                        // Use the Jameson's switch with threshold to tag cells.
                        for (int k = 0; k < interior_dims[2]; k++)
                        {
                            for (int j = 0; j < interior_dims[1]; j++)
                            {
                                for (int i = 0; i < interior_dims[0]; i++)
                                {
                                    // Compute indices.
                                    const int idx = (i + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const int idx_x_L = (i - 1 + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const int idx_x_R = (i + 1 + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const int idx_y_B = (i + d_num_ghosts[0]) +
                                        (j - 1 + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const int idx_y_T = (i + d_num_ghosts[0]) +
                                        (j + 1 + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const int idx_z_B = (i + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k - 1 + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const int idx_z_F = (i + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + 1 + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const int idx_nghost = i + j*interior_dims[0] + k*interior_dims[0]*interior_dims[1];
                                    
                                    const double psi_x = fabs(rho[idx_x_R] - 2*rho[idx] + rho[idx_x_L])/
                                        (rho[idx_x_R] + 2*rho[idx] + rho[idx_x_L] + EPSILON);
                                    
                                    const double psi_y = fabs(rho[idx_y_T] - 2*rho[idx] + rho[idx_y_B])/
                                        (rho[idx_y_T] + 2*rho[idx] + rho[idx_y_B] + EPSILON);
                                    
                                    const double psi_z = fabs(rho[idx_z_F] - 2*rho[idx] + rho[idx_z_B])/
                                        (rho[idx_z_F] + 2*rho[idx] + rho[idx_z_B] + EPSILON);
                                    
                                    const double psi = fmax(fmax(psi_x, psi_y), psi_z);
                                    
                                    if (psi > d_density_Jameson_tol)
                                    {
                                        tag_ptr[idx_nghost] |= 1;
                                    }
                                }
                            }
                        }
                    }
                    
                    break;
                }
                default:
                {
                    TBOX_ERROR(d_object_name
                        << ": "
                        << "d_flow_model '"
                        << d_flow_model
                        << "' not yet implemented."
                        << std::endl);
                }
            }
        }
    }
}
