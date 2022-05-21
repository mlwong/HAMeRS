#include "flow/flow_models/single-species/FlowModelSourceUtilitiesSingleSpecies.hpp"

FlowModelSourceUtilitiesSingleSpecies::FlowModelSourceUtilitiesSingleSpecies(
    const std::string& object_name,
    const std::string& project_name,
    const tbox::Dimension& dim,
    const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
    const int& num_species,
    const HAMERS_SHARED_PTR<tbox::Database>& flow_model_db,
    const HAMERS_SHARED_PTR<EquationOfStateMixingRules> equation_of_state_mixing_rules):
        FlowModelSourceUtilities(
            object_name,
            project_name,
            dim,
            grid_geometry,
            num_species,
            2 + dim.getValue(),
            flow_model_db),
    d_has_gravity(false),
    d_num_subghosts_gruneisen_parameter(-hier::IntVector::getOne(d_dim)),
    d_subghost_box_gruneisen_parameter(hier::Box::getEmptyBox(d_dim)),
    d_subghostcell_dims_gruneisen_parameter(hier::IntVector::getZero(d_dim)),
    d_cell_data_computed_gruneisen_parameter(false),
    d_equation_of_state_mixing_rules(equation_of_state_mixing_rules)
{
    if (d_has_source_terms)
    {
        HAMERS_SHARED_PTR<tbox::Database> source_terms_db;
        
        if (flow_model_db->keyExists("Source_terms"))
        {
            source_terms_db = flow_model_db->getDatabase("Source_terms");
        }
        else if (flow_model_db->keyExists("d_source_terms"))
        {
            source_terms_db = flow_model_db->getDatabase("d_source_terms");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "No key 'Source_terms'/'d_source_terms' found in data for flow model."
                << std::endl);
        }
        
        if (source_terms_db->keyExists("has_gravity"))
        {
            d_has_gravity = source_terms_db->getBool("has_gravity");
            if (d_has_gravity)
            {
                if (source_terms_db->keyExists("gravity"))
                {
                    source_terms_db->getVector("gravity", d_gravity);
                }
                else
                {
                    TBOX_ERROR(d_object_name
                        << ": "
                        << "No key 'gravity' found in data for source terms."
                        << std::endl);
                }
            }
        }
        else if (source_terms_db->keyExists("d_has_gravity"))
        {
            d_has_gravity = source_terms_db->getBool("d_has_gravity");
            if (d_has_gravity)
            {
                if (source_terms_db->keyExists("d_gravity"))
                {
                    source_terms_db->getVector("d_gravity", d_gravity);
                }
                else
                {
                    TBOX_ERROR(d_object_name
                        << ": "
                        << "No key 'd_gravity' found in data for source terms."
                        << std::endl);
                }
            }
        }
        
        if (d_has_gravity)
        {
            if (d_gravity.size() != d_dim.getValue())
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSourceUtilitiesSingleSpecies::FlowModelSourceUtilitiesSingleSpecies\n"
                    << "Size of 'gravity' or 'd_gravity' is not consistent with problem dimension."
                    << std::endl);
            }
        }
    }
}


/*
 * Register the required variables for the computation of source terms in the registered patch.
 */
void
FlowModelSourceUtilitiesSingleSpecies::registerDerivedVariablesForSourceTerms(
    const hier::IntVector& num_subghosts)
{
    if (d_has_source_terms)
    {
        if (d_flow_model.expired())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "The object is not setup yet!"
                << std::endl);
        }
        
        HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
        
        // Check whether a patch is already registered.
        if (!flow_model_tmp->hasRegisteredPatch())
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSourceUtilitiesSingleSpecies::"
                << "registerDerivedVariablesForSourceTerms()\n"
                << "No patch is registered yet."
                << std::endl);
        }
        
        // Check whether all or part of derived cell data is alredy computed.
        if (d_derived_cell_data_computed)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSourceUtilitiesSingleSpecies::registerDerivedVariablesForSourceTerms()\n"
                << "Derived cell data is already computed."
                << std::endl);
        }
        
        if ((num_subghosts < hier::IntVector::getZero(d_dim)) ||
            (num_subghosts > flow_model_tmp->getNumberOfGhostCells()))
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSourceUtilitiesSingleSpecies::registerDerivedVariablesForSourceTerms()\n"
                << "The number of sub-ghost cells of variable is not between zero and number of ghosts of"
                << " conservative variables."
                << std::endl);
        }
        
        setNumberOfSubGhosts(
            num_subghosts,
            "SOURCE_TERMS",
            "SOURCE_TERMS");
        
        if (d_has_gravity)
        {
            /*
             * Register the required derived variables in flow model.
             */
            
            // Nothing needs to be done.
        }
    }
}


/*
 * Register the required variables for the computation of local stable time increment for
 * source terms in the registered patch.
 */
void
FlowModelSourceUtilitiesSingleSpecies::registerDerivedVariablesForSourceTermsStableDt(
    const hier::IntVector& num_subghosts)
{
    if (d_has_source_terms)
    {
        if (d_flow_model.expired())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "The object is not setup yet!"
                << std::endl);
        }
        
        HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
        
        // Check whether a patch is already registered.
        if (!flow_model_tmp->hasRegisteredPatch())
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSourceUtilitiesSingleSpecies::"
                << "registerDerivedVariablesForSourceTermsStableDt()\n"
                << "No patch is registered yet."
                << std::endl);
        }
        
        // Check whether all or part of derived cell data is alredy computed.
        if (d_derived_cell_data_computed)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSourceUtilitiesSingleSpecies::registerDerivedVariablesForSourceTermsStableDt()\n"
                << "Derived cell data is already computed."
                << std::endl);
        }
        
        if ((num_subghosts < hier::IntVector::getZero(d_dim)) ||
            (num_subghosts > flow_model_tmp->getNumberOfGhostCells()))
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSourceUtilitiesSingleSpecies::registerDerivedVariablesForSourceTermsStableDt()\n"
                << "The number of sub-ghost cells of variable is not between zero and number of ghosts of"
                << " conservative variables."
                << std::endl);
        }
        
        setNumberOfSubGhosts(
            num_subghosts,
            "SOURCE_TERMS",
            "SOURCE_TERMS");
        
        if (d_has_gravity)
        {
            /*
             * Register the required derived variables in flow model.
             */
            
            std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
            
            num_subghosts_of_data.insert(
                std::pair<std::string, hier::IntVector>("PRESSURE", num_subghosts));
            
            num_subghosts_of_data.insert(
                std::pair<std::string, hier::IntVector>("SOUND_SPEED", num_subghosts));
            
            flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
            
            setNumberOfSubGhosts(
                num_subghosts,
                "GRUNEISEN_PARAMETER",
                "GRUNEISEN_PARAMETER");
        }
    }
}


/*
 * Allocate memory for cell data of different registered derived variables related to this
 * class in the registered patch.
 */
void
FlowModelSourceUtilitiesSingleSpecies::allocateMemoryForDerivedCellData()
{
    if (d_has_source_terms)
    {
        if (d_flow_model.expired())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "The object is not setup yet!"
                << std::endl);
        }
        
        HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
        
        flow_model_tmp->allocateMemoryForDerivedCellData();
        
        const hier::Patch& patch = flow_model_tmp->getRegisteredPatch();
        const hier::Box interior_box = patch.getBox();
        
        if (d_has_gravity)
        {
            if (!d_cell_data_computed_gruneisen_parameter)
            {
                if (!d_data_gruneisen_parameter)
                {
                    // Create the cell data of Gruneisen parameter.
                    d_data_gruneisen_parameter.reset(new pdat::CellData<double>(
                        interior_box, 1, d_subghostcell_dims_gruneisen_parameter));
                }
            }
            else
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSourceUtilitiesSingleSpecies::allocateMemoryForDerivedCellData()\n"
                    << "Cell data of 'GRUNEISEN_PARAMETER' is aleady computed."
                    << std::endl);
            }
        }
    }
}


/*
 * Clear cell data of different derived variables related to this class in the registered patch.
 */
void
FlowModelSourceUtilitiesSingleSpecies::clearCellData()
{
    d_num_subghosts_source_terms        = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_gruneisen_parameter = -hier::IntVector::getOne(d_dim);
    
    d_subghost_box_source_terms        = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_gruneisen_parameter = hier::Box::getEmptyBox(d_dim);
    
    d_subghostcell_dims_source_terms        = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_gruneisen_parameter = hier::IntVector::getZero(d_dim);
    
    d_data_gruneisen_parameter.reset();
    
    d_cell_data_computed_gruneisen_parameter = false;
    
    d_derived_cell_data_computed = false;
}


/*
 * Compute cell data of different registered derived variables related to this class.
 */
void
FlowModelSourceUtilitiesSingleSpecies::computeDerivedCellData()
{
    if (d_has_source_terms)
    {
        if (d_flow_model.expired())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "The object is not setup yet!"
                << std::endl);
        }
        
        HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
        
        // Check whether a patch is already registered.
        if (!flow_model_tmp->hasRegisteredPatch())
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSourceUtilitiesSingleSpecies::"
                << "computeDerivedCellData()\n"
                << "No patch is registered yet."
                << std::endl);
        }
        
        flow_model_tmp->computeDerivedCellData();
        
        /*
         * Set the boxes and their dimensions for the derived cell variables.
         */
        if (!d_derived_cell_data_computed)
        {
            setDerivedCellVariableGhostBoxes();
        }
        
        // Compute the Gruneisen parameter cell data.
        if (d_num_subghosts_gruneisen_parameter > -hier::IntVector::getOne(d_dim))
        {
            if (!d_cell_data_computed_gruneisen_parameter)
            {
                computeCellDataOfGruneisenParameter();
            }
        }
        
        d_derived_cell_data_computed = true;
    }
}


/*
 * Compute all source terms on a patch.
 */
void
FlowModelSourceUtilitiesSingleSpecies::computeSourceTermsOnPatch(
    const HAMERS_SHARED_PTR<pdat::CellVariable<double> >& variable_source,
    const double time,
    const double dt,
    const int RK_step_number)
{
    NULL_USE(time);
    NULL_USE(RK_step_number);
    
    if (d_has_source_terms)
    {
        if (d_flow_model.expired())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "The object is not setup yet!"
                << std::endl);
        }
        
        HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
        const hier::Patch& patch = flow_model_tmp->getRegisteredPatch();
        const HAMERS_SHARED_PTR<hier::VariableContext> data_context = flow_model_tmp->getDataContext();
        
        // Get the cell data of source.
        HAMERS_SHARED_PTR<pdat::CellData<double> > source(
            HAMERS_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                patch.getPatchData(variable_source, data_context)));
        
        std::vector<double*> S;
        S.reserve(d_num_eqn);
        for (int si = 0; si < d_num_eqn; si++)
        {
            S.push_back(source->getPointer(si));
        }
        
        /*
         * Get the dimension of the interior box.
         */
        
        const hier::Box interior_box = patch.getBox();
        const hier::IntVector interior_dims = interior_box.numberCells();
        
        if (d_has_gravity)
        {
            // Get the cell data of momentum.
            HAMERS_SHARED_PTR<pdat::CellData<double> > data_momentum =
                flow_model_tmp->getCellData("MOMENTUM");
            
            // Get the cell data of density.
            HAMERS_SHARED_PTR<pdat::CellData<double> > data_density =
                flow_model_tmp->getCellData("DENSITY");
            
            /*
             * Get the numbers of ghost cells of momentum and density.
             */
            
            const hier::IntVector num_subghosts_momentum = data_momentum->getGhostCellWidth();
            const hier::IntVector num_subghosts_density  = data_density->getGhostCellWidth();
            
            /*
             * Get the dimensions of the ghost cell box of momentum and density.
             */
            
            const hier::Box subghost_box_momentum = data_momentum->getGhostBox();
            const hier::IntVector subghostcell_dims_momentum = subghost_box_momentum.numberCells();
            
            const hier::Box subghost_box_density = data_density->getGhostBox();
            const hier::IntVector subghostcell_dims_density = subghost_box_density.numberCells();
            
            //Get the pointer to the cell data of density.
            double* rho = data_density->getPointer(0);
            
            if (d_dim == tbox::Dimension(1))
            {
                // Get the pointer to cell data of momentum.
                double* rho_u = data_momentum->getPointer(0);
                
                /*
                 * Compute the source due to gravity.
                 */
                for (int i = -d_num_subghosts_source_terms[0];
                     i < interior_dims[0] + d_num_subghosts_source_terms[0];
                     i++)
                {
                    // Compute the linear indices.
                    const int idx_source   = i + d_num_subghosts_source_terms[0];
                    const int idx_momentum = i + num_subghosts_momentum[0];
                    const int idx_density  = i + num_subghosts_density[0];
                    
                    S[1][idx_source] += dt*rho[idx_density]*d_gravity[0];
                    S[2][idx_source] += dt*rho_u[idx_momentum]*d_gravity[0];
                }
            }
            else if (d_dim == tbox::Dimension(2))
            {
                // Get the pointer to cell data of momentum.
                double* rho_u = data_momentum->getPointer(0);
                double* rho_v = data_momentum->getPointer(1);
                
                /*
                 * Compute the source due to gravity.
                 */
                for (int j = -d_num_subghosts_source_terms[1];
                     j < interior_dims[1] + d_num_subghosts_source_terms[1];
                     j++)
                {
                    for (int i = -d_num_subghosts_source_terms[0];
                         i < interior_dims[0] + d_num_subghosts_source_terms[0];
                         i++)
                    {
                        // Compute the linear indices.
                        const int idx_source = (i + d_num_subghosts_source_terms[0]) +
                            (j + d_num_subghosts_source_terms[1])*d_subghostcell_dims_source_terms[0];
                        
                        const int idx_momentum = (i + num_subghosts_momentum[0]) +
                            (j + num_subghosts_momentum[1])*subghostcell_dims_momentum[0];
                        
                        const int idx_density = (i + num_subghosts_density[0]) +
                            (j + num_subghosts_density[1])*subghostcell_dims_density[0];
                        
                        S[1][idx_source] += dt*rho[idx_density]*d_gravity[0];
                        S[2][idx_source] += dt*rho[idx_density]*d_gravity[1];
                        S[3][idx_source] += dt*
                            (rho_u[idx_momentum]*d_gravity[0] +
                             rho_v[idx_momentum]*d_gravity[1]);
                    }
                }
            }
            else if (d_dim == tbox::Dimension(3))
            {
                // Get the pointer to cell data of momentum.
                double* rho_u = data_momentum->getPointer(0);
                double* rho_v = data_momentum->getPointer(1);
                double* rho_w = data_momentum->getPointer(2);
                
                /*
                 * Compute the source due to gravity.
                 */
                for (int k = -d_num_subghosts_source_terms[2];
                     k < interior_dims[2] + d_num_subghosts_source_terms[2];
                     k++)
                {
                    for (int j = -d_num_subghosts_source_terms[1];
                         j < interior_dims[1] + d_num_subghosts_source_terms[1];
                         j++)
                    {
                        for (int i = -d_num_subghosts_source_terms[0];
                             i < interior_dims[0] + d_num_subghosts_source_terms[0];
                             i++)
                        {
                            const int idx_source = (i + d_num_subghosts_source_terms[0]) +
                                (j + d_num_subghosts_source_terms[1])*d_subghostcell_dims_source_terms[0] +
                                (k + d_num_subghosts_source_terms[2])*d_subghostcell_dims_source_terms[0]*
                                    d_subghostcell_dims_source_terms[1];
                            
                            const int idx_momentum = (i + num_subghosts_momentum[0]) +
                                (j + num_subghosts_momentum[1])*subghostcell_dims_momentum[0] +
                                (k + num_subghosts_momentum[2])*subghostcell_dims_momentum[0]*
                                    subghostcell_dims_momentum[1];
                            
                            const int idx_density = (i + num_subghosts_density[0]) +
                                (j + num_subghosts_density[1])*subghostcell_dims_density[0] +
                                (k + num_subghosts_density[2])*subghostcell_dims_density[0]*
                                    subghostcell_dims_density[1];
                            
                            S[1][idx_source] += dt*rho[idx_density]*d_gravity[0];
                            S[2][idx_source] += dt*rho[idx_density]*d_gravity[1];
                            S[3][idx_source] += dt*rho[idx_density]*d_gravity[2];
                            S[4][idx_source] += dt*
                                (rho_u[idx_momentum]*d_gravity[0] +
                                 rho_v[idx_momentum]*d_gravity[1] +
                                 rho_w[idx_momentum]*d_gravity[2]);
                        }
                    }
                }
            }
        }
        
        flow_model_tmp.reset();
        computeSpecialSourceTermsOnPatch(
            variable_source,
            time,
            dt,
            RK_step_number);
    }
}


/*
 * Get local stable time increment for source terms.
 */
double
FlowModelSourceUtilitiesSingleSpecies::getStableDtOnPatch()
{
    // Create empty box.
    const hier::Box empty_box(d_dim);
    
    double dt_stable_global = tbox::MathUtilities<double>::getMax();
    
    if (d_has_source_terms)
    {
        if (d_flow_model.expired())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "The object is not setup yet!"
                << std::endl);
        }
        
        HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
        const hier::Patch& patch = flow_model_tmp->getRegisteredPatch();
        const HAMERS_SHARED_PTR<hier::VariableContext> data_context = flow_model_tmp->getDataContext();
        
        /*
         * Get the dimension of the interior box.
         */
        
        const hier::Box interior_box = patch.getBox();
        const hier::IntVector interior_dims = interior_box.numberCells();
        
        if (d_has_gravity)
        {
            if (!d_cell_data_computed_gruneisen_parameter)
            {
                computeCellDataOfGruneisenParameter();
            }
            
            // Get the cell data of sound speed.
            HAMERS_SHARED_PTR<pdat::CellData<double> > data_sound_speed =
                flow_model_tmp->getCellData("SOUND_SPEED");
            
            /*
             * Get the numbers of ghost cells of Gruneisen parameter and sound speed.
             */
            
            const hier::IntVector num_subghosts_gruneisen_parameter = d_data_gruneisen_parameter->getGhostCellWidth();
            const hier::IntVector num_subghosts_sound_speed = data_sound_speed->getGhostCellWidth();
            
            /*
             * Get the dimensions of the ghost cell box of Gruneisen parameter and sound speed.
             */
            
            const hier::Box subghost_box_gruneisen_parameter = d_data_gruneisen_parameter->getGhostBox();
            const hier::IntVector subghostcell_dims_gruneisen_parameter = subghost_box_gruneisen_parameter.numberCells();
            
            const hier::Box subghost_box_sound_speed = data_sound_speed->getGhostBox();
            const hier::IntVector subghostcell_dims_sound_speed = subghost_box_sound_speed.numberCells();
            
            // Get the pointer to the cell data of Gruneisen parameter.
            double* Gamma = d_data_gruneisen_parameter->getPointer(0);
            
            // Get the pointer to the cell data of sound speed.
            double* c = data_sound_speed->getPointer(0);
            
            if (d_dim == tbox::Dimension(1))
            {
                const double g_mag = fabs(d_gravity[0]);
                
                /*
                 * Compute the stable time step size.
                 */
                for (int i = -d_num_subghosts_source_terms[0];
                     i < interior_dims[0] + d_num_subghosts_source_terms[0];
                     i++)
                {
                    // Compute the linear indices.
                    const int idx_gruneisen_parameter = i + d_num_subghosts_gruneisen_parameter[0];
                    const int idx_sound_speed         = i + num_subghosts_sound_speed[0];
                    
                    const double dt_stable_local = c[idx_sound_speed]/
                        (g_mag*sqrt(double(2)*Gamma[idx_gruneisen_parameter]*(Gamma[idx_gruneisen_parameter] + double(1))));
                    
                    dt_stable_global = fmin(dt_stable_local, dt_stable_global);
                }
            }
            else if (d_dim == tbox::Dimension(2))
            {
                const double g_mag = sqrt(d_gravity[0]*d_gravity[0] + d_gravity[1]*d_gravity[1]);
                
                /*
                 * Compute the stable time step size.
                 */
                for (int j = -d_num_subghosts_source_terms[1];
                     j < interior_dims[1] + d_num_subghosts_source_terms[1];
                     j++)
                {
                    for (int i = -d_num_subghosts_source_terms[0];
                         i < interior_dims[0] + d_num_subghosts_source_terms[0];
                         i++)
                    {
                        // Compute the linear indices.
                        const int idx_gruneisen_parameter = (i + num_subghosts_gruneisen_parameter[0]) +
                            (j + num_subghosts_gruneisen_parameter[1])*subghostcell_dims_gruneisen_parameter[0];
                        
                        const int idx_sound_speed = (i + num_subghosts_sound_speed[0]) +
                            (j + num_subghosts_sound_speed[1])*subghostcell_dims_sound_speed[0];
                        
                        const double dt_stable_local = c[idx_sound_speed]/
                            (g_mag*sqrt(double(2)*Gamma[idx_gruneisen_parameter]*(Gamma[idx_gruneisen_parameter] + double(1))));
                        
                        dt_stable_global = fmin(dt_stable_local, dt_stable_global);
                    }
                }
            }
            else if (d_dim == tbox::Dimension(3))
            {
                const double g_mag = sqrt(d_gravity[0]*d_gravity[0] + d_gravity[1]*d_gravity[1] + d_gravity[2]*d_gravity[2]);
                
                /*
                 * Compute the stable time step size.
                 */
                for (int k = -d_num_subghosts_source_terms[2];
                     k < interior_dims[2] + d_num_subghosts_source_terms[2];
                     k++)
                {
                    for (int j = -d_num_subghosts_source_terms[1];
                         j < interior_dims[1] + d_num_subghosts_source_terms[1];
                         j++)
                    {
                        for (int i = -d_num_subghosts_source_terms[0];
                             i < interior_dims[0] + d_num_subghosts_source_terms[0];
                             i++)
                        {
                            const int idx_gruneisen_parameter = (i + num_subghosts_gruneisen_parameter[0]) +
                                (j + num_subghosts_gruneisen_parameter[1])*subghostcell_dims_gruneisen_parameter[0] +
                                (k + num_subghosts_gruneisen_parameter[2])*subghostcell_dims_gruneisen_parameter[0]*
                                    subghostcell_dims_gruneisen_parameter[1];
                            
                            const int idx_sound_speed = (i + num_subghosts_sound_speed[0]) +
                                (j + num_subghosts_sound_speed[1])*subghostcell_dims_sound_speed[0] +
                                (k + num_subghosts_sound_speed[2])*subghostcell_dims_sound_speed[0]*
                                    subghostcell_dims_sound_speed[1];
                            
                            const double dt_stable_local = c[idx_sound_speed]/
                                (g_mag*sqrt(double(2)*Gamma[idx_gruneisen_parameter]*(Gamma[idx_gruneisen_parameter] + double(1))));
                            
                            dt_stable_global = fmin(dt_stable_local, dt_stable_global);
                        }
                    }
                }
            }
        }
    }
    
    return dt_stable_global;
}


/*
 * Put the characteristics of this class into the restart database.
 */
void
FlowModelSourceUtilitiesSingleSpecies::putToRestart(
    const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const
{
    putToRestartBase(restart_db);
    
    if (d_has_source_terms)
    {
        HAMERS_SHARED_PTR<tbox::Database> restart_source_terms_db =
            restart_db->putDatabase("d_source_terms");
        
        putToRestartSourceBase(restart_source_terms_db);
        
        restart_source_terms_db->putBool("d_has_gravity", d_has_gravity);
        if (d_has_gravity)
        {
            restart_source_terms_db->putVector("d_gravity", d_gravity);
        }
    }
}


/*
 * Set the number of sub-ghost cells of a variable.
 * This function can be called recursively if the variables are computed recursively.
 */
void
FlowModelSourceUtilitiesSingleSpecies::setNumberOfSubGhosts(
    const hier::IntVector& num_subghosts,
    const std::string& variable_name,
    const std::string& parent_variable_name)
{
    NULL_USE(parent_variable_name);
    
    if (variable_name == "SOURCE_TERMS")
    {
        if (d_num_subghosts_source_terms > -hier::IntVector::getOne(d_dim))
        {
            if (num_subghosts > d_num_subghosts_source_terms)
            {
                d_num_subghosts_source_terms = num_subghosts;
            }
        }
        else
        {
            d_num_subghosts_source_terms = num_subghosts;
        }
    }
    
    if (variable_name == "GRUNEISEN_PARAMETER")
    {
        if (d_num_subghosts_gruneisen_parameter > -hier::IntVector::getOne(d_dim))
        {
            if (num_subghosts > d_num_subghosts_gruneisen_parameter)
            {
                d_num_subghosts_gruneisen_parameter = num_subghosts;
            }
        }
        else
        {
            d_num_subghosts_gruneisen_parameter = num_subghosts;
        }
    }
}


/*
 * Set the ghost boxes of derived cell variables.
 */
void
FlowModelSourceUtilitiesSingleSpecies::setDerivedCellVariableGhostBoxes()
{
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    const hier::Patch& patch = flow_model_tmp->getRegisteredPatch();
    const hier::Box interior_box = patch.getBox();
    
    if (d_num_subghosts_source_terms > -hier::IntVector::getOne(d_dim))
    {
        d_subghost_box_source_terms = interior_box;
        d_subghost_box_source_terms.grow(d_num_subghosts_source_terms);
        d_subghostcell_dims_source_terms = d_subghost_box_source_terms.numberCells();
    }
    
    if (d_num_subghosts_gruneisen_parameter > -hier::IntVector::getOne(d_dim))
    {
        d_subghost_box_gruneisen_parameter = interior_box;
        d_subghost_box_gruneisen_parameter.grow(d_num_subghosts_gruneisen_parameter);
        d_subghostcell_dims_gruneisen_parameter = d_subghost_box_gruneisen_parameter.numberCells();
    }
}


/*
 * Compute the cell data of Gruneisen parameter in the registered patch.
 */
void
FlowModelSourceUtilitiesSingleSpecies::computeCellDataOfGruneisenParameter()
{
    // Create empty box.
    const hier::Box empty_box(d_dim);
    
    if (d_num_subghosts_gruneisen_parameter > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_gruneisen_parameter)
        {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(d_data_gruneisen_parameter);
#endif
            
            if (d_flow_model.expired())
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "The object is not setup yet!"
                    << std::endl);
            }
            
            HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
            
            // Get the cell data of density.
            HAMERS_SHARED_PTR<pdat::CellData<double> > data_density =
                flow_model_tmp->getCellData("DENSITY");
            
            // Get the cell data of pressure.
            HAMERS_SHARED_PTR<pdat::CellData<double> > data_pressure =
                flow_model_tmp->getCellData("PRESSURE");
            
            // Get the thermodynamic properties of the species.
            
            const int num_thermo_properties = d_equation_of_state_mixing_rules->
                getNumberOfSpeciesThermodynamicProperties();
            
            std::vector<double> thermo_properties;
            std::vector<double*> thermo_properties_ptr;
            std::vector<const double*> thermo_properties_const_ptr;
            
            thermo_properties.resize(num_thermo_properties);
            thermo_properties_ptr.reserve(num_thermo_properties);
            thermo_properties_const_ptr.reserve(num_thermo_properties);
            
            for (int ti = 0; ti < num_thermo_properties; ti++)
            {
                thermo_properties_ptr.push_back(&thermo_properties[ti]);
                thermo_properties_const_ptr.push_back(&thermo_properties[ti]);
            }
            
            d_equation_of_state_mixing_rules->getSpeciesThermodynamicProperties(
                thermo_properties_ptr,
                0);
            
            // Compute the Gruneisen parameter.
            d_equation_of_state_mixing_rules->getEquationOfState()->computeGruneisenParameter(
                d_data_gruneisen_parameter,
                data_density,
                data_pressure,
                thermo_properties_const_ptr,
                empty_box);
            
            d_cell_data_computed_gruneisen_parameter = true;
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSourceUtilitiesSingleSpecies::computeCellDataOfGruneisenParameter()\n"
            << "Cell data of 'GRUNEISEN_PARAMETER' is not yet registered."
            << std::endl);
    }
}
