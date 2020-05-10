#include "flow/flow_models/single-species/FlowModelSourceUtilitiesSingleSpecies.hpp"

FlowModelSourceUtilitiesSingleSpecies::FlowModelSourceUtilitiesSingleSpecies(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
    const int& num_species,
    const boost::shared_ptr<tbox::Database>& flow_model_db,
    const boost::shared_ptr<EquationOfStateMixingRules> equation_of_state_mixing_rules):
        FlowModelSourceUtilities(
            object_name,
            dim,
            grid_geometry,
            num_species,
            2 + dim.getValue(),
            flow_model_db),
    d_has_gravity(false),
    d_equation_of_state_mixing_rules(equation_of_state_mixing_rules)
{
    if (d_has_source_terms)
    {
        boost::shared_ptr<tbox::Database> source_terms_db;
        
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
                if (source_terms_db->keyExists("d_has_gravity"))
                {
                    source_terms_db->getVector("d_gravity", d_gravity);
                }
                else
                {
                    TBOX_ERROR(d_object_name
                        << ": "
                        << "No key 'd_has_gravity' found in data for source terms."
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
        
        boost::shared_ptr<FlowModel> flow_model_tmp = d_flow_model.lock();
        const hier::Patch& patch = flow_model_tmp->getRegisteredPatch();
        const hier::Box interior_box = patch.getBox();
        
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
        
        d_num_subghosts_source_terms  = num_subghosts;
        d_subghost_box_source_terms = interior_box;
        d_subghost_box_source_terms.grow(d_num_subghosts_source_terms);
        d_subghostcell_dims_source_terms = d_subghost_box_source_terms.numberCells();
        
        if (d_has_gravity)
        {
            /*
             * Register the required derived variables in flow model.
             */
            
            // Nothing need to be done.
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
        
        // Nothing need to be done.
    }
}


/*
 * Clear cell data of different derived variables related to this class in the registered patch.
 */
void
FlowModelSourceUtilitiesSingleSpecies::clearCellData()
{
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
        
        boost::shared_ptr<FlowModel> flow_model_tmp = d_flow_model.lock();
        
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
        
        d_derived_cell_data_computed = true;
    }
}


/*
 * Compute all source terms on a patch.
 */
void
FlowModelSourceUtilitiesSingleSpecies::computeSourceTermsOnPatch(
    const boost::shared_ptr<pdat::CellVariable<double> >& variable_source,
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
        
        boost::shared_ptr<FlowModel> flow_model_tmp = d_flow_model.lock();
        const hier::Patch& patch = flow_model_tmp->getRegisteredPatch();
        const boost::shared_ptr<hier::VariableContext> data_context = flow_model_tmp->getDataContext();
        
        // Get the cell data of source.
        boost::shared_ptr<pdat::CellData<double> > source(
            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
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
            boost::shared_ptr<pdat::CellData<double> > data_momentum =
                flow_model_tmp->getCellData("MOMENTUM");
            
            // Get the cell data of density.
            boost::shared_ptr<pdat::CellData<double> > data_density =
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
            else if (d_dim == tbox::Dimension(2))
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
    }
}


/*
 * Put the characteristics of this class into the restart database.
 */
void
FlowModelSourceUtilitiesSingleSpecies::putToRestart(
    const boost::shared_ptr<tbox::Database>& restart_db) const
{
    restart_db->putBool("d_has_source_terms", d_has_source_terms);
    
    restart_db->putBool("d_has_gravity", d_has_gravity);
    if (d_has_gravity)
    {
        restart_db->putVector("d_gravity", d_gravity);
    }
}
