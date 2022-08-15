#include "flow/flow_models/five-eqn_Allaire/FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire.hpp"

FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
    const int& num_species,
    const HAMERS_SHARED_PTR<EquationOfShearViscosityMixingRules> equation_of_shear_viscosity_mixing_rules,
    const HAMERS_SHARED_PTR<EquationOfBulkViscosityMixingRules> equation_of_bulk_viscosity_mixing_rules):
        FlowModelDiffusiveFluxUtilities(
            object_name,
            dim,
            grid_geometry,
            num_species,
            dim.getValue() + 2*num_species),
        d_num_subghosts_shear_viscosity(-hier::IntVector::getOne(d_dim)),
        d_num_subghosts_bulk_viscosity(-hier::IntVector::getOne(d_dim)),
        d_subghost_box_shear_viscosity(hier::Box::getEmptyBox(d_dim)),
        d_subghost_box_bulk_viscosity(hier::Box::getEmptyBox(d_dim)),
        d_subghostcell_dims_shear_viscosity(hier::IntVector::getZero(d_dim)),
        d_subghostcell_dims_bulk_viscosity(hier::IntVector::getZero(d_dim)),
        d_cell_data_computed_shear_viscosity(false),
        d_cell_data_computed_bulk_viscosity(false),
        d_equation_of_shear_viscosity_mixing_rules(equation_of_shear_viscosity_mixing_rules),
        d_equation_of_bulk_viscosity_mixing_rules(equation_of_bulk_viscosity_mixing_rules)
{}


/*
 * Register different derived variables related to this class in the registered patch. The
 * derived variables to be registered are given as entries in a map of the variable name to
 * the number of sub-ghost cells required.
 */
void
FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::registerDerivedVariables(
    const std::unordered_map<std::string, hier::IntVector>& num_subghosts_of_data)
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
            << ": FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::"
            << "registerDerivedVariables()\n"
            << "No patch is registered yet."
            << std::endl);
    }
    
    // Check whether all or part of derived cell data is alredy computed.
    if (d_derived_cell_data_computed)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::registerDerivedVariables()\n"
            << "Derived cell data is already computed."
            << std::endl);
    }
    
    for (std::unordered_map<std::string, hier::IntVector>::const_iterator it = num_subghosts_of_data.begin();
         it != num_subghosts_of_data.end();
         it++)
    {
        if ((it->second < hier::IntVector::getZero(d_dim)) ||
            (it->second > flow_model_tmp->getNumberOfGhostCells()))
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::registerDerivedVariables()\n"
                << "The number of sub-ghost cells of variable '"
                << it->first
                << "' is not between zero and number of ghosts of conservative variables."
                << std::endl);
        }
    }
    
    if (num_subghosts_of_data.find("SHEAR_VISCOSITY") != num_subghosts_of_data.end())
    {
        std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data_flow_model;
        
        num_subghosts_of_data_flow_model.insert(std::pair<std::string, hier::IntVector>(
            "MASS_FRACTIONS",
            num_subghosts_of_data.find("SHEAR_VISCOSITY")->second));
        
        num_subghosts_of_data_flow_model.insert(std::pair<std::string, hier::IntVector>(
            "PRESSURE",
            num_subghosts_of_data.find("SHEAR_VISCOSITY")->second));
        
        num_subghosts_of_data_flow_model.insert(std::pair<std::string, hier::IntVector>(
            "SPECIES_TEMPERATURES",
            num_subghosts_of_data.find("SHEAR_VISCOSITY")->second));
        
        flow_model_tmp->registerDerivedVariables(num_subghosts_of_data_flow_model);
        
        setNumberOfSubGhosts(
            num_subghosts_of_data.find("SHEAR_VISCOSITY")->second,
            "SHEAR_VISCOSITY",
            "SHEAR_VISCOSITY");
    }
    
    if (num_subghosts_of_data.find("BULK_VISCOSITY") != num_subghosts_of_data.end())
    {
        std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data_flow_model;
        
        num_subghosts_of_data_flow_model.insert(std::pair<std::string, hier::IntVector>(
            "MASS_FRACTIONS",
            num_subghosts_of_data.find("BULK_VISCOSITY")->second));
        
        num_subghosts_of_data_flow_model.insert(std::pair<std::string, hier::IntVector>(
            "PRESSURE",
            num_subghosts_of_data.find("BULK_VISCOSITY")->second));
        
        num_subghosts_of_data_flow_model.insert(std::pair<std::string, hier::IntVector>(
            "SPECIES_TEMPERATURES",
            num_subghosts_of_data.find("BULK_VISCOSITY")->second));
        
        flow_model_tmp->registerDerivedVariables(num_subghosts_of_data_flow_model);
        
        setNumberOfSubGhosts(
            num_subghosts_of_data.find("BULK_VISCOSITY")->second,
            "BULK_VISCOSITY",
            "BULK_VISCOSITY");
    }
}


/*
 * Register the required variables for the computation of diffusive fluxes in the registered patch.
 */
void
FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::registerDerivedVariablesForDiffusiveFluxes(
    const hier::IntVector& num_subghosts,
    const bool need_side_diffusivities)
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
            << ": FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::"
            << "registerDerivedVariablesForDiffusiveFluxes()\n"
            << "No patch is registered yet."
            << std::endl);
    }
    
    // Check whether all or part of derived cell data is alredy computed.
    if (d_derived_cell_data_computed)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::registerDerivedVariablesForDiffusiveFluxes()\n"
            << "Derived cell data is already computed."
            << std::endl);
    }
    
    if ((num_subghosts < hier::IntVector::getZero(d_dim)) ||
        (num_subghosts > flow_model_tmp->getNumberOfGhostCells()))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::registerDerivedVariablesForDiffusiveFluxes()\n"
            << "The number of sub-ghost cells of variable is not between zero and number of ghosts of conservative variables."
            << std::endl);
    }
    
    /*
     * Register the required derived variables in flow model.
     */
    
    std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
    
    num_subghosts_of_data.insert(
        std::pair<std::string, hier::IntVector>("MASS_FRACTIONS", num_subghosts));
    
    num_subghosts_of_data.insert(
        std::pair<std::string, hier::IntVector>("VELOCITY", num_subghosts));
    
    num_subghosts_of_data.insert(
        std::pair<std::string, hier::IntVector>("PRESSURE", num_subghosts));
    
    num_subghosts_of_data.insert(
        std::pair<std::string, hier::IntVector>("SPECIES_TEMPERATURES", num_subghosts));
    
    flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
    
    setNumberOfSubGhosts(
        num_subghosts,
        "DIFFUSIVITIES",
        "DIFFUSIVITIES");
    
    d_need_side_diffusivities = need_side_diffusivities;
}


/*
 * Allocate memory for cell data of different registered derived variables related to this
 * class in the registered patch.
 */
void
FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::allocateMemoryForDerivedCellData()
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
    
    if (d_num_subghosts_shear_viscosity > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_shear_viscosity)
        {
            if (!d_data_shear_viscosity)
            {
                // Create the cell data of shear viscosity.
                d_data_shear_viscosity.reset(new pdat::CellData<double>(
                    interior_box, 1, d_num_subghosts_shear_viscosity));
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::allocateMemoryForDerivedCellData()\n"
                << "Cell data of 'SHEAR_VISCOSITY' is aleady computed."
                << std::endl);
        }
    }
    
    if (d_num_subghosts_bulk_viscosity > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_bulk_viscosity)
        {
            if (!d_data_bulk_viscosity)
            {
                // Create the cell data of bulk viscosity.
                d_data_bulk_viscosity.reset(new pdat::CellData<double>(
                    interior_box, 1, d_num_subghosts_bulk_viscosity));
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::allocateMemoryForDerivedCellData()\n"
                << "Cell data of 'BULK_VISCOSITY' is aleady computed."
                << std::endl);
        }
    }
    
    if (!d_need_side_diffusivities)
    {
        if (d_num_subghosts_diffusivities > -hier::IntVector::getOne(d_dim))
        {
            if (!d_cell_data_computed_diffusivities)
            {
                if (!d_data_diffusivities)
                {
                    if (d_dim == tbox::Dimension(1))
                    {
                        d_data_diffusivities.reset(new pdat::CellData<double>(
                            interior_box,
                            2,
                            d_num_subghosts_diffusivities));
                    }
                    else if (d_dim == tbox::Dimension(2))
                    {
                        d_data_diffusivities.reset(new pdat::CellData<double>(
                            interior_box,
                            9,
                            d_num_subghosts_diffusivities));
                    }
                    else if (d_dim == tbox::Dimension(3))
                    {
                        d_data_diffusivities.reset(new pdat::CellData<double>(
                            interior_box,
                            12,
                            d_num_subghosts_diffusivities));
                    }
                }
            }
            else
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::allocateMemoryForDerivedCellData()\n"
                    << "Cell data of 'DIFFUSIVITIES' is aleady computed."
                    << std::endl);
            }
        }
    }
}


/*
 * Allocate memory for side data of the diffusivities.
 */
void
FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::allocateMemoryForSideDataOfDiffusiveFluxDiffusivities()
{
    if (!d_need_side_diffusivities)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::allocateMemoryForSideDataOfDiffusiveFluxDiffusivities()\n"
            << "Side data of 'DIFFUSIVITIES' is not needed."
            << std::endl);
    }
    
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
}


/*
 * Clear cell data of different derived variables related to this class in the registered patch.
 */
void
FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::clearCellData()
{
    d_num_subghosts_diffusivities   = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_shear_viscosity = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_bulk_viscosity  = -hier::IntVector::getOne(d_dim);
    
    d_subghost_box_diffusivities   = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_shear_viscosity = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_bulk_viscosity  = hier::Box::getEmptyBox(d_dim);
    
    d_subghostcell_dims_diffusivities   = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_shear_viscosity = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_bulk_viscosity  = hier::IntVector::getZero(d_dim);
    
    d_data_diffusivities.reset();
    d_data_shear_viscosity.reset();
    d_data_bulk_viscosity.reset();
    
    d_side_data_diffusivities.reset();
    
    d_cell_data_computed_diffusivities   = false;
    d_cell_data_computed_shear_viscosity = false;
    d_cell_data_computed_bulk_viscosity  = false;
    
    d_derived_cell_data_computed = false;
    
    d_side_data_computed_diffusivities = false;
}


/*
 * Compute cell data of different registered derived variables related to this class.
 */
void
FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::computeDerivedCellData()
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
            << ": FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::"
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
    
    // Compute the shear viscosity cell data.
    if (d_num_subghosts_shear_viscosity > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_shear_viscosity)
        {
            computeCellDataOfShearViscosity();
        }
    }
    
    // Compute the bulk viscosity cell data.
    if (d_num_subghosts_bulk_viscosity > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_bulk_viscosity)
        {
            computeCellDataOfBulkViscosity();
        }
    }
    
    // Compute the diffusivities cell data.
    if (d_num_subghosts_diffusivities > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_diffusivities)
        {
            computeCellDataOfDiffusivities();
        }
    }
    
    d_derived_cell_data_computed = true;
}


/*
 * Get the cell data of one cell variable related to this class in the registered patch.
 */
HAMERS_SHARED_PTR<pdat::CellData<double> >
FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::getCellData(const std::string& variable_key)
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
            << ": FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::"
            << "getCellData()\n"
            << "No patch is registered yet."
            << std::endl);
    }
    
    HAMERS_SHARED_PTR<pdat::CellData<double> > cell_data;
    
    if (variable_key == "SHEAR_VISCOSITY")
    {
        if (!d_cell_data_computed_shear_viscosity)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::getCellData()\n"
                << "Cell data of 'SHEAR_VISCOSITY' is not registered/computed yet."
                << std::endl);
        }
        cell_data = d_data_shear_viscosity;
    }
    else if (variable_key == "BULK_VISCOSITY")
    {
        if (!d_cell_data_computed_bulk_viscosity)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::getCellData()\n"
                << "Cell data of 'BULK_VISCOSITY' is not registered/computed yet."
                << std::endl);
        }
        cell_data = d_data_bulk_viscosity;
    }
    
    return cell_data;
}


/*
 * Get the cell data of different cell variables related to this class in the registered patch.
 */
std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >
FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::getCellData(
    const std::vector<std::string>& variable_keys)
{
    std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > cell_data(
        static_cast<int>(variable_keys.size()));
    
    for (int vi = 0; static_cast<int>(variable_keys.size()); vi++)
    {
        cell_data[vi] = getCellData(variable_keys[vi]);
    }
    
    return cell_data;
}


/*
 * Get the variables for the derivatives in the diffusive fluxes.
 */
void
FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::getCellDataOfDiffusiveFluxVariablesForDerivative(
    std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& derivative_var_data,
    std::vector<std::vector<int> >& derivative_var_component_idx,
    const DIRECTION::TYPE& flux_direction,
    const DIRECTION::TYPE& derivative_direction)
{
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
    derivative_var_data.resize(d_num_eqn);
    derivative_var_component_idx.resize(d_num_eqn);
    
    // Get the cell data of velocity.
    HAMERS_SHARED_PTR<pdat::CellData<double> > data_velocity =
        flow_model_tmp->getCellData("VELOCITY");
    
    if (d_dim == tbox::Dimension(1))
    {
        switch (flux_direction)
        {
            case DIRECTION::X_DIRECTION:
            {
                switch (derivative_direction)
                {
                    case DIRECTION::X_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            derivative_var_data[si].resize(0);
                            derivative_var_component_idx[si].resize(0);
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[d_num_species].resize(1);
                        derivative_var_component_idx[d_num_species].resize(1);
                        
                        // Variable u.
                        derivative_var_data[d_num_species][0] = data_velocity;
                        derivative_var_component_idx[d_num_species][0] = 0;
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[d_num_species + 1].resize(1);
                        derivative_var_component_idx[d_num_species + 1].resize(1);
                        
                        // Variable u.
                        derivative_var_data[d_num_species + 1][0] = data_velocity;
                        derivative_var_component_idx[d_num_species + 1][0] = 0;
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::"
                            << "getCellDataOfDiffusiveFluxVariablesForDerivative()\n"
                            << "There are only x-direction for one-dimensional problem."
                            << std::endl);
                    }
                }
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::"
                    << "getCellDataOfDiffusiveFluxVariablesForDerivative()\n"
                    << "There are only x-direction for one-dimensional problem."
                    << std::endl);
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        switch (flux_direction)
        {
            case DIRECTION::X_DIRECTION:
            {
                switch (derivative_direction)
                {
                    case DIRECTION::X_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            derivative_var_data[si].resize(0);
                            derivative_var_component_idx[si].resize(0);
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[d_num_species].resize(1);
                        derivative_var_component_idx[d_num_species].resize(1);
                        
                        // Variable u.
                        derivative_var_data[d_num_species][0] = data_velocity;
                        derivative_var_component_idx[d_num_species][0] = 0;
                        
                        derivative_var_data[d_num_species + 1].resize(1);
                        derivative_var_component_idx[d_num_species + 1].resize(1);
                        
                        // Variable v.
                        derivative_var_data[d_num_species + 1][0] = data_velocity;
                        derivative_var_component_idx[d_num_species + 1][0] = 1;
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[d_num_species + 2].resize(2);
                        derivative_var_component_idx[d_num_species + 2].resize(2);
                        
                        // Variable u.
                        derivative_var_data[d_num_species + 2][0] = data_velocity;
                        derivative_var_component_idx[d_num_species + 2][0] = 0;
                        
                        // Variable v.
                        derivative_var_data[d_num_species + 2][1] = data_velocity;
                        derivative_var_component_idx[d_num_species + 2][1] = 1;
                        
                        break;
                    }
                    case DIRECTION::Y_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            derivative_var_data[si].resize(0);
                            derivative_var_component_idx[si].resize(0);
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[d_num_species].resize(1);
                        derivative_var_component_idx[d_num_species].resize(1);
                        
                        // Variable v.
                        derivative_var_data[d_num_species][0] = data_velocity;
                        derivative_var_component_idx[d_num_species][0] = 1;
                        
                        derivative_var_data[d_num_species + 1].resize(1);
                        derivative_var_component_idx[d_num_species + 1].resize(1);
                        
                        // Variable u.
                        derivative_var_data[d_num_species + 1][0] = data_velocity;
                        derivative_var_component_idx[d_num_species + 1][0] = 0;
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[d_num_species + 2].resize(2);
                        derivative_var_component_idx[d_num_species + 2].resize(2);
                        
                        // Variable u.
                        derivative_var_data[d_num_species + 2][0] = data_velocity;
                        derivative_var_component_idx[d_num_species + 2][0] = 0;
                        
                        // Variable v.
                        derivative_var_data[d_num_species + 2][1] = data_velocity;
                        derivative_var_component_idx[d_num_species + 2][1] = 1;
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::"
                            << "getCellDataOfDiffusiveFluxVariablesForDerivative()\n"
                            << "There are only x-direction and y-direction for two-dimensional problem."
                            << std::endl);
                    }
                }
                
                break;
            }
            case DIRECTION::Y_DIRECTION:
            {
                switch (derivative_direction)
                {
                    case DIRECTION::X_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            derivative_var_data[si].resize(0);
                            derivative_var_component_idx[si].resize(0);
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[d_num_species].resize(1);
                        derivative_var_component_idx[d_num_species].resize(1);
                        
                        // Variable v.
                        derivative_var_data[d_num_species][0] = data_velocity;
                        derivative_var_component_idx[d_num_species][0] = 1;
                        
                        derivative_var_data[d_num_species + 1].resize(1);
                        derivative_var_component_idx[d_num_species + 1].resize(1);
                        
                        // Variable u.
                        derivative_var_data[d_num_species + 1][0] = data_velocity;
                        derivative_var_component_idx[d_num_species + 1][0] = 0;
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[d_num_species + 2].resize(2);
                        derivative_var_component_idx[d_num_species + 2].resize(2);
                        
                        // Variable u.
                        derivative_var_data[d_num_species + 2][0] = data_velocity;
                        derivative_var_component_idx[d_num_species + 2][0] = 0;
                        
                        // Variable v.
                        derivative_var_data[d_num_species + 2][1] = data_velocity;
                        derivative_var_component_idx[d_num_species + 2][1] = 1;
                        
                        break;
                    }
                    case DIRECTION::Y_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            derivative_var_data[si].resize(0);
                            derivative_var_component_idx[si].resize(0);
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[d_num_species].resize(1);
                        derivative_var_component_idx[d_num_species].resize(1);
                        
                        // Variable u.
                        derivative_var_data[d_num_species][0] = data_velocity;
                        derivative_var_component_idx[d_num_species][0] = 0;
                        
                        derivative_var_data[d_num_species + 1].resize(1);
                        derivative_var_component_idx[d_num_species + 1].resize(1);
                        
                        // Variable v.
                        derivative_var_data[d_num_species + 1][0] = data_velocity;
                        derivative_var_component_idx[d_num_species + 1][0] = 1;
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[d_num_species + 2].resize(2);
                        derivative_var_component_idx[d_num_species + 2].resize(2);
                        
                        // Variable u.
                        derivative_var_data[d_num_species + 2][0] = data_velocity;
                        derivative_var_component_idx[d_num_species + 2][0] = 0;
                        
                        // Variable v.
                        derivative_var_data[d_num_species + 2][1] = data_velocity;
                        derivative_var_component_idx[d_num_species + 2][1] = 1;
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::"
                            << "getCellDataOfDiffusiveFluxVariablesForDerivative()\n"
                            << "There are only x-direction and y-direction for two-dimensional problem."
                            << std::endl);
                    }
                }
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::"
                    << "getCellDataOfDiffusiveFluxVariablesForDerivative()\n"
                    << "There are only x-direction and y-direction for two-dimensional problem."
                    << std::endl);
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        switch (flux_direction)
        {
            case DIRECTION::X_DIRECTION:
            {
                switch (derivative_direction)
                {
                    case DIRECTION::X_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            derivative_var_data[si].resize(0);
                            derivative_var_component_idx[si].resize(0);
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[d_num_species].resize(1);
                        derivative_var_component_idx[d_num_species].resize(1);
                        
                        // Variable u.
                        derivative_var_data[d_num_species][0] = data_velocity;
                        derivative_var_component_idx[d_num_species][0] = 0;
                        
                        derivative_var_data[d_num_species + 1].resize(1);
                        derivative_var_component_idx[d_num_species + 1].resize(1);
                        
                        // Variable v.
                        derivative_var_data[d_num_species + 1][0] = data_velocity;
                        derivative_var_component_idx[d_num_species + 1][0] = 1;
                        
                        derivative_var_data[d_num_species + 2].resize(1);
                        derivative_var_component_idx[d_num_species + 2].resize(1);
                        
                        // Variable w.
                        derivative_var_data[d_num_species + 2][0] = data_velocity;
                        derivative_var_component_idx[d_num_species + 2][0] = 2;
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[d_num_species + 3].resize(3);
                        derivative_var_component_idx[d_num_species + 3].resize(3);
                        
                        // Variable u.
                        derivative_var_data[d_num_species + 3][0] = data_velocity;
                        derivative_var_component_idx[d_num_species + 3][0] = 0;
                        
                        // Variable v.
                        derivative_var_data[d_num_species + 3][1] = data_velocity;
                        derivative_var_component_idx[d_num_species + 3][1] = 1;
                        
                        // Variable w.
                        derivative_var_data[d_num_species + 3][2] = data_velocity;
                        derivative_var_component_idx[d_num_species + 3][2] = 2;
                        
                        break;
                    }
                    case DIRECTION::Y_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            derivative_var_data[si].resize(0);
                            derivative_var_component_idx[si].resize(0);
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[d_num_species].resize(1);
                        derivative_var_component_idx[d_num_species].resize(1);
                        
                        // Variable v.
                        derivative_var_data[d_num_species][0] = data_velocity;
                        derivative_var_component_idx[d_num_species][0] = 1;
                        
                        derivative_var_data[d_num_species + 1].resize(1);
                        derivative_var_component_idx[d_num_species + 1].resize(1);
                        
                        // Variable u.
                        derivative_var_data[d_num_species + 1][0] = data_velocity;
                        derivative_var_component_idx[d_num_species + 1][0] = 0;
                        
                        derivative_var_data[d_num_species + 2].resize(0);
                        derivative_var_component_idx[d_num_species + 2].resize(0);
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[d_num_species + 3].resize(2);
                        derivative_var_component_idx[d_num_species + 3].resize(2);
                        
                        // Variable u.
                        derivative_var_data[d_num_species + 3][0] = data_velocity;
                        derivative_var_component_idx[d_num_species + 3][0] = 0;
                        
                        // Variable v.
                        derivative_var_data[d_num_species + 3][1] = data_velocity;
                        derivative_var_component_idx[d_num_species + 3][1] = 1;
                        
                        break;
                    }
                    case DIRECTION::Z_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            derivative_var_data[si].resize(0);
                            derivative_var_component_idx[si].resize(0);
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[d_num_species].resize(1);
                        derivative_var_component_idx[d_num_species].resize(1);
                        
                        // Variable w.
                        derivative_var_data[d_num_species][0] = data_velocity;
                        derivative_var_component_idx[d_num_species][0] = 2;
                        
                        derivative_var_data[d_num_species + 1].resize(0);
                        derivative_var_component_idx[d_num_species + 1].resize(0);
                        
                        derivative_var_data[d_num_species + 2].resize(1);
                        derivative_var_component_idx[d_num_species + 2].resize(1);
                        
                        // Variable u.
                        derivative_var_data[d_num_species + 2][0] = data_velocity;
                        derivative_var_component_idx[d_num_species + 2][0] = 0;
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[d_num_species + 3].resize(2);
                        derivative_var_component_idx[d_num_species + 3].resize(2);
                        
                        // Variable u.
                        derivative_var_data[d_num_species + 3][0] = data_velocity;
                        derivative_var_component_idx[d_num_species + 3][0] = 0;
                        
                        // Variable w.
                        derivative_var_data[d_num_species + 3][1] = data_velocity;
                        derivative_var_component_idx[d_num_species + 3][1] = 2;
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::"
                            << "getCellDataOfDiffusiveFluxVariablesForDerivative()\n"
                            << "There are only x-direction, y-direction and z-direction for three-dimensional problem."
                            << std::endl);
                    }
                }
                
                break;
            }
            case DIRECTION::Y_DIRECTION:
            {
                switch (derivative_direction)
                {
                    case DIRECTION::X_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            derivative_var_data[si].resize(0);
                            derivative_var_component_idx[si].resize(0);
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[d_num_species].resize(1);
                        derivative_var_component_idx[d_num_species].resize(1);
                        
                        // Variable v.
                        derivative_var_data[d_num_species][0] = data_velocity;
                        derivative_var_component_idx[d_num_species][0] = 1;
                        
                        derivative_var_data[d_num_species + 1].resize(1);
                        derivative_var_component_idx[d_num_species + 1].resize(1);
                        
                        // Variable u.
                        derivative_var_data[d_num_species + 1][0] = data_velocity;
                        derivative_var_component_idx[d_num_species + 1][0] = 0;
                        
                        derivative_var_data[d_num_species + 2].resize(0);
                        derivative_var_component_idx[d_num_species + 2].resize(0);
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[d_num_species + 3].resize(2);
                        derivative_var_component_idx[d_num_species + 3].resize(2);
                        
                        // Variable u.
                        derivative_var_data[d_num_species + 3][0] = data_velocity;
                        derivative_var_component_idx[d_num_species + 3][0] = 0;
                        
                        // Variable v.
                        derivative_var_data[d_num_species + 3][1] = data_velocity;
                        derivative_var_component_idx[d_num_species + 3][1] = 1;
                        
                        break;
                    }
                    case DIRECTION::Y_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            derivative_var_data[si].resize(0);
                            derivative_var_component_idx[si].resize(0);
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[d_num_species].resize(1);
                        derivative_var_component_idx[d_num_species].resize(1);
                        
                        // Variable u.
                        derivative_var_data[d_num_species][0] = data_velocity;
                        derivative_var_component_idx[d_num_species][0] = 0;
                        
                        derivative_var_data[d_num_species + 1].resize(1);
                        derivative_var_component_idx[d_num_species + 1].resize(1);
                        
                        // Variable v.
                        derivative_var_data[d_num_species + 1][0] = data_velocity;
                        derivative_var_component_idx[d_num_species + 1][0] = 1;
                        
                        derivative_var_data[d_num_species + 2].resize(1);
                        derivative_var_component_idx[d_num_species + 2].resize(1);
                        
                        // Variable w.
                        derivative_var_data[d_num_species + 2][0] = data_velocity;
                        derivative_var_component_idx[d_num_species + 2][0] = 2;
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[d_num_species + 3].resize(3);
                        derivative_var_component_idx[d_num_species + 3].resize(3);
                        
                        // Variable u.
                        derivative_var_data[d_num_species + 3][0] = data_velocity;
                        derivative_var_component_idx[d_num_species + 3][0] = 0;
                        
                        // Variable v.
                        derivative_var_data[d_num_species + 3][1] = data_velocity;
                        derivative_var_component_idx[d_num_species + 3][1] = 1;
                        
                        // Variable w.
                        derivative_var_data[d_num_species + 3][2] = data_velocity;
                        derivative_var_component_idx[d_num_species + 3][2] = 2;
                        
                        break;
                    }
                    case DIRECTION::Z_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            derivative_var_data[si].resize(0);
                            derivative_var_component_idx[si].resize(0);
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[d_num_species].resize(0);
                        derivative_var_component_idx[d_num_species].resize(0);
                        
                        derivative_var_data[d_num_species + 1].resize(1);
                        derivative_var_component_idx[d_num_species + 1].resize(1);
                        
                        // Variable w.
                        derivative_var_data[d_num_species + 1][0] = data_velocity;
                        derivative_var_component_idx[d_num_species + 1][0] = 2;
                        
                        derivative_var_data[d_num_species + 2].resize(1);
                        derivative_var_component_idx[d_num_species + 2].resize(1);
                        
                        // Variable v.
                        derivative_var_data[d_num_species + 2][0] = data_velocity;
                        derivative_var_component_idx[d_num_species + 2][0] = 1;
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[d_num_species + 3].resize(2);
                        derivative_var_component_idx[d_num_species + 3].resize(2);
                        
                        // Variable v.
                        derivative_var_data[d_num_species + 3][0] = data_velocity;
                        derivative_var_component_idx[d_num_species + 3][0] = 1;
                        
                        // Variable w.
                        derivative_var_data[d_num_species + 3][1] = data_velocity;
                        derivative_var_component_idx[d_num_species + 3][1] = 2;
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::"
                            << "getCellDataOfDiffusiveFluxVariablesForDerivative()\n"
                            << "There are only x-direction, y-direction and z-direction for three-dimensional problem."
                            << std::endl);
                    }
                }
                
                break;
            }
            case DIRECTION::Z_DIRECTION:
            {
                switch (derivative_direction)
                {
                    case DIRECTION::X_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            derivative_var_data[si].resize(0);
                            derivative_var_component_idx[si].resize(0);
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[d_num_species].resize(1);
                        derivative_var_component_idx[d_num_species].resize(1);
                        
                        // Variable w.
                        derivative_var_data[d_num_species][0] = data_velocity;
                        derivative_var_component_idx[d_num_species][0] = 2;
                        
                        derivative_var_data[d_num_species + 1].resize(0);
                        derivative_var_component_idx[d_num_species + 1].resize(0);
                        
                        derivative_var_data[d_num_species + 2].resize(1);
                        derivative_var_component_idx[d_num_species + 2].resize(1);
                        
                        // Variable u.
                        derivative_var_data[d_num_species + 2][0] = data_velocity;
                        derivative_var_component_idx[d_num_species + 2][0] = 0;
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[d_num_species + 3].resize(2);
                        derivative_var_component_idx[d_num_species + 3].resize(2);
                        
                        // Variable u.
                        derivative_var_data[d_num_species + 3][0] = data_velocity;
                        derivative_var_component_idx[d_num_species + 3][0] = 0;
                        
                        // Variable w.
                        derivative_var_data[d_num_species + 3][1] = data_velocity;
                        derivative_var_component_idx[d_num_species + 3][1] = 2;
                        
                        break;
                    }
                    case DIRECTION::Y_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            derivative_var_data[si].resize(0);
                            derivative_var_component_idx[si].resize(0);
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[d_num_species].resize(0);
                        derivative_var_component_idx[d_num_species].resize(0);
                        
                        derivative_var_data[d_num_species + 1].resize(1);
                        derivative_var_component_idx[d_num_species + 1].resize(1);
                        
                        // Variable w.
                        derivative_var_data[d_num_species + 1][0] = data_velocity;
                        derivative_var_component_idx[d_num_species + 1][0] = 2;
                        
                        derivative_var_data[d_num_species + 2].resize(1);
                        derivative_var_component_idx[d_num_species + 2].resize(1);
                        
                        // Variable v.
                        derivative_var_data[d_num_species + 2][0] = data_velocity;
                        derivative_var_component_idx[d_num_species + 2][0] = 1;
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[d_num_species + 3].resize(2);
                        derivative_var_component_idx[d_num_species + 3].resize(2);
                        
                        // Variable v.
                        derivative_var_data[d_num_species + 3][0] = data_velocity;
                        derivative_var_component_idx[d_num_species + 3][0] = 1;
                        
                        // Variable w.
                        derivative_var_data[d_num_species + 3][1] = data_velocity;
                        derivative_var_component_idx[d_num_species + 3][1] = 2;
                        
                        break;
                    }
                    case DIRECTION::Z_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            derivative_var_data[si].resize(0);
                            derivative_var_component_idx[si].resize(0);
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[d_num_species].resize(1);
                        derivative_var_component_idx[d_num_species].resize(1);
                        
                        // Variable u.
                        derivative_var_data[d_num_species][0] = data_velocity;
                        derivative_var_component_idx[d_num_species][0] = 0;
                        
                        derivative_var_data[d_num_species + 1].resize(1);
                        derivative_var_component_idx[d_num_species + 1].resize(1);
                        
                        // Variable v.
                        derivative_var_data[d_num_species + 1][0] = data_velocity;
                        derivative_var_component_idx[d_num_species + 1][0] = 1;
                        
                        derivative_var_data[d_num_species + 2].resize(1);
                        derivative_var_component_idx[d_num_species + 2].resize(1);
                        
                        // Variable w.
                        derivative_var_data[d_num_species + 2][0] = data_velocity;
                        derivative_var_component_idx[d_num_species + 2][0] = 2;
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[d_num_species + 3].resize(3);
                        derivative_var_component_idx[d_num_species + 3].resize(3);
                        
                        // Variable u.
                        derivative_var_data[d_num_species + 3][0] = data_velocity;
                        derivative_var_component_idx[d_num_species + 3][0] = 0;
                        
                        // Variable v.
                        derivative_var_data[d_num_species + 3][1] = data_velocity;
                        derivative_var_component_idx[d_num_species + 3][1] = 1;
                        
                        // Variable w.
                        derivative_var_data[d_num_species + 3][2] = data_velocity;
                        derivative_var_component_idx[d_num_species + 3][2] = 2;
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::"
                            << "getCellDataOfDiffusiveFluxVariablesForDerivative()\n"
                            << "There are only x-direction, y-direction and z-direction for three-dimensional problem."
                            << std::endl);
                    }
                }
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::"
                    << "getCellDataOfDiffusiveFluxVariablesForDerivative()\n"
                    << "There are only x-direction, y-direction and z-direction for three-dimensional problem."
                    << std::endl);
            }
        }
    }
}


/*
 * Get the diffusivities in the diffusive flux.
 */
void
FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::getCellDataOfDiffusiveFluxDiffusivities(
    std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& diffusivities_data,
    std::vector<std::vector<int> >& diffusivities_component_idx,
    const DIRECTION::TYPE& flux_direction,
    const DIRECTION::TYPE& derivative_direction)
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
            << ": FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::"
            << "getCellDataOfDiffusiveFluxDiffusivities()\n"
            << "No patch is registered yet."
            << std::endl);
    }
    
    if (!d_cell_data_computed_diffusivities)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::getCellDataOfDiffusiveFluxDiffusivities()\n"
            << "Cell data of 'DIFFUSIVITIES' is not registered/computed yet."
            << std::endl);
    }
    
    diffusivities_data.resize(d_num_eqn);
    diffusivities_component_idx.resize(d_num_eqn);
    
    if (d_dim == tbox::Dimension(1))
    {
        switch (flux_direction)
        {
            case DIRECTION::X_DIRECTION:
            {
                switch (derivative_direction)
                {
                    case DIRECTION::X_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            diffusivities_data[si].resize(0);
                            diffusivities_component_idx[si].resize(0);
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[d_num_species].resize(1);
                        diffusivities_component_idx[d_num_species].resize(1);
                        
                        // -(4/3*mu + mu_v).
                        diffusivities_data[d_num_species][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species][0] = 0;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[d_num_species + 1].resize(1);
                        diffusivities_component_idx[d_num_species + 1].resize(1);
                        
                        // -u*(4/3*mu + mu_v).
                        diffusivities_data[d_num_species + 1][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 1][0] = 1;
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::"
                            << "getCellDataOfDiffusiveFluxDiffusivities()\n"
                            << "There are only x-direction for one-dimensional problem."
                            << std::endl);
                    }
                }
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::getCellDataOfDiffusiveFluxDiffusivities()\n"
                    << "There are only x-direction for one-dimensional problem."
                    << std::endl);
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        switch (flux_direction)
        {
            case DIRECTION::X_DIRECTION:
            {
                switch (derivative_direction)
                {
                    case DIRECTION::X_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            diffusivities_data[si].resize(0);
                            diffusivities_component_idx[si].resize(0);
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[d_num_species].resize(1);
                        diffusivities_component_idx[d_num_species].resize(1);
                        
                        // -(4/3*mu + mu_v).
                        diffusivities_data[d_num_species][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species][0] = 0;
                        
                        diffusivities_data[d_num_species + 1].resize(1);
                        diffusivities_component_idx[d_num_species + 1].resize(1);
                        
                        // -mu.
                        diffusivities_data[d_num_species + 1][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 1][0] = 2;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[d_num_species + 2].resize(2);
                        diffusivities_component_idx[d_num_species + 2].resize(2);
                        
                        // -u*(4/3*mu + mu_v).
                        diffusivities_data[d_num_species + 2][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 2][0] = 3;
                        
                        // -v*mu.
                        diffusivities_data[d_num_species + 2][1] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 2][1] = 8;
                        
                        break;
                    }
                    case DIRECTION::Y_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            diffusivities_data[si].resize(0);
                            diffusivities_component_idx[si].resize(0);
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[d_num_species].resize(1);
                        diffusivities_component_idx[d_num_species].resize(1);
                        
                        // 2/3*mu - mu_v.
                        diffusivities_data[d_num_species][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species][0] = 1;
                        
                        diffusivities_data[d_num_species + 1].resize(1);
                        diffusivities_component_idx[d_num_species + 1].resize(1);
                        
                        // -mu.
                        diffusivities_data[d_num_species + 1][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 1][0] = 2;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[d_num_species + 2].resize(2);
                        diffusivities_component_idx[d_num_species + 2].resize(2);
                        
                        // -v*mu.
                        diffusivities_data[d_num_species + 2][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 2][0] = 8;
                        
                        // u*(2/3*mu - mu_v).
                        diffusivities_data[d_num_species + 2][1] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 2][1] = 5;
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::"
                            << "getCellDataOfDiffusiveFluxDiffusivities()\n"
                            << "There are only x-direction and y-direction for two-dimensional problem."
                            << std::endl);
                    }
                }
                
                break;
            }
            case DIRECTION::Y_DIRECTION:
            {
                switch (derivative_direction)
                {
                    case DIRECTION::X_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            diffusivities_data[si].resize(0);
                            diffusivities_component_idx[si].resize(0);
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[d_num_species].resize(1);
                        diffusivities_component_idx[d_num_species].resize(1);
                        
                        // -mu.
                        diffusivities_data[d_num_species][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species][0] = 2;
                        
                        diffusivities_data[d_num_species + 1].resize(1);
                        diffusivities_component_idx[d_num_species + 1].resize(1);
                        
                        // 2/3*mu - mu_v.
                        diffusivities_data[d_num_species + 1][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 1][0] = 1;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[d_num_species + 2].resize(2);
                        diffusivities_component_idx[d_num_species + 2].resize(2);
                        
                        // v*(2/3*mu - mu_v).
                        diffusivities_data[d_num_species + 2][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 2][0] = 6;
                        
                        // -u*mu.
                        diffusivities_data[d_num_species + 2][1] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 2][1] = 7;
                        
                        break;
                    }
                    case DIRECTION::Y_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            diffusivities_data[si].resize(0);
                            diffusivities_component_idx[si].resize(0);
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[d_num_species].resize(1);
                        diffusivities_component_idx[d_num_species].resize(1);
                        
                        // -mu.
                        diffusivities_data[d_num_species][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species][0] = 2;
                        
                        diffusivities_data[d_num_species + 1].resize(1);
                        diffusivities_component_idx[d_num_species + 1].resize(1);
                        
                        // -(4/3*mu + mu_v).
                        diffusivities_data[d_num_species + 1][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 1][0] = 0;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[d_num_species + 2].resize(2);
                        diffusivities_component_idx[d_num_species + 2].resize(2);
                        
                        // -u*mu.
                        diffusivities_data[d_num_species + 2][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 2][0] = 7;
                        
                        // -v*(4/3*mu + mu_v).
                        diffusivities_data[d_num_species + 2][1] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 2][1] = 4;
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::"
                            << "getCellDataOfDiffusiveFluxDiffusivities()\n"
                            << "There are only x-direction and y-direction for two-dimensional problem."
                            << std::endl);
                    }
                }
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::getCellDataOfDiffusiveFluxDiffusivities()\n"
                    << "There are only x-direction and y-direction for two-dimensional problem."
                    << std::endl);
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        switch (flux_direction)
        {
            case DIRECTION::X_DIRECTION:
            {
                switch (derivative_direction)
                {
                    case DIRECTION::X_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            diffusivities_data[si].resize(0);
                            diffusivities_component_idx[si].resize(0);
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[d_num_species].resize(1);
                        diffusivities_component_idx[d_num_species].resize(1);
                        
                        // -(4/3*mu + mu_v).
                        diffusivities_data[d_num_species][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species][0] = 0;
                        
                        diffusivities_data[d_num_species + 1].resize(1);
                        diffusivities_component_idx[d_num_species + 1].resize(1);
                        
                        // -mu.
                        diffusivities_data[d_num_species + 1][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 1][0] = 2;
                        
                        diffusivities_data[d_num_species + 2].resize(1);
                        diffusivities_component_idx[d_num_species + 2].resize(1);
                        
                        // -mu.
                        diffusivities_data[d_num_species + 2][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 2][0] = 2;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[d_num_species + 3].resize(3);
                        diffusivities_component_idx[d_num_species + 3].resize(3);
                        
                        // -u*(4/3*mu + mu_v).
                        diffusivities_data[d_num_species + 3][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][0] = 3;
                        
                        // -v*mu.
                        diffusivities_data[d_num_species + 3][1] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][1] = 10;
                        
                        // -w*mu.
                        diffusivities_data[d_num_species + 3][2] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][2] = 11;
                        
                        break;
                    }
                    case DIRECTION::Y_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            diffusivities_data[si].resize(0);
                            diffusivities_component_idx[si].resize(0);
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[d_num_species].resize(1);
                        diffusivities_component_idx[d_num_species].resize(1);
                        
                        // 2/3*mu - mu_v.
                        diffusivities_data[d_num_species][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species][0] = 1;
                        
                        diffusivities_data[d_num_species + 1].resize(1);
                        diffusivities_component_idx[d_num_species + 1].resize(1);
                        
                        // -mu.
                        diffusivities_data[d_num_species + 1][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 1][0] = 2;
                        
                        diffusivities_data[d_num_species + 2].resize(0);
                        diffusivities_component_idx[d_num_species + 2].resize(0);
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[d_num_species + 3].resize(2);
                        diffusivities_component_idx[d_num_species + 3].resize(2);
                        
                        // -v*mu.
                        diffusivities_data[d_num_species + 3][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][0] = 10;
                        
                        // u*(2/3*mu - mu_v).
                        diffusivities_data[d_num_species + 3][1] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][1] = 6;
                        
                        break;
                    }
                    case DIRECTION::Z_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            diffusivities_data[si].resize(0);
                            diffusivities_component_idx[si].resize(0);
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[d_num_species].resize(1);
                        diffusivities_component_idx[d_num_species].resize(1);
                        
                        // 2/3*mu - mu_v.
                        diffusivities_data[d_num_species][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species][0] = 1;
                        
                        diffusivities_data[d_num_species + 1].resize(0);
                        diffusivities_component_idx[d_num_species + 1].resize(0);
                        
                        diffusivities_data[d_num_species + 2].resize(1);
                        diffusivities_component_idx[d_num_species + 2].resize(1);
                        
                        // -mu.
                        diffusivities_data[d_num_species + 2][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 2][0] = 2;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[d_num_species + 3].resize(2);
                        diffusivities_component_idx[d_num_species + 3].resize(2);
                        
                        // -w*mu.
                        diffusivities_data[d_num_species + 3][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][0] = 11;
                        
                        // u*(2/3*mu - mu_v).
                        diffusivities_data[d_num_species + 3][1] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][1] = 6;
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::"
                            << "getCellDataOfDiffusiveFluxDiffusivities()\n"
                            << "There are only x-direction, y-direction and z-direction for three-dimensional problem."
                            << std::endl);
                    }
                }
                
                break;
            }
            case DIRECTION::Y_DIRECTION:
            {
                switch (derivative_direction)
                {
                    case DIRECTION::X_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            diffusivities_data[si].resize(0);
                            diffusivities_component_idx[si].resize(0);
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[d_num_species].resize(1);
                        diffusivities_component_idx[d_num_species].resize(1);
                        
                        // -mu.
                        diffusivities_data[d_num_species][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species][0] = 2;
                        
                        diffusivities_data[d_num_species + 1].resize(1);
                        diffusivities_component_idx[d_num_species + 1].resize(1);
                        
                        // 2/3*mu - mu_v.
                        diffusivities_data[d_num_species + 1][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 1][0] = 1;
                        
                        diffusivities_data[d_num_species + 2].resize(0);
                        diffusivities_component_idx[d_num_species + 2].resize(0);
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[d_num_species + 3].resize(2);
                        diffusivities_component_idx[d_num_species + 3].resize(2);
                        
                        // v*(2/3*mu - mu_v).
                        diffusivities_data[d_num_species + 3][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][0] = 7;
                        
                        // -u*mu.
                        diffusivities_data[d_num_species + 3][1] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][1] = 9;
                        
                        break;
                    }
                    case DIRECTION::Y_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            diffusivities_data[si].resize(0);
                            diffusivities_component_idx[si].resize(0);
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[d_num_species].resize(1);
                        diffusivities_component_idx[d_num_species].resize(1);
                        
                        // -mu.
                        diffusivities_data[d_num_species][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species][0] = 2;
                        
                        diffusivities_data[d_num_species + 1].resize(1);
                        diffusivities_component_idx[d_num_species + 1].resize(1);
                        
                        // -(4/3*mu + mu_v).
                        diffusivities_data[d_num_species + 1][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 1][0] = 0;
                        
                        diffusivities_data[d_num_species + 2].resize(1);
                        diffusivities_component_idx[d_num_species + 2].resize(1);
                        
                        // -mu.
                        diffusivities_data[d_num_species + 2][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 2][0] = 2;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[d_num_species + 3].resize(3);
                        diffusivities_component_idx[d_num_species + 3].resize(3);
                        
                        // -u*mu.
                        diffusivities_data[d_num_species + 3][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][0] = 9;
                        
                        // -v*(4/3*mu + mu_v).
                        diffusivities_data[d_num_species + 3][1] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][1] = 4;
                        
                        // -w*mu.
                        diffusivities_data[d_num_species + 3][2] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][2] = 11;
                        
                        break;
                    }
                    case DIRECTION::Z_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            diffusivities_data[si].resize(0);
                            diffusivities_component_idx[si].resize(0);
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[d_num_species].resize(0);
                        diffusivities_component_idx[d_num_species].resize(0);
                        
                        diffusivities_data[d_num_species + 1].resize(1);
                        diffusivities_component_idx[d_num_species + 1].resize(1);
                        
                        // 2/3*(mu - mu_v).
                        diffusivities_data[d_num_species + 1][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 1][0] = 1;
                        
                        diffusivities_data[d_num_species + 2].resize(1);
                        diffusivities_component_idx[d_num_species + 2].resize(1);
                        
                        // -mu.
                        diffusivities_data[d_num_species + 2][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 2][0] = 2;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[d_num_species + 3].resize(2);
                        diffusivities_component_idx[d_num_species + 3].resize(2);
                        
                        // -w*u.
                        diffusivities_data[d_num_species + 3][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][0] = 11;
                        
                        // v*(2/3*mu - mu_v).
                        diffusivities_data[d_num_species + 3][1] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][1] = 7;
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::"
                            << "getCellDataOfDiffusiveFluxDiffusivities()\n"
                            << "There are only x-direction, y-direction and z-direction for three-dimensional problem."
                            << std::endl);
                    }
                }
                
                break;
            }
            case DIRECTION::Z_DIRECTION:
            {
                switch (derivative_direction)
                {
                    case DIRECTION::X_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            diffusivities_data[si].resize(0);
                            diffusivities_component_idx[si].resize(0);
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[d_num_species].resize(1);
                        diffusivities_component_idx[d_num_species].resize(1);
                        
                        // -mu.
                        diffusivities_data[d_num_species][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species][0] = 2;
                        
                        diffusivities_data[d_num_species + 1].resize(0);
                        diffusivities_component_idx[d_num_species + 1].resize(0);
                        
                        diffusivities_data[d_num_species + 2].resize(1);
                        diffusivities_component_idx[d_num_species + 2].resize(1);
                        
                        // 2/3*mu - mu_v.
                        diffusivities_data[d_num_species + 2][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 2][0] = 1;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[d_num_species + 3].resize(2);
                        diffusivities_component_idx[d_num_species + 3].resize(2);
                        
                        // w*(2/3*mu - mu_v).
                        diffusivities_data[d_num_species + 3][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][0] = 8;
                        
                        // -u*mu.
                        diffusivities_data[d_num_species + 3][1] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][1] = 9;
                        
                        break;
                    }
                    case DIRECTION::Y_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            diffusivities_data[si].resize(0);
                            diffusivities_component_idx[si].resize(0);
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[d_num_species].resize(0);
                        diffusivities_component_idx[d_num_species].resize(0);
                        
                        diffusivities_data[d_num_species + 1].resize(1);
                        diffusivities_component_idx[d_num_species + 1].resize(1);
                        
                        // -mu.
                        diffusivities_data[d_num_species + 1][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 1][0] = 2;
                        
                        diffusivities_data[d_num_species + 2].resize(1);
                        diffusivities_component_idx[d_num_species + 2].resize(1);
                        
                        // 2/3*mu - mu_v.
                        diffusivities_data[d_num_species + 2][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 2][0] = 1;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[d_num_species + 3].resize(2);
                        diffusivities_component_idx[d_num_species + 3].resize(2);
                        
                        // w*(2/3*mu - mu_v).
                        diffusivities_data[d_num_species + 3][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][0] = 8;
                        
                        // -v*mu.
                        diffusivities_data[d_num_species + 3][1] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][1] = 10;
                        
                        break;
                    }
                    case DIRECTION::Z_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            diffusivities_data[si].resize(0);
                            diffusivities_component_idx[si].resize(0);
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[d_num_species].resize(1);
                        diffusivities_component_idx[d_num_species].resize(1);
                        
                        // -mu.
                        diffusivities_data[d_num_species][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species][0] = 2;
                        
                        diffusivities_data[d_num_species + 1].resize(1);
                        diffusivities_component_idx[d_num_species + 1].resize(1);
                        
                        // -mu.
                        diffusivities_data[d_num_species + 1][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 1][0] = 2;
                        
                        diffusivities_data[d_num_species + 2].resize(1);
                        diffusivities_component_idx[d_num_species + 2].resize(1);
                        
                        // -(4/3*mu + mu_v).
                        diffusivities_data[d_num_species + 2][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 2][0] = 0;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[d_num_species + 3].resize(3);
                        diffusivities_component_idx[d_num_species + 3].resize(3);
                        
                        // -u*mu.
                        diffusivities_data[d_num_species + 3][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][0] = 9;
                        
                        // -v*mu.
                        diffusivities_data[d_num_species + 3][1] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][1] = 10;
                        
                        // -w*(4/3*mu + mu_v).
                        diffusivities_data[d_num_species + 3][2] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][2] = 5;
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::"
                            << "getCellDataOfDiffusiveFluxDiffusivities()\n"
                            << "There are only x-direction, y-direction and z-direction for three-dimensional problem."
                            << std::endl);
                    }
                }
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::getCellDataOfDiffusiveFluxDiffusivities()\n"
                    << "There are only x-direction, y-direction and z-direction for three-dimensional problem."
                    << std::endl);
            }
        }
    }
}


/*
 * Get the cell data that needs interpolation to midpoints for computing side data of diffusivities in the
 * diffusive flux.
 */
void
FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::getCellDataForInterpolationToSideDataForDiffusiveFluxDiffusivities(
    std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& var_data_for_diffusivities,
    std::vector<int>& var_data_for_diffusivities_component_idx)
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
            << ": FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::"
            << "getCellDataOfDiffusiveFluxDiffusivities()\n"
            << "No patch is registered yet."
            << std::endl);
    }
    
    if (d_dim == tbox::Dimension(1))
    {
    }
    else if (d_dim == tbox::Dimension(2))
    {
    }
    else if (d_dim == tbox::Dimension(3))
    {
    }
}


/*
 * Compute the side data of the diffusivities in the diffusive flux with the interpolated side data.
 */
void
FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::computeSideDataOfDiffusiveFluxDiffusivities(
    const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& var_data_for_diffusivities)
{
    if (d_num_subghosts_diffusivities > -hier::IntVector::getOne(d_dim))
    {
        if (!d_side_data_computed_diffusivities)
        {
            if (d_flow_model.expired())
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "The object is not setup yet!"
                    << std::endl);
            }
            
            if (d_dim == tbox::Dimension(1))
            {
            }
            else if (d_dim == tbox::Dimension(2))
            {
            }
            else if (d_dim == tbox::Dimension(3))
            {
            }
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::computeSideDataOfDiffusiveFluxDiffusivities()\n"
            << "Cell data of 'DIFFUSIVITIES' is not yet registered."
            << std::endl);
    }
}


/*
 * Get the side data of the diffusivities in the diffusive fluxa.
 */
void
FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::getSideDataOfDiffusiveFluxDiffusivities(
    std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > >& diffusivities_data,
    std::vector<std::vector<int> >& diffusivities_component_idx,
    const DIRECTION::TYPE& flux_direction,
    const DIRECTION::TYPE& derivative_direction)
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
            << ": FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::"
            << "getSideDataOfDiffusiveFluxDiffusivities()\n"
            << "No patch is registered yet."
            << std::endl);
    }
    
    if (!d_cell_data_computed_diffusivities)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::getSideDataOfDiffusiveFluxDiffusivities()\n"
            << "Cell data of 'DIFFUSIVITIES' is not registered/computed yet."
            << std::endl);
    }
    
    diffusivities_data.resize(d_num_eqn);
    diffusivities_component_idx.resize(d_num_eqn);
    
    if (d_dim == tbox::Dimension(1))
    {
        switch (flux_direction)
        {
            case DIRECTION::X_DIRECTION:
            {
                switch (derivative_direction)
                {
                    case DIRECTION::X_DIRECTION:
                    {
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::"
                            << "getSideDataOfDiffusiveFluxDiffusivities()\n"
                            << "There are only x-direction for one-dimensional problem."
                            << std::endl);
                    }
                }
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::getSideDataOfDiffusiveFluxDiffusivities()\n"
                    << "There are only x-direction for one-dimensional problem."
                    << std::endl);
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        switch (flux_direction)
        {
            case DIRECTION::X_DIRECTION:
            {
                switch (derivative_direction)
                {
                    case DIRECTION::X_DIRECTION:
                    {
                        break;
                    }
                    case DIRECTION::Y_DIRECTION:
                    {
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::"
                            << "getSideDataOfDiffusiveFluxDiffusivities()\n"
                            << "There are only x-direction and y-direction for two-dimensional problem."
                            << std::endl);
                    }
                }
                
                break;
            }
            case DIRECTION::Y_DIRECTION:
            {
                switch (derivative_direction)
                {
                    case DIRECTION::X_DIRECTION:
                    {
                        break;
                    }
                    case DIRECTION::Y_DIRECTION:
                    {
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::"
                            << "getSideDataOfDiffusiveFluxDiffusivities()\n"
                            << "There are only x-direction and y-direction for two-dimensional problem."
                            << std::endl);
                    }
                }
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::getSideDataOfDiffusiveFluxDiffusivities()\n"
                    << "There are only x-direction and y-direction for two-dimensional problem."
                    << std::endl);
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        switch (flux_direction)
        {
            case DIRECTION::X_DIRECTION:
            {
                switch (derivative_direction)
                {
                    case DIRECTION::X_DIRECTION:
                    {
                        break;
                    }
                    case DIRECTION::Y_DIRECTION:
                    {
                        break;
                    }
                    case DIRECTION::Z_DIRECTION:
                    {
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::"
                            << "getSideDataOfDiffusiveFluxDiffusivities()\n"
                            << "There are only x-direction, y-direction and z-direction for three-dimensional problem."
                            << std::endl);
                    }
                }
                
                break;
            }
            case DIRECTION::Y_DIRECTION:
            {
                switch (derivative_direction)
                {
                    case DIRECTION::X_DIRECTION:
                    {
                        break;
                    }
                    case DIRECTION::Y_DIRECTION:
                    {
                        break;
                    }
                    case DIRECTION::Z_DIRECTION:
                    {
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::"
                            << "getSideDataOfDiffusiveFluxDiffusivities()\n"
                            << "There are only x-direction, y-direction and z-direction for three-dimensional problem."
                            << std::endl);
                    }
                }
                
                break;
            }
            case DIRECTION::Z_DIRECTION:
            {
                switch (derivative_direction)
                {
                    case DIRECTION::X_DIRECTION:
                    {
                        break;
                    }
                    case DIRECTION::Y_DIRECTION:
                    {
                        break;
                    }
                    case DIRECTION::Z_DIRECTION:
                    {
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::"
                            << "getSideDataOfDiffusiveFluxDiffusivities()\n"
                            << "There are only x-direction, y-direction and z-direction for three-dimensional problem."
                            << std::endl);
                    }
                }
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::getSideDataOfDiffusiveFluxDiffusivities()\n"
                    << "There are only x-direction, y-direction and z-direction for three-dimensional problem."
                    << std::endl);
            }
        }
    }
}


/*
 * Set the number of sub-ghost cells of a variable.
 * This function can be called recursively if the variables are computed recursively.
 */
void
FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::setNumberOfSubGhosts(
    const hier::IntVector& num_subghosts,
    const std::string& variable_name,
    const std::string& parent_variable_name)
{
    NULL_USE(parent_variable_name);
    
    if (variable_name == "SHEAR_VISCOSITY")
    {
        if (d_num_subghosts_shear_viscosity > -hier::IntVector::getOne(d_dim))
        {
            if (num_subghosts > d_num_subghosts_shear_viscosity)
            {
                d_num_subghosts_shear_viscosity = num_subghosts;
            }
        }
        else
        {
            d_num_subghosts_shear_viscosity = num_subghosts;
        }
    }
    
    if (variable_name == "BULK_VISCOSITY")
    {
        if (d_num_subghosts_bulk_viscosity > -hier::IntVector::getOne(d_dim))
        {
            if (num_subghosts > d_num_subghosts_bulk_viscosity)
            {
                d_num_subghosts_bulk_viscosity = num_subghosts;
            }
        }
        else
        {
            d_num_subghosts_bulk_viscosity = num_subghosts;
        }
    }
    
    if (variable_name == "DIFFUSIVITIES")
    {
        if (d_num_subghosts_diffusivities > -hier::IntVector::getOne(d_dim))
        {
            if (num_subghosts > d_num_subghosts_diffusivities)
            {
                d_num_subghosts_diffusivities = num_subghosts;
            }
        }
        else
        {
            d_num_subghosts_diffusivities = num_subghosts;
        }
        
        setNumberOfSubGhosts(num_subghosts, "SHEAR_VISCOSITY", parent_variable_name);
        setNumberOfSubGhosts(num_subghosts, "BULK_VISCOSITY", parent_variable_name);
    }
}


/*
 * Set the ghost boxes of derived cell variables.
 */
void
FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::setDerivedCellVariableGhostBoxes()
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
    
    if (d_num_subghosts_shear_viscosity > -hier::IntVector::getOne(d_dim))
    {
        d_subghost_box_shear_viscosity = interior_box;
        d_subghost_box_shear_viscosity.grow(d_num_subghosts_shear_viscosity);
        d_subghostcell_dims_shear_viscosity = d_subghost_box_shear_viscosity.numberCells();
    }
    
    if (d_num_subghosts_bulk_viscosity > -hier::IntVector::getOne(d_dim))
    {
        d_subghost_box_bulk_viscosity = interior_box;
        d_subghost_box_bulk_viscosity.grow(d_num_subghosts_bulk_viscosity);
        d_subghostcell_dims_bulk_viscosity = d_subghost_box_bulk_viscosity.numberCells();
    }
    
    if (d_num_subghosts_diffusivities > -hier::IntVector::getOne(d_dim))
    {
        d_subghost_box_diffusivities = interior_box;
        d_subghost_box_diffusivities.grow(d_num_subghosts_diffusivities);
        d_subghostcell_dims_diffusivities = d_subghost_box_diffusivities.numberCells();
    }
}


/*
 * Compute the cell data of shear viscosity in the registered patch.
 */
void
FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::computeCellDataOfShearViscosity()
{
    // Create empty box.
    const hier::Box empty_box(d_dim);
    
    if (d_num_subghosts_shear_viscosity > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_shear_viscosity)
        {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(d_data_shear_viscosity);
#endif
            
            if (d_flow_model.expired())
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "The object is not setup yet!"
                    << std::endl);
            }
            
            HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
            
            // Get the cell data of volume fractions.
            HAMERS_SHARED_PTR<pdat::CellData<double> > data_volume_fractions =
                flow_model_tmp->getCellData("VOLUME_FRACTIONS");
            
            // Get the cell data of mass fractions.
            HAMERS_SHARED_PTR<pdat::CellData<double> > data_mass_fractions =
                flow_model_tmp->getCellData("MASS_FRACTIONS");
            
            // Get the cell data of pressure.
            HAMERS_SHARED_PTR<pdat::CellData<double> > data_pressure =
                flow_model_tmp->getCellData("PRESSURE");
            
            // Get the cell data of species temperatures.
            std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > data_species_temperatures =
                flow_model_tmp->getSpeciesCellData("SPECIES_TEMPERATURES");
            
            // Compute the shear viscosity field.
            d_equation_of_shear_viscosity_mixing_rules->computeShearViscosity(
                d_data_shear_viscosity,
                data_pressure,
                data_species_temperatures,
                data_mass_fractions,
                data_volume_fractions,
                empty_box);
            
            d_cell_data_computed_shear_viscosity = true;
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::computeCellDataOfShearViscosity()\n"
            << "Cell data of 'SHEAR_VISCOSITY' is not yet registered."
            << std::endl);
    }
}


/*
 * Compute the cell data of bulk viscosity in the registered patch.
 */
void
FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::computeCellDataOfBulkViscosity()
{
    // Create empty box.
    const hier::Box empty_box(d_dim);
    
    if (d_num_subghosts_bulk_viscosity > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_bulk_viscosity)
        {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(d_data_bulk_viscosity);
#endif
            
            if (d_flow_model.expired())
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "The object is not setup yet!"
                    << std::endl);
            }
            
            HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
            
            // Get the cell data of volume fractions.
            HAMERS_SHARED_PTR<pdat::CellData<double> > data_volume_fractions =
                flow_model_tmp->getCellData("VOLUME_FRACTIONS");
            
            // Get the cell data of mass fractions.
            HAMERS_SHARED_PTR<pdat::CellData<double> > data_mass_fractions =
                flow_model_tmp->getCellData("MASS_FRACTIONS");
            
            // Get the cell data of pressure.
            HAMERS_SHARED_PTR<pdat::CellData<double> > data_pressure =
                flow_model_tmp->getCellData("PRESSURE");
            
            // Get the cell data of species temperatures.
            std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > data_species_temperatures =
                flow_model_tmp->getSpeciesCellData("SPECIES_TEMPERATURES");
            
            // Compute the bulk viscosity field.
            d_equation_of_bulk_viscosity_mixing_rules->computeBulkViscosity(
                d_data_bulk_viscosity,
                data_pressure,
                data_species_temperatures,
                data_mass_fractions,
                data_volume_fractions,
                empty_box);
            
            d_cell_data_computed_bulk_viscosity = true;
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::computeCellDataOfBulkViscosity()\n"
            << "Cell data of 'BULK_VISCOSITY' is not yet registered."
            << std::endl);
    }
}


/*
 * Compute the cell data of diffusivities in the registered patch.
 */
void
FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::computeCellDataOfDiffusivities()
{
    if (d_num_subghosts_diffusivities > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_diffusivities)
        {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(d_data_diffusivities);
#endif
            
            if (d_flow_model.expired())
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "The object is not setup yet!"
                    << std::endl);
            }
            
            HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
            const hier::Patch& patch = flow_model_tmp->getRegisteredPatch();
            
            /*
             * Get the dimension of the interior box.
             */
            
            const hier::Box interior_box = patch.getBox();
            const hier::IntVector interior_dims = interior_box.numberCells();
            
            if (!d_cell_data_computed_shear_viscosity)
            {
                computeCellDataOfShearViscosity();
            }
            
            if (!d_cell_data_computed_bulk_viscosity)
            {
                computeCellDataOfBulkViscosity();
            }
            
            // Get the cell data of velocity.
            HAMERS_SHARED_PTR<pdat::CellData<double> > data_velocity =
                flow_model_tmp->getCellData("VELOCITY");
            
            /*
             * Get the number of ghost cells of velocity.
             */
            
            const hier::IntVector num_subghosts_velocity = data_velocity->getGhostCellWidth();
            
            /*
             * Get the dimensions of the ghost cell box of velocity.
             */
            
            const hier::Box subghost_box_velocity = data_velocity->getGhostBox();
            const hier::IntVector subghostcell_dims_velocity = subghost_box_velocity.numberCells();
            
            /*
             * Get the pointers to the cell data of shear viscosity and bulk viscosity.
             */
            
            double* mu    = d_data_shear_viscosity->getPointer(0);
            double* mu_v  = d_data_bulk_viscosity->getPointer(0);
            
            if (d_dim == tbox::Dimension(1))
            {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(d_data_diffusivities->getDepth() == 2);
#endif
                
                // Get the pointer to cell data of velocity and diffusivities.
                double* u = data_velocity->getPointer(0);
                
                std::vector<double*> D_ptr;
                D_ptr.reserve(2);
                
                for (int i = 0; i < 2; i++)
                {
                    D_ptr.push_back(d_data_diffusivities->getPointer(i));
                }
                
                /*
                 * Compute the diffusivities.
                 */
                for (int i = -d_num_subghosts_diffusivities[0];
                     i < interior_dims[0] + d_num_subghosts_diffusivities[0];
                     i++)
                {
                    // Compute the linear indices.
                    const int idx_diffusivities = i + d_num_subghosts_diffusivities[0];
                    const int idx_shear_viscosity = i + d_num_subghosts_shear_viscosity[0];
                    const int idx_bulk_viscosity = i + d_num_subghosts_bulk_viscosity[0];
                    const int idx_velocity = i + num_subghosts_velocity[0];
                    
                    D_ptr[0][idx_diffusivities] =
                        -(double(4)/double(3)*mu[idx_shear_viscosity] + mu_v[idx_bulk_viscosity]);
                    D_ptr[1][idx_diffusivities] =
                        -u[idx_velocity]*(double(4)/double(3)*mu[idx_shear_viscosity] + mu_v[idx_bulk_viscosity]);
                }
            }
            else if (d_dim == tbox::Dimension(2))
            {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(d_data_diffusivities->getDepth() == 9);
#endif
                
                // Get the pointer to cell data of velocity and diffusivities.
                double* u = data_velocity->getPointer(0);
                double* v = data_velocity->getPointer(1);
                
                std::vector<double*> D_ptr;
                D_ptr.reserve(9);
                
                for (int i = 0; i < 9; i++)
                {
                    D_ptr.push_back(d_data_diffusivities->getPointer(i));
                }
                
                /*
                 * Compute the diffusivities.
                 */
                for (int j = -d_num_subghosts_diffusivities[1];
                     j < interior_dims[1] + d_num_subghosts_diffusivities[1];
                     j++)
                {
                    for (int i = -d_num_subghosts_diffusivities[0];
                         i < interior_dims[0] + d_num_subghosts_diffusivities[0];
                         i++)
                    {
                        // Compute the linear indices.
                        const int idx_diffusivities = (i + d_num_subghosts_diffusivities[0]) +
                            (j + d_num_subghosts_diffusivities[1])*d_subghostcell_dims_diffusivities[0];
                        
                        const int idx_shear_viscosity = (i + d_num_subghosts_shear_viscosity[0]) +
                            (j + d_num_subghosts_shear_viscosity[1])*d_subghostcell_dims_shear_viscosity[0];
                        
                        const int idx_bulk_viscosity = (i + d_num_subghosts_bulk_viscosity[0]) +
                            (j + d_num_subghosts_bulk_viscosity[1])*d_subghostcell_dims_bulk_viscosity[0];
                        
                        const int idx_velocity = (i + num_subghosts_velocity[0]) +
                            (j + num_subghosts_velocity[1])*subghostcell_dims_velocity[0];
                        
                        D_ptr[0][idx_diffusivities] =
                            -(double(4)/double(3)*mu[idx_shear_viscosity] + mu_v[idx_bulk_viscosity]);
                        D_ptr[1][idx_diffusivities] =
                            double(2)/double(3)*mu[idx_shear_viscosity] - mu_v[idx_bulk_viscosity];
                        D_ptr[2][idx_diffusivities] =
                            -mu[idx_shear_viscosity];
                        D_ptr[3][idx_diffusivities] =
                            -u[idx_velocity]*(double(4)/double(3)*mu[idx_shear_viscosity] + mu_v[idx_bulk_viscosity]);
                        D_ptr[4][idx_diffusivities] =
                            -v[idx_velocity]*(double(4)/double(3)*mu[idx_shear_viscosity] + mu_v[idx_bulk_viscosity]);
                        D_ptr[5][idx_diffusivities] =
                            u[idx_velocity]*(double(2)/double(3)*mu[idx_shear_viscosity] - mu_v[idx_bulk_viscosity]);
                        D_ptr[6][idx_diffusivities] =
                            v[idx_velocity]*(double(2)/double(3)*mu[idx_shear_viscosity] - mu_v[idx_bulk_viscosity]);
                        D_ptr[7][idx_diffusivities] =
                            -u[idx_velocity]*mu[idx_shear_viscosity];
                        D_ptr[8][idx_diffusivities] =
                            -v[idx_velocity]*mu[idx_shear_viscosity];
                    }
                }
            }
            else if (d_dim == tbox::Dimension(3))
            {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(d_data_diffusivities->getDepth() == 12);
#endif
                
                // Get the pointer to cell data of velocity and diffusivities.
                double* u = data_velocity->getPointer(0);
                double* v = data_velocity->getPointer(1);
                double* w = data_velocity->getPointer(2);
                
                std::vector<double*> D_ptr;
                D_ptr.reserve(12);
                
                for (int i = 0; i < 12; i++)
                {
                    D_ptr.push_back(d_data_diffusivities->getPointer(i));
                }
                
                /*
                 * Compute the diffusivities.
                 */
                for (int k = -d_num_subghosts_diffusivities[2];
                     k < interior_dims[2] + d_num_subghosts_diffusivities[2];
                     k++)
                {
                    for (int j = -d_num_subghosts_diffusivities[1];
                         j < interior_dims[1] + d_num_subghosts_diffusivities[1];
                         j++)
                    {
                        for (int i = -d_num_subghosts_diffusivities[0];
                             i < interior_dims[0] + d_num_subghosts_diffusivities[0];
                             i++)
                        {
                            const int idx_diffusivities = (i + d_num_subghosts_diffusivities[0]) +
                                (j + d_num_subghosts_diffusivities[1])*d_subghostcell_dims_diffusivities[0] +
                                (k + d_num_subghosts_diffusivities[2])*d_subghostcell_dims_diffusivities[0]*
                                    d_subghostcell_dims_diffusivities[1];
                            
                            const int idx_shear_viscosity = (i + d_num_subghosts_shear_viscosity[0]) +
                                (j + d_num_subghosts_shear_viscosity[1])*d_subghostcell_dims_shear_viscosity[0] +
                                (k + d_num_subghosts_shear_viscosity[2])*d_subghostcell_dims_shear_viscosity[0]*
                                    d_subghostcell_dims_shear_viscosity[1];
                            
                            const int idx_bulk_viscosity = (i + d_num_subghosts_bulk_viscosity[0]) +
                                (j + d_num_subghosts_bulk_viscosity[1])*d_subghostcell_dims_bulk_viscosity[0] +
                                (k + d_num_subghosts_bulk_viscosity[2])*d_subghostcell_dims_bulk_viscosity[0]*
                                    d_subghostcell_dims_bulk_viscosity[1];
                            
                            const int idx_velocity = (i + num_subghosts_velocity[0]) +
                                (j + num_subghosts_velocity[1])*subghostcell_dims_velocity[0] +
                                (k + num_subghosts_velocity[2])*subghostcell_dims_velocity[0]*
                                    subghostcell_dims_velocity[1];
                            
                            D_ptr[0][idx_diffusivities] =
                                -(double(4)/double(3)*mu[idx_shear_viscosity] + mu_v[idx_bulk_viscosity]);
                            D_ptr[1][idx_diffusivities] =
                                double(2)/double(3)*mu[idx_shear_viscosity] - mu_v[idx_bulk_viscosity];
                            D_ptr[2][idx_diffusivities] =
                                -mu[idx_shear_viscosity];
                            D_ptr[3][idx_diffusivities] =
                                -u[idx_velocity]*(double(4)/double(3)*mu[idx_shear_viscosity] + mu_v[idx_bulk_viscosity]);
                            D_ptr[4][idx_diffusivities] =
                                -v[idx_velocity]*(double(4)/double(3)*mu[idx_shear_viscosity] + mu_v[idx_bulk_viscosity]);
                            D_ptr[5][idx_diffusivities] =
                                -w[idx_velocity]*(double(4)/double(3)*mu[idx_shear_viscosity] + mu_v[idx_bulk_viscosity]);
                            D_ptr[6][idx_diffusivities] =
                                u[idx_velocity]*(double(2)/double(3)*mu[idx_shear_viscosity] - mu_v[idx_bulk_viscosity]);
                            D_ptr[7][idx_diffusivities] =
                                v[idx_velocity]*(double(2)/double(3)*mu[idx_shear_viscosity] - mu_v[idx_bulk_viscosity]);
                            D_ptr[8][idx_diffusivities] =
                                w[idx_velocity]*(double(2)/double(3)*mu[idx_shear_viscosity] - mu_v[idx_bulk_viscosity]);
                            D_ptr[9][idx_diffusivities] =
                                -u[idx_velocity]*mu[idx_shear_viscosity];
                            D_ptr[10][idx_diffusivities] =
                                -v[idx_velocity]*mu[idx_shear_viscosity];
                            D_ptr[11][idx_diffusivities] =
                                -w[idx_velocity]*mu[idx_shear_viscosity];
                        }
                    }
                }
            }
            
            d_cell_data_computed_diffusivities = true;
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::computeCellDataOfDiffusivities()\n"
            << "Cell data of 'DIFFUSIVITIES' is not yet registered."
            << std::endl);
    }
}
