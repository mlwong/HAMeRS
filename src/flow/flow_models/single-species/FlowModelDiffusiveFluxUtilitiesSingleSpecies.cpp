#include "flow/flow_models/single-species/FlowModelDiffusiveFluxUtilitiesSingleSpecies.hpp"

#include "flow/flow_models/single-species/FlowModelSubgridScaleModelSingleSpecies.hpp"

FlowModelDiffusiveFluxUtilitiesSingleSpecies::FlowModelDiffusiveFluxUtilitiesSingleSpecies(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
    const int& num_species,
    const HAMERS_SHARED_PTR<tbox::Database>& flow_model_db,
    const HAMERS_SHARED_PTR<EquationOfShearViscosityMixingRules> equation_of_shear_viscosity_mixing_rules,
    const HAMERS_SHARED_PTR<EquationOfBulkViscosityMixingRules> equation_of_bulk_viscosity_mixing_rules,
    const HAMERS_SHARED_PTR<EquationOfThermalConductivityMixingRules> equation_of_thermal_conductivity_mixing_rules):
        FlowModelDiffusiveFluxUtilities(
            object_name,
            dim,
            grid_geometry,
            num_species,
            2 + dim.getValue(),
            flow_model_db),
        d_num_subghosts_shear_viscosity(-hier::IntVector::getOne(d_dim)),
        d_num_subghosts_bulk_viscosity(-hier::IntVector::getOne(d_dim)),
        d_num_subghosts_thermal_conductivity(-hier::IntVector::getOne(d_dim)),
        d_subghost_box_shear_viscosity(hier::Box::getEmptyBox(d_dim)),
        d_subghost_box_bulk_viscosity(hier::Box::getEmptyBox(d_dim)),
        d_subghost_box_thermal_conductivity(hier::Box::getEmptyBox(d_dim)),
        d_subghostcell_dims_shear_viscosity(hier::IntVector::getZero(d_dim)),
        d_subghostcell_dims_bulk_viscosity(hier::IntVector::getZero(d_dim)),
        d_subghostcell_dims_thermal_conductivity(hier::IntVector::getZero(d_dim)),
        d_cell_data_computed_shear_viscosity(false),
        d_cell_data_computed_bulk_viscosity(false),
        d_cell_data_computed_thermal_conductivity(false),
        d_equation_of_shear_viscosity_mixing_rules(equation_of_shear_viscosity_mixing_rules),
        d_equation_of_bulk_viscosity_mixing_rules(equation_of_bulk_viscosity_mixing_rules),
        d_equation_of_thermal_conductivity_mixing_rules(equation_of_thermal_conductivity_mixing_rules)
{
    if (d_use_subgrid_scale_model)
    {
        HAMERS_SHARED_PTR<tbox::Database> subgrid_scale_model_db;
        
        if (flow_model_db->keyExists("Subgrid_scale_model"))
        {
            subgrid_scale_model_db = flow_model_db->getDatabase("Subgrid_scale_model");
        }
        else if (flow_model_db->keyExists("d_subgrid_scale_model"))
        {
            subgrid_scale_model_db = flow_model_db->getDatabase("d_subgrid_scale_model");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "No key 'Subgrid_scale_model'/'d_subgrid_scale_model' found in data for flow model."
                << std::endl);
        }
        
        /*
         * Initialize subgrid scale model object.
         */
        d_flow_model_subgrid_scale_model.reset(new FlowModelSubgridScaleModelSingleSpecies(
            "d_flow_model_subgrid_scale_model",
            d_dim,
            d_grid_geometry,
            d_num_species,
            subgrid_scale_model_db));
    }
}


/*
 * Register different derived variables related to this class in the registered patch. The
 * derived variables to be registered are given as entries in a map of the variable name to
 * the number of sub-ghost cells required.
 */
void
FlowModelDiffusiveFluxUtilitiesSingleSpecies::registerDerivedVariables(
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
            << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::"
            << "registerDerivedVariables()\n"
            << "No patch is registered yet."
            << std::endl);
    }
    
    // Check whether all or part of derived cell data is alredy computed.
    if (d_derived_cell_data_computed)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::registerDerivedVariables()\n"
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
                << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::registerDerivedVariables()\n"
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
            "PRESSURE",
            num_subghosts_of_data.find("SHEAR_VISCOSITY")->second));
        
        num_subghosts_of_data_flow_model.insert(std::pair<std::string, hier::IntVector>(
            "TEMPERATURE",
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
            "PRESSURE",
            num_subghosts_of_data.find("BULK_VISCOSITY")->second));
        
        num_subghosts_of_data_flow_model.insert(std::pair<std::string, hier::IntVector>(
            "TEMPERATURE",
            num_subghosts_of_data.find("BULK_VISCOSITY")->second));
        
        flow_model_tmp->registerDerivedVariables(num_subghosts_of_data_flow_model);
        
        setNumberOfSubGhosts(
            num_subghosts_of_data.find("BULK_VISCOSITY")->second,
            "BULK_VISCOSITY",
            "BULK_VISCOSITY");
    }
    
    if (num_subghosts_of_data.find("THERMAL_CONDUCTIVITY") != num_subghosts_of_data.end())
    {
        std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data_flow_model;
        
        num_subghosts_of_data_flow_model.insert(std::pair<std::string, hier::IntVector>(
            "PRESSURE",
            num_subghosts_of_data.find("THERMAL_CONDUCTIVITY")->second));
        
        num_subghosts_of_data_flow_model.insert(std::pair<std::string, hier::IntVector>(
            "TEMPERATURE",
            num_subghosts_of_data.find("THERMAL_CONDUCTIVITY")->second));
        
        flow_model_tmp->registerDerivedVariables(num_subghosts_of_data_flow_model);
        
        setNumberOfSubGhosts(
            num_subghosts_of_data.find("THERMAL_CONDUCTIVITY")->second,
            "THERMAL_CONDUCTIVITY",
            "THERMAL_CONDUCTIVITY");
    }
}


/*
 * Register the required variables for the computation of diffusive fluxes in the registered patch.
 */
void
FlowModelDiffusiveFluxUtilitiesSingleSpecies::registerDerivedVariablesForDiffusiveFluxes(
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
            << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::"
            << "registerDerivedVariablesForDiffusiveFluxes()\n"
            << "No patch is registered yet."
            << std::endl);
    }
    
    // Check whether all or part of derived cell data is alredy computed.
    if (d_derived_cell_data_computed)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::registerDerivedVariablesForDiffusiveFluxes()\n"
            << "Derived cell data is already computed."
            << std::endl);
    }
    
    if ((num_subghosts < hier::IntVector::getZero(d_dim)) ||
        (num_subghosts > flow_model_tmp->getNumberOfGhostCells()))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::registerDerivedVariablesForDiffusiveFluxes()\n"
            << "The number of sub-ghost cells of variable is not between zero and number of ghosts of conservative variables."
            << std::endl);
    }
    
    /*
     * Register the required derived variables in flow model.
     */
    
    std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
    
    num_subghosts_of_data.insert(
        std::pair<std::string, hier::IntVector>("VELOCITY", num_subghosts));
    
    num_subghosts_of_data.insert(
        std::pair<std::string, hier::IntVector>("PRESSURE", num_subghosts));
    
    num_subghosts_of_data.insert(
        std::pair<std::string, hier::IntVector>("TEMPERATURE", num_subghosts));
    
    flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
    
    setNumberOfSubGhosts(
        num_subghosts,
        "DIFFUSIVITIES",
        "DIFFUSIVITIES");
    
    d_need_side_diffusivities = need_side_diffusivities;
    
    if (d_use_subgrid_scale_model)
    {
        const std::vector<std::string>& var_to_register = d_flow_model_subgrid_scale_model->getDerivedVariablesToRegister();
        
        for (int vi = 0; vi < static_cast<int>(var_to_register.size()); vi++)
        {
            num_subghosts_of_data.insert(
                std::pair<std::string, hier::IntVector>(var_to_register[vi], num_subghosts));
        }
    }
}


/*
 * Allocate memory for cell data of different registered derived variables related to this
 * class in the registered patch.
 */
void
FlowModelDiffusiveFluxUtilitiesSingleSpecies::allocateMemoryForDerivedCellData()
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
                d_data_shear_viscosity.reset(new pdat::CellData<Real>(
                    interior_box, 1, d_num_subghosts_shear_viscosity));
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::allocateMemoryForDerivedCellData()\n"
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
                d_data_bulk_viscosity.reset(new pdat::CellData<Real>(
                    interior_box, 1, d_num_subghosts_bulk_viscosity));
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::allocateMemoryForDerivedCellData()\n"
                << "Cell data of 'BULK_VISCOSITY' is aleady computed."
                << std::endl);
        }
    }
    
    if (d_num_subghosts_thermal_conductivity > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_thermal_conductivity)
        {
            if (!d_data_thermal_conductivity)
            {
                d_data_thermal_conductivity.reset(new pdat::CellData<Real>(
                    interior_box, 1, d_num_subghosts_thermal_conductivity));
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::allocateMemoryForDerivedCellData()\n"
                << "Cell data of 'THERMAL_CONDUCTIVITY' is aleady computed."
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
                        d_data_diffusivities.reset(new pdat::CellData<Real>(
                            interior_box, 3, d_num_subghosts_diffusivities));
                    }
                    else if (d_dim == tbox::Dimension(2))
                    {
                        d_data_diffusivities.reset(new pdat::CellData<Real>(
                            interior_box, 10, d_num_subghosts_diffusivities));
                    }
                    else if (d_dim == tbox::Dimension(3))
                    {
                        d_data_diffusivities.reset(new pdat::CellData<Real>(
                            interior_box, 13, d_num_subghosts_diffusivities));
                    }
                }
            }
            else
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::allocateMemoryForDerivedCellData()\n"
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
FlowModelDiffusiveFluxUtilitiesSingleSpecies::allocateMemoryForSideDataOfDiffusiveFluxDiffusivities()
{
    if (!d_need_side_diffusivities)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::allocateMemoryForSideDataOfDiffusiveFluxDiffusivities()\n"
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
    
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
    const hier::Patch& patch = flow_model_tmp->getRegisteredPatch();
    const hier::Box interior_box = patch.getBox();
    
    if (d_num_subghosts_diffusivities > -hier::IntVector::getOne(d_dim))
    {
        if (!d_side_data_diffusivities_computed)
        {
            if (!d_side_data_diffusivities)
            {
                if (d_dim == tbox::Dimension(1))
                {
                    d_side_data_diffusivities.reset(new pdat::SideData<Real>(
                        interior_box, 3, d_num_subghosts_diffusivities));
                }
                else if (d_dim == tbox::Dimension(2))
                {
                    d_side_data_diffusivities.reset(new pdat::SideData<Real>(
                        interior_box, 7, d_num_subghosts_diffusivities));
                }
                else if (d_dim == tbox::Dimension(3))
                {
                    d_side_data_diffusivities.reset(new pdat::SideData<Real>(
                        interior_box, 8, d_num_subghosts_diffusivities));
                }
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::allocateMemoryForSideDataOfDiffusiveFluxDiffusivities()\n"
                << "Side data of 'DIFFUSIVITIES' is aleady computed."
                << std::endl);
        }
    }
}


/*
 * Clear cell and side data of different derived variables related to this class in the registered patch.
 */
void
FlowModelDiffusiveFluxUtilitiesSingleSpecies::clearCellAndSideData()
{
    d_num_subghosts_diffusivities        = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_shear_viscosity      = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_bulk_viscosity       = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_thermal_conductivity = -hier::IntVector::getOne(d_dim);
    
    d_subghost_box_diffusivities        = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_shear_viscosity      = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_bulk_viscosity       = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_thermal_conductivity = hier::Box::getEmptyBox(d_dim);
    
    d_subghostcell_dims_diffusivities        = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_shear_viscosity      = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_bulk_viscosity       = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_thermal_conductivity = hier::IntVector::getZero(d_dim);
    
    d_data_diffusivities.reset();
    d_data_shear_viscosity.reset();
    d_data_bulk_viscosity.reset();
    d_data_thermal_conductivity.reset();
    
    d_side_data_diffusivities.reset();
    
    d_cell_data_computed_diffusivities        = false;
    d_cell_data_computed_shear_viscosity      = false;
    d_cell_data_computed_bulk_viscosity       = false;
    d_cell_data_computed_thermal_conductivity = false;
    
    d_derived_cell_data_computed = false;
    
    d_side_data_diffusivities_computed = false;
}


/*
 * Compute cell data of different registered derived variables related to this class.
 */
void
FlowModelDiffusiveFluxUtilitiesSingleSpecies::computeDerivedCellData()
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
            << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::"
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
    
    // Compute the thermal conductivity cell data.
    if (d_num_subghosts_thermal_conductivity > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_thermal_conductivity)
        {
            computeCellDataOfThermalConductivity();
        }
    }
    
    if (!d_need_side_diffusivities)
    {
        // Compute the diffusivities cell data.
        if (d_num_subghosts_diffusivities > -hier::IntVector::getOne(d_dim))
        {
            if (!d_cell_data_computed_diffusivities)
            {
                computeCellDataOfDiffusivities();
            }
        }
    }
    
    d_derived_cell_data_computed = true;
}


/*
 * Get the cell data of one cell variable related to this class in the registered patch.
 */
HAMERS_SHARED_PTR<pdat::CellData<Real> >
FlowModelDiffusiveFluxUtilitiesSingleSpecies::getCellData(const std::string& variable_key)
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
            << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::"
            << "getCellData()\n"
            << "No patch is registered yet."
            << std::endl);
    }
    
    HAMERS_SHARED_PTR<pdat::CellData<Real> > cell_data;
    
    if (variable_key == "SHEAR_VISCOSITY")
    {
        if (!d_cell_data_computed_shear_viscosity)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::getCellData()\n"
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
                << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::getCellData()\n"
                << "Cell data of 'BULK_VISCOSITY' is not registered/computed yet."
                << std::endl);
        }
        cell_data = d_data_bulk_viscosity;
    }
    else if (variable_key == "THERMAL_CONDUCTIVITY")
    {
        if (!d_cell_data_computed_thermal_conductivity)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::getCellData()\n"
                << "Cell data of 'THERMAL_CONDUCTIVITY' is not registered/computed yet."
                << std::endl);
        }
        cell_data = d_data_thermal_conductivity;
    }
    
    return cell_data;
}


/*
 * Get the cell data of different cell variables related to this class in the registered patch.
 */
std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >
FlowModelDiffusiveFluxUtilitiesSingleSpecies::getCellData(
    const std::vector<std::string>& variable_keys)
{
    std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > cell_data(
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
FlowModelDiffusiveFluxUtilitiesSingleSpecies::getCellDataOfDiffusiveFluxVariablesForDerivative(
    std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > >& derivative_var_data,
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
    HAMERS_SHARED_PTR<pdat::CellData<Real> > data_velocity =
        flow_model_tmp->getCellData("VELOCITY");
    
    // Get the cell data of temperature.
    HAMERS_SHARED_PTR<pdat::CellData<Real> > data_temperature =
        flow_model_tmp->getCellData("TEMPERATURE");
    
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
                         * Mass equation.
                         */
                        
                        derivative_var_data[0].resize(0);
                        derivative_var_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[1].resize(1);
                        derivative_var_component_idx[1].resize(1);
                        
                        // Variable u.
                        derivative_var_data[1][0] = data_velocity;
                        derivative_var_component_idx[1][0] = 0;
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[2].resize(2);
                        derivative_var_component_idx[2].resize(2);
                        
                        // Variable u.
                        derivative_var_data[2][0] = data_velocity;
                        derivative_var_component_idx[2][0] = 0;
                        
                        // Variable T.
                        derivative_var_data[2][1] = data_temperature;
                        derivative_var_component_idx[2][1] = 0;
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::"
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
                    << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::"
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
                         * Mass equation.
                         */
                        
                        derivative_var_data[0].resize(0);
                        derivative_var_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[1].resize(1);
                        derivative_var_component_idx[1].resize(1);
                        
                        // Variable u.
                        derivative_var_data[1][0] = data_velocity;
                        derivative_var_component_idx[1][0] = 0;
                        
                        derivative_var_data[2].resize(1);
                        derivative_var_component_idx[2].resize(1);
                        
                        // Variable v.
                        derivative_var_data[2][0] = data_velocity;
                        derivative_var_component_idx[2][0] = 1;
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[3].resize(3);
                        derivative_var_component_idx[3].resize(3);
                        
                        // Variable u.
                        derivative_var_data[3][0] = data_velocity;
                        derivative_var_component_idx[3][0] = 0;
                        
                        // Variable v.
                        derivative_var_data[3][1] = data_velocity;
                        derivative_var_component_idx[3][1] = 1;
                        
                        // Variable T.
                        derivative_var_data[3][2] = data_temperature;
                        derivative_var_component_idx[3][2] = 0;
                        
                        break;
                    }
                    case DIRECTION::Y_DIRECTION:
                    {
                        /*
                         * Mass equation.
                         */
                        
                        derivative_var_data[0].resize(0);
                        derivative_var_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[1].resize(1);
                        derivative_var_component_idx[1].resize(1);
                        
                        // Variable v.
                        derivative_var_data[1][0] = data_velocity;
                        derivative_var_component_idx[1][0] = 1;
                        
                        derivative_var_data[2].resize(1);
                        derivative_var_component_idx[2].resize(1);
                        
                        // Variable u.
                        derivative_var_data[2][0] = data_velocity;
                        derivative_var_component_idx[2][0] = 0;
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[3].resize(2);
                        derivative_var_component_idx[3].resize(2);
                        
                        // Variable u.
                        derivative_var_data[3][0] = data_velocity;
                        derivative_var_component_idx[3][0] = 0;
                        
                        // Variable v.
                        derivative_var_data[3][1] = data_velocity;
                        derivative_var_component_idx[3][1] = 1;
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::"
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
                         * Mass equation.
                         */
                        
                        derivative_var_data[0].resize(0);
                        derivative_var_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[1].resize(1);
                        derivative_var_component_idx[1].resize(1);
                        
                        // Variable v.
                        derivative_var_data[1][0] = data_velocity;
                        derivative_var_component_idx[1][0] = 1;
                        
                        derivative_var_data[2].resize(1);
                        derivative_var_component_idx[2].resize(1);
                        
                        // Variable u.
                        derivative_var_data[2][0] = data_velocity;
                        derivative_var_component_idx[2][0] = 0;
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[3].resize(2);
                        derivative_var_component_idx[3].resize(2);
                        
                        // Variable u.
                        derivative_var_data[3][0] = data_velocity;
                        derivative_var_component_idx[3][0] = 0;
                        
                        // Variable v.
                        derivative_var_data[3][1] = data_velocity;
                        derivative_var_component_idx[3][1] = 1;
                        
                        break;
                    }
                    case DIRECTION::Y_DIRECTION:
                    {
                        /*
                         * Mass equation.
                         */
                        
                        derivative_var_data[0].resize(0);
                        derivative_var_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[1].resize(1);
                        derivative_var_component_idx[1].resize(1);
                        
                        // Variable u.
                        derivative_var_data[1][0] = data_velocity;
                        derivative_var_component_idx[1][0] = 0;
                        
                        derivative_var_data[2].resize(1);
                        derivative_var_component_idx[2].resize(1);
                        
                        // Variable v.
                        derivative_var_data[2][0] = data_velocity;
                        derivative_var_component_idx[2][0] = 1;
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[3].resize(3);
                        derivative_var_component_idx[3].resize(3);
                        
                        // Variable u.
                        derivative_var_data[3][0] = data_velocity;
                        derivative_var_component_idx[3][0] = 0;
                        
                        // Variable v.
                        derivative_var_data[3][1] = data_velocity;
                        derivative_var_component_idx[3][1] = 1;
                        
                        // Variable T.
                        derivative_var_data[3][2] = data_temperature;
                        derivative_var_component_idx[3][2] = 0;
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::"
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
                    << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::"
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
                         * Mass equation.
                         */
                        
                        derivative_var_data[0].resize(0);
                        derivative_var_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[1].resize(1);
                        derivative_var_component_idx[1].resize(1);
                        
                        // Variable u.
                        derivative_var_data[1][0] = data_velocity;
                        derivative_var_component_idx[1][0] = 0;
                        
                        derivative_var_data[2].resize(1);
                        derivative_var_component_idx[2].resize(1);
                        
                        // Variable v.
                        derivative_var_data[2][0] = data_velocity;
                        derivative_var_component_idx[2][0] = 1;
                        
                        derivative_var_data[3].resize(1);
                        derivative_var_component_idx[3].resize(1);
                        
                        // Variable w.
                        derivative_var_data[3][0] = data_velocity;
                        derivative_var_component_idx[3][0] = 2;
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[4].resize(4);
                        derivative_var_component_idx[4].resize(4);
                        
                        // Variable u.
                        derivative_var_data[4][0] = data_velocity;
                        derivative_var_component_idx[4][0] = 0;
                        
                        // Variable v.
                        derivative_var_data[4][1] = data_velocity;
                        derivative_var_component_idx[4][1] = 1;
                        
                        // Variable w.
                        derivative_var_data[4][2] = data_velocity;
                        derivative_var_component_idx[4][2] = 2;
                        
                        // Variable T.
                        derivative_var_data[4][3] = data_temperature;
                        derivative_var_component_idx[4][3] = 0;
                        
                        break;
                    }
                    case DIRECTION::Y_DIRECTION:
                    {
                        /*
                         * Mass equation.
                         */
                        
                        derivative_var_data[0].resize(0);
                        derivative_var_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[1].resize(1);
                        derivative_var_component_idx[1].resize(1);
                        
                        // Variable v.
                        derivative_var_data[1][0] = data_velocity;
                        derivative_var_component_idx[1][0] = 1;
                        
                        derivative_var_data[2].resize(1);
                        derivative_var_component_idx[2].resize(1);
                        
                        // Variable u.
                        derivative_var_data[2][0] = data_velocity;
                        derivative_var_component_idx[2][0] = 0;
                        
                        derivative_var_data[3].resize(0);
                        derivative_var_component_idx[3].resize(0);
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[4].resize(2);
                        derivative_var_component_idx[4].resize(2);
                        
                        // Variable u.
                        derivative_var_data[4][0] = data_velocity;
                        derivative_var_component_idx[4][0] = 0;
                        
                        // Variable v.
                        derivative_var_data[4][1] = data_velocity;
                        derivative_var_component_idx[4][1] = 1;
                        
                        break;
                    }
                    case DIRECTION::Z_DIRECTION:
                    {
                        /*
                         * Mass equation.
                         */
                        
                        derivative_var_data[0].resize(0);
                        derivative_var_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[1].resize(1);
                        derivative_var_component_idx[1].resize(1);
                        
                        // Variable w.
                        derivative_var_data[1][0] = data_velocity;
                        derivative_var_component_idx[1][0] = 2;
                        
                        derivative_var_data[2].resize(0);
                        derivative_var_component_idx[2].resize(0);
                        
                        derivative_var_data[3].resize(1);
                        derivative_var_component_idx[3].resize(1);
                        
                        // Variable u.
                        derivative_var_data[3][0] = data_velocity;
                        derivative_var_component_idx[3][0] = 0;
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[4].resize(2);
                        derivative_var_component_idx[4].resize(2);
                        
                        // Variable u.
                        derivative_var_data[4][0] = data_velocity;
                        derivative_var_component_idx[4][0] = 0;
                        
                        // Variable w.
                        derivative_var_data[4][1] = data_velocity;
                        derivative_var_component_idx[4][1] = 2;
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::"
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
                         * Mass equation.
                         */
                        
                        derivative_var_data[0].resize(0);
                        derivative_var_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[1].resize(1);
                        derivative_var_component_idx[1].resize(1);
                        
                        // Variable v.
                        derivative_var_data[1][0] = data_velocity;
                        derivative_var_component_idx[1][0] = 1;
                        
                        derivative_var_data[2].resize(1);
                        derivative_var_component_idx[2].resize(1);
                        
                        // Variable u.
                        derivative_var_data[2][0] = data_velocity;
                        derivative_var_component_idx[2][0] = 0;
                        
                        derivative_var_data[3].resize(0);
                        derivative_var_component_idx[3].resize(0);
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[4].resize(2);
                        derivative_var_component_idx[4].resize(2);
                        
                        // Variable u.
                        derivative_var_data[4][0] = data_velocity;
                        derivative_var_component_idx[4][0] = 0;
                        
                        // Variable v.
                        derivative_var_data[4][1] = data_velocity;
                        derivative_var_component_idx[4][1] = 1;
                        
                        break;
                    }
                    case DIRECTION::Y_DIRECTION:
                    {
                        /*
                         * Mass equation.
                         */
                        
                        derivative_var_data[0].resize(0);
                        derivative_var_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[1].resize(1);
                        derivative_var_component_idx[1].resize(1);
                        
                        // Variable u.
                        derivative_var_data[1][0] = data_velocity;
                        derivative_var_component_idx[1][0] = 0;
                        
                        derivative_var_data[2].resize(1);
                        derivative_var_component_idx[2].resize(1);
                        
                        // Variable v.
                        derivative_var_data[2][0] = data_velocity;
                        derivative_var_component_idx[2][0] = 1;
                        
                        derivative_var_data[3].resize(1);
                        derivative_var_component_idx[3].resize(1);
                        
                        // Variable w.
                        derivative_var_data[3][0] = data_velocity;
                        derivative_var_component_idx[3][0] = 2;
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[4].resize(4);
                        derivative_var_component_idx[4].resize(4);
                        
                        // Variable u.
                        derivative_var_data[4][0] = data_velocity;
                        derivative_var_component_idx[4][0] = 0;
                        
                        // Variable v.
                        derivative_var_data[4][1] = data_velocity;
                        derivative_var_component_idx[4][1] = 1;
                        
                        // Variable w.
                        derivative_var_data[4][2] = data_velocity;
                        derivative_var_component_idx[4][2] = 2;
                        
                        // Variable T.
                        derivative_var_data[4][3] = data_temperature;
                        derivative_var_component_idx[4][3] = 0;
                        
                        break;
                    }
                    case DIRECTION::Z_DIRECTION:
                    {
                        /*
                         * Mass equation.
                         */
                        
                        derivative_var_data[0].resize(0);
                        derivative_var_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[1].resize(0);
                        derivative_var_component_idx[1].resize(0);
                        
                        derivative_var_data[2].resize(1);
                        derivative_var_component_idx[2].resize(1);
                        
                        // Variable w.
                        derivative_var_data[2][0] = data_velocity;
                        derivative_var_component_idx[2][0] = 2;
                        
                        derivative_var_data[3].resize(1);
                        derivative_var_component_idx[3].resize(1);
                        
                        // Variable v.
                        derivative_var_data[3][0] = data_velocity;
                        derivative_var_component_idx[3][0] = 1;
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[4].resize(2);
                        derivative_var_component_idx[4].resize(2);
                        
                        // Variable v.
                        derivative_var_data[4][0] = data_velocity;
                        derivative_var_component_idx[4][0] = 1;
                        
                        // Variable w.
                        derivative_var_data[4][1] = data_velocity;
                        derivative_var_component_idx[4][1] = 2;
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::"
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
                         * Mass equation.
                         */
                        
                        derivative_var_data[0].resize(0);
                        derivative_var_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[1].resize(1);
                        derivative_var_component_idx[1].resize(1);
                        
                        // Variable w.
                        derivative_var_data[1][0] = data_velocity;
                        derivative_var_component_idx[1][0] = 2;
                        
                        derivative_var_data[2].resize(0);
                        derivative_var_component_idx[2].resize(0);
                        
                        derivative_var_data[3].resize(1);
                        derivative_var_component_idx[3].resize(1);
                        
                        // Variable u.
                        derivative_var_data[3][0] = data_velocity;
                        derivative_var_component_idx[3][0] = 0;
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[4].resize(2);
                        derivative_var_component_idx[4].resize(2);
                        
                        // Variable u.
                        derivative_var_data[4][0] = data_velocity;
                        derivative_var_component_idx[4][0] = 0;
                        
                        // Variable w.
                        derivative_var_data[4][1] = data_velocity;
                        derivative_var_component_idx[4][1] = 2;
                        
                        break;
                    }
                    case DIRECTION::Y_DIRECTION:
                    {
                        /*
                         * Mass equation.
                         */
                        
                        derivative_var_data[0].resize(0);
                        derivative_var_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[1].resize(0);
                        derivative_var_component_idx[1].resize(0);
                        
                        derivative_var_data[2].resize(1);
                        derivative_var_component_idx[2].resize(1);
                        
                        // Variable w.
                        derivative_var_data[2][0] = data_velocity;
                        derivative_var_component_idx[2][0] = 2;
                        
                        derivative_var_data[3].resize(1);
                        derivative_var_component_idx[3].resize(1);
                        
                        // Variable v.
                        derivative_var_data[3][0] = data_velocity;
                        derivative_var_component_idx[3][0] = 1;
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[4].resize(2);
                        derivative_var_component_idx[4].resize(2);
                        
                        // Variable v.
                        derivative_var_data[4][0] = data_velocity;
                        derivative_var_component_idx[4][0] = 1;
                        
                        // Variable w.
                        derivative_var_data[4][1] = data_velocity;
                        derivative_var_component_idx[4][1] = 2;
                        
                        break;
                    }
                    case DIRECTION::Z_DIRECTION:
                    {
                        /*
                         * Mass equation.
                         */
                        
                        derivative_var_data[0].resize(0);
                        derivative_var_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[1].resize(1);
                        derivative_var_component_idx[1].resize(1);
                        
                        // Variable u.
                        derivative_var_data[1][0] = data_velocity;
                        derivative_var_component_idx[1][0] = 0;
                        
                        derivative_var_data[2].resize(1);
                        derivative_var_component_idx[2].resize(1);
                        
                        // Variable v.
                        derivative_var_data[2][0] = data_velocity;
                        derivative_var_component_idx[2][0] = 1;
                        
                        derivative_var_data[3].resize(1);
                        derivative_var_component_idx[3].resize(1);
                        
                        // Variable w.
                        derivative_var_data[3][0] = data_velocity;
                        derivative_var_component_idx[3][0] = 2;
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[4].resize(4);
                        derivative_var_component_idx[4].resize(4);
                        
                        // Variable u.
                        derivative_var_data[4][0] = data_velocity;
                        derivative_var_component_idx[4][0] = 0;
                        
                        // Variable v.
                        derivative_var_data[4][1] = data_velocity;
                        derivative_var_component_idx[4][1] = 1;
                        
                        // Variable w.
                        derivative_var_data[4][2] = data_velocity;
                        derivative_var_component_idx[4][2] = 2;
                        
                        // Variable T.
                        derivative_var_data[4][3] = data_temperature;
                        derivative_var_component_idx[4][3] = 0;
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::"
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
                    << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::"
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
FlowModelDiffusiveFluxUtilitiesSingleSpecies::getCellDataOfDiffusiveFluxDiffusivities(
    std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > >& diffusivities_data,
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
            << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::"
            << "getCellDataOfDiffusiveFluxDiffusivities()\n"
            << "No patch is registered yet."
            << std::endl);
    }
    
    if (!d_cell_data_computed_diffusivities)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::getCellDataOfDiffusiveFluxDiffusivities()\n"
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
                         * Mass equation.
                         */
                        
                        diffusivities_data[0].resize(0);
                        diffusivities_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[1].resize(1);
                        diffusivities_component_idx[1].resize(1);
                        
                        // -(4/3*mu + mu_v).
                        diffusivities_data[1][0] = d_data_diffusivities;
                        diffusivities_component_idx[1][0] = 0;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[2].resize(2);
                        diffusivities_component_idx[2].resize(2);
                        
                        // -u*(4/3*mu + mu_v).
                        diffusivities_data[2][0] = d_data_diffusivities;
                        diffusivities_component_idx[2][0] = 1;
                        
                        // -kappa.
                        diffusivities_data[2][1] = d_data_diffusivities;
                        diffusivities_component_idx[2][1] = 2;
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::"
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
                    << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::getCellDataOfDiffusiveFluxDiffusivities()\n"
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
                         * Mass equation.
                         */
                        
                        diffusivities_data[0].resize(0);
                        diffusivities_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[1].resize(1);
                        diffusivities_component_idx[1].resize(1);
                        
                        // -(4/3*mu + mu_v).
                        diffusivities_data[1][0] = d_data_diffusivities;
                        diffusivities_component_idx[1][0] = 0;
                        
                        diffusivities_data[2].resize(1);
                        diffusivities_component_idx[2].resize(1);
                        
                        // -mu.
                        diffusivities_data[2][0] = d_data_diffusivities;
                        diffusivities_component_idx[2][0] = 2;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[3].resize(3);
                        diffusivities_component_idx[3].resize(3);
                        
                        // -u*(4/3*mu + mu_v).
                        diffusivities_data[3][0] = d_data_diffusivities;
                        diffusivities_component_idx[3][0] = 3;
                        
                        // -v*mu.
                        diffusivities_data[3][1] = d_data_diffusivities;
                        diffusivities_component_idx[3][1] = 8;
                        
                        // -kappa.
                        diffusivities_data[3][2] = d_data_diffusivities;
                        diffusivities_component_idx[3][2] = 9;
                        
                        break;
                    }
                    case DIRECTION::Y_DIRECTION:
                    {
                        /*
                         * Mass equation.
                         */
                        
                        diffusivities_data[0].resize(0);
                        diffusivities_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[1].resize(1);
                        diffusivities_component_idx[1].resize(1);
                        
                        // 2/3*mu - mu_v.
                        diffusivities_data[1][0] = d_data_diffusivities;
                        diffusivities_component_idx[1][0] = 1;
                        
                        diffusivities_data[2].resize(1);
                        diffusivities_component_idx[2].resize(1);
                        
                        // -mu.
                        diffusivities_data[2][0] = d_data_diffusivities;
                        diffusivities_component_idx[2][0] = 2;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[3].resize(2);
                        diffusivities_component_idx[3].resize(2);
                        
                        // -v*mu.
                        diffusivities_data[3][0] = d_data_diffusivities;
                        diffusivities_component_idx[3][0] = 8;
                        
                        // u*(2/3*mu - mu_v).
                        diffusivities_data[3][1] = d_data_diffusivities;
                        diffusivities_component_idx[3][1] = 5;
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::"
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
                         * Mass equation.
                         */
                        
                        diffusivities_data[0].resize(0);
                        diffusivities_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[1].resize(1);
                        diffusivities_component_idx[1].resize(1);
                        
                        // -mu.
                        diffusivities_data[1][0] = d_data_diffusivities;
                        diffusivities_component_idx[1][0] = 2;
                        
                        diffusivities_data[2].resize(1);
                        diffusivities_component_idx[2].resize(1);
                        
                        // 2/3*mu - mu_v.
                        diffusivities_data[2][0] = d_data_diffusivities;
                        diffusivities_component_idx[2][0] = 1;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[3].resize(2);
                        diffusivities_component_idx[3].resize(2);
                        
                        // v*(2/3*mu - mu_v).
                        diffusivities_data[3][0] = d_data_diffusivities;
                        diffusivities_component_idx[3][0] = 6;
                        
                        // -u*mu.
                        diffusivities_data[3][1] = d_data_diffusivities;
                        diffusivities_component_idx[3][1] = 7;
                        
                        break;
                    }
                    case DIRECTION::Y_DIRECTION:
                    {
                        /*
                         * Mass equation.
                         */
                        
                        diffusivities_data[0].resize(0);
                        diffusivities_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[1].resize(1);
                        diffusivities_component_idx[1].resize(1);
                        
                        // -mu.
                        diffusivities_data[1][0] = d_data_diffusivities;
                        diffusivities_component_idx[1][0] = 2;
                        
                        diffusivities_data[2].resize(1);
                        diffusivities_component_idx[2].resize(1);
                        
                        // -(4/3*mu + mu_v).
                        diffusivities_data[2][0] = d_data_diffusivities;
                        diffusivities_component_idx[2][0] = 0;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[3].resize(3);
                        diffusivities_component_idx[3].resize(3);
                        
                        // -u*mu.
                        diffusivities_data[3][0] = d_data_diffusivities;
                        diffusivities_component_idx[3][0] = 7;
                        
                        // -v*(4/3*mu + mu_v).
                        diffusivities_data[3][1] = d_data_diffusivities;
                        diffusivities_component_idx[3][1] = 4;
                        
                        // -kappa.
                        diffusivities_data[3][2] = d_data_diffusivities;
                        diffusivities_component_idx[3][2] = 9;
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::"
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
                    << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::getCellDataOfDiffusiveFluxDiffusivities()\n"
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
                         * Mass equation.
                         */
                        
                        diffusivities_data[0].resize(0);
                        diffusivities_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[1].resize(1);
                        diffusivities_component_idx[1].resize(1);
                        
                        // -(4/3*mu + mu_v).
                        diffusivities_data[1][0] = d_data_diffusivities;
                        diffusivities_component_idx[1][0] = 0;
                        
                        diffusivities_data[2].resize(1);
                        diffusivities_component_idx[2].resize(1);
                        
                        // -mu.
                        diffusivities_data[2][0] = d_data_diffusivities;
                        diffusivities_component_idx[2][0] = 2;
                        
                        diffusivities_data[3].resize(1);
                        diffusivities_component_idx[3].resize(1);
                        
                        // -mu.
                        diffusivities_data[3][0] = d_data_diffusivities;
                        diffusivities_component_idx[3][0] = 2;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[4].resize(4);
                        diffusivities_component_idx[4].resize(4);
                        
                        // -u*(4/3*mu + mu_v).
                        diffusivities_data[4][0] = d_data_diffusivities;
                        diffusivities_component_idx[4][0] = 3;
                        
                        // -v*mu.
                        diffusivities_data[4][1] = d_data_diffusivities;
                        diffusivities_component_idx[4][1] = 10;
                        
                        // -w*mu.
                        diffusivities_data[4][2] = d_data_diffusivities;
                        diffusivities_component_idx[4][2] = 11;
                        
                        // -kappa.
                        diffusivities_data[4][3] = d_data_diffusivities;
                        diffusivities_component_idx[4][3] = 12;
                        
                        break;
                    }
                    case DIRECTION::Y_DIRECTION:
                    {
                        /*
                         * Mass equation.
                         */
                        
                        diffusivities_data[0].resize(0);
                        diffusivities_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[1].resize(1);
                        diffusivities_component_idx[1].resize(1);
                        
                        // 2/3*mu - mu_v.
                        diffusivities_data[1][0] = d_data_diffusivities;
                        diffusivities_component_idx[1][0] = 1;
                        
                        diffusivities_data[2].resize(1);
                        diffusivities_component_idx[2].resize(1);
                        
                        // -mu.
                        diffusivities_data[2][0] = d_data_diffusivities;
                        diffusivities_component_idx[2][0] = 2;
                        
                        diffusivities_data[3].resize(0);
                        diffusivities_component_idx[3].resize(0);
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[4].resize(2);
                        diffusivities_component_idx[4].resize(2);
                        
                        // -v*mu.
                        diffusivities_data[4][0] = d_data_diffusivities;
                        diffusivities_component_idx[4][0] = 10;
                        
                        // u*(2/3*mu - mu_v).
                        diffusivities_data[4][1] = d_data_diffusivities;
                        diffusivities_component_idx[4][1] = 6;
                        
                        break;
                    }
                    case DIRECTION::Z_DIRECTION:
                    {
                        /*
                         * Mass equation.
                         */
                        
                        diffusivities_data[0].resize(0);
                        diffusivities_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[1].resize(1);
                        diffusivities_component_idx[1].resize(1);
                        
                        // 2/3*mu - mu_v.
                        diffusivities_data[1][0] = d_data_diffusivities;
                        diffusivities_component_idx[1][0] = 1;
                        
                        diffusivities_data[2].resize(0);
                        diffusivities_component_idx[2].resize(0);
                        
                        diffusivities_data[3].resize(1);
                        diffusivities_component_idx[3].resize(1);
                        
                        // -mu.
                        diffusivities_data[3][0] = d_data_diffusivities;
                        diffusivities_component_idx[3][0] = 2;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[4].resize(2);
                        diffusivities_component_idx[4].resize(2);
                        
                        // -w*mu.
                        diffusivities_data[4][0] = d_data_diffusivities;
                        diffusivities_component_idx[4][0] = 11;
                        
                        // u*(2/3*mu - mu_v).
                        diffusivities_data[4][1] = d_data_diffusivities;
                        diffusivities_component_idx[4][1] = 6;
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::"
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
                         * Mass equation.
                         */
                        
                        diffusivities_data[0].resize(0);
                        diffusivities_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[1].resize(1);
                        diffusivities_component_idx[1].resize(1);
                        
                        // -mu.
                        diffusivities_data[1][0] = d_data_diffusivities;
                        diffusivities_component_idx[1][0] = 2;
                        
                        diffusivities_data[2].resize(1);
                        diffusivities_component_idx[2].resize(1);
                        
                        // 2/3*mu - mu_v.
                        diffusivities_data[2][0] = d_data_diffusivities;
                        diffusivities_component_idx[2][0] = 1;
                        
                        diffusivities_data[3].resize(0);
                        diffusivities_component_idx[3].resize(0);
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[4].resize(2);
                        diffusivities_component_idx[4].resize(2);
                        
                        // v*(2/3*mu - mu_v).
                        diffusivities_data[4][0] = d_data_diffusivities;
                        diffusivities_component_idx[4][0] = 7;
                        
                        // -u*mu.
                        diffusivities_data[4][1] = d_data_diffusivities;
                        diffusivities_component_idx[4][1] = 9;
                        
                        break;
                    }
                    case DIRECTION::Y_DIRECTION:
                    {
                        /*
                         * Mass equation.
                         */
                        
                        diffusivities_data[0].resize(0);
                        diffusivities_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[1].resize(1);
                        diffusivities_component_idx[1].resize(1);
                        
                        // -mu.
                        diffusivities_data[1][0] = d_data_diffusivities;
                        diffusivities_component_idx[1][0] = 2;
                        
                        diffusivities_data[2].resize(1);
                        diffusivities_component_idx[2].resize(1);
                        
                        // -(4/3*mu + mu_v).
                        diffusivities_data[2][0] = d_data_diffusivities;
                        diffusivities_component_idx[2][0] = 0;
                        
                        diffusivities_data[3].resize(1);
                        diffusivities_component_idx[3].resize(1);
                        
                        // -mu.
                        diffusivities_data[3][0] = d_data_diffusivities;
                        diffusivities_component_idx[3][0] = 2;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[4].resize(4);
                        diffusivities_component_idx[4].resize(4);
                        
                        // -u*mu.
                        diffusivities_data[4][0] = d_data_diffusivities;
                        diffusivities_component_idx[4][0] = 9;
                        
                        // -v*(4/3*mu + mu_v).
                        diffusivities_data[4][1] = d_data_diffusivities;
                        diffusivities_component_idx[4][1] = 4;
                        
                        // -w*mu.
                        diffusivities_data[4][2] = d_data_diffusivities;
                        diffusivities_component_idx[4][2] = 11;
                        
                        // -kappa.
                        diffusivities_data[4][3] = d_data_diffusivities;
                        diffusivities_component_idx[4][3] = 12;
                        
                        break;
                    }
                    case DIRECTION::Z_DIRECTION:
                    {
                        /*
                         * Mass equation.
                         */
                        
                        diffusivities_data[0].resize(0);
                        diffusivities_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[1].resize(0);
                        diffusivities_component_idx[1].resize(0);
                        
                        diffusivities_data[2].resize(1);
                        diffusivities_component_idx[2].resize(1);
                        
                        // 2/3*mu - mu_v.
                        diffusivities_data[2][0] = d_data_diffusivities;
                        diffusivities_component_idx[2][0] = 1;
                        
                        diffusivities_data[3].resize(1);
                        diffusivities_component_idx[3].resize(1);
                        
                        // -mu.
                        diffusivities_data[3][0] = d_data_diffusivities;
                        diffusivities_component_idx[3][0] = 2;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[4].resize(2);
                        diffusivities_component_idx[4].resize(2);
                        
                        // -w*mu.
                        diffusivities_data[4][0] = d_data_diffusivities;
                        diffusivities_component_idx[4][0] = 11;
                        
                        // v*(2/3*mu - mu_v).
                        diffusivities_data[4][1] = d_data_diffusivities;
                        diffusivities_component_idx[4][1] = 7;
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::"
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
                         * Mass equation.
                         */
                        
                        diffusivities_data[0].resize(0);
                        diffusivities_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[1].resize(1);
                        diffusivities_component_idx[1].resize(1);
                        
                        // -mu.
                        diffusivities_data[1][0] = d_data_diffusivities;
                        diffusivities_component_idx[1][0] = 2;
                        
                        diffusivities_data[2].resize(0);
                        diffusivities_component_idx[2].resize(0);
                        
                        diffusivities_data[3].resize(1);
                        diffusivities_component_idx[3].resize(1);
                        
                        // 2/3*mu - mu_v.
                        diffusivities_data[3][0] = d_data_diffusivities;
                        diffusivities_component_idx[3][0] = 1;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[4].resize(2);
                        diffusivities_component_idx[4].resize(2);
                        
                        // w*(2/3*mu - mu_v).
                        diffusivities_data[4][0] = d_data_diffusivities;
                        diffusivities_component_idx[4][0] = 8;
                        
                        // -u*mu.
                        diffusivities_data[4][1] = d_data_diffusivities;
                        diffusivities_component_idx[4][1] = 9;
                        
                        break;
                    }
                    case DIRECTION::Y_DIRECTION:
                    {
                        /*
                         * Mass equation.
                         */
                        
                        diffusivities_data[0].resize(0);
                        diffusivities_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[1].resize(0);
                        diffusivities_component_idx[1].resize(0);
                        
                        diffusivities_data[2].resize(1);
                        diffusivities_component_idx[2].resize(1);
                        
                        // -mu.
                        diffusivities_data[2][0] = d_data_diffusivities;
                        diffusivities_component_idx[2][0] = 2;
                        
                        diffusivities_data[3].resize(1);
                        diffusivities_component_idx[3].resize(1);
                        
                        // 2/3*mu - mu_v.
                        diffusivities_data[3][0] = d_data_diffusivities;
                        diffusivities_component_idx[3][0] = 1;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[4].resize(2);
                        diffusivities_component_idx[4].resize(2);
                        
                        // w*(2/3*mu - mu_v).
                        diffusivities_data[4][0] = d_data_diffusivities;
                        diffusivities_component_idx[4][0] = 8;
                        
                        // -v*mu.
                        diffusivities_data[4][1] = d_data_diffusivities;
                        diffusivities_component_idx[4][1] = 10;
                        
                        break;
                    }
                    case DIRECTION::Z_DIRECTION:
                    {
                        /*
                         * Mass equation.
                         */
                        
                        diffusivities_data[0].resize(0);
                        diffusivities_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[1].resize(1);
                        diffusivities_component_idx[1].resize(1);
                        
                        // -mu.
                        diffusivities_data[1][0] = d_data_diffusivities;
                        diffusivities_component_idx[1][0] = 2;
                        
                        diffusivities_data[2].resize(1);
                        diffusivities_component_idx[2].resize(1);
                        
                        // -mu.
                        diffusivities_data[2][0] = d_data_diffusivities;
                        diffusivities_component_idx[2][0] = 2;
                        
                        diffusivities_data[3].resize(1);
                        diffusivities_component_idx[3].resize(1);
                        
                        // -(4/3*mu + mu_v).
                        diffusivities_data[3][0] = d_data_diffusivities;
                        diffusivities_component_idx[3][0] = 0;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[4].resize(4);
                        diffusivities_component_idx[4].resize(4);
                        
                        // -u*mu.
                        diffusivities_data[4][0] = d_data_diffusivities;
                        diffusivities_component_idx[4][0] = 9;
                        
                        // -v*mu.
                        diffusivities_data[4][1] = d_data_diffusivities;
                        diffusivities_component_idx[4][1] = 10;
                        
                        // -w*(4/3*mu + mu_v).
                        diffusivities_data[4][2] = d_data_diffusivities;
                        diffusivities_component_idx[4][2] = 5;
                        
                        // -kappa.
                        diffusivities_data[4][3] = d_data_diffusivities;
                        diffusivities_component_idx[4][3] = 12;
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::"
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
                    << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::getCellDataOfDiffusiveFluxDiffusivities()\n"
                    << "There are only x-direction, y-direction and z-direction for three-dimensional problem."
                    << std::endl);
            }
        }
    }
}


/*
 * Get the cell data that needs interpolation to sides for computing side data of diffusivities in the
 * diffusive flux.
 */
void
FlowModelDiffusiveFluxUtilitiesSingleSpecies::getCellDataForInterpolationToSideDataForDiffusiveFluxDiffusivities(
    std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& var_data_for_diffusivities,
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
            << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::"
            << "getCellDataForInterpolationToSideDataForDiffusiveFluxDiffusivities()\n"
            << "No patch is registered yet."
            << std::endl);
    }
    
    if (!d_cell_data_computed_shear_viscosity)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::getCellDataForInterpolationToSideDataForDiffusiveFluxDiffusivities()\n"
            << "Cell data of 'SHEAR_VISCOSITY' is not registered/computed yet."
            << std::endl);
    }
    
    if (!d_cell_data_computed_bulk_viscosity)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::getCellDataForInterpolationToSideDataForDiffusiveFluxDiffusivities()\n"
            << "Cell data of 'BULK_VISCOSITY' is not registered/computed yet."
            << std::endl);
    }
    
    if (!d_cell_data_computed_thermal_conductivity)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::getCellDataForInterpolationToSideDataForDiffusiveFluxDiffusivities()\n"
            << "Cell data of 'THERMAL_CONDUCTIVITY' is not registered/computed yet."
            << std::endl);
    }
    
    var_data_for_diffusivities.resize(3 + d_dim.getValue());
    var_data_for_diffusivities_component_idx.resize(3 + d_dim.getValue());
    
    var_data_for_diffusivities[0] = d_data_shear_viscosity;
    var_data_for_diffusivities_component_idx[0] = 0;
    var_data_for_diffusivities[1] = d_data_bulk_viscosity;
    var_data_for_diffusivities_component_idx[1] = 0;
    var_data_for_diffusivities[2] = d_data_thermal_conductivity;
    var_data_for_diffusivities_component_idx[2] = 0;
    
    // Get the cell data of velocity.
    HAMERS_SHARED_PTR<pdat::CellData<Real> > data_velocity =
        flow_model_tmp->getCellData("VELOCITY");
    
    if (d_dim == tbox::Dimension(1))
    {
        var_data_for_diffusivities[3] = data_velocity;
        var_data_for_diffusivities_component_idx[3] = 0;
    }
    else if (d_dim == tbox::Dimension(2))
    {
        var_data_for_diffusivities[3] = data_velocity;
        var_data_for_diffusivities_component_idx[3] = 0;
        var_data_for_diffusivities[4] = data_velocity;
        var_data_for_diffusivities_component_idx[4] = 1;
    }
    else if (d_dim == tbox::Dimension(3))
    {
        var_data_for_diffusivities[3] = data_velocity;
        var_data_for_diffusivities_component_idx[3] = 0;
        var_data_for_diffusivities[4] = data_velocity;
        var_data_for_diffusivities_component_idx[4] = 1;
        var_data_for_diffusivities[5] = data_velocity;
        var_data_for_diffusivities_component_idx[5] = 2;
    }
    
    if (d_use_subgrid_scale_model)
    {
        std::vector<std::string> var_to_interpolate;
        std::vector<int> var_to_interpolate_component_idx;
        d_flow_model_subgrid_scale_model->getDerivedVariablesForInterpolationToSideData(var_to_interpolate, var_to_interpolate_component_idx);
        
        for (int vi = 0; vi < static_cast<int>(var_to_interpolate.size()); vi++)
        {
            const std::string var_str = var_to_interpolate[vi];
            HAMERS_SHARED_PTR<pdat::CellData<Real> > data_u = flow_model_tmp->getCellData(var_str);
            var_data_for_diffusivities.push_back(data_u);
            var_data_for_diffusivities_component_idx.push_back(var_to_interpolate_component_idx[vi]);
        }
    }
}


/*
 * Compute the side data of the diffusivities in the diffusive flux with the interpolated side data.
 */
void
FlowModelDiffusiveFluxUtilitiesSingleSpecies::computeSideDataOfDiffusiveFluxDiffusivities(
    const std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& var_data_for_diffusivities)
{
    if (d_num_subghosts_diffusivities > -hier::IntVector::getOne(d_dim))
    {
        if (!d_side_data_diffusivities_computed)
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
            
            /*
             * Get the dimension of the interior box.
             */
            
            const hier::Box interior_box = patch.getBox();
            const hier::IntVector interior_dims = interior_box.numberCells();
            
            const hier::IntVector num_ghosts = var_data_for_diffusivities[0]->getGhostCellWidth();
            const hier::IntVector ghostcell_dims = var_data_for_diffusivities[0]->getGhostBox().numberCells();
            
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
            std::vector<std::string> var_to_interpolate;
            std::vector<int> var_to_interpolate_component_idx;
            if (d_use_subgrid_scale_model)
            {
                d_flow_model_subgrid_scale_model->getDerivedVariablesForInterpolationToSideData(var_to_interpolate, var_to_interpolate_component_idx);
            }
            
            TBOX_ASSERT(static_cast<int>(var_data_for_diffusivities.size()) == 3 + d_dim.getValue() + static_cast<int>(var_to_interpolate.size()));
            TBOX_ASSERT(num_ghosts <= d_num_subghosts_diffusivities);
            
            for (int vi = 0; vi < static_cast<int>(var_data_for_diffusivities.size()); vi++)
            {
                TBOX_ASSERT(var_data_for_diffusivities[vi]->getGhostCellWidth() == num_ghosts);
                TBOX_ASSERT(var_data_for_diffusivities[vi]->getGhostBox().contains(interior_box));
            }
#endif
            
            if (d_dim == tbox::Dimension(1))
            {
                Real* mu_x    = var_data_for_diffusivities[0]->getPointer(0, 0);
                Real* mu_v_x  = var_data_for_diffusivities[1]->getPointer(0, 0);
                Real* kappa_x = var_data_for_diffusivities[2]->getPointer(0, 0);
                Real* u_x     = var_data_for_diffusivities[3]->getPointer(0, 0);
                
                /*
                 * Compute the diffusivities for the x-direction derivatives.
                 */
                
                Real* D_00 = d_side_data_diffusivities->getPointer(0, 0);
                Real* D_01 = d_side_data_diffusivities->getPointer(0, 1);
                Real* D_02 = d_side_data_diffusivities->getPointer(0, 2);
                
                for (int i = -num_ghosts[0]; i < interior_dims[0] + num_ghosts[0] + 1; i++)
                {
                    const int idx_diffusivities = i + d_num_subghosts_diffusivities[0];
                    const int idx_var_data      = i + num_ghosts[0];
                    
                    D_00[idx_diffusivities] = -(Real(4)/Real(3)*mu_x[idx_var_data] + mu_v_x[idx_var_data]);
                    D_01[idx_diffusivities] = -u_x[idx_var_data]*(Real(4)/Real(3)*mu_x[idx_var_data] +
                        mu_v_x[idx_var_data]);
                    D_02[idx_diffusivities] = -kappa_x[idx_var_data];
                }
            }
            else if (d_dim == tbox::Dimension(2))
            {
                Real* mu_x    = var_data_for_diffusivities[0]->getPointer(0, 0);
                Real* mu_v_x  = var_data_for_diffusivities[1]->getPointer(0, 0);
                Real* kappa_x = var_data_for_diffusivities[2]->getPointer(0, 0);
                
                Real* mu_y    = var_data_for_diffusivities[0]->getPointer(1, 0);
                Real* mu_v_y  = var_data_for_diffusivities[1]->getPointer(1, 0);
                Real* kappa_y = var_data_for_diffusivities[2]->getPointer(1, 0);
                
                Real* u_x = var_data_for_diffusivities[3]->getPointer(0, 0);
                Real* v_x = var_data_for_diffusivities[4]->getPointer(0, 0);
                
                Real* u_y = var_data_for_diffusivities[3]->getPointer(1, 0);
                Real* v_y = var_data_for_diffusivities[4]->getPointer(1, 0);
                
                /*
                 * Compute the diffusivities for the x-direction derivatives.
                 */
                
                Real* D_00 = d_side_data_diffusivities->getPointer(0, 0);
                Real* D_01 = d_side_data_diffusivities->getPointer(0, 1);
                Real* D_02 = d_side_data_diffusivities->getPointer(0, 2);
                Real* D_03 = d_side_data_diffusivities->getPointer(0, 3);
                Real* D_04 = d_side_data_diffusivities->getPointer(0, 4);
                Real* D_05 = d_side_data_diffusivities->getPointer(0, 5);
                Real* D_06 = d_side_data_diffusivities->getPointer(0, 6);
                
                for (int j = -num_ghosts[1]; j < interior_dims[1] + num_ghosts[1]; j++)
                {
                    for (int i = -num_ghosts[0]; i < interior_dims[0] + num_ghosts[0] + 1; i++)
                    {
                        const int idx_diffusivities = (i + d_num_subghosts_diffusivities[0]) +
                            (j + d_num_subghosts_diffusivities[1])*(d_subghostcell_dims_diffusivities[0] + 1);
                        
                        const int idx_var_data = (i + num_ghosts[0]) +
                            (j + num_ghosts[1])*(ghostcell_dims[0] + 1);
                        
                        D_00[idx_diffusivities] = -(Real(4)/Real(3)*mu_x[idx_var_data] + mu_v_x[idx_var_data]);
                        D_01[idx_diffusivities] = Real(2)/Real(3)*mu_x[idx_var_data] - mu_v_x[idx_var_data];
                        D_02[idx_diffusivities] = -mu_x[idx_var_data];
                        D_03[idx_diffusivities] = -u_x[idx_var_data]*(Real(4)/Real(3)*mu_x[idx_var_data] +
                            mu_v_x[idx_var_data]);
                        D_04[idx_diffusivities] = u_x[idx_var_data]*(Real(2)/Real(3)*mu_x[idx_var_data] -
                            mu_v_x[idx_var_data]);
                        D_05[idx_diffusivities] = -v_x[idx_var_data]*mu_x[idx_var_data];
                        D_06[idx_diffusivities] = -kappa_x[idx_var_data];
                    }
                }
                
                /*
                 * Compute the diffusivities for the y-direction derivatives.
                 */
                
                D_00 = d_side_data_diffusivities->getPointer(1, 0);
                D_01 = d_side_data_diffusivities->getPointer(1, 1);
                D_02 = d_side_data_diffusivities->getPointer(1, 2);
                D_03 = d_side_data_diffusivities->getPointer(1, 3);
                D_04 = d_side_data_diffusivities->getPointer(1, 4);
                D_05 = d_side_data_diffusivities->getPointer(1, 5);
                D_06 = d_side_data_diffusivities->getPointer(1, 6);
                
                for (int j = -num_ghosts[1]; j < interior_dims[1] + num_ghosts[1] + 1; j++)
                {
                    for (int i = -num_ghosts[0]; i < interior_dims[0] + num_ghosts[0]; i++)
                    {
                        const int idx_diffusivities = (i + d_num_subghosts_diffusivities[0]) +
                            (j + d_num_subghosts_diffusivities[1])*d_subghostcell_dims_diffusivities[0];
                        
                        const int idx_var_data = (i + num_ghosts[0]) +
                            (j + num_ghosts[1])*ghostcell_dims[0];
                        
                        D_00[idx_diffusivities] = -(Real(4)/Real(3)*mu_y[idx_var_data] + mu_v_y[idx_var_data]);
                        D_01[idx_diffusivities] = Real(2)/Real(3)*mu_y[idx_var_data] - mu_v_y[idx_var_data];
                        D_02[idx_diffusivities] = -mu_y[idx_var_data];
                        D_03[idx_diffusivities] = -v_y[idx_var_data]*(Real(4)/Real(3)*mu_y[idx_var_data] +
                            mu_v_y[idx_var_data]);
                        D_04[idx_diffusivities] = v_y[idx_var_data]*(Real(2)/Real(3)*mu_y[idx_var_data] -
                            mu_v_y[idx_var_data]);
                        D_05[idx_diffusivities] = -u_y[idx_var_data]*mu_y[idx_var_data];
                        D_06[idx_diffusivities] = -kappa_y[idx_var_data];
                    }
                }
            }
            else if (d_dim == tbox::Dimension(3))
            {
                Real* mu_x    = var_data_for_diffusivities[0]->getPointer(0, 0);
                Real* mu_v_x  = var_data_for_diffusivities[1]->getPointer(0, 0);
                Real* kappa_x = var_data_for_diffusivities[2]->getPointer(0, 0);
                
                Real* mu_y    = var_data_for_diffusivities[0]->getPointer(1, 0);
                Real* mu_v_y  = var_data_for_diffusivities[1]->getPointer(1, 0);
                Real* kappa_y = var_data_for_diffusivities[2]->getPointer(1, 0);
                
                Real* mu_z    = var_data_for_diffusivities[0]->getPointer(2, 0);
                Real* mu_v_z  = var_data_for_diffusivities[1]->getPointer(2, 0);
                Real* kappa_z = var_data_for_diffusivities[2]->getPointer(2, 0);
                
                Real* u_x = var_data_for_diffusivities[3]->getPointer(0, 0);
                Real* v_x = var_data_for_diffusivities[4]->getPointer(0, 0);
                Real* w_x = var_data_for_diffusivities[5]->getPointer(0, 0);
                
                Real* u_y = var_data_for_diffusivities[3]->getPointer(1, 0);
                Real* v_y = var_data_for_diffusivities[4]->getPointer(1, 0);
                Real* w_y = var_data_for_diffusivities[5]->getPointer(1, 0);
                
                Real* u_z = var_data_for_diffusivities[3]->getPointer(2, 0);
                Real* v_z = var_data_for_diffusivities[4]->getPointer(2, 0);
                Real* w_z = var_data_for_diffusivities[5]->getPointer(2, 0);
                
                /*
                 * Compute the diffusivities for the x-direction derivatives.
                 */
                
                Real* D_00 = d_side_data_diffusivities->getPointer(0, 0);
                Real* D_01 = d_side_data_diffusivities->getPointer(0, 1);
                Real* D_02 = d_side_data_diffusivities->getPointer(0, 2);
                Real* D_03 = d_side_data_diffusivities->getPointer(0, 3);
                Real* D_04 = d_side_data_diffusivities->getPointer(0, 4);
                Real* D_05 = d_side_data_diffusivities->getPointer(0, 5);
                Real* D_06 = d_side_data_diffusivities->getPointer(0, 6);
                Real* D_07 = d_side_data_diffusivities->getPointer(0, 7);
                
                for (int k = -num_ghosts[2]; k < interior_dims[2] + num_ghosts[2]; k++)
                {
                    for (int j = -num_ghosts[1]; j < interior_dims[1] + num_ghosts[1]; j++)
                    {
                        for (int i = -num_ghosts[0]; i < interior_dims[0] + num_ghosts[0] + 1; i++)
                        {
                            const int idx_diffusivities = (i + d_num_subghosts_diffusivities[0]) +
                                (j + d_num_subghosts_diffusivities[1])*(d_subghostcell_dims_diffusivities[0] + 1) +
                                (k + d_num_subghosts_diffusivities[2])*(d_subghostcell_dims_diffusivities[0] + 1)*
                                    d_subghostcell_dims_diffusivities[1];
                            
                            const int idx_var_data = (i + num_ghosts[0]) +
                                (j + num_ghosts[1])*(ghostcell_dims[0] + 1) +
                                (k + num_ghosts[2])*(ghostcell_dims[0] + 1)*
                                    ghostcell_dims[1];
                            
                            D_00[idx_diffusivities] = -(Real(4)/Real(3)*mu_x[idx_var_data] + mu_v_x[idx_var_data]);
                            D_01[idx_diffusivities] = Real(2)/Real(3)*mu_x[idx_var_data] - mu_v_x[idx_var_data];
                            D_02[idx_diffusivities] = -mu_x[idx_var_data];
                            D_03[idx_diffusivities] = -u_x[idx_var_data]*(Real(4)/Real(3)*mu_x[idx_var_data] +
                                mu_v_x[idx_var_data]);
                            D_04[idx_diffusivities] = u_x[idx_var_data]*(Real(2)/Real(3)*mu_x[idx_var_data] -
                                mu_v_x[idx_var_data]);
                            D_05[idx_diffusivities] = -v_x[idx_var_data]*mu_x[idx_var_data];
                            D_06[idx_diffusivities] = -w_x[idx_var_data]*mu_x[idx_var_data];
                            D_07[idx_diffusivities] = -kappa_x[idx_var_data];
                        }
                    }
                }
                
                /*
                 * Compute the diffusivities for the y-direction derivatives.
                 */
                
                D_00 = d_side_data_diffusivities->getPointer(1, 0);
                D_01 = d_side_data_diffusivities->getPointer(1, 1);
                D_02 = d_side_data_diffusivities->getPointer(1, 2);
                D_03 = d_side_data_diffusivities->getPointer(1, 3);
                D_04 = d_side_data_diffusivities->getPointer(1, 4);
                D_05 = d_side_data_diffusivities->getPointer(1, 5);
                D_06 = d_side_data_diffusivities->getPointer(1, 6);
                D_07 = d_side_data_diffusivities->getPointer(1, 7);
                
                for (int k = -num_ghosts[2]; k < interior_dims[2] + num_ghosts[2]; k++)
                {
                    for (int j = -num_ghosts[1]; j < interior_dims[1] + num_ghosts[1] + 1; j++)
                    {
                        for (int i = -num_ghosts[0]; i < interior_dims[0] + num_ghosts[0]; i++)
                        {
                            const int idx_diffusivities = (i + d_num_subghosts_diffusivities[0]) +
                                (j + d_num_subghosts_diffusivities[1])*d_subghostcell_dims_diffusivities[0] +
                                (k + d_num_subghosts_diffusivities[2])*d_subghostcell_dims_diffusivities[0]*
                                    (d_subghostcell_dims_diffusivities[1] + 1);
                            
                            const int idx_var_data = (i + num_ghosts[0]) +
                                (j + num_ghosts[1])*ghostcell_dims[0] +
                                (k + num_ghosts[2])*ghostcell_dims[0]*
                                    (ghostcell_dims[1] + 1);
                            
                            D_00[idx_diffusivities] = -(Real(4)/Real(3)*mu_y[idx_var_data] + mu_v_y[idx_var_data]);
                            D_01[idx_diffusivities] = Real(2)/Real(3)*mu_y[idx_var_data] - mu_v_y[idx_var_data];
                            D_02[idx_diffusivities] = -mu_y[idx_var_data];
                            D_03[idx_diffusivities] = -v_y[idx_var_data]*(Real(4)/Real(3)*mu_y[idx_var_data] +
                                mu_v_y[idx_var_data]);
                            D_04[idx_diffusivities] = v_y[idx_var_data]*(Real(2)/Real(3)*mu_y[idx_var_data] -
                                mu_v_y[idx_var_data]);
                            D_05[idx_diffusivities] = -u_y[idx_var_data]*mu_y[idx_var_data];
                            D_06[idx_diffusivities] = -w_y[idx_var_data]*mu_y[idx_var_data];
                            D_07[idx_diffusivities] = -kappa_y[idx_var_data];
                        }
                    }
                }
                
                /*
                 * Compute the diffusivities for the z-direction derivatives.
                 */
                
                D_00 = d_side_data_diffusivities->getPointer(2, 0);
                D_01 = d_side_data_diffusivities->getPointer(2, 1);
                D_02 = d_side_data_diffusivities->getPointer(2, 2);
                D_03 = d_side_data_diffusivities->getPointer(2, 3);
                D_04 = d_side_data_diffusivities->getPointer(2, 4);
                D_05 = d_side_data_diffusivities->getPointer(2, 5);
                D_06 = d_side_data_diffusivities->getPointer(2, 6);
                D_07 = d_side_data_diffusivities->getPointer(2, 7);
                
                for (int k = -num_ghosts[2]; k < interior_dims[2] + num_ghosts[2] + 1; k++)
                {
                    for (int j = -num_ghosts[1]; j < interior_dims[1] + num_ghosts[1]; j++)
                    {
                        for (int i = -num_ghosts[0]; i < interior_dims[0] + num_ghosts[0]; i++)
                        {
                            const int idx_diffusivities = (i + d_num_subghosts_diffusivities[0]) +
                                (j + d_num_subghosts_diffusivities[1])*d_subghostcell_dims_diffusivities[0] +
                                (k + d_num_subghosts_diffusivities[2])*d_subghostcell_dims_diffusivities[0]*
                                    d_subghostcell_dims_diffusivities[1];
                            
                            const int idx_var_data = (i + num_ghosts[0]) +
                                (j + num_ghosts[1])*ghostcell_dims[0] +
                                (k + num_ghosts[2])*ghostcell_dims[0]*
                                    ghostcell_dims[1];
                            
                            D_00[idx_diffusivities] = -(Real(4)/Real(3)*mu_z[idx_var_data] + mu_v_z[idx_var_data]);
                            D_01[idx_diffusivities] = Real(2)/Real(3)*mu_z[idx_var_data] - mu_v_z[idx_var_data];
                            D_02[idx_diffusivities] = -mu_z[idx_var_data];
                            D_03[idx_diffusivities] = -w_z[idx_var_data]*(Real(4)/Real(3)*mu_z[idx_var_data] +
                                mu_v_z[idx_var_data]);
                            D_04[idx_diffusivities] = w_z[idx_var_data]*(Real(2)/Real(3)*mu_z[idx_var_data] -
                                mu_v_z[idx_var_data]);
                            D_05[idx_diffusivities] = -u_z[idx_var_data]*mu_z[idx_var_data];
                            D_06[idx_diffusivities] = -v_z[idx_var_data]*mu_z[idx_var_data];
                            D_07[idx_diffusivities] = -kappa_z[idx_var_data];
                        }
                    }
                }
            }
            
            d_side_data_diffusivities_computed = true;
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::computeSideDataOfDiffusiveFluxDiffusivities()\n"
            << "Cell data of 'DIFFUSIVITIES' is not yet registered."
            << std::endl);
    }
}


/*
 * Get the side data of the diffusivities in the diffusive fluxa.
 */
void
FlowModelDiffusiveFluxUtilitiesSingleSpecies::getSideDataOfDiffusiveFluxDiffusivities(
    std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > >& diffusivities_data,
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
            << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::"
            << "getSideDataOfDiffusiveFluxDiffusivities()\n"
            << "No patch is registered yet."
            << std::endl);
    }
    
    if (!d_side_data_diffusivities_computed)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::getSideDataOfDiffusiveFluxDiffusivities()\n"
            << "Side data of 'DIFFUSIVITIES' is not registered/computed yet."
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
                         * Mass equation.
                         */
                        
                        diffusivities_data[0].resize(0);
                        diffusivities_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[1].resize(1);
                        diffusivities_component_idx[1].resize(1);
                        
                        // -(4/3*mu + mu_v).
                        diffusivities_data[1][0] = d_side_data_diffusivities;
                        diffusivities_component_idx[1][0] = 0;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[2].resize(2);
                        diffusivities_component_idx[2].resize(2);
                        
                        // -u*(4/3*mu + mu_v).
                        diffusivities_data[2][0] = d_side_data_diffusivities;
                        diffusivities_component_idx[2][0] = 1;
                        
                        // -kappa.
                        diffusivities_data[2][1] = d_side_data_diffusivities;
                        diffusivities_component_idx[2][1] = 2;
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::"
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
                    << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::getSideDataOfDiffusiveFluxDiffusivities()\n"
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
                         * Mass equation.
                         */
                        
                        diffusivities_data[0].resize(0);
                        diffusivities_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[1].resize(1);
                        diffusivities_component_idx[1].resize(1);
                        
                        // -(4/3*mu + mu_v).
                        diffusivities_data[1][0] = d_side_data_diffusivities;
                        diffusivities_component_idx[1][0] = 0;
                        
                        diffusivities_data[2].resize(1);
                        diffusivities_component_idx[2].resize(1);
                        
                        // -mu.
                        diffusivities_data[2][0] = d_side_data_diffusivities;
                        diffusivities_component_idx[2][0] = 2;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[3].resize(3);
                        diffusivities_component_idx[3].resize(3);
                        
                        // -u*(4/3*mu + mu_v).
                        diffusivities_data[3][0] = d_side_data_diffusivities;
                        diffusivities_component_idx[3][0] = 3;
                        
                        // -v*mu.
                        diffusivities_data[3][1] = d_side_data_diffusivities;
                        diffusivities_component_idx[3][1] = 5;
                        
                        // -kappa.
                        diffusivities_data[3][2] = d_side_data_diffusivities;
                        diffusivities_component_idx[3][2] = 6;
                        
                        break;
                    }
                    case DIRECTION::Y_DIRECTION:
                    {
                        /*
                         * Mass equation.
                         */
                        
                        diffusivities_data[0].resize(0);
                        diffusivities_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[1].resize(1);
                        diffusivities_component_idx[1].resize(1);
                        
                        // 2/3*mu - mu_v.
                        diffusivities_data[1][0] = d_side_data_diffusivities;
                        diffusivities_component_idx[1][0] = 1;
                        
                        diffusivities_data[2].resize(1);
                        diffusivities_component_idx[2].resize(1);
                        
                        // -mu.
                        diffusivities_data[2][0] = d_side_data_diffusivities;
                        diffusivities_component_idx[2][0] = 2;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[3].resize(2);
                        diffusivities_component_idx[3].resize(2);
                        
                        // -v*mu.
                        diffusivities_data[3][0] = d_side_data_diffusivities;
                        diffusivities_component_idx[3][0] = 5;
                        
                        // u*(2/3*mu - mu_v).
                        diffusivities_data[3][1] = d_side_data_diffusivities;
                        diffusivities_component_idx[3][1] = 4;
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::"
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
                        /*
                         * Mass equation.
                         */
                        
                        diffusivities_data[0].resize(0);
                        diffusivities_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[1].resize(1);
                        diffusivities_component_idx[1].resize(1);
                        
                        // -mu.
                        diffusivities_data[1][0] = d_side_data_diffusivities;
                        diffusivities_component_idx[1][0] = 2;
                        
                        diffusivities_data[2].resize(1);
                        diffusivities_component_idx[2].resize(1);
                        
                        // 2/3*mu - mu_v.
                        diffusivities_data[2][0] = d_side_data_diffusivities;
                        diffusivities_component_idx[2][0] = 1;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[3].resize(2);
                        diffusivities_component_idx[3].resize(2);
                        
                        // v*(2/3*mu - mu_v).
                        diffusivities_data[3][0] = d_side_data_diffusivities;
                        diffusivities_component_idx[3][0] = 4;
                        
                        // -u*mu.
                        diffusivities_data[3][1] = d_side_data_diffusivities;
                        diffusivities_component_idx[3][1] = 5;
                        
                        break;
                    }
                    case DIRECTION::Y_DIRECTION:
                    {
                        /*
                         * Mass equation.
                         */
                        
                        diffusivities_data[0].resize(0);
                        diffusivities_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[1].resize(1);
                        diffusivities_component_idx[1].resize(1);
                        
                        // -mu.
                        diffusivities_data[1][0] = d_side_data_diffusivities;
                        diffusivities_component_idx[1][0] = 2;
                        
                        diffusivities_data[2].resize(1);
                        diffusivities_component_idx[2].resize(1);
                        
                        // -(4/3*mu + mu_v).
                        diffusivities_data[2][0] = d_side_data_diffusivities;
                        diffusivities_component_idx[2][0] = 0;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[3].resize(3);
                        diffusivities_component_idx[3].resize(3);
                        
                        // -u*mu.
                        diffusivities_data[3][0] = d_side_data_diffusivities;
                        diffusivities_component_idx[3][0] = 5;
                        
                        // -v*(4/3*mu + mu_v).
                        diffusivities_data[3][1] = d_side_data_diffusivities;
                        diffusivities_component_idx[3][1] = 3;
                        
                        // -kappa.
                        diffusivities_data[3][2] = d_side_data_diffusivities;
                        diffusivities_component_idx[3][2] = 6;
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::"
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
                    << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::getSideDataOfDiffusiveFluxDiffusivities()\n"
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
                         * Mass equation.
                         */
                        
                        diffusivities_data[0].resize(0);
                        diffusivities_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[1].resize(1);
                        diffusivities_component_idx[1].resize(1);
                        
                        // -(4/3*mu + mu_v).
                        diffusivities_data[1][0] = d_side_data_diffusivities;
                        diffusivities_component_idx[1][0] = 0;
                        
                        diffusivities_data[2].resize(1);
                        diffusivities_component_idx[2].resize(1);
                        
                        // -mu.
                        diffusivities_data[2][0] = d_side_data_diffusivities;
                        diffusivities_component_idx[2][0] = 2;
                        
                        diffusivities_data[3].resize(1);
                        diffusivities_component_idx[3].resize(1);
                        
                        // -mu.
                        diffusivities_data[3][0] = d_side_data_diffusivities;
                        diffusivities_component_idx[3][0] = 2;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[4].resize(4);
                        diffusivities_component_idx[4].resize(4);
                        
                        // -u*(4/3*mu + mu_v).
                        diffusivities_data[4][0] = d_side_data_diffusivities;
                        diffusivities_component_idx[4][0] = 3;
                        
                        // -v*mu.
                        diffusivities_data[4][1] = d_side_data_diffusivities;
                        diffusivities_component_idx[4][1] = 5;
                        
                        // -w*mu.
                        diffusivities_data[4][2] = d_side_data_diffusivities;
                        diffusivities_component_idx[4][2] = 6;
                        
                        // -kappa.
                        diffusivities_data[4][3] = d_side_data_diffusivities;
                        diffusivities_component_idx[4][3] = 7;
                        
                        break;
                    }
                    case DIRECTION::Y_DIRECTION:
                    {
                        /*
                         * Mass equation.
                         */
                        
                        diffusivities_data[0].resize(0);
                        diffusivities_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[1].resize(1);
                        diffusivities_component_idx[1].resize(1);
                        
                        // 2/3*mu - mu_v.
                        diffusivities_data[1][0] = d_side_data_diffusivities;
                        diffusivities_component_idx[1][0] = 1;
                        
                        diffusivities_data[2].resize(1);
                        diffusivities_component_idx[2].resize(1);
                        
                        // -mu.
                        diffusivities_data[2][0] = d_side_data_diffusivities;
                        diffusivities_component_idx[2][0] = 2;
                        
                        diffusivities_data[3].resize(0);
                        diffusivities_component_idx[3].resize(0);
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[4].resize(2);
                        diffusivities_component_idx[4].resize(2);
                        
                        // -v*mu.
                        diffusivities_data[4][0] = d_side_data_diffusivities;
                        diffusivities_component_idx[4][0] = 5;
                        
                        // u*(2/3*mu - mu_v).
                        diffusivities_data[4][1] = d_side_data_diffusivities;
                        diffusivities_component_idx[4][1] = 4;
                        
                        break;
                    }
                    case DIRECTION::Z_DIRECTION:
                    {
                        /*
                         * Mass equation.
                         */
                        
                        diffusivities_data[0].resize(0);
                        diffusivities_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[1].resize(1);
                        diffusivities_component_idx[1].resize(1);
                        
                        // 2/3*mu - mu_v.
                        diffusivities_data[1][0] = d_side_data_diffusivities;
                        diffusivities_component_idx[1][0] = 1;
                        
                        diffusivities_data[2].resize(0);
                        diffusivities_component_idx[2].resize(0);
                        
                        diffusivities_data[3].resize(1);
                        diffusivities_component_idx[3].resize(1);
                        
                        // -mu.
                        diffusivities_data[3][0] = d_side_data_diffusivities;
                        diffusivities_component_idx[3][0] = 2;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[4].resize(2);
                        diffusivities_component_idx[4].resize(2);
                        
                        // -w*mu.
                        diffusivities_data[4][0] = d_side_data_diffusivities;
                        diffusivities_component_idx[4][0] = 6;
                        
                        // u*(2/3*mu - mu_v).
                        diffusivities_data[4][1] = d_side_data_diffusivities;
                        diffusivities_component_idx[4][1] = 4;
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::"
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
                        /*
                         * Mass equation.
                         */
                        
                        diffusivities_data[0].resize(0);
                        diffusivities_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[1].resize(1);
                        diffusivities_component_idx[1].resize(1);
                        
                        // -mu.
                        diffusivities_data[1][0] = d_side_data_diffusivities;
                        diffusivities_component_idx[1][0] = 2;
                        
                        diffusivities_data[2].resize(1);
                        diffusivities_component_idx[2].resize(1);
                        
                        // 2/3*mu - mu_v.
                        diffusivities_data[2][0] = d_side_data_diffusivities;
                        diffusivities_component_idx[2][0] = 1;
                        
                        diffusivities_data[3].resize(0);
                        diffusivities_component_idx[3].resize(0);
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[4].resize(2);
                        diffusivities_component_idx[4].resize(2);
                        
                        // v*(2/3*mu - mu_v).
                        diffusivities_data[4][0] = d_side_data_diffusivities;
                        diffusivities_component_idx[4][0] = 4;
                        
                        // -u*mu.
                        diffusivities_data[4][1] = d_side_data_diffusivities;
                        diffusivities_component_idx[4][1] = 5;
                        
                        break;
                    }
                    case DIRECTION::Y_DIRECTION:
                    {
                        /*
                         * Mass equation.
                         */
                        
                        diffusivities_data[0].resize(0);
                        diffusivities_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[1].resize(1);
                        diffusivities_component_idx[1].resize(1);
                        
                        // -mu.
                        diffusivities_data[1][0] = d_side_data_diffusivities;
                        diffusivities_component_idx[1][0] = 2;
                        
                        diffusivities_data[2].resize(1);
                        diffusivities_component_idx[2].resize(1);
                        
                        // -(4/3*mu + mu_v).
                        diffusivities_data[2][0] = d_side_data_diffusivities;
                        diffusivities_component_idx[2][0] = 0;
                        
                        diffusivities_data[3].resize(1);
                        diffusivities_component_idx[3].resize(1);
                        
                        // -mu.
                        diffusivities_data[3][0] = d_side_data_diffusivities;
                        diffusivities_component_idx[3][0] = 2;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[4].resize(4);
                        diffusivities_component_idx[4].resize(4);
                        
                        // -u*mu.
                        diffusivities_data[4][0] = d_side_data_diffusivities;
                        diffusivities_component_idx[4][0] = 5;
                        
                        // -v*(4/3*mu + mu_v).
                        diffusivities_data[4][1] = d_side_data_diffusivities;
                        diffusivities_component_idx[4][1] = 3;
                        
                        // -w*mu.
                        diffusivities_data[4][2] = d_side_data_diffusivities;
                        diffusivities_component_idx[4][2] = 6;
                        
                        // -kappa.
                        diffusivities_data[4][3] = d_side_data_diffusivities;
                        diffusivities_component_idx[4][3] = 7;
                        
                        break;
                    }
                    case DIRECTION::Z_DIRECTION:
                    {
                        /*
                         * Mass equation.
                         */
                        
                        diffusivities_data[0].resize(0);
                        diffusivities_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[1].resize(0);
                        diffusivities_component_idx[1].resize(0);
                        
                        diffusivities_data[2].resize(1);
                        diffusivities_component_idx[2].resize(1);
                        
                        // 2/3*mu - mu_v.
                        diffusivities_data[2][0] = d_side_data_diffusivities;
                        diffusivities_component_idx[2][0] = 1;
                        
                        diffusivities_data[3].resize(1);
                        diffusivities_component_idx[3].resize(1);
                        
                        // -mu.
                        diffusivities_data[3][0] = d_side_data_diffusivities;
                        diffusivities_component_idx[3][0] = 2;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[4].resize(2);
                        diffusivities_component_idx[4].resize(2);
                        
                        // -w*mu.
                        diffusivities_data[4][0] = d_side_data_diffusivities;
                        diffusivities_component_idx[4][0] = 6;
                        
                        // v*(2/3*mu - mu_v).
                        diffusivities_data[4][1] = d_side_data_diffusivities;
                        diffusivities_component_idx[4][1] = 4;
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::"
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
                        /*
                         * Mass equation.
                         */
                        
                        diffusivities_data[0].resize(0);
                        diffusivities_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[1].resize(1);
                        diffusivities_component_idx[1].resize(1);
                        
                        // -mu.
                        diffusivities_data[1][0] = d_side_data_diffusivities;
                        diffusivities_component_idx[1][0] = 2;
                        
                        diffusivities_data[2].resize(0);
                        diffusivities_component_idx[2].resize(0);
                        
                        diffusivities_data[3].resize(1);
                        diffusivities_component_idx[3].resize(1);
                        
                        // 2/3*mu - mu_v.
                        diffusivities_data[3][0] = d_side_data_diffusivities;
                        diffusivities_component_idx[3][0] = 1;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[4].resize(2);
                        diffusivities_component_idx[4].resize(2);
                        
                        // w*(2/3*mu - mu_v).
                        diffusivities_data[4][0] = d_side_data_diffusivities;
                        diffusivities_component_idx[4][0] = 4;
                        
                        // -u*mu.
                        diffusivities_data[4][1] = d_side_data_diffusivities;
                        diffusivities_component_idx[4][1] = 5;
                        
                        break;
                    }
                    case DIRECTION::Y_DIRECTION:
                    {
                        /*
                         * Mass equation.
                         */
                        
                        diffusivities_data[0].resize(0);
                        diffusivities_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[1].resize(0);
                        diffusivities_component_idx[1].resize(0);
                        
                        diffusivities_data[2].resize(1);
                        diffusivities_component_idx[2].resize(1);
                        
                        // -mu.
                        diffusivities_data[2][0] = d_side_data_diffusivities;
                        diffusivities_component_idx[2][0] = 2;
                        
                        diffusivities_data[3].resize(1);
                        diffusivities_component_idx[3].resize(1);
                        
                        // 2/3*mu - mu_v.
                        diffusivities_data[3][0] = d_side_data_diffusivities;
                        diffusivities_component_idx[3][0] = 1;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[4].resize(2);
                        diffusivities_component_idx[4].resize(2);
                        
                        // w*(2/3*mu - mu_v).
                        diffusivities_data[4][0] = d_side_data_diffusivities;
                        diffusivities_component_idx[4][0] = 4;
                        
                        // -v*mu.
                        diffusivities_data[4][1] = d_side_data_diffusivities;
                        diffusivities_component_idx[4][1] = 6;
                        
                        break;
                    }
                    case DIRECTION::Z_DIRECTION:
                    {
                        /*
                         * Mass equation.
                         */
                        
                        diffusivities_data[0].resize(0);
                        diffusivities_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[1].resize(1);
                        diffusivities_component_idx[1].resize(1);
                        
                        // -mu.
                        diffusivities_data[1][0] = d_side_data_diffusivities;
                        diffusivities_component_idx[1][0] = 2;
                        
                        diffusivities_data[2].resize(1);
                        diffusivities_component_idx[2].resize(1);
                        
                        // -mu.
                        diffusivities_data[2][0] = d_side_data_diffusivities;
                        diffusivities_component_idx[2][0] = 2;
                        
                        diffusivities_data[3].resize(1);
                        diffusivities_component_idx[3].resize(1);
                        
                        // -(4/3*mu + mu_v).
                        diffusivities_data[3][0] = d_side_data_diffusivities;
                        diffusivities_component_idx[3][0] = 0;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[4].resize(4);
                        diffusivities_component_idx[4].resize(4);
                        
                        // -u*mu.
                        diffusivities_data[4][0] = d_side_data_diffusivities;
                        diffusivities_component_idx[4][0] = 5;
                        
                        // -v*mu.
                        diffusivities_data[4][1] = d_side_data_diffusivities;
                        diffusivities_component_idx[4][1] = 6;
                        
                        // -w*(4/3*mu + mu_v).
                        diffusivities_data[4][2] = d_side_data_diffusivities;
                        diffusivities_component_idx[4][2] = 3;
                        
                        // -kappa.
                        diffusivities_data[4][3] = d_side_data_diffusivities;
                        diffusivities_component_idx[4][3] = 7;
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::"
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
                    << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::getSideDataOfDiffusiveFluxDiffusivities()\n"
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
FlowModelDiffusiveFluxUtilitiesSingleSpecies::setNumberOfSubGhosts(
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
    
    if (variable_name == "THERMAL_CONDUCTIVITY")
    {
        if (d_num_subghosts_thermal_conductivity > -hier::IntVector::getOne(d_dim))
        {
            if (num_subghosts > d_num_subghosts_thermal_conductivity)
            {
                d_num_subghosts_thermal_conductivity = num_subghosts;
            }
        }
        else
        {
            d_num_subghosts_thermal_conductivity = num_subghosts;
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
        setNumberOfSubGhosts(num_subghosts, "THERMAL_CONDUCTIVITY", parent_variable_name);
    }
}


/*
 * Set the ghost boxes of derived cell variables.
 */
void
FlowModelDiffusiveFluxUtilitiesSingleSpecies::setDerivedCellVariableGhostBoxes()
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
    
    if (d_num_subghosts_thermal_conductivity > -hier::IntVector::getOne(d_dim))
    {
        d_subghost_box_thermal_conductivity = interior_box;
        d_subghost_box_thermal_conductivity.grow(d_num_subghosts_thermal_conductivity);
        d_subghostcell_dims_thermal_conductivity = d_subghost_box_thermal_conductivity.numberCells();
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
FlowModelDiffusiveFluxUtilitiesSingleSpecies::computeCellDataOfShearViscosity()
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
            
            // Get the cell data of pressure.
            HAMERS_SHARED_PTR<pdat::CellData<Real> > data_pressure =
                flow_model_tmp->getCellData("PRESSURE");
            
            // Get the cell data of temperature.
            HAMERS_SHARED_PTR<pdat::CellData<Real> > data_temperature =
                flow_model_tmp->getCellData("TEMPERATURE");
            
            /*
             * Get the molecular properties of the species for shear viscosity.
             */
            
            std::vector<Real> molecular_properties_shear_viscosity;
            std::vector<Real*> molecular_properties_shear_viscosity_ptr;
            std::vector<const Real*> molecular_properties_shear_viscosity_const_ptr;
            
            const int num_molecular_properties_shear_viscosity = d_equation_of_shear_viscosity_mixing_rules->
                getNumberOfSpeciesMolecularProperties();
            
            molecular_properties_shear_viscosity.resize(num_molecular_properties_shear_viscosity);
            molecular_properties_shear_viscosity_ptr.reserve(num_molecular_properties_shear_viscosity);
            molecular_properties_shear_viscosity_const_ptr.reserve(num_molecular_properties_shear_viscosity);
            
            for (int ti = 0; ti < num_molecular_properties_shear_viscosity; ti++)
            {
                molecular_properties_shear_viscosity_ptr.push_back(&molecular_properties_shear_viscosity[ti]);
                molecular_properties_shear_viscosity_const_ptr.push_back(&molecular_properties_shear_viscosity[ti]);
            }
            
            d_equation_of_shear_viscosity_mixing_rules->getSpeciesMolecularProperties(
                molecular_properties_shear_viscosity_ptr,
                0);
            
            // Compute the shear viscosity field.
            d_equation_of_shear_viscosity_mixing_rules->getEquationOfShearViscosity()->
                computeShearViscosity(
                    d_data_shear_viscosity,
                    data_pressure,
                    data_temperature,
                    molecular_properties_shear_viscosity_const_ptr,
                    empty_box);
            
            d_cell_data_computed_shear_viscosity = true;
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::computeCellDataOfShearViscosity()\n"
            << "Cell data of 'SHEAR_VISCOSITY' is not yet registered."
            << std::endl);
    }
}


/*
 * Compute the cell data of bulk viscosity in the registered patch.
 */
void
FlowModelDiffusiveFluxUtilitiesSingleSpecies::computeCellDataOfBulkViscosity()
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
            
            // Get the cell data of pressure.
            HAMERS_SHARED_PTR<pdat::CellData<Real> > data_pressure =
                flow_model_tmp->getCellData("PRESSURE");
            
            // Get the cell data of temperature.
            HAMERS_SHARED_PTR<pdat::CellData<Real> > data_temperature =
                flow_model_tmp->getCellData("TEMPERATURE");
            
            /*
             * Get the molecular properties of the species for bulk viscosity.
             */
            
            std::vector<Real> molecular_properties_bulk_viscosity;
            std::vector<Real*> molecular_properties_bulk_viscosity_ptr;
            std::vector<const Real*> molecular_properties_bulk_viscosity_const_ptr;
            
            const int num_molecular_properties_bulk_viscosity = d_equation_of_bulk_viscosity_mixing_rules->
                getNumberOfSpeciesMolecularProperties();
            
            molecular_properties_bulk_viscosity.resize(num_molecular_properties_bulk_viscosity);
            molecular_properties_bulk_viscosity_ptr.reserve(num_molecular_properties_bulk_viscosity);
            molecular_properties_bulk_viscosity_const_ptr.reserve(num_molecular_properties_bulk_viscosity);
            
            for (int ti = 0; ti < num_molecular_properties_bulk_viscosity; ti++)
            {
                molecular_properties_bulk_viscosity_ptr.push_back(&molecular_properties_bulk_viscosity[ti]);
                molecular_properties_bulk_viscosity_const_ptr.push_back(&molecular_properties_bulk_viscosity[ti]);
            }
            
            d_equation_of_bulk_viscosity_mixing_rules->getSpeciesMolecularProperties(
                molecular_properties_bulk_viscosity_ptr,
                0);
            
            // Compute the bulk viscosity field.
            d_equation_of_bulk_viscosity_mixing_rules->getEquationOfBulkViscosity()->
                computeBulkViscosity(
                    d_data_bulk_viscosity,
                    data_pressure,
                    data_temperature,
                    molecular_properties_bulk_viscosity_const_ptr,
                    empty_box);
            
            d_cell_data_computed_bulk_viscosity = true;
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::computeCellDataOfBulkViscosity()\n"
            << "Cell data of 'BULK_VISCOSITY' is not yet registered."
            << std::endl);
    }
}


/*
 * Compute the cell data of thermal conductivity in the registered patch.
 */
void
FlowModelDiffusiveFluxUtilitiesSingleSpecies::computeCellDataOfThermalConductivity()
{
    // Create empty box.
    const hier::Box empty_box(d_dim);
    
    if (d_num_subghosts_thermal_conductivity > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_thermal_conductivity)
        {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(d_data_thermal_conductivity);
#endif
            
            if (d_flow_model.expired())
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "The object is not setup yet!"
                    << std::endl);
            }
            
            HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
            
            // Get the cell data of pressure.
            HAMERS_SHARED_PTR<pdat::CellData<Real> > data_pressure =
                flow_model_tmp->getCellData("PRESSURE");
            
            // Get the cell data of temperature.
            HAMERS_SHARED_PTR<pdat::CellData<Real> > data_temperature =
                flow_model_tmp->getCellData("TEMPERATURE");
            
            /*
             * Get the molecular properties of the species for thermal conductivity.
             */
            
            std::vector<Real> molecular_properties_thermal_conductivity;
            std::vector<Real*> molecular_properties_thermal_conductivity_ptr;
            std::vector<const Real*> molecular_properties_thermal_conductivity_const_ptr;
            
            const int num_molecular_properties_thermal_conductivity = d_equation_of_thermal_conductivity_mixing_rules->
                getNumberOfSpeciesMolecularProperties();
            
            molecular_properties_thermal_conductivity.resize(num_molecular_properties_thermal_conductivity);
            molecular_properties_thermal_conductivity_ptr.reserve(num_molecular_properties_thermal_conductivity);
            molecular_properties_thermal_conductivity_const_ptr.reserve(num_molecular_properties_thermal_conductivity);
            
            for (int ti = 0; ti < num_molecular_properties_thermal_conductivity; ti++)
            {
                molecular_properties_thermal_conductivity_ptr.push_back(&molecular_properties_thermal_conductivity[ti]);
                molecular_properties_thermal_conductivity_const_ptr.push_back(&molecular_properties_thermal_conductivity[ti]);
            }
            
            d_equation_of_thermal_conductivity_mixing_rules->getSpeciesMolecularProperties(
                molecular_properties_thermal_conductivity_ptr,
                0);
            
            // Compute the thermal conductivity field.
            d_equation_of_thermal_conductivity_mixing_rules->getEquationOfThermalConductivity()->
                computeThermalConductivity(
                    d_data_thermal_conductivity,
                    data_pressure,
                    data_temperature,
                    molecular_properties_thermal_conductivity_const_ptr,
                    empty_box);
            
            d_cell_data_computed_thermal_conductivity = true;
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::computeCellDataOfThermalConductivity()\n"
            << "Cell data of 'THERMAL_CONDUCTIVITY' is not yet registered."
            << std::endl);
    }
}


/*
 * Compute the cell data of diffusivities in the registered patch.
 */
void
FlowModelDiffusiveFluxUtilitiesSingleSpecies::computeCellDataOfDiffusivities()
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
            
            if (!d_cell_data_computed_thermal_conductivity)
            {
                computeCellDataOfThermalConductivity();
            }
            
            // Get the cell data of velocity.
            HAMERS_SHARED_PTR<pdat::CellData<Real> > data_velocity =
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
             * Get the pointers to the cell data of shear viscosity, bulk viscosity and thermal conductivity.
             */
            
            Real* mu    = d_data_shear_viscosity->getPointer(0);
            Real* mu_v  = d_data_bulk_viscosity->getPointer(0);
            Real* kappa = d_data_thermal_conductivity->getPointer(0);
            
            if (d_dim == tbox::Dimension(1))
            {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(d_data_diffusivities->getDepth() == 3);
#endif
                
                // Get the pointer to cell data of velocity and diffusivities.
                Real* u = data_velocity->getPointer(0);
                
                Real* D_00 = d_data_diffusivities->getPointer(0);
                Real* D_01 = d_data_diffusivities->getPointer(1);
                Real* D_02 = d_data_diffusivities->getPointer(2);
                
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
                    const int idx_thermal_conductivity = i + d_num_subghosts_thermal_conductivity[0];
                    const int idx_velocity = i + num_subghosts_velocity[0];
                    
                    D_00[idx_diffusivities] = -(Real(4)/Real(3)*mu[idx_shear_viscosity] + mu_v[idx_bulk_viscosity]);
                    D_01[idx_diffusivities] = -u[idx_velocity]*(Real(4)/Real(3)*mu[idx_shear_viscosity] +
                        mu_v[idx_bulk_viscosity]);
                    D_02[idx_diffusivities] = -kappa[idx_thermal_conductivity];
                }
            }
            else if (d_dim == tbox::Dimension(2))
            {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(d_data_diffusivities->getDepth() == 10);
#endif
                
                // Get the pointer to cell data of velocity and diffusivities.
                Real* u = data_velocity->getPointer(0);
                Real* v = data_velocity->getPointer(1);
                
                Real* D_00 = d_data_diffusivities->getPointer(0);
                Real* D_01 = d_data_diffusivities->getPointer(1);
                Real* D_02 = d_data_diffusivities->getPointer(2);
                Real* D_03 = d_data_diffusivities->getPointer(3);
                Real* D_04 = d_data_diffusivities->getPointer(4);
                Real* D_05 = d_data_diffusivities->getPointer(5);
                Real* D_06 = d_data_diffusivities->getPointer(6);
                Real* D_07 = d_data_diffusivities->getPointer(7);
                Real* D_08 = d_data_diffusivities->getPointer(8);
                Real* D_09 = d_data_diffusivities->getPointer(9);
                
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
                        
                        const int idx_thermal_conductivity = (i + d_num_subghosts_thermal_conductivity[0]) +
                            (j + d_num_subghosts_thermal_conductivity[1])*d_subghostcell_dims_thermal_conductivity[0];
                        
                        const int idx_velocity = (i + num_subghosts_velocity[0]) +
                            (j + num_subghosts_velocity[1])*subghostcell_dims_velocity[0];
                        
                        D_00[idx_diffusivities] = -(Real(4)/Real(3)*mu[idx_shear_viscosity] + mu_v[idx_bulk_viscosity]);
                        D_01[idx_diffusivities] = Real(2)/Real(3)*mu[idx_shear_viscosity] - mu_v[idx_bulk_viscosity];
                        D_02[idx_diffusivities] = -mu[idx_shear_viscosity];
                        D_03[idx_diffusivities] = -u[idx_velocity]*(Real(4)/Real(3)*mu[idx_shear_viscosity] +
                            mu_v[idx_bulk_viscosity]);
                        D_04[idx_diffusivities] = -v[idx_velocity]*(Real(4)/Real(3)*mu[idx_shear_viscosity] +
                            mu_v[idx_bulk_viscosity]);
                        D_05[idx_diffusivities] = u[idx_velocity]*(Real(2)/Real(3)*mu[idx_shear_viscosity] -
                            mu_v[idx_bulk_viscosity]);
                        D_06[idx_diffusivities] = v[idx_velocity]*(Real(2)/Real(3)*mu[idx_shear_viscosity] -
                            mu_v[idx_bulk_viscosity]);
                        D_07[idx_diffusivities] = -u[idx_velocity]*mu[idx_shear_viscosity];
                        D_08[idx_diffusivities] = -v[idx_velocity]*mu[idx_shear_viscosity];
                        D_09[idx_diffusivities] = -kappa[idx_thermal_conductivity];
                    }
                }
            }
            else if (d_dim == tbox::Dimension(3))
            {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(d_data_diffusivities->getDepth() == 13);
#endif
                
                // Get the pointer to cell data of velocity and diffusivities.
                Real* u = data_velocity->getPointer(0);
                Real* v = data_velocity->getPointer(1);
                Real* w = data_velocity->getPointer(2);
                
                Real* D_00 = d_data_diffusivities->getPointer(0);
                Real* D_01 = d_data_diffusivities->getPointer(1);
                Real* D_02 = d_data_diffusivities->getPointer(2);
                Real* D_03 = d_data_diffusivities->getPointer(3);
                Real* D_04 = d_data_diffusivities->getPointer(4);
                Real* D_05 = d_data_diffusivities->getPointer(5);
                Real* D_06 = d_data_diffusivities->getPointer(6);
                Real* D_07 = d_data_diffusivities->getPointer(7);
                Real* D_08 = d_data_diffusivities->getPointer(8);
                Real* D_09 = d_data_diffusivities->getPointer(9);
                Real* D_10 = d_data_diffusivities->getPointer(10);
                Real* D_11 = d_data_diffusivities->getPointer(11);
                Real* D_12 = d_data_diffusivities->getPointer(12);
                
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
                            
                            const int idx_thermal_conductivity = (i + d_num_subghosts_thermal_conductivity[0]) +
                                (j + d_num_subghosts_thermal_conductivity[1])*d_subghostcell_dims_thermal_conductivity[0] +
                                (k + d_num_subghosts_thermal_conductivity[2])*d_subghostcell_dims_thermal_conductivity[0]*
                                    d_subghostcell_dims_thermal_conductivity[1];
                            
                            const int idx_velocity = (i + num_subghosts_velocity[0]) +
                                (j + num_subghosts_velocity[1])*subghostcell_dims_velocity[0] +
                                (k + num_subghosts_velocity[2])*subghostcell_dims_velocity[0]*
                                    subghostcell_dims_velocity[1];
                            
                            D_00[idx_diffusivities] = -(Real(4)/Real(3)*mu[idx_shear_viscosity] + mu_v[idx_bulk_viscosity]);
                            D_01[idx_diffusivities] = Real(2)/Real(3)*mu[idx_shear_viscosity] - mu_v[idx_bulk_viscosity];
                            D_02[idx_diffusivities] = -mu[idx_shear_viscosity];
                            D_03[idx_diffusivities] = -u[idx_velocity]*(Real(4)/Real(3)*mu[idx_shear_viscosity] +
                                mu_v[idx_bulk_viscosity]);
                            D_04[idx_diffusivities] = -v[idx_velocity]*(Real(4)/Real(3)*mu[idx_shear_viscosity] +
                                mu_v[idx_bulk_viscosity]);
                            D_05[idx_diffusivities] = -w[idx_velocity]*(Real(4)/Real(3)*mu[idx_shear_viscosity] +
                                mu_v[idx_bulk_viscosity]);
                            D_06[idx_diffusivities] = u[idx_velocity]*(Real(2)/Real(3)*mu[idx_shear_viscosity] -
                                mu_v[idx_bulk_viscosity]);
                            D_07[idx_diffusivities] = v[idx_velocity]*(Real(2)/Real(3)*mu[idx_shear_viscosity] -
                                mu_v[idx_bulk_viscosity]);
                            D_08[idx_diffusivities] = w[idx_velocity]*(Real(2)/Real(3)*mu[idx_shear_viscosity] -
                                mu_v[idx_bulk_viscosity]);
                            D_09[idx_diffusivities] = -u[idx_velocity]*mu[idx_shear_viscosity];
                            D_10[idx_diffusivities] = -v[idx_velocity]*mu[idx_shear_viscosity];
                            D_11[idx_diffusivities] = -w[idx_velocity]*mu[idx_shear_viscosity];
                            D_12[idx_diffusivities] = -kappa[idx_thermal_conductivity];
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
            << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::computeCellDataOfDiffusivities()\n"
            << "Cell data of 'DIFFUSIVITIES' is not yet registered."
            << std::endl);
    }
}
