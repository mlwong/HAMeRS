#include "flow/flow_models/single-species/FlowModelDiffusiveFluxUtilitiesSingleSpecies.hpp"

FlowModelDiffusiveFluxUtilitiesSingleSpecies::FlowModelDiffusiveFluxUtilitiesSingleSpecies(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
    const int& num_species,
    const boost::shared_ptr<EquationOfShearViscosityMixingRules> equation_of_shear_viscosity_mixing_rules,
    const boost::shared_ptr<EquationOfBulkViscosityMixingRules> equation_of_bulk_viscosity_mixing_rules,
    const boost::shared_ptr<EquationOfThermalConductivityMixingRules> equation_of_thermal_conductivity_mixing_rules):
        FlowModelDiffusiveFluxUtilities(
            object_name,
            dim,
            grid_geometry,
            num_species,
            2 + dim.getValue()),
        d_num_subghosts_shear_viscosity(-hier::IntVector::getOne(d_dim)),
        d_num_subghosts_bulk_viscosity(-hier::IntVector::getOne(d_dim)),
        d_num_subghosts_thermal_conductivity(-hier::IntVector::getOne(d_dim)),
        d_subghost_box_shear_viscosity(hier::Box::getEmptyBox(d_dim)),
        d_subghost_box_bulk_viscosity(hier::Box::getEmptyBox(d_dim)),
        d_subghost_box_thermal_conductivity(hier::Box::getEmptyBox(d_dim)),
        d_subghostcell_dims_shear_viscosity(hier::IntVector::getZero(d_dim)),
        d_subghostcell_dims_bulk_viscosity(hier::IntVector::getZero(d_dim)),
        d_subghostcell_dims_thermal_conductivity(hier::IntVector::getZero(d_dim)),
        d_cell_data_shear_viscosity_computed(false),
        d_cell_data_bulk_viscosity_computed(false),
        d_cell_data_thermal_conductivity_computed(false),
        d_equation_of_shear_viscosity_mixing_rules(equation_of_shear_viscosity_mixing_rules),
        d_equation_of_bulk_viscosity_mixing_rules(equation_of_bulk_viscosity_mixing_rules),
        d_equation_of_thermal_conductivity_mixing_rules(equation_of_thermal_conductivity_mixing_rules)
{}


/*
 * Register the required variables for the computation of diffusive fluxes in the registered patch.
 */
void
FlowModelDiffusiveFluxUtilitiesSingleSpecies::registerDerivedVariablesForDiffusiveFluxes(
    const hier::IntVector& num_subghosts)
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
    
    /*
     * Get the interior box.
     */
    
    const hier::Box interior_box = patch.getBox();
    
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
    
    /*
     * Set ghost boxes of derived cell variables for this class.
     */
    
    d_num_subghosts_diffusivities = num_subghosts;
    d_subghost_box_diffusivities = interior_box;
    d_subghost_box_diffusivities.grow(d_num_subghosts_diffusivities);
    d_subghostcell_dims_diffusivities = d_subghost_box_diffusivities.numberCells();
    
    d_num_subghosts_shear_viscosity      = d_num_subghosts_diffusivities;
    d_num_subghosts_bulk_viscosity       = d_num_subghosts_diffusivities;
    d_num_subghosts_thermal_conductivity = d_num_subghosts_diffusivities;
    
    d_subghost_box_shear_viscosity      = d_subghost_box_diffusivities;
    d_subghost_box_bulk_viscosity       = d_subghost_box_diffusivities;
    d_subghost_box_thermal_conductivity = d_subghost_box_diffusivities;
    
    d_subghostcell_dims_shear_viscosity      = d_subghostcell_dims_diffusivities;
    d_subghostcell_dims_bulk_viscosity       = d_subghostcell_dims_diffusivities;
    d_subghostcell_dims_thermal_conductivity = d_subghostcell_dims_diffusivities;
}


/*
 * Allocate memory for cell data of all registered derived variables in the registered patch for this class.
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
    
    boost::shared_ptr<FlowModel> flow_model_tmp = d_flow_model.lock();
    const hier::Patch& patch = flow_model_tmp->getRegisteredPatch();
    const hier::Box interior_box = patch.getBox();
    
    if (d_num_subghosts_shear_viscosity > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_shear_viscosity_computed)
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
                << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::allocateMemoryForDerivedCellData()\n"
                << "Cell data of 'SHEAR_VISCOSITY' is aleady computed."
                << std::endl);
        }
    }
    
    if (d_num_subghosts_bulk_viscosity > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_bulk_viscosity_computed)
        {
            if (!d_data_bulk_viscosity)
            {
                d_data_bulk_viscosity.reset(new pdat::CellData<double>(
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
        if (!d_cell_data_thermal_conductivity_computed)
        {
            if (!d_data_thermal_conductivity)
            {
                d_data_thermal_conductivity.reset(new pdat::CellData<double>(
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
    
    if (d_num_subghosts_diffusivities > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_diffusivities_computed)
        {
            if (!d_data_diffusivities)
            {
                if (d_dim == tbox::Dimension(1))
                {
                    d_data_diffusivities.reset(new pdat::CellData<double>(
                        interior_box, 3, d_num_subghosts_diffusivities));
                }
                else if (d_dim == tbox::Dimension(2))
                {
                    d_data_diffusivities.reset(new pdat::CellData<double>(
                        interior_box, 10, d_num_subghosts_diffusivities));
                }
                else if (d_dim == tbox::Dimension(3))
                {
                    d_data_diffusivities.reset(new pdat::CellData<double>(
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


/*
 * Clear cell data of all derived variables in the registered patch for this class.
 */
void
FlowModelDiffusiveFluxUtilitiesSingleSpecies::clearCellData()
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
    
    d_cell_data_diffusivities_computed        = false;
    d_cell_data_shear_viscosity_computed      = false;
    d_cell_data_bulk_viscosity_computed       = false;
    d_cell_data_thermal_conductivity_computed = false;
    
    d_derived_cell_data_computed = false;
}


/*
 * Get the variables for the derivatives in the diffusive fluxes.
 */
void
FlowModelDiffusiveFluxUtilitiesSingleSpecies::getCellDataOfDiffusiveFluxVariablesForDerivative(
    std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > >& derivative_var_data,
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
    
    boost::shared_ptr<FlowModel> flow_model_tmp = d_flow_model.lock();
    
    derivative_var_data.resize(d_num_eqn);
    derivative_var_component_idx.resize(d_num_eqn);
    
    // Get the cell data of velocity.
    boost::shared_ptr<pdat::CellData<double> > data_velocity =
        flow_model_tmp->getCellData("VELOCITY");
    
    // Get the cell data of temperature.
    boost::shared_ptr<pdat::CellData<double> > data_temperature =
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
                            << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::getCellDataOfDiffusiveFluxVariablesForDerivative()\n"
                            << "There are only x-direction for one-dimensional problem."
                            << std::endl);
                    }
                }
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::getCellDataOfDiffusiveFluxVariablesForDerivative()\n"
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
                            << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::getCellDataOfDiffusiveFluxVariablesForDerivative()\n"
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
                            << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::getCellDataOfDiffusiveFluxVariablesForDerivative()\n"
                            << "There are only x-direction and y-direction for two-dimensional problem."
                            << std::endl);
                    }
                }
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::getCellDataOfDiffusiveFluxVariablesForDerivative()\n"
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
                            << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::getCellDataOfDiffusiveFluxVariablesForDerivative()\n"
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
                            << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::getCellDataOfDiffusiveFluxVariablesForDerivative()\n"
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
                            << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::getCellDataOfDiffusiveFluxVariablesForDerivative()\n"
                            << "There are only x-direction, y-direction and z-direction for three-dimensional problem."
                            << std::endl);
                    }
                }
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::getCellDataOfDiffusiveFluxVariablesForDerivative()\n"
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
    std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > >& diffusivities_data,
    std::vector<std::vector<int> >& diffusivities_component_idx,
    const DIRECTION::TYPE& flux_direction,
    const DIRECTION::TYPE& derivative_direction)
{
    // Create empty box.
    const hier::Box empty_box(d_dim);
    
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    boost::shared_ptr<FlowModel> flow_model_tmp = d_flow_model.lock();
    const hier::Patch& patch = flow_model_tmp->getRegisteredPatch();
    
    diffusivities_data.resize(d_num_eqn);
    diffusivities_component_idx.resize(d_num_eqn);
    
    /*
     * Get the dimension of the interior box.
     */
    
    const hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    if (!d_derived_cell_data_computed)
    {
        // Get the cell data of velocity.
        boost::shared_ptr<pdat::CellData<double> > data_velocity =
            flow_model_tmp->getCellData("VELOCITY");
        
        // Get the cell data of pressure.
        boost::shared_ptr<pdat::CellData<double> > data_pressure =
            flow_model_tmp->getCellData("PRESSURE");
        
        // Get the cell data of temperature.
        boost::shared_ptr<pdat::CellData<double> > data_temperature =
            flow_model_tmp->getCellData("TEMPERATURE");
        
        /*
         * Get the molecular properties of the species for shear viscosity.
         */
        
        std::vector<double> molecular_properties_shear_viscosity;
        std::vector<double*> molecular_properties_shear_viscosity_ptr;
        std::vector<const double*> molecular_properties_shear_viscosity_const_ptr;
        
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
        
        /*
         * Get the molecular properties of the species for bulk viscosity.
         */
        
        std::vector<double> molecular_properties_bulk_viscosity;
        std::vector<double*> molecular_properties_bulk_viscosity_ptr;
        std::vector<const double*> molecular_properties_bulk_viscosity_const_ptr;
        
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
        
        /*
         * Get the molecular properties of the species for thermal conductivity.
         */
        
        std::vector<double> molecular_properties_thermal_conductivity;
        std::vector<double*> molecular_properties_thermal_conductivity_ptr;
        std::vector<const double*> molecular_properties_thermal_conductivity_const_ptr;
        
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
        
        if (!d_cell_data_shear_viscosity_computed)
        {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(d_data_shear_viscosity);
#endif
            
            // Compute the shear viscosity field.
            d_equation_of_shear_viscosity_mixing_rules->getEquationOfShearViscosity()->
                computeShearViscosity(
                    d_data_shear_viscosity,
                    data_pressure,
                    data_temperature,
                    molecular_properties_shear_viscosity_const_ptr,
                    empty_box);
            
            d_cell_data_shear_viscosity_computed = true;
        }
        
        if (!d_cell_data_bulk_viscosity_computed)
        {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(d_data_bulk_viscosity);
#endif
            
            // Compute the bulk viscosity field.
            d_equation_of_bulk_viscosity_mixing_rules->getEquationOfBulkViscosity()->
                computeBulkViscosity(
                    d_data_bulk_viscosity,
                    data_pressure,
                    data_temperature,
                    molecular_properties_bulk_viscosity_const_ptr,
                    empty_box);
            
            d_cell_data_bulk_viscosity_computed = true;
        }
        
        if (!d_cell_data_thermal_conductivity_computed)
        {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(d_data_thermal_conductivity);
#endif
            
            // Compute the thermal conductivity field.
            d_equation_of_thermal_conductivity_mixing_rules->getEquationOfThermalConductivity()->
                computeThermalConductivity(
                    d_data_thermal_conductivity,
                    data_pressure,
                    data_temperature,
                    molecular_properties_thermal_conductivity_const_ptr,
                    empty_box);
            
            d_cell_data_thermal_conductivity_computed = true;
        }
        
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
        
        double* mu    = d_data_shear_viscosity->getPointer(0);
        double* mu_v  = d_data_bulk_viscosity->getPointer(0);
        double* kappa = d_data_thermal_conductivity->getPointer(0);
        
        if (!d_cell_data_diffusivities_computed)
        {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(d_data_diffusivities);
#endif
            
            if (d_dim == tbox::Dimension(1))
            {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(d_data_diffusivities->getDepth() == 3);
#endif
                
                // Get the pointer to cell data of velocity and diffusivities.
                double* u = data_velocity->getPointer(0);
                
                double* D_00 = d_data_diffusivities->getPointer(0);
                double* D_01 = d_data_diffusivities->getPointer(1);
                double* D_02 = d_data_diffusivities->getPointer(2);
                
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
                    
                    D_00[idx_diffusivities] = -(double(4)/double(3)*mu[idx_shear_viscosity] + mu_v[idx_bulk_viscosity]);
                    D_01[idx_diffusivities] = -u[idx_velocity]*(double(4)/double(3)*mu[idx_shear_viscosity] +
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
                double* u = data_velocity->getPointer(0);
                double* v = data_velocity->getPointer(1);
                
                double* D_00 = d_data_diffusivities->getPointer(0);
                double* D_01 = d_data_diffusivities->getPointer(1);
                double* D_02 = d_data_diffusivities->getPointer(2);
                double* D_03 = d_data_diffusivities->getPointer(3);
                double* D_04 = d_data_diffusivities->getPointer(4);
                double* D_05 = d_data_diffusivities->getPointer(5);
                double* D_06 = d_data_diffusivities->getPointer(6);
                double* D_07 = d_data_diffusivities->getPointer(7);
                double* D_08 = d_data_diffusivities->getPointer(8);
                double* D_09 = d_data_diffusivities->getPointer(9);
                
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
                        
                        D_00[idx_diffusivities] = -(double(4)/double(3)*mu[idx_shear_viscosity] + mu_v[idx_bulk_viscosity]);
                        D_01[idx_diffusivities] = double(2)/double(3)*mu[idx_shear_viscosity] - mu_v[idx_bulk_viscosity];
                        D_02[idx_diffusivities] = -mu[idx_shear_viscosity];
                        D_03[idx_diffusivities] = -u[idx_velocity]*(double(4)/double(3)*mu[idx_shear_viscosity] +
                            mu_v[idx_bulk_viscosity]);
                        D_04[idx_diffusivities] = -v[idx_velocity]*(double(4)/double(3)*mu[idx_shear_viscosity] +
                            mu_v[idx_bulk_viscosity]);
                        D_05[idx_diffusivities] = u[idx_velocity]*(double(2)/double(3)*mu[idx_shear_viscosity] -
                            mu_v[idx_bulk_viscosity]);
                        D_06[idx_diffusivities] = v[idx_velocity]*(double(2)/double(3)*mu[idx_shear_viscosity] -
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
                double* u = data_velocity->getPointer(0);
                double* v = data_velocity->getPointer(1);
                double* w = data_velocity->getPointer(2);
                
                double* D_00 = d_data_diffusivities->getPointer(0);
                double* D_01 = d_data_diffusivities->getPointer(1);
                double* D_02 = d_data_diffusivities->getPointer(2);
                double* D_03 = d_data_diffusivities->getPointer(3);
                double* D_04 = d_data_diffusivities->getPointer(4);
                double* D_05 = d_data_diffusivities->getPointer(5);
                double* D_06 = d_data_diffusivities->getPointer(6);
                double* D_07 = d_data_diffusivities->getPointer(7);
                double* D_08 = d_data_diffusivities->getPointer(8);
                double* D_09 = d_data_diffusivities->getPointer(9);
                double* D_10 = d_data_diffusivities->getPointer(10);
                double* D_11 = d_data_diffusivities->getPointer(11);
                double* D_12 = d_data_diffusivities->getPointer(12);
                
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
                            
                            D_00[idx_diffusivities] = -(double(4)/double(3)*mu[idx_shear_viscosity] + mu_v[idx_bulk_viscosity]);
                            D_01[idx_diffusivities] = double(2)/double(3)*mu[idx_shear_viscosity] - mu_v[idx_bulk_viscosity];
                            D_02[idx_diffusivities] = -mu[idx_shear_viscosity];
                            D_03[idx_diffusivities] = -u[idx_velocity]*(double(4)/double(3)*mu[idx_shear_viscosity] +
                                mu_v[idx_bulk_viscosity]);
                            D_04[idx_diffusivities] = -v[idx_velocity]*(double(4)/double(3)*mu[idx_shear_viscosity] +
                                mu_v[idx_bulk_viscosity]);
                            D_05[idx_diffusivities] = -w[idx_velocity]*(double(4)/double(3)*mu[idx_shear_viscosity] +
                                mu_v[idx_bulk_viscosity]);
                            D_06[idx_diffusivities] = u[idx_velocity]*(double(2)/double(3)*mu[idx_shear_viscosity] -
                                mu_v[idx_bulk_viscosity]);
                            D_07[idx_diffusivities] = v[idx_velocity]*(double(2)/double(3)*mu[idx_shear_viscosity] -
                                mu_v[idx_bulk_viscosity]);
                            D_08[idx_diffusivities] = w[idx_velocity]*(double(2)/double(3)*mu[idx_shear_viscosity] -
                                mu_v[idx_bulk_viscosity]);
                            D_09[idx_diffusivities] = -u[idx_velocity]*mu[idx_shear_viscosity];
                            D_10[idx_diffusivities] = -v[idx_velocity]*mu[idx_shear_viscosity];
                            D_11[idx_diffusivities] = -w[idx_velocity]*mu[idx_shear_viscosity];
                            D_12[idx_diffusivities] = -kappa[idx_thermal_conductivity];
                        }
                    }
                }
            }
            
            d_cell_data_diffusivities_computed = true;
        }
        
        d_derived_cell_data_computed = true;
    }
    
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
                            << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::getCellDataOfDiffusiveFluxDiffusivities()\n"
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
                            << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::getCellDataOfDiffusiveFluxDiffusivities()\n"
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
                            << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::getCellDataOfDiffusiveFluxDiffusivities()\n"
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
                            << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::getCellDataOfDiffusiveFluxDiffusivities()\n"
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
                        
                        // 2/3*(mu - mu_v).
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
                        
                        // -w*u.
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
                            << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::getCellDataOfDiffusiveFluxDiffusivities()\n"
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
                            << ": FlowModelDiffusiveFluxUtilitiesSingleSpecies::getCellDataOfDiffusiveFluxDiffusivities()\n"
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
