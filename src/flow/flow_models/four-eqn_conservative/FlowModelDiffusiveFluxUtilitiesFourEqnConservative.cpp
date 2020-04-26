#include "flow/flow_models/four-eqn_conservative/FlowModelDiffusiveFluxUtilitiesFourEqnConservative.hpp"

FlowModelDiffusiveFluxUtilitiesFourEqnConservative::FlowModelDiffusiveFluxUtilitiesFourEqnConservative(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
    const int& num_species,
    const boost::shared_ptr<EquationOfMassDiffusivityMixingRules> equation_of_mass_diffusivity_mixing_rules,
    const boost::shared_ptr<EquationOfShearViscosityMixingRules> equation_of_shear_viscosity_mixing_rules,
    const boost::shared_ptr<EquationOfBulkViscosityMixingRules> equation_of_bulk_viscosity_mixing_rules,
    const boost::shared_ptr<EquationOfThermalConductivityMixingRules> equation_of_thermal_conductivity_mixing_rules):
        FlowModelDiffusiveFluxUtilities(
            object_name,
            dim,
            grid_geometry,
            num_species,
            num_species + dim.getValue() + 1),
        d_num_subghosts_mass_diffusivities(-hier::IntVector::getOne(d_dim)),
        d_num_subghosts_shear_viscosity(-hier::IntVector::getOne(d_dim)),
        d_num_subghosts_bulk_viscosity(-hier::IntVector::getOne(d_dim)),
        d_num_subghosts_thermal_conductivity(-hier::IntVector::getOne(d_dim)),
        d_subghost_box_mass_diffusivities(hier::Box::getEmptyBox(dim)),
        d_subghost_box_shear_viscosity(hier::Box::getEmptyBox(dim)),
        d_subghost_box_bulk_viscosity(hier::Box::getEmptyBox(dim)),
        d_subghost_box_thermal_conductivity(hier::Box::getEmptyBox(dim)),
        d_subghostcell_dims_mass_diffusivities(hier::IntVector::getZero(d_dim)),
        d_subghostcell_dims_shear_viscosity(hier::IntVector::getZero(d_dim)),
        d_subghostcell_dims_bulk_viscosity(hier::IntVector::getZero(d_dim)),
        d_subghostcell_dims_thermal_conductivity(hier::IntVector::getZero(d_dim)),
        d_equation_of_mass_diffusivity_mixing_rules(equation_of_mass_diffusivity_mixing_rules),
        d_equation_of_shear_viscosity_mixing_rules(equation_of_shear_viscosity_mixing_rules),
        d_equation_of_bulk_viscosity_mixing_rules(equation_of_bulk_viscosity_mixing_rules),
        d_equation_of_thermal_conductivity_mixing_rules(equation_of_thermal_conductivity_mixing_rules)
{}


/*
 * Register the required variables for the computation of diffusive fluxes in the registered patch.
 */
void
FlowModelDiffusiveFluxUtilitiesFourEqnConservative::registerDerivedVariablesForDiffusiveFluxes(
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
        std::pair<std::string, hier::IntVector>("DENSITY", num_subghosts));
    
    num_subghosts_of_data.insert(
        std::pair<std::string, hier::IntVector>("MASS_FRACTIONS", num_subghosts));
    
    num_subghosts_of_data.insert(
        std::pair<std::string, hier::IntVector>("VELOCITY", num_subghosts));
    
    num_subghosts_of_data.insert(
        std::pair<std::string, hier::IntVector>("PRESSURE", num_subghosts));
    
    num_subghosts_of_data.insert(
        std::pair<std::string, hier::IntVector>("TEMPERATURE", num_subghosts));
    
    num_subghosts_of_data.insert(
        std::pair<std::string, hier::IntVector>("SPECIES_ENTHALPIES", num_subghosts));
    
    flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
    
    /*
     * Set ghost boxes of derived cell variables for this class.
     */
    
    d_num_subghosts_diffusivities = num_subghosts;
    d_subghost_box_diffusivities = interior_box;
    d_subghost_box_diffusivities.grow(d_num_subghosts_diffusivities);
    d_subghostcell_dims_diffusivities = d_subghost_box_diffusivities.numberCells();
    
    d_num_subghosts_mass_diffusivities   = d_num_subghosts_diffusivities;
    d_num_subghosts_shear_viscosity      = d_num_subghosts_diffusivities;
    d_num_subghosts_bulk_viscosity       = d_num_subghosts_diffusivities;
    d_num_subghosts_thermal_conductivity = d_num_subghosts_diffusivities;
    
    d_subghost_box_mass_diffusivities   = d_subghost_box_diffusivities;
    d_subghost_box_shear_viscosity      = d_subghost_box_diffusivities;
    d_subghost_box_bulk_viscosity       = d_subghost_box_diffusivities;
    d_subghost_box_thermal_conductivity = d_subghost_box_diffusivities;
    
    d_subghostcell_dims_mass_diffusivities   = d_subghostcell_dims_diffusivities;
    d_subghostcell_dims_shear_viscosity      = d_subghostcell_dims_diffusivities;
    d_subghostcell_dims_bulk_viscosity       = d_subghostcell_dims_diffusivities;
    d_subghostcell_dims_thermal_conductivity = d_subghostcell_dims_diffusivities;
}


/*
 * The cell data of all derived variables in the patch for this class are cleared.
 */
void
FlowModelDiffusiveFluxUtilitiesFourEqnConservative::clearCellData()
{
    d_num_subghosts_diffusivities        = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_mass_diffusivities   = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_shear_viscosity      = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_bulk_viscosity       = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_thermal_conductivity = -hier::IntVector::getOne(d_dim);
    
    d_subghost_box_diffusivities        = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_mass_diffusivities   = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_shear_viscosity      = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_bulk_viscosity       = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_thermal_conductivity = hier::Box::getEmptyBox(d_dim);
    
    d_subghostcell_dims_diffusivities        = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_mass_diffusivities   = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_shear_viscosity      = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_bulk_viscosity       = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_thermal_conductivity = hier::IntVector::getZero(d_dim);
    
    d_data_diffusivities.reset();
    d_data_mass_diffusivities.reset();
    d_data_shear_viscosity.reset();
    d_data_bulk_viscosity.reset();
    d_data_thermal_conductivity.reset();
    
    d_derived_cell_data_computed = false;
}


/*
 * Get the variables for the derivatives in the diffusive fluxes.
 */
void
FlowModelDiffusiveFluxUtilitiesFourEqnConservative::getCellDataOfDiffusiveFluxVariablesForDerivative(
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
    
    // Get the cell data of mass fractions.
    boost::shared_ptr<pdat::CellData<double> > data_mass_fractions =
        flow_model_tmp->getCellData("MASS_FRACTIONS");
    
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
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            derivative_var_data[si].resize(d_num_species + 1);
                            derivative_var_component_idx[si].resize(d_num_species + 1);
                            
                            derivative_var_data[si][0] = data_mass_fractions;
                            derivative_var_component_idx[si][0] = si;
                            
                            for (int sj = 0; sj < d_num_species; sj++)
                            {
                                derivative_var_data[si][1 + sj] = data_mass_fractions;
                                derivative_var_component_idx[si][1 + sj] = sj;
                            }
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
                        
                        derivative_var_data[d_num_species + 1].resize(
                            2 + d_num_species*(d_num_species + 1));
                        derivative_var_component_idx[d_num_species + 1].resize(
                            2 + d_num_species*(d_num_species + 1));
                        
                        // Variable u.
                        derivative_var_data[d_num_species + 1][0] = data_velocity;
                        derivative_var_component_idx[d_num_species + 1][0] = 0;
                        
                        // Variable T.
                        derivative_var_data[d_num_species + 1][1] = data_temperature;
                        derivative_var_component_idx[d_num_species + 1][1] = 0;
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            derivative_var_data[d_num_species + 1][2 + si*(d_num_species + 1)] =
                                data_mass_fractions;
                            derivative_var_component_idx[d_num_species + 1][2 + si*(d_num_species + 1)] =
                                si;
                            
                            for (int sj = 0; sj < d_num_species; sj++)
                            {
                                derivative_var_data[d_num_species + 1][3 + si*(d_num_species + 1) + sj] =
                                    data_mass_fractions;
                                derivative_var_component_idx[d_num_species + 1][3 + si*(d_num_species + 1) + sj] =
                                    sj;
                            }
                        }
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::getCellDataOfDiffusiveFluxVariablesForDerivative()\n"
                            << "There are only x-direction for one-dimensional problem."
                            << std::endl);
                    }
                }
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::getCellDataOfDiffusiveFluxVariablesForDerivative()\n"
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
                            derivative_var_data[si].resize(d_num_species + 1);
                            derivative_var_component_idx[si].resize(d_num_species + 1);
                            
                            derivative_var_data[si][0] = data_mass_fractions;
                            derivative_var_component_idx[si][0] = si;
                            
                            for (int sj = 0; sj < d_num_species; sj++)
                            {
                                derivative_var_data[si][1 + sj] = data_mass_fractions;
                                derivative_var_component_idx[si][1 + sj] = sj;
                            }
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
                        
                        derivative_var_data[d_num_species + 2].resize(
                            3 + d_num_species*(d_num_species + 1));
                        derivative_var_component_idx[d_num_species + 2].resize(
                            3 + d_num_species*(d_num_species + 1));
                        
                        // Variable u.
                        derivative_var_data[d_num_species + 2][0] = data_velocity;
                        derivative_var_component_idx[d_num_species + 2][0] = 0;
                        
                        // Variable v.
                        derivative_var_data[d_num_species + 2][1] = data_velocity;
                        derivative_var_component_idx[d_num_species + 2][1] = 1;
                        
                        // Variable T.
                        derivative_var_data[d_num_species + 2][2] = data_temperature;
                        derivative_var_component_idx[d_num_species + 2][2] = 0;
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            derivative_var_data[d_num_species + 2][3 + si*(d_num_species + 1)] =
                                data_mass_fractions;
                            derivative_var_component_idx[d_num_species + 2][3 + si*(d_num_species + 1)] =
                                si;
                            
                            for (int sj = 0; sj < d_num_species; sj++)
                            {
                                derivative_var_data[d_num_species + 2][4 + si*(d_num_species + 1) + sj] =
                                    data_mass_fractions;
                                derivative_var_component_idx[d_num_species + 2][4 + si*(d_num_species + 1) + sj] =
                                    sj;
                            }
                        }
                        
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
                            << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::getCellDataOfDiffusiveFluxVariablesForDerivative()\n"
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
                            derivative_var_data[si].resize(d_num_species + 1);
                            derivative_var_component_idx[si].resize(d_num_species + 1);
                            
                            derivative_var_data[si][0] = data_mass_fractions;
                            derivative_var_component_idx[si][0] = si;
                            
                            for (int sj = 0; sj < d_num_species; sj++)
                            {
                                derivative_var_data[si][1 + sj] = data_mass_fractions;
                                derivative_var_component_idx[si][1 + sj] = sj;
                            }
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
                        
                        derivative_var_data[d_num_species + 2].resize(
                            3 + d_num_species*(d_num_species + 1));
                        derivative_var_component_idx[d_num_species + 2].resize(
                            3 + d_num_species*(d_num_species + 1));
                        
                        // Variable u.
                        derivative_var_data[d_num_species + 2][0] = data_velocity;
                        derivative_var_component_idx[d_num_species + 2][0] = 0;
                        
                        // Variable v.
                        derivative_var_data[d_num_species + 2][1] = data_velocity;
                        derivative_var_component_idx[d_num_species + 2][1] = 1;
                        
                        // Variable T.
                        derivative_var_data[d_num_species + 2][2] = data_temperature;
                        derivative_var_component_idx[d_num_species + 2][2] = 0;
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            derivative_var_data[d_num_species + 2][3 + si*(d_num_species + 1)] =
                                data_mass_fractions;
                            derivative_var_component_idx[d_num_species + 2][3 + si*(d_num_species + 1)] =
                                si;
                            
                            for (int sj = 0; sj < d_num_species; sj++)
                            {
                                derivative_var_data[d_num_species + 2][4 + si*(d_num_species + 1) + sj] =
                                    data_mass_fractions;
                                derivative_var_component_idx[d_num_species + 2][4 + si*(d_num_species + 1) + sj] =
                                    sj;
                            }
                        }
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::getCellDataOfDiffusiveFluxVariablesForDerivative()\n"
                            << "There are only x-direction and y-direction for two-dimensional problem."
                            << std::endl);
                    }
                }
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::getCellDataOfDiffusiveFluxVariablesForDerivative()\n"
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
                            derivative_var_data[si].resize(d_num_species + 1);
                            derivative_var_component_idx[si].resize(d_num_species + 1);
                            
                            derivative_var_data[si][0] = data_mass_fractions;
                            derivative_var_component_idx[si][0] = si;
                            
                            for (int sj = 0; sj < d_num_species; sj++)
                            {
                                derivative_var_data[si][1 + sj] = data_mass_fractions;
                                derivative_var_component_idx[si][1 + sj] = sj;
                            }
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
                        
                        derivative_var_data[d_num_species + 3].resize(
                            4 + d_num_species*(d_num_species + 1));
                        derivative_var_component_idx[d_num_species + 3].resize(
                            4 + d_num_species*(d_num_species + 1));
                        
                        // Variable u.
                        derivative_var_data[d_num_species + 3][0] = data_velocity;
                        derivative_var_component_idx[d_num_species + 3][0] = 0;
                        
                        // Variable v.
                        derivative_var_data[d_num_species + 3][1] = data_velocity;
                        derivative_var_component_idx[d_num_species + 3][1] = 1;
                        
                        // Variable w.
                        derivative_var_data[d_num_species + 3][2] = data_velocity;
                        derivative_var_component_idx[d_num_species + 3][2] = 2;
                        
                        // Variable T.
                        derivative_var_data[d_num_species + 3][3] = data_temperature;
                        derivative_var_component_idx[d_num_species + 3][3] = 0;
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            derivative_var_data[d_num_species + 3][4 + si*(d_num_species + 1)] =
                                data_mass_fractions;
                            derivative_var_component_idx[d_num_species + 3][4 + si*(d_num_species + 1)] =
                                si;
                            
                            for (int sj = 0; sj < d_num_species; sj++)
                            {
                                derivative_var_data[d_num_species + 3][5 + si*(d_num_species + 1) + sj] =
                                    data_mass_fractions;
                                derivative_var_component_idx[d_num_species + 3][5 + si*(d_num_species + 1) + sj] =
                                    sj;
                            }
                        }
                        
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
                            << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::getCellDataOfDiffusiveFluxVariablesForDerivative()\n"
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
                            derivative_var_data[si].resize(d_num_species + 1);
                            derivative_var_component_idx[si].resize(d_num_species + 1);
                            
                            derivative_var_data[si][0] = data_mass_fractions;
                            derivative_var_component_idx[si][0] = si;
                            
                            for (int sj = 0; sj < d_num_species; sj++)
                            {
                                derivative_var_data[si][1 + sj] = data_mass_fractions;
                                derivative_var_component_idx[si][1 + sj] = sj;
                            }
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
                        
                        derivative_var_data[d_num_species + 3].resize(
                            4 + d_num_species*(d_num_species + 1));
                        derivative_var_component_idx[d_num_species + 3].resize(
                            4 + d_num_species*(d_num_species + 1));
                        
                        // Variable u.
                        derivative_var_data[d_num_species + 3][0] = data_velocity;
                        derivative_var_component_idx[d_num_species + 3][0] = 0;
                        
                        // Variable v.
                        derivative_var_data[d_num_species + 3][1] = data_velocity;
                        derivative_var_component_idx[d_num_species + 3][1] = 1;
                        
                        // Variable w.
                        derivative_var_data[d_num_species + 3][2] = data_velocity;
                        derivative_var_component_idx[d_num_species + 3][2] = 2;
                        
                        // Variable T.
                        derivative_var_data[d_num_species + 3][3] = data_temperature;
                        derivative_var_component_idx[d_num_species + 3][3] = 0;
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            derivative_var_data[d_num_species + 3][4 + si*(d_num_species + 1)] =
                                data_mass_fractions;
                            derivative_var_component_idx[d_num_species + 3][4 + si*(d_num_species + 1)] =
                                si;
                            
                            for (int sj = 0; sj < d_num_species; sj++)
                            {
                                derivative_var_data[d_num_species + 3][5 + si*(d_num_species + 1) + sj] =
                                    data_mass_fractions;
                                derivative_var_component_idx[d_num_species + 3][5 + si*(d_num_species + 1) + sj] =
                                    sj;
                            }
                        }
                        
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
                            << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::getCellDataOfDiffusiveFluxVariablesForDerivative()\n"
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
                            derivative_var_data[si].resize(d_num_species + 1);
                            derivative_var_component_idx[si].resize(d_num_species + 1);
                            
                            derivative_var_data[si][0] = data_mass_fractions;
                            derivative_var_component_idx[si][0] = si;
                            
                            for (int sj = 0; sj < d_num_species; sj++)
                            {
                                derivative_var_data[si][1 + sj] = data_mass_fractions;
                                derivative_var_component_idx[si][1 + sj] = sj;
                            }
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
                        
                        derivative_var_data[d_num_species + 3].resize(
                            4 + d_num_species*(d_num_species + 1));
                        derivative_var_component_idx[d_num_species + 3].resize(
                            4 + d_num_species*(d_num_species + 1));
                        
                        // Variable u.
                        derivative_var_data[d_num_species + 3][0] = data_velocity;
                        derivative_var_component_idx[d_num_species + 3][0] = 0;
                        
                        // Variable v.
                        derivative_var_data[d_num_species + 3][1] = data_velocity;
                        derivative_var_component_idx[d_num_species + 3][1] = 1;
                        
                        // Variable w.
                        derivative_var_data[d_num_species + 3][2] = data_velocity;
                        derivative_var_component_idx[d_num_species + 3][2] = 2;
                        
                        // Variable T.
                        derivative_var_data[d_num_species + 3][3] = data_temperature;
                        derivative_var_component_idx[d_num_species + 3][3] = 0;
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            derivative_var_data[d_num_species + 3][4 + si*(d_num_species + 1)] =
                                data_mass_fractions;
                            derivative_var_component_idx[d_num_species + 3][4 + si*(d_num_species + 1)] =
                                si;
                            
                            for (int sj = 0; sj < d_num_species; sj++)
                            {
                                derivative_var_data[d_num_species + 3][5 + si*(d_num_species + 1) + sj] =
                                    data_mass_fractions;
                                derivative_var_component_idx[d_num_species + 3][5 + si*(d_num_species + 1) + sj] =
                                    sj;
                            }
                        }
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::getCellDataOfDiffusiveFluxVariablesForDerivative()\n"
                            << "There are only x-direction, y-direction and z-direction for three-dimensional problem."
                            << std::endl);
                    }
                }
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::getCellDataOfDiffusiveFluxVariablesForDerivative()\n"
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
FlowModelDiffusiveFluxUtilitiesFourEqnConservative::getCellDataOfDiffusiveFluxDiffusivities(
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
        // Get the cell data of density.
        boost::shared_ptr<pdat::CellData<double> > data_density =
            flow_model_tmp->getCellData("DENSITY");
        
        // Get the cell data of mass fractions.
        boost::shared_ptr<pdat::CellData<double> > data_mass_fractions =
            flow_model_tmp->getCellData("MASS_FRACTIONS");
        
        // Get the cell data of velocity.
        boost::shared_ptr<pdat::CellData<double> > data_velocity =
            flow_model_tmp->getCellData("VELOCITY");
        
        // Get the cell data of pressure.
        boost::shared_ptr<pdat::CellData<double> > data_pressure =
            flow_model_tmp->getCellData("PRESSURE");
        
        // Get the cell data of temperature.
        boost::shared_ptr<pdat::CellData<double> > data_temperature =
            flow_model_tmp->getCellData("TEMPERATURE");
        
        // Get the cell data of species enthalpies.
        std::vector<boost::shared_ptr<pdat::CellData<double> > > data_species_enthalpies =
            flow_model_tmp->getSpeciesCellData("SPECIES_ENTHALPIES");
        
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
        for (int si = 1; si < d_num_species; si++)
        {
            TBOX_ASSERT(data_species_enthalpies[si]->getBox().isSpatiallyEqual(data_species_enthalpies[0]->getBox()));
            TBOX_ASSERT(data_species_enthalpies[si]->getGhostCellWidth() ==
                data_species_enthalpies[0]->getGhostCellWidth());
        }
#endif
        
        /*
         * Get the numbers of ghost cells.
         */
        
        const hier::IntVector num_subghosts_density = data_density->getGhostCellWidth();
        const hier::IntVector num_subghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
        const hier::IntVector num_subghosts_velocity = data_velocity->getGhostCellWidth();
        const hier::IntVector num_subghosts_species_enthalpies = data_species_enthalpies[0]->getGhostCellWidth();
        
        /*
         * Get the dimensions of the ghost cell boxes.
         */
        
        const hier::Box subghost_box_mass_fractions = data_mass_fractions->getGhostBox();
        const hier::IntVector subghostcell_dims_mass_fractions = subghost_box_mass_fractions.numberCells();
        
        const hier::Box subghost_box_density = data_density->getGhostBox();
        const hier::IntVector subghostcell_dims_density = subghost_box_density.numberCells();
        
        const hier::Box subghost_box_velocity = data_velocity->getGhostBox();
        const hier::IntVector subghostcell_dims_velocity = subghost_box_velocity.numberCells();
        
        const hier::Box subghost_box_species_enthalpies = data_species_enthalpies[0]->getGhostBox();
        const hier::IntVector subghostcell_dims_species_enthalpies = subghost_box_species_enthalpies.numberCells();
        
        /*
         * Create cell data of mass diffusivities, shear viscosity, bulk viscosity and thermal conductivity.
         */
        
        d_data_mass_diffusivities.reset(new pdat::CellData<double>(
            interior_box, d_num_species, d_num_subghosts_mass_diffusivities));
        
        d_data_shear_viscosity.reset(new pdat::CellData<double>(
            interior_box, 1, d_num_subghosts_shear_viscosity));
        
        d_data_bulk_viscosity.reset(new pdat::CellData<double>(
            interior_box, 1, d_num_subghosts_bulk_viscosity));
        
        d_data_thermal_conductivity.reset(new pdat::CellData<double>(
            interior_box, 1, d_num_subghosts_thermal_conductivity));
        
        /*
         * Get the pointers to the cell data of density, mass fractions, mass diffusivities, shear
         * viscosity, bulk viscosity and thermal conductivity.
         */
        
        double* rho = data_density->getPointer(0);
        std::vector<double*> Y;
        Y.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            Y.push_back(data_mass_fractions->getPointer(si));
        }
        std::vector<double*> D;
        D.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            D.push_back(d_data_mass_diffusivities->getPointer(si));
        }
        std::vector<double*> h_i;
        h_i.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            h_i.push_back(data_species_enthalpies[si]->getPointer(0));
        }
        double* mu    = d_data_shear_viscosity->getPointer(0);
        double* mu_v  = d_data_bulk_viscosity->getPointer(0);
        double* kappa = d_data_thermal_conductivity->getPointer(0);
        
        // Compute the mass diffusivity fields.
        d_equation_of_mass_diffusivity_mixing_rules->computeMassDiffusivities(
            d_data_mass_diffusivities,
            data_pressure,
            data_temperature,
            data_mass_fractions,
            empty_box);
        
        // Compute the shear viscosity field.
        d_equation_of_shear_viscosity_mixing_rules->computeShearViscosity(
            d_data_shear_viscosity,
            data_pressure,
            data_temperature,
            data_mass_fractions,
            empty_box);
        
        // Compute the bulk viscosity field.
        d_equation_of_bulk_viscosity_mixing_rules->computeBulkViscosity(
            d_data_bulk_viscosity,
            data_pressure,
            data_temperature,
            data_mass_fractions,
            empty_box);
        
        // Compute the thermal conductivity field.
        d_equation_of_thermal_conductivity_mixing_rules->computeThermalConductivity(
            d_data_thermal_conductivity,
            data_pressure,
            data_temperature,
            data_mass_fractions,
            empty_box);
        
        if (d_dim == tbox::Dimension(1))
        {
            d_data_diffusivities.reset(new pdat::CellData<double>(
                interior_box,
                2*d_num_species*(d_num_species + 1) + 3,
                d_num_subghosts_diffusivities));
            
            // Get the pointer to cell data of velocity and diffusivities.
            double* u = data_velocity->getPointer(0);
            
            std::vector<double*> D_ptr;
            D_ptr.reserve(2*d_num_species*(d_num_species + 1) + 3);
            
            for (int i = 0; i < 2*d_num_species*(d_num_species + 1) + 3; i++)
            {
                D_ptr.push_back(d_data_diffusivities->getPointer(i));
            }
            
            /*
             * Compute the diffusivities.
             */
            
            // Mass equations.
            for (int si = 0; si < d_num_species; si++)
            {
                const int component_idx = si*(d_num_species + 1);
                
                for (int i = -d_num_subghosts_diffusivities[0];
                     i < interior_dims[0] + d_num_subghosts_diffusivities[0];
                     i++)
                {
                    // Compute the linear indices.
                    const int idx_diffusivities = i + d_num_subghosts_diffusivities[0];
                    const int idx_mass_diffusivities = i + d_num_subghosts_mass_diffusivities[0];
                    const int idx_density = i + num_subghosts_density[0];
                    
                    D_ptr[component_idx][idx_diffusivities] =
                        -rho[idx_density]*D[si][idx_mass_diffusivities];
                }
            }
            
            for (int si = 0; si < d_num_species; si++)
            {
                for (int sj = 0; sj < d_num_species; sj++)
                {
                    const int component_idx = si*(d_num_species + 1) + sj + 1;
                    
                    for (int i = -d_num_subghosts_diffusivities[0];
                         i < interior_dims[0] + d_num_subghosts_diffusivities[0];
                         i++)
                    {
                        // Compute the linear indices.
                        const int idx_diffusivities = i + d_num_subghosts_diffusivities[0];
                        const int idx_mass_diffusivities = i + d_num_subghosts_mass_diffusivities[0];
                        const int idx_density = i + num_subghosts_density[0];
                        const int idx_mass_fractions = i + num_subghosts_mass_fractions[0];
                        
                        D_ptr[component_idx][idx_diffusivities] =
                            rho[idx_density]*Y[si][idx_mass_fractions]*D[sj][idx_mass_diffusivities];
                    }
                }
            }
            
            // Momentum and energy equations.
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
                
                D_ptr[d_num_species*(d_num_species + 1)][idx_diffusivities] =
                    -(double(4)/double(3)*mu[idx_shear_viscosity] + mu_v[idx_bulk_viscosity]);
                D_ptr[d_num_species*(d_num_species + 1) + 1][idx_diffusivities] =
                    -u[idx_velocity]*(double(4)/double(3)*mu[idx_shear_viscosity] + mu_v[idx_bulk_viscosity]);
                D_ptr[d_num_species*(d_num_species + 1) + 2][idx_diffusivities] =
                    -kappa[idx_thermal_conductivity];
            }
            
            // Energy equation.
            for (int si = 0; si < d_num_species; si++)
            {
                const int component_idx =
                    d_num_species*(d_num_species + 1) + 3 + si*(d_num_species + 1);
                
                for (int i = -d_num_subghosts_diffusivities[0];
                     i < interior_dims[0] + d_num_subghosts_diffusivities[0];
                     i++)
                {
                    // Compute the linear indices.
                    const int idx_diffusivities = i + d_num_subghosts_diffusivities[0];
                    const int idx_mass_diffusivities = i + d_num_subghosts_mass_diffusivities[0];
                    const int idx_density = i + num_subghosts_density[0];
                    const int idx_species_enthalpies = i + num_subghosts_species_enthalpies[0];
                    
                    D_ptr[component_idx][idx_diffusivities] =
                        -rho[idx_density]*D[si][idx_mass_diffusivities]*
                            h_i[si][idx_species_enthalpies];
                }
            }
            
            for (int si = 0; si < d_num_species; si++)
            {
                for (int sj = 0; sj < d_num_species; sj++)
                {
                    const int component_idx =
                        d_num_species*(d_num_species + 1) + 4 + si*(d_num_species + 1) + sj;
                    
                    for (int i = -d_num_subghosts_diffusivities[0];
                         i < interior_dims[0] + d_num_subghosts_diffusivities[0];
                         i++)
                    {
                        // Compute the linear indices.
                        const int idx_diffusivities = i + d_num_subghosts_diffusivities[0];
                        const int idx_mass_diffusivities = i + d_num_subghosts_mass_diffusivities[0];
                        const int idx_density = i + num_subghosts_density[0];
                        const int idx_mass_fractions = i + num_subghosts_mass_fractions[0];
                        const int idx_species_enthalpies = i + num_subghosts_species_enthalpies[0];
                        
                        D_ptr[component_idx][idx_diffusivities] =
                            rho[idx_density]*Y[si][idx_mass_fractions]*D[sj][idx_mass_diffusivities]*
                                h_i[si][idx_species_enthalpies];
                    }
                }
            }
        }
        else if (d_dim == tbox::Dimension(2))
        {
            d_data_diffusivities.reset(new pdat::CellData<double>(
                interior_box,
                2*d_num_species*(d_num_species + 1) + 10,
                d_num_subghosts_diffusivities));
            
            // Get the pointer to cell data of velocity and diffusivities.
            double* u = data_velocity->getPointer(0);
            double* v = data_velocity->getPointer(1);
            
            std::vector<double*> D_ptr;
            D_ptr.reserve(2*d_num_species*(d_num_species + 1) + 10);
            
            for (int i = 0; i < 2*d_num_species*(d_num_species + 1) + 10; i++)
            {
                D_ptr.push_back(d_data_diffusivities->getPointer(i));
            }
            
            /*
             * Compute the diffusivities.
             */
            
            // Mass equations.
            for (int si = 0; si < d_num_species; si++)
            {
                const int component_idx = si*(d_num_species + 1);
                
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
                        
                        const int idx_mass_diffusivities = (i + d_num_subghosts_mass_diffusivities[0]) +
                            (j + d_num_subghosts_mass_diffusivities[1])*d_subghostcell_dims_mass_diffusivities[0];
                        
                        const int idx_density = (i + num_subghosts_density[0]) +
                            (j + num_subghosts_density[1])*subghostcell_dims_density[0];
                        
                        D_ptr[component_idx][idx_diffusivities] =
                            -rho[idx_density]*D[si][idx_mass_diffusivities];
                    }
                }
            }
            
            for (int si = 0; si < d_num_species; si++)
            {
                for (int sj = 0; sj < d_num_species; sj++)
                {
                    const int component_idx = si*(d_num_species + 1) + sj + 1;
                    
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
                            
                            const int idx_mass_diffusivities = (i + d_num_subghosts_mass_diffusivities[0]) +
                                (j + d_num_subghosts_mass_diffusivities[1])*d_subghostcell_dims_mass_diffusivities[0];
                            
                            const int idx_density = (i + num_subghosts_density[0]) +
                                (j + num_subghosts_density[1])*subghostcell_dims_density[0];
                            
                            const int idx_mass_fractions = (i + num_subghosts_mass_fractions[0]) +
                                (j + num_subghosts_mass_fractions[1])*subghostcell_dims_mass_fractions[0];
                            
                            D_ptr[component_idx][idx_diffusivities] =
                                rho[idx_density]*Y[si][idx_mass_fractions]*D[sj][idx_mass_diffusivities];
                        }
                    }
                }
            }
            
            // Momentum and energy equations.
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
                    
                    D_ptr[d_num_species*(d_num_species + 1)][idx_diffusivities] =
                        -(double(4)/double(3)*mu[idx_shear_viscosity] + mu_v[idx_bulk_viscosity]);
                    D_ptr[d_num_species*(d_num_species + 1) + 1][idx_diffusivities] =
                        double(2)/double(3)*mu[idx_shear_viscosity] - mu_v[idx_bulk_viscosity];
                    D_ptr[d_num_species*(d_num_species + 1) + 2][idx_diffusivities] =
                        -mu[idx_shear_viscosity];
                    D_ptr[d_num_species*(d_num_species + 1) + 3][idx_diffusivities] =
                        -u[idx_velocity]*(double(4)/double(3)*mu[idx_shear_viscosity] + mu_v[idx_bulk_viscosity]);
                    D_ptr[d_num_species*(d_num_species + 1) + 4][idx_diffusivities] =
                        -v[idx_velocity]*(double(4)/double(3)*mu[idx_shear_viscosity] + mu_v[idx_bulk_viscosity]);
                    D_ptr[d_num_species*(d_num_species + 1) + 5][idx_diffusivities] =
                        u[idx_velocity]*(double(2)/double(3)*mu[idx_shear_viscosity] - mu_v[idx_bulk_viscosity]);
                    D_ptr[d_num_species*(d_num_species + 1) + 6][idx_diffusivities] =
                        v[idx_velocity]*(double(2)/double(3)*mu[idx_shear_viscosity] - mu_v[idx_bulk_viscosity]);
                    D_ptr[d_num_species*(d_num_species + 1) + 7][idx_diffusivities] =
                        -u[idx_velocity]*mu[idx_shear_viscosity];
                    D_ptr[d_num_species*(d_num_species + 1) + 8][idx_diffusivities] =
                        -v[idx_velocity]*mu[idx_shear_viscosity];
                    D_ptr[d_num_species*(d_num_species + 1) + 9][idx_diffusivities] =
                        -kappa[idx_thermal_conductivity];
                }
            }
            
            // Energy equation.
            for (int si = 0; si < d_num_species; si++)
            {
                const int component_idx =
                    d_num_species*(d_num_species + 1) + 10 + si*(d_num_species + 1);
                
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
                        
                        const int idx_mass_diffusivities = (i + d_num_subghosts_mass_diffusivities[0]) +
                            (j + d_num_subghosts_mass_diffusivities[1])*d_subghostcell_dims_mass_diffusivities[0];
                        
                        const int idx_density = (i + num_subghosts_density[0]) +
                            (j + num_subghosts_density[1])*subghostcell_dims_density[0];
                        
                        const int idx_species_enthalpies = (i + num_subghosts_species_enthalpies[0]) +
                            (j + num_subghosts_species_enthalpies[1])*subghostcell_dims_species_enthalpies[0];
                        
                        D_ptr[component_idx][idx_diffusivities] =
                            -rho[idx_density]*D[si][idx_mass_diffusivities]*
                                h_i[si][idx_species_enthalpies];
                    }
                }
            }
            
            for (int si = 0; si < d_num_species; si++)
            {
                for (int sj = 0; sj < d_num_species; sj++)
                {
                    const int component_idx =
                        d_num_species*(d_num_species + 1) + 11 + si*(d_num_species + 1) + sj;
                    
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
                            
                            const int idx_mass_diffusivities = (i + d_num_subghosts_mass_diffusivities[0]) +
                                (j + d_num_subghosts_mass_diffusivities[1])*d_subghostcell_dims_mass_diffusivities[0];
                            
                            const int idx_density = (i + num_subghosts_density[0]) +
                                (j + num_subghosts_density[1])*subghostcell_dims_density[0];
                            
                            const int idx_mass_fractions = (i + num_subghosts_mass_fractions[0]) +
                                (j + num_subghosts_mass_fractions[1])*subghostcell_dims_mass_fractions[0];
                            
                            const int idx_species_enthalpies = (i + num_subghosts_species_enthalpies[0]) +
                                (j + num_subghosts_species_enthalpies[1])*subghostcell_dims_species_enthalpies[0];
                            
                            D_ptr[component_idx][idx_diffusivities] =
                                rho[idx_density]*Y[si][idx_mass_fractions]*D[sj][idx_mass_diffusivities]*
                                    h_i[si][idx_species_enthalpies];
                        }
                    }
                }
            }
        }
        else if (d_dim == tbox::Dimension(3))
        {
            d_data_diffusivities.reset(new pdat::CellData<double>(
                interior_box,
                2*d_num_species*(d_num_species + 1) + 13,
                d_num_subghosts_diffusivities));
            
            // Get the pointer to cell data of velocity and diffusivities.
            double* u = data_velocity->getPointer(0);
            double* v = data_velocity->getPointer(1);
            double* w = data_velocity->getPointer(2);
            
            std::vector<double*> D_ptr;
            D_ptr.reserve(2*d_num_species*(d_num_species + 1) + 13);
            
            for (int i = 0; i < 2*d_num_species*(d_num_species + 1) + 13; i++)
            {
                D_ptr.push_back(d_data_diffusivities->getPointer(i));
            }
            
            /*
             * Compute the diffusivities.
             */
            
            // Mass equations.
            for (int si = 0; si < d_num_species; si++)
            {
                const int component_idx = si*(d_num_species + 1);
                
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
                            // Compute the linear indices.
                            const int idx_diffusivities = (i + d_num_subghosts_diffusivities[0]) +
                                (j + d_num_subghosts_diffusivities[1])*d_subghostcell_dims_diffusivities[0] +
                                (k + d_num_subghosts_diffusivities[2])*d_subghostcell_dims_diffusivities[0]*
                                    d_subghostcell_dims_diffusivities[1];
                            
                            const int idx_mass_diffusivities = (i + d_num_subghosts_mass_diffusivities[0]) +
                                (j + d_num_subghosts_mass_diffusivities[1])*d_subghostcell_dims_mass_diffusivities[0] +
                                (k + d_num_subghosts_mass_diffusivities[2])*d_subghostcell_dims_mass_diffusivities[0]*
                                    d_subghostcell_dims_mass_diffusivities[1];
                            
                            const int idx_density = (i + num_subghosts_density[0]) +
                                (j + num_subghosts_density[1])*subghostcell_dims_density[0] +
                                (k + num_subghosts_density[2])*subghostcell_dims_density[0]*
                                    subghostcell_dims_density[1];
                            
                            D_ptr[component_idx][idx_diffusivities] =
                                -rho[idx_density]*D[si][idx_mass_diffusivities];
                        }
                    }
                }
            }
            
            for (int si = 0; si < d_num_species; si++)
            {
                for (int sj = 0; sj < d_num_species; sj++)
                {
                    const int component_idx = si*(d_num_species + 1) + sj + 1;
                    
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
                                // Compute the linear indices.
                                const int idx_diffusivities = (i + d_num_subghosts_diffusivities[0]) +
                                    (j + d_num_subghosts_diffusivities[1])*d_subghostcell_dims_diffusivities[0] +
                                    (k + d_num_subghosts_diffusivities[2])*d_subghostcell_dims_diffusivities[0]*
                                        d_subghostcell_dims_diffusivities[1];
                                
                                const int idx_mass_diffusivities = (i + d_num_subghosts_mass_diffusivities[0]) +
                                    (j + d_num_subghosts_mass_diffusivities[1])*d_subghostcell_dims_mass_diffusivities[0] +
                                    (k + d_num_subghosts_mass_diffusivities[2])*d_subghostcell_dims_mass_diffusivities[0]*
                                        d_subghostcell_dims_mass_diffusivities[1];
                                
                                const int idx_density = (i + num_subghosts_density[0]) +
                                    (j + num_subghosts_density[1])*subghostcell_dims_density[0] +
                                    (k + num_subghosts_density[2])*subghostcell_dims_density[0]*
                                        subghostcell_dims_density[1];
                                
                                const int idx_mass_fractions = (i + num_subghosts_mass_fractions[0]) +
                                    (j + num_subghosts_mass_fractions[1])*subghostcell_dims_mass_fractions[0] +
                                    (k + num_subghosts_mass_fractions[2])*subghostcell_dims_mass_fractions[0]*
                                        subghostcell_dims_mass_fractions[1];
                                
                                D_ptr[component_idx][idx_diffusivities] =
                                    rho[idx_density]*Y[si][idx_mass_fractions]*D[sj][idx_mass_diffusivities];
                            }
                        }
                    }
                }
            }
            
            // Momentum and energy equations.
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
                        // Compute the linear indices.
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
                        
                        D_ptr[d_num_species*(d_num_species + 1)][idx_diffusivities] =
                            -(double(4)/double(3)*mu[idx_shear_viscosity] + mu_v[idx_bulk_viscosity]);
                        D_ptr[d_num_species*(d_num_species + 1) + 1][idx_diffusivities] =
                            double(2)/double(3)*mu[idx_shear_viscosity] - mu_v[idx_bulk_viscosity];
                        D_ptr[d_num_species*(d_num_species + 1) + 2][idx_diffusivities] =
                            -mu[idx_shear_viscosity];
                        D_ptr[d_num_species*(d_num_species + 1) + 3][idx_diffusivities] =
                            -u[idx_velocity]*(double(4)/double(3)*mu[idx_shear_viscosity] + mu_v[idx_bulk_viscosity]);
                        D_ptr[d_num_species*(d_num_species + 1) + 4][idx_diffusivities] =
                            -v[idx_velocity]*(double(4)/double(3)*mu[idx_shear_viscosity] + mu_v[idx_bulk_viscosity]);
                        D_ptr[d_num_species*(d_num_species + 1) + 5][idx_diffusivities] =
                            -w[idx_velocity]*(double(4)/double(3)*mu[idx_shear_viscosity] + mu_v[idx_bulk_viscosity]);
                        D_ptr[d_num_species*(d_num_species + 1) + 6][idx_diffusivities] =
                            u[idx_velocity]*(double(2)/double(3)*mu[idx_shear_viscosity] - mu_v[idx_bulk_viscosity]);
                        D_ptr[d_num_species*(d_num_species + 1) + 7][idx_diffusivities] =
                            v[idx_velocity]*(double(2)/double(3)*mu[idx_shear_viscosity] - mu_v[idx_bulk_viscosity]);
                        D_ptr[d_num_species*(d_num_species + 1) + 8][idx_diffusivities] =
                            w[idx_velocity]*(double(2)/double(3)*mu[idx_shear_viscosity] - mu_v[idx_bulk_viscosity]);
                        D_ptr[d_num_species*(d_num_species + 1) + 9][idx_diffusivities] =
                            -u[idx_velocity]*mu[idx_shear_viscosity];
                        D_ptr[d_num_species*(d_num_species + 1) + 10][idx_diffusivities] =
                            -v[idx_velocity]*mu[idx_shear_viscosity];
                        D_ptr[d_num_species*(d_num_species + 1) + 11][idx_diffusivities] =
                            -w[idx_velocity]*mu[idx_shear_viscosity];
                        D_ptr[d_num_species*(d_num_species + 1) + 12][idx_diffusivities] =
                            -kappa[idx_thermal_conductivity];
                    }
                }
            }
            
            // Energy equation.
            for (int si = 0; si < d_num_species; si++)
            {
                const int component_idx =
                    d_num_species*(d_num_species + 1) + 13 + si*(d_num_species + 1);
                
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
                            // Compute the linear indices.
                            const int idx_diffusivities = (i + d_num_subghosts_diffusivities[0]) +
                                (j + d_num_subghosts_diffusivities[1])*d_subghostcell_dims_diffusivities[0] +
                                (k + d_num_subghosts_diffusivities[2])*d_subghostcell_dims_diffusivities[0]*
                                    d_subghostcell_dims_diffusivities[1];
                            
                            const int idx_mass_diffusivities = (i + d_num_subghosts_mass_diffusivities[0]) +
                                (j + d_num_subghosts_mass_diffusivities[1])*d_subghostcell_dims_mass_diffusivities[0] +
                                (k + d_num_subghosts_mass_diffusivities[2])*d_subghostcell_dims_mass_diffusivities[0]*
                                    d_subghostcell_dims_mass_diffusivities[1];
                            
                            const int idx_density = (i + num_subghosts_density[0]) +
                                (j + num_subghosts_density[1])*subghostcell_dims_density[0] +
                                (k + num_subghosts_density[2])*subghostcell_dims_density[0]*
                                    subghostcell_dims_density[1];
                            
                            const int idx_species_enthalpies = (i + num_subghosts_species_enthalpies[0]) +
                                (j + num_subghosts_species_enthalpies[1])*subghostcell_dims_species_enthalpies[0] +
                                (k + num_subghosts_species_enthalpies[2])*subghostcell_dims_species_enthalpies[0]*
                                    subghostcell_dims_species_enthalpies[1];
                            
                            D_ptr[component_idx][idx_diffusivities] =
                                -rho[idx_density]*D[si][idx_mass_diffusivities]*
                                    h_i[si][idx_species_enthalpies];
                        }
                    }
                }
            }
            
            for (int si = 0; si < d_num_species; si++)
            {
                for (int sj = 0; sj < d_num_species; sj++)
                {
                    const int component_idx =
                        d_num_species*(d_num_species + 1) + 14 + si*(d_num_species + 1) + sj;
                    
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
                                // Compute the linear indices.
                                const int idx_diffusivities = (i + d_num_subghosts_diffusivities[0]) +
                                    (j + d_num_subghosts_diffusivities[1])*d_subghostcell_dims_diffusivities[0] +
                                    (k + d_num_subghosts_diffusivities[2])*d_subghostcell_dims_diffusivities[0]*
                                        d_subghostcell_dims_diffusivities[1];
                                
                                const int idx_mass_diffusivities = (i + d_num_subghosts_mass_diffusivities[0]) +
                                    (j + d_num_subghosts_mass_diffusivities[1])*d_subghostcell_dims_mass_diffusivities[0] +
                                    (k + d_num_subghosts_mass_diffusivities[2])*d_subghostcell_dims_mass_diffusivities[0]*
                                        d_subghostcell_dims_mass_diffusivities[1];
                                
                                const int idx_density = (i + num_subghosts_density[0]) +
                                    (j + num_subghosts_density[1])*subghostcell_dims_density[0] +
                                    (k + num_subghosts_density[2])*subghostcell_dims_density[0]*
                                        subghostcell_dims_density[1];
                                
                                const int idx_mass_fractions = (i + num_subghosts_mass_fractions[0]) +
                                    (j + num_subghosts_mass_fractions[1])*subghostcell_dims_mass_fractions[0] +
                                    (k + num_subghosts_mass_fractions[2])*subghostcell_dims_mass_fractions[0]*
                                        subghostcell_dims_mass_fractions[1];
                                
                                const int idx_species_enthalpies = (i + num_subghosts_species_enthalpies[0]) +
                                    (j + num_subghosts_species_enthalpies[1])*subghostcell_dims_species_enthalpies[0] +
                                    (k + num_subghosts_species_enthalpies[2])*subghostcell_dims_species_enthalpies[0]*
                                        subghostcell_dims_species_enthalpies[1];
                                
                                D_ptr[component_idx][idx_diffusivities] =
                                    rho[idx_density]*Y[si][idx_mass_fractions]*D[sj][idx_mass_diffusivities]*
                                        h_i[si][idx_species_enthalpies];
                            }
                        }
                    }
                }
            }
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
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            diffusivities_data[si].resize(d_num_species + 1);
                            diffusivities_component_idx[si].resize(d_num_species + 1);
                            
                            diffusivities_data[si][0] = d_data_diffusivities;
                            diffusivities_component_idx[si][0] = si*(d_num_species + 1);
                            
                            for (int sj = 0; sj < d_num_species; sj++)
                            {
                                diffusivities_data[si][1 + sj] = d_data_diffusivities;
                                diffusivities_component_idx[si][1 + sj] = si*(d_num_species + 1) + sj + 1;
                            }
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[d_num_species].resize(1);
                        diffusivities_component_idx[d_num_species].resize(1);
                        
                        // -(4/3*mu + mu_v).
                        diffusivities_data[d_num_species][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species][0] =
                            d_num_species*(d_num_species + 1);
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[d_num_species + 1].resize(
                            2 + d_num_species*(d_num_species + 1));
                        diffusivities_component_idx[d_num_species + 1].resize(
                            2 + d_num_species*(d_num_species + 1));
                        
                        // -u*(4/3*mu + mu_v).
                        diffusivities_data[d_num_species + 1][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 1][0] =
                            d_num_species*(d_num_species + 1) + 1;
                        
                        // -kappa.
                        diffusivities_data[d_num_species + 1][1] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 1][1] =
                            d_num_species*(d_num_species + 1) + 2;
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            diffusivities_data[d_num_species + 1][2 + si*(d_num_species + 1)] =
                                d_data_diffusivities;
                            diffusivities_component_idx[d_num_species + 1][2 + si*(d_num_species + 1)] =
                                d_num_species*(d_num_species + 1) + 3 + si*(d_num_species + 1);
                            
                            for (int sj = 0; sj < d_num_species; sj++)
                            {
                                diffusivities_data[d_num_species + 1][3 + si*(d_num_species + 1) + sj] =
                                    d_data_diffusivities;
                                diffusivities_component_idx[d_num_species + 1][3 + si*(d_num_species + 1) + sj] =
                                    d_num_species*(d_num_species + 1) + 4 + si*(d_num_species + 1) + sj;
                            }
                        }
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::getCellDataOfDiffusiveFluxDiffusivities()\n"
                            << "There are only x-direction for one-dimensional problem."
                            << std::endl);
                    }
                }
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::getCellDataOfDiffusiveFluxDiffusivities()\n"
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
                            diffusivities_data[si].resize(d_num_species + 1);
                            diffusivities_component_idx[si].resize(d_num_species + 1);
                            
                            diffusivities_data[si][0] = d_data_diffusivities;
                            diffusivities_component_idx[si][0] = si*(d_num_species + 1);
                            
                            for (int sj = 0; sj < d_num_species; sj++)
                            {
                                diffusivities_data[si][1 + sj] = d_data_diffusivities;
                                diffusivities_component_idx[si][1 + sj] = si*(d_num_species + 1) + sj + 1;
                            }
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[d_num_species].resize(1);
                        diffusivities_component_idx[d_num_species].resize(1);
                        
                        // -(4/3*mu + mu_v).
                        diffusivities_data[d_num_species][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species][0] =
                            d_num_species*(d_num_species + 1);
                        
                        diffusivities_data[d_num_species + 1].resize(1);
                        diffusivities_component_idx[d_num_species + 1].resize(1);
                        
                        // -mu.
                        diffusivities_data[d_num_species + 1][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 1][0] =
                            d_num_species*(d_num_species + 1) + 2;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[d_num_species + 2].resize(
                            3 + d_num_species*(d_num_species + 1));
                        diffusivities_component_idx[d_num_species + 2].resize(
                            3 + d_num_species*(d_num_species + 1));
                        
                        // -u*(4/3*mu + mu_v).
                        diffusivities_data[d_num_species + 2][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 2][0] =
                            d_num_species*(d_num_species + 1) + 3;
                        
                        // -v*mu.
                        diffusivities_data[d_num_species + 2][1] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 2][1] =
                            d_num_species*(d_num_species + 1) + 8;
                        
                        // -kappa.
                        diffusivities_data[d_num_species + 2][2] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 2][2] =
                            d_num_species*(d_num_species + 1) + 9;
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            diffusivities_data[d_num_species + 2][3 + si*(d_num_species + 1)] =
                                d_data_diffusivities;
                            diffusivities_component_idx[d_num_species + 2][3 + si*(d_num_species + 1)] =
                                d_num_species*(d_num_species + 1) + 10 + si*(d_num_species + 1);
                            
                            for (int sj = 0; sj < d_num_species; sj++)
                            {
                                diffusivities_data[d_num_species + 2][4 + si*(d_num_species + 1) + sj] =
                                    d_data_diffusivities;
                                diffusivities_component_idx[d_num_species + 2][4 + si*(d_num_species + 1) + sj] =
                                    d_num_species*(d_num_species + 1) + 11 + si*(d_num_species + 1) + sj;
                            }
                        }
                        
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
                        diffusivities_component_idx[d_num_species][0] =
                            d_num_species*(d_num_species + 1) + 1;
                        
                        diffusivities_data[d_num_species + 1].resize(1);
                        diffusivities_component_idx[d_num_species + 1].resize(1);
                        
                        // -mu.
                        diffusivities_data[d_num_species + 1][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 1][0] =
                            d_num_species*(d_num_species + 1) + 2;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[d_num_species + 2].resize(2);
                        diffusivities_component_idx[d_num_species + 2].resize(2);
                        
                        // -v*mu.
                        diffusivities_data[d_num_species + 2][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 2][0] =
                            d_num_species*(d_num_species + 1) + 8;
                        
                        // u*(2/3*mu - mu_v).
                        diffusivities_data[d_num_species + 2][1] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 2][1] =
                            d_num_species*(d_num_species + 1) + 5;
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::getCellDataOfDiffusiveFluxDiffusivities()\n"
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
                        diffusivities_component_idx[d_num_species][0] =
                            d_num_species*(d_num_species + 1) + 2;
                        
                        diffusivities_data[d_num_species + 1].resize(1);
                        diffusivities_component_idx[d_num_species + 1].resize(1);
                        
                        // 2/3*mu - mu_v.
                        diffusivities_data[d_num_species + 1][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 1][0] =
                            d_num_species*(d_num_species + 1) + 1;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[d_num_species + 2].resize(2);
                        diffusivities_component_idx[d_num_species + 2].resize(2);
                        
                        // v*(2/3*mu - mu_v).
                        diffusivities_data[d_num_species + 2][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 2][0] =
                            d_num_species*(d_num_species + 1) + 6;
                        
                        // -u*mu.
                        diffusivities_data[d_num_species + 2][1] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 2][1] =
                            d_num_species*(d_num_species + 1) + 7;
                        
                        break;
                    }
                    case DIRECTION::Y_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            diffusivities_data[si].resize(d_num_species + 1);
                            diffusivities_component_idx[si].resize(d_num_species + 1);
                            
                            diffusivities_data[si][0] = d_data_diffusivities;
                            diffusivities_component_idx[si][0] = si*(d_num_species + 1);
                            
                            for (int sj = 0; sj < d_num_species; sj++)
                            {
                                diffusivities_data[si][1 + sj] = d_data_diffusivities;
                                diffusivities_component_idx[si][1 + sj] = si*(d_num_species + 1) + sj + 1;
                            }
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[d_num_species].resize(1);
                        diffusivities_component_idx[d_num_species].resize(1);
                        
                        // -mu.
                        diffusivities_data[d_num_species][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species][0] =
                            d_num_species*(d_num_species + 1) + 2;
                        
                        diffusivities_data[d_num_species + 1].resize(1);
                        diffusivities_component_idx[d_num_species + 1].resize(1);
                        
                        // -(4/3*mu + mu_v).
                        diffusivities_data[d_num_species + 1][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 1][0] =
                            d_num_species*(d_num_species + 1);
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[d_num_species + 2].resize(
                            3 + d_num_species*(d_num_species + 1));
                        diffusivities_component_idx[d_num_species + 2].resize(
                            3 + d_num_species*(d_num_species + 1));
                        
                        // -u*mu.
                        diffusivities_data[d_num_species + 2][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 2][0] =
                            d_num_species*(d_num_species + 1) + 7;
                        
                        // -v*(4/3*mu + mu_v).
                        diffusivities_data[d_num_species + 2][1] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 2][1] =
                            d_num_species*(d_num_species + 1) + 4;
                        
                        // -kappa.
                        diffusivities_data[d_num_species + 2][2] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 2][2] =
                            d_num_species*(d_num_species + 1) + 9;
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            diffusivities_data[d_num_species + 2][3 + si*(d_num_species + 1)] =
                                d_data_diffusivities;
                            diffusivities_component_idx[d_num_species + 2][3 + si*(d_num_species + 1)] =
                                d_num_species*(d_num_species + 1) + 10 + si*(d_num_species + 1);
                            
                            for (int sj = 0; sj < d_num_species; sj++)
                            {
                                diffusivities_data[d_num_species + 2][4 + si*(d_num_species + 1) + sj] =
                                    d_data_diffusivities;
                                diffusivities_component_idx[d_num_species + 2][4 + si*(d_num_species + 1) + sj] =
                                    d_num_species*(d_num_species + 1) + 11 + si*(d_num_species + 1) + sj;
                            }
                        }
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::getCellDataOfDiffusiveFluxDiffusivities()\n"
                            << "There are only x-direction and y-direction for two-dimensional problem."
                            << std::endl);
                    }
                }
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::getCellDataOfDiffusiveFluxDiffusivities()\n"
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
                            diffusivities_data[si].resize(d_num_species + 1);
                            diffusivities_component_idx[si].resize(d_num_species + 1);
                            
                            diffusivities_data[si][0] = d_data_diffusivities;
                            diffusivities_component_idx[si][0] = si*(d_num_species + 1);
                            
                            for (int sj = 0; sj < d_num_species; sj++)
                            {
                                diffusivities_data[si][1 + sj] = d_data_diffusivities;
                                diffusivities_component_idx[si][1 + sj] = si*(d_num_species + 1) + sj + 1;
                            }
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[d_num_species].resize(1);
                        diffusivities_component_idx[d_num_species].resize(1);
                        
                        // -(4/3*mu + mu_v).
                        diffusivities_data[d_num_species][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species][0] =
                            d_num_species*(d_num_species + 1);
                        
                        diffusivities_data[d_num_species + 1].resize(1);
                        diffusivities_component_idx[d_num_species + 1].resize(1);
                        
                        // -mu.
                        diffusivities_data[d_num_species + 1][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 1][0] =
                            d_num_species*(d_num_species + 1) + 2;
                        
                        diffusivities_data[d_num_species + 2].resize(1);
                        diffusivities_component_idx[d_num_species + 2].resize(1);
                        
                        // -mu.
                        diffusivities_data[d_num_species + 2][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 2][0] =
                            d_num_species*(d_num_species + 1) + 2;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[d_num_species + 3].resize(
                            4 + d_num_species*(d_num_species + 1));
                        diffusivities_component_idx[d_num_species + 3].resize(
                            4 + d_num_species*(d_num_species + 1));
                        
                        // -u*(4/3*mu + mu_v).
                        diffusivities_data[d_num_species + 3][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][0] =
                            d_num_species*(d_num_species + 1) + 3;
                        
                        // -v*mu.
                        diffusivities_data[d_num_species + 3][1] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][1] =
                            d_num_species*(d_num_species + 1) + 10;
                        
                        // -w*mu.
                        diffusivities_data[d_num_species + 3][2] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][2] =
                            d_num_species*(d_num_species + 1) + 11;
                        
                        // -kappa.
                        diffusivities_data[d_num_species + 3][3] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][3] =
                            d_num_species*(d_num_species + 1) + 12;
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            diffusivities_data[d_num_species + 3][4 + si*(d_num_species + 1)] =
                                d_data_diffusivities;
                            diffusivities_component_idx[d_num_species + 3][4 + si*(d_num_species + 1)] =
                                d_num_species*(d_num_species + 1) + 13 + si*(d_num_species + 1);
                            
                            for (int sj = 0; sj < d_num_species; sj++)
                            {
                                diffusivities_data[d_num_species + 3][5 + si*(d_num_species + 1) + sj] =
                                    d_data_diffusivities;
                                diffusivities_component_idx[d_num_species + 3][5 + si*(d_num_species + 1) + sj] =
                                    d_num_species*(d_num_species + 1) + 14 + si*(d_num_species + 1) + sj;
                            }
                        }
                        
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
                        diffusivities_component_idx[d_num_species][0] =
                            d_num_species*(d_num_species + 1) + 1;
                        
                        diffusivities_data[d_num_species + 1].resize(1);
                        diffusivities_component_idx[d_num_species + 1].resize(1);
                        
                        // -mu.
                        diffusivities_data[d_num_species + 1][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 1][0] =
                            d_num_species*(d_num_species + 1) + 2;
                        
                        diffusivities_data[d_num_species + 2].resize(0);
                        diffusivities_component_idx[d_num_species + 2].resize(0);
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[d_num_species + 3].resize(2);
                        diffusivities_component_idx[d_num_species + 3].resize(2);
                        
                        // -v*mu.
                        diffusivities_data[d_num_species + 3][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][0] =
                            d_num_species*(d_num_species + 1) + 10;
                        
                        // u*(2/3*mu - mu_v).
                        diffusivities_data[d_num_species + 3][1] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][1] =
                            d_num_species*(d_num_species + 1) + 6;
                        
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
                        diffusivities_component_idx[d_num_species][0] =
                            d_num_species*(d_num_species + 1) + 1;
                        
                        diffusivities_data[d_num_species + 1].resize(0);
                        diffusivities_component_idx[d_num_species + 1].resize(0);
                        
                        diffusivities_data[d_num_species + 2].resize(1);
                        diffusivities_component_idx[d_num_species + 2].resize(1);
                        
                        // -mu.
                        diffusivities_data[d_num_species + 2][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 2][0] =
                            d_num_species*(d_num_species + 1) + 2;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[d_num_species + 3].resize(2);
                        diffusivities_component_idx[d_num_species + 3].resize(2);
                        
                        // -w*mu.
                        diffusivities_data[d_num_species + 3][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][0] =
                            d_num_species*(d_num_species + 1) + 11;
                        
                        // u*(2/3*mu - mu_v).
                        diffusivities_data[d_num_species + 3][1] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][1] =
                            d_num_species*(d_num_species + 1) + 6;
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::getCellDataOfDiffusiveFluxDiffusivities()\n"
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
                        diffusivities_component_idx[d_num_species][0] =
                            d_num_species*(d_num_species + 1) + 2;
                        
                        diffusivities_data[d_num_species + 1].resize(1);
                        diffusivities_component_idx[d_num_species + 1].resize(1);
                        
                        // 2/3*mu - mu_v.
                        diffusivities_data[d_num_species + 1][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 1][0] =
                            d_num_species*(d_num_species + 1) + 1;
                        
                        diffusivities_data[d_num_species + 2].resize(0);
                        diffusivities_component_idx[d_num_species + 2].resize(0);
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[d_num_species + 3].resize(2);
                        diffusivities_component_idx[d_num_species + 3].resize(2);
                        
                        // v*(2/3*mu - mu_v).
                        diffusivities_data[d_num_species + 3][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][0] =
                            d_num_species*(d_num_species + 1) + 7;
                        
                        // -u*mu.
                        diffusivities_data[d_num_species + 3][1] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][1] =
                            d_num_species*(d_num_species + 1) + 9;
                        
                        break;
                    }
                    case DIRECTION::Y_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            diffusivities_data[si].resize(d_num_species + 1);
                            diffusivities_component_idx[si].resize(d_num_species + 1);
                            
                            diffusivities_data[si][0] = d_data_diffusivities;
                            diffusivities_component_idx[si][0] = si*(d_num_species + 1);
                            
                            for (int sj = 0; sj < d_num_species; sj++)
                            {
                                diffusivities_data[si][1 + sj] = d_data_diffusivities;
                                diffusivities_component_idx[si][1 + sj] = si*(d_num_species + 1) + sj + 1;
                            }
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[d_num_species].resize(1);
                        diffusivities_component_idx[d_num_species].resize(1);
                        
                        // -mu.
                        diffusivities_data[d_num_species][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species][0] =
                            d_num_species*(d_num_species + 1) + 2;
                        
                        diffusivities_data[d_num_species + 1].resize(1);
                        diffusivities_component_idx[d_num_species + 1].resize(1);
                        
                        // -(4/3*mu + mu_v).
                        diffusivities_data[d_num_species + 1][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 1][0] =
                            d_num_species*(d_num_species + 1);
                        
                        diffusivities_data[d_num_species + 2].resize(1);
                        diffusivities_component_idx[d_num_species + 2].resize(1);
                        
                        // -mu.
                        diffusivities_data[d_num_species + 2][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 2][0] =
                            d_num_species*(d_num_species + 1) + 2;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[d_num_species + 3].resize(
                            4 + d_num_species*(d_num_species + 1));
                        diffusivities_component_idx[d_num_species + 3].resize(
                            4 + d_num_species*(d_num_species + 1));
                        
                        // -u*mu.
                        diffusivities_data[d_num_species + 3][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][0] =
                            d_num_species*(d_num_species + 1) + 9;
                        
                        // -v*(4/3*mu + mu_v).
                        diffusivities_data[d_num_species + 3][1] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][1] =
                            d_num_species*(d_num_species + 1) + 4;
                        
                        // -w*mu.
                        diffusivities_data[d_num_species + 3][2] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][2] =
                            d_num_species*(d_num_species + 1) + 11;
                        
                        // -kappa.
                        diffusivities_data[d_num_species + 3][3] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][3] =
                            d_num_species*(d_num_species + 1) + 12;
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            diffusivities_data[d_num_species + 3][4 + si*(d_num_species + 1)] =
                                d_data_diffusivities;
                            diffusivities_component_idx[d_num_species + 3][4 + si*(d_num_species + 1)] =
                                d_num_species*(d_num_species + 1) + 13 + si*(d_num_species + 1);
                            
                            for (int sj = 0; sj < d_num_species; sj++)
                            {
                                diffusivities_data[d_num_species + 3][5 + si*(d_num_species + 1) + sj] =
                                    d_data_diffusivities;
                                diffusivities_component_idx[d_num_species + 3][5 + si*(d_num_species + 1) + sj] =
                                    d_num_species*(d_num_species + 1) + 14 + si*(d_num_species + 1) + sj;
                            }
                        }
                        
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
                        diffusivities_component_idx[d_num_species + 1][0] =
                            d_num_species*(d_num_species + 1) + 1;
                        
                        diffusivities_data[d_num_species + 2].resize(1);
                        diffusivities_component_idx[d_num_species + 2].resize(1);
                        
                        // -mu.
                        diffusivities_data[d_num_species + 2][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 2][0] =
                            d_num_species*(d_num_species + 1) + 2;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[d_num_species + 3].resize(2);
                        diffusivities_component_idx[d_num_species + 3].resize(2);
                        
                        // -w*u.
                        diffusivities_data[d_num_species + 3][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][0] =
                            d_num_species*(d_num_species + 1) + 11;
                        
                        // v*(2/3*mu - mu_v).
                        diffusivities_data[d_num_species + 3][1] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][1] =
                            d_num_species*(d_num_species + 1) + 7;
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::getCellDataOfDiffusiveFluxDiffusivities()\n"
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
                        diffusivities_component_idx[d_num_species][0] =
                            d_num_species*(d_num_species + 1) + 2;
                        
                        diffusivities_data[d_num_species + 1].resize(0);
                        diffusivities_component_idx[d_num_species + 1].resize(0);
                        
                        diffusivities_data[d_num_species + 2].resize(1);
                        diffusivities_component_idx[d_num_species + 2].resize(1);
                        
                        // 2/3*mu - mu_v.
                        diffusivities_data[d_num_species + 2][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 2][0] =
                            d_num_species*(d_num_species + 1) + 1;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[d_num_species + 3].resize(2);
                        diffusivities_component_idx[d_num_species + 3].resize(2);
                        
                        // w*(2/3*mu - mu_v).
                        diffusivities_data[d_num_species + 3][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][0] =
                            d_num_species*(d_num_species + 1) + 8;
                        
                        // -u*mu.
                        diffusivities_data[d_num_species + 3][1] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][1] =
                            d_num_species*(d_num_species + 1) + 9;
                        
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
                        diffusivities_component_idx[d_num_species + 1][0] =
                            d_num_species*(d_num_species + 1) + 2;
                        
                        diffusivities_data[d_num_species + 2].resize(1);
                        diffusivities_component_idx[d_num_species + 2].resize(1);
                        
                        // 2/3*mu - mu_v.
                        diffusivities_data[d_num_species + 2][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 2][0] =
                            d_num_species*(d_num_species + 1) + 1;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[d_num_species + 3].resize(2);
                        diffusivities_component_idx[d_num_species + 3].resize(2);
                        
                        // w*(2/3*mu - mu_v).
                        diffusivities_data[d_num_species + 3][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][0] =
                            d_num_species*(d_num_species + 1) + 8;
                        
                        // -v*mu.
                        diffusivities_data[d_num_species + 3][1] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][1] =
                            d_num_species*(d_num_species + 1) + 10;
                        
                        break;
                    }
                    case DIRECTION::Z_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            diffusivities_data[si].resize(d_num_species + 1);
                            diffusivities_component_idx[si].resize(d_num_species + 1);
                            
                            diffusivities_data[si][0] = d_data_diffusivities;
                            diffusivities_component_idx[si][0] = si*(d_num_species + 1);
                            
                            for (int sj = 0; sj < d_num_species; sj++)
                            {
                                diffusivities_data[si][1 + sj] = d_data_diffusivities;
                                diffusivities_component_idx[si][1 + sj] = si*(d_num_species + 1) + sj + 1;
                            }
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[d_num_species].resize(1);
                        diffusivities_component_idx[d_num_species].resize(1);
                        
                        // -mu.
                        diffusivities_data[d_num_species][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species][0] =
                            d_num_species*(d_num_species + 1) + 2;
                        
                        diffusivities_data[d_num_species + 1].resize(1);
                        diffusivities_component_idx[d_num_species + 1].resize(1);
                        
                        // -mu.
                        diffusivities_data[d_num_species + 1][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 1][0] =
                            d_num_species*(d_num_species + 1) + 2;
                        
                        diffusivities_data[d_num_species + 2].resize(1);
                        diffusivities_component_idx[d_num_species + 2].resize(1);
                        
                        // -(4/3*mu + mu_v).
                        diffusivities_data[d_num_species + 2][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 2][0] =
                            d_num_species*(d_num_species + 1);
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[d_num_species + 3].resize(
                            4 + d_num_species*(d_num_species + 1));
                        diffusivities_component_idx[d_num_species + 3].resize(
                            4 + d_num_species*(d_num_species + 1));
                        
                        // -u*mu.
                        diffusivities_data[d_num_species + 3][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][0] =
                            d_num_species*(d_num_species + 1) + 9;
                        
                        // -v*mu.
                        diffusivities_data[d_num_species + 3][1] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][1] =
                            d_num_species*(d_num_species + 1) + 10;
                        
                        // -w*(4/3*mu + mu_v).
                        diffusivities_data[d_num_species + 3][2] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][2] =
                            d_num_species*(d_num_species + 1) + 5;
                        
                        // -kappa.
                        diffusivities_data[d_num_species + 3][3] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][3] =
                            d_num_species*(d_num_species + 1) + 12;
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            diffusivities_data[d_num_species + 3][4 + si*(d_num_species + 1)] =
                                d_data_diffusivities;
                            diffusivities_component_idx[d_num_species + 3][4 + si*(d_num_species + 1)] =
                                d_num_species*(d_num_species + 1) + 13 + si*(d_num_species + 1);
                            
                            for (int sj = 0; sj < d_num_species; sj++)
                            {
                                diffusivities_data[d_num_species + 3][5 + si*(d_num_species + 1) + sj] =
                                    d_data_diffusivities;
                                diffusivities_component_idx[d_num_species + 3][5 + si*(d_num_species + 1) + sj] =
                                    d_num_species*(d_num_species + 1) + 14 + si*(d_num_species + 1) + sj;
                            }
                        }
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::getCellDataOfDiffusiveFluxDiffusivities()\n"
                            << "There are only x-direction, y-direction and z-direction for three-dimensional problem."
                            << std::endl);
                    }
                }
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::getCellDataOfDiffusiveFluxDiffusivities()\n"
                    << "There are only x-direction, y-direction and z-direction for three-dimensional problem."
                    << std::endl);
            }
        }
    }
}
