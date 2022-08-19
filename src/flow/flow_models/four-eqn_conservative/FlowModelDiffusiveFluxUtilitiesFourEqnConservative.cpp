#include "flow/flow_models/four-eqn_conservative/FlowModelDiffusiveFluxUtilitiesFourEqnConservative.hpp"

FlowModelDiffusiveFluxUtilitiesFourEqnConservative::FlowModelDiffusiveFluxUtilitiesFourEqnConservative(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
    const int& num_species,
    const HAMERS_SHARED_PTR<EquationOfMassDiffusivityMixingRules> equation_of_mass_diffusivity_mixing_rules,
    const HAMERS_SHARED_PTR<EquationOfShearViscosityMixingRules> equation_of_shear_viscosity_mixing_rules,
    const HAMERS_SHARED_PTR<EquationOfBulkViscosityMixingRules> equation_of_bulk_viscosity_mixing_rules,
    const HAMERS_SHARED_PTR<EquationOfThermalConductivityMixingRules> equation_of_thermal_conductivity_mixing_rules):
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
        d_subghost_box_mass_diffusivities(hier::Box::getEmptyBox(d_dim)),
        d_subghost_box_shear_viscosity(hier::Box::getEmptyBox(d_dim)),
        d_subghost_box_bulk_viscosity(hier::Box::getEmptyBox(d_dim)),
        d_subghost_box_thermal_conductivity(hier::Box::getEmptyBox(d_dim)),
        d_subghostcell_dims_mass_diffusivities(hier::IntVector::getZero(d_dim)),
        d_subghostcell_dims_shear_viscosity(hier::IntVector::getZero(d_dim)),
        d_subghostcell_dims_bulk_viscosity(hier::IntVector::getZero(d_dim)),
        d_subghostcell_dims_thermal_conductivity(hier::IntVector::getZero(d_dim)),
        d_cell_data_computed_mass_diffusivities(false),
        d_cell_data_computed_shear_viscosity(false),
        d_cell_data_computed_bulk_viscosity(false),
        d_cell_data_computed_thermal_conductivity(false),
        d_equation_of_mass_diffusivity_mixing_rules(equation_of_mass_diffusivity_mixing_rules),
        d_equation_of_shear_viscosity_mixing_rules(equation_of_shear_viscosity_mixing_rules),
        d_equation_of_bulk_viscosity_mixing_rules(equation_of_bulk_viscosity_mixing_rules),
        d_equation_of_thermal_conductivity_mixing_rules(equation_of_thermal_conductivity_mixing_rules)
{}


/*
 * Register different derived variables related to this class in the registered patch. The
 * derived variables to be registered are given as entries in a map of the variable name to
 * the number of sub-ghost cells required.
 */
void
FlowModelDiffusiveFluxUtilitiesFourEqnConservative::registerDerivedVariables(
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
            << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::"
            << "registerDerivedVariables()\n"
            << "No patch is registered yet."
            << std::endl);
    }
    
    // Check whether all or part of derived cell data is alredy computed.
    if (d_derived_cell_data_computed)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::registerDerivedVariables()\n"
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
                << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::registerDerivedVariables()\n"
                << "The number of sub-ghost cells of variable '"
                << it->first
                << "' is not between zero and number of ghosts of conservative variables."
                << std::endl);
        }
    }
    
    if (num_subghosts_of_data.find("MASS_DIFFUSIVITIES") != num_subghosts_of_data.end())
    {
        std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data_flow_model;
        
        num_subghosts_of_data_flow_model.insert(std::pair<std::string, hier::IntVector>(
            "MASS_FRACTIONS",
            num_subghosts_of_data.find("MASS_DIFFUSIVITIES")->second));
        
        num_subghosts_of_data_flow_model.insert(std::pair<std::string, hier::IntVector>(
            "PRESSURE",
            num_subghosts_of_data.find("MASS_DIFFUSIVITIES")->second));
        
        num_subghosts_of_data_flow_model.insert(std::pair<std::string, hier::IntVector>(
            "TEMPERATURE",
            num_subghosts_of_data.find("MASS_DIFFUSIVITIES")->second));
        
        flow_model_tmp->registerDerivedVariables(num_subghosts_of_data_flow_model);
        
        setNumberOfSubGhosts(
            num_subghosts_of_data.find("MASS_DIFFUSIVITIES")->second,
            "MASS_DIFFUSIVITIES",
            "MASS_DIFFUSIVITIES");
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
            "MASS_FRACTIONS",
            num_subghosts_of_data.find("BULK_VISCOSITY")->second));
        
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
            "MASS_FRACTIONS",
            num_subghosts_of_data.find("THERMAL_CONDUCTIVITY")->second));
        
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
FlowModelDiffusiveFluxUtilitiesFourEqnConservative::registerDerivedVariablesForDiffusiveFluxes(
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
            << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::"
            << "registerDerivedVariablesForDiffusiveFluxes()\n"
            << "No patch is registered yet."
            << std::endl);
    }
    
    // Check whether all or part of derived cell data is alredy computed.
    if (d_derived_cell_data_computed)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::registerDerivedVariablesForDiffusiveFluxes()\n"
            << "Derived cell data is already computed."
            << std::endl);
    }
    
    if ((num_subghosts < hier::IntVector::getZero(d_dim)) ||
        (num_subghosts > flow_model_tmp->getNumberOfGhostCells()))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::registerDerivedVariablesForDiffusiveFluxes()\n"
            << "The number of sub-ghost cells of variable is not between zero and number of ghosts of conservative variables."
            << std::endl);
    }
    
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
FlowModelDiffusiveFluxUtilitiesFourEqnConservative::allocateMemoryForDerivedCellData()
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
    
    if (d_num_subghosts_mass_diffusivities > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_mass_diffusivities)
        {
            if (!d_data_mass_diffusivities)
            {
                // Create the cell data of mass diffusivities.
                d_data_mass_diffusivities.reset(new pdat::CellData<double>(
                    interior_box, d_num_species, d_num_subghosts_mass_diffusivities));
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::allocateMemoryForDerivedCellData()\n"
                << "Cell data of 'MASS_DIFFUSIVITIES' is aleady computed."
                << std::endl);
        }
    }
    
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
                << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::allocateMemoryForDerivedCellData()\n"
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
                << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::allocateMemoryForDerivedCellData()\n"
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
                // Create the cell data of thermal conductivity.
                d_data_thermal_conductivity.reset(new pdat::CellData<double>(
                    interior_box, 1, d_num_subghosts_thermal_conductivity));
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::allocateMemoryForDerivedCellData()\n"
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
                        d_data_diffusivities.reset(new pdat::CellData<double>(
                            interior_box,
                            2*d_num_species*(d_num_species + 1) + 3,
                            d_num_subghosts_diffusivities));
                    }
                    else if (d_dim == tbox::Dimension(2))
                    {
                        d_data_diffusivities.reset(new pdat::CellData<double>(
                            interior_box,
                            2*d_num_species*(d_num_species + 1) + 10,
                            d_num_subghosts_diffusivities));
                    }
                    else if (d_dim == tbox::Dimension(3))
                    {
                        d_data_diffusivities.reset(new pdat::CellData<double>(
                            interior_box,
                            2*d_num_species*(d_num_species + 1) + 13,
                            d_num_subghosts_diffusivities));
                    }
                }
            }
            else
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::allocateMemoryForDerivedCellData()\n"
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
FlowModelDiffusiveFluxUtilitiesFourEqnConservative::allocateMemoryForSideDataOfDiffusiveFluxDiffusivities()
{
    if (!d_need_side_diffusivities)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::allocateMemoryForSideDataOfDiffusiveFluxDiffusivities()\n"
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
 * Clear cell and side data of different derived variables related to this class in the registered patch.
 */
void
FlowModelDiffusiveFluxUtilitiesFourEqnConservative::clearCellAndSideData()
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
    
    d_side_data_diffusivities.reset();
    
    d_cell_data_computed_diffusivities        = false;
    d_cell_data_computed_mass_diffusivities   = false;
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
FlowModelDiffusiveFluxUtilitiesFourEqnConservative::computeDerivedCellData()
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
            << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::"
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
    
    // Compute the mass diffusivities cell data.
    if (d_num_subghosts_mass_diffusivities > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_mass_diffusivities)
        {
            computeCellDataOfMassDiffusivities();
        }
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
FlowModelDiffusiveFluxUtilitiesFourEqnConservative::getCellData(const std::string& variable_key)
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
            << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::"
            << "getCellData()\n"
            << "No patch is registered yet."
            << std::endl);
    }
    
    HAMERS_SHARED_PTR<pdat::CellData<double> > cell_data;
    
    if (variable_key == "MASS_DIFFUSIVITIES")
    {
        if (!d_cell_data_computed_mass_diffusivities)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::getCellData()\n"
                << "Cell data of 'MASS_DIFFUSIVITIES' is not registered/computed yet."
                << std::endl);
        }
        cell_data = d_data_mass_diffusivities;
    }
    else if (variable_key == "SHEAR_VISCOSITY")
    {
        if (!d_cell_data_computed_shear_viscosity)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::getCellData()\n"
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
                << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::getCellData()\n"
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
                << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::getCellData()\n"
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
std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >
FlowModelDiffusiveFluxUtilitiesFourEqnConservative::getCellData(
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
FlowModelDiffusiveFluxUtilitiesFourEqnConservative::getCellDataOfDiffusiveFluxVariablesForDerivative(
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
    
    // Get the cell data of mass fractions.
    HAMERS_SHARED_PTR<pdat::CellData<double> > data_mass_fractions =
        flow_model_tmp->getCellData("MASS_FRACTIONS");
    
    // Get the cell data of velocity.
    HAMERS_SHARED_PTR<pdat::CellData<double> > data_velocity =
        flow_model_tmp->getCellData("VELOCITY");
    
    // Get the cell data of temperature.
    HAMERS_SHARED_PTR<pdat::CellData<double> > data_temperature =
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
                            << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::"
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
                    << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::"
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
                            << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::"
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
                            << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::"
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
                    << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::"
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
                            << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::"
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
                            << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::"
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
                            << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::"
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
                    << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::"
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
FlowModelDiffusiveFluxUtilitiesFourEqnConservative::getCellDataOfDiffusiveFluxDiffusivities(
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
            << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::"
            << "getCellDataOfDiffusiveFluxDiffusivities()\n"
            << "No patch is registered yet."
            << std::endl);
    }
    
    if (!d_cell_data_computed_diffusivities)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::getCellDataOfDiffusiveFluxDiffusivities()\n"
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
                            << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::"
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
                    << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::"
                    << "getCellDataOfDiffusiveFluxDiffusivities()\n"
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
                            << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::"
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
                            << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::"
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
                    << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::"
                    << "getCellDataOfDiffusiveFluxDiffusivities()\n"
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
                            << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::"
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
                            << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::"
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
                            << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::"
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
                    << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::"
                    << "getCellDataOfDiffusiveFluxDiffusivities()\n"
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
FlowModelDiffusiveFluxUtilitiesFourEqnConservative::getCellDataForInterpolationToSideDataForDiffusiveFluxDiffusivities(
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
            << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::"
            << "getCellDataForInterpolationToSideDataForDiffusiveFluxDiffusivities()\n"
            << "No patch is registered yet."
            << std::endl);
    }
    
    if (!d_cell_data_computed_mass_diffusivities)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::getCellDataForInterpolationToSideDataForDiffusiveFluxDiffusivities()\n"
            << "Cell data of 'MASS_DIFFUSIVITIES' is not registered/computed yet."
            << std::endl);
    }
    
    if (!d_cell_data_computed_shear_viscosity)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::getCellDataForInterpolationToSideDataForDiffusiveFluxDiffusivities()\n"
            << "Cell data of 'SHEAR_VISCOSITY' is not registered/computed yet."
            << std::endl);
    }
    
    if (!d_cell_data_computed_bulk_viscosity)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::getCellDataForInterpolationToSideDataForDiffusiveFluxDiffusivities()\n"
            << "Cell data of 'BULK_VISCOSITY' is not registered/computed yet."
            << std::endl);
    }
    
    if (!d_cell_data_computed_thermal_conductivity)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::getCellDataForInterpolationToSideDataForDiffusiveFluxDiffusivities()\n"
            << "Cell data of 'THERMAL_CONDUCTIVITY' is not registered/computed yet."
            << std::endl);
    }
    
    var_data_for_diffusivities.resize(3*d_num_species + 4 + d_dim.getValue());
    var_data_for_diffusivities_component_idx.resize(3*d_num_species + 4 + d_dim.getValue());
    
    for (int si = 0; si < d_num_species; si++)
    {
        var_data_for_diffusivities[si] = d_data_mass_diffusivities;
        var_data_for_diffusivities_component_idx[si] = si;
    }
    
    var_data_for_diffusivities[d_num_species] = d_data_shear_viscosity;
    var_data_for_diffusivities_component_idx[d_num_species] = 0;
    var_data_for_diffusivities[d_num_species + 1] = d_data_bulk_viscosity;
    var_data_for_diffusivities_component_idx[d_num_species + 1] = 0;
    var_data_for_diffusivities[d_num_species + 2] = d_data_thermal_conductivity;
    var_data_for_diffusivities_component_idx[d_num_species + 2] = 0;
    
    // Get the cell data of density.
    HAMERS_SHARED_PTR<pdat::CellData<double> > data_density =
        flow_model_tmp->getCellData("DENSITY");
    
    // Get the cell data of mass fractions.
    HAMERS_SHARED_PTR<pdat::CellData<double> > data_mass_fractions =
        flow_model_tmp->getCellData("MASS_FRACTIONS");
    
    // Get the cell data of velocity.
    HAMERS_SHARED_PTR<pdat::CellData<double> > data_velocity =
        flow_model_tmp->getCellData("VELOCITY");
    
    // Get the cell data of species enthalpies.
    std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > data_species_enthalpies =
        flow_model_tmp->getSpeciesCellData("SPECIES_ENTHALPIES");
    
    var_data_for_diffusivities[d_num_species + 3] = data_density;
    var_data_for_diffusivities_component_idx[d_num_species + 3] = 0;
    
    for (int si = 0; si < d_num_species; si++)
    {
        var_data_for_diffusivities[d_num_species + 4 + si] = data_mass_fractions;
        var_data_for_diffusivities_component_idx[d_num_species + 4 + si] = si;
    }
    
    if (d_dim == tbox::Dimension(1))
    {
        var_data_for_diffusivities[2*d_num_species + 4] = data_velocity;
        var_data_for_diffusivities_component_idx[2*d_num_species + 4] = 0;
    }
    else if (d_dim == tbox::Dimension(2))
    {
        var_data_for_diffusivities[2*d_num_species + 4] = data_velocity;
        var_data_for_diffusivities_component_idx[2*d_num_species + 4] = 0;
        var_data_for_diffusivities[2*d_num_species + 5] = data_velocity;
        var_data_for_diffusivities_component_idx[2*d_num_species + 5] = 1;
    }
    else if (d_dim == tbox::Dimension(3))
    {
        var_data_for_diffusivities[2*d_num_species + 4] = data_velocity;
        var_data_for_diffusivities_component_idx[2*d_num_species + 4] = 0;
        var_data_for_diffusivities[2*d_num_species + 5] = data_velocity;
        var_data_for_diffusivities_component_idx[2*d_num_species + 5] = 1;
        var_data_for_diffusivities[2*d_num_species + 6] = data_velocity;
        var_data_for_diffusivities_component_idx[2*d_num_species + 6] = 2;
    }
    
    for (int si = 0; si < d_num_species; si++)
    {
        var_data_for_diffusivities[2*d_num_species + 4 + d_dim.getValue() + si] = data_species_enthalpies[si];
        var_data_for_diffusivities_component_idx[2*d_num_species + 4 + d_dim.getValue() + si] = 0;
    }
}


/*
 * Compute the side data of the diffusivities in the diffusive flux with the interpolated side data.
 */
void
FlowModelDiffusiveFluxUtilitiesFourEqnConservative::computeSideDataOfDiffusiveFluxDiffusivities(
    const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& var_data_for_diffusivities)
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
            TBOX_ASSERT(static_cast<int>(var_data_for_diffusivities.size()) == 3*d_num_species + 4 + d_dim.getValue());
            TBOX_ASSERT(num_ghosts <= d_num_subghosts_diffusivities);
            
            for (int vi = 0; vi < static_cast<int>(var_data_for_diffusivities.size()); vi++)
            {
                TBOX_ASSERT(var_data_for_diffusivities[vi] == num_ghosts);
                TBOX_ASSERT(var_data_for_diffusivities[vi]->getGhostBox().contains(interior_box));
            }
#endif
            
            if (d_dim == tbox::Dimension(1))
            {
                std::vector<double*> D_x;
                D_x.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    D_x.push_back(var_data_for_diffusivities[si]->getPointer(0, 0));
                }
                
                double* mu_x    = var_data_for_diffusivities[d_num_species    ]->getPointer(0, 0);
                double* mu_v_x  = var_data_for_diffusivities[d_num_species + 1]->getPointer(0, 0);
                double* kappa_x = var_data_for_diffusivities[d_num_species + 2]->getPointer(0, 0);
                
                double* rho_x = var_data_for_diffusivities[d_num_species + 3]->getPointer(0, 0);
                
                std::vector<double*> Y_x;
                Y_x.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    Y_x.push_back(var_data_for_diffusivities[d_num_species + 4 + si]->getPointer(0, 0));
                }
                
                double* u_x = var_data_for_diffusivities[2*d_num_species + 4]->getPointer(0, 0);
                
                std::vector<double*> h_i_x;
                h_i_x.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    h_i_x.push_back(var_data_for_diffusivities[2*d_num_species + 5 + si]->getPointer(0, 0));
                }
                
                std::vector<double*> D_ptr;
                D_ptr.reserve(2*d_num_species*(d_num_species + 1) + 8);
                
                /*
                 * Compute the diffusivities for the x-direction derivatives.
                 */
                
                for (int i = 0; i < 2*d_num_species*(d_num_species + 1) + 3; i++)
                {
                    D_ptr.push_back(d_side_data_diffusivities->getPointer(0, i));
                }
                
                // Mass equations.
                for (int si = 0; si < d_num_species; si++)
                {
                    const int component_idx = si*(d_num_species + 1);
                    
                    for (int i = -num_ghosts[0]; i < interior_dims[0] + num_ghosts[0] + 1; i++)
                    {
                        const int idx_diffusivities = i + d_num_subghosts_diffusivities[0];
                        const int idx_var_data      = i + num_ghosts[0];
                        
                        D_ptr[component_idx][idx_diffusivities] =
                            -rho_x[idx_var_data]*D_x[si][idx_var_data];
                    }
                }
                
                for (int si = 0; si < d_num_species; si++)
                {
                    for (int sj = 0; sj < d_num_species; sj++)
                    {
                        const int component_idx = si*(d_num_species + 1) + sj + 1;
                        
                        for (int i = -num_ghosts[0]; i < interior_dims[0] + num_ghosts[0] + 1; i++)
                        {
                            const int idx_diffusivities = i + d_num_subghosts_diffusivities[0];
                            const int idx_var_data      = i + num_ghosts[0];
                            
                            D_ptr[component_idx][idx_diffusivities] =
                                rho_x[idx_var_data]*Y_x[si][idx_var_data]*D_x[sj][idx_var_data];
                        }
                    }
                }
                
                // Momentum and energy equations.
                for (int i = -num_ghosts[0]; i < interior_dims[0] + num_ghosts[0] + 1; i++)
                {
                    const int idx_diffusivities = i + d_num_subghosts_diffusivities[0];
                    const int idx_var_data      = i + num_ghosts[0];
                    
                    D_ptr[d_num_species*(d_num_species + 1) + 0][idx_diffusivities] = -(double(4)/double(3)*mu_x[idx_var_data] + mu_v_x[idx_var_data]);
                    D_ptr[d_num_species*(d_num_species + 1) + 1][idx_diffusivities] = -u_x[idx_var_data]*(double(4)/double(3)*mu_x[idx_var_data] +
                        mu_v_x[idx_var_data]);
                    D_ptr[d_num_species*(d_num_species + 1) + 2][idx_diffusivities] = -kappa_x[idx_var_data];
                }
                
                // Energy equation.
                for (int si = 0; si < d_num_species; si++)
                {
                    const int component_idx =
                        d_num_species*(d_num_species + 1) + 3 + si*(d_num_species + 1);
                    
                    for (int i = -num_ghosts[0]; i < interior_dims[0] + num_ghosts[0] + 1; i++)
                    {
                        const int idx_diffusivities = i + d_num_subghosts_diffusivities[0];
                        const int idx_var_data      = i + num_ghosts[0];
                        
                        D_ptr[component_idx][idx_diffusivities] =
                            -rho_x[idx_var_data]*D_x[si][idx_var_data]*
                                h_i_x[si][idx_var_data];
                    }
                }
                
                for (int si = 0; si < d_num_species; si++)
                {
                    for (int sj = 0; sj < d_num_species; sj++)
                    {
                        const int component_idx =
                            d_num_species*(d_num_species + 1) + 4 + si*(d_num_species + 1) + sj;
                        
                        for (int i = -num_ghosts[0]; i < interior_dims[0] + num_ghosts[0] + 1; i++)
                        {
                            const int idx_diffusivities = i + d_num_subghosts_diffusivities[0];
                            const int idx_var_data      = i + num_ghosts[0];
                            
                            D_ptr[component_idx][idx_diffusivities] =
                                rho_x[idx_var_data]*Y_x[si][idx_var_data]*D_x[sj][idx_var_data]*
                                    h_i_x[si][idx_var_data];
                        }
                    }
                }
                
                D_ptr.clear();
            }
            else if (d_dim == tbox::Dimension(2))
            {
                std::vector<double*> D_x;
                std::vector<double*> D_y;
                D_x.reserve(d_num_species);
                D_y.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    D_x.push_back(var_data_for_diffusivities[si]->getPointer(0, 0));
                    D_y.push_back(var_data_for_diffusivities[si]->getPointer(1, 0));
                }
                
                double* mu_x    = var_data_for_diffusivities[d_num_species    ]->getPointer(0, 0);
                double* mu_v_x  = var_data_for_diffusivities[d_num_species + 1]->getPointer(0, 0);
                double* kappa_x = var_data_for_diffusivities[d_num_species + 2]->getPointer(0, 0);
                
                double* mu_y    = var_data_for_diffusivities[d_num_species    ]->getPointer(1, 0);
                double* mu_v_y  = var_data_for_diffusivities[d_num_species + 1]->getPointer(1, 0);
                double* kappa_y = var_data_for_diffusivities[d_num_species + 2]->getPointer(1, 0);
                
                double* rho_x = var_data_for_diffusivities[d_num_species + 3]->getPointer(0, 0);
                double* rho_y = var_data_for_diffusivities[d_num_species + 3]->getPointer(1, 0);
                
                std::vector<double*> Y_x;
                std::vector<double*> Y_y;
                Y_x.reserve(d_num_species);
                Y_y.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    Y_x.push_back(var_data_for_diffusivities[d_num_species + 4 + si]->getPointer(0, 0));
                    Y_y.push_back(var_data_for_diffusivities[d_num_species + 4 + si]->getPointer(1, 0));
                }
                
                double* u_x = var_data_for_diffusivities[2*d_num_species + 4]->getPointer(0, 0);
                double* v_x = var_data_for_diffusivities[2*d_num_species + 5]->getPointer(0, 0);
                
                double* u_y = var_data_for_diffusivities[2*d_num_species + 4]->getPointer(1, 0);
                double* v_y = var_data_for_diffusivities[2*d_num_species + 5]->getPointer(1, 0);
                
                std::vector<double*> h_i_x;
                std::vector<double*> h_i_y;
                h_i_x.reserve(d_num_species);
                h_i_y.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    h_i_x.push_back(var_data_for_diffusivities[2*d_num_species + 6 + si]->getPointer(0, 0));
                    h_i_y.push_back(var_data_for_diffusivities[2*d_num_species + 6 + si]->getPointer(1, 0));
                }
                
                std::vector<double*> D_ptr;
                D_ptr.reserve(2*d_num_species*(d_num_species + 1) + 8);
                
                /*
                 * Compute the diffusivities for the x-direction derivatives.
                 */
                
                for (int i = 0; i < 2*d_num_species*(d_num_species + 1) + 8; i++)
                {
                    D_ptr.push_back(d_side_data_diffusivities->getPointer(0, i));
                }
                
                // Mass equations.
                for (int si = 0; si < d_num_species; si++)
                {
                    const int component_idx = si*(d_num_species + 1);
                    
                    for (int j = -num_ghosts[1]; j < interior_dims[1] + num_ghosts[1]; j++)
                    {
                        for (int i = -num_ghosts[0]; i < interior_dims[0] + num_ghosts[0] + 1; i++)
                        {
                            const int idx_diffusivities = (i + d_num_subghosts_diffusivities[0]) +
                                (j + d_num_subghosts_diffusivities[1])*(d_subghostcell_dims_diffusivities[0] + 1);
                            
                            const int idx_var_data = (i + num_ghosts[0]) +
                                (j + num_ghosts[1])*(ghostcell_dims[0] + 1);
                            
                            D_ptr[component_idx][idx_diffusivities] =
                                -rho_x[idx_var_data]*D_x[si][idx_var_data];
                        }
                    }
                }
                
                for (int si = 0; si < d_num_species; si++)
                {
                    for (int sj = 0; sj < d_num_species; sj++)
                    {
                        const int component_idx = si*(d_num_species + 1) + sj + 1;
                        
                        for (int j = -num_ghosts[1]; j < interior_dims[1] + num_ghosts[1]; j++)
                        {
                            for (int i = -num_ghosts[0]; i < interior_dims[0] + num_ghosts[0] + 1; i++)
                            {
                                const int idx_diffusivities = (i + d_num_subghosts_diffusivities[0]) +
                                    (j + d_num_subghosts_diffusivities[1])*(d_subghostcell_dims_diffusivities[0] + 1);
                                
                                const int idx_var_data = (i + num_ghosts[0]) +
                                    (j + num_ghosts[1])*(ghostcell_dims[0] + 1);
                                
                                D_ptr[component_idx][idx_diffusivities] =
                                    rho_x[idx_var_data]*Y_x[si][idx_var_data]*D_x[sj][idx_var_data];
                            }
                        }
                    }
                }
                
                // Momentum and energy equations.
                for (int j = -num_ghosts[1]; j < interior_dims[1] + num_ghosts[1]; j++)
                {
                    for (int i = -num_ghosts[0]; i < interior_dims[0] + num_ghosts[0] + 1; i++)
                    {
                        const int idx_diffusivities = (i + d_num_subghosts_diffusivities[0]) +
                            (j + d_num_subghosts_diffusivities[1])*(d_subghostcell_dims_diffusivities[0] + 1);
                        
                        const int idx_var_data = (i + num_ghosts[0]) +
                            (j + num_ghosts[1])*(ghostcell_dims[0] + 1);
                        
                        D_ptr[d_num_species*(d_num_species + 1) + 0][idx_diffusivities] = -(double(4)/double(3)*mu_x[idx_var_data] + mu_v_x[idx_var_data]);
                        D_ptr[d_num_species*(d_num_species + 1) + 1][idx_diffusivities] = double(2)/double(3)*mu_x[idx_var_data] - mu_v_x[idx_var_data];
                        D_ptr[d_num_species*(d_num_species + 1) + 2][idx_diffusivities] = -mu_x[idx_var_data];
                        D_ptr[d_num_species*(d_num_species + 1) + 3][idx_diffusivities] = -u_x[idx_var_data]*(double(4)/double(3)*mu_x[idx_var_data] +
                            mu_v_x[idx_var_data]);
                        D_ptr[d_num_species*(d_num_species + 1) + 4][idx_diffusivities] = v_x[idx_var_data]*(double(2)/double(3)*mu_x[idx_var_data] -
                            mu_v_x[idx_var_data]);
                        D_ptr[d_num_species*(d_num_species + 1) + 5][idx_diffusivities] = -u_x[idx_var_data]*mu_x[idx_var_data];
                        D_ptr[d_num_species*(d_num_species + 1) + 6][idx_diffusivities] = -v_x[idx_var_data]*mu_x[idx_var_data];
                        D_ptr[d_num_species*(d_num_species + 1) + 7][idx_diffusivities] = -kappa_x[idx_var_data];
                    }
                }
                
                // Energy equation.
                for (int si = 0; si < d_num_species; si++)
                {
                    const int component_idx =
                        d_num_species*(d_num_species + 1) + 8 + si*(d_num_species + 1);
                    
                    for (int j = -num_ghosts[1]; j < interior_dims[1] + num_ghosts[1]; j++)
                    {
                        for (int i = -num_ghosts[0]; i < interior_dims[0] + num_ghosts[0] + 1; i++)
                        {
                            const int idx_diffusivities = (i + d_num_subghosts_diffusivities[0]) +
                                (j + d_num_subghosts_diffusivities[1])*(d_subghostcell_dims_diffusivities[0] + 1);
                            
                            const int idx_var_data = (i + num_ghosts[0]) +
                                (j + num_ghosts[1])*(ghostcell_dims[0] + 1);
                            
                            D_ptr[component_idx][idx_diffusivities] =
                                -rho_x[idx_var_data]*D_x[si][idx_var_data]*
                                    h_i_x[si][idx_var_data];
                        }
                    }
                }
                
                for (int si = 0; si < d_num_species; si++)
                {
                    for (int sj = 0; sj < d_num_species; sj++)
                    {
                        const int component_idx =
                            d_num_species*(d_num_species + 1) + 9 + si*(d_num_species + 1) + sj;
                        
                        for (int j = -num_ghosts[1]; j < interior_dims[1] + num_ghosts[1]; j++)
                        {
                            for (int i = -num_ghosts[0]; i < interior_dims[0] + num_ghosts[0] + 1; i++)
                            {
                                const int idx_diffusivities = (i + d_num_subghosts_diffusivities[0]) +
                                    (j + d_num_subghosts_diffusivities[1])*(d_subghostcell_dims_diffusivities[0] + 1);
                                
                                const int idx_var_data = (i + num_ghosts[0]) +
                                    (j + num_ghosts[1])*(ghostcell_dims[0] + 1);
                                
                                D_ptr[component_idx][idx_diffusivities] =
                                    rho_x[idx_var_data]*Y_x[si][idx_var_data]*D_x[sj][idx_var_data]*
                                        h_i_x[si][idx_var_data];
                            }
                        }
                    }
                }
                
                D_ptr.clear();
                
                /*
                 * Compute the diffusivities for the y-direction derivatives.
                 */
                
                for (int i = 0; i < 2*d_num_species*(d_num_species + 1) + 8; i++)
                {
                    D_ptr.push_back(d_side_data_diffusivities->getPointer(1, i));
                }
                
                // Mass equations.
                for (int si = 0; si < d_num_species; si++)
                {
                    const int component_idx = si*(d_num_species + 1);
                    
                    for (int j = -num_ghosts[1]; j < interior_dims[1] + num_ghosts[1] + 1; j++)
                    {
                        for (int i = -num_ghosts[0]; i < interior_dims[0] + num_ghosts[0]; i++)
                        {
                            const int idx_diffusivities = (i + d_num_subghosts_diffusivities[0]) +
                                (j + d_num_subghosts_diffusivities[1])*d_subghostcell_dims_diffusivities[0];
                            
                            const int idx_var_data = (i + num_ghosts[0]) +
                                (j + num_ghosts[1])*ghostcell_dims[0];;
                            
                            D_ptr[component_idx][idx_diffusivities] =
                                -rho_y[idx_var_data]*D_y[si][idx_var_data];
                        }
                    }
                }
                
                for (int si = 0; si < d_num_species; si++)
                {
                    for (int sj = 0; sj < d_num_species; sj++)
                    {
                        const int component_idx = si*(d_num_species + 1) + sj + 1;
                        
                        for (int j = -num_ghosts[1]; j < interior_dims[1] + num_ghosts[1] + 1; j++)
                        {
                            for (int i = -num_ghosts[0]; i < interior_dims[0] + num_ghosts[0]; i++)
                            {
                                const int idx_diffusivities = (i + d_num_subghosts_diffusivities[0]) +
                                    (j + d_num_subghosts_diffusivities[1])*d_subghostcell_dims_diffusivities[0];
                                
                                const int idx_var_data = (i + num_ghosts[0]) +
                                    (j + num_ghosts[1])*ghostcell_dims[0];;
                                
                                D_ptr[component_idx][idx_diffusivities] =
                                    rho_y[idx_var_data]*Y_y[si][idx_var_data]*D_y[sj][idx_var_data];
                            }
                        }
                    }
                }
                
                // Momentum and energy equations.
                for (int j = -num_ghosts[1]; j < interior_dims[1] + num_ghosts[1] + 1; j++)
                {
                    for (int i = -num_ghosts[0]; i < interior_dims[0] + num_ghosts[0]; i++)
                    {
                        const int idx_diffusivities = (i + d_num_subghosts_diffusivities[0]) +
                            (j + d_num_subghosts_diffusivities[1])*d_subghostcell_dims_diffusivities[0];
                        
                        const int idx_var_data = (i + num_ghosts[0]) +
                            (j + num_ghosts[1])*ghostcell_dims[0];
                        
                        D_ptr[d_num_species*(d_num_species + 1) + 0][idx_diffusivities] = -(double(4)/double(3)*mu_y[idx_var_data] + mu_v_y[idx_var_data]);
                        D_ptr[d_num_species*(d_num_species + 1) + 1][idx_diffusivities] = double(2)/double(3)*mu_y[idx_var_data] - mu_v_y[idx_var_data];
                        D_ptr[d_num_species*(d_num_species + 1) + 2][idx_diffusivities] = -mu_y[idx_var_data];
                        D_ptr[d_num_species*(d_num_species + 1) + 3][idx_diffusivities] = -v_y[idx_var_data]*(double(4)/double(3)*mu_y[idx_var_data] +
                            mu_v_y[idx_var_data]);
                        D_ptr[d_num_species*(d_num_species + 1) + 4][idx_diffusivities] = u_y[idx_var_data]*(double(2)/double(3)*mu_y[idx_var_data] -
                            mu_v_y[idx_var_data]);
                        D_ptr[d_num_species*(d_num_species + 1) + 5][idx_diffusivities] = -u_y[idx_var_data]*mu_y[idx_var_data];
                        D_ptr[d_num_species*(d_num_species + 1) + 6][idx_diffusivities] = -v_y[idx_var_data]*mu_y[idx_var_data];
                        D_ptr[d_num_species*(d_num_species + 1) + 7][idx_diffusivities] = -kappa_y[idx_var_data];
                    }
                }
                
                // Energy equation.
                for (int si = 0; si < d_num_species; si++)
                {
                    const int component_idx =
                        d_num_species*(d_num_species + 1) + 8 + si*(d_num_species + 1);
                    
                    for (int j = -num_ghosts[1]; j < interior_dims[1] + num_ghosts[1] + 1; j++)
                    {
                        for (int i = -num_ghosts[0]; i < interior_dims[0] + num_ghosts[0]; i++)
                        {
                            const int idx_diffusivities = (i + d_num_subghosts_diffusivities[0]) +
                                (j + d_num_subghosts_diffusivities[1])*d_subghostcell_dims_diffusivities[0];
                            
                            const int idx_var_data = (i + num_ghosts[0]) +
                                (j + num_ghosts[1])*ghostcell_dims[0];
                            
                            D_ptr[component_idx][idx_diffusivities] =
                                -rho_y[idx_var_data]*D_y[si][idx_var_data]*
                                    h_i_y[si][idx_var_data];
                        }
                    }
                }
                
                for (int si = 0; si < d_num_species; si++)
                {
                    for (int sj = 0; sj < d_num_species; sj++)
                    {
                        const int component_idx =
                            d_num_species*(d_num_species + 1) + 9 + si*(d_num_species + 1) + sj;
                        
                        for (int j = -num_ghosts[1]; j < interior_dims[1] + num_ghosts[1] + 1; j++)
                        {
                            for (int i = -num_ghosts[0]; i < interior_dims[0] + num_ghosts[0]; i++)
                            {
                                const int idx_diffusivities = (i + d_num_subghosts_diffusivities[0]) +
                                    (j + d_num_subghosts_diffusivities[1])*d_subghostcell_dims_diffusivities[0];
                                
                                const int idx_var_data = (i + num_ghosts[0]) +
                                    (j + num_ghosts[1])*ghostcell_dims[0];
                                
                                D_ptr[component_idx][idx_diffusivities] =
                                    rho_y[idx_var_data]*Y_y[si][idx_var_data]*D_y[sj][idx_var_data]*
                                        h_i_y[si][idx_var_data];
                            }
                        }
                    }
                }
                
                D_ptr.clear();
            }
            else if (d_dim == tbox::Dimension(3))
            {
                std::vector<double*> D_x;
                std::vector<double*> D_y;
                std::vector<double*> D_z;
                D_x.reserve(d_num_species);
                D_y.reserve(d_num_species);
                D_z.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    D_x.push_back(var_data_for_diffusivities[si]->getPointer(0, 0));
                    D_y.push_back(var_data_for_diffusivities[si]->getPointer(1, 0));
                    D_z.push_back(var_data_for_diffusivities[si]->getPointer(2, 0));
                }
                
                double* mu_x    = var_data_for_diffusivities[d_num_species    ]->getPointer(0, 0);
                double* mu_v_x  = var_data_for_diffusivities[d_num_species + 1]->getPointer(0, 0);
                double* kappa_x = var_data_for_diffusivities[d_num_species + 2]->getPointer(0, 0);
                
                double* mu_y    = var_data_for_diffusivities[d_num_species    ]->getPointer(1, 0);
                double* mu_v_y  = var_data_for_diffusivities[d_num_species + 1]->getPointer(1, 0);
                double* kappa_y = var_data_for_diffusivities[d_num_species + 2]->getPointer(1, 0);
                
                double* mu_z    = var_data_for_diffusivities[d_num_species    ]->getPointer(2, 0);
                double* mu_v_z  = var_data_for_diffusivities[d_num_species + 1]->getPointer(2, 0);
                double* kappa_z = var_data_for_diffusivities[d_num_species + 2]->getPointer(2, 0);
                
                double* rho_x = var_data_for_diffusivities[d_num_species + 3]->getPointer(0, 0);
                double* rho_y = var_data_for_diffusivities[d_num_species + 3]->getPointer(1, 0);
                double* rho_z = var_data_for_diffusivities[d_num_species + 3]->getPointer(2, 0);
                
                std::vector<double*> Y_x;
                std::vector<double*> Y_y;
                std::vector<double*> Y_z;
                Y_x.reserve(d_num_species);
                Y_y.reserve(d_num_species);
                Y_z.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    Y_x.push_back(var_data_for_diffusivities[d_num_species + 4 + si]->getPointer(0, 0));
                    Y_y.push_back(var_data_for_diffusivities[d_num_species + 4 + si]->getPointer(1, 0));
                    Y_z.push_back(var_data_for_diffusivities[d_num_species + 4 + si]->getPointer(2, 0));
                }
                
                double* u_x = var_data_for_diffusivities[2*d_num_species + 4]->getPointer(0, 0);
                double* v_x = var_data_for_diffusivities[2*d_num_species + 5]->getPointer(0, 0);
                double* w_x = var_data_for_diffusivities[2*d_num_species + 6]->getPointer(0, 0);
                
                double* u_y = var_data_for_diffusivities[2*d_num_species + 4]->getPointer(1, 0);
                double* v_y = var_data_for_diffusivities[2*d_num_species + 5]->getPointer(1, 0);
                double* w_y = var_data_for_diffusivities[2*d_num_species + 6]->getPointer(1, 0);
                
                double* u_z = var_data_for_diffusivities[2*d_num_species + 4]->getPointer(2, 0);
                double* v_z = var_data_for_diffusivities[2*d_num_species + 5]->getPointer(2, 0);
                double* w_z = var_data_for_diffusivities[2*d_num_species + 6]->getPointer(2, 0);
                
                std::vector<double*> h_i_x;
                std::vector<double*> h_i_y;
                std::vector<double*> h_i_z;
                h_i_x.reserve(d_num_species);
                h_i_y.reserve(d_num_species);
                h_i_z.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    h_i_x.push_back(var_data_for_diffusivities[2*d_num_species + 7 + si]->getPointer(0, 0));
                    h_i_y.push_back(var_data_for_diffusivities[2*d_num_species + 7 + si]->getPointer(1, 0));
                    h_i_z.push_back(var_data_for_diffusivities[2*d_num_species + 7 + si]->getPointer(2, 0));
                }
                
                
                std::vector<double*> D_ptr;
                D_ptr.reserve(2*d_num_species*(d_num_species + 1) + 10);
                
                /*
                 * Compute the diffusivities for the x-direction derivatives.
                 */
                
                for (int i = 0; i < 2*d_num_species*(d_num_species + 1) + 10; i++)
                {
                    D_ptr.push_back(d_side_data_diffusivities->getPointer(0, i));
                }
                
                // Mass equations.
                for (int si = 0; si < d_num_species; si++)
                {
                    const int component_idx = si*(d_num_species + 1);
                    
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
                                
                                D_ptr[component_idx][idx_diffusivities] =
                                    -rho_x[idx_var_data]*D_x[si][idx_var_data];
                            }
                        }
                    }
                }
                
                for (int si = 0; si < d_num_species; si++)
                {
                    for (int sj = 0; sj < d_num_species; sj++)
                    {
                        const int component_idx = si*(d_num_species + 1) + sj + 1;
                        
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
                                    
                                    D_ptr[component_idx][idx_diffusivities] =
                                        rho_x[idx_var_data]*Y_x[si][idx_var_data]*D_x[sj][idx_var_data];
                                }
                            }
                        }
                    }
                }
                
                // Momentum and energy equations.
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
                            
                            D_ptr[d_num_species*(d_num_species + 1) + 0][idx_diffusivities] = -(double(4)/double(3)*mu_x[idx_var_data] + mu_v_x[idx_var_data]);
                            D_ptr[d_num_species*(d_num_species + 1) + 1][idx_diffusivities] = double(2)/double(3)*mu_x[idx_var_data] - mu_v_x[idx_var_data];
                            D_ptr[d_num_species*(d_num_species + 1) + 2][idx_diffusivities] = -mu_x[idx_var_data];
                            D_ptr[d_num_species*(d_num_species + 1) + 3][idx_diffusivities] = -u_x[idx_var_data]*(double(4)/double(3)*mu_x[idx_var_data] +
                                mu_v_x[idx_var_data]);
                            D_ptr[d_num_species*(d_num_species + 1) + 4][idx_diffusivities] = v_x[idx_var_data]*(double(2)/double(3)*mu_x[idx_var_data] -
                                mu_v_x[idx_var_data]);
                            D_ptr[d_num_species*(d_num_species + 1) + 5][idx_diffusivities] = w_x[idx_var_data]*(double(2)/double(3)*mu_x[idx_var_data] -
                                mu_v_x[idx_var_data]);
                            D_ptr[d_num_species*(d_num_species + 1) + 6][idx_diffusivities] = -u_x[idx_var_data]*mu_x[idx_var_data];
                            D_ptr[d_num_species*(d_num_species + 1) + 7][idx_diffusivities] = -v_x[idx_var_data]*mu_x[idx_var_data];
                            D_ptr[d_num_species*(d_num_species + 1) + 8][idx_diffusivities] = -w_x[idx_var_data]*mu_x[idx_var_data];
                            D_ptr[d_num_species*(d_num_species + 1) + 9][idx_diffusivities] = -kappa_x[idx_var_data];
                        }
                    }
                }
                
                // Energy equation.
                for (int si = 0; si < d_num_species; si++)
                {
                    const int component_idx =
                        d_num_species*(d_num_species + 1) + 10 + si*(d_num_species + 1);
                    
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
                                
                                D_ptr[component_idx][idx_diffusivities] = 
                                    -rho_x[idx_var_data]*D_x[si][idx_var_data]*
                                        h_i_x[si][idx_var_data];
                            }
                        }
                    }
                }
                
                for (int si = 0; si < d_num_species; si++)
                {
                    for (int sj = 0; sj < d_num_species; sj++)
                    {
                        const int component_idx =
                            d_num_species*(d_num_species + 1) + 11 + si*(d_num_species + 1) + sj;
                        
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
                                    
                                    D_ptr[component_idx][idx_diffusivities] =
                                        rho_x[idx_var_data]*Y_x[si][idx_var_data]*D_x[sj][idx_var_data]*
                                            h_i_x[si][idx_var_data];
                                }
                            }
                        }
                    }
                }
                
                D_ptr.clear();
                
                /*
                 * Compute the diffusivities for the y-direction derivatives.
                 */
                
                for (int i = 0; i < 2*d_num_species*(d_num_species + 1) + 10; i++)
                {
                    D_ptr.push_back(d_side_data_diffusivities->getPointer(1, i));
                }
                
                // Mass equations.
                for (int si = 0; si < d_num_species; si++)
                {
                    const int component_idx = si*(d_num_species + 1);
                    
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
                                
                                D_ptr[component_idx][idx_diffusivities] =
                                    -rho_y[idx_var_data]*D_y[si][idx_var_data];
                            }
                        }
                    }
                }
                
                for (int si = 0; si < d_num_species; si++)
                {
                    for (int sj = 0; sj < d_num_species; sj++)
                    {
                        const int component_idx = si*(d_num_species + 1) + sj + 1;
                        
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
                                    
                                    D_ptr[component_idx][idx_diffusivities] =
                                        rho_y[idx_var_data]*Y_y[si][idx_var_data]*D_y[sj][idx_var_data];
                                }
                            }
                        }
                    }
                }
                
                // Momentum and energy equations.
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
                            
                            D_ptr[d_num_species*(d_num_species + 1) + 0][idx_diffusivities] = -(double(4)/double(3)*mu_y[idx_var_data] + mu_v_y[idx_var_data]);
                            D_ptr[d_num_species*(d_num_species + 1) + 1][idx_diffusivities] = double(2)/double(3)*mu_y[idx_var_data] - mu_v_y[idx_var_data];
                            D_ptr[d_num_species*(d_num_species + 1) + 2][idx_diffusivities] = -mu_y[idx_var_data];
                            D_ptr[d_num_species*(d_num_species + 1) + 3][idx_diffusivities] = -v_y[idx_var_data]*(double(4)/double(3)*mu_y[idx_var_data] +
                                mu_v_y[idx_var_data]);
                            D_ptr[d_num_species*(d_num_species + 1) + 4][idx_diffusivities] = u_y[idx_var_data]*(double(2)/double(3)*mu_y[idx_var_data] -
                                mu_v_y[idx_var_data]);
                            D_ptr[d_num_species*(d_num_species + 1) + 5][idx_diffusivities] = w_y[idx_var_data]*(double(2)/double(3)*mu_y[idx_var_data] -
                                mu_v_y[idx_var_data]);
                            D_ptr[d_num_species*(d_num_species + 1) + 6][idx_diffusivities] = -u_y[idx_var_data]*mu_y[idx_var_data];
                            D_ptr[d_num_species*(d_num_species + 1) + 7][idx_diffusivities] = -v_y[idx_var_data]*mu_y[idx_var_data];
                            D_ptr[d_num_species*(d_num_species + 1) + 8][idx_diffusivities] = -w_y[idx_var_data]*mu_y[idx_var_data];
                            D_ptr[d_num_species*(d_num_species + 1) + 9][idx_diffusivities] = -kappa_y[idx_var_data];
                        }
                    }
                }
                
                // Energy equation.
                for (int si = 0; si < d_num_species; si++)
                {
                    const int component_idx =
                        d_num_species*(d_num_species + 1) + 10 + si*(d_num_species + 1);
                    
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
                                
                                D_ptr[component_idx][idx_diffusivities] = 
                                    -rho_y[idx_var_data]*D_y[si][idx_var_data]*
                                        h_i_y[si][idx_var_data];
                            }
                        }
                    }
                }
                
                for (int si = 0; si < d_num_species; si++)
                {
                    for (int sj = 0; sj < d_num_species; sj++)
                    {
                        const int component_idx =
                            d_num_species*(d_num_species + 1) + 11 + si*(d_num_species + 1) + sj;
                        
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
                                    
                                    D_ptr[component_idx][idx_diffusivities] =
                                        rho_y[idx_var_data]*Y_y[si][idx_var_data]*D_y[sj][idx_var_data]*
                                            h_i_y[si][idx_var_data];
                                }
                            }
                        }
                    }
                }
                
                D_ptr.clear();
                
                /*
                 * Compute the diffusivities for the z-direction derivatives.
                 */
                
                for (int i = 0; i < 2*d_num_species*(d_num_species + 1) + 10; i++)
                {
                    D_ptr.push_back(d_side_data_diffusivities->getPointer(2, i));
                }
                
                // Mass equations.
                for (int si = 0; si < d_num_species; si++)
                {
                    const int component_idx = si*(d_num_species + 1);
                    
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
                                
                                D_ptr[component_idx][idx_diffusivities] =
                                    -rho_z[idx_var_data]*D_z[si][idx_var_data];
                            }
                        }
                    }
                }
                
                for (int si = 0; si < d_num_species; si++)
                {
                    for (int sj = 0; sj < d_num_species; sj++)
                    {
                        const int component_idx = si*(d_num_species + 1) + sj + 1;
                        
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
                                    
                                    D_ptr[component_idx][idx_diffusivities] =
                                        rho_z[idx_var_data]*Y_z[si][idx_var_data]*D_z[sj][idx_var_data];
                                }
                            }
                        }
                    }
                }
                
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
                            
                            D_ptr[d_num_species*(d_num_species + 1) + 0][idx_diffusivities] = -(double(4)/double(3)*mu_z[idx_var_data] + mu_v_z[idx_var_data]);
                            D_ptr[d_num_species*(d_num_species + 1) + 1][idx_diffusivities] = double(2)/double(3)*mu_z[idx_var_data] - mu_v_z[idx_var_data];
                            D_ptr[d_num_species*(d_num_species + 1) + 2][idx_diffusivities] = -mu_z[idx_var_data];
                            D_ptr[d_num_species*(d_num_species + 1) + 3][idx_diffusivities] = -w_z[idx_var_data]*(double(4)/double(3)*mu_z[idx_var_data] +
                                mu_v_z[idx_var_data]);
                            D_ptr[d_num_species*(d_num_species + 1) + 4][idx_diffusivities] = u_z[idx_var_data]*(double(2)/double(3)*mu_z[idx_var_data] -
                                mu_v_z[idx_var_data]);
                            D_ptr[d_num_species*(d_num_species + 1) + 5][idx_diffusivities] = v_z[idx_var_data]*(double(2)/double(3)*mu_z[idx_var_data] -
                                mu_v_z[idx_var_data]);
                            D_ptr[d_num_species*(d_num_species + 1) + 6][idx_diffusivities] = -u_z[idx_var_data]*mu_z[idx_var_data];
                            D_ptr[d_num_species*(d_num_species + 1) + 7][idx_diffusivities] = -v_z[idx_var_data]*mu_z[idx_var_data];
                            D_ptr[d_num_species*(d_num_species + 1) + 8][idx_diffusivities] = -w_z[idx_var_data]*mu_z[idx_var_data];
                            D_ptr[d_num_species*(d_num_species + 1) + 9][idx_diffusivities] = -kappa_z[idx_var_data];
                        }
                    }
                }
                
                // Energy equation.
                for (int si = 0; si < d_num_species; si++)
                {
                    const int component_idx =
                        d_num_species*(d_num_species + 1) + 10 + si*(d_num_species + 1);
                    
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
                                
                                D_ptr[component_idx][idx_diffusivities] = 
                                    -rho_z[idx_var_data]*D_z[si][idx_var_data]*
                                        h_i_z[si][idx_var_data];
                            }
                        }
                    }
                }
                
                for (int si = 0; si < d_num_species; si++)
                {
                    for (int sj = 0; sj < d_num_species; sj++)
                    {
                        const int component_idx =
                            d_num_species*(d_num_species + 1) + 11 + si*(d_num_species + 1) + sj;
                        
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
                                    
                                    D_ptr[component_idx][idx_diffusivities] =
                                        rho_z[idx_var_data]*Y_z[si][idx_var_data]*D_z[sj][idx_var_data]*
                                            h_i_z[si][idx_var_data];
                                }
                            }
                        }
                    }
                }
                
                D_ptr.clear();
            }
            
            d_side_data_diffusivities_computed = true;
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::computeSideDataOfDiffusiveFluxDiffusivities()\n"
            << "Cell data of 'DIFFUSIVITIES' is not yet registered."
            << std::endl);
    }
}


/*
 * Get the side data of the diffusivities in the diffusive fluxa.
 */
void
FlowModelDiffusiveFluxUtilitiesFourEqnConservative::getSideDataOfDiffusiveFluxDiffusivities(
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
            << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::"
            << "getSideDataOfDiffusiveFluxDiffusivities()\n"
            << "No patch is registered yet."
            << std::endl);
    }
    
    if (!d_cell_data_computed_diffusivities)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::getSideDataOfDiffusiveFluxDiffusivities()\n"
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
                            << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::"
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
                    << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::getSideDataOfDiffusiveFluxDiffusivities()\n"
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
                            << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::"
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
                            << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::"
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
                    << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::getSideDataOfDiffusiveFluxDiffusivities()\n"
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
                            << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::"
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
                            << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::"
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
                            << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::"
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
                    << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::getSideDataOfDiffusiveFluxDiffusivities()\n"
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
FlowModelDiffusiveFluxUtilitiesFourEqnConservative::setNumberOfSubGhosts(
    const hier::IntVector& num_subghosts,
    const std::string& variable_name,
    const std::string& parent_variable_name)
{
    NULL_USE(parent_variable_name);
    
    if (variable_name == "MASS_DIFFUSIVITIES")
    {
        if (d_num_subghosts_mass_diffusivities > -hier::IntVector::getOne(d_dim))
        {
            if (num_subghosts > d_num_subghosts_mass_diffusivities)
            {
                d_num_subghosts_mass_diffusivities = num_subghosts;
            }
        }
        else
        {
            d_num_subghosts_mass_diffusivities = num_subghosts;
        }
    }
    
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
        
        setNumberOfSubGhosts(num_subghosts, "MASS_DIFFUSIVITIES", parent_variable_name);
        setNumberOfSubGhosts(num_subghosts, "SHEAR_VISCOSITY", parent_variable_name);
        setNumberOfSubGhosts(num_subghosts, "BULK_VISCOSITY", parent_variable_name);
        setNumberOfSubGhosts(num_subghosts, "THERMAL_CONDUCTIVITY", parent_variable_name);
    }
}


/*
 * Set the ghost boxes of derived cell variables.
 */
void
FlowModelDiffusiveFluxUtilitiesFourEqnConservative::setDerivedCellVariableGhostBoxes()
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
    
    if (d_num_subghosts_mass_diffusivities > -hier::IntVector::getOne(d_dim))
    {
        d_subghost_box_mass_diffusivities = interior_box;
        d_subghost_box_mass_diffusivities.grow(d_num_subghosts_mass_diffusivities);
        d_subghostcell_dims_mass_diffusivities = d_subghost_box_mass_diffusivities.numberCells();
    }
    
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
 * Compute the cell data of mass diffusivities in the registered patch.
 */
void
FlowModelDiffusiveFluxUtilitiesFourEqnConservative::computeCellDataOfMassDiffusivities()
{
    // Create empty box.
    const hier::Box empty_box(d_dim);
    
    if (d_num_subghosts_mass_diffusivities > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_mass_diffusivities)
        {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(d_data_mass_diffusivities);
#endif
            
            if (d_flow_model.expired())
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "The object is not setup yet!"
                    << std::endl);
            }
            
            HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
            
            // Get the cell data of mass fractions.
            HAMERS_SHARED_PTR<pdat::CellData<double> > data_mass_fractions =
                flow_model_tmp->getCellData("MASS_FRACTIONS");
            
            // Get the cell data of pressure.
            HAMERS_SHARED_PTR<pdat::CellData<double> > data_pressure =
                flow_model_tmp->getCellData("PRESSURE");
            
            // Get the cell data of temperature.
            HAMERS_SHARED_PTR<pdat::CellData<double> > data_temperature =
                flow_model_tmp->getCellData("TEMPERATURE");
            
            // Compute the mass diffusivity fields.
            d_equation_of_mass_diffusivity_mixing_rules->computeMassDiffusivities(
                d_data_mass_diffusivities,
                data_pressure,
                data_temperature,
                data_mass_fractions,
                empty_box);
            
            d_cell_data_computed_mass_diffusivities = true;
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::computeCellDataOfMassDiffusivities()\n"
            << "Cell data of 'MASS_DIFFUSIVITIES' is not yet registered."
            << std::endl);
    }
}


/*
 * Compute the cell data of shear viscosity in the registered patch.
 */
void
FlowModelDiffusiveFluxUtilitiesFourEqnConservative::computeCellDataOfShearViscosity()
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
            
            // Get the cell data of mass fractions.
            HAMERS_SHARED_PTR<pdat::CellData<double> > data_mass_fractions =
                flow_model_tmp->getCellData("MASS_FRACTIONS");
            
            // Get the cell data of pressure.
            HAMERS_SHARED_PTR<pdat::CellData<double> > data_pressure =
                flow_model_tmp->getCellData("PRESSURE");
            
            // Get the cell data of temperature.
            HAMERS_SHARED_PTR<pdat::CellData<double> > data_temperature =
                flow_model_tmp->getCellData("TEMPERATURE");
            
            // Compute the shear viscosity field.
            d_equation_of_shear_viscosity_mixing_rules->computeShearViscosity(
                d_data_shear_viscosity,
                data_pressure,
                data_temperature,
                data_mass_fractions,
                empty_box);
            
            d_cell_data_computed_shear_viscosity = true;
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::computeCellDataOfShearViscosity()\n"
            << "Cell data of 'SHEAR_VISCOSITY' is not yet registered."
            << std::endl);
    }
}


/*
 * Compute the cell data of bulk viscosity in the registered patch.
 */
void
FlowModelDiffusiveFluxUtilitiesFourEqnConservative::computeCellDataOfBulkViscosity()
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
            
            // Get the cell data of mass fractions.
            HAMERS_SHARED_PTR<pdat::CellData<double> > data_mass_fractions =
                flow_model_tmp->getCellData("MASS_FRACTIONS");
            
            // Get the cell data of pressure.
            HAMERS_SHARED_PTR<pdat::CellData<double> > data_pressure =
                flow_model_tmp->getCellData("PRESSURE");
            
            // Get the cell data of temperature.
            HAMERS_SHARED_PTR<pdat::CellData<double> > data_temperature =
                flow_model_tmp->getCellData("TEMPERATURE");
            
            // Compute the bulk viscosity field.
            d_equation_of_bulk_viscosity_mixing_rules->computeBulkViscosity(
                d_data_bulk_viscosity,
                data_pressure,
                data_temperature,
                data_mass_fractions,
                empty_box);
            
            d_cell_data_computed_bulk_viscosity = true;
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::computeCellDataOfBulkViscosity()\n"
            << "Cell data of 'BULK_VISCOSITY' is not yet registered."
            << std::endl);
    }
}


/*
 * Compute the cell data of thermal conductivity in the registered patch.
 */
void
FlowModelDiffusiveFluxUtilitiesFourEqnConservative::computeCellDataOfThermalConductivity()
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
            
            // Get the cell data of mass fractions.
            HAMERS_SHARED_PTR<pdat::CellData<double> > data_mass_fractions =
                flow_model_tmp->getCellData("MASS_FRACTIONS");
            
            // Get the cell data of pressure.
            HAMERS_SHARED_PTR<pdat::CellData<double> > data_pressure =
                flow_model_tmp->getCellData("PRESSURE");
            
            // Get the cell data of temperature.
            HAMERS_SHARED_PTR<pdat::CellData<double> > data_temperature =
                flow_model_tmp->getCellData("TEMPERATURE");
            
            // Compute the thermal conductivity field.
            d_equation_of_thermal_conductivity_mixing_rules->computeThermalConductivity(
                d_data_thermal_conductivity,
                data_pressure,
                data_temperature,
                data_mass_fractions,
                empty_box);
            
            d_cell_data_computed_thermal_conductivity = true;
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::computeCellDataOfThermalConductivity()\n"
            << "Cell data of 'THERMAL_CONDUCTIVITY' is not yet registered."
            << std::endl);
    }
}


/*
 * Compute the cell data of diffusivities in the registered patch.
 */
void
FlowModelDiffusiveFluxUtilitiesFourEqnConservative::computeCellDataOfDiffusivities()
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
            
            if (!d_cell_data_computed_mass_diffusivities)
            {
                computeCellDataOfMassDiffusivities();
            }
            
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
            
            // Get the cell data of density.
            HAMERS_SHARED_PTR<pdat::CellData<double> > data_density =
                flow_model_tmp->getCellData("DENSITY");
            
            // Get the cell data of mass fractions.
            HAMERS_SHARED_PTR<pdat::CellData<double> > data_mass_fractions =
                flow_model_tmp->getCellData("MASS_FRACTIONS");
            
            // Get the cell data of velocity.
            HAMERS_SHARED_PTR<pdat::CellData<double> > data_velocity =
                flow_model_tmp->getCellData("VELOCITY");
            
            // Get the cell data of species enthalpies.
            std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > data_species_enthalpies =
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
            
            const hier::Box subghost_box_density = data_density->getGhostBox();
            const hier::IntVector subghostcell_dims_density = subghost_box_density.numberCells();
            
            const hier::Box subghost_box_mass_fractions = data_mass_fractions->getGhostBox();
            const hier::IntVector subghostcell_dims_mass_fractions = subghost_box_mass_fractions.numberCells();
            
            const hier::Box subghost_box_velocity = data_velocity->getGhostBox();
            const hier::IntVector subghostcell_dims_velocity = subghost_box_velocity.numberCells();
            
            const hier::Box subghost_box_species_enthalpies = data_species_enthalpies[0]->getGhostBox();
            const hier::IntVector subghostcell_dims_species_enthalpies = subghost_box_species_enthalpies.numberCells();
            
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
            
            if (d_dim == tbox::Dimension(1))
            {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(d_data_diffusivities->getDepth() == 2*d_num_species*(d_num_species + 1) + 3);
#endif
                
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
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(d_data_diffusivities->getDepth() == 2*d_num_species*(d_num_species + 1) + 10);
#endif
                
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
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(d_data_diffusivities->getDepth() == 2*d_num_species*(d_num_species + 1) + 13);
#endif
                
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
            
            d_cell_data_computed_diffusivities = true;
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelDiffusiveFluxUtilitiesFourEqnConservative::computeCellDataOfDiffusivities()\n"
            << "Cell data of 'DIFFUSIVITIES' is not yet registered."
            << std::endl);
    }
}
