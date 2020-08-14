#include "flow/flow_models/four-eqn_conservative/FlowModelStatisticsUtilitiesFourEqnConservative.hpp"

#include "extn/patch_hierarchies/ExtendedFlattenedHierarchy.hpp"

#include <fstream>

boost::shared_ptr<pdat::CellVariable<double> > FlowModelStatisticsUtilitiesFourEqnConservative::s_variable_partial_density_unfiltered;
boost::shared_ptr<pdat::CellVariable<double> > FlowModelStatisticsUtilitiesFourEqnConservative::s_variable_momentum_unfiltered;
boost::shared_ptr<pdat::CellVariable<double> > FlowModelStatisticsUtilitiesFourEqnConservative::s_variable_total_energy_unfiltered;

boost::shared_ptr<pdat::CellVariable<double> > FlowModelStatisticsUtilitiesFourEqnConservative::s_variable_partial_density_filtered;
boost::shared_ptr<pdat::CellVariable<double> > FlowModelStatisticsUtilitiesFourEqnConservative::s_variable_momentum_filtered;
boost::shared_ptr<pdat::CellVariable<double> > FlowModelStatisticsUtilitiesFourEqnConservative::s_variable_total_energy_filtered;

boost::shared_ptr<pdat::CellVariable<double> > FlowModelStatisticsUtilitiesFourEqnConservative::s_variable_pressure_unfiltered;
boost::shared_ptr<pdat::CellVariable<double> > FlowModelStatisticsUtilitiesFourEqnConservative::s_variable_shear_stress_unfiltered;
boost::shared_ptr<pdat::CellVariable<double> > FlowModelStatisticsUtilitiesFourEqnConservative::s_variable_convective_stress_unfiltered;

boost::shared_ptr<pdat::CellVariable<double> > FlowModelStatisticsUtilitiesFourEqnConservative::s_variable_pressure_filtered;
boost::shared_ptr<pdat::CellVariable<double> > FlowModelStatisticsUtilitiesFourEqnConservative::s_variable_shear_stress_filtered;

boost::shared_ptr<pdat::CellVariable<double> > FlowModelStatisticsUtilitiesFourEqnConservative::s_variable_convective_stress_filtered;
boost::shared_ptr<pdat::CellVariable<double> > FlowModelStatisticsUtilitiesFourEqnConservative::s_variable_SFS_stress;

boost::shared_ptr<pdat::CellVariable<double> > FlowModelStatisticsUtilitiesFourEqnConservative::s_variable_velocity_Favre_filtered;
boost::shared_ptr<pdat::CellVariable<double> > FlowModelStatisticsUtilitiesFourEqnConservative::s_variable_specific_volume_Favre_filtered;
boost::shared_ptr<pdat::CellVariable<double> > FlowModelStatisticsUtilitiesFourEqnConservative::s_variable_density_filtered;

FlowModelStatisticsUtilitiesFourEqnConservative::FlowModelStatisticsUtilitiesFourEqnConservative(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
    const int& num_species,
    const boost::shared_ptr<tbox::Database>& flow_model_db,
    const boost::shared_ptr<EquationOfStateMixingRules> equation_of_state_mixing_rules,
    const boost::shared_ptr<EquationOfMassDiffusivityMixingRules> equation_of_mass_diffusivity_mixing_rules,
    const boost::shared_ptr<EquationOfShearViscosityMixingRules> equation_of_shear_viscosity_mixing_rules,
    const boost::shared_ptr<EquationOfBulkViscosityMixingRules> equation_of_bulk_viscosity_mixing_rules,
    const boost::shared_ptr<EquationOfThermalConductivityMixingRules> equation_of_thermal_conductivity_mixing_rules):
        FlowModelStatisticsUtilities(
            object_name,
            dim,
            grid_geometry,
            num_species,
            flow_model_db),
        d_pressure_filtered(false),
        d_shear_stress_filtered(false),
        d_convective_stress_filtered(false),
        d_derivatives_filtered(false),
        d_conservative_variables_filtered(false),
        d_equation_of_state_mixing_rules(equation_of_state_mixing_rules),
        d_equation_of_mass_diffusivity_mixing_rules(equation_of_mass_diffusivity_mixing_rules),
        d_equation_of_shear_viscosity_mixing_rules(equation_of_shear_viscosity_mixing_rules),
        d_equation_of_bulk_viscosity_mixing_rules(equation_of_bulk_viscosity_mixing_rules),
        d_equation_of_thermal_conductivity_mixing_rules(equation_of_thermal_conductivity_mixing_rules)
{
    /*
     * Initialize the variables.
     */
    
    s_variable_partial_density_unfiltered = boost::shared_ptr<pdat::CellVariable<double> > (
        new pdat::CellVariable<double>(d_dim, "partial density unfiltered", d_num_species));
    
    s_variable_momentum_unfiltered = boost::shared_ptr<pdat::CellVariable<double> > (
        new pdat::CellVariable<double>(d_dim, "momentum unfiltered", d_dim.getValue()));
    
    s_variable_total_energy_unfiltered = boost::shared_ptr<pdat::CellVariable<double> > (
        new pdat::CellVariable<double>(d_dim, "total energy unfiltered", 1));
    
    s_variable_partial_density_filtered = boost::shared_ptr<pdat::CellVariable<double> > (
        new pdat::CellVariable<double>(d_dim, "partial density filtered", d_num_species));
    
    s_variable_momentum_filtered = boost::shared_ptr<pdat::CellVariable<double> > (
        new pdat::CellVariable<double>(d_dim, "momentum filtered", d_dim.getValue()));
    
    s_variable_total_energy_filtered = boost::shared_ptr<pdat::CellVariable<double> > (
        new pdat::CellVariable<double>(d_dim, "total energy filtered", 1));
    
    
    s_variable_pressure_unfiltered = boost::shared_ptr<pdat::CellVariable<double> > (
        new pdat::CellVariable<double>(d_dim, "pressure unfiltered", 1));
    
    if (d_dim == tbox::Dimension(1))
    {
        s_variable_shear_stress_unfiltered = boost::shared_ptr<pdat::CellVariable<double> > (
            new pdat::CellVariable<double>(d_dim, "shear stress unfiltered", 1));
    }
    else if (d_dim == tbox::Dimension(2))
    {
        s_variable_shear_stress_unfiltered = boost::shared_ptr<pdat::CellVariable<double> > (
            new pdat::CellVariable<double>(d_dim, "shear stress unfiltered", 3));
    }
    else if (d_dim == tbox::Dimension(3))
    {
        s_variable_shear_stress_unfiltered = boost::shared_ptr<pdat::CellVariable<double> > (
            new pdat::CellVariable<double>(d_dim, "shear stress unfiltered", 6));
    }
    
    if (d_dim == tbox::Dimension(1))
    {
        s_variable_convective_stress_unfiltered = boost::shared_ptr<pdat::CellVariable<double> > (
            new pdat::CellVariable<double>(d_dim, "convective stress unfiltered", 1));
    }
    else if (d_dim == tbox::Dimension(2))
    {
        s_variable_convective_stress_unfiltered = boost::shared_ptr<pdat::CellVariable<double> > (
            new pdat::CellVariable<double>(d_dim, "convective stress unfiltered", 3));
    }
    else if (d_dim == tbox::Dimension(3))
    {
        s_variable_convective_stress_unfiltered = boost::shared_ptr<pdat::CellVariable<double> > (
            new pdat::CellVariable<double>(d_dim, "convective stress unfiltered", 6));
    }
    
    s_variable_pressure_filtered = boost::shared_ptr<pdat::CellVariable<double> > (
        new pdat::CellVariable<double>(d_dim, "pressure filtered", 1));
    
    if (d_dim == tbox::Dimension(1))
    {
        s_variable_shear_stress_filtered = boost::shared_ptr<pdat::CellVariable<double> > (
            new pdat::CellVariable<double>(d_dim, "shear stress filtered", 1));
    }
    else if (d_dim == tbox::Dimension(2))
    {
        s_variable_shear_stress_filtered = boost::shared_ptr<pdat::CellVariable<double> > (
            new pdat::CellVariable<double>(d_dim, "shear stress filtered", 3));
    }
    else if (d_dim == tbox::Dimension(3))
    {
        s_variable_shear_stress_filtered = boost::shared_ptr<pdat::CellVariable<double> > (
            new pdat::CellVariable<double>(d_dim, "shear stress filtered", 6));
    }
    
    if (d_dim == tbox::Dimension(1))
    {
        s_variable_convective_stress_filtered = boost::shared_ptr<pdat::CellVariable<double> > (
            new pdat::CellVariable<double>(d_dim, "convective stress filtered", 1));
        
        s_variable_SFS_stress = boost::shared_ptr<pdat::CellVariable<double> > (
            new pdat::CellVariable<double>(d_dim, "SFS stress", 1));
    }
    else if (d_dim == tbox::Dimension(2))
    {
        s_variable_convective_stress_filtered = boost::shared_ptr<pdat::CellVariable<double> > (
            new pdat::CellVariable<double>(d_dim, "convective stress filtered", 3));
        
        s_variable_SFS_stress = boost::shared_ptr<pdat::CellVariable<double> > (
            new pdat::CellVariable<double>(d_dim, "SFS stress", 3));
    }
    else if (d_dim == tbox::Dimension(3))
    {
        s_variable_convective_stress_filtered = boost::shared_ptr<pdat::CellVariable<double> > (
            new pdat::CellVariable<double>(d_dim, "convective stress filtered", 6));
        
        s_variable_SFS_stress = boost::shared_ptr<pdat::CellVariable<double> > (
            new pdat::CellVariable<double>(d_dim, "SFS stress", 6));
    }
    
    s_variable_velocity_Favre_filtered = boost::shared_ptr<pdat::CellVariable<double> > (
        new pdat::CellVariable<double>(d_dim, "velocity Favre-filtered", d_dim.getValue()));
    
    s_variable_specific_volume_Favre_filtered = boost::shared_ptr<pdat::CellVariable<double> > (
        new pdat::CellVariable<double>(d_dim, "specific volume Favre-filtered", 1));
    
    s_variable_density_filtered = boost::shared_ptr<pdat::CellVariable<double> > (
        new pdat::CellVariable<double>(d_dim, "density filtered", 1));
    
    /*
     * Initialize the filters.
     */
    
    // For debugging.
    // d_filter_x = boost::shared_ptr<FilterNone> (
    //     new FilterNone("d_filter_x", d_dim, DIRECTION::X_DIRECTION));
    
    // if ((d_dim == tbox::Dimension(2)) || (d_dim == tbox::Dimension(3)))
    // {
    //     d_filter_y = boost::shared_ptr<FilterNone> (
    //         new FilterNone("d_filter_y", d_dim, DIRECTION::Y_DIRECTION));
    // }
    
    // if (d_dim == tbox::Dimension(3))
    // {
    //     d_filter_z = boost::shared_ptr<FilterNone> (
    //         new FilterNone("d_filter_z", d_dim, DIRECTION::Z_DIRECTION));
    // }
    
    /*
     * Initialize the filters.
     */
    
    d_filter_x = boost::shared_ptr<FilterTruncatedGaussian> (
        new FilterTruncatedGaussian("d_filter_x", d_dim, DIRECTION::X_DIRECTION));
    
    if ((d_dim == tbox::Dimension(2)) || (d_dim == tbox::Dimension(3)))
    {
        d_filter_y = boost::shared_ptr<FilterTruncatedGaussian> (
            new FilterTruncatedGaussian("d_filter_y", d_dim, DIRECTION::Y_DIRECTION));
    }
    
    if (d_dim == tbox::Dimension(3))
    {
        d_filter_z = boost::shared_ptr<FilterTruncatedGaussian> (
            new FilterTruncatedGaussian("d_filter_z", d_dim, DIRECTION::Z_DIRECTION));
    }
}


/*
 * Register the variables required for computing statistics.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::registerVariables(
    RungeKuttaLevelIntegrator* integrator,
    const hier::IntVector& num_ghosts)
{
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    integrator->registerVariable(
        s_variable_partial_density_unfiltered,
        num_ghosts,
        num_ghosts,
        RungeKuttaLevelIntegrator::STATISTICS,
        d_grid_geometry,
        "NO_COARSEN",
        "NO_REFINE");
    
    integrator->registerVariable(
        s_variable_momentum_unfiltered,
        num_ghosts,
        num_ghosts,
        RungeKuttaLevelIntegrator::STATISTICS,
        d_grid_geometry,
        "NO_COARSEN",
        "NO_REFINE");
    
    integrator->registerVariable(
        s_variable_total_energy_unfiltered,
        num_ghosts,
        num_ghosts,
        RungeKuttaLevelIntegrator::STATISTICS,
        d_grid_geometry,
        "NO_COARSEN",
        "NO_REFINE");
    
    integrator->registerVariable(
        s_variable_partial_density_filtered,
        num_ghosts,
        num_ghosts,
        RungeKuttaLevelIntegrator::STATISTICS,
        d_grid_geometry,
        "CONSERVATIVE_COARSEN",
        "CONSERVATIVE_LINEAR_REFINE");
    
    integrator->registerVariable(
        s_variable_momentum_filtered,
        num_ghosts,
        num_ghosts,
        RungeKuttaLevelIntegrator::STATISTICS,
        d_grid_geometry,
        "CONSERVATIVE_COARSEN",
        "CONSERVATIVE_LINEAR_REFINE");
    
    integrator->registerVariable(
        s_variable_total_energy_filtered,
        num_ghosts,
        num_ghosts,
        RungeKuttaLevelIntegrator::STATISTICS,
        d_grid_geometry,
        "CONSERVATIVE_COARSEN",
        "CONSERVATIVE_LINEAR_REFINE");
    
    integrator->registerVariable(
        s_variable_pressure_unfiltered,
        num_ghosts,
        num_ghosts,
        RungeKuttaLevelIntegrator::STATISTICS,
        d_grid_geometry,
        "NO_COARSEN",
        "NO_REFINE");
    
    integrator->registerVariable(
        s_variable_shear_stress_unfiltered,
        num_ghosts,
        num_ghosts,
        RungeKuttaLevelIntegrator::STATISTICS,
        d_grid_geometry,
        "CONSERVATIVE_COARSEN",
        "LINEAR_REFINE");
    
    integrator->registerVariable(
        s_variable_convective_stress_unfiltered,
        num_ghosts,
        num_ghosts,
        RungeKuttaLevelIntegrator::STATISTICS,
        d_grid_geometry,
        "NO_COARSEN",
        "NO_REFINE");
    
    integrator->registerVariable(
        s_variable_pressure_filtered,
        num_ghosts,
        num_ghosts,
        RungeKuttaLevelIntegrator::STATISTICS,
        d_grid_geometry,
        "CONSERVATIVE_COARSEN",
        "LINEAR_REFINE");
    
    integrator->registerVariable(
        s_variable_shear_stress_filtered,
        num_ghosts,
        num_ghosts,
        RungeKuttaLevelIntegrator::STATISTICS,
        d_grid_geometry,
        "CONSERVATIVE_COARSEN",
        "LINEAR_REFINE");
    
    integrator->registerVariable(
        s_variable_convective_stress_filtered,
        num_ghosts,
        num_ghosts,
        RungeKuttaLevelIntegrator::STATISTICS,
        d_grid_geometry,
        "CONSERVATIVE_COARSEN",
        "LINEAR_REFINE");
    
    integrator->registerVariable(
        s_variable_SFS_stress,
        num_ghosts,
        num_ghosts,
        RungeKuttaLevelIntegrator::STATISTICS,
        d_grid_geometry,
        "NO_COARSEN",
        "LINEAR_REFINE");
    
    integrator->registerVariable(
        s_variable_velocity_Favre_filtered,
        num_ghosts,
        num_ghosts,
        RungeKuttaLevelIntegrator::STATISTICS,
        d_grid_geometry,
        "NO_COARSEN",
        "NO_REFINE");
    
    integrator->registerVariable(
        s_variable_specific_volume_Favre_filtered,
        num_ghosts,
        num_ghosts,
        RungeKuttaLevelIntegrator::STATISTICS,
        d_grid_geometry,
        "NO_COARSEN",
        "NO_REFINE");
    
    integrator->registerVariable(
        s_variable_density_filtered,
        num_ghosts,
        num_ghosts,
        RungeKuttaLevelIntegrator::STATISTICS,
        d_grid_geometry,
        "NO_COARSEN",
        "NO_REFINE");
}


/*
 * Compute the variables required for computing statistics.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::computeVariables(
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
    tbox::pout << "FlowModelStatisticsUtilitiesFourEqnConservative::computeVariables: start" << std::endl;
    
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    // computePressure(
    //     patch_hierarchy,
    //     data_context);
    
    computeShearStress(
        patch_hierarchy,
        data_context);
    
    // computeConvectiveStress(
    //     patch_hierarchy,
    //     data_context);
    
    tbox::pout << "FlowModelStatisticsUtilitiesFourEqnConservative::computeVariables: end" << std::endl;
}


/*
 * Filter the variables required for computing statistics.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::filterVariables(
    const int level,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
    tbox::pout << "FlowModelStatisticsUtilitiesFourEqnConservative::filterVariables: start" << std::endl;
    
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    if (level == num_levels - 1)
    {
        if (!d_pressure_filtered)
        {
            computePressure(
                patch_hierarchy,
                data_context);
        }
        
        if (!d_convective_stress_filtered)
        {
            computeConvectiveStress(
                patch_hierarchy,
                data_context);
        }
    }
    
    filterPressure(
        level,
        patch_hierarchy,
        data_context);
    
    filterShearStress(
        level,
        patch_hierarchy,
        data_context);
    
    filterConvectiveStress(
        level,
        patch_hierarchy,
        data_context);
    
    // Get the flow model.
    boost::shared_ptr<FlowModel> d_flow_model_tmp = d_flow_model.lock();
    
    /*
     * Get the patch level.
     */
    
    boost::shared_ptr<hier::PatchLevel> patch_level(
        patch_hierarchy->getPatchLevel(level));
    
    for (hier::PatchLevel::iterator ip(patch_level->begin());
         ip != patch_level->end();
         ip++)
    {
        const boost::shared_ptr<hier::Patch> patch = *ip;
        
        d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
        
        // Get the conservative variables stored in flow model.
        
        boost::shared_ptr<pdat::CellData<double> > data_partial_density =
            d_flow_model_tmp->getGlobalCellData("PARTIAL_DENSITY");
        
        boost::shared_ptr<pdat::CellData<double> > data_momentum =
            d_flow_model_tmp->getGlobalCellData("MOMENTUM");
        
        boost::shared_ptr<pdat::CellData<double> > data_total_energy =
            d_flow_model_tmp->getGlobalCellData("TOTAL_ENERGY");
        
        // Get the data containers of filtered conservative variables.
        
        boost::shared_ptr<pdat::CellData<double> > data_partial_density_filtered(
            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                patch->getPatchData(s_variable_partial_density_filtered, data_context)));
        
        boost::shared_ptr<pdat::CellData<double> > data_momentum_filtered(
            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                patch->getPatchData(s_variable_momentum_filtered, data_context)));
        
        boost::shared_ptr<pdat::CellData<double> > data_total_energy_filtered(
            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                patch->getPatchData(s_variable_total_energy_filtered, data_context)));
        
        // Copy data to unfiltered variables for temporary storage.
        
        if (!d_conservative_variables_filtered)
        {
            // Get the data containers of unfiltered conservative variables.
            
            boost::shared_ptr<pdat::CellData<double> > data_partial_density_unfiltered(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch->getPatchData(s_variable_partial_density_unfiltered, data_context)));
            
            boost::shared_ptr<pdat::CellData<double> > data_momentum_unfiltered(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch->getPatchData(s_variable_momentum_unfiltered, data_context)));
            
            boost::shared_ptr<pdat::CellData<double> > data_total_energy_unfiltered(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch->getPatchData(s_variable_total_energy_unfiltered, data_context)));
            
            data_partial_density_unfiltered->copy(*data_partial_density);
            data_momentum_unfiltered->copy(*data_momentum);
            data_total_energy_unfiltered->copy(*data_total_energy);
            
            data_partial_density_filtered->copy(*data_partial_density);
            data_momentum_filtered->copy(*data_momentum);
            data_total_energy_filtered->copy(*data_total_energy);
        }
        
        // Apply filter in x-direction to conservative variables.
        
        for (int si = 0; si < d_num_species; si++)
        {
            d_filter_x->applyFilter(
                data_partial_density,
                data_partial_density_filtered,
                si,
                si);
        }
        
        for (int di = 0; di < d_dim.getValue(); di++)
        {
            d_filter_x->applyFilter(
                data_momentum,
                data_momentum_filtered,
                di,
                di);
        }
        
        d_filter_x->applyFilter(
            data_total_energy,
            data_total_energy_filtered,
            0,
            0);
        
        if ((d_dim == tbox::Dimension(2)) || (d_dim == tbox::Dimension(3)))
        {
            // Apply filter in y-direction to conservative variables.
            
            data_partial_density_filtered->copy(*data_partial_density);
            data_momentum_filtered->copy(*data_momentum);
            data_total_energy_filtered->copy(*data_total_energy);
            
            for (int si = 0; si < d_num_species; si++)
            {
                d_filter_y->applyFilter(
                    data_partial_density,
                    data_partial_density_filtered,
                    si,
                    si);
            }
            
            for (int di = 0; di < d_dim.getValue(); di++)
            {
                d_filter_y->applyFilter(
                    data_momentum,
                    data_momentum_filtered,
                    di,
                    di);
            }
            
            d_filter_y->applyFilter(
                data_total_energy,
                data_total_energy_filtered,
                0,
                0);
        }
        
        if (d_dim == tbox::Dimension(3))
        {
            // Apply filter in z-direction to conservative variables.
            
            data_partial_density_filtered->copy(*data_partial_density);
            data_momentum_filtered->copy(*data_momentum);
            data_total_energy_filtered->copy(*data_total_energy);
            
            for (int si = 0; si < d_num_species; si++)
            {
                d_filter_z->applyFilter(
                    data_partial_density,
                    data_partial_density_filtered,
                    si,
                    si);
            }
            
            for (int di = 0; di < d_dim.getValue(); di++)
            {
                d_filter_z->applyFilter(
                    data_momentum,
                    data_momentum_filtered,
                    di,
                    di);
            }
            
            d_filter_z->applyFilter(
                data_total_energy,
                data_total_energy_filtered,
                0,
                0);
        }
        
        data_partial_density_filtered->copy(*data_partial_density);
        data_momentum_filtered->copy(*data_momentum);
        data_total_energy_filtered->copy(*data_total_energy);
            
        /*
         * Unregister the patch and data of all registered derived cell variables in the flow model.
         */
        
        d_flow_model_tmp->unregisterPatch();
    }
    
    if (level == 0)
    {
        d_conservative_variables_filtered = true;
    }
    
    tbox::pout << "FlowModelStatisticsUtilitiesFourEqnConservative::filterVariables: end" << std::endl;
}


/*
 * Output names of statistical quantities to output to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputStatisticalQuantitiesNames(
    const std::string& stat_dump_filename)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    if (mpi.getRank() == 0)
    {
        // Loop over statistical quantities.
        for (int qi = 0; qi < static_cast<int>(d_statistical_quantities.size()); qi++)
        {
            // Get the key of the current variable.
            std::string statistical_quantity_key = d_statistical_quantities[qi];
            
            // DO NOTHING.
        }
    }
}


/*
 * Output statisitcal quantities to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputStatisticalQuantities(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double output_time)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    setConservativeVariablesToFilteredValues(
        patch_hierarchy,
        data_context);
    
    computeFavreFilteredVelocity(
        patch_hierarchy,
        data_context);
    
    computeFavreFilteredSpecificVolume(
        patch_hierarchy,
        data_context);
    
    computeFilteredDensity(
        patch_hierarchy,
        data_context);
    
    computeStressSFS(
        patch_hierarchy,
        data_context);
    
    // Loop over statistical quantities.
    for (int qi = 0; qi < static_cast<int>(d_statistical_quantities.size()); qi++)
    {
        // Get the key of the current variable.
        std::string statistical_quantity_key = d_statistical_quantities[qi];
        
        if (statistical_quantity_key == "DENSITY")
        {
            outputAveragedDensityWithInhomogeneousXDirection(
                "rho_mean.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "PRESSURE")
        {
            outputAveragedPressureWithInhomogeneousXDirection(
                "p_mean.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "VELOCITY_X")
        {
            outputAveragedVelocityXWithInhomogeneousXDirection(
                "u_mean.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "F_VELOCITY_X")
        {
            outputFavreAveragedVelocityXWithInhomogeneousXDirection(
                "u_tilde.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "SPECIFIC_VOLUME")
        {
            outputAveragedSpecificVolumeWithInhomogeneousXDirection(
                "rho_inv_mean.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "DENSITY_VARIANCE")
        {
            outputDensityVarianceWithInhomogeneousXDirection(
                "rho_variance.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "TKE")
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Not implemented!"
                << std::endl);
        }
        else if (statistical_quantity_key == "R11")
        {
            outputReynoldsNormalStressInXDirectionWithInhomogeneousXDirection(
                "R11.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "R22")
        {
            outputReynoldsNormalStressInYDirectionWithInhomogeneousXDirection(
                "R22.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "R33")
        {
            outputReynoldsNormalStressInZDirectionWithInhomogeneousXDirection(
                "R33.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "R12")
        {
            outputReynoldsShearStressInXYDirectionWithInhomogeneousXDirection(
                "R12.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "R13")
        {
            outputReynoldsShearStressInXZDirectionWithInhomogeneousXDirection(
                "R13.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "R23")
        {
            outputReynoldsShearStressInYZDirectionWithInhomogeneousXDirection(
                "R23.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "ra1")
        {
            outputAveragedTurbMassFluxXWithInhomogeneousXDirection(
                "ra1.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "ra2")
        {
            outputAveragedTurbMassFluxYWithInhomogeneousXDirection(
                "ra2.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "ra3")
        {
            outputAveragedTurbMassFluxZWithInhomogeneousXDirection(
                "ra3.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "u_p_u_p")
        {
            outputVelocityComponentInXDirectionVarianceWithInhomogeneousXDirection(
                "u_p_u_p.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "v_p_v_p")
        {
            outputVelocityComponentInYDirectionVarianceWithInhomogeneousXDirection(
                "v_p_v_p.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "w_p_w_p")
        {
            outputVelocityComponentInZDirectionVarianceWithInhomogeneousXDirection(
                "w_p_w_p.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "rho_p_u_p_u_p")
        {
            outputDensityVelocityComponentSquareInXDirectionCorrelationWithInhomogeneousXDirection(
                "rho_p_u_p_u_p.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "rho_p_v_p_v_p")
        {
            outputDensityVelocityComponentSquareInYDirectionCorrelationWithInhomogeneousXDirection(
                "rho_p_v_p_v_p.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "rho_p_w_p_w_p")
        {
            outputDensityVelocityComponentSquareInZDirectionCorrelationWithInhomogeneousXDirection(
                "rho_p_w_p_w_p.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "b")
        {
            outputDensitySpecificVolumeCovarianceWithInhomogeneousXDirection(
                "b.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "ra1_budget_filtered")
        {
            outputBudgetFilteredTurbMassFluxXWithInhomogeneousXDirection(
                "ra1_budget_filtered.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "rb_budget_filtered")
        {
            outputBudgetFilteredDensitySpecificVolumeCovarianceWithInhomogeneousXDirection(
                "rb_budget_filtered.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "rR11_budget_filtered")
        {
            outputBudgetFilteredReynoldsNormalStressInXDirectionWithInhomogeneousXDirection(
                "rR11_budget_filtered.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "rR22_budget_filtered")
        {
            outputBudgetFilteredReynoldsNormalStressInYDirectionWithInhomogeneousXDirection(
                "rR22_budget_filtered.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "rR33_budget_filtered")
        {
            outputBudgetFilteredReynoldsNormalStressInZDirectionWithInhomogeneousXDirection(
                "rR33_budget_filtered.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Unknown statistical quantity key = '"
                << statistical_quantity_key
                << " found."
                << std::endl);
        }
    }
    
    resetConservativeVariablesToUnfilteredValues(
        patch_hierarchy,
        data_context);
    
    d_pressure_filtered = false;
    d_shear_stress_filtered = false;
    d_convective_stress_filtered = false;
    d_conservative_variables_filtered = false;
}


/*
 * Set conservative variables to filtered value.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::setConservativeVariablesToFilteredValues(
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    // Get the flow model.
    boost::shared_ptr<FlowModel> d_flow_model_tmp = d_flow_model.lock();
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    for (int li = 0; li < num_levels; li++)
    {
        /*
         * Get the current patch level.
         */
        
        boost::shared_ptr<hier::PatchLevel> patch_level(
            patch_hierarchy->getPatchLevel(li));
        
        for (hier::PatchLevel::iterator ip(patch_level->begin());
             ip != patch_level->end();
             ip++)
        {
            const boost::shared_ptr<hier::Patch> patch = *ip;
            
            d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
            
            // Get the conservative variables stored in flow model.
            
            boost::shared_ptr<pdat::CellData<double> > data_partial_density =
                d_flow_model_tmp->getGlobalCellData("PARTIAL_DENSITY");
            
            boost::shared_ptr<pdat::CellData<double> > data_momentum =
                d_flow_model_tmp->getGlobalCellData("MOMENTUM");
            
            boost::shared_ptr<pdat::CellData<double> > data_total_energy =
                d_flow_model_tmp->getGlobalCellData("TOTAL_ENERGY");
            
            // Get the data containers of filtered conservative variables.
            
            boost::shared_ptr<pdat::CellData<double> > data_partial_density_filtered(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch->getPatchData(s_variable_partial_density_filtered, data_context)));
            
            boost::shared_ptr<pdat::CellData<double> > data_momentum_filtered(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch->getPatchData(s_variable_momentum_filtered, data_context)));
            
            boost::shared_ptr<pdat::CellData<double> > data_total_energy_filtered(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch->getPatchData(s_variable_total_energy_filtered, data_context)));
            
            data_partial_density->copy(*data_partial_density_filtered);
            data_momentum->copy(*data_momentum_filtered);
            data_total_energy->copy(*data_total_energy_filtered);
            
            /*
             * Unregister the patch and data of all registered derived cell variables in the flow model.
             */
            
            d_flow_model_tmp->unregisterPatch();
        }
    }
}


/*
 * Reset conservative variables to unfiltered value.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::resetConservativeVariablesToUnfilteredValues(
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    // Get the flow model.
    boost::shared_ptr<FlowModel> d_flow_model_tmp = d_flow_model.lock();
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    for (int li = 0; li < num_levels; li++)
    {
        /*
         * Get the current patch level.
         */
        
        boost::shared_ptr<hier::PatchLevel> patch_level(
            patch_hierarchy->getPatchLevel(li));
        
        for (hier::PatchLevel::iterator ip(patch_level->begin());
             ip != patch_level->end();
             ip++)
        {
            const boost::shared_ptr<hier::Patch> patch = *ip;
            
            d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
            
            // Get the conservative variables stored in flow model.
            
            boost::shared_ptr<pdat::CellData<double> > data_partial_density =
                d_flow_model_tmp->getGlobalCellData("PARTIAL_DENSITY");
            
            boost::shared_ptr<pdat::CellData<double> > data_momentum =
                d_flow_model_tmp->getGlobalCellData("MOMENTUM");
            
            boost::shared_ptr<pdat::CellData<double> > data_total_energy =
                d_flow_model_tmp->getGlobalCellData("TOTAL_ENERGY");
            
            // Get the data containers of unfiltered conservative variables.
            
            boost::shared_ptr<pdat::CellData<double> > data_partial_density_unfiltered(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch->getPatchData(s_variable_partial_density_unfiltered, data_context)));
            
            boost::shared_ptr<pdat::CellData<double> > data_momentum_unfiltered(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch->getPatchData(s_variable_momentum_unfiltered, data_context)));
            
            boost::shared_ptr<pdat::CellData<double> > data_total_energy_unfiltered(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch->getPatchData(s_variable_total_energy_unfiltered, data_context)));
            
            // // Get the data containers of filtered conservative variables.
            
            // boost::shared_ptr<pdat::CellData<double> > data_partial_density_filtered(
            //     BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
            //         patch->getPatchData(s_variable_partial_density_filtered, data_context)));
            
            // boost::shared_ptr<pdat::CellData<double> > data_momentum_filtered(
            //     BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
            //         patch->getPatchData(s_variable_momentum_filtered, data_context)));
            
            // boost::shared_ptr<pdat::CellData<double> > data_total_energy_filtered(
            //     BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
            //         patch->getPatchData(s_variable_total_energy_filtered, data_context)));
            
            // Reset conservative variables to unfiltered value.
            
            data_partial_density->copy(*data_partial_density_unfiltered);
            data_momentum->copy(*data_momentum_unfiltered);
            data_total_energy->copy(*data_total_energy_unfiltered);
            
            /*
             * Unregister the patch and data of all registered derived cell variables in the flow model.
             */
            
            d_flow_model_tmp->unregisterPatch();
        }
    }
}


/*
 * Output averaged density with inhomogeneous x-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputAveragedDensityWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double output_time) const
{
    std::vector<double> rho_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        patch_hierarchy,
        data_context);
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_output;
        f_output.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        
        f_output.write((char*)&output_time, sizeof(double));
        f_output.write((char*)&rho_mean[0], sizeof(double)*rho_mean.size());
        
        f_output.close();
    }
}


/*
 * Output averaged pressure with inhomogeneous x-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputAveragedPressureWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double output_time) const
{
    std::vector<double> p_mean = getAveragedQuantityWithInhomogeneousXDirection(
        s_variable_pressure_filtered,
        0,
        patch_hierarchy,
        data_context);
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_output;
        f_output.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        
        f_output.write((char*)&output_time, sizeof(double));
        f_output.write((char*)&p_mean[0], sizeof(double)*p_mean.size());
        
        f_output.close();
    }
}


/*
 * Output averaged velocity component in x-direction with inhomogeneous x-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputAveragedVelocityXWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double output_time) const
{
    std::vector<double> u_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "VELOCITY",
        0,
        patch_hierarchy,
        data_context);
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_output;
        f_output.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        
        f_output.write((char*)&output_time, sizeof(double));
        f_output.write((char*)&u_mean[0], sizeof(double)*u_mean.size());
        
        f_output.close();
    }
}


/*
 * Output Favre-averaged velocity component in x-direction with inhomogeneous x-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputFavreAveragedVelocityXWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double output_time) const
{
    const int finest_level_dim_0 = getRefinedDomainNumberOfPointsX(patch_hierarchy);
    
    std::vector<double> rho_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        patch_hierarchy,
        data_context);
    
    std::vector<double> rho_u_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "MOMENTUM",
        0,
        patch_hierarchy,
        data_context);
    
    std::vector<double> u_tilde(rho_u_mean);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        u_tilde[i] /= rho_mean[i];
    }
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_output;
        f_output.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        
        f_output.write((char*)&output_time, sizeof(double));
        f_output.write((char*)&u_tilde[0], sizeof(double)*u_tilde.size());
        
        f_output.close();
    }
}


/*
 * Output averaged specific volume with inhomogeneous x-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputAveragedSpecificVolumeWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double output_time) const
{
    std::vector<double> v_mean = getAveragedReciprocalOfQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        patch_hierarchy,
        data_context);
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_output;
        f_output.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        
        f_output.write((char*)&output_time, sizeof(double));
        f_output.write((char*)&v_mean[0], sizeof(double)*v_mean.size());
        
        f_output.close();
    }
}



/*
 * Output turbulent mass flux in x-direction with inhomogeneous x-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputAveragedTurbMassFluxXWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double output_time) const
{
    std::vector<double> rho_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        patch_hierarchy,
        data_context);
    
    std::vector<double> u_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "VELOCITY",
        0,
        patch_hierarchy,
        data_context);
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    std::vector<std::vector<double> > averaged_quantities;
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(rho_mean);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(u_mean);
    
    std::vector<double> rho_p_u_p = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_output;
        f_output.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        
        f_output.write((char*)&output_time, sizeof(double));
        f_output.write((char*)&rho_p_u_p[0], sizeof(double)*rho_p_u_p.size());
        
        f_output.close();
    }
}


/*
 * Output turbulent mass flux in y-direction with inhomogeneous x-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputAveragedTurbMassFluxYWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double output_time) const
{
    std::vector<double> rho_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        patch_hierarchy,
        data_context);
    
    std::vector<double> v_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "VELOCITY",
        1,
        patch_hierarchy,
        data_context);
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    std::vector<std::vector<double> > averaged_quantities;
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(rho_mean);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(1);
    averaged_quantities.push_back(v_mean);
    
    std::vector<double> rho_p_v_p = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_output;
        f_output.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        
        f_output.write((char*)&output_time, sizeof(double));
        f_output.write((char*)&rho_p_v_p[0], sizeof(double)*rho_p_v_p.size());
        
        f_output.close();
    }
}


/*
 * Output turbulent mass flux in z-direction with inhomogeneous x-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputAveragedTurbMassFluxZWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double output_time) const
{
    std::vector<double> rho_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        patch_hierarchy,
        data_context);
    
    std::vector<double> w_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "VELOCITY",
        2,
        patch_hierarchy,
        data_context);
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    std::vector<std::vector<double> > averaged_quantities;
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(rho_mean);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(2);
    averaged_quantities.push_back(w_mean);
    
    std::vector<double> rho_p_w_p = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_output;
        f_output.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        
        f_output.write((char*)&output_time, sizeof(double));
        f_output.write((char*)&rho_p_w_p[0], sizeof(double)*rho_p_w_p.size());
        
        f_output.close();
    }
}


/*
 * Output Reynolds normal stress in x-direction with inhomogeneous x-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputReynoldsNormalStressInXDirectionWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double output_time) const
{
    const int finest_level_dim_0 = getRefinedDomainNumberOfPointsX(patch_hierarchy);
    
    std::vector<double> rho_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        patch_hierarchy,
        data_context);
    
    // Compute u_tilde.
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    
    std::vector<double> rho_u_mean = getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        patch_hierarchy,
        data_context);
    
    std::vector<double> u_tilde(rho_u_mean);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        u_tilde[i] /= rho_mean[i];
    }
    
    // Compute R_11.
    
    std::vector<double> zeros(finest_level_dim_0, double(0));
    
    quantity_names.clear();
    component_indices.clear();
    
    std::vector<std::vector<double> > averaged_quantities;
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(zeros);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(u_tilde);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(u_tilde);
    
    std::vector<double> R_11 = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        R_11[i] /= rho_mean[i];
    }
    
    // Output R_11.
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_output;
        f_output.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        
        f_output.write((char*)&output_time, sizeof(double));
        f_output.write((char*)&R_11[0], sizeof(double)*R_11.size());
        
        f_output.close();
    }
}


/*
 * Output Reynolds normal stress in y-direction with inhomogeneous x-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputReynoldsNormalStressInYDirectionWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double output_time) const
{
    const int finest_level_dim_0 = getRefinedDomainNumberOfPointsX(patch_hierarchy);
    
    std::vector<double> rho_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        patch_hierarchy,
        data_context);
    
    // Compute v_tilde.
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(1);
    
    std::vector<double> rho_v_mean = getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        patch_hierarchy,
        data_context);
    
    std::vector<double> v_tilde(rho_v_mean);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        v_tilde[i] /= rho_mean[i];
    }
    
    // Compute R_22.
    
    std::vector<double> zeros(finest_level_dim_0, double(0));
    
    quantity_names.clear();
    component_indices.clear();
    
    std::vector<std::vector<double> > averaged_quantities;
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(zeros);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(1);
    averaged_quantities.push_back(v_tilde);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(1);
    averaged_quantities.push_back(v_tilde);
    
    std::vector<double> R_22 = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        R_22[i] /= rho_mean[i];
    }
    
    // Output R_22.
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_output;
        f_output.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        
        f_output.write((char*)&output_time, sizeof(double));
        f_output.write((char*)&R_22[0], sizeof(double)*R_22.size());
        
        f_output.close();
    }
}


/*
 * Output Reynolds normal stress in z-direction with inhomogeneous x-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputReynoldsNormalStressInZDirectionWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double output_time) const
{
    const int finest_level_dim_0 = getRefinedDomainNumberOfPointsX(patch_hierarchy);
    
    std::vector<double> rho_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        patch_hierarchy,
        data_context);
    
    // Compute w_tilde.
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(2);
    
    std::vector<double> rho_w_mean = getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        patch_hierarchy,
        data_context);
    
    std::vector<double> w_tilde(rho_w_mean);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        w_tilde[i] /= rho_mean[i];
    }
    
    // Compute R_33.
    
    std::vector<double> zeros(finest_level_dim_0, double(0));
    
    quantity_names.clear();
    component_indices.clear();
    
    std::vector<std::vector<double> > averaged_quantities;
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(zeros);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(2);
    averaged_quantities.push_back(w_tilde);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(2);
    averaged_quantities.push_back(w_tilde);
    
    std::vector<double> R_33 = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        R_33[i] /= rho_mean[i];
    }
    
    // Output R_33.
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_output;
        f_output.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        
        f_output.write((char*)&output_time, sizeof(double));
        f_output.write((char*)&R_33[0], sizeof(double)*R_33.size());
        
        f_output.close();
    }
}


/*
 * Output Reynolds shear stress in xy-direction with inhomogeneous x-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputReynoldsShearStressInXYDirectionWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double output_time) const
{
    const int finest_level_dim_0 = getRefinedDomainNumberOfPointsX(patch_hierarchy);
    
    std::vector<double> rho_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        patch_hierarchy,
        data_context);
    
    // Compute u_tilde.
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    
    std::vector<double> rho_u_mean = getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        patch_hierarchy,
        data_context);
    
    std::vector<double> u_tilde(rho_u_mean);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        u_tilde[i] /= rho_mean[i];
    }
    
    // Compute v_tilde.
    
    quantity_names.clear();
    component_indices.clear();
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(1);
    
    std::vector<double> rho_v_mean = getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        patch_hierarchy,
        data_context);
    
    std::vector<double> v_tilde(rho_v_mean);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        v_tilde[i] /= rho_mean[i];
    }
    
    // Compute R_12.
    
    std::vector<double> zeros(finest_level_dim_0, double(0));
    
    quantity_names.clear();
    component_indices.clear();
    
    std::vector<std::vector<double> > averaged_quantities;
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(zeros);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(u_tilde);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(1);
    averaged_quantities.push_back(v_tilde);
    
    std::vector<double> R_12 = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        R_12[i] /= rho_mean[i];
    }
    
    // Output R_12.
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_output;
        f_output.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        
        f_output.write((char*)&output_time, sizeof(double));
        f_output.write((char*)&R_12[0], sizeof(double)*R_12.size());
        
        f_output.close();
    }
}


/*
 * Output Reynolds shear stress in xz-direction with inhomogeneous x-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputReynoldsShearStressInXZDirectionWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double output_time) const
{
    const int finest_level_dim_0 = getRefinedDomainNumberOfPointsX(patch_hierarchy);
    
    std::vector<double> rho_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        patch_hierarchy,
        data_context);
    
    // Compute u_tilde.
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    
    std::vector<double> rho_u_mean = getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        patch_hierarchy,
        data_context);
    
    std::vector<double> u_tilde(rho_u_mean);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        u_tilde[i] /= rho_mean[i];
    }
    
    // Compute w_tilde.
    
    quantity_names.clear();
    component_indices.clear();
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(2);
    
    std::vector<double> rho_w_mean = getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        patch_hierarchy,
        data_context);
    
    std::vector<double> w_tilde(rho_w_mean);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        w_tilde[i] /= rho_mean[i];
    }
    
    // Compute R_13.
    
    std::vector<double> zeros(finest_level_dim_0, double(0));
    
    quantity_names.clear();
    component_indices.clear();
    
    std::vector<std::vector<double> > averaged_quantities;
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(zeros);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(u_tilde);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(2);
    averaged_quantities.push_back(w_tilde);
    
    std::vector<double> R_13 = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        R_13[i] /= rho_mean[i];
    }
    
    // Output R_13.
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_output;
        f_output.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        
        f_output.write((char*)&output_time, sizeof(double));
        f_output.write((char*)&R_13[0], sizeof(double)*R_13.size());
        
        f_output.close();
    }
}


/*
 * Output Reynolds shear stress in yz-direction with inhomogeneous x-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputReynoldsShearStressInYZDirectionWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double output_time) const
{
    const int finest_level_dim_0 = getRefinedDomainNumberOfPointsX(patch_hierarchy);
    
    std::vector<double> rho_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        patch_hierarchy,
        data_context);
    
    // Compute v_tilde.
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(1);
    
    std::vector<double> rho_v_mean = getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        patch_hierarchy,
        data_context);
    
    std::vector<double> v_tilde(rho_v_mean);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        v_tilde[i] /= rho_mean[i];
    }
    
    // Compute w_tilde.
    
    quantity_names.clear();
    component_indices.clear();
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(2);
    
    std::vector<double> rho_w_mean = getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        patch_hierarchy,
        data_context);
    
    std::vector<double> w_tilde(rho_w_mean);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        w_tilde[i] /= rho_mean[i];
    }
    
    // Compute R_23.
    
    std::vector<double> zeros(finest_level_dim_0, double(0));
    
    quantity_names.clear();
    component_indices.clear();
    
    std::vector<std::vector<double> > averaged_quantities;
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(zeros);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(1);
    averaged_quantities.push_back(v_tilde);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(2);
    averaged_quantities.push_back(w_tilde);
    
    std::vector<double> R_23 = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        R_23[i] /= rho_mean[i];
    }
    
    // Output R_23.
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_output;
        f_output.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        
        f_output.write((char*)&output_time, sizeof(double));
        f_output.write((char*)&R_23[0], sizeof(double)*R_23.size());
        
        f_output.close();
    }
}


/*
 * Output variance of velocity component in x-direction with inhomogeneous x-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputVelocityComponentInXDirectionVarianceWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double output_time) const
{
    std::vector<double> u_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "VELOCITY",
        0,
        patch_hierarchy,
        data_context);
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    std::vector<std::vector<double> > averaged_quantities;
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(u_mean);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(u_mean);
    
    std::vector<double> u_p_u_p = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_output;
        f_output.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        
        f_output.write((char*)&output_time, sizeof(double));
        f_output.write((char*)&u_p_u_p[0], sizeof(double)*u_p_u_p.size());
        
        f_output.close();
    }
}


/*
 * Output variance of velocity component in y-direction with inhomogeneous x-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputVelocityComponentInYDirectionVarianceWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double output_time) const
{
    std::vector<double> v_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "VELOCITY",
        1,
        patch_hierarchy,
        data_context);
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    std::vector<std::vector<double> > averaged_quantities;
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(1);
    averaged_quantities.push_back(v_mean);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(1);
    averaged_quantities.push_back(v_mean);
    
    std::vector<double> v_p_v_p = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_output;
        f_output.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        
        f_output.write((char*)&output_time, sizeof(double));
        f_output.write((char*)&v_p_v_p[0], sizeof(double)*v_p_v_p.size());
        
        f_output.close();
    }
}


/*
 * Output variance of velocity component in z-direction with inhomogeneous x-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputVelocityComponentInZDirectionVarianceWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double output_time) const
{
    std::vector<double> w_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "VELOCITY",
        2,
        patch_hierarchy,
        data_context);
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    std::vector<std::vector<double> > averaged_quantities;
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(2);
    averaged_quantities.push_back(w_mean);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(2);
    averaged_quantities.push_back(w_mean);
    
    std::vector<double> w_p_w_p = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_output;
        f_output.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        
        f_output.write((char*)&output_time, sizeof(double));
        f_output.write((char*)&w_p_w_p[0], sizeof(double)*w_p_w_p.size());
        
        f_output.close();
    }
}


/*
 * Output correlation of density and velocity component square in x-direction with inhomogeneous x-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::
outputDensityVelocityComponentSquareInXDirectionCorrelationWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double output_time) const
{
    std::vector<double> rho_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        patch_hierarchy,
        data_context);
    
    std::vector<double> u_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "VELOCITY",
        0,
        patch_hierarchy,
        data_context);
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    std::vector<std::vector<double> > averaged_quantities;
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(rho_mean);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(u_mean);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(u_mean);
    
    std::vector<double> rho_p_u_p_u_p = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_output;
        f_output.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        
        f_output.write((char*)&output_time, sizeof(double));
        f_output.write((char*)&rho_p_u_p_u_p[0], sizeof(double)*rho_p_u_p_u_p.size());
        
        f_output.close();
    }
}


/*
 * Output correlation of density and velocity component square in y-direction with inhomogeneous x-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::
outputDensityVelocityComponentSquareInYDirectionCorrelationWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double output_time) const
{
    std::vector<double> rho_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        patch_hierarchy,
        data_context);
    
    std::vector<double> v_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "VELOCITY",
        1,
        patch_hierarchy,
        data_context);
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    std::vector<std::vector<double> > averaged_quantities;
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(rho_mean);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(1);
    averaged_quantities.push_back(v_mean);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(1);
    averaged_quantities.push_back(v_mean);
    
    std::vector<double> rho_p_v_p_v_p = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_output;
        f_output.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        
        f_output.write((char*)&output_time, sizeof(double));
        f_output.write((char*)&rho_p_v_p_v_p[0], sizeof(double)*rho_p_v_p_v_p.size());
        
        f_output.close();
    }
}


/*
 * Output correlation of density and velocity component square in z-direction with inhomogeneous x-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::
outputDensityVelocityComponentSquareInZDirectionCorrelationWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double output_time) const
{
    std::vector<double> rho_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        patch_hierarchy,
        data_context);
    
    std::vector<double> w_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "VELOCITY",
        2,
        patch_hierarchy,
        data_context);
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    std::vector<std::vector<double> > averaged_quantities;
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(rho_mean);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(2);
    averaged_quantities.push_back(w_mean);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(2);
    averaged_quantities.push_back(w_mean);
    
    std::vector<double> rho_p_w_p_w_p = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_output;
        f_output.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        
        f_output.write((char*)&output_time, sizeof(double));
        f_output.write((char*)&rho_p_w_p_w_p[0], sizeof(double)*rho_p_w_p_w_p.size());
        
        f_output.close();
    }
}



/*
 * Output density variance with inhomogeneous x-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputDensityVarianceWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double output_time) const
{
    std::vector<double> rho_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        patch_hierarchy,
        data_context);
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    std::vector<std::vector<double> > averaged_quantities;
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(rho_mean);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(rho_mean);
    
    std::vector<double> rho_p_rho_p = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_output;
        f_output.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        
        f_output.write((char*)&output_time, sizeof(double));
        f_output.write((char*)&rho_p_rho_p[0], sizeof(double)*rho_p_rho_p.size());
        
        f_output.close();
    }
}


/*
 * Output density specific volume covariance with inhomogeneous x-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputDensitySpecificVolumeCovarianceWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double output_time) const
{
    const int finest_level_dim_0 = getRefinedDomainNumberOfPointsX(patch_hierarchy);
    
    std::vector<double> rho_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        patch_hierarchy,
        data_context);
    
    std::vector<double> v_mean = getAveragedReciprocalOfQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        patch_hierarchy,
        data_context);
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    std::vector<bool> use_reciprocal;
    std::vector<std::vector<double> > averaged_quantities;
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_reciprocal.push_back(false);
    averaged_quantities.push_back(rho_mean);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_reciprocal.push_back(true);
    averaged_quantities.push_back(v_mean);
    
    std::vector<double> rho_p_v_p = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        use_reciprocal,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    std::vector<double> b(rho_p_v_p);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        b[i] = -b[i];
    }
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_output;
        f_output.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        
        f_output.write((char*)&output_time, sizeof(double));
        f_output.write((char*)&b[0], sizeof(double)*b.size());
        
        f_output.close();
    }
}


/**
 ** Helper functions.
 **/


/*
 * Get number of points in the x-direction of the refined domain.
 */
int
FlowModelStatisticsUtilitiesFourEqnConservative::getRefinedDomainNumberOfPointsX(
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy) const
{
    /*
     * Get the refinement ratio from the finest level to the coarest level.
     */
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    hier::IntVector ratioFinestLevelToCoarestLevel =
        patch_hierarchy->getRatioToCoarserLevel(num_levels - 1);
    for (int li = num_levels - 2; li > 0 ; li--)
    {
        ratioFinestLevelToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(li);
    }
    
    /*
     * Get the number of cells of physical domain refined to the finest level.
     */
    
    const hier::BoxContainer& physical_domain = d_grid_geometry->getPhysicalDomain();
    const hier::Box& physical_domain_box = physical_domain.front();
    const hier::IntVector& physical_domain_dims = physical_domain_box.numberCells();
    const hier::IntVector finest_level_dims = physical_domain_dims*ratioFinestLevelToCoarestLevel;
    
    const int finest_level_dim_0 = finest_level_dims[0];
    return finest_level_dim_0;
}


/*
 * Get grid spacing in the x-direction of the refined domain.
 */
double
FlowModelStatisticsUtilitiesFourEqnConservative::getRefinedDomainGridSpacingX(
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy) const
{
    /*
     * Get the refinement ratio from the finest level to the coarest level.
     */
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    hier::IntVector ratioFinestLevelToCoarestLevel =
        patch_hierarchy->getRatioToCoarserLevel(num_levels - 1);
    for (int li = num_levels - 2; li > 0 ; li--)
    {
        ratioFinestLevelToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(li);
    }
    
    const double* dx = d_grid_geometry->getDx();
    
    return dx[0]/ratioFinestLevelToCoarestLevel[0];
}


/*
 * Compute the one-dimensional derivative given a vector.
 */
std::vector<double> FlowModelStatisticsUtilitiesFourEqnConservative::computeDerivativeOfVector1D(
    const std::vector<double> quantity_vector,
    const double dx) const
{
    const int vector_length = quantity_vector.size();
    
    std::vector<double> derivative;
    derivative.resize(vector_length);
    
    const double* u = quantity_vector.data();
    double* dudx = derivative.data();
    
    // HARD CODE TO BE SIXTH ORDER CENTRAL SCHEME FOR DIFFERENTIATION.
    for (int i = 3; i < vector_length - 3; i++)
    {
        // Compute linear indices.
        const int idx     = i;
        
        const int idx_LLL = i - 3;
        const int idx_LL  = i - 2;
        const int idx_L   = i - 1;
        const int idx_R   = i + 1;
        const int idx_RR  = i + 2;
        const int idx_RRR = i + 3;
        
        dudx[idx] = (double(1)/double(60)*(u[idx_RRR] - u[idx_LLL])
            - double(3)/double(20)*(u[idx_RR] - u[idx_LL])
            + double(3)/double(4)*(u[idx_R] - u[idx_L]))/dx;
    }
    
    for (int i = 0; i < 3; i++)
    {
        dudx[i]                     = dudx[3];
        dudx[vector_length - i - 1] = dudx[vector_length - 4];
    }
    
    return derivative;
}


/*
 * Compute averaged value with only x direction as inhomogeneous direction.
 */
std::vector<double>
FlowModelStatisticsUtilitiesFourEqnConservative::getAveragedQuantityWithInhomogeneousXDirection(
    const std::string quantity_name,
    const int component_idx,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context) const
{
    std::vector<double> averaged_quantity;
    
    boost::shared_ptr<FlowModel> d_flow_model_tmp = d_flow_model.lock();
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Get the refinement ratio from the finest level to the coarest level.
     */
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    hier::IntVector ratioFinestLevelToCoarestLevel =
        patch_hierarchy->getRatioToCoarserLevel(num_levels - 1);
    for (int li = num_levels - 2; li > 0 ; li--)
    {
        ratioFinestLevelToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(li);
    }
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    boost::shared_ptr<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *patch_hierarchy,
            0,
            num_levels - 1));
    
    /*
     * Get the number of cells of physical domain refined to the finest level.
     */
    
    const hier::BoxContainer& physical_domain = d_grid_geometry->getPhysicalDomain();
    const hier::Box& physical_domain_box = physical_domain.front();
    const hier::IntVector& physical_domain_dims = physical_domain_box.numberCells();
    const hier::IntVector finest_level_dims = physical_domain_dims*ratioFinestLevelToCoarestLevel;
    
    /*
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        double* u_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_quantity.resize(finest_level_dim_0);
        double* u_avg_global = averaged_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_avg_local[i] = double(0);
            u_avg_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_quantity =
                    d_flow_model_tmp->getGlobalCellData(quantity_name);
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        /*
                         * Compute the index of the data point and count how many times the data is repeated.
                         */
                        
                        const hier::Index idx_pt(tbox::Dimension(1), idx_lo_0 + i);
                        
                        int n_overlapped = 1;
                        
                        for (hier::BoxContainer::BoxContainerConstIterator iob(
                                patch_overlapped_visible_boxes.begin());
                             iob != patch_overlapped_visible_boxes.end();
                             iob++)
                        {
                            const hier::Box& patch_overlapped_visible_box = *iob;
                            
                            if (patch_overlapped_visible_box.contains(idx_pt))
                            {
                                n_overlapped++;
                            }
                        }
                        
                        /*
                         * Compute the linear index and the data to add.
                         */
                        
                        const int idx = relative_idx_lo_0 + i + num_ghosts_0_quantity;
                        
                        const double value_to_add = u[idx]/((double) n_overlapped);
                        
                        for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                        {
                            const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                            
                            u_avg_local[idx_fine] += value_to_add;
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        mpi.Allreduce(
            u_avg_local,
            u_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_avg_local);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        
        double* u_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_quantity.resize(finest_level_dim_0);
        double* u_avg_global = averaged_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_avg_local[i] = double(0);
            u_avg_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_quantity =
                    d_flow_model_tmp->getGlobalCellData(quantity_name);
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                
                const double weight = dx[1]/L_y;
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            /*
                             * Compute the index of the data point and count how many times the data is repeated.
                             */
                            
                            const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j);
                            
                            int n_overlapped = 1;
                            
                            for (hier::BoxContainer::BoxContainerConstIterator iob(
                                    patch_overlapped_visible_boxes.begin());
                                 iob != patch_overlapped_visible_boxes.end();
                                 iob++)
                            {
                                const hier::Box& patch_overlapped_visible_box = *iob;
                                
                                if (patch_overlapped_visible_box.contains(idx_pt))
                                {
                                    n_overlapped++;
                                }
                            }
                            
                            /*
                             * Compute the linear index and the data to add.
                             */
                            
                            const int idx = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                            
                            const double value_to_add = u[idx]*weight/((double) n_overlapped);
                            
                            for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                            {
                                const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                
                                u_avg_local[idx_fine] += value_to_add;
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        mpi.Allreduce(
            u_avg_local,
            u_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_avg_local);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        const double L_z = x_hi[2] - x_lo[2];
        
        double* u_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_quantity.resize(finest_level_dim_0);
        double* u_avg_global = averaged_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_avg_local[i] = double(0);
            u_avg_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_quantity =
                    d_flow_model_tmp->getGlobalCellData(quantity_name);
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int num_ghosts_2_quantity = num_ghosts_quantity[2];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                const int ghostcell_dim_1_quantity = ghostcell_dims_quantity[1];
                
                const double weight = (dx[1]*dx[2])/(L_y*L_z);
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    const int interior_dim_2 = interior_dims[2];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int idx_lo_2 = index_lo[2];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    const int relative_idx_lo_2 = relative_index_lo[2];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                /*
                                 * Compute the index of the data point and count how many times the data is repeated.
                                 */
                                
                                const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j, idx_lo_2 + k);
                                
                                int n_overlapped = 1;
                                
                                for (hier::BoxContainer::BoxContainerConstIterator iob(
                                        patch_overlapped_visible_boxes.begin());
                                     iob != patch_overlapped_visible_boxes.end();
                                     iob++)
                                {
                                    const hier::Box& patch_overlapped_visible_box = *iob;
                                    
                                    if (patch_overlapped_visible_box.contains(idx_pt))
                                    {
                                        n_overlapped++;
                                    }
                                }
                                
                                /*
                                 * Compute the linear index and the data to add.
                                 */
                                
                                const int idx = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                    (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                        ghostcell_dim_1_quantity;
                                
                                const double value_to_add = u[idx]*weight/((double) n_overlapped);
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    u_avg_local[idx_fine] += value_to_add;
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        mpi.Allreduce(
            u_avg_local,
            u_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_avg_local);
    }
    
    return averaged_quantity;
}


/*
 * Compute averaged reciprocal of value with only x direction as inhomogeneous direction.
 */
std::vector<double>
FlowModelStatisticsUtilitiesFourEqnConservative::getAveragedReciprocalOfQuantityWithInhomogeneousXDirection(
    const std::string quantity_name,
    const int component_idx,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context) const
{
    std::vector<double> averaged_reciprocal_quantity;
    
    boost::shared_ptr<FlowModel> d_flow_model_tmp = d_flow_model.lock();
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Get the refinement ratio from the finest level to the coarest level.
     */
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    hier::IntVector ratioFinestLevelToCoarestLevel =
        patch_hierarchy->getRatioToCoarserLevel(num_levels - 1);
    for (int li = num_levels - 2; li > 0 ; li--)
    {
        ratioFinestLevelToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(li);
    }
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    boost::shared_ptr<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *patch_hierarchy,
            0,
            num_levels - 1));
    
    /*
     * Get the number of cells of physical domain refined to the finest level.
     */
    
    const hier::BoxContainer& physical_domain = d_grid_geometry->getPhysicalDomain();
    const hier::Box& physical_domain_box = physical_domain.front();
    const hier::IntVector& physical_domain_dims = physical_domain_box.numberCells();
    const hier::IntVector finest_level_dims = physical_domain_dims*ratioFinestLevelToCoarestLevel;
    
    /*
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        double* u_inv_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_reciprocal_quantity.resize(finest_level_dim_0);
        double* u_inv_avg_global = averaged_reciprocal_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_inv_avg_local[i] = double(0);
            u_inv_avg_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_quantity =
                    d_flow_model_tmp->getGlobalCellData(quantity_name);
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        /*
                         * Compute the index of the data point and count how many times the data is repeated.
                         */
                        
                        const hier::Index idx_pt(tbox::Dimension(1), idx_lo_0 + i);
                        
                        int n_overlapped = 1;
                        
                        for (hier::BoxContainer::BoxContainerConstIterator iob(
                                patch_overlapped_visible_boxes.begin());
                             iob != patch_overlapped_visible_boxes.end();
                             iob++)
                        {
                            const hier::Box& patch_overlapped_visible_box = *iob;
                            
                            if (patch_overlapped_visible_box.contains(idx_pt))
                            {
                                n_overlapped++;
                            }
                        }
                        
                        /*
                         * Compute the linear index and the data to add.
                         */
                        
                        const int idx = relative_idx_lo_0 + i + num_ghosts_0_quantity;
                        
                        const double value_to_add = (double(1)/u[idx])/((double) n_overlapped);
                        
                        for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                        {
                            const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                            
                            u_inv_avg_local[idx_fine] += value_to_add;
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average of reciprocal.
         */
        
        mpi.Allreduce(
            u_inv_avg_local,
            u_inv_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_inv_avg_local);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        
        double* u_inv_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_reciprocal_quantity.resize(finest_level_dim_0);
        double* u_inv_avg_global = averaged_reciprocal_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_inv_avg_local[i] = double(0);
            u_inv_avg_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_quantity =
                    d_flow_model_tmp->getGlobalCellData(quantity_name);
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                
                const double weight = dx[1]/L_y;
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            /*
                             * Compute the index of the data point and count how many times the data is repeated.
                             */
                            
                            const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j);
                            
                            int n_overlapped = 1;
                            
                            for (hier::BoxContainer::BoxContainerConstIterator iob(
                                    patch_overlapped_visible_boxes.begin());
                                 iob != patch_overlapped_visible_boxes.end();
                                 iob++)
                            {
                                const hier::Box& patch_overlapped_visible_box = *iob;
                                
                                if (patch_overlapped_visible_box.contains(idx_pt))
                                {
                                    n_overlapped++;
                                }
                            }
                            
                            /*
                             * Compute the linear index and the data to add.
                             */
                            
                            const int idx = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                            
                            const double value_to_add = (double(1)/u[idx])*weight/((double) n_overlapped);
                            
                            for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                            {
                                const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                
                                u_inv_avg_local[idx_fine] += value_to_add;
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average of reciprocal.
         */
        
        mpi.Allreduce(
            u_inv_avg_local,
            u_inv_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_inv_avg_local);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        const double L_z = x_hi[2] - x_lo[2];
        
        double* u_inv_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_reciprocal_quantity.resize(finest_level_dim_0);
        double* u_inv_avg_global = averaged_reciprocal_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_inv_avg_local[i] = double(0);
            u_inv_avg_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_quantity =
                    d_flow_model_tmp->getGlobalCellData(quantity_name);
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int num_ghosts_2_quantity = num_ghosts_quantity[2];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                const int ghostcell_dim_1_quantity = ghostcell_dims_quantity[1];
                
                const double weight = (dx[1]*dx[2])/(L_y*L_z);
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    const int interior_dim_2 = interior_dims[2];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int idx_lo_2 = index_lo[2];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    const int relative_idx_lo_2 = relative_index_lo[2];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                /*
                                 * Compute the index of the data point and count how many times the data is repeated.
                                 */
                                
                                const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j, idx_lo_2 + k);
                                
                                int n_overlapped = 1;
                                
                                for (hier::BoxContainer::BoxContainerConstIterator iob(
                                        patch_overlapped_visible_boxes.begin());
                                     iob != patch_overlapped_visible_boxes.end();
                                     iob++)
                                {
                                    const hier::Box& patch_overlapped_visible_box = *iob;
                                    
                                    if (patch_overlapped_visible_box.contains(idx_pt))
                                    {
                                        n_overlapped++;
                                    }
                                }
                                
                                /*
                                 * Compute the linear index and the data to add.
                                 */
                                
                                const int idx = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                    (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                        ghostcell_dim_1_quantity;
                                
                                const double value_to_add = (double(1)/u[idx])*weight/((double) n_overlapped);
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    u_inv_avg_local[idx_fine] += value_to_add;
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average of reciprocal.
         */
        
        mpi.Allreduce(
            u_inv_avg_local,
            u_inv_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_inv_avg_local);
    }
    
    return averaged_reciprocal_quantity;
}


/*
 * Compute averaged derivative of value with only x direction as inhomogeneous direction.
 */
std::vector<double>
FlowModelStatisticsUtilitiesFourEqnConservative::getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
    const std::string quantity_name,
    const int component_idx,
    const int derivative_direction,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context) const
{
    std::vector<double> averaged_derivative_quantity;
    
    boost::shared_ptr<FlowModel> d_flow_model_tmp = d_flow_model.lock();
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Get the refinement ratio from the finest level to the coarest level.
     */
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    hier::IntVector ratioFinestLevelToCoarestLevel =
        patch_hierarchy->getRatioToCoarserLevel(num_levels - 1);
    for (int li = num_levels - 2; li > 0 ; li--)
    {
        ratioFinestLevelToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(li);
    }
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    boost::shared_ptr<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *patch_hierarchy,
            0,
            num_levels - 1));
    
    /*
     * Get the number of cells of physical domain refined to the finest level.
     */
    
    const hier::BoxContainer& physical_domain = d_grid_geometry->getPhysicalDomain();
    const hier::Box& physical_domain_box = physical_domain.front();
    const hier::IntVector& physical_domain_dims = physical_domain_box.numberCells();
    const hier::IntVector finest_level_dims = physical_domain_dims*ratioFinestLevelToCoarestLevel;
    
    /*
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        double* u_der_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_derivative_quantity.resize(finest_level_dim_0);
        double* u_der_avg_global = averaged_derivative_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_der_avg_local[i] = double(0);
            u_der_avg_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                // HARD CODE TO BE SIXTH ORDER CENTRAL SCHEME FOR DIFFERENTIATION.
                TBOX_ASSERT(num_ghosts >= hier::IntVector::getOne(d_dim)*3);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, num_ghosts));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_quantity =
                    d_flow_model_tmp->getGlobalCellData(quantity_name);
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        /*
                         * Compute the index of the data point and count how many times the data is repeated.
                         */
                        
                        const hier::Index idx_pt(tbox::Dimension(1), idx_lo_0 + i);
                        
                        int n_overlapped = 1;
                        
                        for (hier::BoxContainer::BoxContainerConstIterator iob(
                                patch_overlapped_visible_boxes.begin());
                             iob != patch_overlapped_visible_boxes.end();
                             iob++)
                        {
                            const hier::Box& patch_overlapped_visible_box = *iob;
                            
                            if (patch_overlapped_visible_box.contains(idx_pt))
                            {
                                n_overlapped++;
                            }
                        }
                        
                        /*
                         * Compute the derivative.
                         */
                        
                        double value_to_add = double(0);
                        
                        if (derivative_direction == 0)
                        {
                            // Compute linear indices.
                            const int idx_LLL = relative_idx_lo_0 + (i - 3) + num_ghosts_0_quantity;
                            const int idx_LL  = relative_idx_lo_0 + (i - 2) + num_ghosts_0_quantity;
                            const int idx_L   = relative_idx_lo_0 + (i - 1) + num_ghosts_0_quantity;
                            const int idx_R   = relative_idx_lo_0 + (i + 1) + num_ghosts_0_quantity;
                            const int idx_RR  = relative_idx_lo_0 + (i + 2) + num_ghosts_0_quantity;
                            const int idx_RRR = relative_idx_lo_0 + (i + 3) + num_ghosts_0_quantity;
                            
                            const double dudx = (double(1)/double(60)*(u[idx_RRR] - u[idx_LLL])
                                - double(3)/double(20)*(u[idx_RR] - u[idx_LL])
                                + double(3)/double(4)*(u[idx_R] - u[idx_L]))/dx[0];
                            
                            value_to_add = dudx/((double) n_overlapped);
                        }
                        else
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "Cannot take derivative for one-dimensional problem!\n"
                                << "derivative_direction = " << derivative_direction << " given!\n"
                                << std::endl);
                        }
                        
                        /*
                         * Add the data.
                         */
                        
                        for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                        {
                            const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                            
                            u_der_avg_local[idx_fine] += value_to_add;
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average of derivative.
         */
        
        mpi.Allreduce(
            u_der_avg_local,
            u_der_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_der_avg_local);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        
        double* u_der_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_derivative_quantity.resize(finest_level_dim_0);
        double* u_der_avg_global = averaged_derivative_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_der_avg_local[i] = double(0);
            u_der_avg_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                // HARD CODE TO BE SIXTH ORDER CENTRAL SCHEME FOR DIFFERENTIATION.
                TBOX_ASSERT(num_ghosts >= hier::IntVector::getOne(d_dim)*3);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, num_ghosts));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_quantity =
                    d_flow_model_tmp->getGlobalCellData(quantity_name);
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                
                const double weight = dx[1]/L_y;
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            /*
                             * Compute the index of the data point and count how many times the data is repeated.
                             */
                            
                            const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j);
                            
                            int n_overlapped = 1;
                            
                            for (hier::BoxContainer::BoxContainerConstIterator iob(
                                    patch_overlapped_visible_boxes.begin());
                                 iob != patch_overlapped_visible_boxes.end();
                                 iob++)
                            {
                                const hier::Box& patch_overlapped_visible_box = *iob;
                                
                                if (patch_overlapped_visible_box.contains(idx_pt))
                                {
                                    n_overlapped++;
                                }
                            }
                            
                            /*
                             * Compute the derivative.
                             */
                            
                            double value_to_add = double(0);
                            
                            if (derivative_direction == 0)
                            {
                                const int idx_LLL = (relative_idx_lo_0 + (i - 3) + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const int idx_LL  = (relative_idx_lo_0 + (i - 2) + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const int idx_L   = (relative_idx_lo_0 + (i - 1) + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const int idx_R   = (relative_idx_lo_0 + (i + 1) + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const int idx_RR  = (relative_idx_lo_0 + (i + 2) + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const int idx_RRR = (relative_idx_lo_0 + (i + 3) + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const double dudx = (double(1)/double(60)*(u[idx_RRR] - u[idx_LLL])
                                    - double(3)/double(20)*(u[idx_RR] - u[idx_LL])
                                    + double(3)/double(4)*(u[idx_R] - u[idx_L]))/dx[0];
                                
                                value_to_add = dudx*weight/((double) n_overlapped);
                            }
                            else if (derivative_direction == 1)
                            {
                                const int idx_BBB = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + (j - 3) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const int idx_BB  = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + (j - 2) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const int idx_B   = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + (j - 1) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const int idx_T   = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + (j + 1) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const int idx_TT  = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + (j + 2) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const int idx_TTT = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + (j + 3) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const double dudy = (double(1)/double(60)*(u[idx_TTT] - u[idx_BBB])
                                    - double(3)/double(20)*(u[idx_TT] - u[idx_BB])
                                    + double(3)/double(4)*(u[idx_T] - u[idx_B]))/dx[1];
                                
                                value_to_add = dudy*weight/((double) n_overlapped);
                            }
                            else
                            {
                                TBOX_ERROR(d_object_name
                                    << ": "
                                    << "Cannot take derivative for two-dimensional problem!\n"
                                    << "derivative_direction = " << derivative_direction << " given!\n"
                                    << std::endl);
                            }
                            
                            /*
                             * Add the data.
                             */
                            
                            for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                            {
                                const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                
                                u_der_avg_local[idx_fine] += value_to_add;
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average of derivative.
         */
        
        mpi.Allreduce(
            u_der_avg_local,
            u_der_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_der_avg_local);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        const double L_z = x_hi[2] - x_lo[2];
        
        double* u_der_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_derivative_quantity.resize(finest_level_dim_0);
        double* u_der_avg_global = averaged_derivative_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_der_avg_local[i] = double(0);
            u_der_avg_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                // HARD CODE TO BE SIXTH ORDER CENTRAL SCHEME FOR DIFFERENTIATION.
                TBOX_ASSERT(num_ghosts >= hier::IntVector::getOne(d_dim)*3);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, num_ghosts));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_quantity =
                    d_flow_model_tmp->getGlobalCellData(quantity_name);
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int num_ghosts_2_quantity = num_ghosts_quantity[2];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                const int ghostcell_dim_1_quantity = ghostcell_dims_quantity[1];
                
                const double weight = (dx[1]*dx[2])/(L_y*L_z);
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    const int interior_dim_2 = interior_dims[2];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int idx_lo_2 = index_lo[2];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    const int relative_idx_lo_2 = relative_index_lo[2];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                /*
                                 * Compute the index of the data point and count how many times the data is repeated.
                                 */
                                
                                const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j, idx_lo_2 + k);
                                
                                int n_overlapped = 1;
                                
                                for (hier::BoxContainer::BoxContainerConstIterator iob(
                                        patch_overlapped_visible_boxes.begin());
                                     iob != patch_overlapped_visible_boxes.end();
                                     iob++)
                                {
                                    const hier::Box& patch_overlapped_visible_box = *iob;
                                    
                                    if (patch_overlapped_visible_box.contains(idx_pt))
                                    {
                                        n_overlapped++;
                                    }
                                }
                                
                                /*
                                 * Compute the derivative.
                                 */
                                
                                double value_to_add = double(0);
                                
                                if (derivative_direction == 0)
                                {
                                    const int idx_LLL = (relative_idx_lo_0 + (i - 3) + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_LL  = (relative_idx_lo_0 + (i - 2) + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_L   = (relative_idx_lo_0 + (i - 1) + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_R   = (relative_idx_lo_0 + (i + 1) + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_RR  = (relative_idx_lo_0 + (i + 2) + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_RRR = (relative_idx_lo_0 + (i + 3) + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const double dudx = (double(1)/double(60)*(u[idx_RRR] - u[idx_LLL])
                                        - double(3)/double(20)*(u[idx_RR] - u[idx_LL])
                                        + double(3)/double(4)*(u[idx_R] - u[idx_L]))/dx[0];
                                    
                                    value_to_add = dudx*weight/((double) n_overlapped);
                                }
                                else if (derivative_direction == 1)
                                {
                                    const int idx_BBB = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + (j - 3) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_BB  = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + (j - 2) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_B   = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + (j - 1) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_T   = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + (j + 1) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_TT  = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + (j + 2) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_TTT = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + (j + 3) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const double dudy = (double(1)/double(60)*(u[idx_TTT] - u[idx_BBB])
                                        - double(3)/double(20)*(u[idx_TT] - u[idx_BB])
                                        + double(3)/double(4)*(u[idx_T] - u[idx_B]))/dx[1];
                                    
                                    value_to_add = dudy*weight/((double) n_overlapped);
                                }
                                else if (derivative_direction == 2)
                                {
                                    const int idx_BBB = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + (k - 3) + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_BB  = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + (k - 2) + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_B   = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + (k - 1) + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_F   = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + (k + 1) + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_FF  = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + (k + 2) + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_FFF = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + (k + 3) + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const double dudz = (double(1)/double(60)*(u[idx_FFF] - u[idx_BBB])
                                        - double(3)/double(20)*(u[idx_FF] - u[idx_BB])
                                        + double(3)/double(4)*(u[idx_F] - u[idx_B]))/dx[2];
                                    
                                    value_to_add = dudz*weight/((double) n_overlapped);
                                }
                                else
                                {
                                    TBOX_ERROR(d_object_name
                                        << ": "
                                        << "Cannot take derivative for three-dimensional problem!\n"
                                        << "derivative_direction = " << derivative_direction << " given!\n"
                                        << std::endl);
                                }
                                
                                /*
                                 * Add the data.
                                 */
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    u_der_avg_local[idx_fine] += value_to_add;
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average of derivative.
         */
        
        mpi.Allreduce(
            u_der_avg_local,
            u_der_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_der_avg_local);
    }
    
    return averaged_derivative_quantity;
}


/*
 * Compute averaged derivative of reciprocal of value with only x direction as inhomogeneous direction.
 */
std::vector<double>
FlowModelStatisticsUtilitiesFourEqnConservative::
getAveragedDerivativeOfReciprocalOfQuantityWithInhomogeneousXDirection(
    const std::string quantity_name,
    const int component_idx,
    const int derivative_direction,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context) const
{
    std::vector<double> averaged_derivative_reciprocal_quantity;
    
    boost::shared_ptr<FlowModel> d_flow_model_tmp = d_flow_model.lock();
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Get the refinement ratio from the finest level to the coarest level.
     */
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    hier::IntVector ratioFinestLevelToCoarestLevel =
        patch_hierarchy->getRatioToCoarserLevel(num_levels - 1);
    for (int li = num_levels - 2; li > 0 ; li--)
    {
        ratioFinestLevelToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(li);
    }
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    boost::shared_ptr<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *patch_hierarchy,
            0,
            num_levels - 1));
    
    /*
     * Get the number of cells of physical domain refined to the finest level.
     */
    
    const hier::BoxContainer& physical_domain = d_grid_geometry->getPhysicalDomain();
    const hier::Box& physical_domain_box = physical_domain.front();
    const hier::IntVector& physical_domain_dims = physical_domain_box.numberCells();
    const hier::IntVector finest_level_dims = physical_domain_dims*ratioFinestLevelToCoarestLevel;
    
    /*
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        double* u_inv_der_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_derivative_reciprocal_quantity.resize(finest_level_dim_0);
        double* u_inv_der_avg_global = averaged_derivative_reciprocal_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_inv_der_avg_local[i] = double(0);
            u_inv_der_avg_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                // HARD CODE TO BE SIXTH ORDER CENTRAL SCHEME FOR DIFFERENTIATION.
                TBOX_ASSERT(num_ghosts >= hier::IntVector::getOne(d_dim)*3);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, num_ghosts));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_quantity =
                    d_flow_model_tmp->getGlobalCellData(quantity_name);
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        /*
                         * Compute the index of the data point and count how many times the data is repeated.
                         */
                        
                        const hier::Index idx_pt(tbox::Dimension(1), idx_lo_0 + i);
                        
                        int n_overlapped = 1;
                        
                        for (hier::BoxContainer::BoxContainerConstIterator iob(
                                patch_overlapped_visible_boxes.begin());
                             iob != patch_overlapped_visible_boxes.end();
                             iob++)
                        {
                            const hier::Box& patch_overlapped_visible_box = *iob;
                            
                            if (patch_overlapped_visible_box.contains(idx_pt))
                            {
                                n_overlapped++;
                            }
                        }
                        
                        /*
                         * Compute the derivative.
                         */
                        
                        double value_to_add = double(0);
                        
                        if (derivative_direction == 0)
                        {
                            // Compute linear indices.
                            const int idx_LLL = relative_idx_lo_0 + (i - 3) + num_ghosts_0_quantity;
                            const int idx_LL  = relative_idx_lo_0 + (i - 2) + num_ghosts_0_quantity;
                            const int idx_L   = relative_idx_lo_0 + (i - 1) + num_ghosts_0_quantity;
                            const int idx_R   = relative_idx_lo_0 + (i + 1) + num_ghosts_0_quantity;
                            const int idx_RR  = relative_idx_lo_0 + (i + 2) + num_ghosts_0_quantity;
                            const int idx_RRR = relative_idx_lo_0 + (i + 3) + num_ghosts_0_quantity;
                            
                            const double dudx = (double(1)/double(60)*(double(1)/u[idx_RRR] - double(1)/u[idx_LLL])
                                - double(3)/double(20)*(double(1)/u[idx_RR] - double(1)/u[idx_LL])
                                + double(3)/double(4)*(double(1)/u[idx_R] - double(1)/u[idx_L]))/dx[0];
                            
                            value_to_add = dudx/((double) n_overlapped);
                        }
                        else
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "Cannot take derivative for one-dimensional problem!\n"
                                << "derivative_direction = " << derivative_direction << " given!\n"
                                << std::endl);
                        }
                        
                        /*
                         * Add the data.
                         */
                        
                        for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                        {
                            const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                            
                            u_inv_der_avg_local[idx_fine] += value_to_add;
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average of derivative of reciprocal.
         */
        
        mpi.Allreduce(
            u_inv_der_avg_local,
            u_inv_der_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_inv_der_avg_local);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        
        double* u_inv_der_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_derivative_reciprocal_quantity.resize(finest_level_dim_0);
        double* u_inv_der_avg_global = averaged_derivative_reciprocal_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_inv_der_avg_local[i] = double(0);
            u_inv_der_avg_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                // HARD CODE TO BE SIXTH ORDER CENTRAL SCHEME FOR DIFFERENTIATION.
                TBOX_ASSERT(num_ghosts >= hier::IntVector::getOne(d_dim)*3);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, num_ghosts));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_quantity =
                    d_flow_model_tmp->getGlobalCellData(quantity_name);
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                
                const double weight = dx[1]/L_y;
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            /*
                             * Compute the index of the data point and count how many times the data is repeated.
                             */
                            
                            const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j);
                            
                            int n_overlapped = 1;
                            
                            for (hier::BoxContainer::BoxContainerConstIterator iob(
                                    patch_overlapped_visible_boxes.begin());
                                 iob != patch_overlapped_visible_boxes.end();
                                 iob++)
                            {
                                const hier::Box& patch_overlapped_visible_box = *iob;
                                
                                if (patch_overlapped_visible_box.contains(idx_pt))
                                {
                                    n_overlapped++;
                                }
                            }
                            
                            /*
                             * Compute the derivative.
                             */
                            
                            double value_to_add = double(0);
                            
                            if (derivative_direction == 0)
                            {
                                const int idx_LLL = (relative_idx_lo_0 + (i - 3) + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const int idx_LL  = (relative_idx_lo_0 + (i - 2) + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const int idx_L   = (relative_idx_lo_0 + (i - 1) + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const int idx_R   = (relative_idx_lo_0 + (i + 1) + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const int idx_RR  = (relative_idx_lo_0 + (i + 2) + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const int idx_RRR = (relative_idx_lo_0 + (i + 3) + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const double dudx = (double(1)/double(60)*(double(1)/u[idx_RRR] - double(1)/u[idx_LLL])
                                    - double(3)/double(20)*(double(1)/u[idx_RR] - double(1)/u[idx_LL])
                                    + double(3)/double(4)*(double(1)/u[idx_R] - double(1)/u[idx_L]))/dx[0];
                                
                                value_to_add = dudx*weight/((double) n_overlapped);
                            }
                            else if (derivative_direction == 1)
                            {
                                const int idx_BBB = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + (j - 3) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const int idx_BB  = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + (j - 2) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const int idx_B   = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + (j - 1) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const int idx_T   = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + (j + 1) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const int idx_TT  = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + (j + 2) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const int idx_TTT = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + (j + 3) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const double dudy = (double(1)/double(60)*(double(1)/u[idx_TTT] - double(1)/u[idx_BBB])
                                    - double(3)/double(20)*(double(1)/u[idx_TT] - double(1)/u[idx_BB])
                                    + double(3)/double(4)*(double(1)/u[idx_T] - double(1)/u[idx_B]))/dx[1];
                                
                                value_to_add = dudy*weight/((double) n_overlapped);
                            }
                            else
                            {
                                TBOX_ERROR(d_object_name
                                    << ": "
                                    << "Cannot take derivative for two-dimensional problem!\n"
                                    << "derivative_direction = " << derivative_direction << " given!\n"
                                    << std::endl);
                            }
                            
                            /*
                             * Add the data.
                             */
                            
                            for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                            {
                                const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                
                                u_inv_der_avg_local[idx_fine] += value_to_add;
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average of derivative of reciprocal.
         */
        
        mpi.Allreduce(
            u_inv_der_avg_local,
            u_inv_der_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_inv_der_avg_local);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        const double L_z = x_hi[2] - x_lo[2];
        
        double* u_inv_der_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_derivative_reciprocal_quantity.resize(finest_level_dim_0);
        double* u_inv_der_avg_global = averaged_derivative_reciprocal_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_inv_der_avg_local[i] = double(0);
            u_inv_der_avg_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                // HARD CODE TO BE SIXTH ORDER CENTRAL SCHEME FOR DIFFERENTIATION.
                TBOX_ASSERT(num_ghosts >= hier::IntVector::getOne(d_dim)*3);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, num_ghosts));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_quantity =
                    d_flow_model_tmp->getGlobalCellData(quantity_name);
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int num_ghosts_2_quantity = num_ghosts_quantity[2];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                const int ghostcell_dim_1_quantity = ghostcell_dims_quantity[1];
                
                const double weight = (dx[1]*dx[2])/(L_y*L_z);
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    const int interior_dim_2 = interior_dims[2];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int idx_lo_2 = index_lo[2];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    const int relative_idx_lo_2 = relative_index_lo[2];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                /*
                                 * Compute the index of the data point and count how many times the data is repeated.
                                 */
                                
                                const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j, idx_lo_2 + k);
                                
                                int n_overlapped = 1;
                                
                                for (hier::BoxContainer::BoxContainerConstIterator iob(
                                        patch_overlapped_visible_boxes.begin());
                                     iob != patch_overlapped_visible_boxes.end();
                                     iob++)
                                {
                                    const hier::Box& patch_overlapped_visible_box = *iob;
                                    
                                    if (patch_overlapped_visible_box.contains(idx_pt))
                                    {
                                        n_overlapped++;
                                    }
                                }
                                
                                /*
                                 * Compute the derivative.
                                 */
                                
                                double value_to_add = double(0);
                                
                                if (derivative_direction == 0)
                                {
                                    const int idx_LLL = (relative_idx_lo_0 + (i - 3) + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_LL  = (relative_idx_lo_0 + (i - 2) + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_L   = (relative_idx_lo_0 + (i - 1) + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_R   = (relative_idx_lo_0 + (i + 1) + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_RR  = (relative_idx_lo_0 + (i + 2) + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_RRR = (relative_idx_lo_0 + (i + 3) + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const double dudx = (double(1)/double(60)*(double(1)/u[idx_RRR] - double(1)/u[idx_LLL])
                                        - double(3)/double(20)*(double(1)/u[idx_RR] - double(1)/u[idx_LL])
                                        + double(3)/double(4)*(double(1)/u[idx_R] - double(1)/u[idx_L]))/dx[0];
                                    
                                    value_to_add = dudx*weight/((double) n_overlapped);
                                }
                                else if (derivative_direction == 1)
                                {
                                    const int idx_BBB = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + (j - 3) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_BB  = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + (j - 2) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_B   = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + (j - 1) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_T   = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + (j + 1) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_TT  = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + (j + 2) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_TTT = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + (j + 3) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const double dudy = (double(1)/double(60)*(double(1)/u[idx_TTT] - double(1)/u[idx_BBB])
                                        - double(3)/double(20)*(double(1)/u[idx_TT] - double(1)/u[idx_BB])
                                        + double(3)/double(4)*(double(1)/u[idx_T] - double(1)/u[idx_B]))/dx[1];
                                    
                                    value_to_add = dudy*weight/((double) n_overlapped);
                                }
                                else if (derivative_direction == 2)
                                {
                                    const int idx_BBB = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + (k - 3) + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_BB  = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + (k - 2) + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_B   = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + (k - 1) + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_F   = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + (k + 1) + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_FF  = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + (k + 2) + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_FFF = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + (k + 3) + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const double dudz = (double(1)/double(60)*(double(1)/u[idx_FFF] - double(1)/u[idx_BBB])
                                        - double(3)/double(20)*(double(1)/u[idx_FF] - double(1)/u[idx_BB])
                                        + double(3)/double(4)*(double(1)/u[idx_F] - double(1)/u[idx_B]))/dx[2];
                                    
                                    value_to_add = dudz*weight/((double) n_overlapped);
                                }
                                else
                                {
                                    TBOX_ERROR(d_object_name
                                        << ": "
                                        << "Cannot take derivative for three-dimensional problem!\n"
                                        << "derivative_direction = " << derivative_direction << " given!\n"
                                        << std::endl);
                                }
                                
                                /*
                                 * Add the data.
                                 */
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    u_inv_der_avg_local[idx_fine] += value_to_add;
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average of derivative of reciprocal.
         */
        
        mpi.Allreduce(
            u_inv_der_avg_local,
            u_inv_der_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_inv_der_avg_local);
    }
    
    return averaged_derivative_reciprocal_quantity;
}


/*
 * Compute averaged derivative of value (on product of variables) with only x direction as inhomogeneous direction.
 */
std::vector<double>
FlowModelStatisticsUtilitiesFourEqnConservative::getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
    const std::vector<std::string>& quantity_names,
    const std::vector<int>& component_indices,
    const int derivative_direction,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context) const
{
    int num_quantities = static_cast<int>(quantity_names.size());
    
    TBOX_ASSERT(static_cast<int>(component_indices.size()) == num_quantities);
    
    std::vector<double> averaged_derivative_quantity;
    
    boost::shared_ptr<FlowModel> d_flow_model_tmp = d_flow_model.lock();
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Get the refinement ratio from the finest level to the coarest level.
     */
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    hier::IntVector ratioFinestLevelToCoarestLevel =
        patch_hierarchy->getRatioToCoarserLevel(num_levels - 1);
    for (int li = num_levels - 2; li > 0 ; li--)
    {
        ratioFinestLevelToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(li);
    }
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    boost::shared_ptr<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *patch_hierarchy,
            0,
            num_levels - 1));
    
    /*
     * Get the number of cells of physical domain refined to the finest level.
     */
    
    const hier::BoxContainer& physical_domain = d_grid_geometry->getPhysicalDomain();
    const hier::Box& physical_domain_box = physical_domain.front();
    const hier::IntVector& physical_domain_dims = physical_domain_box.numberCells();
    const hier::IntVector finest_level_dims = physical_domain_dims*ratioFinestLevelToCoarestLevel;
    
    /*
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        double* u_der_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_derivative_quantity.resize(finest_level_dim_0);
        double* u_der_avg_global = averaged_derivative_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_der_avg_local[i] = double(0);
            u_der_avg_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and the quantities in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                // HARD CODE TO BE SIXTH ORDER CENTRAL SCHEME FOR DIFFERENTIATION.
                TBOX_ASSERT(num_ghosts >= hier::IntVector::getOne(d_dim)*3);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], num_ghosts));
                }
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<boost::shared_ptr<pdat::CellData<double> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = d_flow_model_tmp->getGlobalCellData(quantity_names[qi]);
                }
                
                std::vector<double*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                }
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        /*
                         * Compute the index of the data point and count how many times the data is repeated.
                         */
                        
                        const hier::Index idx_pt(tbox::Dimension(1), idx_lo_0 + i);
                        
                        int n_overlapped = 1;
                        
                        for (hier::BoxContainer::BoxContainerConstIterator iob(
                                patch_overlapped_visible_boxes.begin());
                             iob != patch_overlapped_visible_boxes.end();
                             iob++)
                        {
                            const hier::Box& patch_overlapped_visible_box = *iob;
                            
                            if (patch_overlapped_visible_box.contains(idx_pt))
                            {
                                n_overlapped++;
                            }
                        }
                        
                        /*
                         * Compute the derivative.
                         */
                        
                        double value_to_add = double(0);
                        
                        if (derivative_direction == 0)
                        {
                            double u_LLL = double(1);
                            double u_LL  = double(1);
                            double u_L   = double(1);
                            double u_R   = double(1);
                            double u_RR  = double(1);
                            double u_RRR = double(1);
                            
                            for (int qi = 0; qi < num_quantities; qi++)
                            {
                                // Compute linear indices.
                                const int idx_LLL = relative_idx_lo_0 + (i - 3) + num_ghosts_0_u_qi[qi];
                                const int idx_LL  = relative_idx_lo_0 + (i - 2) + num_ghosts_0_u_qi[qi];
                                const int idx_L   = relative_idx_lo_0 + (i - 1) + num_ghosts_0_u_qi[qi];
                                const int idx_R   = relative_idx_lo_0 + (i + 1) + num_ghosts_0_u_qi[qi];
                                const int idx_RR  = relative_idx_lo_0 + (i + 2) + num_ghosts_0_u_qi[qi];
                                const int idx_RRR = relative_idx_lo_0 + (i + 3) + num_ghosts_0_u_qi[qi];
                                
                                u_LLL *= u_qi[qi][idx_LLL];
                                u_LL  *= u_qi[qi][idx_LL ];
                                u_L   *= u_qi[qi][idx_L  ];
                                u_R   *= u_qi[qi][idx_R  ];
                                u_RR  *= u_qi[qi][idx_RR ];
                                u_RRR *= u_qi[qi][idx_RRR];
                            }
                            
                            const double dudx = (double(1)/double(60)*(u_RRR - u_LLL)
                                - double(3)/double(20)*(u_RR - u_LL)
                                + double(3)/double(4)*(u_R - u_L))/dx[0];
                            
                            value_to_add = dudx/((double) n_overlapped);
                        }
                        else
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "Cannot take derivative for one-dimensional problem!\n"
                                << "derivative_direction = " << derivative_direction << " given!\n"
                                << std::endl);
                        }
                        
                        /*
                         * Add the data.
                         */
                        
                        for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                        {
                            const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                            
                            u_der_avg_local[idx_fine] += value_to_add;
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average of derivative.
         */
        
        mpi.Allreduce(
            u_der_avg_local,
            u_der_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_der_avg_local);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        
        double* u_der_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_derivative_quantity.resize(finest_level_dim_0);
        double* u_der_avg_global = averaged_derivative_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_der_avg_local[i] = double(0);
            u_der_avg_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and the quantities in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                // HARD CODE TO BE SIXTH ORDER CENTRAL SCHEME FOR DIFFERENTIATION.
                TBOX_ASSERT(num_ghosts >= hier::IntVector::getOne(d_dim)*3);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], num_ghosts));
                }
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<boost::shared_ptr<pdat::CellData<double> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = d_flow_model_tmp->getGlobalCellData(quantity_names[qi]);
                }
                
                std::vector<double*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                }
                
                const double weight = dx[1]/L_y;
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            /*
                             * Compute the index of the data point and count how many times the data is repeated.
                             */
                            
                            const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j);
                            
                            int n_overlapped = 1;
                            
                            for (hier::BoxContainer::BoxContainerConstIterator iob(
                                    patch_overlapped_visible_boxes.begin());
                                 iob != patch_overlapped_visible_boxes.end();
                                 iob++)
                            {
                                const hier::Box& patch_overlapped_visible_box = *iob;
                                
                                if (patch_overlapped_visible_box.contains(idx_pt))
                                {
                                    n_overlapped++;
                                }
                            }
                            
                            /*
                             * Compute the derivative.
                             */
                            
                            double value_to_add = double(0);
                            
                            if (derivative_direction == 0)
                            {
                                double u_LLL = double(1);
                                double u_LL  = double(1);
                                double u_L   = double(1);
                                double u_R   = double(1);
                                double u_RR  = double(1);
                                double u_RRR = double(1);
                                
                                for (int qi = 0; qi < num_quantities; qi++)
                                {
                                    // Compute linear indices.
                                    const int idx_LLL = (relative_idx_lo_0 + (i - 3) + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    const int idx_LL  = (relative_idx_lo_0 + (i - 2) + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    const int idx_L   = (relative_idx_lo_0 + (i - 1) + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    const int idx_R   = (relative_idx_lo_0 + (i + 1) + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    const int idx_RR  = (relative_idx_lo_0 + (i + 2) + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    const int idx_RRR = (relative_idx_lo_0 + (i + 3) + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    u_LLL *= u_qi[qi][idx_LLL];
                                    u_LL  *= u_qi[qi][idx_LL ];
                                    u_L   *= u_qi[qi][idx_L  ];
                                    u_R   *= u_qi[qi][idx_R  ];
                                    u_RR  *= u_qi[qi][idx_RR ];
                                    u_RRR *= u_qi[qi][idx_RRR];
                                }
                                
                                const double dudx = (double(1)/double(60)*(u_RRR - u_LLL)
                                    - double(3)/double(20)*(u_RR - u_LL)
                                    + double(3)/double(4)*(u_R - u_L))/dx[0];
                                
                                value_to_add = dudx*weight/((double) n_overlapped);
                            }
                            else if (derivative_direction == 1)
                            {
                                double u_BBB = double(1);
                                double u_BB  = double(1);
                                double u_B   = double(1);
                                double u_T   = double(1);
                                double u_TT  = double(1);
                                double u_TTT = double(1);
                                
                                for (int qi = 0; qi < num_quantities; qi++)
                                {
                                    // Compute linear indices.
                                    const int idx_BBB = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + (j - 3) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    const int idx_BB  = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + (j - 2) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    const int idx_B   = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + (j - 1) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    const int idx_T   = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + (j + 1) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    const int idx_TT  = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + (j + 2) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    const int idx_TTT = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + (j + 3) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    u_BBB *= u_qi[qi][idx_BBB];
                                    u_BB  *= u_qi[qi][idx_BB ];
                                    u_B   *= u_qi[qi][idx_B  ];
                                    u_T   *= u_qi[qi][idx_T  ];
                                    u_TT  *= u_qi[qi][idx_TT ];
                                    u_TTT *= u_qi[qi][idx_TTT];
                                }
                                
                                const double dudy = (double(1)/double(60)*(u_TTT - u_BBB)
                                    - double(3)/double(20)*(u_TT - u_BB)
                                    + double(3)/double(4)*(u_T - u_B))/dx[1];
                                
                                value_to_add = dudy*weight/((double) n_overlapped);
                            }
                            else
                            {
                                TBOX_ERROR(d_object_name
                                    << ": "
                                    << "Cannot take derivative for two-dimensional problem!\n"
                                    << "derivative_direction = " << derivative_direction << " given!\n"
                                    << std::endl);
                            }
                            
                            /*
                             * Add the data.
                             */
                            
                            for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                            {
                                const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                
                                u_der_avg_local[idx_fine] += value_to_add;
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average of derivative.
         */
        
        mpi.Allreduce(
            u_der_avg_local,
            u_der_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_der_avg_local);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        const double L_z = x_hi[2] - x_lo[2];
        
        double* u_der_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_derivative_quantity.resize(finest_level_dim_0);
        double* u_der_avg_global = averaged_derivative_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_der_avg_local[i] = double(0);
            u_der_avg_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and the quantities in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                // HARD CODE TO BE SIXTH ORDER CENTRAL SCHEME FOR DIFFERENTIATION.
                TBOX_ASSERT(num_ghosts >= hier::IntVector::getOne(d_dim)*3);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], num_ghosts));
                }
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<boost::shared_ptr<pdat::CellData<double> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = d_flow_model_tmp->getGlobalCellData(quantity_names[qi]);
                }
                
                std::vector<double*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> num_ghosts_2_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                std::vector<int> ghostcell_dim_1_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                num_ghosts_2_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                ghostcell_dim_1_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    num_ghosts_2_u_qi.push_back(num_ghosts_u_qi[2]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                    ghostcell_dim_1_u_qi.push_back(ghostcell_dims_u_qi[1]);
                }
                
                const double weight = (dx[1]*dx[2])/(L_y*L_z);
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    const int interior_dim_2 = interior_dims[2];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int idx_lo_2 = index_lo[2];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    const int relative_idx_lo_2 = relative_index_lo[2];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                /*
                                 * Compute the index of the data point and count how many times the data is repeated.
                                 */
                                
                                const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j, idx_lo_2 + k);
                                
                                int n_overlapped = 1;
                                
                                for (hier::BoxContainer::BoxContainerConstIterator iob(
                                        patch_overlapped_visible_boxes.begin());
                                     iob != patch_overlapped_visible_boxes.end();
                                     iob++)
                                {
                                    const hier::Box& patch_overlapped_visible_box = *iob;
                                    
                                    if (patch_overlapped_visible_box.contains(idx_pt))
                                    {
                                        n_overlapped++;
                                    }
                                }
                                
                                /*
                                 * Compute the derivative.
                                 */
                                
                                double value_to_add = double(0);
                                
                                if (derivative_direction == 0)
                                {
                                    double u_LLL = double(1);
                                    double u_LL  = double(1);
                                    double u_L   = double(1);
                                    double u_R   = double(1);
                                    double u_RR  = double(1);
                                    double u_RRR = double(1);
                                    
                                    for (int qi = 0; qi < num_quantities; qi++)
                                    {
                                        // Compute linear indices.
                                        const int idx_LLL = (relative_idx_lo_0 + (i - 3) + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_LL  = (relative_idx_lo_0 + (i - 2) + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_L   = (relative_idx_lo_0 + (i - 1) + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_R   = (relative_idx_lo_0 + (i + 1) + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_RR  = (relative_idx_lo_0 + (i + 2) + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_RRR = (relative_idx_lo_0 + (i + 3) + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        u_LLL *= u_qi[qi][idx_LLL];
                                        u_LL  *= u_qi[qi][idx_LL ];
                                        u_L   *= u_qi[qi][idx_L  ];
                                        u_R   *= u_qi[qi][idx_R  ];
                                        u_RR  *= u_qi[qi][idx_RR ];
                                        u_RRR *= u_qi[qi][idx_RRR];
                                    }
                                    
                                    const double dudx = (double(1)/double(60)*(u_RRR - u_LLL)
                                        - double(3)/double(20)*(u_RR - u_LL)
                                        + double(3)/double(4)*(u_R - u_L))/dx[0];
                                    
                                    value_to_add = dudx*weight/((double) n_overlapped);
                                }
                                else if (derivative_direction == 1)
                                {
                                    double u_BBB = double(1);
                                    double u_BB  = double(1);
                                    double u_B   = double(1);
                                    double u_T   = double(1);
                                    double u_TT  = double(1);
                                    double u_TTT = double(1);
                                    
                                    for (int qi = 0; qi < num_quantities; qi++)
                                    {
                                        // Compute linear indices.
                                        const int idx_BBB = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + (j - 3) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_BB  = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + (j - 2) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_B   = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + (j - 1) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_T   = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + (j + 1) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_TT  = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + (j + 2) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_TTT = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + (j + 3) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        u_BBB *= u_qi[qi][idx_BBB];
                                        u_BB  *= u_qi[qi][idx_BB ];
                                        u_B   *= u_qi[qi][idx_B  ];
                                        u_T   *= u_qi[qi][idx_T  ];
                                        u_TT  *= u_qi[qi][idx_TT ];
                                        u_TTT *= u_qi[qi][idx_TTT];
                                    }
                                    
                                    const double dudy = (double(1)/double(60)*(u_TTT - u_BBB)
                                        - double(3)/double(20)*(u_TT - u_BB)
                                        + double(3)/double(4)*(u_T - u_B))/dx[1];
                                    
                                    value_to_add = dudy*weight/((double) n_overlapped);
                                }
                                else if (derivative_direction == 2)
                                {
                                    double u_BBB = double(1);
                                    double u_BB  = double(1);
                                    double u_B   = double(1);
                                    double u_F   = double(1);
                                    double u_FF  = double(1);
                                    double u_FFF = double(1);
                                    
                                    for (int qi = 0; qi < num_quantities; qi++)
                                    {
                                        // Compute linear indices.
                                        const int idx_BBB = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + (k - 3) + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_BB  = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + (k - 2) + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_B   = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + (k - 1) + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_F   = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + (k + 1) + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_FF  = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + (k + 2) + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_FFF = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + (k + 3) + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        u_BBB *= u_qi[qi][idx_BBB];
                                        u_BB  *= u_qi[qi][idx_BB ];
                                        u_B   *= u_qi[qi][idx_B  ];
                                        u_F   *= u_qi[qi][idx_F  ];
                                        u_FF  *= u_qi[qi][idx_FF ];
                                        u_FFF *= u_qi[qi][idx_FFF];
                                    }
                                    
                                    const double dudz = (double(1)/double(60)*(u_FFF - u_BBB)
                                        - double(3)/double(20)*(u_FF - u_BB)
                                        + double(3)/double(4)*(u_F - u_B))/dx[2];
                                    
                                    value_to_add = dudz*weight/((double) n_overlapped);
                                }
                                else
                                {
                                    TBOX_ERROR(d_object_name
                                        << ": "
                                        << "Cannot take derivative for three-dimensional problem!\n"
                                        << "derivative_direction = " << derivative_direction << " given!\n"
                                        << std::endl);
                                }
                                
                                /*
                                 * Add the data.
                                 */
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    u_der_avg_local[idx_fine] += value_to_add;
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average of derivative.
         */
        
        mpi.Allreduce(
            u_der_avg_local,
            u_der_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_der_avg_local);
    }
    
    return averaged_derivative_quantity;
}


/*
 * Compute averaged derivative of value (on product of variables) with only x direction as inhomogeneous direction.
 */
std::vector<double>
FlowModelStatisticsUtilitiesFourEqnConservative::getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
    const std::vector<std::string>& quantity_names,
    const std::vector<int>& component_indices,
    const std::vector<bool>& use_reciprocal,
    const int derivative_direction,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context) const
{
    int num_quantities = static_cast<int>(quantity_names.size());
    
    TBOX_ASSERT(static_cast<int>(component_indices.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(use_reciprocal.size()) == num_quantities);
    
    std::vector<double> averaged_derivative_quantity;
    
    boost::shared_ptr<FlowModel> d_flow_model_tmp = d_flow_model.lock();
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Get the refinement ratio from the finest level to the coarest level.
     */
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    hier::IntVector ratioFinestLevelToCoarestLevel =
        patch_hierarchy->getRatioToCoarserLevel(num_levels - 1);
    for (int li = num_levels - 2; li > 0 ; li--)
    {
        ratioFinestLevelToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(li);
    }
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    boost::shared_ptr<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *patch_hierarchy,
            0,
            num_levels - 1));
    
    /*
     * Get the number of cells of physical domain refined to the finest level.
     */
    
    const hier::BoxContainer& physical_domain = d_grid_geometry->getPhysicalDomain();
    const hier::Box& physical_domain_box = physical_domain.front();
    const hier::IntVector& physical_domain_dims = physical_domain_box.numberCells();
    const hier::IntVector finest_level_dims = physical_domain_dims*ratioFinestLevelToCoarestLevel;
    
    /*
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        double* u_der_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_derivative_quantity.resize(finest_level_dim_0);
        double* u_der_avg_global = averaged_derivative_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_der_avg_local[i] = double(0);
            u_der_avg_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and the quantities in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                // HARD CODE TO BE SIXTH ORDER CENTRAL SCHEME FOR DIFFERENTIATION.
                TBOX_ASSERT(num_ghosts >= hier::IntVector::getOne(d_dim)*3);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], num_ghosts));
                }
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<boost::shared_ptr<pdat::CellData<double> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = d_flow_model_tmp->getGlobalCellData(quantity_names[qi]);
                }
                
                std::vector<double*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                }
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        /*
                         * Compute the index of the data point and count how many times the data is repeated.
                         */
                        
                        const hier::Index idx_pt(tbox::Dimension(1), idx_lo_0 + i);
                        
                        int n_overlapped = 1;
                        
                        for (hier::BoxContainer::BoxContainerConstIterator iob(
                                patch_overlapped_visible_boxes.begin());
                             iob != patch_overlapped_visible_boxes.end();
                             iob++)
                        {
                            const hier::Box& patch_overlapped_visible_box = *iob;
                            
                            if (patch_overlapped_visible_box.contains(idx_pt))
                            {
                                n_overlapped++;
                            }
                        }
                        
                        /*
                         * Compute the derivative.
                         */
                        
                        double value_to_add = double(0);
                        
                        if (derivative_direction == 0)
                        {
                            double u_LLL = double(1);
                            double u_LL  = double(1);
                            double u_L   = double(1);
                            double u_R   = double(1);
                            double u_RR  = double(1);
                            double u_RRR = double(1);
                            
                            for (int qi = 0; qi < num_quantities; qi++)
                            {
                                // Compute linear indices.
                                const int idx_LLL = relative_idx_lo_0 + (i - 3) + num_ghosts_0_u_qi[qi];
                                const int idx_LL  = relative_idx_lo_0 + (i - 2) + num_ghosts_0_u_qi[qi];
                                const int idx_L   = relative_idx_lo_0 + (i - 1) + num_ghosts_0_u_qi[qi];
                                const int idx_R   = relative_idx_lo_0 + (i + 1) + num_ghosts_0_u_qi[qi];
                                const int idx_RR  = relative_idx_lo_0 + (i + 2) + num_ghosts_0_u_qi[qi];
                                const int idx_RRR = relative_idx_lo_0 + (i + 3) + num_ghosts_0_u_qi[qi];
                                
                                if (use_reciprocal[qi])
                                {
                                    u_LLL /= u_qi[qi][idx_LLL];
                                    u_LL  /= u_qi[qi][idx_LL ];
                                    u_L   /= u_qi[qi][idx_L  ];
                                    u_R   /= u_qi[qi][idx_R  ];
                                    u_RR  /= u_qi[qi][idx_RR ];
                                    u_RRR /= u_qi[qi][idx_RRR];
                                }
                                else
                                {
                                    u_LLL *= u_qi[qi][idx_LLL];
                                    u_LL  *= u_qi[qi][idx_LL ];
                                    u_L   *= u_qi[qi][idx_L  ];
                                    u_R   *= u_qi[qi][idx_R  ];
                                    u_RR  *= u_qi[qi][idx_RR ];
                                    u_RRR *= u_qi[qi][idx_RRR];
                                }
                            }
                            
                            const double dudx = (double(1)/double(60)*(u_RRR - u_LLL)
                                - double(3)/double(20)*(u_RR - u_LL)
                                + double(3)/double(4)*(u_R - u_L))/dx[0];
                            
                            value_to_add = dudx/((double) n_overlapped);
                        }
                        else
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "Cannot take derivative for one-dimensional problem!\n"
                                << "derivative_direction = " << derivative_direction << " given!\n"
                                << std::endl);
                        }
                        
                        /*
                         * Add the data.
                         */
                        
                        for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                        {
                            const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                            
                            u_der_avg_local[idx_fine] += value_to_add;
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average of derivative.
         */
        
        mpi.Allreduce(
            u_der_avg_local,
            u_der_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_der_avg_local);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        
        double* u_der_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_derivative_quantity.resize(finest_level_dim_0);
        double* u_der_avg_global = averaged_derivative_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_der_avg_local[i] = double(0);
            u_der_avg_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and the quantities in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                // HARD CODE TO BE SIXTH ORDER CENTRAL SCHEME FOR DIFFERENTIATION.
                TBOX_ASSERT(num_ghosts >= hier::IntVector::getOne(d_dim)*3);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], num_ghosts));
                }
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<boost::shared_ptr<pdat::CellData<double> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = d_flow_model_tmp->getGlobalCellData(quantity_names[qi]);
                }
                
                std::vector<double*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                }
                
                const double weight = dx[1]/L_y;
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            /*
                             * Compute the index of the data point and count how many times the data is repeated.
                             */
                            
                            const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j);
                            
                            int n_overlapped = 1;
                            
                            for (hier::BoxContainer::BoxContainerConstIterator iob(
                                    patch_overlapped_visible_boxes.begin());
                                 iob != patch_overlapped_visible_boxes.end();
                                 iob++)
                            {
                                const hier::Box& patch_overlapped_visible_box = *iob;
                                
                                if (patch_overlapped_visible_box.contains(idx_pt))
                                {
                                    n_overlapped++;
                                }
                            }
                            
                            /*
                             * Compute the derivative.
                             */
                            
                            double value_to_add = double(0);
                            
                            if (derivative_direction == 0)
                            {
                                double u_LLL = double(1);
                                double u_LL  = double(1);
                                double u_L   = double(1);
                                double u_R   = double(1);
                                double u_RR  = double(1);
                                double u_RRR = double(1);
                                
                                for (int qi = 0; qi < num_quantities; qi++)
                                {
                                    // Compute linear indices.
                                    const int idx_LLL = (relative_idx_lo_0 + (i - 3) + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    const int idx_LL  = (relative_idx_lo_0 + (i - 2) + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    const int idx_L   = (relative_idx_lo_0 + (i - 1) + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    const int idx_R   = (relative_idx_lo_0 + (i + 1) + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    const int idx_RR  = (relative_idx_lo_0 + (i + 2) + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    const int idx_RRR = (relative_idx_lo_0 + (i + 3) + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    if (use_reciprocal[qi])
                                    {
                                        u_LLL /= u_qi[qi][idx_LLL];
                                        u_LL  /= u_qi[qi][idx_LL ];
                                        u_L   /= u_qi[qi][idx_L  ];
                                        u_R   /= u_qi[qi][idx_R  ];
                                        u_RR  /= u_qi[qi][idx_RR ];
                                        u_RRR /= u_qi[qi][idx_RRR];
                                    }
                                    else
                                    {
                                        u_LLL *= u_qi[qi][idx_LLL];
                                        u_LL  *= u_qi[qi][idx_LL ];
                                        u_L   *= u_qi[qi][idx_L  ];
                                        u_R   *= u_qi[qi][idx_R  ];
                                        u_RR  *= u_qi[qi][idx_RR ];
                                        u_RRR *= u_qi[qi][idx_RRR];
                                    }
                                }
                                
                                const double dudx = (double(1)/double(60)*(u_RRR - u_LLL)
                                    - double(3)/double(20)*(u_RR - u_LL)
                                    + double(3)/double(4)*(u_R - u_L))/dx[0];
                                
                                value_to_add = dudx*weight/((double) n_overlapped);
                            }
                            else if (derivative_direction == 1)
                            {
                                double u_BBB = double(1);
                                double u_BB  = double(1);
                                double u_B   = double(1);
                                double u_T   = double(1);
                                double u_TT  = double(1);
                                double u_TTT = double(1);
                                
                                for (int qi = 0; qi < num_quantities; qi++)
                                {
                                    // Compute linear indices.
                                    const int idx_BBB = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + (j - 3) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    const int idx_BB  = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + (j - 2) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    const int idx_B   = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + (j - 1) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    const int idx_T   = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + (j + 1) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    const int idx_TT  = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + (j + 2) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    const int idx_TTT = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + (j + 3) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    if (use_reciprocal[qi])
                                    {
                                        u_BBB /= u_qi[qi][idx_BBB];
                                        u_BB  /= u_qi[qi][idx_BB ];
                                        u_B   /= u_qi[qi][idx_B  ];
                                        u_T   /= u_qi[qi][idx_T  ];
                                        u_TT  /= u_qi[qi][idx_TT ];
                                        u_TTT /= u_qi[qi][idx_TTT];
                                    }
                                    else
                                    {
                                        u_BBB *= u_qi[qi][idx_BBB];
                                        u_BB  *= u_qi[qi][idx_BB ];
                                        u_B   *= u_qi[qi][idx_B  ];
                                        u_T   *= u_qi[qi][idx_T  ];
                                        u_TT  *= u_qi[qi][idx_TT ];
                                        u_TTT *= u_qi[qi][idx_TTT];
                                    }
                                }
                                
                                const double dudy = (double(1)/double(60)*(u_TTT - u_BBB)
                                    - double(3)/double(20)*(u_TT - u_BB)
                                    + double(3)/double(4)*(u_T - u_B))/dx[1];
                                
                                value_to_add = dudy*weight/((double) n_overlapped);
                            }
                            else
                            {
                                TBOX_ERROR(d_object_name
                                    << ": "
                                    << "Cannot take derivative for two-dimensional problem!\n"
                                    << "derivative_direction = " << derivative_direction << " given!\n"
                                    << std::endl);
                            }
                            
                            /*
                             * Add the data.
                             */
                            
                            for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                            {
                                const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                
                                u_der_avg_local[idx_fine] += value_to_add;
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average of derivative.
         */
        
        mpi.Allreduce(
            u_der_avg_local,
            u_der_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_der_avg_local);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        const double L_z = x_hi[2] - x_lo[2];
        
        double* u_der_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_derivative_quantity.resize(finest_level_dim_0);
        double* u_der_avg_global = averaged_derivative_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_der_avg_local[i] = double(0);
            u_der_avg_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and the quantities in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                // HARD CODE TO BE SIXTH ORDER CENTRAL SCHEME FOR DIFFERENTIATION.
                TBOX_ASSERT(num_ghosts >= hier::IntVector::getOne(d_dim)*3);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], num_ghosts));
                }
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<boost::shared_ptr<pdat::CellData<double> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = d_flow_model_tmp->getGlobalCellData(quantity_names[qi]);
                }
                
                std::vector<double*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> num_ghosts_2_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                std::vector<int> ghostcell_dim_1_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                num_ghosts_2_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                ghostcell_dim_1_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    num_ghosts_2_u_qi.push_back(num_ghosts_u_qi[2]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                    ghostcell_dim_1_u_qi.push_back(ghostcell_dims_u_qi[1]);
                }
                
                const double weight = (dx[1]*dx[2])/(L_y*L_z);
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    const int interior_dim_2 = interior_dims[2];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int idx_lo_2 = index_lo[2];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    const int relative_idx_lo_2 = relative_index_lo[2];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                /*
                                 * Compute the index of the data point and count how many times the data is repeated.
                                 */
                                
                                const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j, idx_lo_2 + k);
                                
                                int n_overlapped = 1;
                                
                                for (hier::BoxContainer::BoxContainerConstIterator iob(
                                        patch_overlapped_visible_boxes.begin());
                                     iob != patch_overlapped_visible_boxes.end();
                                     iob++)
                                {
                                    const hier::Box& patch_overlapped_visible_box = *iob;
                                    
                                    if (patch_overlapped_visible_box.contains(idx_pt))
                                    {
                                        n_overlapped++;
                                    }
                                }
                                
                                /*
                                 * Compute the derivative.
                                 */
                                
                                double value_to_add = double(0);
                                
                                if (derivative_direction == 0)
                                {
                                    double u_LLL = double(1);
                                    double u_LL  = double(1);
                                    double u_L   = double(1);
                                    double u_R   = double(1);
                                    double u_RR  = double(1);
                                    double u_RRR = double(1);
                                    
                                    for (int qi = 0; qi < num_quantities; qi++)
                                    {
                                        // Compute linear indices.
                                        const int idx_LLL = (relative_idx_lo_0 + (i - 3) + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_LL  = (relative_idx_lo_0 + (i - 2) + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_L   = (relative_idx_lo_0 + (i - 1) + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_R   = (relative_idx_lo_0 + (i + 1) + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_RR  = (relative_idx_lo_0 + (i + 2) + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_RRR = (relative_idx_lo_0 + (i + 3) + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        if (use_reciprocal[qi])
                                        {
                                            u_LLL /= u_qi[qi][idx_LLL];
                                            u_LL  /= u_qi[qi][idx_LL ];
                                            u_L   /= u_qi[qi][idx_L  ];
                                            u_R   /= u_qi[qi][idx_R  ];
                                            u_RR  /= u_qi[qi][idx_RR ];
                                            u_RRR /= u_qi[qi][idx_RRR];
                                        }
                                        else
                                        {
                                            u_LLL *= u_qi[qi][idx_LLL];
                                            u_LL  *= u_qi[qi][idx_LL ];
                                            u_L   *= u_qi[qi][idx_L  ];
                                            u_R   *= u_qi[qi][idx_R  ];
                                            u_RR  *= u_qi[qi][idx_RR ];
                                            u_RRR *= u_qi[qi][idx_RRR];
                                        }
                                    }
                                    
                                    const double dudx = (double(1)/double(60)*(u_RRR - u_LLL)
                                        - double(3)/double(20)*(u_RR - u_LL)
                                        + double(3)/double(4)*(u_R - u_L))/dx[0];
                                    
                                    value_to_add = dudx*weight/((double) n_overlapped);
                                }
                                else if (derivative_direction == 1)
                                {
                                    double u_BBB = double(1);
                                    double u_BB  = double(1);
                                    double u_B   = double(1);
                                    double u_T   = double(1);
                                    double u_TT  = double(1);
                                    double u_TTT = double(1);
                                    
                                    for (int qi = 0; qi < num_quantities; qi++)
                                    {
                                        // Compute linear indices.
                                        const int idx_BBB = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + (j - 3) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_BB  = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + (j - 2) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_B   = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + (j - 1) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_T   = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + (j + 1) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_TT  = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + (j + 2) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_TTT = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + (j + 3) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        if (use_reciprocal[qi])
                                        {
                                            u_BBB /= u_qi[qi][idx_BBB];
                                            u_BB  /= u_qi[qi][idx_BB ];
                                            u_B   /= u_qi[qi][idx_B  ];
                                            u_T   /= u_qi[qi][idx_T  ];
                                            u_TT  /= u_qi[qi][idx_TT ];
                                            u_TTT /= u_qi[qi][idx_TTT];
                                        }
                                        else
                                        {
                                            u_BBB *= u_qi[qi][idx_BBB];
                                            u_BB  *= u_qi[qi][idx_BB ];
                                            u_B   *= u_qi[qi][idx_B  ];
                                            u_T   *= u_qi[qi][idx_T  ];
                                            u_TT  *= u_qi[qi][idx_TT ];
                                            u_TTT *= u_qi[qi][idx_TTT];
                                        }
                                    }
                                    
                                    const double dudy = (double(1)/double(60)*(u_TTT - u_BBB)
                                        - double(3)/double(20)*(u_TT - u_BB)
                                        + double(3)/double(4)*(u_T - u_B))/dx[1];
                                    
                                    value_to_add = dudy*weight/((double) n_overlapped);
                                }
                                else if (derivative_direction == 2)
                                {
                                    double u_BBB = double(1);
                                    double u_BB  = double(1);
                                    double u_B   = double(1);
                                    double u_F   = double(1);
                                    double u_FF  = double(1);
                                    double u_FFF = double(1);
                                    
                                    for (int qi = 0; qi < num_quantities; qi++)
                                    {
                                        // Compute linear indices.
                                        const int idx_BBB = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + (k - 3) + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_BB  = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + (k - 2) + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_B   = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + (k - 1) + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_F   = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + (k + 1) + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_FF  = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + (k + 2) + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_FFF = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + (k + 3) + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        if (use_reciprocal[qi])
                                        {
                                            u_BBB /= u_qi[qi][idx_BBB];
                                            u_BB  /= u_qi[qi][idx_BB ];
                                            u_B   /= u_qi[qi][idx_B  ];
                                            u_F   /= u_qi[qi][idx_F  ];
                                            u_FF  /= u_qi[qi][idx_FF ];
                                            u_FFF /= u_qi[qi][idx_FFF];
                                        }
                                        else
                                        {
                                            u_BBB *= u_qi[qi][idx_BBB];
                                            u_BB  *= u_qi[qi][idx_BB ];
                                            u_B   *= u_qi[qi][idx_B  ];
                                            u_F   *= u_qi[qi][idx_F  ];
                                            u_FF  *= u_qi[qi][idx_FF ];
                                            u_FFF *= u_qi[qi][idx_FFF];
                                        }
                                    }
                                    
                                    const double dudz = (double(1)/double(60)*(u_FFF - u_BBB)
                                        - double(3)/double(20)*(u_FF - u_BB)
                                        + double(3)/double(4)*(u_F - u_B))/dx[2];
                                    
                                    value_to_add = dudz*weight/((double) n_overlapped);
                                }
                                else
                                {
                                    TBOX_ERROR(d_object_name
                                        << ": "
                                        << "Cannot take derivative for three-dimensional problem!\n"
                                        << "derivative_direction = " << derivative_direction << " given!\n"
                                        << std::endl);
                                }
                                
                                /*
                                 * Add the data.
                                 */
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    u_der_avg_local[idx_fine] += value_to_add;
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average of derivative.
         */
        
        mpi.Allreduce(
            u_der_avg_local,
            u_der_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_der_avg_local);
    }
    
    return averaged_derivative_quantity;
}


/*
 * Compute averaged value (on product of variables) with only x direction as inhomogeneous direction.
 */
std::vector<double>
FlowModelStatisticsUtilitiesFourEqnConservative::getAveragedQuantityWithInhomogeneousXDirection(
    const std::vector<std::string>& quantity_names,
    const std::vector<int>& component_indices,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context) const
{
    int num_quantities = static_cast<int>(quantity_names.size());
    
    TBOX_ASSERT(static_cast<int>(component_indices.size()) == num_quantities);
    
    std::vector<double> averaged_quantity;
    
    boost::shared_ptr<FlowModel> d_flow_model_tmp = d_flow_model.lock();
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Get the refinement ratio from the finest level to the coarest level.
     */
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    hier::IntVector ratioFinestLevelToCoarestLevel =
        patch_hierarchy->getRatioToCoarserLevel(num_levels - 1);
    for (int li = num_levels - 2; li > 0 ; li--)
    {
        ratioFinestLevelToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(li);
    }
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    boost::shared_ptr<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *patch_hierarchy,
            0,
            num_levels - 1));
    
    /*
     * Get the number of cells of physical domain refined to the finest level.
     */
    
    const hier::BoxContainer& physical_domain = d_grid_geometry->getPhysicalDomain();
    const hier::Box& physical_domain_box = physical_domain.front();
    const hier::IntVector& physical_domain_dims = physical_domain_box.numberCells();
    const hier::IntVector finest_level_dims = physical_domain_dims*ratioFinestLevelToCoarestLevel;
    
    /*
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        double* avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_quantity.resize(finest_level_dim_0);
        double* avg_global = averaged_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            avg_local[i] = double(0);
            avg_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                /*
                 * Register the patch and data in the flow model and compute the corresponding
                 * correlation.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
                }
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<boost::shared_ptr<pdat::CellData<double> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = d_flow_model_tmp->getGlobalCellData(quantity_names[qi]);
                }
                
                std::vector<double*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                }
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        /*
                         * Compute the index of the data point and count how many times the data is repeated.
                         */
                        
                        const hier::Index idx_pt(tbox::Dimension(1), idx_lo_0 + i);
                        
                        int n_overlapped = 1;
                        
                        for (hier::BoxContainer::BoxContainerConstIterator iob(
                                patch_overlapped_visible_boxes.begin());
                             iob != patch_overlapped_visible_boxes.end();
                             iob++)
                        {
                            const hier::Box& patch_overlapped_visible_box = *iob;
                            
                            if (patch_overlapped_visible_box.contains(idx_pt))
                            {
                                n_overlapped++;
                            }
                        }
                        
                        /*
                         * Compute the linear indices and the data to add.
                         */
                        
                        for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                        {
                            const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                            
                            const int idx_q0 = relative_idx_lo_0 + i + num_ghosts_0_u_qi[0];
                            
                            double avg = u_qi[0][idx_q0];
                            
                            for (int qi = 1; qi < num_quantities; qi++)
                            {
                                const int idx_qi = relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi];
                                
                                avg *= u_qi[qi][idx_qi];
                            }
                            
                            avg_local[idx_fine] += (avg/((double) n_overlapped));
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        mpi.Allreduce(
            avg_local,
            avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(avg_local);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        
        double* avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_quantity.resize(finest_level_dim_0);
        double* avg_global = averaged_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            avg_local[i] = double(0);
            avg_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and data in the flow model and compute the corresponding
                 * correlation.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
                }
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<boost::shared_ptr<pdat::CellData<double> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = d_flow_model_tmp->getGlobalCellData(quantity_names[qi]);
                }
                
                std::vector<double*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                }
                
                const double weight = dx[1]/L_y;
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            /*
                             * Compute the index of the data point and count how many times the data is repeated.
                             */
                            
                            const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j);
                            
                            int n_overlapped = 1;
                            
                            for (hier::BoxContainer::BoxContainerConstIterator iob(
                                    patch_overlapped_visible_boxes.begin());
                                 iob != patch_overlapped_visible_boxes.end();
                                 iob++)
                            {
                                const hier::Box& patch_overlapped_visible_box = *iob;
                                
                                if (patch_overlapped_visible_box.contains(idx_pt))
                                {
                                    n_overlapped++;
                                }
                            }
                            
                            /*
                             * Compute the linear indices and the data to add.
                             */
                            
                            for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                            {
                                const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                
                                const int idx_q0 = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[0]) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0];
                                
                                double avg = u_qi[0][idx_q0];
                                
                                for (int qi = 1; qi < num_quantities; qi++)
                                {
                                    const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    avg *= u_qi[qi][idx_qi];
                                }
                                
                                avg_local[idx_fine] += (avg*weight/((double) n_overlapped));
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        mpi.Allreduce(
            avg_local,
            avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(avg_local);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        const double L_z = x_hi[2] - x_lo[2];
        
        double* avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_quantity.resize(finest_level_dim_0);
        double* avg_global = averaged_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            avg_local[i] = double(0);
            avg_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and data in the flow model and compute the corresponding
                 * correlation.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
                }
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<boost::shared_ptr<pdat::CellData<double> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = d_flow_model_tmp->getGlobalCellData(quantity_names[qi]);
                }
                
                std::vector<double*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> num_ghosts_2_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                std::vector<int> ghostcell_dim_1_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                num_ghosts_2_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                ghostcell_dim_1_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    num_ghosts_2_u_qi.push_back(num_ghosts_u_qi[2]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                    ghostcell_dim_1_u_qi.push_back(ghostcell_dims_u_qi[1]);
                }
                
                const double weight = dx[1]*dx[2]/(L_y*L_z);
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    const int interior_dim_2 = interior_dims[2];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int idx_lo_2 = index_lo[2];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    const int relative_idx_lo_2 = relative_index_lo[2];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                /*
                                 * Compute the index of the data point and count how many times the data is repeated.
                                 */
                                
                                const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j, idx_lo_2 + k);
                                
                                int n_overlapped = 1;
                                
                                for (hier::BoxContainer::BoxContainerConstIterator iob(
                                        patch_overlapped_visible_boxes.begin());
                                     iob != patch_overlapped_visible_boxes.end();
                                     iob++)
                                {
                                    const hier::Box& patch_overlapped_visible_box = *iob;
                                    
                                    if (patch_overlapped_visible_box.contains(idx_pt))
                                    {
                                        n_overlapped++;
                                    }
                                }
                                
                                /*
                                 * Compute the linear index and the data to add.
                                 */
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    const int idx_q0 = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[0]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0] +
                                        (relative_idx_lo_2 + k + num_ghosts_2_u_qi[0])*ghostcell_dim_0_u_qi[0]*
                                            ghostcell_dim_1_u_qi[0];
                                    
                                    double avg = u_qi[0][idx_q0];
                                    
                                    for (int qi = 1; qi < num_quantities; qi++)
                                    {
                                        const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        avg *= u_qi[qi][idx_qi];
                                    }
                                    
                                    avg_local[idx_fine] += (avg*weight/((double) n_overlapped));
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        mpi.Allreduce(
            avg_local,
            avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(avg_local);
    }
    
    return averaged_quantity;
}


/*
 * Compute averaged value (on product of variables) with only x direction as inhomogeneous direction.
 */
std::vector<double>
FlowModelStatisticsUtilitiesFourEqnConservative::getAveragedQuantityWithInhomogeneousXDirection(
    const std::vector<std::string>& quantity_names,
    const std::vector<int>& component_indices,
    const std::vector<bool>& use_reciprocal,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context) const
{
    int num_quantities = static_cast<int>(quantity_names.size());
    
    TBOX_ASSERT(static_cast<int>(component_indices.size()) == num_quantities);
    
    std::vector<double> averaged_quantity;
    
    boost::shared_ptr<FlowModel> d_flow_model_tmp = d_flow_model.lock();
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Get the refinement ratio from the finest level to the coarest level.
     */
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    hier::IntVector ratioFinestLevelToCoarestLevel =
        patch_hierarchy->getRatioToCoarserLevel(num_levels - 1);
    for (int li = num_levels - 2; li > 0 ; li--)
    {
        ratioFinestLevelToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(li);
    }
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    boost::shared_ptr<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *patch_hierarchy,
            0,
            num_levels - 1));
    
    /*
     * Get the number of cells of physical domain refined to the finest level.
     */
    
    const hier::BoxContainer& physical_domain = d_grid_geometry->getPhysicalDomain();
    const hier::Box& physical_domain_box = physical_domain.front();
    const hier::IntVector& physical_domain_dims = physical_domain_box.numberCells();
    const hier::IntVector finest_level_dims = physical_domain_dims*ratioFinestLevelToCoarestLevel;
    
    /*
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        double* avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_quantity.resize(finest_level_dim_0);
        double* avg_global = averaged_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            avg_local[i] = double(0);
            avg_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                /*
                 * Register the patch and data in the flow model and compute the corresponding
                 * correlation.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
                }
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<boost::shared_ptr<pdat::CellData<double> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = d_flow_model_tmp->getGlobalCellData(quantity_names[qi]);
                }
                
                std::vector<double*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                }
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        /*
                         * Compute the index of the data point and count how many times the data is repeated.
                         */
                        
                        const hier::Index idx_pt(tbox::Dimension(1), idx_lo_0 + i);
                        
                        int n_overlapped = 1;
                        
                        for (hier::BoxContainer::BoxContainerConstIterator iob(
                                patch_overlapped_visible_boxes.begin());
                             iob != patch_overlapped_visible_boxes.end();
                             iob++)
                        {
                            const hier::Box& patch_overlapped_visible_box = *iob;
                            
                            if (patch_overlapped_visible_box.contains(idx_pt))
                            {
                                n_overlapped++;
                            }
                        }
                        
                        /*
                         * Compute the linear indices and the data to add.
                         */
                        
                        for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                        {
                            const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                            
                            double avg = double(1);
                            
                            for (int qi = 0; qi < num_quantities; qi++)
                            {
                                const int idx_qi = relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi];
                                
                                if (use_reciprocal[qi])
                                {
                                    avg /= u_qi[qi][idx_qi];
                                }
                                else
                                {
                                    avg *= u_qi[qi][idx_qi];
                                }
                            }
                            
                            avg_local[idx_fine] += (avg/((double) n_overlapped));
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        mpi.Allreduce(
            avg_local,
            avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(avg_local);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        
        double* avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_quantity.resize(finest_level_dim_0);
        double* avg_global = averaged_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            avg_local[i] = double(0);
            avg_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and data in the flow model and compute the corresponding
                 * correlation.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
                }
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<boost::shared_ptr<pdat::CellData<double> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = d_flow_model_tmp->getGlobalCellData(quantity_names[qi]);
                }
                
                std::vector<double*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                }
                
                const double weight = dx[1]/L_y;
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            /*
                             * Compute the index of the data point and count how many times the data is repeated.
                             */
                            
                            const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j);
                            
                            int n_overlapped = 1;
                            
                            for (hier::BoxContainer::BoxContainerConstIterator iob(
                                    patch_overlapped_visible_boxes.begin());
                                 iob != patch_overlapped_visible_boxes.end();
                                 iob++)
                            {
                                const hier::Box& patch_overlapped_visible_box = *iob;
                                
                                if (patch_overlapped_visible_box.contains(idx_pt))
                                {
                                    n_overlapped++;
                                }
                            }
                            
                            /*
                             * Compute the linear indices and the data to add.
                             */
                            
                            for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                            {
                                const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                
                                double avg = double(1);
                                
                                for (int qi = 0; qi < num_quantities; qi++)
                                {
                                    const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    if (use_reciprocal[qi])
                                    {
                                        avg /= u_qi[qi][idx_qi];
                                    }
                                    else
                                    {
                                        avg *= u_qi[qi][idx_qi];
                                    }
                                }
                                
                                avg_local[idx_fine] += (avg*weight/((double) n_overlapped));
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        mpi.Allreduce(
            avg_local,
            avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(avg_local);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        const double L_z = x_hi[2] - x_lo[2];
        
        double* avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_quantity.resize(finest_level_dim_0);
        double* avg_global = averaged_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            avg_local[i] = double(0);
            avg_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and data in the flow model and compute the corresponding
                 * correlation.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
                }
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<boost::shared_ptr<pdat::CellData<double> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = d_flow_model_tmp->getGlobalCellData(quantity_names[qi]);
                }
                
                std::vector<double*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> num_ghosts_2_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                std::vector<int> ghostcell_dim_1_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                num_ghosts_2_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                ghostcell_dim_1_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    num_ghosts_2_u_qi.push_back(num_ghosts_u_qi[2]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                    ghostcell_dim_1_u_qi.push_back(ghostcell_dims_u_qi[1]);
                }
                
                const double weight = dx[1]*dx[2]/(L_y*L_z);
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    const int interior_dim_2 = interior_dims[2];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int idx_lo_2 = index_lo[2];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    const int relative_idx_lo_2 = relative_index_lo[2];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                /*
                                 * Compute the index of the data point and count how many times the data is repeated.
                                 */
                                
                                const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j, idx_lo_2 + k);
                                
                                int n_overlapped = 1;
                                
                                for (hier::BoxContainer::BoxContainerConstIterator iob(
                                        patch_overlapped_visible_boxes.begin());
                                     iob != patch_overlapped_visible_boxes.end();
                                     iob++)
                                {
                                    const hier::Box& patch_overlapped_visible_box = *iob;
                                    
                                    if (patch_overlapped_visible_box.contains(idx_pt))
                                    {
                                        n_overlapped++;
                                    }
                                }
                                
                                /*
                                 * Compute the linear index and the data to add.
                                 */
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    double avg = double(1);
                                    
                                    for (int qi = 0; qi < num_quantities; qi++)
                                    {
                                        const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        if (use_reciprocal[qi])
                                        {
                                            avg /= u_qi[qi][idx_qi];
                                        }
                                        else
                                        {
                                            avg *= u_qi[qi][idx_qi];
                                        }
                                    }
                                    
                                    avg_local[idx_fine] += (avg*weight/((double) n_overlapped));
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        mpi.Allreduce(
            avg_local,
            avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(avg_local);
    }
    
    return averaged_quantity;
}


/*
 * Compute correlation with only x direction as inhomogeneous direction.
 */
std::vector<double>
FlowModelStatisticsUtilitiesFourEqnConservative::getQuantityCorrelationWithInhomogeneousXDirection(
    const std::vector<std::string>& quantity_names,
    const std::vector<int>& component_indices,
    const std::vector<std::vector<double> >& averaged_quantities,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context) const
{
    int num_quantities = static_cast<int>(quantity_names.size());
    
    TBOX_ASSERT(static_cast<int>(component_indices.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(averaged_quantities.size()) == num_quantities);
    
    std::vector<double> correlation;
    
    boost::shared_ptr<FlowModel> d_flow_model_tmp = d_flow_model.lock();
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Get the refinement ratio from the finest level to the coarest level.
     */
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    hier::IntVector ratioFinestLevelToCoarestLevel =
        patch_hierarchy->getRatioToCoarserLevel(num_levels - 1);
    for (int li = num_levels - 2; li > 0 ; li--)
    {
        ratioFinestLevelToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(li);
    }
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    boost::shared_ptr<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *patch_hierarchy,
            0,
            num_levels - 1));
    
    /*
     * Get the number of cells of physical domain refined to the finest level.
     */
    
    const hier::BoxContainer& physical_domain = d_grid_geometry->getPhysicalDomain();
    const hier::Box& physical_domain_box = physical_domain.front();
    const hier::IntVector& physical_domain_dims = physical_domain_box.numberCells();
    const hier::IntVector finest_level_dims = physical_domain_dims*ratioFinestLevelToCoarestLevel;
    
    /*
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << "FlowModelStatisticsUtilitiesFourEqnConservative::\n"
            << "getQuantityCorrelationWithInhomogeneousXDirection:\n"
            << "Not implemented for 1D problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        
        std::vector<const double*> u_qi_avg_global;
        u_qi_avg_global.reserve(num_quantities);
        for (int qi = 0; qi < num_quantities; qi++)
        {
            u_qi_avg_global.push_back(averaged_quantities[qi].data());
        }
        
        double* corr_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        correlation.resize(finest_level_dim_0);
        double* corr_global = correlation.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            corr_local[i] = double(0);
            corr_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                // Get the dimensions of box that covers the interior of patch.
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and data in the flow model and compute the corresponding
                 * correlation.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
                }
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<boost::shared_ptr<pdat::CellData<double> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = d_flow_model_tmp->getGlobalCellData(quantity_names[qi]);
                }
                
                std::vector<double*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                }
                
                const double weight = dx[1]/L_y;
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            /*
                             * Compute the index of the data point and count how many times the data is repeated.
                             */
                            
                            const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j);
                            
                            int n_overlapped = 1;
                            
                            for (hier::BoxContainer::BoxContainerConstIterator iob(
                                    patch_overlapped_visible_boxes.begin());
                                 iob != patch_overlapped_visible_boxes.end();
                                 iob++)
                            {
                                const hier::Box& patch_overlapped_visible_box = *iob;
                                
                                if (patch_overlapped_visible_box.contains(idx_pt))
                                {
                                    n_overlapped++;
                                }
                            }
                            
                            /*
                             * Compute the linear indices and the data to add.
                             */
                            
                            for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                            {
                                const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                
                                const int idx_q0 = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[0]) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0];
                                
                                double corr = u_qi[0][idx_q0] - u_qi_avg_global[0][idx_fine];
                                
                                for (int qi = 1; qi < num_quantities; qi++)
                                {
                                    const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    corr *= (u_qi[qi][idx_qi] - u_qi_avg_global[qi][idx_fine]);
                                }
                                
                                corr_local[idx_fine] += (corr*weight/((double) n_overlapped));
                            }
                            
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        mpi.Allreduce(
            corr_local,
            corr_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(corr_local);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        const double L_z = x_hi[2] - x_lo[2];
        
        std::vector<const double*> u_qi_avg_global;
        u_qi_avg_global.reserve(num_quantities);
        for (int qi = 0; qi < num_quantities; qi++)
        {
            u_qi_avg_global.push_back(averaged_quantities[qi].data());
        }
        
        double* corr_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        correlation.resize(finest_level_dim_0);
        double* corr_global = correlation.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            corr_local[i] = double(0);
            corr_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                // Get the dimensions of box that covers the interior of patch.
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and data in the flow model and compute the corresponding
                 * correlation.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
                }
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<boost::shared_ptr<pdat::CellData<double> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = d_flow_model_tmp->getGlobalCellData(quantity_names[qi]);
                }
                
                std::vector<double*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> num_ghosts_2_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                std::vector<int> ghostcell_dim_1_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                num_ghosts_2_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                ghostcell_dim_1_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    num_ghosts_2_u_qi.push_back(num_ghosts_u_qi[2]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                    ghostcell_dim_1_u_qi.push_back(ghostcell_dims_u_qi[1]);
                }
                
                const double weight = dx[1]*dx[2]/(L_y*L_z);
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    const int interior_dim_2 = interior_dims[2];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int idx_lo_2 = index_lo[2];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    const int relative_idx_lo_2 = relative_index_lo[2];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                /*
                                 * Compute the index of the data point and count how many times the data is repeated.
                                 */
                                
                                const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j, idx_lo_2 + k);
                                
                                int n_overlapped = 1;
                                
                                for (hier::BoxContainer::BoxContainerConstIterator iob(
                                        patch_overlapped_visible_boxes.begin());
                                     iob != patch_overlapped_visible_boxes.end();
                                     iob++)
                                {
                                    const hier::Box& patch_overlapped_visible_box = *iob;
                                    
                                    if (patch_overlapped_visible_box.contains(idx_pt))
                                    {
                                        n_overlapped++;
                                    }
                                }
                                
                                /*
                                 * Compute the linear index and the data to add.
                                 */
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    const int idx_q0 = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[0]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0] +
                                        (relative_idx_lo_2 + k + num_ghosts_2_u_qi[0])*ghostcell_dim_0_u_qi[0]*
                                            ghostcell_dim_1_u_qi[0];
                                    
                                    double corr = u_qi[0][idx_q0] - u_qi_avg_global[0][idx_fine];
                                    
                                    for (int qi = 1; qi < num_quantities; qi++)
                                    {
                                        const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        corr *= (u_qi[qi][idx_qi] - u_qi_avg_global[qi][idx_fine]);
                                    }
                                    
                                    corr_local[idx_fine] += (corr*weight/((double) n_overlapped));
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        mpi.Allreduce(
            corr_local,
            corr_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(corr_local);
    }
    
    return correlation;
}


/*
 * Compute correlation with only x direction as inhomogeneous direction.
 */
std::vector<double>
FlowModelStatisticsUtilitiesFourEqnConservative::getQuantityCorrelationWithInhomogeneousXDirection(
    const std::vector<std::string>& quantity_names,
    const std::vector<int>& component_indices,
    const std::vector<bool>& use_reciprocal,
    const std::vector<std::vector<double> >& averaged_quantities,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context) const
{
    int num_quantities = static_cast<int>(quantity_names.size());
    
    TBOX_ASSERT(static_cast<int>(component_indices.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(use_reciprocal.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(averaged_quantities.size()) == num_quantities);
    
    std::vector<double> correlation;
    
    boost::shared_ptr<FlowModel> d_flow_model_tmp = d_flow_model.lock();
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Get the refinement ratio from the finest level to the coarest level.
     */
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    hier::IntVector ratioFinestLevelToCoarestLevel =
        patch_hierarchy->getRatioToCoarserLevel(num_levels - 1);
    for (int li = num_levels - 2; li > 0 ; li--)
    {
        ratioFinestLevelToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(li);
    }
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    boost::shared_ptr<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *patch_hierarchy,
            0,
            num_levels - 1));
    
    /*
     * Get the number of cells of physical domain refined to the finest level.
     */
    
    const hier::BoxContainer& physical_domain = d_grid_geometry->getPhysicalDomain();
    const hier::Box& physical_domain_box = physical_domain.front();
    const hier::IntVector& physical_domain_dims = physical_domain_box.numberCells();
    const hier::IntVector finest_level_dims = physical_domain_dims*ratioFinestLevelToCoarestLevel;
    
    /*
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << "FlowModelStatisticsUtilitiesFourEqnConservative::\n"
            << "getQuantityCorrelationWithInhomogeneousXDirection:\n"
            << "Not implemented for 1D problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        
        std::vector<const double*> u_qi_avg_global;
        u_qi_avg_global.reserve(num_quantities);
        for (int qi = 0; qi < num_quantities; qi++)
        {
            u_qi_avg_global.push_back(averaged_quantities[qi].data());
        }
        
        double* corr_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        correlation.resize(finest_level_dim_0);
        double* corr_global = correlation.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            corr_local[i] = double(0);
            corr_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                // Get the dimensions of box that covers the interior of patch.
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and data in the flow model and compute the corresponding
                 * correlation.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
                }
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<boost::shared_ptr<pdat::CellData<double> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = d_flow_model_tmp->getGlobalCellData(quantity_names[qi]);
                }
                
                std::vector<double*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                }
                
                const double weight = dx[1]/L_y;
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            /*
                             * Compute the index of the data point and count how many times the data is repeated.
                             */
                            
                            const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j);
                            
                            int n_overlapped = 1;
                            
                            for (hier::BoxContainer::BoxContainerConstIterator iob(
                                    patch_overlapped_visible_boxes.begin());
                                 iob != patch_overlapped_visible_boxes.end();
                                 iob++)
                            {
                                const hier::Box& patch_overlapped_visible_box = *iob;
                                
                                if (patch_overlapped_visible_box.contains(idx_pt))
                                {
                                    n_overlapped++;
                                }
                            }
                            
                            /*
                             * Compute the linear indices and the data to add.
                             */
                            
                            for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                            {
                                const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                
                                const int idx_q0 = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[0]) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0];
                                
                                double corr = double(0);
                                if (use_reciprocal[0])
                                {
                                    corr = double(1)/(u_qi[0][idx_q0]) - u_qi_avg_global[0][idx_fine];
                                }
                                else
                                {
                                    corr = u_qi[0][idx_q0] - u_qi_avg_global[0][idx_fine];
                                }
                                
                                for (int qi = 1; qi < num_quantities; qi++)
                                {
                                    const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    if (use_reciprocal[qi])
                                    {
                                        corr *= (double(1)/(u_qi[qi][idx_qi]) - u_qi_avg_global[qi][idx_fine]);
                                    }
                                    else
                                    {
                                        corr *= (u_qi[qi][idx_qi] - u_qi_avg_global[qi][idx_fine]);
                                    }
                                }
                                
                                corr_local[idx_fine] += (corr*weight/((double) n_overlapped));
                            }
                            
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        mpi.Allreduce(
            corr_local,
            corr_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(corr_local);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        const double L_z = x_hi[2] - x_lo[2];
        
        std::vector<const double*> u_qi_avg_global;
        u_qi_avg_global.reserve(num_quantities);
        for (int qi = 0; qi < num_quantities; qi++)
        {
            u_qi_avg_global.push_back(averaged_quantities[qi].data());
        }
        
        double* corr_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        correlation.resize(finest_level_dim_0);
        double* corr_global = correlation.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            corr_local[i] = double(0);
            corr_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                // Get the dimensions of box that covers the interior of patch.
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and data in the flow model and compute the corresponding
                 * correlation.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
                }
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<boost::shared_ptr<pdat::CellData<double> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = d_flow_model_tmp->getGlobalCellData(quantity_names[qi]);
                }
                
                std::vector<double*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> num_ghosts_2_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                std::vector<int> ghostcell_dim_1_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                num_ghosts_2_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                ghostcell_dim_1_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    num_ghosts_2_u_qi.push_back(num_ghosts_u_qi[2]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                    ghostcell_dim_1_u_qi.push_back(ghostcell_dims_u_qi[1]);
                }
                
                const double weight = dx[1]*dx[2]/(L_y*L_z);
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    const int interior_dim_2 = interior_dims[2];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int idx_lo_2 = index_lo[2];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    const int relative_idx_lo_2 = relative_index_lo[2];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                /*
                                 * Compute the index of the data point and count how many times the data is repeated.
                                 */
                                
                                const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j, idx_lo_2 + k);
                                
                                int n_overlapped = 1;
                                
                                for (hier::BoxContainer::BoxContainerConstIterator iob(
                                        patch_overlapped_visible_boxes.begin());
                                     iob != patch_overlapped_visible_boxes.end();
                                     iob++)
                                {
                                    const hier::Box& patch_overlapped_visible_box = *iob;
                                    
                                    if (patch_overlapped_visible_box.contains(idx_pt))
                                    {
                                        n_overlapped++;
                                    }
                                }
                                
                                /*
                                 * Compute the linear index and the data to add.
                                 */
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    const int idx_q0 = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[0]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0] +
                                        (relative_idx_lo_2 + k + num_ghosts_2_u_qi[0])*ghostcell_dim_0_u_qi[0]*
                                            ghostcell_dim_1_u_qi[0];
                                    
                                    double corr = double(0);
                                    if (use_reciprocal[0])
                                    {
                                        corr = double(1)/(u_qi[0][idx_q0]) - u_qi_avg_global[0][idx_fine];
                                    }
                                    else
                                    {
                                        corr = u_qi[0][idx_q0] - u_qi_avg_global[0][idx_fine];
                                    }
                                    
                                    for (int qi = 1; qi < num_quantities; qi++)
                                    {
                                        const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        if (use_reciprocal[qi])
                                        {
                                            corr *= (double(1)/(u_qi[qi][idx_qi]) - u_qi_avg_global[qi][idx_fine]);
                                        }
                                        else
                                        {
                                            corr *= (u_qi[qi][idx_qi] - u_qi_avg_global[qi][idx_fine]);
                                        }
                                    }
                                    
                                    corr_local[idx_fine] += (corr*weight/((double) n_overlapped));
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        mpi.Allreduce(
            corr_local,
            corr_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(corr_local);
    }
    
    return correlation;
}


/*
 * Compute correlation with only x direction as inhomogeneous direction.
 */
std::vector<double>
FlowModelStatisticsUtilitiesFourEqnConservative::getQuantityCorrelationWithInhomogeneousXDirection(
    const std::vector<std::string>& quantity_names,
    const std::vector<int>& component_indices,
    const std::vector<bool>& use_reciprocal,
    const std::vector<int>& derivative_directions,
    const std::vector<std::vector<double> >& averaged_quantities,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context) const
{
    int num_quantities = static_cast<int>(quantity_names.size());
    
    TBOX_ASSERT(static_cast<int>(component_indices.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(use_reciprocal.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(derivative_directions.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(averaged_quantities.size()) == num_quantities);
    
    std::vector<double> correlation;
    
    boost::shared_ptr<FlowModel> d_flow_model_tmp = d_flow_model.lock();
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Get the refinement ratio from the finest level to the coarest level.
     */
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    hier::IntVector ratioFinestLevelToCoarestLevel =
        patch_hierarchy->getRatioToCoarserLevel(num_levels - 1);
    for (int li = num_levels - 2; li > 0 ; li--)
    {
        ratioFinestLevelToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(li);
    }
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    boost::shared_ptr<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *patch_hierarchy,
            0,
            num_levels - 1));
    
    /*
     * Get the number of cells of physical domain refined to the finest level.
     */
    
    const hier::BoxContainer& physical_domain = d_grid_geometry->getPhysicalDomain();
    const hier::Box& physical_domain_box = physical_domain.front();
    const hier::IntVector& physical_domain_dims = physical_domain_box.numberCells();
    const hier::IntVector finest_level_dims = physical_domain_dims*ratioFinestLevelToCoarestLevel;
    
    /*
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << "FlowModelStatisticsUtilitiesFourEqnConservative::\n"
            << "getQuantityCorrelationWithInhomogeneousXDirection:\n"
            << "Not implemented for 1D problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        
        std::vector<const double*> u_qi_avg_global;
        u_qi_avg_global.reserve(num_quantities);
        for (int qi = 0; qi < num_quantities; qi++)
        {
            u_qi_avg_global.push_back(averaged_quantities[qi].data());
        }
        
        double* corr_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        correlation.resize(finest_level_dim_0);
        double* corr_global = correlation.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            corr_local[i] = double(0);
            corr_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                // Get the dimensions of box that covers the interior of patch.
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and data in the flow model and compute the corresponding
                 * correlation.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                // HARD CODE TO BE SIXTH ORDER CENTRAL SCHEME FOR DIFFERENTIATION.
                TBOX_ASSERT(num_ghosts >= hier::IntVector::getOne(d_dim)*3);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    std::unordered_map<std::string, hier::IntVector>::const_iterator quantity_it =
                        num_subghosts_of_data.find(quantity_names[qi]);
                    
                    if (quantity_it == num_subghosts_of_data.end())
                    {
                        if (derivative_directions[qi] >= 0)
                        {
                            num_subghosts_of_data.insert(
                                std::pair<std::string, hier::IntVector>(quantity_names[qi], num_ghosts));
                        }
                        else
                        {
                            num_subghosts_of_data.insert(
                                std::pair<std::string, hier::IntVector>(
                                    quantity_names[qi], hier::IntVector::getZero(d_dim)));
                        }
                    }
                    else
                    {
                        if (derivative_directions[qi] >= 0 && quantity_it->second < num_ghosts)
                        {
                            num_subghosts_of_data.erase(quantity_it);
                            
                            num_subghosts_of_data.insert(
                                std::pair<std::string, hier::IntVector>(quantity_names[qi], num_ghosts));
                        }
                    }
                }
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<boost::shared_ptr<pdat::CellData<double> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = d_flow_model_tmp->getGlobalCellData(quantity_names[qi]);
                }
                
                std::vector<double*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                }
                
                const double weight = dx[1]/L_y;
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            /*
                             * Compute the index of the data point and count how many times the data is repeated.
                             */
                            
                            const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j);
                            
                            int n_overlapped = 1;
                            
                            for (hier::BoxContainer::BoxContainerConstIterator iob(
                                    patch_overlapped_visible_boxes.begin());
                                 iob != patch_overlapped_visible_boxes.end();
                                 iob++)
                            {
                                const hier::Box& patch_overlapped_visible_box = *iob;
                                
                                if (patch_overlapped_visible_box.contains(idx_pt))
                                {
                                    n_overlapped++;
                                }
                            }
                            
                            /*
                             * Perform operations on quantities.
                             */
                            
                            double u_values[num_quantities];
                            
                            for (int qi = 0; qi < num_quantities; qi++)
                            {
                                if (derivative_directions[qi] == -1)
                                {
                                    const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    if (use_reciprocal[qi])
                                    {
                                        u_values[qi] = double(1)/(u_qi[qi][idx_qi]);
                                    }
                                    else
                                    {
                                        u_values[qi] = u_qi[qi][idx_qi];
                                    }
                                }
                                else if (derivative_directions[qi] == 0)
                                {
                                    const int idx_LLL = (relative_idx_lo_0 + (i - 3) + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    const int idx_LL  = (relative_idx_lo_0 + (i - 2) + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    const int idx_L   = (relative_idx_lo_0 + (i - 1) + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    const int idx_R   = (relative_idx_lo_0 + (i + 1) + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    const int idx_RR  = (relative_idx_lo_0 + (i + 2) + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    const int idx_RRR = (relative_idx_lo_0 + (i + 3) + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    if (use_reciprocal[qi])
                                    {
                                        u_values[qi] = (double(1)/double(60)*(double(1)/u_qi[qi][idx_RRR] - double(1)/u_qi[qi][idx_LLL])
                                            - double(3)/double(20)*(double(1)/u_qi[qi][idx_RR] - double(1)/u_qi[qi][idx_LL])
                                            + double(3)/double(4)*(double(1)/u_qi[qi][idx_R] - double(1)/u_qi[qi][idx_L]))/dx[0];
                                    }
                                    else
                                    {
                                        u_values[qi] = (double(1)/double(60)*(u_qi[qi][idx_RRR] - u_qi[qi][idx_LLL])
                                            - double(3)/double(20)*(u_qi[qi][idx_RR] - u_qi[qi][idx_LL])
                                            + double(3)/double(4)*(u_qi[qi][idx_R] - u_qi[qi][idx_L]))/dx[0];
                                    }
                                }
                                else if (derivative_directions[qi] == 1)
                                {
                                    const int idx_BBB = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + (j - 3) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    const int idx_BB  = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + (j - 2) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    const int idx_B   = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + (j - 1) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    const int idx_T   = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + (j + 1) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    const int idx_TT  = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + (j + 2) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    const int idx_TTT = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + (j + 3) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    if (use_reciprocal[qi])
                                    {
                                        u_values[qi] = (double(1)/double(60)*(double(1)/u_qi[qi][idx_TTT] - double(1)/u_qi[qi][idx_BBB])
                                            - double(3)/double(20)*(double(1)/u_qi[qi][idx_TT] - double(1)/u_qi[qi][idx_BB])
                                            + double(3)/double(4)*(double(1)/u_qi[qi][idx_T] - double(1)/u_qi[qi][idx_B]))/dx[1];
                                    }
                                    else
                                    {
                                        u_values[qi] = (double(1)/double(60)*(u_qi[qi][idx_TTT] - u_qi[qi][idx_BBB])
                                            - double(3)/double(20)*(u_qi[qi][idx_TT] - u_qi[qi][idx_BB])
                                            + double(3)/double(4)*(u_qi[qi][idx_T] - u_qi[qi][idx_B]))/dx[1];
                                    }
                                }
                                else
                                {
                                    TBOX_ERROR(d_object_name
                                        << ": "
                                        << "Cannot take derivative for two-dimensional problem!\n"
                                        << "derivative_direction = " << derivative_directions[qi] << " given!\n"
                                        << std::endl);
                                }
                            }
                            
                            /*
                             * Compute the correlation.
                             */
                            
                            for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                            {
                                const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                
                                double corr = double(1);
                                for (int qi = 0; qi < num_quantities; qi++)
                                {
                                    corr *= (u_values[qi] - u_qi_avg_global[qi][idx_fine]);
                                }
                                
                                corr_local[idx_fine] += (corr*weight/((double) n_overlapped));
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        mpi.Allreduce(
            corr_local,
            corr_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(corr_local);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        const double L_z = x_hi[2] - x_lo[2];
        
        std::vector<const double*> u_qi_avg_global;
        u_qi_avg_global.reserve(num_quantities);
        for (int qi = 0; qi < num_quantities; qi++)
        {
            u_qi_avg_global.push_back(averaged_quantities[qi].data());
        }
        
        double* corr_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        correlation.resize(finest_level_dim_0);
        double* corr_global = correlation.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            corr_local[i] = double(0);
            corr_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                // Get the dimensions of box that covers the interior of patch.
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and data in the flow model and compute the corresponding
                 * correlation.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                // HARD CODE TO BE SIXTH ORDER CENTRAL SCHEME FOR DIFFERENTIATION.
                TBOX_ASSERT(num_ghosts >= hier::IntVector::getOne(d_dim)*3);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    std::unordered_map<std::string, hier::IntVector>::const_iterator quantity_it =
                        num_subghosts_of_data.find(quantity_names[qi]);
                    
                    if (quantity_it == num_subghosts_of_data.end())
                    {
                        if (derivative_directions[qi] >= 0)
                        {
                            num_subghosts_of_data.insert(
                                std::pair<std::string, hier::IntVector>(quantity_names[qi], num_ghosts));
                        }
                        else
                        {
                            num_subghosts_of_data.insert(
                                std::pair<std::string, hier::IntVector>(
                                    quantity_names[qi], hier::IntVector::getZero(d_dim)));
                        }
                    }
                    else
                    {
                        if (derivative_directions[qi] >= 0 && quantity_it->second < num_ghosts)
                        {
                            num_subghosts_of_data.erase(quantity_it);
                            
                            num_subghosts_of_data.insert(
                                std::pair<std::string, hier::IntVector>(quantity_names[qi], num_ghosts));
                        }
                    }
                }
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<boost::shared_ptr<pdat::CellData<double> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = d_flow_model_tmp->getGlobalCellData(quantity_names[qi]);
                }
                
                std::vector<double*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> num_ghosts_2_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                std::vector<int> ghostcell_dim_1_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                num_ghosts_2_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                ghostcell_dim_1_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    num_ghosts_2_u_qi.push_back(num_ghosts_u_qi[2]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                    ghostcell_dim_1_u_qi.push_back(ghostcell_dims_u_qi[1]);
                }
                
                const double weight = dx[1]*dx[2]/(L_y*L_z);
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    const int interior_dim_2 = interior_dims[2];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int idx_lo_2 = index_lo[2];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    const int relative_idx_lo_2 = relative_index_lo[2];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                /*
                                 * Compute the index of the data point and count how many times the data is repeated.
                                 */
                                
                                const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j, idx_lo_2 + k);
                                
                                int n_overlapped = 1;
                                
                                for (hier::BoxContainer::BoxContainerConstIterator iob(
                                        patch_overlapped_visible_boxes.begin());
                                     iob != patch_overlapped_visible_boxes.end();
                                     iob++)
                                {
                                    const hier::Box& patch_overlapped_visible_box = *iob;
                                    
                                    if (patch_overlapped_visible_box.contains(idx_pt))
                                    {
                                        n_overlapped++;
                                    }
                                }
                                
                                /*
                                 * Perform operations on quantities.
                                 */
                                
                                double u_values[num_quantities];
                                
                                for (int qi = 0; qi < num_quantities; qi++)
                                {
                                    if (derivative_directions[qi] == -1)
                                    {
                                        const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        if (use_reciprocal[qi])
                                        {
                                            u_values[qi] = double(1)/(u_qi[qi][idx_qi]);
                                        }
                                        else
                                        {
                                            u_values[qi] = u_qi[qi][idx_qi];
                                        }
                                    }
                                    else if (derivative_directions[qi] == 0)
                                    {
                                        const int idx_LLL = (relative_idx_lo_0 + (i - 3) + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_LL  = (relative_idx_lo_0 + (i - 2) + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_L   = (relative_idx_lo_0 + (i - 1) + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_R   = (relative_idx_lo_0 + (i + 1) + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_RR  = (relative_idx_lo_0 + (i + 2) + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_RRR = (relative_idx_lo_0 + (i + 3) + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        if (use_reciprocal[qi])
                                        {
                                            u_values[qi] = (double(1)/double(60)*(double(1)/u_qi[qi][idx_RRR] - double(1)/u_qi[qi][idx_LLL])
                                                - double(3)/double(20)*(double(1)/u_qi[qi][idx_RR] - double(1)/u_qi[qi][idx_LL])
                                                + double(3)/double(4)*(double(1)/u_qi[qi][idx_R] - double(1)/u_qi[qi][idx_L]))/dx[0];
                                        }
                                        else
                                        {
                                            u_values[qi] = (double(1)/double(60)*(u_qi[qi][idx_RRR] - u_qi[qi][idx_LLL])
                                                - double(3)/double(20)*(u_qi[qi][idx_RR] - u_qi[qi][idx_LL])
                                                + double(3)/double(4)*(u_qi[qi][idx_R] - u_qi[qi][idx_L]))/dx[0];
                                        }
                                    }
                                    else if (derivative_directions[qi] == 1)
                                    {
                                        const int idx_BBB = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + (j - 3) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_BB  = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + (j - 2) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_B   = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + (j - 1) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_T   = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + (j + 1) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_TT  = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + (j + 2) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_TTT = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + (j + 3) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        if (use_reciprocal[qi])
                                        {
                                            u_values[qi] = (double(1)/double(60)*(double(1)/u_qi[qi][idx_TTT] - double(1)/u_qi[qi][idx_BBB])
                                                - double(3)/double(20)*(double(1)/u_qi[qi][idx_TT] - double(1)/u_qi[qi][idx_BB])
                                                + double(3)/double(4)*(double(1)/u_qi[qi][idx_T] - double(1)/u_qi[qi][idx_B]))/dx[1];
                                        }
                                        else
                                        {
                                            u_values[qi] = (double(1)/double(60)*(u_qi[qi][idx_TTT] - u_qi[qi][idx_BBB])
                                                - double(3)/double(20)*(u_qi[qi][idx_TT] - u_qi[qi][idx_BB])
                                                + double(3)/double(4)*(u_qi[qi][idx_T] - u_qi[qi][idx_B]))/dx[1];
                                        }
                                    }
                                    else if (derivative_directions[qi] == 2)
                                    {
                                        const int idx_BBB = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + (k - 3) + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_BB  = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + (k - 2) + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_B   = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + (k - 1) + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_F   = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + (k + 1) + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_FF  = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + (k + 2) + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_FFF = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + (k + 3) + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        if (use_reciprocal[qi])
                                        {
                                            u_values[qi] = (double(1)/double(60)*(double(1)/u_qi[qi][idx_FFF] - double(1)/u_qi[qi][idx_BBB])
                                                - double(3)/double(20)*(double(1)/u_qi[qi][idx_FF] - double(1)/u_qi[qi][idx_BB])
                                                + double(3)/double(4)*(double(1)/u_qi[qi][idx_F] - double(1)/u_qi[qi][idx_B]))/dx[2];
                                        }
                                        else
                                        {
                                            u_values[qi] = (double(1)/double(60)*(u_qi[qi][idx_FFF] - u_qi[qi][idx_BBB])
                                                - double(3)/double(20)*(u_qi[qi][idx_FF] - u_qi[qi][idx_BB])
                                                + double(3)/double(4)*(u_qi[qi][idx_F] - u_qi[qi][idx_B]))/dx[2];
                                        }
                                    }
                                    else
                                    {
                                        TBOX_ERROR(d_object_name
                                            << ": "
                                            << "Cannot take derivative for three-dimensional problem!\n"
                                            << "derivative_direction = " << derivative_directions[qi] << " given!\n"
                                            << std::endl);
                                    }
                                }
                                
                                /*
                                 * Compute the correlation.
                                 */
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    double corr = double(1);
                                    for (int qi = 0; qi < num_quantities; qi++)
                                    {
                                        corr *= (u_values[qi] - u_qi_avg_global[qi][idx_fine]);
                                    }
                                    
                                    corr_local[idx_fine] += (corr*weight/((double) n_overlapped));
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        mpi.Allreduce(
            corr_local,
            corr_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(corr_local);
    }
    
    return correlation;
}


/**************************************************************************************************
 * For budgets on filtered variables.
 *************************************************************************************************/


/*
 * Compute pressure with only x direction as inhomogeneous direction.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::computePressure(
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
    boost::shared_ptr<FlowModel> d_flow_model_tmp = d_flow_model.lock();
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    for (int li = 0; li < num_levels; li++)
    {
        /*
         * Get the current patch level.
         */
        
        boost::shared_ptr<hier::PatchLevel> patch_level(
            patch_hierarchy->getPatchLevel(li));
        
        for (hier::PatchLevel::iterator ip(patch_level->begin());
             ip != patch_level->end();
             ip++)
        {
            const boost::shared_ptr<hier::Patch> patch = *ip;
            
            // Get the unfiltered pressure cell data.
            
            boost::shared_ptr<pdat::CellData<double> > data_pressure_unfiltered(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch->getPatchData(s_variable_pressure_unfiltered, data_context)));
            
            // Get the unfiltered pressure cell data.
            
            d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
            
            hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
            
            std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
            
            num_subghosts_of_data.insert(
                std::pair<std::string, hier::IntVector>("PRESSURE", num_ghosts));
            
            d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
            
            d_flow_model_tmp->computeGlobalDerivedCellData();
            
            boost::shared_ptr<pdat::CellData<double> > data_pressure =
                d_flow_model_tmp->getGlobalCellData("PRESSURE");
            
            data_pressure_unfiltered->copy(*data_pressure);
            
            /*
             * Unregister the patch and data of all registered derived cell variables in the flow model.
             */
            
            d_flow_model_tmp->unregisterPatch();
        }
    }
}


/*
 * Compute shear stress with only x direction as inhomogeneous direction.
 */
void 
FlowModelStatisticsUtilitiesFourEqnConservative::computeShearStress(
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
    boost::shared_ptr<FlowModel> d_flow_model_tmp = d_flow_model.lock();
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    if (d_dim == tbox::Dimension(1))
    {
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                // Get the shear stress cell data.
                
                boost::shared_ptr<pdat::CellData<double> > data_shear_stress(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch->getPatchData(s_variable_shear_stress_unfiltered, data_context)));
                
                // Get the pointer to shear stress data.
                double* tau_11 = data_shear_stress->getPointer(0);
                
                // Get the number of ghost cells of the shear stress cell data.
                const hier::IntVector num_ghosts_shear_stress = data_shear_stress->getGhostCellWidth();
                const int num_ghosts_0_shear_stress = num_ghosts_shear_stress[0];
                
                // Get the box that covers the interior of patch.
                const hier::Box interior_box = data_shear_stress->getBox();
                const hier::IntVector interior_dims = interior_box.numberCells();
                
                const int interior_dim_0 = interior_dims[0];
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                // HARD CODE TO BE SIXTH ORDER CENTRAL SCHEME FOR DIFFERENTIATION.
                TBOX_ASSERT(num_ghosts >= hier::IntVector::getOne(d_dim)*3);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", num_ghosts));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTION", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("PRESSURE", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("TEMPERATURE", hier::IntVector::getZero(d_dim)));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_velocity =
                    d_flow_model_tmp->getGlobalCellData("VELOCITY");
                
                boost::shared_ptr<pdat::CellData<double> > data_mass_fraction =
                    d_flow_model_tmp->getGlobalCellData("MASS_FRACTION");
                
                boost::shared_ptr<pdat::CellData<double> > data_pressure =
                    d_flow_model_tmp->getGlobalCellData("PRESSURE");
                
                boost::shared_ptr<pdat::CellData<double> > data_temperature =
                    d_flow_model_tmp->getGlobalCellData("TEMPERATURE");
                
                double* u = data_velocity->getPointer(0);
                std::vector<double*> Y;
                Y.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    Y.push_back(data_mass_fraction->getPointer(si));
                }
                double* p = data_pressure->getPointer(0);
                double* T = data_temperature->getPointer(0);
                
                const hier::IntVector num_ghosts_velocity = data_velocity->getGhostCellWidth();
                const hier::IntVector num_ghosts_mass_fraction = data_mass_fraction->getGhostCellWidth();
                const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
                const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
                
                const int num_ghosts_0_velocity = num_ghosts_velocity[0];
                const int num_ghosts_0_mass_fraction = num_ghosts_mass_fraction[0];
                const int num_ghosts_0_pressure = num_ghosts_pressure[0];
                const int num_ghosts_0_temperature = num_ghosts_temperature[0];
                
                for (int i = 0; i < interior_dim_0; i++)
                {
                    /*
                     * Compute the shear stress.
                     */
                    
                    // Compute the linear index of shear stress.
                    const int idx_shear_stress = i + num_ghosts_0_shear_stress;
                    
                    // Compute the linear indices of mass fraction, pressure and temperature.
                    const int idx_mass_fraction = i + num_ghosts_0_mass_fraction;
                    const int idx_pressure = i + num_ghosts_0_pressure;
                    const int idx_temperature = i + num_ghosts_0_temperature;
                    
                    /*
                     * Compute component 0.
                     */
                    
                    std::vector<const double*> Y_ptr;
                    Y_ptr.resize(d_num_species);
                    for (int si = 0; si < d_num_species; si++)
                    {
                        Y_ptr[si] = &Y[si][idx_mass_fraction];
                    }
                    const double mu = d_equation_of_shear_viscosity_mixing_rules->
                        getShearViscosity(
                            &p[idx_pressure],
                            &T[idx_temperature],
                            Y_ptr);
                    const double mu_v = d_equation_of_bulk_viscosity_mixing_rules->
                        getBulkViscosity(
                            &p[idx_pressure],
                            &T[idx_temperature],
                            Y_ptr);
                    
                    const int idx_vel_x_LLL = (i - 3) + num_ghosts_0_velocity;
                    const int idx_vel_x_LL  = (i - 2) + num_ghosts_0_velocity;
                    const int idx_vel_x_L   = (i - 1) + num_ghosts_0_velocity;
                    const int idx_vel_x_R   = (i + 1) + num_ghosts_0_velocity;
                    const int idx_vel_x_RR  = (i + 2) + num_ghosts_0_velocity;
                    const int idx_vel_x_RRR = (i + 3) + num_ghosts_0_velocity;
                    
                    const double dudx = (double(1)/double(60)*(u[idx_vel_x_RRR] - u[idx_vel_x_LLL])
                        - double(3)/double(20)*(u[idx_vel_x_RR] - u[idx_vel_x_LL])
                        + double(3)/double(4)*(u[idx_vel_x_R] - u[idx_vel_x_L]))/dx[0];
                    
                    tau_11[idx_shear_stress] = (double(4)/double(3)*mu + mu_v)*dudx;
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                // Get the shear stress cell data.
                
                boost::shared_ptr<pdat::CellData<double> > data_shear_stress(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch->getPatchData(s_variable_shear_stress_unfiltered, data_context)));
                
                // Get the pointers to shear stress data.
                double* tau_11 = data_shear_stress->getPointer(0);
                double* tau_12 = data_shear_stress->getPointer(1);
                double* tau_22 = data_shear_stress->getPointer(2);
                
                // Get the number of ghost cells of the shear stress cell data.
                const hier::IntVector num_ghosts_shear_stress = data_shear_stress->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_shear_stress = data_shear_stress->getGhostBox().numberCells();
                
                const int num_ghosts_0_shear_stress = num_ghosts_shear_stress[0];
                const int num_ghosts_1_shear_stress = num_ghosts_shear_stress[1];
                const int ghostcell_dim_0_shear_stress = ghostcell_dims_shear_stress[0];
                
                // Get the box that covers the interior of patch.
                const hier::Box interior_box = data_shear_stress->getBox();
                const hier::IntVector interior_dims = interior_box.numberCells();
                
                const int interior_dim_0 = interior_dims[0];
                const int interior_dim_1 = interior_dims[1];
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                // HARD CODE TO BE SIXTH ORDER CENTRAL SCHEME FOR DIFFERENTIATION.
                TBOX_ASSERT(num_ghosts >= hier::IntVector::getOne(d_dim)*3);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", num_ghosts));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTION", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("PRESSURE", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("TEMPERATURE", hier::IntVector::getZero(d_dim)));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_velocity =
                    d_flow_model_tmp->getGlobalCellData("VELOCITY");
                
                boost::shared_ptr<pdat::CellData<double> > data_mass_fraction =
                    d_flow_model_tmp->getGlobalCellData("MASS_FRACTION");
                
                boost::shared_ptr<pdat::CellData<double> > data_pressure =
                    d_flow_model_tmp->getGlobalCellData("PRESSURE");
                
                boost::shared_ptr<pdat::CellData<double> > data_temperature =
                    d_flow_model_tmp->getGlobalCellData("TEMPERATURE");
                
                double* u = data_velocity->getPointer(0);
                double* v = data_velocity->getPointer(1);
                std::vector<double*> Y;
                Y.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    Y.push_back(data_mass_fraction->getPointer(si));
                }
                double* p = data_pressure->getPointer(0);
                double* T = data_temperature->getPointer(0);
                
                const hier::IntVector num_ghosts_velocity = data_velocity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_velocity = data_velocity->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_mass_fraction = data_mass_fraction->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mass_fraction = data_mass_fraction->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_pressure = data_pressure->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_temperature = data_temperature->getGhostBox().numberCells();
                
                const int num_ghosts_0_velocity = num_ghosts_velocity[0];
                const int num_ghosts_1_velocity = num_ghosts_velocity[1];
                const int ghostcell_dim_0_velocity = ghostcell_dims_velocity[0];
                
                const int num_ghosts_0_mass_fraction = num_ghosts_mass_fraction[0];
                const int num_ghosts_1_mass_fraction = num_ghosts_mass_fraction[1];
                const int ghostcell_dim_0_mass_fraction = ghostcell_dims_mass_fraction[0];
                
                const int num_ghosts_0_pressure = num_ghosts_pressure[0];
                const int num_ghosts_1_pressure = num_ghosts_pressure[1];
                const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
                
                const int num_ghosts_0_temperature = num_ghosts_temperature[0];
                const int num_ghosts_1_temperature = num_ghosts_temperature[1];
                const int ghostcell_dim_0_temperature = ghostcell_dims_temperature[0];
                
                for (int j = 0; j < interior_dim_1; j++)
                {
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        /*
                         * Compute the shear stress.
                         */
                            
                        // Compute the linear indices.
                        const int idx_shear_stress = i + num_ghosts_0_shear_stress +
                            (j + num_ghosts_1_shear_stress)*ghostcell_dim_0_shear_stress;
                        
                        const int idx_mass_fraction = (i + num_ghosts_0_mass_fraction) +
                            (j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                        
                        const int idx_pressure = (i + num_ghosts_0_pressure) +
                            (j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure;
                        
                        const int idx_temperature = (i + num_ghosts_0_temperature) +
                            (j + num_ghosts_1_temperature)*ghostcell_dim_0_temperature;
                        
                        const int idx_vel_x_LLL = ((i - 3) + num_ghosts_0_velocity) +
                            (j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                        
                        const int idx_vel_x_LL  = ((i - 2) + num_ghosts_0_velocity) +
                            (j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                        
                        const int idx_vel_x_L   = ((i - 1) + num_ghosts_0_velocity) +
                            (j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                        
                        const int idx_vel_x_R   = ((i + 1) + num_ghosts_0_velocity) +
                            (j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                        
                        const int idx_vel_x_RR  = ((i + 2) + num_ghosts_0_velocity) +
                            (j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                        
                        const int idx_vel_x_RRR = ((i + 3) + num_ghosts_0_velocity) +
                            (j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                        
                        const int idx_vel_y_BBB = (i + num_ghosts_0_velocity) +
                            ((j - 3) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                        
                        const int idx_vel_y_BB  = (i + num_ghosts_0_velocity) +
                            ((j - 2) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                        
                        const int idx_vel_y_B   = (i + num_ghosts_0_velocity) +
                            ((j - 1) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                        
                        const int idx_vel_y_T   = (i + num_ghosts_0_velocity) +
                            ((j + 1) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                        
                        const int idx_vel_y_TT  = (i + num_ghosts_0_velocity) +
                            ((j + 2) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                        
                        const int idx_vel_y_TTT = (i + num_ghosts_0_velocity) +
                            ((j + 3) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                        
                        std::vector<const double*> Y_ptr;
                        Y_ptr.resize(d_num_species);
                        for (int si = 0; si < d_num_species; si++)
                        {
                            Y_ptr[si] = &Y[si][idx_mass_fraction];
                        }
                        const double mu = d_equation_of_shear_viscosity_mixing_rules->
                            getShearViscosity(
                                &p[idx_pressure],
                                &T[idx_temperature],
                                Y_ptr);
                        const double mu_v = d_equation_of_bulk_viscosity_mixing_rules->
                            getBulkViscosity(
                                &p[idx_pressure],
                                &T[idx_temperature],
                                Y_ptr);
                        
                        const double dudx = (double(1)/double(60)*(u[idx_vel_x_RRR] - u[idx_vel_x_LLL])
                            - double(3)/double(20)*(u[idx_vel_x_RR] - u[idx_vel_x_LL])
                            + double(3)/double(4)*(u[idx_vel_x_R] - u[idx_vel_x_L]))/dx[0];
                        
                        const double dudy = (double(1)/double(60)*(u[idx_vel_y_TTT] - u[idx_vel_y_BBB])
                            - double(3)/double(20)*(u[idx_vel_y_TT] - u[idx_vel_y_BB])
                            + double(3)/double(4)*(u[idx_vel_y_T] - u[idx_vel_y_B]))/dx[1];
                        
                        const double dvdx = (double(1)/double(60)*(v[idx_vel_x_RRR] - v[idx_vel_x_LLL])
                            - double(3)/double(20)*(v[idx_vel_x_RR] - v[idx_vel_x_LL])
                            + double(3)/double(4)*(v[idx_vel_x_R] - v[idx_vel_x_L]))/dx[0];
                        
                        const double dvdy = (double(1)/double(60)*(v[idx_vel_y_TTT] - v[idx_vel_y_BBB])
                            - double(3)/double(20)*(v[idx_vel_y_TT] - v[idx_vel_y_BB])
                            + double(3)/double(4)*(v[idx_vel_y_T] - v[idx_vel_y_B]))/dx[1];
                        
                        tau_11[idx_shear_stress] = (double(4)/double(3)*mu + mu_v)*dudx - (double(2)/double(3)*mu - mu_v)*dvdy;
                        tau_12[idx_shear_stress] = mu*(dudy + dvdx);
                        tau_22[idx_shear_stress] = (double(4)/double(3)*mu + mu_v)*dvdy - (double(2)/double(3)*mu - mu_v)*dudx;
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                // Get the shear stress cell data.
                
                boost::shared_ptr<pdat::CellData<double> > data_shear_stress(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch->getPatchData(s_variable_shear_stress_unfiltered, data_context)));
                
                // Get the pointers to shear stress data.
                double* tau_11 = data_shear_stress->getPointer(0);
                double* tau_12 = data_shear_stress->getPointer(1);
                double* tau_13 = data_shear_stress->getPointer(2);
                double* tau_22 = data_shear_stress->getPointer(3);
                double* tau_23 = data_shear_stress->getPointer(4);
                double* tau_33 = data_shear_stress->getPointer(5);
                
                // Get the number of ghost cells of the shear stress cell data.
                const hier::IntVector num_ghosts_shear_stress = data_shear_stress->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_shear_stress = data_shear_stress->getGhostBox().numberCells();
                
                const int num_ghosts_0_shear_stress = num_ghosts_shear_stress[0];
                const int num_ghosts_1_shear_stress = num_ghosts_shear_stress[1];
                const int num_ghosts_2_shear_stress = num_ghosts_shear_stress[2];
                const int ghostcell_dim_0_shear_stress = ghostcell_dims_shear_stress[0];
                const int ghostcell_dim_1_shear_stress = ghostcell_dims_shear_stress[1];
                
                // Get the box that covers the interior of patch.
                const hier::Box interior_box = data_shear_stress->getBox();
                const hier::IntVector interior_dims = interior_box.numberCells();
                
                const int interior_dim_0 = interior_dims[0];
                const int interior_dim_1 = interior_dims[1];
                const int interior_dim_2 = interior_dims[2];
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                // HARD CODE TO BE SIXTH ORDER CENTRAL SCHEME FOR DIFFERENTIATION.
                TBOX_ASSERT(num_ghosts >= hier::IntVector::getOne(d_dim)*3);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", num_ghosts));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTION", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("PRESSURE", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("TEMPERATURE", hier::IntVector::getZero(d_dim)));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_velocity =
                    d_flow_model_tmp->getGlobalCellData("VELOCITY");
                
                boost::shared_ptr<pdat::CellData<double> > data_mass_fraction =
                    d_flow_model_tmp->getGlobalCellData("MASS_FRACTION");
                
                boost::shared_ptr<pdat::CellData<double> > data_pressure =
                    d_flow_model_tmp->getGlobalCellData("PRESSURE");
                
                boost::shared_ptr<pdat::CellData<double> > data_temperature =
                    d_flow_model_tmp->getGlobalCellData("TEMPERATURE");
                
                double* u = data_velocity->getPointer(0);
                double* v = data_velocity->getPointer(1);
                double* w = data_velocity->getPointer(2);
                std::vector<double*> Y;
                Y.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    Y.push_back(data_mass_fraction->getPointer(si));
                }
                double* p = data_pressure->getPointer(0);
                double* T = data_temperature->getPointer(0);
                
                const hier::IntVector num_ghosts_velocity = data_velocity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_velocity = data_velocity->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_mass_fraction = data_mass_fraction->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mass_fraction = data_mass_fraction->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_pressure = data_pressure->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_temperature = data_temperature->getGhostBox().numberCells();
                
                const int num_ghosts_0_velocity = num_ghosts_velocity[0];
                const int num_ghosts_1_velocity = num_ghosts_velocity[1];
                const int num_ghosts_2_velocity = num_ghosts_velocity[2];
                const int ghostcell_dim_0_velocity = ghostcell_dims_velocity[0];
                const int ghostcell_dim_1_velocity = ghostcell_dims_velocity[1];
                
                const int num_ghosts_0_mass_fraction = num_ghosts_mass_fraction[0];
                const int num_ghosts_1_mass_fraction = num_ghosts_mass_fraction[1];
                const int num_ghosts_2_mass_fraction = num_ghosts_mass_fraction[2];
                const int ghostcell_dim_0_mass_fraction = ghostcell_dims_mass_fraction[0];
                const int ghostcell_dim_1_mass_fraction = ghostcell_dims_mass_fraction[1];
                
                const int num_ghosts_0_pressure = num_ghosts_pressure[0];
                const int num_ghosts_1_pressure = num_ghosts_pressure[1];
                const int num_ghosts_2_pressure = num_ghosts_pressure[2];
                const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
                const int ghostcell_dim_1_pressure = ghostcell_dims_pressure[1];
                
                const int num_ghosts_0_temperature = num_ghosts_temperature[0];
                const int num_ghosts_1_temperature = num_ghosts_temperature[1];
                const int num_ghosts_2_temperature = num_ghosts_temperature[2];
                const int ghostcell_dim_0_temperature = ghostcell_dims_temperature[0];
                const int ghostcell_dim_1_temperature = ghostcell_dims_temperature[1];
                
                for (int k = 0; k < interior_dim_2; k++)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            /*
                             * Compute the shear stress.
                             */
                                
                            // Compute the linear indices.
                            const int idx_shear_stress = (i + num_ghosts_0_shear_stress) +
                                (j + num_ghosts_1_shear_stress)*ghostcell_dim_0_shear_stress +
                                (k + num_ghosts_2_shear_stress)*ghostcell_dim_0_shear_stress*
                                    ghostcell_dim_1_shear_stress;
                            
                            const int idx_mass_fraction = (i + num_ghosts_0_mass_fraction) +
                                (j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                (k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                    ghostcell_dim_1_mass_fraction;
                            
                            const int idx_pressure = (i + num_ghosts_0_pressure) +
                                (j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure +
                                (k + num_ghosts_2_pressure)*ghostcell_dim_0_pressure*
                                    ghostcell_dim_1_pressure;
                            
                            const int idx_temperature = (i + num_ghosts_0_temperature) +
                                (j + num_ghosts_1_temperature)*ghostcell_dim_0_temperature +
                                (k + num_ghosts_2_temperature)*ghostcell_dim_0_temperature*
                                    ghostcell_dim_1_temperature;
                            
                            const int idx_vel_x_LLL = ((i - 3) + num_ghosts_0_velocity) +
                                (j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                (k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                    ghostcell_dim_1_velocity;
                            
                            const int idx_vel_x_LL  = ((i - 2) + num_ghosts_0_velocity) +
                                (j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                (k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                    ghostcell_dim_1_velocity;
                            
                            const int idx_vel_x_L   = ((i - 1) + num_ghosts_0_velocity) +
                                (j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                (k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                    ghostcell_dim_1_velocity;
                            
                            const int idx_vel_x_R   = ((i + 1) + num_ghosts_0_velocity) +
                                (j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                (k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                    ghostcell_dim_1_velocity;
                            
                            const int idx_vel_x_RR  = ((i + 2) + num_ghosts_0_velocity) +
                                (j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                (k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                    ghostcell_dim_1_velocity;
                            
                            const int idx_vel_x_RRR = ((i + 3) + num_ghosts_0_velocity) +
                                (j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                (k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                    ghostcell_dim_1_velocity;
                            
                            const int idx_vel_y_BBB = (i + num_ghosts_0_velocity) +
                                ((j - 3) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                (k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                    ghostcell_dim_1_velocity;
                            
                            const int idx_vel_y_BB  = (i + num_ghosts_0_velocity) +
                                ((j - 2) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                (k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                    ghostcell_dim_1_velocity;
                            
                            const int idx_vel_y_B   = (i + num_ghosts_0_velocity) +
                                ((j - 1) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                (k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                    ghostcell_dim_1_velocity;
                            
                            const int idx_vel_y_T   = (i + num_ghosts_0_velocity) +
                                ((j + 1) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                (k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                    ghostcell_dim_1_velocity;
                            
                            const int idx_vel_y_TT  = (i + num_ghosts_0_velocity) +
                                ((j + 2) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                (k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                    ghostcell_dim_1_velocity;
                            
                            const int idx_vel_y_TTT = (i + num_ghosts_0_velocity) +
                                ((j + 3) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                (k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                    ghostcell_dim_1_velocity;
                            
                            const int idx_vel_z_BBB = (i + num_ghosts_0_velocity) +
                                (j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                ((k - 3) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                    ghostcell_dim_1_velocity;
                            
                            const int idx_vel_z_BB  = (i + num_ghosts_0_velocity) +
                                (j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                ((k - 2) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                    ghostcell_dim_1_velocity;
                            
                            const int idx_vel_z_B   = (i + num_ghosts_0_velocity) +
                                (j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                ((k - 1) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                    ghostcell_dim_1_velocity;
                            
                            const int idx_vel_z_F   = (i + num_ghosts_0_velocity) +
                                (j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                ((k + 1) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                    ghostcell_dim_1_velocity;
                            
                            const int idx_vel_z_FF  = (i + num_ghosts_0_velocity) +
                                (j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                ((k + 2) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                    ghostcell_dim_1_velocity;
                            
                            const int idx_vel_z_FFF = (i + num_ghosts_0_velocity) +
                                (j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                ((k + 3) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                    ghostcell_dim_1_velocity;
                            
                            std::vector<const double*> Y_ptr;
                            Y_ptr.resize(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_ptr[si] = &Y[si][idx_mass_fraction];
                            }
                            const double mu = d_equation_of_shear_viscosity_mixing_rules->
                                getShearViscosity(
                                    &p[idx_pressure],
                                    &T[idx_temperature],
                                    Y_ptr);
                            const double mu_v = d_equation_of_bulk_viscosity_mixing_rules->
                                getBulkViscosity(
                                    &p[idx_pressure],
                                    &T[idx_temperature],
                                    Y_ptr);
                            
                            const double dudx = (double(1)/double(60)*(u[idx_vel_x_RRR] - u[idx_vel_x_LLL])
                                - double(3)/double(20)*(u[idx_vel_x_RR] - u[idx_vel_x_LL])
                                + double(3)/double(4)*(u[idx_vel_x_R] - u[idx_vel_x_L]))/dx[0];
                            
                            const double dudy = (double(1)/double(60)*(u[idx_vel_y_TTT] - u[idx_vel_y_BBB])
                                - double(3)/double(20)*(u[idx_vel_y_TT] - u[idx_vel_y_BB])
                                + double(3)/double(4)*(u[idx_vel_y_T] - u[idx_vel_y_B]))/dx[1];
                            
                            const double dudz = (double(1)/double(60)*(u[idx_vel_z_FFF] - u[idx_vel_z_BBB])
                                - double(3)/double(20)*(u[idx_vel_z_FF] - u[idx_vel_z_BB])
                                + double(3)/double(4)*(u[idx_vel_z_F] - u[idx_vel_z_B]))/dx[2];
                            
                            const double dvdx = (double(1)/double(60)*(v[idx_vel_x_RRR] - v[idx_vel_x_LLL])
                                - double(3)/double(20)*(v[idx_vel_x_RR] - v[idx_vel_x_LL])
                                + double(3)/double(4)*(v[idx_vel_x_R] - v[idx_vel_x_L]))/dx[0];
                            
                            const double dvdy = (double(1)/double(60)*(v[idx_vel_y_TTT] - v[idx_vel_y_BBB])
                                - double(3)/double(20)*(v[idx_vel_y_TT] - v[idx_vel_y_BB])
                                + double(3)/double(4)*(v[idx_vel_y_T] - v[idx_vel_y_B]))/dx[1];
                            
                            const double dvdz = (double(1)/double(60)*(v[idx_vel_z_FFF] - v[idx_vel_z_BBB])
                                - double(3)/double(20)*(v[idx_vel_z_FF] - v[idx_vel_z_BB])
                                + double(3)/double(4)*(v[idx_vel_z_F] - v[idx_vel_z_B]))/dx[2];
                            
                            const double dwdx = (double(1)/double(60)*(w[idx_vel_x_RRR] - w[idx_vel_x_LLL])
                                - double(3)/double(20)*(w[idx_vel_x_RR] - w[idx_vel_x_LL])
                                + double(3)/double(4)*(w[idx_vel_x_R] - w[idx_vel_x_L]))/dx[0];
                            
                            const double dwdy = (double(1)/double(60)*(w[idx_vel_y_TTT] - w[idx_vel_y_BBB])
                                - double(3)/double(20)*(w[idx_vel_y_TT] - w[idx_vel_y_BB])
                                + double(3)/double(4)*(w[idx_vel_y_T] - w[idx_vel_y_B]))/dx[1];
                            
                            const double dwdz = (double(1)/double(60)*(w[idx_vel_z_FFF] - w[idx_vel_z_BBB])
                                - double(3)/double(20)*(w[idx_vel_z_FF] - w[idx_vel_z_BB])
                                + double(3)/double(4)*(w[idx_vel_z_F] - w[idx_vel_z_B]))/dx[2];
                            
                            tau_11[idx_shear_stress] = (double(4)/double(3)*mu + mu_v)*dudx - (double(2)/double(3)*mu - mu_v)*(dvdy + dwdz);
                            tau_12[idx_shear_stress] = mu*(dudy + dvdx);
                            tau_13[idx_shear_stress] = mu*(dudz + dwdx);
                            tau_22[idx_shear_stress] = (double(4)/double(3)*mu + mu_v)*dvdy - (double(2)/double(3)*mu - mu_v)*(dudx + dwdz);
                            tau_23[idx_shear_stress] = mu*(dvdz + dwdy);
                            tau_33[idx_shear_stress] = (double(4)/double(3)*mu + mu_v)*dwdz - (double(2)/double(3)*mu - mu_v)*(dudx + dvdy);
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
    }
}


/*
 * Compute convective stress with only x direction as inhomogeneous direction.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::computeConvectiveStress(
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
    boost::shared_ptr<FlowModel> d_flow_model_tmp = d_flow_model.lock();
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    if (d_dim == tbox::Dimension(1))
    {
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                // Get the unfiltered convective stress cell data.
                boost::shared_ptr<pdat::CellData<double> > data_convective_stress(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch->getPatchData(s_variable_convective_stress_unfiltered, data_context)));
                
                // Get the pointer to convective stress data.
                double* F_conv_11 = data_convective_stress->getPointer(0);
                
                // Get the number of ghost cells of the convective stress cell data.
                const hier::IntVector num_ghosts_convective_stress = data_convective_stress->getGhostCellWidth();
                const int num_ghosts_0_convective_stress = num_ghosts_convective_stress[0];
                
                // Get the box that covers the interior of patch.
                const hier::Box interior_box = data_convective_stress->getBox();
                const hier::IntVector interior_dims = interior_box.numberCells();
                
                const int interior_dim_0 = interior_dims[0];
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("DENSITY", num_ghosts));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", num_ghosts));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_density =
                    d_flow_model_tmp->getGlobalCellData("DENSITY");
                
                boost::shared_ptr<pdat::CellData<double> > data_velocity =
                    d_flow_model_tmp->getGlobalCellData("VELOCITY");
                
                double* rho = data_density->getPointer(0);
                double* u   = data_velocity->getPointer(0);
                
                const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
                const hier::IntVector num_ghosts_velocity = data_velocity->getGhostCellWidth();
                
                const int num_ghosts_0_density = num_ghosts_density[0];
                const int num_ghosts_0_velocity = num_ghosts_velocity[0];
                
                for (int i = -num_ghosts_0_convective_stress; i < interior_dim_0 + num_ghosts_0_convective_stress; i++)
                {
                    /*
                     * Compute the convective stress.
                     */
                        
                    // Compute the linear indices.
                    const int idx_convective_stress = i + num_ghosts_0_convective_stress;
                    const int idx_density    = i + num_ghosts_0_density;
                    const int idx_velocity   = i + num_ghosts_0_velocity;
                    
                    F_conv_11[idx_convective_stress] = rho[idx_density]*u[idx_velocity]*u[idx_velocity];
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                // Get the unfiltered convective stress cell data.
                boost::shared_ptr<pdat::CellData<double> > data_convective_stress(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch->getPatchData(s_variable_convective_stress_unfiltered, data_context)));
                
                // Get the pointer to convective stress data.
                double* F_conv_11 = data_convective_stress->getPointer(0);
                double* F_conv_12 = data_convective_stress->getPointer(1);
                double* F_conv_22 = data_convective_stress->getPointer(2);
                
                // Get the number of ghost cells of the convective stress cell data.
                const hier::IntVector num_ghosts_convective_stress = data_convective_stress->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_convective_stress = data_convective_stress->getGhostBox().numberCells();
                
                const int num_ghosts_0_convective_stress = num_ghosts_convective_stress[0];
                const int num_ghosts_1_convective_stress = num_ghosts_convective_stress[1];
                const int ghostcell_dim_0_convective_stress = ghostcell_dims_convective_stress[0];
                
                // Get the box that covers the interior of patch.
                const hier::Box interior_box = data_convective_stress->getBox();
                const hier::IntVector interior_dims = interior_box.numberCells();
                
                const int interior_dim_0 = interior_dims[0];
                const int interior_dim_1 = interior_dims[1];
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("DENSITY", num_ghosts));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", num_ghosts));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_density =
                    d_flow_model_tmp->getGlobalCellData("DENSITY");
                
                boost::shared_ptr<pdat::CellData<double> > data_velocity =
                    d_flow_model_tmp->getGlobalCellData("VELOCITY");
                
                double* rho = data_density->getPointer(0);
                double* u   = data_velocity->getPointer(0);
                double* v   = data_velocity->getPointer(1);
                
                const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_density = data_density->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_velocity = data_velocity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_velocity = data_velocity->getGhostBox().numberCells();
                
                const int num_ghosts_0_density = num_ghosts_density[0];
                const int num_ghosts_1_density = num_ghosts_density[1];
                const int ghostcell_dim_0_density = ghostcell_dims_density[0];
                
                const int num_ghosts_0_velocity = num_ghosts_velocity[0];
                const int num_ghosts_1_velocity = num_ghosts_velocity[1];
                const int ghostcell_dim_0_velocity = ghostcell_dims_velocity[0];
                
                for (int j = -num_ghosts_1_convective_stress; j < interior_dim_1 + num_ghosts_1_convective_stress; j++)
                {
                    for (int i = -num_ghosts_0_convective_stress; i < interior_dim_0 + num_ghosts_0_convective_stress; i++)
                    {
                        /*
                         * Compute the convective stress.
                         */
                            
                        // Compute the linear indices.
                        const int idx_convective_stress = i + num_ghosts_0_convective_stress +
                            (j + num_ghosts_1_convective_stress)*ghostcell_dim_0_convective_stress;
                        
                        const int idx_density = i + num_ghosts_0_density +
                            (j + num_ghosts_1_density)*ghostcell_dim_0_density;
                        
                        const int idx_velocity = i + num_ghosts_0_velocity +
                            (j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                        
                        F_conv_11[idx_convective_stress] = rho[idx_density]*u[idx_velocity]*u[idx_velocity];
                        F_conv_12[idx_convective_stress] = rho[idx_density]*u[idx_velocity]*v[idx_velocity];
                        F_conv_22[idx_convective_stress] = rho[idx_density]*v[idx_velocity]*v[idx_velocity];
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                // Get the unfiltered convective stress cell data.
                boost::shared_ptr<pdat::CellData<double> > data_convective_stress(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch->getPatchData(s_variable_convective_stress_unfiltered, data_context)));
                
                // Get the pointers to convective stress data.
                double* F_conv_11 = data_convective_stress->getPointer(0);
                double* F_conv_12 = data_convective_stress->getPointer(1);
                double* F_conv_13 = data_convective_stress->getPointer(2);
                double* F_conv_22 = data_convective_stress->getPointer(3);
                double* F_conv_23 = data_convective_stress->getPointer(4);
                double* F_conv_33 = data_convective_stress->getPointer(5);
                
                // Get the number of ghost cells of the convective stress cell data.
                const hier::IntVector num_ghosts_convective_stress = data_convective_stress->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_convective_stress = data_convective_stress->getGhostBox().numberCells();
                
                const int num_ghosts_0_convective_stress = num_ghosts_convective_stress[0];
                const int num_ghosts_1_convective_stress = num_ghosts_convective_stress[1];
                const int num_ghosts_2_convective_stress = num_ghosts_convective_stress[2];
                const int ghostcell_dim_0_convective_stress = ghostcell_dims_convective_stress[0];
                const int ghostcell_dim_1_convective_stress = ghostcell_dims_convective_stress[1];
                
                // Get the box that covers the interior of patch.
                const hier::Box interior_box = data_convective_stress->getBox();
                const hier::IntVector interior_dims = interior_box.numberCells();
                
                const int interior_dim_0 = interior_dims[0];
                const int interior_dim_1 = interior_dims[1];
                const int interior_dim_2 = interior_dims[2];
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("DENSITY", num_ghosts));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", num_ghosts));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_density =
                    d_flow_model_tmp->getGlobalCellData("DENSITY");
                
                boost::shared_ptr<pdat::CellData<double> > data_velocity =
                    d_flow_model_tmp->getGlobalCellData("VELOCITY");
                
                double* rho = data_density->getPointer(0);
                double* u   = data_velocity->getPointer(0);
                double* v   = data_velocity->getPointer(1);
                double* w   = data_velocity->getPointer(2);
                
                const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_density = data_density->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_velocity = data_velocity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_velocity = data_velocity->getGhostBox().numberCells();
                
                const int num_ghosts_0_density = num_ghosts_density[0];
                const int num_ghosts_1_density = num_ghosts_density[1];
                const int num_ghosts_2_density = num_ghosts_density[2];
                const int ghostcell_dim_0_density = ghostcell_dims_density[0];
                const int ghostcell_dim_1_density = ghostcell_dims_density[1];
                
                const int num_ghosts_0_velocity = num_ghosts_velocity[0];
                const int num_ghosts_1_velocity = num_ghosts_velocity[1];
                const int num_ghosts_2_velocity = num_ghosts_velocity[2];
                const int ghostcell_dim_0_velocity = ghostcell_dims_velocity[0];
                const int ghostcell_dim_1_velocity = ghostcell_dims_velocity[1];
                
                for (int k = -num_ghosts_2_convective_stress; k < interior_dim_2 + num_ghosts_2_convective_stress; k++)
                {
                    for (int j = -num_ghosts_1_convective_stress; j < interior_dim_1 + num_ghosts_1_convective_stress; j++)
                    {
                        for (int i = -num_ghosts_0_convective_stress; i < interior_dim_0 + num_ghosts_0_convective_stress; i++)
                        {
                            /*
                             * Compute the convective stress.
                             */
                                
                            // Compute the linear indices.
                            const int idx_convective_stress = (i + num_ghosts_0_convective_stress) +
                                (j + num_ghosts_1_convective_stress)*ghostcell_dim_0_convective_stress +
                                (k + num_ghosts_2_convective_stress)*ghostcell_dim_0_convective_stress*
                                    ghostcell_dim_1_convective_stress;
                            
                            const int idx_density = (i + num_ghosts_0_density) +
                                (j + num_ghosts_1_density)*ghostcell_dim_0_density +
                                (k + num_ghosts_2_density)*ghostcell_dim_0_density*
                                    ghostcell_dim_1_density;
                            
                            const int idx_velocity = (i + num_ghosts_0_velocity) +
                                (j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                (k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                    ghostcell_dim_1_velocity;
                            
                            F_conv_11[idx_convective_stress] = rho[idx_density]*u[idx_velocity]*u[idx_velocity];
                            F_conv_12[idx_convective_stress] = rho[idx_density]*u[idx_velocity]*v[idx_velocity];
                            F_conv_13[idx_convective_stress] = rho[idx_density]*u[idx_velocity]*w[idx_velocity];
                            F_conv_22[idx_convective_stress] = rho[idx_density]*v[idx_velocity]*v[idx_velocity];
                            F_conv_23[idx_convective_stress] = rho[idx_density]*v[idx_velocity]*w[idx_velocity];
                            F_conv_33[idx_convective_stress] = rho[idx_density]*w[idx_velocity]*w[idx_velocity];
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
    }
}


/*
 * Filter pressure.
 */
void 
FlowModelStatisticsUtilitiesFourEqnConservative::filterPressure(
    const int level,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    TBOX_ASSERT(level < num_levels);
    
    /*
     * Get the patch level.
     */
    
    boost::shared_ptr<hier::PatchLevel> patch_level(
        patch_hierarchy->getPatchLevel(level));
    
    for (hier::PatchLevel::iterator ip(patch_level->begin());
         ip != patch_level->end();
         ip++)
    {
        const boost::shared_ptr<hier::Patch> patch = *ip;
        
        // Get the unfiltered and filtered pressure cell data.
        
        boost::shared_ptr<pdat::CellData<double> > data_pressure_unfiltered(
            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                patch->getPatchData(s_variable_pressure_unfiltered, data_context)));
        
        boost::shared_ptr<pdat::CellData<double> > data_pressure_filtered(
            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                patch->getPatchData(s_variable_pressure_filtered, data_context)));
        
        if (d_pressure_filtered)
        {
            data_pressure_unfiltered->copy(*data_pressure_filtered);
        }
        
        // Apply filter in x-direction.
        
        d_filter_x->applyFilter(
            data_pressure_filtered,
            data_pressure_unfiltered,
            0,
            0);
        
        if ((d_dim == tbox::Dimension(2)) || (d_dim == tbox::Dimension(3)))
        {
            // Apply filter in y-direction.
            
            data_pressure_unfiltered->copy(*data_pressure_filtered);
            
            d_filter_y->applyFilter(
                data_pressure_filtered,
                data_pressure_unfiltered,
                0,
                0);
        }
        
        if (d_dim == tbox::Dimension(3))
        {
            // Apply filter in z-direction.
            
            data_pressure_unfiltered->copy(*data_pressure_filtered);
            
            d_filter_z->applyFilter(
                data_pressure_filtered,
                data_pressure_unfiltered,
                0,
                0);
        }
    }
    
    if (level == 0)
    {
        d_pressure_filtered = true;
    }
}


/*
 * Filter shear stress.
 */
void 
FlowModelStatisticsUtilitiesFourEqnConservative::filterShearStress(
    const int level,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    TBOX_ASSERT(level < num_levels);
    
    /*
     * Get the patch level.
     */
    
    boost::shared_ptr<hier::PatchLevel> patch_level(
        patch_hierarchy->getPatchLevel(level));
    
    for (hier::PatchLevel::iterator ip(patch_level->begin());
         ip != patch_level->end();
         ip++)
    {
        const boost::shared_ptr<hier::Patch> patch = *ip;
        
        // Get the unfiltered and filtered shear stress cell data.
        
        boost::shared_ptr<pdat::CellData<double> > data_shear_stress_unfiltered(
            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                patch->getPatchData(s_variable_shear_stress_unfiltered, data_context)));
        
        boost::shared_ptr<pdat::CellData<double> > data_shear_stress_filtered(
            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                patch->getPatchData(s_variable_shear_stress_filtered, data_context)));
        
        if (d_shear_stress_filtered)
        {
            data_shear_stress_unfiltered->copy(*data_shear_stress_filtered);
        }
        
        if (d_dim == tbox::Dimension(1))
        {
            for (int di = 0; di < 1; di++)
            {
                // Apply filter in x-direction.
                
                d_filter_x->applyFilter(
                    data_shear_stress_filtered,
                    data_shear_stress_unfiltered,
                    di,
                    di);
            }
        }
        else if (d_dim == tbox::Dimension(2))
        {
            // Apply filter in x-direction.
            
            for (int di = 0; di < 3; di++)
            {
                d_filter_x->applyFilter(
                    data_shear_stress_filtered,
                    data_shear_stress_unfiltered,
                    di,
                    di);
            }
            
            // Apply filter in y-direction.
            
            data_shear_stress_unfiltered->copy(*data_shear_stress_filtered);
            
            for (int di = 0; di < 3; di++)
            {
                d_filter_y->applyFilter(
                    data_shear_stress_filtered,
                    data_shear_stress_unfiltered,
                    di,
                    di);
            }
        }
        else if (d_dim == tbox::Dimension(3))
        {
            // Apply filter in x-direction.
            
            for (int di = 0; di < 6; di++)
            {
                d_filter_x->applyFilter(
                    data_shear_stress_filtered,
                    data_shear_stress_unfiltered,
                    di,
                    di);
            }
            
            // Apply filter in y-direction.
            
            data_shear_stress_unfiltered->copy(*data_shear_stress_filtered);
            
            for (int di = 0; di < 6; di++)
            {
                d_filter_y->applyFilter(
                    data_shear_stress_filtered,
                    data_shear_stress_unfiltered,
                    di,
                    di);
            }
            
            // Apply filter in z-direction.
            
            data_shear_stress_unfiltered->copy(*data_shear_stress_filtered);
            
            for (int di = 0; di < 6; di++)
            {
                d_filter_z->applyFilter(
                    data_shear_stress_filtered,
                    data_shear_stress_unfiltered,
                    di,
                    di);
            }
        }
    }
    
    if (level == 0)
    {
        d_shear_stress_filtered = true;
    }
}


/*
 * Filter convective stress.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::filterConvectiveStress(
    const int level,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    TBOX_ASSERT(level < num_levels);
    
    /*
     * Get the patch level.
     */
    
    boost::shared_ptr<hier::PatchLevel> patch_level(
        patch_hierarchy->getPatchLevel(level));
    
    for (hier::PatchLevel::iterator ip(patch_level->begin());
         ip != patch_level->end();
         ip++)
    {
        const boost::shared_ptr<hier::Patch> patch = *ip;
        
        // Get the unfiltered and filtered convective stress cell data.
        
        boost::shared_ptr<pdat::CellData<double> > data_convective_stress_unfiltered(
            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                patch->getPatchData(s_variable_convective_stress_unfiltered, data_context)));
        
        boost::shared_ptr<pdat::CellData<double> > data_convective_stress_filtered(
            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                patch->getPatchData(s_variable_convective_stress_filtered, data_context)));
        
        if (d_convective_stress_filtered)
        {
            data_convective_stress_unfiltered->copy(*data_convective_stress_filtered);
        }
        
        if (d_dim == tbox::Dimension(1))
        {
            for (int di = 0; di < 1; di++)
            {
                // Apply filter in x-direction.
                
                d_filter_x->applyFilter(
                    data_convective_stress_filtered,
                    data_convective_stress_unfiltered,
                    di,
                    di);
            }
        }
        else if (d_dim == tbox::Dimension(2))
        {
            for (int di = 0; di < 3; di++)
            {
                // Apply filter in x-direction.
                
                d_filter_x->applyFilter(
                    data_convective_stress_filtered,
                    data_convective_stress_unfiltered,
                    di,
                    di);
            }
            
            // Apply filter in y-direction.
            
            data_convective_stress_unfiltered->copy(*data_convective_stress_filtered);
            
            for (int di = 0; di < 3; di++)
            {
                d_filter_y->applyFilter(
                    data_convective_stress_filtered,
                    data_convective_stress_unfiltered,
                    di,
                    di);
            }
        }
        else if (d_dim == tbox::Dimension(3))
        {
            for (int di = 0; di < 6; di++)
            {
                // Apply filter in x-direction.
                
                d_filter_x->applyFilter(
                    data_convective_stress_filtered,
                    data_convective_stress_unfiltered,
                    di,
                    di);
            }
            
            // Apply filter in y-direction.
            
            data_convective_stress_unfiltered->copy(*data_convective_stress_filtered);
            
            for (int di = 0; di < 6; di++)
            {
                d_filter_y->applyFilter(
                    data_convective_stress_filtered,
                    data_convective_stress_unfiltered,
                    di,
                    di);
            }
            
            // Apply filter in z-direction.
            
            data_convective_stress_unfiltered->copy(*data_convective_stress_filtered);
            
            for (int di = 0; di < 6; di++)
            {
                d_filter_z->applyFilter(
                    data_convective_stress_filtered,
                    data_convective_stress_unfiltered,
                    di,
                    di);
            }
        }
    }
    
    if (level == 0)
    {
        d_convective_stress_filtered = true;
    }
}


/*
 * Compute SFS stress.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::computeStressSFS(
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
    boost::shared_ptr<FlowModel> d_flow_model_tmp = d_flow_model.lock();
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    if (d_dim == tbox::Dimension(1))
    {
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                // Get the SFS stress cell data.
                
                boost::shared_ptr<pdat::CellData<double> > data_SFS_stress(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch->getPatchData(s_variable_SFS_stress, data_context)));
                
                boost::shared_ptr<pdat::CellData<double> > data_convective_stress_filtered(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch->getPatchData(s_variable_convective_stress_filtered, data_context)));
                
                // Get the pointer to SFS stress data.
                double* tau_SFS_11 = data_SFS_stress->getPointer(0);
                
                // Get the number of ghost cells of the SFS stress cell data.
                const hier::IntVector num_ghosts_SFS_stress = data_SFS_stress->getGhostCellWidth();
                const int num_ghosts_0_SFS_stress = num_ghosts_SFS_stress[0];
                
                // Get the box that covers the interior of patch.
                const hier::Box interior_box = data_SFS_stress->getBox();
                const hier::IntVector interior_dims = interior_box.numberCells();
                
                const int interior_dim_0 = interior_dims[0];
                
                // Get the pointer to filtered convective stress data.
                double* F_conv_filtered_11 = data_convective_stress_filtered->getPointer(0);
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("DENSITY", num_ghosts));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", num_ghosts));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_density =
                    d_flow_model_tmp->getGlobalCellData("DENSITY");
                
                boost::shared_ptr<pdat::CellData<double> > data_velocity =
                    d_flow_model_tmp->getGlobalCellData("VELOCITY");
                
                double* rho = data_density->getPointer(0);
                double* u   = data_velocity->getPointer(0);
                
                const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
                const hier::IntVector num_ghosts_velocity = data_velocity->getGhostCellWidth();
                
                const int num_ghosts_0_density = num_ghosts_density[0];
                const int num_ghosts_0_velocity = num_ghosts_velocity[0];
                
                for (int i = -num_ghosts_0_SFS_stress; i < interior_dim_0 + num_ghosts_0_SFS_stress; i++)
                {
                    /*
                        * Compute the convective stress.
                        */
                        
                    // Compute the linear indices.
                    const int idx_SFS_stress = i + num_ghosts_0_SFS_stress;
                    const int idx_density    = i + num_ghosts_0_density;
                    const int idx_velocity   = i + num_ghosts_0_velocity;
                    
                    tau_SFS_11[idx_SFS_stress] = F_conv_filtered_11[idx_SFS_stress] - rho[idx_density]*u[idx_velocity]*u[idx_velocity];
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                // Get the SFS stress cell data.
                
                boost::shared_ptr<pdat::CellData<double> > data_SFS_stress(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch->getPatchData(s_variable_SFS_stress, data_context)));
                
                boost::shared_ptr<pdat::CellData<double> > data_convective_stress_filtered(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch->getPatchData(s_variable_convective_stress_filtered, data_context)));
                
                // Get the pointers to SFS stress data.
                double* tau_SFS_11 = data_SFS_stress->getPointer(0);
                double* tau_SFS_12 = data_SFS_stress->getPointer(1);
                double* tau_SFS_22 = data_SFS_stress->getPointer(2);
                
                // Get the number of ghost cells of the SFS stress cell data.
                const hier::IntVector num_ghosts_SFS_stress = data_SFS_stress->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_SFS_stress = data_SFS_stress->getGhostBox().numberCells();
                
                const int num_ghosts_0_SFS_stress = num_ghosts_SFS_stress[0];
                const int num_ghosts_1_SFS_stress = num_ghosts_SFS_stress[1];
                const int ghostcell_dim_0_SFS_stress = ghostcell_dims_SFS_stress[0];
                
                // Get the box that covers the interior of patch.
                const hier::Box interior_box = data_SFS_stress->getBox();
                const hier::IntVector interior_dims = interior_box.numberCells();
                
                const int interior_dim_0 = interior_dims[0];
                const int interior_dim_1 = interior_dims[1];
                
                // Get the pointers to filtered convective stress data.
                double* F_conv_filtered_11 = data_convective_stress_filtered->getPointer(0);
                double* F_conv_filtered_12 = data_convective_stress_filtered->getPointer(1);
                double* F_conv_filtered_22 = data_convective_stress_filtered->getPointer(2);
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("DENSITY", num_ghosts));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", num_ghosts));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_density =
                    d_flow_model_tmp->getGlobalCellData("DENSITY");
                
                boost::shared_ptr<pdat::CellData<double> > data_velocity =
                    d_flow_model_tmp->getGlobalCellData("VELOCITY");
                
                double* rho = data_density->getPointer(0);
                double* u   = data_velocity->getPointer(0);
                double* v   = data_velocity->getPointer(1);
                
                const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_density = data_density->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_velocity = data_velocity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_velocity = data_velocity->getGhostBox().numberCells();
                
                const int num_ghosts_0_density = num_ghosts_density[0];
                const int num_ghosts_1_density = num_ghosts_density[1];
                const int ghostcell_dim_0_density = ghostcell_dims_density[0];
                
                const int num_ghosts_0_velocity = num_ghosts_velocity[0];
                const int num_ghosts_1_velocity = num_ghosts_velocity[1];
                const int ghostcell_dim_0_velocity = ghostcell_dims_velocity[0];
                
                for (int j = -num_ghosts_1_SFS_stress; j < interior_dim_1 + num_ghosts_1_SFS_stress; j++)
                {
                    for (int i = -num_ghosts_0_SFS_stress; i < interior_dim_0 + num_ghosts_0_SFS_stress; i++)
                    {
                        /*
                         * Compute the convective stress.
                         */
                            
                        // Compute the linear indices.
                        const int idx_SFS_stress = i + num_ghosts_0_SFS_stress +
                            (j + num_ghosts_1_SFS_stress)*ghostcell_dim_0_SFS_stress;
                        
                        const int idx_density = i + num_ghosts_0_density +
                            (j + num_ghosts_1_density)*ghostcell_dim_0_density;
                        
                        const int idx_velocity = i + num_ghosts_0_velocity +
                            (j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                        
                        tau_SFS_11[idx_SFS_stress] = F_conv_filtered_11[idx_SFS_stress] - rho[idx_density]*u[idx_velocity]*u[idx_velocity];
                        tau_SFS_12[idx_SFS_stress] = F_conv_filtered_12[idx_SFS_stress] - rho[idx_density]*u[idx_velocity]*v[idx_velocity];
                        tau_SFS_22[idx_SFS_stress] = F_conv_filtered_22[idx_SFS_stress] - rho[idx_density]*v[idx_velocity]*v[idx_velocity];
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                // Get the SFS stress cell data.
                
                boost::shared_ptr<pdat::CellData<double> > data_SFS_stress(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch->getPatchData(s_variable_SFS_stress, data_context)));
                
                boost::shared_ptr<pdat::CellData<double> > data_convective_stress_filtered(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch->getPatchData(s_variable_convective_stress_filtered, data_context)));
                
                // Get the pointers to SFS stress data.
                double* tau_SFS_11 = data_SFS_stress->getPointer(0);
                double* tau_SFS_12 = data_SFS_stress->getPointer(1);
                double* tau_SFS_13 = data_SFS_stress->getPointer(2);
                double* tau_SFS_22 = data_SFS_stress->getPointer(3);
                double* tau_SFS_23 = data_SFS_stress->getPointer(4);
                double* tau_SFS_33 = data_SFS_stress->getPointer(5);
                
                // Get the number of ghost cells of the SFS stress cell data.
                const hier::IntVector num_ghosts_SFS_stress = data_SFS_stress->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_SFS_stress = data_SFS_stress->getGhostBox().numberCells();
                
                const int num_ghosts_0_SFS_stress = num_ghosts_SFS_stress[0];
                const int num_ghosts_1_SFS_stress = num_ghosts_SFS_stress[1];
                const int num_ghosts_2_SFS_stress = num_ghosts_SFS_stress[2];
                const int ghostcell_dim_0_SFS_stress = ghostcell_dims_SFS_stress[0];
                const int ghostcell_dim_1_SFS_stress = ghostcell_dims_SFS_stress[1];
                
                // Get the box that covers the interior of patch.
                const hier::Box interior_box = data_SFS_stress->getBox();
                const hier::IntVector interior_dims = interior_box.numberCells();
                
                const int interior_dim_0 = interior_dims[0];
                const int interior_dim_1 = interior_dims[1];
                const int interior_dim_2 = interior_dims[2];
                
                // Get the pointers to filtered convective stress data.
                double* F_conv_filtered_11 = data_convective_stress_filtered->getPointer(0);
                double* F_conv_filtered_12 = data_convective_stress_filtered->getPointer(1);
                double* F_conv_filtered_13 = data_convective_stress_filtered->getPointer(2);
                double* F_conv_filtered_22 = data_convective_stress_filtered->getPointer(3);
                double* F_conv_filtered_23 = data_convective_stress_filtered->getPointer(4);
                double* F_conv_filtered_33 = data_convective_stress_filtered->getPointer(5);
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("DENSITY", num_ghosts));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", num_ghosts));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_density =
                    d_flow_model_tmp->getGlobalCellData("DENSITY");
                
                boost::shared_ptr<pdat::CellData<double> > data_velocity =
                    d_flow_model_tmp->getGlobalCellData("VELOCITY");
                
                double* rho = data_density->getPointer(0);
                double* u   = data_velocity->getPointer(0);
                double* v   = data_velocity->getPointer(1);
                double* w   = data_velocity->getPointer(2);
                
                const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_density = data_density->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_velocity = data_velocity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_velocity = data_velocity->getGhostBox().numberCells();
                
                const int num_ghosts_0_density = num_ghosts_density[0];
                const int num_ghosts_1_density = num_ghosts_density[1];
                const int num_ghosts_2_density = num_ghosts_density[2];
                const int ghostcell_dim_0_density = ghostcell_dims_density[0];
                const int ghostcell_dim_1_density = ghostcell_dims_density[1];
                
                const int num_ghosts_0_velocity = num_ghosts_velocity[0];
                const int num_ghosts_1_velocity = num_ghosts_velocity[1];
                const int num_ghosts_2_velocity = num_ghosts_velocity[2];
                const int ghostcell_dim_0_velocity = ghostcell_dims_velocity[0];
                const int ghostcell_dim_1_velocity = ghostcell_dims_velocity[1];
                
                for (int k = -num_ghosts_2_SFS_stress; k < interior_dim_2 + num_ghosts_2_SFS_stress; k++)
                {
                    for (int j = -num_ghosts_1_SFS_stress; j < interior_dim_1 + num_ghosts_1_SFS_stress; j++)
                    {
                        for (int i = -num_ghosts_0_SFS_stress; i < interior_dim_0 + num_ghosts_0_SFS_stress; i++)
                        {
                            /*
                             * Compute the convective stress.
                             */
                                
                            // Compute the linear indices.
                            const int idx_SFS_stress = (i + num_ghosts_0_SFS_stress) +
                                (j + num_ghosts_1_SFS_stress)*ghostcell_dim_0_SFS_stress +
                                (k + num_ghosts_2_SFS_stress)*ghostcell_dim_0_SFS_stress*
                                    ghostcell_dim_1_SFS_stress;
                            
                            const int idx_density = (i + num_ghosts_0_density) +
                                (j + num_ghosts_1_density)*ghostcell_dim_0_density +
                                (k + num_ghosts_2_density)*ghostcell_dim_0_density*
                                    ghostcell_dim_1_density;
                            
                            const int idx_velocity = (i + num_ghosts_0_velocity) +
                                (j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                (k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                    ghostcell_dim_1_velocity;
                            
                            tau_SFS_11[idx_SFS_stress] = F_conv_filtered_11[idx_SFS_stress] - rho[idx_density]*u[idx_velocity]*u[idx_velocity];
                            tau_SFS_12[idx_SFS_stress] = F_conv_filtered_12[idx_SFS_stress] - rho[idx_density]*u[idx_velocity]*v[idx_velocity];
                            tau_SFS_13[idx_SFS_stress] = F_conv_filtered_13[idx_SFS_stress] - rho[idx_density]*u[idx_velocity]*w[idx_velocity];
                            tau_SFS_22[idx_SFS_stress] = F_conv_filtered_22[idx_SFS_stress] - rho[idx_density]*v[idx_velocity]*v[idx_velocity];
                            tau_SFS_23[idx_SFS_stress] = F_conv_filtered_23[idx_SFS_stress] - rho[idx_density]*v[idx_velocity]*w[idx_velocity];
                            tau_SFS_33[idx_SFS_stress] = F_conv_filtered_33[idx_SFS_stress] - rho[idx_density]*w[idx_velocity]*w[idx_velocity];
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
    }
}


/*
 * Compute Favre-filtered velocity.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::computeFavreFilteredVelocity(
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
    boost::shared_ptr<FlowModel> d_flow_model_tmp = d_flow_model.lock();
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    for (int li = 0; li < num_levels; li++)
    {
        /*
         * Get the current patch level.
         */
        
        boost::shared_ptr<hier::PatchLevel> patch_level(
            patch_hierarchy->getPatchLevel(li));
        
        for (hier::PatchLevel::iterator ip(patch_level->begin());
             ip != patch_level->end();
             ip++)
        {
            const boost::shared_ptr<hier::Patch> patch = *ip;
            
            // Get the Favre-filtered velocity cell data.
            
            boost::shared_ptr<pdat::CellData<double> > data_velocity_Favre_filtered(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch->getPatchData(s_variable_velocity_Favre_filtered, data_context)));
            
            /*
             * Register the patch and the quantity in the flow model and compute the
             * corresponding cell data.
             */
            
            d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
            
            hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
            
            std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
            
            num_subghosts_of_data.insert(
                std::pair<std::string, hier::IntVector>("VELOCITY", num_ghosts));
            
            d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
            
            d_flow_model_tmp->computeGlobalDerivedCellData();
            
            boost::shared_ptr<pdat::CellData<double> > data_velocity =
                d_flow_model_tmp->getGlobalCellData("VELOCITY");
            
            data_velocity_Favre_filtered->copy(*data_velocity);
            
            /*
             * Unregister the patch and data of all registered derived cell variables in the flow model.
             */
            
            d_flow_model_tmp->unregisterPatch();
        }
    }
}


/*
 * Compute Favre-filtered specific volume.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::computeFavreFilteredSpecificVolume(
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
    boost::shared_ptr<FlowModel> d_flow_model_tmp = d_flow_model.lock();
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    if (d_dim == tbox::Dimension(1))
    {
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                // Get the Favre-filtered specific volume cell data.
                
                boost::shared_ptr<pdat::CellData<double> > data_specific_volume(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch->getPatchData(s_variable_specific_volume_Favre_filtered, data_context)));
                
                // Get the pointer to specific volume cell data.
                double* rho_inv = data_specific_volume->getPointer(0);
                
                // Get the number of ghost cells of the specific volume cell data.
                const hier::IntVector num_ghosts_specific_volume = data_specific_volume->getGhostCellWidth();
                const int num_ghosts_0_specific_volume = num_ghosts_specific_volume[0];
                
                // Get the box that covers the interior of patch.
                const hier::Box interior_box = data_specific_volume->getBox();
                const hier::IntVector interior_dims = interior_box.numberCells();
                
                const int interior_dim_0 = interior_dims[0];
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("DENSITY", num_ghosts));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_density =
                    d_flow_model_tmp->getGlobalCellData("DENSITY");
                
                double* rho = data_density->getPointer(0);
                
                const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
                
                const int num_ghosts_0_density = num_ghosts_density[0];
                
                for (int i = -num_ghosts_0_specific_volume; i < interior_dim_0 + num_ghosts_0_specific_volume; i++)
                {
                    /*
                     * Compute the specific volume.
                     */
                        
                    // Compute the linear indices.
                    const int idx_specific_volume = i + num_ghosts_0_specific_volume;
                    const int idx_density = i + num_ghosts_0_density;
                    
                    rho_inv[idx_specific_volume] = double(1)/rho[idx_density];
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                // Get the Favre-filtered specific volume cell data.
                
                boost::shared_ptr<pdat::CellData<double> > data_specific_volume(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch->getPatchData(s_variable_specific_volume_Favre_filtered, data_context)));
                
                // Get the pointer to specific volume cell data.
                double* rho_inv = data_specific_volume->getPointer(0);
                
                // Get the number of ghost cells of the specific volume cell data.
                const hier::IntVector num_ghosts_specific_volume = data_specific_volume->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_specific_volume = data_specific_volume->getGhostBox().numberCells();
                
                const int num_ghosts_0_specific_volume = num_ghosts_specific_volume[0];
                const int num_ghosts_1_specific_volume = num_ghosts_specific_volume[1];
                const int ghostcell_dim_0_specific_volume = ghostcell_dims_specific_volume[0];
                
                // Get the box that covers the interior of patch.
                const hier::Box interior_box = data_specific_volume->getBox();
                const hier::IntVector interior_dims = interior_box.numberCells();
                
                const int interior_dim_0 = interior_dims[0];
                const int interior_dim_1 = interior_dims[1];
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("DENSITY", num_ghosts));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_density =
                    d_flow_model_tmp->getGlobalCellData("DENSITY");
                
                double* rho = data_density->getPointer(0);
                
                const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_density = data_density->getGhostBox().numberCells();
                
                const int num_ghosts_0_density = num_ghosts_density[0];
                const int num_ghosts_1_density = num_ghosts_density[1];
                const int ghostcell_dim_0_density = ghostcell_dims_density[0];
                
                for (int j = -num_ghosts_1_specific_volume; j < interior_dim_1 + num_ghosts_1_specific_volume; j++)
                {
                    for (int i = -num_ghosts_0_specific_volume; i < interior_dim_0 + num_ghosts_0_specific_volume; i++)
                    {
                        /*
                         * Compute the specific volume.
                         */
                            
                        // Compute the linear indices.
                        const int idx_specific_volume = i + num_ghosts_0_specific_volume +
                            (j + num_ghosts_1_specific_volume)*ghostcell_dim_0_specific_volume;
                        
                        const int idx_density = i + num_ghosts_0_density +
                            (j + num_ghosts_1_density)*ghostcell_dim_0_density;
                        
                        rho_inv[idx_specific_volume] = double(1)/rho[idx_density];
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                // Get the Favre-filtered specific volume cell data.
                
                boost::shared_ptr<pdat::CellData<double> > data_specific_volume(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch->getPatchData(s_variable_specific_volume_Favre_filtered, data_context)));
                
                // Get the pointer to specific volume cell data.
                double* rho_inv = data_specific_volume->getPointer(0);
                
                // Get the number of ghost cells of the specific volume cell data.
                const hier::IntVector num_ghosts_specific_volume = data_specific_volume->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_specific_volume = data_specific_volume->getGhostBox().numberCells();
                
                const int num_ghosts_0_specific_volume = num_ghosts_specific_volume[0];
                const int num_ghosts_1_specific_volume = num_ghosts_specific_volume[1];
                const int num_ghosts_2_specific_volume = num_ghosts_specific_volume[2];
                const int ghostcell_dim_0_specific_volume = ghostcell_dims_specific_volume[0];
                const int ghostcell_dim_1_specific_volume = ghostcell_dims_specific_volume[1];
                
                // Get the box that covers the interior of patch.
                const hier::Box interior_box = data_specific_volume->getBox();
                const hier::IntVector interior_dims = interior_box.numberCells();
                
                const int interior_dim_0 = interior_dims[0];
                const int interior_dim_1 = interior_dims[1];
                const int interior_dim_2 = interior_dims[2];
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("DENSITY", num_ghosts));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_density =
                    d_flow_model_tmp->getGlobalCellData("DENSITY");
                
                double* rho = data_density->getPointer(0);
                
                const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_density = data_density->getGhostBox().numberCells();
                
                const int num_ghosts_0_density = num_ghosts_density[0];
                const int num_ghosts_1_density = num_ghosts_density[1];
                const int num_ghosts_2_density = num_ghosts_density[2];
                const int ghostcell_dim_0_density = ghostcell_dims_density[0];
                const int ghostcell_dim_1_density = ghostcell_dims_density[1];
                
                for (int k = -num_ghosts_2_specific_volume; k < interior_dim_2 + num_ghosts_2_specific_volume; k++)
                {
                    for (int j = -num_ghosts_1_specific_volume; j < interior_dim_1 + num_ghosts_1_specific_volume; j++)
                    {
                        for (int i = -num_ghosts_0_specific_volume; i < interior_dim_0 + num_ghosts_0_specific_volume; i++)
                        {
                            /*
                             * Compute the specific volume.
                             */
                                
                            // Compute the linear indices.
                            const int idx_specific_volume = (i + num_ghosts_0_specific_volume) +
                                (j + num_ghosts_1_specific_volume)*ghostcell_dim_0_specific_volume +
                                (k + num_ghosts_2_specific_volume)*ghostcell_dim_0_specific_volume*
                                    ghostcell_dim_1_specific_volume;
                            
                            const int idx_density = (i + num_ghosts_0_density) +
                                (j + num_ghosts_1_density)*ghostcell_dim_0_density +
                                (k + num_ghosts_2_density)*ghostcell_dim_0_density*
                                    ghostcell_dim_1_density;
                            
                            rho_inv[idx_specific_volume] = double(1)/rho[idx_density];
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
    }
}


/*
 * Compute filtered density.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::computeFilteredDensity(
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
    boost::shared_ptr<FlowModel> d_flow_model_tmp = d_flow_model.lock();
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    for (int li = 0; li < num_levels; li++)
    {
        /*
         * Get the current patch level.
         */
        
        boost::shared_ptr<hier::PatchLevel> patch_level(
            patch_hierarchy->getPatchLevel(li));
        
        for (hier::PatchLevel::iterator ip(patch_level->begin());
             ip != patch_level->end();
             ip++)
        {
            const boost::shared_ptr<hier::Patch> patch = *ip;
            
            // Get the filtered density cell data.
            
            boost::shared_ptr<pdat::CellData<double> > data_density_filtered(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch->getPatchData(s_variable_density_filtered, data_context)));
            
            /*
             * Register the patch and the quantity in the flow model and compute the
             * corresponding cell data.
             */
            
            d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
            
            hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
            
            std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
            
            num_subghosts_of_data.insert(
                std::pair<std::string, hier::IntVector>("DENSITY", num_ghosts));
            
            d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
            
            d_flow_model_tmp->computeGlobalDerivedCellData();
            
            /*
             * Get the pointers to data inside the flow model.
             */
            
            boost::shared_ptr<pdat::CellData<double> > data_density =
                d_flow_model_tmp->getGlobalCellData("DENSITY");
            
            data_density_filtered->copy(*data_density);
            
            /*
             * Unregister the patch and data of all registered derived cell variables in the flow model.
             */
            
            d_flow_model_tmp->unregisterPatch();
        }
    }
}


/*
 * Compute averaged value with only x direction as inhomogeneous direction.
 */
std::vector<double>
FlowModelStatisticsUtilitiesFourEqnConservative::getAveragedQuantityWithInhomogeneousXDirection(
    boost::shared_ptr<pdat::CellVariable<double> >& variable_quantity,
    const int component_idx,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context) const
{
    std::vector<double> averaged_quantity;
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Get the refinement ratio from the finest level to the coarest level.
     */
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    hier::IntVector ratioFinestLevelToCoarestLevel =
        patch_hierarchy->getRatioToCoarserLevel(num_levels - 1);
    for (int li = num_levels - 2; li > 0 ; li--)
    {
        ratioFinestLevelToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(li);
    }
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    boost::shared_ptr<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *patch_hierarchy,
            0,
            num_levels - 1));
    
    /*
     * Get the number of cells of physical domain refined to the finest level.
     */
    
    const hier::BoxContainer& physical_domain = d_grid_geometry->getPhysicalDomain();
    const hier::Box& physical_domain_box = physical_domain.front();
    const hier::IntVector& physical_domain_dims = physical_domain_box.numberCells();
    const hier::IntVector finest_level_dims = physical_domain_dims*ratioFinestLevelToCoarestLevel;
    
    /*
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        double* u_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_quantity.resize(finest_level_dim_0);
        double* u_avg_global = averaged_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_avg_local[i] = double(0);
            u_avg_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                boost::shared_ptr<pdat::CellData<double> > data_quantity(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch->getPatchData(variable_quantity, data_context)));
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        /*
                         * Compute the index of the data point and count how many times the data is repeated.
                         */
                        
                        const hier::Index idx_pt(tbox::Dimension(1), idx_lo_0 + i);
                        
                        int n_overlapped = 1;
                        
                        for (hier::BoxContainer::BoxContainerConstIterator iob(
                                patch_overlapped_visible_boxes.begin());
                             iob != patch_overlapped_visible_boxes.end();
                             iob++)
                        {
                            const hier::Box& patch_overlapped_visible_box = *iob;
                            
                            if (patch_overlapped_visible_box.contains(idx_pt))
                            {
                                n_overlapped++;
                            }
                        }
                        
                        /*
                         * Compute the linear index and the data to add.
                         */
                        
                        const int idx = relative_idx_lo_0 + i + num_ghosts_0_quantity;
                        
                        const double value_to_add = u[idx]/((double) n_overlapped);
                        
                        for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                        {
                            const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                            
                            u_avg_local[idx_fine] += value_to_add;
                        }
                    }
                }
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        mpi.Allreduce(
            u_avg_local,
            u_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_avg_local);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        
        double* u_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_quantity.resize(finest_level_dim_0);
        double* u_avg_global = averaged_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_avg_local[i] = double(0);
            u_avg_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                boost::shared_ptr<pdat::CellData<double> > data_quantity(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch->getPatchData(variable_quantity, data_context)));
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                
                const double weight = dx[1]/L_y;
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            /*
                             * Compute the index of the data point and count how many times the data is repeated.
                             */
                            
                            const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j);
                            
                            int n_overlapped = 1;
                            
                            for (hier::BoxContainer::BoxContainerConstIterator iob(
                                    patch_overlapped_visible_boxes.begin());
                                 iob != patch_overlapped_visible_boxes.end();
                                 iob++)
                            {
                                const hier::Box& patch_overlapped_visible_box = *iob;
                                
                                if (patch_overlapped_visible_box.contains(idx_pt))
                                {
                                    n_overlapped++;
                                }
                            }
                            
                            /*
                             * Compute the linear index and the data to add.
                             */
                            
                            const int idx = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                            
                            const double value_to_add = u[idx]*weight/((double) n_overlapped);
                            
                            for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                            {
                                const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                
                                u_avg_local[idx_fine] += value_to_add;
                            }
                        }
                    }
                }
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        mpi.Allreduce(
            u_avg_local,
            u_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_avg_local);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        const double L_z = x_hi[2] - x_lo[2];
        
        double* u_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_quantity.resize(finest_level_dim_0);
        double* u_avg_global = averaged_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_avg_local[i] = double(0);
            u_avg_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                boost::shared_ptr<pdat::CellData<double> > data_quantity(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch->getPatchData(variable_quantity, data_context)));
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int num_ghosts_2_quantity = num_ghosts_quantity[2];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                const int ghostcell_dim_1_quantity = ghostcell_dims_quantity[1];
                
                const double weight = (dx[1]*dx[2])/(L_y*L_z);
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    const int interior_dim_2 = interior_dims[2];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int idx_lo_2 = index_lo[2];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    const int relative_idx_lo_2 = relative_index_lo[2];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                /*
                                 * Compute the index of the data point and count how many times the data is repeated.
                                 */
                                
                                const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j, idx_lo_2 + k);
                                
                                int n_overlapped = 1;
                                
                                for (hier::BoxContainer::BoxContainerConstIterator iob(
                                        patch_overlapped_visible_boxes.begin());
                                     iob != patch_overlapped_visible_boxes.end();
                                     iob++)
                                {
                                    const hier::Box& patch_overlapped_visible_box = *iob;
                                    
                                    if (patch_overlapped_visible_box.contains(idx_pt))
                                    {
                                        n_overlapped++;
                                    }
                                }
                                
                                /*
                                 * Compute the linear index and the data to add.
                                 */
                                
                                const int idx = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                    (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                        ghostcell_dim_1_quantity;
                                
                                const double value_to_add = u[idx]*weight/((double) n_overlapped);
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    u_avg_local[idx_fine] += value_to_add;
                                }
                            }
                        }
                    }
                }
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        mpi.Allreduce(
            u_avg_local,
            u_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_avg_local);
    }
    
    return averaged_quantity;
}


/*
 * Compute averaged value (on product of variables) with only x direction as inhomogeneous direction.
 */
std::vector<double>
FlowModelStatisticsUtilitiesFourEqnConservative::getAveragedQuantityWithInhomogeneousXDirection(
    std::vector<boost::shared_ptr<pdat::CellVariable<double> > >& variable_quantities,
    const std::vector<int>& component_indices,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context) const
{
    int num_quantities = static_cast<int>(variable_quantities.size());
    
    TBOX_ASSERT(static_cast<int>(component_indices.size()) == num_quantities);
    
    std::vector<double> averaged_quantity;
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Get the refinement ratio from the finest level to the coarest level.
     */
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    hier::IntVector ratioFinestLevelToCoarestLevel =
        patch_hierarchy->getRatioToCoarserLevel(num_levels - 1);
    for (int li = num_levels - 2; li > 0 ; li--)
    {
        ratioFinestLevelToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(li);
    }
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    boost::shared_ptr<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *patch_hierarchy,
            0,
            num_levels - 1));
    
    /*
     * Get the number of cells of physical domain refined to the finest level.
     */
    
    const hier::BoxContainer& physical_domain = d_grid_geometry->getPhysicalDomain();
    const hier::Box& physical_domain_box = physical_domain.front();
    const hier::IntVector& physical_domain_dims = physical_domain_box.numberCells();
    const hier::IntVector finest_level_dims = physical_domain_dims*ratioFinestLevelToCoarestLevel;
    
    /*
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        double* avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_quantity.resize(finest_level_dim_0);
        double* avg_global = averaged_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            avg_local[i] = double(0);
            avg_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<boost::shared_ptr<pdat::CellData<double> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch->getPatchData(variable_quantities[qi], data_context));
                }
                
                std::vector<double*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                }
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        /*
                         * Compute the index of the data point and count how many times the data is repeated.
                         */
                        
                        const hier::Index idx_pt(tbox::Dimension(1), idx_lo_0 + i);
                        
                        int n_overlapped = 1;
                        
                        for (hier::BoxContainer::BoxContainerConstIterator iob(
                                patch_overlapped_visible_boxes.begin());
                             iob != patch_overlapped_visible_boxes.end();
                             iob++)
                        {
                            const hier::Box& patch_overlapped_visible_box = *iob;
                            
                            if (patch_overlapped_visible_box.contains(idx_pt))
                            {
                                n_overlapped++;
                            }
                        }
                        
                        /*
                         * Compute the linear indices and the data to add.
                         */
                        
                        for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                        {
                            const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                            
                            const int idx_q0 = relative_idx_lo_0 + i + num_ghosts_0_u_qi[0];
                            
                            double avg = u_qi[0][idx_q0];
                            
                            for (int qi = 1; qi < num_quantities; qi++)
                            {
                                const int idx_qi = relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi];
                                
                                avg *= u_qi[qi][idx_qi];
                            }
                            
                            avg_local[idx_fine] += (avg/((double) n_overlapped));
                        }
                    }
                }
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        mpi.Allreduce(
            avg_local,
            avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(avg_local);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        
        double* avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_quantity.resize(finest_level_dim_0);
        double* avg_global = averaged_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            avg_local[i] = double(0);
            avg_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<boost::shared_ptr<pdat::CellData<double> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch->getPatchData(variable_quantities[qi], data_context));
                }
                
                std::vector<double*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                }
                
                const double weight = dx[1]/L_y;
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            /*
                             * Compute the index of the data point and count how many times the data is repeated.
                             */
                            
                            const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j);
                            
                            int n_overlapped = 1;
                            
                            for (hier::BoxContainer::BoxContainerConstIterator iob(
                                    patch_overlapped_visible_boxes.begin());
                                 iob != patch_overlapped_visible_boxes.end();
                                 iob++)
                            {
                                const hier::Box& patch_overlapped_visible_box = *iob;
                                
                                if (patch_overlapped_visible_box.contains(idx_pt))
                                {
                                    n_overlapped++;
                                }
                            }
                            
                            /*
                             * Compute the linear indices and the data to add.
                             */
                            
                            for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                            {
                                const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                
                                const int idx_q0 = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[0]) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0];
                                
                                double avg = u_qi[0][idx_q0];
                                
                                for (int qi = 1; qi < num_quantities; qi++)
                                {
                                    const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    avg *= u_qi[qi][idx_qi];
                                }
                                
                                avg_local[idx_fine] += (avg*weight/((double) n_overlapped));
                            }
                        }
                    }
                }
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        mpi.Allreduce(
            avg_local,
            avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(avg_local);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        const double L_z = x_hi[2] - x_lo[2];
        
        double* avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_quantity.resize(finest_level_dim_0);
        double* avg_global = averaged_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            avg_local[i] = double(0);
            avg_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<boost::shared_ptr<pdat::CellData<double> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch->getPatchData(variable_quantities[qi], data_context));
                }
                
                std::vector<double*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> num_ghosts_2_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                std::vector<int> ghostcell_dim_1_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                num_ghosts_2_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                ghostcell_dim_1_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    num_ghosts_2_u_qi.push_back(num_ghosts_u_qi[2]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                    ghostcell_dim_1_u_qi.push_back(ghostcell_dims_u_qi[1]);
                }
                
                const double weight = dx[1]*dx[2]/(L_y*L_z);
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    const int interior_dim_2 = interior_dims[2];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int idx_lo_2 = index_lo[2];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    const int relative_idx_lo_2 = relative_index_lo[2];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                /*
                                 * Compute the index of the data point and count how many times the data is repeated.
                                 */
                                
                                const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j, idx_lo_2 + k);
                                
                                int n_overlapped = 1;
                                
                                for (hier::BoxContainer::BoxContainerConstIterator iob(
                                        patch_overlapped_visible_boxes.begin());
                                     iob != patch_overlapped_visible_boxes.end();
                                     iob++)
                                {
                                    const hier::Box& patch_overlapped_visible_box = *iob;
                                    
                                    if (patch_overlapped_visible_box.contains(idx_pt))
                                    {
                                        n_overlapped++;
                                    }
                                }
                                
                                /*
                                 * Compute the linear index and the data to add.
                                 */
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    const int idx_q0 = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[0]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0] +
                                        (relative_idx_lo_2 + k + num_ghosts_2_u_qi[0])*ghostcell_dim_0_u_qi[0]*
                                            ghostcell_dim_1_u_qi[0];
                                    
                                    double avg = u_qi[0][idx_q0];
                                    
                                    for (int qi = 1; qi < num_quantities; qi++)
                                    {
                                        const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        avg *= u_qi[qi][idx_qi];
                                    }
                                    
                                    avg_local[idx_fine] += (avg*weight/((double) n_overlapped));
                                }
                            }
                        }
                    }
                }
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        mpi.Allreduce(
            avg_local,
            avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(avg_local);
    }
    
    return averaged_quantity;
}


/*
 * Compute averaged derivative of value with only x direction as inhomogeneous direction.
 */
std::vector<double>
FlowModelStatisticsUtilitiesFourEqnConservative::getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
    boost::shared_ptr<pdat::CellVariable<double> >& variable_quantity,
    const int component_idx,
    const int derivative_direction,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context) const
{
    std::vector<double> averaged_derivative_quantity;
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Get the refinement ratio from the finest level to the coarest level.
     */
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    hier::IntVector ratioFinestLevelToCoarestLevel =
        patch_hierarchy->getRatioToCoarserLevel(num_levels - 1);
    for (int li = num_levels - 2; li > 0 ; li--)
    {
        ratioFinestLevelToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(li);
    }
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    boost::shared_ptr<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *patch_hierarchy,
            0,
            num_levels - 1));
    
    /*
     * Get the number of cells of physical domain refined to the finest level.
     */
    
    const hier::BoxContainer& physical_domain = d_grid_geometry->getPhysicalDomain();
    const hier::Box& physical_domain_box = physical_domain.front();
    const hier::IntVector& physical_domain_dims = physical_domain_box.numberCells();
    const hier::IntVector finest_level_dims = physical_domain_dims*ratioFinestLevelToCoarestLevel;
    
    /*
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        double* u_der_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_derivative_quantity.resize(finest_level_dim_0);
        double* u_der_avg_global = averaged_derivative_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_der_avg_local[i] = double(0);
            u_der_avg_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                boost::shared_ptr<pdat::CellData<double> > data_quantity(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch->getPatchData(variable_quantity, data_context)));
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        /*
                         * Compute the index of the data point and count how many times the data is repeated.
                         */
                        
                        const hier::Index idx_pt(tbox::Dimension(1), idx_lo_0 + i);
                        
                        int n_overlapped = 1;
                        
                        for (hier::BoxContainer::BoxContainerConstIterator iob(
                                patch_overlapped_visible_boxes.begin());
                             iob != patch_overlapped_visible_boxes.end();
                             iob++)
                        {
                            const hier::Box& patch_overlapped_visible_box = *iob;
                            
                            if (patch_overlapped_visible_box.contains(idx_pt))
                            {
                                n_overlapped++;
                            }
                        }
                        
                        /*
                         * Compute the derivative.
                         */
                        
                        double value_to_add = double(0);
                        
                        if (derivative_direction == 0)
                        {
                            // Compute linear indices.
                            const int idx_LLL = relative_idx_lo_0 + (i - 3) + num_ghosts_0_quantity;
                            const int idx_LL  = relative_idx_lo_0 + (i - 2) + num_ghosts_0_quantity;
                            const int idx_L   = relative_idx_lo_0 + (i - 1) + num_ghosts_0_quantity;
                            const int idx_R   = relative_idx_lo_0 + (i + 1) + num_ghosts_0_quantity;
                            const int idx_RR  = relative_idx_lo_0 + (i + 2) + num_ghosts_0_quantity;
                            const int idx_RRR = relative_idx_lo_0 + (i + 3) + num_ghosts_0_quantity;
                            
                            const double dudx = (double(1)/double(60)*(u[idx_RRR] - u[idx_LLL])
                                - double(3)/double(20)*(u[idx_RR] - u[idx_LL])
                                + double(3)/double(4)*(u[idx_R] - u[idx_L]))/dx[0];
                            
                            value_to_add = dudx/((double) n_overlapped);
                        }
                        else
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "Cannot take derivative for one-dimensional problem!\n"
                                << "derivative_direction = " << derivative_direction << " given!\n"
                                << std::endl);
                        }
                        
                        /*
                         * Add the data.
                         */
                        
                        for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                        {
                            const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                            
                            u_der_avg_local[idx_fine] += value_to_add;
                        }
                    }
                }
            }
        }
        
        /*
         * Reduction to get the global average of derivative.
         */
        
        mpi.Allreduce(
            u_der_avg_local,
            u_der_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_der_avg_local);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        
        double* u_der_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_derivative_quantity.resize(finest_level_dim_0);
        double* u_der_avg_global = averaged_derivative_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_der_avg_local[i] = double(0);
            u_der_avg_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                boost::shared_ptr<pdat::CellData<double> > data_quantity(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch->getPatchData(variable_quantity, data_context)));
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                
                const double weight = dx[1]/L_y;
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            /*
                             * Compute the index of the data point and count how many times the data is repeated.
                             */
                            
                            const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j);
                            
                            int n_overlapped = 1;
                            
                            for (hier::BoxContainer::BoxContainerConstIterator iob(
                                    patch_overlapped_visible_boxes.begin());
                                 iob != patch_overlapped_visible_boxes.end();
                                 iob++)
                            {
                                const hier::Box& patch_overlapped_visible_box = *iob;
                                
                                if (patch_overlapped_visible_box.contains(idx_pt))
                                {
                                    n_overlapped++;
                                }
                            }
                            
                            /*
                             * Compute the derivative.
                             */
                            
                            double value_to_add = double(0);
                            
                            if (derivative_direction == 0)
                            {
                                const int idx_LLL = (relative_idx_lo_0 + (i - 3) + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const int idx_LL  = (relative_idx_lo_0 + (i - 2) + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const int idx_L   = (relative_idx_lo_0 + (i - 1) + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const int idx_R   = (relative_idx_lo_0 + (i + 1) + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const int idx_RR  = (relative_idx_lo_0 + (i + 2) + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const int idx_RRR = (relative_idx_lo_0 + (i + 3) + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const double dudx = (double(1)/double(60)*(u[idx_RRR] - u[idx_LLL])
                                    - double(3)/double(20)*(u[idx_RR] - u[idx_LL])
                                    + double(3)/double(4)*(u[idx_R] - u[idx_L]))/dx[0];
                                
                                value_to_add = dudx*weight/((double) n_overlapped);
                            }
                            else if (derivative_direction == 1)
                            {
                                const int idx_BBB = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + (j - 3) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const int idx_BB  = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + (j - 2) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const int idx_B   = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + (j - 1) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const int idx_T   = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + (j + 1) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const int idx_TT  = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + (j + 2) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const int idx_TTT = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + (j + 3) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const double dudy = (double(1)/double(60)*(u[idx_TTT] - u[idx_BBB])
                                    - double(3)/double(20)*(u[idx_TT] - u[idx_BB])
                                    + double(3)/double(4)*(u[idx_T] - u[idx_B]))/dx[1];
                                
                                value_to_add = dudy*weight/((double) n_overlapped);
                            }
                            else
                            {
                                TBOX_ERROR(d_object_name
                                    << ": "
                                    << "Cannot take derivative for two-dimensional problem!\n"
                                    << "derivative_direction = " << derivative_direction << " given!\n"
                                    << std::endl);
                            }
                            
                            /*
                             * Add the data.
                             */
                            
                            for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                            {
                                const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                
                                u_der_avg_local[idx_fine] += value_to_add;
                            }
                        }
                    }
                }
            }
        }
        
        /*
         * Reduction to get the global average of derivative.
         */
        
        mpi.Allreduce(
            u_der_avg_local,
            u_der_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_der_avg_local);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        const double L_z = x_hi[2] - x_lo[2];
        
        double* u_der_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_derivative_quantity.resize(finest_level_dim_0);
        double* u_der_avg_global = averaged_derivative_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_der_avg_local[i] = double(0);
            u_der_avg_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                boost::shared_ptr<pdat::CellData<double> > data_quantity(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch->getPatchData(variable_quantity, data_context)));
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int num_ghosts_2_quantity = num_ghosts_quantity[2];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                const int ghostcell_dim_1_quantity = ghostcell_dims_quantity[1];
                
                const double weight = (dx[1]*dx[2])/(L_y*L_z);
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    const int interior_dim_2 = interior_dims[2];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int idx_lo_2 = index_lo[2];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    const int relative_idx_lo_2 = relative_index_lo[2];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                /*
                                 * Compute the index of the data point and count how many times the data is repeated.
                                 */
                                
                                const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j, idx_lo_2 + k);
                                
                                int n_overlapped = 1;
                                
                                for (hier::BoxContainer::BoxContainerConstIterator iob(
                                        patch_overlapped_visible_boxes.begin());
                                     iob != patch_overlapped_visible_boxes.end();
                                     iob++)
                                {
                                    const hier::Box& patch_overlapped_visible_box = *iob;
                                    
                                    if (patch_overlapped_visible_box.contains(idx_pt))
                                    {
                                        n_overlapped++;
                                    }
                                }
                                
                                /*
                                 * Compute the derivative.
                                 */
                                
                                double value_to_add = double(0);
                                
                                if (derivative_direction == 0)
                                {
                                    const int idx_LLL = (relative_idx_lo_0 + (i - 3) + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_LL  = (relative_idx_lo_0 + (i - 2) + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_L   = (relative_idx_lo_0 + (i - 1) + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_R   = (relative_idx_lo_0 + (i + 1) + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_RR  = (relative_idx_lo_0 + (i + 2) + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_RRR = (relative_idx_lo_0 + (i + 3) + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const double dudx = (double(1)/double(60)*(u[idx_RRR] - u[idx_LLL])
                                        - double(3)/double(20)*(u[idx_RR] - u[idx_LL])
                                        + double(3)/double(4)*(u[idx_R] - u[idx_L]))/dx[0];
                                    
                                    value_to_add = dudx*weight/((double) n_overlapped);
                                }
                                else if (derivative_direction == 1)
                                {
                                    const int idx_BBB = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + (j - 3) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_BB  = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + (j - 2) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_B   = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + (j - 1) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_T   = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + (j + 1) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_TT  = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + (j + 2) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_TTT = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + (j + 3) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const double dudy = (double(1)/double(60)*(u[idx_TTT] - u[idx_BBB])
                                        - double(3)/double(20)*(u[idx_TT] - u[idx_BB])
                                        + double(3)/double(4)*(u[idx_T] - u[idx_B]))/dx[1];
                                    
                                    value_to_add = dudy*weight/((double) n_overlapped);
                                }
                                else if (derivative_direction == 2)
                                {
                                    const int idx_BBB = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + (k - 3) + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_BB  = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + (k - 2) + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_B   = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + (k - 1) + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_F   = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + (k + 1) + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_FF  = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + (k + 2) + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_FFF = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + (k + 3) + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const double dudz = (double(1)/double(60)*(u[idx_FFF] - u[idx_BBB])
                                        - double(3)/double(20)*(u[idx_FF] - u[idx_BB])
                                        + double(3)/double(4)*(u[idx_F] - u[idx_B]))/dx[2];
                                    
                                    value_to_add = dudz*weight/((double) n_overlapped);
                                }
                                else
                                {
                                    TBOX_ERROR(d_object_name
                                        << ": "
                                        << "Cannot take derivative for three-dimensional problem!\n"
                                        << "derivative_direction = " << derivative_direction << " given!\n"
                                        << std::endl);
                                }
                                
                                /*
                                 * Add the data.
                                 */
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    u_der_avg_local[idx_fine] += value_to_add;
                                }
                            }
                        }
                    }
                }
            }
        }
        
        /*
         * Reduction to get the global average of derivative.
         */
        
        mpi.Allreduce(
            u_der_avg_local,
            u_der_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_der_avg_local);
    }
    
    return averaged_derivative_quantity;
}


/*
 * Compute averaged derivative of value (on product of variables) with only x direction as inhomogeneous direction.
 */
std::vector<double>
FlowModelStatisticsUtilitiesFourEqnConservative::getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
    std::vector<boost::shared_ptr<pdat::CellVariable<double> > >& variable_quantities,
    const std::vector<int>& component_indices,
    const int derivative_direction,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context) const
{
    int num_quantities = static_cast<int>(variable_quantities.size());
    
    TBOX_ASSERT(static_cast<int>(component_indices.size()) == num_quantities);
    
    std::vector<double> averaged_derivative_quantity;
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Get the refinement ratio from the finest level to the coarest level.
     */
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    hier::IntVector ratioFinestLevelToCoarestLevel =
        patch_hierarchy->getRatioToCoarserLevel(num_levels - 1);
    for (int li = num_levels - 2; li > 0 ; li--)
    {
        ratioFinestLevelToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(li);
    }
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    boost::shared_ptr<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *patch_hierarchy,
            0,
            num_levels - 1));
    
    /*
     * Get the number of cells of physical domain refined to the finest level.
     */
    
    const hier::BoxContainer& physical_domain = d_grid_geometry->getPhysicalDomain();
    const hier::Box& physical_domain_box = physical_domain.front();
    const hier::IntVector& physical_domain_dims = physical_domain_box.numberCells();
    const hier::IntVector finest_level_dims = physical_domain_dims*ratioFinestLevelToCoarestLevel;
    
    /*
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        double* u_der_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_derivative_quantity.resize(finest_level_dim_0);
        double* u_der_avg_global = averaged_derivative_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_der_avg_local[i] = double(0);
            u_der_avg_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<boost::shared_ptr<pdat::CellData<double> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch->getPatchData(variable_quantities[qi], data_context));
                }
                
                std::vector<double*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                }
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        /*
                         * Compute the index of the data point and count how many times the data is repeated.
                         */
                        
                        const hier::Index idx_pt(tbox::Dimension(1), idx_lo_0 + i);
                        
                        int n_overlapped = 1;
                        
                        for (hier::BoxContainer::BoxContainerConstIterator iob(
                                patch_overlapped_visible_boxes.begin());
                             iob != patch_overlapped_visible_boxes.end();
                             iob++)
                        {
                            const hier::Box& patch_overlapped_visible_box = *iob;
                            
                            if (patch_overlapped_visible_box.contains(idx_pt))
                            {
                                n_overlapped++;
                            }
                        }
                        
                        /*
                         * Compute the derivative.
                         */
                        
                        double value_to_add = double(0);
                        
                        if (derivative_direction == 0)
                        {
                            double u_LLL = double(1);
                            double u_LL  = double(1);
                            double u_L   = double(1);
                            double u_R   = double(1);
                            double u_RR  = double(1);
                            double u_RRR = double(1);
                            
                            for (int qi = 0; qi < num_quantities; qi++)
                            {
                                // Compute linear indices.
                                const int idx_LLL = relative_idx_lo_0 + (i - 3) + num_ghosts_0_u_qi[qi];
                                const int idx_LL  = relative_idx_lo_0 + (i - 2) + num_ghosts_0_u_qi[qi];
                                const int idx_L   = relative_idx_lo_0 + (i - 1) + num_ghosts_0_u_qi[qi];
                                const int idx_R   = relative_idx_lo_0 + (i + 1) + num_ghosts_0_u_qi[qi];
                                const int idx_RR  = relative_idx_lo_0 + (i + 2) + num_ghosts_0_u_qi[qi];
                                const int idx_RRR = relative_idx_lo_0 + (i + 3) + num_ghosts_0_u_qi[qi];
                                
                                u_LLL *= u_qi[qi][idx_LLL];
                                u_LL  *= u_qi[qi][idx_LL ];
                                u_L   *= u_qi[qi][idx_L  ];
                                u_R   *= u_qi[qi][idx_R  ];
                                u_RR  *= u_qi[qi][idx_RR ];
                                u_RRR *= u_qi[qi][idx_RRR];
                            }
                            
                            const double dudx = (double(1)/double(60)*(u_RRR - u_LLL)
                                - double(3)/double(20)*(u_RR - u_LL)
                                + double(3)/double(4)*(u_R - u_L))/dx[0];
                            
                            value_to_add = dudx/((double) n_overlapped);
                        }
                        else
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "Cannot take derivative for one-dimensional problem!\n"
                                << "derivative_direction = " << derivative_direction << " given!\n"
                                << std::endl);
                        }
                        
                        /*
                         * Add the data.
                         */
                        
                        for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                        {
                            const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                            
                            u_der_avg_local[idx_fine] += value_to_add;
                        }
                    }
                }
            }
        }
        
        /*
         * Reduction to get the global average of derivative.
         */
        
        mpi.Allreduce(
            u_der_avg_local,
            u_der_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_der_avg_local);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        
        double* u_der_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_derivative_quantity.resize(finest_level_dim_0);
        double* u_der_avg_global = averaged_derivative_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_der_avg_local[i] = double(0);
            u_der_avg_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<boost::shared_ptr<pdat::CellData<double> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch->getPatchData(variable_quantities[qi], data_context));
                }
                
                std::vector<double*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                }
                
                const double weight = dx[1]/L_y;
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            /*
                             * Compute the index of the data point and count how many times the data is repeated.
                             */
                            
                            const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j);
                            
                            int n_overlapped = 1;
                            
                            for (hier::BoxContainer::BoxContainerConstIterator iob(
                                    patch_overlapped_visible_boxes.begin());
                                 iob != patch_overlapped_visible_boxes.end();
                                 iob++)
                            {
                                const hier::Box& patch_overlapped_visible_box = *iob;
                                
                                if (patch_overlapped_visible_box.contains(idx_pt))
                                {
                                    n_overlapped++;
                                }
                            }
                            
                            /*
                             * Compute the derivative.
                             */
                            
                            double value_to_add = double(0);
                            
                            if (derivative_direction == 0)
                            {
                                double u_LLL = double(1);
                                double u_LL  = double(1);
                                double u_L   = double(1);
                                double u_R   = double(1);
                                double u_RR  = double(1);
                                double u_RRR = double(1);
                                
                                for (int qi = 0; qi < num_quantities; qi++)
                                {
                                    // Compute linear indices.
                                    const int idx_LLL = (relative_idx_lo_0 + (i - 3) + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    const int idx_LL  = (relative_idx_lo_0 + (i - 2) + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    const int idx_L   = (relative_idx_lo_0 + (i - 1) + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    const int idx_R   = (relative_idx_lo_0 + (i + 1) + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    const int idx_RR  = (relative_idx_lo_0 + (i + 2) + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    const int idx_RRR = (relative_idx_lo_0 + (i + 3) + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    u_LLL *= u_qi[qi][idx_LLL];
                                    u_LL  *= u_qi[qi][idx_LL ];
                                    u_L   *= u_qi[qi][idx_L  ];
                                    u_R   *= u_qi[qi][idx_R  ];
                                    u_RR  *= u_qi[qi][idx_RR ];
                                    u_RRR *= u_qi[qi][idx_RRR];
                                }
                                
                                const double dudx = (double(1)/double(60)*(u_RRR - u_LLL)
                                    - double(3)/double(20)*(u_RR - u_LL)
                                    + double(3)/double(4)*(u_R - u_L))/dx[0];
                                
                                value_to_add = dudx*weight/((double) n_overlapped);
                            }
                            else if (derivative_direction == 1)
                            {
                                double u_BBB = double(1);
                                double u_BB  = double(1);
                                double u_B   = double(1);
                                double u_T   = double(1);
                                double u_TT  = double(1);
                                double u_TTT = double(1);
                                
                                for (int qi = 0; qi < num_quantities; qi++)
                                {
                                    // Compute linear indices.
                                    const int idx_BBB = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + (j - 3) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    const int idx_BB  = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + (j - 2) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    const int idx_B   = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + (j - 1) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    const int idx_T   = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + (j + 1) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    const int idx_TT  = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + (j + 2) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    const int idx_TTT = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + (j + 3) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    u_BBB *= u_qi[qi][idx_BBB];
                                    u_BB  *= u_qi[qi][idx_BB ];
                                    u_B   *= u_qi[qi][idx_B  ];
                                    u_T   *= u_qi[qi][idx_T  ];
                                    u_TT  *= u_qi[qi][idx_TT ];
                                    u_TTT *= u_qi[qi][idx_TTT];
                                }
                                
                                const double dudy = (double(1)/double(60)*(u_TTT - u_BBB)
                                    - double(3)/double(20)*(u_TT - u_BB)
                                    + double(3)/double(4)*(u_T - u_B))/dx[1];
                                
                                value_to_add = dudy*weight/((double) n_overlapped);
                            }
                            else
                            {
                                TBOX_ERROR(d_object_name
                                    << ": "
                                    << "Cannot take derivative for two-dimensional problem!\n"
                                    << "derivative_direction = " << derivative_direction << " given!\n"
                                    << std::endl);
                            }
                            
                            /*
                             * Add the data.
                             */
                            
                            for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                            {
                                const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                
                                u_der_avg_local[idx_fine] += value_to_add;
                            }
                        }
                    }
                }
            }
        }
        
        /*
         * Reduction to get the global average of derivative.
         */
        
        mpi.Allreduce(
            u_der_avg_local,
            u_der_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_der_avg_local);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        const double L_z = x_hi[2] - x_lo[2];
        
        double* u_der_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_derivative_quantity.resize(finest_level_dim_0);
        double* u_der_avg_global = averaged_derivative_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_der_avg_local[i] = double(0);
            u_der_avg_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<boost::shared_ptr<pdat::CellData<double> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch->getPatchData(variable_quantities[qi], data_context));
                }
                
                std::vector<double*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> num_ghosts_2_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                std::vector<int> ghostcell_dim_1_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                num_ghosts_2_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                ghostcell_dim_1_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    num_ghosts_2_u_qi.push_back(num_ghosts_u_qi[2]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                    ghostcell_dim_1_u_qi.push_back(ghostcell_dims_u_qi[1]);
                }
                
                const double weight = (dx[1]*dx[2])/(L_y*L_z);
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    const int interior_dim_2 = interior_dims[2];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int idx_lo_2 = index_lo[2];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    const int relative_idx_lo_2 = relative_index_lo[2];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                /*
                                 * Compute the index of the data point and count how many times the data is repeated.
                                 */
                                
                                const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j, idx_lo_2 + k);
                                
                                int n_overlapped = 1;
                                
                                for (hier::BoxContainer::BoxContainerConstIterator iob(
                                        patch_overlapped_visible_boxes.begin());
                                     iob != patch_overlapped_visible_boxes.end();
                                     iob++)
                                {
                                    const hier::Box& patch_overlapped_visible_box = *iob;
                                    
                                    if (patch_overlapped_visible_box.contains(idx_pt))
                                    {
                                        n_overlapped++;
                                    }
                                }
                                
                                /*
                                 * Compute the derivative.
                                 */
                                
                                double value_to_add = double(0);
                                
                                if (derivative_direction == 0)
                                {
                                    double u_LLL = double(1);
                                    double u_LL  = double(1);
                                    double u_L   = double(1);
                                    double u_R   = double(1);
                                    double u_RR  = double(1);
                                    double u_RRR = double(1);
                                    
                                    for (int qi = 0; qi < num_quantities; qi++)
                                    {
                                        // Compute linear indices.
                                        const int idx_LLL = (relative_idx_lo_0 + (i - 3) + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_LL  = (relative_idx_lo_0 + (i - 2) + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_L   = (relative_idx_lo_0 + (i - 1) + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_R   = (relative_idx_lo_0 + (i + 1) + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_RR  = (relative_idx_lo_0 + (i + 2) + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_RRR = (relative_idx_lo_0 + (i + 3) + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        u_LLL *= u_qi[qi][idx_LLL];
                                        u_LL  *= u_qi[qi][idx_LL ];
                                        u_L   *= u_qi[qi][idx_L  ];
                                        u_R   *= u_qi[qi][idx_R  ];
                                        u_RR  *= u_qi[qi][idx_RR ];
                                        u_RRR *= u_qi[qi][idx_RRR];
                                    }
                                    
                                    const double dudx = (double(1)/double(60)*(u_RRR - u_LLL)
                                        - double(3)/double(20)*(u_RR - u_LL)
                                        + double(3)/double(4)*(u_R - u_L))/dx[0];
                                    
                                    value_to_add = dudx*weight/((double) n_overlapped);
                                }
                                else if (derivative_direction == 1)
                                {
                                    double u_BBB = double(1);
                                    double u_BB  = double(1);
                                    double u_B   = double(1);
                                    double u_T   = double(1);
                                    double u_TT  = double(1);
                                    double u_TTT = double(1);
                                    
                                    for (int qi = 0; qi < num_quantities; qi++)
                                    {
                                        // Compute linear indices.
                                        const int idx_BBB = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + (j - 3) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_BB  = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + (j - 2) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_B   = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + (j - 1) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_T   = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + (j + 1) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_TT  = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + (j + 2) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_TTT = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + (j + 3) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        u_BBB *= u_qi[qi][idx_BBB];
                                        u_BB  *= u_qi[qi][idx_BB ];
                                        u_B   *= u_qi[qi][idx_B  ];
                                        u_T   *= u_qi[qi][idx_T  ];
                                        u_TT  *= u_qi[qi][idx_TT ];
                                        u_TTT *= u_qi[qi][idx_TTT];
                                    }
                                    
                                    const double dudy = (double(1)/double(60)*(u_TTT - u_BBB)
                                        - double(3)/double(20)*(u_TT - u_BB)
                                        + double(3)/double(4)*(u_T - u_B))/dx[1];
                                    
                                    value_to_add = dudy*weight/((double) n_overlapped);
                                }
                                else if (derivative_direction == 2)
                                {
                                    double u_BBB = double(1);
                                    double u_BB  = double(1);
                                    double u_B   = double(1);
                                    double u_F   = double(1);
                                    double u_FF  = double(1);
                                    double u_FFF = double(1);
                                    
                                    for (int qi = 0; qi < num_quantities; qi++)
                                    {
                                        // Compute linear indices.
                                        const int idx_BBB = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + (k - 3) + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_BB  = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + (k - 2) + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_B   = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + (k - 1) + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_F   = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + (k + 1) + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_FF  = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + (k + 2) + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_FFF = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + (k + 3) + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        u_BBB *= u_qi[qi][idx_BBB];
                                        u_BB  *= u_qi[qi][idx_BB ];
                                        u_B   *= u_qi[qi][idx_B  ];
                                        u_F   *= u_qi[qi][idx_F  ];
                                        u_FF  *= u_qi[qi][idx_FF ];
                                        u_FFF *= u_qi[qi][idx_FFF];
                                    }
                                    
                                    const double dudz = (double(1)/double(60)*(u_FFF - u_BBB)
                                        - double(3)/double(20)*(u_FF - u_BB)
                                        + double(3)/double(4)*(u_F - u_B))/dx[2];
                                    
                                    value_to_add = dudz*weight/((double) n_overlapped);
                                }
                                else
                                {
                                    TBOX_ERROR(d_object_name
                                        << ": "
                                        << "Cannot take derivative for three-dimensional problem!\n"
                                        << "derivative_direction = " << derivative_direction << " given!\n"
                                        << std::endl);
                                }
                                
                                /*
                                 * Add the data.
                                 */
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    u_der_avg_local[idx_fine] += value_to_add;
                                }
                            }
                        }
                    }
                }
            }
        }
        
        /*
         * Reduction to get the global average of derivative.
         */
        
        mpi.Allreduce(
            u_der_avg_local,
            u_der_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_der_avg_local);
    }
    
    return averaged_derivative_quantity;
}


/*
 * Compute correlation with only x direction as inhomogeneous direction.
 */
std::vector<double>
FlowModelStatisticsUtilitiesFourEqnConservative::getQuantityCorrelationWithInhomogeneousXDirection(
    std::vector<boost::shared_ptr<pdat::CellVariable<double> > >& variable_quantities,
    const std::vector<int>& component_indices,
    const std::vector<std::vector<double> >& averaged_quantities,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context) const
{
    int num_quantities = static_cast<int>(variable_quantities.size());
    
    TBOX_ASSERT(static_cast<int>(component_indices.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(averaged_quantities.size()) == num_quantities);
    
    std::vector<double> correlation;
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Get the refinement ratio from the finest level to the coarest level.
     */
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    hier::IntVector ratioFinestLevelToCoarestLevel =
        patch_hierarchy->getRatioToCoarserLevel(num_levels - 1);
    for (int li = num_levels - 2; li > 0 ; li--)
    {
        ratioFinestLevelToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(li);
    }
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    boost::shared_ptr<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *patch_hierarchy,
            0,
            num_levels - 1));
    
    /*
     * Get the number of cells of physical domain refined to the finest level.
     */
    
    const hier::BoxContainer& physical_domain = d_grid_geometry->getPhysicalDomain();
    const hier::Box& physical_domain_box = physical_domain.front();
    const hier::IntVector& physical_domain_dims = physical_domain_box.numberCells();
    const hier::IntVector finest_level_dims = physical_domain_dims*ratioFinestLevelToCoarestLevel;
    
    /*
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << "FlowModelStatisticsUtilitiesFourEqnConservative::\n"
            << "getQuantityCorrelationWithInhomogeneousXDirection:\n"
            << "Not implemented for 1D problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        
        std::vector<const double*> u_qi_avg_global;
        u_qi_avg_global.reserve(num_quantities);
        for (int qi = 0; qi < num_quantities; qi++)
        {
            u_qi_avg_global.push_back(averaged_quantities[qi].data());
        }
        
        double* corr_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        correlation.resize(finest_level_dim_0);
        double* corr_global = correlation.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            corr_local[i] = double(0);
            corr_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                // Get the dimensions of box that covers the interior of patch.
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<boost::shared_ptr<pdat::CellData<double> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch->getPatchData(variable_quantities[qi], data_context));
                }
                
                std::vector<double*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                }
                
                const double weight = dx[1]/L_y;
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            /*
                             * Compute the index of the data point and count how many times the data is repeated.
                             */
                            
                            const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j);
                            
                            int n_overlapped = 1;
                            
                            for (hier::BoxContainer::BoxContainerConstIterator iob(
                                    patch_overlapped_visible_boxes.begin());
                                 iob != patch_overlapped_visible_boxes.end();
                                 iob++)
                            {
                                const hier::Box& patch_overlapped_visible_box = *iob;
                                
                                if (patch_overlapped_visible_box.contains(idx_pt))
                                {
                                    n_overlapped++;
                                }
                            }
                            
                            /*
                             * Compute the linear indices and the data to add.
                             */
                            
                            for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                            {
                                const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                
                                const int idx_q0 = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[0]) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0];
                                
                                double corr = u_qi[0][idx_q0] - u_qi_avg_global[0][idx_fine];
                                
                                for (int qi = 1; qi < num_quantities; qi++)
                                {
                                    const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    corr *= (u_qi[qi][idx_qi] - u_qi_avg_global[qi][idx_fine]);
                                }
                                
                                corr_local[idx_fine] += (corr*weight/((double) n_overlapped));
                            }
                            
                        }
                    }
                }
            }
        }
        
        mpi.Allreduce(
            corr_local,
            corr_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(corr_local);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        const double L_z = x_hi[2] - x_lo[2];
        
        std::vector<const double*> u_qi_avg_global;
        u_qi_avg_global.reserve(num_quantities);
        for (int qi = 0; qi < num_quantities; qi++)
        {
            u_qi_avg_global.push_back(averaged_quantities[qi].data());
        }
        
        double* corr_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        correlation.resize(finest_level_dim_0);
        double* corr_global = correlation.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            corr_local[i] = double(0);
            corr_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                // Get the dimensions of box that covers the interior of patch.
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<boost::shared_ptr<pdat::CellData<double> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch->getPatchData(variable_quantities[qi], data_context));
                }
                
                std::vector<double*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> num_ghosts_2_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                std::vector<int> ghostcell_dim_1_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                num_ghosts_2_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                ghostcell_dim_1_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    num_ghosts_2_u_qi.push_back(num_ghosts_u_qi[2]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                    ghostcell_dim_1_u_qi.push_back(ghostcell_dims_u_qi[1]);
                }
                
                const double weight = dx[1]*dx[2]/(L_y*L_z);
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    const int interior_dim_2 = interior_dims[2];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int idx_lo_2 = index_lo[2];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    const int relative_idx_lo_2 = relative_index_lo[2];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                /*
                                 * Compute the index of the data point and count how many times the data is repeated.
                                 */
                                
                                const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j, idx_lo_2 + k);
                                
                                int n_overlapped = 1;
                                
                                for (hier::BoxContainer::BoxContainerConstIterator iob(
                                        patch_overlapped_visible_boxes.begin());
                                     iob != patch_overlapped_visible_boxes.end();
                                     iob++)
                                {
                                    const hier::Box& patch_overlapped_visible_box = *iob;
                                    
                                    if (patch_overlapped_visible_box.contains(idx_pt))
                                    {
                                        n_overlapped++;
                                    }
                                }
                                
                                /*
                                 * Compute the linear index and the data to add.
                                 */
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    const int idx_q0 = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[0]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0] +
                                        (relative_idx_lo_2 + k + num_ghosts_2_u_qi[0])*ghostcell_dim_0_u_qi[0]*
                                            ghostcell_dim_1_u_qi[0];
                                    
                                    double corr = u_qi[0][idx_q0] - u_qi_avg_global[0][idx_fine];
                                    
                                    for (int qi = 1; qi < num_quantities; qi++)
                                    {
                                        const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        corr *= (u_qi[qi][idx_qi] - u_qi_avg_global[qi][idx_fine]);
                                    }
                                    
                                    corr_local[idx_fine] += (corr*weight/((double) n_overlapped));
                                }
                            }
                        }
                    }
                }
            }
        }
        
        mpi.Allreduce(
            corr_local,
            corr_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(corr_local);
    }
    
    return correlation;
}


/*
 * Compute correlation with only x direction as inhomogeneous direction.
 */
std::vector<double>
FlowModelStatisticsUtilitiesFourEqnConservative::getQuantityCorrelationWithInhomogeneousXDirection(
    std::vector<boost::shared_ptr<pdat::CellVariable<double> > >& variable_quantities,
    const std::vector<int>& component_indices,
    const std::vector<bool>& use_reciprocal,
    const std::vector<std::vector<double> >& averaged_quantities,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context) const
{
    int num_quantities = static_cast<int>(variable_quantities.size());
    
    TBOX_ASSERT(static_cast<int>(component_indices.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(use_reciprocal.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(averaged_quantities.size()) == num_quantities);
    
    std::vector<double> correlation;
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Get the refinement ratio from the finest level to the coarest level.
     */
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    hier::IntVector ratioFinestLevelToCoarestLevel =
        patch_hierarchy->getRatioToCoarserLevel(num_levels - 1);
    for (int li = num_levels - 2; li > 0 ; li--)
    {
        ratioFinestLevelToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(li);
    }
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    boost::shared_ptr<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *patch_hierarchy,
            0,
            num_levels - 1));
    
    /*
     * Get the number of cells of physical domain refined to the finest level.
     */
    
    const hier::BoxContainer& physical_domain = d_grid_geometry->getPhysicalDomain();
    const hier::Box& physical_domain_box = physical_domain.front();
    const hier::IntVector& physical_domain_dims = physical_domain_box.numberCells();
    const hier::IntVector finest_level_dims = physical_domain_dims*ratioFinestLevelToCoarestLevel;
    
    /*
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << "FlowModelStatisticsUtilitiesFourEqnConservative::\n"
            << "getQuantityCorrelationWithInhomogeneousXDirection:\n"
            << "Not implemented for 1D problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        
        std::vector<const double*> u_qi_avg_global;
        u_qi_avg_global.reserve(num_quantities);
        for (int qi = 0; qi < num_quantities; qi++)
        {
            u_qi_avg_global.push_back(averaged_quantities[qi].data());
        }
        
        double* corr_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        correlation.resize(finest_level_dim_0);
        double* corr_global = correlation.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            corr_local[i] = double(0);
            corr_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                // Get the dimensions of box that covers the interior of patch.
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<boost::shared_ptr<pdat::CellData<double> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch->getPatchData(variable_quantities[qi], data_context));
                }
                
                std::vector<double*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                }
                
                const double weight = dx[1]/L_y;
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            /*
                             * Compute the index of the data point and count how many times the data is repeated.
                             */
                            
                            const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j);
                            
                            int n_overlapped = 1;
                            
                            for (hier::BoxContainer::BoxContainerConstIterator iob(
                                    patch_overlapped_visible_boxes.begin());
                                 iob != patch_overlapped_visible_boxes.end();
                                 iob++)
                            {
                                const hier::Box& patch_overlapped_visible_box = *iob;
                                
                                if (patch_overlapped_visible_box.contains(idx_pt))
                                {
                                    n_overlapped++;
                                }
                            }
                            
                            /*
                             * Compute the linear indices and the data to add.
                             */
                            
                            for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                            {
                                const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                
                                const int idx_q0 = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[0]) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0];
                                
                                double corr = double(0);
                                if (use_reciprocal[0])
                                {
                                    corr = double(1)/(u_qi[0][idx_q0]) - u_qi_avg_global[0][idx_fine];
                                }
                                else
                                {
                                    corr = u_qi[0][idx_q0] - u_qi_avg_global[0][idx_fine];
                                }
                                
                                for (int qi = 1; qi < num_quantities; qi++)
                                {
                                    const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    if (use_reciprocal[qi])
                                    {
                                        corr *= (double(1)/(u_qi[qi][idx_qi]) - u_qi_avg_global[qi][idx_fine]);
                                    }
                                    else
                                    {
                                        corr *= (u_qi[qi][idx_qi] - u_qi_avg_global[qi][idx_fine]);
                                    }
                                }
                                
                                corr_local[idx_fine] += (corr*weight/((double) n_overlapped));
                            }
                            
                        }
                    }
                }
            }
        }
        
        mpi.Allreduce(
            corr_local,
            corr_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(corr_local);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        const double L_z = x_hi[2] - x_lo[2];
        
        std::vector<const double*> u_qi_avg_global;
        u_qi_avg_global.reserve(num_quantities);
        for (int qi = 0; qi < num_quantities; qi++)
        {
            u_qi_avg_global.push_back(averaged_quantities[qi].data());
        }
        
        double* corr_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        correlation.resize(finest_level_dim_0);
        double* corr_global = correlation.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            corr_local[i] = double(0);
            corr_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                // Get the dimensions of box that covers the interior of patch.
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<boost::shared_ptr<pdat::CellData<double> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch->getPatchData(variable_quantities[qi], data_context));
                }
                
                std::vector<double*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> num_ghosts_2_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                std::vector<int> ghostcell_dim_1_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                num_ghosts_2_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                ghostcell_dim_1_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    num_ghosts_2_u_qi.push_back(num_ghosts_u_qi[2]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                    ghostcell_dim_1_u_qi.push_back(ghostcell_dims_u_qi[1]);
                }
                
                const double weight = dx[1]*dx[2]/(L_y*L_z);
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    const int interior_dim_2 = interior_dims[2];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int idx_lo_2 = index_lo[2];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    const int relative_idx_lo_2 = relative_index_lo[2];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                /*
                                 * Compute the index of the data point and count how many times the data is repeated.
                                 */
                                
                                const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j, idx_lo_2 + k);
                                
                                int n_overlapped = 1;
                                
                                for (hier::BoxContainer::BoxContainerConstIterator iob(
                                        patch_overlapped_visible_boxes.begin());
                                     iob != patch_overlapped_visible_boxes.end();
                                     iob++)
                                {
                                    const hier::Box& patch_overlapped_visible_box = *iob;
                                    
                                    if (patch_overlapped_visible_box.contains(idx_pt))
                                    {
                                        n_overlapped++;
                                    }
                                }
                                
                                /*
                                 * Compute the linear index and the data to add.
                                 */
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    const int idx_q0 = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[0]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0] +
                                        (relative_idx_lo_2 + k + num_ghosts_2_u_qi[0])*ghostcell_dim_0_u_qi[0]*
                                            ghostcell_dim_1_u_qi[0];
                                    
                                    double corr = double(0);
                                    if (use_reciprocal[0])
                                    {
                                        corr = double(1)/(u_qi[0][idx_q0]) - u_qi_avg_global[0][idx_fine];
                                    }
                                    else
                                    {
                                        corr = u_qi[0][idx_q0] - u_qi_avg_global[0][idx_fine];
                                    }
                                    
                                    for (int qi = 1; qi < num_quantities; qi++)
                                    {
                                        const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        if (use_reciprocal[qi])
                                        {
                                            corr *= (double(1)/(u_qi[qi][idx_qi]) - u_qi_avg_global[qi][idx_fine]);
                                        }
                                        else
                                        {
                                            corr *= (u_qi[qi][idx_qi] - u_qi_avg_global[qi][idx_fine]);
                                        }
                                    }
                                    
                                    corr_local[idx_fine] += (corr*weight/((double) n_overlapped));
                                }
                            }
                        }
                    }
                }
            }
        }
        
        mpi.Allreduce(
            corr_local,
            corr_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(corr_local);
    }
    
    return correlation;
}


/*
 * Compute correlation with only x direction as inhomogeneous direction.
 */
std::vector<double>
FlowModelStatisticsUtilitiesFourEqnConservative::getQuantityCorrelationWithInhomogeneousXDirection(
    std::vector<boost::shared_ptr<pdat::CellVariable<double> > >& variable_quantities,
    const std::vector<int>& component_indices,
    const std::vector<bool>& use_reciprocal,
    const std::vector<int>& derivative_directions,
    const std::vector<std::vector<double> >& averaged_quantities,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context) const
{
    int num_quantities = static_cast<int>(variable_quantities.size());
    
    TBOX_ASSERT(static_cast<int>(component_indices.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(use_reciprocal.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(derivative_directions.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(averaged_quantities.size()) == num_quantities);
    
    std::vector<double> correlation;
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Get the refinement ratio from the finest level to the coarest level.
     */
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    hier::IntVector ratioFinestLevelToCoarestLevel =
        patch_hierarchy->getRatioToCoarserLevel(num_levels - 1);
    for (int li = num_levels - 2; li > 0 ; li--)
    {
        ratioFinestLevelToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(li);
    }
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    boost::shared_ptr<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *patch_hierarchy,
            0,
            num_levels - 1));
    
    /*
     * Get the number of cells of physical domain refined to the finest level.
     */
    
    const hier::BoxContainer& physical_domain = d_grid_geometry->getPhysicalDomain();
    const hier::Box& physical_domain_box = physical_domain.front();
    const hier::IntVector& physical_domain_dims = physical_domain_box.numberCells();
    const hier::IntVector finest_level_dims = physical_domain_dims*ratioFinestLevelToCoarestLevel;
    
    /*
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << "FlowModelStatisticsUtilitiesFourEqnConservative::\n"
            << "getQuantityCorrelationWithInhomogeneousXDirection:\n"
            << "Not implemented for 1D problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        
        std::vector<const double*> u_qi_avg_global;
        u_qi_avg_global.reserve(num_quantities);
        for (int qi = 0; qi < num_quantities; qi++)
        {
            u_qi_avg_global.push_back(averaged_quantities[qi].data());
        }
        
        double* corr_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        correlation.resize(finest_level_dim_0);
        double* corr_global = correlation.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            corr_local[i] = double(0);
            corr_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                // Get the dimensions of box that covers the interior of patch.
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<boost::shared_ptr<pdat::CellData<double> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch->getPatchData(variable_quantities[qi], data_context));
                }
                
                std::vector<double*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                }
                
                const double weight = dx[1]/L_y;
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            /*
                             * Compute the index of the data point and count how many times the data is repeated.
                             */
                            
                            const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j);
                            
                            int n_overlapped = 1;
                            
                            for (hier::BoxContainer::BoxContainerConstIterator iob(
                                    patch_overlapped_visible_boxes.begin());
                                 iob != patch_overlapped_visible_boxes.end();
                                 iob++)
                            {
                                const hier::Box& patch_overlapped_visible_box = *iob;
                                
                                if (patch_overlapped_visible_box.contains(idx_pt))
                                {
                                    n_overlapped++;
                                }
                            }
                            
                            /*
                             * Perform operations on quantities.
                             */
                            
                            double u_values[num_quantities];
                            
                            for (int qi = 0; qi < num_quantities; qi++)
                            {
                                if (derivative_directions[qi] == -1)
                                {
                                    const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    if (use_reciprocal[qi])
                                    {
                                        u_values[qi] = double(1)/(u_qi[qi][idx_qi]);
                                    }
                                    else
                                    {
                                        u_values[qi] = u_qi[qi][idx_qi];
                                    }
                                }
                                else if (derivative_directions[qi] == 0)
                                {
                                    const int idx_LLL = (relative_idx_lo_0 + (i - 3) + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    const int idx_LL  = (relative_idx_lo_0 + (i - 2) + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    const int idx_L   = (relative_idx_lo_0 + (i - 1) + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    const int idx_R   = (relative_idx_lo_0 + (i + 1) + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    const int idx_RR  = (relative_idx_lo_0 + (i + 2) + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    const int idx_RRR = (relative_idx_lo_0 + (i + 3) + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    if (use_reciprocal[qi])
                                    {
                                        u_values[qi] = (double(1)/double(60)*(double(1)/u_qi[qi][idx_RRR] - double(1)/u_qi[qi][idx_LLL])
                                            - double(3)/double(20)*(double(1)/u_qi[qi][idx_RR] - double(1)/u_qi[qi][idx_LL])
                                            + double(3)/double(4)*(double(1)/u_qi[qi][idx_R] - double(1)/u_qi[qi][idx_L]))/dx[0];
                                    }
                                    else
                                    {
                                        u_values[qi] = (double(1)/double(60)*(u_qi[qi][idx_RRR] - u_qi[qi][idx_LLL])
                                            - double(3)/double(20)*(u_qi[qi][idx_RR] - u_qi[qi][idx_LL])
                                            + double(3)/double(4)*(u_qi[qi][idx_R] - u_qi[qi][idx_L]))/dx[0];
                                    }
                                }
                                else if (derivative_directions[qi] == 1)
                                {
                                    const int idx_BBB = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + (j - 3) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    const int idx_BB  = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + (j - 2) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    const int idx_B   = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + (j - 1) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    const int idx_T   = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + (j + 1) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    const int idx_TT  = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + (j + 2) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    const int idx_TTT = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + (j + 3) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    if (use_reciprocal[qi])
                                    {
                                        u_values[qi] = (double(1)/double(60)*(double(1)/u_qi[qi][idx_TTT] - double(1)/u_qi[qi][idx_BBB])
                                            - double(3)/double(20)*(double(1)/u_qi[qi][idx_TT] - double(1)/u_qi[qi][idx_BB])
                                            + double(3)/double(4)*(double(1)/u_qi[qi][idx_T] - double(1)/u_qi[qi][idx_B]))/dx[1];
                                    }
                                    else
                                    {
                                        u_values[qi] = (double(1)/double(60)*(u_qi[qi][idx_TTT] - u_qi[qi][idx_BBB])
                                            - double(3)/double(20)*(u_qi[qi][idx_TT] - u_qi[qi][idx_BB])
                                            + double(3)/double(4)*(u_qi[qi][idx_T] - u_qi[qi][idx_B]))/dx[1];
                                    }
                                }
                                else
                                {
                                    TBOX_ERROR(d_object_name
                                        << ": "
                                        << "Cannot take derivative for two-dimensional problem!\n"
                                        << "derivative_direction = " << derivative_directions[qi] << " given!\n"
                                        << std::endl);
                                }
                            }
                            
                            /*
                             * Compute the correlation.
                             */
                            
                            for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                            {
                                const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                
                                double corr = double(1);
                                for (int qi = 0; qi < num_quantities; qi++)
                                {
                                    corr *= (u_values[qi] - u_qi_avg_global[qi][idx_fine]);
                                }
                                
                                corr_local[idx_fine] += (corr*weight/((double) n_overlapped));
                            }
                        }
                    }
                }
            }
        }
        
        mpi.Allreduce(
            corr_local,
            corr_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(corr_local);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        const double L_z = x_hi[2] - x_lo[2];
        
        std::vector<const double*> u_qi_avg_global;
        u_qi_avg_global.reserve(num_quantities);
        for (int qi = 0; qi < num_quantities; qi++)
        {
            u_qi_avg_global.push_back(averaged_quantities[qi].data());
        }
        
        double* corr_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        correlation.resize(finest_level_dim_0);
        double* corr_global = correlation.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            corr_local[i] = double(0);
            corr_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                // Get the dimensions of box that covers the interior of patch.
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<boost::shared_ptr<pdat::CellData<double> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch->getPatchData(variable_quantities[qi], data_context));
                }
                
                std::vector<double*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> num_ghosts_2_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                std::vector<int> ghostcell_dim_1_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                num_ghosts_2_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                ghostcell_dim_1_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    num_ghosts_2_u_qi.push_back(num_ghosts_u_qi[2]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                    ghostcell_dim_1_u_qi.push_back(ghostcell_dims_u_qi[1]);
                }
                
                const double weight = dx[1]*dx[2]/(L_y*L_z);
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    const int interior_dim_2 = interior_dims[2];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int idx_lo_2 = index_lo[2];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    const int relative_idx_lo_2 = relative_index_lo[2];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                /*
                                 * Compute the index of the data point and count how many times the data is repeated.
                                 */
                                
                                const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j, idx_lo_2 + k);
                                
                                int n_overlapped = 1;
                                
                                for (hier::BoxContainer::BoxContainerConstIterator iob(
                                        patch_overlapped_visible_boxes.begin());
                                     iob != patch_overlapped_visible_boxes.end();
                                     iob++)
                                {
                                    const hier::Box& patch_overlapped_visible_box = *iob;
                                    
                                    if (patch_overlapped_visible_box.contains(idx_pt))
                                    {
                                        n_overlapped++;
                                    }
                                }
                                
                                /*
                                 * Perform operations on quantities.
                                 */
                                
                                double u_values[num_quantities];
                                
                                for (int qi = 0; qi < num_quantities; qi++)
                                {
                                    if (derivative_directions[qi] == -1)
                                    {
                                        const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        if (use_reciprocal[qi])
                                        {
                                            u_values[qi] = double(1)/(u_qi[qi][idx_qi]);
                                        }
                                        else
                                        {
                                            u_values[qi] = u_qi[qi][idx_qi];
                                        }
                                    }
                                    else if (derivative_directions[qi] == 0)
                                    {
                                        const int idx_LLL = (relative_idx_lo_0 + (i - 3) + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_LL  = (relative_idx_lo_0 + (i - 2) + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_L   = (relative_idx_lo_0 + (i - 1) + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_R   = (relative_idx_lo_0 + (i + 1) + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_RR  = (relative_idx_lo_0 + (i + 2) + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_RRR = (relative_idx_lo_0 + (i + 3) + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        if (use_reciprocal[qi])
                                        {
                                            u_values[qi] = (double(1)/double(60)*(double(1)/u_qi[qi][idx_RRR] - double(1)/u_qi[qi][idx_LLL])
                                                - double(3)/double(20)*(double(1)/u_qi[qi][idx_RR] - double(1)/u_qi[qi][idx_LL])
                                                + double(3)/double(4)*(double(1)/u_qi[qi][idx_R] - double(1)/u_qi[qi][idx_L]))/dx[0];
                                        }
                                        else
                                        {
                                            u_values[qi] = (double(1)/double(60)*(u_qi[qi][idx_RRR] - u_qi[qi][idx_LLL])
                                                - double(3)/double(20)*(u_qi[qi][idx_RR] - u_qi[qi][idx_LL])
                                                + double(3)/double(4)*(u_qi[qi][idx_R] - u_qi[qi][idx_L]))/dx[0];
                                        }
                                    }
                                    else if (derivative_directions[qi] == 1)
                                    {
                                        const int idx_BBB = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + (j - 3) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_BB  = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + (j - 2) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_B   = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + (j - 1) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_T   = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + (j + 1) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_TT  = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + (j + 2) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_TTT = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + (j + 3) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        if (use_reciprocal[qi])
                                        {
                                            u_values[qi] = (double(1)/double(60)*(double(1)/u_qi[qi][idx_TTT] - double(1)/u_qi[qi][idx_BBB])
                                                - double(3)/double(20)*(double(1)/u_qi[qi][idx_TT] - double(1)/u_qi[qi][idx_BB])
                                                + double(3)/double(4)*(double(1)/u_qi[qi][idx_T] - double(1)/u_qi[qi][idx_B]))/dx[1];
                                        }
                                        else
                                        {
                                            u_values[qi] = (double(1)/double(60)*(u_qi[qi][idx_TTT] - u_qi[qi][idx_BBB])
                                                - double(3)/double(20)*(u_qi[qi][idx_TT] - u_qi[qi][idx_BB])
                                                + double(3)/double(4)*(u_qi[qi][idx_T] - u_qi[qi][idx_B]))/dx[1];
                                        }
                                    }
                                    else if (derivative_directions[qi] == 2)
                                    {
                                        const int idx_BBB = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + (k - 3) + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_BB  = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + (k - 2) + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_B   = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + (k - 1) + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_F   = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + (k + 1) + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_FF  = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + (k + 2) + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        const int idx_FFF = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + (k + 3) + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        if (use_reciprocal[qi])
                                        {
                                            u_values[qi] = (double(1)/double(60)*(double(1)/u_qi[qi][idx_FFF] - double(1)/u_qi[qi][idx_BBB])
                                                - double(3)/double(20)*(double(1)/u_qi[qi][idx_FF] - double(1)/u_qi[qi][idx_BB])
                                                + double(3)/double(4)*(double(1)/u_qi[qi][idx_F] - double(1)/u_qi[qi][idx_B]))/dx[2];
                                        }
                                        else
                                        {
                                            u_values[qi] = (double(1)/double(60)*(u_qi[qi][idx_FFF] - u_qi[qi][idx_BBB])
                                                - double(3)/double(20)*(u_qi[qi][idx_FF] - u_qi[qi][idx_BB])
                                                + double(3)/double(4)*(u_qi[qi][idx_F] - u_qi[qi][idx_B]))/dx[2];
                                        }
                                    }
                                    else
                                    {
                                        TBOX_ERROR(d_object_name
                                            << ": "
                                            << "Cannot take derivative for three-dimensional problem!\n"
                                            << "derivative_direction = " << derivative_directions[qi] << " given!\n"
                                            << std::endl);
                                    }
                                }
                                
                                /*
                                 * Compute the correlation.
                                 */
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    double corr = double(1);
                                    for (int qi = 0; qi < num_quantities; qi++)
                                    {
                                        corr *= (u_values[qi] - u_qi_avg_global[qi][idx_fine]);
                                    }
                                    
                                    corr_local[idx_fine] += (corr*weight/((double) n_overlapped));
                                }
                            }
                        }
                    }
                }
            }
        }
        
        mpi.Allreduce(
            corr_local,
            corr_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(corr_local);
    }
    
    return correlation;
}


/*
 * Output budget of turbulent mass flux in x-direction with inhomogeneous x-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::
outputBudgetFilteredTurbMassFluxXWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double output_time) const
{
    tbox::pout << "FlowModelStatisticsUtilitiesFourEqnConservative::"
        << "outputBudgetFilteredTurbMassFluxXWithInhomogeneousXDirection: start" << std::endl;
    
    const int finest_level_dim_0 = getRefinedDomainNumberOfPointsX(patch_hierarchy);
    
    const double dx = getRefinedDomainGridSpacingX(patch_hierarchy);
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    std::vector<std::vector<double> > averaged_quantities;
    std::vector<bool> use_reciprocal;
    std::vector<int> derivative_directions;
    
    std::vector<boost::shared_ptr<pdat::CellVariable<double> > > variable_quantities;
    
    /*
     * Compute rho_a1.
     */
    
    std::vector<double> rho_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        patch_hierarchy,
        data_context);
    
    std::vector<double> u_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "VELOCITY",
        0,
        patch_hierarchy,
        data_context);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(rho_mean);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(u_mean);
    
    std::vector<double> rho_p_u_p = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    averaged_quantities.clear();
    
    std::vector<double> a_1(rho_p_u_p);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        a_1[i] /= rho_mean[i];
    }
    
    /*
     * Compute u_tilde.
     */
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    
    std::vector<double> rho_u_mean = getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        patch_hierarchy,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    
    std::vector<double> u_tilde(rho_u_mean);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        u_tilde[i] /= rho_mean[i];
    }
    
    /*
     * Compute term II.
     */
    
    std::vector<double> drho_u_dx_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        "MOMENTUM",
        0,
        0,
        patch_hierarchy,
        data_context);
    
    std::vector<double> drho_dx_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        0,
        patch_hierarchy,
        data_context);
    
    std::vector<double> du_dx_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        "VELOCITY",
        0,
        0,
        patch_hierarchy,
        data_context);
    
    std::vector<double> da1_dx(finest_level_dim_0, double(0));
    std::vector<double> d_rho_u_tilde_a1_dx(finest_level_dim_0, double(0));
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        da1_dx[i] = -rho_p_u_p[i]/(rho_mean[i]*rho_mean[i])*drho_dx_mean[i] +
            double(1)/rho_mean[i]*(drho_u_dx_mean[i] - u_mean[i]*drho_dx_mean[i]) - du_dx_mean[i];
        
        d_rho_u_tilde_a1_dx[i] = rho_u_mean[i]*da1_dx[i] + a_1[i]*drho_u_dx_mean[i];
    }
    
    /*
     * Compute term II in moving frame of mixing layer.
     */
    
    std::vector<double> rho_a1_a1(a_1);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        rho_a1_a1[i] *= rho_p_u_p[i];
    }
    
    std::vector<double> d_rho_a1_a1_dx = computeDerivativeOfVector1D(
        rho_a1_a1,
        dx);
    
    /*
     * Compute term III(1).
     */
    
    std::vector<double> rho_inv_mean = getAveragedReciprocalOfQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        patch_hierarchy,
        data_context);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_reciprocal.push_back(false);
    averaged_quantities.push_back(rho_mean);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_reciprocal.push_back(true);
    averaged_quantities.push_back(rho_inv_mean);
    
    std::vector<double> b = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        use_reciprocal,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        b[i] = -b[i];
    }
    
    quantity_names.clear();
    component_indices.clear();
    use_reciprocal.clear();
    averaged_quantities.clear();
    
    std::vector<double> dp_dx_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        s_variable_pressure_filtered,
        0,
        0,
        patch_hierarchy,
        data_context);
    
    std::vector<double> b_dp_dx(dp_dx_mean);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        b_dp_dx[i] *= b[i];
    }
    
    /*
     * Compute term III(2).
     */
    
    std::vector<double> dtau_11_dx_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        s_variable_shear_stress_filtered,
        0,
        0,
        patch_hierarchy,
        data_context);
    
    std::vector<double> m_b_dtau_11_dx(dtau_11_dx_mean);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        m_b_dtau_11_dx[i] *= (-b[i]);
    }
    
    /*
     * Compute term III(3).
     */
    
    std::vector<double> dtau_SFS_11_dx_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        s_variable_SFS_stress,
        0,
        0,
        patch_hierarchy,
        data_context);
    
    std::vector<double> b_dtau_SFS_11_dx(dtau_SFS_11_dx_mean);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        b_dtau_SFS_11_dx[i] *= (b[i]);
    }
    
    /*
     * Compute term III(4).
     */
    
    // Compute R_11.
    
    std::vector<double> zeros(finest_level_dim_0, double(0));
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(zeros);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(u_tilde);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(u_tilde);
    
    std::vector<double> R_11 = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    averaged_quantities.clear();
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        R_11[i] /= rho_mean[i];
    }
    
    std::vector<double> m_R_11_drho_dx(drho_dx_mean);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        m_R_11_drho_dx[i] *= (-R_11[i]);
    }
    
    /*
     * Compute term IV(1).
     */
    
    std::vector<double> rho_da_1_sq_dx(finest_level_dim_0, double(0));
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        rho_da_1_sq_dx[i] = double(2)*rho_mean[i]*a_1[i]*da1_dx[i];
    }
    
    /*
     * Compute term IV(2).
     */
    
    std::vector<double> m_rho_a_1_du_dx(du_dx_mean);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        m_rho_a_1_du_dx[i] *= (-rho_mean[i]*a_1[i]);
    }
    
    /*
     * Compute term V.
     */
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    
    std::vector<double> u_u_mean = getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        patch_hierarchy,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();

    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(rho_mean);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(u_mean);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(u_mean);
    
    std::vector<double> rho_p_u_p_u_p = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    averaged_quantities.clear();
    
    quantity_names.push_back("MOMENTUM");
    component_indices.push_back(0);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    
    std::vector<double> drho_u_u_dx_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        0,
        patch_hierarchy,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    
    std::vector<double> du_u_dx_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        0,
        patch_hierarchy,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    
    std::vector<double> m_rho_drho_p_u_p_sq_over_rho_dx(finest_level_dim_0, double(0));
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        m_rho_drho_p_u_p_sq_over_rho_dx[i] = -rho_mean[i]*(
            -rho_p_u_p_u_p[i]/(rho_mean[i]*rho_mean[i])*drho_dx_mean[i] +
            double(1)/rho_mean[i]*(
                drho_u_u_dx_mean[i] - double(2)*rho_u_mean[i]*du_dx_mean[i] - double(2)*u_mean[i]*drho_u_dx_mean[i] +
                double(2)*u_mean[i]*u_mean[i]*drho_dx_mean[i] + double(4)*rho_mean[i]*u_mean[i]*du_dx_mean[i] -
                rho_mean[i]*du_u_dx_mean[i] - u_u_mean[i]*drho_dx_mean[i]
            )
        );
    }
    
    // Old implementation.
    // 
    // quantity_names.push_back("DENSITY");
    // component_indices.push_back(0);
    // averaged_quantities.push_back(rho_mean);
    // 
    // quantity_names.push_back("VELOCITY");
    // component_indices.push_back(0);
    // averaged_quantities.push_back(u_mean);
    // 
    // quantity_names.push_back("VELOCITY");
    // component_indices.push_back(0);
    // averaged_quantities.push_back(u_mean);
    // 
    // std::vector<double> rho_p_u_p_sq_over_rho = getQuantityCorrelationWithInhomogeneousXDirection(
    //     quantity_names,
    //     component_indices,
    //     averaged_quantities,
    //     patch_hierarchy,
    //     data_context);
    // 
    // quantity_names.clear();
    // component_indices.clear();
    // averaged_quantities.clear();
    // 
    // for (int i = 0; i < finest_level_dim_0; i++)
    // {
    //     rho_p_u_p_sq_over_rho[i] /= rho_mean[i];
    // }
    // 
    // std::vector<double> drho_p_u_p_sq_over_rho_dx = computeDerivativeOfVector1D(
    //     rho_p_u_p_sq_over_rho,
    //     dx);
    // 
    // std::vector<double> m_rho_drho_p_u_p_sq_over_rho_dx(drho_p_u_p_sq_over_rho_dx);
    // 
    // for (int i = 0; i < finest_level_dim_0; i++)
    // {
    //     m_rho_drho_p_u_p_sq_over_rho_dx[i] *= (-rho_mean[i]);
    // }
    
    /*
     * Compute term VI(1).
     */
    
    variable_quantities.push_back(s_variable_specific_volume_Favre_filtered);
    component_indices.push_back(0);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(-1);
    averaged_quantities.push_back(rho_inv_mean);
    
    variable_quantities.push_back(s_variable_pressure_filtered);
    component_indices.push_back(0);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(0);
    averaged_quantities.push_back(dp_dx_mean);
    
    std::vector<double> rho_rho_inv_p_dp_dx_p = getQuantityCorrelationWithInhomogeneousXDirection(
        variable_quantities,
        component_indices,
        use_reciprocal,
        derivative_directions,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    variable_quantities.clear();
    component_indices.clear();
    use_reciprocal.clear();
    derivative_directions.clear();
    averaged_quantities.clear();
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        rho_rho_inv_p_dp_dx_p[i] *= rho_mean[i];
    }
    
    /*
     * Compute term VI(2).
     */
    
    variable_quantities.push_back(s_variable_specific_volume_Favre_filtered);
    component_indices.push_back(0);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(-1);
    averaged_quantities.push_back(rho_inv_mean);
    
    variable_quantities.push_back(s_variable_shear_stress_filtered);
    component_indices.push_back(0);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(0);
    averaged_quantities.push_back(dtau_11_dx_mean);
    
    std::vector<double> rho_inv_p_dtau_11_dx_p =
        getQuantityCorrelationWithInhomogeneousXDirection(
            variable_quantities,
            component_indices,
            use_reciprocal,
            derivative_directions,
            averaged_quantities,
            patch_hierarchy,
            data_context) ;
    
    variable_quantities.clear();
    component_indices.clear();
    use_reciprocal.clear();
    derivative_directions.clear();
    averaged_quantities.clear();
    
    std::vector<double> dtau_12_dy_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        s_variable_shear_stress_filtered,
        1,
        1,
        patch_hierarchy,
        data_context);
    
    variable_quantities.push_back(s_variable_specific_volume_Favre_filtered);
    component_indices.push_back(0);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(-1);
    averaged_quantities.push_back(rho_inv_mean);
    
    variable_quantities.push_back(s_variable_shear_stress_filtered);
    component_indices.push_back(1);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(1);
    averaged_quantities.push_back(dtau_12_dy_mean);
    
    std::vector<double> rho_inv_p_dtau_12_dy_p =
        getQuantityCorrelationWithInhomogeneousXDirection(
            variable_quantities,
            component_indices,
            use_reciprocal,
            derivative_directions,
            averaged_quantities,
            patch_hierarchy,
            data_context) ;
    
    variable_quantities.clear();
    component_indices.clear();
    use_reciprocal.clear();
    derivative_directions.clear();
    averaged_quantities.clear();
    
    std::vector<double> rho_inv_p_dtau_13_dz_p;
    if (d_dim == tbox::Dimension(3))
    {
        std::vector<double> dtau_13_dz_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
            s_variable_shear_stress_filtered,
            2,
            2,
            patch_hierarchy,
            data_context);
        
        variable_quantities.push_back(s_variable_specific_volume_Favre_filtered);
        component_indices.push_back(0);
        use_reciprocal.push_back(false);
        derivative_directions.push_back(-1);
        averaged_quantities.push_back(rho_inv_mean);
        
        variable_quantities.push_back(s_variable_shear_stress_filtered);
        component_indices.push_back(2);
        use_reciprocal.push_back(false);
        derivative_directions.push_back(2);
        averaged_quantities.push_back(dtau_13_dz_mean);
        
        rho_inv_p_dtau_13_dz_p =
            getQuantityCorrelationWithInhomogeneousXDirection(
                variable_quantities,
                component_indices,
                use_reciprocal,
                derivative_directions,
                averaged_quantities,
                patch_hierarchy,
                data_context) ;
        
        variable_quantities.clear();
        component_indices.clear();
        use_reciprocal.clear();
        derivative_directions.clear();
        averaged_quantities.clear();
    }
    
    std::vector<double> m_rho_rho_inv_p_dtau_ij_dx_p(finest_level_dim_0, double(0));
    
    if (d_dim == tbox::Dimension(2))
    {
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            m_rho_rho_inv_p_dtau_ij_dx_p[i] = -rho_mean[i]*(rho_inv_p_dtau_11_dx_p[i] + rho_inv_p_dtau_12_dy_p[i]);
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            m_rho_rho_inv_p_dtau_ij_dx_p[i] = -rho_mean[i]*(rho_inv_p_dtau_11_dx_p[i] + rho_inv_p_dtau_12_dy_p[i] +
                rho_inv_p_dtau_13_dz_p[i]);
        }
    }
    
    /*
     * Compute term VI(3).
     */
    
    variable_quantities.push_back(s_variable_specific_volume_Favre_filtered);
    component_indices.push_back(0);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(-1);
    averaged_quantities.push_back(rho_inv_mean);
    
    variable_quantities.push_back(s_variable_SFS_stress);
    component_indices.push_back(0);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(0);
    averaged_quantities.push_back(dtau_SFS_11_dx_mean);
    
    std::vector<double> rho_inv_p_dtau_SFS_11_dx_p =
        getQuantityCorrelationWithInhomogeneousXDirection(
            variable_quantities,
            component_indices,
            use_reciprocal,
            derivative_directions,
            averaged_quantities,
            patch_hierarchy,
            data_context) ;
    
    variable_quantities.clear();
    component_indices.clear();
    use_reciprocal.clear();
    derivative_directions.clear();
    averaged_quantities.clear();
    
    std::vector<double> dtau_SFS_12_dy_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        s_variable_SFS_stress,
        1,
        1,
        patch_hierarchy,
        data_context);
    
    variable_quantities.push_back(s_variable_specific_volume_Favre_filtered);
    component_indices.push_back(0);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(-1);
    averaged_quantities.push_back(rho_inv_mean);
    
    variable_quantities.push_back(s_variable_SFS_stress);
    component_indices.push_back(1);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(1);
    averaged_quantities.push_back(dtau_SFS_12_dy_mean);
    
    std::vector<double> rho_inv_p_dtau_SFS_12_dy_p =
        getQuantityCorrelationWithInhomogeneousXDirection(
            variable_quantities,
            component_indices,
            use_reciprocal,
            derivative_directions,
            averaged_quantities,
            patch_hierarchy,
            data_context) ;
    
    variable_quantities.clear();
    component_indices.clear();
    use_reciprocal.clear();
    derivative_directions.clear();
    averaged_quantities.clear();
    
    std::vector<double> rho_inv_p_dtau_SFS_13_dz_p;
    if (d_dim == tbox::Dimension(3))
    {
        std::vector<double> dtau_SFS_13_dz_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
            s_variable_SFS_stress,
            2,
            2,
            patch_hierarchy,
            data_context);
        
        variable_quantities.push_back(s_variable_specific_volume_Favre_filtered);
        component_indices.push_back(0);
        use_reciprocal.push_back(false);
        derivative_directions.push_back(-1);
        averaged_quantities.push_back(rho_inv_mean);
        
        variable_quantities.push_back(s_variable_SFS_stress);
        component_indices.push_back(2);
        use_reciprocal.push_back(false);
        derivative_directions.push_back(2);
        averaged_quantities.push_back(dtau_SFS_13_dz_mean);
        
        rho_inv_p_dtau_SFS_13_dz_p =
            getQuantityCorrelationWithInhomogeneousXDirection(
                variable_quantities,
                component_indices,
                use_reciprocal,
                derivative_directions,
                averaged_quantities,
                patch_hierarchy,
                data_context) ;
        
        variable_quantities.clear();
        component_indices.clear();
        use_reciprocal.clear();
        derivative_directions.clear();
        averaged_quantities.clear();
    }
    
    std::vector<double> rho_rho_inv_p_dtau_SFS_ij_dx_p(finest_level_dim_0, double(0));
    
    if (d_dim == tbox::Dimension(2))
    {
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            rho_rho_inv_p_dtau_SFS_ij_dx_p[i] = rho_mean[i]*(rho_inv_p_dtau_SFS_11_dx_p[i] + rho_inv_p_dtau_SFS_12_dy_p[i]);
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            rho_rho_inv_p_dtau_SFS_ij_dx_p[i] = rho_mean[i]*(rho_inv_p_dtau_SFS_11_dx_p[i] + rho_inv_p_dtau_SFS_12_dy_p[i] +
                rho_inv_p_dtau_SFS_13_dz_p[i]);
        }
    }
    
    /*
     * Compute term VI(4).
     */
    
    std::vector<double> dv_dy_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        "VELOCITY",
        1,
        1,
        patch_hierarchy,
        data_context);
    
    std::vector<double> dw_dz_mean;
    if (d_dim == tbox::Dimension(3))
    {
        dw_dz_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
            "VELOCITY",
            2,
            2,
            patch_hierarchy,
            data_context);
    }
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(-1);
    averaged_quantities.push_back(u_mean);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(0);
    averaged_quantities.push_back(du_dx_mean);
    
    std::vector<double> epsilon_a1_1 = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        use_reciprocal,
        derivative_directions,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    use_reciprocal.clear();
    derivative_directions.clear();
    averaged_quantities.clear();
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(-1);
    averaged_quantities.push_back(u_mean);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(1);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(1);
    averaged_quantities.push_back(dv_dy_mean);
    
    std::vector<double> epsilon_a1_2 = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        use_reciprocal,
        derivative_directions,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    use_reciprocal.clear();
    derivative_directions.clear();
    averaged_quantities.clear();
    
    std::vector<double> epsilon_a1_3;
    if (d_dim == tbox::Dimension(3))
    {
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(0);
        use_reciprocal.push_back(false);
        derivative_directions.push_back(-1);
        averaged_quantities.push_back(u_mean);
        
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(2);
        use_reciprocal.push_back(false);
        derivative_directions.push_back(2);
        averaged_quantities.push_back(dw_dz_mean);
        
        epsilon_a1_3 = getQuantityCorrelationWithInhomogeneousXDirection(
            quantity_names,
            component_indices,
            use_reciprocal,
            derivative_directions,
            averaged_quantities,
            patch_hierarchy,
            data_context);
        
        quantity_names.clear();
        component_indices.clear();
        use_reciprocal.clear();
        derivative_directions.clear();
        averaged_quantities.clear();
    }
    
    std::vector<double> rho_epsilon_a1(finest_level_dim_0, double(0));
    
    if (d_dim == tbox::Dimension(2))
    {
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            rho_epsilon_a1[i] = -rho_mean[i]*(epsilon_a1_1[i] + epsilon_a1_2[i]);
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            rho_epsilon_a1[i] = -rho_mean[i]*(epsilon_a1_1[i] + epsilon_a1_2[i] + epsilon_a1_3[i]);
        }
    }
    
    /*
     * Output budget.
     */
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_output;
        f_output.open(stat_dump_filename.c_str(), std::ios_base::app | std::ios::out | std::ios::binary);
        
        if (!f_output.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output turbulent mass flux budget!"
                << std::endl);
        }
        
        f_output.write((char*)&output_time, sizeof(double));
        f_output.write((char*)&rho_p_u_p[0], sizeof(double)*rho_p_u_p.size());
        // Term II.
        f_output.write((char*)&d_rho_u_tilde_a1_dx[0], sizeof(double)*d_rho_u_tilde_a1_dx.size());
        
        // Term III(1).
        f_output.write((char*)&b_dp_dx[0], sizeof(double)*b_dp_dx.size());
        // Term III(2).
        f_output.write((char*)&m_b_dtau_11_dx[0], sizeof(double)*m_b_dtau_11_dx.size());
        // Term III(2).
        f_output.write((char*)&b_dtau_SFS_11_dx[0], sizeof(double)*b_dtau_SFS_11_dx.size());
        // Term III(4).
        f_output.write((char*)&m_R_11_drho_dx[0], sizeof(double)*m_R_11_drho_dx.size());
        
        // Term IV(1).
        f_output.write((char*)&rho_da_1_sq_dx[0], sizeof(double)*rho_da_1_sq_dx.size());
        // Term IV(2).
        f_output.write((char*)&m_rho_a_1_du_dx[0], sizeof(double)*m_rho_a_1_du_dx.size());
        
        // Term V.
        f_output.write((char*)&m_rho_drho_p_u_p_sq_over_rho_dx[0], sizeof(double)*m_rho_drho_p_u_p_sq_over_rho_dx.size());
        
        // Term VI(1).
        f_output.write((char*)&rho_rho_inv_p_dp_dx_p[0], sizeof(double)*rho_rho_inv_p_dp_dx_p.size());
        // Term VI(2).
        f_output.write((char*)&m_rho_rho_inv_p_dtau_ij_dx_p[0], sizeof(double)*m_rho_rho_inv_p_dtau_ij_dx_p.size());
        // Term VI(3).
        f_output.write((char*)&rho_rho_inv_p_dtau_SFS_ij_dx_p[0], sizeof(double)*rho_rho_inv_p_dtau_SFS_ij_dx_p.size());
        // Term VI(4).
        f_output.write((char*)&rho_epsilon_a1[0], sizeof(double)*rho_epsilon_a1.size());
        
        // Term II in moving frame of mixing layer.
        f_output.write((char*)&d_rho_a1_a1_dx[0], sizeof(double)*d_rho_a1_a1_dx.size());
        
        f_output.close();
    }
    
    tbox::pout << "FlowModelStatisticsUtilitiesFourEqnConservative::"
        << "outputBudgetFilteredTurbMassFluxXWithInhomogeneousXDirection: end" << std::endl;
}


/*
 * Output budget of density specific volume covariance with inhomogeneous x-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::
outputBudgetFilteredDensitySpecificVolumeCovarianceWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double output_time) const
{
    tbox::pout << "FlowModelStatisticsUtilitiesFourEqnConservative::"
        << "outputBudgetFilteredDensitySpecificVolumeCovarianceWithInhomogeneousXDirection: start" << std::endl;
    
    const int finest_level_dim_0 = getRefinedDomainNumberOfPointsX(patch_hierarchy);
    
    const double dx = getRefinedDomainGridSpacingX(patch_hierarchy);
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    std::vector<std::vector<double> > averaged_quantities;
    std::vector<bool> use_reciprocal;
    std::vector<int> derivative_directions;
    
    /*
     * Compute rho_mean.
     */
    
    std::vector<double> rho_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        patch_hierarchy,
        data_context);
    
    /*
     * Compute u_mean.
     */
    
    std::vector<double> u_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "VELOCITY",
        0,
        patch_hierarchy,
        data_context);
    
    /*
     * Compute b.
     */
    
    std::vector<double> rho_inv_mean = getAveragedReciprocalOfQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        patch_hierarchy,
        data_context);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_reciprocal.push_back(false);
    averaged_quantities.push_back(rho_mean);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_reciprocal.push_back(true);
    averaged_quantities.push_back(rho_inv_mean);
    
    std::vector<double> b = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        use_reciprocal,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        b[i] = -b[i];
    }
    
    quantity_names.clear();
    component_indices.clear();
    use_reciprocal.clear();
    averaged_quantities.clear();
    
    /*
     * Compute rho_b.
     */
    
    std::vector<double> rho_b(b);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        rho_b[i] *= rho_mean[i];
    }
    
    /*
     * Compute u_tilde.
     */
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    
    std::vector<double> rho_u_mean = getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        patch_hierarchy,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    
    std::vector<double> u_tilde(rho_u_mean);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        u_tilde[i] /= rho_mean[i];
    }
    
    /*
     * Compute term II.
     */
    
   std::vector<double> drho_u_dx_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        "MOMENTUM",
        0,
        0,
        patch_hierarchy,
        data_context);
    
    std::vector<double> drho_dx_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        0,
        patch_hierarchy,
        data_context);
    
    std::vector<double> drho_inv_dx_mean = getAveragedDerivativeOfReciprocalOfQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        0,
        patch_hierarchy,
        data_context);
    
    std::vector<double> db_dx(finest_level_dim_0, double(0));
    std::vector<double> d_rho_u_tilde_b_dx(finest_level_dim_0, double(0));
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        db_dx[i] = rho_inv_mean[i]*drho_dx_mean[i] + rho_mean[i]*drho_inv_dx_mean[i];
        
        d_rho_u_tilde_b_dx[i] = rho_u_mean[i]*db_dx[i] + b[i]*drho_u_dx_mean[i];
    }
    
    // Old implementation.
    // 
    // std::vector<double> rho_u_tilde_b(u_tilde);
    // 
    // for (int i = 0; i < finest_level_dim_0; i++)
    // {
    //     rho_u_tilde_b[i] *= rho_b[i];
    // }
    // 
    // std::vector<double> d_rho_u_tilde_b_dx = computeDerivativeOfVector1D(
    //     rho_u_tilde_b,
    //     dx);
    
    /*
     * Compute term II in moving frame of mixing layer.
     */
    
    // Compute a1.
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(rho_mean);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(u_mean);
    
    std::vector<double> rho_p_u_p = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    averaged_quantities.clear();
    
    std::vector<double> a_1(rho_p_u_p);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        a_1[i] /= rho_mean[i];
    }
    
    std::vector<double> rho_a1_b(a_1);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        rho_a1_b[i] *= rho_b[i];
    }
    
    std::vector<double> d_rho_a1_b_dx = computeDerivativeOfVector1D(
        rho_a1_b,
        dx);
    
    /*
     * Compute term III.
     */
    
    std::vector<double> production(drho_dx_mean);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        production[i] *= (-double(2)*(b[i] + double(1))*a_1[i]);
    }
    
    /*
     * Compute term IV.
     */
    
    // std::vector<double> db_dx = computeDerivativeOfVector1D(
    //     b,
    //     dx);
    
    std::vector<double> redistribution(db_dx);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        redistribution[i] *= (double(2)*rho_mean[i]*a_1[i]);
    }
    
    /*
     * Compute term V.
     */
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_reciprocal.push_back(false);
    averaged_quantities.push_back(rho_mean);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_reciprocal.push_back(true);
    averaged_quantities.push_back(rho_inv_mean);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    use_reciprocal.push_back(false);
    averaged_quantities.push_back(u_mean);
    
    std::vector<double> rho_p_rho_inv_p_u_p = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        use_reciprocal,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    use_reciprocal.clear();
    averaged_quantities.clear();
    
   std::vector<double> du_dx_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        "VELOCITY",
        0,
        0,
        patch_hierarchy,
        data_context);
    
    quantity_names.push_back("VELOCITY");
    use_reciprocal.push_back(false);
    component_indices.push_back(0);
    
    quantity_names.push_back("DENSITY");
    use_reciprocal.push_back(true);
    component_indices.push_back(0);
    
    std::vector<double> du_rho_inv_dx_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        use_reciprocal,
        0,
        patch_hierarchy,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    use_reciprocal.clear();
    
    quantity_names.push_back("VELOCITY");
    use_reciprocal.push_back(false);
    component_indices.push_back(0);
    
    quantity_names.push_back("DENSITY");
    use_reciprocal.push_back(true);
    component_indices.push_back(0);
    
    std::vector<double> u_rho_inv_mean = getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        use_reciprocal,
        patch_hierarchy,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    use_reciprocal.clear();
    
    std::vector<double> turb_transport(finest_level_dim_0, double(0));
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        turb_transport[i] = (rho_mean[i]*rho_mean[i])*(
            -rho_p_rho_inv_p_u_p[i]/(rho_mean[i]*rho_mean[i])*drho_dx_mean[i] +
            double(1)/rho_mean[i]*(
                -rho_mean[i]*du_rho_inv_dx_mean[i] - u_rho_inv_mean[i]*drho_dx_mean[i] -
                rho_u_mean[i]*drho_inv_dx_mean[i] - rho_inv_mean[i]*drho_u_dx_mean[i] +
                double(2)*rho_mean[i]*rho_inv_mean[i]*du_dx_mean[i] + double(2)*rho_inv_mean[i]*u_mean[i]*drho_dx_mean[i] +
                double(2)*rho_mean[i]*u_mean[i]*drho_inv_dx_mean[i]
            )
        );
    }
    
    // Old implementation.
    // 
    // quantity_names.push_back("DENSITY");
    // component_indices.push_back(0);
    // use_reciprocal.push_back(false);
    // averaged_quantities.push_back(rho_mean);
    // 
    // quantity_names.push_back("DENSITY");
    // component_indices.push_back(0);
    // use_reciprocal.push_back(true);
    // averaged_quantities.push_back(rho_inv_mean);
    // 
    // quantity_names.push_back("VELOCITY");
    // component_indices.push_back(0);
    // use_reciprocal.push_back(false);
    // averaged_quantities.push_back(u_mean);
    // 
    // std::vector<double> rho_p_rho_inv_p_u_p_over_rho = getQuantityCorrelationWithInhomogeneousXDirection(
    //     quantity_names,
    //     component_indices,
    //     use_reciprocal,
    //     averaged_quantities,
    //     patch_hierarchy,
    //     data_context);
    // 
    // quantity_names.clear();
    // component_indices.clear();
    // use_reciprocal.clear();
    // averaged_quantities.clear();
    // 
    // for (int i = 0; i < finest_level_dim_0; i++)
    // {
    //     rho_p_rho_inv_p_u_p_over_rho[i] /= rho_mean[i];
    // }
    // 
    // std::vector<double> drho_p_rho_inv_p_u_p_over_rho_dx = computeDerivativeOfVector1D(
    //     rho_p_rho_inv_p_u_p_over_rho,
    //     dx);
    // 
    // std::vector<double> turb_transport(drho_p_rho_inv_p_u_p_over_rho_dx);
    // 
    // for (int i = 0; i < finest_level_dim_0; i++)
    // {
    //     turb_transport[i] *= (rho_mean[i]*rho_mean[i]);
    // }
    
    /*
     * Compute term VI.
     */
    
    // std::vector<double> du_dx_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
    //     "VELOCITY",
    //     0,
    //     0,
    //     patch_hierarchy,
    //     data_context);
    
    std::vector<double> dv_dy_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        "VELOCITY",
        1,
        1,
        patch_hierarchy,
        data_context);
    
    std::vector<double> dw_dz_mean;
    if (d_dim == tbox::Dimension(3))
    {
        dw_dz_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
            "VELOCITY",
            2,
            2,
            patch_hierarchy,
            data_context);
    }
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_reciprocal.push_back(true);
    derivative_directions.push_back(-1);
    averaged_quantities.push_back(rho_inv_mean);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(0);
    averaged_quantities.push_back(du_dx_mean);
    
    std::vector<double> epsilon_b_1 = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        use_reciprocal,
        derivative_directions,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    use_reciprocal.clear();
    derivative_directions.clear();
    averaged_quantities.clear();
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_reciprocal.push_back(true);
    derivative_directions.push_back(-1);
    averaged_quantities.push_back(rho_inv_mean);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(1);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(1);
    averaged_quantities.push_back(dv_dy_mean);
    
    std::vector<double> epsilon_b_2 = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        use_reciprocal,
        derivative_directions,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    use_reciprocal.clear();
    derivative_directions.clear();
    averaged_quantities.clear();
    
    std::vector<double> epsilon_b_3;
    if (d_dim == tbox::Dimension(3))
    {
        quantity_names.push_back("DENSITY");
        component_indices.push_back(0);
        use_reciprocal.push_back(true);
        derivative_directions.push_back(-1);
        averaged_quantities.push_back(rho_inv_mean);
        
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(2);
        use_reciprocal.push_back(false);
        derivative_directions.push_back(2);
        averaged_quantities.push_back(dw_dz_mean);
        
        epsilon_b_3 = getQuantityCorrelationWithInhomogeneousXDirection(
            quantity_names,
            component_indices,
            use_reciprocal,
            derivative_directions,
            averaged_quantities,
            patch_hierarchy,
            data_context);
        
        quantity_names.clear();
        component_indices.clear();
        use_reciprocal.clear();
        derivative_directions.clear();
        averaged_quantities.clear();
    }
    
    std::vector<double> rho_epsilon_b(finest_level_dim_0, double(0));
    
    if (d_dim == tbox::Dimension(2))
    {
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            rho_epsilon_b[i] = double(2)*rho_mean[i]*rho_mean[i]*(epsilon_b_1[i] + epsilon_b_2[i]);
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            rho_epsilon_b[i] = double(2)*rho_mean[i]*rho_mean[i]*(epsilon_b_1[i] + epsilon_b_2[i] + epsilon_b_3[i]);
        }
    }
    
    /*
     * Output budget.
     */
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_output;
        f_output.open(stat_dump_filename.c_str(), std::ios_base::app | std::ios::out | std::ios::binary);
        
        if (!f_output.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output density specific volume covariance budget!"
                << std::endl);
        }
        
        f_output.write((char*)&output_time, sizeof(double));
        
        f_output.write((char*)&rho_b[0], sizeof(double)*rho_b.size());
        
        // Term II.
        f_output.write((char*)&d_rho_u_tilde_b_dx[0], sizeof(double)*d_rho_u_tilde_b_dx.size());
        
        // Term III.
        f_output.write((char*)&production[0], sizeof(double)*production.size());
        
        // Term IV.
        f_output.write((char*)&redistribution[0], sizeof(double)*redistribution.size());
        
        // Term V.
        f_output.write((char*)&turb_transport[0], sizeof(double)*turb_transport.size());
        
        // Term VI.
        f_output.write((char*)&rho_epsilon_b[0], sizeof(double)*rho_epsilon_b.size());
        
        // Term II in moving frame of mixing layer.
        f_output.write((char*)&d_rho_a1_b_dx[0], sizeof(double)*d_rho_a1_b_dx.size());
        
        // Mean density profile.
        f_output.write((char*)&rho_mean[0], sizeof(double)*rho_mean.size());
        
        f_output.close();
    }
    
    tbox::pout << "FlowModelStatisticsUtilitiesFourEqnConservative::"
        << "outputBudgetFilteredDensitySpecificVolumeCovarianceWithInhomogeneousXDirection: end" << std::endl;
}


/*
 * Output budget of Reynolds normal stress in x-direction with inhomogeneous x-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::
outputBudgetFilteredReynoldsNormalStressInXDirectionWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double output_time) const
{
    tbox::pout << "FlowModelStatisticsUtilitiesFourEqnConservative::"
        << "outputBudgetFilteredReynoldsNormalStressInXDirectionWithInhomogeneousXDirection: start" << std::endl;
    
    const int finest_level_dim_0 = getRefinedDomainNumberOfPointsX(patch_hierarchy);
    
    const double dx = getRefinedDomainGridSpacingX(patch_hierarchy);
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    std::vector<std::vector<double> > averaged_quantities;
    std::vector<bool> use_reciprocal;
    std::vector<int> derivative_directions;
    
    std::vector<boost::shared_ptr<pdat::CellVariable<double> > > variable_quantities;
    
    std::vector<double> zeros(finest_level_dim_0, double(0));
    
    /*
     * Compute rho_mean.
     */
    
    std::vector<double> rho_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        patch_hierarchy,
        data_context);
    
    /*
     * Compute u_mean.
     */
    
    std::vector<double> u_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "VELOCITY",
        0,
        patch_hierarchy,
        data_context);
    
    /*
     * Compute p_mean.
     */
    
    std::vector<double> p_mean = getAveragedQuantityWithInhomogeneousXDirection(
        s_variable_pressure_filtered,
        0,
        patch_hierarchy,
        data_context);
    
    /*
     * Compute u_tilde.
     */
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    
    std::vector<double> rho_u_mean = getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        patch_hierarchy,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    
    std::vector<double> u_tilde(rho_u_mean);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        u_tilde[i] /= rho_mean[i];
    }
    
    /*
     * Compute a1.
     */
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(rho_mean);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(u_mean);
    
    std::vector<double> rho_p_u_p = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    averaged_quantities.clear();
    
    std::vector<double> a1(rho_p_u_p);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        a1[i] /= rho_mean[i];
    }
    
    /*
     * Compute R11.
     */
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(zeros);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(u_tilde);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(u_tilde);
    
    std::vector<double> rho_R11 = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    averaged_quantities.clear();
    
    std::vector<double> R11(rho_R11);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        R11[i] /= rho_mean[i];
    }
    
    /*
     * Compute term II.
     */
    
    tbox::pout << "FlowModelStatisticsUtilitiesFourEqnConservative::"
        << "outputBudgetFilteredReynoldsNormalStressInXDirectionWithInhomogeneousXDirection: Compute term II" << std::endl;
    
    std::vector<double> drho_u_dx_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        "MOMENTUM",
        0,
        0,
        patch_hierarchy,
        data_context);
    
    std::vector<double> drho_dx_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        0,
        patch_hierarchy,
        data_context);
    
    variable_quantities.push_back(s_variable_velocity_Favre_filtered);
    component_indices.push_back(0);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(-1);
    averaged_quantities.push_back(zeros);
    
    variable_quantities.push_back(s_variable_momentum_filtered);
    component_indices.push_back(0);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(0);
    averaged_quantities.push_back(zeros);
    
    std::vector<double> u_drho_u_dx_mean = getQuantityCorrelationWithInhomogeneousXDirection(
        variable_quantities,
        component_indices,
        use_reciprocal,
        derivative_directions,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    variable_quantities.clear();
    component_indices.clear();
    use_reciprocal.clear();
    derivative_directions.clear();
    averaged_quantities.clear();
    
    
    variable_quantities.push_back(s_variable_velocity_Favre_filtered);
    component_indices.push_back(0);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(-1);
    averaged_quantities.push_back(zeros);
    
    variable_quantities.push_back(s_variable_velocity_Favre_filtered);
    component_indices.push_back(0);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(-1);
    averaged_quantities.push_back(zeros);
    
    variable_quantities.push_back(s_variable_density_filtered);
    component_indices.push_back(0);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(0);
    averaged_quantities.push_back(zeros);
    
    std::vector<double> u_sq_drho_dx_mean =  getQuantityCorrelationWithInhomogeneousXDirection(
        variable_quantities,
        component_indices,
        use_reciprocal,
        derivative_directions,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    variable_quantities.clear();
    component_indices.clear();
    use_reciprocal.clear();
    derivative_directions.clear();
    averaged_quantities.clear();
    
    std::vector<double> drho_u_u_dx_mean(finest_level_dim_0, double(0));
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        drho_u_u_dx_mean[i] = double(2)*u_drho_u_dx_mean[i] - u_sq_drho_dx_mean[i];
    }
    
    std::vector<double> d_rho_u_tilde_R11_dx(finest_level_dim_0, double(0));
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        const double dR_11_tilde_dx = -(rho_R11[i]/(rho_mean[i]*rho_mean[i]))*drho_dx_mean[i] + 
            double(1)/rho_mean[i]*(drho_u_u_dx_mean[i] - double(2)*u_tilde[i]*drho_u_dx_mean[i] + 
            u_tilde[i]*u_tilde[i]*drho_dx_mean[i]);
        
        d_rho_u_tilde_R11_dx[i] = rho_u_mean[i]*dR_11_tilde_dx + R11[i]*drho_u_dx_mean[i];
    }
    
    /*
     * Compute term II in moving frame of mixing layer.
     */
    
    tbox::pout << "FlowModelStatisticsUtilitiesFourEqnConservative::"
        << "outputBudgetFilteredReynoldsNormalStressInXDirectionWithInhomogeneousXDirection: Compute term II in MF" << std::endl;
    
    std::vector<double> rho_a1_R11(rho_R11);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        rho_a1_R11[i] *= a1[i];
    }
    
    std::vector<double> d_rho_a1_R11_dx = computeDerivativeOfVector1D(
        rho_a1_R11,
        dx);
    
    /*
     * Compute term III(1).
     */
    
    tbox::pout << "FlowModelStatisticsUtilitiesFourEqnConservative::"
        << "outputBudgetFilteredReynoldsNormalStressInXDirectionWithInhomogeneousXDirection: Compute term III(1)" << std::endl;
    
    std::vector<double> dp_dx_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        s_variable_pressure_filtered,
        0,
        0,
        patch_hierarchy,
        data_context);
    
    std::vector<double> two_a1_dp_dx(dp_dx_mean);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        two_a1_dp_dx[i] *= (double(2)*a1[i]);
    }
    
    /*
     * Compute term III(2).
     */
    
    tbox::pout << "FlowModelStatisticsUtilitiesFourEqnConservative::"
        << "outputBudgetFilteredReynoldsNormalStressInXDirectionWithInhomogeneousXDirection: Compute term III(2)" << std::endl;
    
    std::vector<double> dtau_11_dx_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        s_variable_shear_stress_filtered,
        0,
        0,
        patch_hierarchy,
        data_context);
    
    std::vector<double> m_2a1_dtau_11_dx(dtau_11_dx_mean);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        m_2a1_dtau_11_dx[i] *= (-double(2)*a1[i]);
    }
    
    /*
     * Compute term III(3).
     */
    
    tbox::pout << "FlowModelStatisticsUtilitiesFourEqnConservative::"
        << "outputBudgetFilteredReynoldsNormalStressInXDirectionWithInhomogeneousXDirection: Compute term III(3)" << std::endl;
    
    std::vector<double> dtau_SFS_11_dx_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        s_variable_SFS_stress,
        0,
        0,
        patch_hierarchy,
        data_context);
    
    std::vector<double> two_a1_dtau_SFS_11_dx(dtau_SFS_11_dx_mean);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        two_a1_dtau_SFS_11_dx[i] *= (double(2)*a1[i]);
    }
    
    /*
     * Compute term III(4).
     */
    
    tbox::pout << "FlowModelStatisticsUtilitiesFourEqnConservative::"
        << "outputBudgetFilteredReynoldsNormalStressInXDirectionWithInhomogeneousXDirection: Compute term III(4)" << std::endl;
    
    std::vector<double> du_tilde_dx(finest_level_dim_0, double(0));
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        du_tilde_dx[i] = drho_u_dx_mean[i]/rho_mean[i] - rho_u_mean[i]/(rho_mean[i]*rho_mean[i])*drho_dx_mean[i];
    }
    
    std::vector<double> m_2rho_R11_du_tilde_dx(du_tilde_dx);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        m_2rho_R11_du_tilde_dx[i] *= (-double(2)*rho_R11[i]);
    }
    
    /*
     * Compute term IV(1).
     */
    
    tbox::pout << "FlowModelStatisticsUtilitiesFourEqnConservative::"
        << "outputBudgetFilteredReynoldsNormalStressInXDirectionWithInhomogeneousXDirection: Compute term IV(1)" << std::endl;
    
    // quantity_names.push_back("MOMENTUM");
    variable_quantities.push_back(s_variable_momentum_filtered);
    component_indices.push_back(0);
    
    // quantity_names.push_back("VELOCITY");
    variable_quantities.push_back(s_variable_velocity_Favre_filtered);
    component_indices.push_back(0);
    
    // quantity_names.push_back("VELOCITY");
    variable_quantities.push_back(s_variable_velocity_Favre_filtered);
    component_indices.push_back(0);
    
    std::vector<double> drho_u_u_u_dx_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        variable_quantities,
        component_indices,
        0,
        patch_hierarchy,
        data_context);
    
    // quantity_names.clear();
    variable_quantities.clear();
    component_indices.clear();
    
    // quantity_names.push_back("MOMENTUM");
    variable_quantities.push_back(s_variable_momentum_filtered);
    component_indices.push_back(0);
    
    // quantity_names.push_back("VELOCITY");
    variable_quantities.push_back(s_variable_velocity_Favre_filtered);
    component_indices.push_back(0);
    
    std::vector<double> rho_u_u_mean = getAveragedQuantityWithInhomogeneousXDirection(
        variable_quantities,
        component_indices,
        patch_hierarchy,
        data_context);
    
    // quantity_names.clear();
    variable_quantities.clear();
    component_indices.clear();
    
    std::vector<double> m_drho_u_pp_u_pp_u_pp_dx(finest_level_dim_0, double(0));
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        m_drho_u_pp_u_pp_u_pp_dx[i] = -(drho_u_u_u_dx_mean[i] - double(2)*rho_u_u_mean[i]*du_tilde_dx[i] -
            double(2)*u_tilde[i]*drho_u_u_dx_mean[i] + u_tilde[i]*u_tilde[i]*drho_u_dx_mean[i] +
            double(2)*rho_u_mean[i]*u_tilde[i]*du_tilde_dx[i] - d_rho_u_tilde_R11_dx[i]
            );
    }
    
    /*
     * Compute term IV(2).
     */
    
    tbox::pout << "FlowModelStatisticsUtilitiesFourEqnConservative::"
        << "outputBudgetFilteredReynoldsNormalStressInXDirectionWithInhomogeneousXDirection: Compute term IV(2)" << std::endl;
    
    std::vector<double> du_dx_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        "VELOCITY",
        0,
        0,
        patch_hierarchy,
        data_context);
    
    // quantity_names.push_back("PRESSURE");
    variable_quantities.push_back(s_variable_pressure_filtered);
    component_indices.push_back(0);
    
    // quantity_names.push_back("VELOCITY");
    variable_quantities.push_back(s_variable_velocity_Favre_filtered);
    component_indices.push_back(0);
    
    std::vector<double> du_pressure_dx_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        variable_quantities,
        component_indices,
        0,
        patch_hierarchy,
        data_context);
    
    // quantity_names.clear();
    variable_quantities.clear();
    component_indices.clear();
    
    std::vector<double> m_2du_p_p_p_dx(finest_level_dim_0, double(0));
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        m_2du_p_p_p_dx[i] = -double(2)*(du_pressure_dx_mean[i] - u_mean[i]*dp_dx_mean[i] - p_mean[i]*du_dx_mean[i]);
    }
    
    /*
     * Compute term IV(3).
     */
    
    tbox::pout << "FlowModelStatisticsUtilitiesFourEqnConservative::"
        << "outputBudgetFilteredReynoldsNormalStressInXDirectionWithInhomogeneousXDirection: Compute term IV(3)" << std::endl;
    
    std::vector<double> tau11_mean = getAveragedQuantityWithInhomogeneousXDirection(
        s_variable_shear_stress_filtered,
        0,
        patch_hierarchy,
        data_context);
    
    variable_quantities.push_back(s_variable_shear_stress_filtered);
    component_indices.push_back(0);
    
    variable_quantities.push_back(s_variable_velocity_Favre_filtered);
    component_indices.push_back(0);
    
    std::vector<double> du_tau_11_dx_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        variable_quantities,
        component_indices,
        0,
        patch_hierarchy,
        data_context);
    
    variable_quantities.clear();
    component_indices.clear();
    
    std::vector<double> two_du_p_tau11_p_dx(finest_level_dim_0, double(0));
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        two_du_p_tau11_p_dx[i] = double(2)*(du_tau_11_dx_mean[i] - u_mean[i]*dtau_11_dx_mean[i] - tau11_mean[i]*du_dx_mean[i]);
    }
    
    /*
     * Compute term IV(4).
     */
    
    tbox::pout << "FlowModelStatisticsUtilitiesFourEqnConservative::"
        << "outputBudgetFilteredReynoldsNormalStressInXDirectionWithInhomogeneousXDirection: Compute term IV(4)" << std::endl;
    
    std::vector<double> tau_SFS_11_mean = getAveragedQuantityWithInhomogeneousXDirection(
        s_variable_SFS_stress,
        0,
        patch_hierarchy,
        data_context);
    
    variable_quantities.push_back(s_variable_SFS_stress);
    component_indices.push_back(0);
    
    variable_quantities.push_back(s_variable_velocity_Favre_filtered);
    component_indices.push_back(0);
    
    std::vector<double> du_tau_SFS_11_dx_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        variable_quantities,
        component_indices,
        0,
        patch_hierarchy,
        data_context);
    
    variable_quantities.clear();
    component_indices.clear();
    
    std::vector<double> m_two_du_p_tau_SFS_11_p_dx(finest_level_dim_0, double(0));
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        m_two_du_p_tau_SFS_11_p_dx[i] = -double(2)*(du_tau_SFS_11_dx_mean[i] - u_mean[i]*dtau_SFS_11_dx_mean[i] - tau_SFS_11_mean[i]*du_dx_mean[i]);
    }
    
    /*
     * Compute term V.
     */
    
    tbox::pout << "FlowModelStatisticsUtilitiesFourEqnConservative::"
        << "outputBudgetFilteredReynoldsNormalStressInXDirectionWithInhomogeneousXDirection: Compute term V" << std::endl;
    
    variable_quantities.push_back(s_variable_pressure_filtered);
    component_indices.push_back(0);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(-1);
    averaged_quantities.push_back(p_mean);
    
    variable_quantities.push_back(s_variable_velocity_Favre_filtered);
    component_indices.push_back(0);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(0);
    averaged_quantities.push_back(du_dx_mean);
    
    std::vector<double> two_p_p_du_dx_p = getQuantityCorrelationWithInhomogeneousXDirection(
        variable_quantities,
        component_indices,
        use_reciprocal,
        derivative_directions,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    variable_quantities.clear();
    component_indices.clear();
    use_reciprocal.clear();
    derivative_directions.clear();
    averaged_quantities.clear();
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        two_p_p_du_dx_p[i] *= double(2);
    }
    
    /*
     * Compute term VI(1).
     */
    
    tbox::pout << "FlowModelStatisticsUtilitiesFourEqnConservative::"
        << "outputBudgetFilteredReynoldsNormalStressInXDirectionWithInhomogeneousXDirection: Compute term VI(1)" << std::endl;
    
    variable_quantities.push_back(s_variable_velocity_Favre_filtered);
    component_indices.push_back(0);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(0);
    averaged_quantities.push_back(du_dx_mean);
    
    variable_quantities.push_back(s_variable_shear_stress_filtered);
    component_indices.push_back(0);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(-1);
    averaged_quantities.push_back(tau11_mean);
    
    std::vector<double> tau_11_p_du_dx_p =
        getQuantityCorrelationWithInhomogeneousXDirection(
            variable_quantities,
            component_indices,
            use_reciprocal,
            derivative_directions,
            averaged_quantities,
            patch_hierarchy,
            data_context) ;
    
    variable_quantities.clear();
    component_indices.clear();
    use_reciprocal.clear();
    derivative_directions.clear();
    averaged_quantities.clear();
    
    std::vector<double> du_dy_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        "VELOCITY",
        0,
        1,
        patch_hierarchy,
        data_context);
    
    std::vector<double> tau12_mean = getAveragedQuantityWithInhomogeneousXDirection(
        s_variable_shear_stress_filtered,
        1,
        patch_hierarchy,
        data_context);
    
    variable_quantities.push_back(s_variable_velocity_Favre_filtered);
    component_indices.push_back(0);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(1);
    averaged_quantities.push_back(du_dy_mean);
    
    variable_quantities.push_back(s_variable_shear_stress_filtered);
    component_indices.push_back(1);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(-1);
    averaged_quantities.push_back(tau12_mean);
    
    std::vector<double> tau_12_p_du_dy_p =
        getQuantityCorrelationWithInhomogeneousXDirection(
            variable_quantities,
            component_indices,
            use_reciprocal,
            derivative_directions,
            averaged_quantities,
            patch_hierarchy,
            data_context) ;
    
    variable_quantities.clear();
    component_indices.clear();
    use_reciprocal.clear();
    derivative_directions.clear();
    averaged_quantities.clear();
    
    std::vector<double> tau_13_p_du_dz_p;
    if (d_dim == tbox::Dimension(3))
    {
        std::vector<double> du_dz_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
            "VELOCITY",
            0,
            2,
            patch_hierarchy,
            data_context);
        
        std::vector<double> tau13_mean = getAveragedQuantityWithInhomogeneousXDirection(
            s_variable_shear_stress_filtered,
            2,
            patch_hierarchy,
            data_context);
        
        variable_quantities.push_back(s_variable_velocity_Favre_filtered);
        component_indices.push_back(0);
        use_reciprocal.push_back(false);
        derivative_directions.push_back(2);
        averaged_quantities.push_back(du_dz_mean);
        
        variable_quantities.push_back(s_variable_shear_stress_filtered);
        component_indices.push_back(2);
        use_reciprocal.push_back(false);
        derivative_directions.push_back(-1);
        averaged_quantities.push_back(tau13_mean);
        
        tau_13_p_du_dz_p =
            getQuantityCorrelationWithInhomogeneousXDirection(
                variable_quantities,
                component_indices,
                use_reciprocal,
                derivative_directions,
                averaged_quantities,
                patch_hierarchy,
                data_context);
        
        variable_quantities.clear();
        component_indices.clear();
        use_reciprocal.clear();
        derivative_directions.clear();
        averaged_quantities.clear();
    }
    
    std::vector<double> m_2tau_1i_p_du_p_dxi(finest_level_dim_0, double(0));
    
    if (d_dim == tbox::Dimension(2))
    {
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            m_2tau_1i_p_du_p_dxi[i] = -double(2)*(tau_11_p_du_dx_p[i] + tau_12_p_du_dy_p[i]);
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            m_2tau_1i_p_du_p_dxi[i] = -double(2)*(tau_11_p_du_dx_p[i] + tau_12_p_du_dy_p[i] + tau_13_p_du_dz_p[i]);
        }
    }
    
    /*
     * Compute term VI(2).
     */
    
    tbox::pout << "FlowModelStatisticsUtilitiesFourEqnConservative::"
        << "outputBudgetFilteredReynoldsNormalStressInXDirectionWithInhomogeneousXDirection: Compute term VI(2)" << std::endl;
    
    variable_quantities.push_back(s_variable_velocity_Favre_filtered);
    component_indices.push_back(0);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(0);
    averaged_quantities.push_back(du_dx_mean);
    
    variable_quantities.push_back(s_variable_SFS_stress);
    component_indices.push_back(0);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(-1);
    averaged_quantities.push_back(tau_SFS_11_mean);
    
    std::vector<double> tau_SFS_11_p_du_dx_p =
        getQuantityCorrelationWithInhomogeneousXDirection(
            variable_quantities,
            component_indices,
            use_reciprocal,
            derivative_directions,
            averaged_quantities,
            patch_hierarchy,
            data_context) ;
    
    variable_quantities.clear();
    component_indices.clear();
    use_reciprocal.clear();
    derivative_directions.clear();
    averaged_quantities.clear();
    
    std::vector<double> tau_SFS_12_mean = getAveragedQuantityWithInhomogeneousXDirection(
        s_variable_SFS_stress,
        1,
        patch_hierarchy,
        data_context);
    
    variable_quantities.push_back(s_variable_velocity_Favre_filtered);
    component_indices.push_back(0);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(1);
    averaged_quantities.push_back(du_dy_mean);
    
    variable_quantities.push_back(s_variable_SFS_stress);
    component_indices.push_back(1);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(-1);
    averaged_quantities.push_back(tau_SFS_12_mean);
    
    std::vector<double> tau_SFS_12_p_du_dy_p =
        getQuantityCorrelationWithInhomogeneousXDirection(
            variable_quantities,
            component_indices,
            use_reciprocal,
            derivative_directions,
            averaged_quantities,
            patch_hierarchy,
            data_context);
    
    variable_quantities.clear();
    component_indices.clear();
    use_reciprocal.clear();
    derivative_directions.clear();
    averaged_quantities.clear();
    
    std::vector<double> tau_SFS_13_p_du_dz_p;
    if (d_dim == tbox::Dimension(3))
    {
        std::vector<double> du_dz_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
            "VELOCITY",
            0,
            2,
            patch_hierarchy,
            data_context);
        
        std::vector<double> tau_SFS_13_mean = getAveragedQuantityWithInhomogeneousXDirection(
            s_variable_SFS_stress,
            2,
            patch_hierarchy,
            data_context);
        
        variable_quantities.push_back(s_variable_velocity_Favre_filtered);
        component_indices.push_back(0);
        use_reciprocal.push_back(false);
        derivative_directions.push_back(2);
        averaged_quantities.push_back(du_dz_mean);
        
        variable_quantities.push_back(s_variable_SFS_stress);
        component_indices.push_back(2);
        use_reciprocal.push_back(false);
        derivative_directions.push_back(-1);
        averaged_quantities.push_back(tau_SFS_13_mean);
        
        tau_SFS_13_p_du_dz_p =
            getQuantityCorrelationWithInhomogeneousXDirection(
                variable_quantities,
                component_indices,
                use_reciprocal,
                derivative_directions,
                averaged_quantities,
                patch_hierarchy,
                data_context);
        
        variable_quantities.clear();
        component_indices.clear();
        use_reciprocal.clear();
        derivative_directions.clear();
        averaged_quantities.clear();
    }
    
    std::vector<double> two_tau_SFS_1i_p_du_p_dxi(finest_level_dim_0, double(0));
    
    if (d_dim == tbox::Dimension(2))
    {
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            two_tau_SFS_1i_p_du_p_dxi[i] = double(2)*(tau_SFS_11_p_du_dx_p[i] + tau_SFS_12_p_du_dy_p[i]);
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            two_tau_SFS_1i_p_du_p_dxi[i] = double(2)*(tau_SFS_11_p_du_dx_p[i] + tau_SFS_12_p_du_dy_p[i] + tau_SFS_13_p_du_dz_p[i]);
        }
    }
    
    /*
     * Output budget.
     */
    
    tbox::pout << "FlowModelStatisticsUtilitiesFourEqnConservative::"
        << "outputBudgetFilteredReynoldsNormalStressInXDirectionWithInhomogeneousXDirection: Output budget" << std::endl;
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_output;
        f_output.open(stat_dump_filename.c_str(), std::ios_base::app | std::ios::out | std::ios::binary);
        
        if (!f_output.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output R11 budget!"
                << std::endl);
        }
        
        f_output.write((char*)&output_time, sizeof(double));
        f_output.write((char*)&rho_R11[0], sizeof(double)*rho_R11.size());
        // Term II.
        f_output.write((char*)&d_rho_u_tilde_R11_dx[0], sizeof(double)*d_rho_u_tilde_R11_dx.size());
        
        // Term III(1).
        f_output.write((char*)&two_a1_dp_dx[0], sizeof(double)*two_a1_dp_dx.size());
        // Term III(2).
        f_output.write((char*)&m_2a1_dtau_11_dx[0], sizeof(double)*m_2a1_dtau_11_dx.size());
        // Term III(3).
        f_output.write((char*)&two_a1_dtau_SFS_11_dx[0], sizeof(double)*two_a1_dtau_SFS_11_dx.size());
        // Term III(4).
        f_output.write((char*)&m_2rho_R11_du_tilde_dx[0], sizeof(double)*m_2rho_R11_du_tilde_dx.size());
        
        // Term IV(1).
        f_output.write((char*)&m_drho_u_pp_u_pp_u_pp_dx[0], sizeof(double)*m_drho_u_pp_u_pp_u_pp_dx.size());
        // Term IV(2).
        f_output.write((char*)&m_2du_p_p_p_dx[0], sizeof(double)*m_2du_p_p_p_dx.size());
        // Term IV(3).
        f_output.write((char*)&two_du_p_tau11_p_dx[0], sizeof(double)*two_du_p_tau11_p_dx.size());
        // Term IV(4).
        f_output.write((char*)&m_two_du_p_tau_SFS_11_p_dx[0], sizeof(double)*m_two_du_p_tau_SFS_11_p_dx.size());
        
        // Term V.
        f_output.write((char*)&two_p_p_du_dx_p[0], sizeof(double)*two_p_p_du_dx_p.size());
        
        // Term VI(1).
        f_output.write((char*)&m_2tau_1i_p_du_p_dxi[0], sizeof(double)*m_2tau_1i_p_du_p_dxi.size());
        // Term VI(2).
        f_output.write((char*)&two_tau_SFS_1i_p_du_p_dxi[0], sizeof(double)*two_tau_SFS_1i_p_du_p_dxi.size());
        
        // Term II in moving frame of mixing layer.
        f_output.write((char*)&d_rho_a1_R11_dx[0], sizeof(double)*d_rho_a1_R11_dx.size());
        
        // SFS stress profile.
        f_output.write((char*)&tau_SFS_11_mean[0], sizeof(double)*tau_SFS_11_mean.size());
        
        f_output.close();
    }
    
    tbox::pout << "FlowModelStatisticsUtilitiesFourEqnConservative::"
        << "outputBudgetFilteredReynoldsNormalStressInXDirectionWithInhomogeneousXDirection: end" << std::endl;
}


/*
 * Output budget of Reynolds normal stress in y-direction with inhomogeneous x-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::
outputBudgetFilteredReynoldsNormalStressInYDirectionWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double output_time) const
{
    tbox::pout << "FlowModelStatisticsUtilitiesFourEqnConservative::"
        << "outputBudgetFilteredReynoldsNormalStressInYDirectionWithInhomogeneousXDirection: start" << std::endl;
    
    const int finest_level_dim_0 = getRefinedDomainNumberOfPointsX(patch_hierarchy);
    
    const double dx = getRefinedDomainGridSpacingX(patch_hierarchy);
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    std::vector<std::vector<double> > averaged_quantities;
    std::vector<bool> use_reciprocal;
    std::vector<int> derivative_directions;
    
    std::vector<boost::shared_ptr<pdat::CellVariable<double> > > variable_quantities;
    
    std::vector<double> zeros(finest_level_dim_0, double(0));
    
    /*
     * Compute rho_mean.
     */
    
    std::vector<double> rho_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        patch_hierarchy,
        data_context);
    
    /*
     * Compute u_mean.
     */
    
    std::vector<double> u_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "VELOCITY",
        0,
        patch_hierarchy,
        data_context);
    
    /*
     * Compute v_mean.
     */
    
    std::vector<double> v_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "VELOCITY",
        1,
        patch_hierarchy,
        data_context);
    
    /*
     * Compute p_mean.
     */
    
    std::vector<double> p_mean = getAveragedQuantityWithInhomogeneousXDirection(
        s_variable_pressure_filtered,
        0,
        patch_hierarchy,
        data_context);
    
    /*
     * Compute u_tilde.
     */
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    
    std::vector<double> rho_u_mean = getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        patch_hierarchy,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    
    std::vector<double> u_tilde(rho_u_mean);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        u_tilde[i] /= rho_mean[i];
    }
    
    /*
     * Compute v_tilde.
     */
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(1);
    
    std::vector<double> rho_v_mean = getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        patch_hierarchy,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    
    std::vector<double> v_tilde(rho_v_mean);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        v_tilde[i] /= rho_mean[i];
    }
    
    /*
     * Compute a1.
     */
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(rho_mean);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(u_mean);
    
    std::vector<double> rho_p_u_p = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    averaged_quantities.clear();
    
    std::vector<double> a1(rho_p_u_p);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        a1[i] /= rho_mean[i];
    }
    
    /*
     * Compute R22.
     */
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(zeros);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(1);
    averaged_quantities.push_back(v_tilde);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(1);
    averaged_quantities.push_back(v_tilde);
    
    std::vector<double> rho_R22 = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    averaged_quantities.clear();
    
    std::vector<double> R22(rho_R22);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        R22[i] /= rho_mean[i];
    }
    
    /*
     * Compute term II.
     */
    
    tbox::pout << "FlowModelStatisticsUtilitiesFourEqnConservative::"
        << "outputBudgetFilteredReynoldsNormalStressInYDirectionWithInhomogeneousXDirection: Compute term II" << std::endl;
    
    std::vector<double> drho_u_dx_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        "MOMENTUM",
        0,
        0,
        patch_hierarchy,
        data_context);
    
    std::vector<double> drho_v_dx_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        "MOMENTUM",
        1,
        0,
        patch_hierarchy,
        data_context);
    
    std::vector<double> drho_dx_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        0,
        patch_hierarchy,
        data_context);
    
    // quantity_names.push_back("MOMENTUM");
    variable_quantities.push_back(s_variable_momentum_filtered);
    component_indices.push_back(1);
    
    // quantity_names.push_back("VELOCITY");
    variable_quantities.push_back(s_variable_velocity_Favre_filtered);
    component_indices.push_back(1);
    
    std::vector<double> drho_v_v_dx_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        variable_quantities,
        component_indices,
        0,
        patch_hierarchy,
        data_context);
    
    // quantity_names.clear();
    variable_quantities.clear();
    component_indices.clear();
    
    std::vector<double> d_rho_u_tilde_R22_dx(finest_level_dim_0, double(0));
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        const double dR_22_tilde_dx = -(rho_R22[i]/(rho_mean[i]*rho_mean[i]))*drho_dx_mean[i] + 
            double(1)/rho_mean[i]*(drho_v_v_dx_mean[i] - double(2)*v_tilde[i]*drho_v_dx_mean[i] + 
            v_tilde[i]*v_tilde[i]*drho_dx_mean[i]);
        
        d_rho_u_tilde_R22_dx[i] = rho_u_mean[i]*dR_22_tilde_dx + R22[i]*drho_u_dx_mean[i];
    }
    
    /*
     * Compute term II in moving frame of mixing layer.
     */
    
    tbox::pout << "FlowModelStatisticsUtilitiesFourEqnConservative::"
        << "outputBudgetFilteredReynoldsNormalStressInYDirectionWithInhomogeneousXDirection: Compute term II in MF" << std::endl;
    
    std::vector<double> rho_a1_R22(rho_R22);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        rho_a1_R22[i] *= a1[i];
    }
    
    std::vector<double> d_rho_a1_R22_dx = computeDerivativeOfVector1D(
        rho_a1_R22,
        dx);
    
    /*
     * Compute term IV(1).
     */
    
    tbox::pout << "FlowModelStatisticsUtilitiesFourEqnConservative::"
        << "outputBudgetFilteredReynoldsNormalStressInYDirectionWithInhomogeneousXDirection: Compute term IV(1)" << std::endl;
    
    std::vector<double> dv_tilde_dx(finest_level_dim_0, double(0));
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        dv_tilde_dx[i] = drho_v_dx_mean[i]/rho_mean[i] - rho_v_mean[i]/(rho_mean[i]*rho_mean[i])*drho_dx_mean[i];
    }
    
    // quantity_names.push_back("MOMENTUM");
    variable_quantities.push_back(s_variable_momentum_filtered);
    component_indices.push_back(0);
    
    // quantity_names.push_back("VELOCITY");
    variable_quantities.push_back(s_variable_velocity_Favre_filtered);
    component_indices.push_back(1);
    
    std::vector<double> rho_u_v_mean = getAveragedQuantityWithInhomogeneousXDirection(
        variable_quantities,
        component_indices,
        patch_hierarchy,
        data_context);
    
    // quantity_names.clear();
    variable_quantities.clear();
    component_indices.clear();
    
    // quantity_names.push_back("MOMENTUM");
    variable_quantities.push_back(s_variable_momentum_filtered);
    component_indices.push_back(0);
    
    // quantity_names.push_back("VELOCITY");
    variable_quantities.push_back(s_variable_velocity_Favre_filtered);
    component_indices.push_back(1);
    
    // quantity_names.push_back("VELOCITY");
    variable_quantities.push_back(s_variable_velocity_Favre_filtered);
    component_indices.push_back(1);
    
    std::vector<double> drho_u_v_v_dx_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        variable_quantities,
        component_indices,
        0,
        patch_hierarchy,
        data_context);
    
    // quantity_names.clear();
    variable_quantities.clear();
    component_indices.clear();
    
    // quantity_names.push_back("MOMENTUM");
    variable_quantities.push_back(s_variable_momentum_filtered);
    component_indices.push_back(0);
    
    // quantity_names.push_back("VELOCITY");
    variable_quantities.push_back(s_variable_velocity_Favre_filtered);
    component_indices.push_back(1);
    
    std::vector<double> drho_u_v_dx_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        variable_quantities,
        component_indices,
        0,
        patch_hierarchy,
        data_context);
    
    // quantity_names.clear();
    variable_quantities.clear();
    component_indices.clear();
    
    std::vector<double> m_drho_v_pp_v_pp_u_pp_dx(finest_level_dim_0, double(0));
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        m_drho_v_pp_v_pp_u_pp_dx[i] = -(drho_u_v_v_dx_mean[i] - double(2)*rho_u_v_mean[i]*dv_tilde_dx[i] -
            double(2)*v_tilde[i]*drho_u_v_dx_mean[i] + v_tilde[i]*v_tilde[i]*drho_u_dx_mean[i] +
            double(2)*rho_u_mean[i]*v_tilde[i]*dv_tilde_dx[i] - d_rho_u_tilde_R22_dx[i]
            );
    }
    
    /*
     * Compute term IV(2).
     */
    
    tbox::pout << "FlowModelStatisticsUtilitiesFourEqnConservative::"
        << "outputBudgetFilteredReynoldsNormalStressInYDirectionWithInhomogeneousXDirection: Compute term IV(2)" << std::endl;
    
    std::vector<double> tau21_mean = getAveragedQuantityWithInhomogeneousXDirection(
        s_variable_shear_stress_filtered,
        1,
        patch_hierarchy,
        data_context);
    
    variable_quantities.push_back(s_variable_velocity_Favre_filtered);
    component_indices.push_back(1);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(-1);
    averaged_quantities.push_back(v_mean);
    
    variable_quantities.push_back(s_variable_shear_stress_filtered);
    component_indices.push_back(1);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(-1);
    averaged_quantities.push_back(tau21_mean);
    
    std::vector<double> v_p_tau21_p = getQuantityCorrelationWithInhomogeneousXDirection(
        variable_quantities,
        component_indices,
        use_reciprocal,
        derivative_directions,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    variable_quantities.clear();
    component_indices.clear();
    use_reciprocal.clear();
    derivative_directions.clear();
    averaged_quantities.clear();
    
    std::vector<double> two_dv_p_tau21_p_dx = computeDerivativeOfVector1D(
        v_p_tau21_p,
        dx);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        two_dv_p_tau21_p_dx[i] *= double(2);
    }
    
    /*
     * Compute term IV(3).
     */
    
    tbox::pout << "FlowModelStatisticsUtilitiesFourEqnConservative::"
        << "outputBudgetFilteredReynoldsNormalStressInYDirectionWithInhomogeneousXDirection: Compute term IV(3)" << std::endl;
    
    std::vector<double> tau_SFS_21_mean = getAveragedQuantityWithInhomogeneousXDirection(
        s_variable_SFS_stress,
        1,
        patch_hierarchy,
        data_context);
    
    variable_quantities.push_back(s_variable_velocity_Favre_filtered);
    component_indices.push_back(1);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(-1);
    averaged_quantities.push_back(v_mean);
    
    variable_quantities.push_back(s_variable_SFS_stress);
    component_indices.push_back(1);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(-1);
    averaged_quantities.push_back(tau_SFS_21_mean);
    
    std::vector<double> v_p_tau_SFS_21_p = getQuantityCorrelationWithInhomogeneousXDirection(
        variable_quantities,
        component_indices,
        use_reciprocal,
        derivative_directions,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    variable_quantities.clear();
    component_indices.clear();
    use_reciprocal.clear();
    derivative_directions.clear();
    averaged_quantities.clear();
    
    std::vector<double> m_2dv_p_tau_SFS_21_p_dx = computeDerivativeOfVector1D(
        v_p_tau_SFS_21_p,
        dx);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        m_2dv_p_tau_SFS_21_p_dx[i] *= (-double(2));
    }
    
    /*
     * Compute term V.
     */
    
    tbox::pout << "FlowModelStatisticsUtilitiesFourEqnConservative::"
        << "outputBudgetFilteredReynoldsNormalStressInYDirectionWithInhomogeneousXDirection: Compute term V" << std::endl;
    
    std::vector<double> dv_dy_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        "VELOCITY",
        1,
        1,
        patch_hierarchy,
        data_context);
    
    variable_quantities.push_back(s_variable_pressure_filtered);
    component_indices.push_back(0);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(-1);
    averaged_quantities.push_back(p_mean);
    
    variable_quantities.push_back(s_variable_velocity_Favre_filtered);
    component_indices.push_back(1);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(1);
    averaged_quantities.push_back(dv_dy_mean);
    
    std::vector<double> two_p_p_dv_dy_p = getQuantityCorrelationWithInhomogeneousXDirection(
        variable_quantities,
        component_indices,
        use_reciprocal,
        derivative_directions,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    variable_quantities.clear();
    component_indices.clear();
    use_reciprocal.clear();
    derivative_directions.clear();
    averaged_quantities.clear();
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        two_p_p_dv_dy_p[i] *= double(2);
    }
    
    /*
     * Compute term VI(1).
     */
    
    tbox::pout << "FlowModelStatisticsUtilitiesFourEqnConservative::"
        << "outputBudgetFilteredReynoldsNormalStressInYDirectionWithInhomogeneousXDirection: Compute term VI(1)" << std::endl;
    
    std::vector<double> dv_dx_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        "VELOCITY",
        1,
        0,
        patch_hierarchy,
        data_context);
    
    variable_quantities.push_back(s_variable_velocity_Favre_filtered);
    component_indices.push_back(1);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(0);
    averaged_quantities.push_back(dv_dx_mean);
    
    variable_quantities.push_back(s_variable_shear_stress_filtered);
    component_indices.push_back(1);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(-1);
    averaged_quantities.push_back(tau21_mean);
    
    std::vector<double> tau_21_p_dv_dx_p =
        getQuantityCorrelationWithInhomogeneousXDirection(
            variable_quantities,
            component_indices,
            use_reciprocal,
            derivative_directions,
            averaged_quantities,
            patch_hierarchy,
            data_context) ;
    
    variable_quantities.clear();
    component_indices.clear();
    use_reciprocal.clear();
    derivative_directions.clear();
    averaged_quantities.clear();
    
    std::vector<double> tau22_mean;
    
    if (d_dim == tbox::Dimension(2))
    {
        tau22_mean = getAveragedQuantityWithInhomogeneousXDirection(
            s_variable_shear_stress_filtered,
            2,
            patch_hierarchy,
            data_context);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        tau22_mean = getAveragedQuantityWithInhomogeneousXDirection(
            s_variable_shear_stress_filtered,
            3,
            patch_hierarchy,
            data_context);
    }
    
    variable_quantities.push_back(s_variable_velocity_Favre_filtered);
    component_indices.push_back(1);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(1);
    averaged_quantities.push_back(dv_dy_mean);
    
    if (d_dim == tbox::Dimension(2))
    {
        variable_quantities.push_back(s_variable_shear_stress_filtered);
        component_indices.push_back(2);
        use_reciprocal.push_back(false);
        derivative_directions.push_back(-1);
        averaged_quantities.push_back(tau22_mean);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        variable_quantities.push_back(s_variable_shear_stress_filtered);
        component_indices.push_back(3);
        use_reciprocal.push_back(false);
        derivative_directions.push_back(-1);
        averaged_quantities.push_back(tau22_mean);
    }
    
    std::vector<double> tau_22_p_dv_dy_p =
        getQuantityCorrelationWithInhomogeneousXDirection(
            variable_quantities,
            component_indices,
            use_reciprocal,
            derivative_directions,
            averaged_quantities,
            patch_hierarchy,
            data_context) ;
    
    variable_quantities.clear();
    component_indices.clear();
    use_reciprocal.clear();
    derivative_directions.clear();
    averaged_quantities.clear();
    
    std::vector<double> tau_23_p_dv_dz_p;
    if (d_dim == tbox::Dimension(3))
    {
        std::vector<double> dv_dz_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
            "VELOCITY",
            1,
            2,
            patch_hierarchy,
            data_context);
        
        std::vector<double> tau23_mean = getAveragedQuantityWithInhomogeneousXDirection(
            s_variable_shear_stress_filtered,
            4,
            patch_hierarchy,
            data_context);
        
        variable_quantities.push_back(s_variable_velocity_Favre_filtered);
        component_indices.push_back(1);
        use_reciprocal.push_back(false);
        derivative_directions.push_back(2);
        averaged_quantities.push_back(dv_dz_mean);
        
        variable_quantities.push_back(s_variable_shear_stress_filtered);
        component_indices.push_back(4);
        use_reciprocal.push_back(false);
        derivative_directions.push_back(-1);
        averaged_quantities.push_back(tau23_mean);
        
        tau_23_p_dv_dz_p =
            getQuantityCorrelationWithInhomogeneousXDirection(
                variable_quantities,
                component_indices,
                use_reciprocal,
                derivative_directions,
                averaged_quantities,
                patch_hierarchy,
                data_context);
        
        variable_quantities.clear();
        component_indices.clear();
        use_reciprocal.clear();
        derivative_directions.clear();
        averaged_quantities.clear();
    }
    
    std::vector<double> m_2tau_2i_p_dv_p_dxi(finest_level_dim_0, double(0));
    
    if (d_dim == tbox::Dimension(2))
    {
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            m_2tau_2i_p_dv_p_dxi[i] = -double(2)*(tau_21_p_dv_dx_p[i] + tau_22_p_dv_dy_p[i]);
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            m_2tau_2i_p_dv_p_dxi[i] = -double(2)*(tau_21_p_dv_dx_p[i] + tau_22_p_dv_dy_p[i] + tau_23_p_dv_dz_p[i]);
        }
    }
    
    /*
     * Compute term VI(2).
     */
    
    tbox::pout << "FlowModelStatisticsUtilitiesFourEqnConservative::"
        << "outputBudgetFilteredReynoldsNormalStressInYDirectionWithInhomogeneousXDirection: Compute term VI(2)" << std::endl;
    
    variable_quantities.push_back(s_variable_velocity_Favre_filtered);
    component_indices.push_back(1);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(0);
    averaged_quantities.push_back(dv_dx_mean);
    
    variable_quantities.push_back(s_variable_SFS_stress);
    component_indices.push_back(1);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(-1);
    averaged_quantities.push_back(tau_SFS_21_mean);
    
    std::vector<double> tau_SFS_21_p_dv_dx_p =
        getQuantityCorrelationWithInhomogeneousXDirection(
            variable_quantities,
            component_indices,
            use_reciprocal,
            derivative_directions,
            averaged_quantities,
            patch_hierarchy,
            data_context) ;
    
    variable_quantities.clear();
    component_indices.clear();
    use_reciprocal.clear();
    derivative_directions.clear();
    averaged_quantities.clear();
    
    std::vector<double> tau_SFS_22_mean;
    
    if (d_dim == tbox::Dimension(2))
    {
        tau_SFS_22_mean = getAveragedQuantityWithInhomogeneousXDirection(
            s_variable_SFS_stress,
            2,
            patch_hierarchy,
            data_context);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        tau_SFS_22_mean = getAveragedQuantityWithInhomogeneousXDirection(
            s_variable_SFS_stress,
            3,
            patch_hierarchy,
            data_context);
    }
    
    variable_quantities.push_back(s_variable_velocity_Favre_filtered);
    component_indices.push_back(1);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(1);
    averaged_quantities.push_back(dv_dy_mean);
    
    if (d_dim == tbox::Dimension(2))
    {
      variable_quantities.push_back(s_variable_SFS_stress);
        component_indices.push_back(2);
        use_reciprocal.push_back(false);
        derivative_directions.push_back(-1);
        averaged_quantities.push_back(tau_SFS_22_mean);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        variable_quantities.push_back(s_variable_SFS_stress);
        component_indices.push_back(3);
        use_reciprocal.push_back(false);
        derivative_directions.push_back(-1);
        averaged_quantities.push_back(tau_SFS_22_mean);
    }
    
    std::vector<double> tau_SFS_22_p_dv_dy_p =
        getQuantityCorrelationWithInhomogeneousXDirection(
            variable_quantities,
            component_indices,
            use_reciprocal,
            derivative_directions,
            averaged_quantities,
            patch_hierarchy,
            data_context) ;
    
    variable_quantities.clear();
    component_indices.clear();
    use_reciprocal.clear();
    derivative_directions.clear();
    averaged_quantities.clear();
    
    std::vector<double> tau_SFS_23_p_dv_dz_p;
    if (d_dim == tbox::Dimension(3))
    {
        std::vector<double> dv_dz_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
            "VELOCITY",
            1,
            2,
            patch_hierarchy,
            data_context);
        
        std::vector<double> tau_SFS_23_mean = getAveragedQuantityWithInhomogeneousXDirection(
            s_variable_SFS_stress,
            4,
            patch_hierarchy,
            data_context);
        
        variable_quantities.push_back(s_variable_velocity_Favre_filtered);
        component_indices.push_back(1);
        use_reciprocal.push_back(false);
        derivative_directions.push_back(2);
        averaged_quantities.push_back(dv_dz_mean);
        
        variable_quantities.push_back(s_variable_SFS_stress);
        component_indices.push_back(4);
        use_reciprocal.push_back(false);
        derivative_directions.push_back(-1);
        averaged_quantities.push_back(tau_SFS_23_mean);
        
        tau_SFS_23_p_dv_dz_p =
            getQuantityCorrelationWithInhomogeneousXDirection(
                variable_quantities,
                component_indices,
                use_reciprocal,
                derivative_directions,
                averaged_quantities,
                patch_hierarchy,
                data_context);
        
        variable_quantities.clear();
        component_indices.clear();
        use_reciprocal.clear();
        derivative_directions.clear();
        averaged_quantities.clear();
    }
    
    std::vector<double> two_tau_SFS_2i_p_dv_p_dxi(finest_level_dim_0, double(0));
    
    if (d_dim == tbox::Dimension(2))
    {
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            two_tau_SFS_2i_p_dv_p_dxi[i] = double(2)*(tau_SFS_21_p_dv_dx_p[i] + tau_SFS_22_p_dv_dy_p[i]);
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            two_tau_SFS_2i_p_dv_p_dxi[i] = double(2)*(tau_SFS_21_p_dv_dx_p[i] + tau_SFS_22_p_dv_dy_p[i] + tau_SFS_23_p_dv_dz_p[i]);
        }
    }
    
    /*
     * Output budget.
     */
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_output;
        f_output.open(stat_dump_filename.c_str(), std::ios_base::app | std::ios::out | std::ios::binary);
        
        if (!f_output.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output R22 budget!"
                << std::endl);
        }
        
        f_output.write((char*)&output_time, sizeof(double));
        f_output.write((char*)&rho_R22[0], sizeof(double)*rho_R22.size());
        // Term II.
        f_output.write((char*)&d_rho_u_tilde_R22_dx[0], sizeof(double)*d_rho_u_tilde_R22_dx.size());
        
        // Term IV(1).
        f_output.write((char*)&m_drho_v_pp_v_pp_u_pp_dx[0], sizeof(double)*m_drho_v_pp_v_pp_u_pp_dx.size());
        // Term IV(2).
        f_output.write((char*)&two_dv_p_tau21_p_dx[0], sizeof(double)*two_dv_p_tau21_p_dx.size());
        // Term IV(3).
        f_output.write((char*)&m_2dv_p_tau_SFS_21_p_dx[0], sizeof(double)*m_2dv_p_tau_SFS_21_p_dx.size());
        
        // Term V.
        f_output.write((char*)&two_p_p_dv_dy_p[0], sizeof(double)*two_p_p_dv_dy_p.size());
        
        // Term VI(1).
        f_output.write((char*)&m_2tau_2i_p_dv_p_dxi[0], sizeof(double)*m_2tau_2i_p_dv_p_dxi.size());
        // Term VI(2).
        f_output.write((char*)&two_tau_SFS_2i_p_dv_p_dxi[0], sizeof(double)*two_tau_SFS_2i_p_dv_p_dxi.size());
        
        // Term II in moving frame of mixing layer.
        f_output.write((char*)&d_rho_a1_R22_dx[0], sizeof(double)*d_rho_a1_R22_dx.size());
        
        // SFS stress profile.
        f_output.write((char*)&tau_SFS_22_mean[0], sizeof(double)*tau_SFS_22_mean.size());
        
        f_output.close();
    }
}


/*
 * Output budget of Reynolds normal stress in z-direction with inhomogeneous x-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::
outputBudgetFilteredReynoldsNormalStressInZDirectionWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double output_time) const
{
    tbox::pout << "FlowModelStatisticsUtilitiesFourEqnConservative::"
        << "outputBudgetFilteredReynoldsNormalStressInZDirectionWithInhomogeneousXDirection: start" << std::endl;
    
    const int finest_level_dim_0 = getRefinedDomainNumberOfPointsX(patch_hierarchy);
    
    const double dx = getRefinedDomainGridSpacingX(patch_hierarchy);
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    std::vector<std::vector<double> > averaged_quantities;
    std::vector<bool> use_reciprocal;
    std::vector<int> derivative_directions;
    
    std::vector<boost::shared_ptr<pdat::CellVariable<double> > > variable_quantities;
    
    std::vector<double> zeros(finest_level_dim_0, double(0));
    
    /*
     * Compute rho_mean.
     */
    
    std::vector<double> rho_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        patch_hierarchy,
        data_context);
    
    /*
     * Compute u_mean.
     */
    
    std::vector<double> u_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "VELOCITY",
        0,
        patch_hierarchy,
        data_context);
    
    /*
     * Compute w_mean.
     */
    
    std::vector<double> w_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "VELOCITY",
        2,
        patch_hierarchy,
        data_context);
    
    /*
     * Compute p_mean.
     */
    
    std::vector<double> p_mean = getAveragedQuantityWithInhomogeneousXDirection(
        s_variable_pressure_filtered,
        0,
        patch_hierarchy,
        data_context);
    
    /*
     * Compute u_tilde.
     */
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    
    std::vector<double> rho_u_mean = getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        patch_hierarchy,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    
    std::vector<double> u_tilde(rho_u_mean);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        u_tilde[i] /= rho_mean[i];
    }
    
    /*
     * Compute w_tilde.
     */
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(2);
    
    std::vector<double> rho_w_mean = getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        patch_hierarchy,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    
    std::vector<double> w_tilde(rho_w_mean);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        w_tilde[i] /= rho_mean[i];
    }
    
    /*
     * Compute a1.
     */
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(rho_mean);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(u_mean);
    
    std::vector<double> rho_p_u_p = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    averaged_quantities.clear();
    
    std::vector<double> a1(rho_p_u_p);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        a1[i] /= rho_mean[i];
    }
    
    /*
     * Compute R33.
     */
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(zeros);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(2);
    averaged_quantities.push_back(w_tilde);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(2);
    averaged_quantities.push_back(w_tilde);
    
    std::vector<double> rho_R33 = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    averaged_quantities.clear();
    
    std::vector<double> R33(rho_R33);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        R33[i] /= rho_mean[i];
    }
    
    /*
     * Compute term II.
     */
    
    tbox::pout << "FlowModelStatisticsUtilitiesFourEqnConservative::"
        << "outputBudgetFilteredReynoldsNormalStressInZDirectionWithInhomogeneousXDirection: Compute term II" << std::endl;
    
    std::vector<double> drho_u_dx_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        "MOMENTUM",
        0,
        0,
        patch_hierarchy,
        data_context);
    
    std::vector<double> drho_w_dx_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        "MOMENTUM",
        2,
        0,
        patch_hierarchy,
        data_context);
    
    std::vector<double> drho_dx_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        0,
        patch_hierarchy,
        data_context);
    
    // quantity_names.push_back("MOMENTUM");
    variable_quantities.push_back(s_variable_momentum_filtered);
    component_indices.push_back(2);
    
    // quantity_names.push_back("VELOCITY");
    variable_quantities.push_back(s_variable_velocity_Favre_filtered);
    component_indices.push_back(2);
    
    std::vector<double> drho_w_w_dx_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        variable_quantities,
        component_indices,
        0,
        patch_hierarchy,
        data_context);
    
    // quantity_names.clear();
    variable_quantities.clear();
    component_indices.clear();
    
    std::vector<double> d_rho_u_tilde_R33_dx(finest_level_dim_0, double(0));
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        const double dR_33_tilde_dx = -(rho_R33[i]/(rho_mean[i]*rho_mean[i]))*drho_dx_mean[i] + 
            double(1)/rho_mean[i]*(drho_w_w_dx_mean[i] - double(2)*w_tilde[i]*drho_w_dx_mean[i] + 
            w_tilde[i]*w_tilde[i]*drho_dx_mean[i]);
        
        d_rho_u_tilde_R33_dx[i] = rho_u_mean[i]*dR_33_tilde_dx + R33[i]*drho_u_dx_mean[i];
    }
    
    /*
     * Compute term II in moving frame of mixing layer.
     */
    
    tbox::pout << "FlowModelStatisticsUtilitiesFourEqnConservative::"
        << "outputBudgetFilteredReynoldsNormalStressInZDirectionWithInhomogeneousXDirection: Compute term II in MF" << std::endl;
    
    std::vector<double> rho_a1_R33(rho_R33);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        rho_a1_R33[i] *= a1[i];
    }
    
    std::vector<double> d_rho_a1_R33_dx = computeDerivativeOfVector1D(
        rho_a1_R33,
        dx);
    
    /*
     * Compute term IV(1).
     */
    
    tbox::pout << "FlowModelStatisticsUtilitiesFourEqnConservative::"
        << "outputBudgetFilteredReynoldsNormalStressInZDirectionWithInhomogeneousXDirection: Compute term IV(1)" << std::endl;
    
    std::vector<double> dw_tilde_dx(finest_level_dim_0, double(0));
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        dw_tilde_dx[i] = drho_w_dx_mean[i]/rho_mean[i] - rho_w_mean[i]/(rho_mean[i]*rho_mean[i])*drho_dx_mean[i];
    }
    
    // quantity_names.push_back("MOMENTUM");
    variable_quantities.push_back(s_variable_momentum_filtered);
    component_indices.push_back(0);
    
    // quantity_names.push_back("VELOCITY");
    variable_quantities.push_back(s_variable_velocity_Favre_filtered);
    component_indices.push_back(2);
    
    std::vector<double> rho_u_w_mean = getAveragedQuantityWithInhomogeneousXDirection(
        variable_quantities,
        component_indices,
        patch_hierarchy,
        data_context);
    
    // quantity_names.clear();
    variable_quantities.clear();
    component_indices.clear();
    
    // quantity_names.push_back("MOMENTUM");
    variable_quantities.push_back(s_variable_momentum_filtered);
    component_indices.push_back(0);
    
    // quantity_names.push_back("VELOCITY");
    variable_quantities.push_back(s_variable_velocity_Favre_filtered);
    component_indices.push_back(2);
    
    // quantity_names.push_back("VELOCITY");
    variable_quantities.push_back(s_variable_velocity_Favre_filtered);
    component_indices.push_back(2);
    
    std::vector<double> drho_u_w_w_dx_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        variable_quantities,
        component_indices,
        0,
        patch_hierarchy,
        data_context);
    
    // quantity_names.clear();
    variable_quantities.clear();
    component_indices.clear();
    
    // quantity_names.push_back("MOMENTUM");
    variable_quantities.push_back(s_variable_momentum_filtered);
    component_indices.push_back(0);
    
    // quantity_names.push_back("VELOCITY");
    variable_quantities.push_back(s_variable_velocity_Favre_filtered);
    component_indices.push_back(2);
    
    std::vector<double> drho_u_w_dx_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        variable_quantities,
        component_indices,
        0,
        patch_hierarchy,
        data_context);
    
    // quantity_names.clear();
    variable_quantities.clear();
    component_indices.clear();
    
    std::vector<double> m_drho_w_pp_w_pp_u_pp_dx(finest_level_dim_0, double(0));
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        m_drho_w_pp_w_pp_u_pp_dx[i] = -(drho_u_w_w_dx_mean[i] - double(2)*rho_u_w_mean[i]*dw_tilde_dx[i] -
            double(2)*w_tilde[i]*drho_u_w_dx_mean[i] + w_tilde[i]*w_tilde[i]*drho_u_dx_mean[i] +
            double(2)*rho_u_mean[i]*w_tilde[i]*dw_tilde_dx[i] - d_rho_u_tilde_R33_dx[i]
            );
    }
    
    /*
     * Compute term IV(2).
     */
    
    tbox::pout << "FlowModelStatisticsUtilitiesFourEqnConservative::"
        << "outputBudgetFilteredReynoldsNormalStressInZDirectionWithInhomogeneousXDirection: Compute term IV(2)" << std::endl;
    
    std::vector<double> tau31_mean = getAveragedQuantityWithInhomogeneousXDirection(
        s_variable_shear_stress_filtered,
        2,
        patch_hierarchy,
        data_context);
    
    variable_quantities.push_back(s_variable_velocity_Favre_filtered);
    component_indices.push_back(2);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(-1);
    averaged_quantities.push_back(w_mean);
    
    variable_quantities.push_back(s_variable_shear_stress_filtered);
    component_indices.push_back(2);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(-1);
    averaged_quantities.push_back(tau31_mean);
    
    std::vector<double> w_p_tau31_p = getQuantityCorrelationWithInhomogeneousXDirection(
        variable_quantities,
        component_indices,
        use_reciprocal,
        derivative_directions,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    variable_quantities.clear();
    component_indices.clear();
    use_reciprocal.clear();
    derivative_directions.clear();
    averaged_quantities.clear();
    
    std::vector<double> two_dw_p_tau31_p_dx = computeDerivativeOfVector1D(
        w_p_tau31_p,
        dx);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        two_dw_p_tau31_p_dx[i] *= double(2);
    }
    
    /*
     * Compute term IV(3).
     */
    
    tbox::pout << "FlowModelStatisticsUtilitiesFourEqnConservative::"
        << "outputBudgetFilteredReynoldsNormalStressInZDirectionWithInhomogeneousXDirection: Compute term IV(3)" << std::endl;
    
    std::vector<double> tau_SFS_31_mean = getAveragedQuantityWithInhomogeneousXDirection(
        s_variable_SFS_stress,
        2,
        patch_hierarchy,
        data_context);
    
    variable_quantities.push_back(s_variable_velocity_Favre_filtered);
    component_indices.push_back(2);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(-1);
    averaged_quantities.push_back(w_mean);
    
    variable_quantities.push_back(s_variable_SFS_stress);
    component_indices.push_back(2);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(-1);
    averaged_quantities.push_back(tau_SFS_31_mean);
    
    std::vector<double> w_p_tau_SFS_31_p = getQuantityCorrelationWithInhomogeneousXDirection(
        variable_quantities,
        component_indices,
        use_reciprocal,
        derivative_directions,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    variable_quantities.clear();
    component_indices.clear();
    use_reciprocal.clear();
    derivative_directions.clear();
    averaged_quantities.clear();
    
    std::vector<double> m_2dw_p_tau_SFS_31_p_dx = computeDerivativeOfVector1D(
        w_p_tau_SFS_31_p,
        dx);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        m_2dw_p_tau_SFS_31_p_dx[i] *= (-double(2));
    }
    
    /*
     * Compute term V.
     */
    
    tbox::pout << "FlowModelStatisticsUtilitiesFourEqnConservative::"
        << "outputBudgetFilteredReynoldsNormalStressInZDirectionWithInhomogeneousXDirection: Compute term V" << std::endl;
    
    std::vector<double> dw_dz_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        "VELOCITY",
        2,
        2,
        patch_hierarchy,
        data_context);
    
    variable_quantities.push_back(s_variable_pressure_filtered);
    component_indices.push_back(0);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(-1);
    averaged_quantities.push_back(p_mean);
    
    variable_quantities.push_back(s_variable_velocity_Favre_filtered);
    component_indices.push_back(2);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(2);
    averaged_quantities.push_back(dw_dz_mean);
    
    std::vector<double> two_p_p_dw_dz_p = getQuantityCorrelationWithInhomogeneousXDirection(
        variable_quantities,
        component_indices,
        use_reciprocal,
        derivative_directions,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    variable_quantities.clear();
    component_indices.clear();
    use_reciprocal.clear();
    derivative_directions.clear();
    averaged_quantities.clear();
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        two_p_p_dw_dz_p[i] *= double(2);
    }
    
    /*
     * Compute term VI(1).
     */
    
    tbox::pout << "FlowModelStatisticsUtilitiesFourEqnConservative::"
        << "outputBudgetFilteredReynoldsNormalStressInZDirectionWithInhomogeneousXDirection: Compute term VI(1)" << std::endl;
    
    std::vector<double> dw_dx_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        "VELOCITY",
        2,
        0,
        patch_hierarchy,
        data_context);
    
    variable_quantities.push_back(s_variable_velocity_Favre_filtered);
    component_indices.push_back(2);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(0);
    averaged_quantities.push_back(dw_dx_mean);
    
    variable_quantities.push_back(s_variable_shear_stress_filtered);
    component_indices.push_back(2);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(-1);
    averaged_quantities.push_back(tau31_mean);
    
    std::vector<double> tau_31_p_dw_dx_p =
        getQuantityCorrelationWithInhomogeneousXDirection(
            variable_quantities,
            component_indices,
            use_reciprocal,
            derivative_directions,
            averaged_quantities,
            patch_hierarchy,
            data_context) ;
    
    variable_quantities.clear();
    component_indices.clear();
    use_reciprocal.clear();
    derivative_directions.clear();
    averaged_quantities.clear();
    
    std::vector<double> dw_dy_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        "VELOCITY",
        2,
        1,
        patch_hierarchy,
        data_context);
    
    std::vector<double> tau32_mean = getAveragedQuantityWithInhomogeneousXDirection(
        s_variable_shear_stress_filtered,
        4,
        patch_hierarchy,
        data_context);
    
    variable_quantities.push_back(s_variable_velocity_Favre_filtered);
    component_indices.push_back(2);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(1);
    averaged_quantities.push_back(dw_dy_mean);
    
    variable_quantities.push_back(s_variable_shear_stress_filtered);
    component_indices.push_back(4);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(-1);
    averaged_quantities.push_back(tau32_mean);
    
    std::vector<double> tau_32_p_dw_dy_p =
        getQuantityCorrelationWithInhomogeneousXDirection(
            variable_quantities,
            component_indices,
            use_reciprocal,
            derivative_directions,
            averaged_quantities,
            patch_hierarchy,
            data_context) ;
    
    variable_quantities.clear();
    component_indices.clear();
    use_reciprocal.clear();
    derivative_directions.clear();
    averaged_quantities.clear();
    
    std::vector<double> tau33_mean = getAveragedQuantityWithInhomogeneousXDirection(
        s_variable_shear_stress_filtered,
        5,
        patch_hierarchy,
        data_context);
    
    variable_quantities.push_back(s_variable_velocity_Favre_filtered);
    component_indices.push_back(2);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(2);
    averaged_quantities.push_back(dw_dz_mean);
    
    variable_quantities.push_back(s_variable_shear_stress_filtered);
    component_indices.push_back(5);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(-1);
    averaged_quantities.push_back(tau33_mean);
    
    std::vector<double> tau_33_p_dw_dz_p =
        getQuantityCorrelationWithInhomogeneousXDirection(
            variable_quantities,
            component_indices,
            use_reciprocal,
            derivative_directions,
            averaged_quantities,
            patch_hierarchy,
            data_context);
    
    variable_quantities.clear();
    component_indices.clear();
    use_reciprocal.clear();
    derivative_directions.clear();
    averaged_quantities.clear();
    
    std::vector<double> m_2tau_3i_p_dw_p_dxi(finest_level_dim_0, double(0));
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        m_2tau_3i_p_dw_p_dxi[i] = -double(2)*(tau_31_p_dw_dx_p[i] + tau_32_p_dw_dy_p[i] + tau_33_p_dw_dz_p[i]);
    }
    
    /*
     * Compute term VI(2).
     */
    
    tbox::pout << "FlowModelStatisticsUtilitiesFourEqnConservative::"
        << "outputBudgetFilteredReynoldsNormalStressInZDirectionWithInhomogeneousXDirection: Compute term VI(2)" << std::endl;
    
    variable_quantities.push_back(s_variable_velocity_Favre_filtered);
    component_indices.push_back(2);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(0);
    averaged_quantities.push_back(dw_dx_mean);
    
    variable_quantities.push_back(s_variable_SFS_stress);
    component_indices.push_back(2);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(-1);
    averaged_quantities.push_back(tau_SFS_31_mean);
    
    std::vector<double> tau_SFS_31_p_dw_dx_p =
        getQuantityCorrelationWithInhomogeneousXDirection(
            variable_quantities,
            component_indices,
            use_reciprocal,
            derivative_directions,
            averaged_quantities,
            patch_hierarchy,
            data_context) ;
    
    variable_quantities.clear();
    component_indices.clear();
    use_reciprocal.clear();
    derivative_directions.clear();
    averaged_quantities.clear();
    
    std::vector<double> tau_SFS_32_mean = getAveragedQuantityWithInhomogeneousXDirection(
        s_variable_SFS_stress,
        4,
        patch_hierarchy,
        data_context);
    
    variable_quantities.push_back(s_variable_velocity_Favre_filtered);
    component_indices.push_back(2);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(1);
    averaged_quantities.push_back(dw_dy_mean);
    
    variable_quantities.push_back(s_variable_SFS_stress);
    component_indices.push_back(4);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(-1);
    averaged_quantities.push_back(tau_SFS_32_mean);
    
    std::vector<double> tau_SFS_32_p_dw_dy_p =
        getQuantityCorrelationWithInhomogeneousXDirection(
            variable_quantities,
            component_indices,
            use_reciprocal,
            derivative_directions,
            averaged_quantities,
            patch_hierarchy,
            data_context) ;
    
    variable_quantities.clear();
    component_indices.clear();
    use_reciprocal.clear();
    derivative_directions.clear();
    averaged_quantities.clear();
    
    std::vector<double> tau_SFS_33_mean = getAveragedQuantityWithInhomogeneousXDirection(
        s_variable_SFS_stress,
        5,
        patch_hierarchy,
        data_context);
    
    variable_quantities.push_back(s_variable_velocity_Favre_filtered);
    component_indices.push_back(2);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(2);
    averaged_quantities.push_back(dw_dz_mean);
    
    variable_quantities.push_back(s_variable_SFS_stress);
    component_indices.push_back(5);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(-1);
    averaged_quantities.push_back(tau_SFS_33_mean);
    
    std::vector<double> tau_SFS_33_p_dw_dz_p =
        getQuantityCorrelationWithInhomogeneousXDirection(
            variable_quantities,
            component_indices,
            use_reciprocal,
            derivative_directions,
            averaged_quantities,
            patch_hierarchy,
            data_context);
    
    variable_quantities.clear();
    component_indices.clear();
    use_reciprocal.clear();
    derivative_directions.clear();
    averaged_quantities.clear();
    
    std::vector<double> two_tau_SFS_3i_p_dw_p_dxi(finest_level_dim_0, double(0));
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        two_tau_SFS_3i_p_dw_p_dxi[i] = double(2)*(tau_SFS_31_p_dw_dx_p[i] + tau_SFS_32_p_dw_dy_p[i] + tau_SFS_33_p_dw_dz_p[i]);
    }
    
    /*
     * Output budget.
     */
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_output;
        f_output.open(stat_dump_filename.c_str(), std::ios_base::app | std::ios::out | std::ios::binary);
        
        if (!f_output.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output R33 budget!"
                << std::endl);
        }
        
        f_output.write((char*)&output_time, sizeof(double));
        f_output.write((char*)&rho_R33[0], sizeof(double)*rho_R33.size());
        // Term II.
        f_output.write((char*)&d_rho_u_tilde_R33_dx[0], sizeof(double)*d_rho_u_tilde_R33_dx.size());
        
        // Term IV(1).
        f_output.write((char*)&m_drho_w_pp_w_pp_u_pp_dx[0], sizeof(double)*m_drho_w_pp_w_pp_u_pp_dx.size());
        // Term IV(2).
        f_output.write((char*)&two_dw_p_tau31_p_dx[0], sizeof(double)*two_dw_p_tau31_p_dx.size());
        // Term IV(3).
        f_output.write((char*)&m_2dw_p_tau_SFS_31_p_dx[0], sizeof(double)*m_2dw_p_tau_SFS_31_p_dx.size());
        
        // Term V.
        f_output.write((char*)&two_p_p_dw_dz_p[0], sizeof(double)*two_p_p_dw_dz_p.size());
        
        // Term VI(1).
        f_output.write((char*)&m_2tau_3i_p_dw_p_dxi[0], sizeof(double)*m_2tau_3i_p_dw_p_dxi.size());
        // Term VI(1).
        f_output.write((char*)&two_tau_SFS_3i_p_dw_p_dxi[0], sizeof(double)*two_tau_SFS_3i_p_dw_p_dxi.size());
        
        // Term II in moving frame of mixing layer.
        f_output.write((char*)&d_rho_a1_R33_dx[0], sizeof(double)*d_rho_a1_R33_dx.size());
        
        // SFS stress profile.
        f_output.write((char*)&tau_SFS_33_mean[0], sizeof(double)*tau_SFS_33_mean.size());
        
        f_output.close();
    }
}
