#include "flow/refinement_taggers/gradient_tagger/GradientTagger.hpp"

#include "boost/lexical_cast.hpp"

// #define PLOTTING_GRADIENT_TAGGER

#define EPSILON 1e-40

GradientTagger::GradientTagger(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
    const int& num_species,
    const boost::shared_ptr<FlowModel>& flow_model,
    const boost::shared_ptr<tbox::Database>& gradient_tagger_db):
        d_object_name(object_name),
        d_dim(dim),
        d_grid_geometry(grid_geometry),
        d_num_gradient_ghosts(hier::IntVector::getZero(d_dim)),
        d_num_species(num_species),
        d_flow_model(flow_model)
{
    if (gradient_tagger_db != nullptr)
    {
        std::vector<std::string> sensor_keys = gradient_tagger_db->getAllKeys();
        
        const int num_keys = static_cast<int>(sensor_keys.size());
        
        if (gradient_tagger_db->keyExists("gradient_sensors"))
        {
            d_gradient_sensors = gradient_tagger_db->getStringVector("gradient_sensors");
        }
        else if (gradient_tagger_db->keyExists("d_gradient_sensors"))
        {
            d_gradient_sensors = gradient_tagger_db->getStringVector("d_gradient_sensors");
        }
        else
        {
            TBOX_WARNING(d_object_name
                << ": "
                << "No key 'gradient_sensors'/'d_gradient_sensors' found in data for"
                << " Gradient_tagger. No refinement with gradient sensors will occur."
                << std::endl);
        }
        
        std::vector<std::string> sensor_keys_defined(num_keys);
        int sensor_keys_count = 0;
        boost::shared_ptr<tbox::Database> sensor_db;
        for (int i = 0; i < num_keys; i++)
        {
            std::string sensor_key = sensor_keys[i];
            sensor_db.reset();
            
            if (!((sensor_key == "gradient_sensors") || (sensor_key == "d_gradient_sensors")))
            {
                if (!((sensor_key == "JAMESON_GRADIENT")))
                {
                    TBOX_ERROR(d_object_name
                        << ": "
                        << "Unknown sensor: "
                        << sensor_key
                        << "\nin input."
                        << std::endl);
                }
                else
                {
                    sensor_db = gradient_tagger_db->getDatabase(sensor_key);
                    sensor_keys_defined[sensor_keys_count] = sensor_key;
                    sensor_keys_count++;
                }
                
                if (sensor_db && sensor_key == "JAMESON_GRADIENT")
                {
                    d_gradient_sensor_Jameson.reset(new GradientSensorJameson(
                        "Jameson gradient sensor",
                        d_dim));
                    
                    if (sensor_db->keyExists("Jameson_gradient_variables"))
                    {
                        d_Jameson_gradient_variables = sensor_db->getStringVector("Jameson_gradient_variables");
                    }
                    else if (sensor_db->keyExists("d_Jameson_gradient_variables"))
                    {
                        d_Jameson_gradient_variables = sensor_db->getStringVector("d_Jameson_gradient_variables");
                    }
                    else
                    {
                        TBOX_ERROR(d_object_name
                            << ": "
                            << "No key 'Jameson_gradient_variables'/'d_Jameson_gradient_variables' found in data for "
                            << sensor_key
                            << "."
                            << std::endl);
                    }
                    
                    for (int vi = 0; vi < static_cast<int>(d_Jameson_gradient_variables.size()); vi++)
                    {
                        std::string variable_key = d_Jameson_gradient_variables[vi];
                        
                        if (!((variable_key == "DENSITY") ||
                              (variable_key == "TOTAL_ENERGY") ||
                              (variable_key == "PRESSURE") ||
                              (variable_key == "DILATATION") ||
                              (variable_key == "ENSTROPHY") ||
                              (variable_key == "MASS_FRACTION")))
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "Unknown/unsupported variable: "
                                << variable_key
                                << "\nin input."
                                << std::endl);
                        }
                    }
                    
                    if (sensor_db->keyExists("Jameson_gradient_tol"))
                    {
                        d_Jameson_gradient_tol = sensor_db->getDoubleVector("Jameson_gradient_tol");
                    }
                    else if (sensor_db->keyExists("d_Jameson_gradient_tol"))
                    {
                        d_Jameson_gradient_tol = sensor_db->getDoubleVector("d_Jameson_gradient_tol");
                    }
                    else
                    {
                        TBOX_ERROR(d_object_name
                            << ": "
                            << "No key 'Jameson_gradient_tol'/'d_Jameson_gradient_tol' found in data for "
                            << sensor_key
                            << "."
                            << std::endl);
                    }
                    
                    if (static_cast<int>(d_Jameson_gradient_variables.size()) != static_cast<int>(d_Jameson_gradient_tol.size()))
                    {
                        TBOX_ERROR(d_object_name
                            << ": "
                            << "The numbers of variables and tolerances provided don't match for "
                            << sensor_key
                            << "."
                            << std::endl);
                    }
                }
            }
        } // Loop over sensors.
        
        /*
         * Check that input is found for each string identifier in key list.
         */
        for (int ki = 0;
             ki < static_cast<int>(d_gradient_sensors.size());
             ki++)
        {
            std::string use_key = d_gradient_sensors[ki];
            bool key_found = false;
            for (int k1 = 0; k1 < sensor_keys_count; k1++)
            {
                std::string def_key = sensor_keys_defined[k1];
                if (def_key == use_key)
                {
                    key_found = true;
                }
            }
            
            if (!key_found)
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "No input found for specified gradient sensor: "
                    << d_gradient_sensors[ki]
                    << "."
                    << std::endl);
            }
        }
    }
    else
    {
        TBOX_WARNING(d_object_name
            << ": "
            << "Key data 'Gradient_tagger' not found in input/restart database."
            << " No refinement with gradient sensors will occur."
            << std::endl);
    }
    
    /*
     * Compute the number of ghost cells requried by the gradient detector.
     */
    
    for (int si = 0;
         si < static_cast<int>(d_gradient_sensors.size());
         si++)
    {
        std::string sensor_key = d_gradient_sensors[si];
        
        if (sensor_key == "JAMESON_GRADIENT")
        {
            d_num_gradient_ghosts = hier::IntVector::max(
                d_num_gradient_ghosts,
                d_gradient_sensor_Jameson->getGradientSensorNumberOfGhostCells());
        }
    }
}


/*
 * Register the temporary variables used in gradient tagger class.
 */
void
GradientTagger::registerGradientTaggerVariables(
    RungeKuttaLevelIntegrator* integrator)
{
    for (int si = 0;
             si < static_cast<int>(d_gradient_sensors.size());
             si++)
    {
        std::string sensor_key = d_gradient_sensors[si];
        
        if (sensor_key == "JAMESON_GRADIENT")
        {
            for (int vi = 0; vi < static_cast<int>(d_Jameson_gradient_variables.size()); vi++)
            {
                // Get the key of the current variable.
                std::string variable_key = d_Jameson_gradient_variables[vi];
                
                if (variable_key == "DENSITY")
                {
                    d_Jameson_density_gradient =
                        boost::shared_ptr<pdat::CellVariable<double> > (
                            new pdat::CellVariable<double>(d_dim, "Jameson density gradient", 1));
                }
                else if (variable_key == "TOTAL_ENERGY")
                {
                    d_Jameson_total_energy_gradient =
                        boost::shared_ptr<pdat::CellVariable<double> > (
                            new pdat::CellVariable<double>(d_dim, "Jameson total_energy gradient", 1));
                }
                else if (variable_key == "PRESSURE")
                {
                    d_Jameson_pressure_gradient =
                        boost::shared_ptr<pdat::CellVariable<double> > (
                            new pdat::CellVariable<double>(d_dim, "Jameson pressure gradient", 1));
                }
                else if (variable_key == "DILATATION")
                {
                    d_Jameson_pressure_gradient =
                        boost::shared_ptr<pdat::CellVariable<double> > (
                            new pdat::CellVariable<double>(d_dim, "Jameson dilatation gradient", 1));
                }
                else if (variable_key == "ENSTROPHY")
                {
                    d_Jameson_enstrophy_gradient =
                        boost::shared_ptr<pdat::CellVariable<double> > (
                            new pdat::CellVariable<double>(d_dim, "Jameson enstrophy gradient", 1));
                }
                else if (variable_key == "MASS_FRACTION")
                {
                    d_Jameson_mass_fraction_gradient.reserve(d_num_species);
                    
                    for (int spi = 0; spi < d_num_species; spi++)
                    {
                        d_Jameson_mass_fraction_gradient.push_back(
                            boost::shared_ptr<pdat::CellVariable<double> > (
                                new pdat::CellVariable<double>(d_dim, "Jameson mass fraction " +
                                    tbox::Utilities::intToString(spi) + " gradient", 1)));
                    }
                }
                else
                {
                    TBOX_ERROR(d_object_name
                        << ": "
                        << "Unknown/unsupported variable: "
                        << variable_key
                        << "\nin input."
                        << std::endl);
                }
            }
            
            for (int vi = 0; vi < static_cast<int>(d_Jameson_gradient_variables.size()); vi++)
            {
                // Get the key of the current variable.
                std::string variable_key = d_Jameson_gradient_variables[vi];
                
                if (variable_key == "DENSITY")
                {
                    integrator->registerVariable(
                        d_Jameson_density_gradient,
                        d_num_gradient_ghosts,
                        RungeKuttaLevelIntegrator::TIME_DEP,
                            d_grid_geometry,
                            "CONSERVATIVE_COARSEN",
                            "CONSERVATIVE_LINEAR_REFINE");
                }
                else if (variable_key == "TOTAL_ENERGY")
                {
                    integrator->registerVariable(
                        d_Jameson_total_energy_gradient,
                        d_num_gradient_ghosts,
                        RungeKuttaLevelIntegrator::TIME_DEP,
                            d_grid_geometry,
                            "CONSERVATIVE_COARSEN",
                            "CONSERVATIVE_LINEAR_REFINE");
                }
                else if (variable_key == "PRESSURE")
                {
                    integrator->registerVariable(
                        d_Jameson_pressure_gradient,
                        d_num_gradient_ghosts,
                        RungeKuttaLevelIntegrator::TIME_DEP,
                            d_grid_geometry,
                            "CONSERVATIVE_COARSEN",
                            "CONSERVATIVE_LINEAR_REFINE");
                }
                else if (variable_key == "DILATATION")
                {
                    integrator->registerVariable(
                        d_Jameson_dilatation_gradient,
                        d_num_gradient_ghosts,
                        RungeKuttaLevelIntegrator::TIME_DEP,
                            d_grid_geometry,
                            "CONSERVATIVE_COARSEN",
                            "CONSERVATIVE_LINEAR_REFINE");
                }
                else if (variable_key == "ENSTROPHY")
                {
                    integrator->registerVariable(
                        d_Jameson_enstrophy_gradient,
                        d_num_gradient_ghosts,
                        RungeKuttaLevelIntegrator::TIME_DEP,
                            d_grid_geometry,
                            "CONSERVATIVE_COARSEN",
                            "CONSERVATIVE_LINEAR_REFINE");
                }
                else if (variable_key == "MASS_FRACTION")
                {
                    for (int spi = 0; spi < d_num_species; spi++)
                    {
                        integrator->registerVariable(
                            d_Jameson_mass_fraction_gradient[spi],
                            d_num_gradient_ghosts,
                            RungeKuttaLevelIntegrator::TIME_DEP,
                                d_grid_geometry,
                                "CONSERVATIVE_COARSEN",
                                "CONSERVATIVE_LINEAR_REFINE");
                    }
                }
                else
                {
                    TBOX_ERROR(d_object_name
                        << ": "
                        << "Unknown/unsupported variable: "
                        << variable_key
                        << "\nin input."
                        << std::endl);
                }
            }
        }
    }
}


/*
 * Register the plotting quantities.
 */
void
GradientTagger::registerPlotQuantities(
    const boost::shared_ptr<appu::VisItDataWriter>& visit_writer,
    const boost::shared_ptr<hier::VariableContext>& plot_context)
{
#ifdef PLOTTING_GRADIENT_TAGGER
    hier::VariableDatabase* vardb = hier::VariableDatabase::getDatabase();
    
    for (int si = 0;
             si < static_cast<int>(d_gradient_sensors.size());
             si++)
    {
        std::string sensor_key = d_gradient_sensors[si];
        
        if (sensor_key == "JAMESON_GRADIENT")
        {
            for (int vi = 0; vi < static_cast<int>(d_Jameson_gradient_variables.size()); vi++)
            {
                // Get the key of the current variable.
                std::string variable_key = d_Jameson_gradient_variables[vi];
                
                if (variable_key == "DENSITY")
                {
                    visit_writer->registerPlotQuantity(
                        "Jameson density gradient",
                        "SCALAR",
                        vardb->mapVariableAndContextToIndex(
                           d_Jameson_density_gradient,
                           plot_context));
                }
                else if (variable_key == "TOTAL_ENERGY")
                {
                    visit_writer->registerPlotQuantity(
                        "Jameson total energy gradient",
                        "SCALAR",
                        vardb->mapVariableAndContextToIndex(
                           d_Jameson_total_energy_gradient,
                           plot_context));
                }
                else if (variable_key == "PRESSURE")
                {
                    visit_writer->registerPlotQuantity(
                        "Jameson pressure gradient",
                        "SCALAR",
                        vardb->mapVariableAndContextToIndex(
                           d_Jameson_pressure_gradient,
                           plot_context));
                }
                else if (variable_key == "DILATATION")
                {
                    visit_writer->registerPlotQuantity(
                        "Jameson dilatation gradient",
                        "SCALAR",
                        vardb->mapVariableAndContextToIndex(
                           d_Jameson_dilatation_gradient,
                           plot_context));
                }
                else if (variable_key == "ENSTROPHY")
                {
                    visit_writer->registerPlotQuantity(
                        "Jameson enstrophy gradient",
                        "SCALAR",
                        vardb->mapVariableAndContextToIndex(
                           d_Jameson_enstrophy_gradient,
                           plot_context));
                }
                else if (variable_key == "MASS_FRACTION")
                {
                    for (int spi = 0; spi < d_num_species; spi++)
                    {
                        visit_writer->registerPlotQuantity(
                            "Jameson mass fraction " +
                                tbox::Utilities::intToString(spi) + " gradient",
                            "SCALAR",
                            vardb->mapVariableAndContextToIndex(
                               d_Jameson_mass_fraction_gradient[spi],
                               plot_context));
                    }
                }
            }
        }
    }
#endif
}


/*
 * Print all characteristics of the gradient tagger class.
 */
void
GradientTagger::printClassData(std::ostream& os) const
{
    os << "\nPrint GradientTagger object..."
       << std::endl;
    
    os << std::endl;
    
    os << "d_gradient_sensors = ";
    for (int si = 0; si < static_cast<int>(d_gradient_sensors.size()); si++)
    {
        os << "\"" << d_gradient_sensors[si] << "\"";
        
        if (si < static_cast<int>(d_gradient_sensors.size()) - 1)
        {
            os << ", ";
        }
    }
    os << std::endl;
    
    if (d_gradient_sensor_Jameson != nullptr)
    {
        os << std::endl;
        os << "d_Jameson_gradient_variables = ";
        for (int vi = 0; vi < static_cast<int>(d_Jameson_gradient_variables.size()); vi++)
        {
            os << "\"" << d_Jameson_gradient_variables[vi] << "\"";
            
            if (vi < static_cast<int>(d_Jameson_gradient_variables.size()) - 1)
            {
                os << ", ";
            }
        }
        os << std::endl;
        os << "d_Jameson_gradient_tol = ";
        for (int ti = 0; ti < static_cast<int>(d_Jameson_gradient_tol.size()); ti++)
        {
            os << "\"" << d_Jameson_gradient_tol[ti] << "\"";
            
            if (ti < static_cast<int>(d_Jameson_gradient_tol.size()) - 1)
            {
                os << ", ";
            }
        }
        os << std::endl;
    }
}


/*
 * Put the characteristics of the gradient tagger class into the restart
 * database.
 */
void
GradientTagger::putToRestart(
    const boost::shared_ptr<tbox::Database>& restart_db) const
{
    if (static_cast<int>(d_gradient_sensors.size()) > 0)
    {
        restart_db->putStringVector("d_gradient_sensors", d_gradient_sensors);
    }
    
    for (int si = 0; si < static_cast<int>(d_gradient_sensors.size()); si++)
    {
        if (d_gradient_sensors[si] == "JAMESON_GRADIENT")
        {
            boost::shared_ptr<tbox::Database> sensor_db =
                restart_db->putDatabase("JAMESON_GRADIENT");
            
            sensor_db->putStringVector("d_Jameson_gradient_variables",
                d_Jameson_gradient_variables);
            
            sensor_db->putDoubleVector("d_Jameson_gradient_tol",
                d_Jameson_gradient_tol);
        }
    }
}


/*
 * Tag cells for refinement using gradient sensors.
 */
void
GradientTagger::tagCells(
   hier::Patch& patch,
   boost::shared_ptr<pdat::CellData<int> > tags,
   const boost::shared_ptr<hier::VariableContext>& data_context)
{
    // Loop over gradient sensors chosen.
    for (int si = 0;
             si < static_cast<int>(d_gradient_sensors.size());
             si++)
    {
        std::string sensor_key = d_gradient_sensors[si];
        
        if (sensor_key == "JAMESON_GRADIENT")
        {
            // Looop over variables chosen.
            for (int vi = 0; vi < static_cast<int>(d_Jameson_gradient_variables.size()); vi++)
            {
                // Get the key of the current variable.
                std::string variable_key = d_Jameson_gradient_variables[vi];
                
                // Ge the tolerance for the current variable.
                double tol = d_Jameson_gradient_tol[vi];
                
                if (variable_key == "DENSITY")
                {
                    /*
                     * Register the patch and density in the flow model and compute the corresponding cell data.
                     */
                    
                    d_flow_model->registerPatchWithDataContext(patch, data_context);
                    
                    std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                    
                    num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("DENSITY", d_num_gradient_ghosts));
                    
                    d_flow_model->registerDerivedCellVariable(num_subghosts_of_data);
                    
                    d_flow_model->computeGlobalDerivedCellData();
                    
                    /*
                     * Get the pointer to density data inside the flow model.
                     */
                    
                    boost::shared_ptr<pdat::CellData<double> > data_density =
                        d_flow_model->getGlobalCellData("DENSITY");
                    
                    // Get the cell data of the density gradient.
                    boost::shared_ptr<pdat::CellData<double> > gradient(
                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                            patch.getPatchData(d_Jameson_density_gradient, data_context)));
                    
                    // Compute the gradient.
                    d_gradient_sensor_Jameson->computeGradient(patch, data_density, gradient);
                    
                    // Tag the cells.
                    tagCellsWithGradientSensor(
                        patch,
                        tags,
                        gradient,
                        tol,
                        sensor_key);
                    
                    /*
                     * Unregister the patch and data of all registered derived cell variables in the flow model.
                     */
                    
                    d_flow_model->unregisterPatch();
                    
                }
                else if (variable_key == "TOTAL_ENERGY")
                {
                    /*
                     * Register the patch and total energy in the flow model and compute the corresponding cell data.
                     */
                    
                    d_flow_model->registerPatchWithDataContext(patch, data_context);
                    
                    std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                    
                    num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("TOTAL_ENERGY", d_num_gradient_ghosts));
                    
                    d_flow_model->registerDerivedCellVariable(num_subghosts_of_data);
                    
                    d_flow_model->computeGlobalDerivedCellData();
                    
                    /*
                     * Get the pointer to total energy data inside the flow model.
                     */
                    
                    boost::shared_ptr<pdat::CellData<double> > data_total_energy =
                        d_flow_model->getGlobalCellData("TOTAL_ENERGY");
                    
                    // Get the cell data of the total energy gradient.
                    boost::shared_ptr<pdat::CellData<double> > gradient(
                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                            patch.getPatchData(d_Jameson_total_energy_gradient, data_context)));
                    
                    // Compute the gradient.
                    d_gradient_sensor_Jameson->computeGradient(patch, data_total_energy, gradient);
                    
                    // Tag the cells.
                    tagCellsWithGradientSensor(
                        patch,
                        tags,
                        gradient,
                        tol,
                        sensor_key);
                    
                    /*
                     * Unregister the patch and data of all registered derived cell variables in the flow model.
                     */
                    
                    d_flow_model->unregisterPatch();
                    
                }
                else if (variable_key == "PRESSURE")
                {
                    /*
                     * Register the patch and pressure in the flow model and compute the corresponding cell data.
                     */
                    
                    d_flow_model->registerPatchWithDataContext(patch, data_context);
                    
                    std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                    
                    num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("PRESSURE", d_num_gradient_ghosts));
                    
                    d_flow_model->registerDerivedCellVariable(num_subghosts_of_data);
                    
                    d_flow_model->computeGlobalDerivedCellData();
                    
                    /*
                     * Get the pointer to pressure data inside the flow model.
                     */
                    
                    boost::shared_ptr<pdat::CellData<double> > data_pressure =
                        d_flow_model->getGlobalCellData("PRESSURE");
                    
                    // Get the cell data of the pressure gradient.
                    boost::shared_ptr<pdat::CellData<double> > gradient(
                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                            patch.getPatchData(d_Jameson_pressure_gradient, data_context)));
                    
                    // Compute the gradient.
                    d_gradient_sensor_Jameson->computeGradient(patch, data_pressure, gradient);
                    
                    // Tag the cells.
                    tagCellsWithGradientSensor(
                        patch,
                        tags,
                        gradient,
                        tol,
                        sensor_key);
                    
                    /*
                     * Unregister the patch and data of all registered derived cell variables in the flow model.
                     */
                    
                    d_flow_model->unregisterPatch();
                    
                }
                else if (variable_key == "DILATATION")
                {
                    /*
                     * Register the patch and dilatation in the flow model and compute the corresponding cell data.
                     */
                    
                    d_flow_model->registerPatchWithDataContext(patch, data_context);
                    
                    std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                    
                    num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("DILATATION", d_num_gradient_ghosts));
                    
                    d_flow_model->registerDerivedCellVariable(num_subghosts_of_data);
                    
                    d_flow_model->computeGlobalDerivedCellData();
                    
                    /*
                     * Get the pointer to dilatation data inside the flow model.
                     */
                    
                    boost::shared_ptr<pdat::CellData<double> > data_dilatation =
                        d_flow_model->getGlobalCellData("DILATATION");
                    
                    // Get the cell data of the dilatation gradient.
                    boost::shared_ptr<pdat::CellData<double> > gradient(
                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                            patch.getPatchData(d_Jameson_dilatation_gradient, data_context)));
                    
                    // Compute the gradient.
                    d_gradient_sensor_Jameson->computeGradient(patch, data_dilatation, gradient);
                    
                    // Tag the cells.
                    tagCellsWithGradientSensor(
                        patch,
                        tags,
                        gradient,
                        tol,
                        sensor_key);
                    
                    /*
                     * Unregister the patch and data of all registered derived cell variables in the flow model.
                     */
                    
                    d_flow_model->unregisterPatch();
                    
                }
                else if (variable_key == "ENSTROPHY")
                {
                    /*
                     * Register the patch and enstrophy in the flow model and compute the corresponding cell data.
                     */
                    
                    d_flow_model->registerPatchWithDataContext(patch, data_context);
                    
                    std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                    
                    num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("ENSTROPHY", d_num_gradient_ghosts));
                    
                    d_flow_model->registerDerivedCellVariable(num_subghosts_of_data);
                    
                    d_flow_model->computeGlobalDerivedCellData();
                    
                    /*
                     * Get the pointer to enstrophy data inside the flow model.
                     */
                    
                    boost::shared_ptr<pdat::CellData<double> > data_enstrophy =
                        d_flow_model->getGlobalCellData("ENSTROPHY");
                    
                    // Get the cell data of the enstrophy gradient.
                    boost::shared_ptr<pdat::CellData<double> > gradient(
                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                            patch.getPatchData(d_Jameson_enstrophy_gradient, data_context)));
                    
                    // Compute the gradient.
                    d_gradient_sensor_Jameson->computeGradient(patch, data_enstrophy, gradient);
                    
                    // Tag the cells.
                    tagCellsWithGradientSensor(
                        patch,
                        tags,
                        gradient,
                        tol,
                        sensor_key);
                    
                    /*
                     * Unregister the patch and data of all registered derived cell variables in the flow model.
                     */
                    
                    d_flow_model->unregisterPatch();
                    
                }
                else if (variable_key == "MASS_FRACTION")
                {
                    /*
                     * Register the patch and dilatation in the flow model and compute the corresponding cell data.
                     */
                    
                    d_flow_model->registerPatchWithDataContext(patch, data_context);
                    
                    std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                    
                    num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("MASS_FRACTION", d_num_gradient_ghosts));
                    
                    d_flow_model->registerDerivedCellVariable(num_subghosts_of_data);
                    
                    d_flow_model->computeGlobalDerivedCellData();
                    
                    /*
                     * Get the pointer to mass fraction data inside the flow model.
                     */
                    
                    boost::shared_ptr<pdat::CellData<double> > data_mass_fraction =
                        d_flow_model->getGlobalCellData("MASS_FRACTION");
                    
                    for (int spi = 0; spi < d_num_species; spi++)
                    {
                        // Get the cell data of the dilatation gradient.
                        boost::shared_ptr<pdat::CellData<double> > gradient(
                            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                patch.getPatchData(d_Jameson_mass_fraction_gradient[spi], data_context)));
                        
                        // Compute the gradient.
                        d_gradient_sensor_Jameson->computeGradient(patch, data_mass_fraction, gradient, spi);
                        
                        // Tag the cells.
                        tagCellsWithGradientSensor(
                            patch,
                            tags,
                            gradient,
                            tol,
                            sensor_key);
                    }
                    
                    /*
                     * Unregister the patch and data of all registered derived cell variables in the flow model.
                     */
                    
                    d_flow_model->unregisterPatch();
                }
                else
                {
                    TBOX_ERROR(d_object_name
                        << ": "
                        << "Unknown/unsupported variable: "
                        << variable_key
                        << "\nin input."
                        << std::endl);
                }
            } // Looop over variables chosen.
        }
    } // Loop over gradient sensors chosen.
}


/*
 * Tag cells using gradient sensor.
 */
void
GradientTagger::tagCellsWithGradientSensor(
    hier::Patch& patch,
    boost::shared_ptr<pdat::CellData<int> > tags,
    boost::shared_ptr<pdat::CellData<double> > gradient,
    double& tol,
    std::string& sensor_key)
{
    TBOX_ASSERT(tags->getGhostCellWidth() == hier::IntVector::getZero(d_dim));
    
    const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(patch_geom);
#endif
    
    // Get the pointer of the tags.
    int* tag_ptr  = tags->getPointer(0);
    
    // Get the pointer to the gradient.
    double* psi = gradient->getPointer(0);
    
    // Get the number of ghost cells of gradient.
    hier::IntVector num_ghosts = gradient->getGhostCellWidth();
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = gradient->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of patch plus
    // ghost cells.
    const hier::Box ghost_box = gradient->getGhostBox();
    const hier::IntVector ghostcell_dims = ghost_box.numberCells();
    
    if (sensor_key == "JAMESON_GRADIENT")
    {
        if (d_dim == tbox::Dimension(1))
        {
            for (int i = 0; i < interior_dims[0]; i++)
            {
                // Compute indices.
                const int idx = i + num_ghosts[0];
                const int idx_nghost = i;
                
                if (psi[idx] > tol)
                {
                    tag_ptr[idx_nghost] |= 1;
                }
            }
        }
        else if (d_dim == tbox::Dimension(2))
        {
            for (int j = 0; j < interior_dims[1]; j++)
            {
                for (int i = 0; i < interior_dims[0]; i++)
                {
                    // Compute indices.
                    const int idx = (i + num_ghosts[0]) +
                        (j + num_ghosts[1])*ghostcell_dims[0];
                    
                    const int idx_nghost = i + j*interior_dims[0];
                    
                    if (psi[idx] > tol)
                    {
                        tag_ptr[idx_nghost] |= 1;
                    }
                }
            }
        }
        else if (d_dim == tbox::Dimension(3))
        {
            for (int k = 0; k < interior_dims[2]; k++)
            {
                for (int j = 0; j < interior_dims[1]; j++)
                {
                    for (int i = 0; i < interior_dims[0]; i++)
                    {
                        // Compute indices.
                        const int idx = (i + num_ghosts[0]) +
                            (j + num_ghosts[1])*ghostcell_dims[0] +
                            (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                        
                        const int idx_nghost = i + j*interior_dims[0] + k*interior_dims[0]*interior_dims[1];
                        
                        if (psi[idx] > tol)
                        {
                            tag_ptr[idx_nghost] |= 1;
                        }
                    }
                }
            }
        }
    }
}
