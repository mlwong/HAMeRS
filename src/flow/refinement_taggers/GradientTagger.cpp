#include "flow/refinement_taggers/GradientTagger.hpp"

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
                if (!((sensor_key == "FIRST_DERIVATIVE") ||
                      (sensor_key == "SECOND_DERIVATIVE") ||
                      (sensor_key == "JAMESON_GRADIENT")))
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
                
                if (sensor_db && sensor_key == "FIRST_DERIVATIVE")
                {
                    d_gradient_sensor_first_derivative.reset(new GradientSensorFirstDerivative(
                        "first derivative gradient sensor",
                        d_dim));
                    
                    if (sensor_db->keyExists("first_derivative_variables"))
                    {
                        d_first_derivative_variables = sensor_db->getStringVector("first_derivative_variables");
                    }
                    else if (sensor_db->keyExists("d_first_derivative_variables"))
                    {
                        d_first_derivative_variables = sensor_db->getStringVector("d_first_derivative_variables");
                    }
                    else
                    {
                        TBOX_ERROR(d_object_name
                            << ": "
                            << "No key 'first_derivative_variables'/'d_first_derivative_variables' found in database for "
                            << sensor_key
                            << "."
                            << std::endl);
                    }
                    
                    for (int vi = 0; vi < static_cast<int>(d_first_derivative_variables.size()); vi++)
                    {
                        std::string variable_key = d_first_derivative_variables[vi];
                        
                        if (!((variable_key == "DENSITY") ||
                              (variable_key == "TOTAL_ENERGY") ||
                              (variable_key == "PRESSURE")))
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "Unknown/unsupported variable: "
                                << variable_key
                                << "\nin input."
                                << std::endl);
                        }
                    }
                    
                    int num_true = 0;
                    
                    /*
                     * Get the settings for derivative global tolerances.
                     */
                    
                    if (sensor_db->keyExists("first_derivative_uses_global_tol"))
                    {
                        d_first_derivative_uses_global_tol = sensor_db->getBoolVector("first_derivative_uses_global_tol");
                    }
                    else if (sensor_db->keyExists("d_first_derivative_uses_global_tol"))
                    {
                        d_first_derivative_uses_global_tol = sensor_db->getBoolVector("d_first_derivative_uses_global_tol");
                    }
                    else
                    {
                        TBOX_ERROR(d_object_name
                            << ": "
                            << "No key 'first_derivative_uses_global_tol'/'d_first_derivative_uses_global_tol'"
                            << " found in database for "
                            << sensor_key
                            << "."
                            << std::endl);
                    }
                    
                    if (static_cast<int>(d_first_derivative_variables.size()) !=
                        static_cast<int>(d_first_derivative_uses_global_tol.size()))
                    {
                        TBOX_ERROR(d_object_name
                            << ": "
                            << "The numbers of variables and switches for global tolerances provided don't match"
                            << " in database for "
                            << sensor_key
                            << "."
                            << std::endl);
                    }
                    
                    num_true = std::count(d_first_derivative_uses_global_tol.begin(),
                        d_first_derivative_uses_global_tol.end(),
                        true);
                    
                    if (num_true > 0)
                    {
                        if (sensor_db->keyExists("first_derivative_global_tol"))
                        {
                            d_first_derivative_global_tol = sensor_db->getDoubleVector("first_derivative_global_tol");
                        }
                        else if (sensor_db->keyExists("d_first_derivative_global_tol"))
                        {
                            d_first_derivative_global_tol = sensor_db->getDoubleVector("d_first_derivative_global_tol");
                        }
                        else
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "No key 'first_derivative_global_tol'/'d_first_derivative_global_tol' found in database for "
                                << sensor_key
                                << "."
                                << std::endl);
                        }
                        
                        if (num_true != static_cast<int>(d_first_derivative_global_tol.size()))
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "The number of variables that use global tolerances and number of"
                                << " global tolerances provided don't match"
                                << " in database for "
                                << sensor_key
                                << "."
                                << std::endl);
                        }
                    }
                    
                    /*
                     * Get the settings for derivative local tolerances.
                     */
                    
                    if (sensor_db->keyExists("first_derivative_uses_local_tol"))
                    {
                        d_first_derivative_uses_local_tol = sensor_db->getBoolVector("first_derivative_uses_local_tol");
                    }
                    else if (sensor_db->keyExists("d_first_derivative_uses_local_tol"))
                    {
                        d_first_derivative_uses_local_tol = sensor_db->getBoolVector("d_first_derivative_uses_local_tol");
                    }
                    else
                    {
                        TBOX_ERROR(d_object_name
                            << ": "
                            << "No key 'first_derivative_uses_local_tol'/'d_first_derivative_uses_local_tol'"
                            << " found in database for "
                            << sensor_key
                            << "."
                            << std::endl);
                    }
                    
                    if (static_cast<int>(d_first_derivative_variables.size()) !=
                        static_cast<int>(d_first_derivative_uses_local_tol.size()))
                    {
                        TBOX_ERROR(d_object_name
                            << ": "
                            << "The numbers of variables and switches for local tolerances provided don't match"
                            << " in database for "
                            << sensor_key
                            << "."
                            << std::endl);
                    }
                    
                    num_true = std::count(d_first_derivative_uses_local_tol.begin(),
                        d_first_derivative_uses_local_tol.end(),
                        true);
                    
                    if (num_true > 0)
                    {
                        if (sensor_db->keyExists("first_derivative_local_tol"))
                        {
                            d_first_derivative_local_tol = sensor_db->getDoubleVector("first_derivative_local_tol");
                        }
                        else if (sensor_db->keyExists("d_first_derivative_local_tol"))
                        {
                            d_first_derivative_local_tol = sensor_db->getDoubleVector("d_first_derivative_local_tol");
                        }
                        else
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "No key 'first_derivative_local_tol'/'d_first_derivative_local_tol' found in database for "
                                << sensor_key
                                << "."
                                << std::endl);
                        }
                        
                        if (num_true != static_cast<int>(d_first_derivative_local_tol.size()))
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "The number of variables that use local tolerances and number of"
                                << " local tolerances provided don't match"
                                << " in database for "
                                << sensor_key
                                << "."
                                << std::endl);
                        }
                    }
                }
                else if (sensor_db && sensor_key == "SECOND_DERIVATIVE")
                {
                    d_gradient_sensor_second_derivative.reset(new GradientSensorSecondDerivative(
                        "second derivative gradient sensor",
                        d_dim));
                    
                    if (sensor_db->keyExists("second_derivative_variables"))
                    {
                        d_second_derivative_variables = sensor_db->getStringVector("second_derivative_variables");
                    }
                    else if (sensor_db->keyExists("d_second_derivative_variables"))
                    {
                        d_second_derivative_variables = sensor_db->getStringVector("d_second_derivative_variables");
                    }
                    else
                    {
                        TBOX_ERROR(d_object_name
                            << ": "
                            << "No key 'second_derivative_variables'/'d_second_derivative_variables' found in database for "
                            << sensor_key
                            << "."
                            << std::endl);
                    }
                    
                    for (int vi = 0; vi < static_cast<int>(d_second_derivative_variables.size()); vi++)
                    {
                        std::string variable_key = d_second_derivative_variables[vi];
                        
                        if (!((variable_key == "DENSITY") ||
                              (variable_key == "TOTAL_ENERGY") ||
                              (variable_key == "PRESSURE")))
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "Unknown/unsupported variable: "
                                << variable_key
                                << "\nin input."
                                << std::endl);
                        }
                    }
                    
                    int num_true = 0;
                    
                    /*
                     * Get the settings for derivative global tolerances.
                     */
                    
                    if (sensor_db->keyExists("second_derivative_uses_global_tol"))
                    {
                        d_second_derivative_uses_global_tol = sensor_db->getBoolVector("second_derivative_uses_global_tol");
                    }
                    else if (sensor_db->keyExists("d_second_derivative_uses_global_tol"))
                    {
                        d_second_derivative_uses_global_tol = sensor_db->getBoolVector("d_second_derivative_uses_global_tol");
                    }
                    else
                    {
                        TBOX_ERROR(d_object_name
                            << ": "
                            << "No key 'second_derivative_uses_global_tol'/'d_second_derivative_uses_global_tol'"
                            << " found in database for "
                            << sensor_key
                            << "."
                            << std::endl);
                    }
                    
                    if (static_cast<int>(d_second_derivative_variables.size()) !=
                        static_cast<int>(d_second_derivative_uses_global_tol.size()))
                    {
                        TBOX_ERROR(d_object_name
                            << ": "
                            << "The numbers of variables and switches for global tolerances provided don't match"
                            << " in database for "
                            << sensor_key
                            << "."
                            << std::endl);
                    }
                    
                    num_true = std::count(d_second_derivative_uses_global_tol.begin(),
                        d_second_derivative_uses_global_tol.end(),
                        true);
                    
                    if (num_true > 0)
                    {
                        if (sensor_db->keyExists("second_derivative_global_tol"))
                        {
                            d_second_derivative_global_tol = sensor_db->getDoubleVector("second_derivative_global_tol");
                        }
                        else if (sensor_db->keyExists("d_second_derivative_global_tol"))
                        {
                            d_second_derivative_global_tol = sensor_db->getDoubleVector("d_second_derivative_global_tol");
                        }
                        else
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "No key 'second_derivative_global_tol'/'d_second_derivative_global_tol' found in database for "
                                << sensor_key
                                << "."
                                << std::endl);
                        }
                        
                        if (num_true != static_cast<int>(d_second_derivative_global_tol.size()))
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "The number of variables that use global tolerances and number of"
                                << " global tolerances provided don't match"
                                << " in database for "
                                << sensor_key
                                << "."
                                << std::endl);
                        }
                    }
                    
                    /*
                     * Get the settings for derivative local tolerances.
                     */
                    
                    if (sensor_db->keyExists("second_derivative_uses_local_tol"))
                    {
                        d_second_derivative_uses_local_tol = sensor_db->getBoolVector("second_derivative_uses_local_tol");
                    }
                    else if (sensor_db->keyExists("d_second_derivative_uses_local_tol"))
                    {
                        d_second_derivative_uses_local_tol = sensor_db->getBoolVector("d_second_derivative_uses_local_tol");
                    }
                    else
                    {
                        TBOX_ERROR(d_object_name
                            << ": "
                            << "No key 'second_derivative_uses_local_tol'/'d_second_derivative_uses_local_tol'"
                            << " found in database for "
                            << sensor_key
                            << "."
                            << std::endl);
                    }
                    
                    if (static_cast<int>(d_second_derivative_variables.size()) !=
                        static_cast<int>(d_second_derivative_uses_local_tol.size()))
                    {
                        TBOX_ERROR(d_object_name
                            << ": "
                            << "The numbers of variables and switches for local tolerances provided don't match"
                            << " in database for "
                            << sensor_key
                            << "."
                            << std::endl);
                    }
                    
                    num_true = std::count(d_second_derivative_uses_local_tol.begin(),
                        d_second_derivative_uses_local_tol.end(),
                        true);
                    
                    if (num_true > 0)
                    {
                        if (sensor_db->keyExists("second_derivative_local_tol"))
                        {
                            d_second_derivative_local_tol = sensor_db->getDoubleVector("second_derivative_local_tol");
                        }
                        else if (sensor_db->keyExists("d_second_derivative_local_tol"))
                        {
                            d_second_derivative_local_tol = sensor_db->getDoubleVector("d_second_derivative_local_tol");
                        }
                        else
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "No key 'second_derivative_local_tol'/'d_second_derivative_local_tol' found in database for "
                                << sensor_key
                                << "."
                                << std::endl);
                        }
                        
                        if (num_true != static_cast<int>(d_second_derivative_local_tol.size()))
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "The number of variables that use local tolerances and number of"
                                << " local tolerances provided don't match"
                                << " in database for "
                                << sensor_key
                                << "."
                                << std::endl);
                        }
                    }
                }
                else if (sensor_db && sensor_key == "JAMESON_GRADIENT")
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
                              (variable_key == "PRESSURE")))
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
        
        if (sensor_key == "FIRST_DERIVATIVE")
        {
            d_num_gradient_ghosts = hier::IntVector::max(
                d_num_gradient_ghosts,
                d_gradient_sensor_first_derivative->getGradientSensorNumberOfGhostCells());
        }
        else if (sensor_key == "SECOND_DERIVATIVE")
        {
            d_num_gradient_ghosts = hier::IntVector::max(
                d_num_gradient_ghosts,
                d_gradient_sensor_second_derivative->getGradientSensorNumberOfGhostCells());
        }
        else if (sensor_key == "JAMESON_GRADIENT")
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
        
        if (sensor_key == "FIRST_DERIVATIVE")
        {
            for (int vi = 0; vi < static_cast<int>(d_first_derivative_variables.size()); vi++)
            {
                // Get the key of the current variable.
                std::string variable_key = d_first_derivative_variables[vi];
                
                if (variable_key == "DENSITY")
                {
                    d_first_derivative_density =
                        boost::shared_ptr<pdat::CellVariable<double> > (
                            new pdat::CellVariable<double>(d_dim, "first derivative of density", 1));
                    
                    if (d_first_derivative_uses_global_tol[vi])
                    {
                        d_first_derivative_max_density = 0.0;
                    }
                    
                    if (d_first_derivative_uses_local_tol[vi])
                    {
                        d_first_derivative_local_mean_density =
                            boost::shared_ptr<pdat::CellVariable<double> > (
                                new pdat::CellVariable<double>(d_dim, "first derivative local mean of density", 1));
                    }
                }
                else if (variable_key == "TOTAL_ENERGY")
                {
                    d_first_derivative_total_energy =
                        boost::shared_ptr<pdat::CellVariable<double> > (
                            new pdat::CellVariable<double>(d_dim, "first derivative of total energy", 1));
                    
                    if (d_first_derivative_uses_global_tol[vi])
                    {
                        d_first_derivative_max_total_energy = 0.0;
                    }
                    
                    if (d_first_derivative_uses_local_tol[vi])
                    {
                        d_first_derivative_local_mean_total_energy =
                            boost::shared_ptr<pdat::CellVariable<double> > (
                                new pdat::CellVariable<double>(d_dim, "first derivative local mean of total energy", 1));
                    }
                }
                else if (variable_key == "PRESSURE")
                {
                    d_first_derivative_pressure =
                        boost::shared_ptr<pdat::CellVariable<double> > (
                            new pdat::CellVariable<double>(d_dim, "first derivative of pressure", 1));
                    
                    if (d_first_derivative_uses_global_tol[vi])
                    {
                        d_first_derivative_max_pressure = 0.0;
                    }
                    
                    if (d_first_derivative_uses_local_tol[vi])
                    {
                        d_first_derivative_local_mean_pressure =
                            boost::shared_ptr<pdat::CellVariable<double> > (
                                new pdat::CellVariable<double>(d_dim, "first derivative local mean of pressure", 1));
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
            
            for (int vi = 0; vi < static_cast<int>(d_first_derivative_variables.size()); vi++)
            {
                // Get the key of the current variable.
                std::string variable_key = d_first_derivative_variables[vi];
                
                if (variable_key == "DENSITY")
                {
                    integrator->registerVariable(
                        d_first_derivative_density,
                        d_num_gradient_ghosts,
                        RungeKuttaLevelIntegrator::TIME_DEP,
                            d_grid_geometry,
                            "NO_COARSEN",
                            "NO_REFINE");
                }
                else if (variable_key == "TOTAL_ENERGY")
                {
                    integrator->registerVariable(
                        d_first_derivative_total_energy,
                        d_num_gradient_ghosts,
                        RungeKuttaLevelIntegrator::TIME_DEP,
                            d_grid_geometry,
                            "NO_COARSEN",
                            "NO_REFINE");
                }
                else if (variable_key == "PRESSURE")
                {
                    integrator->registerVariable(
                        d_first_derivative_pressure,
                        d_num_gradient_ghosts,
                        RungeKuttaLevelIntegrator::TIME_DEP,
                            d_grid_geometry,
                            "NO_COARSEN",
                            "NO_REFINE");
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
        else if (sensor_key == "SECOND_DERIVATIVE")
        {
            for (int vi = 0; vi < static_cast<int>(d_second_derivative_variables.size()); vi++)
            {
                // Get the key of the current variable.
                std::string variable_key = d_second_derivative_variables[vi];
                
                if (variable_key == "DENSITY")
                {
                    d_second_derivative_density =
                        boost::shared_ptr<pdat::CellVariable<double> > (
                            new pdat::CellVariable<double>(d_dim, "second derivative of density", 1));
                    
                    if (d_second_derivative_uses_global_tol[vi])
                    {
                        d_second_derivative_max_density = 0.0;
                    }
                    
                    if (d_second_derivative_uses_local_tol[vi])
                    {
                        d_second_derivative_local_mean_density =
                            boost::shared_ptr<pdat::CellVariable<double> > (
                                new pdat::CellVariable<double>(d_dim, "second derivative local mean of density", 1));
                    }
                }
                else if (variable_key == "TOTAL_ENERGY")
                {
                    d_second_derivative_total_energy =
                        boost::shared_ptr<pdat::CellVariable<double> > (
                            new pdat::CellVariable<double>(d_dim, "second derivative of total energy", 1));
                    
                    if (d_second_derivative_uses_global_tol[vi])
                    {
                        d_second_derivative_max_total_energy = 0.0;
                    }
                    
                    if (d_second_derivative_uses_local_tol[vi])
                    {
                        d_second_derivative_local_mean_total_energy =
                            boost::shared_ptr<pdat::CellVariable<double> > (
                                new pdat::CellVariable<double>(d_dim, "second derivative local mean of total energy", 1));
                    }
                }
                else if (variable_key == "PRESSURE")
                {
                    d_second_derivative_pressure =
                        boost::shared_ptr<pdat::CellVariable<double> > (
                            new pdat::CellVariable<double>(d_dim, "second derivative of pressure", 1));
                    
                    if (d_second_derivative_uses_global_tol[vi])
                    {
                        d_second_derivative_max_pressure = 0.0;
                    }
                    
                    if (d_second_derivative_uses_local_tol[vi])
                    {
                        d_second_derivative_local_mean_pressure =
                            boost::shared_ptr<pdat::CellVariable<double> > (
                                new pdat::CellVariable<double>(d_dim, "second derivative local mean of pressure", 1));
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
            
            for (int vi = 0; vi < static_cast<int>(d_second_derivative_variables.size()); vi++)
            {
                // Get the key of the current variable.
                std::string variable_key = d_second_derivative_variables[vi];
                
                if (variable_key == "DENSITY")
                {
                    integrator->registerVariable(
                        d_second_derivative_density,
                        d_num_gradient_ghosts,
                        RungeKuttaLevelIntegrator::TIME_DEP,
                            d_grid_geometry,
                            "NO_COARSEN",
                            "NO_REFINE");
                }
                else if (variable_key == "TOTAL_ENERGY")
                {
                    integrator->registerVariable(
                        d_second_derivative_total_energy,
                        d_num_gradient_ghosts,
                        RungeKuttaLevelIntegrator::TIME_DEP,
                            d_grid_geometry,
                            "NO_COARSEN",
                            "NO_REFINE");
                }
                else if (variable_key == "PRESSURE")
                {
                    integrator->registerVariable(
                        d_second_derivative_pressure,
                        d_num_gradient_ghosts,
                        RungeKuttaLevelIntegrator::TIME_DEP,
                            d_grid_geometry,
                            "NO_COARSEN",
                            "NO_REFINE");
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
        else if (sensor_key == "JAMESON_GRADIENT")
        {
            for (int vi = 0; vi < static_cast<int>(d_Jameson_gradient_variables.size()); vi++)
            {
                // Get the key of the current variable.
                std::string variable_key = d_Jameson_gradient_variables[vi];
                
                if (variable_key == "DENSITY")
                {
                    d_Jameson_gradient_density =
                        boost::shared_ptr<pdat::CellVariable<double> > (
                            new pdat::CellVariable<double>(d_dim, "Jameson density gradient", 1));
                }
                else if (variable_key == "TOTAL_ENERGY")
                {
                    d_Jameson_gradient_total_energy =
                        boost::shared_ptr<pdat::CellVariable<double> > (
                            new pdat::CellVariable<double>(d_dim, "Jameson total_energy gradient", 1));
                }
                else if (variable_key == "PRESSURE")
                {
                    d_Jameson_gradient_pressure =
                        boost::shared_ptr<pdat::CellVariable<double> > (
                            new pdat::CellVariable<double>(d_dim, "Jameson pressure gradient", 1));
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
                        d_Jameson_gradient_density,
                        d_num_gradient_ghosts,
                        RungeKuttaLevelIntegrator::TIME_DEP,
                            d_grid_geometry,
                            "NO_COARSEN",
                            "NO_REFINE");
                }
                else if (variable_key == "TOTAL_ENERGY")
                {
                    integrator->registerVariable(
                        d_Jameson_gradient_total_energy,
                        d_num_gradient_ghosts,
                        RungeKuttaLevelIntegrator::TIME_DEP,
                            d_grid_geometry,
                            "NO_COARSEN",
                            "NO_REFINE");
                }
                else if (variable_key == "PRESSURE")
                {
                    integrator->registerVariable(
                        d_Jameson_gradient_pressure,
                        d_num_gradient_ghosts,
                        RungeKuttaLevelIntegrator::TIME_DEP,
                            d_grid_geometry,
                            "NO_COARSEN",
                            "NO_REFINE");
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
        
        if (sensor_key == "FIRST_DERIVATIVE")
        {
            for (int vi = 0; vi < static_cast<int>(d_first_derivative_variables.size()); vi++)
            {
                // Get the key of the current variable.
                std::string variable_key = d_first_derivative_variables[vi];
                
                if (variable_key == "DENSITY")
                {
                    visit_writer->registerPlotQuantity(
                        "first derivative of density",
                        "SCALAR",
                        vardb->mapVariableAndContextToIndex(
                           d_first_derivative_density,
                           plot_context));
                }
                else if (variable_key == "TOTAL_ENERGY")
                {
                    visit_writer->registerPlotQuantity(
                        "first derivative of total energy",
                        "SCALAR",
                        vardb->mapVariableAndContextToIndex(
                           d_first_derivative_total_energy,
                           plot_context));
                }
                else if (variable_key == "PRESSURE")
                {
                    visit_writer->registerPlotQuantity(
                        "first derivative of pressure",
                        "SCALAR",
                        vardb->mapVariableAndContextToIndex(
                           d_first_derivative_pressure,
                           plot_context));
                }
            }
        }
        else if (sensor_key == "SECOND_DERIVATIVE")
        {
            for (int vi = 0; vi < static_cast<int>(d_second_derivative_variables.size()); vi++)
            {
                // Get the key of the current variable.
                std::string variable_key = d_second_derivative_variables[vi];
                
                if (variable_key == "DENSITY")
                {
                    visit_writer->registerPlotQuantity(
                        "second derivative of density",
                        "SCALAR",
                        vardb->mapVariableAndContextToIndex(
                           d_second_derivative_density,
                           plot_context));
                }
                else if (variable_key == "TOTAL_ENERGY")
                {
                    visit_writer->registerPlotQuantity(
                        "second derivative of total energy",
                        "SCALAR",
                        vardb->mapVariableAndContextToIndex(
                           d_second_derivative_total_energy,
                           plot_context));
                }
                else if (variable_key == "PRESSURE")
                {
                    visit_writer->registerPlotQuantity(
                        "second derivative of pressure",
                        "SCALAR",
                        vardb->mapVariableAndContextToIndex(
                           d_second_derivative_pressure,
                           plot_context));
                }
            }
        }
        else if (sensor_key == "JAMESON_GRADIENT")
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
                           d_Jameson_gradient_density,
                           plot_context));
                }
                else if (variable_key == "TOTAL_ENERGY")
                {
                    visit_writer->registerPlotQuantity(
                        "Jameson total energy gradient",
                        "SCALAR",
                        vardb->mapVariableAndContextToIndex(
                           d_Jameson_gradient_total_energy,
                           plot_context));
                }
                else if (variable_key == "PRESSURE")
                {
                    visit_writer->registerPlotQuantity(
                        "Jameson pressure gradient",
                        "SCALAR",
                        vardb->mapVariableAndContextToIndex(
                           d_Jameson_gradient_pressure,
                           plot_context));
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
    
    if (d_gradient_sensor_first_derivative != nullptr)
    {
        os << std::endl;
        os << "d_first_derivative_variables = ";
        for (int vi = 0; vi < static_cast<int>(d_first_derivative_variables.size()); vi++)
        {
            os << "\"" << d_first_derivative_variables[vi] << "\"";
            
            if (vi < static_cast<int>(d_first_derivative_variables.size()) - 1)
            {
                os << ", ";
            }
        }
        os << std::endl;
        
        os << "d_first_derivative_uses_global_tol = ";
        for (int ti = 0; ti < static_cast<int>(d_first_derivative_uses_global_tol.size()); ti++)
        {
            os << "\"" << d_first_derivative_uses_global_tol[ti] << "\"";
            
            if (ti < static_cast<int>(d_first_derivative_uses_global_tol.size()) - 1)
            {
                os << ", ";
            }
        }
        os << std::endl;
        
        os << "d_first_derivative_global_tol = ";
        for (int ti = 0; ti < static_cast<int>(d_first_derivative_global_tol.size()); ti++)
        {
            os << "\"" << d_first_derivative_global_tol[ti] << "\"";
            
            if (ti < static_cast<int>(d_first_derivative_global_tol.size()) - 1)
            {
                os << ", ";
            }
        }
        os << std::endl;
        
        os << "d_first_derivative_uses_local_tol = ";
        for (int ti = 0; ti < static_cast<int>(d_first_derivative_uses_local_tol.size()); ti++)
        {
            os << "\"" << d_first_derivative_uses_local_tol[ti] << "\"";
            
            if (ti < static_cast<int>(d_first_derivative_uses_local_tol.size()) - 1)
            {
                os << ", ";
            }
        }
        os << std::endl;
        
        os << "d_first_derivative_local_tol = ";
        for (int ti = 0; ti < static_cast<int>(d_first_derivative_local_tol.size()); ti++)
        {
            os << "\"" << d_first_derivative_local_tol[ti] << "\"";
            
            if (ti < static_cast<int>(d_first_derivative_local_tol.size()) - 1)
            {
                os << ", ";
            }
        }
        os << std::endl;
    }
    
    if (d_gradient_sensor_second_derivative != nullptr)
    {
        os << std::endl;
        os << "d_second_derivative_variables = ";
        for (int vi = 0; vi < static_cast<int>(d_second_derivative_variables.size()); vi++)
        {
            os << "\"" << d_second_derivative_variables[vi] << "\"";
            
            if (vi < static_cast<int>(d_second_derivative_variables.size()) - 1)
            {
                os << ", ";
            }
        }
        os << std::endl;
        
        os << "d_second_derivative_uses_global_tol = ";
        for (int ti = 0; ti < static_cast<int>(d_second_derivative_uses_global_tol.size()); ti++)
        {
            os << "\"" << d_second_derivative_uses_global_tol[ti] << "\"";
            
            if (ti < static_cast<int>(d_second_derivative_uses_global_tol.size()) - 1)
            {
                os << ", ";
            }
        }
        os << std::endl;
        
        os << "d_second_derivative_global_tol = ";
        for (int ti = 0; ti < static_cast<int>(d_second_derivative_global_tol.size()); ti++)
        {
            os << "\"" << d_second_derivative_global_tol[ti] << "\"";
            
            if (ti < static_cast<int>(d_second_derivative_global_tol.size()) - 1)
            {
                os << ", ";
            }
        }
        os << std::endl;
        
        os << "d_second_derivative_uses_local_tol = ";
        for (int ti = 0; ti < static_cast<int>(d_second_derivative_uses_local_tol.size()); ti++)
        {
            os << "\"" << d_second_derivative_uses_local_tol[ti] << "\"";
            
            if (ti < static_cast<int>(d_second_derivative_uses_local_tol.size()) - 1)
            {
                os << ", ";
            }
        }
        os << std::endl;
        
        os << "d_second_derivative_local_tol = ";
        for (int ti = 0; ti < static_cast<int>(d_second_derivative_local_tol.size()); ti++)
        {
            os << "\"" << d_second_derivative_local_tol[ti] << "\"";
            
            if (ti < static_cast<int>(d_second_derivative_local_tol.size()) - 1)
            {
                os << ", ";
            }
        }
        os << std::endl;
    }
    
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
        if (d_gradient_sensors[si] == "FIRST_DERIVATIVE")
        {
            boost::shared_ptr<tbox::Database> sensor_db =
                restart_db->putDatabase("FIRST_DERIVATIVE");
            
            sensor_db->putStringVector("d_first_derivative_variables",
                d_first_derivative_variables);
            
            sensor_db->putBoolVector("d_first_derivative_uses_global_tol",
                d_first_derivative_uses_global_tol);
            
            sensor_db->putBoolVector("d_first_derivative_uses_local_tol",
                d_first_derivative_uses_local_tol);
            
            int num_true = 0;
            
            num_true = std::count(d_first_derivative_uses_global_tol.begin(),
                d_first_derivative_uses_global_tol.end(),
                true);
            if (num_true > 0)
            {
                sensor_db->putDoubleVector("d_first_derivative_global_tol",
                    d_first_derivative_global_tol);
            }
            
            num_true = std::count(d_first_derivative_uses_local_tol.begin(),
                d_first_derivative_uses_local_tol.end(),
                true);
            if (num_true > 0)
            {
                sensor_db->putDoubleVector("d_first_derivative_local_tol",
                    d_first_derivative_local_tol);
            }
        }
        
        if (d_gradient_sensors[si] == "SECOND_DERIVATIVE")
        {
            boost::shared_ptr<tbox::Database> sensor_db =
                restart_db->putDatabase("SECOND_DERIVATIVE");
            
            sensor_db->putStringVector("d_second_derivative_variables",
                d_second_derivative_variables);
            
            sensor_db->putBoolVector("d_second_derivative_uses_global_tol",
                d_second_derivative_uses_global_tol);
            
            sensor_db->putBoolVector("d_second_derivative_uses_local_tol",
                d_second_derivative_uses_local_tol);
            
            int num_true = 0;
            
            num_true = std::count(d_second_derivative_uses_global_tol.begin(),
                d_second_derivative_uses_global_tol.end(),
                true);
            if (num_true > 0)
            {
                sensor_db->putDoubleVector("d_second_derivative_global_tol",
                    d_second_derivative_global_tol);
            }
            
            num_true = std::count(d_second_derivative_uses_local_tol.begin(),
                d_second_derivative_uses_local_tol.end(),
                true);
            if (num_true > 0)
            {
                sensor_db->putDoubleVector("d_second_derivative_local_tol",
                    d_second_derivative_local_tol);
            }
        }
        
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
 * Compute values of gradient sensors.
 */
void
GradientTagger::computeGradientSensorValues(
    hier::Patch& patch,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
    // Loop over gradient sensors chosen.
    for (int si = 0;
             si < static_cast<int>(d_gradient_sensors.size());
             si++)
    {
        std::string sensor_key = d_gradient_sensors[si];
        
        if (sensor_key == "FIRST_DERIVATIVE")
        {
            for (int vi = 0; vi < static_cast<int>(d_first_derivative_variables.size()); vi++)
            {
                // Get the key of the current variable.
                std::string variable_key = d_first_derivative_variables[vi];
                
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
                    
                    // Get the pointer to density data inside the flow model.
                    boost::shared_ptr<pdat::CellData<double> > data_density =
                        d_flow_model->getGlobalCellData("DENSITY");
                    
                    // Get the cell data of the derivative.
                    boost::shared_ptr<pdat::CellData<double> > derivative(
                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                            patch.getPatchData(d_first_derivative_density, data_context)));
                    
                    // Compute the derivative.
                    if (d_first_derivative_uses_local_tol[vi])
                    {
                        // Get the cell data of local mean.
                        boost::shared_ptr<pdat::CellData<double> > variable_local_mean(
                            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                patch.getPatchData(d_first_derivative_local_mean_density, data_context)));
                        
                        d_gradient_sensor_first_derivative->computeGradientWithVariableLocalMean(
                            patch,
                            data_density,
                            derivative,
                            variable_local_mean);
                    }
                    else
                    {
                        d_gradient_sensor_first_derivative->computeGradient(
                            patch,
                            data_density,
                            derivative);
                    }
                    
                    /*
                     * Unregister the patch and data of all registered derived cell variables in the flow model.
                     */
                    
                    d_flow_model->unregisterPatch();
                    
                }
                else if (variable_key ==  "TOTAL_ENERGY")
                {
                    /*
                     * Register the patch and total energy in the flow model and compute the corresponding cell data.
                     */
                    
                    d_flow_model->registerPatchWithDataContext(patch, data_context);
                    
                    std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                    
                    num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("TOTAL_ENERGY", d_num_gradient_ghosts));
                    
                    d_flow_model->registerDerivedCellVariable(num_subghosts_of_data);
                    
                    d_flow_model->computeGlobalDerivedCellData();
                    
                    // Get the pointer to total energy data inside the flow model.
                    boost::shared_ptr<pdat::CellData<double> > data_total_energy =
                        d_flow_model->getGlobalCellData("TOTAL_ENERGY");
                    
                    // Get the cell data of the derivative.
                    boost::shared_ptr<pdat::CellData<double> > derivative(
                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                            patch.getPatchData(d_first_derivative_total_energy, data_context)));
                    
                    // Compute the derivative.
                    if (d_first_derivative_uses_local_tol[vi])
                    {
                        // Get the cell data of local mean.
                        boost::shared_ptr<pdat::CellData<double> > variable_local_mean(
                            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                patch.getPatchData(d_first_derivative_local_mean_total_energy, data_context)));
                        
                        d_gradient_sensor_first_derivative->computeGradientWithVariableLocalMean(
                            patch,
                            data_total_energy,
                            derivative,
                            variable_local_mean);
                    }
                    else
                    {
                        d_gradient_sensor_first_derivative->computeGradient(
                            patch,
                            data_total_energy,
                            derivative);
                    }
                    
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
                    
                    // Get the pointer to pressure data inside the flow model.
                    boost::shared_ptr<pdat::CellData<double> > data_pressure =
                        d_flow_model->getGlobalCellData("PRESSURE");
                    
                    // Get the cell data of the derivative.
                    boost::shared_ptr<pdat::CellData<double> > derivative(
                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                            patch.getPatchData(d_first_derivative_pressure, data_context)));
                    
                    // Compute the derivative.
                    if (d_first_derivative_uses_local_tol[vi])
                    {
                        // Get the local mean.
                        boost::shared_ptr<pdat::CellData<double> > variable_local_mean(
                            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                patch.getPatchData(d_first_derivative_local_mean_pressure, data_context)));
                        
                        d_gradient_sensor_first_derivative->computeGradientWithVariableLocalMean(
                            patch,
                            data_pressure,
                            derivative,
                            variable_local_mean);
                    }
                    else
                    {
                        d_gradient_sensor_first_derivative->computeGradient(
                            patch,
                            data_pressure,
                            derivative);
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
        else if (sensor_key == "SECOND_DERIVATIVE")
        {
            for (int vi = 0; vi < static_cast<int>(d_second_derivative_variables.size()); vi++)
            {
                // Get the key of the current variable.
                std::string variable_key = d_second_derivative_variables[vi];
                
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
                    
                    // Get the pointer to density data inside the flow model.
                    boost::shared_ptr<pdat::CellData<double> > data_density =
                        d_flow_model->getGlobalCellData("DENSITY");
                    
                    // Get the cell data of the derivative.
                    boost::shared_ptr<pdat::CellData<double> > derivative(
                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                            patch.getPatchData(d_second_derivative_density, data_context)));
                    
                    // Compute the derivative.
                    if (d_second_derivative_uses_local_tol[vi])
                    {
                        // Get the cell data of local mean.
                        boost::shared_ptr<pdat::CellData<double> > variable_local_mean(
                            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                patch.getPatchData(d_second_derivative_local_mean_density, data_context)));
                        
                        d_gradient_sensor_second_derivative->computeGradientWithVariableLocalMean(
                            patch,
                            data_density,
                            derivative,
                            variable_local_mean);
                    }
                    else
                    {
                        d_gradient_sensor_second_derivative->computeGradient(
                            patch,
                            data_density,
                            derivative);
                    }
                    
                    /*
                     * Unregister the patch and data of all registered derived cell variables in the flow model.
                     */
                    
                    d_flow_model->unregisterPatch();
                    
                }
                else if (variable_key ==  "TOTAL_ENERGY")
                {
                    /*
                     * Register the patch and total energy in the flow model and compute the corresponding cell data.
                     */
                    
                    d_flow_model->registerPatchWithDataContext(patch, data_context);
                    
                    std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                    
                    num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("TOTAL_ENERGY", d_num_gradient_ghosts));
                    
                    d_flow_model->registerDerivedCellVariable(num_subghosts_of_data);
                    
                    d_flow_model->computeGlobalDerivedCellData();
                    
                    // Get the pointer to total energy data inside the flow model.
                    boost::shared_ptr<pdat::CellData<double> > data_total_energy =
                        d_flow_model->getGlobalCellData("TOTAL_ENERGY");
                    
                    // Get the cell data of the derivative.
                    boost::shared_ptr<pdat::CellData<double> > derivative(
                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                            patch.getPatchData(d_second_derivative_total_energy, data_context)));
                    
                    // Compute the derivative.
                    if (d_second_derivative_uses_local_tol[vi])
                    {
                        // Get the cell data of local mean.
                        boost::shared_ptr<pdat::CellData<double> > variable_local_mean(
                            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                patch.getPatchData(d_second_derivative_local_mean_total_energy, data_context)));
                        
                        d_gradient_sensor_second_derivative->computeGradientWithVariableLocalMean(
                            patch,
                            data_total_energy,
                            derivative,
                            variable_local_mean);
                    }
                    else
                    {
                        d_gradient_sensor_second_derivative->computeGradient(
                            patch,
                            data_total_energy,
                            derivative);
                    }
                    
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
                    
                    // Get the pointer to pressure data inside the flow model.
                    boost::shared_ptr<pdat::CellData<double> > data_pressure =
                        d_flow_model->getGlobalCellData("PRESSURE");
                    
                    // Get the cell data of the derivative.
                    boost::shared_ptr<pdat::CellData<double> > derivative(
                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                            patch.getPatchData(d_second_derivative_pressure, data_context)));
                    
                    // Compute the derivative.
                    if (d_second_derivative_uses_local_tol[vi])
                    {
                        // Get the local mean.
                        boost::shared_ptr<pdat::CellData<double> > variable_local_mean(
                            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                patch.getPatchData(d_second_derivative_local_mean_pressure, data_context)));
                        
                        d_gradient_sensor_second_derivative->computeGradientWithVariableLocalMean(
                            patch,
                            data_pressure,
                            derivative,
                            variable_local_mean);
                    }
                    else
                    {
                        d_gradient_sensor_second_derivative->computeGradient(
                            patch,
                            data_pressure,
                            derivative);
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
    } // Loop over multiresolution sensors chosen.
}


/*
 * Get the statistics of the sensor values that are required by the
 * gradient sensors at a given patch level.
 */
void
GradientTagger::getSensorValueStatistics(
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const int level_number,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    math::HierarchyCellDataOpsReal<double> cell_double_operator(patch_hierarchy, level_number, level_number);
    
    hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();
    
    for (int si = 0;
             si < static_cast<int>(d_gradient_sensors.size());
             si++)
    {
        std::string sensor_key = d_gradient_sensors[si];
        
        if (sensor_key == "FIRST_DERIVATIVE")
        {
            for (int vi = 0; vi < static_cast<int>(d_first_derivative_variables.size()); vi++)
            {
                // Only get the statistics if global tolerances is used.
                if (d_first_derivative_uses_global_tol[vi])
                {
                    // Get the key of the current variable.
                    std::string variable_key = d_first_derivative_variables[vi];
                    
                    if (variable_key == "DENSITY")
                    {
                        const int w_rho_id = variable_db->mapVariableAndContextToIndex(
                            d_first_derivative_density,
                            data_context);
                        
                        double w_max_rho_local = cell_double_operator.max(w_rho_id);
                        d_first_derivative_max_density = 0.0;
                        
                        mpi.Allreduce(
                            &w_max_rho_local,
                            &d_first_derivative_max_density,
                            1,
                            MPI_DOUBLE,
                            MPI_MAX);
                    }
                    else if (variable_key == "TOTAL_ENERGY")
                    {
                        const int w_E_id = variable_db->mapVariableAndContextToIndex(
                            d_first_derivative_total_energy,
                            data_context);
                        
                        double w_max_E_local = cell_double_operator.max(w_E_id);
                        d_first_derivative_max_total_energy = 0.0;
                        
                        mpi.Allreduce(
                            &w_max_E_local,
                            &d_first_derivative_max_total_energy,
                            1,
                            MPI_DOUBLE,
                            MPI_MAX);
                    }
                    else if (variable_key == "PRESSURE")
                    {
                        const int w_p_id = variable_db->mapVariableAndContextToIndex(
                            d_first_derivative_pressure,
                            data_context);
                        
                        double w_max_p_local = cell_double_operator.max(w_p_id);
                        d_first_derivative_max_pressure = 0.0;
                        
                        mpi.Allreduce(
                            &w_max_p_local,
                            &d_first_derivative_max_pressure,
                            1,
                            MPI_DOUBLE,
                            MPI_MAX);
                    }
                }
            }
        }
        else if (sensor_key == "SECOND_DERIVATIVE")
        {
            for (int vi = 0; vi < static_cast<int>(d_second_derivative_variables.size()); vi++)
            {
                // Only get the statistics if global tolerances is used.
                if (d_second_derivative_uses_global_tol[vi])
                {
                    // Get the key of the current variable.
                    std::string variable_key = d_second_derivative_variables[vi];
                    
                    if (variable_key == "DENSITY")
                    {
                        const int w_rho_id = variable_db->mapVariableAndContextToIndex(
                            d_second_derivative_density,
                            data_context);
                        
                        double w_max_rho_local = cell_double_operator.max(w_rho_id);
                        d_second_derivative_max_density = 0.0;
                        
                        mpi.Allreduce(
                            &w_max_rho_local,
                            &d_second_derivative_max_density,
                            1,
                            MPI_DOUBLE,
                            MPI_MAX);
                    }
                    else if (variable_key == "TOTAL_ENERGY")
                    {
                        const int w_E_id = variable_db->mapVariableAndContextToIndex(
                            d_second_derivative_total_energy,
                            data_context);
                        
                        double w_max_E_local = cell_double_operator.max(w_E_id);
                        d_second_derivative_max_total_energy = 0.0;
                        
                        mpi.Allreduce(
                            &w_max_E_local,
                            &d_second_derivative_max_total_energy,
                            1,
                            MPI_DOUBLE,
                            MPI_MAX);
                    }
                    else if (variable_key == "PRESSURE")
                    {
                        const int w_p_id = variable_db->mapVariableAndContextToIndex(
                            d_second_derivative_pressure,
                            data_context);
                        
                        double w_max_p_local = cell_double_operator.max(w_p_id);
                        d_second_derivative_max_pressure = 0.0;
                        
                        mpi.Allreduce(
                            &w_max_p_local,
                            &d_second_derivative_max_pressure,
                            1,
                            MPI_DOUBLE,
                            MPI_MAX);
                    }
                }
            }
        }
    }
}


/*
 * Tag cells for refinement using gradient sensors.
 */
void
GradientTagger::tagCells(
   hier::Patch& patch,
   const boost::shared_ptr<pdat::CellData<int> >& tags,
   const boost::shared_ptr<hier::VariableContext>& data_context)
{
    // Loop over gradient sensors chosen.
    for (int si = 0;
             si < static_cast<int>(d_gradient_sensors.size());
             si++)
    {
        std::string sensor_key = d_gradient_sensors[si];
        
        if (sensor_key == "FIRST_DERIVATIVE")
        {
            int count_global_tol = 0;
            int count_local_tol = 0;
            
            // Looop over variables chosen.
            for (int vi = 0; vi < static_cast<int>(d_first_derivative_variables.size()); vi++)
            {
                // Get the key of the current variable.
                std::string variable_key = d_first_derivative_variables[vi];
                
                // Get tolerance settings for the current variable.
                bool uses_global_tol = d_first_derivative_uses_global_tol[vi];
                bool uses_local_tol = d_first_derivative_uses_local_tol[vi];
                
                double global_tol = 0.0;
                double local_tol = 0.0;
                
                if (uses_global_tol)
                {
                    global_tol = d_first_derivative_global_tol[count_global_tol];
                    count_global_tol++;
                }
                
                if (uses_local_tol)
                {
                    local_tol = d_first_derivative_local_tol[count_local_tol];
                    count_local_tol++;
                }
                
                /*
                 * Get the derivative and the statistics of the derivative at different levels
                 * and tag cells.
                 */
                
                if (variable_key == "DENSITY")
                {
                    // Get the cell data of the derivative.
                    boost::shared_ptr<pdat::CellData<double> > derivative(
                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                            patch.getPatchData(d_first_derivative_density, data_context)));
                    
                    boost::shared_ptr<pdat::CellData<double> > variable_local_mean;
                    if (uses_local_tol)
                    {
                        // Get the cell data of local mean.
                        boost::shared_ptr<pdat::CellData<double> > variable_local_mean(
                            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                patch.getPatchData(d_first_derivative_local_mean_density, data_context)));
                    }
                    
                    tagCellsWithDerivativeSensor(
                        patch,
                        tags,
                        derivative,
                        d_first_derivative_max_density,
                        variable_local_mean,
                        uses_global_tol,
                        uses_local_tol,
                        global_tol,
                        local_tol);
                }
                else if (variable_key == "TOTAL_ENERGY")
                {
                    // Get the cell data of the derivative.
                    boost::shared_ptr<pdat::CellData<double> > derivative(
                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                            patch.getPatchData(d_first_derivative_total_energy, data_context)));
                    
                    boost::shared_ptr<pdat::CellData<double> > variable_local_mean;
                    if (uses_local_tol)
                    {
                        // Get the cell data of local mean.
                        boost::shared_ptr<pdat::CellData<double> > variable_local_mean(
                            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                patch.getPatchData(d_first_derivative_local_mean_total_energy, data_context)));
                    }
                    
                    tagCellsWithDerivativeSensor(
                        patch,
                        tags,
                        derivative,
                        d_first_derivative_max_total_energy,
                        variable_local_mean,
                        uses_global_tol,
                        uses_local_tol,
                        global_tol,
                        local_tol);
                }
                else if (variable_key == "PRESSURE")
                {
                    // Get the cell data of the derivative.
                    boost::shared_ptr<pdat::CellData<double> > derivative(
                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                            patch.getPatchData(d_first_derivative_pressure, data_context)));
                    
                    boost::shared_ptr<pdat::CellData<double> > variable_local_mean;
                    if (uses_local_tol)
                    {
                        // Get the cell data of local mean.
                        boost::shared_ptr<pdat::CellData<double> > variable_local_mean(
                            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                patch.getPatchData(d_first_derivative_local_mean_pressure, data_context)));
                    }
                    
                    tagCellsWithDerivativeSensor(
                        patch,
                        tags,
                        derivative,
                        d_first_derivative_max_pressure,
                        variable_local_mean,
                        uses_global_tol,
                        uses_local_tol,
                        global_tol,
                        local_tol);
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
        else if (sensor_key == "SECOND_DERIVATIVE")
        {
            int count_global_tol = 0;
            int count_local_tol = 0;
            
            // Looop over variables chosen.
            for (int vi = 0; vi < static_cast<int>(d_second_derivative_variables.size()); vi++)
            {
                // Get the key of the current variable.
                std::string variable_key = d_second_derivative_variables[vi];
                
                // Get tolerance settings for the current variable.
                bool uses_global_tol = d_second_derivative_uses_global_tol[vi];
                bool uses_local_tol = d_second_derivative_uses_local_tol[vi];
                
                double global_tol = 0.0;
                double local_tol = 0.0;
                
                if (uses_global_tol)
                {
                    global_tol = d_second_derivative_global_tol[count_global_tol];
                    count_global_tol++;
                }
                
                if (uses_local_tol)
                {
                    local_tol = d_second_derivative_local_tol[count_local_tol];
                    count_local_tol++;
                }
                
                /*
                 * Get the derivative and the statistics of the derivative at different levels
                 * and tag cells.
                 */
                
                if (variable_key == "DENSITY")
                {
                    // Get the cell data of the derivative.
                    boost::shared_ptr<pdat::CellData<double> > derivative(
                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                            patch.getPatchData(d_second_derivative_density, data_context)));
                    
                    boost::shared_ptr<pdat::CellData<double> > variable_local_mean;
                    if (uses_local_tol)
                    {
                        // Get the cell data of local mean.
                        boost::shared_ptr<pdat::CellData<double> > variable_local_mean(
                            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                patch.getPatchData(d_second_derivative_local_mean_density, data_context)));
                    }
                    
                    tagCellsWithDerivativeSensor(
                        patch,
                        tags,
                        derivative,
                        d_second_derivative_max_density,
                        variable_local_mean,
                        uses_global_tol,
                        uses_local_tol,
                        global_tol,
                        local_tol);
                }
                else if (variable_key == "TOTAL_ENERGY")
                {
                    // Get the cell data of the derivative.
                    boost::shared_ptr<pdat::CellData<double> > derivative(
                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                            patch.getPatchData(d_second_derivative_total_energy, data_context)));
                    
                    boost::shared_ptr<pdat::CellData<double> > variable_local_mean;
                    if (uses_local_tol)
                    {
                        // Get the cell data of local mean.
                        boost::shared_ptr<pdat::CellData<double> > variable_local_mean(
                            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                patch.getPatchData(d_second_derivative_local_mean_total_energy, data_context)));
                    }
                    
                    tagCellsWithDerivativeSensor(
                        patch,
                        tags,
                        derivative,
                        d_second_derivative_max_total_energy,
                        variable_local_mean,
                        uses_global_tol,
                        uses_local_tol,
                        global_tol,
                        local_tol);
                }
                else if (variable_key == "PRESSURE")
                {
                    // Get the cell data of the derivative.
                    boost::shared_ptr<pdat::CellData<double> > derivative(
                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                            patch.getPatchData(d_second_derivative_pressure, data_context)));
                    
                    boost::shared_ptr<pdat::CellData<double> > variable_local_mean;
                    if (uses_local_tol)
                    {
                        // Get the cell data of local mean.
                        boost::shared_ptr<pdat::CellData<double> > variable_local_mean(
                            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                patch.getPatchData(d_second_derivative_local_mean_pressure, data_context)));
                    }
                    
                    tagCellsWithDerivativeSensor(
                        patch,
                        tags,
                        derivative,
                        d_second_derivative_max_pressure,
                        variable_local_mean,
                        uses_global_tol,
                        uses_local_tol,
                        global_tol,
                        local_tol);
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
        else if (sensor_key == "JAMESON_GRADIENT")
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
                            patch.getPatchData(d_Jameson_gradient_density, data_context)));
                    
                    // Compute the gradient.
                    d_gradient_sensor_Jameson->computeGradient(patch, data_density, gradient);
                    
                    // Tag the cells.
                    tagCellsWithGradientSensor(
                        patch,
                        sensor_key,
                        tags,
                        gradient,
                        tol);
                    
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
                            patch.getPatchData(d_Jameson_gradient_total_energy, data_context)));
                    
                    // Compute the gradient.
                    d_gradient_sensor_Jameson->computeGradient(patch, data_total_energy, gradient);
                    
                    // Tag the cells.
                    tagCellsWithGradientSensor(
                        patch,
                        sensor_key,
                        tags,
                        gradient,
                        tol);
                    
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
                            patch.getPatchData(d_Jameson_gradient_pressure, data_context)));
                    
                    // Compute the gradient.
                    d_gradient_sensor_Jameson->computeGradient(patch, data_pressure, gradient);
                    
                    // Tag the cells.
                    tagCellsWithGradientSensor(
                        patch,
                        sensor_key,
                        tags,
                        gradient,
                        tol);
                    
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
 * Tag cells using value of gradient sensor.
 */
void
GradientTagger::tagCellsWithGradientSensor(
    hier::Patch& patch,
    const std::string& sensor_key,
    const boost::shared_ptr<pdat::CellData<int> >& tags,
    const boost::shared_ptr<pdat::CellData<double> >& gradient,
    const double tol)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(tags->getGhostCellWidth() == hier::IntVector::getZero(d_dim));
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the number of ghost cells and dimensions of box that covers interior of patch plus
    // ghost cells.
    const hier::IntVector num_ghosts_gradient_tagger = gradient->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_gradient_tagger = gradient->getGhostBox().numberCells();
    
    // Get the pointer of the tags.
    int* tag_ptr  = tags->getPointer(0);
    
    // Get the pointer to the data.
    double* psi = gradient->getPointer(0);
    
    if (sensor_key == "JAMESON_GRADIENT")
    {
        if (d_dim == tbox::Dimension(1))
        {
            const int interior_dim_0 = interior_dims[0];
            
            const int num_ghosts_0_gradient_tagger = num_ghosts_gradient_tagger[0];
            
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute indices.
                const int idx = i + num_ghosts_0_gradient_tagger;
                const int idx_nghost = i;
                
                if (psi[idx] > tol)
                {
                    tag_ptr[idx_nghost] |= 1;
                }
            }
        }
        else if (d_dim == tbox::Dimension(2))
        {
            const int interior_dim_0 = interior_dims[0];
            const int interior_dim_1 = interior_dims[1];
            
            const int num_ghosts_0_gradient_tagger = num_ghosts_gradient_tagger[0];
            const int num_ghosts_1_gradient_tagger = num_ghosts_gradient_tagger[1];
            const int ghostcell_dim_0_gradient_tagger = ghostcell_dims_gradient_tagger[0];
            
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute indices.
                    const int idx = (i + num_ghosts_0_gradient_tagger) +
                        (j + num_ghosts_1_gradient_tagger)*ghostcell_dim_0_gradient_tagger;
                    
                    const int idx_nghost = i +
                        j*interior_dim_0;
                    
                    if (psi[idx] > tol)
                    {
                        tag_ptr[idx_nghost] |= 1;
                    }
                }
            }
        }
        else if (d_dim == tbox::Dimension(3))
        {
            const int interior_dim_0 = interior_dims[0];
            const int interior_dim_1 = interior_dims[1];
            const int interior_dim_2 = interior_dims[2];
            
            const int num_ghosts_0_gradient_tagger = num_ghosts_gradient_tagger[0];
            const int num_ghosts_1_gradient_tagger = num_ghosts_gradient_tagger[1];
            const int num_ghosts_2_gradient_tagger = num_ghosts_gradient_tagger[2];
            const int ghostcell_dim_0_gradient_tagger = ghostcell_dims_gradient_tagger[0];
            const int ghostcell_dim_1_gradient_tagger = ghostcell_dims_gradient_tagger[1];
            
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute indices.
                        const int idx = (i + num_ghosts_0_gradient_tagger) +
                            (j + num_ghosts_1_gradient_tagger)*ghostcell_dim_0_gradient_tagger +
                            (k + num_ghosts_2_gradient_tagger)*ghostcell_dim_0_gradient_tagger*
                                ghostcell_dim_1_gradient_tagger;
                        
                        const int idx_nghost = i +
                            j*interior_dim_0 +
                            k*interior_dim_0*interior_dim_1;
                        
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


/*
 * Tag cells using value of derivative sensor.
 */
void
GradientTagger::tagCellsWithDerivativeSensor(
    hier::Patch& patch,
    const boost::shared_ptr<pdat::CellData<int> >& tags,
    const boost::shared_ptr<pdat::CellData<double> >& derivative,
    const double& derivative_max,
    const boost::shared_ptr<pdat::CellData<double> >& variable_local_mean,
    const bool& uses_global_tol,
    const bool& uses_local_tol,
    const double& global_tol,
    const double& local_tol)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(tags->getGhostCellWidth() == hier::IntVector::getZero(d_dim));
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the number of ghost cells and dimensions of box that covers interior of patch plus
    // ghost cells.
    const hier::IntVector num_ghosts_derivative = derivative->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_derivative = derivative->getGhostBox().numberCells();
    
    // Allocate temporary patch data.
    boost::shared_ptr<pdat::CellData<int> > tags_gradient_tagger(
        new pdat::CellData<int>(interior_box, d_dim.getValue(), hier::IntVector::getZero(d_dim)));
    
    tags_gradient_tagger->fillAll(1);
    
    // Get the pointers to the tags.
    int* tag_ptr_gradient_tagger = tags_gradient_tagger->getPointer(0);
    int* tag_ptr = tags->getPointer(0);
    
    // Get the pointers to the derivative.
    double* w = derivative->getPointer(0);
    
    // Get the pointers to the variable local means.
    double* u_mean = variable_local_mean->getPointer(0);
    
    if (d_dim == tbox::Dimension(1))
    {
        const int interior_dim_0 = interior_dims[0];
        
        const int num_ghosts_0_derivative = num_ghosts_derivative[0];
        
        if (uses_global_tol)
        {
            // Allocate temporary patch data.
            boost::shared_ptr<pdat::CellData<int> > tags_global_tol(
                new pdat::CellData<int>(interior_box, d_dim.getValue(), hier::IntVector::getZero(d_dim)));
            
            tags_global_tol->fillAll(0);
            
            // Get the pointer to the tags.
            int* tag_ptr_global_tol = tags_global_tol->getPointer(0);
            
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx = i + num_ghosts_0_derivative;
                const int idx_nghost = i;
                
                if (w[idx]/(derivative_max + EPSILON) > global_tol)
                {
                    tag_ptr_global_tol[idx_nghost] = 1;
                }
            }
            
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear index.
                const int idx_nghost = i;
                
                tag_ptr_gradient_tagger[idx_nghost] &= tag_ptr_global_tol[idx_nghost];
            }
        }
        
        if (uses_local_tol)
        {
            // Allocate temporary patch data.
            boost::shared_ptr<pdat::CellData<int> > tags_local_tol(
                new pdat::CellData<int>(interior_box, d_dim.getValue(), hier::IntVector::getZero(d_dim)));
            
            tags_local_tol->fillAll(0);
            
            // Get the pointer to the tags.
            int* tag_ptr_local_tol = tags_local_tol->getPointer(0);
            
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx = i + num_ghosts_0_derivative;
                const int idx_nghost = i;
                
                if (w[idx]/(u_mean[idx] + EPSILON) > local_tol)
                {
                    tag_ptr_local_tol[idx_nghost] = 1;
                }
            }
            
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear index.
                const int idx_nghost = i;
                
                tag_ptr_gradient_tagger[idx_nghost] &= tag_ptr_local_tol[idx_nghost];
            }
        }
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = 0; i < interior_dim_0; i++)
        {
            // Compute the linear index.
            const int idx_nghost = i;
            
            tag_ptr[idx_nghost] |= tag_ptr_gradient_tagger[idx_nghost];
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        
        const int num_ghosts_0_derivative = num_ghosts_derivative[0];
        const int num_ghosts_1_derivative = num_ghosts_derivative[1];
        const int ghostcell_dim_0_derivative = ghostcell_dims_derivative[0];
        
        if (uses_global_tol)
        {
            // Allocate temporary patch data.
            boost::shared_ptr<pdat::CellData<int> > tags_global_tol(
                new pdat::CellData<int>(interior_box, d_dim.getValue(), hier::IntVector::getZero(d_dim)));
            
            tags_global_tol->fillAll(0);
            
            // Get the pointer to the tags.
            int* tag_ptr_global_tol = tags_global_tol->getPointer(0);
            
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx = (i + num_ghosts_0_derivative) +
                        (j + num_ghosts_1_derivative)*ghostcell_dim_0_derivative;
                    
                    const int idx_nghost = i +
                        j*interior_dim_0;
                    
                    if (w[idx]/(derivative_max + EPSILON) > global_tol)
                    {
                        tag_ptr_global_tol[idx_nghost] = 1;
                    }
                }
            }
            
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_nghost = i +
                        j*interior_dim_0;
                    
                    tag_ptr_gradient_tagger[idx_nghost] &= tag_ptr_global_tol[idx_nghost];
                }
            }
        }
        
        if (uses_local_tol)
        {
            // Allocate temporary patch data.
            boost::shared_ptr<pdat::CellData<int> > tags_local_tol(
                new pdat::CellData<int>(interior_box, d_dim.getValue(), hier::IntVector::getZero(d_dim)));
            
            tags_local_tol->fillAll(0);
            
            // Get the pointer to the tags.
            int* tag_ptr_local_tol = tags_local_tol->getPointer(0);
            
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx = (i + num_ghosts_0_derivative) +
                        (j + num_ghosts_1_derivative)*ghostcell_dim_0_derivative;
                    
                    const int idx_nghost = i +
                        j*interior_dim_0;
                    
                    if (w[idx]/(u_mean[idx] + EPSILON) > local_tol)
                    {
                        tag_ptr_local_tol[idx_nghost] = 1;
                    }
                }
            }
            
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_nghost = i +
                        j*interior_dim_0;
                    
                    tag_ptr_gradient_tagger[idx_nghost] &= tag_ptr_local_tol[idx_nghost];
                }
            }
        }
        
        for (int j = 0; j < interior_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear index.
                const int idx_nghost = i +
                    j*interior_dim_0;
                
                tag_ptr[idx_nghost] |= tag_ptr_gradient_tagger[idx_nghost];
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        const int interior_dim_2 = interior_dims[2];
        
        const int num_ghosts_0_derivative = num_ghosts_derivative[0];
        const int num_ghosts_1_derivative = num_ghosts_derivative[1];
        const int num_ghosts_2_derivative = num_ghosts_derivative[2];
        const int ghostcell_dim_0_derivative = ghostcell_dims_derivative[0];
        const int ghostcell_dim_1_derivative = ghostcell_dims_derivative[1];
        
        if (uses_global_tol)
        {
            // Allocate temporary patch data.
            boost::shared_ptr<pdat::CellData<int> > tags_global_tol(
                new pdat::CellData<int>(interior_box, d_dim.getValue(), hier::IntVector::getZero(d_dim)));
            
            tags_global_tol->fillAll(0);
            
            // Get the pointer to the tags.
            int* tag_ptr_global_tol = tags_global_tol->getPointer(0);
            
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx = (i + num_ghosts_0_derivative) +
                            (j + num_ghosts_1_derivative)*ghostcell_dim_0_derivative +
                            (k + num_ghosts_2_derivative)*ghostcell_dim_0_derivative*
                                ghostcell_dim_1_derivative;
                        
                        const int idx_nghost = i +
                            j*interior_dim_0 +
                            k*interior_dim_0*
                                interior_dim_1;
                        
                        if (w[idx]/(derivative_max + EPSILON) > global_tol)
                        {
                            tag_ptr_global_tol[idx_nghost] = 1;
                        }
                    }
                }
            }
            
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx_nghost = i +
                            j*interior_dim_0 +
                            k*interior_dim_0*
                                interior_dim_1;
                        
                        tag_ptr_gradient_tagger[idx_nghost] &= tag_ptr_global_tol[idx_nghost];
                    }
                }
            }
        }
        
        if (uses_local_tol)
        {
            // Allocate temporary patch data.
            boost::shared_ptr<pdat::CellData<int> > tags_local_tol(
                new pdat::CellData<int>(interior_box, d_dim.getValue(), hier::IntVector::getZero(d_dim)));
            
            tags_local_tol->fillAll(0);
            
            // Get the pointer to the tags.
            int* tag_ptr_local_tol = tags_local_tol->getPointer(0);
            
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx = (i + num_ghosts_0_derivative) +
                            (j + num_ghosts_1_derivative)*ghostcell_dim_0_derivative +
                            (k + num_ghosts_2_derivative)*ghostcell_dim_0_derivative*
                                ghostcell_dim_1_derivative;
                        
                        const int idx_nghost = i +
                            j*interior_dim_0 +
                            k*interior_dim_0*
                                interior_dim_1;
                        
                        if (w[idx]/(u_mean[idx] + EPSILON) > local_tol)
                        {
                            tag_ptr_local_tol[idx_nghost] = 1;
                        }
                    }
                }
            }
            
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx_nghost = i +
                            j*interior_dim_0 +
                            k*interior_dim_0*
                                interior_dim_1;
                        
                        tag_ptr_gradient_tagger[idx_nghost] &= tag_ptr_local_tol[idx_nghost];
                    }
                }
            }
        }
        
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_nghost = i +
                            j*interior_dim_0 +
                            k*interior_dim_0*
                                interior_dim_1;
                    
                    tag_ptr[idx_nghost] |= tag_ptr_gradient_tagger[idx_nghost];
                }
            }
        }
    }
}
