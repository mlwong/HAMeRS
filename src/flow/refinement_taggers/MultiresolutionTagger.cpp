#include "flow/refinement_taggers/MultiresolutionTagger.hpp"

#include <algorithm>

// #define HAMERS_PLOTTING_MULTIRESOLUTION_TAGGER

#define EPSILON HAMERS_EPSILON

MultiresolutionTagger::MultiresolutionTagger(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
    const HAMERS_SHARED_PTR<FlowModel>& flow_model,
    const HAMERS_SHARED_PTR<tbox::Database>& multiresolution_tagger_db):
        d_object_name(object_name),
        d_dim(dim),
        d_grid_geometry(grid_geometry),
        d_num_multiresolution_ghosts(hier::IntVector::getZero(d_dim)),
        d_flow_model(flow_model)
{
    if (multiresolution_tagger_db != nullptr)
    {
        std::vector<std::string> sensor_keys = multiresolution_tagger_db->getAllKeys();
        
        const int num_keys = static_cast<int>(sensor_keys.size());
        
        if (multiresolution_tagger_db->keyExists("multiresolution_sensors"))
        {
            d_multiresolution_sensors =
                multiresolution_tagger_db->getStringVector("multiresolution_sensors");
        }
        else if (multiresolution_tagger_db->keyExists("d_multiresolution_sensors"))
        {
            d_multiresolution_sensors =
                multiresolution_tagger_db->getStringVector("d_multiresolution_sensors");
        }
        else
        {
            TBOX_WARNING(d_object_name
                << ": "
                << "No key 'multiresolution_sensors'/"
                << "'d_multiresolution_sensors' found in data for"
                << " Multiresolution_tagger. No refinement with multiresolution sensors will occur."
                << std::endl);
        }
        
        /*
         * Loop over the multiresolution sensors chosen.
         */
        
        std::vector<std::string> sensor_keys_defined(num_keys);
        int sensor_keys_count = 0;
        HAMERS_SHARED_PTR<tbox::Database> sensor_db;
        for (int i = 0; i < num_keys; i++)
        {
            std::string sensor_key = sensor_keys[i];
            sensor_db.reset();
            
            if (!((sensor_key == "multiresolution_sensors") ||
                  (sensor_key == "d_multiresolution_sensors")))
            {
                if (!(sensor_key == "HARTEN_WAVELET"))
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
                    sensor_db = multiresolution_tagger_db->getDatabase(sensor_key);
                    sensor_keys_defined[sensor_keys_count] = sensor_key;
                    sensor_keys_count++;
                }
                
                if (sensor_db && sensor_key == "HARTEN_WAVELET")
                {
                    // Get the number of wavelet levels.
                    if (sensor_db->keyExists("Harten_wavelet_num_level"))
                    {
                        d_Harten_wavelet_num_level =
                            sensor_db->getInteger("Harten_wavelet_num_level");
                    }
                    else if (sensor_db->keyExists("d_Harten_wavelet_num_level"))
                    {
                        d_Harten_wavelet_num_level =
                            sensor_db->getInteger("d_Harten_wavelet_num_level");
                    }
                    else
                    {
                        TBOX_ERROR(d_object_name
                            << ": "
                            << "No key 'Harten_wavelet_num_level'/"
                            << "'d_Harten_wavelet_num_level' found in data for "
                            << sensor_key
                            << "."
                            << std::endl);
                    }
                    
                    if (d_Harten_wavelet_num_level < 2 || d_Harten_wavelet_num_level > 4)
                    {
                        TBOX_ERROR(d_object_name
                            << ": "
                            << "Only number of wavelet levels between 2 and 4 is allowed. \n"
                            << "Number of Harten wavelet levels = "
                            << d_Harten_wavelet_num_level
                            << " is provided."
                            << std::endl);
                    }
                    
                    // Get the number of vanishing moments.
                    if (sensor_db->keyExists("Harten_wavelet_num_vanishing_moments"))
                    {
                        d_Harten_wavelet_num_vanishing_moments =
                            sensor_db->getInteger("Harten_wavelet_num_vanishing_moments");
                    }
                    else if (sensor_db->keyExists("d_Harten_wavelet_num_vanishing_moments"))
                    {
                        d_Harten_wavelet_num_vanishing_moments =
                            sensor_db->getInteger("d_Harten_wavelet_num_vanishing_moments");
                    }
                    else
                    {
                        TBOX_ERROR(d_object_name
                            << ": "
                            << "No key 'Harten_wavelet_num_vanishing_moments'/"
                            << "'d_Harten_wavelet_num_vanishing_moments' found in data for "
                            << sensor_key
                            << "."
                            << std::endl);
                    }
                    
                    d_wavelet_transfrom_Harten.reset(new WaveletTransformHarten(
                        "Harten wavelet transform",
                        d_dim,
                        d_Harten_wavelet_num_level,
                        d_Harten_wavelet_num_vanishing_moments));
                    
                    // Get the variables for the wavelet sensor to apply on.
                    if (sensor_db->keyExists("Harten_wavelet_variables"))
                    {
                        d_Harten_wavelet_variables =
                            sensor_db->getStringVector("Harten_wavelet_variables");
                    }
                    else if (sensor_db->keyExists("d_Harten_wavelet_variables"))
                    {
                        d_Harten_wavelet_variables =
                            sensor_db->getStringVector("d_Harten_wavelet_variables");
                    }
                    else
                    {
                        TBOX_ERROR(d_object_name
                            << ": "
                            << "No key 'Harten_wavelet_variables'/"
                            << "'d_Harten_wavelet_variables' found in database for "
                            << sensor_key
                            << "."
                            << std::endl);
                    }
                    
                    for (int vi = 0; vi < static_cast<int>(d_Harten_wavelet_variables.size()); vi++)
                    {
                        std::string variable_key = d_Harten_wavelet_variables[vi];
                        
                        if (!((variable_key == "DENSITY") ||
                              (variable_key == "TOTAL_ENERGY") ||
                              (variable_key == "PRESSURE")))
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "Unknown/unsupported variable '"
                                << variable_key
                                << "'in input."
                                << std::endl);
                        }
                    }
                    
                    int num_true = 0;
                    
                    /*
                     * Get the settings for wavelet global tolerances.
                     */
                    
                    if (sensor_db->keyExists("Harten_wavelet_uses_global_tol"))
                    {
                        d_Harten_wavelet_uses_global_tol =
                            sensor_db->getBoolVector("Harten_wavelet_uses_global_tol");
                    }
                    else if (sensor_db->keyExists("d_Harten_wavelet_uses_global_tol"))
                    {
                        d_Harten_wavelet_uses_global_tol =
                            sensor_db->getBoolVector("d_Harten_wavelet_uses_global_tol");
                    }
                    else
                    {
                        TBOX_ERROR(d_object_name
                            << ": "
                            << "No key 'Harten_wavelet_uses_global_tol'/"
                            << "'d_Harten_wavelet_uses_global_tol'"
                            << " found in database for "
                            << sensor_key
                            << "."
                            << std::endl);
                    }
                    
                    if (static_cast<int>(d_Harten_wavelet_variables.size()) !=
                        static_cast<int>(d_Harten_wavelet_uses_global_tol.size()))
                    {
                        TBOX_ERROR(d_object_name
                            << ": "
                            << "The numbers of variables and switches for global tolerances"
                            << " provided don't match"
                            << " in database for "
                            << sensor_key
                            << "."
                            << std::endl);
                    }
                    
                    num_true = std::count(d_Harten_wavelet_uses_global_tol.begin(),
                        d_Harten_wavelet_uses_global_tol.end(),
                        true);
                    
                    if (num_true > 0)
                    {
                        if (sensor_db->keyExists("Harten_wavelet_global_tol"))
                        {
                            d_Harten_wavelet_global_tol =
                                sensor_db->getDoubleVector("Harten_wavelet_global_tol");
                        }
                        else if (sensor_db->keyExists("d_Harten_wavelet_global_tol"))
                        {
                            d_Harten_wavelet_global_tol =
                                sensor_db->getDoubleVector("d_Harten_wavelet_global_tol");
                        }
                        else
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "No key 'Harten_wavelet_global_tol'/"
                                << "'d_Harten_wavelet_global_tol' found in database for "
                                << sensor_key
                                << "."
                                << std::endl);
                        }
                        
                        if (num_true != static_cast<int>(d_Harten_wavelet_global_tol.size()))
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "The number of variables that use global tolerances and number"
                                << " of global tolerances provided don't match"
                                << " in database for "
                                << sensor_key
                                << "."
                                << std::endl);
                        }
                    }
                    
                    /*
                     * Get the settings for wavelet local tolerances.
                     */
                    
                    if (sensor_db->keyExists("Harten_wavelet_uses_local_tol"))
                    {
                        d_Harten_wavelet_uses_local_tol =
                            sensor_db->getBoolVector("Harten_wavelet_uses_local_tol");
                    }
                    else if (sensor_db->keyExists("d_Harten_wavelet_uses_local_tol"))
                    {
                        d_Harten_wavelet_uses_local_tol =
                            sensor_db->getBoolVector("d_Harten_wavelet_uses_local_tol");
                    }
                    else
                    {
                        TBOX_ERROR(d_object_name
                            << ": "
                            << "No key 'Harten_wavelet_uses_local_tol'/"
                            << "'d_Harten_wavelet_uses_local_tol'"
                            << " found in database for "
                            << sensor_key
                            << "."
                            << std::endl);
                    }
                    
                    if (static_cast<int>(d_Harten_wavelet_variables.size()) !=
                        static_cast<int>(d_Harten_wavelet_uses_local_tol.size()))
                    {
                        TBOX_ERROR(d_object_name
                            << ": "
                            << "The numbers of variables and switches for local tolerances"
                            << " provided don't match"
                            << " in database for "
                            << sensor_key
                            << "."
                            << std::endl);
                    }
                    
                    num_true = std::count(d_Harten_wavelet_uses_local_tol.begin(),
                        d_Harten_wavelet_uses_local_tol.end(),
                        true);
                    
                    if (num_true > 0)
                    {
                        if (sensor_db->keyExists("Harten_wavelet_local_tol"))
                        {
                            d_Harten_wavelet_local_tol =
                                sensor_db->getDoubleVector("Harten_wavelet_local_tol");
                        }
                        else if (sensor_db->keyExists("d_Harten_wavelet_local_tol"))
                        {
                            d_Harten_wavelet_local_tol =
                                sensor_db->getDoubleVector("d_Harten_wavelet_local_tol");
                        }
                        else
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "No key 'Harten_wavelet_local_tol'/"
                                << "'d_Harten_wavelet_local_tol' found in database for "
                                << sensor_key
                                << "."
                                << std::endl);
                        }
                        
                        if (num_true != static_cast<int>(d_Harten_wavelet_local_tol.size()))
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "The number of variables that use local tolerances and number"
                                << " of local tolerances provided don't match"
                                << " in database for "
                                << sensor_key
                                << "."
                                << std::endl);
                        }
                    }
                    
                    /*
                     * Get the settings for wavelet Lipschitz's tolerances.
                     */
                    
                    if (sensor_db->keyExists("Harten_wavelet_uses_alpha_tol"))
                    {
                        d_Harten_wavelet_uses_alpha_tol =
                            sensor_db->getBoolVector("Harten_wavelet_uses_alpha_tol");
                    }
                    else if (sensor_db->keyExists("d_Harten_wavelet_uses_alpha_tol"))
                    {
                        d_Harten_wavelet_uses_alpha_tol =
                            sensor_db->getBoolVector("d_Harten_wavelet_uses_alpha_tol");
                    }
                    else
                    {
                        TBOX_ERROR(d_object_name
                            << ": "
                            << "No key 'Harten_wavelet_uses_alpha_tol'/"
                            << "'d_Harten_wavelet_uses_alpha_tol'"
                            << " found in database for "
                            << sensor_key
                            << "."
                            << std::endl);
                    }
                    
                    if (static_cast<int>(d_Harten_wavelet_variables.size()) !=
                        static_cast<int>(d_Harten_wavelet_uses_alpha_tol.size()))
                    {
                        TBOX_ERROR(d_object_name
                            << ": "
                            << "The numbers of variables and switches for Lipschitz's tolerances"
                            << " provided don't match"
                            << " in database for "
                            << sensor_key
                            << "."
                            << std::endl);
                    }
                    
                    num_true = std::count(d_Harten_wavelet_uses_alpha_tol.begin(),
                        d_Harten_wavelet_uses_alpha_tol.end(),
                        true);
                    
                    if (num_true > 0)
                    {
                        if (sensor_db->keyExists("Harten_wavelet_alpha_tol"))
                        {
                            d_Harten_wavelet_alpha_tol =
                                sensor_db->getDoubleVector("Harten_wavelet_alpha_tol");
                        }
                        else if (sensor_db->keyExists("d_Harten_wavelet_alpha_tol"))
                        {
                            d_Harten_wavelet_alpha_tol =
                                sensor_db->getDoubleVector("d_Harten_wavelet_alpha_tol");
                        }
                        else
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "No key 'Harten_wavelet_alpha_tol'/"
                                << "'d_Harten_wavelet_alpha_tol' found in database for "
                                << sensor_key
                                << "."
                                << std::endl);
                        }
                        
                        if (num_true != static_cast<int>(d_Harten_wavelet_alpha_tol.size()))
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "The number of variables that use Lipschitz's tolerances and number"
                                << " of Lipschitz's tolerances provided don't match"
                                << " in database for "
                                << sensor_key
                                << "."
                                << std::endl);
                        }
                    }
                }
            }
        } // Loop over sensors.
        
        /*
         * Check that input is found for each string identifier in key list.
         */
        for (int ki = 0;
             ki < static_cast<int>(d_multiresolution_sensors.size());
             ki++)
        {
            std::string use_key = d_multiresolution_sensors[ki];
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
                    << "No input found for specified mulitresolution sensor: "
                    << d_multiresolution_sensors[ki]
                    << "."
                    << std::endl);
            }
        }
        
    }
    else
    {
        TBOX_WARNING(d_object_name
            << ": "
            << "Key data 'Multiresolution_tagger' not found in input/restart database."
            << " No refinement with multiresolution sensors will occur."
            << std::endl);
    }
    
    /*
     * Compute the number of ghost cells required by the multiresolution detecor.
     */
    
    if (d_wavelet_transfrom_Harten != nullptr)
    {
        d_num_multiresolution_ghosts = hier::IntVector::max(
            d_num_multiresolution_ghosts,
            d_wavelet_transfrom_Harten->getWaveletTransformNumberOfGhostCells());
    }
}


/*
 * Register the temporary variables used in multiresolution tagger class.
 */
void
MultiresolutionTagger::registerMultiresolutionTaggerVariables(
    RungeKuttaLevelIntegrator* integrator)
{
    for (int si = 0; si < static_cast<int>(d_multiresolution_sensors.size()); si++)
    {
        std::string sensor_key = d_multiresolution_sensors[si];
        
        if (sensor_key == "HARTEN_WAVELET")
        {
            for (int vi = 0; vi < static_cast<int>(d_Harten_wavelet_variables.size()); vi++)
            {
                // Get the key of the current variable.
                std::string variable_key = d_Harten_wavelet_variables[vi];
                
                if (variable_key == "DENSITY")
                {
                    for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                    {
                        d_Harten_wavelet_coeffs_density.push_back(
                            HAMERS_MAKE_SHARED<pdat::CellVariable<double> >(
                                d_dim,
                                "Harten wavelet coefficient of density at level " +
                                    std::to_string(li),
                                1));
                        
                        if (d_Harten_wavelet_uses_global_tol[vi])
                        {
                            d_Harten_wavelet_coeffs_maxs_density.push_back(0.0);
                        }
                        
                        if (d_Harten_wavelet_uses_local_tol[vi])
                        {
                            d_Harten_local_means_density.push_back(
                                HAMERS_MAKE_SHARED<pdat::CellVariable<double> >(
                                    d_dim,
                                    "Harten local mean of density at level " +
                                        std::to_string(li),
                                    1));
                        }
                    }
                    
                    if (d_Harten_wavelet_uses_alpha_tol[vi])
                    {
                        d_Harten_Lipschitz_exponent_density =
                            HAMERS_SHARED_PTR<pdat::CellVariable<double> > (
                                new pdat::CellVariable<double>(
                                    d_dim,
                                    "Harten Lipschitz's exponent of density",
                                    1));
                    }
                }
                else if (variable_key == "TOTAL_ENERGY")
                {
                    for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                    {
                        d_Harten_wavelet_coeffs_total_energy.push_back(
                            HAMERS_MAKE_SHARED<pdat::CellVariable<double> >(
                                d_dim,
                                "Harten wavelet coefficient of total energy at level " +
                                    std::to_string(li),
                                1));
                        
                        if (d_Harten_wavelet_uses_global_tol[vi])
                        {
                            d_Harten_wavelet_coeffs_maxs_total_energy.push_back(0.0);
                        }
                        
                        if (d_Harten_wavelet_uses_local_tol[vi])
                        {
                            d_Harten_local_means_total_energy.push_back(
                                HAMERS_MAKE_SHARED<pdat::CellVariable<double> >(
                                    d_dim,
                                    "Harten local mean of total energy at level " +
                                        std::to_string(li),
                                    1));
                        }
                    }
                    
                    if (d_Harten_wavelet_uses_alpha_tol[vi])
                    {
                        d_Harten_Lipschitz_exponent_total_energy =
                            HAMERS_SHARED_PTR<pdat::CellVariable<double> > (
                                new pdat::CellVariable<double>(
                                    d_dim,
                                    "Harten Lipschitz's exponent of total energy",
                                    1));
                    }
                }
                else if (variable_key == "PRESSURE")
                {
                    for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                    {
                        d_Harten_wavelet_coeffs_pressure.push_back(
                            HAMERS_MAKE_SHARED<pdat::CellVariable<double> >(
                                d_dim,
                                "Harten wavelet coefficient of pressure at level " +
                                    std::to_string(li),
                                1));
                        
                        if (d_Harten_wavelet_uses_global_tol[vi])
                        {
                            d_Harten_wavelet_coeffs_maxs_pressure.push_back(0.0);
                        }
                        
                        if (d_Harten_wavelet_uses_local_tol[vi])
                        {
                            d_Harten_local_means_pressure.push_back(
                                HAMERS_MAKE_SHARED<pdat::CellVariable<double> >(
                                    d_dim,
                                    "Harten local mean of pressure at level " +
                                        std::to_string(li),
                                    1));
                        }
                    }
                    
                    if (d_Harten_wavelet_uses_alpha_tol[vi])
                    {
                        d_Harten_Lipschitz_exponent_pressure =
                            HAMERS_SHARED_PTR<pdat::CellVariable<double> > (
                                new pdat::CellVariable<double>(
                                    d_dim,
                                    "Harten Lipschitz's exponent of pressure",
                                    1));
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
            
            for (int vi = 0; vi < static_cast<int>(d_Harten_wavelet_variables.size()); vi++)
            {
                // Get the key of the current variable.
                std::string variable_key = d_Harten_wavelet_variables[vi];
                
                if (variable_key == "DENSITY")
                {
                    for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                    {
                        integrator->registerVariable(
                            d_Harten_wavelet_coeffs_density[li],
                            d_num_multiresolution_ghosts,
                            d_num_multiresolution_ghosts,
                            RungeKuttaLevelIntegrator::TIME_DEP,
                                d_grid_geometry,
                                "NO_COARSEN",
                                "NO_REFINE");
                        
                        if (d_Harten_wavelet_uses_local_tol[vi])
                        {
                            integrator->registerVariable(
                                d_Harten_local_means_density[li],
                                d_num_multiresolution_ghosts,
                                d_num_multiresolution_ghosts,
                                RungeKuttaLevelIntegrator::TEMPORARY,
                                    d_grid_geometry,
                                    "NO_COARSEN",
                                    "NO_REFINE");
                        }
                    }
                    
                    if (d_Harten_wavelet_uses_alpha_tol[vi])
                    {
                        integrator->registerVariable(
                            d_Harten_Lipschitz_exponent_density,
                            d_num_multiresolution_ghosts,
                            d_num_multiresolution_ghosts,
                            RungeKuttaLevelIntegrator::TIME_DEP,
                                d_grid_geometry,
                                "NO_COARSEN",
                                "NO_REFINE");
                    }
                }
                else if (variable_key == "TOTAL_ENERGY")
                {
                    for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                    {
                        integrator->registerVariable(
                            d_Harten_wavelet_coeffs_total_energy[li],
                            d_num_multiresolution_ghosts,
                            d_num_multiresolution_ghosts,
                            RungeKuttaLevelIntegrator::TIME_DEP,
                            d_grid_geometry,
                            "NO_COARSEN",
                            "NO_REFINE");
                        
                        if (d_Harten_wavelet_uses_local_tol[vi])
                        {
                            integrator->registerVariable(
                                d_Harten_local_means_total_energy[li],
                                d_num_multiresolution_ghosts,
                                d_num_multiresolution_ghosts,
                                RungeKuttaLevelIntegrator::TEMPORARY,
                                    d_grid_geometry,
                                    "NO_COARSEN",
                                    "NO_REFINE");
                        }
                    }
                    
                    if (d_Harten_wavelet_uses_alpha_tol[vi])
                    {
                        integrator->registerVariable(
                            d_Harten_Lipschitz_exponent_total_energy,
                            d_num_multiresolution_ghosts,
                            d_num_multiresolution_ghosts,
                            RungeKuttaLevelIntegrator::TIME_DEP,
                                d_grid_geometry,
                                "NO_COARSEN",
                                "NO_REFINE");
                    }
                }
                else if (variable_key == "PRESSURE")
                {
                    for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                    {
                        integrator->registerVariable(
                            d_Harten_wavelet_coeffs_pressure[li],
                            d_num_multiresolution_ghosts,
                            d_num_multiresolution_ghosts,
                            RungeKuttaLevelIntegrator::TIME_DEP,
                            d_grid_geometry,
                            "NO_COARSEN",
                            "NO_REFINE");
                        
                        if (d_Harten_wavelet_uses_local_tol[vi])
                        {
                            integrator->registerVariable(
                                d_Harten_local_means_pressure[li],
                                d_num_multiresolution_ghosts,
                                d_num_multiresolution_ghosts,
                                RungeKuttaLevelIntegrator::TEMPORARY,
                                    d_grid_geometry,
                                    "NO_COARSEN",
                                    "NO_REFINE");
                        }
                    }
                    
                    if (d_Harten_wavelet_uses_alpha_tol[vi])
                    {
                        integrator->registerVariable(
                            d_Harten_Lipschitz_exponent_pressure,
                            d_num_multiresolution_ghosts,
                            d_num_multiresolution_ghosts,
                            RungeKuttaLevelIntegrator::TIME_DEP,
                                d_grid_geometry,
                                "NO_COARSEN",
                                "NO_REFINE");
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
MultiresolutionTagger::registerPlotQuantities(
    const HAMERS_SHARED_PTR<ExtendedVisItDataWriter>& visit_writer,
    const HAMERS_SHARED_PTR<hier::VariableContext>& plot_context)
{
#ifdef HAMERS_PLOTTING_MULTIRESOLUTION_TAGGER
    hier::VariableDatabase* vardb = hier::VariableDatabase::getDatabase();
    
    for (int si = 0; si < static_cast<int>(d_multiresolution_sensors.size()); si++)
    {
        std::string sensor_key = d_multiresolution_sensors[si];
        
        if (sensor_key == "HARTEN_WAVELET")
        {
            for (int vi = 0; vi < static_cast<int>(d_Harten_wavelet_variables.size()); vi++)
            {
                // Get the key of the current variable.
                std::string variable_key = d_Harten_wavelet_variables[vi];
                
                if (variable_key == "DENSITY")
                {
                    for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                    {
                        visit_writer->registerPlotQuantity(
                            "Harten wavelet coefficient of density at level " +
                                std::to_string(li),
                            "SCALAR",
                            vardb->mapVariableAndContextToIndex(
                               d_Harten_wavelet_coeffs_density[li],
                               plot_context));
                    }
                    
                    if (d_Harten_wavelet_uses_alpha_tol)
                    {
                        visit_writer->registerPlotQuantity(
                            "Harten Lipschitz's exponent of density",
                            "SCALAR",
                            vardb->mapVariableAndContextToIndex(
                               d_Harten_Lipschitz_exponent_density,
                               plot_context));
                    }
                }
                else if (variable_key == "TOTAL_ENERGY")
                {
                    for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                    {
                        visit_writer->registerPlotQuantity(
                            "Harten wavelet coefficient of total energy at level " +
                                std::to_string(li),
                            "SCALAR",
                            vardb->mapVariableAndContextToIndex(
                               d_Harten_wavelet_coeffs_total_energy[li],
                               plot_context));
                    }
                    
                    if (d_Harten_wavelet_uses_alpha_tol)
                    {
                        visit_writer->registerPlotQuantity(
                            "Harten Lipschitz's exponent of total energy",
                            "SCALAR",
                            vardb->mapVariableAndContextToIndex(
                               d_Harten_Lipschitz_exponent_total_energy,
                               plot_context));
                    }
                }
                else if (variable_key == "PRESSURE")
                {
                    for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                    {
                        visit_writer->registerPlotQuantity(
                            "Harten wavelet coefficient of pressure at level " +
                                std::to_string(li),
                            "SCALAR",
                            vardb->mapVariableAndContextToIndex(
                               d_Harten_wavelet_coeffs_pressure[li],
                               plot_context));
                    }
                    
                    if (d_Harten_wavelet_uses_alpha_tol)
                    {
                        visit_writer->registerPlotQuantity(
                            "Harten Lipschitz's exponent of pressure",
                            "SCALAR",
                            vardb->mapVariableAndContextToIndex(
                               d_Harten_Lipschitz_exponent_pressure,
                               plot_context));
                    }
                }
            }
        }
    }
#endif
}


/*
 * Print all characteristics of the multiresolution tagger class.
 */
void
MultiresolutionTagger::printClassData(std::ostream& os) const
{
    os << "\nPrint MultiresolutionTagger object..."
       << std::endl;
    
    os << std::endl;
    
    os << "d_multiresolution_sensors = ";
    for (int si = 0; si < static_cast<int>(d_multiresolution_sensors.size()); si++)
    {
        os << "\"" << d_multiresolution_sensors[si] << "\"";
        
        if (si < static_cast<int>(d_multiresolution_sensors.size()) - 1)
        {
            os << ", ";
        }
    }
    os << std::endl;
    
    if (d_wavelet_transfrom_Harten != nullptr)
    {
        os << std::endl;
        os << "d_Harten_wavelet_variables = ";
        for (int vi = 0; vi < static_cast<int>(d_Harten_wavelet_variables.size()); vi++)
        {
            os << "\"" << d_Harten_wavelet_variables[vi] << "\"";
            
            if (vi < static_cast<int>(d_Harten_wavelet_variables.size()) - 1)
            {
                os << ", ";
            }
        }
        os << std::endl;
        
        os << "d_Harten_wavelet_uses_global_tol = ";
        for (int ti = 0; ti < static_cast<int>(d_Harten_wavelet_uses_global_tol.size()); ti++)
        {
            os << "\"" << d_Harten_wavelet_uses_global_tol[ti] << "\"";
            
            if (ti < static_cast<int>(d_Harten_wavelet_uses_global_tol.size()) - 1)
            {
                os << ", ";
            }
        }
        os << std::endl;
        
        os << "d_Harten_wavelet_global_tol = ";
        for (int ti = 0; ti < static_cast<int>(d_Harten_wavelet_global_tol.size()); ti++)
        {
            os << "\"" << d_Harten_wavelet_global_tol[ti] << "\"";
            
            if (ti < static_cast<int>(d_Harten_wavelet_global_tol.size()) - 1)
            {
                os << ", ";
            }
        }
        os << std::endl;
        
        os << "d_Harten_wavelet_uses_local_tol = ";
        for (int ti = 0; ti < static_cast<int>(d_Harten_wavelet_uses_local_tol.size()); ti++)
        {
            os << "\"" << d_Harten_wavelet_uses_local_tol[ti] << "\"";
            
            if (ti < static_cast<int>(d_Harten_wavelet_uses_local_tol.size()) - 1)
            {
                os << ", ";
            }
        }
        os << std::endl;
        
        os << "d_Harten_wavelet_local_tol = ";
        for (int ti = 0; ti < static_cast<int>(d_Harten_wavelet_local_tol.size()); ti++)
        {
            os << "\"" << d_Harten_wavelet_local_tol[ti] << "\"";
            
            if (ti < static_cast<int>(d_Harten_wavelet_local_tol.size()) - 1)
            {
                os << ", ";
            }
        }
        os << std::endl;
        
        os << "d_uses_global_d_Harten_wavelet_uses_alpha_toltol_up = ";
        for (int ti = 0; ti < static_cast<int>(d_Harten_wavelet_uses_alpha_tol.size()); ti++)
        {
            os << "\"" << d_Harten_wavelet_uses_alpha_tol[ti] << "\"";
            
            if (ti < static_cast<int>(d_Harten_wavelet_uses_alpha_tol.size()) - 1)
            {
                os << ", ";
            }
        }
        os << std::endl;
        
        os << "d_Harten_wavelet_alpha_tol = ";
        for (int ti = 0; ti < static_cast<int>(d_Harten_wavelet_alpha_tol.size()); ti++)
        {
            os << "\"" << d_Harten_wavelet_alpha_tol[ti] << "\"";
            
            if (ti < static_cast<int>(d_Harten_wavelet_alpha_tol.size()) - 1)
            {
                os << ", ";
            }
        }
        os << std::endl;
    }
}


/*
 * Put the characteristics of the multiresolution tagger into the restart
 * database.
 */
void
MultiresolutionTagger::putToRestart(
    const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const
{
    if (static_cast<int>(d_multiresolution_sensors.size()) > 0)
    {
        restart_db->putStringVector("d_multiresolution_sensors", d_multiresolution_sensors);
    }
    
    for (int si = 0; si < static_cast<int>(d_multiresolution_sensors.size()); si++)
    {
        if (d_multiresolution_sensors[si] == "HARTEN_WAVELET")
        {
            HAMERS_SHARED_PTR<tbox::Database> sensor_db =
                restart_db->putDatabase("HARTEN_WAVELET");
            
            sensor_db->putInteger("d_Harten_wavelet_num_level",
                d_Harten_wavelet_num_level);
            
            sensor_db->putInteger("d_Harten_wavelet_num_vanishing_moments",
                d_Harten_wavelet_num_vanishing_moments);
            
            sensor_db->putStringVector("d_Harten_wavelet_variables",
                d_Harten_wavelet_variables);
            
            sensor_db->putBoolVector("d_Harten_wavelet_uses_global_tol",
                d_Harten_wavelet_uses_global_tol);
            
            sensor_db->putBoolVector("d_Harten_wavelet_uses_local_tol",
                d_Harten_wavelet_uses_local_tol);
            
            sensor_db->putBoolVector("d_Harten_wavelet_uses_alpha_tol",
                d_Harten_wavelet_uses_alpha_tol);
            
            int num_true = 0;
            
            num_true = std::count(d_Harten_wavelet_uses_global_tol.begin(),
                d_Harten_wavelet_uses_global_tol.end(),
                true);
            if (num_true > 0)
            {
                sensor_db->putDoubleVector("d_Harten_wavelet_global_tol",
                    d_Harten_wavelet_global_tol);
            }
            
            num_true = std::count(d_Harten_wavelet_uses_local_tol.begin(),
                d_Harten_wavelet_uses_local_tol.end(),
                true);
            if (num_true > 0)
            {
                sensor_db->putDoubleVector("d_Harten_wavelet_local_tol",
                    d_Harten_wavelet_local_tol);
            }
            
            num_true = std::count(d_Harten_wavelet_uses_alpha_tol.begin(),
                d_Harten_wavelet_uses_alpha_tol.end(),
                true);
            if (num_true > 0)
            {
                sensor_db->putDoubleVector("d_Harten_wavelet_alpha_tol",
                    d_Harten_wavelet_alpha_tol);
            }
        }
    }
}


/*
 * Compute values of multiresolution sensors on a  patch.
 */
void
MultiresolutionTagger::computeMultiresolutionSensorValuesOnPatch(
    hier::Patch& patch,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    // Loop over multiresolution sensors chosen.
    for (int si = 0; si < static_cast<int>(d_multiresolution_sensors.size()); si++)
    {
        std::string sensor_key = d_multiresolution_sensors[si];
        
        if (sensor_key == "HARTEN_WAVELET")
        {
            for (int vi = 0; vi < static_cast<int>(d_Harten_wavelet_variables.size()); vi++)
            {
                // Get the key of the current variable.
                std::string variable_key = d_Harten_wavelet_variables[vi];
                
                if (variable_key == "DENSITY")
                {
                    /*
                     * Register the patch and density in the flow model and compute the
                     * corresponding cell data.
                     */
                    
                    d_flow_model->registerPatchWithDataContext(patch, data_context);
                    
                    std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                    
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(
                            "DENSITY", d_num_multiresolution_ghosts));
                    
                    d_flow_model->registerDerivedVariables(num_subghosts_of_data);
                    
                    d_flow_model->allocateMemoryForDerivedCellData();
                    
                    d_flow_model->computeDerivedCellData();
                    
                    // Get the pointer to density data inside the flow model.
                    HAMERS_SHARED_PTR<pdat::CellData<double> > data_density = d_flow_model->getCellData("DENSITY");
                    
                    // Get the wavelet coefficients.
                    std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > wavelet_coeffs;
                    for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                    {
                        wavelet_coeffs.push_back(
                            HAMERS_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                                patch.getPatchData(
                                    d_Harten_wavelet_coeffs_density[li],
                                    data_context)));
                    }
                    
                    // Compute the wavelet coefficients.
                    if (d_Harten_wavelet_uses_local_tol[vi])
                    {
                        // Get the local means.
                        std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > variable_local_means;
                        for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                        {
                            variable_local_means.push_back(
                                HAMERS_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                                    patch.getPatchData(
                                        d_Harten_local_means_density[li],
                                        data_context)));
                        }
                        
                        d_wavelet_transfrom_Harten->computeWaveletCoefficientsWithVariableLocalMeans(
                            wavelet_coeffs,
                            variable_local_means,
                            data_density,
                            patch);
                    }
                    else
                    {
                        d_wavelet_transfrom_Harten->computeWaveletCoefficients(
                            wavelet_coeffs,
                            data_density,
                            patch);
                    }
                    
                    /*
                     * Unregister the patch and data of all registered derived cell variables in
                     * the flow model.
                     */
                    
                    d_flow_model->unregisterPatch();
                    
                }
                else if (variable_key ==  "TOTAL_ENERGY")
                {
                    /*
                     * Register the patch and total energy in the flow model and compute the
                     * corresponding cell data.
                     */
                    
                    d_flow_model->registerPatchWithDataContext(patch, data_context);
                    
                    std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                    
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(
                            "TOTAL_ENERGY", d_num_multiresolution_ghosts));
                    
                    d_flow_model->registerDerivedVariables(num_subghosts_of_data);
                    
                    d_flow_model->allocateMemoryForDerivedCellData();
                    
                    d_flow_model->computeDerivedCellData();
                    
                    // Get the pointer to total energy data inside the flow model.
                    HAMERS_SHARED_PTR<pdat::CellData<double> > data_total_energy =
                        d_flow_model->getCellData("TOTAL_ENERGY");
                    
                    // Get the wavelet coefficients.
                    std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > wavelet_coeffs;
                    for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                    {
                        wavelet_coeffs.push_back(
                            HAMERS_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                                patch.getPatchData(
                                    d_Harten_wavelet_coeffs_total_energy[li],
                                    data_context)));
                    }
                    
                    // Compute the wavelet coefficients.
                    if (d_Harten_wavelet_uses_local_tol[vi])
                    {
                        // Get the local means.
                        std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > variable_local_means;
                        for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                        {
                            variable_local_means.push_back(
                                HAMERS_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                                    patch.getPatchData(
                                        d_Harten_local_means_total_energy[li],
                                        data_context)));
                        }
                        
                        d_wavelet_transfrom_Harten->computeWaveletCoefficientsWithVariableLocalMeans(
                            wavelet_coeffs,
                            variable_local_means,
                            data_total_energy,
                            patch);
                    }
                    else
                    {
                        d_wavelet_transfrom_Harten->computeWaveletCoefficients(
                            wavelet_coeffs,
                            data_total_energy,
                            patch);
                    }
                    
                    /*
                     * Unregister the patch and data of all registered derived cell variables in
                     * the flow model.
                     */
                    
                    d_flow_model->unregisterPatch();
                    
                }
                else if (variable_key == "PRESSURE")
                {
                    /*
                     * Register the patch and pressure in the flow model and compute the
                     * corresponding cell data.
                     */
                    
                    d_flow_model->registerPatchWithDataContext(patch, data_context);
                    
                    std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                    
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(
                            "PRESSURE", d_num_multiresolution_ghosts));
                    
                    d_flow_model->registerDerivedVariables(num_subghosts_of_data);
                    
                    d_flow_model->allocateMemoryForDerivedCellData();
                    
                    d_flow_model->computeDerivedCellData();
                    
                    // Get the pointer to pressure data inside the flow model.
                    HAMERS_SHARED_PTR<pdat::CellData<double> > data_pressure = d_flow_model->getCellData("PRESSURE");
                    
                    // Get the wavelet coefficients.
                    std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > wavelet_coeffs;
                    for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                    {
                        wavelet_coeffs.push_back(
                                HAMERS_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                                    patch.getPatchData(
                                        d_Harten_wavelet_coeffs_pressure[li],
                                        data_context)));
                    }
                    
                    // Compute the wavelet coefficients.
                    if (d_Harten_wavelet_uses_local_tol[vi])
                    {
                        // Get the local means.
                        std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > variable_local_means;
                        for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                        {
                            variable_local_means.push_back(
                                    HAMERS_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                                        patch.getPatchData(
                                            d_Harten_local_means_pressure[li],
                                            data_context)));
                        }
                        
                        d_wavelet_transfrom_Harten->computeWaveletCoefficientsWithVariableLocalMeans(
                            wavelet_coeffs,
                            variable_local_means,
                            data_pressure,
                            patch);
                    }
                    else
                    {
                        d_wavelet_transfrom_Harten->computeWaveletCoefficients(
                            wavelet_coeffs,
                            data_pressure,
                            patch);
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
 * multiresolution sensors at a given patch level.
 */
void
MultiresolutionTagger::getSensorValueStatistics(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const int level_number,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    math::HierarchyCellDataOpsReal<double> cell_double_operator(patch_hierarchy, level_number, level_number);
    
    hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();
    
    for (int si = 0; si < static_cast<int>(d_multiresolution_sensors.size()); si++)
    {
        std::string sensor_key = d_multiresolution_sensors[si];
        
        if (sensor_key == "HARTEN_WAVELET")
        {
            for (int vi = 0; vi < static_cast<int>(d_Harten_wavelet_variables.size()); vi++)
            {
                // Only get the statistics if global tolerances is used.
                if (d_Harten_wavelet_uses_global_tol[vi])
                {
                    // Get the key of the current variable.
                    std::string variable_key = d_Harten_wavelet_variables[vi];
                    
                    if (variable_key == "DENSITY")
                    {
                        for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                        {
                            const int w_rho_id = variable_db->mapVariableAndContextToIndex(
                                d_Harten_wavelet_coeffs_density[li],
                                data_context);
                            
                            double w_max_rho_local = cell_double_operator.max(w_rho_id);
                            d_Harten_wavelet_coeffs_maxs_density[li] = 0.0;
                            
                            mpi.Allreduce(
                                &w_max_rho_local,
                                &d_Harten_wavelet_coeffs_maxs_density[li],
                                1,
                                MPI_DOUBLE,
                                MPI_MAX);
                        }
                    }
                    else if (variable_key == "TOTAL_ENERGY")
                    {
                        for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                        {
                            const int w_E_id = variable_db->mapVariableAndContextToIndex(
                                d_Harten_wavelet_coeffs_total_energy[li],
                                data_context);
                            
                            double w_E_max_local = cell_double_operator.max(w_E_id);
                            d_Harten_wavelet_coeffs_maxs_total_energy[li] = 0.0;
                            
                            mpi.Allreduce(
                                &w_E_max_local,
                                &d_Harten_wavelet_coeffs_maxs_total_energy[li],
                                1,
                                MPI_DOUBLE,
                                MPI_MAX);
                        }
                    }
                    else if (variable_key == "PRESSURE")
                    {
                        for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                        {
                            const int w_p_id = variable_db->mapVariableAndContextToIndex(
                                d_Harten_wavelet_coeffs_pressure[li],
                                data_context);
                            
                            double w_p_max_local = cell_double_operator.max(w_p_id);
                            d_Harten_wavelet_coeffs_maxs_pressure[li] = 0.0;
                            
                            mpi.Allreduce(
                                &w_p_max_local,
                                &d_Harten_wavelet_coeffs_maxs_pressure[li],
                                1,
                                MPI_DOUBLE,
                                MPI_MAX);
                        }
                    }
                }
            }
        }
    }
}


/*
 * Tag cells on a patch for refinement using multiresolution sensors.
 */
void
MultiresolutionTagger::tagCellsOnPatch(
   hier::Patch& patch,
   const HAMERS_SHARED_PTR<pdat::CellData<int> >& tags,
   const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    for (int si = 0; si < static_cast<int>(d_multiresolution_sensors.size()); si++)
    {
        std::string sensor_key = d_multiresolution_sensors[si];
        
        if (sensor_key == "HARTEN_WAVELET")
        {
            int count_global_tol = 0;
            int count_local_tol = 0;
            int count_alpha_tol = 0;
            
            // Loop over the variables for wavelet analysis.
            for (int vi = 0; vi < static_cast<int>(d_Harten_wavelet_variables.size()); vi++)
            {
                // Get the key of the current variable.
                std::string variable_key = d_Harten_wavelet_variables[vi];
                
                // Get tolerance settings for the current variable.
                bool uses_global_tol = d_Harten_wavelet_uses_global_tol[vi];
                bool uses_local_tol = d_Harten_wavelet_uses_local_tol[vi];
                bool uses_alpha_tol = d_Harten_wavelet_uses_alpha_tol[vi];
                
                double global_tol = 0.0;
                double local_tol = 0.0;
                double alpha_tol = 0.0;
                
                if (uses_global_tol)
                {
                    global_tol = d_Harten_wavelet_global_tol[count_global_tol];
                    count_global_tol++;
                }
                
                if (uses_local_tol)
                {
                    local_tol = d_Harten_wavelet_local_tol[count_local_tol];
                    count_local_tol++;
                }
                
                if (uses_alpha_tol)
                {
                    alpha_tol = d_Harten_wavelet_alpha_tol[count_alpha_tol];
                    count_alpha_tol++;
                }
                
                /*
                 * Get the wavelet coefficients and the statistics of the wavelet coefficients
                 * at different levels and tag cells.
                 */
                
                if (variable_key == "DENSITY")
                {
                    std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > wavelet_coeffs;
                    
                    for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                    {
                        wavelet_coeffs.push_back(
                            HAMERS_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                                patch.getPatchData(
                                    d_Harten_wavelet_coeffs_density[li],
                                    data_context)));
                    }
                    
                    std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > variable_local_means;
                    if (uses_local_tol)
                    {
                        // Get the local means.
                        for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                        {
                            variable_local_means.push_back(
                                HAMERS_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                                    patch.getPatchData(
                                        d_Harten_local_means_density[li],
                                        data_context)));
                        }
                    }
                    
                    HAMERS_SHARED_PTR<pdat::CellData<double> > Lipschitz_exponent;
                    if (uses_alpha_tol)
                    {
                        Lipschitz_exponent = HAMERS_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                            patch.getPatchData(
                                d_Harten_Lipschitz_exponent_density,
                                data_context));
                    }
                    
                    tagCellsOnPatchWithWaveletSensor(
                        patch,
                        tags,
                        wavelet_coeffs,
                        d_Harten_wavelet_coeffs_maxs_density,
                        variable_local_means,
                        Lipschitz_exponent,
                        sensor_key,
                        uses_global_tol,
                        uses_local_tol,
                        uses_alpha_tol,
                        global_tol,
                        local_tol,
                        alpha_tol);
                }
                else if (variable_key == "TOTAL_ENERGY")
                {
                    std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > wavelet_coeffs;
                    
                    for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                    {
                        wavelet_coeffs.push_back(
                                HAMERS_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                                    patch.getPatchData(
                                        d_Harten_wavelet_coeffs_total_energy[li],
                                        data_context)));
                    }
                    
                    std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > variable_local_means;
                    if (uses_local_tol)
                    {
                        // Get the local means.
                        for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                        {
                            variable_local_means.push_back(
                                    HAMERS_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                                        patch.getPatchData(
                                            d_Harten_local_means_total_energy[li],
                                            data_context)));
                        }
                    }
                    
                    HAMERS_SHARED_PTR<pdat::CellData<double> > Lipschitz_exponent;
                    if (uses_alpha_tol)
                    {
                        Lipschitz_exponent = HAMERS_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                            patch.getPatchData(
                                d_Harten_Lipschitz_exponent_total_energy,
                                data_context));
                    }
                    
                    tagCellsOnPatchWithWaveletSensor(
                        patch,
                        tags,
                        wavelet_coeffs,
                        d_Harten_wavelet_coeffs_maxs_total_energy,
                        variable_local_means,
                        Lipschitz_exponent,
                        sensor_key,
                        uses_global_tol,
                        uses_local_tol,
                        uses_alpha_tol,
                        global_tol,
                        local_tol,
                        alpha_tol);
                }
                else if (variable_key == "PRESSURE")
                {
                    std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > wavelet_coeffs;
                    
                    for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                    {
                        wavelet_coeffs.push_back(
                                HAMERS_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                                    patch.getPatchData(
                                        d_Harten_wavelet_coeffs_pressure[li],
                                        data_context)));
                    }
                    
                    std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > variable_local_means;
                    if (uses_local_tol)
                    {
                        // Get the local means.
                        for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                        {
                            variable_local_means.push_back(
                                    HAMERS_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                                        patch.getPatchData(
                                            d_Harten_local_means_pressure[li],
                                            data_context)));
                        }
                    }
                    
                    HAMERS_SHARED_PTR<pdat::CellData<double> > Lipschitz_exponent;
                    if (uses_alpha_tol)
                    {
                        Lipschitz_exponent = HAMERS_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                            patch.getPatchData(
                                d_Harten_Lipschitz_exponent_pressure,
                                data_context));
                    }
                    
                    tagCellsOnPatchWithWaveletSensor(
                        patch,
                        tags,
                        wavelet_coeffs,
                        d_Harten_wavelet_coeffs_maxs_pressure,
                        variable_local_means,
                        Lipschitz_exponent,
                        sensor_key,
                        uses_global_tol,
                        uses_local_tol,
                        uses_alpha_tol,
                        global_tol,
                        local_tol,
                        alpha_tol);
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
            } // Loop over variables.
        }
    }
}


/*
 * Compute the Lipschitz's exponent on a patch. There are two steps:
 * 1. Find the maximum wavelet coefficients in the domain of dependence.
 * 2. Compute Lipschitz's exponent.
 */
void
MultiresolutionTagger::computeLipschitzExponentOnPatch(
    hier::Patch& patch,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& Lipschitz_exponent,
    const std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& wavelet_coeffs,
    const std::string& sensor_key)
{
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    if (sensor_key == "HARTEN_WAVELET")
    {
        /*
         * 1. Find the maximum wavelet coefficients in domain of dependence.
         */
        
        // Get the pointers to the wavelet coefficients.
        std::vector<double*> w;
        for (int li = 0; li < d_Harten_wavelet_num_level; li++)
        {
            w.push_back(wavelet_coeffs[li]->getPointer(0));
        }
        
        // Get the number of ghost cells and dimensions of box that covers interior of patch plus
        // ghost cells.
        const hier::IntVector num_ghosts_wavelet_coeffs = wavelet_coeffs[0]->getGhostCellWidth();
        const hier::IntVector ghostcell_dims_wavelet_coeffs = wavelet_coeffs[0]->getGhostBox().numberCells();
        
        // Create a vector of maximum wavelet coefficients in domain of dependence.
        std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > wavelet_coeffs_local_max;
        for (int li = 0; li < d_Harten_wavelet_num_level; li++)
        {
            wavelet_coeffs_local_max.push_back(HAMERS_MAKE_SHARED<pdat::CellData<double> >(
                interior_box, 1, num_ghosts_wavelet_coeffs));
        }
        
        // Get the pointers to the maximum wavelet coefficients in the domain of dependence.
        std::vector<double*> r;
        for (int li = 0; li < d_Harten_wavelet_num_level; li++)
        {
            r.push_back(wavelet_coeffs_local_max[li]->getPointer(0));
        }
        
        // Get the stencil width of the wavelet transform.
        const int p = d_wavelet_transfrom_Harten->getLeftStencilWidth();
        const int q = d_wavelet_transfrom_Harten->getRightStencilWidth();
        
        for (int li = 0; li < d_Harten_wavelet_num_level; li++)
        {
            if (d_dim == tbox::Dimension(1))
            {
                const int interior_dim_0 = interior_dims[0];
                
                const int num_ghosts_0_wavelet_coeffs = num_ghosts_wavelet_coeffs[0];
                
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute linear index.
                    const int idx = i + num_ghosts_0_wavelet_coeffs;
                    
                    // Find the maximum wavelet coefficient over the domain of dependence.
                    r[li][idx] = 0.0;
                    for (int ii = -p*pow(2, li+1); ii <= q*pow(2, li+1); ii++)
                    {
                        // Compute the index.
                        const int idx_s = i + ii + num_ghosts_0_wavelet_coeffs;
                        
                        r[li][idx] = fmax(r[li][idx], w[li][idx_s]);
                    }
                }
            }
            else if (d_dim == tbox::Dimension(2))
            {
                const int interior_dim_0 = interior_dims[0];
                const int interior_dim_1 = interior_dims[1];
                
                const int num_ghosts_0_wavelet_coeffs = num_ghosts_wavelet_coeffs[0];
                const int num_ghosts_1_wavelet_coeffs = num_ghosts_wavelet_coeffs[1];
                const int ghostcell_dim_0_wavelet_coeffs = ghostcell_dims_wavelet_coeffs[0];
                
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute linear index.
                        const int idx = (i + num_ghosts_0_wavelet_coeffs) +
                            (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs;
                        
                        // Find the maximum wavelet coefficient over the domain of dependence.
                        r[li][idx] = 0.0;
                        for (int ii = -p*pow(2, li+1); ii <= q*pow(2, li+1); ii++)
                        {
                            // Compute the index.
                            const int idx_s = (i + ii + num_ghosts_0_wavelet_coeffs) +
                                (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs;
                            
                            r[li][idx] = fmax(r[li][idx], w[li][idx_s]);
                        }
                        for (int jj = -p*pow(2, li+1); jj <= q*pow(2, li+1); jj++)
                        {
                            // Compute the index.
                            const int idx_s = (i + num_ghosts_0_wavelet_coeffs) +
                                (j + jj + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs;
                            
                            r[li][idx] = fmax(r[li][idx], w[li][idx_s]);
                        }
                    }
                }
            }
            else if (d_dim == tbox::Dimension(3))
            {
                const int interior_dim_0 = interior_dims[0];
                const int interior_dim_1 = interior_dims[1];
                const int interior_dim_2 = interior_dims[2];
                
                const int num_ghosts_0_wavelet_coeffs = num_ghosts_wavelet_coeffs[0];
                const int num_ghosts_1_wavelet_coeffs = num_ghosts_wavelet_coeffs[1];
                const int num_ghosts_2_wavelet_coeffs = num_ghosts_wavelet_coeffs[2];
                const int ghostcell_dim_0_wavelet_coeffs = ghostcell_dims_wavelet_coeffs[0];
                const int ghostcell_dim_1_wavelet_coeffs = ghostcell_dims_wavelet_coeffs[1];
                
                for (int k = 0; k < interior_dim_2; k++)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute linear index.
                            const int idx = (i + num_ghosts_0_wavelet_coeffs) +
                                (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                (k + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                    ghostcell_dim_1_wavelet_coeffs;
                            
                            // Find the maximum wavelet coefficient over the domain of dependence.
                            r[li][idx] = 0.0;
                            for (int ii = -p*pow(2, li+1); ii <= q*pow(2, li+1); ii++)
                            {
                                // Compute the index.
                                const int idx_s = (i + ii + num_ghosts_0_wavelet_coeffs) +
                                    (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                    (k + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                        ghostcell_dim_1_wavelet_coeffs;
                                
                                r[li][idx] = fmax(r[li][idx], w[li][idx_s]);
                            }
                            for (int jj = -p*pow(2, li+1); jj <= q*pow(2, li+1); jj++)
                            {
                                // Compute the index.
                                const int idx_s = (i + num_ghosts_0_wavelet_coeffs) +
                                    (j + jj + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                    (k + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                        ghostcell_dim_1_wavelet_coeffs;
                                
                                r[li][idx] = fmax(r[li][idx], w[li][idx_s]);
                            }
                            for (int kk = -p*pow(2, li+1); kk <=q*pow(2, li+1); kk++)
                            {
                                // Compute the index.
                                const int idx_s = (i + num_ghosts_0_wavelet_coeffs) +
                                    (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                    (k + kk + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                        ghostcell_dim_1_wavelet_coeffs;
                                
                                r[li][idx] = fmax(r[li][idx], w[li][idx_s]);
                            }
                        }
                    }
                }
            }
        }
        
        /*
         * 2. Compute Lipschitz's exponent.
         */
        
        double* alpha = Lipschitz_exponent->getPointer(0);
        
        switch (d_Harten_wavelet_num_level)
        {
            case 2:
            {
                if (d_dim == tbox::Dimension(1))
                {
                    const int interior_dim_0 = interior_dims[0];
                    
                    const int num_ghosts_0_wavelet_coeffs = num_ghosts_wavelet_coeffs[0];
                    
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute linear index.
                        const int idx = i + num_ghosts_0_wavelet_coeffs;
                        
                        if ((r[0][idx] > 1.0e-8) &&
                            (r[1][idx] > 1.0e-8))
                        {
                            alpha[idx] = fmin(
                                log2(r[1][idx]/r[0][idx]),
                                (double) d_Harten_wavelet_num_vanishing_moments);
                        }
                        else
                        {
                            alpha[idx] = (double) d_Harten_wavelet_num_vanishing_moments;
                        }
                    }
                }
                else if (d_dim == tbox::Dimension(2))
                {
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    
                    const int num_ghosts_0_wavelet_coeffs = num_ghosts_wavelet_coeffs[0];
                    const int num_ghosts_1_wavelet_coeffs = num_ghosts_wavelet_coeffs[1];
                    const int ghostcell_dim_0_wavelet_coeffs = ghostcell_dims_wavelet_coeffs[0];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute linear index.
                            const int idx = (i + num_ghosts_0_wavelet_coeffs) +
                                (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs;
                            
                            if ((r[0][idx] > 1.0e-8) &&
                                (r[1][idx] > 1.0e-8))
                            {
                                alpha[idx] = fmin(
                                    log2(r[1][idx]/r[0][idx]),
                                    (double) d_Harten_wavelet_num_vanishing_moments);
                            }
                            else
                            {
                                alpha[idx] = (double) d_Harten_wavelet_num_vanishing_moments;
                            }
                        }
                    }
                }
                else if (d_dim == tbox::Dimension(3))
                {
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    const int interior_dim_2 = interior_dims[2];
                    
                    const int num_ghosts_0_wavelet_coeffs = num_ghosts_wavelet_coeffs[0];
                    const int num_ghosts_1_wavelet_coeffs = num_ghosts_wavelet_coeffs[1];
                    const int num_ghosts_2_wavelet_coeffs = num_ghosts_wavelet_coeffs[2];
                    const int ghostcell_dim_0_wavelet_coeffs = ghostcell_dims_wavelet_coeffs[0];
                    const int ghostcell_dim_1_wavelet_coeffs = ghostcell_dims_wavelet_coeffs[1];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
#ifdef HAMERS_ENABLE_SIMD
                            #pragma omp simd
#endif
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                // Compute linear index.
                                const int idx = (i + num_ghosts_0_wavelet_coeffs) +
                                    (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                    (k + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                        ghostcell_dim_1_wavelet_coeffs;
                                
                                if ((r[0][idx] > 1.0e-8) &&
                                    (r[1][idx] > 1.0e-8))
                                {
                                    alpha[idx] = fmin(
                                        log2(r[1][idx]/r[0][idx]),
                                        (double) d_Harten_wavelet_num_vanishing_moments);
                                }
                                else
                                {
                                    alpha[idx] = (double) d_Harten_wavelet_num_vanishing_moments;
                                }
                            }
                        }
                    }
                }
                
                break;
            }
            case 3:
            {
                if (d_dim == tbox::Dimension(1))
                {
                    const int interior_dim_0 = interior_dims[0];
                    
                    const int num_ghosts_0_wavelet_coeffs = num_ghosts_wavelet_coeffs[0];
                    
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute linear index.
                        const int idx = i + num_ghosts_0_wavelet_coeffs;
                        
                        if ((r[0][idx] > 1.0e-8) &&
                            (r[1][idx] > 1.0e-8) &&
                            (r[2][idx] > 1.0e-8))
                        {
                            alpha[idx] = fmin(
                                0.5*log2(r[2][idx]/r[0][idx]),
                                (double) d_Harten_wavelet_num_vanishing_moments);
                        }
                        else
                        {
                            alpha[idx] = (double) d_Harten_wavelet_num_vanishing_moments;
                        }
                    }
                }
                else if (d_dim == tbox::Dimension(2))
                {
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    
                    const int num_ghosts_0_wavelet_coeffs = num_ghosts_wavelet_coeffs[0];
                    const int num_ghosts_1_wavelet_coeffs = num_ghosts_wavelet_coeffs[1];
                    const int ghostcell_dim_0_wavelet_coeffs = ghostcell_dims_wavelet_coeffs[0];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute linear index.
                            const int idx = (i + num_ghosts_0_wavelet_coeffs) +
                                (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs;
                            
                            if ((r[0][idx] > 1.0e-8) &&
                                (r[1][idx] > 1.0e-8) &&
                                (r[2][idx] > 1.0e-8))
                            {
                                alpha[idx] = fmin(
                                    0.5*log2(r[2][idx]/r[0][idx]),
                                    (double) d_Harten_wavelet_num_vanishing_moments);
                            }
                            else
                            {
                                alpha[idx] = (double) d_Harten_wavelet_num_vanishing_moments;
                            }
                        }
                    }
                }
                else if (d_dim == tbox::Dimension(3))
                {
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    const int interior_dim_2 = interior_dims[2];
                    
                    const int num_ghosts_0_wavelet_coeffs = num_ghosts_wavelet_coeffs[0];
                    const int num_ghosts_1_wavelet_coeffs = num_ghosts_wavelet_coeffs[1];
                    const int num_ghosts_2_wavelet_coeffs = num_ghosts_wavelet_coeffs[2];
                    const int ghostcell_dim_0_wavelet_coeffs = ghostcell_dims_wavelet_coeffs[0];
                    const int ghostcell_dim_1_wavelet_coeffs = ghostcell_dims_wavelet_coeffs[1];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
#ifdef HAMERS_ENABLE_SIMD
                            #pragma omp simd
#endif
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                // Compute linear index.
                                const int idx = (i + num_ghosts_0_wavelet_coeffs) +
                                    (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                    (k + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                        ghostcell_dim_1_wavelet_coeffs;
                                
                                if ((r[0][idx] > 1.0e-8) &&
                                    (r[1][idx] > 1.0e-8) &&
                                    (r[2][idx] > 1.0e-8))
                                {
                                    alpha[idx] = fmin(
                                        0.5*log2(r[2][idx]/r[0][idx]),
                                        (double) d_Harten_wavelet_num_vanishing_moments);
                                }
                                else
                                {
                                    alpha[idx] = (double) d_Harten_wavelet_num_vanishing_moments;
                                }
                            }
                        }
                    }
                }
                
                break;
            }
            case 4:
            {
                if (d_dim == tbox::Dimension(1))
                {
                    const int interior_dim_0 = interior_dims[0];
                    
                    const int num_ghosts_0_wavelet_coeffs = num_ghosts_wavelet_coeffs[0];
                    
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute linear index.
                        const int idx = i + num_ghosts_0_wavelet_coeffs;
                        
                        if ((r[0][idx] > 1.0e-8) &&
                            (r[1][idx] > 1.0e-8) &&
                            (r[2][idx] > 1.0e-8) &&
                            (r[3][idx] > 1.0e-8))
                        {
                            alpha[idx] = fmin(
                                (3.0/10.0*log2(r[3][idx]/r[0][idx]) +
                                    1.0/10.0*log2(r[2][idx]/r[1][idx])),
                                (double) d_Harten_wavelet_num_vanishing_moments);
                        }
                        else
                        {
                            alpha[idx] = (double) d_Harten_wavelet_num_vanishing_moments;
                        }
                    }
                }
                else if (d_dim == tbox::Dimension(2))
                {
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    
                    const int num_ghosts_0_wavelet_coeffs = num_ghosts_wavelet_coeffs[0];
                    const int num_ghosts_1_wavelet_coeffs = num_ghosts_wavelet_coeffs[1];
                    const int ghostcell_dim_0_wavelet_coeffs = ghostcell_dims_wavelet_coeffs[0];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute linear index.
                            const int idx = (i + num_ghosts_0_wavelet_coeffs) +
                                (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs;
                            
                            if ((r[0][idx] > 1.0e-8) &&
                                (r[1][idx] > 1.0e-8) &&
                                (r[2][idx] > 1.0e-8) &&
                                (r[3][idx] > 1.0e-8))
                            {
                                alpha[idx] = fmin(
                                    (3.0/10.0*log2(r[3][idx]/r[0][idx]) +
                                        1.0/10.0*log2(r[2][idx]/r[1][idx])),
                                    (double) d_Harten_wavelet_num_vanishing_moments);
                            }
                            else
                            {
                                alpha[idx] = (double) d_Harten_wavelet_num_vanishing_moments;
                            }
                        }
                    }
                }
                else if (d_dim == tbox::Dimension(3))
                {
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    const int interior_dim_2 = interior_dims[2];
                    
                    const int num_ghosts_0_wavelet_coeffs = num_ghosts_wavelet_coeffs[0];
                    const int num_ghosts_1_wavelet_coeffs = num_ghosts_wavelet_coeffs[1];
                    const int num_ghosts_2_wavelet_coeffs = num_ghosts_wavelet_coeffs[2];
                    const int ghostcell_dim_0_wavelet_coeffs = ghostcell_dims_wavelet_coeffs[0];
                    const int ghostcell_dim_1_wavelet_coeffs = ghostcell_dims_wavelet_coeffs[1];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
#ifdef HAMERS_ENABLE_SIMD
                            #pragma omp simd
#endif
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                // Compute linear index.
                                const int idx = (i + num_ghosts_0_wavelet_coeffs) +
                                    (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                    (k + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                        ghostcell_dim_1_wavelet_coeffs;
                                
                                if ((r[0][idx] > 1.0e-8) &&
                                    (r[1][idx] > 1.0e-8) &&
                                    (r[2][idx] > 1.0e-8) &&
                                    (r[3][idx] > 1.0e-8))
                                {
                                    alpha[idx] = fmin(
                                        (3.0/10.0*log2(r[3][idx]/r[0][idx]) +
                                            1.0/10.0*log2(r[2][idx]/r[1][idx])),
                                        (double) d_Harten_wavelet_num_vanishing_moments);
                                }
                                else
                                {
                                    alpha[idx] = (double) d_Harten_wavelet_num_vanishing_moments;
                                }
                            }
                        }
                    }
                }
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "d_Harten_wavelet_num_level = "
                    << d_Harten_wavelet_num_level
                    << " not yet implemented."
                    << std::endl);
            }
        }
    }
}


/*
 * Tag cells on a patch using wavelet sensor with the combination of three possible criteria:
 * 1. When ratio between wavelet coefficient and global maximum at any level is greater than the tolerance.
 * 2. When ratio between wavelet coefficient and local mean at any level is greater than the tolerance.
 * 3. When the Lipschitz's exponent is smaller than the tolerance.
 */
void
MultiresolutionTagger::tagCellsOnPatchWithWaveletSensor(
    hier::Patch& patch,
    const HAMERS_SHARED_PTR<pdat::CellData<int> >& tags,
    const std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& wavelet_coeffs,
    const std::vector<double>& wavelet_coeffs_maxs,
    const std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& variable_local_means,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& Lipschitz_exponent,
    const std::string& sensor_key,
    const bool uses_global_tol,
    const bool uses_local_tol,
    const bool uses_alpha_tol,
    const double global_tol,
    const double local_tol,
    const double alpha_tol)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(tags->getGhostCellWidth() == hier::IntVector::getZero(d_dim));
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the number of ghost cells and dimensions of box that covers interior of patch plus
    // ghost cells.
    const hier::IntVector num_ghosts_wavelet_coeffs = wavelet_coeffs[0]->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_wavelet_coeffs = wavelet_coeffs[0]->getGhostBox().numberCells();
    
    // Allocate temporary patch data.
    HAMERS_SHARED_PTR<pdat::CellData<int> > tags_multiresolution_tagger(
        new pdat::CellData<int>(interior_box, d_dim.getValue(), hier::IntVector::getZero(d_dim)));
    
    tags_multiresolution_tagger->fillAll(1);
    
    // Get the pointers to the tags.
    int* tag_ptr_multiresolution_tagger = tags_multiresolution_tagger->getPointer(0);
    int* tag_ptr = tags->getPointer(0);
    
    // Get the pointers to the wavelet coefficients.
    std::vector<double*> w;
    for (int li = 0; li < d_Harten_wavelet_num_level; li++)
    {
        w.push_back(wavelet_coeffs[li]->getPointer(0));
    }
    
    // Get the pointers to the variable local means.
    std::vector<double*> u_mean;
    if (uses_local_tol)
    {
        for (int li = 0; li < d_Harten_wavelet_num_level; li++)
        {
            u_mean.push_back(variable_local_means[li]->getPointer(0));
        }
    }
    
    // Declare pointer to the Lipschitz's exponent.
    double* alpha = NULL;    
    
    if (uses_alpha_tol)
    {
        alpha = Lipschitz_exponent->getPointer(0);
        
        computeLipschitzExponentOnPatch(
            patch,
            Lipschitz_exponent,
            wavelet_coeffs,
            sensor_key);
    }
    
    if (d_dim == tbox::Dimension(1))
    {
        const int interior_dim_0 = interior_dims[0];
        
        const int num_ghosts_0_wavelet_coeffs = num_ghosts_wavelet_coeffs[0];
        
        if (uses_global_tol)
        {
            // Allocate temporary patch data.
            HAMERS_SHARED_PTR<pdat::CellData<int> > tags_global_tol(
                new pdat::CellData<int>(interior_box, d_dim.getValue(), hier::IntVector::getZero(d_dim)));
            
            tags_global_tol->fillAll(0);
            
            // Get the pointer to the tags.
            int* tag_ptr_global_tol = tags_global_tol->getPointer(0);
            
            for (int li = 0; li < d_Harten_wavelet_num_level; li++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx = i + num_ghosts_0_wavelet_coeffs;
                    const int idx_nghost = i;
                    
                    if (w[li][idx]/(wavelet_coeffs_maxs[li] + EPSILON) > global_tol)
                    {
                        tag_ptr_global_tol[idx_nghost] = 1;
                    }
                }
            }
            
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear index.
                const int idx_nghost = i;
                
                tag_ptr_multiresolution_tagger[idx_nghost] &= tag_ptr_global_tol[idx_nghost];
            }
        }
        
        if (uses_local_tol)
        {
            // Allocate temporary patch data.
            HAMERS_SHARED_PTR<pdat::CellData<int> > tags_local_tol(
                new pdat::CellData<int>(interior_box, d_dim.getValue(), hier::IntVector::getZero(d_dim)));
            
            tags_local_tol->fillAll(0);
            
            // Get the pointer to the tags.
            int* tag_ptr_local_tol = tags_local_tol->getPointer(0);
            
            for (int li = 0; li < d_Harten_wavelet_num_level; li++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx = i + num_ghosts_0_wavelet_coeffs;
                    const int idx_nghost = i;
                    
                    if (w[li][idx]/(u_mean[li][idx] + EPSILON) > local_tol)
                    {
                        tag_ptr_local_tol[idx_nghost] = 1;
                    }
                }
            }
            
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear index.
                const int idx_nghost = i;
                
                tag_ptr_multiresolution_tagger[idx_nghost] &= tag_ptr_local_tol[idx_nghost];
            }
        }
        
        if (uses_alpha_tol)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx = i + num_ghosts_0_wavelet_coeffs;
                const int idx_nghost = i;
                
                if (alpha[idx] < alpha_tol)
                {
                    tag_ptr_multiresolution_tagger[idx_nghost] &= 1;
                }
                else
                {
                    tag_ptr_multiresolution_tagger[idx_nghost] &= 0;
                }
            }
        }
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = 0; i < interior_dim_0; i++)
        {
            // Compute the linear index.
            const int idx_nghost = i;
            
            tag_ptr[idx_nghost] |= tag_ptr_multiresolution_tagger[idx_nghost];
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        
        const int num_ghosts_0_wavelet_coeffs = num_ghosts_wavelet_coeffs[0];
        const int num_ghosts_1_wavelet_coeffs = num_ghosts_wavelet_coeffs[1];
        const int ghostcell_dim_0_wavelet_coeffs = ghostcell_dims_wavelet_coeffs[0];
        
        if (uses_global_tol)
        {
            // Allocate temporary patch data.
            HAMERS_SHARED_PTR<pdat::CellData<int> > tags_global_tol(
                new pdat::CellData<int>(interior_box, d_dim.getValue(), hier::IntVector::getZero(d_dim)));
            
            tags_global_tol->fillAll(0);
            
            // Get the pointer to the tags.
            int* tag_ptr_global_tol = tags_global_tol->getPointer(0);
            
            for (int li = 0; li < d_Harten_wavelet_num_level; li++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx = (i + num_ghosts_0_wavelet_coeffs) +
                            (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs;
                        
                        const int idx_nghost = i + j*interior_dim_0;
                        
                        if (w[li][idx]/(wavelet_coeffs_maxs[li] + EPSILON) > global_tol)
                        {
                            tag_ptr_global_tol[idx_nghost] = 1;
                        }
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
                    const int idx_nghost = i + j*interior_dim_0;
                    
                    tag_ptr_multiresolution_tagger[idx_nghost] &= tag_ptr_global_tol[idx_nghost];
                }
            }
        }
        
        if (uses_local_tol)
        {
            // Allocate temporary patch data.
            HAMERS_SHARED_PTR<pdat::CellData<int> > tags_local_tol(
                new pdat::CellData<int>(interior_box, d_dim.getValue(), hier::IntVector::getZero(d_dim)));
            
            tags_local_tol->fillAll(0);
            
            // Get the pointer to the tags.
            int* tag_ptr_local_tol = tags_local_tol->getPointer(0);
            
            for (int li = 0; li < d_Harten_wavelet_num_level; li++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx = (i + num_ghosts_0_wavelet_coeffs) +
                            (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs;
                        
                        const int idx_nghost = i + j*interior_dim_0;
                        
                        if (w[li][idx]/(u_mean[li][idx] + EPSILON) > local_tol)
                        {
                            tag_ptr_local_tol[idx_nghost] = 1;
                        }
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
                    const int idx_nghost = i + j*interior_dim_0;
                    
                    tag_ptr_multiresolution_tagger[idx_nghost] &= tag_ptr_local_tol[idx_nghost];
                }
            }
        }
        
        if (uses_alpha_tol)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx = (i + num_ghosts_0_wavelet_coeffs) +
                        (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs;
                    
                    const int idx_nghost = i + j*interior_dim_0;
                    
                    if (alpha[idx] < alpha_tol)
                    {
                        tag_ptr_multiresolution_tagger[idx_nghost] &= 1;
                    }
                    else
                    {
                        tag_ptr_multiresolution_tagger[idx_nghost] &= 0;
                    }
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
                const int idx_nghost = i + j*interior_dim_0;
                
                tag_ptr[idx_nghost] |= tag_ptr_multiresolution_tagger[idx_nghost];
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        const int interior_dim_2 = interior_dims[2];
        
        const int num_ghosts_0_wavelet_coeffs = num_ghosts_wavelet_coeffs[0];
        const int num_ghosts_1_wavelet_coeffs = num_ghosts_wavelet_coeffs[1];
        const int num_ghosts_2_wavelet_coeffs = num_ghosts_wavelet_coeffs[2];
        const int ghostcell_dim_0_wavelet_coeffs = ghostcell_dims_wavelet_coeffs[0];
        const int ghostcell_dim_1_wavelet_coeffs = ghostcell_dims_wavelet_coeffs[1];
        
        if (uses_global_tol)
        {
            // Allocate temporary patch data.
            HAMERS_SHARED_PTR<pdat::CellData<int> > tags_global_tol(
                new pdat::CellData<int>(interior_box, d_dim.getValue(), hier::IntVector::getZero(d_dim)));
            
            tags_global_tol->fillAll(0);
            
            // Get the pointer to the tags.
            int* tag_ptr_global_tol = tags_global_tol->getPointer(0);
            
            for (int li = 0; li < d_Harten_wavelet_num_level; li++)
            {
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
                            const int idx = (i + num_ghosts_0_wavelet_coeffs) +
                                (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                (k + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                    ghostcell_dim_1_wavelet_coeffs;
                            
                            const int idx_nghost = i + j*interior_dim_0 + k*interior_dim_0*interior_dim_1;
                            
                            if (w[li][idx]/(wavelet_coeffs_maxs[li] + EPSILON) > global_tol)
                            {
                                tag_ptr_global_tol[idx_nghost] = 1;
                            }
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
                        const int idx_nghost = i + j*interior_dim_0 + k*interior_dim_0*interior_dim_1;
                        
                        tag_ptr_multiresolution_tagger[idx_nghost] &= tag_ptr_global_tol[idx_nghost];
                    }
                }
            }
        }
        
        if (uses_local_tol)
        {
            // Allocate temporary patch data.
            HAMERS_SHARED_PTR<pdat::CellData<int> > tags_local_tol(
                new pdat::CellData<int>(interior_box, d_dim.getValue(), hier::IntVector::getZero(d_dim)));
            
            tags_local_tol->fillAll(0);
            
            // Get the pointer to the tags.
            int* tag_ptr_local_tol = tags_local_tol->getPointer(0);
            
            for (int li = 0; li < d_Harten_wavelet_num_level; li++)
            {
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
                            const int idx = (i + num_ghosts_0_wavelet_coeffs) +
                                (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                (k + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                    ghostcell_dim_1_wavelet_coeffs;
                            
                            const int idx_nghost = i + j*interior_dim_0 + k*interior_dim_0*interior_dim_1;
                            
                            if (w[li][idx]/(u_mean[li][idx] + EPSILON) > local_tol)
                            {
                                tag_ptr_local_tol[idx_nghost] = 1;
                            }
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
                        const int idx_nghost = i + j*interior_dim_0 + k*interior_dim_0*interior_dim_1;
                        
                        tag_ptr_multiresolution_tagger[idx_nghost] &= tag_ptr_local_tol[idx_nghost];
                    }
                }
            }
        }
        
        if (uses_alpha_tol)
        {
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
                        const int idx = (i + num_ghosts_0_wavelet_coeffs) +
                            (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                            (k + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                ghostcell_dim_1_wavelet_coeffs;
                        
                        const int idx_nghost = i + j*interior_dim_0 + k*interior_dim_0*interior_dim_1;
                        
                        if (alpha[idx] < alpha_tol)
                        {
                            tag_ptr_multiresolution_tagger[idx_nghost] &= 1;
                        }
                        else
                        {
                            tag_ptr_multiresolution_tagger[idx_nghost] &= 0;
                        }
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
                    const int idx_nghost = i + j*interior_dim_0 + k*interior_dim_0*interior_dim_1;
                    
                    tag_ptr[idx_nghost] |= tag_ptr_multiresolution_tagger[idx_nghost];
                }
            }
        }
    }
}
