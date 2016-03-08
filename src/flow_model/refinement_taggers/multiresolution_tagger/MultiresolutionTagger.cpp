#include "flow_model/refinement_taggers/multiresolution_tagger/MultiresolutionTagger.hpp"

#include "boost/lexical_cast.hpp"

#define PLOTTING_MULTIRESOLUTION_TAGGER

#define EPSILON 1e-40

MultiresolutionTagger::MultiresolutionTagger(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
    const FLOW_MODEL& flow_model,
    const int& num_species,
    const boost::shared_ptr<EquationOfState>& equation_of_state,
    const boost::shared_ptr<tbox::Database>& multiresolution_tagger_db):
        d_object_name(object_name),
        d_dim(dim),
        d_grid_geometry(grid_geometry),
        d_num_ghosts(hier::IntVector::getZero(d_dim)),
        d_flow_model(flow_model),
        d_num_species(num_species),
        d_equation_of_state(equation_of_state),
        d_variables_set(false),
        d_num_ghosts_set(false),
        d_Harten_wavelet_uses_global_tol(false),
        d_Harten_wavelet_uses_local_tol(false),
        d_Harten_wavelet_uses_alpha_tol(false)
{
    if (multiresolution_tagger_db != nullptr)
    {
        std::vector<std::string> sensor_keys = multiresolution_tagger_db->getAllKeys();
        
        const int num_keys = static_cast<int>(sensor_keys.size());
        
        if (multiresolution_tagger_db->keyExists("multiresolution_sensors"))
        {
            d_multiresolution_sensors = multiresolution_tagger_db->getStringVector("multiresolution_sensors");
        }
        else if (multiresolution_tagger_db->keyExists("d_multiresolution_sensors"))
        {
            d_multiresolution_sensors = multiresolution_tagger_db->getStringVector("d_multiresolution_sensors");
        }
        else
        {
            TBOX_WARNING(d_object_name
                << ": "
                << "No key 'multiresolution_sensors'/'d_multiresolution_sensors' found in data for"
                << " Multiresolution_tagger. No refinement with multiresolution sensors will occur."
                << std::endl);
        }
        
        /*
         * Loop over the multiresolution sensors chosen.
         */
        
        std::vector<std::string> sensor_keys_defined(num_keys);
        int sensor_keys_count = 0;
        boost::shared_ptr<tbox::Database> sensor_db;
        for (int i = 0; i < num_keys; i++)
        {
            std::string sensor_key = sensor_keys[i];
            sensor_db.reset();
            
            if (!((sensor_key == "multiresolution_sensors") || (sensor_key == "d_multiresolution_sensors")))
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
                        d_Harten_wavelet_num_level = sensor_db->getInteger("Harten_wavelet_num_level");
                    }
                    else if (sensor_db->keyExists("d_Harten_wavelet_num_level"))
                    {
                        d_Harten_wavelet_num_level = sensor_db->getInteger("d_Harten_wavelet_num_level");
                    }
                    else
                    {
                        TBOX_ERROR(d_object_name
                            << ": "
                            << "No key 'Harten_wavelet_num_level'/'d_Harten_wavelet_num_level' found in data for "
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
                        d_Harten_wavelet_num_vanishing_moments = sensor_db->getInteger("Harten_wavelet_num_vanishing_moments");
                    }
                    else if (sensor_db->keyExists("d_Harten_wavelet_num_vanishing_moments"))
                    {
                        d_Harten_wavelet_num_vanishing_moments = sensor_db->getInteger("d_Harten_wavelet_num_vanishing_moments");
                    }
                    else
                    {
                        TBOX_ERROR(d_object_name
                            << ": "
                            << "No key 'Harten_wavelet_num_vanishing_moments'/'d_Harten_wavelet_num_vanishing_moments' found in data for "
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
                        d_Harten_wavelet_variables = sensor_db->getStringVector("Harten_wavelet_variables");
                    }
                    else if (sensor_db->keyExists("d_Harten_wavelet_variables"))
                    {
                        d_Harten_wavelet_variables = sensor_db->getStringVector("d_Harten_wavelet_variables");
                    }
                    else
                    {
                        TBOX_ERROR(d_object_name
                            << ": "
                            << "No key 'Harten_wavelet_variables'/'d_Harten_wavelet_variables' found in data for "
                            << sensor_key
                            << "."
                            << std::endl);
                    }
                    
                    for (int vi = 0; vi < static_cast<int>(d_Harten_wavelet_variables.size()); vi++)
                    {
                        std::string variable_key = d_Harten_wavelet_variables[vi];
                        
                        switch (d_flow_model)
                        {
                            case SINGLE_SPECIES:
                            {
                                if (!((variable_key == "DENSITY") ||
                                      (variable_key == "TOTAL_ENERGY") ||
                                      (variable_key == "PRESSURE") ||
                                      (variable_key == "ENSTROPHY")))
                                {
                                    TBOX_ERROR(d_object_name
                                        << ": "
                                        << "Unknown/unsupported variable: "
                                        << variable_key
                                        << "\nin input."
                                        << std::endl);
                                }
                                
                                break;
                            }
                            case FOUR_EQN_SHYUE:
                            {
                                if (!((variable_key == "DENSITY") ||
                                      (variable_key == "TOTAL_ENERGY") ||
                                      (variable_key == "PRESSURE") ||
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
                                
                                break;
                            }
                            case FIVE_EQN_ALLAIRE:
                            {
                                if (!((variable_key == "DENSITY") ||
                                      (variable_key == "TOTAL_ENERGY") ||
                                      (variable_key == "PRESSURE") ||
                                      (variable_key == "ENSTROPHY") ||
                                      (variable_key == "MASS_FRACTION") ||
                                      (variable_key == "VOLUME_FRACTION")))
                                {
                                    TBOX_ERROR(d_object_name
                                        << ": "
                                        << "Unknown/unsupported variable: "
                                        << variable_key
                                        << "\nin input."
                                        << std::endl);
                                }
                                
                                break;
                            }
                            default:
                            {
                                TBOX_ERROR(d_object_name
                                    << ": "
                                    << "d_flow_model '"
                                    << d_flow_model
                                    << "' not yet implemented."
                                    << std::endl);
                            }
                        }
                    }
                    
                    // Get the settings for wavelet global filter.
                    if (sensor_db->keyExists("Harten_wavelet_uses_global_tol"))
                    {
                        d_Harten_wavelet_uses_global_tol = sensor_db->getBool("Harten_wavelet_uses_global_tol");
                    }
                    else if (sensor_db->keyExists("d_Harten_wavelet_uses_global_tol"))
                    {
                        d_Harten_wavelet_uses_global_tol = sensor_db->getBool("d_Harten_wavelet_uses_global_tol");
                    }
                    else
                    {
                        TBOX_ERROR(d_object_name
                            << ": "
                            << "No key 'Harten_wavelet_uses_global_tol'/'d_Harten_wavelet_uses_global_tol'"
                            << " found in data for "
                            << sensor_key
                            << "."
                            << std::endl);
                    }
                    
                    if (d_Harten_wavelet_uses_global_tol)
                    {
                        if (sensor_db->keyExists("Harten_wavelet_global_tol"))
                        {
                            d_Harten_wavelet_global_tol = sensor_db->getDoubleVector("Harten_wavelet_global_tol");
                        }
                        else if (sensor_db->keyExists("d_Harten_wavelet_global_tol"))
                        {
                            d_Harten_wavelet_global_tol = sensor_db->getDoubleVector("d_Harten_wavelet_global_tol");
                        }
                        else
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "No key 'Harten_wavelet_global_tol'/'d_Harten_wavelet_global_tol' found in data for "
                                << sensor_key
                                << "."
                                << std::endl);
                        }
                        
                        if (static_cast<int>(d_Harten_wavelet_variables.size()) != static_cast<int>(d_Harten_wavelet_global_tol.size()))
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "The numbers of variables and global tolerances provided don't match for "
                                << sensor_key
                                << "."
                                << std::endl);
                        }
                    }
                    
                    // Get the settings for wavelet local filter.
                    if (sensor_db->keyExists("Harten_wavelet_uses_local_tol"))
                    {
                        d_Harten_wavelet_uses_local_tol = sensor_db->getBool("Harten_wavelet_uses_local_tol");
                    }
                    else if (sensor_db->keyExists("d_Harten_wavelet_uses_local_tol"))
                    {
                        d_Harten_wavelet_uses_local_tol = sensor_db->getBool("d_Harten_wavelet_uses_local_tol");
                    }
                    else
                    {
                        TBOX_ERROR(d_object_name
                            << ": "
                            << "No key 'Harten_wavelet_uses_local_tol'/'d_Harten_wavelet_uses_local_tol'"
                            << " found in data for "
                            << sensor_key
                            << "."
                            << std::endl);
                    }
                    
                    if (d_Harten_wavelet_uses_local_tol)
                    {
                        if (sensor_db->keyExists("Harten_wavelet_local_tol"))
                        {
                            d_Harten_wavelet_local_tol = sensor_db->getDoubleVector("Harten_wavelet_local_tol");
                        }
                        else if (sensor_db->keyExists("d_Harten_wavelet_local_tol"))
                        {
                            d_Harten_wavelet_local_tol = sensor_db->getDoubleVector("d_Harten_wavelet_local_tol");
                        }
                        else
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "No key 'Harten_wavelet_local_tol'/'d_Harten_wavelet_local_tol' found in data for "
                                << sensor_key
                                << "."
                                << std::endl);
                        }
                        
                        if (static_cast<int>(d_Harten_wavelet_variables.size()) != static_cast<int>(d_Harten_wavelet_local_tol.size()))
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "The numbers of variables and local tolerances provided don't match for "
                                << sensor_key
                                << "."
                                << std::endl);
                        }
                    }
                    
                    // Get the settings for wavelet Lipschitz's filter.
                    if (sensor_db->keyExists("Harten_wavelet_uses_alpha_tol"))
                    {
                        d_Harten_wavelet_uses_alpha_tol = sensor_db->getBool("Harten_wavelet_uses_alpha_tol");
                    }
                    else if (sensor_db->keyExists("d_Harten_wavelet_uses_alpha_tol"))
                    {
                        d_Harten_wavelet_uses_alpha_tol = sensor_db->getBool("d_Harten_wavelet_uses_alpha_tol");
                    }
                    else
                    {
                        TBOX_ERROR(d_object_name
                            << ": "
                            << "No key 'Harten_wavelet_uses_alpha_tol'/'d_Harten_wavelet_uses_alpha_tol'"
                            << " found in data for "
                            << sensor_key
                            << "."
                            << std::endl);
                    }
                    
                    if (d_Harten_wavelet_uses_alpha_tol)
                    {
                        if (sensor_db->keyExists("Harten_wavelet_alpha_tol"))
                        {
                            d_Harten_wavelet_alpha_tol = sensor_db->getDoubleVector("Harten_wavelet_alpha_tol");
                        }
                        else if (sensor_db->keyExists("d_Harten_wavelet_alpha_tol"))
                        {
                            d_Harten_wavelet_alpha_tol = sensor_db->getDoubleVector("d_Harten_wavelet_alpha_tol");
                        }
                        else
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "No key 'Harten_wavelet_alpha_tol'/'d_Harten_wavelet_alpha_tol' found in data for "
                                << sensor_key
                                << "."
                                << std::endl);
                        }
                        
                        if (static_cast<int>(d_Harten_wavelet_variables.size()) != static_cast<int>(d_Harten_wavelet_alpha_tol.size()))
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "The numbers of variables and Lipschitz's tolerances provided don't match for "
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
}


/*
 * Register the temporary variables used in multiresolution tagger class.
 */
void
MultiresolutionTagger::registerMultiresolutionTaggerVariables(
    RungeKuttaLevelIntegrator* integrator)
{
    if (d_variables_set == true)
    {
        if (d_num_ghosts_set == true)
        {
            for (int si = 0;
                     si < static_cast<int>(d_multiresolution_sensors.size());
                     si++)
            {
                std::string sensor_key = d_multiresolution_sensors[si];
                
                if (sensor_key == "HARTEN_WAVELET")
                {
                    d_Harten_Lipschitz_exponent =
                        boost::shared_ptr<pdat::CellVariable<double> > (
                            new pdat::CellVariable<double>(d_dim, "Harten Lipschitz's exponent", 1));
                    
                    for (int vi = 0; vi < static_cast<int>(d_Harten_wavelet_variables.size()); vi++)
                    {
                        // Get the key of the current variable.
                        std::string variable_key = d_Harten_wavelet_variables[vi];
                        
                        if (variable_key == "DENSITY")
                        {
                            for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                            {
                                d_Harten_density_wavelet_coeffs.push_back(
                                    boost::make_shared<pdat::CellVariable<double> >(
                                        d_dim,
                                        "Harten density wavelet coefficients at level " + boost::lexical_cast<std::string>(li),
                                        1));
                                
                                if (d_Harten_wavelet_uses_global_tol)
                                {
                                    d_Harten_density_wavelet_coeffs_maxs.push_back(0.0);
                                }
                                
                                if (d_Harten_wavelet_uses_local_tol)
                                {
                                    d_Harten_density_local_means.push_back(
                                        boost::make_shared<pdat::CellVariable<double> >(
                                            d_dim,
                                            "Harten density local means at level " + boost::lexical_cast<std::string>(li),
                                            1));
                                }
                            }
                        }
                        else if (variable_key == "TOTAL_ENERGY")
                        {
                            for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                            {
                                d_Harten_total_energy_wavelet_coeffs.push_back(
                                    boost::make_shared<pdat::CellVariable<double> >(
                                        d_dim,
                                        "Harten total energy wavelet coefficients at level " + boost::lexical_cast<std::string>(li),
                                        1));
                                
                                if (d_Harten_wavelet_uses_global_tol)
                                {
                                    d_Harten_total_energy_wavelet_coeffs_maxs.push_back(0.0);
                                }
                                
                                if (d_Harten_wavelet_uses_local_tol)
                                {
                                    d_Harten_total_energy_local_means.push_back(
                                        boost::make_shared<pdat::CellVariable<double> >(
                                            d_dim,
                                            "Harten total energy local means at level " + boost::lexical_cast<std::string>(li),
                                            1));
                                }
                            }
                        }
                        else if (variable_key == "PRESSURE")
                        {
                            for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                            {
                                d_Harten_pressure_wavelet_coeffs.push_back(
                                    boost::make_shared<pdat::CellVariable<double> >(
                                        d_dim,
                                        "Harten pressure wavelet coefficients at level " + boost::lexical_cast<std::string>(li),
                                        1));
                                
                                if (d_Harten_wavelet_uses_global_tol)
                                {
                                    d_Harten_pressure_wavelet_coeffs_maxs.push_back(0.0);
                                }
                                
                                if (d_Harten_wavelet_uses_local_tol)
                                {
                                    d_Harten_pressure_local_means.push_back(
                                        boost::make_shared<pdat::CellVariable<double> >(
                                            d_dim,
                                            "Harten pressure local means at level " + boost::lexical_cast<std::string>(li),
                                            1));
                                }
                            }
                        }
                        else if (variable_key == "ENSTROPHY")
                        {
                            for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                            {
                                d_Harten_enstrophy_wavelet_coeffs.push_back(
                                    boost::make_shared<pdat::CellVariable<double> >(
                                        d_dim,
                                        "Harten enstrophy wavelet coefficients at level " + boost::lexical_cast<std::string>(li),
                                        1));
                                
                                if (d_Harten_wavelet_uses_global_tol)
                                {
                                    d_Harten_enstrophy_wavelet_coeffs_maxs.push_back(0.0);
                                }
                                
                                if (d_Harten_wavelet_uses_local_tol)
                                {
                                    d_Harten_enstrophy_local_means.push_back(
                                        boost::make_shared<pdat::CellVariable<double> >(
                                            d_dim,
                                            "Harten enstrophy local means at level " + boost::lexical_cast<std::string>(li),
                                            1));
                                }
                            }
                        }
                        else if (variable_key == "MASS_FRACTION")
                        {
                            switch (d_flow_model)
                            {
                                case SINGLE_SPECIES:
                                {
                                    TBOX_ERROR(d_object_name
                                    << ": '"
                                    << variable_key
                                    << "' not supported for '"
                                    << d_flow_model
                                    << "' flow model."
                                    << std::endl);
                                
                                    break;
                                }
                                case FOUR_EQN_SHYUE:
                                {
                                    for (int di = 0; di < d_mass_fraction->getDepth(); di++)
                                    {
                                        for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                                        {
                                            d_Harten_mass_fraction_wavelet_coeffs.push_back(
                                                boost::make_shared<pdat::CellVariable<double> >(
                                                    d_dim,
                                                    "Harten mass fraction " + boost::lexical_cast<std::string>(di) +
                                                        " wavelet coefficients at level " +
                                                        boost::lexical_cast<std::string>(li),
                                                    1));
                                            
                                            if (d_Harten_wavelet_uses_global_tol)
                                            {
                                                d_Harten_mass_fraction_wavelet_coeffs_maxs.push_back(0.0);
                                            }
                                            
                                            if (d_Harten_wavelet_uses_local_tol)
                                            {
                                                d_Harten_mass_fraction_local_means.push_back(
                                                    boost::make_shared<pdat::CellVariable<double> >(
                                                        d_dim,
                                                        "Harten mass fraction " + boost::lexical_cast<std::string>(di) +
                                                            " local means at level " +
                                                            boost::lexical_cast<std::string>(li),
                                                        1));
                                            }
                                        }
                                    }
                                    
                                    break;
                                }
                                case FIVE_EQN_ALLAIRE:
                                {
                                    for (int di = 0; di < d_num_species; di++)
                                    {
                                        for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                                        {
                                            d_Harten_mass_fraction_wavelet_coeffs.push_back(
                                                boost::make_shared<pdat::CellVariable<double> >(
                                                    d_dim,
                                                    "Harten mass fraction " + boost::lexical_cast<std::string>(di) +
                                                        " wavelet coefficients at level " +
                                                        boost::lexical_cast<std::string>(li),
                                                    1));
                                            
                                            if (d_Harten_wavelet_uses_global_tol)
                                            {
                                                d_Harten_mass_fraction_wavelet_coeffs_maxs.push_back(0.0);
                                            }
                                            
                                            if (d_Harten_wavelet_uses_local_tol)
                                            {
                                                d_Harten_mass_fraction_local_means.push_back(
                                                    boost::make_shared<pdat::CellVariable<double> >(
                                                        d_dim,
                                                        "Harten mass fraction " + boost::lexical_cast<std::string>(di) +
                                                            " local means at level " +
                                                            boost::lexical_cast<std::string>(li),
                                                        1));
                                            }
                                        }
                                    }
                                    
                                    break;
                                }
                                default:
                                {
                                    TBOX_ERROR(d_object_name
                                        << ": "
                                        << "d_flow_model '"
                                        << d_flow_model
                                        << "' not yet implemented."
                                        << std::endl);
                                }
                            }
                        }
                        else if (variable_key == "VOLUME_FRACTION")
                        {
                            switch (d_flow_model)
                            {
                                case SINGLE_SPECIES: case FOUR_EQN_SHYUE:
                                {
                                    TBOX_ERROR(d_object_name
                                        << ": '"
                                        << variable_key
                                        << "' not supported for '"
                                        << d_flow_model
                                        << "' flow model."
                                        << std::endl);
                                    
                                    break;
                                }
                                case FIVE_EQN_ALLAIRE:
                                {
                                    for (int di = 0; di < d_volume_fraction->getDepth(); di++)
                                    {
                                        for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                                        {
                                            d_Harten_volume_fraction_wavelet_coeffs.push_back(
                                                boost::make_shared<pdat::CellVariable<double> >(
                                                    d_dim,
                                                    "Harten volume fraction " + boost::lexical_cast<std::string>(di) +
                                                        " wavelet coefficients at level " +
                                                        boost::lexical_cast<std::string>(li),
                                                    1));
                                            
                                            if (d_Harten_wavelet_uses_global_tol)
                                            {
                                                d_Harten_volume_fraction_wavelet_coeffs_maxs.push_back(0.0);
                                            }
                                            
                                            if (d_Harten_wavelet_uses_local_tol)
                                            {
                                                d_Harten_volume_fraction_local_means.push_back(
                                                    boost::make_shared<pdat::CellVariable<double> >(
                                                        d_dim,
                                                        "Harten volume fraction " + boost::lexical_cast<std::string>(di) +
                                                            " local means at level " +
                                                            boost::lexical_cast<std::string>(li),
                                                        1));
                                            }
                                        }
                                    }
                                    
                                    break;
                                }
                                default:
                                {
                                    TBOX_ERROR(d_object_name
                                        << ": "
                                        << "d_flow_model '"
                                        << d_flow_model
                                        << "' not yet implemented."
                                        << std::endl);
                                }
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
                    
                    integrator->registerVariable(
                        d_Harten_Lipschitz_exponent,
                        d_num_ghosts,
                        RungeKuttaLevelIntegrator::TIME_DEP,
                            d_grid_geometry,
                            "CONSERVATIVE_COARSEN",
                            "CONSERVATIVE_LINEAR_REFINE");
                    
                    for (int vi = 0; vi < static_cast<int>(d_Harten_wavelet_variables.size()); vi++)
                    {
                        // Get the key of the current variable.
                        std::string variable_key = d_Harten_wavelet_variables[vi];
                        
                        if (variable_key == "DENSITY")
                        {
                            for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                            {
                                integrator->registerVariable(
                                    d_Harten_density_wavelet_coeffs[li],
                                    d_num_ghosts,
                                    RungeKuttaLevelIntegrator::TIME_DEP,
                                        d_grid_geometry,
                                        "CONSERVATIVE_COARSEN",
                                        "CONSERVATIVE_LINEAR_REFINE");
                                
                                if (d_Harten_wavelet_uses_local_tol)
                                {
                                    integrator->registerVariable(
                                        d_Harten_density_local_means[li],
                                        d_num_ghosts,
                                        RungeKuttaLevelIntegrator::TEMPORARY,
                                            d_grid_geometry,
                                            "NO_COARSEN",
                                            "NO_REFINE");
                                }
                            }
                        }
                        else if (variable_key == "TOTAL_ENERGY")
                        {
                            for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                            {
                                integrator->registerVariable(
                                    d_Harten_total_energy_wavelet_coeffs[li],
                                    d_num_ghosts,
                                    RungeKuttaLevelIntegrator::TIME_DEP,
                                    d_grid_geometry,
                                    "CONSERVATIVE_COARSEN",
                                    "CONSERVATIVE_LINEAR_REFINE");
                                
                                if (d_Harten_wavelet_uses_local_tol)
                                {
                                    integrator->registerVariable(
                                        d_Harten_total_energy_local_means[li],
                                        d_num_ghosts,
                                        RungeKuttaLevelIntegrator::TEMPORARY,
                                            d_grid_geometry,
                                            "NO_COARSEN",
                                            "NO_REFINE");
                                }
                            }
                        }
                        else if (variable_key == "PRESSURE")
                        {
                            for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                            {
                                integrator->registerVariable(
                                    d_Harten_pressure_wavelet_coeffs[li],
                                    d_num_ghosts,
                                    RungeKuttaLevelIntegrator::TIME_DEP,
                                    d_grid_geometry,
                                    "CONSERVATIVE_COARSEN",
                                    "CONSERVATIVE_LINEAR_REFINE");
                                
                                if (d_Harten_wavelet_uses_local_tol)
                                {
                                    integrator->registerVariable(
                                        d_Harten_pressure_local_means[li],
                                        d_num_ghosts,
                                        RungeKuttaLevelIntegrator::TEMPORARY,
                                            d_grid_geometry,
                                            "NO_COARSEN",
                                            "NO_REFINE");
                                }
                            }
                        }
                        else if (variable_key == "ENSTROPHY")
                        {
                            for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                            {
                                integrator->registerVariable(
                                    d_Harten_enstrophy_wavelet_coeffs[li],
                                    d_num_ghosts,
                                    RungeKuttaLevelIntegrator::TIME_DEP,
                                    d_grid_geometry,
                                    "CONSERVATIVE_COARSEN",
                                    "CONSERVATIVE_LINEAR_REFINE");
                                
                                if (d_Harten_wavelet_uses_local_tol)
                                {
                                    integrator->registerVariable(
                                        d_Harten_enstrophy_local_means[li],
                                        d_num_ghosts,
                                        RungeKuttaLevelIntegrator::TEMPORARY,
                                            d_grid_geometry,
                                            "NO_COARSEN",
                                            "NO_REFINE");
                                }
                            }
                        }
                        else if (variable_key == "MASS_FRACTION")
                        {
                            switch (d_flow_model)
                            {
                                case SINGLE_SPECIES:
                                {
                                    TBOX_ERROR(d_object_name
                                        << ": '"
                                        << variable_key
                                        << "' not supported for '"
                                        << d_flow_model
                                        << "' flow model."
                                        << std::endl);
                                    
                                    break;
                                }
                                case FOUR_EQN_SHYUE:
                                {
                                    for (int di = 0; di < d_mass_fraction->getDepth(); di++)
                                    {
                                        for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                                        {
                                            integrator->registerVariable(
                                                d_Harten_mass_fraction_wavelet_coeffs[di*d_Harten_wavelet_num_level + li],
                                                d_num_ghosts,
                                                RungeKuttaLevelIntegrator::TIME_DEP,
                                                d_grid_geometry,
                                                "CONSERVATIVE_COARSEN",
                                                "CONSERVATIVE_LINEAR_REFINE");
                                            
                                            if (d_Harten_wavelet_uses_local_tol)
                                            {
                                                integrator->registerVariable(
                                                    d_Harten_mass_fraction_local_means[di*d_Harten_wavelet_num_level + li],
                                                    d_num_ghosts,
                                                    RungeKuttaLevelIntegrator::TEMPORARY,
                                                        d_grid_geometry,
                                                        "NO_COARSEN",
                                                        "NO_REFINE");
                                            }
                                        }
                                    }
                                    
                                    break;
                                }
                                case FIVE_EQN_ALLAIRE:
                                {
                                    for (int di = 0; di < d_num_species; di++)
                                    {
                                        for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                                        {
                                            integrator->registerVariable(
                                                d_Harten_mass_fraction_wavelet_coeffs[di*d_Harten_wavelet_num_level + li],
                                                d_num_ghosts,
                                                RungeKuttaLevelIntegrator::TIME_DEP,
                                                d_grid_geometry,
                                                "CONSERVATIVE_COARSEN",
                                                "CONSERVATIVE_LINEAR_REFINE");
                                            
                                            if (d_Harten_wavelet_uses_local_tol)
                                            {
                                                integrator->registerVariable(
                                                    d_Harten_mass_fraction_local_means[di*d_Harten_wavelet_num_level + li],
                                                    d_num_ghosts,
                                                    RungeKuttaLevelIntegrator::TEMPORARY,
                                                        d_grid_geometry,
                                                        "NO_COARSEN",
                                                        "NO_REFINE");
                                            }
                                        }
                                    }
                                    
                                    break;
                                }
                                default:
                                {
                                    TBOX_ERROR(d_object_name
                                        << ": "
                                        << "d_flow_model '"
                                        << d_flow_model
                                        << "' not yet implemented."
                                        << std::endl);
                                }
                            }
                        }
                        else if (variable_key == "VOLUME_FRACTION")
                        {
                            switch (d_flow_model)
                            {
                                case SINGLE_SPECIES: case FOUR_EQN_SHYUE:
                                {
                                    TBOX_ERROR(d_object_name
                                        << ": '"
                                        << variable_key
                                        << "' not supported for '"
                                        << d_flow_model
                                        << "' flow model."
                                        << std::endl);
                                    
                                    break;
                                }
                                case FIVE_EQN_ALLAIRE:
                                {
                                    for (int di = 0; di < d_volume_fraction->getDepth(); di++)
                                    {
                                        for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                                        {
                                            integrator->registerVariable(
                                                d_Harten_volume_fraction_wavelet_coeffs[di*d_Harten_wavelet_num_level + li],
                                                d_num_ghosts,
                                                RungeKuttaLevelIntegrator::TIME_DEP,
                                                d_grid_geometry,
                                                "CONSERVATIVE_COARSEN",
                                                "CONSERVATIVE_LINEAR_REFINE");
                                            
                                            if (d_Harten_wavelet_uses_local_tol)
                                            {
                                                integrator->registerVariable(
                                                    d_Harten_volume_fraction_local_means[di*d_Harten_wavelet_num_level + li],
                                                    d_num_ghosts,
                                                    RungeKuttaLevelIntegrator::TEMPORARY,
                                                        d_grid_geometry,
                                                        "NO_COARSEN",
                                                        "NO_REFINE");
                                            }
                                        }
                                    }
                                    
                                    break;
                                }
                                default:
                                {
                                    TBOX_ERROR(d_object_name
                                        << ": "
                                        << "d_flow_model '"
                                        << d_flow_model
                                        << "' not yet implemented."
                                        << std::endl);
                                }
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
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Number of ghost cells is not set yet."
                << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Variables are not set yet."
            << std::endl);
    }
}


/*
 * Register the plotting quantities.
 */
void
MultiresolutionTagger::registerPlotQuantities(
    RungeKuttaLevelIntegrator* integrator,
    const boost::shared_ptr<appu::VisItDataWriter>& visit_writer,
    const boost::shared_ptr<hier::VariableContext>& plot_context)
{
#ifdef PLOTTING_MULTIRESOLUTION_TAGGER
    if (d_variables_set == true)
    {
        if (d_num_ghosts_set == true)
        {
            hier::VariableDatabase* vardb = hier::VariableDatabase::getDatabase();
            
            for (int si = 0;
                     si < static_cast<int>(d_multiresolution_sensors.size());
                     si++)
            {
                std::string sensor_key = d_multiresolution_sensors[si];
                
                if (sensor_key == "HARTEN_WAVELET")
                {
                    visit_writer->registerPlotQuantity(
                        "Harten Lipschitz's exponent",
                        "SCALAR",
                        vardb->mapVariableAndContextToIndex(
                           d_Harten_Lipschitz_exponent,
                           plot_context));
                    
                    for (int vi = 0; vi < static_cast<int>(d_Harten_wavelet_variables.size()); vi++)
                    {
                        // Get the key of the current variable.
                        std::string variable_key = d_Harten_wavelet_variables[vi];
                        
                        if (variable_key == "DENSITY")
                        {
                            for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                            {
                                visit_writer->registerPlotQuantity(
                                    "Harten density wavelet coefficients at level " + boost::lexical_cast<std::string>(li),
                                    "SCALAR",
                                    vardb->mapVariableAndContextToIndex(
                                       d_Harten_density_wavelet_coeffs[li],
                                       plot_context));
                            }
                        }
                        else if (variable_key == "TOTAL_ENERGY")
                        {
                            for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                            {
                                visit_writer->registerPlotQuantity(
                                    "Harten total energy wavelet coefficients at level " + boost::lexical_cast<std::string>(li),
                                    "SCALAR",
                                    vardb->mapVariableAndContextToIndex(
                                       d_Harten_total_energy_wavelet_coeffs[li],
                                       plot_context));
                            }
                        }
                        else if (variable_key == "PRESSURE")
                        {
                            for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                            {
                                visit_writer->registerPlotQuantity(
                                    "Harten pressure wavelet coefficients at level " + boost::lexical_cast<std::string>(li),
                                    "SCALAR",
                                    vardb->mapVariableAndContextToIndex(
                                       d_Harten_pressure_wavelet_coeffs[li],
                                       plot_context));
                            }
                        }
                        else if (variable_key == "ENSTROPHY")
                        {
                            for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                            {
                                visit_writer->registerPlotQuantity(
                                    "Harten enstrophy wavelet coefficients at level " + boost::lexical_cast<std::string>(li),
                                    "SCALAR",
                                    vardb->mapVariableAndContextToIndex(
                                       d_Harten_enstrophy_wavelet_coeffs[li],
                                       plot_context));
                            }
                        }
                        else if (variable_key == "MASS_FRACTION")
                        {
                            switch (d_flow_model)
                            {
                                case SINGLE_SPECIES:
                                {
                                    TBOX_ERROR(d_object_name
                                        << ": '"
                                        << variable_key
                                        << "' not supported for '"
                                        << d_flow_model
                                        << "' flow model."
                                        << std::endl);
                                    
                                    break;
                                }
                                case FOUR_EQN_SHYUE:
                                {
                                    for (int di = 0; di < d_mass_fraction->getDepth(); di++)
                                    {
                                        for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                                        {
                                            visit_writer->registerPlotQuantity(
                                                "Harten mass fraction " + boost::lexical_cast<std::string>(di) +
                                                    " wavelet coefficients at level " + boost::lexical_cast<std::string>(li),
                                                "SCALAR",
                                                vardb->mapVariableAndContextToIndex(
                                                   d_Harten_mass_fraction_wavelet_coeffs[di*d_Harten_wavelet_num_level + li],
                                                   plot_context));
                                        }
                                    }
                                    
                                    break;
                                }
                                case FIVE_EQN_ALLAIRE:
                                {
                                    for (int di = 0; di < d_num_species; di++)
                                    {
                                        for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                                        {
                                            visit_writer->registerPlotQuantity(
                                                "Harten mass fraction " + boost::lexical_cast<std::string>(di) +
                                                    " wavelet coefficients at level " + boost::lexical_cast<std::string>(li),
                                                "SCALAR",
                                                vardb->mapVariableAndContextToIndex(
                                                   d_Harten_mass_fraction_wavelet_coeffs[di*d_Harten_wavelet_num_level + li],
                                                   plot_context));
                                        }
                                    }
                                    
                                    break;
                                }
                                default:
                                {
                                    TBOX_ERROR(d_object_name
                                        << ": "
                                        << "d_flow_model '"
                                        << d_flow_model
                                        << "' not yet implemented."
                                        << std::endl);
                                }
                            }
                        }
                        else if (variable_key == "VOLUME_FRACTION")
                        {
                            switch (d_flow_model)
                            {
                                case SINGLE_SPECIES: case FOUR_EQN_SHYUE:
                                {
                                    TBOX_ERROR(d_object_name
                                        << ": '"
                                        << variable_key
                                        << "' not supported for '"
                                        << d_flow_model
                                        << "' flow model."
                                        << std::endl);
                                    
                                    break;
                                }
                                case FIVE_EQN_ALLAIRE:
                                {
                                    for (int di = 0; di < d_volume_fraction->getDepth(); di++)
                                    {
                                        for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                                        {
                                            visit_writer->registerPlotQuantity(
                                                    "Harten volume fraction " + boost::lexical_cast<std::string>(di) +
                                                        " wavelet coefficients at level " + boost::lexical_cast<std::string>(li),
                                                    "SCALAR",
                                                    vardb->mapVariableAndContextToIndex(
                                                       d_Harten_volume_fraction_wavelet_coeffs[di*d_Harten_wavelet_num_level + li],
                                                       plot_context));
                                        }
                                    }
                                }
                                default:
                                {
                                    TBOX_ERROR(d_object_name
                                        << ": "
                                        << "d_flow_model '"
                                        << d_flow_model
                                        << "' not yet implemented."
                                        << std::endl);
                                }
                            }
                        }
                    }
                }
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Number of ghost cells is not set yet."
                << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Variables are not set yet."
            << std::endl);
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
    
    os << "d_variables_set = "
       << d_variables_set
       << std::endl;
    os << "d_num_ghosts_set = "
       << d_num_ghosts_set
       << std::endl;
    
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
    const boost::shared_ptr<tbox::Database>& restart_db) const
{
    if (static_cast<int>(d_multiresolution_sensors.size()) > 0)
    {
        restart_db->putStringVector("d_multiresolution_sensors", d_multiresolution_sensors);
    }
    
    for (int si = 0; si < static_cast<int>(d_multiresolution_sensors.size()); si++)
    {
        if (d_multiresolution_sensors[si] == "HARTEN_WAVELET")
        {
            boost::shared_ptr<tbox::Database> sensor_db =
                restart_db->putDatabase("HARTEN_WAVELET");
            
            sensor_db->putInteger("d_Harten_wavelet_num_level",
                d_Harten_wavelet_num_level);
            
            sensor_db->putInteger("d_Harten_wavelet_num_vanishing_moments",
                d_Harten_wavelet_num_vanishing_moments);
            
            sensor_db->putStringVector("d_Harten_wavelet_variables",
                d_Harten_wavelet_variables);
            
            sensor_db->putBool("d_Harten_wavelet_uses_global_tol",
                d_Harten_wavelet_uses_global_tol);
            
            sensor_db->putBool("d_Harten_wavelet_uses_local_tol",
                d_Harten_wavelet_uses_local_tol);
            
            sensor_db->putBool("d_Harten_wavelet_uses_alpha_tol",
                d_Harten_wavelet_uses_alpha_tol);
            
            if (d_Harten_wavelet_uses_global_tol)
            {
                sensor_db->putDoubleVector("d_Harten_wavelet_global_tol",
                    d_Harten_wavelet_global_tol);
            }
            
            if (d_Harten_wavelet_uses_local_tol)
            {
                sensor_db->putDoubleVector("d_Harten_wavelet_local_tol",
                    d_Harten_wavelet_local_tol);
            }
            
            if (d_Harten_wavelet_uses_alpha_tol)
            {
                sensor_db->putDoubleVector("d_Harten_wavelet_alpha_tol",
                    d_Harten_wavelet_alpha_tol);
            }
        }
    }
}


/*
 * Compute values of multiresolution sensors.
 */
void
MultiresolutionTagger::computeMultiresolutionSensorValues(
    hier::Patch& patch,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
    if (d_variables_set == true)
    {
        if (d_num_ghosts_set == true)
        {
            const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                    patch.getPatchGeometry()));
            
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(patch_geom);
#endif
            
            const double* dx = patch_geom->getDx();
            
            // Get the dimensions of box that covers the interior of Patch.
            hier::Box dummy_box = patch.getBox();
            const hier::Box interior_box = dummy_box;
            const hier::IntVector interior_dims = interior_box.numberCells();
            
            // Get the dimensions of box that covers interior of Patch plus
            // ghost cells.
            dummy_box.grow(d_num_ghosts);
            const hier::Box ghost_box = dummy_box;
            const hier::IntVector ghostcell_dims = ghost_box.numberCells();
            
            for (int si = 0;
                     si < static_cast<int>(d_multiresolution_sensors.size());
                     si++)
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
                            switch (d_flow_model)
                            {
                                case SINGLE_SPECIES: case FOUR_EQN_SHYUE:
                                {
                                    // Get the cell data of density.
                                    boost::shared_ptr<pdat::CellData<double> > density(
                                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                            patch.getPatchData(d_density, data_context)));
                                    
                                    // Get the wavelet coefficients.
                                    std::vector<boost::shared_ptr<pdat::CellData<double> > > wavelet_coeffs;
                                    for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                                    {
                                        wavelet_coeffs.push_back(BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                            patch.getPatchData(d_Harten_density_wavelet_coeffs[li], data_context)));
                                    }
                                    
                                    // Compute the wavelet coefficients.
                                    if (d_Harten_wavelet_uses_local_tol)
                                    {
                                        // Get the local means.
                                        std::vector<boost::shared_ptr<pdat::CellData<double> > > variable_local_means;
                                        for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                                        {
                                            variable_local_means.push_back(BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                                patch.getPatchData(d_Harten_density_local_means[li], data_context)));
                                        }
                                        
                                        d_wavelet_transfrom_Harten->computeWaveletCoefficientsWithVariableLocalMeans(
                                            patch,
                                            density,
                                            wavelet_coeffs,
                                            variable_local_means);
                                    }
                                    else
                                    {
                                        d_wavelet_transfrom_Harten->computeWaveletCoefficients(
                                            patch,
                                            density,
                                            wavelet_coeffs);
                                    }
                                    
                                    break;
                                }
                                case FIVE_EQN_ALLAIRE:
                                {
                                    // Get the cell data of partial density.
                                    boost::shared_ptr<pdat::CellData<double> > partial_density(
                                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                            patch.getPatchData(d_partial_density, data_context)));
                                    
                                    // Allocate temporary cell data for density.
                                    boost::shared_ptr<pdat::CellData<double> > density(
                                        new pdat::CellData<double>(interior_box, 1, d_num_ghosts));
                                    
                                    // Get the pointers to partial density and density.
                                    std::vector<double*> Z_rho;
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        Z_rho.push_back(partial_density->getPointer(si));
                                    } 
                                    double* rho   = density->getPointer(0);
                                    
                                    // Compute the field of density.
                                    if (d_dim == tbox::Dimension(1))
                                    {
                                        // NOT YET IMPLEMENTED
                                    }
                                    else if (d_dim == tbox::Dimension(2))
                                    {
                                        for (int j = -d_num_ghosts[1]; j < interior_dims[1] + d_num_ghosts[1]; j++)
                                        {
                                            for (int i = -d_num_ghosts[0]; i < interior_dims[0] + d_num_ghosts[0]; i++)
                                            {
                                                // Compute index into linear data array.
                                                const int idx = (i + d_num_ghosts[0]) +
                                                    (j + d_num_ghosts[1])*ghostcell_dims[0];
                                                
                                                std::vector<const double*> Z_rho_ptr;
                                                for (int si = 0; si < d_num_species; si++)
                                                {
                                                    Z_rho_ptr.push_back(&Z_rho[si][idx]);
                                                }
                                                
                                                rho[idx] = d_equation_of_state->
                                                    getTotalDensity(
                                                        Z_rho_ptr);
                                            }
                                        }
                                    }
                                    else if (d_dim == tbox::Dimension(3))
                                    {
                                        for (int k = -d_num_ghosts[2]; k < interior_dims[2] + d_num_ghosts[2]; k++)
                                        {
                                            for (int j = -d_num_ghosts[1]; j < interior_dims[1] + d_num_ghosts[1]; j++)
                                            {
                                                for (int i = -d_num_ghosts[0]; i < interior_dims[0] + d_num_ghosts[0]; i++)
                                                {
                                                    // Compute index into linear data array.
                                                    const int idx = (i + d_num_ghosts[0]) +
                                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                                    
                                                    std::vector<const double*> Z_rho_ptr;
                                                    for (int si = 0; si < d_num_species; si++)
                                                    {
                                                        Z_rho_ptr.push_back(&Z_rho[si][idx]);
                                                    }
                                                    
                                                    rho[idx] = d_equation_of_state->
                                                        getTotalDensity(
                                                            Z_rho_ptr);
                                                }
                                            }
                                        }
                                    }
                                    
                                    // Get the wavelet coefficients.
                                    std::vector<boost::shared_ptr<pdat::CellData<double> > > wavelet_coeffs;
                                    for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                                    {
                                        wavelet_coeffs.push_back(BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                            patch.getPatchData(d_Harten_density_wavelet_coeffs[li], data_context)));
                                    }
                                    
                                    // Compute the wavelet coefficients.
                                    if (d_Harten_wavelet_uses_local_tol)
                                    {
                                        // Get the local means.
                                        std::vector<boost::shared_ptr<pdat::CellData<double> > > variable_local_means;
                                        for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                                        {
                                            variable_local_means.push_back(BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                                patch.getPatchData(d_Harten_density_local_means[li], data_context)));
                                        }
                                        
                                        d_wavelet_transfrom_Harten->computeWaveletCoefficientsWithVariableLocalMeans(
                                            patch,
                                            density,
                                            wavelet_coeffs,
                                            variable_local_means);
                                    }
                                    else
                                    {
                                        d_wavelet_transfrom_Harten->computeWaveletCoefficients(
                                            patch,
                                            density,
                                            wavelet_coeffs);
                                    }
                                    
                                    break;
                                }
                                default:
                                {
                                    TBOX_ERROR(d_object_name
                                        << ": "
                                        << "d_flow_model '"
                                        << d_flow_model
                                        << "' not yet implemented."
                                        << std::endl);
                                }
                            }
                        }
                        else if (variable_key ==  "TOTAL_ENERGY")
                        {
                            switch (d_flow_model)
                            {
                                case SINGLE_SPECIES: case FOUR_EQN_SHYUE: case FIVE_EQN_ALLAIRE:
                                {
                                    // Get the cell data of total energy.
                                    boost::shared_ptr<pdat::CellData<double> > total_energy(
                                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                            patch.getPatchData(d_total_energy, data_context)));
                                    
                                    // Get the wavelet coefficients.
                                    std::vector<boost::shared_ptr<pdat::CellData<double> > > wavelet_coeffs;
                                    for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                                    {
                                        wavelet_coeffs.push_back(BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                            patch.getPatchData(d_Harten_total_energy_wavelet_coeffs[li], data_context)));
                                    }
                                    
                                    // Compute the wavelet coefficients.
                                    if (d_Harten_wavelet_uses_local_tol)
                                    {
                                        // Get the local means.
                                        std::vector<boost::shared_ptr<pdat::CellData<double> > > variable_local_means;
                                        for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                                        {
                                            variable_local_means.push_back(BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                                patch.getPatchData(d_Harten_total_energy_local_means[li], data_context)));
                                        }
                                        
                                        d_wavelet_transfrom_Harten->computeWaveletCoefficientsWithVariableLocalMeans(
                                            patch,
                                            total_energy,
                                            wavelet_coeffs,
                                            variable_local_means);
                                    }
                                    else
                                    {
                                        d_wavelet_transfrom_Harten->computeWaveletCoefficients(
                                            patch,
                                            total_energy,
                                            wavelet_coeffs);
                                    }
                                    
                                    break;
                                }
                                default:
                                {
                                    TBOX_ERROR(d_object_name
                                        << ": "
                                        << "d_flow_model '"
                                        << d_flow_model
                                        << "' not yet implemented."
                                        << std::endl);
                                }
                            }
                        }
                        else if (variable_key == "PRESSURE")
                        {
                            switch (d_flow_model)
                            {
                                case SINGLE_SPECIES: case FOUR_EQN_SHYUE: case FIVE_EQN_ALLAIRE:
                                {
                                    // Get the cell data of density, momentum and total energy.
                                    boost::shared_ptr<pdat::CellData<double> > density(
                                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                            patch.getPatchData(d_density, data_context)));
                                    
                                    boost::shared_ptr<pdat::CellData<double> > momentum(
                                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                            patch.getPatchData(d_momentum, data_context)));
                                    
                                    boost::shared_ptr<pdat::CellData<double> > total_energy(
                                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                            patch.getPatchData(d_total_energy, data_context)));
                                    
                                    // Allocate temporary cell data for pressure.
                                    boost::shared_ptr<pdat::CellData<double> > pressure(
                                        new pdat::CellData<double>(interior_box, 1, d_num_ghosts));
                                    
                                    // Get the pointers to density, total energy and pressure.
                                    double* rho   = density->getPointer(0);
                                    double* E     = total_energy->getPointer(0);
                                    double* p     = pressure->getPointer(0);
                                    
                                    // Compute the field of pressure.
                                    if (d_dim == tbox::Dimension(1))
                                    {
                                        // NOT YET IMPLEMENTED
                                    }
                                    else if (d_dim == tbox::Dimension(2))
                                    {
                                        // Get the pointers to momentum components.
                                        double* rho_u = momentum->getPointer(0);
                                        double* rho_v = momentum->getPointer(1);
                                        
                                        for (int j = -d_num_ghosts[1]; j < interior_dims[1] + d_num_ghosts[1]; j++)
                                        {
                                            for (int i = -d_num_ghosts[0]; i < interior_dims[0] + d_num_ghosts[0]; i++)
                                            {
                                                // Compute index into linear data array.
                                                const int idx = (i + d_num_ghosts[0]) +
                                                    (j + d_num_ghosts[1])*ghostcell_dims[0];
                                                
                                                std::vector<const double*> m_ptr;
                                                m_ptr.push_back(&rho_u[idx]);
                                                m_ptr.push_back(&rho_v[idx]);
                                                
                                                p[idx] = d_equation_of_state->getPressure(
                                                    &rho[idx],
                                                    m_ptr,
                                                    &E[idx]);
                                            }
                                        }
                                    }
                                    else if (d_dim == tbox::Dimension(3))
                                    {
                                        // Get the pointers to momentum components.
                                        double* rho_u = momentum->getPointer(0);
                                        double* rho_v = momentum->getPointer(1);
                                        double* rho_w = momentum->getPointer(2);
                                        
                                        for (int k = -d_num_ghosts[2]; k < interior_dims[2] + d_num_ghosts[2]; k++)
                                        {
                                            for (int j = -d_num_ghosts[1]; j < interior_dims[1] + d_num_ghosts[1]; j++)
                                            {
                                                for (int i = -d_num_ghosts[0]; i < interior_dims[0] + d_num_ghosts[0]; i++)
                                                {
                                                    // Compute index into linear data array.
                                                    const int idx = (i + d_num_ghosts[0]) +
                                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                                    
                                                    std::vector<const double*> m_ptr;
                                                    m_ptr.push_back(&rho_u[idx]);
                                                    m_ptr.push_back(&rho_v[idx]);
                                                    m_ptr.push_back(&rho_w[idx]);
                                                    
                                                    p[idx] = d_equation_of_state->getPressure(
                                                        &rho[idx],
                                                        m_ptr,
                                                        &E[idx]);
                                                }
                                            }
                                        }
                                    }
                                    
                                    // Get the wavelet coefficients.
                                    std::vector<boost::shared_ptr<pdat::CellData<double> > > wavelet_coeffs;
                                    for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                                    {
                                        wavelet_coeffs.push_back(BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                            patch.getPatchData(d_Harten_pressure_wavelet_coeffs[li], data_context)));
                                    }
                                    
                                    // Compute the wavelet coefficients.
                                    if (d_Harten_wavelet_uses_local_tol)
                                    {
                                        // Get the local means.
                                        std::vector<boost::shared_ptr<pdat::CellData<double> > > variable_local_means;
                                        for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                                        {
                                            variable_local_means.push_back(BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                                patch.getPatchData(d_Harten_pressure_local_means[li], data_context)));
                                        }
                                        
                                        d_wavelet_transfrom_Harten->computeWaveletCoefficientsWithVariableLocalMeans(
                                            patch,
                                            pressure,
                                            wavelet_coeffs,
                                            variable_local_means);
                                    }
                                    else
                                    {
                                        d_wavelet_transfrom_Harten->computeWaveletCoefficients(
                                            patch,
                                            pressure,
                                            wavelet_coeffs);
                                    }
                                    
                                    break;
                                }
                                default:
                                {
                                    TBOX_ERROR(d_object_name
                                        << ": "
                                        << "d_flow_model '"
                                        << d_flow_model
                                        << "' not yet implemented."
                                        << std::endl);
                                }
                            }
                        }
                        else if (variable_key == "ENSTROPHY")
                        {
                            switch (d_flow_model)
                            {
                                case SINGLE_SPECIES: case FOUR_EQN_SHYUE: case FIVE_EQN_ALLAIRE:
                                {
                                    // Get the cell data of density, momentum and total energy.
                                    boost::shared_ptr<pdat::CellData<double> > density(
                                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                            patch.getPatchData(d_density, data_context)));
                                    
                                    boost::shared_ptr<pdat::CellData<double> > momentum(
                                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                            patch.getPatchData(d_momentum, data_context)));
                                    
                                    // Allocate temporary cell data for velocity.
                                    boost::shared_ptr<pdat::CellData<double> > velocity(
                                        new pdat::CellData<double>(interior_box, d_dim.getValue(), d_num_ghosts));
                                    
                                    // Allocate temporary cell data for enstrophy.
                                    boost::shared_ptr<pdat::CellData<double> > enstrophy(
                                        new pdat::CellData<double>(interior_box, 1, d_num_ghosts));
                                    
                                    // Get the pointer to density.
                                    double* rho   = density->getPointer(0);
                                    
                                    // Get the pointer to enstrophy.
                                    double* Omega = enstrophy->getPointer(0);
                                    
                                    // Compute the enstrophy.
                                    if (d_dim == tbox::Dimension(1))
                                    {
                                        // NOT YET IMPLEMENTED
                                    }
                                    else if (d_dim == tbox::Dimension(2))
                                    {
                                        // Get the pointers to momentum components.
                                        double* rho_u = momentum->getPointer(0);
                                        double* rho_v = momentum->getPointer(1);
                                        
                                        // Get the pointers to velocity components.
                                        double* u = velocity->getPointer(0);
                                        double* v = velocity->getPointer(1);
                                        
                                        // Compute the velocity.
                                        for (int j = -d_num_ghosts[1]; j < interior_dims[1] + d_num_ghosts[1]; j++)
                                        {
                                            for (int i = -d_num_ghosts[0]; i < interior_dims[0] + d_num_ghosts[0]; i++)
                                            {
                                                // Compute index into linear data array.
                                                const int idx = (i + d_num_ghosts[0]) +
                                                    (j + d_num_ghosts[1])*ghostcell_dims[0];
                                                
                                                u[idx] = rho_u[idx]/rho[idx];
                                                v[idx] = rho_v[idx]/rho[idx];
                                            }
                                        }
                                        
                                        // Compute the enstrophy.
                                        for (int j = -d_num_ghosts[1]; j < interior_dims[1] + d_num_ghosts[1]; j++)
                                        {
                                            for (int i = -d_num_ghosts[0]; i < interior_dims[0] + d_num_ghosts[0]; i++)
                                            {
                                                // Compute indices.
                                                const int idx = (i + d_num_ghosts[0]) +
                                                    (j + d_num_ghosts[1])*ghostcell_dims[0];
                                                
                                                const int idx_x_L = (i - 1 + d_num_ghosts[0]) +
                                                    (j + d_num_ghosts[1])*ghostcell_dims[0];
                                                
                                                const int idx_x_R = (i + 1 + d_num_ghosts[0]) +
                                                    (j + d_num_ghosts[1])*ghostcell_dims[0];
                                                
                                                const int idx_y_B = (i + d_num_ghosts[0]) +
                                                    (j - 1 + d_num_ghosts[1])*ghostcell_dims[0];
                                                
                                                const int idx_y_T = (i + d_num_ghosts[0]) +
                                                    (j + 1 + d_num_ghosts[1])*ghostcell_dims[0];
                                                
                                                double dudy, dvdx;
                                                
                                                if (i == -d_num_ghosts[0])
                                                {
                                                    dvdx = (v[idx_x_R] - v[idx])/(dx[0]);
                                                }
                                                else if (i == interior_dims[0] + d_num_ghosts[0] - 1)
                                                {
                                                    dvdx = (v[idx] - v[idx_x_L])/(dx[0]);
                                                }
                                                else
                                                {
                                                    dvdx = (v[idx_x_R] - v[idx_x_L])/(2*dx[0]);
                                                }
                                                
                                                if (j == -d_num_ghosts[1])
                                                {
                                                    dudy = (u[idx_y_T] - u[idx])/(dx[1]);
                                                }
                                                else if (j == interior_dims[1] + d_num_ghosts[1] - 1)
                                                {
                                                    dudy = (u[idx] - u[idx_y_B])/(dx[1]);
                                                }
                                                else
                                                {
                                                    dudy = (u[idx_y_T] - u[idx_y_B])/(2*dx[1]);
                                                }
                                                
                                                Omega[idx] = pow(dvdx - dudy, 2);
                                            }
                                        }
                                    }
                                    else if (d_dim == tbox::Dimension(3))
                                    {
                                        // Get the pointers to momentum components.
                                        double* rho_u = momentum->getPointer(0);
                                        double* rho_v = momentum->getPointer(1);
                                        double* rho_w = momentum->getPointer(2);
                                        
                                        // Get the pointers to velocity components.
                                        double* u = velocity->getPointer(0);
                                        double* v = velocity->getPointer(1);
                                        double* w = velocity->getPointer(2);
                                        
                                        // Compute the velocity.
                                        for (int k = -d_num_ghosts[2]; k < interior_dims[2] + d_num_ghosts[2]; k++)
                                        {
                                            for (int j = -d_num_ghosts[1]; j < interior_dims[1] + d_num_ghosts[1]; j++)
                                            {
                                                for (int i = -d_num_ghosts[0]; i < interior_dims[0] + d_num_ghosts[0]; i++)
                                                {
                                                    // Compute index into linear data array.
                                                    const int idx = (i + d_num_ghosts[0]) +
                                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                                    
                                                    u[idx] = rho_u[idx]/rho[idx];
                                                    v[idx] = rho_v[idx]/rho[idx];
                                                    w[idx] = rho_w[idx]/rho[idx];
                                                }
                                            }
                                        }
                                        
                                        // Compute the enstrophy.
                                        for (int k = -d_num_ghosts[2]; k < interior_dims[2] + d_num_ghosts[2]; k++)
                                        {
                                            for (int j = -d_num_ghosts[1]; j < interior_dims[1] + d_num_ghosts[1]; j++)
                                            {
                                                for (int i = -d_num_ghosts[0]; i < interior_dims[0] + d_num_ghosts[0]; i++)
                                                {
                                                    // Compute indices.
                                                    const int idx = (i + d_num_ghosts[0]) +
                                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                                    
                                                    const int idx_x_L = (i - 1 + d_num_ghosts[0]) +
                                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                                    
                                                    const int idx_x_R = (i + 1 + d_num_ghosts[0]) +
                                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                                    
                                                    const int idx_y_B = (i + d_num_ghosts[0]) +
                                                        (j - 1 + d_num_ghosts[1])*ghostcell_dims[0] +
                                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                                    
                                                    const int idx_y_T = (i + d_num_ghosts[0]) +
                                                        (j + 1 + d_num_ghosts[1])*ghostcell_dims[0] +
                                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                                    
                                                    const int idx_z_B = (i + d_num_ghosts[0]) +
                                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                                        (k - 1 + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                                    
                                                    const int idx_z_F = (i + d_num_ghosts[0]) +
                                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                                        (k + 1 + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                                    
                                                    double dudy, dudz, dvdx, dvdz, dwdx, dwdy;
                                                    
                                                    if (i == -d_num_ghosts[0])
                                                    {
                                                        dvdx = (v[idx_x_R] - v[idx])/(dx[0]);
                                                        dwdx = (w[idx_x_R] - w[idx])/(dx[0]);
                                                    }
                                                    else if (i == interior_dims[0] + d_num_ghosts[0] - 1)
                                                    {
                                                        dvdx = (v[idx] - v[idx_x_L])/(dx[0]);
                                                        dwdx = (w[idx] - w[idx_x_L])/(dx[0]);
                                                    }
                                                    else
                                                    {
                                                        dvdx = (v[idx_x_R] - v[idx_x_L])/(2*dx[0]);
                                                        dwdx = (w[idx_x_R] - w[idx_x_L])/(2*dx[0]);
                                                    }
                                                    
                                                    if (j == -d_num_ghosts[1])
                                                    {
                                                        dudy = (u[idx_y_T] - u[idx])/(dx[1]);
                                                        dwdy = (w[idx_y_T] - w[idx])/(dx[1]);
                                                    }
                                                    else if (j == interior_dims[1] + d_num_ghosts[1] - 1)
                                                    {
                                                        dudy = (u[idx] - u[idx_y_B])/(dx[1]);
                                                        dwdy = (w[idx] - w[idx_y_B])/(dx[1]);
                                                    }
                                                    else
                                                    {
                                                        dudy = (u[idx_y_T] - u[idx_y_B])/(2*dx[1]);
                                                        dwdy = (w[idx_y_T] - w[idx_y_B])/(2*dx[1]);
                                                    }
                                                    
                                                    if (k == -d_num_ghosts[2])
                                                    {
                                                        dudz = (u[idx_z_F] - u[idx])/(dx[2]);
                                                        dvdz = (v[idx_z_F] - v[idx])/(dx[2]);
                                                    }
                                                    else if (k == interior_dims[2] + d_num_ghosts[2] - 1)
                                                    {
                                                        dudz = (u[idx] - u[idx_z_B])/(dx[2]);
                                                        dvdz = (v[idx] - v[idx_z_B])/(dx[2]);
                                                    }
                                                    else
                                                    {
                                                        dudz = (u[idx_z_F] - u[idx_z_B])/(2*dx[2]);
                                                        dvdz = (v[idx_z_F] - v[idx_z_B])/(2*dx[2]);
                                                    }
                                                    
                                                    Omega[idx] = pow(dwdy - dvdz, 2) +
                                                        pow(dudz - dwdx, 2) +
                                                        pow(dvdx - dudy, 2);
                                                }
                                            }
                                        }
                                    }
                                    
                                    // Get the wavelet coefficients.
                                    std::vector<boost::shared_ptr<pdat::CellData<double> > > wavelet_coeffs;
                                    for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                                    {
                                        wavelet_coeffs.push_back(BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                            patch.getPatchData(d_Harten_enstrophy_wavelet_coeffs[li], data_context)));
                                    }
                                    
                                    // Compute the wavelet coefficients.
                                    if (d_Harten_wavelet_uses_local_tol)
                                    {
                                        // Get the local means.
                                        std::vector<boost::shared_ptr<pdat::CellData<double> > > variable_local_means;
                                        for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                                        {
                                            variable_local_means.push_back(BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                                patch.getPatchData(d_Harten_enstrophy_local_means[li], data_context)));
                                        }
                                        
                                        d_wavelet_transfrom_Harten->computeWaveletCoefficientsWithVariableLocalMeans(
                                            patch,
                                            enstrophy,
                                            wavelet_coeffs,
                                            variable_local_means);
                                    }
                                    else
                                    {
                                        d_wavelet_transfrom_Harten->computeWaveletCoefficients(
                                            patch,
                                            enstrophy,
                                            wavelet_coeffs);
                                    }
                                    
                                    break;
                                }
                                default:
                                {
                                    TBOX_ERROR(d_object_name
                                        << ": "
                                        << "d_flow_model '"
                                        << d_flow_model
                                        << "' not yet implemented."
                                        << std::endl);
                                }
                            }
                        }
                        else if (variable_key == "MASS_FRACTION")
                        {
                            switch (d_flow_model)
                            {
                                case SINGLE_SPECIES:
                                {
                                    TBOX_ERROR(d_object_name
                                        << ": '"
                                        << variable_key
                                        << "' not supported for '"
                                        << d_flow_model
                                        << "' flow model."
                                        << std::endl);
                                    
                                    break;
                                }
                                case FOUR_EQN_SHYUE:
                                {
                                    // Get the cell data of mass fraction.
                                    boost::shared_ptr<pdat::CellData<double> > mass_fraction(
                                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                            patch.getPatchData(d_mass_fraction, data_context)));
                                    
                                    for (int di = 0; di < d_mass_fraction->getDepth(); di++)
                                    {
                                        // Get the wavelet coefficients.
                                        std::vector<boost::shared_ptr<pdat::CellData<double> > > wavelet_coeffs;
                                        for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                                        {
                                            wavelet_coeffs.push_back(BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                                patch.getPatchData(d_Harten_mass_fraction_wavelet_coeffs[di*d_Harten_wavelet_num_level + li],
                                                    data_context)));
                                        }
                                        
                                        // Compute the wavelet coefficients.
                                        if (d_Harten_wavelet_uses_local_tol)
                                        {
                                            // Get the local means.
                                            std::vector<boost::shared_ptr<pdat::CellData<double> > > variable_local_means;
                                            for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                                            {
                                                variable_local_means.push_back(BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                                    patch.getPatchData(
                                                        d_Harten_mass_fraction_local_means[di*d_Harten_wavelet_num_level + li],
                                                        data_context)));
                                            }
                                            
                                            d_wavelet_transfrom_Harten->computeWaveletCoefficientsWithVariableLocalMeans(
                                                patch,
                                                mass_fraction,
                                                wavelet_coeffs,
                                                variable_local_means,
                                                di);
                                        }
                                        else
                                        {
                                            d_wavelet_transfrom_Harten->computeWaveletCoefficients(
                                                patch,
                                                mass_fraction,
                                                wavelet_coeffs,
                                                di);
                                        }
                                    }
                                    
                                    break;
                                }
                                case FIVE_EQN_ALLAIRE:
                                {
                                    // Get the cell data of partial density.
                                    boost::shared_ptr<pdat::CellData<double> > partial_density(
                                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                            patch.getPatchData(d_partial_density, data_context)));
                                    
                                    // Allocate temporary cell data for mass fraction.
                                    boost::shared_ptr<pdat::CellData<double> > mass_fraction(
                                        new pdat::CellData<double>(interior_box, d_num_species, d_num_ghosts));
                                    
                                    // Get the pointers to partial density and mass fraction.
                                    std::vector<double*> Z_rho;
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        Z_rho.push_back(partial_density->getPointer(si));
                                    }
                                    std::vector<double*> Y;
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        Y.push_back(mass_fraction->getPointer(si));
                                    }
                                    
                                    // Compute the field of pressure.
                                    if (d_dim == tbox::Dimension(1))
                                    {
                                        // NOT YET IMPLEMENTED
                                    }
                                    else if (d_dim == tbox::Dimension(2))
                                    {
                                        for (int j = -d_num_ghosts[1]; j < interior_dims[1] + d_num_ghosts[1]; j++)
                                        {
                                            for (int i = -d_num_ghosts[0]; i < interior_dims[0] + d_num_ghosts[0]; i++)
                                            {
                                                // Compute index into linear data array.
                                                const int idx = (i + d_num_ghosts[0]) +
                                                    (j + d_num_ghosts[1])*ghostcell_dims[0];
                                                
                                                std::vector<const double*> Z_rho_ptr;
                                                for (int si = 0; si < d_num_species; si++)
                                                {
                                                    Z_rho_ptr.push_back(&Z_rho[si][idx]);
                                                }
                                                
                                                const double rho = d_equation_of_state->
                                                    getTotalDensity(
                                                        Z_rho_ptr);
                                                
                                                for (int si = 0; si < d_num_species; si++)
                                                {
                                                    Y[si][idx] = Z_rho[si][idx]/rho;
                                                }
                                            }
                                        }
                                    }
                                    else if (d_dim == tbox::Dimension(3))
                                    {
                                        for (int k = -d_num_ghosts[2]; k < interior_dims[2] + d_num_ghosts[2]; k++)
                                        {
                                            for (int j = -d_num_ghosts[1]; j < interior_dims[1] + d_num_ghosts[1]; j++)
                                            {
                                                for (int i = -d_num_ghosts[0]; i < interior_dims[0] + d_num_ghosts[0]; i++)
                                                {
                                                    // Compute index into linear data array.
                                                    const int idx = (i + d_num_ghosts[0]) +
                                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                                    
                                                    std::vector<const double*> Z_rho_ptr;
                                                    for (int si = 0; si < d_num_species; si++)
                                                    {
                                                        Z_rho_ptr.push_back(&Z_rho[si][idx]);
                                                    }
                                                    
                                                    const double rho = d_equation_of_state->
                                                        getTotalDensity(
                                                            Z_rho_ptr);
                                                    
                                                    for (int si = 0; si < d_num_species; si++)
                                                    {
                                                        Y[si][idx] = Z_rho[si][idx]/rho;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    
                                    for (int di = 0; di < d_num_species; di++)
                                    {
                                        // Get the wavelet coefficients.
                                        std::vector<boost::shared_ptr<pdat::CellData<double> > > wavelet_coeffs;
                                        for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                                        {
                                            wavelet_coeffs.push_back(BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                                patch.getPatchData(d_Harten_mass_fraction_wavelet_coeffs[di*d_Harten_wavelet_num_level + li],
                                                    data_context)));
                                        }
                                        
                                        // Compute the wavelet coefficients.
                                        if (d_Harten_wavelet_uses_local_tol)
                                        {
                                            // Get the local means.
                                            std::vector<boost::shared_ptr<pdat::CellData<double> > > variable_local_means;
                                            for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                                            {
                                                variable_local_means.push_back(BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                                    patch.getPatchData(
                                                        d_Harten_mass_fraction_local_means[di*d_Harten_wavelet_num_level + li],
                                                        data_context)));
                                            }
                                            
                                            d_wavelet_transfrom_Harten->computeWaveletCoefficientsWithVariableLocalMeans(
                                                patch,
                                                mass_fraction,
                                                wavelet_coeffs,
                                                variable_local_means,
                                                di);
                                        }
                                        else
                                        {
                                            d_wavelet_transfrom_Harten->computeWaveletCoefficients(
                                                patch,
                                                mass_fraction,
                                                wavelet_coeffs,
                                                di);
                                        }
                                    }
                                    
                                    break;
                                }
                                default:
                                {
                                    TBOX_ERROR(d_object_name
                                        << ": "
                                        << "d_flow_model '"
                                        << d_flow_model
                                        << "' not yet implemented."
                                        << std::endl);
                                }
                            }
                        }
                        else if (variable_key == "VOLUME_FRACTION")
                        {
                            switch (d_flow_model)
                            {
                                case SINGLE_SPECIES: case FOUR_EQN_SHYUE:
                                {
                                    TBOX_ERROR(d_object_name
                                        << ": '"
                                        << variable_key
                                        << "' not supported for '"
                                        << d_flow_model
                                        << "' flow model."
                                        << std::endl);
                                    
                                    break;
                                }
                                case FIVE_EQN_ALLAIRE:
                                {
                                    // Get the cell data of volume fraction.
                                    boost::shared_ptr<pdat::CellData<double> > volume_fraction(
                                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                            patch.getPatchData(d_volume_fraction, data_context)));
                                    
                                    for (int di = 0; di < d_volume_fraction->getDepth(); di++)
                                    {
                                        // Get the wavelet coefficients.
                                        std::vector<boost::shared_ptr<pdat::CellData<double> > > wavelet_coeffs;
                                        for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                                        {
                                            wavelet_coeffs.push_back(BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                                patch.getPatchData(d_Harten_volume_fraction_wavelet_coeffs[di*d_Harten_wavelet_num_level + li],
                                                    data_context)));
                                        }
                                        
                                        // Compute the wavelet coefficients.
                                        if (d_Harten_wavelet_uses_local_tol)
                                        {
                                            // Get the local means.
                                            std::vector<boost::shared_ptr<pdat::CellData<double> > > variable_local_means;
                                            for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                                            {
                                                variable_local_means.push_back(BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                                    patch.getPatchData(
                                                        d_Harten_volume_fraction_local_means[di*d_Harten_wavelet_num_level + li],
                                                        data_context)));
                                            }
                                            
                                            d_wavelet_transfrom_Harten->computeWaveletCoefficientsWithVariableLocalMeans(
                                                patch,
                                                volume_fraction,
                                                wavelet_coeffs,
                                                variable_local_means,
                                                di);
                                        }
                                        else
                                        {
                                            d_wavelet_transfrom_Harten->computeWaveletCoefficients(
                                                patch,
                                                volume_fraction,
                                                wavelet_coeffs,
                                                di);
                                        }
                                    }
                                    
                                    break;
                                }
                                default:
                                {
                                    TBOX_ERROR(d_object_name
                                        << ": "
                                        << "d_flow_model '"
                                        << d_flow_model
                                        << "' not yet implemented."
                                        << std::endl);
                                }
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
                    } // Loop over variables.
                }
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Number of ghost cells is not set yet."
                << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Variables are not set yet."
            << std::endl);
    }
}


/*
 * Get the statistics of the sensor values that are required by the
 * multiresolution sensors at a given patch level.
 */
void
MultiresolutionTagger::getSensorValueStatistics(
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const int level_number,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    math::HierarchyCellDataOpsReal<double> cell_double_operator(patch_hierarchy, level_number, level_number);
    
    hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();
    
    if (d_variables_set == true)
    {
        if (d_num_ghosts_set == true)
        {
            for (int si = 0;
                     si < static_cast<int>(d_multiresolution_sensors.size());
                     si++)
            {
                std::string sensor_key = d_multiresolution_sensors[si];
                
                if (sensor_key == "HARTEN_WAVELET")
                {
                    // Only get the statistics if global filter is used.
                    if (d_Harten_wavelet_uses_global_tol)
                    {
                        for (int vi = 0; vi < static_cast<int>(d_Harten_wavelet_variables.size()); vi++)
                        {
                            // Get the key of the current variable.
                            std::string variable_key = d_Harten_wavelet_variables[vi];
                            
                            if (variable_key == "DENSITY")
                            {
                                for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                                {
                                    const int rho_w_id = variable_db->mapVariableAndContextToIndex(
                                        d_Harten_density_wavelet_coeffs[li],
                                        data_context);
                                    
                                    double rho_w_max_local = cell_double_operator.max(rho_w_id);
                                    d_Harten_density_wavelet_coeffs_maxs[li] = 0.0;
                                    
                                    mpi.Allreduce(
                                        &rho_w_max_local,
                                        &d_Harten_density_wavelet_coeffs_maxs[li],
                                        1,
                                        MPI_DOUBLE,
                                        MPI_MAX);
                                }
                            }
                            else if (variable_key == "TOTAL_ENERGY")
                            {
                                for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                                {
                                    const int E_w_id = variable_db->mapVariableAndContextToIndex(
                                        d_Harten_total_energy_wavelet_coeffs[li],
                                        data_context);
                                    
                                    double E_w_max_local = cell_double_operator.max(E_w_id);
                                    d_Harten_total_energy_wavelet_coeffs_maxs[li] = 0.0;
                                    
                                    mpi.Allreduce(
                                        &E_w_max_local,
                                        &d_Harten_total_energy_wavelet_coeffs_maxs[li],
                                        1,
                                        MPI_DOUBLE,
                                        MPI_MAX);
                                }
                            }
                            else if (variable_key == "PRESSURE")
                            {
                                for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                                {
                                    const int p_w_id = variable_db->mapVariableAndContextToIndex(
                                        d_Harten_pressure_wavelet_coeffs[li],
                                        data_context);
                                    
                                    double p_w_max_local = cell_double_operator.max(p_w_id);
                                    d_Harten_pressure_wavelet_coeffs_maxs[li] = 0.0;
                                    
                                    mpi.Allreduce(
                                        &p_w_max_local,
                                        &d_Harten_pressure_wavelet_coeffs_maxs[li],
                                        1,
                                        MPI_DOUBLE,
                                        MPI_MAX);
                                }
                            }
                            else if (variable_key == "ENSTROPHY")
                            {
                                for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                                {
                                    const int Omega_w_id = variable_db->mapVariableAndContextToIndex(
                                        d_Harten_enstrophy_wavelet_coeffs[li],
                                        data_context);
                                    
                                    double Omega_w_max_local = cell_double_operator.max(Omega_w_id);
                                    d_Harten_enstrophy_wavelet_coeffs_maxs[li] = 0.0;
                                    
                                    mpi.Allreduce(
                                        &Omega_w_max_local,
                                        &d_Harten_enstrophy_wavelet_coeffs_maxs[li],
                                        1,
                                        MPI_DOUBLE,
                                        MPI_MAX);
                                }
                            }
                            else if (variable_key == "MASS_FRACTION")
                            {
                                switch (d_flow_model)
                                {
                                    case SINGLE_SPECIES:
                                    {
                                        TBOX_ERROR(d_object_name
                                            << ": '"
                                            << variable_key
                                            << "' not supported for '"
                                            << d_flow_model
                                            << "' flow model."
                                            << std::endl);
                                        
                                        break;
                                    }
                                    case FOUR_EQN_SHYUE:
                                    {
                                        for (int di = 0; di < d_mass_fraction->getDepth(); di++)
                                        {
                                            for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                                            {
                                                const int Y_w_id = variable_db->mapVariableAndContextToIndex(
                                                    d_Harten_mass_fraction_wavelet_coeffs[di*d_Harten_wavelet_num_level + li],
                                                    data_context);
                                                
                                                double Y_w_max_local = cell_double_operator.max(Y_w_id);
                                                d_Harten_mass_fraction_wavelet_coeffs_maxs[di*d_Harten_wavelet_num_level + li] = 0.0;
                                                
                                                mpi.Allreduce(
                                                    &Y_w_max_local,
                                                    &d_Harten_mass_fraction_wavelet_coeffs_maxs[di*d_Harten_wavelet_num_level + li],
                                                    1,
                                                    MPI_DOUBLE,
                                                    MPI_MAX);
                                            }
                                        }
                                        
                                        break;
                                    }
                                    case FIVE_EQN_ALLAIRE:
                                    {
                                        for (int di = 0; di < d_num_species; di++)
                                        {
                                            for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                                            {
                                                const int Y_w_id = variable_db->mapVariableAndContextToIndex(
                                                    d_Harten_mass_fraction_wavelet_coeffs[di*d_Harten_wavelet_num_level + li],
                                                    data_context);
                                                
                                                double Y_w_max_local = cell_double_operator.max(Y_w_id);
                                                d_Harten_mass_fraction_wavelet_coeffs_maxs[di*d_Harten_wavelet_num_level + li] = 0.0;
                                                
                                                mpi.Allreduce(
                                                    &Y_w_max_local,
                                                    &d_Harten_mass_fraction_wavelet_coeffs_maxs[di*d_Harten_wavelet_num_level + li],
                                                    1,
                                                    MPI_DOUBLE,
                                                    MPI_MAX);
                                            }
                                        }
                                        
                                        break;
                                    }
                                    default:
                                    {
                                        TBOX_ERROR(d_object_name
                                            << ": "
                                            << "d_flow_model '"
                                            << d_flow_model
                                            << "' not yet implemented."
                                            << std::endl);
                                    }
                                }
                            }
                            else if (variable_key == "VOLUME_FRACTION")
                            {
                                switch (d_flow_model)
                                {
                                    case SINGLE_SPECIES: case FOUR_EQN_SHYUE:
                                    {
                                        TBOX_ERROR(d_object_name
                                            << ": '"
                                            << variable_key
                                            << "' not supported for '"
                                            << d_flow_model
                                            << "' flow model."
                                            << std::endl);
                                        
                                        break;
                                    }
                                    case FIVE_EQN_ALLAIRE:
                                    {
                                        for (int di = 0; di < d_volume_fraction->getDepth(); di++)
                                        {
                                            for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                                            {
                                                const int Z_w_id = variable_db->mapVariableAndContextToIndex(
                                                    d_Harten_volume_fraction_wavelet_coeffs[di*d_Harten_wavelet_num_level + li],
                                                    data_context);
                                                
                                                double Z_w_max_local = cell_double_operator.max(Z_w_id);
                                                d_Harten_volume_fraction_wavelet_coeffs_maxs[di*d_Harten_wavelet_num_level + li] = 0.0;
                                                
                                                mpi.Allreduce(
                                                    &Z_w_max_local,
                                                    &d_Harten_volume_fraction_wavelet_coeffs_maxs[di*d_Harten_wavelet_num_level + li],
                                                    1,
                                                    MPI_DOUBLE,
                                                    MPI_MAX);
                                            }
                                        }
                                        
                                        break;
                                    }
                                    default:
                                    {
                                        TBOX_ERROR(d_object_name
                                            << ": "
                                            << "d_flow_model '"
                                            << d_flow_model
                                            << "' not yet implemented."
                                            << std::endl);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Number of ghost cells is not set yet."
                << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Variables are not set yet."
            << std::endl);
    }
}


/*
 * Tag cells for refinement using multiresolution sensors.
 */
void
MultiresolutionTagger::tagCells(
   hier::Patch& patch,
   boost::shared_ptr<pdat::CellData<int> > tags,
   const boost::shared_ptr<hier::VariableContext>& data_context)
{
    if (d_variables_set == true)
    {
        if (d_num_ghosts_set == true)
        {
            for (int si = 0;
                     si < static_cast<int>(d_multiresolution_sensors.size());
                     si++)
            {
                std::string sensor_key = d_multiresolution_sensors[si];
                
                if (sensor_key == "HARTEN_WAVELET")
                {
                    boost::shared_ptr<pdat::CellData<double> > Lipschitz_exponent(
                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                            patch.getPatchData(d_Harten_Lipschitz_exponent,
                            data_context)));
                    
                    // Loop over the variables for wavelet analysis.
                    for (int vi = 0; vi < static_cast<int>(d_Harten_wavelet_variables.size()); vi++)
                    {
                        // Get the key of the current variable.
                        std::string variable_key = d_Harten_wavelet_variables[vi];
                        
                        /*
                         * Get the wavelet coefficients and the statistics of the wavelet coefficients
                         * at different levels and tag cells using Lipschitz's exponent.
                         */
                        
                        double global_tol = 0.0;
                        double local_tol = 0.0;
                        double alpha_tol = 0.0;
                        
                        if (d_Harten_wavelet_uses_global_tol)
                        {
                            global_tol = d_Harten_wavelet_global_tol[vi];
                        }
                        if (d_Harten_wavelet_uses_local_tol)
                        {
                            local_tol = d_Harten_wavelet_local_tol[vi];
                        }
                        if (d_Harten_wavelet_uses_alpha_tol)
                        {
                            alpha_tol = d_Harten_wavelet_alpha_tol[vi];
                        }
                        
                        if (variable_key == "DENSITY")
                        {
                            std::vector<boost::shared_ptr<pdat::CellData<double> > > wavelet_coeffs;
                            
                            for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                            {
                                wavelet_coeffs.push_back(BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                    patch.getPatchData(d_Harten_density_wavelet_coeffs[li], data_context)));
                            }
                            
                            std::vector<boost::shared_ptr<pdat::CellData<double> > > variable_local_means;
                            if (d_Harten_wavelet_uses_local_tol)
                            {
                                // Get the local means.
                                for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                                {
                                    variable_local_means.push_back(BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                        patch.getPatchData(d_Harten_density_local_means[li], data_context)));
                                }
                            }
                            
                            tagCellsWithWaveletSensor(
                                patch,
                                tags,
                                wavelet_coeffs,
                                d_Harten_density_wavelet_coeffs_maxs,
                                variable_local_means,
                                Lipschitz_exponent,
                                global_tol,
                                local_tol,
                                alpha_tol,
                                sensor_key);
                        }
                        else if (variable_key == "TOTAL_ENERGY")
                        {
                            std::vector<boost::shared_ptr<pdat::CellData<double> > > wavelet_coeffs;
                            
                            for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                            {
                                wavelet_coeffs.push_back(BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                    patch.getPatchData(d_Harten_total_energy_wavelet_coeffs[li], data_context)));
                            }
                            
                            std::vector<boost::shared_ptr<pdat::CellData<double> > > variable_local_means;
                            if (d_Harten_wavelet_uses_local_tol)
                            {
                                // Get the local means.
                                for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                                {
                                    variable_local_means.push_back(BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                        patch.getPatchData(d_Harten_total_energy_local_means[li], data_context)));
                                }
                            }
                            
                            tagCellsWithWaveletSensor(
                                patch,
                                tags,
                                wavelet_coeffs,
                                d_Harten_total_energy_wavelet_coeffs_maxs,
                                variable_local_means,
                                Lipschitz_exponent,
                                global_tol,
                                local_tol,
                                alpha_tol,
                                sensor_key);
                        }
                        else if (variable_key == "PRESSURE")
                        {
                            std::vector<boost::shared_ptr<pdat::CellData<double> > > wavelet_coeffs;
                            
                            for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                            {
                                wavelet_coeffs.push_back(BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                    patch.getPatchData(d_Harten_pressure_wavelet_coeffs[li], data_context)));
                            }
                            
                            std::vector<boost::shared_ptr<pdat::CellData<double> > > variable_local_means;
                            if (d_Harten_wavelet_uses_local_tol)
                            {
                                // Get the local means.
                                for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                                {
                                    variable_local_means.push_back(BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                        patch.getPatchData(d_Harten_pressure_local_means[li], data_context)));
                                }
                            }
                            
                            tagCellsWithWaveletSensor(
                                patch,
                                tags,
                                wavelet_coeffs,
                                d_Harten_pressure_wavelet_coeffs_maxs,
                                variable_local_means,
                                Lipschitz_exponent,
                                global_tol,
                                local_tol,
                                alpha_tol,
                                sensor_key);
                        }
                        else if (variable_key == "ENSTROPHY")
                        {
                            std::vector<boost::shared_ptr<pdat::CellData<double> > > wavelet_coeffs;
                            
                            for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                            {
                                wavelet_coeffs.push_back(BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                    patch.getPatchData(d_Harten_enstrophy_wavelet_coeffs[li], data_context)));
                            }
                            
                            std::vector<boost::shared_ptr<pdat::CellData<double> > > variable_local_means;
                            if (d_Harten_wavelet_uses_local_tol)
                            {
                                // Get the local means.
                                for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                                {
                                    variable_local_means.push_back(BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                        patch.getPatchData(d_Harten_enstrophy_local_means[li], data_context)));
                                }
                            }
                            
                            tagCellsWithWaveletSensor(
                                patch,
                                tags,
                                wavelet_coeffs,
                                d_Harten_enstrophy_wavelet_coeffs_maxs,
                                variable_local_means,
                                Lipschitz_exponent,
                                global_tol,
                                local_tol,
                                alpha_tol,
                                sensor_key);
                        }
                        else if (variable_key == "MASS_FRACTION")
                        {
                            switch (d_flow_model)
                            {
                                case SINGLE_SPECIES:
                                {
                                    TBOX_ERROR(d_object_name
                                        << ": '"
                                        << variable_key
                                        << "' not supported for '"
                                        << d_flow_model
                                        << "' flow model."
                                        << std::endl);
                                    
                                    break;
                                }
                                case FOUR_EQN_SHYUE:
                                {
                                    for (int di = 0; di < d_mass_fraction->getDepth(); di++)
                                    {
                                        std::vector<boost::shared_ptr<pdat::CellData<double> > > wavelet_coeffs;
                                        std::vector<double> wavelet_coeffs_maxs;
                                        std::vector<boost::shared_ptr<pdat::CellData<double> > > variable_local_means;
                                        
                                        for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                                        {
                                            wavelet_coeffs.push_back(BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                                patch.getPatchData(
                                                    d_Harten_mass_fraction_wavelet_coeffs[di*d_Harten_wavelet_num_level + li],
                                                    data_context)));
                                        }
                                        
                                        if (d_Harten_wavelet_uses_global_tol)
                                        {
                                            for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                                            {
                                                wavelet_coeffs_maxs.push_back(
                                                    d_Harten_mass_fraction_wavelet_coeffs_maxs[
                                                        di*d_Harten_wavelet_num_level + li]);
                                            }
                                        }
                                        
                                        if (d_Harten_wavelet_uses_local_tol)
                                        {
                                            for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                                            {
                                                variable_local_means.push_back(BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                                    patch.getPatchData(
                                                        d_Harten_mass_fraction_local_means[di*d_Harten_wavelet_num_level + li],
                                                        data_context)));
                                            }
                                        }
                                        
                                        tagCellsWithWaveletSensor(
                                            patch,
                                            tags,
                                            wavelet_coeffs,
                                            wavelet_coeffs_maxs,
                                            variable_local_means,
                                            Lipschitz_exponent,
                                            global_tol,
                                            local_tol,
                                            alpha_tol,
                                            sensor_key);
                                    }
                                    
                                    break;
                                }
                                case FIVE_EQN_ALLAIRE:
                                {
                                    for (int di = 0; di < d_num_species; di++)
                                    {
                                        std::vector<boost::shared_ptr<pdat::CellData<double> > > wavelet_coeffs;
                                        std::vector<double> wavelet_coeffs_maxs;
                                        std::vector<boost::shared_ptr<pdat::CellData<double> > > variable_local_means;
                                        
                                        if (d_Harten_wavelet_uses_global_tol)
                                        {
                                            for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                                            {
                                                wavelet_coeffs_maxs.push_back(
                                                    d_Harten_mass_fraction_wavelet_coeffs_maxs[
                                                        di*d_Harten_wavelet_num_level + li]);
                                            }
                                        }
                                        
                                        if (d_Harten_wavelet_uses_local_tol)
                                        {
                                            for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                                            {
                                                variable_local_means.push_back(BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                                    patch.getPatchData(
                                                        d_Harten_mass_fraction_local_means[di*d_Harten_wavelet_num_level + li],
                                                        data_context)));
                                            }
                                        }
                                        
                                        tagCellsWithWaveletSensor(
                                            patch,
                                            tags,
                                            wavelet_coeffs,
                                            wavelet_coeffs_maxs,
                                            variable_local_means,
                                            Lipschitz_exponent,
                                            global_tol,
                                            local_tol,
                                            alpha_tol,
                                            sensor_key);
                                    }
                                    
                                    break;
                                }
                                default:
                                {
                                    TBOX_ERROR(d_object_name
                                        << ": "
                                        << "d_flow_model '"
                                        << d_flow_model
                                        << "' not yet implemented."
                                        << std::endl);
                                }
                            }
                        }
                        else if (variable_key == "VOLUME_FRACTION")
                        {
                            switch (d_flow_model)
                            {
                                case SINGLE_SPECIES: case FOUR_EQN_SHYUE:
                                {
                                    TBOX_ERROR(d_object_name
                                        << ": '"
                                        << variable_key
                                        << "' not supported for '"
                                        << d_flow_model
                                        << "' flow model."
                                        << std::endl);
                                    
                                    break;
                                }
                                case FIVE_EQN_ALLAIRE:
                                {
                                    for (int di = 0; di < d_volume_fraction->getDepth(); di++)
                                    {
                                        std::vector<boost::shared_ptr<pdat::CellData<double> > > wavelet_coeffs;
                                        std::vector<double> wavelet_coeffs_maxs;
                                        std::vector<boost::shared_ptr<pdat::CellData<double> > > variable_local_means;
                                        
                                        for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                                        {
                                            wavelet_coeffs.push_back(BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                                patch.getPatchData(
                                                    d_Harten_volume_fraction_wavelet_coeffs[di*d_Harten_wavelet_num_level + li],
                                                    data_context)));
                                        }
                                        
                                        if (d_Harten_wavelet_uses_global_tol)
                                        {
                                            for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                                            {
                                                wavelet_coeffs_maxs.push_back(
                                                    d_Harten_volume_fraction_wavelet_coeffs_maxs[
                                                        di*d_Harten_wavelet_num_level + li]);
                                            }
                                        }
                                        
                                        if (d_Harten_wavelet_uses_local_tol)
                                        {
                                            for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                                            {
                                                variable_local_means.push_back(BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                                    patch.getPatchData(
                                                        d_Harten_volume_fraction_local_means[di*d_Harten_wavelet_num_level + li],
                                                        data_context)));
                                            }
                                        }
                                        
                                        tagCellsWithWaveletSensor(
                                            patch,
                                            tags,
                                            wavelet_coeffs,
                                            wavelet_coeffs_maxs,
                                            variable_local_means,
                                            Lipschitz_exponent,
                                            global_tol,
                                            local_tol,
                                            alpha_tol,
                                            sensor_key);
                                    }
                                    
                                    break;
                                }
                                default:
                                {
                                    TBOX_ERROR(d_object_name
                                        << ": "
                                        << "d_flow_model '"
                                        << d_flow_model
                                        << "' not yet implemented."
                                        << std::endl);
                                }
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
                    } // Loop over variables.
                }
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Number of ghost cells is not set yet."
                << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Variables are not set yet."
            << std::endl);
    }
}


/*
 * Compute the Lipschitz's exponent. There are two steps:
 * 1. Find the maximum wavelet coefficients in the domain of dependence.
 * 2. Compute Lipschitz's exponent.
 */
void
MultiresolutionTagger::computeLipschitzExponent(
    hier::Patch& patch,
    std::vector<boost::shared_ptr<pdat::CellData<double> > >& wavelet_coeffs,
    boost::shared_ptr<pdat::CellData<double> > Lipschitz_exponent,
    std::string& sensor_key)
{
    const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(patch_geom);
#endif
    
    // Get the dimensions of box that covers the interior of Patch.
    hier::Box dummy_box = patch.getBox();
    const hier::Box interior_box = dummy_box;
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of Patch plus
    // ghost cells.
    dummy_box.grow(d_num_ghosts);
    const hier::Box ghost_box = dummy_box;
    const hier::IntVector ghostcell_dims = ghost_box.numberCells();
    
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
        
        // Create a vector of maximum wavelet coefficients in domain of dependence.
        std::vector<boost::shared_ptr<pdat::CellData<double> > > wavelet_coeffs_local_max;
        for (int li = 0; li < d_Harten_wavelet_num_level; li++)
        {
            wavelet_coeffs_local_max.push_back(boost::make_shared<pdat::CellData<double> >(
                interior_box, 1, d_num_ghosts));
        }
        
        // Get the pointers to the maximum wavelet coefficients in the domain of dependence.
        std::vector<double*> r;
        for (int li = 0; li < d_Harten_wavelet_num_level; li++)
        {
            r.push_back(wavelet_coeffs_local_max[li]->getPointer(0));
        }
        
        // Get the stencil width of the wavelet transform.
        int p = d_wavelet_transfrom_Harten->getLeftStencilWidth();
        int q = d_wavelet_transfrom_Harten->getRightStencilWidth();
        
        for (int li = 0; li < d_Harten_wavelet_num_level; li++)
        {
            if (d_dim == tbox::Dimension(1))
            {
                // NOT YET IMPLEMENTED
            }
            else if (d_dim == tbox::Dimension(2))
            {
                for (int j = 0; j < interior_dims[1]; j++)
                {
                    for (int i = 0; i < interior_dims[0]; i++)
                    {
                        // Compute index into linear data array.
                        const int idx = (i + d_num_ghosts[0]) +
                            (j + d_num_ghosts[1])*ghostcell_dims[0];
                        
                        // Find the maximum wavelet coefficient over the domain
                        // of dependence.
                        r[li][idx] = 0.0;
                        for (int ii = -p*pow(2,li+1); ii <= q*pow(2,li+1); ii++)
                        {
                            // Compute the index.
                            const int idx_s = (i + ii + d_num_ghosts[0]) +
                                (j + d_num_ghosts[1])*ghostcell_dims[0];
                            
                            r[li][idx] = fmax(r[li][idx], w[li][idx_s]);
                        }
                        for (int jj = -p*pow(2,li+1); jj <= q*pow(2,li+1); jj++)
                        {
                            // Compute the index.
                            const int idx_s = (i + d_num_ghosts[0]) +
                                (j + jj + d_num_ghosts[1])*ghostcell_dims[0];
                            
                            r[li][idx] = fmax(r[li][idx], w[li][idx_s]);
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
                            // Compute index into linear data array.
                            const int idx = (i + d_num_ghosts[0]) +
                                (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            
                            // Find the maximum wavelet coefficient over the domain
                            // of dependence.
                            r[li][idx] = 0.0;
                            for (int ii = -p*pow(2,li+1); ii <= q*pow(2,li+1); ii++)
                            {
                                // Compute the index.
                                const int idx_s = (i + ii + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                    (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                
                                r[li][idx] = fmax(r[li][idx], w[li][idx_s]);
                            }
                            for (int jj = -p*pow(2,li+1); jj <= q*pow(2,li+1); jj++)
                            {
                                // Compute the index.
                                const int idx_s = (i + d_num_ghosts[0]) +
                                    (j + jj + d_num_ghosts[1])*ghostcell_dims[0] +
                                    (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                
                                r[li][idx] = fmax(r[li][idx], w[li][idx_s]);
                            }
                            for (int kk = -p*pow(2,li+1); kk <=q*pow(2,li+1); kk++)
                            {
                                // Compute the index.
                                const int idx_s = (i + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                    (k + kk + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                
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
                    // NOT YET IMPLEMENTED
                }
                else if (d_dim == tbox::Dimension(2))
                {
                    for (int j = 0; j < interior_dims[1]; j++)
                    {
                        for (int i = 0; i < interior_dims[0]; i++)
                        {
                            // Compute indices.
                            const int idx = (i + d_num_ghosts[0]) +
                                (j + d_num_ghosts[1])*ghostcell_dims[0];
                            
                            if ((r[0][idx] > 1.0e-8) &&
                                (r[1][idx] > 1.0e-8))
                            {
                                alpha[idx] = log2(r[1][idx]/r[0][idx]);
                                
                                alpha[idx] = fmin(
                                    alpha[idx],
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
                    for (int k = 0; k < interior_dims[2]; k++)
                    {
                        for (int j = 0; j < interior_dims[1]; j++)
                        {
                            for (int i = 0; i < interior_dims[0]; i++)
                            {
                                // Compute indices.
                                const int idx = (i + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                    (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                
                                if ((r[0][idx] > 1.0e-8) &&
                                    (r[1][idx] > 1.0e-8))
                                {
                                    alpha[idx] = log2(r[1][idx]/r[0][idx]);
                                    
                                    alpha[idx] = fmin(
                                        alpha[idx],
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
                    // NOT YET IMPLEMENTED
                }
                else if (d_dim == tbox::Dimension(2))
                {
                    for (int j = 0; j < interior_dims[1]; j++)
                    {
                        for (int i = 0; i < interior_dims[0]; i++)
                        {
                            // Compute indices.
                            const int idx = (i + d_num_ghosts[0]) +
                                (j + d_num_ghosts[1])*ghostcell_dims[0];
                            
                            if ((r[0][idx] > 1.0e-8) &&
                                (r[1][idx] > 1.0e-8) &&
                                (r[2][idx] > 1.0e-8))
                            {
                                alpha[idx] = 0.5*log2(r[2][idx]/r[0][idx]);
                                
                                alpha[idx] = fmin(
                                    alpha[idx],
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
                    for (int k = 0; k < interior_dims[2]; k++)
                    {
                        for (int j = 0; j < interior_dims[1]; j++)
                        {
                            for (int i = 0; i < interior_dims[0]; i++)
                            {
                                // Compute indices.
                                const int idx = (i + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                    (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                
                                if ((r[0][idx] > 1.0e-8) &&
                                    (r[1][idx] > 1.0e-8) &&
                                    (r[2][idx] > 1.0e-8))
                                {
                                    alpha[idx] = 0.5*log2(r[2][idx]/r[0][idx]);
                                    
                                    alpha[idx] = fmin(
                                        alpha[idx],
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
                    // NOT YET IMPLEMENTED
                }
                else if (d_dim == tbox::Dimension(2))
                {
                    for (int j = 0; j < interior_dims[1]; j++)
                    {
                        for (int i = 0; i < interior_dims[0]; i++)
                        {
                            // Compute indices.
                            const int idx = (i + d_num_ghosts[0]) +
                                (j + d_num_ghosts[1])*ghostcell_dims[0];
                            
                            if ((r[0][idx] > 1.0e-8) &&
                                (r[1][idx] > 1.0e-8) &&
                                (r[2][idx] > 1.0e-8) &&
                                (r[3][idx] > 1.0e-8))
                            {
                                alpha[idx] = (3.0*log2(r[3][idx]/r[0][idx]) +
                                    log2(r[2][idx]/r[1][idx]))/10.0;
                                
                                alpha[idx] = fmin(
                                    alpha[idx],
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
                    for (int k = 0; k < interior_dims[2]; k++)
                    {
                        for (int j = 0; j < interior_dims[1]; j++)
                        {
                            for (int i = 0; i < interior_dims[0]; i++)
                            {
                                // Compute indices.
                                const int idx = (i + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                    (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                
                                if ((r[0][idx] > 1.0e-8) &&
                                    (r[1][idx] > 1.0e-8) &&
                                    (r[2][idx] > 1.0e-8) &&
                                    (r[3][idx] > 1.0e-8))
                                {
                                    alpha[idx] = (3*log2(r[3][idx]/r[0][idx]) +
                                        log2(r[2][idx]/r[1][idx]))/10;
                                    
                                    alpha[idx] = fmin(
                                        alpha[idx],
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
 * Tag cells using wavelet sensor with the combination of three possible criteria:
 * 1. When ratio between wavelet coefficient and global maximum at any level is greater than the tolerance.
 * 2. When ratio between wavelet coefficient and local mean at any level is greater than the tolerance.
 * 3. When the Lipschitz's exponent is smaller than the tolerance.
 */
void
MultiresolutionTagger::tagCellsWithWaveletSensor(
    hier::Patch& patch,
    boost::shared_ptr<pdat::CellData<int> > tags,
    std::vector<boost::shared_ptr<pdat::CellData<double> > >& wavelet_coeffs,
    std::vector<double>& wavelet_coeffs_maxs,
    std::vector<boost::shared_ptr<pdat::CellData<double> > >& variable_local_means,
    boost::shared_ptr<pdat::CellData<double> > Lipschitz_exponent,
    double& global_tol,
    double& local_tol,
    double& alpha_tol,
    std::string& sensor_key)
{
    const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(patch_geom);
#endif
    
    // Get the dimensions of box that covers the interior of Patch.
    hier::Box dummy_box = patch.getBox();
    const hier::Box interior_box = dummy_box;
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of Patch plus
    // ghost cells.
    dummy_box.grow(d_num_ghosts);
    const hier::Box ghost_box = dummy_box;
    const hier::IntVector ghostcell_dims = ghost_box.numberCells();
    
    // Get the pointer of the tags.
    int* tag_ptr  = tags->getPointer(0);
    
    std::vector<boost::shared_ptr<pdat::CellData<double> > > wavelet_coeffs_local_mean;
    
    // Get the pointers to the wavelet coefficients.
    std::vector<double*> w;
    for (int li = 0; li < d_Harten_wavelet_num_level; li++)
    {
        w.push_back(wavelet_coeffs[li]->getPointer(0));
    }
    
    // Get the pointers to the variable local means.
    std::vector<double*> u_mean;
    if (d_Harten_wavelet_uses_local_tol)
    {
        for (int li = 0; li < d_Harten_wavelet_num_level; li++)
        {
            u_mean.push_back(variable_local_means[li]->getPointer(0));
        }
    }
    
    // Declare pointer to the Lipschitz's exponent.
    double* alpha = NULL;    
    
    if (sensor_key == "HARTEN_WAVELET")
    {
        if (d_Harten_wavelet_uses_alpha_tol)
        {
            alpha = Lipschitz_exponent->getPointer(0);
            
            computeLipschitzExponent(
                patch,
                wavelet_coeffs,
                Lipschitz_exponent,
                sensor_key);
        }
        
        if (d_dim == tbox::Dimension(1))
        {
            // NOT YET IMPLEMENTED
        }
        else if (d_dim == tbox::Dimension(2))
        {
            for (int j = 0; j < interior_dims[1]; j++)
            {
                for (int i = 0; i < interior_dims[0]; i++)
                {
                    // Compute indices.
                    const int idx = (i + d_num_ghosts[0]) +
                        (j + d_num_ghosts[1])*ghostcell_dims[0];
                    
                    const int idx_nghost = i + j*interior_dims[0];
                    
                    int tag_cell            = 1;
                    int tag_cell_global_tol = 0;
                    int tag_cell_local_tol  = 0;
                    int tag_cell_alpha_tol  = 0;
                    
                    if (d_Harten_wavelet_uses_global_tol)
                    {
                        /*
                        for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                        {
                            if (w[li][idx]/(wavelet_coeffs_maxs[li] + EPSILON) >= global_tol)
                            {
                                tag_cell_global_tol = 1;
                                
                                break;
                            }
                        }
                        */
                        
                        if (w[0][idx]/(wavelet_coeffs_maxs[0] + EPSILON) >= global_tol)
                        {
                            tag_cell_global_tol = 1;
                        }
                        
                        tag_cell &= tag_cell_global_tol;
                    }
                    
                    if (d_Harten_wavelet_uses_local_tol)
                    {
                        /*
                        for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                        {
                            if (w[li][idx]/(u_mean[li][idx] + EPSILON) >= local_tol)
                            {
                                tag_cell_local_tol = 1;
                                
                                break;
                            }
                        }
                        */
                        
                        if (w[0][idx]/(u_mean[0][idx] + EPSILON) >= local_tol)
                        {
                            tag_cell_local_tol = 1;
                        }
alpha[idx] = w[0][idx]/(u_mean[0][idx] + EPSILON);
                        tag_cell &= tag_cell_local_tol;
                    }
                    
                    if (d_Harten_wavelet_uses_alpha_tol)
                    {
                        if (alpha[idx] <= alpha_tol)
                        {
                            tag_cell_alpha_tol = 1;
                        }
                        
                        tag_cell &= tag_cell_alpha_tol;
                    }
                    
                    tag_ptr[idx_nghost] |= tag_cell;
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
                        const int idx = (i + d_num_ghosts[0]) +
                            (j + d_num_ghosts[1])*ghostcell_dims[0] +
                            (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                        
                        const int idx_nghost = i + j*interior_dims[0] + k*interior_dims[0]*interior_dims[1];
                        
                        int tag_cell            = 1;
                        int tag_cell_global_tol = 0;
                        int tag_cell_local_tol  = 0;
                        int tag_cell_alpha_tol  = 0;
                        
                        if (d_Harten_wavelet_uses_global_tol)
                        {
                            /*
                            for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                            {
                                if (w[li][idx]/(wavelet_coeffs_maxs[li] + EPSILON) >= global_tol)
                                {
                                    tag_cell_global_tol = 1;
                                    
                                    break;
                                }
                            }
                            */
                            
                            if (w[0][idx]/(wavelet_coeffs_maxs[0] + EPSILON) >= global_tol)
                            {
                                tag_cell_global_tol = 1;
                            }
                            
                            tag_cell &= tag_cell_global_tol;
                        }
                        
                        if (d_Harten_wavelet_uses_local_tol)
                        {
                            /*
                            for (int li = 0; li < d_Harten_wavelet_num_level; li++)
                            {
                                if (w[li][idx]/(u_mean[li][idx] + EPSILON) >= local_tol)
                                {
                                    tag_cell_local_tol = 1;
                                    
                                    break;
                                }
                            }
                            */
                            
                            if (w[0][idx]/(u_mean[0][idx] + EPSILON) >= local_tol)
                            {
                                tag_cell_local_tol = 1;
                            }
                            
                            tag_cell &= tag_cell_local_tol;
                        }
                        
                        if (d_Harten_wavelet_uses_alpha_tol)
                        {
                            if (alpha[idx] <= alpha_tol)
                            {
                                tag_cell_alpha_tol = 1;
                            }
                            
                            tag_cell &= tag_cell_alpha_tol;
                        }
                        
                        tag_ptr[idx_nghost] |= tag_cell;
                    }
                }
            }
        }
    }
}
