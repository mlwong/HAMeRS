#include "flow_model/refinement_taggers/gradient_tagger/GradientTagger.hpp"

#include "boost/lexical_cast.hpp"

#define PLOTTING_GRADIENT_TAGGER

#define EPSILON 1e-40

GradientTagger::GradientTagger(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
    const FLOW_MODEL& flow_model,
    const int& num_species,
    const boost::shared_ptr<EquationOfState>& equation_of_state,
    const boost::shared_ptr<tbox::Database>& gradient_tagger_db):
        d_object_name(object_name),
        d_dim(dim),
        d_grid_geometry(grid_geometry),
        d_num_ghosts(hier::IntVector::getZero(d_dim)),
        d_flow_model(flow_model),
        d_num_species(num_species),
        d_equation_of_state(equation_of_state),
        d_variables_set(false),
        d_num_ghosts_set(false)
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
}


/*
 * Register the temporary variables used in gradient tagger class.
 */
void
GradientTagger::registerGradientTaggerVariables(
    RungeKuttaLevelIntegrator* integrator)
{
    if (d_variables_set == true)
    {
        if (d_num_ghosts_set == true)
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
                        else if (variable_key == "ENSTROPHY")
                        {
                            d_Jameson_enstrophy_gradient =
                                boost::shared_ptr<pdat::CellVariable<double> > (
                                    new pdat::CellVariable<double>(d_dim, "Jameson enstrophy gradient", 1));
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
                                        d_Jameson_mass_fraction_gradient.push_back(
                                            boost::make_shared<pdat::CellVariable<double> >(
                                                d_dim,
                                                "Jameson mass fraction " + boost::lexical_cast<std::string>(di) +
                                                " gradient",
                                                1));
                                    }
                                    
                                    break;
                                }
                                case FIVE_EQN_ALLAIRE:
                                {
                                    for (int di = 0; di < d_num_species; di++)
                                    {
                                        d_Jameson_mass_fraction_gradient.push_back(
                                            boost::make_shared<pdat::CellVariable<double> >(
                                                d_dim,
                                                "Jameson mass fraction " + boost::lexical_cast<std::string>(di) +
                                                " gradient",
                                                1));
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
                                        d_Jameson_volume_fraction_gradient.push_back(
                                            boost::make_shared<pdat::CellVariable<double> >(
                                                d_dim,
                                                "Jameson volume fraction " + boost::lexical_cast<std::string>(di) +
                                                " gradient",
                                                1));
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
                    
                    for (int vi = 0; vi < static_cast<int>(d_Jameson_gradient_variables.size()); vi++)
                    {
                        // Get the key of the current variable.
                        std::string variable_key = d_Jameson_gradient_variables[vi];
                        
                        if (variable_key == "DENSITY")
                        {
                            integrator->registerVariable(
                                d_Jameson_density_gradient,
                                d_num_ghosts,
                                RungeKuttaLevelIntegrator::TIME_DEP,
                                    d_grid_geometry,
                                    "CONSERVATIVE_COARSEN",
                                    "CONSERVATIVE_LINEAR_REFINE");
                        }
                        else if (variable_key == "TOTAL_ENERGY")
                        {
                            integrator->registerVariable(
                                d_Jameson_total_energy_gradient,
                                d_num_ghosts,
                                RungeKuttaLevelIntegrator::TIME_DEP,
                                    d_grid_geometry,
                                    "CONSERVATIVE_COARSEN",
                                    "CONSERVATIVE_LINEAR_REFINE");
                        }
                        else if (variable_key == "PRESSURE")
                        {
                            integrator->registerVariable(
                                d_Jameson_pressure_gradient,
                                d_num_ghosts,
                                RungeKuttaLevelIntegrator::TIME_DEP,
                                    d_grid_geometry,
                                    "CONSERVATIVE_COARSEN",
                                    "CONSERVATIVE_LINEAR_REFINE");
                        }
                        else if (variable_key == "ENSTROPHY")
                        {
                            integrator->registerVariable(
                                d_Jameson_enstrophy_gradient,
                                d_num_ghosts,
                                RungeKuttaLevelIntegrator::TIME_DEP,
                                    d_grid_geometry,
                                    "CONSERVATIVE_COARSEN",
                                    "CONSERVATIVE_LINEAR_REFINE");
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
                                        integrator->registerVariable(
                                            d_Jameson_mass_fraction_gradient[di],
                                            d_num_ghosts,
                                            RungeKuttaLevelIntegrator::TIME_DEP,
                                            d_grid_geometry,
                                            "CONSERVATIVE_COARSEN",
                                            "CONSERVATIVE_LINEAR_REFINE");
                                    }
                                        
                                    break;
                                }
                                case FIVE_EQN_ALLAIRE:
                                {
                                    for (int di = 0; di < d_num_species; di++)
                                    {
                                        integrator->registerVariable(
                                            d_Jameson_mass_fraction_gradient[di],
                                            d_num_ghosts,
                                            RungeKuttaLevelIntegrator::TIME_DEP,
                                            d_grid_geometry,
                                            "CONSERVATIVE_COARSEN",
                                            "CONSERVATIVE_LINEAR_REFINE");
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
                                        integrator->registerVariable(
                                            d_Jameson_volume_fraction_gradient[di],
                                            d_num_ghosts,
                                            RungeKuttaLevelIntegrator::TIME_DEP,
                                            d_grid_geometry,
                                            "CONSERVATIVE_COARSEN",
                                            "CONSERVATIVE_LINEAR_REFINE");
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
GradientTagger::registerPlotQuantities(
    RungeKuttaLevelIntegrator* integrator,
    const boost::shared_ptr<appu::VisItDataWriter>& visit_writer,
    const boost::shared_ptr<hier::VariableContext>& plot_context)
{
#ifdef PLOTTING_GRADIENT_TAGGER
    if (d_variables_set == true)
    {
        if (d_num_ghosts_set == true)
        {
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
                                        visit_writer->registerPlotQuantity(
                                            "Jameson mass fraction " + boost::lexical_cast<std::string>(di) +
                                                " gradient",
                                            "SCALAR",
                                            vardb->mapVariableAndContextToIndex(
                                               d_Jameson_mass_fraction_gradient[di],
                                               plot_context));
                                    }
                                    
                                    break;
                                }
                                case FIVE_EQN_ALLAIRE:
                                {
                                    for (int di = 0; di < d_num_species; di++)
                                    {
                                        visit_writer->registerPlotQuantity(
                                            "Jameson mass fraction " + boost::lexical_cast<std::string>(di) +
                                                " gradient",
                                            "SCALAR",
                                            vardb->mapVariableAndContextToIndex(
                                               d_Jameson_mass_fraction_gradient[di],
                                               plot_context));
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
                                        visit_writer->registerPlotQuantity(
                                            "Jameson volume fraction " + boost::lexical_cast<std::string>(di) +
                                                " gradient",
                                            "SCALAR",
                                            vardb->mapVariableAndContextToIndex(
                                               d_Jameson_volume_fraction_gradient[di],
                                               plot_context));
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
 * Print all characteristics of the gradient tagger class.
 */
void
GradientTagger::printClassData(std::ostream& os) const
{
    os << "\nPrint GradientTagger object..."
       << std::endl;
    
    os << std::endl;
    
    os << "d_variables_set = "
       << d_variables_set
       << std::endl;
    os << "d_num_ghosts_set = "
       << d_num_ghosts_set
       << std::endl;
    
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
                            switch (d_flow_model)
                            {
                                case SINGLE_SPECIES: case FOUR_EQN_SHYUE:
                                {
                                    // Get the cell data of the density gradient.
                                    boost::shared_ptr<pdat::CellData<double> > gradient(
                                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                            patch.getPatchData(d_Jameson_density_gradient, data_context)));
                                    
                                    // Get the cell data of density.
                                    boost::shared_ptr<pdat::CellData<double> > density(
                                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                            patch.getPatchData(d_density, data_context)));
                                    
                                    // Compute the gradient.
                                    d_gradient_sensor_Jameson->computeGradient(patch, density, gradient);
                                    
                                    tagCellsWithGradientSensor(
                                        patch,
                                        tags,
                                        gradient,
                                        tol,
                                        sensor_key);
                                    
                                    break;
                                }
                                case FIVE_EQN_ALLAIRE:
                                {
                                    // Get the cell data of the density gradient.
                                    boost::shared_ptr<pdat::CellData<double> > gradient(
                                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                            patch.getPatchData(d_Jameson_density_gradient, data_context)));
                                    
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
                                    
                                    // Compute the gradient.
                                    d_gradient_sensor_Jameson->computeGradient(patch, density, gradient);
                                    
                                    tagCellsWithGradientSensor(
                                        patch,
                                        tags,
                                        gradient,
                                        tol,
                                        sensor_key);
                                    
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
                        else if (variable_key == "TOTAL_ENERGY")
                        {
                            switch (d_flow_model)
                            {
                                case SINGLE_SPECIES: case FOUR_EQN_SHYUE: case FIVE_EQN_ALLAIRE:
                                {
                                    // Get the cell data of the total energy gradient.
                                    boost::shared_ptr<pdat::CellData<double> > gradient(
                                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                            patch.getPatchData(d_Jameson_total_energy_gradient, data_context)));
                                    
                                    // Get the cell data of total energy.
                                    boost::shared_ptr<pdat::CellData<double> > total_energy(
                                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                            patch.getPatchData(d_total_energy, data_context)));
                                    
                                    // Compute the gradient.
                                    d_gradient_sensor_Jameson->computeGradient(patch, total_energy, gradient);
                                    
                                    tagCellsWithGradientSensor(
                                        patch,
                                        tags,
                                        gradient,
                                        tol,
                                        sensor_key);
                                    
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
                                    // Get the cell data of the pressure gradient.
                                    boost::shared_ptr<pdat::CellData<double> > gradient(
                                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                            patch.getPatchData(d_Jameson_pressure_gradient, data_context)));
                                    
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
                                    
                                    // Compute the gradient.
                                    d_gradient_sensor_Jameson->computeGradient(patch, pressure, gradient);
                                    
                                    tagCellsWithGradientSensor(
                                        patch,
                                        tags,
                                        gradient,
                                        tol,
                                        sensor_key);
                                    
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
                                    // Get the cell data of the enstrophy gradient.
                                    boost::shared_ptr<pdat::CellData<double> > gradient(
                                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                            patch.getPatchData(d_Jameson_enstrophy_gradient, data_context)));
                                    
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
                                    
                                    // Compute the enstrophy and gradients.
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
                                    
                                    // Compute the gradient.
                                    d_gradient_sensor_Jameson->computeGradient(patch, enstrophy, gradient);
                                    
                                    tagCellsWithGradientSensor(
                                        patch,
                                        tags,
                                        gradient,
                                        tol,
                                        sensor_key);
                                    
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
                                    for (int di = 0; di < d_mass_fraction->getDepth(); di++)
                                    {
                                        // Get the cell data of the mass fraction gradient.
                                        boost::shared_ptr<pdat::CellData<double> > gradient(
                                            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                                patch.getPatchData(d_Jameson_mass_fraction_gradient[di], data_context)));
                                        
                                        // Get the cell data of mass fraction.
                                        boost::shared_ptr<pdat::CellData<double> > mass_fraction(
                                            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                                patch.getPatchData(d_mass_fraction, data_context)));
                                        
                                        // Compute the gradient.
                                        d_gradient_sensor_Jameson->computeGradient(patch, mass_fraction, gradient, di);
                                        
                                        tagCellsWithGradientSensor(
                                            patch,
                                            tags,
                                            gradient,
                                            tol,
                                            sensor_key);
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
                                        // Get the cell data of the mass fraction gradient.
                                        boost::shared_ptr<pdat::CellData<double> > gradient(
                                            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                                patch.getPatchData(d_Jameson_mass_fraction_gradient[di], data_context)));
                                        
                                        // Compute the gradient.
                                        d_gradient_sensor_Jameson->computeGradient(patch, mass_fraction, gradient, di);
                                        
                                        tagCellsWithGradientSensor(
                                            patch,
                                            tags,
                                            gradient,
                                            tol,
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
                                        // Get the cell data of the volume fraction gradient.
                                        boost::shared_ptr<pdat::CellData<double> > gradient(
                                            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                                patch.getPatchData(d_Jameson_volume_fraction_gradient[di], data_context)));
                                        
                                        // Get the cell data of volume fraction.
                                        boost::shared_ptr<pdat::CellData<double> > volume_fraction(
                                            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                                patch.getPatchData(d_volume_fraction, data_context)));
                                        
                                        // Compute the gradient.
                                        d_gradient_sensor_Jameson->computeGradient(patch, volume_fraction, gradient, di);
                                        
                                        tagCellsWithGradientSensor(
                                            patch,
                                            tags,
                                            gradient,
                                            tol,
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
                    } // Looop over variables chosen.
                }
            } // Loop over gradient sensors chosen.
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
    
    // Get the pointer to the gradient.
    double* psi = gradient->getPointer(0);
    
    if (sensor_key == "JAMESON_GRADIENT")
    {
        if (d_dim == tbox::Dimension(1))
        {
            for (int i = 0; i < interior_dims[0]; i++)
            {
                // Compute indices.
                const int idx = i + d_num_ghosts[0];
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
                    const int idx = (i + d_num_ghosts[0]) +
                        (j + d_num_ghosts[1])*ghostcell_dims[0];
                    
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
                        const int idx = (i + d_num_ghosts[0]) +
                            (j + d_num_ghosts[1])*ghostcell_dims[0] +
                            (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                        
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
