#include "flow/refinement_taggers/ValueTagger.hpp"

#include <algorithm>
#include "boost/lexical_cast.hpp"

// #define PLOTTING_VALUE_TAGGER

#define EPSILON 1e-40

ValueTagger::ValueTagger(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
    const int& num_species,
    const boost::shared_ptr<FlowModel>& flow_model,
    const boost::shared_ptr<tbox::Database>& value_tagger_db):
        d_object_name(object_name),
        d_dim(dim),
        d_grid_geometry(grid_geometry),
        d_num_value_ghosts(hier::IntVector::getZero(d_dim)),
        d_num_species(num_species),
        d_flow_model(flow_model),
        d_value_tagger_density_max(0.0),
        d_value_tagger_total_energy_max(0.0),
        d_value_tagger_pressure_max(0.0),
        d_value_tagger_dilatation_max(0.0),
        d_value_tagger_enstrophy_max(0.0)
{
    if (value_tagger_db != nullptr)
    {
        if (value_tagger_db->keyExists("variables"))
        {
            d_variables = value_tagger_db->getStringVector("variables");
        }
        else if (value_tagger_db->keyExists("d_variables"))
        {
            d_variables = value_tagger_db->getStringVector("d_variables");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "No key 'variables'/'d_variables' found in database of value tagger."
                << std::endl);
        }
        
        // Loop over variables chosen.
        for (int vi = 0; vi < static_cast<int>(d_variables.size()); vi++)
        {
            std::string variable_key = d_variables[vi];
            
            if (!((variable_key == "DENSITY") ||
                  (variable_key == "TOTAL_ENERGY") ||
                  (variable_key == "PRESSURE") ||
                  (variable_key == "DILATATION") ||
                  (variable_key == "ENSTROPHY") ||
                  (variable_key == "MASS_FRACTION")))
            {
                TBOX_ERROR(d_object_name
                    << ":"
                    << "Unknown/unsupported variable: "
                    << variable_key
                    << "\nin input."
                    << std::endl);
            }
        }
        
        int num_true = 0;
        
        /*
         * Get the settings for global upper tolerances.
         */
        
        if (value_tagger_db->keyExists("uses_global_tol_up"))
        {
            d_uses_global_tol_up = value_tagger_db->getBoolVector("uses_global_tol_up");
        }
        else if (value_tagger_db->keyExists("d_uses_global_tol_up"))
        {
            d_uses_global_tol_up = value_tagger_db->getBoolVector("d_uses_global_tol_up");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "No key 'uses_global_tol_up'/'d_uses_global_tol_up"
                << " found in database of value tagger."
                << std::endl);
        }
        
        if (static_cast<int>(d_variables.size()) != static_cast<int>(d_uses_global_tol_up.size()))
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "The numbers of variables and switches for global upper tolerances provided don't match"
                << " in database of value tagger."
                << std::endl);
        }
        
        num_true = std::count(d_uses_global_tol_up.begin(), d_uses_global_tol_up.end(), true);
        
        if (num_true > 0)
        {
            if (value_tagger_db->keyExists("global_tol_up"))
            {
                d_global_tol_up = value_tagger_db->getDoubleVector("global_tol_up");
            }
            else if (value_tagger_db->keyExists("d_global_tol_up"))
            {
                d_global_tol_up = value_tagger_db->getDoubleVector("d_global_tol_up");
            }
            else
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "No key 'global_tol_up'/'d_global_tol_up' found in database of value tagger."
                    << std::endl);
            }
            
            if (num_true != static_cast<int>(d_global_tol_up.size()))
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "The number of variables that use global upper tolerances and number of"
                    << " global upper tolerances provided don't match in database of value tagger."
                    << std::endl);
            }
        }
        
        /*
         * Get the settings for global lower tolerances.
         */
        
        if (value_tagger_db->keyExists("uses_global_tol_lo"))
        {
            d_uses_global_tol_lo = value_tagger_db->getBoolVector("uses_global_tol_lo");
        }
        else if (value_tagger_db->keyExists("d_uses_global_tol_lo"))
        {
            d_uses_global_tol_lo = value_tagger_db->getBoolVector("d_uses_global_tol_lo");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "No key 'uses_global_tol_lo'/'d_uses_global_tol_lo"
                << " found in database of value tagger."
                << std::endl);
        }
        
        if (static_cast<int>(d_variables.size()) != static_cast<int>(d_uses_global_tol_lo.size()))
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "The numbers of variables and switches for global lower tolerances provided don't match"
                << " in database of value tagger."
                << std::endl);
        }
        
        num_true = std::count(d_uses_global_tol_lo.begin(), d_uses_global_tol_lo.end(), true);
        
        if (num_true > 0)
        {
            if (value_tagger_db->keyExists("global_tol_lo"))
            {
                d_global_tol_lo = value_tagger_db->getDoubleVector("global_tol_lo");
            }
            else if (value_tagger_db->keyExists("d_global_tol_lo"))
            {
                d_global_tol_lo = value_tagger_db->getDoubleVector("d_global_tol_lo");
            }
            else
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "No key 'global_tol_lo'/'d_global_tol_lo' found in database of value tagger."
                    << std::endl);
            }
            
            if (num_true != static_cast<int>(d_global_tol_lo.size()))
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "The number of variables that use global lower tolerances and number of"
                    << " global lower tolerances provided don't match in database of value tagger."
                    << std::endl);
            }
        }
        
        /*
         * Get the settings for local upper tolerances.
         */
        
        if (value_tagger_db->keyExists("uses_local_tol_up"))
        {
            d_uses_local_tol_up = value_tagger_db->getBoolVector("uses_local_tol_up");
        }
        else if (value_tagger_db->keyExists("d_uses_local_tol_up"))
        {
            d_uses_local_tol_up = value_tagger_db->getBoolVector("d_uses_local_tol_up");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "No key 'uses_local_tol_up'/'d_uses_local_tol_up"
                << " found in database of value tagger."
                << std::endl);
        }
        
        if (static_cast<int>(d_variables.size()) != static_cast<int>(d_uses_local_tol_up.size()))
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "The numbers of variables and switches for local upper tolerances provided don't match"
                << " in database of value tagger."
                << std::endl);
        }
        
        num_true = std::count(d_uses_local_tol_up.begin(), d_uses_local_tol_up.end(), true);
        
        if (num_true > 0)
        {
            if (value_tagger_db->keyExists("local_tol_up"))
            {
                d_local_tol_up = value_tagger_db->getDoubleVector("local_tol_up");
            }
            else if (value_tagger_db->keyExists("d_local_tol_up"))
            {
                d_local_tol_up = value_tagger_db->getDoubleVector("d_local_tol_up");
            }
            else
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "No key 'local_tol_up'/'d_local_tol_up' found in database of value tagger."
                    << std::endl);
            }
            
            if (num_true != static_cast<int>(d_local_tol_up.size()))
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "The number of variables that use local upper tolerances and number of"
                    << " local upper tolerances provided don't match in database of value tagger."
                    << std::endl);
            }
        }
        
        /*
         * Get the settings for local lower tolerances.
         */
        
        if (value_tagger_db->keyExists("uses_local_tol_lo"))
        {
            d_uses_local_tol_lo = value_tagger_db->getBoolVector("uses_local_tol_lo");
        }
        else if (value_tagger_db->keyExists("d_uses_local_tol_lo"))
        {
            d_uses_local_tol_lo = value_tagger_db->getBoolVector("d_uses_local_tol_lo");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "No key 'uses_local_tol_lo'/'d_uses_local_tol_lo"
                << " found in database of value tagger."
                << std::endl);
        }
        
        if (static_cast<int>(d_variables.size()) != static_cast<int>(d_uses_local_tol_lo.size()))
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "The numbers of variables and switches for local lower tolerances provided don't match"
                << " in database of value tagger."
                << std::endl);
        }
        
        num_true = std::count(d_uses_local_tol_lo.begin(), d_uses_local_tol_lo.end(), true);
        
        if (num_true > 0)
        {
            if (value_tagger_db->keyExists("local_tol_lo"))
            {
                d_local_tol_lo = value_tagger_db->getDoubleVector("local_tol_lo");
            }
            else if (value_tagger_db->keyExists("d_local_tol_lo"))
            {
                d_local_tol_lo = value_tagger_db->getDoubleVector("d_local_tol_lo");
            }
            else
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "No key 'local_tol_lo'/'d_local_tol_lo' found in database of value tagger."
                    << std::endl);
            }
            
            if (num_true != static_cast<int>(d_local_tol_lo.size()))
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "The number of variables that use local lower tolerances and number of"
                    << " local lower tolerances provided don't match in database of value tagger."
                    << std::endl);
            }
        }
    }
    else
    {
        TBOX_WARNING(d_object_name
            << ": "
            << "Key data 'Value_tagger' not found in input/restart database."
            << " No refinement with value tagger will occur."
            << std::endl);
    }
}


/*
 * Register the temporary variables used in value tagger class.
 */
void
ValueTagger::registerValueTaggerVariables(
    RungeKuttaLevelIntegrator* integrator)
{
    // Loop over variables chosen.
    for (int vi = 0; vi < static_cast<int>(d_variables.size()); vi++)
    {
        // Get the key of the current variable.
        std::string variable_key = d_variables[vi];
        
        if (variable_key == "DENSITY")
        {
            d_value_tagger_variable_density = boost::make_shared<pdat::CellVariable<double> >(
                d_dim,
                "Value tagger density",
                1);
        }
        else if (variable_key == "TOTAL_ENERGY")
        {
            d_value_tagger_variable_total_energy = boost::make_shared<pdat::CellVariable<double> >(
                d_dim,
                "Value tagger total energy",
                1);
        }
        else if (variable_key == "PRESSURE")
        {
            d_value_tagger_variable_pressure = boost::make_shared<pdat::CellVariable<double> >(
                d_dim,
                "Value tagger pressure",
                1);
        }
        else if (variable_key == "DILATATION")
        {
            d_value_tagger_variable_dilatation = boost::make_shared<pdat::CellVariable<double> >(
                d_dim,
                "Value tagger dilatation",
                1);
        }
        else if (variable_key == "ENSTROPHY")
        {
            d_value_tagger_variable_enstrophy = boost::make_shared<pdat::CellVariable<double> >(
                d_dim,
                "Value tagger enstrophy",
                1);
        }
        else if (variable_key == "MASS_FRACTION")
        {
            d_value_tagger_variable_mass_fraction.reserve(d_num_species);
            
            for (int si = 0; si < d_num_species; si++)
            {
                d_value_tagger_variable_mass_fraction.push_back(boost::make_shared<pdat::CellVariable<double> >(
                    d_dim,
                    "Value tagger mass fraction "  + boost::lexical_cast<std::string>(si),
                    1));
            }
            
            if (d_uses_global_tol_up[vi] || d_uses_global_tol_lo[vi])
            {
                d_value_tagger_mass_fraction_max.reserve(d_num_species);
                
                for (int si = 0; si < d_num_species; si++)
                {
                    d_value_tagger_mass_fraction_max.push_back(0.0);
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
    
    // Loop over variables chosen.
    for (int vi = 0; vi < static_cast<int>(d_variables.size()); vi++)
    {
        // Get the key of the current variable.
        std::string variable_key = d_variables[vi];
        
        if (variable_key == "DENSITY")
        {
            integrator->registerVariable(
                d_value_tagger_variable_density,
                d_num_value_ghosts,
                RungeKuttaLevelIntegrator::TEMPORARY,
                    d_grid_geometry,
                    "NO_COARSEN",
                    "NO_REFINE");
        }
        else if (variable_key == "TOTAL_ENERGY")
        {
            integrator->registerVariable(
                d_value_tagger_variable_total_energy,
                d_num_value_ghosts,
                RungeKuttaLevelIntegrator::TEMPORARY,
                    d_grid_geometry,
                    "NO_COARSEN",
                    "NO_REFINE");
        }
        else if (variable_key == "PRESSURE")
        {
            integrator->registerVariable(
                d_value_tagger_variable_pressure,
                d_num_value_ghosts,
                RungeKuttaLevelIntegrator::TEMPORARY,
                    d_grid_geometry,
                    "NO_COARSEN",
                    "NO_REFINE");
        }
        else if (variable_key == "DILATATION")
        {
            integrator->registerVariable(
                d_value_tagger_variable_dilatation,
                d_num_value_ghosts,
                RungeKuttaLevelIntegrator::TEMPORARY,
                    d_grid_geometry,
                    "NO_COARSEN",
                    "NO_REFINE");
        }
        else if (variable_key == "ENSTROPHY")
        {
            integrator->registerVariable(
                d_value_tagger_variable_enstrophy,
                d_num_value_ghosts,
                RungeKuttaLevelIntegrator::TEMPORARY,
                    d_grid_geometry,
                    "NO_COARSEN",
                    "NO_REFINE");
        }
        else if (variable_key == "MASS_FRACTION")
        {
            for (int si = 0; si < d_num_species; si++)
            {
                integrator->registerVariable(
                    d_value_tagger_variable_mass_fraction[si],
                    d_num_value_ghosts,
                    RungeKuttaLevelIntegrator::TEMPORARY,
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


/*
 * Register the plotting quantities.
 */
void
ValueTagger::registerPlotQuantities(
    const boost::shared_ptr<appu::VisItDataWriter>& visit_writer,
    const boost::shared_ptr<hier::VariableContext>& plot_context)
{
#ifdef PLOTTING_VALUE_TAGGER
#endif
}


/*
 * Print all characteristics of the value tagger class.
 */
void
ValueTagger::printClassData(std::ostream& os) const
{
    os << "\nPrint ValueTagger object..."
       << std::endl;
    
    os << std::endl;
    
    os << "d_variables = ";
    for (int vi = 0; vi < static_cast<int>(d_variables.size()); vi++)
    {
        os << "\"" << d_variables[vi] << "\"";
        
        if (vi < static_cast<int>(d_variables.size()) - 1)
        {
            os << ", ";
        }
    }
    os << std::endl;
    
    os << "d_uses_global_tol_up = ";
    for (int ti = 0; ti < static_cast<int>(d_uses_global_tol_up.size()); ti++)
    {
        os << "\"" << d_uses_global_tol_up[ti] << "\"";
        
        if (ti < static_cast<int>(d_uses_global_tol_up.size()) - 1)
        {
            os << ", ";
        }
    }
    os << std::endl;
    
    os << "d_uses_global_tol_lo = ";
    for (int ti = 0; ti < static_cast<int>(d_uses_global_tol_lo.size()); ti++)
    {
        os << "\"" << d_uses_global_tol_lo[ti] << "\"";
        
        if (ti < static_cast<int>(d_uses_global_tol_lo.size()) - 1)
        {
            os << ", ";
        }
    }
    os << std::endl;
    
    os << "d_global_tol_up = ";
    for (int ti = 0; ti < static_cast<int>(d_global_tol_up.size()); ti++)
    {
        os << "\"" << d_global_tol_up[ti] << "\"";
        
        if (ti < static_cast<int>(d_global_tol_up.size()) - 1)
        {
            os << ", ";
        }
    }
    os << std::endl;
    
    os << "d_global_tol_lo = ";
    for (int ti = 0; ti < static_cast<int>(d_global_tol_lo.size()); ti++)
    {
        os << "\"" << d_global_tol_lo[ti] << "\"";
        
        if (ti < static_cast<int>(d_global_tol_lo.size()) - 1)
        {
            os << ", ";
        }
    }
    os << std::endl;
    
    os << "d_uses_local_tol_up = ";
    for (int ti = 0; ti < static_cast<int>(d_uses_local_tol_up.size()); ti++)
    {
        os << "\"" << d_uses_local_tol_up[ti] << "\"";
        
        if (ti < static_cast<int>(d_uses_local_tol_up.size()) - 1)
        {
            os << ", ";
        }
    }
    os << std::endl;
    
    os << "d_uses_local_tol_lo = ";
    for (int ti = 0; ti < static_cast<int>(d_uses_local_tol_lo.size()); ti++)
    {
        os << "\"" << d_uses_local_tol_lo[ti] << "\"";
        
        if (ti < static_cast<int>(d_uses_local_tol_lo.size()) - 1)
        {
            os << ", ";
        }
    }
    os << std::endl;
    
    os << "d_local_tol_up = ";
    for (int ti = 0; ti < static_cast<int>(d_local_tol_up.size()); ti++)
    {
        os << "\"" << d_local_tol_up[ti] << "\"";
        
        if (ti < static_cast<int>(d_local_tol_up.size()) - 1)
        {
            os << ", ";
        }
    }
    os << std::endl;
    
    os << "d_local_tol_lo = ";
    for (int ti = 0; ti < static_cast<int>(d_local_tol_lo.size()); ti++)
    {
        os << "\"" << d_local_tol_lo[ti] << "\"";
        
        if (ti < static_cast<int>(d_local_tol_lo.size()) - 1)
        {
            os << ", ";
        }
    }
    os << std::endl;
}


/*
 * Put the characteristics of the value tagger class into the restart
 * database.
 */
void
ValueTagger::putToRestart(
    const boost::shared_ptr<tbox::Database>& restart_db) const
{
    restart_db->putStringVector("d_variables",
        d_variables);
    
    restart_db->putBoolVector("d_uses_global_tol_up", d_uses_global_tol_up);
    restart_db->putBoolVector("d_uses_global_tol_lo", d_uses_global_tol_lo);
    restart_db->putBoolVector("d_uses_local_tol_up", d_uses_local_tol_up);
    restart_db->putBoolVector("d_uses_local_tol_lo", d_uses_local_tol_lo);
    
    int num_true = 0;
    
    num_true = std::count(d_uses_global_tol_up.begin(), d_uses_global_tol_up.end(), true);
    if (num_true > 0)
    {
        restart_db->putDoubleVector("d_global_tol_up", d_global_tol_up);
    }
    
    num_true = std::count(d_uses_global_tol_lo.begin(), d_uses_global_tol_lo.end(), true);
    if (num_true > 0)
    {
        restart_db->putDoubleVector("d_global_tol_lo", d_global_tol_lo);
    }
    
    num_true = std::count(d_uses_local_tol_up.begin(), d_uses_local_tol_up.end(), true);
    if (num_true > 0)
    {
        restart_db->putDoubleVector("d_local_tol_up", d_local_tol_up);
    }
    
    num_true = std::count(d_uses_local_tol_lo.begin(), d_uses_local_tol_lo.end(), true);
    if (num_true > 0)
    {
        restart_db->putDoubleVector("d_local_tol_lo", d_local_tol_lo);
    }
}


/*
 * Compute values for value tagger.
 */
void
ValueTagger::computeValueTaggerValues(
    hier::Patch& patch,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
    // Loop over variables chosen.
    for (int vi = 0; vi < static_cast<int>(d_variables.size()); vi++)
    {
        // Get the key of the current variable.
        std::string variable_key = d_variables[vi];
        
        if (variable_key == "DENSITY")
        {
            /*
             * Register the patch and density in the flow model and compute the corresponding cell data.
             */
            
            d_flow_model->registerPatchWithDataContext(patch, data_context);
            
            std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
            
            num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("DENSITY", d_num_value_ghosts));
            
            d_flow_model->registerDerivedCellVariable(num_subghosts_of_data);
            
            d_flow_model->computeGlobalDerivedCellData();
            
            /*
             * Get the pointer to density data inside the flow model.
             */
            
            boost::shared_ptr<pdat::CellData<double> > flow_model_data_density =
                d_flow_model->getGlobalCellData("DENSITY");
            
            /*
             * Transfer data from flow model to the class variable.
             */
            
            transferDataToClassVariable(
                patch,
                data_context,
                flow_model_data_density,
                d_value_tagger_variable_density,
                0);
            
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
            
            num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("TOTAL_ENERGY", d_num_value_ghosts));
            
            d_flow_model->registerDerivedCellVariable(num_subghosts_of_data);
            
            d_flow_model->computeGlobalDerivedCellData();
            
            /*
             * Get the pointer to total energy data inside the flow model.
             */
            
            boost::shared_ptr<pdat::CellData<double> > flow_model_data_total_energy =
                d_flow_model->getGlobalCellData("TOTAL_ENERGY");
            
            /*
             * Transfer data from flow model to the class variable.
             */
            
            transferDataToClassVariable(
                patch,
                data_context,
                flow_model_data_total_energy,
                d_value_tagger_variable_total_energy,
                0);
            
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
            
            num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("PRESSURE", d_num_value_ghosts));
            
            d_flow_model->registerDerivedCellVariable(num_subghosts_of_data);
            
            d_flow_model->computeGlobalDerivedCellData();
            
            /*
             * Get the pointer to pressure data inside the flow model.
             */
            
            boost::shared_ptr<pdat::CellData<double> > flow_model_data_pressure =
                d_flow_model->getGlobalCellData("PRESSURE");
            
            /*
             * Transfer data from flow model to the class variable.
             */
            
            transferDataToClassVariable(
                patch,
                data_context,
                flow_model_data_pressure,
                d_value_tagger_variable_pressure,
                0);
            
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
            
            num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("DILATATION", d_num_value_ghosts));
            
            d_flow_model->registerDerivedCellVariable(num_subghosts_of_data);
            
            d_flow_model->computeGlobalDerivedCellData();
            
            /*
             * Get the pointer to dilatation data inside the flow model.
             */
            
            boost::shared_ptr<pdat::CellData<double> > flow_model_data_dilatation =
                d_flow_model->getGlobalCellData("DILATATION");
            
            /*
             * Transfer data from flow model to the class variable.
             */
            
            transferDataToClassVariable(
                patch,
                data_context,
                flow_model_data_dilatation,
                d_value_tagger_variable_dilatation,
                0);
            
            /*
             * Unregister the patch and data of all registered derived cell variables in the flow model.
             */
            
            d_flow_model->unregisterPatch();
        }
        else if (variable_key == "ENSTROPHY")
        {
            /*
             * Register the patch and dilatation in the flow model and compute the corresponding cell data.
             */
            
            d_flow_model->registerPatchWithDataContext(patch, data_context);
            
            std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
            
            num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("ENSTROPHY", d_num_value_ghosts));
            
            d_flow_model->registerDerivedCellVariable(num_subghosts_of_data);
            
            d_flow_model->computeGlobalDerivedCellData();
            
            /*
             * Get the pointer to enstrophy data inside the flow model.
             */
            
            boost::shared_ptr<pdat::CellData<double> > flow_model_data_enstrophy =
                d_flow_model->getGlobalCellData("ENSTROPHY");
            
            /*
             * Transfer data from flow model to the class variable.
             */
            
            transferDataToClassVariable(
                patch,
                data_context,
                flow_model_data_enstrophy,
                d_value_tagger_variable_enstrophy,
                0);
            
            /*
             * Unregister the patch and data of all registered derived cell variables in the flow model.
             */
            
            d_flow_model->unregisterPatch();
        }
        else if (variable_key == "MASS_FRACTION")
        {
            /*
             * Register the patch and mass fraction in the flow model and compute the corresponding cell data.
             */
            
            d_flow_model->registerPatchWithDataContext(patch, data_context);
            
            std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
            
            num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("MASS_FRACTION", d_num_value_ghosts));
            
            d_flow_model->registerDerivedCellVariable(num_subghosts_of_data);
            
            d_flow_model->computeGlobalDerivedCellData();
            
            /*
             * Get the pointer to mass fraction data inside the flow model.
             */
            
            boost::shared_ptr<pdat::CellData<double> > flow_model_data_mass_fraction =
                d_flow_model->getGlobalCellData("MASS_FRACTION");
            
            /*
             * Transfer data from flow model to the class variable.
             */
            
            for (int si = 0; si < d_num_species; si++)
            {
                transferDataToClassVariable(
                    patch,
                    data_context,
                    flow_model_data_mass_fraction,
                    d_value_tagger_variable_mass_fraction[si],
                    si);
            }
            
            /*
             * Unregister the patch and data of all registered derived cell variables in the flow model.
             */
            
            d_flow_model->unregisterPatch();
        }
    }
}


/*
 * Get the statistics of values that are required by the value tagger.
 */
void
ValueTagger::getValueStatistics(
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const int level_number,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    math::HierarchyCellDataOpsReal<double> cell_double_operator(patch_hierarchy, level_number, level_number);
    
    hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();
    
    // Loop over variables chosen.
    for (int vi = 0; vi < static_cast<int>(d_variables.size()); vi++)
    {
        if (d_uses_global_tol_up[vi] || d_uses_global_tol_lo[vi])
        {
            // Get the key of the current variable.
            std::string variable_key = d_variables[vi];
            
            if (variable_key == "DENSITY")
            {
                const int rho_id = variable_db->mapVariableAndContextToIndex(
                    d_value_tagger_variable_density,
                    data_context);
                
                double rho_max_local = cell_double_operator.max(rho_id);
                d_value_tagger_density_max = 0.0;
                
                mpi.Allreduce(
                    &rho_max_local,
                    &d_value_tagger_density_max,
                    1,
                    MPI_DOUBLE,
                    MPI_MAX);
            }
            else if (variable_key == "TOTAL_ENERGY")
            {
                const int E_id = variable_db->mapVariableAndContextToIndex(
                    d_value_tagger_variable_total_energy,
                    data_context);
                
                double E_max_local = cell_double_operator.max(E_id);
                d_value_tagger_total_energy_max = 0.0;
                
                mpi.Allreduce(
                    &E_max_local,
                    &d_value_tagger_total_energy_max,
                    1,
                    MPI_DOUBLE,
                    MPI_MAX);
            }
            else if (variable_key == "PRESSURE")
            {
                const int p_id = variable_db->mapVariableAndContextToIndex(
                    d_value_tagger_variable_pressure,
                    data_context);
                
                double p_max_local = cell_double_operator.max(p_id);
                d_value_tagger_pressure_max = 0.0;
                
                mpi.Allreduce(
                    &p_max_local,
                    &d_value_tagger_pressure_max,
                    1,
                    MPI_DOUBLE,
                    MPI_MAX);
            }
            else if (variable_key == "DILATATION")
            {
                const int theta_id = variable_db->mapVariableAndContextToIndex(
                    d_value_tagger_variable_dilatation,
                    data_context);
                
                double theta_max_local = cell_double_operator.max(theta_id);
                d_value_tagger_dilatation_max = 0.0;
                
                mpi.Allreduce(
                    &theta_max_local,
                    &d_value_tagger_dilatation_max,
                    1,
                    MPI_DOUBLE,
                    MPI_MAX);
            }
            else if (variable_key == "ENSTROPHY")
            {
                const int Omega_id = variable_db->mapVariableAndContextToIndex(
                    d_value_tagger_variable_enstrophy,
                    data_context);
                
                double Omega_max_local = cell_double_operator.max(Omega_id);
                d_value_tagger_enstrophy_max = 0.0;
                
                mpi.Allreduce(
                    &Omega_max_local,
                    &d_value_tagger_enstrophy_max,
                    1,
                    MPI_DOUBLE,
                    MPI_MAX);
            }
            else if (variable_key == "MASS_FRACTION")
            {
                for (int si = 0; si < d_num_species; si++)
                {
                    const int Y_id = variable_db->mapVariableAndContextToIndex(
                        d_value_tagger_variable_mass_fraction[si],
                        data_context);
                    
                    double Y_max_local = cell_double_operator.max(Y_id);
                    d_value_tagger_mass_fraction_max[si] = 0.0;
                    
                    mpi.Allreduce(
                        &Y_max_local,
                        &d_value_tagger_mass_fraction_max[si],
                        1,
                        MPI_DOUBLE,
                        MPI_MAX);
                }
            }
        }
    }
}


/*
 * Tag cells for refinement using value tagger.
 */
void
ValueTagger::tagCells(
   hier::Patch& patch,
   boost::shared_ptr<pdat::CellData<int> > tags,
   const boost::shared_ptr<hier::VariableContext>& data_context)
{
    int count_global_tol_up = 0;
    int count_global_tol_lo = 0;
    int count_local_tol_up = 0;
    int count_local_tol_lo = 0;
    
    // Loop over variables chosen.
    for (int vi = 0; vi < static_cast<int>(d_variables.size()); vi++)
    {
        // Get the key of the current variable.
        std::string variable_key = d_variables[vi];
        
        bool uses_global_tol_up = d_uses_global_tol_up[vi];
        bool uses_global_tol_lo = d_uses_global_tol_lo[vi];
        bool uses_local_tol_up = d_uses_local_tol_up[vi];
        bool uses_local_tol_lo = d_uses_local_tol_lo[vi];
        
        double global_tol_up = 0.0;
        double global_tol_lo = 0.0;
        double local_tol_up = 0.0;
        double local_tol_lo = 0.0;
        
        if (uses_global_tol_up)
        {
            global_tol_up = d_global_tol_up[count_global_tol_up];
            count_global_tol_up++;
        }
        
        if (uses_global_tol_lo)
        {
            global_tol_lo = d_global_tol_lo[count_global_tol_lo];
            count_global_tol_lo++;
        }
        
        if (uses_local_tol_up)
        {
            local_tol_up = d_local_tol_up[count_local_tol_up];
            count_local_tol_up++;
        }
        
        if (uses_local_tol_lo)
        {
            local_tol_lo = d_local_tol_lo[count_local_tol_lo];
            count_local_tol_lo++;
        }
        
        if (variable_key == "DENSITY")
        {
            tagCellsWithValue(patch,
                data_context,
                tags,
                d_value_tagger_variable_density,
                d_value_tagger_density_max,
                uses_global_tol_up,
                uses_global_tol_lo,
                uses_local_tol_up,
                uses_local_tol_lo,
                global_tol_up,
                global_tol_lo,
                local_tol_up,
                local_tol_lo);
        }
        else if (variable_key == "TOTAL_ENERGY")
        {
            tagCellsWithValue(patch,
                data_context,
                tags,
                d_value_tagger_variable_total_energy,
                d_value_tagger_total_energy_max,
                uses_global_tol_up,
                uses_global_tol_lo,
                uses_local_tol_up,
                uses_local_tol_lo,
                global_tol_up,
                global_tol_lo,
                local_tol_up,
                local_tol_lo);
        }
        else if (variable_key == "PRESSURE")
        {
            tagCellsWithValue(patch,
                data_context,
                tags,
                d_value_tagger_variable_pressure,
                d_value_tagger_pressure_max,
                uses_global_tol_up,
                uses_global_tol_lo,
                uses_local_tol_up,
                uses_local_tol_lo,
                global_tol_up,
                global_tol_lo,
                local_tol_up,
                local_tol_lo);
        }
        else if (variable_key == "DILATATION")
        {
            tagCellsWithValue(patch,
                data_context,
                tags,
                d_value_tagger_variable_dilatation,
                d_value_tagger_dilatation_max,
                uses_global_tol_up,
                uses_global_tol_lo,
                uses_local_tol_up,
                uses_local_tol_lo,
                global_tol_up,
                global_tol_lo,
                local_tol_up,
                local_tol_lo);
        }
        else if (variable_key == "ENSTROPHY")
        {
            tagCellsWithValue(patch,
                data_context,
                tags,
                d_value_tagger_variable_enstrophy,
                d_value_tagger_enstrophy_max,
                uses_global_tol_up,
                uses_global_tol_lo,
                uses_local_tol_up,
                uses_local_tol_lo,
                global_tol_up,
                global_tol_lo,
                local_tol_up,
                local_tol_lo);
        }
        else if (variable_key == "MASS_FRACTION")
        {
            for (int si = 0; si < d_num_species; si++)
            {
                tagCellsWithValue(patch,
                    data_context,
                    tags,
                    d_value_tagger_variable_mass_fraction[si],
                    d_value_tagger_mass_fraction_max[si],
                    uses_global_tol_up,
                    uses_global_tol_lo,
                    uses_local_tol_up,
                    uses_local_tol_lo,
                    global_tol_up,
                    global_tol_lo,
                    local_tol_up,
                    local_tol_lo);
            }
        }
    }
}


/*
 * Tag cells for refinement using data values.
 */
void
ValueTagger::tagCellsWithValue(
    hier::Patch& patch,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const boost::shared_ptr<pdat::CellData<int> >& tags,
    const boost::shared_ptr<pdat::CellVariable<double> >& variable_value_tagger,
    const double value_max,
    const bool uses_global_tol_up,
    const bool uses_global_tol_lo,
    const bool uses_local_tol_up,
    const bool uses_local_tol_lo,
    const double global_tol_up,
    const double global_tol_lo,
    const double local_tol_up,
    const double local_tol_lo)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(tags->getGhostCellWidth() == hier::IntVector::getZero(d_dim));
#endif
    
    boost::shared_ptr<pdat::CellData<double> > data_value_tagger(
        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
            patch.getPatchData(variable_value_tagger, data_context)));
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the number of ghost cells and dimensions of boxes that cover interior of patch plus
    // ghost cells.
    const hier::IntVector num_ghosts_value_tagger = data_value_tagger->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_value_tagger = data_value_tagger->getGhostBox().numberCells();
    
    boost::shared_ptr<pdat::CellData<int> > tags_value_tagger(
        new pdat::CellData<int>(interior_box, d_dim.getValue(), hier::IntVector::getZero(d_dim)));
    
    tags_value_tagger->fillAll(1);
    
    // Get the pointer to the tags.
    int* tag_ptr_value_tagger = tags_value_tagger->getPointer(0);
    int* tag_ptr  = tags->getPointer(0);
    
    // Get the pointer to the data.
    double* u = data_value_tagger->getPointer(0);
    
    if (d_dim == tbox::Dimension(1))
    {
        const int interior_dim_0 = interior_dims[0];
        
        const int num_ghosts_0_value_tagger = num_ghosts_value_tagger[0];
        
        if (uses_global_tol_up)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx = i + num_ghosts_0_value_tagger;
                const int idx_nghost = i;
                
                if (u[idx]/(value_max + EPSILON) <= global_tol_up)
                {
                    tag_ptr_value_tagger[idx_nghost] &= 1;
                }
                else
                {
                    tag_ptr_value_tagger[idx_nghost] &= 0;
                }
            }
        }
        
        if (uses_global_tol_lo)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx = i + num_ghosts_0_value_tagger;
                const int idx_nghost = i;
                
                if (u[idx]/(value_max + EPSILON) >= global_tol_lo)
                {
                    tag_ptr_value_tagger[idx_nghost] &= 1;
                }
                else
                {
                    tag_ptr_value_tagger[idx_nghost] &= 0;
                }
            }
        }
        
        if (uses_local_tol_up)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx = i + num_ghosts_0_value_tagger;
                const int idx_nghost = i;
                
                if (u[idx] <= local_tol_up)
                {
                    tag_ptr_value_tagger[idx_nghost] &= 1;
                }
                else
                {
                    tag_ptr_value_tagger[idx_nghost] &= 0;
                }
            }
        }
        
        if (uses_local_tol_lo)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx = i + num_ghosts_0_value_tagger;
                const int idx_nghost = i;
                
                if (u[idx] >= local_tol_lo)
                {
                    tag_ptr_value_tagger[idx_nghost] &= 1;
                }
                else
                {
                    tag_ptr_value_tagger[idx_nghost] &= 0;
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
            
            tag_ptr[idx_nghost] |= tag_ptr_value_tagger[idx_nghost];
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        
        const int num_ghosts_0_value_tagger = num_ghosts_value_tagger[0];
        const int num_ghosts_1_value_tagger = num_ghosts_value_tagger[1];
        const int ghostcell_dim_0_value_tagger = ghostcell_dims_value_tagger[0];
        
        for (int j = 0; j < interior_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx = (i + num_ghosts_0_value_tagger) +
                    (j + num_ghosts_1_value_tagger)*ghostcell_dim_0_value_tagger;
                
                const int idx_nghost = i +
                    j*interior_dim_0;
                
                int tag_cell = 1;
                
                if (uses_global_tol_up)
                {
                    if (u[idx]/(value_max + EPSILON) <= global_tol_up)
                    {
                        tag_cell &= 1;
                    }
                    else
                    {
                        tag_cell &= 0;
                    }
                }
                
                if (uses_global_tol_lo)
                {
                    if (u[idx]/(value_max + EPSILON) >= global_tol_lo)
                    {
                        tag_cell &= 1;
                    }
                    else
                    {
                        tag_cell &= 0;
                    }
                }
                
                if (uses_local_tol_up)
                {
                    if (u[idx] <= local_tol_up)
                    {
                        tag_cell &= 1;
                    }
                    else
                    {
                        tag_cell &= 0;
                    }
                }
                
                if (uses_local_tol_lo)
                {
                    if (u[idx] >= local_tol_lo)
                    {
                        tag_cell &= 1;
                    }
                    else
                    {
                        tag_cell &= 0;
                    }
                }
                
                tag_ptr[idx_nghost] |= tag_cell;
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        const int interior_dim_2 = interior_dims[2];
        
        const int num_ghosts_0_value_tagger = num_ghosts_value_tagger[0];
        const int num_ghosts_1_value_tagger = num_ghosts_value_tagger[1];
        const int num_ghosts_2_value_tagger = num_ghosts_value_tagger[2];
        const int ghostcell_dim_0_value_tagger = ghostcell_dims_value_tagger[0];
        const int ghostcell_dim_1_value_tagger = ghostcell_dims_value_tagger[1];
        
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
                    const int idx = (i + num_ghosts_0_value_tagger) +
                        (j + num_ghosts_1_value_tagger)*ghostcell_dim_0_value_tagger +
                        (k + num_ghosts_2_value_tagger)*ghostcell_dim_0_value_tagger*
                            ghostcell_dim_1_value_tagger;
                    
                    const int idx_nghost = i +
                        j*interior_dim_0 +
                        k*interior_dim_0*interior_dim_1;
                    
                    int tag_cell = 1;
                    
                    if (uses_global_tol_up)
                    {
                        if (u[idx]/(value_max + EPSILON) <= global_tol_up)
                        {
                            tag_cell &= 1;
                        }
                        else
                        {
                            tag_cell &= 0;
                        }
                    }
                    
                    if (uses_global_tol_lo)
                    {
                        if (u[idx]/(value_max + EPSILON) >= global_tol_lo)
                        {
                            tag_cell &= 1;
                        }
                        else
                        {
                            tag_cell &= 0;
                        }
                    }
                    
                    if (uses_local_tol_up)
                    {
                        if (u[idx] <= local_tol_up)
                        {
                            tag_cell &= 1;
                        }
                        else
                        {
                            tag_cell &= 0;
                        }
                    }
                    
                    if (uses_local_tol_lo)
                    {
                        if (u[idx] >= local_tol_lo)
                        {
                            tag_cell &= 1;
                        }
                        else
                        {
                            tag_cell &= 0;
                        }
                    }
                    
                    tag_ptr[idx_nghost] |= tag_cell;
                }
            }
        }
    }
}


/*
 * Transfer data input to data in class variable.
 */
void
ValueTagger::transferDataToClassVariable(
    hier::Patch& patch,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const boost::shared_ptr<pdat::CellData<double> >& data_input,
    const boost::shared_ptr<pdat::CellVariable<double> >& variable_value_tagger,
    const int depth)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_input->getDepth() > depth);
#endif
    
    boost::shared_ptr<pdat::CellData<double> > data_value_tagger(
        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
            patch.getPatchData(variable_value_tagger, data_context)));
    
    // Get the snumber of ghost cells and dimensions of box that covers interior of patch plus
    // ghost cells.
    const hier::IntVector num_ghosts_input = data_input->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_input = data_input->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_value_tagger = data_value_tagger->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_value_tagger = data_value_tagger->getGhostBox().numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(num_ghosts_input >= num_ghosts_value_tagger);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the pointer to the data.
    double* u_input = data_input->getPointer(depth);
    double* u_value_tagger = data_value_tagger->getPointer(0);
    
    if (d_dim == tbox::Dimension(1))
    {
        const int interior_dim_0 = interior_dims[0];
        
        const int num_ghosts_0_input = num_ghosts_input[0];
        const int num_ghosts_0_value_tagger = num_ghosts_value_tagger[0];
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = -num_ghosts_0_value_tagger;
             i < interior_dim_0 + num_ghosts_0_value_tagger;
             i++)
        {
            // Compute the linear indices.
            const int idx_input = i + num_ghosts_0_input;
            const int idx_value_tagger = i + num_ghosts_0_value_tagger;
            
            u_value_tagger[idx_value_tagger] = u_input[idx_input];
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        
        const int num_ghosts_0_input = num_ghosts_input[0];
        const int num_ghosts_1_input = num_ghosts_input[1];
        const int ghostcell_dim_0_input = ghostcell_dims_input[0];
        
        const int num_ghosts_0_value_tagger = num_ghosts_value_tagger[0];
        const int num_ghosts_1_value_tagger = num_ghosts_value_tagger[1];
        const int ghostcell_dim_0_value_tagger = ghostcell_dims_value_tagger[0];
        
        for (int j = -num_ghosts_1_value_tagger;
             j < interior_dim_1 + num_ghosts_1_value_tagger;
             j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = -num_ghosts_0_value_tagger;
                 i < interior_dim_0 + num_ghosts_0_value_tagger;
                 i++)
            {
                // Compute the linear indices.
                const int idx_input = (i + num_ghosts_0_input) +
                    (j + num_ghosts_1_input)*ghostcell_dim_0_input;
                
                const int idx_value_tagger = (i + num_ghosts_0_value_tagger) +
                    (j + num_ghosts_1_value_tagger)*ghostcell_dim_0_value_tagger;
                
                u_value_tagger[idx_value_tagger] = u_input[idx_input];
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        const int interior_dim_2 = interior_dims[2];
        
        const int num_ghosts_0_input = num_ghosts_input[0];
        const int num_ghosts_1_input = num_ghosts_input[1];
        const int num_ghosts_2_input = num_ghosts_input[2];
        const int ghostcell_dim_0_input = ghostcell_dims_input[0];
        const int ghostcell_dim_1_input = ghostcell_dims_input[1];
        
        const int num_ghosts_0_value_tagger = num_ghosts_value_tagger[0];
        const int num_ghosts_1_value_tagger = num_ghosts_value_tagger[1];
        const int num_ghosts_2_value_tagger = num_ghosts_value_tagger[2];
        const int ghostcell_dim_0_value_tagger = ghostcell_dims_value_tagger[0];
        const int ghostcell_dim_1_value_tagger = ghostcell_dims_value_tagger[1];
        
        for (int k = -num_ghosts_2_value_tagger;
             k < interior_dim_2 + num_ghosts_2_value_tagger;
             k++)
        {
            for (int j = -num_ghosts_1_value_tagger;
                 j < interior_dim_1 + num_ghosts_1_value_tagger;
                 j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = -num_ghosts_0_value_tagger;
                     i < interior_dim_0 + num_ghosts_0_value_tagger;
                     i++)
                {
                    // Compute the linear indices.
                    const int idx_input = (i + num_ghosts_0_input) +
                        (j + num_ghosts_1_input)*ghostcell_dim_0_input +
                        (k + num_ghosts_2_input)*ghostcell_dim_0_input*
                            ghostcell_dim_1_input;
                    
                    const int idx_value_tagger = (i + num_ghosts_0_value_tagger) +
                        (j + num_ghosts_1_value_tagger)*ghostcell_dim_0_value_tagger +
                        (k + num_ghosts_2_value_tagger)*ghostcell_dim_0_value_tagger*
                            ghostcell_dim_1_value_tagger;
                    
                    u_value_tagger[idx_value_tagger] = u_input[idx_input];
                }
            }
        }
    }
}
