#include "flow/refinement_taggers/ValueTagger.hpp"

#include <algorithm>

// #define HAMERS_PLOTTING_VALUE_TAGGER

#define EPSILON HAMERS_EPSILON

ValueTagger::ValueTagger(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
    const HAMERS_SHARED_PTR<FlowModel>& flow_model,
    const HAMERS_SHARED_PTR<tbox::Database>& value_tagger_db):
        d_object_name(object_name),
        d_dim(dim),
        d_grid_geometry(grid_geometry),
        d_num_value_ghosts(hier::IntVector::getZero(d_dim)),
        d_flow_model(flow_model),
        d_num_ghosts_derivative(1),
        d_value_tagger_max_density(Real(0)),
        d_value_tagger_max_total_energy(Real(0)),
        d_value_tagger_max_pressure(Real(0)),
        d_value_tagger_max_dilatation(Real(0)),
        d_value_tagger_max_enstrophy(Real(0))
{
    if (value_tagger_db != nullptr)
    {
        if (value_tagger_db->keyExists("num_ghosts_derivative"))
        {
            d_num_ghosts_derivative = value_tagger_db->getInteger("num_ghosts_derivative");
        }
        else if (value_tagger_db->keyExists("d_num_ghosts_derivative"))
        {
            d_num_ghosts_derivative = value_tagger_db->getInteger("d_num_ghosts_derivative");
        }
        
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
                  (variable_key == "MASS_FRACTION") || (variable_key == "MASS_FRACTIONS")))
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
                << "The numbers of variables and switches for global upper tolerances"
                << " provided don't match"
                << " in database of value tagger."
                << std::endl);
        }
        
        num_true = int(std::count(d_uses_global_tol_up.begin(), d_uses_global_tol_up.end(), true));
        
        if (num_true > 0)
        {
            if (value_tagger_db->keyExists("global_tol_up"))
            {
                d_global_tol_up = value_tagger_db->getRealVector("global_tol_up");
            }
            else if (value_tagger_db->keyExists("d_global_tol_up"))
            {
                d_global_tol_up = value_tagger_db->getRealVector("d_global_tol_up");
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
        
        num_true = int(std::count(d_uses_global_tol_lo.begin(), d_uses_global_tol_lo.end(), true));
        
        if (num_true > 0)
        {
            if (value_tagger_db->keyExists("global_tol_lo"))
            {
                d_global_tol_lo = value_tagger_db->getRealVector("global_tol_lo");
            }
            else if (value_tagger_db->keyExists("d_global_tol_lo"))
            {
                d_global_tol_lo = value_tagger_db->getRealVector("d_global_tol_lo");
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
        
        num_true = int(std::count(d_uses_local_tol_up.begin(), d_uses_local_tol_up.end(), true));
        
        if (num_true > 0)
        {
            if (value_tagger_db->keyExists("local_tol_up"))
            {
                d_local_tol_up = value_tagger_db->getRealVector("local_tol_up");
            }
            else if (value_tagger_db->keyExists("d_local_tol_up"))
            {
                d_local_tol_up = value_tagger_db->getRealVector("d_local_tol_up");
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
        
        num_true = int(std::count(d_uses_local_tol_lo.begin(), d_uses_local_tol_lo.end(), true));
        
        if (num_true > 0)
        {
            if (value_tagger_db->keyExists("local_tol_lo"))
            {
                d_local_tol_lo = value_tagger_db->getRealVector("local_tol_lo");
            }
            else if (value_tagger_db->keyExists("d_local_tol_lo"))
            {
                d_local_tol_lo = value_tagger_db->getRealVector("d_local_tol_lo");
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
        
        if ((std::find(d_variables.begin(), d_variables.end(), "DILATATION") != d_variables.end()) ||
            (std::find(d_variables.begin(), d_variables.end(), "ENSTROPHY") != d_variables.end()))
        {
            d_num_value_ghosts = hier::IntVector::getOne(d_dim)*d_num_ghosts_derivative;
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
            d_value_tagger_variable_density = HAMERS_MAKE_SHARED<pdat::CellVariable<Real> >(
                d_dim,
                "Value tagger density",
                1);
        }
        else if (variable_key == "TOTAL_ENERGY")
        {
            d_value_tagger_variable_total_energy = HAMERS_MAKE_SHARED<pdat::CellVariable<Real> >(
                d_dim,
                "Value tagger total energy",
                1);
        }
        else if (variable_key == "PRESSURE")
        {
            d_value_tagger_variable_pressure = HAMERS_MAKE_SHARED<pdat::CellVariable<Real> >(
                d_dim,
                "Value tagger pressure",
                1);
        }
        else if (variable_key == "DILATATION")
        {
            d_value_tagger_variable_dilatation = HAMERS_MAKE_SHARED<pdat::CellVariable<Real> >(
                d_dim,
                "Value tagger dilatation",
                1);
        }
        else if (variable_key == "ENSTROPHY")
        {
            d_value_tagger_variable_enstrophy = HAMERS_MAKE_SHARED<pdat::CellVariable<Real> >(
                d_dim,
                "Value tagger enstrophy",
                1);
        }
        else if (variable_key == "MASS_FRACTION" || variable_key == "MASS_FRACTIONS")
        {
            const int num_species = d_flow_model->getNumberOfSpecies();
            
            d_value_tagger_variable_mass_fractions.reserve(num_species);
            
            for (int si = 0; si < num_species; si++)
            {
                d_value_tagger_variable_mass_fractions.push_back(
                    HAMERS_MAKE_SHARED<pdat::CellVariable<Real> >(
                        d_dim,
                        "Value tagger mass fraction " + std::to_string(si),
                        1));
            }
            
            d_value_tagger_max_mass_fractions.reserve(num_species);
            
            for (int si = 0; si < num_species; si++)
            {
                d_value_tagger_max_mass_fractions.push_back(Real(0));
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
                hier::IntVector::getZero(d_dim),
                hier::IntVector::getZero(d_dim),
                RungeKuttaLevelIntegrator::TEMPORARY,
                    d_grid_geometry,
                    "NO_COARSEN",
                    "NO_REFINE");
        }
        else if (variable_key == "TOTAL_ENERGY")
        {
            integrator->registerVariable(
                d_value_tagger_variable_total_energy,
                hier::IntVector::getZero(d_dim),
                hier::IntVector::getZero(d_dim),
                RungeKuttaLevelIntegrator::TEMPORARY,
                    d_grid_geometry,
                    "NO_COARSEN",
                    "NO_REFINE");
        }
        else if (variable_key == "PRESSURE")
        {
            integrator->registerVariable(
                d_value_tagger_variable_pressure,
                hier::IntVector::getZero(d_dim),
                hier::IntVector::getZero(d_dim),
                RungeKuttaLevelIntegrator::TEMPORARY,
                    d_grid_geometry,
                    "NO_COARSEN",
                    "NO_REFINE");
        }
        else if (variable_key == "DILATATION")
        {
            integrator->registerVariable(
                d_value_tagger_variable_dilatation,
                hier::IntVector::getZero(d_dim),
                hier::IntVector::getZero(d_dim),
                RungeKuttaLevelIntegrator::TEMPORARY,
                    d_grid_geometry,
                    "NO_COARSEN",
                    "NO_REFINE");
        }
        else if (variable_key == "ENSTROPHY")
        {
            integrator->registerVariable(
                d_value_tagger_variable_enstrophy,
                hier::IntVector::getZero(d_dim),
                hier::IntVector::getZero(d_dim),
                RungeKuttaLevelIntegrator::TEMPORARY,
                    d_grid_geometry,
                    "NO_COARSEN",
                    "NO_REFINE");
        }
        else if (variable_key == "MASS_FRACTION" || variable_key == "MASS_FRACTIONS")
        {
            for (int si = 0; si < d_flow_model->getNumberOfSpecies(); si++)
            {
                integrator->registerVariable(
                    d_value_tagger_variable_mass_fractions[si],
                    hier::IntVector::getZero(d_dim),
                    hier::IntVector::getZero(d_dim),
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
    const HAMERS_SHARED_PTR<ExtendedVisItDataWriter>& visit_writer,
    const HAMERS_SHARED_PTR<hier::VariableContext>& plot_context)
{
#ifdef HAMERS_PLOTTING_VALUE_TAGGER
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
    const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const
{
    restart_db->putInteger("d_num_ghosts_derivative", d_num_ghosts_derivative);
    
    restart_db->putStringVector("d_variables", d_variables);
    
    restart_db->putBoolVector("d_uses_global_tol_up", d_uses_global_tol_up);
    restart_db->putBoolVector("d_uses_global_tol_lo", d_uses_global_tol_lo);
    restart_db->putBoolVector("d_uses_local_tol_up", d_uses_local_tol_up);
    restart_db->putBoolVector("d_uses_local_tol_lo", d_uses_local_tol_lo);
    
    int num_true = 0;
    
    num_true = int(std::count(d_uses_global_tol_up.begin(), d_uses_global_tol_up.end(), true));
    if (num_true > 0)
    {
        restart_db->putRealVector("d_global_tol_up", d_global_tol_up);
    }
    
    num_true = int(std::count(d_uses_global_tol_lo.begin(), d_uses_global_tol_lo.end(), true));
    if (num_true > 0)
    {
        restart_db->putRealVector("d_global_tol_lo", d_global_tol_lo);
    }
    
    num_true = int(std::count(d_uses_local_tol_up.begin(), d_uses_local_tol_up.end(), true));
    if (num_true > 0)
    {
        restart_db->putRealVector("d_local_tol_up", d_local_tol_up);
    }
    
    num_true = int(std::count(d_uses_local_tol_lo.begin(), d_uses_local_tol_lo.end(), true));
    if (num_true > 0)
    {
        restart_db->putRealVector("d_local_tol_lo", d_local_tol_lo);
    }
}


/*
 * Compute values on a patch for value tagger.
 */
void
ValueTagger::computeValueTaggerValuesOnPatch(
    hier::Patch& patch,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
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
            
            num_subghosts_of_data.insert(
                std::pair<std::string, hier::IntVector>("DENSITY", hier::IntVector::getZero(d_dim)));
            
            d_flow_model->registerDerivedVariables(num_subghosts_of_data);
            
            d_flow_model->allocateMemoryForDerivedCellData();
            
            d_flow_model->computeDerivedCellData();
            
            /*
             * Get the pointer to density data inside the flow model.
             */
            
            HAMERS_SHARED_PTR<pdat::CellData<Real> > flow_model_data_density = d_flow_model->getCellData("DENSITY");
            
            /*
             * Transfer data from flow model to the class variable.
             */
            
            transferDataOnPatchToClassVariable(
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
            
            num_subghosts_of_data.insert(
                std::pair<std::string, hier::IntVector>("TOTAL_ENERGY", hier::IntVector::getZero(d_dim)));
            
            d_flow_model->registerDerivedVariables(num_subghosts_of_data);
            
            d_flow_model->allocateMemoryForDerivedCellData();
            
            d_flow_model->computeDerivedCellData();
            
            /*
             * Get the pointer to total energy data inside the flow model.
             */
            
            HAMERS_SHARED_PTR<pdat::CellData<Real> > flow_model_data_total_energy =
                d_flow_model->getCellData("TOTAL_ENERGY");
            
            /*
             * Transfer data from flow model to the class variable.
             */
            
            transferDataOnPatchToClassVariable(
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
            
            num_subghosts_of_data.insert(
                std::pair<std::string, hier::IntVector>("PRESSURE", hier::IntVector::getZero(d_dim)));
            
            d_flow_model->registerDerivedVariables(num_subghosts_of_data);
            
            d_flow_model->allocateMemoryForDerivedCellData();
            
            d_flow_model->computeDerivedCellData();
            
            /*
             * Get the pointer to pressure data inside the flow model.
             */
            
            HAMERS_SHARED_PTR<pdat::CellData<Real> > flow_model_data_pressure = d_flow_model->getCellData("PRESSURE");
            
            /*
             * Transfer data from flow model to the class variable.
             */
            
            transferDataOnPatchToClassVariable(
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
             * Register the patch and velocity in the flow model and compute the corresponding cell
             * data of dilatation.
             */
            
            d_flow_model->registerPatchWithDataContext(patch, data_context);
            
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
            const hier::IntVector& num_ghosts = d_flow_model->getNumberOfGhostCells();
            TBOX_ASSERT(num_ghosts >= hier::IntVector::getOne(d_dim)*d_num_ghosts_derivative);
            NULL_USE(num_ghosts);
#endif
            
            std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
            
            num_subghosts_of_data.insert(
                std::pair<std::string, hier::IntVector>("VELOCITY", d_num_value_ghosts));
            
            d_flow_model->registerDerivedVariables(num_subghosts_of_data);
            
            d_flow_model->allocateMemoryForDerivedCellData();
            
            d_flow_model->computeDerivedCellData();
            
            // Get the cell data.
            HAMERS_SHARED_PTR<pdat::CellData<Real> > dilatation(
                HAMERS_SHARED_PTR_CAST<pdat::CellData<Real>, hier::PatchData>(
                    patch.getPatchData(d_value_tagger_variable_dilatation, data_context)));
            
            HAMERS_SHARED_PTR<pdat::CellData<Real> > velocity = d_flow_model->getCellData("VELOCITY");
            
            // Get the dimensions of box that covers the interior of patch.
            const hier::Box interior_box = patch.getBox();
            const hier::IntVector interior_dims = interior_box.numberCells();
            
            // Get the numbers of ghost cells and the dimensions of the ghost cell box.
            const hier::IntVector num_ghosts_dilatation = dilatation->getGhostCellWidth();
            const hier::IntVector ghostcell_dims_dilatation =
                dilatation->getGhostBox().numberCells();
    
            // Initialize cell data for velocity derivatives.
            HAMERS_SHARED_PTR<pdat::CellData<Real> > velocity_derivatives(
                new pdat::CellData<Real>(interior_box, d_dim.getValue(), hier::IntVector::getZero(d_dim)));
            
            // Get the grid spacing.
            const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                    patch.getPatchGeometry()));
            
            const double* const dx = patch_geom->getDx();
            
            // Get the pointer to the cell data of dilatation.
            Real* theta = dilatation->getPointer(0);
            
            if (d_dim == tbox::Dimension(1))
            {
                /*
                 * Get the number of cells in each dimension and number of ghost cells.
                 */
                
                const int interior_dim_0 = interior_dims[0];
                
                const int num_ghosts_0_dilatation = num_ghosts_dilatation[0];
                
                HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_x(
                    new DerivativeFirstOrder(
                        "first order derivative in x-direction",
                        d_dim, DIRECTION::X_DIRECTION,
                        d_num_ghosts_derivative));
                
                // Compute dudx.
                derivative_first_order_x->computeDerivative(
                    velocity_derivatives,
                    velocity,
                    Real(dx[0]),
                    0,
                    0);
                
                // Get the pointer to the cell data of velocity derivative.
                Real* dudx = velocity_derivatives->getPointer(0);
                
                HAMERS_PRAGMA_SIMD
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx = i;
                    const int idx_dilatation = i + num_ghosts_0_dilatation;
                    
                    theta[idx_dilatation] = dudx[idx];
                }
            }
            else if (d_dim == tbox::Dimension(2))
            {
                /*
                 * Get the numbers of cells in each dimension and numbers of ghost cells.
                 */
                
                const int interior_dim_0 = interior_dims[0];
                const int interior_dim_1 = interior_dims[1];
                
                const int num_ghosts_0_dilatation = num_ghosts_dilatation[0];
                const int num_ghosts_1_dilatation = num_ghosts_dilatation[1];
                const int ghostcell_dim_0_dilatation = ghostcell_dims_dilatation[0];
                
                HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_x(
                    new DerivativeFirstOrder(
                        "first order derivative in x-direction",
                        d_dim, DIRECTION::X_DIRECTION,
                        d_num_ghosts_derivative));
                
                HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_y(
                    new DerivativeFirstOrder(
                        "first order derivative in y-direction",
                        d_dim, DIRECTION::Y_DIRECTION,
                        d_num_ghosts_derivative));
                
                // Compute dudx.
                derivative_first_order_x->computeDerivative(
                    velocity_derivatives,
                    velocity,
                    Real(dx[0]),
                    0,
                    0);
                
                // Compute dvdy.
                derivative_first_order_y->computeDerivative(
                    velocity_derivatives,
                    velocity,
                    Real(dx[1]),
                    1,
                    1);
                
                // Get the pointers to the cell data of velocity derivatives.
                Real* dudx = velocity_derivatives->getPointer(0);
                Real* dvdy = velocity_derivatives->getPointer(1);
                
                // Compute the dilatation.
                for (int j = 0; j < interior_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx = i + j*interior_dim_0;
                        
                        const int idx_dilatation = (i + num_ghosts_0_dilatation) +
                            (j + num_ghosts_1_dilatation)*ghostcell_dim_0_dilatation;
                        
                        theta[idx_dilatation] = dudx[idx] + dvdy[idx];
                    }
                }
            }
            else if (d_dim == tbox::Dimension(3))
            {
                /*
                 * Get the numbers of cells in each dimension and numbers of ghost cells.
                 */
                
                const int interior_dim_0 = interior_dims[0];
                const int interior_dim_1 = interior_dims[1];
                const int interior_dim_2 = interior_dims[2];
                
                const int num_ghosts_0_dilatation = num_ghosts_dilatation[0];
                const int num_ghosts_1_dilatation = num_ghosts_dilatation[1];
                const int num_ghosts_2_dilatation = num_ghosts_dilatation[2];
                const int ghostcell_dim_0_dilatation = ghostcell_dims_dilatation[0];
                const int ghostcell_dim_1_dilatation = ghostcell_dims_dilatation[1];
                
                HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_x(
                    new DerivativeFirstOrder(
                        "first order derivative in x-direction",
                        d_dim, DIRECTION::X_DIRECTION,
                        d_num_ghosts_derivative));
                
                HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_y(
                    new DerivativeFirstOrder(
                        "first order derivative in y-direction",
                        d_dim, DIRECTION::Y_DIRECTION,
                        d_num_ghosts_derivative));
                
                HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_z(
                    new DerivativeFirstOrder(
                        "first order derivative in z-direction",
                        d_dim, DIRECTION::Z_DIRECTION,
                        d_num_ghosts_derivative));
                
                // Compute dudx.
                derivative_first_order_x->computeDerivative(
                    velocity_derivatives,
                    velocity,
                    Real(dx[0]),
                    0,
                    0);
                
                // Compute dvdy.
                derivative_first_order_y->computeDerivative(
                    velocity_derivatives,
                    velocity,
                    Real(dx[1]),
                    1,
                    1);
                
                // Compute dwdz.
                derivative_first_order_z->computeDerivative(
                    velocity_derivatives,
                    velocity,
                    Real(dx[2]),
                    2,
                    2);
                
                // Get the pointers to the cell data of velocity derivatives.
                Real* dudx = velocity_derivatives->getPointer(0);
                Real* dvdy = velocity_derivatives->getPointer(1);
                Real* dwdz = velocity_derivatives->getPointer(2);
                
                // Compute the dilatation.
                for (int k = 0; k < interior_dim_2; k++)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx = i + j*interior_dim_0 + k*interior_dim_0*interior_dim_1;
                            
                            const int idx_dilatation = (i + num_ghosts_0_dilatation) +
                                (j + num_ghosts_1_dilatation)*ghostcell_dim_0_dilatation +
                                (k + num_ghosts_2_dilatation)*ghostcell_dim_0_dilatation*
                                    ghostcell_dim_1_dilatation;
                            
                            theta[idx_dilatation] = dudx[idx] + dvdy[idx] + dwdz[idx];
                        }
                    }
                }
            }
            
            /*
             * Unregister the patch and data of all registered derived cell variables in the flow model.
             */
            
            d_flow_model->unregisterPatch();
        }
        else if (variable_key == "ENSTROPHY")
        {
            /*
             * Register the patch and velocity in the flow model and compute the corresponding cell
             * data of enstrophy.
             */
            
            d_flow_model->registerPatchWithDataContext(patch, data_context);
            
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
            const hier::IntVector& num_ghosts = d_flow_model->getNumberOfGhostCells();
            TBOX_ASSERT(num_ghosts >= hier::IntVector::getOne(d_dim)*d_num_ghosts_derivative);
            NULL_USE(num_ghosts);
#endif
            
            std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
            
            num_subghosts_of_data.insert(
                std::pair<std::string, hier::IntVector>("VELOCITY", d_num_value_ghosts));
            
            d_flow_model->registerDerivedVariables(num_subghosts_of_data);
            
            d_flow_model->allocateMemoryForDerivedCellData();
            
            d_flow_model->computeDerivedCellData();
            
            // Get the cell data.
            HAMERS_SHARED_PTR<pdat::CellData<Real> > enstrophy(
                HAMERS_SHARED_PTR_CAST<pdat::CellData<Real>, hier::PatchData>(
                    patch.getPatchData(d_value_tagger_variable_enstrophy, data_context)));
            
            HAMERS_SHARED_PTR<pdat::CellData<Real> > velocity = d_flow_model->getCellData("VELOCITY");
            
            // Get the dimensions of box that covers the interior of patch.
            const hier::Box interior_box = patch.getBox();
            const hier::IntVector interior_dims = interior_box.numberCells();
            
            // Get the numbers of ghost cells and the dimensions of the ghost cell box.
            const hier::IntVector num_ghosts_enstrophy = enstrophy->getGhostCellWidth();
            const hier::IntVector ghostcell_dims_enstrophy =
                enstrophy->getGhostBox().numberCells();
    
            // Get the grid spacing.
            const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                    patch.getPatchGeometry()));
            
            const double* const dx = patch_geom->getDx();
            
            // Get the pointer to the cell data of enstrophy.
            Real* Omega = enstrophy->getPointer(0);
            
            if (d_dim == tbox::Dimension(1))
            {
                TBOX_ERROR(d_object_name
                    << ": ValueTagger::computeValueTaggerValuesOnPatch()\n"
                    << "Enstrophy cannot be found for one-dimensional flow."
                    << std::endl);
            }
            else if (d_dim == tbox::Dimension(2))
            {
                /*
                 * Get the numbers of cells in each dimension and numbers of ghost cells.
                 */
                
                const int interior_dim_0 = interior_dims[0];
                const int interior_dim_1 = interior_dims[1];
                
                const int num_ghosts_0_enstrophy = num_ghosts_enstrophy[0];
                const int num_ghosts_1_enstrophy = num_ghosts_enstrophy[1];
                const int ghostcell_dim_0_enstrophy = ghostcell_dims_enstrophy[0];
                
                // Initialize cell data for velocity derivatives.
                HAMERS_SHARED_PTR<pdat::CellData<Real> > velocity_derivatives(
                    new pdat::CellData<Real>(interior_box, d_dim.getValue()*d_dim.getValue() - d_dim.getValue(),
                        hier::IntVector::getZero(d_dim)));
                
                HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_x(
                    new DerivativeFirstOrder(
                        "first order derivative in x-direction",
                        d_dim, DIRECTION::X_DIRECTION,
                        d_num_ghosts_derivative));
                
                HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_y(
                    new DerivativeFirstOrder(
                        "first order derivative in y-direction",
                        d_dim, DIRECTION::Y_DIRECTION,
                        d_num_ghosts_derivative));
                
                // Compute dudy.
                derivative_first_order_y->computeDerivative(
                    velocity_derivatives,
                    velocity,
                    Real(dx[1]),
                    0,
                    0);
                
                // Compute dvdx.
                derivative_first_order_x->computeDerivative(
                    velocity_derivatives,
                    velocity,
                    Real(dx[0]),
                    1,
                    1);
                
                // Get the pointers to the cell data of velocity derivatives.
                Real* dudy = velocity_derivatives->getPointer(0);
                Real* dvdx = velocity_derivatives->getPointer(1);
                
                // Compute the enstrophy.
                for (int j = 0; j < interior_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx = i + j*interior_dim_0;
                        
                        const int idx_enstrophy = (i + num_ghosts_0_enstrophy) +
                            (j + num_ghosts_1_enstrophy)*ghostcell_dim_0_enstrophy;
                        
                        const Real omega = dvdx[idx] - dudy[idx];
                        Omega[idx_enstrophy] = omega*omega;
                    }
                }
            }
            else if (d_dim == tbox::Dimension(3))
            {
                /*
                 * Get the dimension.
                 */
                
                const int interior_dim_0 = interior_dims[0];
                const int interior_dim_1 = interior_dims[1];
                const int interior_dim_2 = interior_dims[2];
                
                const int num_ghosts_0_enstrophy = num_ghosts_enstrophy[0];
                const int num_ghosts_1_enstrophy = num_ghosts_enstrophy[1];
                const int num_ghosts_2_enstrophy = num_ghosts_enstrophy[2];
                const int ghostcell_dim_0_enstrophy = ghostcell_dims_enstrophy[0];
                const int ghostcell_dim_1_enstrophy = ghostcell_dims_enstrophy[1];
                
                // Initialize cell data for velocity derivatives.
                HAMERS_SHARED_PTR<pdat::CellData<Real> > velocity_derivatives(
                    new pdat::CellData<Real>(interior_box, d_dim.getValue()*d_dim.getValue() - d_dim.getValue(),
                        hier::IntVector::getZero(d_dim)));
                
                HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_x(
                    new DerivativeFirstOrder(
                        "first order derivative in x-direction",
                        d_dim, DIRECTION::X_DIRECTION,
                        d_num_ghosts_derivative));
                
                HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_y(
                    new DerivativeFirstOrder(
                        "first order derivative in y-direction",
                        d_dim, DIRECTION::Y_DIRECTION,
                        d_num_ghosts_derivative));
                
                HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_z(
                    new DerivativeFirstOrder(
                        "first order derivative in z-direction",
                        d_dim, DIRECTION::Z_DIRECTION,
                        d_num_ghosts_derivative));
                
                // Compute dudy.
                derivative_first_order_y->computeDerivative(
                    velocity_derivatives,
                    velocity,
                    Real(dx[1]),
                    0,
                    0);
                
                // Compute dudz.
                derivative_first_order_z->computeDerivative(
                    velocity_derivatives,
                    velocity,
                    Real(dx[2]),
                    1,
                    0);
                
                // Compute dvdx.
                derivative_first_order_x->computeDerivative(
                    velocity_derivatives,
                    velocity,
                    Real(dx[0]),
                    2,
                    1);
                
                // Compute dvdz.
                derivative_first_order_z->computeDerivative(
                    velocity_derivatives,
                    velocity,
                    Real(dx[2]),
                    3,
                    1);
                
                // Compute dwdx.
                derivative_first_order_x->computeDerivative(
                    velocity_derivatives,
                    velocity,
                    Real(dx[0]),
                    4,
                    2);
                
                // Compute dwdy.
                derivative_first_order_y->computeDerivative(
                    velocity_derivatives,
                    velocity,
                    Real(dx[1]),
                    5,
                    2);
                
                // Get the pointers to the cell data of velocity derivatives.
                Real* dudy = velocity_derivatives->getPointer(0);
                Real* dudz = velocity_derivatives->getPointer(1);
                Real* dvdx = velocity_derivatives->getPointer(2);
                Real* dvdz = velocity_derivatives->getPointer(3);
                Real* dwdx = velocity_derivatives->getPointer(4);
                Real* dwdy = velocity_derivatives->getPointer(5);
                
                // Compute the enstrophy.
                for (int k = 0; k < interior_dim_2; k++)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear index.
                            const int idx = i + j*interior_dim_0 + k*interior_dim_0*interior_dim_1;
                            
                            const int idx_enstrophy = (i + num_ghosts_0_enstrophy) +
                                (j + num_ghosts_1_enstrophy)*ghostcell_dim_0_enstrophy +
                                (k + num_ghosts_2_enstrophy)*ghostcell_dim_0_enstrophy*
                                    ghostcell_dim_1_enstrophy;
                            
                            const Real omega_x = dwdy[idx] - dvdz[idx];
                            const Real omega_y = dudz[idx] - dwdx[idx];
                            const Real omega_z = dvdx[idx] - dudy[idx];
                            
                            Omega[idx_enstrophy] = omega_x*omega_x + omega_y*omega_y + omega_z*omega_z;
                        }
                    }
                }
            }
            
            /*
             * Unregister the patch and data of all registered derived cell variables in the flow model.
             */
            
            d_flow_model->unregisterPatch();
        }
        else if (variable_key == "MASS_FRACTION" || variable_key == "MASS_FRACTIONS")
        {
            /*
             * Register the patch and mass fractions in the flow model and compute the corresponding cell data.
             */
            
            d_flow_model->registerPatchWithDataContext(patch, data_context);
            
            std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
            
            num_subghosts_of_data.insert(
                std::pair<std::string, hier::IntVector>("MASS_FRACTIONS", hier::IntVector::getZero(d_dim)));
            
            d_flow_model->registerDerivedVariables(num_subghosts_of_data);
            
            d_flow_model->allocateMemoryForDerivedCellData();
            
            d_flow_model->computeDerivedCellData();
            
            /*
             * Get the pointer to mass fraction data inside the flow model.
             */
            
            HAMERS_SHARED_PTR<pdat::CellData<Real> > flow_model_data_mass_fractions =
                d_flow_model->getCellData("MASS_FRACTIONS");
            
            /*
             * Transfer data from flow model to the class variable.
             */
            
            for (int si = 0; si < d_flow_model->getNumberOfSpecies(); si++)
            {
                transferDataOnPatchToClassVariable(
                    patch,
                    data_context,
                    flow_model_data_mass_fractions,
                    d_value_tagger_variable_mass_fractions[si],
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
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const int level_number,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    math::HierarchyCellDataOpsReal<Real> cell_double_operator(patch_hierarchy, level_number, level_number);
    
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
                
                Real rho_max_local = cell_double_operator.max(rho_id);
                d_value_tagger_max_density = Real(0);
                
                mpi.Allreduce(
                    &rho_max_local,
                    &d_value_tagger_max_density,
                    1,
                    MPI_DOUBLE,
                    MPI_MAX);
            }
            else if (variable_key == "TOTAL_ENERGY")
            {
                const int E_id = variable_db->mapVariableAndContextToIndex(
                    d_value_tagger_variable_total_energy,
                    data_context);
                
                Real E_max_local = cell_double_operator.max(E_id);
                d_value_tagger_max_total_energy = Real(0);
                
                mpi.Allreduce(
                    &E_max_local,
                    &d_value_tagger_max_total_energy,
                    1,
                    MPI_DOUBLE,
                    MPI_MAX);
            }
            else if (variable_key == "PRESSURE")
            {
                const int p_id = variable_db->mapVariableAndContextToIndex(
                    d_value_tagger_variable_pressure,
                    data_context);
                
                Real p_max_local = cell_double_operator.max(p_id);
                d_value_tagger_max_pressure = Real(0);
                
                mpi.Allreduce(
                    &p_max_local,
                    &d_value_tagger_max_pressure,
                    1,
                    MPI_DOUBLE,
                    MPI_MAX);
            }
            else if (variable_key == "DILATATION")
            {
                const int theta_id = variable_db->mapVariableAndContextToIndex(
                    d_value_tagger_variable_dilatation,
                    data_context);
                
                Real theta_max_local = cell_double_operator.max(theta_id);
                d_value_tagger_max_dilatation = Real(0);
                
                mpi.Allreduce(
                    &theta_max_local,
                    &d_value_tagger_max_dilatation,
                    1,
                    MPI_DOUBLE,
                    MPI_MAX);
            }
            else if (variable_key == "ENSTROPHY")
            {
                const int Omega_id = variable_db->mapVariableAndContextToIndex(
                    d_value_tagger_variable_enstrophy,
                    data_context);
                
                Real Omega_max_local = cell_double_operator.max(Omega_id);
                d_value_tagger_max_enstrophy = Real(0);
                
                mpi.Allreduce(
                    &Omega_max_local,
                    &d_value_tagger_max_enstrophy,
                    1,
                    MPI_DOUBLE,
                    MPI_MAX);
            }
            else if (variable_key == "MASS_FRACTION" || variable_key == "MASS_FRACTIONS")
            {
                for (int si = 0; si < d_flow_model->getNumberOfSpecies(); si++)
                {
                    const int Y_id = variable_db->mapVariableAndContextToIndex(
                        d_value_tagger_variable_mass_fractions[si],
                        data_context);
                    
                    Real Y_max_local = cell_double_operator.max(Y_id);
                    d_value_tagger_max_mass_fractions[si] = Real(0);
                    
                    mpi.Allreduce(
                        &Y_max_local,
                        &d_value_tagger_max_mass_fractions[si],
                        1,
                        MPI_DOUBLE,
                        MPI_MAX);
                }
            }
        }
    }
}


/*
 * Tag cells on a patch for refinement using value tagger.
 */
void
ValueTagger::tagCellsOnPatch(
   hier::Patch& patch,
   const HAMERS_SHARED_PTR<pdat::CellData<int> >& tags,
   const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
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
        
        Real global_tol_up = Real(0);
        Real global_tol_lo = Real(0);
        Real local_tol_up = Real(0);
        Real local_tol_lo = Real(0);
        
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
            tagCellsOnPatchWithValue(
                patch,
                data_context,
                tags,
                d_value_tagger_variable_density,
                d_value_tagger_max_density,
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
            tagCellsOnPatchWithValue(
                patch,
                data_context,
                tags,
                d_value_tagger_variable_total_energy,
                d_value_tagger_max_total_energy,
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
            tagCellsOnPatchWithValue(
                patch,
                data_context,
                tags,
                d_value_tagger_variable_pressure,
                d_value_tagger_max_pressure,
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
            tagCellsOnPatchWithValue(
                patch,
                data_context,
                tags,
                d_value_tagger_variable_dilatation,
                d_value_tagger_max_dilatation,
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
            tagCellsOnPatchWithValue(
                patch,
                data_context,
                tags,
                d_value_tagger_variable_enstrophy,
                d_value_tagger_max_enstrophy,
                uses_global_tol_up,
                uses_global_tol_lo,
                uses_local_tol_up,
                uses_local_tol_lo,
                global_tol_up,
                global_tol_lo,
                local_tol_up,
                local_tol_lo);
        }
        else if (variable_key == "MASS_FRACTION" || variable_key == "MASS_FRACTIONS")
        {
            for (int si = 0; si < d_flow_model->getNumberOfSpecies(); si++)
            {
                tagCellsOnPatchWithValue(
                    patch,
                    data_context,
                    tags,
                    d_value_tagger_variable_mass_fractions[si],
                    d_value_tagger_max_mass_fractions[si],
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
 * Tag cells on a patch for refinement using data values.
 */
void
ValueTagger::tagCellsOnPatchWithValue(
    hier::Patch& patch,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const HAMERS_SHARED_PTR<pdat::CellData<int> >& tags,
    const HAMERS_SHARED_PTR<pdat::CellVariable<Real> >& variable_value_tagger,
    const Real value_max,
    const bool uses_global_tol_up,
    const bool uses_global_tol_lo,
    const bool uses_local_tol_up,
    const bool uses_local_tol_lo,
    const Real global_tol_up,
    const Real global_tol_lo,
    const Real local_tol_up,
    const Real local_tol_lo)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(tags->getGhostCellWidth() == hier::IntVector::getZero(d_dim));
#endif
    
    HAMERS_SHARED_PTR<pdat::CellData<Real> > data_value_tagger(
        HAMERS_SHARED_PTR_CAST<pdat::CellData<Real>, hier::PatchData>(
            patch.getPatchData(variable_value_tagger, data_context)));
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the number of ghost cells and dimensions of box that covers interior of patch plus
    // ghost cells.
    const hier::IntVector num_ghosts_value_tagger = data_value_tagger->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_value_tagger = data_value_tagger->getGhostBox().numberCells();
    
    // Allocate temporary patch data.
    HAMERS_SHARED_PTR<pdat::CellData<int> > tags_value_tagger(
        new pdat::CellData<int>(interior_box, d_dim.getValue(), hier::IntVector::getZero(d_dim)));
    
    tags_value_tagger->fillAll(1);
    
    // Get the pointers to the tags.
    int* tag_ptr_value_tagger = tags_value_tagger->getPointer(0);
    int* tag_ptr = tags->getPointer(0);
    
    // Get the pointer to the data.
    Real* u = data_value_tagger->getPointer(0);
    
    if (d_dim == tbox::Dimension(1))
    {
        const int interior_dim_0 = interior_dims[0];
        
        const int num_ghosts_0_value_tagger = num_ghosts_value_tagger[0];
        
        if (uses_global_tol_up)
        {
            HAMERS_PRAGMA_SIMD
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
            HAMERS_PRAGMA_SIMD
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
            HAMERS_PRAGMA_SIMD
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
            HAMERS_PRAGMA_SIMD
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
        
        HAMERS_PRAGMA_SIMD
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
        
        if (uses_global_tol_up)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx = (i + num_ghosts_0_value_tagger) +
                        (j + num_ghosts_1_value_tagger)*ghostcell_dim_0_value_tagger;
                    
                    const int idx_nghost = i +
                        j*interior_dim_0;
                    
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
        }
        
        if (uses_global_tol_lo)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx = (i + num_ghosts_0_value_tagger) +
                        (j + num_ghosts_1_value_tagger)*ghostcell_dim_0_value_tagger;
                    
                    const int idx_nghost = i +
                        j*interior_dim_0;
                    
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
        }
        
        if (uses_local_tol_up)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx = (i + num_ghosts_0_value_tagger) +
                        (j + num_ghosts_1_value_tagger)*ghostcell_dim_0_value_tagger;
                    
                    const int idx_nghost = i +
                        j*interior_dim_0;
                    
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
        }
        
        if (uses_local_tol_lo)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx = (i + num_ghosts_0_value_tagger) +
                        (j + num_ghosts_1_value_tagger)*ghostcell_dim_0_value_tagger;
                    
                    const int idx_nghost = i +
                        j*interior_dim_0;
                    
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
        }
        
        for (int j = 0; j < interior_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear index.
                const int idx_nghost = i +
                    j*interior_dim_0;
                
                tag_ptr[idx_nghost] |= tag_ptr_value_tagger[idx_nghost];
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
        
        if (uses_global_tol_up)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
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
            }
        }
        
        if (uses_global_tol_lo)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
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
            }
        }
        
        if (uses_local_tol_up)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
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
            }
        }
        
        if (uses_local_tol_lo)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
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
            }
        }
        
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_nghost = i +
                        j*interior_dim_0 +
                        k*interior_dim_0*interior_dim_1;
                    
                    tag_ptr[idx_nghost] |= tag_ptr_value_tagger[idx_nghost];
                }
            }
        }
    }
}


/*
 * Transfer data input on a patch to data in class variable.
 */
void
ValueTagger::transferDataOnPatchToClassVariable(
    hier::Patch& patch,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_input,
    const HAMERS_SHARED_PTR<pdat::CellVariable<Real> >& variable_value_tagger,
    const int depth)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_input->getDepth() > depth);
#endif
    
    HAMERS_SHARED_PTR<pdat::CellData<Real> > data_value_tagger(
        HAMERS_SHARED_PTR_CAST<pdat::CellData<Real>, hier::PatchData>(
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
    Real* u_input = data_input->getPointer(depth);
    Real* u_value_tagger = data_value_tagger->getPointer(0);
    
    if (d_dim == tbox::Dimension(1))
    {
        const int interior_dim_0 = interior_dims[0];
        
        const int num_ghosts_0_input = num_ghosts_input[0];
        const int num_ghosts_0_value_tagger = num_ghosts_value_tagger[0];
        
        HAMERS_PRAGMA_SIMD
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
            HAMERS_PRAGMA_SIMD
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
                HAMERS_PRAGMA_SIMD
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
