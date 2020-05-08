#include "flow/flow_models/single-species/FlowModelSourceUtilitiesSingleSpecies.hpp"

FlowModelSourceUtilitiesSingleSpecies::FlowModelSourceUtilitiesSingleSpecies(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
    const int& num_species,
    const boost::shared_ptr<tbox::Database>& flow_model_db,
    const boost::shared_ptr<EquationOfStateMixingRules> equation_of_state_mixing_rules):
        FlowModelSourceUtilities(
            object_name,
            dim,
            grid_geometry,
            num_species,
            2 + dim.getValue(),
            flow_model_db),
    d_has_gravity(false),
    d_equation_of_state_mixing_rules(equation_of_state_mixing_rules)
{
    if (d_has_source_terms)
    {
        boost::shared_ptr<tbox::Database> source_terms_db;
        
        if (flow_model_db->keyExists("Source_terms"))
        {
            source_terms_db = flow_model_db->getDatabase("Source_terms");
        }
        else if (flow_model_db->keyExists("d_source_terms"))
        {
            source_terms_db = flow_model_db->getDatabase("d_source_terms");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "No key 'Source_terms'/'d_source_terms' found in data for flow model"
                << std::endl);
        }
        
        if (source_terms_db->keyExists("has_gravity"))
        {
            d_has_gravity = source_terms_db->getBool("has_gravity");
            if (d_has_gravity)
            {
                source_terms_db->getVector("gravity", d_gravity);
            }
        }
        else if (source_terms_db->keyExists("d_has_gravity"))
        {
            d_has_gravity = source_terms_db->getBool("d_has_gravity");
            if (d_has_gravity)
            {
                source_terms_db->getVector("d_gravity", d_gravity);
            }
        }
        
        if (d_has_gravity)
        {
            if (d_gravity.size() != d_dim.getValue())
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSourceUtilitiesSingleSpecies::FlowModelSourceUtilitiesSingleSpecies\n"
                    << "Size of 'gravity' or 'd_gravity' is not consistent with problem dimension."
                    << std::endl);
            }
        }
    }
}


/*
 * Register the required variables for the computation of source terms in the registered patch.
 */
void
FlowModelSourceUtilitiesSingleSpecies::registerDerivedVariablesForSource(
    const hier::IntVector& num_subghosts)
{
    NULL_USE(num_subghosts);
    
    TBOX_ERROR(d_object_name
        << ": FlowModelSourceUtilitiesSingleSpecies::registerDerivedVariablesForSource()\n"
        << "Function is not yet implemented!"
        << std::endl);
}


/*
 * Allocate memory for cell data of different registered derived variables related to this
 * class in the registered patch.
 */
void
FlowModelSourceUtilitiesSingleSpecies::allocateMemoryForDerivedCellData()
{
    TBOX_ERROR(d_object_name
        << ": FlowModelSourceUtilitiesSingleSpecies::allocateMemoryForDerivedCellData()\n"
        << "Function is not yet implemented!"
        << std::endl);
}


/*
 * Clear cell data of different derived variables related to this class in the registered patch.
 */
void
FlowModelSourceUtilitiesSingleSpecies::clearCellData()
{
}


/*
 * Compute cell data of different registered derived variables related to this class.
 */
void
FlowModelSourceUtilitiesSingleSpecies::computeDerivedCellData()
{
    TBOX_ERROR(d_object_name
        << ": FlowModelSourceUtilitiesSingleSpecies::computeDerivedCellData()\n"
        << "Function is not yet implemented!"
        << std::endl);
}


/*
 * Compute the source on a patch.
 */
void
FlowModelSourceUtilitiesSingleSpecies::computeSourceOnPatch(
    const boost::shared_ptr<pdat::CellVariable<double> >& variable_source,
    const double time,
    const double dt,
    const int RK_step_number)
{
    NULL_USE(variable_source);
    NULL_USE(time);
    NULL_USE(dt);
    NULL_USE(RK_step_number);
}


/*
 * Put the characteristics of this class into the restart database.
 */
void
FlowModelSourceUtilitiesSingleSpecies::putToRestart(
    const boost::shared_ptr<tbox::Database>& restart_db) const
{
    restart_db->putBool("d_has_source_terms", d_has_source_terms);
    
    restart_db->putBool("d_has_gravity", d_has_gravity);
    if (d_has_gravity)
    {
        restart_db->putVector("d_gravity", d_gravity);
    }
}
