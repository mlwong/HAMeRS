#include "flow/flow_models/FlowModelSourceUtilities.hpp"

FlowModelSourceUtilities::FlowModelSourceUtilities(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
    const int& num_species,
    const int& num_eqn,
    const boost::shared_ptr<tbox::Database>& flow_model_db):
        d_object_name(object_name),
        d_dim(dim),
        d_grid_geometry(grid_geometry),
        d_num_species(num_species),
        d_num_eqn(num_eqn),
        d_has_source_terms(false),
        d_derived_cell_data_computed(false),
        d_num_subghosts_source_terms(-hier::IntVector::getOne(d_dim)),
        d_subghost_box_source_terms(hier::Box::getEmptyBox(dim)),
        d_subghostcell_dims_source_terms(hier::IntVector::getZero(d_dim))
{
    if (flow_model_db->keyExists("has_source_terms"))
    {
        d_has_source_terms = flow_model_db->getBool("has_source_terms");
    }
    else if (flow_model_db->keyExists("d_has_source_terms"))
    {
        d_has_source_terms = flow_model_db->getBool("d_has_source_terms");
    }
}


/*
 * Check whether there are any source terms.
 */
bool
FlowModelSourceUtilities::hasSourceTerms() const
{
    return d_has_source_terms;
}


/*
 * Register the required variables for the computation of source terms in the registered patch.
 */
void
FlowModelSourceUtilities::registerDerivedVariablesForSourceTerms(
    const hier::IntVector& num_subghosts)
{
    NULL_USE(num_subghosts);
    
    TBOX_ERROR(d_object_name
        << ": FlowModelSourceUtilities::registerDerivedVariablesForSourceTerms()\n"
        << "Function is not yet implemented!"
        << std::endl);
}


/*
 * Allocate memory for cell data of different registered derived variables related to this
 * class in the registered patch.
 */
void
FlowModelSourceUtilities::allocateMemoryForDerivedCellData()
{
    TBOX_ERROR(d_object_name
        << ": FlowModelSourceUtilities::allocateMemoryForDerivedCellData()\n"
        << "Function is not yet implemented!"
        << std::endl);
}


/*
 * Clear cell data of different derived variables related to this class in the registered patch.
 */
void
FlowModelSourceUtilities::clearCellData()
{
}


/*
 * Compute cell data of different registered derived variables related to this class.
 */
void
FlowModelSourceUtilities::computeDerivedCellData()
{
    TBOX_ERROR(d_object_name
        << ": FlowModelSourceUtilities::computeDerivedCellData()\n"
        << "Function is not yet implemented!"
        << std::endl);
}


/*
 * Compute all source terms on a patch.
 */
void
FlowModelSourceUtilities::computeSourceTermsOnPatch(
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
FlowModelSourceUtilities::putToRestart(
    const boost::shared_ptr<tbox::Database>& restart_db) const
{
    restart_db->putBool("d_has_source_terms", d_has_source_terms);
}
