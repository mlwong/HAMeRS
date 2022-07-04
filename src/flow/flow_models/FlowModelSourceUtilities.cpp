#include "flow/flow_models/FlowModelSourceUtilities.hpp"

FlowModelSourceUtilities::FlowModelSourceUtilities(
    const std::string& object_name,
    const std::string& project_name,
    const tbox::Dimension& dim,
    const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
    const int& num_species,
    const int& num_eqn,
    const HAMERS_SHARED_PTR<tbox::Database>& flow_model_db):
        d_object_name(object_name),
        d_project_name(project_name),
        d_dim(dim),
        d_grid_geometry(grid_geometry),
        d_num_species(num_species),
        d_num_eqn(num_eqn),
        d_has_source_terms(false),
        d_has_special_source_terms(false),
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
    
    if (d_has_source_terms)
    {
        HAMERS_SHARED_PTR<tbox::Database> source_terms_db;
        
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
                << "No key 'Source_terms'/'d_source_terms' found in data for flow model."
                << std::endl);
        }
        
        if (source_terms_db->keyExists("has_special_source_terms"))
        {
            d_has_special_source_terms = source_terms_db->getBool("has_special_source_terms");
        }
        else if (source_terms_db->keyExists("d_has_special_source_terms"))
        {
            d_has_special_source_terms = source_terms_db->getBool("d_has_special_source_terms");
        }
        
        if (d_has_special_source_terms)
        {
            d_special_source_terms.reset(new FlowModelSpecialSourceTerms(
                "d_special_source_terms",
                d_project_name,
                d_dim,
                d_grid_geometry,
                source_terms_db,
                d_num_eqn));
        }
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
 * Register the required variables for the computation of local stable time increment for
 * source terms in the registered patch.
 */
void
FlowModelSourceUtilities::registerDerivedVariablesForSourceTermsStableDt(
    const hier::IntVector& num_subghosts)
{
    NULL_USE(num_subghosts);
    
    TBOX_ERROR(d_object_name
        << ": FlowModelSourceUtilities::registerDerivedVariablesForSourceTermsStableDt()\n"
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
    const HAMERS_SHARED_PTR<pdat::CellVariable<double> >& variable_source,
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
 * Compute special source terms.
 */
void
FlowModelSourceUtilities::computeSpecialSourceTermsOnPatch(
    const HAMERS_SHARED_PTR<pdat::CellVariable<double> >& variable_source,
    const double time,
    const double dt,
    const int RK_step_number)
{
    if (d_has_special_source_terms)
    {
        HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
        const hier::Patch& patch = flow_model_tmp->getRegisteredPatch();
        const HAMERS_SHARED_PTR<hier::VariableContext> data_context = flow_model_tmp->getDataContext();
        
        const std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > conservative_var_data =
            flow_model_tmp->getCellDataOfConservativeVariables();
        
        // Get the cell data of source.
        HAMERS_SHARED_PTR<pdat::CellData<double> > source(
            HAMERS_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                patch.getPatchData(variable_source, data_context)));
        
        d_special_source_terms->computeSpecialSourceTermsOnPatch(
            source,
            patch,
            conservative_var_data,
            time,
            dt,
            RK_step_number);
    }
}


/*
 * Get local stable time increment for source terms.
 */
double
FlowModelSourceUtilities::getStableDtOnPatch()
{
    TBOX_ERROR(d_object_name
        << ": FlowModelSourceUtilities::getStableDtOnPatch()\n"
        << "Function is not yet implemented!"
        << std::endl);
    
    return double(0);
}


/*
 * Put the characteristics of this class into the restart database.
 */
void
FlowModelSourceUtilities::putToRestart(
    const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const
{
    putToRestartBase(restart_db);
    
    if (d_has_source_terms)
    {
        HAMERS_SHARED_PTR<tbox::Database> restart_source_terms_db =
            restart_db->putDatabase("d_source_terms");
        
        putToRestartSourceBase(restart_source_terms_db);
    }
}


/*
 * Put the characteristics of base class into the restart database.
 */
void
FlowModelSourceUtilities::putToRestartBase(
    const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const
{
    restart_db->putBool("d_has_source_terms", d_has_source_terms);
}


/*
 * Put the characteristics of base class into the restart source database.
 */
void
FlowModelSourceUtilities::putToRestartSourceBase(
    const HAMERS_SHARED_PTR<tbox::Database>& restart_source_terms_db) const
{
    restart_source_terms_db->putBool("d_has_special_source_terms", d_has_special_source_terms);
    
    if (d_has_special_source_terms)
    {
        d_special_source_terms->putToRestart(restart_source_terms_db);
    }
}
