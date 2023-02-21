#include "flow/flow_models/FlowModel.hpp"

FlowModel::FlowModel(
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
        d_flow_model_db(flow_model_db),
        d_num_ghosts(-hier::IntVector::getOne(d_dim)),
        d_patch(nullptr),
        d_interior_box(hier::Box::getEmptyBox(d_dim)),
        d_ghost_box(hier::Box::getEmptyBox(d_dim)),
        d_interior_dims(hier::IntVector::getZero(d_dim)),
        d_ghostcell_dims(hier::IntVector::getZero(d_dim)),
        d_subdomain_box(hier::Box::getEmptyBox(d_dim)),
        d_derived_cell_data_computed(false)
{
}


/*
 * Put the characteristics of the flow model class into the restart database.
 */
void
FlowModel::putToRestart(
        const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const
{
    putToRestartBase(restart_db);
}


/*
 * Check whether a patch is registered or not.
 */
bool
FlowModel::hasRegisteredPatch() const
{
    if (d_patch == nullptr)
    {
        return false;
    }
    
    return true;
}


/*
 * Get registered patch.
 */
const hier::Patch&
FlowModel::getRegisteredPatch() const
{
    if (d_patch == nullptr)
    {
        TBOX_ERROR(d_object_name
        << ": FlowModel::getRegisteredPatch()\n"
        << "Patch is not yet registered!"
        << std::endl);
    }
    
    return *d_patch;
}


/*
 * Return HAMERS_SHARED_PTR to patch data context.
 */
const HAMERS_SHARED_PTR<hier::VariableContext>&
FlowModel::getDataContext() const
{
    if (d_patch == nullptr)
    {
        TBOX_ERROR(d_object_name
        << ": FlowModel::getDataContext()\n"
        << "Patch is not yet registered!"
        << std::endl);
    }
    
   return d_data_context;
}


/*
 * Get sub-domain box.
 */
const hier::Box&
FlowModel::getSubdomainBox() const
{
    if (d_patch == nullptr)
    {
        TBOX_ERROR(d_object_name
        << ": FlowModel::getSubdomainBox()\n"
        << "Patch is not yet registered!"
        << std::endl);
    }
    
    return d_subdomain_box;
}


/*
 * Set sub-domain box.
 */
void
FlowModel::setSubdomainBox(const hier::Box& subdomain_box)
{
    if (d_patch == nullptr)
    {
        TBOX_ERROR(d_object_name
        << ": FlowModel::getSubdomainBox()\n"
        << "Patch is not yet registered!"
        << std::endl);
    }
    
    TBOX_ASSERT(d_ghost_box.contains(subdomain_box));
    
    d_subdomain_box = subdomain_box;
}


/*
 * Setup the Riemann solver object.
 */
void
FlowModel::setupRiemannSolver()
{
    d_flow_model_riemann_solver->setFlowModel(shared_from_this());
}


/*
 * Setup the basic utilties object.
 */
void
FlowModel::setupBasicUtilities()
{
    d_flow_model_basic_utilities->setFlowModel(shared_from_this());
}


/*
 * Setup the diffusive flux utilties object.
 */
void
FlowModel::setupDiffusiveFluxUtilities()
{
    d_flow_model_diffusive_flux_utilities->setFlowModel(shared_from_this());
}


/*
 * Setup the source utilties object.
 */
void
FlowModel::setupSourceUtilities()
{
    d_flow_model_source_utilities->setFlowModel(shared_from_this());
}


/*
 * Setup the immersed boundary method object.
 */
void
FlowModel::setupImmersedBoundaryMethod()
{
    d_flow_model_immersed_boundary_method->setFlowModel(shared_from_this());
}


/*
 * Setup the monitoring statistics utilties object.
 */
void
FlowModel::setupMonitoringStatisticsUtilities()
{
    d_flow_model_monitoring_statistics_utilities->setFlowModel(shared_from_this());
}


/*
 * Setup the statistics utilties object.
 */
void
FlowModel::setupStatisticsUtilities()
{
    d_flow_model_statistics_utilities->setFlowModel(shared_from_this());
}


/*
 * Put the characteristics of the base flow model class into the restart database.
 */
void
FlowModel::putToRestartBase(
    const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const
{
    NULL_USE(restart_db);
}
