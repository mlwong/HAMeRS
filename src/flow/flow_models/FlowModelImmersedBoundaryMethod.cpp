#include "flow/flow_models/FlowModelImmersedBoundaryMethod.hpp"

HAMERS_SHARED_PTR<pdat::CellVariable<int> > FlowModelImmersedBoundaryMethod::s_variable_mask;
HAMERS_SHARED_PTR<pdat::CellVariable<Real> > FlowModelImmersedBoundaryMethod::s_variable_wall_distance;
HAMERS_SHARED_PTR<pdat::CellVariable<Real> > FlowModelImmersedBoundaryMethod::s_variable_surface_normal;

FlowModelImmersedBoundaryMethod::FlowModelImmersedBoundaryMethod(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
    const int& num_species,
    const int& num_eqn,
    const HAMERS_SHARED_PTR<ImmersedBoundaries>& immersed_boundaries,
    const HAMERS_SHARED_PTR<tbox::Database>& immersed_boundary_method_db):
        d_object_name(object_name),
        d_dim(dim),
        d_grid_geometry(grid_geometry),
        d_num_species(num_species),
        d_num_eqn(num_eqn),
        d_immersed_boundaries(immersed_boundaries)
{
    /*
     * Initialize the variables.
     */
    
    s_variable_mask = HAMERS_SHARED_PTR<pdat::CellVariable<int> > (
        new pdat::CellVariable<int>(d_dim, "immersed boundary mask", 1));
    
    s_variable_wall_distance = HAMERS_SHARED_PTR<pdat::CellVariable<Real> > (
        new pdat::CellVariable<Real>(d_dim, "wall distance", 1));
    
    s_variable_surface_normal = HAMERS_SHARED_PTR<pdat::CellVariable<Real> > (
        new pdat::CellVariable<Real>(d_dim, "surface_normal", dim.getValue()));
}


/*
 * Register the immersed boundary method variables.
 */
void
FlowModelImmersedBoundaryMethod::registerImmersedBoundaryMethodVariables(
    RungeKuttaLevelIntegrator* integrator,
    const hier::IntVector& num_ghosts,
    const hier::IntVector& num_ghosts_intermediate)
{
    if (s_variable_mask == nullptr)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelImmersedBoundaryMethod::registerVariables()\n"
            << "The mask variable is not yet initialized."
            << std::endl);
    }
    
    integrator->registerVariable(
        s_variable_mask,
        num_ghosts,
        num_ghosts_intermediate,
        RungeKuttaLevelIntegrator::NO_FILL,
        d_grid_geometry,
        "NO_COARSEN",
        "NO_REFINE");
    
    integrator->registerVariable(
        s_variable_wall_distance,
        num_ghosts,
        num_ghosts_intermediate,
        RungeKuttaLevelIntegrator::NO_FILL,
        d_grid_geometry,
        "NO_COARSEN",
        "NO_REFINE");
    
    integrator->registerVariable(
        s_variable_surface_normal,
        num_ghosts,
        num_ghosts_intermediate,
        RungeKuttaLevelIntegrator::NO_FILL,
        d_grid_geometry,
        "NO_COARSEN",
        "NO_REFINE");
}


/*
 * Register the plotting quantities.
 */
#ifdef HAVE_HDF5
void
FlowModelImmersedBoundaryMethod::registerPlotQuantities(
    const HAMERS_SHARED_PTR<ExtendedVisItDataWriter>& visit_writer,
    const HAMERS_SHARED_PTR<hier::VariableContext>& plot_context)
{
    hier::VariableDatabase* vardb = hier::VariableDatabase::getDatabase();
    
    visit_writer->registerPlotQuantity(
        "IB mask",
        "SCALAR",
        vardb->mapVariableAndContextToIndex(
           s_variable_mask,
           plot_context));
    
    visit_writer->registerPlotQuantity(
        "wall distance",
        "SCALAR",
        vardb->mapVariableAndContextToIndex(
           s_variable_wall_distance,
           plot_context));
    
    if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
    {
        visit_writer->registerPlotQuantity(
            "surface normal",
            "VECTOR",
            vardb->mapVariableAndContextToIndex(
               s_variable_surface_normal,
               plot_context));
    }
}
#endif


/*
 * Set the immersed boundary method variables.
 */
void
FlowModelImmersedBoundaryMethod::setImmersedBoundaryMethodVariables(
    const hier::Box& domain,
    const double data_time,
    const bool initial_time,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    const hier::Patch& patch = flow_model_tmp->getRegisteredPatch();
    
    const HAMERS_SHARED_PTR<pdat::CellData<int> > data_mask(
        HAMERS_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
            patch.getPatchData(s_variable_mask, data_context)));
    
    const HAMERS_SHARED_PTR<pdat::CellData<Real> > data_wall_distance(
        HAMERS_SHARED_PTR_CAST<pdat::CellData<Real>, hier::PatchData>(
            patch.getPatchData(s_variable_wall_distance, data_context)));
    
    const HAMERS_SHARED_PTR<pdat::CellData<Real> > data_surface_normal(
        HAMERS_SHARED_PTR_CAST<pdat::CellData<Real>, hier::PatchData>(
            patch.getPatchData(s_variable_surface_normal, data_context)));
    
    d_immersed_boundaries->setImmersedBoundaryVariablesOnPatch(
        patch,
        domain,
        data_mask,
        data_wall_distance,
        data_surface_normal,
        data_time,
        initial_time);
}


/*
 * Get the cell data of the immersed boundary mask in the registered patch.
 */
HAMERS_SHARED_PTR<pdat::CellData<int> >
FlowModelImmersedBoundaryMethod::getCellDataOfImmersedBoundaryMask(
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    const hier::Patch& patch = flow_model_tmp->getRegisteredPatch();
    
    // Get the cell data of the registered variable density.
    HAMERS_SHARED_PTR<pdat::CellData<int> > data_mask(
        HAMERS_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
            patch.getPatchData(s_variable_mask, data_context)));
    
    return data_mask;
}
