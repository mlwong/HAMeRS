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
    const HAMERS_SHARED_PTR<tbox::Database>& immersed_boundary_method_db,
    const HAMERS_SHARED_PTR<EquationOfStateMixingRules>& equation_of_state_mixing_rules):
        d_object_name(object_name),
        d_dim(dim),
        d_grid_geometry(grid_geometry),
        d_num_IBM_ghosts(hier::IntVector::getZero(d_dim)),
        d_num_species(num_species),
        d_num_eqn(num_eqn),
        d_bc_type_velocity(VELOCITY_IBC::NONE),
        d_bc_type_temperature(TEMPERATURE_IBC::NONE),
        d_immersed_boundaries(immersed_boundaries),
        d_equation_of_state_mixing_rules(equation_of_state_mixing_rules)
{
    /*
     * Hard-code the additional number of cells required by the immersed boundary method to be 3.
     * 2D: ceiling of sqrt(2) + 1
     * 3D: ceiling of sqrt(3) + 1
     */
    
    d_num_IBM_ghosts = hier::IntVector::getOne(d_dim)*3;
    
    /*
     * Initialize the variables.
     */
    
    s_variable_mask = HAMERS_SHARED_PTR<pdat::CellVariable<int> > (
        new pdat::CellVariable<int>(d_dim, "immersed boundary mask", 1));
    
    s_variable_wall_distance = HAMERS_SHARED_PTR<pdat::CellVariable<Real> > (
        new pdat::CellVariable<Real>(d_dim, "wall distance", 1));
    
    s_variable_surface_normal = HAMERS_SHARED_PTR<pdat::CellVariable<Real> > (
        new pdat::CellVariable<Real>(d_dim, "surface_normal", dim.getValue()));
    
    /*
     * Read the immersed boundary conditions.
     */
    
    if (immersed_boundary_method_db->keyExists("bc_type_velocity"))
    {
        const std::string bc_type_velocity_str = immersed_boundary_method_db->getString("bc_type_velocity");
        
        if (bc_type_velocity_str == "SLIP")
        {
            d_bc_type_velocity = VELOCITY_IBC::SLIP;
        }
        else if (bc_type_velocity_str == "NO_SLIP")
        {
            d_bc_type_velocity = VELOCITY_IBC::NO_SLIP;
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelImmersedBoundaryMethod::FlowModelImmersedBoundaryMethod()\n"
                << "Unknown 'bc_type_velocity' entry from input database: " << d_bc_type_velocity
                << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelImmersedBoundaryMethod::FlowModelImmersedBoundaryMethod()\n"
            << "Required 'bc_type_velocity' entry from input database missing."
            << std::endl);
    }
    
    if (immersed_boundary_method_db->keyExists("bc_type_temperature"))
    {
        const std::string bc_type_temperature_str = immersed_boundary_method_db->getString("bc_type_temperature");
        
        if (bc_type_temperature_str == "ADIABATIC")
        {
            d_bc_type_temperature = TEMPERATURE_IBC::ADIABATIC;
        }
        // ISOTHERMAL BC NOT IMPLEMENTED YET!
        else if (bc_type_temperature_str == "ISOTHERMAL")
        {
            d_bc_type_temperature = TEMPERATURE_IBC::ISOTHERMAL;
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelImmersedBoundaryMethod::FlowModelImmersedBoundaryMethod()\n"
                << "Unknown 'bc_type_temperature' entry from input database: " << d_bc_type_temperature
                << std::endl);
        }
    }
    else
    {
        // Default is adiabatic wall.
        // d_bc_type_temperature = TEMPERATURE_IBC::ADIABATIC;
        
        TBOX_ERROR(d_object_name
            << ": FlowModelImmersedBoundaryMethod::FlowModelImmersedBoundaryMethod()\n"
            << "Required 'bc_type_temperature' entry from input database missing."
            << std::endl);
    }
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
        data_time,
        initial_time,
        domain,
        data_mask,
        data_wall_distance,
        data_surface_normal);
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


/*
 * Set the immersed boundary method ghost cells for the cell data of conservative variables.
 */
void
FlowModelImmersedBoundaryMethod::setConservativeVariablesCellDataImmersedBoundaryGhosts(
    const hier::Box& domain,
    const double data_time,
    const bool initial_time,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context_IB)
{
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
    // Check whether a patch is already registered.
    if (!flow_model_tmp->hasRegisteredPatch())
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelImmersedBoundaryMethod::"
            << "setConservativeVariablesCellDataImmersedBoundaryGhosts()\n"
            << "No patch is registered yet."
            << std::endl);
    }
    
    const hier::Patch& patch = flow_model_tmp->getRegisteredPatch();
    
    const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > conservative_var_data =
        flow_model_tmp->getCellDataOfConservativeVariables();
    
    const HAMERS_SHARED_PTR<pdat::CellData<int> > data_mask(
        HAMERS_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
            patch.getPatchData(s_variable_mask, data_context_IB)));
    
    const HAMERS_SHARED_PTR<pdat::CellData<Real> > data_wall_distance(
        HAMERS_SHARED_PTR_CAST<pdat::CellData<Real>, hier::PatchData>(
            patch.getPatchData(s_variable_wall_distance, data_context_IB)));
    
    const HAMERS_SHARED_PTR<pdat::CellData<Real> > data_surface_normal(
        HAMERS_SHARED_PTR_CAST<pdat::CellData<Real>, hier::PatchData>(
            patch.getPatchData(s_variable_surface_normal, data_context_IB)));
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_cons_var = conservative_var_data[0]->getGhostBox();
    const hier::IntVector ghostcell_dims_cons_var = ghost_box_cons_var.numberCells();
    
    const hier::Box ghost_box_IB = data_mask->getGhostBox();
    const hier::IntVector ghostcell_dims_IB = ghost_box_IB.numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_cons_var(d_dim);
    hier::IntVector offset_IB(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_cons_var = conservative_var_data[0]->getGhostCellWidth();
        const hier::IntVector num_ghosts_IB = data_mask->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = conservative_var_data[0]->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        for (int vi = 0; vi < static_cast<int>(conservative_var_data.size()); vi++)
        {
            TBOX_ASSERT(num_ghosts_cons_var == conservative_var_data[vi]->getGhostCellWidth());
            TBOX_ASSERT(conservative_var_data[vi]->getBox().isSpatiallyEqual(interior_box));
        }
        
        TBOX_ASSERT(num_ghosts_IB == data_wall_distance->getGhostCellWidth());
        TBOX_ASSERT(num_ghosts_IB == data_surface_normal->getGhostCellWidth());
        
        TBOX_ASSERT(data_mask->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_wall_distance->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_surface_normal->getBox().isSpatiallyEqual(interior_box));
        
        TBOX_ASSERT(num_ghosts_cons_var >= d_num_IBM_ghosts);
#endif
        
        // const hier::IntVector num_ghosts_domain = num_ghosts_cons_var - d_num_IBM_ghosts;
        // domain_lo = -num_ghosts_domain;
        // domain_dims = interior_box.numberCells() + num_ghosts_domain*2;
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = interior_box.numberCells();
        
        offset_cons_var = num_ghosts_cons_var;
        offset_IB = num_ghosts_IB;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        for (int vi = 0; vi < static_cast<int>(conservative_var_data.size()); vi++)
        {
            TBOX_ASSERT(conservative_var_data[vi].contains(domain));
        }
        
        TBOX_ASSERT(data_mask.contains(domain));
        TBOX_ASSERT(data_wall_distance.contains(domain));
        TBOX_ASSERT(data_surface_normal.contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_cons_var = domain.lower() - ghost_box_cons_var.lower();
        offset_IB = domain.lower() - ghost_box_IB.lower();
    }
    
    setConservativeVariablesCellDataImmersedBoundaryGhosts(
        patch,
        data_time,
        initial_time,
        conservative_var_data,
        data_mask,
        data_wall_distance,
        data_surface_normal,
        offset_cons_var,
        offset_IB,
        ghostcell_dims_cons_var,
        ghostcell_dims_IB,
        domain_lo,
        domain_dims);
}
