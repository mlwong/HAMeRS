#include "flow/flow_models/FlowModelImmersedBoundaryMethod.hpp"

#include "TECIO.h"

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
    
    const SurfaceTriangulation& surface_triangulation = d_immersed_boundaries->getSurfaceTriangulation();
    const int num_nodes = static_cast<int>(surface_triangulation.nodes.size());
    
    if (num_nodes > 0)
    {
        d_surface_triangulation_dx_grid.assign(num_nodes, std::numeric_limits<double>::max());
        d_surface_triangulation_coor_ip_1.assign(num_nodes, {0.0, 0.0, 0.0});
        d_surface_triangulation_coor_ip_2.assign(num_nodes, {0.0, 0.0, 0.0});
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


/*
 * Output the surface triangulation.
 */
void
FlowModelImmersedBoundaryMethod::writeSurfaceTriangulationWithData(const std::string& file_name) const
{
#ifdef HAMERS_USE_TECIO
    const SurfaceTriangulation& surface_triangulation = d_immersed_boundaries->getSurfaceTriangulation();
    
    if (surface_triangulation.nodes.size() == 0)
    {
        TBOX_WARNING(d_object_name
            << ": FlowModelImmersedBoundaryMethod::writeSurfaceTriangulationWithData()\n"
            << "The surface triangulation is empty."
            << " No surface file will be written."
            << std::endl);
        return;
    }
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    if (mpi.getRank() == 0)
    {
        const std::vector<std::array<double, 3> >& nodes = surface_triangulation.nodes;
        const std::vector<std::array<int, 3> >& connectivities = surface_triangulation.connectivities;
        const std::vector<std::array<double, 3> >& normal_nodes = surface_triangulation.normal_nodes;
        
        // Check that size of nodes and normal_nodes are the same.
        if (nodes.size() != normal_nodes.size())
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelImmersedBoundaryMethod::writeSurfaceTriangulationWithData()\n"
                << "The size of nodes and normal_nodes are not the same."
                << std::endl);
        }
        
        INTEGER4 num_nodes = static_cast<INTEGER4>(nodes.size());
        INTEGER4 num_centroids = static_cast<INTEGER4>(connectivities.size());
        
        const std::string file_name_full = file_name + ".plt";
        
        INTEGER4 file_format = 0; // 0 == PLT, 1 == SZPLT
        INTEGER4 file_type = 0; // FULL = 0, GRID = 1, SOLUTION = 2
        INTEGER4 debug = 1;
        INTEGER4 v_is_double = 1; // float = 0, double = 1
        INTEGER4 d_is_double = 1; // float = 0, double = 1
        
        std::vector<std::string> variable_names = {
            "x",
            "y",
            "z",
            "node_normal_x",
            "node_normal_y",
            "node_normal_z",
            "d_surface_triangulation_dx_grid"
        };
        
        std::string variable_name_string = "";
        for (int i = 0; i < static_cast<int>(variable_names.size()); ++i)
        {
            variable_name_string += variable_names[i];
            if (i < static_cast<int>(variable_names.size()) - 1)
            {
                variable_name_string += " ";
            }
        }
        
        /*
         * Open the file and write the tecplot datafile  header information
         */
        INTEGER4 i = TECINI142((char*)"DATASET",
            (char*) variable_name_string.c_str(), // NOTE: Make sure and change valueLocation above if this changes.
            (char*) file_name_full.c_str(),
            (char*) ".",
            &file_format,
            &file_type,
            &debug,
            &v_is_double);
        
        INTEGER4 zone_type = 2; // FETRIANGLE
        INTEGER4 num_faces = 1; // Not used.
        INTEGER4 i_cell_max = 0; // Not used.
        INTEGER4 j_cell_max = 0; // Not used.
        INTEGER4 k_cell_max = 0; // Not used.
        double sol_time = 0.0;
        INTEGER4 strand_id = 0; // Static zone.
        INTEGER4 parent_zn = 0; // No parent.
        INTEGER4 is_block  = 1; // Block format.
        INTEGER4 n_fconns  = 0; // Not used.
        INTEGER4 f_nmode   = 0; // Not used.
        INTEGER4 shr_conn  = 0; // Not used.
        
        // cell-centered: 0, nodal: 1
        const std::vector<int> valueLocation(static_cast<int>(variable_names.size()), 1);
        
        /*
         * Write the zone header information.
         */
        i = TECZNE142((char*)"Zone",
            &zone_type,
            &num_nodes,
            &num_centroids,
            &num_faces,
            &i_cell_max,
            &j_cell_max,
            &k_cell_max,
            &sol_time,
            &strand_id,
            &parent_zn,
            &is_block,
            &n_fconns,
            &f_nmode,
            0,                     /* TotalNumFaceNodes */
            0,                     /* NumConnectedBoundaryFaces */
            0,                     /* TotalNumBoundaryConnections */
            NULL,                  /* PassiveVarList */
            valueLocation.data(),  /* ValueLocation = Nodal */
            NULL,                  /* SharVarFromZone */
            &shr_conn);
        
        /*
         * Write out the field data.
         */
        
        std::vector<int> connectivity_array;
        std::vector<double> x, y, z;
        std::vector<double> node_normal_x, node_normal_y, node_normal_z;
        
        for (int ni = 0; ni < num_nodes; ni++)
        {
            const std::array<double, 3>& node = nodes[ni];
            const std::array<double, 3>& normal_node = normal_nodes[ni];
            
            x.push_back(double(node[0]));
            y.push_back(double(node[1]));
            z.push_back(double(node[2]));
            node_normal_x.push_back(double(normal_node[0]));
            node_normal_y.push_back(double(normal_node[1]));
            node_normal_z.push_back(double(normal_node[2]));
        }
        
        for (const auto& conn : connectivities)
        {
            for (int i = 0; i < 3; ++i)
            {
                connectivity_array.push_back(conn[i]);
            }
        }
        INTEGER4 connectivity_count = static_cast<INTEGER4>(connectivities.size())*3;
        
        i = TECDAT142(&num_nodes, x.data(),                               &d_is_double);
        i = TECDAT142(&num_nodes, y.data(),                               &d_is_double);
        i = TECDAT142(&num_nodes, z.data(),                               &d_is_double);
        i = TECDAT142(&num_nodes, node_normal_x.data(),                   &d_is_double);
        i = TECDAT142(&num_nodes, node_normal_y.data(),                   &d_is_double);
        i = TECDAT142(&num_nodes, node_normal_z.data(),                   &d_is_double);
        i = TECDAT142(&num_nodes, d_surface_triangulation_dx_grid.data(), &d_is_double);
        
        i = TECNODE142(&connectivity_count, connectivity_array.data());
         
        i = TECEND142();
    }
#endif
}



/*
 * Compute the data on the surface triangulation.
 */
void
FlowModelImmersedBoundaryMethod::computeSurfaceTriangulationDataBase(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    // if (d_flow_model.expired())
    // {
    //     TBOX_ERROR(d_object_name
    //         << ": "
    //         << "The object is not setup yet!"
    //         << std::endl);
    // }
    
    // HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
    const SurfaceTriangulation& surface_triangulation = d_immersed_boundaries->getSurfaceTriangulation();
    const std::vector<std::array<double, 3> >& nodes = surface_triangulation.nodes;
    
    if (nodes.empty())
    {
        return;
    }
    
    const int num_nodes = static_cast<int>(nodes.size());
    d_surface_triangulation_dx_grid.assign(num_nodes, std::numeric_limits<double>::max());
    std::vector<double> dx_grid_local(num_nodes, std::numeric_limits<double>::max());
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    if (d_dim == tbox::Dimension(1))
    {
        // Do nothing for now.
    }
    else if (d_dim == tbox::Dimension(2))
    {
        // Do nothing for now.
    }
    else if (d_dim == tbox::Dimension(3))
    {
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                // Make sure dx is isotropic.
                const double dx_min = std::min(dx[0], std::min(dx[1], dx[2]));
                const double dx_max = std::max(dx[0], std::max(dx[1], dx[2]));
                const double dx_ratio = dx_max/dx_min;
                if (std::abs(dx_max - dx_min) > 10.0*std::numeric_limits<double>::epsilon())
                {
                    TBOX_ERROR(d_object_name
                        << ": FlowModelImmersedBoundaryMethod::computeSurfaceTriangulationData()\n"
                        << "The grid spacing is not isotropic."
                        << std::endl);
                }
                
                const double* const patch_xlo = patch_geom->getXLower();
                const double* const patch_xhi = patch_geom->getXUpper();
                
                for (int ni = 0; ni < num_nodes; ++ni)
                {
                    const std::array<double, 3>& node = nodes[ni];
                    
                    if ((node[0] >= patch_xlo[0] && node[0] <= patch_xhi[0]) &&
                        (node[1] >= patch_xlo[1] && node[1] <= patch_xhi[1]) &&
                        (node[2] >= patch_xlo[2] && node[2] <= patch_xhi[2]))
                    {
                        dx_grid_local[ni] = std::min(dx_grid_local[ni], dx[0]);
                    }
                }
            }
        }
        
        const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
        
        mpi.Allreduce(
            &dx_grid_local[0],
            &d_surface_triangulation_dx_grid[0],
            num_nodes,
            MPI_DOUBLE,
            MPI_MIN);
        
        if (mpi.getRank() == 0)
        {
            // Make sure d_surface_triangulation_dx_grid is uniform.
            double dx_grid_min = std::numeric_limits<double>::max();
            double dx_grid_max = std::numeric_limits<double>::min();
            for (int ni = 0; ni < num_nodes; ++ni)
            {
                dx_grid_min = std::min(dx_grid_min, d_surface_triangulation_dx_grid[ni]);
                dx_grid_max = std::max(dx_grid_max, d_surface_triangulation_dx_grid[ni]);
            }
            if (std::abs(dx_grid_max - dx_grid_min) > 10.0*std::numeric_limits<double>::epsilon())
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelImmersedBoundaryMethod::computeSurfaceTriangulationData()\n"
                    << "The surfce triangulation is not in the same grid level."
                    << std::endl);
            }
        }
        
        const std::vector<std::array<double, 3> >& normal_nodes = surface_triangulation.normal_nodes;
        
        const double c_ip_1 = sqrt(3.0);
        const double c_ip_2 = 2.0;
        
        for (int ni = 0; ni < num_nodes; ++ni)
        {
            const std::array<double, 3>& node = nodes[ni];
            const std::array<double, 3>& normal_node = normal_nodes[ni];
            const double dx_grid = d_surface_triangulation_dx_grid[ni];
            
            d_surface_triangulation_coor_ip_1[ni][0] = node[0] + normal_node[0]*c_ip_1*dx_grid;
            d_surface_triangulation_coor_ip_1[ni][1] = node[1] + normal_node[1]*c_ip_1*dx_grid;
            d_surface_triangulation_coor_ip_1[ni][2] = node[2] + normal_node[2]*c_ip_1*dx_grid;
            
            d_surface_triangulation_coor_ip_2[ni][0] = node[0] + normal_node[0]*c_ip_2*dx_grid;
            d_surface_triangulation_coor_ip_2[ni][1] = node[1] + normal_node[1]*c_ip_2*dx_grid;
            d_surface_triangulation_coor_ip_2[ni][2] = node[2] + normal_node[2]*c_ip_2*dx_grid;
        }
        
        std::vector<double> ip_1_dx_grid(num_nodes, std::numeric_limits<Real>::max());
        std::vector<double> ip_2_dx_grid(num_nodes, std::numeric_limits<Real>::max());
        std::vector<double> ip_1_dx_grid_local(num_nodes, std::numeric_limits<Real>::max());
        std::vector<double> ip_2_dx_grid_local(num_nodes, std::numeric_limits<Real>::max());
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                const double* const patch_xlo = patch_geom->getXLower();
                const double* const patch_xhi = patch_geom->getXUpper();
                
                for (int ni = 0; ni < num_nodes; ++ni)
                {
                    const std::array<Real, 3>& coor_ip_1 = d_surface_triangulation_coor_ip_1[ni];
                    const std::array<Real, 3>& coor_ip_2 = d_surface_triangulation_coor_ip_2[ni];
                    
                    if ((coor_ip_1[0] >= patch_xlo[0] && coor_ip_1[0] <= patch_xhi[0]) &&
                        (coor_ip_1[1] >= patch_xlo[1] && coor_ip_1[1] <= patch_xhi[1]) &&
                        (coor_ip_1[2] >= patch_xlo[2] && coor_ip_1[2] <= patch_xhi[2]))
                    {
                        ip_1_dx_grid_local[ni] = std::min(ip_1_dx_grid_local[ni], dx[0]);
                    }
                    if ((coor_ip_2[0] >= patch_xlo[0] && coor_ip_2[0] <= patch_xhi[0]) &&
                        (coor_ip_2[1] >= patch_xlo[1] && coor_ip_2[1] <= patch_xhi[1]) &&
                        (coor_ip_2[2] >= patch_xlo[2] && coor_ip_2[2] <= patch_xhi[2]))
                    {
                        ip_2_dx_grid_local[ni] = std::min(ip_2_dx_grid_local[ni], dx[0]);
                    }
                }
            }
        }
        
        mpi.Allreduce(
            &ip_1_dx_grid_local[0],
            &ip_1_dx_grid[0],
            num_nodes,
            MPI_DOUBLE,
            MPI_MIN);
        
        mpi.Allreduce(
            &ip_2_dx_grid_local[0],
            &ip_2_dx_grid[0],
            num_nodes,
            MPI_DOUBLE,
            MPI_MIN);
        
        if (mpi.getRank() == 0)
        {
            // Make sure d_surface_triangulation_dx_grid is uniform.
            double ip_1_dx_grid_min = std::numeric_limits<double>::max();
            double ip_1_dx_grid_max = std::numeric_limits<double>::min();
            double ip_2_dx_grid_min = std::numeric_limits<double>::max();
            double ip_2_dx_grid_max = std::numeric_limits<double>::min();
            for (int ni = 0; ni < num_nodes; ++ni)
            {
                ip_1_dx_grid_min = std::min(ip_1_dx_grid_min, ip_1_dx_grid[ni]);
                ip_1_dx_grid_max = std::max(ip_1_dx_grid_max, ip_1_dx_grid[ni]);
                
                ip_2_dx_grid_min = std::min(ip_2_dx_grid_min, ip_2_dx_grid[ni]);
                ip_2_dx_grid_max = std::max(ip_2_dx_grid_max, ip_2_dx_grid[ni]);
            }
            if (std::abs(ip_1_dx_grid_max - ip_1_dx_grid_min) > 10.0*std::numeric_limits<double>::epsilon())
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelImmersedBoundaryMethod::computeSurfaceTriangulationData()\n"
                    << "The first image points are not in the same grid level."
                    << std::endl);
            }
            if (std::abs(ip_2_dx_grid_max - ip_2_dx_grid_min) > 10.0*std::numeric_limits<double>::epsilon())
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelImmersedBoundaryMethod::computeSurfaceTriangulationData()\n"
                    << "The second image points are not in the same grid level."
                    << std::endl);
            }
        }
    }
}