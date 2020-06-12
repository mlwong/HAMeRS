#include "apps/Navier-Stokes/NavierStokesBoundaryConditions.hpp"

// routines for managing boundary data
#include "util/basic_boundary_conditions/BasicCartesianBoundaryUtilities1.hpp"
#include "util/basic_boundary_conditions/BasicCartesianBoundaryUtilities2.hpp"
#include "util/basic_boundary_conditions/BasicCartesianBoundaryUtilities3.hpp"

//integer constant for debugging improperly set boundary data
#define BOGUS_BDRY_DATA (-9999)

NavierStokesBoundaryConditions::NavierStokesBoundaryConditions(
    const std::string& object_name,
    const std::string& project_name,
    const tbox::Dimension& dim,
    const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
    const FLOW_MODEL::TYPE& flow_model_type,
    const boost::shared_ptr<FlowModel>& flow_model,
    const boost::shared_ptr<tbox::Database>& boundary_conditions_db,
    const bool& boundary_conditions_db_is_from_restart):
        d_object_name(object_name),
        d_project_name(project_name),
        d_dim(dim),
        d_grid_geometry(grid_geometry),
        d_flow_model_type(flow_model_type),
        d_flow_model(flow_model)
{
    /*
     * Defaults for boundary conditions. Set to bogus values for error checking.
     */
    setDefaultBoundaryConditions();
    
    const hier::IntVector &one_vec = hier::IntVector::getOne(d_dim);
    hier::IntVector periodic = d_grid_geometry->getPeriodicShift(one_vec);
    int num_per_dirs = 0;
    for (int di = 0; di < d_dim.getValue(); di++)
    {
        if (periodic(di))
        {
            num_per_dirs++;
        }
    }
    
    /*
     * Get the boost::shared_ptr to flow model boundary utilities object from d_flow_model.
     */
    
    const boost::shared_ptr<FlowModelBoundaryUtilities> flow_model_boundary_utilities =
        d_flow_model->getFlowModelBoundaryUtilities();
    
    if (num_per_dirs < d_dim.getValue())
    {
        if (boundary_conditions_db_is_from_restart)
        {
            std::vector<boost::shared_ptr<pdat::CellVariable<double> > > conservative_var =
                d_flow_model->getConservativeVariables();
            
            d_master_bdry_node_conds = boundary_conditions_db->getIntegerVector("d_master_bdry_node_conds");
            
            if (d_dim == tbox::Dimension(1))
            {
                for (int vi = 0; vi < static_cast<int>(conservative_var.size()); vi++)
                {
                    d_bdry_node_conservative_var[vi] = boundary_conditions_db->getDoubleVector(
                        "d_bdry_node_conservative_var[" + tbox::Utilities::intToString(vi) + "]");
                }
            }
            else if (d_dim == tbox::Dimension(2))
            {
                d_master_bdry_edge_conds = boundary_conditions_db->getIntegerVector("d_master_bdry_edge_conds");
                
                for (int vi = 0; vi < static_cast<int>(conservative_var.size()); vi++)
                {
                    d_bdry_edge_conservative_var[vi] = boundary_conditions_db->getDoubleVector(
                        "d_bdry_edge_conservative_var[" + tbox::Utilities::intToString(vi) + "]");
                }
            }
            else if (d_dim == tbox::Dimension(3))
            {
                d_master_bdry_edge_conds = boundary_conditions_db->getIntegerVector("d_master_bdry_edge_conds");
                d_master_bdry_face_conds = boundary_conditions_db->getIntegerVector("d_master_bdry_face_conds");
                
                for (int vi = 0; vi < static_cast<int>(conservative_var.size()); vi++)
                {
                    d_bdry_face_conservative_var[vi] = boundary_conditions_db->getDoubleVector(
                        "d_bdry_face_conservative_var[" + tbox::Utilities::intToString(vi) + "]");
                }
            }
        }
        else
        {
            /*
             * Get the boundary conditions from the input database.
             */
            if (d_dim == tbox::Dimension(1))
            {
                std::vector<int> node_locs;
                node_locs.reserve(NUM_1D_NODES);
                for (int ni = 0; ni < NUM_1D_NODES; ni++)
                {
                    node_locs.push_back(ni);
                }
                
                flow_model_boundary_utilities->
                    getFromInput1d(
                        boundary_conditions_db,
                        node_locs,
                        d_master_bdry_node_conds,
                        periodic);
                
                BasicCartesianBoundaryUtilities1::getFromInput(
                    this,
                    boundary_conditions_db,
                    node_locs,
                    d_master_bdry_node_conds,
                    periodic);
            }
            if (d_dim == tbox::Dimension(2))
            {
                std::vector<int> edge_locs;
                edge_locs.reserve(NUM_2D_EDGES);
                for (int ei = 0; ei < NUM_2D_EDGES; ei++)
                {
                    edge_locs.push_back(ei);
                }
                
                std::vector<int> node_locs;
                node_locs.reserve(NUM_2D_NODES);
                for (int ni = 0; ni < NUM_2D_NODES; ni++)
                {
                    node_locs.push_back(ni);
                }
                
                flow_model_boundary_utilities->
                    getFromInput2d(
                        boundary_conditions_db,
                        edge_locs,
                        node_locs,
                        d_master_bdry_edge_conds,
                        d_master_bdry_node_conds,
                        periodic);
                
                BasicCartesianBoundaryUtilities2::getFromInput(
                    this,
                    boundary_conditions_db,
                    edge_locs,
                    node_locs,
                    d_master_bdry_edge_conds,
                    d_master_bdry_node_conds,
                    periodic);
            }
            else if (d_dim == tbox::Dimension(3))
            {
                std::vector<int> face_locs;
                face_locs.reserve(NUM_3D_FACES);
                for (int fi = 0; fi < NUM_3D_FACES; fi++)
                {
                    face_locs.push_back(fi);
                }
                
                std::vector<int> edge_locs;
                edge_locs.reserve(NUM_3D_EDGES);
                for (int ei = 0; ei < NUM_3D_EDGES; ei++)
                {
                    edge_locs.push_back(ei);
                }
                
                std::vector<int> node_locs;
                node_locs.reserve(NUM_3D_NODES);
                for (int ni = 0; ni < NUM_3D_NODES; ni++)
                {
                    node_locs.push_back(ni);
                }
                
                flow_model_boundary_utilities->
                    getFromInput3d(
                        boundary_conditions_db,
                        face_locs,
                        edge_locs,
                        node_locs,
                        d_master_bdry_face_conds,
                        d_master_bdry_edge_conds,
                        d_master_bdry_node_conds,
                        periodic);
                
                BasicCartesianBoundaryUtilities3::getFromInput(
                    this,
                    boundary_conditions_db,
                    face_locs,
                    edge_locs,
                    node_locs,
                    d_master_bdry_face_conds,
                    d_master_bdry_edge_conds,
                    d_master_bdry_node_conds,
                    periodic);
            }
        }
    } // if (num_per_dirs < d_dim.getValue())
    
    /*
     * Postprocess boundary data from input/restart values.
     */
    if (d_dim == tbox::Dimension(1))
    {
        for (int i = 0; i < NUM_1D_NODES; i++)
        {
            d_scalar_bdry_node_conds[i] = d_master_bdry_node_conds[i];
            d_vector_bdry_node_conds[i] = d_master_bdry_node_conds[i];
            
            if (d_master_bdry_node_conds[i] == BDRY_COND::BASIC::REFLECT)
            {
                d_scalar_bdry_node_conds[i] = BDRY_COND::BASIC::SYMMETRY;
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const boost::shared_ptr<FlowModelBoundaryUtilities> flow_model_boundary_utilities =
            d_flow_model->getFlowModelBoundaryUtilities();
        
        for (int i = 0; i < NUM_2D_EDGES; i++)
        {
            d_scalar_bdry_edge_conds[i] = d_master_bdry_edge_conds[i];
            d_vector_bdry_edge_conds[i] = d_master_bdry_edge_conds[i];
            
            if (d_master_bdry_edge_conds[i] == BDRY_COND::BASIC::REFLECT)
            {
                d_scalar_bdry_edge_conds[i] = BDRY_COND::BASIC::SYMMETRY;
            }
        }
        
        for (int i = 0; i < NUM_2D_NODES; i++)
        {
            d_scalar_bdry_node_conds[i] = d_master_bdry_node_conds[i];
            d_vector_bdry_node_conds[i] = d_master_bdry_node_conds[i];
            
            if (d_master_bdry_node_conds[i] == BDRY_COND::BASIC::XREFLECT)
            {
                d_scalar_bdry_node_conds[i] = BDRY_COND::BASIC::XSYMMETRY;
            }
            if (d_master_bdry_node_conds[i] == BDRY_COND::BASIC::YREFLECT)
            {
                d_scalar_bdry_node_conds[i] = BDRY_COND::BASIC::YSYMMETRY;
            }
            
            if (d_master_bdry_node_conds[i] != BOGUS_BDRY_DATA)
            {
                d_node_bdry_edge[i] =
                    flow_model_boundary_utilities->getEdgeLocationForNodeBdry(
                        i, d_master_bdry_node_conds[i]);
                
                if (d_node_bdry_edge[i] == -1)
                {
                    d_node_bdry_edge[i] =
                        BasicCartesianBoundaryUtilities2::getEdgeLocationForNodeBdry(
                            i, d_master_bdry_node_conds[i]);
                }
                
                if (d_node_bdry_edge[i] == -1)
                {
                    TBOX_ERROR("NavierStokesBoundaryConditions::NavierStokesBoundaryConditions()\n"
                        << "Node boundary condition type = '"
                        << d_master_bdry_node_conds[i] << "' and \n"
                        << "node location = '" << i
                        << "' passed are inconsistent."
                        << std::endl);
                }
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const boost::shared_ptr<FlowModelBoundaryUtilities> flow_model_boundary_utilities =
            d_flow_model->getFlowModelBoundaryUtilities();
        
        for (int i = 0; i < NUM_3D_FACES; i++)
        {
            d_scalar_bdry_face_conds[i] = d_master_bdry_face_conds[i];
            d_vector_bdry_face_conds[i] = d_master_bdry_face_conds[i];
            
            if (d_master_bdry_face_conds[i] == BDRY_COND::BASIC::REFLECT)
            {
                d_scalar_bdry_face_conds[i] = BDRY_COND::BASIC::SYMMETRY;
            }
        }
        
        for (int i = 0; i < NUM_3D_EDGES; i++)
        {
            d_scalar_bdry_edge_conds[i] = d_master_bdry_edge_conds[i];
            d_vector_bdry_edge_conds[i] = d_master_bdry_edge_conds[i];
            
            if (d_master_bdry_edge_conds[i] == BDRY_COND::BASIC::XREFLECT)
            {
                d_scalar_bdry_edge_conds[i] = BDRY_COND::BASIC::XSYMMETRY;
            }
            if (d_master_bdry_edge_conds[i] == BDRY_COND::BASIC::YREFLECT)
            {
                d_scalar_bdry_edge_conds[i] = BDRY_COND::BASIC::YSYMMETRY;
            }
            if (d_master_bdry_edge_conds[i] == BDRY_COND::BASIC::ZREFLECT)
            {
                d_scalar_bdry_edge_conds[i] = BDRY_COND::BASIC::ZSYMMETRY;
            }
            
            if (d_master_bdry_edge_conds[i] != BOGUS_BDRY_DATA)
            {
                d_edge_bdry_face[i] =
                    flow_model_boundary_utilities->getFaceLocationForEdgeBdry(
                        i, d_master_bdry_edge_conds[i]);
                
                if (d_edge_bdry_face[i] == -1)
                {
                    d_edge_bdry_face[i] =
                        BasicCartesianBoundaryUtilities3::getFaceLocationForEdgeBdry(
                            i, d_master_bdry_edge_conds[i]);
                }
                
                if (d_edge_bdry_face[i] == -1)
                {
                    TBOX_ERROR("NavierStokesBoundaryConditions::NavierStokesBoundaryConditions()\n"
                        << "Edge boundary condition type = '"
                        << d_master_bdry_edge_conds[i] << "' and \n"
                        << "edge location = '" << i
                        << "' passed are inconsistent."
                        << std::endl);
                }
            }
        }
        
        for (int i = 0; i < NUM_3D_NODES; i++)
        {
            d_scalar_bdry_node_conds[i] = d_master_bdry_node_conds[i];
            d_vector_bdry_node_conds[i] = d_master_bdry_node_conds[i];
            
            if (d_master_bdry_node_conds[i] == BDRY_COND::BASIC::XREFLECT)
            {
                d_scalar_bdry_node_conds[i] = BDRY_COND::BASIC::XSYMMETRY;
            }
            if (d_master_bdry_node_conds[i] == BDRY_COND::BASIC::YREFLECT)
            {
                d_scalar_bdry_node_conds[i] = BDRY_COND::BASIC::YSYMMETRY;
            }
            if (d_master_bdry_node_conds[i] == BDRY_COND::BASIC::ZREFLECT)
            {
                d_scalar_bdry_node_conds[i] = BDRY_COND::BASIC::ZSYMMETRY;
            }
            
            if (d_master_bdry_node_conds[i] != BOGUS_BDRY_DATA)
            {
                d_node_bdry_face[i] =
                    flow_model_boundary_utilities->getFaceLocationForNodeBdry(
                        i, d_master_bdry_node_conds[i]);
                
                if (d_node_bdry_face[i] == -1)
                {
                    d_node_bdry_face[i] =
                        BasicCartesianBoundaryUtilities3::getFaceLocationForNodeBdry(
                            i, d_master_bdry_node_conds[i]);
                }
                
                if (d_node_bdry_face[i] == -1)
                {
                    TBOX_ERROR("NavierStokesBoundaryConditions::NavierStokesBoundaryConditions()\n"
                        << "Node boundary condition type = '"
                        << d_master_bdry_node_conds[i] << "' and \n"
                        << "node location = '" << i
                        << "' passed are inconsistent." << std::endl);
                }
            }
        }
    }
    
    d_Navier_Stokes_special_boundary_conditions.reset(new NavierStokesSpecialBoundaryConditions(
        "d_Navier_Stokes_special_boundary_conditions",
        d_project_name,
        d_dim,
        d_grid_geometry,
        d_flow_model_type,
        d_flow_model));
}


/*
 * Print all characteristics of the boundary conditions class.
 */
void
NavierStokesBoundaryConditions::printClassData(std::ostream& os) const
{
    std::vector<boost::shared_ptr<pdat::CellVariable<double> > > conservative_var =
        d_flow_model->getConservativeVariables();
    
    os << "\nPrint NavierStokesBoundaryConditions object..."
       << std::endl;
    
    os << std::endl;
    
    if (d_dim == tbox::Dimension(1))
    {
        for (int j = 0; j < static_cast<int>(d_master_bdry_node_conds.size()); j++)
        {
            os << "d_master_bdry_node_conds["
               << j
               << "] = "
               << static_cast<BDRY_COND::BASIC::TYPE>(d_master_bdry_node_conds[j])
               << std::endl;
            os << "d_scalar_bdry_node_conds["
               << j
               << "] = "
               << static_cast<BDRY_COND::BASIC::TYPE>(d_scalar_bdry_node_conds[j])
               << std::endl;
            os << "d_vector_bdry_node_conds["
               << j
               << "] = "
               << static_cast<BDRY_COND::BASIC::TYPE>(d_vector_bdry_node_conds[j])
               << std::endl;
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        for (int j = 0; j < static_cast<int>(d_master_bdry_node_conds.size()); j++)
        {
             os << "d_master_bdry_node_conds["
                << j
                << "] = "
                << static_cast<BDRY_COND::BASIC::TYPE>(d_master_bdry_node_conds[j])
                << std::endl;
             os << "d_scalar_bdry_node_conds["
                << j
                << "] = "
                << static_cast<BDRY_COND::BASIC::TYPE>(d_scalar_bdry_node_conds[j])
                << std::endl;
             os << "d_vector_bdry_node_conds["
                << j
                << "] = "
                << static_cast<BDRY_COND::BASIC::TYPE>(d_vector_bdry_node_conds[j])
                << std::endl;
             os << "d_node_bdry_edge["
                << j
                << "] = "
                << static_cast<NODE_BDRY_LOC_2D::TYPE>(d_node_bdry_edge[j])
                << std::endl;
        }
        
        os << std::endl;
        
        for (int j = 0; j < static_cast<int>(d_master_bdry_edge_conds.size()); j++)
        {
            os << "d_master_bdry_edge_conds["
               << j
               << "] = "
               << static_cast<BDRY_COND::BASIC::TYPE>(d_master_bdry_edge_conds[j])
               << std::endl;
            
            os << "d_scalar_bdry_edge_conds["
               << j
               << "] = "
               << static_cast<BDRY_COND::BASIC::TYPE>(d_scalar_bdry_edge_conds[j])
               << std::endl;
            
            os << "d_vector_bdry_edge_conds["
               << j
               << "] = "
               << static_cast<BDRY_COND::BASIC::TYPE>(d_vector_bdry_edge_conds[j])
               << std::endl;
            
            if (d_master_bdry_edge_conds[j] == BDRY_COND::BASIC::DIRICHLET)
            {
                for (int vi = 0; vi < static_cast<int>(conservative_var.size()); vi++)
                {
                    const int var_depth = conservative_var[vi]->getDepth();
                    
                    os << "d_bdry_edge_conservative_var["
                       << vi << "]["
                       << j << "] = "
                       << d_bdry_edge_conservative_var[vi][j*var_depth + 0];
                    
                    for (int di = 1; di < var_depth; di++)
                    {
                        os << ", "
                           << d_bdry_edge_conservative_var[vi][j*var_depth + di];
                    }
                    
                    os << std::endl;
                }
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        for (int j = 0; j < static_cast<int>(d_master_bdry_node_conds.size()); j++)
        {
            os << "d_master_bdry_node_conds["
               << j
               << "] = "
               << static_cast<BDRY_COND::BASIC::TYPE>(d_master_bdry_node_conds[j])
               << std::endl;
            os << "d_scalar_bdry_node_conds["
               << j
               << "] = "
               << static_cast<BDRY_COND::BASIC::TYPE>(d_scalar_bdry_node_conds[j])
               << std::endl;
            os << "d_vector_bdry_node_conds["
               << j
               << "] = "
               << static_cast<BDRY_COND::BASIC::TYPE>(d_vector_bdry_node_conds[j])
               << std::endl;
            os << "d_node_bdry_face["
               << j
               << "] = "
               << static_cast<NODE_BDRY_LOC_3D::TYPE>(d_node_bdry_face[j])
               << std::endl;
        }
        
        os << std::endl;
        
        for (int j = 0; j < static_cast<int>(d_master_bdry_edge_conds.size()); j++)
        {
            os << "d_master_bdry_edge_conds["
               << j
               << "] = "
               << static_cast<BDRY_COND::BASIC::TYPE>(d_master_bdry_edge_conds[j])
               << std::endl;
            os << "d_scalar_bdry_edge_conds["
               << j
               << "] = "
               << static_cast<BDRY_COND::BASIC::TYPE>(d_scalar_bdry_edge_conds[j])
               << std::endl;
            os << "d_vector_bdry_edge_conds["
               << j
               << "] = "
               << static_cast<BDRY_COND::BASIC::TYPE>(d_vector_bdry_edge_conds[j])
               << std::endl;
            os << "d_edge_bdry_face["
               << j
               << "] = "
               << static_cast<EDGE_BDRY_LOC_3D::TYPE>(d_edge_bdry_face[j])
               << std::endl;
        }
        
        os << std::endl;
        
        
        for (int j = 0; j < static_cast<int>(d_master_bdry_face_conds.size()); j++)
        {
            os << "d_master_bdry_face_conds["
               << j
               << "] = "
               << static_cast<BDRY_COND::BASIC::TYPE>(d_master_bdry_face_conds[j])
               << std::endl;
            os << "d_scalar_bdry_face_conds["
               << j
               << "] = "
               << static_cast<BDRY_COND::BASIC::TYPE>(d_scalar_bdry_face_conds[j])
               << std::endl;
            os << "d_vector_bdry_face_conds["
               << j
               << "] = "
               << static_cast<BDRY_COND::BASIC::TYPE>(d_vector_bdry_face_conds[j])
               << std::endl;
            
            if (d_master_bdry_face_conds[j] == BDRY_COND::BASIC::DIRICHLET)
            {
                for (int vi = 0; vi < static_cast<int>(conservative_var.size()); vi++)
                {
                    const int var_depth = conservative_var[vi]->getDepth();
                    
                    os << "d_bdry_face_conservative_var["
                       << vi << "]["
                       << j << "] = "
                       << d_bdry_face_conservative_var[vi][j*var_depth + 0];
                    
                    for (int di = 1; di < var_depth; di++)
                    {
                        os << ", "
                           << d_bdry_face_conservative_var[vi][j*var_depth + di];
                    }
                    
                    os << std::endl;
                }
            }
        }
    }
}


/*
 * Put the characteristics of the boundary conditions class into the restart
 * database.
 */
void
NavierStokesBoundaryConditions::putToRestart(
    const boost::shared_ptr<tbox::Database>& restart_db) const
{
    std::vector<boost::shared_ptr<pdat::CellVariable<double> > > conservative_var =
        d_flow_model->getConservativeVariables();
    
    restart_db->putIntegerVector("d_master_bdry_node_conds", d_master_bdry_node_conds);
    
    if (d_dim == tbox::Dimension(1))
    {
        for (int vi = 0; vi < static_cast<int>(conservative_var.size()); vi++)
        {
            restart_db->putDoubleVector(
                "d_bdry_node_conservative_var[" + tbox::Utilities::intToString(vi) + "]",
                d_bdry_node_conservative_var[vi]);
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        restart_db->putIntegerVector("d_master_bdry_edge_conds", d_master_bdry_edge_conds);
        
        for (int vi = 0; vi < static_cast<int>(conservative_var.size()); vi++)
        {
            restart_db->putDoubleVector(
                "d_bdry_edge_conservative_var[" + tbox::Utilities::intToString(vi) + "]",
                d_bdry_edge_conservative_var[vi]);
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        restart_db->putIntegerVector("d_master_bdry_edge_conds",
                                     d_master_bdry_edge_conds);
        
        restart_db->putIntegerVector("d_master_bdry_face_conds",
                                     d_master_bdry_face_conds);
        
        for (int vi = 0; vi < static_cast<int>(conservative_var.size()); vi++)
        {
            restart_db->putDoubleVector(
                "d_bdry_face_conservative_var[" + tbox::Utilities::intToString(vi) + "]",
                d_bdry_face_conservative_var[vi]);
        }
    }
}


/*
 * This routine is a concrete implementation of the virtual function
 * in the base class BoundaryUtilityStrategy.  It reads DIRICHLET
 * boundary state values from the given database with the
 * given name string idenifier.  The integer location index
 * indicates the face (in 3D) or edge (in 2D) to which the boundary
 * condition applies.
 */
void
NavierStokesBoundaryConditions::readDirichletBoundaryDataEntry(
    const boost::shared_ptr<tbox::Database>& db,
    std::string& db_name,
    int bdry_location_index)
{
    TBOX_ASSERT(db);
    TBOX_ASSERT(!db_name.empty());
    
    // Get the primitve data at the boundaries.
    std::vector<double> V = readPrimitiveDataEntry(db, db_name);
    
    // Get a vector pointers to the primitive data.
    std::vector<const double*> V_ptr;
    V_ptr.reserve(V.size());
    for (int ei = 0; ei < static_cast<int>(V.size()); ei++)
    {
        V_ptr.push_back(&V[ei]);
    }
    
    // Create an uninitialized vector of conservative data at the boundaries and get a vector
    // of pointers to the data.
    std::vector<double> Q(V.size());
    std::vector<double*> Q_ptr;
    Q_ptr.reserve(Q.size());
    for (int ei = 0; ei < static_cast<int>(Q.size()); ei++)
    {
        Q_ptr.push_back(&Q[ei]);
    }
    
    // Convert the primitive boundary data to conservative boundary data.
    d_flow_model->setupBasicUtilities();
    boost::shared_ptr<FlowModelBasicUtilities> basic_utilities = d_flow_model->getFlowModelBasicUtilities();
    basic_utilities->convertPrimitiveVariablesToConservativeVariables(V_ptr, Q_ptr);
    
    std::vector<boost::shared_ptr<pdat::CellVariable<double> > > conservative_var =
        d_flow_model->getConservativeVariables();
        
    if (d_dim == tbox::Dimension(1))
    {
        int count = 0;
        for (int vi = 0; vi < static_cast<int>(conservative_var.size()); vi++)
        {
            const int var_depth = conservative_var[vi]->getDepth();
            
            for (int di = 0; di < var_depth; di++)
            {
                d_bdry_node_conservative_var[vi][bdry_location_index*var_depth + di] = Q[count];
                
                count++;
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        int count = 0;
        for (int vi = 0; vi < static_cast<int>(conservative_var.size()); vi++)
        {
            const int var_depth = conservative_var[vi]->getDepth();
            
            for (int di = 0; di < var_depth; di++)
            {
                d_bdry_edge_conservative_var[vi][bdry_location_index*var_depth + di] = Q[count];
                
                count++;
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        int count = 0;
        for (int vi = 0; vi < static_cast<int>(conservative_var.size()); vi++)
        {
            const int var_depth = conservative_var[vi]->getDepth();
            
            for (int di = 0; di < var_depth; di++)
            {
                d_bdry_face_conservative_var[vi][bdry_location_index*var_depth + di] = Q[count];
                
                count++;
            }
        }
    }
}


/*
 * This routine is a concrete implementation of the virtual function
 * in the base class BoundaryUtilityStrategy.  It is a blank implementation
 * for the purposes of this class.
 */
void
NavierStokesBoundaryConditions::readNeumannBoundaryDataEntry(
    const boost::shared_ptr<tbox::Database>& db,
    std::string& db_name,
    int bdry_location_index)
{
    NULL_USE(db);
    NULL_USE(db_name);
    NULL_USE(bdry_location_index);
}


/*
 * Set the data in ghost cells corresponding to physical boundary
 * conditions. Specific boundary conditions are determined by
 * information specified in input file and numerical routines.
 */
void
NavierStokesBoundaryConditions::setPhysicalBoundaryConditions(
    hier::Patch& patch,
    const double fill_time,
    const hier::IntVector& ghost_width_to_fill,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
    NULL_USE(fill_time);
    
    /*
     * Get the conservative variables.
     */
    
    std::vector<boost::shared_ptr<pdat::CellVariable<double> > > conservative_var =
        d_flow_model->getConservativeVariables();
    
    std::vector<std::string> conservative_var_types =
        d_flow_model->getVariableTypesOfConservativeVariables();
    
    std::vector<std::string> conservative_var_names =
        d_flow_model->getNamesOfConservativeVariables();
    
    d_flow_model->registerPatchWithDataContext(patch, data_context);
    
    std::vector<boost::shared_ptr<pdat::CellData<double> > > conservative_var_data =
        d_flow_model->getCellDataOfConservativeVariables();
    
    /*
     * Get the boost::shared_ptr to flow model boundary utilities object from d_flow_model.
     */
    
    const boost::shared_ptr<FlowModelBoundaryUtilities> flow_model_boundary_utilities =
        d_flow_model->getFlowModelBoundaryUtilities();
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Set boundary conditions for cells corresponding to patch nodes.
         */
        
        std::vector<int> node_locs;
        node_locs.reserve(NUM_1D_NODES);
        for (int ni = 0; ni < NUM_1D_NODES; ni++)
        {
            node_locs.push_back(ni);
        }
        
        flow_model_boundary_utilities->
            fill1dNodeBoundaryData(
                conservative_var_data,
                patch,
                node_locs,
                d_vector_bdry_node_conds,
                d_bdry_node_conservative_var,
                ghost_width_to_fill);
        
        for (int vi = 0; vi < static_cast<int>(conservative_var.size()); vi++)
        {
            if (conservative_var_types[vi] == "SCALAR")
            {
                BasicCartesianBoundaryUtilities1::fillNodeBoundaryData(
                    conservative_var_names[vi],
                    conservative_var_data[vi],
                    patch,
                    node_locs,
                    d_scalar_bdry_node_conds,
                    d_bdry_node_conservative_var[vi],
                    ghost_width_to_fill);
            }
            else if (conservative_var_types[vi] == "VECTOR")
            {
                BasicCartesianBoundaryUtilities1::fillNodeBoundaryData(
                    conservative_var_names[vi],
                    conservative_var_data[vi],
                    patch,
                    node_locs,
                    d_vector_bdry_node_conds,
                    d_bdry_node_conservative_var[vi],
                    ghost_width_to_fill);
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Set boundary conditions for cells corresponding to patch edges.
         */
        
        std::vector<int> edge_locs;
        edge_locs.reserve(NUM_2D_EDGES);
        for (int ei = 0; ei < NUM_2D_EDGES; ei++)
        {
            edge_locs.push_back(ei);
        }
        
        flow_model_boundary_utilities->
            fill2dEdgeBoundaryData(
                conservative_var_data,
                patch,
                edge_locs,
                d_vector_bdry_edge_conds,
                d_bdry_edge_conservative_var,
                ghost_width_to_fill);
        
        for (int vi = 0; vi < static_cast<int>(conservative_var.size()); vi++)
        {
            if (conservative_var_types[vi] == "SCALAR")
            {
                BasicCartesianBoundaryUtilities2::fillEdgeBoundaryData(
                    conservative_var_names[vi],
                    conservative_var_data[vi],
                    patch,
                    edge_locs,
                    d_scalar_bdry_edge_conds,
                    d_bdry_edge_conservative_var[vi],
                    ghost_width_to_fill);
            }
            else if (conservative_var_types[vi] == "VECTOR")
            {
                BasicCartesianBoundaryUtilities2::fillEdgeBoundaryData(
                    conservative_var_names[vi],
                    conservative_var_data[vi],
                    patch,
                    edge_locs,
                    d_vector_bdry_edge_conds,
                    d_bdry_edge_conservative_var[vi],
                    ghost_width_to_fill);
            }
        }
        
        /*
         *  Set boundary conditions for cells corresponding to patch nodes.
         */
        
        std::vector<int> node_locs;
        node_locs.reserve(NUM_2D_NODES);
        for (int ni = 0; ni < NUM_2D_NODES; ni++)
        {
            node_locs.push_back(ni);
        }
        
        flow_model_boundary_utilities->
            fill2dNodeBoundaryData(
                conservative_var_data,
                patch,
                node_locs,
                d_vector_bdry_node_conds,
                d_bdry_edge_conservative_var,
                ghost_width_to_fill);
        
        for (int vi = 0; vi < static_cast<int>(conservative_var.size()); vi++)
        {
            if (conservative_var_types[vi] == "SCALAR")
            {
                BasicCartesianBoundaryUtilities2::fillNodeBoundaryData(
                    conservative_var_names[vi],
                    conservative_var_data[vi],
                    patch,
                    node_locs,
                    d_scalar_bdry_node_conds,
                    d_bdry_edge_conservative_var[vi],
                    ghost_width_to_fill);
            }
            else if (conservative_var_types[vi] == "VECTOR")
            {
                BasicCartesianBoundaryUtilities2::fillNodeBoundaryData(
                    conservative_var_names[vi],
                    conservative_var_data[vi],
                    patch,
                    node_locs,
                    d_vector_bdry_node_conds,
                    d_bdry_edge_conservative_var[vi],
                    ghost_width_to_fill);
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         *  Set boundary conditions for cells corresponding to patch faces.
         */
        
        std::vector<int> face_locs;
        face_locs.reserve(NUM_3D_FACES);
        for (int fi = 0; fi < NUM_3D_FACES; fi++)
        {
            face_locs.push_back(fi);
        }
        
        flow_model_boundary_utilities->
            fill3dFaceBoundaryData(
                conservative_var_data,
                patch,
                face_locs,
                d_vector_bdry_face_conds,
                d_bdry_face_conservative_var,
                ghost_width_to_fill);
        
        for (int vi = 0; vi < static_cast<int>(conservative_var.size()); vi++)
        {
            if (conservative_var_types[vi] == "SCALAR")
            {
                BasicCartesianBoundaryUtilities3::fillFaceBoundaryData(
                    conservative_var_names[vi],
                    conservative_var_data[vi],
                    patch,
                    face_locs,
                    d_scalar_bdry_face_conds,
                    d_bdry_face_conservative_var[vi],
                    ghost_width_to_fill);
            }
            else if (conservative_var_types[vi] == "VECTOR")
            {
                BasicCartesianBoundaryUtilities3::fillFaceBoundaryData(
                    conservative_var_names[vi],
                    conservative_var_data[vi],
                    patch,
                    face_locs,
                    d_vector_bdry_face_conds,
                    d_bdry_face_conservative_var[vi],
                    ghost_width_to_fill);
            }
        }
        
        /*
         * Set boundary conditions for cells corresponding to patch edges.
         */
        
        std::vector<int> edge_locs;
        edge_locs.reserve(NUM_3D_EDGES);
        for (int ei = 0; ei < NUM_3D_EDGES; ei++)
        {
            edge_locs.push_back(ei);
        }
        
        flow_model_boundary_utilities->
            fill3dEdgeBoundaryData(
                conservative_var_data,
                patch,
                edge_locs,
                d_vector_bdry_edge_conds,
                d_bdry_face_conservative_var,
                ghost_width_to_fill);
        
        for (int vi = 0; vi < static_cast<int>(conservative_var.size()); vi++)
        {
            if (conservative_var_types[vi] == "SCALAR")
            {
                BasicCartesianBoundaryUtilities3::fillEdgeBoundaryData(
                    conservative_var_names[vi],
                    conservative_var_data[vi],
                    patch,
                    edge_locs,
                    d_scalar_bdry_edge_conds,
                    d_bdry_face_conservative_var[vi],
                    ghost_width_to_fill);
            }
            else if (conservative_var_types[vi] == "VECTOR")
            {
                BasicCartesianBoundaryUtilities3::fillEdgeBoundaryData(
                    conservative_var_names[vi],
                    conservative_var_data[vi],
                    patch,
                    edge_locs,
                    d_vector_bdry_edge_conds,
                    d_bdry_face_conservative_var[vi],
                    ghost_width_to_fill);
            }
        }
        
        /*
         *  Set boundary conditions for cells corresponding to patch nodes.
         */
        
        std::vector<int> node_locs;
        node_locs.reserve(NUM_3D_NODES);
        for (int ni = 0; ni < NUM_3D_NODES; ni++)
        {
            node_locs.push_back(ni);
        }
        
        flow_model_boundary_utilities->
            fill3dNodeBoundaryData(
                conservative_var_data,
                patch,
                node_locs,
                d_vector_bdry_node_conds,
                d_bdry_face_conservative_var,
                ghost_width_to_fill);
        
        for (int vi = 0; vi < static_cast<int>(conservative_var.size()); vi++)
        {
            if (conservative_var_types[vi] == "SCALAR")
            {
                BasicCartesianBoundaryUtilities3::fillNodeBoundaryData(
                    conservative_var_names[vi],
                    conservative_var_data[vi],
                    patch,
                    node_locs,
                    d_scalar_bdry_node_conds,
                    d_bdry_face_conservative_var[vi],
                    ghost_width_to_fill);
            }
            else if (conservative_var_types[vi] == "VECTOR")
            {
                BasicCartesianBoundaryUtilities3::fillNodeBoundaryData(
                    conservative_var_names[vi],
                    conservative_var_data[vi],
                    patch,
                    node_locs,
                    d_vector_bdry_node_conds,
                    d_bdry_face_conservative_var[vi],
                    ghost_width_to_fill);
            }
        }
    }

    d_Navier_Stokes_special_boundary_conditions->setSpecialBoundaryConditions(
        patch,
        conservative_var_data,
        fill_time,
        ghost_width_to_fill);
    
    d_flow_model->unregisterPatch();
}


std::vector<double>
NavierStokesBoundaryConditions::readPrimitiveDataEntry(
    boost::shared_ptr<tbox::Database> db,
    const std::string& db_name)
{
    TBOX_ASSERT(db);
    TBOX_ASSERT(!db_name.empty());
    
    std::vector<double> data_primitive_var;
    
    std::vector<std::string> primitive_var_names = d_flow_model->getNamesOfPrimitiveVariables(true);
    
    for (int vi = 0; vi < static_cast<int>(primitive_var_names.size()); vi++)
    {
        if (db->keyExists(primitive_var_names[vi]))
        {
            std::vector<double> vector_primitive_var_data = db->getDoubleVector(primitive_var_names[vi]);
            
            const int var_depth = static_cast<int>(vector_primitive_var_data.size());
            
            for (int di = 0; di < var_depth; di++)
            {
                data_primitive_var.push_back(vector_primitive_var_data[di]);
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": '"
                << primitive_var_names[vi]
                << "' entry missing from '"
                << db_name
                << "' input database."
                << std::endl);
        }
    }
    
    return data_primitive_var;
}


/*
 * Set defaults for boundary conditions. Set to bogus values
 * for error checking.
 */
void
NavierStokesBoundaryConditions::setDefaultBoundaryConditions()
{
    std::vector<boost::shared_ptr<pdat::CellVariable<double> > > conservative_var =
        d_flow_model->getConservativeVariables();
    
    if (d_dim == tbox::Dimension(1))
    {
        d_bdry_node_conservative_var.resize(conservative_var.size());
        
        d_master_bdry_node_conds.resize(NUM_1D_NODES);
        d_scalar_bdry_node_conds.resize(NUM_1D_NODES);
        d_vector_bdry_node_conds.resize(NUM_1D_NODES);
        for (int ni = 0; ni < NUM_1D_NODES; ni++)
        {
            d_master_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
            d_scalar_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
            d_vector_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
        }
        
        for (int vi = 0; vi < static_cast<int>(conservative_var.size()); vi++)
        {
            d_bdry_node_conservative_var[vi].resize(NUM_1D_NODES*conservative_var[vi]->getDepth());
            
            tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_node_conservative_var[vi]);
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        d_bdry_edge_conservative_var.resize(conservative_var.size());
        
        d_master_bdry_edge_conds.resize(NUM_2D_EDGES);
        d_scalar_bdry_edge_conds.resize(NUM_2D_EDGES);
        d_vector_bdry_edge_conds.resize(NUM_2D_EDGES);
        for (int ei = 0; ei < NUM_2D_EDGES; ei++)
        {
            d_master_bdry_edge_conds[ei] = BOGUS_BDRY_DATA;
            d_scalar_bdry_edge_conds[ei] = BOGUS_BDRY_DATA;
            d_vector_bdry_edge_conds[ei] = BOGUS_BDRY_DATA;
        }
        
        d_master_bdry_node_conds.resize(NUM_2D_NODES);
        d_scalar_bdry_node_conds.resize(NUM_2D_NODES);
        d_vector_bdry_node_conds.resize(NUM_2D_NODES);
        d_node_bdry_edge.resize(NUM_2D_NODES);
        for (int ni = 0; ni < NUM_2D_NODES; ni++)
        {
            d_master_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
            d_scalar_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
            d_vector_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
            d_node_bdry_edge[ni] = BOGUS_BDRY_DATA;
        }
        
        for (int vi = 0; vi < static_cast<int>(conservative_var.size()); vi++)
        {
            d_bdry_edge_conservative_var[vi].resize(NUM_2D_EDGES*conservative_var[vi]->getDepth());
            
            tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_edge_conservative_var[vi]);
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        d_bdry_face_conservative_var.resize(conservative_var.size());
        
        d_master_bdry_face_conds.resize(NUM_3D_FACES);
        d_scalar_bdry_face_conds.resize(NUM_3D_FACES);
        d_vector_bdry_face_conds.resize(NUM_3D_FACES);
        for (int fi = 0; fi < NUM_3D_FACES; fi++) {
           d_master_bdry_face_conds[fi] = BOGUS_BDRY_DATA;
           d_scalar_bdry_face_conds[fi] = BOGUS_BDRY_DATA;
           d_vector_bdry_face_conds[fi] = BOGUS_BDRY_DATA;
        }
     
        d_master_bdry_edge_conds.resize(NUM_3D_EDGES);
        d_scalar_bdry_edge_conds.resize(NUM_3D_EDGES);
        d_vector_bdry_edge_conds.resize(NUM_3D_EDGES);
        d_edge_bdry_face.resize(NUM_3D_EDGES);
        for (int ei = 0; ei < NUM_3D_EDGES; ei++)
        {
            d_master_bdry_edge_conds[ei] = BOGUS_BDRY_DATA;
            d_scalar_bdry_edge_conds[ei] = BOGUS_BDRY_DATA;
            d_vector_bdry_edge_conds[ei] = BOGUS_BDRY_DATA;
            d_edge_bdry_face[ei] = BOGUS_BDRY_DATA;
        }
     
        d_master_bdry_node_conds.resize(NUM_3D_NODES);
        d_scalar_bdry_node_conds.resize(NUM_3D_NODES);
        d_vector_bdry_node_conds.resize(NUM_3D_NODES);
        d_node_bdry_face.resize(NUM_3D_NODES);
        for (int ni = 0; ni < NUM_3D_NODES; ni++)
        {
            d_master_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
            d_scalar_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
            d_vector_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
            d_node_bdry_face[ni] = BOGUS_BDRY_DATA;
        }
        
        for (int vi = 0; vi < static_cast<int>(conservative_var.size()); vi++)
        {
            d_bdry_face_conservative_var[vi].resize(NUM_3D_FACES*conservative_var[vi]->getDepth());
            
            tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_face_conservative_var[vi]);
        }
    }
}
