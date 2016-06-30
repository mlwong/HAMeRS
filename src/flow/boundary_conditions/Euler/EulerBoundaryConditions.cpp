#include "flow/boundary_conditions/Euler/EulerBoundaryConditions.hpp"

//integer constant for debugging improperly set boundary data
#define BOGUS_BDRY_DATA (-9999)

// routines for managing boundary data
#include "util/basic_boundary_conditions/CartesianBoundaryUtilities2.hpp"
#include "util/basic_boundary_conditions/CartesianBoundaryUtilities3.hpp"

EulerBoundaryConditions::EulerBoundaryConditions(
    const std::string& object_name,
    const std::string& project_name,
    const tbox::Dimension& dim,
    const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
    const int& num_species,
    const boost::shared_ptr<FlowModel>& flow_model,
    const boost::shared_ptr<tbox::Database>& boundary_conditions_db,
    const bool& boundary_conditions_db_is_from_restart):
        d_object_name(object_name),
        d_project_name(project_name),
        d_dim(dim),
        d_grid_geometry(grid_geometry),
        d_num_species(num_species),
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
    
    if (num_per_dirs < d_dim.getValue())
    {
        if (boundary_conditions_db_is_from_restart)
        {
            std::vector<boost::shared_ptr<pdat::CellVariable<double> > > conservative_var =
                d_flow_model->getConservativeVariables();
            
            d_master_bdry_node_conds = boundary_conditions_db->getIntegerVector("d_master_bdry_node_conds");
            
            if (d_dim == tbox::Dimension(1))
            {
                // NOT YET IMPLEMENTED
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
                    d_bdry_edge_conservative_var[vi] = boundary_conditions_db->getDoubleVector(
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
                // NOT YET IMPLEMENTED
            }
            if (d_dim == tbox::Dimension(2))
            {
                CartesianBoundaryUtilities2::getFromInput(
                    this,
                    boundary_conditions_db,
                    d_master_bdry_edge_conds,
                    d_master_bdry_node_conds,
                    periodic);
            }
            else if (d_dim == tbox::Dimension(3))
            {
                CartesianBoundaryUtilities3::getFromInput(
                    this,
                    boundary_conditions_db,
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
        // NOT YET IMPLEMENTED
    }
    else if (d_dim == tbox::Dimension(2))
    {
        for (int i = 0; i < NUM_2D_EDGES; i++)
        {
            d_scalar_bdry_edge_conds[i] = d_master_bdry_edge_conds[i];
            d_vector_bdry_edge_conds[i] = d_master_bdry_edge_conds[i];
            
            if (d_master_bdry_edge_conds[i] == BdryCond::Basic::REFLECT)
            {
                d_scalar_bdry_edge_conds[i] = BdryCond::Basic::SYMMETRY;
            }
        }
        
        for (int i = 0; i < NUM_2D_NODES; i++)
        {
            d_scalar_bdry_node_conds[i] = d_master_bdry_node_conds[i];
            d_vector_bdry_node_conds[i] = d_master_bdry_node_conds[i];
            
            if (d_master_bdry_node_conds[i] == BdryCond::Basic::XREFLECT)
            {
                d_scalar_bdry_node_conds[i] = BdryCond::Basic::XSYMMETRY;
            }
            if (d_master_bdry_node_conds[i] == BdryCond::Basic::YREFLECT)
            {
                d_scalar_bdry_node_conds[i] = BdryCond::Basic::YSYMMETRY;
            }
            
            if (d_master_bdry_node_conds[i] != BOGUS_BDRY_DATA)
            {
                d_node_bdry_edge[i] =
                    CartesianBoundaryUtilities2::getEdgeLocationForNodeBdry(
                        i, d_master_bdry_node_conds[i]);
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        for (int i = 0; i < NUM_3D_FACES; i++)
        {
            d_scalar_bdry_face_conds[i] = d_master_bdry_face_conds[i];
            d_vector_bdry_face_conds[i] = d_master_bdry_face_conds[i];
            
            if (d_master_bdry_face_conds[i] == BdryCond::Basic::REFLECT)
            {
                d_scalar_bdry_face_conds[i] = BdryCond::Basic::SYMMETRY;
            }
        }
        
        for (int i = 0; i < NUM_3D_EDGES; i++)
        {
            d_scalar_bdry_edge_conds[i] = d_master_bdry_edge_conds[i];
            d_vector_bdry_edge_conds[i] = d_master_bdry_edge_conds[i];
            
            if (d_master_bdry_edge_conds[i] == BdryCond::Basic::XREFLECT)
            {
                d_scalar_bdry_edge_conds[i] = BdryCond::Basic::XSYMMETRY;
            }
            if (d_master_bdry_edge_conds[i] == BdryCond::Basic::YREFLECT)
            {
                d_scalar_bdry_edge_conds[i] = BdryCond::Basic::YSYMMETRY;
            }
            if (d_master_bdry_edge_conds[i] == BdryCond::Basic::ZREFLECT)
            {
                d_scalar_bdry_edge_conds[i] = BdryCond::Basic::ZSYMMETRY;
            }
            
            if (d_master_bdry_edge_conds[i] != BOGUS_BDRY_DATA)
            {
                d_edge_bdry_face[i] =
                    CartesianBoundaryUtilities3::getFaceLocationForEdgeBdry(
                        i, d_master_bdry_edge_conds[i]);
            }
        }
        
        for (int i = 0; i < NUM_3D_NODES; i++)
        {
            d_scalar_bdry_node_conds[i] = d_master_bdry_node_conds[i];
            d_vector_bdry_node_conds[i] = d_master_bdry_node_conds[i];
            
            if (d_master_bdry_node_conds[i] == BdryCond::Basic::XREFLECT)
            {
                d_scalar_bdry_node_conds[i] = BdryCond::Basic::XSYMMETRY;
            }
            if (d_master_bdry_node_conds[i] == BdryCond::Basic::YREFLECT)
            {
                d_scalar_bdry_node_conds[i] = BdryCond::Basic::YSYMMETRY;
            }
            if (d_master_bdry_node_conds[i] == BdryCond::Basic::ZREFLECT)
            {
                d_scalar_bdry_node_conds[i] = BdryCond::Basic::ZSYMMETRY;
            }
            
            if (d_master_bdry_node_conds[i] != BOGUS_BDRY_DATA)
            {
                d_node_bdry_face[i] =
                    CartesianBoundaryUtilities3::getFaceLocationForNodeBdry(
                        i, d_master_bdry_node_conds[i]);
            }
        }
    }
}


/*
 * Print all characteristics of the boundary conditions class.
 */
void
EulerBoundaryConditions::printClassData(std::ostream& os) const
{
    std::vector<boost::shared_ptr<pdat::CellVariable<double> > > conservative_var =
        d_flow_model->getConservativeVariables();
    
    os << "\nPrint EulerBoundaryConditions object..."
       << std::endl;
    
    os << std::endl;
    
    if (d_dim == tbox::Dimension(1))
    {
        // NOT YET IMPLEMENTED
    }
    else if (d_dim == tbox::Dimension(2))
    {
        for (int j = 0; j < static_cast<int>(d_master_bdry_node_conds.size()); j++)
        {
             os << "d_master_bdry_node_conds["
                << j
                << "] = "
                << static_cast<BdryCond::Basic::Type>(d_master_bdry_node_conds[j])
                << std::endl;
             os << "d_scalar_bdry_node_conds["
                << j
                << "] = "
                << static_cast<BdryCond::Basic::Type>(d_scalar_bdry_node_conds[j])
                << std::endl;
             os << "d_vector_bdry_node_conds["
                << j
                << "] = "
                << static_cast<BdryCond::Basic::Type>(d_vector_bdry_node_conds[j])
                << std::endl;
             os << "d_node_bdry_edge["
                << j
                << "] = "
                << static_cast<NodeBdyLoc2D::Type>(d_node_bdry_edge[j])
                << std::endl;
        }
        
        os << std::endl;
        
        for (int j = 0; j < static_cast<int>(d_master_bdry_edge_conds.size()); j++)
        {
            os << "d_master_bdry_edge_conds["
               << j
               << "] = "
               << static_cast<BdryCond::Basic::Type>(d_master_bdry_edge_conds[j])
               << std::endl;
            
            os << "d_scalar_bdry_edge_conds["
               << j
               << "] = "
               << static_cast<BdryCond::Basic::Type>(d_scalar_bdry_edge_conds[j])
               << std::endl;
            
            os << "d_vector_bdry_edge_conds["
               << j
               << "] = "
               << static_cast<BdryCond::Basic::Type>(d_vector_bdry_edge_conds[j])
               << std::endl;
            
            if (d_master_bdry_edge_conds[j] == BdryCond::Basic::DIRICHLET)
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
               << static_cast<BdryCond::Basic::Type>(d_master_bdry_node_conds[j])
               << std::endl;
            os << "d_scalar_bdry_node_conds["
               << j
               << "] = "
               << static_cast<BdryCond::Basic::Type>(d_scalar_bdry_node_conds[j])
               << std::endl;
            os << "d_vector_bdry_node_conds["
               << j
               << "] = "
               << static_cast<BdryCond::Basic::Type>(d_vector_bdry_node_conds[j])
               << std::endl;
            os << "d_node_bdry_face["
               << j
               << "] = "
               << static_cast<NodeBdyLoc3D::Type>(d_node_bdry_face[j])
               << std::endl;
        }
        
        os << std::endl;
        
        for (int j = 0; j < static_cast<int>(d_master_bdry_edge_conds.size()); j++)
        {
            os << "d_master_bdry_edge_conds["
               << j
               << "] = "
               << static_cast<BdryCond::Basic::Type>(d_master_bdry_edge_conds[j])
               << std::endl;
            os << "d_scalar_bdry_edge_conds["
               << j
               << "] = "
               << static_cast<BdryCond::Basic::Type>(d_scalar_bdry_edge_conds[j])
               << std::endl;
            os << "d_vector_bdry_edge_conds["
               << j
               << "] = "
               << static_cast<BdryCond::Basic::Type>(d_vector_bdry_edge_conds[j])
               << std::endl;
            os << "d_edge_bdry_face["
               << j
               << "] = "
               << static_cast<EdgeBdyLoc3D::Type>(d_edge_bdry_face[j])
               << std::endl;
        }
        
        os << std::endl;
        
        
        for (int j = 0; j < static_cast<int>(d_master_bdry_face_conds.size()); j++)
        {
            os << "d_master_bdry_face_conds["
               << j
               << "] = "
               << static_cast<BdryCond::Basic::Type>(d_master_bdry_face_conds[j])
               << std::endl;
            os << "d_scalar_bdry_face_conds["
               << j
               << "] = "
               << static_cast<BdryCond::Basic::Type>(d_scalar_bdry_face_conds[j])
               << std::endl;
            os << "d_vector_bdry_face_conds["
               << j
               << "] = "
               << static_cast<BdryCond::Basic::Type>(d_vector_bdry_face_conds[j])
               << std::endl;
            
            if (d_master_bdry_face_conds[j] == BdryCond::Basic::DIRICHLET)
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
EulerBoundaryConditions::putToRestart(
    const boost::shared_ptr<tbox::Database>& restart_db) const
{
    std::vector<boost::shared_ptr<pdat::CellVariable<double> > > conservative_var =
        d_flow_model->getConservativeVariables();
    
    restart_db->putIntegerVector("d_master_bdry_node_conds", d_master_bdry_node_conds);
    
    if (d_dim == tbox::Dimension(1))
    {
        // NOT YET IMPLEMENTED
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
EulerBoundaryConditions::readDirichletBoundaryDataEntry(
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
    d_flow_model->convertLocalCellDataPointersPrimitiveVariablesToConservativeVariables(V_ptr, Q_ptr);
    
    std::vector<boost::shared_ptr<pdat::CellVariable<double> > > conservative_var =
        d_flow_model->getConservativeVariables();
        
    if (d_dim == tbox::Dimension(1))
    {
        // Not YET IMPLEMENTED
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
EulerBoundaryConditions::readNeumannBoundaryDataEntry(
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
EulerBoundaryConditions::setPhysicalBoundaryConditions(
    hier::Patch& patch,
    const double fill_time,
    const hier::IntVector& ghost_width_to_fill,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
    NULL_USE(fill_time);
    
    std::vector<boost::shared_ptr<pdat::CellVariable<double> > > conservative_var =
        d_flow_model->getConservativeVariables();
    
    std::vector<std::string> conservative_var_types =
        d_flow_model->getVariableTypesOfConservativeVariables();
    
    std::vector<std::string> conservative_var_names =
        d_flow_model->getNamesOfConservativeVariables();
    
    d_flow_model->registerPatchWithDataContext(patch, data_context);
    
    std::vector<boost::shared_ptr<pdat::CellData<double> > > conservative_var_data =
        d_flow_model->getGlobalCellDataConservativeVariables();
    
    if (d_dim == tbox::Dimension(1))
    {
        // NOT YET IMPLEMENTED
    }
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Set boundary conditions for cells corresponding to patch edges.
         */
        
        for (int vi = 0; vi < static_cast<int>(conservative_var.size()); vi++)
        {
            if (conservative_var_types[vi] == "SCALAR")
            {
                /*
                appu::CartesianBoundaryUtilities2::fillEdgeBoundaryData(
                    conservative_var_names[vi],
                    conservative_var_data[vi],
                    patch,
                    ghost_width_to_fill,
                    d_scalar_bdry_edge_conds,
                    d_bdry_edge_conservative_var[vi]);
                */
                
                CartesianBoundaryUtilities2::fillEdgeBoundaryData(
                    conservative_var_names[vi],
                    conservative_var_data[vi],
                    patch,
                    d_scalar_bdry_edge_conds,
                    d_bdry_edge_conservative_var[vi],
                    ghost_width_to_fill);
            }
            else if (conservative_var_types[vi] == "VECTOR")
            {
                /*
                appu::CartesianBoundaryUtilities2::fillEdgeBoundaryData(
                    conservative_var_names[vi],
                    conservative_var_data[vi],
                    patch,
                    ghost_width_to_fill,
                    d_vector_bdry_edge_conds,
                    d_bdry_edge_conservative_var[vi]);
                */
                
                CartesianBoundaryUtilities2::fillEdgeBoundaryData(
                    conservative_var_names[vi],
                    conservative_var_data[vi],
                    patch,
                    d_vector_bdry_edge_conds,
                    d_bdry_edge_conservative_var[vi],
                    ghost_width_to_fill);
            }
        }
        
        /*
         *  Set boundary conditions for cells corresponding to patch nodes.
         */
        
        for (int vi = 0; vi < static_cast<int>(conservative_var.size()); vi++)
        {
            if (conservative_var_types[vi] == "SCALAR")
            {
                /*
                appu::CartesianBoundaryUtilities2::fillNodeBoundaryData(
                    conservative_var_names[vi],
                    conservative_var_data[vi],
                    patch,
                    ghost_width_to_fill,
                    d_scalar_bdry_node_conds,
                    d_bdry_edge_conservative_var[vi]);
                */
                
                CartesianBoundaryUtilities2::fillNodeBoundaryData(
                    conservative_var_names[vi],
                    conservative_var_data[vi],
                    patch,
                    d_scalar_bdry_node_conds,
                    d_bdry_edge_conservative_var[vi],
                    ghost_width_to_fill);
            }
            else if (conservative_var_types[vi] == "VECTOR")
            {
                /*
                appu::CartesianBoundaryUtilities2::fillNodeBoundaryData(
                    conservative_var_names[vi],
                    conservative_var_data[vi],
                    patch,
                    ghost_width_to_fill,
                    d_vector_bdry_node_conds,
                    d_bdry_edge_conservative_var[vi]);
                */
                
                CartesianBoundaryUtilities2::fillNodeBoundaryData(
                    conservative_var_names[vi],
                    conservative_var_data[vi],
                    patch,
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
        
        for (int vi = 0; vi < static_cast<int>(conservative_var.size()); vi++)
        {
            if (conservative_var_types[vi] == "SCALAR")
            {
                /*
                appu::CartesianBoundaryUtilities3::fillFaceBoundaryData(
                    conservative_var_names[vi],
                    conservative_var_data[vi],
                    patch,
                    ghost_width_to_fill,
                    d_scalar_bdry_face_conds,
                    d_bdry_face_conservative_var[vi]);
                */
                
                CartesianBoundaryUtilities3::fillFaceBoundaryData(
                    conservative_var_names[vi],
                    conservative_var_data[vi],
                    patch,
                    d_scalar_bdry_face_conds,
                    d_bdry_face_conservative_var[vi],
                    ghost_width_to_fill);
            }
            else if (conservative_var_types[vi] == "VECTOR")
            {
                /*
                appu::CartesianBoundaryUtilities3::fillFaceBoundaryData(
                    conservative_var_names[vi],
                    conservative_var_data[vi],
                    patch,
                    ghost_width_to_fill,
                    d_vector_bdry_face_conds,
                    d_bdry_face_conservative_var[vi]);
                */
                
                CartesianBoundaryUtilities3::fillFaceBoundaryData(
                    conservative_var_names[vi],
                    conservative_var_data[vi],
                    patch,
                    d_vector_bdry_face_conds,
                    d_bdry_face_conservative_var[vi],
                    ghost_width_to_fill);
            }
        }
        
        /*
         * Set boundary conditions for cells corresponding to patch edges.
         */
        
        for (int vi = 0; vi < static_cast<int>(conservative_var.size()); vi++)
        {
            if (conservative_var_types[vi] == "SCALAR")
            {
                /*
                appu::CartesianBoundaryUtilities3::fillEdgeBoundaryData(
                    conservative_var_names[vi],
                    conservative_var_data[vi],
                    patch,
                    ghost_width_to_fill,
                    d_scalar_bdry_edge_conds,
                    d_bdry_face_conservative_var[vi]);
                */
                
                CartesianBoundaryUtilities3::fillEdgeBoundaryData(
                    conservative_var_names[vi],
                    conservative_var_data[vi],
                    patch,
                    d_scalar_bdry_edge_conds,
                    d_bdry_face_conservative_var[vi],
                    ghost_width_to_fill);
            }
            else if (conservative_var_types[vi] == "VECTOR")
            {
                /*
                appu::CartesianBoundaryUtilities3::fillEdgeBoundaryData(
                    conservative_var_names[vi],
                    conservative_var_data[vi],
                    patch,
                    ghost_width_to_fill,
                    d_vector_bdry_edge_conds,
                    d_bdry_face_conservative_var[vi]);
                */
                
                CartesianBoundaryUtilities3::fillEdgeBoundaryData(
                    conservative_var_names[vi],
                    conservative_var_data[vi],
                    patch,
                    d_vector_bdry_edge_conds,
                    d_bdry_face_conservative_var[vi],
                    ghost_width_to_fill);
            }
        }
        
        /*
         *  Set boundary conditions for cells corresponding to patch nodes.
         */
        
        for (int vi = 0; vi < static_cast<int>(conservative_var.size()); vi++)
        {
            if (conservative_var_types[vi] == "SCALAR")
            {
                /*
                appu::CartesianBoundaryUtilities3::fillNodeBoundaryData(
                    conservative_var_names[vi],
                    conservative_var_data[vi],
                    patch,
                    ghost_width_to_fill,
                    d_scalar_bdry_node_conds,
                    d_bdry_face_conservative_var[vi]);
                */
                
                CartesianBoundaryUtilities3::fillNodeBoundaryData(
                    conservative_var_names[vi],
                    conservative_var_data[vi],
                    patch,
                    d_scalar_bdry_node_conds,
                    d_bdry_face_conservative_var[vi],
                    ghost_width_to_fill);
            }
            else if (conservative_var_types[vi] == "VECTOR")
            {
                /*
                appu::CartesianBoundaryUtilities3::fillNodeBoundaryData(
                    conservative_var_names[vi],
                    conservative_var_data[vi],
                    patch,
                    ghost_width_to_fill,
                    d_vector_bdry_node_conds,
                    d_bdry_face_conservative_var[vi]);
                */
                
                CartesianBoundaryUtilities3::fillNodeBoundaryData(
                    conservative_var_names[vi],
                    conservative_var_data[vi],
                    patch,
                    d_vector_bdry_node_conds,
                    d_bdry_face_conservative_var[vi],
                    ghost_width_to_fill);
            }
        }
    }
    
    /*
     * Set the boundary conditions for specific problems.
     */
    if (d_project_name == "2D double-Mach reflection")
    {
        TBOX_ASSERT(d_dim == tbox::Dimension(2));
        
        // Get the number of ghost cells of gradient.
        hier::IntVector num_ghosts = conservative_var_data[0]->getGhostCellWidth();
        
        // Get the dimensions of box that covers the interior of patch.
        const hier::Box interior_box = patch.getBox();
        const hier::IntVector interior_dims = interior_box.numberCells();
        
        // Get the dimensions of box that covers interior of patch plus
        // ghost cells.
        const hier::Box ghost_box = conservative_var_data[0]->getGhostBox();
        const hier::IntVector ghostcell_dims = ghost_box.numberCells();
        
        const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
            BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                patch.getPatchGeometry()));
        
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(patch_geom);
#endif
        
        if (patch_geom->getTouchesRegularBoundary(1, 0) ||
            patch_geom->getTouchesRegularBoundary(1, 1))
        {
            const double* const dx = patch_geom->getDx();
            const double* const patch_xlo = patch_geom->getXLower();
            
            boost::shared_ptr<pdat::CellData<double> > density = conservative_var_data[0];
            boost::shared_ptr<pdat::CellData<double> > momentum = conservative_var_data[1];
            boost::shared_ptr<pdat::CellData<double> > total_energy = conservative_var_data[2];
            
            double* rho   = density->getPointer(0);
            double* rho_u = momentum->getPointer(0);
            double* rho_v = momentum->getPointer(1);
            double* E     = total_energy->getPointer(0);
            
            const double x_0 = 1.0/6.0;
            
            const double gamma = 1.4;
            
            const double rho_post_shock = 8.0;
            const double u_post_shock = 8.25*cos(M_PI/6.0);
            const double v_post_shock = -8.25*sin(M_PI/6.0);
            const double p_post_shock = 116.5;
            
            const double rho_pre_shock = 1.4;
            const double u_pre_shock = 0.0;
            const double v_pre_shock = 0.0;
            const double p_pre_shock = 1.0;
            
            const double rho_u_post_shock = rho_post_shock*u_post_shock;
            const double rho_v_post_shock = rho_post_shock*v_post_shock;
            
            const double rho_u_pre_shock = rho_pre_shock*u_pre_shock;
            const double rho_v_pre_shock = rho_pre_shock*v_pre_shock;
            
            const double E_pre_shock = p_pre_shock/(gamma - 1.0) +
                0.5*rho_pre_shock*(u_pre_shock*u_pre_shock + v_pre_shock*v_pre_shock);
            
            const double E_post_shock = p_post_shock/(gamma - 1.0) +
                0.5*rho_post_shock*(u_post_shock*u_post_shock + v_post_shock*v_post_shock);
            
            /*
             * Update the bottom boundary conditions.
             */
            
            if (patch_geom->getTouchesRegularBoundary(1, 0))
            {
                for (int i = 0; i < interior_dims[0]; i++)
                {
                    for (int j = -ghost_width_to_fill[1];
                         j < 0;
                         j++)
                    {
                        const int idx_cell = (i + num_ghosts[0]) +
                            (j + num_ghosts[1])*ghostcell_dims[0];
                        
                        // Compute the coordinates.
                        double x[2];
                        x[0] = patch_xlo[0] + (i + 0.5)*dx[0];
                        x[1] = patch_xlo[1] + (j + 0.5)*dx[1];
                        
                        if (x[0] < x_0)
                        {
                            rho[idx_cell] = rho_post_shock;
                            rho_u[idx_cell] = rho_u_post_shock;
                            rho_v[idx_cell] = rho_v_post_shock;
                            E[idx_cell] = E_post_shock;
                        }
                        else
                        {
                            const int idx_mirror_cell = (i + num_ghosts[0]) +
                                (-j + num_ghosts[1] - 1)*ghostcell_dims[0];
                            
                            rho[idx_cell] = rho[idx_mirror_cell];
                            rho_u[idx_cell] = rho_u[idx_mirror_cell];
                            rho_v[idx_cell] = -rho_v[idx_mirror_cell];
                            E[idx_cell] = E[idx_mirror_cell];
                        }
                    }
                }
            }
            
            /*
             * Update the top boundary conditions.
             */
            
            if (patch_geom->getTouchesRegularBoundary(1, 1))
            {
                const double x_s = x_0 + (1.0 + 20.0*fill_time)/sqrt(3.0);
                
                for (int i = 0; i < interior_dims[0]; i++)
                {
                    for (int j = interior_dims[1];
                         j < interior_dims[1] + ghost_width_to_fill[1];
                         j++)
                    {
                        const int idx_cell = (i + num_ghosts[0]) +
                            (j + num_ghosts[1])*ghostcell_dims[0];
                        
                        // Compute the coordinates.
                        double x[2];
                        x[0] = patch_xlo[0] + (i + 0.5)*dx[0];
                        x[1] = patch_xlo[1] + (j + 0.5)*dx[1];
                        
                        if (x[0] >= x_s)
                        {
                            rho[idx_cell] = rho_pre_shock;
                            rho_u[idx_cell] = rho_u_pre_shock;
                            rho_v[idx_cell] = rho_v_pre_shock;
                            E[idx_cell] = E_pre_shock;
                        }
                        else
                        {
                            rho[idx_cell] = rho_post_shock;
                            rho_u[idx_cell] = rho_u_post_shock;
                            rho_v[idx_cell] = rho_v_post_shock;
                            E[idx_cell] = E_post_shock;
                        }
                    }
                }
            }
        }
    }
    
    d_flow_model->unregisterPatch();
}


std::vector<double>
EulerBoundaryConditions::readPrimitiveDataEntry(
    boost::shared_ptr<tbox::Database> db,
    const std::string& db_name)
{
    TBOX_ASSERT(db);
    TBOX_ASSERT(!db_name.empty());
    
    std::vector<double> data_primitive_var;
    
    std::vector<std::string> primitive_var_names = d_flow_model->getNamesOfPrimitiveVariables();
    
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
EulerBoundaryConditions::setDefaultBoundaryConditions()
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
