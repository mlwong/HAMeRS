#include "flow_model/boundary_conditions/Euler/EulerBoundaryConditions.hpp"

//integer constant for debugging improperly set boundary data
#define BOGUS_BDRY_DATA (-9999)

// routines for managing boundary data
#include "SAMRAI/appu/CartesianBoundaryUtilities2.h"
#include "SAMRAI/appu/CartesianBoundaryUtilities3.h"

EulerBoundaryConditions::EulerBoundaryConditions(
    const std::string& object_name,
    const std::string& project_name,
    const tbox::Dimension& dim,
    const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
    const FLOW_MODEL& flow_model,
    const int& num_species,
    const boost::shared_ptr<EquationOfState>& equation_of_state,
    const boost::shared_ptr<tbox::Database>& boundary_conditions_db,
    const bool& boundary_conditions_db_is_from_restart):
        d_object_name(object_name),
        d_project_name(project_name),
        d_dim(dim),
        d_grid_geometry(grid_geometry),
        d_num_ghosts(hier::IntVector::getZero(d_dim)),
        d_flow_model(flow_model),
        d_num_species(num_species),
        d_equation_of_state(equation_of_state),
        d_variables_set(false),
        d_num_ghosts_set(false)
{
    /*
     * Defaults for boundary conditions. Set to bogus values
     * for error checking.
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
            d_master_bdry_node_conds = boundary_conditions_db->getIntegerVector("d_master_bdry_node_conds");
            
            if (d_dim == tbox::Dimension(1))
            {
                // NOT YET IMPLEMENTED.
                
                /*
                switch (d_flow_model)
                {
                    case SINGLE_SPECIES:
                    {
                        d_bdry_node_density = boundary_conditions_db->getDoubleVector("d_bdry_node_density");
                        d_bdry_node_momentum = boundary_conditions_db->getDoubleVector("d_bdry_node_momentum");
                        d_bdry_node_total_energy = boundary_conditions_db->getDoubleVector("d_bdry_node_total_energy");
                        
                        break;
                    }
                    case FOUR_EQN_SHYUE:
                    {
                        d_bdry_node_density = boundary_conditions_db->getDoubleVector("d_bdry_node_density");
                        d_bdry_node_momentum = boundary_conditions_db->getDoubleVector("d_bdry_node_momentum");
                        d_bdry_node_total_energy = boundary_conditions_db->getDoubleVector("d_bdry_node_total_energy");
                        d_bdry_node_mass_fraction = boundary_conditions_db->getDoubleVector("d_bdry_node_mass_fraction");
                        
                        break;
                    }
                    case FIVE_EQN_ALLAIRE:
                    {
                        d_bdry_node_partial_density = boundary_conditions_db->getDoubleVector("d_bdry_node_partial_density");
                        d_bdry_node_momentum = boundary_conditions_db->getDoubleVector("d_bdry_node_momentum");
                        d_bdry_node_total_energy = boundary_conditions_db->getDoubleVector("d_bdry_node_total_energy");
                        d_bdry_node_volume_fraction = boundary_conditions_db->getDoubleVector("d_bdry_node_volume_fraction");
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": "
                            << "d_flow_model '"
                            << d_flow_model
                            << "' not yet implemented."
                            << std::endl);
                    }
                }
                */
            }
            else if (d_dim == tbox::Dimension(2))
            {
                d_master_bdry_edge_conds = boundary_conditions_db->getIntegerVector("d_master_bdry_edge_conds");
                
                switch (d_flow_model)
                {
                    case SINGLE_SPECIES:
                    {
                        d_bdry_edge_density = boundary_conditions_db->getDoubleVector("d_bdry_edge_density");
                        d_bdry_edge_momentum = boundary_conditions_db->getDoubleVector("d_bdry_edge_momentum");
                        d_bdry_edge_total_energy = boundary_conditions_db->getDoubleVector("d_bdry_edge_total_energy");
                        
                        break;
                    }
                    case FOUR_EQN_SHYUE:
                    {
                        d_bdry_edge_density = boundary_conditions_db->getDoubleVector("d_bdry_edge_density");
                        d_bdry_edge_momentum = boundary_conditions_db->getDoubleVector("d_bdry_edge_momentum");
                        d_bdry_edge_total_energy = boundary_conditions_db->getDoubleVector("d_bdry_edge_total_energy");
                        d_bdry_edge_mass_fraction = boundary_conditions_db->getDoubleVector("d_bdry_edge_mass_fraction");
                        
                        break;
                    }
                    case FIVE_EQN_ALLAIRE:
                    {
                        d_bdry_edge_partial_density = boundary_conditions_db->getDoubleVector("d_bdry_edge_partial_density");
                        d_bdry_edge_momentum = boundary_conditions_db->getDoubleVector("d_bdry_edge_momentum");
                        d_bdry_edge_total_energy = boundary_conditions_db->getDoubleVector("d_bdry_edge_total_energy");
                        d_bdry_edge_volume_fraction = boundary_conditions_db->getDoubleVector("d_bdry_edge_volume_fraction");
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": "
                            << "d_flow_model '"
                            << d_flow_model
                            << "' not yet implemented."
                            << std::endl);
                    }
                }
            }
            else if (d_dim == tbox::Dimension(3))
            {
                d_master_bdry_edge_conds = boundary_conditions_db->getIntegerVector("d_master_bdry_edge_conds");
                d_master_bdry_face_conds = boundary_conditions_db->getIntegerVector("d_master_bdry_face_conds");
                
                switch (d_flow_model)
                {
                    case SINGLE_SPECIES:
                    {
                        d_bdry_face_density = boundary_conditions_db->getDoubleVector("d_bdry_face_density");
                        d_bdry_face_momentum = boundary_conditions_db->getDoubleVector("d_bdry_face_momentum");
                        d_bdry_face_total_energy = boundary_conditions_db->getDoubleVector("d_bdry_face_total_energy");
                        
                        break;
                    }
                    case FOUR_EQN_SHYUE:
                    {
                        d_bdry_face_density = boundary_conditions_db->getDoubleVector("d_bdry_face_density");
                        d_bdry_face_momentum = boundary_conditions_db->getDoubleVector("d_bdry_face_momentum");
                        d_bdry_face_total_energy = boundary_conditions_db->getDoubleVector("d_bdry_face_total_energy");
                        d_bdry_face_mass_fraction = boundary_conditions_db->getDoubleVector("d_bdry_face_mass_fraction");
                        
                        break;
                    }
                    case FIVE_EQN_ALLAIRE:
                    {
                        d_bdry_face_partial_density = boundary_conditions_db->getDoubleVector("d_bdry_face_partial_density");
                        d_bdry_face_momentum = boundary_conditions_db->getDoubleVector("d_bdry_face_momentum");
                        d_bdry_face_total_energy = boundary_conditions_db->getDoubleVector("d_bdry_face_total_energy");
                        d_bdry_face_volume_fraction = boundary_conditions_db->getDoubleVector("d_bdry_face_volume_fraction");
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": "
                            << "d_flow_model '"
                            << d_flow_model
                            << "' not yet implemented."
                            << std::endl);
                    }
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
                appu::CartesianBoundaryUtilities2::getFromInput(
                    this,
                    boundary_conditions_db,
                    d_master_bdry_edge_conds,
                    d_master_bdry_node_conds,
                    periodic);
            }
            else if (d_dim == tbox::Dimension(3))
            {
                appu::CartesianBoundaryUtilities3::getFromInput(
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
            
            if (d_master_bdry_edge_conds[i] == BdryCond::REFLECT)
            {
                d_scalar_bdry_edge_conds[i] = BdryCond::FLOW;
            }
        }
        
        for (int i = 0; i < NUM_2D_NODES; i++)
        {
            d_scalar_bdry_node_conds[i] = d_master_bdry_node_conds[i];
            d_vector_bdry_node_conds[i] = d_master_bdry_node_conds[i];
            
            if (d_master_bdry_node_conds[i] == BdryCond::XREFLECT)
            {
                d_scalar_bdry_node_conds[i] = BdryCond::XFLOW;
            }
            if (d_master_bdry_node_conds[i] == BdryCond::YREFLECT)
            {
                d_scalar_bdry_node_conds[i] = BdryCond::YFLOW;
            }
            
            if (d_master_bdry_node_conds[i] != BOGUS_BDRY_DATA)
            {
                d_node_bdry_edge[i] =
                    appu::CartesianBoundaryUtilities2::getEdgeLocationForNodeBdry(
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
            
            if (d_master_bdry_face_conds[i] == BdryCond::REFLECT)
            {
                d_scalar_bdry_face_conds[i] = BdryCond::FLOW;
            }
        }
        
        for (int i = 0; i < NUM_3D_EDGES; i++)
        {
            d_scalar_bdry_edge_conds[i] = d_master_bdry_edge_conds[i];
            d_vector_bdry_edge_conds[i] = d_master_bdry_edge_conds[i];
            
            if (d_master_bdry_edge_conds[i] == BdryCond::XREFLECT)
            {
                d_scalar_bdry_edge_conds[i] = BdryCond::XFLOW;
            }
            if (d_master_bdry_edge_conds[i] == BdryCond::YREFLECT)
            {
                d_scalar_bdry_edge_conds[i] = BdryCond::YFLOW;
            }
            if (d_master_bdry_edge_conds[i] == BdryCond::ZREFLECT)
            {
                d_scalar_bdry_edge_conds[i] = BdryCond::ZFLOW;
            }
            
            if (d_master_bdry_edge_conds[i] != BOGUS_BDRY_DATA)
            {
                d_edge_bdry_face[i] =
                    appu::CartesianBoundaryUtilities3::getFaceLocationForEdgeBdry(
                        i, d_master_bdry_edge_conds[i]);
            }
        }
        
        for (int i = 0; i < NUM_3D_NODES; i++)
        {
            d_scalar_bdry_node_conds[i] = d_master_bdry_node_conds[i];
            d_vector_bdry_node_conds[i] = d_master_bdry_node_conds[i];
            
            if (d_master_bdry_node_conds[i] == BdryCond::XREFLECT)
            {
                d_scalar_bdry_node_conds[i] = BdryCond::XFLOW;
            }
            if (d_master_bdry_node_conds[i] == BdryCond::YREFLECT)
            {
                d_scalar_bdry_node_conds[i] = BdryCond::YFLOW;
            }
            if (d_master_bdry_node_conds[i] == BdryCond::ZREFLECT)
            {
                d_scalar_bdry_node_conds[i] = BdryCond::ZFLOW;
            }
            
            if (d_master_bdry_node_conds[i] != BOGUS_BDRY_DATA)
            {
                d_node_bdry_face[i] =
                    appu::CartesianBoundaryUtilities3::getFaceLocationForNodeBdry(
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
    os << "\nPrint EulerBoundaryConditions object..."
       << std::endl;
    
    os << std::endl;
    
    os << "d_variables_set = "
       << d_variables_set
       << std::endl;
    os << "d_num_ghosts_set = "
       << d_num_ghosts_set
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
                << static_cast<BdryCond::Type>(d_master_bdry_node_conds[j])
                << std::endl;
             os << "d_scalar_bdry_node_conds["
                << j
                << "] = "
                << static_cast<BdryCond::Type>(d_scalar_bdry_node_conds[j])
                << std::endl;
             os << "d_vector_bdry_node_conds["
                << j
                << "] = "
                << static_cast<BdryCond::Type>(d_vector_bdry_node_conds[j])
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
               << static_cast<BdryCond::Type>(d_master_bdry_edge_conds[j])
               << std::endl;
            
            os << "d_scalar_bdry_edge_conds["
               << j
               << "] = "
               << static_cast<BdryCond::Type>(d_scalar_bdry_edge_conds[j])
               << std::endl;
            
            os << "d_vector_bdry_edge_conds["
               << j
               << "] = "
               << static_cast<BdryCond::Type>(d_vector_bdry_edge_conds[j])
               << std::endl;
            
            if (d_master_bdry_edge_conds[j] == BdryCond::DIRICHLET)
            {
                switch (d_flow_model)
                {
                    case SINGLE_SPECIES:
                    {
                        os << "d_bdry_edge_density["
                           << j
                           << "] = "
                           << d_bdry_edge_density[j]
                           << std::endl;
                        
                        os << "d_bdry_edge_momentum["
                           << j << "] = "
                           << d_bdry_edge_momentum[j*d_dim.getValue() + 0]
                           << " , "
                           << d_bdry_edge_momentum[j*d_dim.getValue() + 1]
                           << std::endl;
                        
                        os << "d_bdry_edge_total_energy["
                           << j
                           << "] = "
                           << d_bdry_edge_total_energy[j]
                           << std::endl;
                        
                        break;
                    }
                    case FOUR_EQN_SHYUE:
                    {
                        os << "d_bdry_edge_density["
                           << j
                           << "] = "
                           << d_bdry_edge_density[j]
                           << std::endl;
                        
                        os << "d_bdry_edge_momentum["
                           << j
                           << "] = "
                           << d_bdry_edge_momentum[j*d_dim.getValue() + 0]
                           << " , "
                           << d_bdry_edge_momentum[j*d_dim.getValue() + 1]
                           << std::endl;
                        
                        os << "d_bdry_edge_total_energy["
                           << j
                           << "] = "
                           << d_bdry_edge_total_energy[j]
                           << std::endl;
                        
                        os << "d_bdry_edge_mass_fraction["
                           << j
                           << "] = "
                           << d_bdry_edge_mass_fraction[j*d_num_species + 0];
                        for (int si = 1; si < d_num_species - 1; si++)
                        {
                            os << " , "
                               << d_bdry_edge_mass_fraction[j*d_num_species + si];
                        }
                        os << std::endl;
                        
                        break;
                    }
                    case FIVE_EQN_ALLAIRE:
                    {
                        os << "d_bdry_edge_partial_density["
                           << j
                           << "] = "
                           << d_bdry_edge_partial_density[j*d_num_species + 0];
                        for (int si = 1; si < d_num_species; si++)
                        {
                            os << " , "
                               << d_bdry_edge_partial_density[j*d_num_species + si];
                        }
                        os << std::endl;
                        
                        os << "d_bdry_edge_momentum["
                           << j << "] = "
                           << d_bdry_edge_momentum[j*d_dim.getValue() + 0]
                           << " , "
                           << d_bdry_edge_momentum[j*d_dim.getValue() + 1]
                           << std::endl;
                        
                        os << "d_bdry_edge_total_energy["
                           << j
                           << "] = "
                           << d_bdry_edge_total_energy[j]
                           << std::endl;
                        
                        os << "d_bdry_edge_volume_fraction["
                           << j
                           << "] = "
                           << d_bdry_edge_volume_fraction[j*d_num_species + 0];
                        for (int si = 1; si < d_num_species - 1; si++)
                        {
                            os << " , "
                               << d_bdry_edge_volume_fraction[j*d_num_species + si];
                        }
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": "
                            << "d_flow_model '"
                            << d_flow_model
                            << "' not yet implemented."
                            << std::endl);
                    }
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
               << static_cast<BdryCond::Type>(d_master_bdry_node_conds[j])
               << std::endl;
            os << "d_scalar_bdry_node_conds["
               << j
               << "] = "
               << static_cast<BdryCond::Type>(d_scalar_bdry_node_conds[j])
               << std::endl;
            os << "d_vector_bdry_node_conds["
               << j
               << "] = "
               << static_cast<BdryCond::Type>(d_vector_bdry_node_conds[j])
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
               << static_cast<BdryCond::Type>(d_master_bdry_edge_conds[j])
               << std::endl;
            os << "d_scalar_bdry_edge_conds["
               << j
               << "] = "
               << static_cast<BdryCond::Type>(d_scalar_bdry_edge_conds[j])
               << std::endl;
            os << "d_vector_bdry_edge_conds["
               << j
               << "] = "
               << static_cast<BdryCond::Type>(d_vector_bdry_edge_conds[j])
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
               << static_cast<BdryCond::Type>(d_master_bdry_face_conds[j])
               << std::endl;
            os << "d_scalar_bdry_face_conds["
               << j
               << "] = "
               << static_cast<BdryCond::Type>(d_scalar_bdry_face_conds[j])
               << std::endl;
            os << "d_vector_bdry_face_conds["
               << j
               << "] = "
               << static_cast<BdryCond::Type>(d_vector_bdry_face_conds[j])
               << std::endl;
            
            if (d_master_bdry_face_conds[j] == BdryCond::DIRICHLET)
            {
                switch (d_flow_model)
                {
                    case SINGLE_SPECIES:
                    {
                        os << "d_bdry_face_density["
                           << j
                           << "] = "
                           << d_bdry_face_density[j]
                           << std::endl;
                        
                        os << "d_bdry_face_momentum["
                           << j << "] = "
                           << d_bdry_face_momentum[j*d_dim.getValue() + 0]
                           << " , "
                           << d_bdry_face_momentum[j*d_dim.getValue() + 1]
                           << " , "
                           << d_bdry_face_momentum[j*d_dim.getValue() + 2]
                           << std::endl;
                        
                        os << "d_bdry_face_total_energy["
                           << j
                           << "] = "
                           << d_bdry_face_total_energy[j]
                           << std::endl;
                        
                        break;
                    }
                    case FOUR_EQN_SHYUE:
                    {
                        os << "d_bdry_face_density["
                           << j
                           << "] = "
                           << d_bdry_face_density[j]
                           << std::endl;
                        
                        os << "d_bdry_face_momentum["
                           << j
                           << "] = "
                           << d_bdry_face_momentum[j*d_dim.getValue() + 0]
                           << " , "
                           << d_bdry_face_momentum[j*d_dim.getValue() + 1]
                           << " , "
                           << d_bdry_face_momentum[j*d_dim.getValue() + 2]
                           << std::endl;
                        
                        os << "d_bdry_face_total_energy["
                           << j
                           << "] = "
                           << d_bdry_face_total_energy[j]
                           << std::endl;
                        
                        os << "d_bdry_face_mass_fraction["
                           << j
                           << "] = "
                           << d_bdry_face_mass_fraction[j*d_num_species + 0];
                        for (int si = 1; si < d_num_species - 1; si++)
                        {
                            os << " , "
                               << d_bdry_face_mass_fraction[j*d_num_species + si];
                        }
                        os << std::endl;
                        
                        break;
                    }
                    case FIVE_EQN_ALLAIRE:
                    {
                        os << "d_bdry_face_partial_density["
                           << j
                           << "] = "
                           << d_bdry_face_partial_density[j*d_num_species + 0];
                        for (int si = 1; si < d_num_species; si++)
                        {
                            os << " , "
                               << d_bdry_face_partial_density[j*d_num_species + si];
                        }
                        os << std::endl;
                        
                        os << "d_bdry_face_momentum["
                           << j << "] = "
                           << d_bdry_face_momentum[j*d_dim.getValue() + 0]
                           << " , "
                           << d_bdry_face_momentum[j*d_dim.getValue() + 1]
                           << " , "
                           << d_bdry_face_momentum[j*d_dim.getValue() + 2]
                           << std::endl;
                        
                        os << "d_bdry_face_total_energy["
                           << j
                           << "] = "
                           << d_bdry_face_total_energy[j]
                           << std::endl;
                        
                        os << "d_bdry_face_volume_fraction["
                           << j
                           << "] = "
                           << d_bdry_face_volume_fraction[j*d_num_species + 0];
                        for (int si = 1; si < d_num_species - 1; si++)
                        {
                            os << " , "
                               << d_bdry_face_volume_fraction[j*d_num_species + si];
                        }
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": "
                            << "d_flow_model '"
                            << d_flow_model
                            << "' not yet implemented."
                            << std::endl);
                    }
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
    restart_db->putIntegerVector("d_master_bdry_node_conds",
        d_master_bdry_node_conds);
    
    if (d_dim == tbox::Dimension(1))
    {
        switch (d_flow_model)
        {
            case SINGLE_SPECIES:
            {
                restart_db->putDoubleVector("d_bdry_node_density",
                                            d_bdry_node_density);
                
                restart_db->putDoubleVector("d_bdry_node_momentum",
                                            d_bdry_node_momentum);
                
                restart_db->putDoubleVector("d_bdry_node_total_energy",
                                            d_bdry_node_total_energy);
                
                break;
            }
            case FOUR_EQN_SHYUE:
            {
                restart_db->putDoubleVector("d_bdry_node_density",
                                            d_bdry_node_density);
                
                restart_db->putDoubleVector("d_bdry_node_momentum",
                                            d_bdry_node_momentum);
                
                restart_db->putDoubleVector("d_bdry_node_total_energy",
                                            d_bdry_node_total_energy);
                
                restart_db->putDoubleVector("d_bdry_node_mass_fraction",
                                            d_bdry_node_mass_fraction);
                
                break;
            }
            case FIVE_EQN_ALLAIRE:
            {
                restart_db->putDoubleVector("d_bdry_node_partial_density",
                                            d_bdry_node_partial_density);
                
                restart_db->putDoubleVector("d_bdry_node_momentum",
                                            d_bdry_node_momentum);
                
                restart_db->putDoubleVector("d_bdry_node_total_energy",
                                            d_bdry_node_total_energy);
                
                restart_db->putDoubleVector("d_bdry_node_volume_fraction",
                                            d_bdry_node_volume_fraction);
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "d_flow_model '"
                    << d_flow_model
                    << "' not yet implemented."
                    << std::endl);
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        restart_db->putIntegerVector("d_master_bdry_edge_conds",
                                     d_master_bdry_edge_conds);
        
        switch (d_flow_model)
        {
            case SINGLE_SPECIES:
            {
                restart_db->putDoubleVector("d_bdry_edge_density",
                                            d_bdry_edge_density);
                
                restart_db->putDoubleVector("d_bdry_edge_momentum",
                                            d_bdry_edge_momentum);
                
                restart_db->putDoubleVector("d_bdry_edge_total_energy",
                                            d_bdry_edge_total_energy);
                
                break;
            }
            case FOUR_EQN_SHYUE:
            {
                restart_db->putDoubleVector("d_bdry_edge_density",
                                            d_bdry_edge_density);
                
                restart_db->putDoubleVector("d_bdry_edge_momentum",
                                            d_bdry_edge_momentum);
                
                restart_db->putDoubleVector("d_bdry_edge_total_energy",
                                            d_bdry_edge_total_energy);
                
                restart_db->putDoubleVector("d_bdry_edge_mass_fraction",
                                            d_bdry_edge_mass_fraction);
                
                break;
            }
            case FIVE_EQN_ALLAIRE:
            {
                restart_db->putDoubleVector("d_bdry_edge_partial_density",
                                            d_bdry_edge_partial_density);
                
                restart_db->putDoubleVector("d_bdry_edge_momentum",
                                            d_bdry_edge_momentum);
                
                restart_db->putDoubleVector("d_bdry_edge_total_energy",
                                            d_bdry_edge_total_energy);
                
                restart_db->putDoubleVector("d_bdry_edge_volume_fraction",
                                            d_bdry_edge_volume_fraction);
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "d_flow_model '"
                    << d_flow_model
                    << "' not yet implemented."
                    << std::endl);
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        restart_db->putIntegerVector("d_master_bdry_edge_conds",
                                     d_master_bdry_edge_conds);
        
        restart_db->putIntegerVector("d_master_bdry_face_conds",
                                     d_master_bdry_face_conds);
        
        switch (d_flow_model)
        {
            case SINGLE_SPECIES:
            {
                restart_db->putDoubleVector("d_bdry_face_density",
                                            d_bdry_face_density);
                
                restart_db->putDoubleVector("d_bdry_face_momentum",
                                            d_bdry_face_momentum);
                
                restart_db->putDoubleVector("d_bdry_face_total_energy",
                                            d_bdry_face_total_energy);
                
                break;
            }
            case FOUR_EQN_SHYUE:
            {
                restart_db->putDoubleVector("d_bdry_face_density",
                                            d_bdry_face_density);
                
                restart_db->putDoubleVector("d_bdry_face_momentum",
                                            d_bdry_face_momentum);
                
                restart_db->putDoubleVector("d_bdry_face_total_energy",
                                            d_bdry_face_total_energy);
                
                restart_db->putDoubleVector("d_bdry_face_mass_fraction",
                                            d_bdry_face_mass_fraction);
                
                break;
            }
            case FIVE_EQN_ALLAIRE:
            {
                restart_db->putDoubleVector("d_bdry_face_partial_density",
                                            d_bdry_face_partial_density);
                
                restart_db->putDoubleVector("d_bdry_face_momentum",
                                            d_bdry_face_momentum);
                
                restart_db->putDoubleVector("d_bdry_face_total_energy",
                                            d_bdry_face_total_energy);
                
                restart_db->putDoubleVector("d_bdry_face_volume_fraction",
                                            d_bdry_face_volume_fraction);
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "d_flow_model '"
                    << d_flow_model
                    << "' not yet implemented."
                    << std::endl);
            }
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
    
    if (d_dim == tbox::Dimension(1))
    {
        // NOT YET IMPLEMENTED
    }
    else if (d_dim == tbox::Dimension(2))
    {
        switch (d_flow_model)
        {
            case SINGLE_SPECIES:
            {
                readStateDataEntryForSingleSpecies(
                    db,
                    db_name,
                    bdry_location_index,
                    d_bdry_edge_density,
                    d_bdry_edge_momentum,
                    d_bdry_edge_total_energy);
                
                break;
            }
            case FOUR_EQN_SHYUE:
            {
                readStateDataEntryForFourEqnShyue(
                    db,
                    db_name,
                    bdry_location_index,
                    d_bdry_edge_density,
                    d_bdry_edge_momentum,
                    d_bdry_edge_total_energy,
                    d_bdry_edge_mass_fraction);
                
                break;
            }
            case FIVE_EQN_ALLAIRE:
            {
                readStateDataEntryForFiveEqnAllaire(
                    db,
                    db_name,
                    bdry_location_index,
                    d_bdry_edge_partial_density,
                    d_bdry_edge_momentum,
                    d_bdry_edge_total_energy,
                    d_bdry_edge_volume_fraction);
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "d_flow_model '"
                    << d_flow_model
                    << "' not yet implemented."
                    << std::endl);
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        switch (d_flow_model)
        {
            case SINGLE_SPECIES:
            {
                readStateDataEntryForSingleSpecies(
                    db,
                    db_name,
                    bdry_location_index,
                    d_bdry_face_density,
                    d_bdry_face_momentum,
                    d_bdry_face_total_energy);
                
                break;
            }
            case FOUR_EQN_SHYUE:
            {
                readStateDataEntryForFourEqnShyue(
                    db,
                    db_name,
                    bdry_location_index,
                    d_bdry_face_density,
                    d_bdry_face_momentum,
                    d_bdry_face_total_energy,
                    d_bdry_face_mass_fraction);
                
                break;
            }
            case FIVE_EQN_ALLAIRE:
            {
                readStateDataEntryForFiveEqnAllaire(
                    db,
                    db_name,
                    bdry_location_index,
                    d_bdry_face_partial_density,
                    d_bdry_face_momentum,
                    d_bdry_face_total_energy,
                    d_bdry_face_volume_fraction);
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "d_flow_model '"
                    << d_flow_model
                    << "' not yet implemented."
                    << std::endl);
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
    
    if (d_variables_set == true)
    {
        if (d_num_ghosts_set == true)
        {
            switch (d_flow_model)
            {
                case SINGLE_SPECIES:
                {
                    boost::shared_ptr<pdat::CellData<double> > density(
                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                            patch.getPatchData(d_density, data_context)));
                    
                    boost::shared_ptr<pdat::CellData<double> > momentum(
                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                            patch.getPatchData(d_momentum, data_context)));
                    
                    boost::shared_ptr<pdat::CellData<double> > total_energy(
                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                            patch.getPatchData(d_total_energy, data_context)));
                    
#ifdef DEBUG_CHECK_ASSERTIONS
                    TBOX_ASSERT(density);
                    TBOX_ASSERT(momentum);
                    TBOX_ASSERT(total_energy);
                    
                    TBOX_ASSERT(density->getGhostCellWidth() == d_num_ghosts);
                    TBOX_ASSERT(momentum->getGhostCellWidth() == d_num_ghosts);
                    TBOX_ASSERT(total_energy->getGhostCellWidth() == d_num_ghosts);
#endif
                    
                    if (d_dim == tbox::Dimension(1))
                    {
                        // NOT YET IMPLEMENTED
                    }
                    else if (d_dim == tbox::Dimension(2))
                    {
                        /*
                         * Set boundary conditions for cells corresponding to patch edges.
                         */
                        
                        appu::CartesianBoundaryUtilities2::fillEdgeBoundaryData(
                            "density",
                            density,
                            patch,
                            ghost_width_to_fill,
                            d_scalar_bdry_edge_conds,
                            d_bdry_edge_density);
                        
                        appu::CartesianBoundaryUtilities2::fillEdgeBoundaryData(
                            "momentum",
                            momentum,
                            patch,
                            ghost_width_to_fill,
                            d_vector_bdry_edge_conds,
                            d_bdry_edge_momentum);
                        
                        appu::CartesianBoundaryUtilities2::fillEdgeBoundaryData(
                            "total energy",
                            total_energy,
                            patch,
                            ghost_width_to_fill,
                            d_scalar_bdry_edge_conds,
                            d_bdry_edge_total_energy);
                        
                        /*
                         *  Set boundary conditions for cells corresponding to patch nodes.
                         */
                        
                        appu::CartesianBoundaryUtilities2::fillNodeBoundaryData(
                            "density",
                            density,
                            patch,
                            ghost_width_to_fill,
                            d_scalar_bdry_node_conds,
                            d_bdry_edge_density);
                        
                        appu::CartesianBoundaryUtilities2::fillNodeBoundaryData(
                            "momentum",
                            momentum,
                            patch,
                            ghost_width_to_fill,
                            d_vector_bdry_node_conds,
                            d_bdry_edge_momentum);
                        
                        appu::CartesianBoundaryUtilities2::fillNodeBoundaryData(
                            "total energy",
                            total_energy,
                            patch,
                            ghost_width_to_fill,
                            d_scalar_bdry_node_conds,
                            d_bdry_edge_total_energy);
                    }
                    else if (d_dim == tbox::Dimension(3))
                    {
                        /*
                         *  Set boundary conditions for cells corresponding to patch faces.
                         */
                        
                        appu::CartesianBoundaryUtilities3::fillFaceBoundaryData(
                            "density",
                            density,
                            patch,
                            ghost_width_to_fill,
                            d_scalar_bdry_face_conds,
                            d_bdry_face_density);
                        
                        appu::CartesianBoundaryUtilities3::fillFaceBoundaryData(
                            "momentum",
                            momentum,
                            patch,
                            ghost_width_to_fill,
                            d_vector_bdry_face_conds,
                            d_bdry_face_momentum);
                        
                        appu::CartesianBoundaryUtilities3::fillFaceBoundaryData(
                            "total energy",
                            total_energy,
                            patch,
                            ghost_width_to_fill,
                            d_scalar_bdry_face_conds,
                            d_bdry_face_total_energy);
                        
                        /*
                         *  Set boundary conditions for cells corresponding to patch edges.
                         */
                        
                        appu::CartesianBoundaryUtilities3::fillEdgeBoundaryData(
                            "density",
                            density,
                            patch,
                            ghost_width_to_fill,
                            d_scalar_bdry_edge_conds,
                            d_bdry_face_density);
                        
                        appu::CartesianBoundaryUtilities3::fillEdgeBoundaryData(
                            "momentum",
                            momentum,
                            patch,
                            ghost_width_to_fill,
                            d_vector_bdry_edge_conds,
                            d_bdry_face_momentum);
                        
                        appu::CartesianBoundaryUtilities3::fillEdgeBoundaryData(
                            "total energy",
                            total_energy,
                            patch,
                            ghost_width_to_fill,
                            d_scalar_bdry_edge_conds,
                            d_bdry_face_total_energy);
                        
                        /*
                         *  Set boundary conditions for cells corresponding to patch nodes.
                         */
                        
                        appu::CartesianBoundaryUtilities3::fillNodeBoundaryData(
                            "density",
                            density,
                            patch,
                            ghost_width_to_fill,
                            d_scalar_bdry_node_conds,
                            d_bdry_face_density);
                        
                        appu::CartesianBoundaryUtilities3::fillNodeBoundaryData(
                            "momentum",
                            momentum,
                            patch,
                            ghost_width_to_fill,
                            d_vector_bdry_node_conds,
                            d_bdry_face_momentum);
                        
                        appu::CartesianBoundaryUtilities3::fillNodeBoundaryData(
                            "total energy",
                            total_energy,
                            patch,
                            ghost_width_to_fill,
                            d_scalar_bdry_node_conds,
                            d_bdry_face_total_energy);
                    }
                    
                    break;
                }
                case FOUR_EQN_SHYUE:
                {
                    boost::shared_ptr<pdat::CellData<double> > density(
                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                            patch.getPatchData(d_density, data_context)));
                    
                    boost::shared_ptr<pdat::CellData<double> > momentum(
                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                            patch.getPatchData(d_momentum, data_context)));
                    
                    boost::shared_ptr<pdat::CellData<double> > total_energy(
                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                            patch.getPatchData(d_total_energy, data_context)));
                    
                    boost::shared_ptr<pdat::CellData<double> > mass_fraction(
                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                            patch.getPatchData(d_mass_fraction, data_context)));
                    
#ifdef DEBUG_CHECK_ASSERTIONS
                    TBOX_ASSERT(density);
                    TBOX_ASSERT(momentum);
                    TBOX_ASSERT(total_energy);
                    TBOX_ASSERT(mass_fraction);
                    
                    TBOX_ASSERT(density->getGhostCellWidth() == d_num_ghosts);
                    TBOX_ASSERT(momentum->getGhostCellWidth() == d_num_ghosts);
                    TBOX_ASSERT(total_energy->getGhostCellWidth() == d_num_ghosts);
                    TBOX_ASSERT(mass_fraction->getGhostCellWidth() == d_num_ghosts);
#endif
                    
                    if (d_dim == tbox::Dimension(1))
                    {
                        // NOT YET IMPLEMENTED
                    }
                    else if (d_dim == tbox::Dimension(2))
                    {
                        /*
                         * Set boundary conditions for cells corresponding to patch edges.
                         */
                        
                        appu::CartesianBoundaryUtilities2::fillEdgeBoundaryData(
                            "density",
                            density,
                            patch,
                            ghost_width_to_fill,
                            d_scalar_bdry_edge_conds,
                            d_bdry_edge_density);
                        
                        appu::CartesianBoundaryUtilities2::fillEdgeBoundaryData(
                            "momentum",
                            momentum,
                            patch,
                            ghost_width_to_fill,
                            d_vector_bdry_edge_conds,
                            d_bdry_edge_momentum);
                        
                        appu::CartesianBoundaryUtilities2::fillEdgeBoundaryData(
                            "total energy",
                            total_energy,
                            patch,
                            ghost_width_to_fill,
                            d_scalar_bdry_edge_conds,
                            d_bdry_edge_total_energy);
                        
                        appu::CartesianBoundaryUtilities2::fillEdgeBoundaryData(
                            "mass fraction",
                            mass_fraction,
                            patch,
                            ghost_width_to_fill,
                            d_scalar_bdry_edge_conds,
                            d_bdry_edge_mass_fraction);
                        
                        /*
                         *  Set boundary conditions for cells corresponding to patch nodes.
                         */
                        
                        appu::CartesianBoundaryUtilities2::fillNodeBoundaryData(
                            "density",
                            density,
                            patch,
                            ghost_width_to_fill,
                            d_scalar_bdry_node_conds,
                            d_bdry_edge_density);
                        
                        appu::CartesianBoundaryUtilities2::fillNodeBoundaryData(
                            "momentum",
                            momentum,
                            patch,
                            ghost_width_to_fill,
                            d_vector_bdry_node_conds,
                            d_bdry_edge_momentum);
                        
                        appu::CartesianBoundaryUtilities2::fillNodeBoundaryData(
                            "total energy",
                            total_energy,
                            patch,
                            ghost_width_to_fill,
                            d_scalar_bdry_node_conds,
                            d_bdry_edge_total_energy);
                        
                        appu::CartesianBoundaryUtilities2::fillNodeBoundaryData(
                            "mass fraction",
                            mass_fraction,
                            patch,
                            ghost_width_to_fill,
                            d_scalar_bdry_node_conds,
                            d_bdry_edge_mass_fraction);
                    }
                    else if (d_dim == tbox::Dimension(3))
                    {
                        /*
                         *  Set boundary conditions for cells corresponding to patch faces.
                         */
                        
                        appu::CartesianBoundaryUtilities3::fillFaceBoundaryData(
                            "density",
                            density,
                            patch,
                            ghost_width_to_fill,
                            d_scalar_bdry_face_conds,
                            d_bdry_face_density);
                        
                        appu::CartesianBoundaryUtilities3::fillFaceBoundaryData(
                            "momentum",
                            momentum,
                            patch,
                            ghost_width_to_fill,
                            d_vector_bdry_face_conds,
                            d_bdry_face_momentum);
                        
                        appu::CartesianBoundaryUtilities3::fillFaceBoundaryData(
                            "total energy",
                            total_energy,
                            patch,
                            ghost_width_to_fill,
                            d_scalar_bdry_face_conds,
                            d_bdry_face_total_energy);
                        
                        appu::CartesianBoundaryUtilities3::fillFaceBoundaryData(
                            "mass fraction",
                            mass_fraction,
                            patch,
                            ghost_width_to_fill,
                            d_scalar_bdry_face_conds,
                            d_bdry_face_mass_fraction);
                        
                        /*
                         *  Set boundary conditions for cells corresponding to patch edges.
                         */
                        
                        appu::CartesianBoundaryUtilities3::fillEdgeBoundaryData(
                            "density",
                            density,
                            patch,
                            ghost_width_to_fill,
                            d_scalar_bdry_edge_conds,
                            d_bdry_face_density);
                        
                        appu::CartesianBoundaryUtilities3::fillEdgeBoundaryData(
                            "momentum",
                            momentum,
                            patch,
                            ghost_width_to_fill,
                            d_vector_bdry_edge_conds,
                            d_bdry_face_momentum);
                        
                        appu::CartesianBoundaryUtilities3::fillEdgeBoundaryData(
                            "total energy",
                            total_energy,
                            patch,
                            ghost_width_to_fill,
                            d_scalar_bdry_edge_conds,
                            d_bdry_face_total_energy);
                        
                        appu::CartesianBoundaryUtilities3::fillEdgeBoundaryData(
                            "mass fraction",
                            mass_fraction,
                            patch,
                            ghost_width_to_fill,
                            d_scalar_bdry_edge_conds,
                            d_bdry_face_mass_fraction);
                        
                        /*
                         *  Set boundary conditions for cells corresponding to patch nodes.
                         */
                        
                        appu::CartesianBoundaryUtilities3::fillNodeBoundaryData(
                            "density",
                            density,
                            patch,
                            ghost_width_to_fill,
                            d_scalar_bdry_node_conds,
                            d_bdry_face_density);
                        
                        appu::CartesianBoundaryUtilities3::fillNodeBoundaryData(
                            "momentum",
                            momentum,
                            patch,
                            ghost_width_to_fill,
                            d_vector_bdry_node_conds,
                            d_bdry_face_momentum);
                        
                        appu::CartesianBoundaryUtilities3::fillNodeBoundaryData(
                            "total energy",
                            total_energy,
                            patch,
                            ghost_width_to_fill,
                            d_scalar_bdry_node_conds,
                            d_bdry_face_total_energy);
                        
                        appu::CartesianBoundaryUtilities3::fillNodeBoundaryData(
                            "mass fraction",
                            mass_fraction,
                            patch,
                            ghost_width_to_fill,
                            d_scalar_bdry_node_conds,
                            d_bdry_face_mass_fraction);
                    }
                    
                    break;
                }
                case FIVE_EQN_ALLAIRE:
                {
                    boost::shared_ptr<pdat::CellData<double> > partial_density(
                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                            patch.getPatchData(d_partial_density, data_context)));
                    
                    boost::shared_ptr<pdat::CellData<double> > momentum(
                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                            patch.getPatchData(d_momentum, data_context)));
                    
                    boost::shared_ptr<pdat::CellData<double> > total_energy(
                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                            patch.getPatchData(d_total_energy, data_context)));
                    
                    boost::shared_ptr<pdat::CellData<double> > volume_fraction(
                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                            patch.getPatchData(d_volume_fraction, data_context)));
                    
#ifdef DEBUG_CHECK_ASSERTIONS
                    TBOX_ASSERT(partial_density);
                    TBOX_ASSERT(momentum);
                    TBOX_ASSERT(total_energy);
                    TBOX_ASSERT(volume_fraction);
                    
                    TBOX_ASSERT(partial_density->getGhostCellWidth() == d_num_ghosts);
                    TBOX_ASSERT(momentum->getGhostCellWidth() == d_num_ghosts);
                    TBOX_ASSERT(total_energy->getGhostCellWidth() == d_num_ghosts);
                    TBOX_ASSERT(volume_fraction->getGhostCellWidth() == d_num_ghosts);
#endif
                    
                    if (d_dim == tbox::Dimension(1))
                    {
                        // NOT YET IMPLEMENTED
                    }
                    else if (d_dim == tbox::Dimension(2))
                    {
                        /*
                         * Set boundary conditions for cells corresponding to patch edges.
                         */
                        
                        appu::CartesianBoundaryUtilities2::fillEdgeBoundaryData(
                            "partial density",
                            partial_density,
                            patch,
                            ghost_width_to_fill,
                            d_scalar_bdry_edge_conds,
                            d_bdry_edge_partial_density);
                        
                        appu::CartesianBoundaryUtilities2::fillEdgeBoundaryData(
                            "momentum",
                            momentum,
                            patch,
                            ghost_width_to_fill,
                            d_vector_bdry_edge_conds,
                            d_bdry_edge_momentum);
                        
                        appu::CartesianBoundaryUtilities2::fillEdgeBoundaryData(
                            "total energy",
                            total_energy,
                            patch,
                            ghost_width_to_fill,
                            d_scalar_bdry_edge_conds,
                            d_bdry_edge_total_energy);
                        
                        appu::CartesianBoundaryUtilities2::fillEdgeBoundaryData(
                            "volume fraction",
                            volume_fraction,
                            patch,
                            ghost_width_to_fill,
                            d_scalar_bdry_edge_conds,
                            d_bdry_edge_volume_fraction);
                        
                        /*
                         *  Set boundary conditions for cells corresponding to patch nodes.
                         */
                        
                        appu::CartesianBoundaryUtilities2::fillNodeBoundaryData(
                            "partial density",
                            partial_density,
                            patch,
                            ghost_width_to_fill,
                            d_scalar_bdry_node_conds,
                            d_bdry_edge_partial_density);
                        
                        appu::CartesianBoundaryUtilities2::fillNodeBoundaryData(
                            "momentum",
                            momentum,
                            patch,
                            ghost_width_to_fill,
                            d_vector_bdry_node_conds,
                            d_bdry_edge_momentum);
                        
                        appu::CartesianBoundaryUtilities2::fillNodeBoundaryData(
                            "total energy",
                            total_energy,
                            patch,
                            ghost_width_to_fill,
                            d_scalar_bdry_node_conds,
                            d_bdry_edge_total_energy);
                        
                        appu::CartesianBoundaryUtilities2::fillNodeBoundaryData(
                            "volume fraction",
                            volume_fraction,
                            patch,
                            ghost_width_to_fill,
                            d_scalar_bdry_node_conds,
                            d_bdry_edge_volume_fraction);
                    }
                    else if (d_dim == tbox::Dimension(3))
                    {
                        /*
                         *  Set boundary conditions for cells corresponding to patch faces.
                         */
                        
                        appu::CartesianBoundaryUtilities3::fillFaceBoundaryData(
                            "partial density",
                            partial_density,
                            patch,
                            ghost_width_to_fill,
                            d_scalar_bdry_face_conds,
                            d_bdry_face_partial_density);
                        
                        appu::CartesianBoundaryUtilities3::fillFaceBoundaryData(
                            "momentum",
                            momentum,
                            patch,
                            ghost_width_to_fill,
                            d_vector_bdry_face_conds,
                            d_bdry_face_momentum);
                        
                        appu::CartesianBoundaryUtilities3::fillFaceBoundaryData(
                            "total energy",
                            total_energy,
                            patch,
                            ghost_width_to_fill,
                            d_scalar_bdry_face_conds,
                            d_bdry_face_total_energy);
                        
                        appu::CartesianBoundaryUtilities3::fillFaceBoundaryData(
                            "volume fraction",
                            volume_fraction,
                            patch,
                            ghost_width_to_fill,
                            d_scalar_bdry_face_conds,
                            d_bdry_face_volume_fraction);
                        
                        /*
                         *  Set boundary conditions for cells corresponding to patch edges.
                         */
                        
                        appu::CartesianBoundaryUtilities3::fillEdgeBoundaryData(
                            "partial density",
                            partial_density,
                            patch,
                            ghost_width_to_fill,
                            d_scalar_bdry_edge_conds,
                            d_bdry_face_partial_density);
                        
                        appu::CartesianBoundaryUtilities3::fillEdgeBoundaryData(
                            "momentum",
                            momentum,
                            patch,
                            ghost_width_to_fill,
                            d_vector_bdry_edge_conds,
                            d_bdry_face_momentum);
                        
                        appu::CartesianBoundaryUtilities3::fillEdgeBoundaryData(
                            "total energy",
                            total_energy,
                            patch,
                            ghost_width_to_fill,
                            d_scalar_bdry_edge_conds,
                            d_bdry_face_total_energy);
                        
                        appu::CartesianBoundaryUtilities3::fillEdgeBoundaryData(
                            "volume fraction",
                            volume_fraction,
                            patch,
                            ghost_width_to_fill,
                            d_scalar_bdry_edge_conds,
                            d_bdry_face_volume_fraction);
                        
                        /*
                         *  Set boundary conditions for cells corresponding to patch nodes.
                         */
                        
                        appu::CartesianBoundaryUtilities3::fillNodeBoundaryData(
                            "partial density",
                            partial_density,
                            patch,
                            ghost_width_to_fill,
                            d_scalar_bdry_node_conds,
                            d_bdry_face_partial_density);
                        
                        appu::CartesianBoundaryUtilities3::fillNodeBoundaryData(
                            "momentum",
                            momentum,
                            patch,
                            ghost_width_to_fill,
                            d_vector_bdry_node_conds,
                            d_bdry_face_momentum);
                        
                        appu::CartesianBoundaryUtilities3::fillNodeBoundaryData(
                            "total energy",
                            total_energy,
                            patch,
                            ghost_width_to_fill,
                            d_scalar_bdry_node_conds,
                            d_bdry_face_total_energy);
                        
                        appu::CartesianBoundaryUtilities3::fillNodeBoundaryData(
                            "volume fraction",
                            volume_fraction,
                            patch,
                            ghost_width_to_fill,
                            d_scalar_bdry_node_conds,
                            d_bdry_face_volume_fraction);
                    }
                    
                    break;
                }
                default:
                {
                    TBOX_ERROR(d_object_name
                        << ": "
                        << "d_flow_model '"
                        << d_flow_model
                        << "' not yet implemented."
                        << std::endl);
                }
            }
            
            if (d_project_name == "2D double-Mach reflection")
            {
                TBOX_ASSERT(d_dim == tbox::Dimension(2));
                TBOX_ASSERT(d_flow_model == SINGLE_SPECIES);
                
                // Get the dimensions of box that covers the interior of patch.
                hier::Box dummy_box = patch.getBox();
                const hier::Box interior_box = dummy_box;
                const hier::IntVector interior_dims = interior_box.numberCells();
                
                // Get the dimensions of box that covers interior of patch plus
                // ghost cells.
                dummy_box.grow(d_num_ghosts);
                const hier::Box ghost_box = dummy_box;
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
                    
                    boost::shared_ptr<pdat::CellData<double> > density(
                            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                                patch.getPatchData(d_density, data_context)));
                    
                    boost::shared_ptr<pdat::CellData<double> > momentum(
                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                            patch.getPatchData(d_momentum, data_context)));
                    
                    boost::shared_ptr<pdat::CellData<double> > total_energy(
                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                            patch.getPatchData(d_total_energy, data_context)));
                    
#ifdef DEBUG_CHECK_ASSERTIONS
                    TBOX_ASSERT(density);
                    TBOX_ASSERT(momentum);
                    TBOX_ASSERT(total_energy);
                    
                    TBOX_ASSERT(density->getGhostCellWidth() == d_num_ghosts);
                    TBOX_ASSERT(momentum->getGhostCellWidth() == d_num_ghosts);
                    TBOX_ASSERT(total_energy->getGhostCellWidth() == d_num_ghosts);
#endif
                    
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
                                const int idx_cell = (i + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0];
                                
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
                                    const int idx_mirror_cell = (i + d_num_ghosts[0]) +
                                        (-j + d_num_ghosts[1] - 1)*ghostcell_dims[0];
                                    
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
                                const int idx_cell = (i + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0];
                                
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
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Number of ghost cells is not set yet."
                << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Variables are not set yet."
            << std::endl);
    }
}


void
EulerBoundaryConditions::readStateDataEntryForSingleSpecies(
    boost::shared_ptr<tbox::Database> db,
    const std::string& db_name,
    int array_indx,
    std::vector<double>& density,
    std::vector<double>& momentum,
    std::vector<double>& total_energy)
{
    TBOX_ASSERT(db);
    TBOX_ASSERT(!db_name.empty());
    TBOX_ASSERT(array_indx >= 0);
    TBOX_ASSERT(static_cast<int>(density.size()) > array_indx);
    TBOX_ASSERT(static_cast<int>(momentum.size()) > array_indx*d_dim.getValue());
    TBOX_ASSERT(static_cast<int>(total_energy.size()) > array_indx);
    
    if (db->keyExists("density"))
    {
        density[array_indx] = db->getDouble("density");
    }
    else
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "'density' entry missing from "
                   << db_name
                   << " input database."
                   << std::endl);
    }
    
    if (db->keyExists("momentum"))
    {
        std::vector<double> tmp_m = db->getDoubleVector("momentum");
        if (static_cast<int>(tmp_m.size()) < d_dim.getValue())
        {
            TBOX_ERROR(d_object_name
                       << ": "
                       << "Insufficient number of 'momentum' values"
                       << " given in "
                       << db_name
                       << " input database."
                       << std::endl);
        }
        for (int di = 0; di < d_dim.getValue(); di++)
        {
            momentum[array_indx*d_dim.getValue() + di] = tmp_m[di];
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "'momentum' entry missing from "
                   << db_name
                   << " input database."
                   << std::endl);
    }
   
    if (db->keyExists("total_energy"))
    {
        total_energy[array_indx] = db->getDouble("total_energy");
    }
    else
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "'total_energy' entry missing from "
                   << db_name
                   << " input database."
                   << std::endl);
    }
}


void
EulerBoundaryConditions::readStateDataEntryForFourEqnShyue(
    boost::shared_ptr<tbox::Database> db,
    const std::string& db_name,
    int array_indx,
    std::vector<double>& density,
    std::vector<double>& momentum,
    std::vector<double>& total_energy,
    std::vector<double>& mass_fraction)
{
    TBOX_ASSERT(db);
    TBOX_ASSERT(!db_name.empty());
    TBOX_ASSERT(array_indx >= 0);
    TBOX_ASSERT(static_cast<int>(density.size()) > array_indx);
    TBOX_ASSERT(static_cast<int>(momentum.size()) > array_indx*d_dim.getValue());
    TBOX_ASSERT(static_cast<int>(total_energy.size()) > array_indx);
    TBOX_ASSERT(static_cast<int>(mass_fraction.size()) > array_indx*d_num_species);
    
    if (db->keyExists("density"))
    {
        density[array_indx] = db->getDouble("density");
    }
    else
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "'density' entry missing from "
                   << db_name
                   << " input database."
                   << std::endl);
    }
    
    if (db->keyExists("momentum"))
    {
        std::vector<double> tmp_m = db->getDoubleVector("momentum");
        if (static_cast<int>(tmp_m.size()) < d_dim.getValue())
        {
            TBOX_ERROR(d_object_name
                       << ": "
                       << "Insufficient number of 'momentum' values"
                       << " given in "
                       << db_name
                       << " input database."
                       << std::endl);
        }
        for (int di = 0; di < d_dim.getValue(); di++)
        {
            momentum[array_indx*d_dim.getValue() + di] = tmp_m[di];
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "'momentum' entry missing from "
                   << db_name
                   << " input database."
                   << std::endl);
    }
   
    if (db->keyExists("total_energy"))
    {
        total_energy[array_indx] = db->getDouble("total_energy");
    }
    else
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "'total_energy' entry missing from "
                   << db_name
                   << " input database."
                   << std::endl);
    }
    
    if (db->keyExists("mass_fraction"))
    {
        std::vector<double> tmp_Y = db->getDoubleVector("mass_fraction");
        if (static_cast<int>(tmp_Y.size()) < d_num_species - 1)
        {
            TBOX_ERROR(d_object_name
                       << ": "
                       << "Insufficient number of 'mass_fraction' values"
                       << " given in "
                       << db_name
                       << " input database."
                       << std::endl);
        }
        double Y_last = 1.0;
        for (int si = 0; si < d_num_species - 1; si++)
        {
            mass_fraction[array_indx*d_num_species + si] = tmp_Y[si];
            Y_last -= tmp_Y[si];
        }
        mass_fraction[(array_indx + 1)*d_num_species - 1] = Y_last;
    }
    else
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "'mass_fraction' entry missing from "
                   << db_name
                   << " input database."
                   << std::endl);
    }
}


void
EulerBoundaryConditions::readStateDataEntryForFiveEqnAllaire(
    boost::shared_ptr<tbox::Database> db,
    const std::string& db_name,
    int array_indx,
    std::vector<double>& partial_density,
    std::vector<double>& momentum,
    std::vector<double>& total_energy,
    std::vector<double>& volume_fraction)
{
    TBOX_ASSERT(db);
    TBOX_ASSERT(!db_name.empty());
    TBOX_ASSERT(array_indx >= 0);
    TBOX_ASSERT(static_cast<int>(partial_density.size()) > array_indx*d_num_species);
    TBOX_ASSERT(static_cast<int>(momentum.size()) > array_indx*d_dim.getValue());
    TBOX_ASSERT(static_cast<int>(total_energy.size()) > array_indx);
    TBOX_ASSERT(static_cast<int>(volume_fraction.size()) > array_indx*d_num_species);
    
    if (db->keyExists("partial_density"))
    {
        std::vector<double> tmp_Z_rho = db->getDoubleVector("partial_density");
        if (static_cast<int>(tmp_Z_rho.size()) < d_num_species)
        {
            TBOX_ERROR(d_object_name
                       << ": "
                       << "Insufficient number of 'partial_density' values"
                       << " given in "
                       << db_name
                       << " input database."
                       << std::endl);
        }
        for (int si = 0; si < d_num_species; si++)
        {
            partial_density[array_indx*d_num_species + si] = tmp_Z_rho[si];
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "'partial_density' entry missing from "
                   << db_name
                   << " input database."
                   << std::endl);
    }
    
    if (db->keyExists("momentum"))
    {
        std::vector<double> tmp_m = db->getDoubleVector("momentum");
        if (static_cast<int>(tmp_m.size()) < d_dim.getValue())
        {
            TBOX_ERROR(d_object_name
                       << ": "
                       << "Insufficient number of 'momentum' values"
                       << " given in "
                       << db_name
                       << " input database."
                       << std::endl);
        }
        for (int di = 0; di < d_dim.getValue(); di++)
        {
            momentum[array_indx*d_dim.getValue() + di] = tmp_m[di];
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "'momentum' entry missing from "
                   << db_name
                   << " input database."
                   << std::endl);
    }
   
    if (db->keyExists("total_energy"))
    {
        total_energy[array_indx] = db->getDouble("total_energy");
    }
    else
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "'total_energy' entry missing from "
                   << db_name
                   << " input database."
                   << std::endl);
    }
    
    if (db->keyExists("volume_fraction"))
    {
        std::vector<double> tmp_Z = db->getDoubleVector("volume_fraction");
        if (static_cast<int>(tmp_Z.size()) < d_num_species - 1)
        {
            TBOX_ERROR(d_object_name
                       << ": "
                       << "Insufficient number of 'volume_fraction' values"
                       << " given in "
                       << db_name
                       << " input database."
                       << std::endl);
        }
        double Z_last = 1.0;
        for (int si = 0; si < d_num_species - 1; si++)
        {
            volume_fraction[array_indx*d_num_species + si] = tmp_Z[si];
            Z_last -= tmp_Z[si];
        }
        volume_fraction[(array_indx + 1)*d_num_species - 1] = Z_last;
    }
    else
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "'volume_fraction' entry missing from "
                   << db_name
                   << " input database."
                   << std::endl);
    }
}


/*
 * Set defaults for boundary conditions. Set to bogus values
 * for error checking.
 */
void
EulerBoundaryConditions::setDefaultBoundaryConditions()
{
    if (d_dim == tbox::Dimension(1))
    {
        d_master_bdry_node_conds.resize(NUM_1D_NODES);
        d_scalar_bdry_node_conds.resize(NUM_1D_NODES);
        d_vector_bdry_node_conds.resize(NUM_1D_NODES);
        for (int ni = 0; ni < NUM_1D_NODES; ni++)
        {
            d_master_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
            d_scalar_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
            d_vector_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
        }
        
        switch (d_flow_model)
        {
            case SINGLE_SPECIES:
            {
                d_bdry_node_density.resize(NUM_1D_NODES);
                d_bdry_node_momentum.resize(NUM_1D_NODES*d_dim.getValue());
                d_bdry_node_total_energy.resize(NUM_1D_NODES);
                
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_node_density);
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_node_momentum);
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_node_total_energy);
                
                break;
            }
            case FOUR_EQN_SHYUE:
            {
                d_bdry_node_density.resize(NUM_1D_NODES);
                d_bdry_node_momentum.resize(NUM_1D_NODES*d_dim.getValue());
                d_bdry_node_total_energy.resize(NUM_1D_NODES);
                d_bdry_node_mass_fraction.resize(NUM_1D_NODES*d_num_species);
                
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_node_density);
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_node_momentum);
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_node_total_energy);
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_node_mass_fraction);
                
                break;
            }
            case FIVE_EQN_ALLAIRE:
            {
                d_bdry_node_partial_density.resize(NUM_1D_NODES*d_num_species);
                d_bdry_node_momentum.resize(NUM_1D_NODES*d_dim.getValue());
                d_bdry_node_total_energy.resize(NUM_1D_NODES);
                d_bdry_node_volume_fraction.resize(NUM_1D_NODES*d_num_species);
                
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_node_partial_density);
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_node_momentum);
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_node_total_energy);
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_node_volume_fraction);
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "d_flow_model '"
                    << d_flow_model
                    << "' not yet implemented."
                    << std::endl);
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
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
        
        switch (d_flow_model)
        {
            case SINGLE_SPECIES:
            {
                d_bdry_edge_density.resize(NUM_2D_EDGES);
                d_bdry_edge_momentum.resize(NUM_2D_EDGES*d_dim.getValue());
                d_bdry_edge_total_energy.resize(NUM_2D_EDGES);
                
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_edge_density);
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_edge_momentum);
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_edge_total_energy);
                
                break;
            }
            case FOUR_EQN_SHYUE:
            {
                d_bdry_edge_density.resize(NUM_2D_EDGES);
                d_bdry_edge_momentum.resize(NUM_2D_EDGES*d_dim.getValue());
                d_bdry_edge_total_energy.resize(NUM_2D_EDGES);
                d_bdry_edge_mass_fraction.resize(NUM_2D_NODES*d_num_species);
                
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_edge_density);
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_edge_momentum);
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_edge_total_energy);
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_edge_mass_fraction);
                
                break;
            }
            case FIVE_EQN_ALLAIRE:
            {
                d_bdry_edge_partial_density.resize(NUM_2D_EDGES*d_num_species);
                d_bdry_edge_momentum.resize(NUM_2D_EDGES*d_dim.getValue());
                d_bdry_edge_total_energy.resize(NUM_2D_EDGES);
                d_bdry_edge_volume_fraction.resize(NUM_2D_EDGES*d_num_species);
                
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_edge_density);
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_edge_momentum);
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_edge_total_energy);
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_edge_volume_fraction);
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "d_flow_model '"
                    << d_flow_model
                    << "' not yet implemented."
                    << std::endl);
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
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
        
        switch (d_flow_model)
        {
            case SINGLE_SPECIES:
            {
                d_bdry_face_density.resize(NUM_3D_FACES);
                d_bdry_face_momentum.resize(NUM_3D_FACES*d_dim.getValue());
                d_bdry_face_total_energy.resize(NUM_3D_FACES);
                
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_face_density);
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_face_momentum);
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_face_total_energy);
                
                break;
            }
            case FOUR_EQN_SHYUE:
            {
                d_bdry_face_density.resize(NUM_3D_FACES);
                d_bdry_face_momentum.resize(NUM_3D_FACES*d_dim.getValue());
                d_bdry_face_total_energy.resize(NUM_3D_FACES);
                d_bdry_face_mass_fraction.resize(NUM_3D_FACES*d_num_species);
                
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_face_density);
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_face_momentum);
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_face_total_energy);
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_face_mass_fraction);
                
                break;
            }
            case FIVE_EQN_ALLAIRE:
            {
                d_bdry_face_partial_density.resize(NUM_3D_FACES*d_num_species);
                d_bdry_face_momentum.resize(NUM_3D_FACES*d_dim.getValue());
                d_bdry_face_total_energy.resize(NUM_3D_FACES);
                d_bdry_face_volume_fraction.resize(NUM_3D_FACES*d_num_species);
                
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_edge_partial_density);
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_edge_momentum);
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_edge_total_energy);
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_edge_volume_fraction);
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "d_flow_model '"
                    << d_flow_model
                    << "' not yet implemented."
                    << std::endl);
            }
        }
    }
}
