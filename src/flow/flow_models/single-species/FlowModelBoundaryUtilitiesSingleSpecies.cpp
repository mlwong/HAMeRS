#include "flow/flow_models/single-species/FlowModelBoundaryUtilitiesSingleSpecies.hpp"

#include "util/basic_boundary_conditions/CartesianBoundaryDefines.hpp"

#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"

#include <algorithm>

// Integer constant for debugging improperly set boundary data.
#define BOGUS_BDRY_LOC (-9999)

FlowModelBoundaryUtilitiesSingleSpecies::FlowModelBoundaryUtilitiesSingleSpecies(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const int& num_species,
    const int& num_eqn,
    const HAMERS_SHARED_PTR<EquationOfStateMixingRules>& equation_of_state_mixing_rules):
        FlowModelBoundaryUtilities(
            object_name,
            dim,
            num_species,
            num_eqn,
            equation_of_state_mixing_rules)
{
    std::vector<Real*> thermo_properties_ptr;
    
    const int num_thermo_properties = d_equation_of_state_mixing_rules->
        getNumberOfSpeciesThermodynamicProperties();
    
    thermo_properties_ptr.reserve(num_thermo_properties);
    d_thermo_properties.resize(num_thermo_properties);
    
    for (int ti = 0; ti < num_thermo_properties; ti++)
    {
        thermo_properties_ptr.push_back(&d_thermo_properties[ti]);
    }
    
    d_equation_of_state_mixing_rules->getSpeciesThermodynamicProperties(
        thermo_properties_ptr,
        0);
    
    /*
     * Set default sizes for containers of boundary conditions.
     */
    
    if (d_dim == tbox::Dimension(1))
    {
        d_bdry_node_adiabatic_no_slip_vel.resize(NUM_1D_NODES);
        
        d_bdry_node_isothermal_no_slip_T.resize(NUM_1D_NODES);
        d_bdry_node_isothermal_no_slip_vel.resize(NUM_1D_NODES);
        
        d_bdry_node_nonreflecting_outflow_p_t.resize(NUM_1D_NODES);
        d_bdry_node_nonreflecting_outflow_sigma.resize(NUM_1D_NODES);
        d_bdry_node_nonreflecting_outflow_beta.resize(NUM_1D_NODES);
        d_bdry_node_nonreflecting_outflow_length_char.resize(NUM_1D_NODES);

    }
    else if (d_dim == tbox::Dimension(2))
    {
        d_bdry_edge_adiabatic_no_slip_vel.resize(NUM_2D_EDGES*2);
        
        d_bdry_edge_isothermal_no_slip_T.resize(NUM_2D_EDGES);
        d_bdry_edge_isothermal_no_slip_vel.resize(NUM_2D_EDGES*2);
        
        d_bdry_edge_nonreflecting_outflow_p_t.resize(NUM_2D_EDGES);
        d_bdry_edge_nonreflecting_outflow_sigma.resize(NUM_2D_EDGES);
        d_bdry_edge_nonreflecting_outflow_beta.resize(NUM_2D_EDGES);
        d_bdry_edge_nonreflecting_outflow_length_char.resize(NUM_2D_EDGES);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        d_bdry_face_adiabatic_no_slip_vel.resize(NUM_3D_FACES*3);
        
        d_bdry_face_isothermal_no_slip_T.resize(NUM_3D_FACES);
        d_bdry_face_isothermal_no_slip_vel.resize(NUM_3D_FACES*3);
        
        d_bdry_face_nonreflecting_outflow_p_t.resize(NUM_3D_FACES);
        d_bdry_face_nonreflecting_outflow_sigma.resize(NUM_3D_FACES);
        d_bdry_face_nonreflecting_outflow_beta.resize(NUM_3D_FACES);
        d_bdry_face_nonreflecting_outflow_length_char.resize(NUM_3D_FACES);
    }
}


/*
 * Function to read 1d boundary data from input database.
 * Node locations that have boundary conditions identified are removed from the container.
 */
void
FlowModelBoundaryUtilitiesSingleSpecies::getFromInput1d(
    const HAMERS_SHARED_PTR<tbox::Database>& input_db,
    std::vector<int>& node_locs,
    std::vector<int>& node_conds,
    const hier::IntVector& periodic)
{
    TBOX_DIM_ASSERT(periodic.getDim() == tbox::Dimension(1));
    TBOX_ASSERT(static_cast<int>(node_locs.size()) <= NUM_1D_NODES);
    if (static_cast<int>(node_locs.size()) > 0)
    {
        TBOX_ASSERT(*min_element(node_locs.begin(), node_locs.end()) >= 0);
        TBOX_ASSERT(*max_element(node_locs.begin(), node_locs.end()) < NUM_1D_NODES);
    }
    TBOX_ASSERT(static_cast<int>(node_conds.size()) == NUM_1D_NODES);
    
    if (!input_db)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBoundaryUtilitiesSingleSpecies::getFromInput1d()\n"
            << "no input database supplied"
            << std::endl);
    }
    
    if (static_cast<int>(node_locs.size()) > 0)
    {
        read1dBdryNodes(
            input_db,
            node_locs,
            node_conds,
            periodic);
    }
}


/*
 * Function to read 2d boundary data from input database.
 * Node and edge locations that have boundary conditions identified are removed from the
 * containers.
 */
void
FlowModelBoundaryUtilitiesSingleSpecies::getFromInput2d(
    const HAMERS_SHARED_PTR<tbox::Database>& input_db,
    std::vector<int>& edge_locs,
    std::vector<int>& node_locs,
    std::vector<int>& edge_conds,
    std::vector<int>& node_conds,
    const hier::IntVector& periodic)
{
    TBOX_DIM_ASSERT(periodic.getDim() == tbox::Dimension(2));
    TBOX_ASSERT(static_cast<int>(edge_locs.size()) <= NUM_2D_EDGES);
    TBOX_ASSERT(static_cast<int>(node_locs.size()) <= NUM_2D_NODES);
    if (static_cast<int>(edge_locs.size()) > 0)
    {
        TBOX_ASSERT(*min_element(edge_locs.begin(), edge_locs.end()) >= 0);
        TBOX_ASSERT(*max_element(edge_locs.begin(), edge_locs.end()) < NUM_2D_EDGES);
    }
    if (static_cast<int>(node_locs.size()) > 0)
    {
        TBOX_ASSERT(*min_element(node_locs.begin(), node_locs.end()) >= 0);
        TBOX_ASSERT(*max_element(node_locs.begin(), node_locs.end()) < NUM_2D_NODES);
    }
    TBOX_ASSERT(static_cast<int>(edge_conds.size()) == NUM_2D_EDGES);
    TBOX_ASSERT(static_cast<int>(node_conds.size()) == NUM_2D_NODES);
    
    if (!input_db)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBoundaryUtilitiesSingleSpecies::getFromInput2d()\n"
            << "no input database supplied"
            << std::endl);
    }
    
    if (static_cast<int>(edge_locs.size()) > 0)
    {
        read2dBdryEdges(
            input_db,
            edge_locs,
            edge_conds,
            periodic);
    }
    
    if (static_cast<int>(node_locs.size()) > 0)
    {
        read2dBdryNodes(
            input_db,
            node_locs,
            edge_conds,
            node_conds,
            periodic);
    }
}


/*
 * Function to read 3d boundary data from input database.
 * Node, edge and face locations that have boundary conditions identified are removed from
 * the containers.
 */
void
FlowModelBoundaryUtilitiesSingleSpecies::getFromInput3d(
    const HAMERS_SHARED_PTR<tbox::Database>& input_db,
    std::vector<int>& face_locs,
    std::vector<int>& edge_locs,
    std::vector<int>& node_locs,
    std::vector<int>& face_conds,
    std::vector<int>& edge_conds,
    std::vector<int>& node_conds,
    const hier::IntVector& periodic)
{
    TBOX_DIM_ASSERT(periodic.getDim() == tbox::Dimension(3));
    TBOX_ASSERT(static_cast<int>(face_locs.size()) <= NUM_3D_FACES);
    TBOX_ASSERT(static_cast<int>(edge_locs.size()) <= NUM_3D_EDGES);
    TBOX_ASSERT(static_cast<int>(node_locs.size()) <= NUM_3D_NODES);
    if (static_cast<int>(face_locs.size()) > 0)
    {
        TBOX_ASSERT(*min_element(face_locs.begin(), face_locs.end()) >= 0);
        TBOX_ASSERT(*max_element(face_locs.begin(), face_locs.end()) < NUM_3D_FACES);
    }
    if (static_cast<int>(edge_locs.size()) > 0)
    {
        TBOX_ASSERT(*min_element(edge_locs.begin(), edge_locs.end()) >= 0);
        TBOX_ASSERT(*max_element(edge_locs.begin(), edge_locs.end()) < NUM_3D_EDGES);
    }
    if (static_cast<int>(node_locs.size()) > 0)
    {
        TBOX_ASSERT(*min_element(node_locs.begin(), node_locs.end()) >= 0);
        TBOX_ASSERT(*max_element(node_locs.begin(), node_locs.end()) < NUM_3D_NODES);
    }
    TBOX_ASSERT(static_cast<int>(face_conds.size()) == NUM_3D_FACES);
    TBOX_ASSERT(static_cast<int>(edge_conds.size()) == NUM_3D_EDGES);
    TBOX_ASSERT(static_cast<int>(node_conds.size()) == NUM_3D_NODES);
    
    if (!input_db)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBoundaryUtilitiesSingleSpecies::getFromInput3d()\n"
            << "no input database supplied"
            << std::endl);
    }
    
    if (static_cast<int>(face_locs.size()) > 0)
    {
        read3dBdryFaces(
            input_db,
            face_locs,
            face_conds,
            periodic);
    }
    
    if (static_cast<int>(edge_locs.size()) > 0)
    {
        read3dBdryEdges(
            input_db,
            edge_locs,
            face_conds,
            edge_conds,
            periodic);
    }
    
    if (static_cast<int>(node_locs.size()) > 0)
    {
        read3dBdryNodes(
            input_db,
            node_locs,
            face_conds,
            node_conds,
            periodic);
    }
}


/*
 * Function that returns the integer edge boundary location
 * corresponding to the given node location and node boundary
 * condition.
 *
 * If the node boundary condition type or node location are unknown,
 * or the boundary condition type is inconsistent with the node location,
 * an error code (-1) is returned.
 */
int
FlowModelBoundaryUtilitiesSingleSpecies::getEdgeLocationForNodeBdry(
    int node_loc,
    int node_btype)
{
    int ret_edge = -1;
    
    switch (node_btype)
    {
        case BDRY_COND::FLOW_MODEL::XADIABATIC_NO_SLIP:
        case BDRY_COND::FLOW_MODEL::XISOTHERMAL_NO_SLIP:
        {
            if (node_loc == NODE_BDRY_LOC_2D::XLO_YLO ||
                node_loc == NODE_BDRY_LOC_2D::XLO_YHI)
            {
                ret_edge = BDRY_LOC::XLO;
            }
            else if (node_loc == NODE_BDRY_LOC_2D::XHI_YLO ||
                     node_loc == NODE_BDRY_LOC_2D::XHI_YHI)
            {
                ret_edge = BDRY_LOC::XHI;
            }
            break;
        }
        case BDRY_COND::FLOW_MODEL::YADIABATIC_NO_SLIP:
        case BDRY_COND::FLOW_MODEL::YISOTHERMAL_NO_SLIP:
        {
            if (node_loc == NODE_BDRY_LOC_2D::XLO_YLO ||
                node_loc == NODE_BDRY_LOC_2D::XHI_YLO)
            {
                ret_edge = BDRY_LOC::YLO;
            }
            else if (node_loc == NODE_BDRY_LOC_2D::XLO_YHI ||
                     node_loc == NODE_BDRY_LOC_2D::XHI_YHI)
            {
                ret_edge = BDRY_LOC::YHI;
            }
            break;
        }
        default:
        {
        }
    }
    
    return ret_edge;
}


/*
 * Function that returns the integer face boundary location
 * corresponding to the given edge location and edge boundary
 * condition.
 *
 * If the edge boundary condition type or edge location are unknown,
 * or the boundary condition type is inconsistent with the edge location,
 * an error code (-1) is returned.
 */
int
FlowModelBoundaryUtilitiesSingleSpecies::getFaceLocationForEdgeBdry(
    int edge_loc,
    int edge_btype)
{
    int ret_face = -1;
    
    switch (edge_btype)
    {
        case BDRY_COND::FLOW_MODEL::XADIABATIC_NO_SLIP:
        case BDRY_COND::FLOW_MODEL::XISOTHERMAL_NO_SLIP:
        {
            if (edge_loc == EDGE_BDRY_LOC_3D::XLO_ZLO ||
                edge_loc == EDGE_BDRY_LOC_3D::XLO_ZHI ||
                edge_loc == EDGE_BDRY_LOC_3D::XLO_YLO ||
                edge_loc == EDGE_BDRY_LOC_3D::XLO_YHI)
            {
                ret_face = BDRY_LOC::XLO;
            }
            else if (edge_loc == EDGE_BDRY_LOC_3D::XHI_ZLO ||
                     edge_loc == EDGE_BDRY_LOC_3D::XHI_ZHI ||
                     edge_loc == EDGE_BDRY_LOC_3D::XHI_YLO ||
                     edge_loc == EDGE_BDRY_LOC_3D::XHI_YHI)
            {
                ret_face = BDRY_LOC::XHI;
            }
            break;
        }
        case BDRY_COND::FLOW_MODEL::YADIABATIC_NO_SLIP:
        case BDRY_COND::FLOW_MODEL::YISOTHERMAL_NO_SLIP:
        {
            if (edge_loc == EDGE_BDRY_LOC_3D::YLO_ZLO ||
                edge_loc == EDGE_BDRY_LOC_3D::YLO_ZHI ||
                edge_loc == EDGE_BDRY_LOC_3D::XLO_YLO ||
                edge_loc == EDGE_BDRY_LOC_3D::XHI_YLO)
            {
                ret_face = BDRY_LOC::YLO;
            }
            else if (edge_loc == EDGE_BDRY_LOC_3D::YHI_ZLO ||
                     edge_loc == EDGE_BDRY_LOC_3D::YHI_ZHI ||
                     edge_loc == EDGE_BDRY_LOC_3D::XLO_YHI ||
                     edge_loc == EDGE_BDRY_LOC_3D::XHI_YHI)
            {
                ret_face = BDRY_LOC::YHI;
            }
            break;
        }
        case BDRY_COND::FLOW_MODEL::ZADIABATIC_NO_SLIP:
        case BDRY_COND::FLOW_MODEL::ZISOTHERMAL_NO_SLIP:
        {
           if (edge_loc == EDGE_BDRY_LOC_3D::YLO_ZLO ||
               edge_loc == EDGE_BDRY_LOC_3D::YHI_ZLO ||
               edge_loc == EDGE_BDRY_LOC_3D::XLO_ZLO ||
               edge_loc == EDGE_BDRY_LOC_3D::XHI_ZLO)
            {
                ret_face = BDRY_LOC::ZLO;
            }
            else if (edge_loc == EDGE_BDRY_LOC_3D::YLO_ZHI ||
                     edge_loc == EDGE_BDRY_LOC_3D::YHI_ZHI ||
                     edge_loc == EDGE_BDRY_LOC_3D::XLO_ZHI ||
                     edge_loc == EDGE_BDRY_LOC_3D::XHI_ZHI)
            {
                ret_face = BDRY_LOC::ZHI;
            }
            break;
        }
        default:
        {
        }
    }
    
    return ret_face;
}


/*
 * Function that returns the integer face boundary location
 * corresponding to the given node location and node boundary
 * condition.
 *
 * If the node boundary condition type or node location are unknown,
 * or the boundary condition type is inconsistent with the node location,
 * an error code (-1) is returned.
 */
int
FlowModelBoundaryUtilitiesSingleSpecies::getFaceLocationForNodeBdry(
    int node_loc,
    int node_btype)
{
    int ret_face = -1;
    
    switch (node_btype)
    {
        case BDRY_COND::FLOW_MODEL::XADIABATIC_NO_SLIP:
        case BDRY_COND::FLOW_MODEL::XISOTHERMAL_NO_SLIP:
        {
            if (node_loc == NODE_BDRY_LOC_3D::XLO_YLO_ZLO ||
                node_loc == NODE_BDRY_LOC_3D::XLO_YHI_ZLO ||
                node_loc == NODE_BDRY_LOC_3D::XLO_YLO_ZHI ||
                node_loc == NODE_BDRY_LOC_3D::XLO_YHI_ZHI)
            {
                ret_face = BDRY_LOC::XLO;
            }
            else if (node_loc == NODE_BDRY_LOC_3D::XHI_YLO_ZLO ||
                     node_loc == NODE_BDRY_LOC_3D::XHI_YHI_ZLO ||
                     node_loc == NODE_BDRY_LOC_3D::XHI_YLO_ZHI ||
                     node_loc == NODE_BDRY_LOC_3D::XHI_YHI_ZHI)
            {
                ret_face = BDRY_LOC::XHI;
            }
            break;
        }
        case BDRY_COND::FLOW_MODEL::YADIABATIC_NO_SLIP:
        case BDRY_COND::FLOW_MODEL::YISOTHERMAL_NO_SLIP:
        {
            if (node_loc == NODE_BDRY_LOC_3D::XLO_YLO_ZLO ||
                node_loc == NODE_BDRY_LOC_3D::XHI_YLO_ZLO ||
                node_loc == NODE_BDRY_LOC_3D::XLO_YLO_ZHI ||
                node_loc == NODE_BDRY_LOC_3D::XHI_YLO_ZHI)
            {
                ret_face = BDRY_LOC::YLO;
            }
            else if (node_loc == NODE_BDRY_LOC_3D::XLO_YHI_ZLO ||
                     node_loc == NODE_BDRY_LOC_3D::XHI_YHI_ZLO ||
                     node_loc == NODE_BDRY_LOC_3D::XLO_YHI_ZHI ||
                     node_loc == NODE_BDRY_LOC_3D::XHI_YHI_ZHI)
            {
                ret_face = BDRY_LOC::YHI;
            }
            
            break;
        }
        case BDRY_COND::FLOW_MODEL::ZADIABATIC_NO_SLIP:
        case BDRY_COND::FLOW_MODEL::ZISOTHERMAL_NO_SLIP:
        {
            if (node_loc == NODE_BDRY_LOC_3D::XLO_YLO_ZLO ||
                node_loc == NODE_BDRY_LOC_3D::XHI_YLO_ZLO ||
                node_loc == NODE_BDRY_LOC_3D::XLO_YHI_ZLO ||
                node_loc == NODE_BDRY_LOC_3D::XHI_YHI_ZLO)
            {
                ret_face = BDRY_LOC::ZLO;
            }
            else if (node_loc == NODE_BDRY_LOC_3D::XLO_YLO_ZHI ||
                     node_loc == NODE_BDRY_LOC_3D::XHI_YLO_ZHI ||
                     node_loc == NODE_BDRY_LOC_3D::XLO_YHI_ZHI ||
                     node_loc == NODE_BDRY_LOC_3D::XHI_YHI_ZHI)
            {
                ret_face = BDRY_LOC::ZHI;
            }
            break;
        }
        default:
        {
        }
    }
    
    return ret_face;
}


/*
 * Function to fill 1d node boundary values for a patch.
 * Node locations that have boundary conditions identified are removed from the container.
 */
void
FlowModelBoundaryUtilitiesSingleSpecies::fill1dNodeBoundaryData(
    const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& conservative_var_data,
    const hier::Patch& patch,
    std::vector<int>& bdry_node_locs,
    const std::vector<int>& bdry_node_conds,
    const std::vector<std::vector<Real> >& bdry_node_values,
    const hier::IntVector& ghost_width_to_fill)
{
    TBOX_ASSERT(static_cast<int>(conservative_var_data.size()) == 3);
    for (int vi = 0; vi < static_cast<int>(conservative_var_data.size()); vi++)
    {
        TBOX_ASSERT(conservative_var_data[vi]);
    }
    TBOX_ASSERT(static_cast<int>(bdry_node_locs.size()) <= NUM_1D_NODES);
    if (static_cast<int>(bdry_node_locs.size()) > 0)
    {
        TBOX_ASSERT(*min_element(bdry_node_locs.begin(), bdry_node_locs.end()) >= 0);
        TBOX_ASSERT(*max_element(bdry_node_locs.begin(), bdry_node_locs.end()) < NUM_1D_NODES);
    }
    TBOX_ASSERT(static_cast<int>(bdry_node_conds.size()) == NUM_1D_NODES);
    TBOX_ASSERT(static_cast<int>(bdry_node_values.size()) == 3);
    for (int vi = 0; vi < static_cast<int>(bdry_node_values.size()); vi++)
    {
        TBOX_ASSERT(static_cast<int>(bdry_node_values[vi].size()) ==
                    NUM_1D_NODES*(conservative_var_data[vi]->getDepth()));
    }
    
    TBOX_DIM_ASSERT(ghost_width_to_fill.getDim() == tbox::Dimension(1));
    
    for (int vi = 0; vi < static_cast<int>(conservative_var_data.size()); vi++)
    {
        TBOX_ASSERT_OBJDIM_EQUALITY3(*conservative_var_data[vi], patch, ghost_width_to_fill);
    }
    
    const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
        HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    TBOX_ASSERT(patch_geom);
    
    hier::IntVector num_ghosts(conservative_var_data[0]->getGhostCellWidth());
    for (int vi = 1; vi < static_cast<int>(conservative_var_data.size()); vi++)
    {
        num_ghosts = hier::IntVector::min(
            num_ghosts,
            conservative_var_data[vi]->getGhostCellWidth());
    }
    
    /*
     * Determine the ghost cell width to fill.
     */
    
    hier::IntVector gcw_to_fill(tbox::Dimension(1));
    
    // If the ghost fill width is not used, it is set to the ghost cell width of the data.
    if (ghost_width_to_fill == -hier::IntVector::getOne(tbox::Dimension(1)))
    {
        gcw_to_fill = num_ghosts;
    }
    else
    {
        gcw_to_fill = hier::IntVector::min(
            num_ghosts,
            ghost_width_to_fill);
    }
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box& interior_box(patch.getBox());
    hier::Index interior_box_lo_idx(interior_box.lower());
    hier::Index interior_box_hi_idx(interior_box.upper());
    
    /*
     * Offset the indices.
     */
    
    interior_box_lo_idx = interior_box_lo_idx - interior_box.lower();
    interior_box_hi_idx = interior_box_hi_idx - interior_box.lower();
    
    const std::vector<hier::BoundaryBox>& node_bdry =
        patch_geom->getCodimensionBoundaries(BDRY::NODE1D);
    
    for (int ni = 0; ni < static_cast<int>(node_bdry.size()); ni++)
    {
        TBOX_ASSERT(node_bdry[ni].getBoundaryType() == BDRY::NODE1D);
        
        int node_loc = node_bdry[ni].getLocationIndex();
        
        if (std::find(bdry_node_locs.begin(), bdry_node_locs.end(), node_loc) !=
            bdry_node_locs.end())
        {
            hier::Box fill_box(patch_geom->getBoundaryFillBox(
                node_bdry[ni],
                interior_box,
                gcw_to_fill));
            
            hier::Index fill_box_lo_idx(fill_box.lower());
            hier::Index fill_box_hi_idx(fill_box.upper());
            
            /*
             * Offset the indices.
             */
            
            fill_box_lo_idx = fill_box_lo_idx - interior_box.lower();
            fill_box_hi_idx = fill_box_hi_idx - interior_box.lower();
            
            if ((bdry_node_conds[node_loc] == BDRY_COND::FLOW_MODEL::ADIABATIC_NO_SLIP) ||
                (bdry_node_conds[node_loc] == BDRY_COND::FLOW_MODEL::ISOTHERMAL_NO_SLIP))
            {
                /*
                 * Get the pointers to the conservative variables.
                 * Get the numbers of ghost cells and the dimensions of the ghost cell boxes of
                 * the conservative variables.
                 */
                
                std::vector<Real*> Q;
                Q.reserve(d_num_eqn);
                
                std::vector<hier::IntVector> num_subghosts_conservative_var;
                std::vector<hier::IntVector> subghostcell_dims_conservative_var;
                
                num_subghosts_conservative_var.reserve(3);
                subghostcell_dims_conservative_var.reserve(3);
                
                int count_eqn = 0;
                
                for (int vi = 0; vi < static_cast<int>(conservative_var_data.size()); vi++)
                {
                    int depth = conservative_var_data[vi]->getDepth();
                    
                    for (int di = 0; di < depth; di++)
                    {
                        // If the last element of the conservative variable vector is not in the
                        // system of equations, ignore it.
                        if (count_eqn >= d_num_eqn)
                            break;
                        
                        Q.push_back(conservative_var_data[vi]->getPointer(di));
                        
                        count_eqn++;
                    }
                    
                    num_subghosts_conservative_var.push_back(
                        conservative_var_data[vi]->getGhostCellWidth());
                    subghostcell_dims_conservative_var.push_back(
                        conservative_var_data[vi]->getGhostBox().numberCells());
                }
                
                // Get the thermodynamic properties of the species.
                std::vector<const Real*> thermo_properties_ptr;
                thermo_properties_ptr.reserve(static_cast<int> (d_thermo_properties.size()));
                for (int ti = 0; ti < static_cast<int> (d_thermo_properties.size()); ti++)
                {
                    thermo_properties_ptr.push_back(&d_thermo_properties[ti]);
                }
                
                if (bdry_node_conds[node_loc] == BDRY_COND::FLOW_MODEL::ADIABATIC_NO_SLIP)
                {
                    for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                    {
                        const int idx_cell_rho = i + num_subghosts_conservative_var[0][0];
                        const int idx_cell_mom = i + num_subghosts_conservative_var[1][0];
                        const int idx_cell_E = i + num_subghosts_conservative_var[2][0];
                        
                        int idx_cell_pivot_rho = idx_cell_rho;
                        int idx_cell_pivot_mom = idx_cell_mom;
                        int idx_cell_pivot_E = idx_cell_E;
                        
                        if (node_loc == BDRY_LOC::XLO)
                        {
                            idx_cell_pivot_rho = interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                num_subghosts_conservative_var[0][0];
                            
                            idx_cell_pivot_mom = interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                num_subghosts_conservative_var[1][0];
                            
                            idx_cell_pivot_E = interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                num_subghosts_conservative_var[2][0];
                            
                        }
                        else if (node_loc == BDRY_LOC::XHI)
                        {
                            idx_cell_pivot_rho = interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                num_subghosts_conservative_var[0][0];
                            
                            idx_cell_pivot_mom = interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                num_subghosts_conservative_var[1][0];
                            
                            idx_cell_pivot_E = interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                num_subghosts_conservative_var[2][0];
                        }
                        
                        /*
                         * Set the values for density and momentum.
                         */
                        
                        Q[0][idx_cell_rho] = Q[0][idx_cell_pivot_rho];
                        Q[1][idx_cell_mom] = -Q[1][idx_cell_pivot_mom] +
                            Real(2)*Q[0][idx_cell_pivot_rho]*d_bdry_node_adiabatic_no_slip_vel[node_loc];
                        
                        /*
                         * Set the values for total internal energy.
                         */
                        
                        Real epsilon_pivot = (Q[2][idx_cell_pivot_E] -
                            Real(1)/Real(2)*Q[1][idx_cell_pivot_mom]*Q[1][idx_cell_pivot_mom]/
                            Q[0][idx_cell_pivot_rho])/Q[0][idx_cell_pivot_rho];
                        
                        Real p_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                            getPressure(
                                &Q[0][idx_cell_pivot_rho],
                                &epsilon_pivot,
                                thermo_properties_ptr);
                        
                        Real T_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                            getTemperature(
                                &Q[0][idx_cell_pivot_rho],
                                &p_pivot,
                                thermo_properties_ptr);
                        
                        Real T = T_pivot;
                        
                        Real epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
                            getInternalEnergyFromTemperature(
                                &Q[0][idx_cell_rho],
                                &T,
                                thermo_properties_ptr);
                        
                        Real E = Q[0][idx_cell_rho]*epsilon +
                            Real(1)/Real(2)*Q[1][idx_cell_mom]*Q[1][idx_cell_mom]/Q[0][idx_cell_rho];
                        
                        Q[2][idx_cell_E] = E;
                    }
                    
                    // Remove node locations that have boundary conditions identified.
                    bdry_node_locs.erase(std::remove(bdry_node_locs.begin(), bdry_node_locs.end(), node_loc),
                        bdry_node_locs.end());
                }
                else if (bdry_node_conds[node_loc] == BDRY_COND::FLOW_MODEL::ISOTHERMAL_NO_SLIP)
                {
                    for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                    {
                        const int idx_cell_rho = i + num_subghosts_conservative_var[0][0];
                        const int idx_cell_mom = i + num_subghosts_conservative_var[1][0];
                        const int idx_cell_E   = i + num_subghosts_conservative_var[2][0];
                        
                        int idx_cell_pivot_rho = idx_cell_rho;
                        int idx_cell_pivot_mom = idx_cell_mom;
                        int idx_cell_pivot_E = idx_cell_E;
                        
                        if (node_loc == BDRY_LOC::XLO)
                        {
                            idx_cell_pivot_rho = interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                num_subghosts_conservative_var[0][0];
                            
                            idx_cell_pivot_mom = interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                num_subghosts_conservative_var[1][0];
                            
                            idx_cell_pivot_E = interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                num_subghosts_conservative_var[2][0];
                        }
                        else if (node_loc == BDRY_LOC::XHI)
                        {
                            idx_cell_pivot_rho = interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                num_subghosts_conservative_var[0][0];
                            
                            idx_cell_pivot_mom = interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                num_subghosts_conservative_var[1][0];
                            
                            idx_cell_pivot_E = interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                num_subghosts_conservative_var[2][0];
                        }
                        
                        /*
                         * Set the values for density, momentum and total internal energy.
                         */
                        
                        Real epsilon_pivot = (Q[2][idx_cell_pivot_E] -
                            Real(1)/Real(2)*Q[1][idx_cell_pivot_mom]*Q[1][idx_cell_pivot_mom]/
                            Q[0][idx_cell_pivot_rho])/Q[0][idx_cell_pivot_rho];
                        
                        Real p_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                            getPressure(
                                &Q[0][idx_cell_pivot_rho],
                                &epsilon_pivot,
                                thermo_properties_ptr);
                        
                        Real p = p_pivot;
                        
                        Real T_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                            getTemperature(
                                &Q[0][idx_cell_pivot_rho],
                                &p_pivot,
                                thermo_properties_ptr);
                        
                        Real T = -T_pivot + Real(2)*d_bdry_node_isothermal_no_slip_T[node_loc];
                        
                        Real rho = d_equation_of_state_mixing_rules->getEquationOfState()->
                            getDensity(
                                &p,
                                &T,
                                thermo_properties_ptr);
                        
                        Real u = -Q[1][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                            Real(2)*d_bdry_node_isothermal_no_slip_vel[node_loc];
                        
                        Q[0][idx_cell_rho] = rho;
                        Q[1][idx_cell_mom] = rho*u;
                        
                        Real epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
                            getInternalEnergyFromTemperature(
                                &rho,
                                &T,
                                thermo_properties_ptr);
                        
                        Real E = rho*epsilon + Real(1)/Real(2)*Q[1][idx_cell_mom]*Q[1][idx_cell_mom]/rho;
                        
                        Q[2][idx_cell_E] = E;
                    }
                    
                    // Remove node locations that have boundary conditions identified.
                    bdry_node_locs.erase(std::remove(bdry_node_locs.begin(), bdry_node_locs.end(), node_loc),
                        bdry_node_locs.end());
                }
            }
        }
    }
    
    for (int ni = 0; ni < static_cast<int>(node_bdry.size()); ni++)
    {
        TBOX_ASSERT(node_bdry[ni].getBoundaryType() == BDRY::NODE1D);
        
        int node_loc = node_bdry[ni].getLocationIndex();
        
        if (std::find(bdry_node_locs.begin(), bdry_node_locs.end(), node_loc) !=
            bdry_node_locs.end())
        {
            TBOX_ERROR("FlowModelBoundaryUtilitiesSingleSpecies::fill1dNodeBoundaryData()\n"
                << "Invalid node boundary condition!\n"
                << "node_loc = '" << node_loc << "'." << std::endl
                << "bdry_node_conds[node_loc] = '" << bdry_node_conds[node_loc] << "'."
                << std::endl);
        }
    }
}


/*
 * Function to fill 2d edge boundary values for a patch.
 * Edge locations that have boundary conditions identified are removed from the container.
 */
void
FlowModelBoundaryUtilitiesSingleSpecies::fill2dEdgeBoundaryData(
    const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& conservative_var_data,
    const hier::Patch& patch,
    std::vector<int>& bdry_edge_locs,
    const std::vector<int>& bdry_edge_conds,
    const std::vector<std::vector<Real> >& bdry_edge_values,
    const hier::IntVector& ghost_width_to_fill)
{
    TBOX_ASSERT(static_cast<int>(conservative_var_data.size()) == 3);
    for (int vi = 0; vi < static_cast<int>(conservative_var_data.size()); vi++)
    {
        TBOX_ASSERT(conservative_var_data[vi]);
    }
    TBOX_ASSERT(static_cast<int>(bdry_edge_locs.size()) <= NUM_2D_EDGES);
    if (static_cast<int>(bdry_edge_locs.size()) > 0)
    {
        TBOX_ASSERT(*min_element(bdry_edge_locs.begin(), bdry_edge_locs.end()) >= 0);
        TBOX_ASSERT(*max_element(bdry_edge_locs.begin(), bdry_edge_locs.end()) < NUM_2D_EDGES);
    }
    TBOX_ASSERT(static_cast<int>(bdry_edge_conds.size()) == NUM_2D_EDGES);
    TBOX_ASSERT(static_cast<int>(bdry_edge_values.size()) == 3);
    for (int vi = 0; vi < static_cast<int>(bdry_edge_values.size()); vi++)
    {
        TBOX_ASSERT(static_cast<int>(bdry_edge_values[vi].size()) ==
                    NUM_2D_EDGES*(conservative_var_data[vi]->getDepth()));
    }
    
    TBOX_DIM_ASSERT(ghost_width_to_fill.getDim() == tbox::Dimension(2));
    
    for (int vi = 0; vi < static_cast<int>(conservative_var_data.size()); vi++)
    {
        TBOX_ASSERT_OBJDIM_EQUALITY3(*conservative_var_data[vi], patch, ghost_width_to_fill);
    }
    
    const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
        HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    TBOX_ASSERT(patch_geom);
    
    hier::IntVector num_ghosts(conservative_var_data[0]->getGhostCellWidth());
    for (int vi = 1; vi < static_cast<int>(conservative_var_data.size()); vi++)
    {
        num_ghosts = hier::IntVector::min(
            num_ghosts,
            conservative_var_data[vi]->getGhostCellWidth());
    }
    
    /*
     * Determine the ghost cell width to fill.
     */
    
    hier::IntVector gcw_to_fill(tbox::Dimension(2));
    
    // If the ghost fill width is not used, it is set to the ghost cell width of the data.
    if (ghost_width_to_fill == -hier::IntVector::getOne(tbox::Dimension(2)))
    {
        gcw_to_fill = num_ghosts;
    }
    else
    {
        gcw_to_fill = hier::IntVector::min(
            num_ghosts,
            ghost_width_to_fill);
    }
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box& interior_box(patch.getBox());
    hier::Index interior_box_lo_idx(interior_box.lower());
    hier::Index interior_box_hi_idx(interior_box.upper());
    
    /*
     * Offset the indices.
     */
    
    interior_box_lo_idx = interior_box_lo_idx - interior_box.lower();
    interior_box_hi_idx = interior_box_hi_idx - interior_box.lower();
    
    const std::vector<hier::BoundaryBox>& edge_bdry =
        patch_geom->getCodimensionBoundaries(BDRY::EDGE2D);
    
    /*
     * 1. Fill ghost cells for edges with regular boundary conditions.
     * 2. Fill ghost cells for edges with non-reflecting boundary conditions.
     *    (assuming non-reflecting boundaries don't touch each other)
     */
    
    /*
     * 1. Fill ghost cells for edges with regular boundary conditions.
     */
    for (int ei = 0; ei < static_cast<int>(edge_bdry.size()); ei++)
    {
        TBOX_ASSERT(edge_bdry[ei].getBoundaryType() == BDRY::EDGE2D);
        
        int edge_loc = edge_bdry[ei].getLocationIndex();
        
        if (std::find(bdry_edge_locs.begin(), bdry_edge_locs.end(), edge_loc) !=
            bdry_edge_locs.end())
        {
            hier::Box fill_box(patch_geom->getBoundaryFillBox(
                edge_bdry[ei],
                interior_box,
                gcw_to_fill));
            
            hier::Index fill_box_lo_idx(fill_box.lower());
            hier::Index fill_box_hi_idx(fill_box.upper());
            
            /*
             * Offset the indices.
             */
            
            fill_box_lo_idx = fill_box_lo_idx - interior_box.lower();
            fill_box_hi_idx = fill_box_hi_idx - interior_box.lower();
            
            if ((bdry_edge_conds[edge_loc] == BDRY_COND::FLOW_MODEL::ADIABATIC_NO_SLIP) ||
                (bdry_edge_conds[edge_loc] == BDRY_COND::FLOW_MODEL::ISOTHERMAL_NO_SLIP))
            {
                /*
                 * Get the pointers to the conservative variables.
                 * Get the numbers of ghost cells and the dimensions of the ghost cell boxes of
                 * the conservative variables.
                 */
                
                std::vector<Real*> Q;
                Q.reserve(d_num_eqn);
                
                std::vector<hier::IntVector> num_subghosts_conservative_var;
                std::vector<hier::IntVector> subghostcell_dims_conservative_var;
                
                num_subghosts_conservative_var.reserve(3);
                subghostcell_dims_conservative_var.reserve(3);
                
                int count_eqn = 0;
                
                for (int vi = 0; vi < static_cast<int>(conservative_var_data.size()); vi++)
                {
                    int depth = conservative_var_data[vi]->getDepth();
                    
                    for (int di = 0; di < depth; di++)
                    {
                        // If the last element of the conservative variable vector is not in the
                        // system of equations, ignore it.
                        if (count_eqn >= d_num_eqn)
                            break;
                        
                        Q.push_back(conservative_var_data[vi]->getPointer(di));
                        
                        count_eqn++;
                    }
                    
                    num_subghosts_conservative_var.push_back(
                        conservative_var_data[vi]->getGhostCellWidth());
                    subghostcell_dims_conservative_var.push_back(
                        conservative_var_data[vi]->getGhostBox().numberCells());
                }
                
                // Get the thermodynamic properties of the species.
                std::vector<const Real*> thermo_properties_ptr;
                thermo_properties_ptr.reserve(static_cast<int> (d_thermo_properties.size()));
                for (int ti = 0; ti < static_cast<int> (d_thermo_properties.size()); ti++)
                {
                    thermo_properties_ptr.push_back(&d_thermo_properties[ti]);
                }
                
                if (bdry_edge_conds[edge_loc] == BDRY_COND::FLOW_MODEL::ADIABATIC_NO_SLIP)
                {
                    for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                    {
                        for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                        {
                            const int idx_cell_rho = (i + num_subghosts_conservative_var[0][0]) +
                                (j + num_subghosts_conservative_var[0][1])*
                                    subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_mom = (i + num_subghosts_conservative_var[1][0]) +
                                (j + num_subghosts_conservative_var[1][1])*
                                    subghostcell_dims_conservative_var[1][0];
                            
                            const int idx_cell_E = (i + num_subghosts_conservative_var[2][0]) +
                                (j + num_subghosts_conservative_var[2][1])*
                                    subghostcell_dims_conservative_var[2][0];
                            
                            int idx_cell_pivot_rho = idx_cell_rho;
                            int idx_cell_pivot_mom = idx_cell_mom;
                            int idx_cell_pivot_E = idx_cell_E;
                            
                            if (edge_loc == BDRY_LOC::XLO)
                            {
                                idx_cell_pivot_rho = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                        num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0];
                                
                                idx_cell_pivot_mom = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                        num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0];
                                
                                idx_cell_pivot_E = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                        num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0];
                            }
                            else if (edge_loc == BDRY_LOC::XHI)
                            {
                                idx_cell_pivot_rho = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                        num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0];
                                
                                idx_cell_pivot_mom = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                        num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0];
                                
                                idx_cell_pivot_E = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                        num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0];
                            }
                            else if (edge_loc == BDRY_LOC::YLO)
                            {
                                idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                    (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) +
                                        num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0];
                                
                                idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                    (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) +
                                        num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0];
                                
                                idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                    (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) +
                                        num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0];
                            }
                            else if (edge_loc == BDRY_LOC::YHI)
                            {
                                idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                    (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) +
                                        num_subghosts_conservative_var[0][1])*
                                            subghostcell_dims_conservative_var[0][0];
                                
                                idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                    (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) +
                                        num_subghosts_conservative_var[1][1])*
                                            subghostcell_dims_conservative_var[1][0];
                                
                                idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                    (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) +
                                        num_subghosts_conservative_var[2][1])*
                                            subghostcell_dims_conservative_var[2][0];
                            }
                            
                            /*
                             * Set the values for density and momentum.
                             */
                            
                            Q[0][idx_cell_rho] = Q[0][idx_cell_pivot_rho];
                            Q[1][idx_cell_mom] = -Q[1][idx_cell_pivot_mom] +
                                Real(2)*Q[0][idx_cell_pivot_rho]*d_bdry_edge_adiabatic_no_slip_vel[edge_loc*2];
                            Q[2][idx_cell_mom] = -Q[2][idx_cell_pivot_mom] +
                                Real(2)*Q[0][idx_cell_pivot_rho]*d_bdry_edge_adiabatic_no_slip_vel[edge_loc*2 + 1];
                            
                            /*
                             * Set the values for total internal energy.
                             */
                            
                            Real epsilon_pivot = (Q[3][idx_cell_pivot_E] -
                                Real(1)/Real(2)*(Q[1][idx_cell_pivot_mom]*Q[1][idx_cell_pivot_mom] +
                                Q[2][idx_cell_pivot_mom]*Q[2][idx_cell_pivot_mom])/
                                Q[0][idx_cell_pivot_rho])/Q[0][idx_cell_pivot_rho];
                            
                            Real p_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getPressure(
                                    &Q[0][idx_cell_pivot_rho],
                                    &epsilon_pivot,
                                    thermo_properties_ptr);
                            
                            Real T_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getTemperature(
                                    &Q[0][idx_cell_pivot_rho],
                                    &p_pivot,
                                    thermo_properties_ptr);
                            
                            Real T = T_pivot;
                            
                            Real epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getInternalEnergyFromTemperature(
                                    &Q[0][idx_cell_rho],
                                    &T,
                                    thermo_properties_ptr);
                            
                            Real E = Q[0][idx_cell_rho]*epsilon +
                                Real(1)/Real(2)*(Q[1][idx_cell_mom]*Q[1][idx_cell_mom] + Q[2][idx_cell_mom]*Q[2][idx_cell_mom])/
                                Q[0][idx_cell_rho];
                            
                            Q[3][idx_cell_E] = E;
                        }
                    }
                    
                    // Remove edge locations that have boundary conditions identified.
                    bdry_edge_locs.erase(std::remove(bdry_edge_locs.begin(), bdry_edge_locs.end(), edge_loc),
                        bdry_edge_locs.end());
                }
                else if (bdry_edge_conds[edge_loc] == BDRY_COND::FLOW_MODEL::ISOTHERMAL_NO_SLIP)
                {
                    for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                    {
                        for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                        {
                            const int idx_cell_rho = (i + num_subghosts_conservative_var[0][0]) +
                                (j + num_subghosts_conservative_var[0][1])*
                                    subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_mom = (i + num_subghosts_conservative_var[1][0]) +
                                (j + num_subghosts_conservative_var[1][1])*
                                    subghostcell_dims_conservative_var[1][0];
                            
                            const int idx_cell_E = (i + num_subghosts_conservative_var[2][0]) +
                                (j + num_subghosts_conservative_var[2][1])*
                                    subghostcell_dims_conservative_var[2][0];
                            
                            int idx_cell_pivot_rho = idx_cell_rho;
                            int idx_cell_pivot_mom = idx_cell_mom;
                            int idx_cell_pivot_E = idx_cell_E;
                            
                            if (edge_loc == BDRY_LOC::XLO)
                            {
                                idx_cell_pivot_rho = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                        num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0];
                                
                                idx_cell_pivot_mom = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                        num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0];
                                
                                idx_cell_pivot_E = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                        num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0];
                            }
                            else if (edge_loc == BDRY_LOC::XHI)
                            {
                                idx_cell_pivot_rho = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                        num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0];
                                
                                idx_cell_pivot_mom = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                        num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0];
                                
                                idx_cell_pivot_E = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                        num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0];
                            }
                            else if (edge_loc == BDRY_LOC::YLO)
                            {
                                idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                    (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) +
                                        num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0];
                                
                                idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                    (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) +
                                        num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0];
                                
                                idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                    (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) +
                                        num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0];
                            }
                            else if (edge_loc == BDRY_LOC::YHI)
                            {
                                idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                    (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) +
                                        num_subghosts_conservative_var[0][1])*
                                            subghostcell_dims_conservative_var[0][0];
                                
                                idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                    (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) +
                                        num_subghosts_conservative_var[1][1])*
                                            subghostcell_dims_conservative_var[1][0];
                                
                                idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                    (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) +
                                        num_subghosts_conservative_var[2][1])*
                                            subghostcell_dims_conservative_var[2][0];
                            }
                            
                            /*
                             * Set the values for density, momentum and total internal energy.
                             */
                            
                            Real epsilon_pivot = (Q[3][idx_cell_pivot_E] -
                                Real(1)/Real(2)*(Q[1][idx_cell_pivot_mom]*Q[1][idx_cell_pivot_mom] +
                                Q[2][idx_cell_pivot_mom]*Q[2][idx_cell_pivot_mom])/
                                Q[0][idx_cell_pivot_rho])/Q[0][idx_cell_pivot_rho];
                            
                            Real p_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getPressure(
                                    &Q[0][idx_cell_pivot_rho],
                                    &epsilon_pivot,
                                    thermo_properties_ptr);
                            
                            Real p = p_pivot;
                            
                            Real T_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getTemperature(
                                    &Q[0][idx_cell_pivot_rho],
                                    &p_pivot,
                                    thermo_properties_ptr);
                            
                            Real T = -T_pivot + Real(2)*d_bdry_edge_isothermal_no_slip_T[edge_loc];
                            
                            Real rho = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getDensity(
                                    &p,
                                    &T,
                                    thermo_properties_ptr);
                            
                            Real u = -Q[1][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                                Real(2)*d_bdry_edge_isothermal_no_slip_vel[edge_loc*2];
                            Real v = -Q[2][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                                Real(2)*d_bdry_edge_isothermal_no_slip_vel[edge_loc*2 + 1];
                            
                            Q[0][idx_cell_rho] = rho;
                            Q[1][idx_cell_mom] = rho*u;
                            Q[2][idx_cell_mom] = rho*v;
                            
                            Real epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getInternalEnergyFromTemperature(
                                    &rho,
                                    &T,
                                    thermo_properties_ptr);
                            
                            Real E = rho*epsilon +
                                Real(1)/Real(2)*(Q[1][idx_cell_mom]*Q[1][idx_cell_mom] + Q[2][idx_cell_mom]*Q[2][idx_cell_mom])/
                                rho;
                            
                            Q[3][idx_cell_E] = E;
                        }
                    }
                    
                    // Remove edge locations that have boundary conditions identified.
                    bdry_edge_locs.erase(std::remove(bdry_edge_locs.begin(), bdry_edge_locs.end(), edge_loc),
                        bdry_edge_locs.end());
                }
            }
        }
    }
    
    /*
     * 2. Fill ghost cells for edges with non-reflecting boundary conditions.
     */
    for (int ei = 0; ei < static_cast<int>(edge_bdry.size()); ei++)
    {
        TBOX_ASSERT(edge_bdry[ei].getBoundaryType() == BDRY::EDGE2D);
        
        int edge_loc = edge_bdry[ei].getLocationIndex();
        
        if (std::find(bdry_edge_locs.begin(), bdry_edge_locs.end(), edge_loc) !=
            bdry_edge_locs.end())
        {
            if (bdry_edge_conds[edge_loc] == BDRY_COND::FLOW_MODEL::NONREFLECTING_OUTFLOW)
            {
                // Follow the method in
                // Motheau, Emmanuel, Ann Almgren, and John B. Bell.
                // "Navier–stokes characteristic boundary conditions using ghost cells."
                // AIAA Journal (2017): 3399-3408.
                
                hier::Box fill_box(patch_geom->getBoundaryFillBox(
                    edge_bdry[ei],
                    interior_box,
                    gcw_to_fill));
            
                hier::Index fill_box_lo_idx(fill_box.lower());
                hier::Index fill_box_hi_idx(fill_box.upper());
            
                /*
                 * Offset the indices.
                 */
            
                fill_box_lo_idx = fill_box_lo_idx - interior_box.lower();
                fill_box_hi_idx = fill_box_hi_idx - interior_box.lower();
                
                /*
                 * Get the pointers to the conservative variables.
                 * Get the numbers of ghost cells and the dimensions of the ghost cell boxes of
                 * the conservative variables.
                 */
                
                std::vector<Real*> Q;
                Q.reserve(d_num_eqn);
                
                std::vector<hier::IntVector> num_subghosts_conservative_var;
                std::vector<hier::IntVector> subghostcell_dims_conservative_var;
                
                num_subghosts_conservative_var.reserve(3);
                subghostcell_dims_conservative_var.reserve(3);
                
                int count_eqn = 0;
                
                for (int vi = 0; vi < static_cast<int>(conservative_var_data.size()); vi++)
                {
                    int depth = conservative_var_data[vi]->getDepth();
                    
                    for (int di = 0; di < depth; di++)
                    {
                        // If the last element of the conservative variable vector is not in the
                        // system of equations, ignore it.
                        if (count_eqn >= d_num_eqn)
                            break;
                        
                        Q.push_back(conservative_var_data[vi]->getPointer(di));
                        
                        count_eqn++;
                    }
                    
                    num_subghosts_conservative_var.push_back(
                        conservative_var_data[vi]->getGhostCellWidth());
                    subghostcell_dims_conservative_var.push_back(
                        conservative_var_data[vi]->getGhostBox().numberCells());
                }
                
                // Get the thermodynamic properties of the species.
                std::vector<const Real*> thermo_properties_ptr;
                thermo_properties_ptr.reserve(static_cast<int> (d_thermo_properties.size()));
                for (int ti = 0; ti < static_cast<int> (d_thermo_properties.size()); ti++)
                {
                    thermo_properties_ptr.push_back(&d_thermo_properties[ti]);
                }
                
                if (edge_loc == BDRY_LOC::XLO)
                {
                    const int num_ghosts_to_fill = fill_box_hi_idx[0] - fill_box_lo_idx[0] + 1;
                    TBOX_ASSERT(fill_box_hi_idx[0] == interior_box_lo_idx[0] - 1);
                    if (num_ghosts_to_fill > 6)
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelBoundaryUtilitiesSingleSpecies::fill2dEdgeBoundaryData()\n"
                            << "Non-reflecting outflow BC doesn't support more than six ghost cells yet!");
                    }
                    
                    for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                    {
                        // Get the grid spacing.
                        const double* const dx = patch_geom->getDx();
                        
                        const int idx_cell_rho_x_R = (interior_box_lo_idx[0] + num_subghosts_conservative_var[0][0]) +
                            (j + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0];
                        
                        const int idx_cell_rho_x_RR = (interior_box_lo_idx[0] + 1 + num_subghosts_conservative_var[0][0]) +
                            (j + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0];
                        
                        const int idx_cell_rho_x_RRR = (interior_box_lo_idx[0] + 2 + num_subghosts_conservative_var[0][0]) +
                            (j + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0];
                        
                        const int idx_cell_mom_x_R = (interior_box_lo_idx[0] + num_subghosts_conservative_var[1][0]) +
                            (j + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0];
                        
                        const int idx_cell_mom_x_RR = (interior_box_lo_idx[0] + 1 + num_subghosts_conservative_var[1][0]) +
                            (j + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0];
                        
                        const int idx_cell_mom_x_RRR = (interior_box_lo_idx[0] + 2 + num_subghosts_conservative_var[1][0]) +
                            (j + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0];
                        
                        const int idx_cell_E_x_R = (interior_box_lo_idx[0] + num_subghosts_conservative_var[2][0]) +
                            (j + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0];
                        
                        const int idx_cell_E_x_RR = (interior_box_lo_idx[0] + 1 + num_subghosts_conservative_var[2][0]) +
                            (j + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0];
                        
                        const int idx_cell_E_x_RRR = (interior_box_lo_idx[0] + 2 + num_subghosts_conservative_var[2][0]) +
                            (j + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0];
                        
                        const Real& rho_x_R   = Q[0][idx_cell_rho_x_R];
                        const Real& rho_x_RR  = Q[0][idx_cell_rho_x_RR];
                        const Real& rho_x_RRR = Q[0][idx_cell_rho_x_RRR];
                        
                        const Real u_x_R   = Q[1][idx_cell_mom_x_R]/rho_x_R;
                        const Real u_x_RR  = Q[1][idx_cell_mom_x_RR]/rho_x_RR;
                        const Real u_x_RRR = Q[1][idx_cell_mom_x_RRR]/rho_x_RRR;
                        
                        const Real v_x_R   = Q[2][idx_cell_mom_x_R]/rho_x_R;
                        const Real v_x_RR  = Q[2][idx_cell_mom_x_RR]/rho_x_RR;
                        const Real v_x_RRR = Q[2][idx_cell_mom_x_RRR]/rho_x_RRR;
                        
                        const Real half = Real(1)/Real(2);
                        const Real epsilon_x_R   = Q[3][idx_cell_E_x_R]/rho_x_R - half*(u_x_R*u_x_R + v_x_R*v_x_R);
                        const Real epsilon_x_RR  = Q[3][idx_cell_E_x_RR]/rho_x_RR - half*(u_x_RR*u_x_RR + v_x_RR*v_x_RR);
                        const Real epsilon_x_RRR = Q[3][idx_cell_E_x_RRR]/rho_x_RRR - half*(u_x_RRR*u_x_RRR + v_x_RRR*v_x_RRR);
                        
                        const Real p_x_R = d_equation_of_state_mixing_rules->getEquationOfState()->
                            getPressure(
                                &rho_x_R,
                                &epsilon_x_R,
                                thermo_properties_ptr);
                        
                        const Real p_x_RR = d_equation_of_state_mixing_rules->getEquationOfState()->
                            getPressure(
                                &rho_x_RR,
                                &epsilon_x_RR,
                                thermo_properties_ptr);
                        
                        const Real p_x_RRR = d_equation_of_state_mixing_rules->getEquationOfState()->
                            getPressure(
                                &rho_x_RRR,
                                &epsilon_x_RRR,
                                thermo_properties_ptr);
                        
                        /*
                         * Compute derivatives in x-direction.
                         */
                        
                        const Real drho_dx = -(rho_x_RRR - Real(4)*rho_x_RR + Real(3)*rho_x_R)/(Real(2)*Real(dx[0]));
                        const Real du_dx   = -(u_x_RRR - Real(4)*u_x_RR + Real(3)*u_x_R)/(Real(2)*Real(dx[0]));
                        const Real dv_dx   = -(v_x_RRR - Real(4)*v_x_RR + Real(3)*v_x_R)/(Real(2)*Real(dx[0]));
                        const Real dp_dx   = -(p_x_RRR - Real(4)*p_x_RR + Real(3)*p_x_R)/(Real(2)*Real(dx[0]));
                        
                        /*
                         * Compute derivatives in y-direction.
                         */
                        
                        Real du_dy = Real(0);
                        Real dv_dy = Real(0);
                        Real dp_dy = Real(0);
                        
                        if ((j + num_subghosts_conservative_var[0][1] == 0) ||
                            (j + num_subghosts_conservative_var[1][1] == 0) ||
                            (j + num_subghosts_conservative_var[2][1] == 0))
                        // if (((patch_geom->getTouchesRegularBoundary(1, 0)) && (j == interior_box_lo_idx[1])) ||
                        //     ((j + num_subghosts_conservative_var[0][1] == 0) ||
                        //      (j + num_subghosts_conservative_var[1][1] == 0) ||
                        //      (j + num_subghosts_conservative_var[2][1] == 0)))
                        {
                            // Patch is touching periodic boundary.
                            
                            const int idx_cell_rho_y_T = (interior_box_lo_idx[0] + num_subghosts_conservative_var[0][0]) +
                                (j + 1 + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_mom_y_T = (interior_box_lo_idx[0] + num_subghosts_conservative_var[1][0]) +
                                (j + 1 + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0];
                            
                            const int idx_cell_E_y_T = (interior_box_lo_idx[0] + num_subghosts_conservative_var[2][0]) +
                                (j + 1 + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0];
                            
                            const Real& rho_y_T = Q[0][idx_cell_rho_y_T];
                            const Real u_y_T = Q[1][idx_cell_mom_y_T]/rho_y_T;
                            const Real v_y_T = Q[2][idx_cell_mom_y_T]/rho_y_T;
                            const Real epsilon_y_T = Q[3][idx_cell_E_y_T]/rho_y_T - half*(u_y_T*u_y_T + v_y_T*v_y_T);
                            
                            const Real p_y_T = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getPressure(
                                    &rho_y_T,
                                    &epsilon_y_T,
                                    thermo_properties_ptr);
                            
                            // One-sided derivatives.
                            du_dy = (u_y_T - u_x_R)/Real(dx[1]);
                            dv_dy = (v_y_T - v_x_R)/Real(dx[1]);
                            dp_dy = (p_y_T - p_x_R)/Real(dx[1]);
                        }
                        else if ((j + num_subghosts_conservative_var[0][1] + 1 == subghostcell_dims_conservative_var[0][1]) ||
                                 (j + num_subghosts_conservative_var[1][1] + 1 == subghostcell_dims_conservative_var[1][1]) ||
                                 (j + num_subghosts_conservative_var[2][1] + 1 == subghostcell_dims_conservative_var[2][1]))
                        // else if (((patch_geom->getTouchesRegularBoundary(1, 1)) && (j == interior_box_hi_idx[1])) ||
                        //          ((j + num_subghosts_conservative_var[0][1] + 1 == subghostcell_dims_conservative_var[0][1]) ||
                        //           (j + num_subghosts_conservative_var[1][1] + 1 == subghostcell_dims_conservative_var[1][1]) ||
                        //           (j + num_subghosts_conservative_var[2][1] + 1 == subghostcell_dims_conservative_var[2][1])))
                        {
                            // Patch is touching periodic boundary.
                            
                            const int idx_cell_rho_y_B = (interior_box_lo_idx[0] + num_subghosts_conservative_var[0][0]) +
                                (j - 1 + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_mom_y_B = (interior_box_lo_idx[0] + num_subghosts_conservative_var[1][0]) +
                                (j - 1 + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0];
                            
                            const int idx_cell_E_y_B = (interior_box_lo_idx[0] + num_subghosts_conservative_var[2][0]) +
                                (j - 1 + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0];
                            
                            const Real& rho_y_B = Q[0][idx_cell_rho_y_B];
                            const Real u_y_B = Q[1][idx_cell_mom_y_B]/rho_y_B;
                            const Real v_y_B = Q[2][idx_cell_mom_y_B]/rho_y_B;
                            const Real epsilon_y_B = Q[3][idx_cell_E_y_B]/rho_y_B - half*(u_y_B*u_y_B + v_y_B*v_y_B);
                            
                            const Real p_y_B = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getPressure(
                                    &rho_y_B,
                                    &epsilon_y_B,
                                    thermo_properties_ptr);
                            
                            // One-sided derivatives.
                            du_dy = (u_x_R - u_y_B)/Real(dx[1]);
                            dv_dy = (v_x_R - v_y_B)/Real(dx[1]);
                            dp_dy = (p_x_R - p_y_B)/Real(dx[1]);
                        }
                        else
                        {
                            const int idx_cell_rho_y_B = (interior_box_lo_idx[0] + num_subghosts_conservative_var[0][0]) +
                                (j - 1 + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_rho_y_T = (interior_box_lo_idx[0] + num_subghosts_conservative_var[0][0]) +
                                (j + 1 + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_mom_y_B = (interior_box_lo_idx[0] + num_subghosts_conservative_var[1][0]) +
                                (j - 1 + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0];
                            
                            const int idx_cell_mom_y_T = (interior_box_lo_idx[0] + num_subghosts_conservative_var[1][0]) +
                                (j + 1 + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0];
                            
                            const int idx_cell_E_y_B = (interior_box_lo_idx[0] + num_subghosts_conservative_var[2][0]) +
                                (j - 1 + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0];
                            
                            const int idx_cell_E_y_T = (interior_box_lo_idx[0] + num_subghosts_conservative_var[2][0]) +
                                (j + 1 + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0];
                            
                            const Real& rho_y_B = Q[0][idx_cell_rho_y_B];
                            const Real& rho_y_T = Q[0][idx_cell_rho_y_T];
                            
                            const Real u_y_B = Q[1][idx_cell_mom_y_B]/rho_y_B;
                            const Real u_y_T = Q[1][idx_cell_mom_y_T]/rho_y_T;
                            
                            const Real v_y_B = Q[2][idx_cell_mom_y_B]/rho_y_B;
                            const Real v_y_T = Q[2][idx_cell_mom_y_T]/rho_y_T;
                            
                            const Real epsilon_y_B = Q[3][idx_cell_E_y_B]/rho_y_B - half*(u_y_B*u_y_B + v_y_B*v_y_B);
                            const Real epsilon_y_T = Q[3][idx_cell_E_y_T]/rho_y_T - half*(u_y_T*u_y_T + v_y_T*v_y_T);
                            
                            const Real p_y_B = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getPressure(
                                    &rho_y_B,
                                    &epsilon_y_B,
                                    thermo_properties_ptr);
                            
                            const Real p_y_T = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getPressure(
                                    &rho_y_T,
                                    &epsilon_y_T,
                                    thermo_properties_ptr);
                            
                            // Central derivatives.
                            du_dy = (u_y_T - u_y_B)/(Real(2)*Real(dx[1]));
                            dv_dy = (v_y_T - v_y_B)/(Real(2)*Real(dx[1]));
                            dp_dy = (p_y_T - p_y_B)/(Real(2)*Real(dx[1]));
                        }
                        
                        const Real c_x_R = d_equation_of_state_mixing_rules->getEquationOfState()->
                            getSoundSpeed(
                                &rho_x_R,
                                &p_x_R,
                                thermo_properties_ptr);
                        
                        const Real lambda_4 = u_x_R + c_x_R;
                        
                        // Compute vector Lambda^(-1) * L.
                        
                        Real Lambda_inv_L[4];
                        
                        const Real& p_t         = d_bdry_edge_nonreflecting_outflow_p_t[edge_loc];
                        const Real& sigma       = d_bdry_edge_nonreflecting_outflow_sigma[edge_loc];
                        const Real& beta        = d_bdry_edge_nonreflecting_outflow_beta[edge_loc];
                        const Real& length_char = d_bdry_edge_nonreflecting_outflow_length_char[edge_loc];
                        
                        const Real T_4 = v_x_R*(dp_dy + rho_x_R*c_x_R*du_dy) + rho_x_R*c_x_R*c_x_R*dv_dy;
                        
                        const Real M_sq = (u_x_R*u_x_R + v_x_R*v_x_R)/(c_x_R*c_x_R);
                        const Real K = sigma*c_x_R*(Real(1) - M_sq)/length_char;
                        
                        Lambda_inv_L[0] = dp_dx - rho_x_R*c_x_R*du_dx;
                        Lambda_inv_L[1] = c_x_R*c_x_R*drho_dx - dp_dx;
                        Lambda_inv_L[2] = dv_dx;
                        Lambda_inv_L[3] = (Real(1)/lambda_4)*(K*(p_x_R - p_t) - (Real(1) - beta)*T_4);
                        
                        // Compute dV_dx.
                        
                        const Real c_sq_inv  = Real(1)/(c_x_R*c_x_R);
                        const Real rho_c_inv = Real(1)/(rho_x_R*c_x_R);
                        
                        Real dV_dx[4];
                        
                        dV_dx[0] = half*c_sq_inv*(Lambda_inv_L[0] + Lambda_inv_L[3]) + c_sq_inv*Lambda_inv_L[1];
                        dV_dx[1] = half*rho_c_inv*(-Lambda_inv_L[0] + Lambda_inv_L[3]);
                        dV_dx[2] = Lambda_inv_L[2];
                        dV_dx[3] = half*(Lambda_inv_L[0] + Lambda_inv_L[3]);
                        
                        Real V_ghost[4*num_ghosts_to_fill];
                        
                        for (int i = num_ghosts_to_fill - 1; i >= 0; i--)
                        {
                            const int idx_cell_rho = (i + fill_box_lo_idx[0] + num_subghosts_conservative_var[0][0]) +
                                (j + num_subghosts_conservative_var[0][1])*
                                    subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_mom = (i + fill_box_lo_idx[0] + num_subghosts_conservative_var[1][0]) +
                                (j + num_subghosts_conservative_var[1][1])*
                                    subghostcell_dims_conservative_var[1][0];
                            
                            const int idx_cell_E = (i + fill_box_lo_idx[0] + num_subghosts_conservative_var[2][0]) +
                                (j + num_subghosts_conservative_var[2][1])*
                                    subghostcell_dims_conservative_var[2][0];
                            
                            if (i == num_ghosts_to_fill - 1)
                            {
                                V_ghost[i*4 + 0] = rho_x_RR - Real(2)*Real(dx[0])*dV_dx[0];
                                V_ghost[i*4 + 1] = u_x_RR   - Real(2)*Real(dx[0])*dV_dx[1];
                                V_ghost[i*4 + 2] = v_x_RR   - Real(2)*Real(dx[0])*dV_dx[2];
                                V_ghost[i*4 + 3] = p_x_RR   - Real(2)*Real(dx[0])*dV_dx[3];
                            }
                            else if (i == num_ghosts_to_fill - 2)
                            {
                                V_ghost[i*4 + 0] = -Real(2)*rho_x_RR - Real(3)*rho_x_R +
                                    Real(6)*V_ghost[(i + 1)*4 + 0] + Real(6)*Real(dx[0])*dV_dx[0];
                                
                                V_ghost[i*4 + 1] = -Real(2)*u_x_RR - Real(3)*u_x_R +
                                    Real(6)*V_ghost[(i + 1)*4 + 1] + Real(6)*Real(dx[0])*dV_dx[1];
                                
                                V_ghost[i*4 + 2] = -Real(2)*v_x_RR - Real(3)*v_x_R +
                                    Real(6)*V_ghost[(i + 1)*4 + 2] + Real(6)*Real(dx[0])*dV_dx[2];
                                
                                V_ghost[i*4 + 3] = -Real(2)*p_x_RR - Real(3)*p_x_R +
                                    Real(6)*V_ghost[(i + 1)*4 + 3] + Real(6)*Real(dx[0])*dV_dx[3];
                            }
                            else if (i == num_ghosts_to_fill - 3)
                            {
                                V_ghost[i*4 + 0] = Real(3)*rho_x_RR + Real(10)*rho_x_R -
                                    Real(18)*V_ghost[(i + 2)*4 + 0] + Real(6)*V_ghost[(i + 1)*4 + 0] -
                                    Real(12)*Real(dx[0])*dV_dx[0];
                                
                                V_ghost[i*4 + 1] = Real(3)*u_x_RR + Real(10)*u_x_R -
                                    Real(18)*V_ghost[(i + 2)*4 + 1] + Real(6)*V_ghost[(i + 1)*4 + 1] -
                                    Real(12)*Real(dx[0])*dV_dx[1];
                                
                                V_ghost[i*4 + 2] = Real(3)*v_x_RR + Real(10)*v_x_R -
                                    Real(18)*V_ghost[(i + 2)*4 + 2] + Real(6)*V_ghost[(i + 1)*4 + 2] -
                                    Real(12)*Real(dx[0])*dV_dx[2];
                                
                                V_ghost[i*4 + 3] = Real(3)*p_x_RR + Real(10)*p_x_R -
                                    Real(18)*V_ghost[(i + 2)*4 + 3] + Real(6)*V_ghost[(i + 1)*4 + 3] -
                                    Real(12)*Real(dx[0])*dV_dx[3];
                            }
                            else if (i == num_ghosts_to_fill - 4)
                            {
                                V_ghost[i*4 + 0] = -Real(4)*rho_x_RR - Real(65)/Real(3)*rho_x_R +
                                    Real(40)*V_ghost[(i + 3)*4 + 0] - Real(20)*V_ghost[(i + 2)*4 + 0] +
                                    Real(20)/Real(3)*V_ghost[(i + 1)*4 + 0] + Real(20)*Real(dx[0])*dV_dx[0];
                                
                                V_ghost[i*4 + 1] = -Real(4)*u_x_RR - Real(65)/Real(3)*u_x_R +
                                    Real(40)*V_ghost[(i + 3)*4 + 1] - Real(20)*V_ghost[(i + 2)*4 + 1] +
                                    Real(20)/Real(3)*V_ghost[(i + 1)*4 + 1] + Real(20)*Real(dx[0])*dV_dx[1];
                                
                                V_ghost[i*4 + 2] = -Real(4)*v_x_RR - Real(65)/Real(3)*v_x_R +
                                    Real(40)*V_ghost[(i + 3)*4 + 2] - Real(20)*V_ghost[(i + 2)*4 + 2] +
                                    Real(20)/Real(3)*V_ghost[(i + 1)*4 + 2] + Real(20)*Real(dx[0])*dV_dx[2];
                                
                                V_ghost[i*4 + 3] = -Real(4)*p_x_RR - Real(65)/Real(3)*p_x_R +
                                    Real(40)*V_ghost[(i + 3)*4 + 3] - Real(20)*V_ghost[(i + 2)*4 + 3] +
                                    Real(20)/Real(3)*V_ghost[(i + 1)*4 + 3] + Real(20)*Real(dx[0])*dV_dx[3];
                            }
                            else if (i == num_ghosts_to_fill - 5)
                            {
                                V_ghost[i*4 + 0] = Real(5)*rho_x_RR + Real(77)/Real(2)*rho_x_R -
                                    Real(75)*V_ghost[(i + 4)*4 + 0] + Real(50)*V_ghost[(i + 3)*4 + 0] -
                                    Real(25)*V_ghost[(i + 2)*4 + 0] + Real(15)/Real(2)*V_ghost[(i + 1)*4 + 0] -
                                    Real(30)*Real(dx[0])*dV_dx[0];
                                
                                V_ghost[i*4 + 1] = Real(5)*u_x_RR + Real(77)/Real(2)*u_x_R -
                                    Real(75)*V_ghost[(i + 4)*4 + 1] + Real(50)*V_ghost[(i + 3)*4 + 1] -
                                    Real(25)*V_ghost[(i + 2)*4 + 1] + Real(15)/Real(2)*V_ghost[(i + 1)*4 + 1] -
                                    Real(30)*Real(dx[0])*dV_dx[1];
                                
                                V_ghost[i*4 + 2] = Real(5)*v_x_RR + Real(77)/Real(2)*v_x_R -
                                    Real(75)*V_ghost[(i + 4)*4 + 2] + Real(50)*V_ghost[(i + 3)*4 + 2] -
                                    Real(25)*V_ghost[(i + 2)*4 + 2] + Real(15)/Real(2)*V_ghost[(i + 1)*4 + 2] -
                                    Real(30)*Real(dx[0])*dV_dx[2];
                                
                                V_ghost[i*4 + 3] = Real(5)*p_x_RR + Real(77)/Real(2)*p_x_R -
                                    Real(75)*V_ghost[(i + 4)*4 + 3] + Real(50)*V_ghost[(i + 3)*4 + 3] -
                                    Real(25)*V_ghost[(i + 2)*4 + 3] + Real(15)/Real(2)*V_ghost[(i + 1)*4 + 3] -
                                    Real(30)*Real(dx[0])*dV_dx[3];
                            }
                            else if (i == num_ghosts_to_fill - 6)
                            {
                                V_ghost[i*4 + 0] = -Real(6)*rho_x_RR - Real(609)/Real(10)*rho_x_R +
                                    Real(126)*V_ghost[(i + 5)*4 + 0] - Real(105)*V_ghost[(i + 4)*4 + 0] +
                                    Real(70)*V_ghost[(i + 3)*4 + 0] - Real(63)/Real(2)*V_ghost[(i + 2)*4 + 0] +
                                    Real(42)/Real(5)*V_ghost[(i + 1)*4 + 0] + Real(42)*Real(dx[0])*dV_dx[0];
                                
                                V_ghost[i*4 + 1] = -Real(6)*u_x_RR - Real(609)/Real(10)*u_x_R +
                                    Real(126)*V_ghost[(i + 5)*4 + 1] - Real(105)*V_ghost[(i + 4)*4 + 1] +
                                    Real(70)*V_ghost[(i + 3)*4 + 1] - Real(63)/Real(2)*V_ghost[(i + 2)*4 + 1] +
                                    Real(42)/Real(5)*V_ghost[(i + 1)*4 + 1] + Real(42)*Real(dx[0])*dV_dx[1];
                                
                                V_ghost[i*4 + 2] = -Real(6)*v_x_RR - Real(609)/Real(10)*v_x_R +
                                    Real(126)*V_ghost[(i + 5)*4 + 2] - Real(105)*V_ghost[(i + 4)*4 + 2] +
                                    Real(70)*V_ghost[(i + 3)*4 + 2] - Real(63)/Real(2)*V_ghost[(i + 2)*4 + 2] +
                                    Real(42)/Real(5)*V_ghost[(i + 1)*4 + 2] + Real(42)*Real(dx[0])*dV_dx[2];
                                
                                V_ghost[i*4 + 3] = -Real(6)*p_x_RR - Real(609)/Real(10)*p_x_R +
                                    Real(126)*V_ghost[(i + 5)*4 + 3] - Real(105)*V_ghost[(i + 4)*4 + 3] +
                                    Real(70)*V_ghost[(i + 3)*4 + 3] - Real(63)/Real(2)*V_ghost[(i + 2)*4 + 3] +
                                    Real(42)/Real(5)*V_ghost[(i + 1)*4 + 3] + Real(42)*Real(dx[0])*dV_dx[3];
                            }
                            
                            Q[0][idx_cell_rho] = V_ghost[i*4 + 0];
                            Q[1][idx_cell_mom] = V_ghost[i*4 + 0]*V_ghost[i*4 + 1];
                            Q[2][idx_cell_mom] = V_ghost[i*4 + 0]*V_ghost[i*4 + 2];
                            
                            const Real epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getInternalEnergy(
                                    &V_ghost[i*4 + 0],
                                    &V_ghost[i*4 + 3],
                                    thermo_properties_ptr);
                            
                            const Real E = V_ghost[i*4 + 0]*epsilon +
                                half*(Q[1][idx_cell_mom]*Q[1][idx_cell_mom] + Q[2][idx_cell_mom]*Q[2][idx_cell_mom])/
                                    V_ghost[i*4 + 0];
                            
                            Q[3][idx_cell_E] = E;
                        }
                    }
                }
                else if (edge_loc == BDRY_LOC::XHI)
                {
                    const int num_ghosts_to_fill = fill_box_hi_idx[0] - fill_box_lo_idx[0] + 1;
                    TBOX_ASSERT(fill_box_lo_idx[0] == interior_box_hi_idx[0] + 1);
                    if (num_ghosts_to_fill > 6)
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelBoundaryUtilitiesSingleSpecies::fill2dEdgeBoundaryData()\n"
                            << "Non-reflecting outflow BC doesn't support more than six ghost cells yet!");
                    }
                    
                    for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                    {
                        // Get the grid spacing.
                        const double* const dx = patch_geom->getDx();
                        
                        // Compute one-sided derivatives in x-direction.
                        const int idx_cell_rho_x_L = (interior_box_hi_idx[0] + num_subghosts_conservative_var[0][0]) +
                            (j + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0];
                        
                        const int idx_cell_rho_x_LL = (interior_box_hi_idx[0] - 1 + num_subghosts_conservative_var[0][0]) +
                            (j + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0];
                        
                        const int idx_cell_rho_x_LLL = (interior_box_hi_idx[0] - 2 + num_subghosts_conservative_var[0][0]) +
                            (j + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0];
                        
                        const int idx_cell_mom_x_L = (interior_box_hi_idx[0] + num_subghosts_conservative_var[1][0]) +
                            (j + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0];
                        
                        const int idx_cell_mom_x_LL = (interior_box_hi_idx[0] - 1 + num_subghosts_conservative_var[1][0]) +
                            (j + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0];
                        
                        const int idx_cell_mom_x_LLL = (interior_box_hi_idx[0] - 2 + num_subghosts_conservative_var[1][0]) +
                            (j + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0];
                        
                        const int idx_cell_E_x_L = (interior_box_hi_idx[0] + num_subghosts_conservative_var[2][0]) +
                            (j + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0];
                        
                        const int idx_cell_E_x_LL = (interior_box_hi_idx[0] - 1 + num_subghosts_conservative_var[2][0]) +
                            (j + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0];
                        
                        const int idx_cell_E_x_LLL = (interior_box_hi_idx[0] - 2 + num_subghosts_conservative_var[2][0]) +
                            (j + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0];
                        
                        const Real& rho_x_L   = Q[0][idx_cell_rho_x_L];
                        const Real& rho_x_LL  = Q[0][idx_cell_rho_x_LL];
                        const Real& rho_x_LLL = Q[0][idx_cell_rho_x_LLL];
                        
                        const Real u_x_L   = Q[1][idx_cell_mom_x_L]/rho_x_L;
                        const Real u_x_LL  = Q[1][idx_cell_mom_x_LL]/rho_x_LL;
                        const Real u_x_LLL = Q[1][idx_cell_mom_x_LLL]/rho_x_LLL;
                        
                        const Real v_x_L   = Q[2][idx_cell_mom_x_L]/rho_x_L;
                        const Real v_x_LL  = Q[2][idx_cell_mom_x_LL]/rho_x_LL;
                        const Real v_x_LLL = Q[2][idx_cell_mom_x_LLL]/rho_x_LLL;
                       
                        const Real half = Real(1)/Real(2);
                        const Real epsilon_x_L   = Q[3][idx_cell_E_x_L]/rho_x_L - half*(u_x_L*u_x_L + v_x_L*v_x_L);
                        const Real epsilon_x_LL  = Q[3][idx_cell_E_x_LL]/rho_x_LL - half*(u_x_LL*u_x_LL + v_x_LL*v_x_LL);
                        const Real epsilon_x_LLL = Q[3][idx_cell_E_x_LLL]/rho_x_LLL - half*(u_x_LLL*u_x_LLL + v_x_LLL*v_x_LLL);
                        
                        const Real p_x_L = d_equation_of_state_mixing_rules->getEquationOfState()->
                            getPressure(
                                &rho_x_L,
                                &epsilon_x_L,
                                thermo_properties_ptr);
                        
                        const Real p_x_LL = d_equation_of_state_mixing_rules->getEquationOfState()->
                            getPressure(
                                &rho_x_LL,
                                &epsilon_x_LL,
                                thermo_properties_ptr);
                        
                        const Real p_x_LLL = d_equation_of_state_mixing_rules->getEquationOfState()->
                            getPressure(
                                &rho_x_LLL,
                                &epsilon_x_LLL,
                                thermo_properties_ptr);
                        
                        /*
                         * Compute derivatives in x-direction.
                         */
                        
                        const Real drho_dx = (rho_x_LLL - Real(4)*rho_x_LL + Real(3)*rho_x_L)/(Real(2)*Real(dx[0]));
                        const Real du_dx   = (u_x_LLL - Real(4)*u_x_LL + Real(3)*u_x_L)/(Real(2)*Real(dx[0]));
                        const Real dv_dx   = (v_x_LLL - Real(4)*v_x_LL + Real(3)*v_x_L)/(Real(2)*Real(dx[0]));
                        const Real dp_dx   = (p_x_LLL - Real(4)*p_x_LL + Real(3)*p_x_L)/(Real(2)*Real(dx[0]));
                        
                        /*
                         * Compute derivatives in y-direction.
                         */
                        
                        Real du_dy = Real(0);
                        Real dv_dy = Real(0);
                        Real dp_dy = Real(0);
                        
                        if ((j + num_subghosts_conservative_var[0][1] == 0) ||
                            (j + num_subghosts_conservative_var[1][1] == 0) ||
                            (j + num_subghosts_conservative_var[2][1] == 0))
                        // if (((patch_geom->getTouchesRegularBoundary(1, 0)) && (j == interior_box_lo_idx[1])) ||
                        //     ((j + num_subghosts_conservative_var[0][1] == 0) ||
                        //      (j + num_subghosts_conservative_var[1][1] == 0) ||
                        //      (j + num_subghosts_conservative_var[2][1] == 0)))
                        {
                            // Patch is touching periodic boundary.
                            
                            const int idx_cell_rho_y_T = (interior_box_hi_idx[0] + num_subghosts_conservative_var[0][0]) +
                                (j + 1 + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_mom_y_T = (interior_box_hi_idx[0] + num_subghosts_conservative_var[1][0]) +
                                (j + 1 + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0];
                            
                            const int idx_cell_E_y_T = (interior_box_hi_idx[0] + num_subghosts_conservative_var[2][0]) +
                                (j + 1 + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0];
                            
                            const Real& rho_y_T = Q[0][idx_cell_rho_y_T];
                            const Real u_y_T = Q[1][idx_cell_mom_y_T]/rho_y_T;
                            const Real v_y_T = Q[2][idx_cell_mom_y_T]/rho_y_T;
                            const Real epsilon_y_T = Q[3][idx_cell_E_y_T]/rho_y_T - half*(u_y_T*u_y_T + v_y_T*v_y_T);
                            
                            const Real p_y_T = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getPressure(
                                    &rho_y_T,
                                    &epsilon_y_T,
                                    thermo_properties_ptr);
                            
                            // One-sided derivatives.
                            du_dy = (u_y_T - u_x_L)/Real(dx[1]);
                            dv_dy = (v_y_T - v_x_L)/Real(dx[1]);
                            dp_dy = (p_y_T - p_x_L)/Real(dx[1]);
                        }
                        else if ((j + num_subghosts_conservative_var[0][1] + 1 == subghostcell_dims_conservative_var[0][1]) ||
                                 (j + num_subghosts_conservative_var[1][1] + 1 == subghostcell_dims_conservative_var[1][1]) ||
                                 (j + num_subghosts_conservative_var[2][1] + 1 == subghostcell_dims_conservative_var[2][1]))
                        // else if (((patch_geom->getTouchesRegularBoundary(1, 1)) && (j == interior_box_hi_idx[1])) ||
                        //          ((j + num_subghosts_conservative_var[0][1] + 1 == subghostcell_dims_conservative_var[0][1]) ||
                        //           (j + num_subghosts_conservative_var[1][1] + 1 == subghostcell_dims_conservative_var[1][1]) ||
                        //           (j + num_subghosts_conservative_var[2][1] + 1 == subghostcell_dims_conservative_var[2][1])))
                        {
                            // Patch is touching periodic boundary.
                            
                            const int idx_cell_rho_y_B = (interior_box_hi_idx[0] + num_subghosts_conservative_var[0][0]) +
                                (j - 1 + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_mom_y_B = (interior_box_hi_idx[0] + num_subghosts_conservative_var[1][0]) +
                                (j - 1 + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0];
                            
                            const int idx_cell_E_y_B = (interior_box_hi_idx[0] + num_subghosts_conservative_var[2][0]) +
                                (j - 1 + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0];
                            
                            const Real& rho_y_B = Q[0][idx_cell_rho_y_B];
                            const Real u_y_B = Q[1][idx_cell_mom_y_B]/rho_y_B;
                            const Real v_y_B = Q[2][idx_cell_mom_y_B]/rho_y_B;
                            const Real epsilon_y_B = Q[3][idx_cell_E_y_B]/rho_y_B - half*(u_y_B*u_y_B + v_y_B*v_y_B);
                            
                            const Real p_y_B = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getPressure(
                                    &rho_y_B,
                                    &epsilon_y_B,
                                    thermo_properties_ptr);
                            
                            // One-sided derivatives.
                            du_dy = (u_x_L - u_y_B)/Real(dx[1]);
                            dv_dy = (v_x_L - v_y_B)/Real(dx[1]);
                            dp_dy = (p_x_L - p_y_B)/Real(dx[1]);
                        }
                        else
                        {
                            const int idx_cell_rho_y_B = (interior_box_hi_idx[0] + num_subghosts_conservative_var[0][0]) +
                                (j - 1 + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_rho_y_T = (interior_box_hi_idx[0] + num_subghosts_conservative_var[0][0]) +
                                (j + 1 + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_mom_y_B = (interior_box_hi_idx[0] + num_subghosts_conservative_var[1][0]) +
                                (j - 1 + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0];
                            
                            const int idx_cell_mom_y_T = (interior_box_hi_idx[0] + num_subghosts_conservative_var[1][0]) +
                                (j + 1 + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0];
                            
                            const int idx_cell_E_y_B = (interior_box_hi_idx[0] + num_subghosts_conservative_var[2][0]) +
                                (j - 1 + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0];
                            
                            const int idx_cell_E_y_T = (interior_box_hi_idx[0] + num_subghosts_conservative_var[2][0]) +
                                (j + 1 + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0];
                            
                            const Real& rho_y_B = Q[0][idx_cell_rho_y_B];
                            const Real& rho_y_T = Q[0][idx_cell_rho_y_T];
                            
                            const Real u_y_B = Q[1][idx_cell_mom_y_B]/rho_y_B;
                            const Real u_y_T = Q[1][idx_cell_mom_y_T]/rho_y_T;
                            
                            const Real v_y_B = Q[2][idx_cell_mom_y_B]/rho_y_B;
                            const Real v_y_T = Q[2][idx_cell_mom_y_T]/rho_y_T;
                            
                            const Real epsilon_y_B = Q[3][idx_cell_E_y_B]/rho_y_B - half*(u_y_B*u_y_B + v_y_B*v_y_B);
                            const Real epsilon_y_T = Q[3][idx_cell_E_y_T]/rho_y_T - half*(u_y_T*u_y_T + v_y_T*v_y_T);
                            
                            const Real p_y_B = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getPressure(
                                    &rho_y_B,
                                    &epsilon_y_B,
                                    thermo_properties_ptr);
                            
                            const Real p_y_T = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getPressure(
                                    &rho_y_T,
                                    &epsilon_y_T,
                                    thermo_properties_ptr);
                            
                            // Central derivatives.
                            du_dy = (u_y_T - u_y_B)/(Real(2)*Real(dx[1]));
                            dv_dy = (v_y_T - v_y_B)/(Real(2)*Real(dx[1]));
                            dp_dy = (p_y_T - p_y_B)/(Real(2)*Real(dx[1]));
                        }
                        
                        // Compute wave speed (u - c) at the boundary.
                        
                        const Real c_x_L = d_equation_of_state_mixing_rules->getEquationOfState()->
                            getSoundSpeed(
                                &rho_x_L,
                                &p_x_L,
                                thermo_properties_ptr);
                        
                        const Real lambda_1 = u_x_L - c_x_L;
                        
                        // Compute vector Lambda^(-1) * L.
                        
                        Real Lambda_inv_L[4];
                        
                        const Real& p_t         = d_bdry_edge_nonreflecting_outflow_p_t[edge_loc];
                        const Real& sigma       = d_bdry_edge_nonreflecting_outflow_sigma[edge_loc];
                        const Real& beta        = d_bdry_edge_nonreflecting_outflow_beta[edge_loc];
                        const Real& length_char = d_bdry_edge_nonreflecting_outflow_length_char[edge_loc];
                        
                        const Real T_1 = v_x_L*(dp_dy - rho_x_L*c_x_L*du_dy) + rho_x_L*c_x_L*c_x_L*dv_dy;
                        
                        const Real M_sq = (u_x_L*u_x_L + v_x_L*v_x_L)/(c_x_L*c_x_L);
                        const Real K = sigma*c_x_L*(Real(1) - M_sq)/length_char;
                        
                        Lambda_inv_L[0] = (Real(1)/lambda_1)*(K*(p_x_L - p_t) - (Real(1) - beta)*T_1);
                        Lambda_inv_L[1] = c_x_L*c_x_L*drho_dx - dp_dx;
                        Lambda_inv_L[2] = dv_dx;
                        Lambda_inv_L[3] = dp_dx + rho_x_L*c_x_L*du_dx;
                        
                        // Compute dV_dx.
                        
                        const Real c_sq_inv  = Real(1)/(c_x_L*c_x_L);
                        const Real rho_c_inv = Real(1)/(rho_x_L*c_x_L);
                        
                        Real dV_dx[4];
                        
                        dV_dx[0] = half*c_sq_inv*(Lambda_inv_L[0] + Lambda_inv_L[3]) + c_sq_inv*Lambda_inv_L[1];
                        dV_dx[1] = half*rho_c_inv*(-Lambda_inv_L[0] + Lambda_inv_L[3]);
                        dV_dx[2] = Lambda_inv_L[2];
                        dV_dx[3] = half*(Lambda_inv_L[0] + Lambda_inv_L[3]);
                        
                        Real V_ghost[4*num_ghosts_to_fill];
                        
                        for (int i = 0; i < num_ghosts_to_fill; i++)
                        {
                            const int idx_cell_rho = (i + fill_box_lo_idx[0] + num_subghosts_conservative_var[0][0]) +
                                (j + num_subghosts_conservative_var[0][1])*
                                    subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_mom = (i + fill_box_lo_idx[0] + num_subghosts_conservative_var[1][0]) +
                                (j + num_subghosts_conservative_var[1][1])*
                                    subghostcell_dims_conservative_var[1][0];
                            
                            const int idx_cell_E = (i + fill_box_lo_idx[0] + num_subghosts_conservative_var[2][0]) +
                                (j + num_subghosts_conservative_var[2][1])*
                                    subghostcell_dims_conservative_var[2][0];
                            
                            if (i == 0)
                            {
                                V_ghost[i*4 + 0] = rho_x_LL + Real(2)*Real(dx[0])*dV_dx[0];
                                V_ghost[i*4 + 1] = u_x_LL   + Real(2)*Real(dx[0])*dV_dx[1];
                                V_ghost[i*4 + 2] = v_x_LL   + Real(2)*Real(dx[0])*dV_dx[2];
                                V_ghost[i*4 + 3] = p_x_LL   + Real(2)*Real(dx[0])*dV_dx[3];
                            }
                            else if (i == 1)
                            {
                                V_ghost[i*4 + 0] = -Real(2)*rho_x_LL - Real(3)*rho_x_L +
                                    Real(6)*V_ghost[(i - 1)*4 + 0] - Real(6)*Real(dx[0])*dV_dx[0];
                                
                                V_ghost[i*4 + 1] = -Real(2)*u_x_LL - Real(3)*u_x_L +
                                    Real(6)*V_ghost[(i - 1)*4 + 1] - Real(6)*Real(dx[0])*dV_dx[1];
                                
                                V_ghost[i*4 + 2] = -Real(2)*v_x_LL - Real(3)*v_x_L +
                                    Real(6)*V_ghost[(i - 1)*4 + 2] - Real(6)*Real(dx[0])*dV_dx[2];
                                
                                V_ghost[i*4 + 3] = -Real(2)*p_x_LL - Real(3)*p_x_L +
                                    Real(6)*V_ghost[(i - 1)*4 + 3] - Real(6)*Real(dx[0])*dV_dx[3];
                            }
                            else if (i == 2)
                            {
                                V_ghost[i*4 + 0] = Real(3)*rho_x_LL + Real(10)*rho_x_L -
                                    Real(18)*V_ghost[(i - 2)*4 + 0] + Real(6)*V_ghost[(i - 1)*4 + 0] +
                                    Real(12)*Real(dx[0])*dV_dx[0];
                                
                                V_ghost[i*4 + 1] = Real(3)*u_x_LL + Real(10)*u_x_L -
                                    Real(18)*V_ghost[(i - 2)*4 + 1] + Real(6)*V_ghost[(i - 1)*4 + 1] +
                                    Real(12)*Real(dx[0])*dV_dx[1];
                                
                                V_ghost[i*4 + 2] = Real(3)*v_x_LL + Real(10)*v_x_L -
                                    Real(18)*V_ghost[(i - 2)*4 + 2] + Real(6)*V_ghost[(i - 1)*4 + 2] +
                                    Real(12)*Real(dx[0])*dV_dx[2];
                                
                                V_ghost[i*4 + 3] = Real(3)*p_x_LL + Real(10)*p_x_L -
                                    Real(18)*V_ghost[(i - 2)*4 + 3] + Real(6)*V_ghost[(i - 1)*4 + 3] +
                                    Real(12)*Real(dx[0])*dV_dx[3];
                            }
                            else if (i == 3)
                            {
                                V_ghost[i*4 + 0] = -Real(4)*rho_x_LL - Real(65)/Real(3)*rho_x_L +
                                    Real(40)*V_ghost[(i - 3)*4 + 0] - Real(20)*V_ghost[(i - 2)*4 + 0] +
                                    Real(20)/Real(3)*V_ghost[(i - 1)*4 + 0] - Real(20)*Real(dx[0])*dV_dx[0];
                                
                                V_ghost[i*4 + 1] = -Real(4)*u_x_LL - Real(65)/Real(3)*u_x_L +
                                    Real(40)*V_ghost[(i - 3)*4 + 1] - Real(20)*V_ghost[(i - 2)*4 + 1] +
                                    Real(20)/Real(3)*V_ghost[(i - 1)*4 + 1] - Real(20)*Real(dx[0])*dV_dx[1];
                                
                                V_ghost[i*4 + 2] = -Real(4)*v_x_LL - Real(65)/Real(3)*v_x_L +
                                    Real(40)*V_ghost[(i - 3)*4 + 2] - Real(20)*V_ghost[(i - 2)*4 + 2] +
                                    Real(20)/Real(3)*V_ghost[(i - 1)*4 + 2] - Real(20)*Real(dx[0])*dV_dx[2];
                                
                                V_ghost[i*4 + 3] = -Real(4)*p_x_LL - Real(65)/Real(3)*p_x_L +
                                    Real(40)*V_ghost[(i - 3)*4 + 3] - Real(20)*V_ghost[(i - 2)*4 + 3] +
                                    Real(20)/Real(3)*V_ghost[(i - 1)*4 + 3] - Real(20)*Real(dx[0])*dV_dx[3];
                            }
                            else if (i == 4)
                            {
                                V_ghost[i*4 + 0] = Real(5)*rho_x_LL + Real(77)/Real(2)*rho_x_L -
                                    Real(75)*V_ghost[(i - 4)*4 + 0] + Real(50)*V_ghost[(i - 3)*4 + 0] -
                                    Real(25)*V_ghost[(i - 2)*4 + 0] + Real(15)/Real(2)*V_ghost[(i - 1)*4 + 0] +
                                    Real(30)*Real(dx[0])*dV_dx[0];
                                
                                V_ghost[i*4 + 1] = Real(5)*u_x_LL + Real(77)/Real(2)*u_x_L -
                                    Real(75)*V_ghost[(i - 4)*4 + 1] + Real(50)*V_ghost[(i - 3)*4 + 1] -
                                    Real(25)*V_ghost[(i - 2)*4 + 1] + Real(15)/Real(2)*V_ghost[(i - 1)*4 + 1] +
                                    Real(30)*Real(dx[0])*dV_dx[1];
                                
                                V_ghost[i*4 + 2] = Real(5)*v_x_LL + Real(77)/Real(2)*v_x_L -
                                    Real(75)*V_ghost[(i - 4)*4 + 2] + Real(50)*V_ghost[(i - 3)*4 + 2] -
                                    Real(25)*V_ghost[(i - 2)*4 + 2] + Real(15)/Real(2)*V_ghost[(i - 1)*4 + 2] +
                                    Real(30)*Real(dx[0])*dV_dx[2];
                                
                                V_ghost[i*4 + 3] = Real(5)*p_x_LL + Real(77)/Real(2)*p_x_L -
                                    Real(75)*V_ghost[(i - 4)*4 + 3] + Real(50)*V_ghost[(i - 3)*4 + 3] -
                                    Real(25)*V_ghost[(i - 2)*4 + 3] + Real(15)/Real(2)*V_ghost[(i - 1)*4 + 3] +
                                    Real(30)*Real(dx[0])*dV_dx[3];
                            }
                            else if (i == 5)
                            {
                                V_ghost[i*4 + 0] = -Real(6)*rho_x_LL - Real(609)/Real(10)*rho_x_L +
                                    Real(126)*V_ghost[(i - 5)*4 + 0] - Real(105)*V_ghost[(i - 4)*4 + 0] +
                                    Real(70)*V_ghost[(i - 3)*4 + 0] - Real(63)/Real(2)*V_ghost[(i - 2)*4 + 0] +
                                    Real(42)/Real(5)*V_ghost[(i - 1)*4 + 0] - Real(42)*Real(dx[0])*dV_dx[0];
                                
                                V_ghost[i*4 + 1] = -Real(6)*u_x_LL - Real(609)/Real(10)*u_x_L +
                                    Real(126)*V_ghost[(i - 5)*4 + 1] - Real(105)*V_ghost[(i - 4)*4 + 1] +
                                    Real(70)*V_ghost[(i - 3)*4 + 1] - Real(63)/Real(2)*V_ghost[(i - 2)*4 + 1] +
                                    Real(42)/Real(5)*V_ghost[(i - 1)*4 + 1] - Real(42)*Real(dx[0])*dV_dx[1];
                                
                                V_ghost[i*4 + 2] = -Real(6)*v_x_LL - Real(609)/Real(10)*v_x_L +
                                    Real(126)*V_ghost[(i - 5)*4 + 2] - Real(105)*V_ghost[(i - 4)*4 + 2] +
                                    Real(70)*V_ghost[(i - 3)*4 + 2] - Real(63)/Real(2)*V_ghost[(i - 2)*4 + 2] +
                                    Real(42)/Real(5)*V_ghost[(i - 1)*4 + 2] - Real(42)*Real(dx[0])*dV_dx[2];
                                
                                V_ghost[i*4 + 3] = -Real(6)*p_x_LL - Real(609)/Real(10)*p_x_L +
                                    Real(126)*V_ghost[(i - 5)*4 + 3] - Real(105)*V_ghost[(i - 4)*4 + 3] +
                                    Real(70)*V_ghost[(i - 3)*4 + 3] - Real(63)/Real(2)*V_ghost[(i - 2)*4 + 3] +
                                    Real(42)/Real(5)*V_ghost[(i - 1)*4 + 3] - Real(42)*Real(dx[0])*dV_dx[3];
                            }
                            
                            Q[0][idx_cell_rho] = V_ghost[i*4 + 0];
                            Q[1][idx_cell_mom] = V_ghost[i*4 + 0]*V_ghost[i*4 + 1];
                            Q[2][idx_cell_mom] = V_ghost[i*4 + 0]*V_ghost[i*4 + 2];
                            
                            const Real epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getInternalEnergy(
                                    &V_ghost[i*4 + 0],
                                    &V_ghost[i*4 + 3],
                                    thermo_properties_ptr);
                            
                            const Real E = V_ghost[i*4 + 0]*epsilon +
                                half*(Q[1][idx_cell_mom]*Q[1][idx_cell_mom] + Q[2][idx_cell_mom]*Q[2][idx_cell_mom])/
                                    V_ghost[i*4 + 0];
                            
                            Q[3][idx_cell_E] = E;
                        }
                    }
                }
                else if (edge_loc == BDRY_LOC::YLO)
                {
                    const int num_ghosts_to_fill = fill_box_hi_idx[1] - fill_box_lo_idx[1] + 1;
                    TBOX_ASSERT(fill_box_hi_idx[1] == interior_box_lo_idx[1] - 1);
                    if (num_ghosts_to_fill > 6)
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelBoundaryUtilitiesSingleSpecies::fill2dEdgeBoundaryData()\n"
                            << "Non-reflecting outflow BC doesn't support more than six ghost cells yet!");
                    }
                    
                    for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                    {
                        // Get the grid spacing.
                        const double* const dx = patch_geom->getDx();
                        
                        const int idx_cell_rho_y_T = (i + num_subghosts_conservative_var[0][0]) +
                            (interior_box_lo_idx[1] + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0];
                        
                        const int idx_cell_rho_y_TT = (i + num_subghosts_conservative_var[0][0]) +
                            (interior_box_lo_idx[1] + 1 + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0];
                        
                        const int idx_cell_rho_y_TTT = (i + num_subghosts_conservative_var[0][0])+
                            (interior_box_lo_idx[1] + 2 + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0];
                        
                        const int idx_cell_mom_y_T = (i + num_subghosts_conservative_var[1][0]) +
                            (interior_box_lo_idx[1] + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0];
                        
                        const int idx_cell_mom_y_TT = (i + num_subghosts_conservative_var[1][0]) +
                            (interior_box_lo_idx[1] + 1 + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0];
                        
                        const int idx_cell_mom_y_TTT = (i + num_subghosts_conservative_var[1][0]) +
                            (interior_box_lo_idx[1] + 2 + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0];
                        
                        const int idx_cell_E_y_T = (i + num_subghosts_conservative_var[2][0]) +
                            (interior_box_lo_idx[1] + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0];
                        
                        const int idx_cell_E_y_TT = (i + num_subghosts_conservative_var[2][0]) +
                            (interior_box_lo_idx[1] + 1 + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0];
                        
                        const int idx_cell_E_y_TTT = (i + num_subghosts_conservative_var[2][0]) +
                            (interior_box_lo_idx[1] + 2 + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0];
                        
                        const Real& rho_y_T   = Q[0][idx_cell_rho_y_T];
                        const Real& rho_y_TT  = Q[0][idx_cell_rho_y_TT];
                        const Real& rho_y_TTT = Q[0][idx_cell_rho_y_TTT];
                        
                        const Real u_y_T   = Q[1][idx_cell_mom_y_T]/rho_y_T;
                        const Real u_y_TT  = Q[1][idx_cell_mom_y_TT]/rho_y_TT;
                        const Real u_y_TTT = Q[1][idx_cell_mom_y_TTT]/rho_y_TTT;
                        
                        const Real v_y_T   = Q[2][idx_cell_mom_y_T]/rho_y_T;
                        const Real v_y_TT  = Q[2][idx_cell_mom_y_TT]/rho_y_TT;
                        const Real v_y_TTT = Q[2][idx_cell_mom_y_TTT]/rho_y_TTT;
                        
                        const Real half = Real(1)/Real(2);
                        const Real epsilon_y_T   = Q[3][idx_cell_E_y_T]/rho_y_T - half*(u_y_T*u_y_T + v_y_T*v_y_T);
                        const Real epsilon_y_TT  = Q[3][idx_cell_E_y_TT]/rho_y_TT - half*(u_y_TT*u_y_TT + v_y_TT*v_y_TT);
                        const Real epsilon_y_TTT = Q[3][idx_cell_E_y_TTT]/rho_y_TTT- half*(u_y_TTT*u_y_TTT + v_y_TTT*v_y_TTT);
                        
                        const Real p_y_T = d_equation_of_state_mixing_rules->getEquationOfState()->
                            getPressure(
                                &rho_y_T,
                                &epsilon_y_T,
                                thermo_properties_ptr);
                        
                        const Real p_y_TT = d_equation_of_state_mixing_rules->getEquationOfState()->
                            getPressure(
                                &rho_y_TT,
                                &epsilon_y_TT,
                                thermo_properties_ptr);
                        
                        const Real p_y_TTT = d_equation_of_state_mixing_rules->getEquationOfState()->
                            getPressure(
                                &rho_y_TTT,
                                &epsilon_y_TTT,
                                thermo_properties_ptr);

                        /*
                         * Compute derivatives in y-direction.
                         */
                        
                        const Real drho_dy = -(rho_y_TTT - Real(4)*rho_y_TT + Real(3)*rho_y_T)/(Real(2)*Real(dx[1]));
                        const Real du_dy   = -(u_y_TTT - Real(4)*u_y_TT + Real(3)*u_y_T)/(Real(2)*Real(dx[1]));
                        const Real dv_dy   = -(v_y_TTT - Real(4)*v_y_TT + Real(3)*v_y_T)/(Real(2)*Real(dx[1]));
                        const Real dp_dy   = -(p_y_TTT - Real(4)*p_y_TT + Real(3)*p_y_T)/(Real(2)*Real(dx[1]));
                        
                        /*
                         * Compute derivatives in x-direction.
                         */
                        
                        Real du_dx = Real(0);
                        Real dv_dx = Real(0);
                        Real dp_dx = Real(0);
                        
                        if ((i + num_subghosts_conservative_var[0][0] == 0) ||
                            (i + num_subghosts_conservative_var[1][0] == 0) ||
                            (i + num_subghosts_conservative_var[2][0] == 0))
                        // if (((patch_geom->getTouchesRegularBoundary(0, 0)) && (i == interior_box_lo_idx[0])) ||
                        //     ((i + num_subghosts_conservative_var[0][0] == 0) ||
                        //      (i + num_subghosts_conservative_var[1][0] == 0) ||
                        //      (i + num_subghosts_conservative_var[2][0] == 0)))
                        {
                            // Patch is touching periodic boundary.
                            
                            const int idx_cell_rho_x_R = (i + 1 + num_subghosts_conservative_var[0][0]) +
                                (interior_box_lo_idx[1] + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_mom_x_R = (i + 1 + num_subghosts_conservative_var[1][0]) +
                                (interior_box_lo_idx[1] + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0];
                            
                            const int idx_cell_E_x_R = (i + 1 + num_subghosts_conservative_var[2][0]) +
                                (interior_box_lo_idx[1] + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0];
                            
                            const Real& rho_x_R = Q[0][idx_cell_rho_x_R];
                            const Real u_x_R = Q[1][idx_cell_mom_x_R]/rho_x_R;
                            const Real v_x_R = Q[2][idx_cell_mom_x_R]/rho_x_R;
                            const Real epsilon_x_R = Q[3][idx_cell_E_x_R]/rho_x_R - half*(u_x_R*u_x_R + v_x_R*v_x_R);
                            
                            const Real p_x_R = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getPressure(
                                    &rho_x_R,
                                    &epsilon_x_R,
                                    thermo_properties_ptr);
                            
                            // One-sided derivatives.
                            du_dx = (u_x_R - u_y_T)/Real(dx[0]);
                            dv_dx = (v_x_R - v_y_T)/Real(dx[0]);
                            dp_dx = (p_x_R - p_y_T)/Real(dx[0]);
                        }
                        else if ((i + num_subghosts_conservative_var[0][0] + 1 == subghostcell_dims_conservative_var[0][0]) ||
                                 (i + num_subghosts_conservative_var[1][0] + 1 == subghostcell_dims_conservative_var[1][0]) ||
                                 (i + num_subghosts_conservative_var[2][0] + 1 == subghostcell_dims_conservative_var[2][0]))
                        // else if (((patch_geom->getTouchesRegularBoundary(0, 1)) && (i == interior_box_hi_idx[0])) ||
                        //          ((i + num_subghosts_conservative_var[0][0] + 1 == subghostcell_dims_conservative_var[0][0]) ||
                        //           (i + num_subghosts_conservative_var[1][0] + 1 == subghostcell_dims_conservative_var[1][0]) ||
                        //           (i + num_subghosts_conservative_var[2][0] + 1 == subghostcell_dims_conservative_var[2][0])))
                        {
                            // Patch is touching periodic boundary.
                            
                            const int idx_cell_rho_x_L = (i - 1 + num_subghosts_conservative_var[0][0]) +
                                (interior_box_lo_idx[1] + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_mom_x_L = (i - 1 + num_subghosts_conservative_var[1][0]) +
                                (interior_box_lo_idx[1] + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0];
                            
                            const int idx_cell_E_x_L = (i - 1 + num_subghosts_conservative_var[2][0]) +
                                (interior_box_lo_idx[1] + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[1][0];
                            
                            const Real& rho_x_L = Q[0][idx_cell_rho_x_L];
                            const Real u_x_L = Q[1][idx_cell_mom_x_L]/rho_x_L;
                            const Real v_x_L = Q[2][idx_cell_mom_x_L]/rho_x_L;
                            const Real epsilon_x_L = Q[3][idx_cell_E_x_L]/rho_x_L - half*(u_x_L*u_x_L + v_x_L*v_x_L);
                            
                            const Real p_x_L = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getPressure(
                                    &rho_x_L,
                                    &epsilon_x_L,
                                    thermo_properties_ptr);
                            
                            // One-sided derivatives.
                            du_dx = (u_y_T - u_x_L)/Real(dx[0]);
                            dv_dx = (v_y_T - v_x_L)/Real(dx[0]);
                            dp_dx = (p_y_T - p_x_L)/Real(dx[0]);
                        }
                        else
                        {
                            const int idx_cell_rho_x_L = (i - 1 + num_subghosts_conservative_var[0][0]) +
                                (interior_box_lo_idx[1] + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_rho_x_R = (i + 1 + num_subghosts_conservative_var[0][0]) +
                                (interior_box_lo_idx[1] + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_mom_x_L = (i - 1 + num_subghosts_conservative_var[1][0]) +
                                (interior_box_lo_idx[1] + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0];
                            
                            const int idx_cell_mom_x_R = (i + 1 + num_subghosts_conservative_var[1][0]) +
                                (interior_box_lo_idx[1] + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0];
                            
                            const int idx_cell_E_x_L = (i + 1 + num_subghosts_conservative_var[1][0]) +
                                (interior_box_lo_idx[1] + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0];
                            
                            const int idx_cell_E_x_R = (i + 1 + num_subghosts_conservative_var[2][0]) +
                                (interior_box_lo_idx[1] + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0];
                            
                            const Real& rho_x_L = Q[0][idx_cell_rho_x_L];
                            const Real& rho_x_R = Q[0][idx_cell_rho_x_R];
                            
                            const Real u_x_L = Q[1][idx_cell_mom_x_L]/rho_x_L;
                            const Real u_x_R = Q[1][idx_cell_mom_x_R]/rho_x_R;
                            
                            const Real v_x_L = Q[2][idx_cell_mom_x_L]/rho_x_L;
                            const Real v_x_R = Q[2][idx_cell_mom_x_R]/rho_x_R;
                            
                            const Real epsilon_x_L = Q[3][idx_cell_E_x_L]/rho_x_L - half*(u_x_L*u_x_L + v_x_L*v_x_L);
                            const Real epsilon_x_R = Q[3][idx_cell_E_x_R]/rho_x_R - half*(u_x_R*u_x_R + v_x_R*v_x_R);
                            
                            const Real p_x_L = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getPressure(
                                    &rho_x_L,
                                    &epsilon_x_L,
                                    thermo_properties_ptr);
                            
                            const Real p_x_R = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getPressure(
                                    &rho_x_R,
                                    &epsilon_x_R,
                                    thermo_properties_ptr);
                            
                            // Central derivatives.
                            du_dx = (u_x_R - u_x_L)/(Real(2)*Real(dx[0]));
                            dv_dx = (v_x_R - v_x_L)/(Real(2)*Real(dx[0]));
                            dp_dx = (p_x_R - p_x_L)/(Real(2)*Real(dx[0]));
                        }
                        
                        // Compute wave speed (v + c) at the boundary.
                        
                        const Real c_y_T = d_equation_of_state_mixing_rules->getEquationOfState()->
                            getSoundSpeed(
                                &rho_y_T,
                                &p_y_T,
                                thermo_properties_ptr);
                        
                        const Real lambda_4 = v_y_T + c_y_T;
                        
                        // Compute vector Lambda^(-1) * L.
                        
                        Real Lambda_inv_L[4];
                        
                        const Real& p_t         = d_bdry_edge_nonreflecting_outflow_p_t[edge_loc];
                        const Real& sigma       = d_bdry_edge_nonreflecting_outflow_sigma[edge_loc];
                        const Real& beta        = d_bdry_edge_nonreflecting_outflow_beta[edge_loc];
                        const Real& length_char = d_bdry_edge_nonreflecting_outflow_length_char[edge_loc];
                        
                        const Real T_4 = u_y_T*(dp_dx + rho_y_T*c_y_T*dv_dx) + rho_y_T*c_y_T*c_y_T*du_dx;
                        
                        const Real M_sq = (v_y_T*v_y_T + u_y_T*u_y_T)/(c_y_T*c_y_T);
                        const Real K = sigma*c_y_T*(Real(1) - M_sq)/length_char;
                        
                        Lambda_inv_L[0] = dp_dy - rho_y_T*c_y_T*dv_dy;
                        Lambda_inv_L[1] = du_dy;
                        Lambda_inv_L[2] = c_y_T*c_y_T*drho_dy - dp_dy;
                        Lambda_inv_L[3] = (Real(1)/lambda_4)*(K*(p_y_T - p_t) - (Real(1) - beta)*T_4);
                        
                        // Compute dV_dy.
                        
                        const Real c_sq_inv  = Real(1)/(c_y_T*c_y_T);
                        const Real rho_c_inv = Real(1)/(rho_y_T*c_y_T);
                        
                        Real dV_dy[4];
                        
                        dV_dy[0] = half*c_sq_inv*(Lambda_inv_L[0] + Lambda_inv_L[3]) + c_sq_inv*Lambda_inv_L[2];
                        dV_dy[1] = Lambda_inv_L[1];
                        dV_dy[2] = half*rho_c_inv*(-Lambda_inv_L[0] + Lambda_inv_L[3]);
                        dV_dy[3] = half*(Lambda_inv_L[0] + Lambda_inv_L[3]);
                        
                        Real V_ghost[4*num_ghosts_to_fill];
                        
                        for (int j = num_ghosts_to_fill - 1; j >= 0; j--)
                        {
                            const int idx_cell_rho = (i + num_subghosts_conservative_var[0][0]) +
                                (j + fill_box_lo_idx[1] + num_subghosts_conservative_var[0][1])*
                                    subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_mom = (i + num_subghosts_conservative_var[1][0]) +
                                (j + fill_box_lo_idx[1] + num_subghosts_conservative_var[1][1])*
                                    subghostcell_dims_conservative_var[1][0];
                            
                            const int idx_cell_E = (i + num_subghosts_conservative_var[2][0]) +
                                (j + fill_box_lo_idx[1] + num_subghosts_conservative_var[2][1])*
                                    subghostcell_dims_conservative_var[2][0];
                            
                            if (j == num_ghosts_to_fill - 1)
                            {
                                V_ghost[j*4 + 0] = rho_y_TT - Real(2)*Real(dx[1])*dV_dy[0];
                                V_ghost[j*4 + 1] = u_y_TT   - Real(2)*Real(dx[1])*dV_dy[1];
                                V_ghost[j*4 + 2] = v_y_TT   - Real(2)*Real(dx[1])*dV_dy[2];
                                V_ghost[j*4 + 3] = p_y_TT   - Real(2)*Real(dx[1])*dV_dy[3];
                            }
                            else if (j == num_ghosts_to_fill - 2)
                            {
                                V_ghost[j*4 + 0] = -Real(2)*rho_y_TT - Real(3)*rho_y_T +
                                    Real(6)*V_ghost[(j + 1)*4 + 0] + Real(6)*Real(dx[1])*dV_dy[0];
                                
                                V_ghost[j*4 + 1] = -Real(2)*u_y_TT - Real(3)*u_y_T +
                                    Real(6)*V_ghost[(j + 1)*4 + 1] + Real(6)*Real(dx[1])*dV_dy[1];
                                
                                V_ghost[j*4 + 2] = -Real(2)*v_y_TT - Real(3)*v_y_T +
                                    Real(6)*V_ghost[(j + 1)*4 + 2] + Real(6)*Real(dx[1])*dV_dy[2];
                                
                                V_ghost[j*4 + 3] = -Real(2)*p_y_TT - Real(3)*p_y_T +
                                    Real(6)*V_ghost[(j + 1)*4 + 3] + Real(6)*Real(dx[1])*dV_dy[3];
                            }
                            else if (j == num_ghosts_to_fill - 3)
                            {
                                V_ghost[j*4 + 0] = Real(3)*rho_y_TT + Real(10)*rho_y_T -
                                    Real(18)*V_ghost[(j + 2)*4 + 0] + Real(6)*V_ghost[(j + 1)*4 + 0] -
                                    Real(12)*Real(dx[1])*dV_dy[0];
                                
                                V_ghost[j*4 + 1] = Real(3)*u_y_TT + Real(10)*u_y_T -
                                    Real(18)*V_ghost[(j + 2)*4 + 1] + Real(6)*V_ghost[(j + 1)*4 + 1] -
                                    Real(12)*Real(dx[1])*dV_dy[1];
                                
                                V_ghost[j*4 + 2] = Real(3)*v_y_TT + Real(10)*v_y_T -
                                    Real(18)*V_ghost[(j + 2)*4 + 2] + Real(6)*V_ghost[(j + 1)*4 + 2] -
                                    Real(12)*Real(dx[1])*dV_dy[2];
                                
                                V_ghost[j*4 + 3] = Real(3)*p_y_TT + Real(10)*p_y_T -
                                    Real(18)*V_ghost[(j + 2)*4 + 3] + Real(6)*V_ghost[(j + 1)*4 + 3] -
                                    Real(12)*Real(dx[1])*dV_dy[3];
                            }
                            else if (j == num_ghosts_to_fill - 4)
                            {
                                V_ghost[j*4 + 0] = -Real(4)*rho_y_TT - Real(65)/Real(3)*rho_y_T +
                                    Real(40)*V_ghost[(j + 3)*4 + 0] - Real(20)*V_ghost[(j + 2)*4 + 0] +
                                    Real(20)/Real(3)*V_ghost[(j + 1)*4 + 0] + Real(20)*Real(dx[1])*dV_dy[0];
                                
                                V_ghost[j*4 + 1] = -Real(4)*u_y_TT - Real(65)/Real(3)*u_y_T +
                                    Real(40)*V_ghost[(j + 3)*4 + 1] - Real(20)*V_ghost[(j + 2)*4 + 1] +
                                    Real(20)/Real(3)*V_ghost[(j + 1)*4 + 1] + Real(20)*Real(dx[1])*dV_dy[1];
                                
                                V_ghost[j*4 + 2] = -Real(4)*v_y_TT - Real(65)/Real(3)*v_y_T +
                                    Real(40)*V_ghost[(j + 3)*4 + 2] - Real(20)*V_ghost[(j + 2)*4 + 2] +
                                    Real(20)/Real(3)*V_ghost[(j + 1)*4 + 2] + Real(20)*Real(dx[1])*dV_dy[2];
                                
                                V_ghost[j*4 + 3] = -Real(4)*p_y_TT - Real(65)/Real(3)*p_y_T +
                                    Real(40)*V_ghost[(j + 3)*4 + 3] - Real(20)*V_ghost[(j + 2)*4 + 3] +
                                    Real(20)/Real(3)*V_ghost[(j + 1)*4 + 3] + Real(20)*Real(dx[1])*dV_dy[3];
                            }
                            else if (j == num_ghosts_to_fill - 5)
                            {
                                V_ghost[j*4 + 0] = Real(5)*rho_y_TT + Real(77)/Real(2)*rho_y_T -
                                    Real(75)*V_ghost[(j + 4)*4 + 0] + Real(50)*V_ghost[(j + 3)*4 + 0] -
                                    Real(25)*V_ghost[(j + 2)*4 + 0] + Real(15)/Real(2)*V_ghost[(j + 1)*4 + 0] -
                                    Real(30)*Real(dx[1])*dV_dy[0];
                                
                                V_ghost[j*4 + 1] = Real(5)*u_y_TT + Real(77)/Real(2)*u_y_T -
                                    Real(75)*V_ghost[(j + 4)*4 + 1] + Real(50)*V_ghost[(j + 3)*4 + 1] -
                                    Real(25)*V_ghost[(j + 2)*4 + 1] + Real(15)/Real(2)*V_ghost[(j + 1)*4 + 1] -
                                    Real(30)*Real(dx[1])*dV_dy[1];
                                
                                V_ghost[j*4 + 2] = Real(5)*v_y_TT + Real(77)/Real(2)*v_y_T -
                                    Real(75)*V_ghost[(j + 4)*4 + 2] + Real(50)*V_ghost[(j + 3)*4 + 2] -
                                    Real(25)*V_ghost[(j + 2)*4 + 2] + Real(15)/Real(2)*V_ghost[(j + 1)*4 + 2] -
                                    Real(30)*Real(dx[1])*dV_dy[2];
                                
                                V_ghost[j*4 + 3] = Real(5)*p_y_TT + Real(77)/Real(2)*p_y_T -
                                    Real(75)*V_ghost[(j + 4)*4 + 3] + Real(50)*V_ghost[(j + 3)*4 + 3] -
                                    Real(25)*V_ghost[(j + 2)*4 + 3] + Real(15)/Real(2)*V_ghost[(j + 1)*4 + 3] -
                                    Real(30)*Real(dx[1])*dV_dy[3];
                            }
                            else if (j == num_ghosts_to_fill - 6)
                            {
                                V_ghost[j*4 + 0] = -Real(6)*rho_y_TT - Real(609)/Real(10)*rho_y_T +
                                    Real(126)*V_ghost[(j + 5)*4 + 0] - Real(105)*V_ghost[(j + 4)*4 + 0] +
                                    Real(70)*V_ghost[(j + 3)*4 + 0] - Real(63)/Real(2)*V_ghost[(j + 2)*4 + 0] +
                                    Real(42)/Real(5)*V_ghost[(j + 1)*4 + 0] + Real(42)*Real(dx[1])*dV_dy[0];
                                
                                V_ghost[j*4 + 1] = -Real(6)*u_y_TT - Real(609)/Real(10)*u_y_T +
                                    Real(126)*V_ghost[(j + 5)*4 + 1] - Real(105)*V_ghost[(j + 4)*4 + 1] +
                                    Real(70)*V_ghost[(j + 3)*4 + 1] - Real(63)/Real(2)*V_ghost[(j + 2)*4 + 1] +
                                    Real(42)/Real(5)*V_ghost[(j + 1)*4 + 1] + Real(42)*Real(dx[1])*dV_dy[1];
                                
                                V_ghost[j*4 + 2] = -Real(6)*v_y_TT - Real(609)/Real(10)*v_y_T +
                                    Real(126)*V_ghost[(j + 5)*4 + 2] - Real(105)*V_ghost[(j + 4)*4 + 2] +
                                    Real(70)*V_ghost[(j + 3)*4 + 2] - Real(63)/Real(2)*V_ghost[(j + 2)*4 + 2] +
                                    Real(42)/Real(5)*V_ghost[(j + 1)*4 + 2] + Real(42)*Real(dx[1])*dV_dy[2];
                                
                                V_ghost[j*4 + 3] = -Real(6)*p_y_TT - Real(609)/Real(10)*p_y_T +
                                    Real(126)*V_ghost[(j + 5)*4 + 3] - Real(105)*V_ghost[(j + 4)*4 + 3] +
                                    Real(70)*V_ghost[(j + 3)*4 + 3] - Real(63)/Real(2)*V_ghost[(j + 2)*4 + 3] +
                                    Real(42)/Real(5)*V_ghost[(j + 1)*4 + 3] + Real(42)*Real(dx[1])*dV_dy[3];
                            }
                            
                            Q[0][idx_cell_rho] = V_ghost[j*4 + 0];
                            Q[1][idx_cell_mom] = V_ghost[j*4 + 0]*V_ghost[j*4 + 1];
                            Q[2][idx_cell_mom] = V_ghost[j*4 + 0]*V_ghost[j*4 + 2];
                            
                            const Real epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getInternalEnergy(
                                    &V_ghost[j*4 + 0],
                                    &V_ghost[j*4 + 3],
                                    thermo_properties_ptr);
                            
                            const Real E = V_ghost[j*4 + 0]*epsilon +
                                half*(Q[1][idx_cell_mom]*Q[1][idx_cell_mom] + Q[2][idx_cell_mom]*Q[2][idx_cell_mom])/
                                    V_ghost[j*4 + 0];
                            
                            Q[3][idx_cell_E] = E;
                        }
                    }
                }
                else if (edge_loc == BDRY_LOC::YHI)
                {
                    const int num_ghosts_to_fill = fill_box_hi_idx[1] - fill_box_lo_idx[1] + 1;
                    TBOX_ASSERT(fill_box_lo_idx[1] == interior_box_hi_idx[1] + 1);
                    if (num_ghosts_to_fill > 6)
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelBoundaryUtilitiesSingleSpecies::fill2dEdgeBoundaryData()\n"
                            << "Non-reflecting outflow BC doesn't support more than six ghost cells yet!");
                    }
                    
                    for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                    {
                        // Get the grid spacing.
                        const double* const dx = patch_geom->getDx();
                        
                        // Compute one-sided derivatives in y-direction.
                        const int idx_cell_rho_y_B = (i + num_subghosts_conservative_var[0][0]) +
                            (interior_box_hi_idx[1] + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0];
                        
                        const int idx_cell_rho_y_BB = (i + num_subghosts_conservative_var[0][0]) +
                            (interior_box_hi_idx[1] - 1 + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0];
                        
                        const int idx_cell_rho_y_BBB = (i + num_subghosts_conservative_var[0][0]) +
                            (interior_box_hi_idx[1] - 2 + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0];
                        
                        const int idx_cell_mom_y_B = (i + num_subghosts_conservative_var[1][0]) +
                            (interior_box_hi_idx[1] + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0];
                        
                        const int idx_cell_mom_y_BB = (i + num_subghosts_conservative_var[1][0]) +
                            (interior_box_hi_idx[1] - 1 + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0];
                        
                        const int idx_cell_mom_y_BBB = (i + num_subghosts_conservative_var[1][0]) +
                            (interior_box_hi_idx[1] - 2 + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0];
                        
                        const int idx_cell_E_y_B = (i + num_subghosts_conservative_var[2][0]) +
                            (interior_box_hi_idx[1] + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0];
                        
                        const int idx_cell_E_y_BB = (i + num_subghosts_conservative_var[2][0]) +
                            (interior_box_hi_idx[1] - 1 + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0];
                        
                        const int idx_cell_E_y_BBB = (i + num_subghosts_conservative_var[2][0]) + 
                            (interior_box_hi_idx[1] - 2 + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0];
                        
                        const Real& rho_y_B   = Q[0][idx_cell_rho_y_B];
                        const Real& rho_y_BB  = Q[0][idx_cell_rho_y_BB];
                        const Real& rho_y_BBB = Q[0][idx_cell_rho_y_BBB];
                        
                        const Real u_y_B   = Q[1][idx_cell_mom_y_B]/rho_y_B;
                        const Real u_y_BB  = Q[1][idx_cell_mom_y_BB]/rho_y_BB;
                        const Real u_y_BBB = Q[1][idx_cell_mom_y_BBB]/rho_y_BBB;
                        
                        const Real v_y_B   = Q[2][idx_cell_mom_y_B]/rho_y_B;
                        const Real v_y_BB  = Q[2][idx_cell_mom_y_BB]/rho_y_BB;
                        const Real v_y_BBB = Q[2][idx_cell_mom_y_BBB]/rho_y_BBB;
                       
                        const Real half = Real(1)/Real(2);
                        const Real epsilon_y_B   = Q[3][idx_cell_E_y_B]/rho_y_B - half*(u_y_B*u_y_B + v_y_B*v_y_B);
                        const Real epsilon_y_BB  = Q[3][idx_cell_E_y_BB]/rho_y_BB - half*(u_y_BB*u_y_BB + v_y_BB*v_y_BB);
                        const Real epsilon_y_BBB = Q[3][idx_cell_E_y_BBB]/rho_y_BBB - half*(u_y_BBB*u_y_BBB + v_y_BBB*v_y_BBB);
                        
                        const Real p_y_B = d_equation_of_state_mixing_rules->getEquationOfState()->
                            getPressure(
                                &rho_y_B,
                                &epsilon_y_B,
                                thermo_properties_ptr);
                        
                        const Real p_y_BB = d_equation_of_state_mixing_rules->getEquationOfState()->
                            getPressure(
                                &rho_y_BB,
                                &epsilon_y_BB,
                                thermo_properties_ptr);
                        
                        const Real p_y_BBB = d_equation_of_state_mixing_rules->getEquationOfState()->
                            getPressure(
                                &rho_y_BBB,
                                &epsilon_y_BBB,
                                thermo_properties_ptr);
                        
                        /*
                         * Compute derivatives in y-direction.
                         */
                        
                        const Real drho_dy = (rho_y_BBB - Real(4)*rho_y_BB + Real(3)*rho_y_B)/(Real(2)*Real(dx[1]));
                        const Real du_dy   = (u_y_BBB - Real(4)*u_y_BB + Real(3)*u_y_B)/(Real(2)*Real(dx[1]));
                        const Real dv_dy   = (v_y_BBB - Real(4)*v_y_BB + Real(3)*v_y_B)/(Real(2)*Real(dx[1]));
                        const Real dp_dy   = (p_y_BBB - Real(4)*p_y_BB + Real(3)*p_y_B)/(Real(2)*Real(dx[1]));
                        
                        /*
                         * Compute derivatives in x-direction.
                         */
                        
                        Real du_dx = Real(0);
                        Real dv_dx = Real(0);
                        Real dp_dx = Real(0);
                        
                        if ((i + num_subghosts_conservative_var[0][0] == 0) ||
                            (i + num_subghosts_conservative_var[1][0] == 0) ||
                            (i + num_subghosts_conservative_var[2][0] == 0))
                        // if (((patch_geom->getTouchesRegularBoundary(0, 0)) && (i == interior_box_lo_idx[0])) ||
                        //     ((i + num_subghosts_conservative_var[0][0] == 0) ||
                        //      (i + num_subghosts_conservative_var[1][0] == 0) ||
                        //      (i + num_subghosts_conservative_var[2][0] == 0)))
                        {
                            // Patch is touching periodic boundary.
                            
                            const int idx_cell_rho_x_R = (i + 1 + num_subghosts_conservative_var[0][0]) +
                                (interior_box_hi_idx[1] + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_mom_x_R = (i + 1 + num_subghosts_conservative_var[1][0]) +
                                (interior_box_hi_idx[1] + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0];
                            
                            const int idx_cell_E_x_R = (i + 1 + num_subghosts_conservative_var[2][0]) +
                                (interior_box_hi_idx[1] + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0];
                            
                            const Real& rho_x_R = Q[0][idx_cell_rho_x_R];
                            const Real u_x_R = Q[1][idx_cell_mom_x_R]/rho_x_R;
                            const Real v_x_R = Q[2][idx_cell_mom_x_R]/rho_x_R;
                            const Real epsilon_x_R = Q[3][idx_cell_E_x_R]/rho_x_R - half*(u_x_R*u_x_R + v_x_R*v_x_R);
                            
                            const Real p_x_R = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getPressure(
                                    &rho_x_R,
                                    &epsilon_x_R,
                                    thermo_properties_ptr);
                            
                            // One-sided derivatives.
                            du_dx = (u_x_R - u_y_B)/Real(dx[0]);
                            dv_dx = (v_x_R - v_y_B)/Real(dx[0]);
                            dp_dx = (p_x_R - p_y_B)/Real(dx[0]);
                        }
                        else if ((i + num_subghosts_conservative_var[0][0] + 1 == subghostcell_dims_conservative_var[0][0]) ||
                                 (i + num_subghosts_conservative_var[1][0] + 1 == subghostcell_dims_conservative_var[1][0]) ||
                                 (i + num_subghosts_conservative_var[2][0] + 1 == subghostcell_dims_conservative_var[2][0]))
                        // else if (((patch_geom->getTouchesRegularBoundary(0, 1)) && (i == interior_box_hi_idx[0])) ||
                        //          ((i + num_subghosts_conservative_var[0][0] + 1 == subghostcell_dims_conservative_var[0][0]) ||
                        //           (i + num_subghosts_conservative_var[1][0] + 1 == subghostcell_dims_conservative_var[1][0]) ||
                        //           (i + num_subghosts_conservative_var[2][0] + 1 == subghostcell_dims_conservative_var[2][0])))
                        {
                            // Patch is touching periodic boundary.
                            
                            const int idx_cell_rho_x_L = (i - 1 + num_subghosts_conservative_var[0][0]) +
                                (interior_box_hi_idx[1] + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_mom_x_L = (i - 1 + num_subghosts_conservative_var[1][0]) +
                                (interior_box_hi_idx[1] + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0];
                            
                            const int idx_cell_E_x_L = (i - 1 + num_subghosts_conservative_var[2][0]) +
                                (interior_box_hi_idx[1] + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0];
                            
                            const Real& rho_x_L = Q[0][idx_cell_rho_x_L];
                            const Real u_x_L = Q[1][idx_cell_mom_x_L]/rho_x_L;
                            const Real v_x_L = Q[2][idx_cell_mom_x_L]/rho_x_L;
                            const Real epsilon_x_L = Q[3][idx_cell_E_x_L]/rho_x_L - half*(u_x_L*u_x_L + v_x_L*v_x_L);
                            
                            const Real p_x_L = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getPressure(
                                    &rho_x_L,
                                    &epsilon_x_L,
                                    thermo_properties_ptr);
                            
                            // One-sided derivatives.
                            du_dx = (u_y_B - u_x_L)/Real(dx[0]);
                            dv_dx = (v_y_B - v_x_L)/Real(dx[0]);
                            dp_dx = (p_y_B - p_x_L)/Real(dx[0]);
                        }
                        else
                        {
                            const int idx_cell_rho_x_L = (i - 1 + num_subghosts_conservative_var[0][0]) +
                                (interior_box_hi_idx[1] + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_rho_x_R = (i + 1 + num_subghosts_conservative_var[0][0]) +
                                (interior_box_hi_idx[1] + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_mom_x_L = (i - 1 + num_subghosts_conservative_var[1][0]) +
                                (interior_box_hi_idx[1] + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0];
                            
                            const int idx_cell_mom_x_R = (i + 1 + num_subghosts_conservative_var[1][0]) +
                                (interior_box_hi_idx[1] + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0];
                            
                            const int idx_cell_E_x_L = (i - 1 + num_subghosts_conservative_var[2][0]) +
                                (interior_box_hi_idx[1] + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0];
                            
                            const int idx_cell_E_x_R = (i + 1 + num_subghosts_conservative_var[2][0]) +
                                (interior_box_hi_idx[1] + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0];
                            
                            const Real& rho_x_L = Q[0][idx_cell_rho_x_L];
                            const Real& rho_x_R = Q[0][idx_cell_rho_x_R];
                            
                            const Real u_x_L = Q[1][idx_cell_mom_x_L]/rho_x_L;
                            const Real u_x_R = Q[1][idx_cell_mom_x_R]/rho_x_R;
                            
                            const Real v_x_L = Q[2][idx_cell_mom_x_L]/rho_x_L;
                            const Real v_x_R = Q[2][idx_cell_mom_x_R]/rho_x_R;
                            
                            const Real epsilon_x_L = Q[3][idx_cell_E_x_L]/rho_x_L - half*(u_x_L*u_x_L + v_x_L*v_x_L);
                            const Real epsilon_x_R = Q[3][idx_cell_E_x_R]/rho_x_R - half*(u_x_R*u_x_R + v_x_R*v_x_R);
                            
                            const Real p_x_L = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getPressure(
                                    &rho_x_L,
                                    &epsilon_x_L,
                                    thermo_properties_ptr);
                            
                            const Real p_x_R = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getPressure(
                                    &rho_x_R,
                                    &epsilon_x_R,
                                    thermo_properties_ptr);
                            
                            // Central derivatives.
                            du_dx = (u_x_R - u_x_L)/(Real(2)*Real(dx[0]));
                            dv_dx = (v_x_R - v_x_L)/(Real(2)*Real(dx[0]));
                            dp_dx = (p_x_R - p_x_L)/(Real(2)*Real(dx[0]));
                        }
                        
                        // Compute wave speed (v - c) at the boundary.
                        
                        const Real c_y_B = d_equation_of_state_mixing_rules->getEquationOfState()->
                            getSoundSpeed(
                                &rho_y_B,
                                &p_y_B,
                                thermo_properties_ptr);
                        
                        const Real lambda_1 = v_y_B - c_y_B;
                        
                        // Compute vector Lambda^(-1) * L.
                        
                        Real Lambda_inv_L[4];
                        
                        const Real& p_t         = d_bdry_edge_nonreflecting_outflow_p_t[edge_loc];
                        const Real& sigma       = d_bdry_edge_nonreflecting_outflow_sigma[edge_loc];
                        const Real& beta        = d_bdry_edge_nonreflecting_outflow_beta[edge_loc];
                        const Real& length_char = d_bdry_edge_nonreflecting_outflow_length_char[edge_loc];
                        
                        const Real T_1 = u_y_B*(dp_dx - rho_y_B*c_y_B*dv_dx) + rho_y_B*c_y_B*c_y_B*du_dx;
                        
                        const Real M_sq = (v_y_B*v_y_B + u_y_B*u_y_B)/(c_y_B*c_y_B);
                        const Real K = sigma*c_y_B*(Real(1) - M_sq)/length_char;
                        
                        Lambda_inv_L[0] = (Real(1)/lambda_1)*(K*(p_y_B - p_t) - (Real(1) - beta)*T_1);
                        Lambda_inv_L[1] = du_dy;
                        Lambda_inv_L[2] = c_y_B*c_y_B*drho_dy - dp_dy;
                        Lambda_inv_L[3] = dp_dy + rho_y_B*c_y_B*dv_dy;
                        
                        // Compute dV_dy.
                        
                        const Real c_sq_inv  = Real(1)/(c_y_B*c_y_B);
                        const Real rho_c_inv = Real(1)/(rho_y_B*c_y_B);
                        
                        Real dV_dy[4];
                        
                        dV_dy[0] = half*c_sq_inv*(Lambda_inv_L[0] + Lambda_inv_L[3]) + c_sq_inv*Lambda_inv_L[2];
                        dV_dy[1] = Lambda_inv_L[1];
                        dV_dy[2] = half*rho_c_inv*(-Lambda_inv_L[0] + Lambda_inv_L[3]);
                        dV_dy[3] = half*(Lambda_inv_L[0] + Lambda_inv_L[3]);
                        
                        Real V_ghost[4*num_ghosts_to_fill];
                        
                        for (int j = 0; j < num_ghosts_to_fill; j++)
                        {
                            const int idx_cell_rho = (i + num_subghosts_conservative_var[0][0]) +
                                (j + fill_box_lo_idx[1] + num_subghosts_conservative_var[0][1])*
                                    subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_mom = (i + num_subghosts_conservative_var[1][0]) +
                                (j + fill_box_lo_idx[1] + num_subghosts_conservative_var[1][1])*
                                    subghostcell_dims_conservative_var[1][0];
                            
                            const int idx_cell_E = (i + num_subghosts_conservative_var[2][0]) +
                                (j + fill_box_lo_idx[1] + num_subghosts_conservative_var[2][1])*
                                    subghostcell_dims_conservative_var[2][0];
                            
                            if (j == 0)
                            {
                                V_ghost[j*4 + 0] = rho_y_BB + Real(2)*Real(dx[1])*dV_dy[0];
                                V_ghost[j*4 + 1] = u_y_BB   + Real(2)*Real(dx[1])*dV_dy[1];
                                V_ghost[j*4 + 2] = v_y_BB   + Real(2)*Real(dx[1])*dV_dy[2];
                                V_ghost[j*4 + 3] = p_y_BB   + Real(2)*Real(dx[1])*dV_dy[3];
                            }
                            else if (j == 1)
                            {
                                V_ghost[j*4 + 0] = -Real(2)*rho_y_BB - Real(3)*rho_y_B +
                                    Real(6)*V_ghost[(j - 1)*4 + 0] - Real(6)*Real(dx[1])*dV_dy[0];
                                
                                V_ghost[j*4 + 1] = -Real(2)*u_y_BB - Real(3)*u_y_B +
                                    Real(6)*V_ghost[(j - 1)*4 + 1] - Real(6)*Real(dx[1])*dV_dy[1];
                                
                                V_ghost[j*4 + 2] = -Real(2)*v_y_BB - Real(3)*v_y_B +
                                    Real(6)*V_ghost[(j - 1)*4 + 2] - Real(6)*Real(dx[1])*dV_dy[2];
                                
                                V_ghost[j*4 + 3] = -Real(2)*p_y_BB - Real(3)*p_y_B +
                                    Real(6)*V_ghost[(j - 1)*4 + 3] - Real(6)*Real(dx[1])*dV_dy[3];
                            }
                            else if (j == 2)
                            {
                                V_ghost[j*4 + 0] = Real(3)*rho_y_BB + Real(10)*rho_y_B -
                                    Real(18)*V_ghost[(j - 2)*4 + 0] + Real(6)*V_ghost[(j - 1)*4 + 0] +
                                    Real(12)*Real(dx[1])*dV_dy[0];
                                
                                V_ghost[j*4 + 1] = Real(3)*u_y_BB + Real(10)*u_y_B -
                                    Real(18)*V_ghost[(j - 2)*4 + 1] + Real(6)*V_ghost[(j - 1)*4 + 1] +
                                    Real(12)*Real(dx[1])*dV_dy[1];
                                
                                V_ghost[j*4 + 2] = Real(3)*v_y_BB + Real(10)*v_y_B -
                                    Real(18)*V_ghost[(j - 2)*4 + 2] + Real(6)*V_ghost[(j - 1)*4 + 2] +
                                    Real(12)*Real(dx[1])*dV_dy[2];
                                
                                V_ghost[j*4 + 3] = Real(3)*p_y_BB + Real(10)*p_y_B -
                                    Real(18)*V_ghost[(j - 2)*4 + 3] + Real(6)*V_ghost[(j - 1)*4 + 3] +
                                    Real(12)*Real(dx[1])*dV_dy[3];
                            }
                            else if (j == 3)
                            {
                                V_ghost[j*4 + 0] = -Real(4)*rho_y_BB - Real(65)/Real(3)*rho_y_B +
                                    Real(40)*V_ghost[(j - 3)*4 + 0] - Real(20)*V_ghost[(j - 2)*4 + 0] +
                                    Real(20)/Real(3)*V_ghost[(j - 1)*4 + 0] - Real(20)*Real(dx[1])*dV_dy[0];
                                
                                V_ghost[j*4 + 1] = -Real(4)*u_y_BB - Real(65)/Real(3)*u_y_B +
                                    Real(40)*V_ghost[(j - 3)*4 + 1] - Real(20)*V_ghost[(j - 2)*4 + 1] +
                                    Real(20)/Real(3)*V_ghost[(j - 1)*4 + 1] - Real(20)*Real(dx[1])*dV_dy[1];
                                
                                V_ghost[j*4 + 2] = -Real(4)*v_y_BB - Real(65)/Real(3)*v_y_B +
                                    Real(40)*V_ghost[(j - 3)*4 + 2] - Real(20)*V_ghost[(j - 2)*4 + 2] +
                                    Real(20)/Real(3)*V_ghost[(j - 1)*4 + 2] - Real(20)*Real(dx[1])*dV_dy[2];
                                
                                V_ghost[j*4 + 3] = -Real(4)*p_y_BB - Real(65)/Real(3)*p_y_B +
                                    Real(40)*V_ghost[(j - 3)*4 + 3] - Real(20)*V_ghost[(j - 2)*4 + 3] +
                                    Real(20)/Real(3)*V_ghost[(j - 1)*4 + 3] - Real(20)*Real(dx[1])*dV_dy[3];
                            }
                            else if (j == 4)
                            {
                                V_ghost[j*4 + 0] = Real(5)*rho_y_BB + Real(77)/Real(2)*rho_y_B -
                                    Real(75)*V_ghost[(j - 4)*4 + 0] + Real(50)*V_ghost[(j - 3)*4 + 0] -
                                    Real(25)*V_ghost[(j - 2)*4 + 0] + Real(15)/Real(2)*V_ghost[(j - 1)*4 + 0] +
                                    Real(30)*Real(dx[1])*dV_dy[0];
                                
                                V_ghost[j*4 + 1] = Real(5)*u_y_BB + Real(77)/Real(2)*u_y_B -
                                    Real(75)*V_ghost[(j - 4)*4 + 1] + Real(50)*V_ghost[(j - 3)*4 + 1] -
                                    Real(25)*V_ghost[(j - 2)*4 + 1] + Real(15)/Real(2)*V_ghost[(j - 1)*4 + 1] +
                                    Real(30)*Real(dx[1])*dV_dy[1];
                                
                                V_ghost[j*4 + 2] = Real(5)*v_y_BB + Real(77)/Real(2)*v_y_B -
                                    Real(75)*V_ghost[(j - 4)*4 + 2] + Real(50)*V_ghost[(j - 3)*4 + 2] -
                                    Real(25)*V_ghost[(j - 2)*4 + 2] + Real(15)/Real(2)*V_ghost[(j - 1)*4 + 2] +
                                    Real(30)*Real(dx[1])*dV_dy[2];
                                
                                V_ghost[j*4 + 3] = Real(5)*p_y_BB + Real(77)/Real(2)*p_y_B -
                                    Real(75)*V_ghost[(j - 4)*4 + 3] + Real(50)*V_ghost[(j - 3)*4 + 3] -
                                    Real(25)*V_ghost[(j - 2)*4 + 3] + Real(15)/Real(2)*V_ghost[(j - 1)*4 + 3] +
                                    Real(30)*Real(dx[1])*dV_dy[3];
                            }
                            else if (j == 5)
                            {
                                V_ghost[j*4 + 0] = -Real(6)*rho_y_BB - Real(609)/Real(10)*rho_y_B +
                                    Real(126)*V_ghost[(j - 5)*4 + 0] - Real(105)*V_ghost[(j - 4)*4 + 0] +
                                    Real(70)*V_ghost[(j - 3)*4 + 0] - Real(63)/Real(2)*V_ghost[(j - 2)*4 + 0] +
                                    Real(42)/Real(5)*V_ghost[(j - 1)*4 + 0] - Real(42)*Real(dx[1])*dV_dy[0];
                                
                                V_ghost[j*4 + 1] = -Real(6)*u_y_BB - Real(609)/Real(10)*u_y_B +
                                    Real(126)*V_ghost[(j - 5)*4 + 1] - Real(105)*V_ghost[(j - 4)*4 + 1] +
                                    Real(70)*V_ghost[(j - 3)*4 + 1] - Real(63)/Real(2)*V_ghost[(j - 2)*4 + 1] +
                                    Real(42)/Real(5)*V_ghost[(j - 1)*4 + 1] - Real(42)*Real(dx[1])*dV_dy[1];
                                
                                V_ghost[j*4 + 2] = -Real(6)*v_y_BB - Real(609)/Real(10)*v_y_B +
                                    Real(126)*V_ghost[(j - 5)*4 + 2] - Real(105)*V_ghost[(j - 4)*4 + 2] +
                                    Real(70)*V_ghost[(j - 3)*4 + 2] - Real(63)/Real(2)*V_ghost[(j - 2)*4 + 2] +
                                    Real(42)/Real(5)*V_ghost[(j - 1)*4 + 2] - Real(42)*Real(dx[1])*dV_dy[2];
                                
                                V_ghost[j*4 + 3] = -Real(6)*p_y_BB - Real(609)/Real(10)*p_y_B +
                                    Real(126)*V_ghost[(j - 5)*4 + 3] - Real(105)*V_ghost[(j - 4)*4 + 3] +
                                    Real(70)*V_ghost[(j - 3)*4 + 3] - Real(63)/Real(2)*V_ghost[(j - 2)*4 + 3] +
                                    Real(42)/Real(5)*V_ghost[(j - 1)*4 + 3] - Real(42)*Real(dx[1])*dV_dy[3];
                            }
                            
                            Q[0][idx_cell_rho] = V_ghost[j*4 + 0];
                            Q[1][idx_cell_mom] = V_ghost[j*4 + 0]*V_ghost[j*4 + 1];
                            Q[2][idx_cell_mom] = V_ghost[j*4 + 0]*V_ghost[j*4 + 2];
                            
                            const Real epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getInternalEnergy(
                                    &V_ghost[j*4 + 0],
                                    &V_ghost[j*4 + 3],
                                    thermo_properties_ptr);
                            
                            const Real E = V_ghost[j*4 + 0]*epsilon +
                                half*(Q[1][idx_cell_mom]*Q[1][idx_cell_mom] + Q[2][idx_cell_mom]*Q[2][idx_cell_mom])/
                                    V_ghost[j*4 + 0];
                            
                            Q[3][idx_cell_E] = E;
                        }
                    }
                }
                
                // Remove edge locations that have boundary conditions identified.
                bdry_edge_locs.erase(std::remove(bdry_edge_locs.begin(), bdry_edge_locs.end(), edge_loc),
                    bdry_edge_locs.end());
            }
        }
    }
    
    for (int ei = 0; ei < static_cast<int>(edge_bdry.size()); ei++)
    {
        TBOX_ASSERT(edge_bdry[ei].getBoundaryType() == BDRY::EDGE2D);
        
        int edge_loc = edge_bdry[ei].getLocationIndex();
        
        if (std::find(bdry_edge_locs.begin(), bdry_edge_locs.end(), edge_loc) !=
            bdry_edge_locs.end())
        {
            TBOX_ERROR("FlowModelBoundaryUtilitiesSingleSpecies::fill2dEdgeBoundaryData()\n"
                << "Invalid edge boundary condition!\n"
                << "edge_loc = '" << edge_loc << "'." << std::endl
                << "bdry_edge_conds[edge_loc] = '" << bdry_edge_conds[edge_loc] << "'."
                << std::endl);
        }
    }
}


/*
 * Function to fill 2d node boundary values for a patch.
 * Node locations that have boundary conditions identified are removed from the container.
 */
void
FlowModelBoundaryUtilitiesSingleSpecies::fill2dNodeBoundaryData(
    const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& conservative_var_data,
    const hier::Patch& patch,
    std::vector<int>& bdry_node_locs,
    const std::vector<int>& bdry_node_conds,
    const std::vector<std::vector<Real> >& bdry_edge_values,
    const hier::IntVector& ghost_width_to_fill)
{
    TBOX_ASSERT(static_cast<int>(conservative_var_data.size()) == 3);
    for (int vi = 0; vi < static_cast<int>(conservative_var_data.size()); vi++)
    {
        TBOX_ASSERT(conservative_var_data[vi]);
    }
    TBOX_ASSERT(static_cast<int>(bdry_node_locs.size()) <= NUM_2D_NODES);
    if (static_cast<int>(bdry_node_locs.size()) > 0)
    {
        TBOX_ASSERT(*min_element(bdry_node_locs.begin(), bdry_node_locs.end()) >= 0);
        TBOX_ASSERT(*max_element(bdry_node_locs.begin(), bdry_node_locs.end()) < NUM_2D_NODES);
    }
    TBOX_ASSERT(static_cast<int>(bdry_node_conds.size()) == NUM_2D_NODES);
    TBOX_ASSERT(static_cast<int>(bdry_edge_values.size()) == 3);
    for (int vi = 0; vi < static_cast<int>(bdry_edge_values.size()); vi++)
    {
        TBOX_ASSERT(static_cast<int>(bdry_edge_values[vi].size()) ==
                    NUM_2D_EDGES*(conservative_var_data[vi]->getDepth()));
    }
    
    TBOX_DIM_ASSERT(ghost_width_to_fill.getDim() == tbox::Dimension(2));
    
    for (int vi = 0; vi < static_cast<int>(conservative_var_data.size()); vi++)
    {
        TBOX_ASSERT_OBJDIM_EQUALITY3(*conservative_var_data[vi], patch, ghost_width_to_fill);
    }
    
    const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
        HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    TBOX_ASSERT(patch_geom);
    
    hier::IntVector num_ghosts(conservative_var_data[0]->getGhostCellWidth());
    for (int vi = 1; vi < static_cast<int>(conservative_var_data.size()); vi++)
    {
        num_ghosts = hier::IntVector::min(
            num_ghosts,
            conservative_var_data[vi]->getGhostCellWidth());
    }
    
    /*
     * Determine the ghost cell width to fill.
     */
    
    hier::IntVector gcw_to_fill(tbox::Dimension(2));
    
    // If the ghost fill width is not used, it is set to the ghost cell width of the data.
    if (ghost_width_to_fill == -hier::IntVector::getOne(tbox::Dimension(2)))
    {
        gcw_to_fill = num_ghosts;
    }
    else
    {
        gcw_to_fill = hier::IntVector::min(
            num_ghosts,
            ghost_width_to_fill);
    }
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box& interior_box(patch.getBox());
    hier::Index interior_box_lo_idx(interior_box.lower());
    hier::Index interior_box_hi_idx(interior_box.upper());
    
    /*
     * Offset the indices.
     */
    
    interior_box_lo_idx = interior_box_lo_idx - interior_box.lower();
    interior_box_hi_idx = interior_box_hi_idx - interior_box.lower();
    
    const std::vector<hier::BoundaryBox>& node_bdry =
        patch_geom->getCodimensionBoundaries(BDRY::NODE2D);
    
    for (int ni = 0; ni < static_cast<int>(node_bdry.size()); ni++)
    {
        TBOX_ASSERT(node_bdry[ni].getBoundaryType() == BDRY::NODE2D);
        
        int node_loc = node_bdry[ni].getLocationIndex();
        
        if (std::find(bdry_node_locs.begin(), bdry_node_locs.end(), node_loc) !=
            bdry_node_locs.end())
        {
            hier::Box fill_box(patch_geom->getBoundaryFillBox(
                node_bdry[ni],
                interior_box,
                gcw_to_fill));
            
            hier::Index fill_box_lo_idx(fill_box.lower());
            hier::Index fill_box_hi_idx(fill_box.upper());
            
            /*
             * Offset the indices.
             */
            
            fill_box_lo_idx = fill_box_lo_idx - interior_box.lower();
            fill_box_hi_idx = fill_box_hi_idx - interior_box.lower();
            
            int edge_loc_0 = -1;
            int edge_loc_1 = -1;
            
            switch (node_loc)
            {
                case NODE_BDRY_LOC_2D::XLO_YLO:
                {
                    edge_loc_0 = BDRY_LOC::XLO;
                    edge_loc_1 = BDRY_LOC::YLO;
                    
                    break;
                }
                case NODE_BDRY_LOC_2D::XHI_YLO:
                {
                    edge_loc_0 = BDRY_LOC::XHI;
                    edge_loc_1 = BDRY_LOC::YLO;
                    
                    break;
                }
                case NODE_BDRY_LOC_2D::XLO_YHI:
                {
                    edge_loc_0 = BDRY_LOC::XLO;
                    edge_loc_1 = BDRY_LOC::YHI;
                    
                    break;
                }
                case NODE_BDRY_LOC_2D::XHI_YHI:
                {
                    edge_loc_0 = BDRY_LOC::XHI;
                    edge_loc_1 = BDRY_LOC::YHI;
                    
                    break;
                }
            }
            
            if ((bdry_node_conds[node_loc] == BDRY_COND::FLOW_MODEL::XADIABATIC_NO_SLIP) ||
                (bdry_node_conds[node_loc] == BDRY_COND::FLOW_MODEL::YADIABATIC_NO_SLIP) ||
                (bdry_node_conds[node_loc] == BDRY_COND::FLOW_MODEL::XISOTHERMAL_NO_SLIP) ||
                (bdry_node_conds[node_loc] == BDRY_COND::FLOW_MODEL::YISOTHERMAL_NO_SLIP))
            {
                /*
                 * Get the pointers to the conservative variables.
                 * Get the numbers of ghost cells and the dimensions of the ghost cell boxes of
                 * the conservative variables.
                 */
                
                std::vector<Real*> Q;
                Q.reserve(d_num_eqn);
                
                std::vector<hier::IntVector> num_subghosts_conservative_var;
                std::vector<hier::IntVector> subghostcell_dims_conservative_var;
                
                num_subghosts_conservative_var.reserve(3);
                subghostcell_dims_conservative_var.reserve(3);
                
                int count_eqn = 0;
                
                for (int vi = 0; vi < static_cast<int>(conservative_var_data.size()); vi++)
                {
                    int depth = conservative_var_data[vi]->getDepth();
                    
                    for (int di = 0; di < depth; di++)
                    {
                        // If the last element of the conservative variable vector is not in the
                        // system of equations, ignore it.
                        if (count_eqn >= d_num_eqn)
                            break;
                        
                        Q.push_back(conservative_var_data[vi]->getPointer(di));
                        
                        count_eqn++;
                    }
                    
                    num_subghosts_conservative_var.push_back(
                        conservative_var_data[vi]->getGhostCellWidth());
                    subghostcell_dims_conservative_var.push_back(
                        conservative_var_data[vi]->getGhostBox().numberCells());
                }
                
                // Get the thermodynamic properties of the species.
                std::vector<const Real*> thermo_properties_ptr;
                thermo_properties_ptr.reserve(static_cast<int> (d_thermo_properties.size()));
                for (int ti = 0; ti < static_cast<int> (d_thermo_properties.size()); ti++)
                {
                    thermo_properties_ptr.push_back(&d_thermo_properties[ti]);
                }
                
                if (bdry_node_conds[node_loc] == BDRY_COND::FLOW_MODEL::XADIABATIC_NO_SLIP)
                {
                    for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                    {
                        for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                        {
                            const int idx_cell_rho = (i + num_subghosts_conservative_var[0][0]) +
                                (j + num_subghosts_conservative_var[0][1])*
                                    subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_mom = (i + num_subghosts_conservative_var[1][0]) +
                                (j + num_subghosts_conservative_var[1][1])*
                                    subghostcell_dims_conservative_var[1][0];
                            
                            const int idx_cell_E = (i + num_subghosts_conservative_var[2][0]) +
                                (j + num_subghosts_conservative_var[2][1])*
                                    subghostcell_dims_conservative_var[2][0];
                            
                            int idx_cell_pivot_rho = idx_cell_rho;
                            int idx_cell_pivot_mom = idx_cell_mom;
                            int idx_cell_pivot_E = idx_cell_E;
                            
                            if (edge_loc_0 == BDRY_LOC::XLO)
                            {
                                idx_cell_pivot_rho = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                        num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0];
                                
                                idx_cell_pivot_mom = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                        num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0];
                                
                                idx_cell_pivot_E = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                        num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0];
                            }
                            else if (edge_loc_0 == BDRY_LOC::XHI)
                            {
                                idx_cell_pivot_rho = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                        num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0];
                                
                                idx_cell_pivot_mom = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                        num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0];
                                
                                idx_cell_pivot_E = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                        num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0];
                            }
                            
                            /*
                             * Set the values for density and momentum.
                             */
                            
                            Q[0][idx_cell_rho] = Q[0][idx_cell_pivot_rho];
                            Q[1][idx_cell_mom] = -Q[1][idx_cell_pivot_mom] +
                                Real(2)*Q[0][idx_cell_pivot_rho]*d_bdry_edge_adiabatic_no_slip_vel[edge_loc_0*2];
                            Q[2][idx_cell_mom] = -Q[2][idx_cell_pivot_mom] +
                                Real(2)*Q[0][idx_cell_pivot_rho]*d_bdry_edge_adiabatic_no_slip_vel[edge_loc_0*2 + 1];
                            
                            /*
                             * Set the values for total internal energy.
                             */
                            
                            Real epsilon_pivot = (Q[3][idx_cell_pivot_E] -
                                Real(1)/Real(2)*(Q[1][idx_cell_pivot_mom]*Q[1][idx_cell_pivot_mom] +
                                Q[2][idx_cell_pivot_mom]*Q[2][idx_cell_pivot_mom])/
                                Q[0][idx_cell_pivot_rho])/Q[0][idx_cell_pivot_rho];
                            
                            Real p_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getPressure(
                                    &Q[0][idx_cell_pivot_rho],
                                    &epsilon_pivot,
                                    thermo_properties_ptr);
                            
                            Real T_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getTemperature(
                                    &Q[0][idx_cell_pivot_rho],
                                    &p_pivot,
                                    thermo_properties_ptr);
                            
                            Real T = T_pivot;
                            
                            Real epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getInternalEnergyFromTemperature(
                                    &Q[0][idx_cell_rho],
                                    &T,
                                    thermo_properties_ptr);
                            
                            Real E = Q[0][idx_cell_rho]*epsilon +
                                Real(1)/Real(2)*(Q[1][idx_cell_mom]*Q[1][idx_cell_mom] + Q[2][idx_cell_mom]*Q[2][idx_cell_mom])/
                                Q[0][idx_cell_rho];
                            
                            Q[3][idx_cell_E] = E;
                        }
                    }
                    
                    // Remove node locations that have boundary conditions identified.
                    bdry_node_locs.erase(std::remove(bdry_node_locs.begin(), bdry_node_locs.end(), node_loc),
                        bdry_node_locs.end());
                }
                else if (bdry_node_conds[node_loc] == BDRY_COND::FLOW_MODEL::YADIABATIC_NO_SLIP)
                {
                    for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                    {
                        for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                        {
                            const int idx_cell_rho = (i + num_subghosts_conservative_var[0][0]) +
                                (j + num_subghosts_conservative_var[0][1])*
                                    subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_mom = (i + num_subghosts_conservative_var[1][0]) +
                                (j + num_subghosts_conservative_var[1][1])*
                                    subghostcell_dims_conservative_var[1][0];
                            
                            const int idx_cell_E = (i + num_subghosts_conservative_var[2][0]) +
                                (j + num_subghosts_conservative_var[2][1])*
                                    subghostcell_dims_conservative_var[2][0];
                            
                            int idx_cell_pivot_rho = idx_cell_rho;
                            int idx_cell_pivot_mom = idx_cell_mom;
                            int idx_cell_pivot_E = idx_cell_E;
                            
                            if (edge_loc_1 == BDRY_LOC::YLO)
                            {
                                idx_cell_pivot_rho = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                        num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0];
                                
                                idx_cell_pivot_mom = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                        num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0];
                                
                                idx_cell_pivot_E = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                        num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0];
                            }
                            else if (edge_loc_1 == BDRY_LOC::YHI)
                            {
                                idx_cell_pivot_rho = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                        num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0];
                                
                                idx_cell_pivot_mom = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                        num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0];
                                
                                idx_cell_pivot_E = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                        num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0];
                            }
                            
                            /*
                             * Set the values for density and momentum.
                             */
                            
                            Q[0][idx_cell_rho] = Q[0][idx_cell_pivot_rho];
                            Q[1][idx_cell_mom] = -Q[1][idx_cell_pivot_mom] +
                                Real(2)*Q[0][idx_cell_pivot_rho]*d_bdry_edge_adiabatic_no_slip_vel[edge_loc_1*2];
                            Q[2][idx_cell_mom] = -Q[2][idx_cell_pivot_mom] +
                                Real(2)*Q[0][idx_cell_pivot_rho]*d_bdry_edge_adiabatic_no_slip_vel[edge_loc_1*2 + 1];
                            
                            /*
                             * Set the values for total internal energy.
                             */
                            
                            Real epsilon_pivot = (Q[3][idx_cell_pivot_E] -
                                Real(1)/Real(2)*(Q[1][idx_cell_pivot_mom]*Q[1][idx_cell_pivot_mom] +
                                Q[2][idx_cell_pivot_mom]*Q[2][idx_cell_pivot_mom])/
                                Q[0][idx_cell_pivot_rho])/Q[0][idx_cell_pivot_rho];
                            
                            Real p_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getPressure(
                                    &Q[0][idx_cell_pivot_rho],
                                    &epsilon_pivot,
                                    thermo_properties_ptr);
                            
                            Real T_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getTemperature(
                                    &Q[0][idx_cell_pivot_rho],
                                    &p_pivot,
                                    thermo_properties_ptr);
                            
                            Real T = T_pivot;
                            
                            Real epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getInternalEnergyFromTemperature(
                                    &Q[0][idx_cell_rho],
                                    &T,
                                    thermo_properties_ptr);
                            
                            Real E = Q[0][idx_cell_rho]*epsilon +
                                Real(1)/Real(2)*(Q[1][idx_cell_mom]*Q[1][idx_cell_mom] + Q[2][idx_cell_mom]*Q[2][idx_cell_mom])/
                                Q[0][idx_cell_rho];
                            
                            Q[3][idx_cell_E] = E;
                        }
                    }
                    
                    // Remove node locations that have boundary conditions identified.
                    bdry_node_locs.erase(std::remove(bdry_node_locs.begin(), bdry_node_locs.end(), node_loc),
                        bdry_node_locs.end());
                }
                else if (bdry_node_conds[node_loc] == BDRY_COND::FLOW_MODEL::XISOTHERMAL_NO_SLIP)
                {
                    for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                    {
                        for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                        {
                            const int idx_cell_rho = (i + num_subghosts_conservative_var[0][0]) +
                                (j + num_subghosts_conservative_var[0][1])*
                                    subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_mom = (i + num_subghosts_conservative_var[1][0]) +
                                (j + num_subghosts_conservative_var[1][1])*
                                    subghostcell_dims_conservative_var[1][0];
                            
                            const int idx_cell_E = (i + num_subghosts_conservative_var[2][0]) +
                                (j + num_subghosts_conservative_var[2][1])*
                                    subghostcell_dims_conservative_var[2][0];
                            
                            int idx_cell_pivot_rho = idx_cell_rho;
                            int idx_cell_pivot_mom = idx_cell_mom;
                            int idx_cell_pivot_E = idx_cell_E;
                            
                            if (edge_loc_0 == BDRY_LOC::XLO)
                            {
                                idx_cell_pivot_rho = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                        num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0];
                                
                                idx_cell_pivot_mom = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                        num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0];
                                
                                idx_cell_pivot_E = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                        num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0];
                            }
                            else if (edge_loc_0 == BDRY_LOC::XHI)
                            {
                                idx_cell_pivot_rho = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                        num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0];
                                
                                idx_cell_pivot_mom = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                        num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0];
                                
                                idx_cell_pivot_E = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                        num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0];
                            }
                            
                            /*
                             * Set the values for density, momentum and total internal energy.
                             */
                            
                            Real epsilon_pivot = (Q[3][idx_cell_pivot_E] -
                                Real(1)/Real(2)*(Q[1][idx_cell_pivot_mom]*Q[1][idx_cell_pivot_mom] +
                                Q[2][idx_cell_pivot_mom]*Q[2][idx_cell_pivot_mom])/
                                Q[0][idx_cell_pivot_rho])/Q[0][idx_cell_pivot_rho];
                            
                            Real p_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getPressure(
                                    &Q[0][idx_cell_pivot_rho],
                                    &epsilon_pivot,
                                    thermo_properties_ptr);
                            
                            Real p = p_pivot;
                            
                            Real T_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getTemperature(
                                    &Q[0][idx_cell_pivot_rho],
                                    &p_pivot,
                                    thermo_properties_ptr);
                            
                            Real T = -T_pivot + Real(2)*d_bdry_edge_isothermal_no_slip_T[edge_loc_0];
                            
                            Real rho = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getDensity(
                                    &p,
                                    &T,
                                    thermo_properties_ptr);
                            
                            Real u = -Q[1][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                                Real(2)*d_bdry_edge_isothermal_no_slip_vel[edge_loc_0*2];
                            Real v = -Q[2][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                                Real(2)*d_bdry_edge_isothermal_no_slip_vel[edge_loc_0*2 + 1];
                            
                            Q[0][idx_cell_rho] = rho;
                            Q[1][idx_cell_mom] = rho*u;
                            Q[2][idx_cell_mom] = rho*v;
                            
                            Real epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getInternalEnergyFromTemperature(
                                    &rho,
                                    &T,
                                    thermo_properties_ptr);
                            
                            Real E = rho*epsilon +
                                Real(1)/Real(2)*(Q[1][idx_cell_mom]*Q[1][idx_cell_mom] + Q[2][idx_cell_mom]*Q[2][idx_cell_mom])/
                                rho;
                            
                            Q[3][idx_cell_E] = E;
                        }
                    }
                    
                    // Remove node locations that have boundary conditions identified.
                    bdry_node_locs.erase(std::remove(bdry_node_locs.begin(), bdry_node_locs.end(), node_loc),
                        bdry_node_locs.end());
                }
                else if (bdry_node_conds[node_loc] == BDRY_COND::FLOW_MODEL::YISOTHERMAL_NO_SLIP)
                {
                    for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                    {
                        for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                        {
                            const int idx_cell_rho = (i + num_subghosts_conservative_var[0][0]) +
                                (j + num_subghosts_conservative_var[0][1])*
                                    subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_mom = (i + num_subghosts_conservative_var[1][0]) +
                                (j + num_subghosts_conservative_var[1][1])*
                                    subghostcell_dims_conservative_var[1][0];
                            
                            const int idx_cell_E = (i + num_subghosts_conservative_var[2][0]) +
                                (j + num_subghosts_conservative_var[2][1])*
                                    subghostcell_dims_conservative_var[2][0];
                            
                            int idx_cell_pivot_rho = idx_cell_rho;
                            int idx_cell_pivot_mom = idx_cell_mom;
                            int idx_cell_pivot_E = idx_cell_E;
                            
                            if (edge_loc_1 == BDRY_LOC::YLO)
                            {
                                idx_cell_pivot_rho = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                        num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0];
                                
                                idx_cell_pivot_mom = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                        num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0];
                                
                                idx_cell_pivot_E = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                        num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0];
                            }
                            else if (edge_loc_1 == BDRY_LOC::YHI)
                            {
                                idx_cell_pivot_rho = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                        num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0];
                                
                                idx_cell_pivot_mom = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                        num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0];
                                
                                idx_cell_pivot_E = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                        num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0];
                            }
                            
                            /*
                             * Set the values for density, momentum and total internal energy.
                             */
                            
                            Real epsilon_pivot = (Q[3][idx_cell_pivot_E] -
                                Real(1)/Real(2)*(Q[1][idx_cell_pivot_mom]*Q[1][idx_cell_pivot_mom] +
                                Q[2][idx_cell_pivot_mom]*Q[2][idx_cell_pivot_mom])/
                                Q[0][idx_cell_pivot_rho])/Q[0][idx_cell_pivot_rho];
                            
                            Real p_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getPressure(
                                    &Q[0][idx_cell_pivot_rho],
                                    &epsilon_pivot,
                                    thermo_properties_ptr);
                            
                            Real p = p_pivot;
                            
                            Real T_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getTemperature(
                                    &Q[0][idx_cell_pivot_rho],
                                    &p_pivot,
                                    thermo_properties_ptr);
                            
                            Real T = -T_pivot + Real(2)*d_bdry_edge_isothermal_no_slip_T[edge_loc_1];
                            
                            Real rho = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getDensity(
                                    &p,
                                    &T,
                                    thermo_properties_ptr);
                            
                            Real u = -Q[1][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                                Real(2)*d_bdry_edge_isothermal_no_slip_vel[edge_loc_1*2];
                            Real v = -Q[2][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                                Real(2)*d_bdry_edge_isothermal_no_slip_vel[edge_loc_1*2 + 1];
                            
                            Q[0][idx_cell_rho] = rho;
                            Q[1][idx_cell_mom] = rho*u;
                            Q[2][idx_cell_mom] = rho*v;
                            
                            Real epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getInternalEnergyFromTemperature(
                                    &rho,
                                    &T,
                                    thermo_properties_ptr);
                            
                            Real E = rho*epsilon +
                                Real(1)/Real(2)*(Q[1][idx_cell_mom]*Q[1][idx_cell_mom] + Q[2][idx_cell_mom]*Q[2][idx_cell_mom])/
                                rho;
                            
                            Q[3][idx_cell_E] = E;
                        }
                    }
                    
                    // Remove node locations that have boundary conditions identified.
                    bdry_node_locs.erase(std::remove(bdry_node_locs.begin(), bdry_node_locs.end(), node_loc),
                        bdry_node_locs.end());
                }
            }
        }
    }
    
    for (int ni = 0; ni < static_cast<int>(node_bdry.size()); ni++)
    {
        TBOX_ASSERT(node_bdry[ni].getBoundaryType() == BDRY::NODE2D);
        
        int node_loc = node_bdry[ni].getLocationIndex();
        
        if (std::find(bdry_node_locs.begin(), bdry_node_locs.end(), node_loc) !=
            bdry_node_locs.end())
        {
            TBOX_ERROR("FlowModelBoundaryUtilitiesSingleSpecies::fill2dNodeBoundaryData\n"
                << "Invalid node boundary condition!\n"
                << "node_loc = '" << node_loc << "'." << std::endl
                << "bdry_node_conds[node_loc] = '" << bdry_node_conds[node_loc] << "'."
                << std::endl);
        }
    }
}


/*
 * Function to fill 3d face boundary values for a patch.
 * Face locations that have boundary conditions identified are removed from the container.
 */
void
FlowModelBoundaryUtilitiesSingleSpecies::fill3dFaceBoundaryData(
    const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& conservative_var_data,
    const hier::Patch& patch,
    std::vector<int>& bdry_face_locs,
    const std::vector<int>& bdry_face_conds,
    const std::vector<std::vector<Real> >& bdry_face_values,
    const hier::IntVector& ghost_width_to_fill)
{
    TBOX_ASSERT(static_cast<int>(conservative_var_data.size()) == 3);
    for (int vi = 0; vi < static_cast<int>(conservative_var_data.size()); vi++)
    {
        TBOX_ASSERT(conservative_var_data[vi]);
    }
    TBOX_ASSERT(static_cast<int>(bdry_face_locs.size()) <= NUM_3D_FACES);
    if (static_cast<int>(bdry_face_locs.size()) > 0)
    {
        TBOX_ASSERT(*min_element(bdry_face_locs.begin(), bdry_face_locs.end()) >= 0);
        TBOX_ASSERT(*max_element(bdry_face_locs.begin(), bdry_face_locs.end()) < NUM_3D_FACES);
    }
    TBOX_ASSERT(static_cast<int>(bdry_face_conds.size()) == NUM_3D_FACES);
    TBOX_ASSERT(static_cast<int>(bdry_face_values.size()) == 3);
    for (int vi = 0; vi < static_cast<int>(bdry_face_values.size()); vi++)
    {
        TBOX_ASSERT(static_cast<int>(bdry_face_values[vi].size()) ==
                    NUM_3D_FACES*(conservative_var_data[vi]->getDepth()));
    }
    
    TBOX_DIM_ASSERT(ghost_width_to_fill.getDim() == tbox::Dimension(3));
    
    for (int vi = 0; vi < static_cast<int>(conservative_var_data.size()); vi++)
    {
        TBOX_ASSERT_OBJDIM_EQUALITY3(*conservative_var_data[vi], patch, ghost_width_to_fill);
    }
    
    const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
        HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    TBOX_ASSERT(patch_geom);
    
    hier::IntVector num_ghosts(conservative_var_data[0]->getGhostCellWidth());
    for (int vi = 1; vi < static_cast<int>(conservative_var_data.size()); vi++)
    {
        num_ghosts = hier::IntVector::min(
            num_ghosts,
            conservative_var_data[vi]->getGhostCellWidth());
    }
    
    /*
     * Determine the ghost cell width to fill.
     */
    
    hier::IntVector gcw_to_fill(tbox::Dimension(3));
    
    // If the ghost fill width is not used, it is set to the ghost cell width of the data.
    if (ghost_width_to_fill == -hier::IntVector::getOne(tbox::Dimension(3)))
    {
        gcw_to_fill = num_ghosts;
    }
    else
    {
        gcw_to_fill = hier::IntVector::min(
            num_ghosts,
            ghost_width_to_fill);
    }
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box& interior_box(patch.getBox());
    hier::Index interior_box_lo_idx(interior_box.lower());
    hier::Index interior_box_hi_idx(interior_box.upper());
    
    /*
     * Offset the indices.
     */
    
    interior_box_lo_idx = interior_box_lo_idx - interior_box.lower();
    interior_box_hi_idx = interior_box_hi_idx - interior_box.lower();
    
    const std::vector<hier::BoundaryBox>& face_bdry =
        patch_geom->getCodimensionBoundaries(BDRY::FACE3D);
    
    /*
     * 1. Fill ghost cells for faces with regular boundary conditions.
     * 2. Fill ghost cells for faces with non-reflecting boundary conditions.
     *    (assuming non-reflecting boundaries don't touch each other)
     */
    
    /*
     * 1. Fill ghost cells for faces with regular boundary conditions.
     */
    for (int fi = 0; fi < static_cast<int>(face_bdry.size()); fi++)
    {
        TBOX_ASSERT(face_bdry[fi].getBoundaryType() == BDRY::FACE3D);
        
        int face_loc = face_bdry[fi].getLocationIndex();
        
        if (std::find(bdry_face_locs.begin(), bdry_face_locs.end(), face_loc) !=
            bdry_face_locs.end())
        {
            hier::Box fill_box(patch_geom->getBoundaryFillBox(
                face_bdry[fi],
                interior_box,
                gcw_to_fill));
            
            hier::Index fill_box_lo_idx(fill_box.lower());
            hier::Index fill_box_hi_idx(fill_box.upper());
            
            /*
             * Offset the indices.
             */
            
            fill_box_lo_idx = fill_box_lo_idx - interior_box.lower();
            fill_box_hi_idx = fill_box_hi_idx - interior_box.lower();
            
            if ((bdry_face_conds[face_loc] == BDRY_COND::FLOW_MODEL::ADIABATIC_NO_SLIP) ||
                (bdry_face_conds[face_loc] == BDRY_COND::FLOW_MODEL::ISOTHERMAL_NO_SLIP))
            {
                /*
                 * Get the pointers to the conservative variables.
                 * Get the numbers of ghost cells and the dimensions of the ghost cell boxes of
                 * the conservative variables.
                 */
                
                std::vector<Real*> Q;
                Q.reserve(d_num_eqn);
                
                std::vector<hier::IntVector> num_subghosts_conservative_var;
                std::vector<hier::IntVector> subghostcell_dims_conservative_var;
                
                num_subghosts_conservative_var.reserve(3);
                subghostcell_dims_conservative_var.reserve(3);
                
                int count_eqn = 0;
                
                for (int vi = 0; vi < static_cast<int>(conservative_var_data.size()); vi++)
                {
                    int depth = conservative_var_data[vi]->getDepth();
                    
                    for (int di = 0; di < depth; di++)
                    {
                        // If the last element of the conservative variable vector is not in the
                        // system of equations, ignore it.
                        if (count_eqn >= d_num_eqn)
                            break;
                        
                        Q.push_back(conservative_var_data[vi]->getPointer(di));
                        
                        count_eqn++;
                    }
                    
                    num_subghosts_conservative_var.push_back(
                        conservative_var_data[vi]->getGhostCellWidth());
                    subghostcell_dims_conservative_var.push_back(
                        conservative_var_data[vi]->getGhostBox().numberCells());
                }
                
                // Get the thermodynamic properties of the species.
                std::vector<const Real*> thermo_properties_ptr;
                thermo_properties_ptr.reserve(static_cast<int> (d_thermo_properties.size()));
                for (int ti = 0; ti < static_cast<int> (d_thermo_properties.size()); ti++)
                {
                    thermo_properties_ptr.push_back(&d_thermo_properties[ti]);
                }
                
                if (bdry_face_conds[face_loc] == BDRY_COND::FLOW_MODEL::ADIABATIC_NO_SLIP)
                {
                    for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                    {
                        for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                        {
                            for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                            {
                                const int idx_cell_rho = (i + num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0] +
                                    (k + num_subghosts_conservative_var[0][2])*
                                        subghostcell_dims_conservative_var[0][0]*
                                            subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom = (i + num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0] +
                                    (k + num_subghosts_conservative_var[1][2])*
                                        subghostcell_dims_conservative_var[1][0]*
                                            subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E = (i + num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0] +
                                    (k + num_subghosts_conservative_var[2][2])*
                                        subghostcell_dims_conservative_var[2][0]*
                                            subghostcell_dims_conservative_var[2][1];
                                
                                int idx_cell_pivot_rho = idx_cell_rho;
                                int idx_cell_pivot_mom = idx_cell_mom;
                                int idx_cell_pivot_E = idx_cell_E;
                                
                                if (face_loc == BDRY_LOC::XLO)
                                {
                                    idx_cell_pivot_rho = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                            num_subghosts_conservative_var[0][0]) +
                                        (j + num_subghosts_conservative_var[0][1])*
                                            subghostcell_dims_conservative_var[0][0] +
                                        (k + num_subghosts_conservative_var[0][2])*
                                            subghostcell_dims_conservative_var[0][0]*
                                                subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                            num_subghosts_conservative_var[1][0]) +
                                        (j + num_subghosts_conservative_var[1][1])*
                                            subghostcell_dims_conservative_var[1][0] +
                                        (k + num_subghosts_conservative_var[1][2])*
                                            subghostcell_dims_conservative_var[1][0]*
                                                subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                            num_subghosts_conservative_var[2][0]) +
                                        (j + num_subghosts_conservative_var[2][1])*
                                            subghostcell_dims_conservative_var[2][0] +
                                        (k + num_subghosts_conservative_var[2][2])*
                                            subghostcell_dims_conservative_var[2][0]*
                                                subghostcell_dims_conservative_var[2][1];
                                }
                                else if (face_loc == BDRY_LOC::XHI)
                                {
                                    idx_cell_pivot_rho = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                            num_subghosts_conservative_var[0][0]) +
                                        (j + num_subghosts_conservative_var[0][1])*
                                            subghostcell_dims_conservative_var[0][0] +
                                        (k + num_subghosts_conservative_var[0][2])*
                                            subghostcell_dims_conservative_var[0][0]*
                                                subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                            num_subghosts_conservative_var[1][0]) +
                                        (j + num_subghosts_conservative_var[1][1])*
                                            subghostcell_dims_conservative_var[1][0] +
                                        (k + num_subghosts_conservative_var[1][2])*
                                            subghostcell_dims_conservative_var[1][0]*
                                                subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                            num_subghosts_conservative_var[2][0]) +
                                        (j + num_subghosts_conservative_var[2][1])*
                                            subghostcell_dims_conservative_var[2][0] +
                                        (k + num_subghosts_conservative_var[2][2])*
                                            subghostcell_dims_conservative_var[2][0]*
                                                subghostcell_dims_conservative_var[2][1];
                                }
                                else if (face_loc == BDRY_LOC::YLO)
                                {
                                    idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                        (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) +
                                            num_subghosts_conservative_var[0][1])*
                                                subghostcell_dims_conservative_var[0][0] +
                                        (k + num_subghosts_conservative_var[0][2])*
                                            subghostcell_dims_conservative_var[0][0]*
                                                subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                        (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) +
                                            num_subghosts_conservative_var[1][1])*
                                                subghostcell_dims_conservative_var[1][0] +
                                        (k + num_subghosts_conservative_var[1][2])*
                                            subghostcell_dims_conservative_var[1][0]*
                                                subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                        (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) +
                                            num_subghosts_conservative_var[2][1])*
                                                subghostcell_dims_conservative_var[2][0] +
                                        (k + num_subghosts_conservative_var[2][2])*
                                            subghostcell_dims_conservative_var[2][0]*
                                                subghostcell_dims_conservative_var[2][1];
                                }
                                else if (face_loc == BDRY_LOC::YHI)
                                {
                                    idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                        (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) +
                                            num_subghosts_conservative_var[0][1])*
                                                subghostcell_dims_conservative_var[0][0] +
                                        (k + num_subghosts_conservative_var[0][2])*
                                            subghostcell_dims_conservative_var[0][0]*
                                                subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                        (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) +
                                            num_subghosts_conservative_var[1][1])*
                                                subghostcell_dims_conservative_var[1][0] +
                                        (k + num_subghosts_conservative_var[1][2])*
                                            subghostcell_dims_conservative_var[1][0]*
                                                subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                        (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) +
                                            num_subghosts_conservative_var[2][1])*
                                                subghostcell_dims_conservative_var[2][0] +
                                        (k + num_subghosts_conservative_var[2][2])*
                                            subghostcell_dims_conservative_var[2][0]*
                                                subghostcell_dims_conservative_var[2][1];
                                }
                                else if (face_loc == BDRY_LOC::ZLO)
                                {
                                    idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                        (j + num_subghosts_conservative_var[0][1])*
                                            subghostcell_dims_conservative_var[0][0] +
                                        (interior_box_lo_idx[2] + (fill_box_hi_idx[2] - k) +
                                            num_subghosts_conservative_var[0][2])*
                                                subghostcell_dims_conservative_var[0][0]*
                                                    subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                        (j + num_subghosts_conservative_var[1][1])*
                                            subghostcell_dims_conservative_var[1][0] +
                                        (interior_box_lo_idx[2] + (fill_box_hi_idx[2] - k) +
                                            num_subghosts_conservative_var[1][2])*
                                                subghostcell_dims_conservative_var[1][0]*
                                                    subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                        (j + num_subghosts_conservative_var[2][1])*
                                            subghostcell_dims_conservative_var[2][0] +
                                        (interior_box_lo_idx[2] + (fill_box_hi_idx[2] - k) +
                                            num_subghosts_conservative_var[2][2])*
                                                subghostcell_dims_conservative_var[2][0]*
                                                    subghostcell_dims_conservative_var[2][1];
                                }
                                else if (face_loc == BDRY_LOC::ZHI)
                                {
                                    idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                        (j + num_subghosts_conservative_var[0][1])*
                                            subghostcell_dims_conservative_var[0][0] +
                                        (interior_box_hi_idx[2] - (k - fill_box_lo_idx[2]) +
                                            num_subghosts_conservative_var[0][2])*
                                                subghostcell_dims_conservative_var[0][0]*
                                                    subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                        (j + num_subghosts_conservative_var[1][1])*
                                            subghostcell_dims_conservative_var[1][0] +
                                        (interior_box_hi_idx[2] - (k - fill_box_lo_idx[2]) +
                                            num_subghosts_conservative_var[1][2])*
                                                subghostcell_dims_conservative_var[1][0]*
                                                    subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                        (j + num_subghosts_conservative_var[2][1])*
                                            subghostcell_dims_conservative_var[2][0] +
                                        (interior_box_hi_idx[2] - (k - fill_box_lo_idx[2]) +
                                            num_subghosts_conservative_var[2][2])*
                                                subghostcell_dims_conservative_var[2][0]*
                                                    subghostcell_dims_conservative_var[2][1];
                                }
                                
                                /*
                                 * Set the values for density and momentum.
                                 */
                                
                                Q[0][idx_cell_rho] = Q[0][idx_cell_pivot_rho];
                                Q[1][idx_cell_mom] = -Q[1][idx_cell_pivot_mom] +
                                    Real(2)*Q[0][idx_cell_pivot_rho]*d_bdry_face_adiabatic_no_slip_vel[face_loc*3];
                                Q[2][idx_cell_mom] = -Q[2][idx_cell_pivot_mom] +
                                    Real(2)*Q[0][idx_cell_pivot_rho]*d_bdry_face_adiabatic_no_slip_vel[face_loc*3 + 1];
                                Q[3][idx_cell_mom] = -Q[3][idx_cell_pivot_mom] +
                                    Real(2)*Q[0][idx_cell_pivot_rho]*d_bdry_face_adiabatic_no_slip_vel[face_loc*3 + 2];
                                
                                /*
                                 * Set the values for total internal energy.
                                 */
                                
                                Real epsilon_pivot = (Q[4][idx_cell_pivot_E] -
                                    Real(1)/Real(2)*(Q[1][idx_cell_pivot_mom]*Q[1][idx_cell_pivot_mom] +
                                    Q[2][idx_cell_pivot_mom]*Q[2][idx_cell_pivot_mom] +
                                    Q[3][idx_cell_pivot_mom]*Q[3][idx_cell_pivot_mom])/
                                    Q[0][idx_cell_pivot_rho])/Q[0][idx_cell_pivot_rho];
                                
                                Real p_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &Q[0][idx_cell_pivot_rho],
                                        &epsilon_pivot,
                                        thermo_properties_ptr);
                                
                                Real T_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getTemperature(
                                        &Q[0][idx_cell_pivot_rho],
                                        &p_pivot,
                                        thermo_properties_ptr);
                                
                                Real T = T_pivot;
                                
                                Real epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getInternalEnergyFromTemperature(
                                        &Q[0][idx_cell_rho],
                                        &T,
                                        thermo_properties_ptr);
                                
                                Real E = Q[0][idx_cell_rho]*epsilon +
                                    Real(1)/Real(2)*(Q[1][idx_cell_mom]*Q[1][idx_cell_mom] + Q[2][idx_cell_mom]*Q[2][idx_cell_mom] +
                                    Q[3][idx_cell_mom]*Q[3][idx_cell_mom])/
                                    Q[0][idx_cell_rho];
                                
                                Q[4][idx_cell_E] = E;
                            }
                        }
                    }
                    
                    // Remove face locations that have boundary conditions identified.
                    bdry_face_locs.erase(std::remove(bdry_face_locs.begin(), bdry_face_locs.end(), face_loc),
                        bdry_face_locs.end());
                }
                else if (bdry_face_conds[face_loc] == BDRY_COND::FLOW_MODEL::ISOTHERMAL_NO_SLIP)
                {
                    for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                    {
                        for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                        {
                            for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                            {
                                const int idx_cell_rho = (i + num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0] +
                                    (k + num_subghosts_conservative_var[0][2])*
                                        subghostcell_dims_conservative_var[0][0]*
                                            subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom = (i + num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0] +
                                    (k + num_subghosts_conservative_var[1][2])*
                                        subghostcell_dims_conservative_var[1][0]*
                                            subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E = (i + num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0] +
                                    (k + num_subghosts_conservative_var[2][2])*
                                        subghostcell_dims_conservative_var[2][0]*
                                            subghostcell_dims_conservative_var[2][1];
                                
                                int idx_cell_pivot_rho = idx_cell_rho;
                                int idx_cell_pivot_mom = idx_cell_mom;
                                int idx_cell_pivot_E = idx_cell_E;
                                
                                if (face_loc == BDRY_LOC::XLO)
                                {
                                    idx_cell_pivot_rho = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                            num_subghosts_conservative_var[0][0]) +
                                        (j + num_subghosts_conservative_var[0][1])*
                                            subghostcell_dims_conservative_var[0][0] +
                                        (k + num_subghosts_conservative_var[0][2])*
                                            subghostcell_dims_conservative_var[0][0]*
                                                subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                            num_subghosts_conservative_var[1][0]) +
                                        (j + num_subghosts_conservative_var[1][1])*
                                            subghostcell_dims_conservative_var[1][0] +
                                        (k + num_subghosts_conservative_var[1][2])*
                                            subghostcell_dims_conservative_var[1][0]*
                                                subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                            num_subghosts_conservative_var[2][0]) +
                                        (j + num_subghosts_conservative_var[2][1])*
                                            subghostcell_dims_conservative_var[2][0] +
                                        (k + num_subghosts_conservative_var[2][2])*
                                            subghostcell_dims_conservative_var[2][0]*
                                                subghostcell_dims_conservative_var[2][1];
                                }
                                else if (face_loc == BDRY_LOC::XHI)
                                {
                                    idx_cell_pivot_rho = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                            num_subghosts_conservative_var[0][0]) +
                                        (j + num_subghosts_conservative_var[0][1])*
                                            subghostcell_dims_conservative_var[0][0] +
                                        (k + num_subghosts_conservative_var[0][2])*
                                            subghostcell_dims_conservative_var[0][0]*
                                                subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                            num_subghosts_conservative_var[1][0]) +
                                        (j + num_subghosts_conservative_var[1][1])*
                                            subghostcell_dims_conservative_var[1][0] +
                                        (k + num_subghosts_conservative_var[1][2])*
                                            subghostcell_dims_conservative_var[1][0]*
                                                subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                            num_subghosts_conservative_var[2][0]) +
                                        (j + num_subghosts_conservative_var[2][1])*
                                            subghostcell_dims_conservative_var[2][0] +
                                        (k + num_subghosts_conservative_var[2][2])*
                                            subghostcell_dims_conservative_var[2][0]*
                                                subghostcell_dims_conservative_var[2][1];
                                }
                                else if (face_loc == BDRY_LOC::YLO)
                                {
                                    idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                        (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) +
                                            num_subghosts_conservative_var[0][1])*
                                                subghostcell_dims_conservative_var[0][0] +
                                        (k + num_subghosts_conservative_var[0][2])*
                                            subghostcell_dims_conservative_var[0][0]*
                                                subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                        (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) +
                                            num_subghosts_conservative_var[1][1])*
                                                subghostcell_dims_conservative_var[1][0] +
                                        (k + num_subghosts_conservative_var[1][2])*
                                            subghostcell_dims_conservative_var[1][0]*
                                                subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                        (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) +
                                            num_subghosts_conservative_var[2][1])*
                                                subghostcell_dims_conservative_var[2][0] +
                                        (k + num_subghosts_conservative_var[2][2])*
                                            subghostcell_dims_conservative_var[2][0]*
                                                subghostcell_dims_conservative_var[2][1];
                                }
                                else if (face_loc == BDRY_LOC::YHI)
                                {
                                    idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                        (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) +
                                            num_subghosts_conservative_var[0][1])*
                                                subghostcell_dims_conservative_var[0][0] +
                                        (k + num_subghosts_conservative_var[0][2])*
                                            subghostcell_dims_conservative_var[0][0]*
                                                subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                        (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) +
                                            num_subghosts_conservative_var[1][1])*
                                                subghostcell_dims_conservative_var[1][0] +
                                        (k + num_subghosts_conservative_var[1][2])*
                                            subghostcell_dims_conservative_var[1][0]*
                                                subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                        (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) +
                                            num_subghosts_conservative_var[2][1])*
                                                subghostcell_dims_conservative_var[2][0] +
                                        (k + num_subghosts_conservative_var[2][2])*
                                            subghostcell_dims_conservative_var[2][0]*
                                                subghostcell_dims_conservative_var[2][1];
                                }
                                else if (face_loc == BDRY_LOC::ZLO)
                                {
                                    idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                        (j + num_subghosts_conservative_var[0][1])*
                                            subghostcell_dims_conservative_var[0][0] +
                                        (interior_box_lo_idx[2] + (fill_box_hi_idx[2] - k) +
                                            num_subghosts_conservative_var[0][2])*
                                                subghostcell_dims_conservative_var[0][0]*
                                                    subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                        (j + num_subghosts_conservative_var[1][1])*
                                            subghostcell_dims_conservative_var[1][0] +
                                        (interior_box_lo_idx[2] + (fill_box_hi_idx[2] - k) +
                                            num_subghosts_conservative_var[1][2])*
                                                subghostcell_dims_conservative_var[1][0]*
                                                    subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                        (j + num_subghosts_conservative_var[2][1])*
                                            subghostcell_dims_conservative_var[2][0] +
                                        (interior_box_lo_idx[2] + (fill_box_hi_idx[2] - k) +
                                            num_subghosts_conservative_var[2][2])*
                                                subghostcell_dims_conservative_var[2][0]*
                                                    subghostcell_dims_conservative_var[2][1];
                                }
                                else if (face_loc == BDRY_LOC::ZHI)
                                {
                                    idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                        (j + num_subghosts_conservative_var[0][1])*
                                            subghostcell_dims_conservative_var[0][0] +
                                        (interior_box_hi_idx[2] - (k - fill_box_lo_idx[2]) +
                                            num_subghosts_conservative_var[0][2])*
                                                subghostcell_dims_conservative_var[0][0]*
                                                    subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                        (j + num_subghosts_conservative_var[1][1])*
                                            subghostcell_dims_conservative_var[1][0] +
                                        (interior_box_hi_idx[2] - (k - fill_box_lo_idx[2]) +
                                            num_subghosts_conservative_var[1][2])*
                                                subghostcell_dims_conservative_var[1][0]*
                                                    subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                        (j + num_subghosts_conservative_var[2][1])*
                                            subghostcell_dims_conservative_var[2][0] +
                                        (interior_box_hi_idx[2] - (k - fill_box_lo_idx[2]) +
                                            num_subghosts_conservative_var[2][2])*
                                                subghostcell_dims_conservative_var[2][0]*
                                                    subghostcell_dims_conservative_var[2][1];
                                }
                                
                                /*
                                 * Set the values for density, momentum and total internal energy.
                                 */
                                
                                Real epsilon_pivot = (Q[4][idx_cell_pivot_E] -
                                    Real(1)/Real(2)*(Q[1][idx_cell_pivot_mom]*Q[1][idx_cell_pivot_mom] +
                                    Q[2][idx_cell_pivot_mom]*Q[2][idx_cell_pivot_mom] +
                                    Q[3][idx_cell_pivot_mom]*Q[3][idx_cell_pivot_mom])/
                                    Q[0][idx_cell_pivot_rho])/Q[0][idx_cell_pivot_rho];
                                
                                Real p_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &Q[0][idx_cell_pivot_rho],
                                        &epsilon_pivot,
                                        thermo_properties_ptr);
                                
                                Real p = p_pivot;
                                
                                Real T_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getTemperature(
                                        &Q[0][idx_cell_pivot_rho],
                                        &p_pivot,
                                        thermo_properties_ptr);
                                
                                Real T = -T_pivot + Real(2)*d_bdry_face_isothermal_no_slip_T[face_loc];
                                
                                Real rho = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getDensity(
                                        &p,
                                        &T,
                                        thermo_properties_ptr);
                                
                                Real u = -Q[1][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                                    Real(2)*d_bdry_face_isothermal_no_slip_vel[face_loc*3];
                                Real v = -Q[2][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                                    Real(2)*d_bdry_face_isothermal_no_slip_vel[face_loc*3 + 1];
                                Real w = -Q[3][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                                    Real(2)*d_bdry_face_isothermal_no_slip_vel[face_loc*3 + 2];
                                
                                Q[0][idx_cell_rho] = rho;
                                Q[1][idx_cell_mom] = rho*u;
                                Q[2][idx_cell_mom] = rho*v;
                                Q[3][idx_cell_mom] = rho*w;
                                
                                Real epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getInternalEnergyFromTemperature(
                                        &rho,
                                        &T,
                                        thermo_properties_ptr);
                                
                                Real E = rho*epsilon +
                                    Real(1)/Real(2)*(Q[1][idx_cell_mom]*Q[1][idx_cell_mom] + Q[2][idx_cell_mom]*Q[2][idx_cell_mom] +
                                    Q[3][idx_cell_mom]*Q[3][idx_cell_mom])/rho;
                                
                                Q[4][idx_cell_E] = E;
                            }
                        }
                    }
                    
                    // Remove face locations that have boundary conditions identified.
                    bdry_face_locs.erase(std::remove(bdry_face_locs.begin(), bdry_face_locs.end(), face_loc),
                        bdry_face_locs.end());
                }
            }
        }
    }
    
    /*
     * 2. Fill ghost cells for faces with non-reflecting boundary conditions.
     */
    for (int fi = 0; fi < static_cast<int>(face_bdry.size()); fi++)
    {
        TBOX_ASSERT(face_bdry[fi].getBoundaryType() == BDRY::FACE3D);
        
        int face_loc = face_bdry[fi].getLocationIndex();
        
        if (std::find(bdry_face_locs.begin(), bdry_face_locs.end(), face_loc) !=
            bdry_face_locs.end())
        {
            if (bdry_face_conds[face_loc] == BDRY_COND::FLOW_MODEL::NONREFLECTING_OUTFLOW)
            {
                // Follow the method in
                // Motheau, Emmanuel, Ann Almgren, and John B. Bell.
                // "Navier–stokes characteristic boundary conditions using ghost cells."
                // AIAA Journal (2017): 3399-3408.
                
                hier::Box fill_box(patch_geom->getBoundaryFillBox(
                    face_bdry[fi],
                    interior_box,
                    gcw_to_fill));
                
                hier::Index fill_box_lo_idx(fill_box.lower());
                hier::Index fill_box_hi_idx(fill_box.upper());
                
                /*
                 * Offset the indices.
                 */
                
                fill_box_lo_idx = fill_box_lo_idx - interior_box.lower();
                fill_box_hi_idx = fill_box_hi_idx - interior_box.lower();
                
                /*
                 * Get the pointers to the conservative variables.
                 * Get the numbers of ghost cells and the dimensions of the ghost cell boxes of
                 * the conservative variables.
                 */
                
                std::vector<Real*> Q;
                Q.reserve(d_num_eqn);
                
                std::vector<hier::IntVector> num_subghosts_conservative_var;
                std::vector<hier::IntVector> subghostcell_dims_conservative_var;
                
                num_subghosts_conservative_var.reserve(3);
                subghostcell_dims_conservative_var.reserve(3);
                
                int count_eqn = 0;
                
                for (int vi = 0; vi < static_cast<int>(conservative_var_data.size()); vi++)
                {
                    int depth = conservative_var_data[vi]->getDepth();
                    
                    for (int di = 0; di < depth; di++)
                    {
                        // If the last element of the conservative variable vector is not in the
                        // system of equations, ignore it.
                        if (count_eqn >= d_num_eqn)
                            break;
                        
                        Q.push_back(conservative_var_data[vi]->getPointer(di));
                        
                        count_eqn++;
                    }
                    
                    num_subghosts_conservative_var.push_back(
                        conservative_var_data[vi]->getGhostCellWidth());
                    subghostcell_dims_conservative_var.push_back(
                        conservative_var_data[vi]->getGhostBox().numberCells());
                }
                
                // Get the thermodynamic properties of the species.
                std::vector<const Real*> thermo_properties_ptr;
                thermo_properties_ptr.reserve(static_cast<int> (d_thermo_properties.size()));
                for (int ti = 0; ti < static_cast<int> (d_thermo_properties.size()); ti++)
                {
                    thermo_properties_ptr.push_back(&d_thermo_properties[ti]);
                }
                
                if (face_loc == BDRY_LOC::XLO)
                {
                    const int num_ghosts_to_fill = fill_box_hi_idx[0] - fill_box_lo_idx[0] + 1;
                    TBOX_ASSERT(fill_box_hi_idx[0] == interior_box_lo_idx[0] - 1);
                    if (num_ghosts_to_fill > 6)
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelBoundaryUtilitiesSingleSpecies::fill3dFaceBoundaryData()\n"
                            << "Non-reflecting outflow BC doesn't support more than six ghost cells yet!");
                    }
                    
                    for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                    {
                        for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                        {
                            // Get the grid spacing.
                            const double* const dx = patch_geom->getDx();
                            
                            const int idx_cell_rho_x_R = (interior_box_lo_idx[0] + num_subghosts_conservative_var[0][0]) +
                                (j + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                (k + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                    subghostcell_dims_conservative_var[0][1];
                            
                            const int idx_cell_rho_x_RR = (interior_box_lo_idx[0] + 1 + num_subghosts_conservative_var[0][0]) +
                                (j + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                (k + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                    subghostcell_dims_conservative_var[0][1];
                            
                            const int idx_cell_rho_x_RRR = (interior_box_lo_idx[0] + 2 + num_subghosts_conservative_var[0][0]) +
                                (j + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                (k + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                    subghostcell_dims_conservative_var[0][1];
                            
                            const int idx_cell_mom_x_R = (interior_box_lo_idx[0] + num_subghosts_conservative_var[1][0]) +
                                (j + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                (k + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                    subghostcell_dims_conservative_var[1][1];
                            
                            const int idx_cell_mom_x_RR = (interior_box_lo_idx[0] + 1 + num_subghosts_conservative_var[1][0]) +
                                (j + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                (k + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                    subghostcell_dims_conservative_var[1][1];
                            
                            const int idx_cell_mom_x_RRR = (interior_box_lo_idx[0] + 2 + num_subghosts_conservative_var[1][0]) +
                                (j + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                (k + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                    subghostcell_dims_conservative_var[1][1];
                            
                            const int idx_cell_E_x_R = (interior_box_lo_idx[0] + num_subghosts_conservative_var[2][0]) +
                                (j + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                (k + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                    subghostcell_dims_conservative_var[2][1];
                            
                            const int idx_cell_E_x_RR = (interior_box_lo_idx[0] + 1 + num_subghosts_conservative_var[2][0]) +
                                (j + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                (k + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                    subghostcell_dims_conservative_var[2][1];
                            
                            const int idx_cell_E_x_RRR = (interior_box_lo_idx[0] + 2 + num_subghosts_conservative_var[2][0]) +
                                (j + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                (k + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                    subghostcell_dims_conservative_var[2][1];
                            
                            const Real& rho_x_R   = Q[0][idx_cell_rho_x_R];
                            const Real& rho_x_RR  = Q[0][idx_cell_rho_x_RR];
                            const Real& rho_x_RRR = Q[0][idx_cell_rho_x_RRR];
                            
                            const Real u_x_R   = Q[1][idx_cell_mom_x_R]/rho_x_R;
                            const Real u_x_RR  = Q[1][idx_cell_mom_x_RR]/rho_x_RR;
                            const Real u_x_RRR = Q[1][idx_cell_mom_x_RRR]/rho_x_RRR;
                            
                            const Real v_x_R   = Q[2][idx_cell_mom_x_R]/rho_x_R;
                            const Real v_x_RR  = Q[2][idx_cell_mom_x_RR]/rho_x_RR;
                            const Real v_x_RRR = Q[2][idx_cell_mom_x_RRR]/rho_x_RRR;
                            
                            const Real w_x_R   = Q[3][idx_cell_mom_x_R]/rho_x_R;
                            const Real w_x_RR  = Q[3][idx_cell_mom_x_RR]/rho_x_RR;
                            const Real w_x_RRR = Q[3][idx_cell_mom_x_RRR]/rho_x_RRR;
                            
                            const Real half = Real(1)/Real(2);
                            const Real epsilon_x_R   = Q[4][idx_cell_E_x_R]/rho_x_R - half*(u_x_R*u_x_R + v_x_R*v_x_R + w_x_R*w_x_R);
                            const Real epsilon_x_RR  = Q[4][idx_cell_E_x_RR]/rho_x_RR - half*(u_x_RR*u_x_RR + v_x_RR*v_x_RR + w_x_RR*w_x_RR);
                            const Real epsilon_x_RRR = Q[4][idx_cell_E_x_RRR]/rho_x_RRR - half*(u_x_RRR*u_x_RRR + v_x_RRR*v_x_RRR + w_x_RRR*w_x_RRR);
                            
                            const Real p_x_R = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getPressure(
                                    &rho_x_R,
                                    &epsilon_x_R,
                                    thermo_properties_ptr);
                            
                            const Real p_x_RR = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getPressure(
                                    &rho_x_RR,
                                    &epsilon_x_RR,
                                    thermo_properties_ptr);
                            
                            const Real p_x_RRR = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getPressure(
                                    &rho_x_RRR,
                                    &epsilon_x_RRR,
                                    thermo_properties_ptr);
                            
                            /*
                             * Compute derivatives in x-direction.
                             */
                            
                            const Real drho_dx = -(rho_x_RRR - Real(4)*rho_x_RR + Real(3)*rho_x_R)/(Real(2)*Real(dx[0]));
                            const Real du_dx   = -(u_x_RRR - Real(4)*u_x_RR + Real(3)*u_x_R)/(Real(2)*Real(dx[0]));
                            const Real dv_dx   = -(v_x_RRR - Real(4)*v_x_RR + Real(3)*v_x_R)/(Real(2)*Real(dx[0]));
                            const Real dw_dx   = -(w_x_RRR - Real(4)*w_x_RR + Real(3)*w_x_R)/(Real(2)*Real(dx[0]));
                            const Real dp_dx   = -(p_x_RRR - Real(4)*p_x_RR + Real(3)*p_x_R)/(Real(2)*Real(dx[0]));
                            
                            /*
                             * Compute derivatives in y-direction.
                             */
                            
                            Real du_dy = Real(0);
                            Real dv_dy = Real(0);
                            // Real dw_dy = Real(0);
                            Real dp_dy = Real(0);
                            
                            if ((j + num_subghosts_conservative_var[0][1] == 0) ||
                                (j + num_subghosts_conservative_var[1][1] == 0) ||
                                (j + num_subghosts_conservative_var[2][1] == 0))
                            // if (((patch_geom->getTouchesRegularBoundary(1, 0)) && (j == interior_box_lo_idx[1])) ||
                            //     ((j + num_subghosts_conservative_var[0][1] == 0) ||
                            //      (j + num_subghosts_conservative_var[1][1] == 0) ||
                            //      (j + num_subghosts_conservative_var[2][1] == 0)))
                            {
                                // Patch is touching periodic boundary.
                                
                                const int idx_cell_rho_y_T = (interior_box_lo_idx[0] + num_subghosts_conservative_var[0][0]) +
                                    (j + 1 + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (k + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom_y_T = (interior_box_lo_idx[0] + num_subghosts_conservative_var[1][0]) +
                                    (j + 1 + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                    (k + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E_y_T = (interior_box_lo_idx[0] + num_subghosts_conservative_var[2][0]) +
                                    (j + 1 + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                    (k + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                const Real& rho_y_T = Q[0][idx_cell_rho_y_T];
                                const Real u_y_T = Q[1][idx_cell_mom_y_T]/rho_y_T;
                                const Real v_y_T = Q[2][idx_cell_mom_y_T]/rho_y_T;
                                const Real w_y_T = Q[3][idx_cell_mom_y_T]/rho_y_T;
                                const Real epsilon_y_T = Q[4][idx_cell_E_y_T]/rho_y_T - half*(u_y_T*u_y_T + v_y_T*v_y_T + w_y_T*w_y_T);
                                
                                const Real p_y_T = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &rho_y_T,
                                        &epsilon_y_T,
                                        thermo_properties_ptr);
                                
                                // One-sided derivatives.
                                du_dy = (u_y_T - u_x_R)/Real(dx[1]);
                                dv_dy = (v_y_T - v_x_R)/Real(dx[1]);
                                // dw_dy = (w_y_T - w_x_R)/Real(dx[1]);
                                dp_dy = (p_y_T - p_x_R)/Real(dx[1]);
                            }
                            else if ((j + num_subghosts_conservative_var[0][1] + 1 == subghostcell_dims_conservative_var[0][1]) ||
                                     (j + num_subghosts_conservative_var[1][1] + 1 == subghostcell_dims_conservative_var[1][1]) ||
                                     (j + num_subghosts_conservative_var[2][1] + 1 == subghostcell_dims_conservative_var[2][1]))
                            // else if (((patch_geom->getTouchesRegularBoundary(1, 1)) && (j == interior_box_hi_idx[1])) ||
                            //          ((j + num_subghosts_conservative_var[0][1] + 1 == subghostcell_dims_conservative_var[0][1]) ||
                            //           (j + num_subghosts_conservative_var[1][1] + 1 == subghostcell_dims_conservative_var[1][1]) ||
                            //           (j + num_subghosts_conservative_var[2][1] + 1 == subghostcell_dims_conservative_var[2][1])))
                            {
                                // Patch is touching periodic boundary.
                                
                                const int idx_cell_rho_y_B = (interior_box_lo_idx[0] + num_subghosts_conservative_var[0][0]) +
                                    (j - 1 + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (k + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom_y_B = (interior_box_lo_idx[0] + num_subghosts_conservative_var[1][0]) +
                                    (j - 1 + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                    (k + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E_y_B = (interior_box_lo_idx[0] + num_subghosts_conservative_var[2][0]) +
                                    (j - 1 + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                    (k + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                const Real& rho_y_B = Q[0][idx_cell_rho_y_B];
                                const Real u_y_B = Q[1][idx_cell_mom_y_B]/rho_y_B;
                                const Real v_y_B = Q[2][idx_cell_mom_y_B]/rho_y_B;
                                const Real w_y_B = Q[3][idx_cell_mom_y_B]/rho_y_B;
                                const Real epsilon_y_B = Q[4][idx_cell_E_y_B]/rho_y_B - half*(u_y_B*u_y_B + v_y_B*v_y_B + w_y_B*w_y_B);
                                
                                const Real p_y_B = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &rho_y_B,
                                        &epsilon_y_B,
                                        thermo_properties_ptr);
                                
                                // One-sided derivatives.
                                du_dy = (u_x_R - u_y_B)/Real(dx[1]);
                                dv_dy = (v_x_R - v_y_B)/Real(dx[1]);
                                // dw_dy = (w_x_R - w_y_B)/Real(dx[1]);
                                dp_dy = (p_x_R - p_y_B)/Real(dx[1]);
                            }
                            else
                            {
                                const int idx_cell_rho_y_B = (interior_box_lo_idx[0] + num_subghosts_conservative_var[0][0]) +
                                    (j - 1 + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (k + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_rho_y_T = (interior_box_lo_idx[0] + num_subghosts_conservative_var[0][0]) +
                                    (j + 1 + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (k + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom_y_B = (interior_box_lo_idx[0] + num_subghosts_conservative_var[1][0]) +
                                    (j - 1 + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                    (k + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_mom_y_T = (interior_box_lo_idx[0] + num_subghosts_conservative_var[1][0]) +
                                    (j + 1 + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                    (k + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E_y_B = (interior_box_lo_idx[0] + num_subghosts_conservative_var[2][0]) +
                                    (j - 1 + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                    (k + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                const int idx_cell_E_y_T = (interior_box_lo_idx[0] + num_subghosts_conservative_var[2][0]) +
                                    (j + 1 + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                    (k + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                const Real& rho_y_B = Q[0][idx_cell_rho_y_B];
                                const Real& rho_y_T = Q[0][idx_cell_rho_y_T];
                                
                                const Real u_y_B = Q[1][idx_cell_mom_y_B]/rho_y_B;
                                const Real u_y_T = Q[1][idx_cell_mom_y_T]/rho_y_T;
                                
                                const Real v_y_B = Q[2][idx_cell_mom_y_B]/rho_y_B;
                                const Real v_y_T = Q[2][idx_cell_mom_y_T]/rho_y_T;
                                
                                const Real w_y_B = Q[3][idx_cell_mom_y_B]/rho_y_B;
                                const Real w_y_T = Q[3][idx_cell_mom_y_T]/rho_y_T;
                                
                                const Real epsilon_y_B = Q[4][idx_cell_E_y_B]/rho_y_B - half*(u_y_B*u_y_B + v_y_B*v_y_B + w_y_B*w_y_B);
                                const Real epsilon_y_T = Q[4][idx_cell_E_y_T]/rho_y_T - half*(u_y_T*u_y_T + v_y_T*v_y_T + w_y_T*w_y_T);
                                
                                const Real p_y_B = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &rho_y_B,
                                        &epsilon_y_B,
                                        thermo_properties_ptr);
                                
                                const Real p_y_T = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &rho_y_T,
                                        &epsilon_y_T,
                                        thermo_properties_ptr);
                                
                                // Central derivatives.
                                du_dy = (u_y_T - u_y_B)/(Real(2)*Real(dx[1]));
                                dv_dy = (v_y_T - v_y_B)/(Real(2)*Real(dx[1]));
                                // dw_dy = (w_y_T - w_y_B)/(Real(2)*Real(dx[1]));
                                dp_dy = (p_y_T - p_y_B)/(Real(2)*Real(dx[1]));
                            }
                            
                            /*
                             * Compute derivatives in z-direction.
                             */
                            
                            Real du_dz = Real(0);
                            // Real dv_dz = Real(0);
                            Real dw_dz = Real(0);
                            Real dp_dz = Real(0);
                            
                            if ((k + num_subghosts_conservative_var[0][2] == 0) ||
                                (k + num_subghosts_conservative_var[1][2] == 0) ||
                                (k + num_subghosts_conservative_var[2][2] == 0))
                            // if (((patch_geom->getTouchesRegularBoundary(2, 0)) && (k == interior_box_lo_idx[2])) ||
                            //     ((k + num_subghosts_conservative_var[0][2] == 0) ||
                            //      (k + num_subghosts_conservative_var[1][2] == 0) ||
                            //      (k + num_subghosts_conservative_var[2][2] == 0)))
                            {
                                // Patch is touching periodic boundary.
                                
                                const int idx_cell_rho_z_F = (interior_box_lo_idx[0] + num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (k + 1 + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom_z_F = (interior_box_lo_idx[0] + num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                    (k + 1 + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E_z_F = (interior_box_lo_idx[0] + num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                    (k + 1 + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                const Real& rho_z_F = Q[0][idx_cell_rho_z_F];
                                const Real u_z_F = Q[1][idx_cell_mom_z_F]/rho_z_F;
                                const Real v_z_F = Q[2][idx_cell_mom_z_F]/rho_z_F;
                                const Real w_z_F = Q[3][idx_cell_mom_z_F]/rho_z_F;
                                const Real epsilon_z_F = Q[4][idx_cell_E_z_F]/rho_z_F - half*(u_z_F*u_z_F + v_z_F*v_z_F + w_z_F*w_z_F);
                                
                                const Real p_z_F = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &rho_z_F,
                                        &epsilon_z_F,
                                        thermo_properties_ptr);
                                
                                // One-sided derivatives.
                                du_dz = (u_z_F - u_x_R)/Real(dx[2]);
                                // dv_dz = (v_z_F - v_x_R)/Real(dx[2]);
                                dw_dz = (w_z_F - w_x_R)/Real(dx[2]);
                                dp_dz = (p_z_F - p_x_R)/Real(dx[2]);
                            }
                            else if ((k + num_subghosts_conservative_var[0][2] + 1 == subghostcell_dims_conservative_var[0][2]) ||
                                     (k + num_subghosts_conservative_var[1][2] + 1 == subghostcell_dims_conservative_var[1][2]) ||
                                     (k + num_subghosts_conservative_var[2][2] + 1 == subghostcell_dims_conservative_var[2][2]))
                            // else if (((patch_geom->getTouchesRegularBoundary(2, 1)) && (k == interior_box_hi_idx[2])) ||
                            //          ((k + num_subghosts_conservative_var[0][2] + 1 == subghostcell_dims_conservative_var[0][2]) ||
                            //           (k + num_subghosts_conservative_var[1][2] + 1 == subghostcell_dims_conservative_var[1][2]) ||
                            //           (k + num_subghosts_conservative_var[2][2] + 1 == subghostcell_dims_conservative_var[2][2])))
                            {
                                // Patch is touching periodic boundary.
                                
                                const int idx_cell_rho_z_B = (interior_box_lo_idx[0] + num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (k - 1 + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom_z_B = (interior_box_lo_idx[0] + num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                    (k - 1 + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E_z_B = (interior_box_lo_idx[0] + num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                    (k - 1 + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                const Real& rho_z_B = Q[0][idx_cell_rho_z_B];
                                const Real u_z_B = Q[1][idx_cell_mom_z_B]/rho_z_B;
                                const Real v_z_B = Q[2][idx_cell_mom_z_B]/rho_z_B;
                                const Real w_z_B = Q[3][idx_cell_mom_z_B]/rho_z_B;
                                const Real epsilon_z_B = Q[4][idx_cell_E_z_B]/rho_z_B - half*(u_z_B*u_z_B + v_z_B*v_z_B + w_z_B*w_z_B);
                                
                                const Real p_z_B = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &rho_z_B,
                                        &epsilon_z_B,
                                        thermo_properties_ptr);
                                
                                // One-sided derivatives.
                                du_dz = (u_x_R - u_z_B)/Real(dx[2]);
                                // dv_dz = (v_x_R - v_z_B)/Real(dx[2]);
                                dw_dz = (w_x_R - w_z_B)/Real(dx[2]);
                                dp_dz = (p_x_R - p_z_B)/Real(dx[2]);
                            }
                            else
                            {
                                const int idx_cell_rho_z_B = (interior_box_lo_idx[0] + num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (k - 1 + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_rho_z_F = (interior_box_lo_idx[0] + num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (k + 1 + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom_z_B = (interior_box_lo_idx[0] + num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                    (k - 1 + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_mom_z_F = (interior_box_lo_idx[0] + num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                    (k + 1 + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E_z_B = (interior_box_lo_idx[0] + num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                    (k - 1 + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                const int idx_cell_E_z_F = (interior_box_lo_idx[0] + num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                    (k + 1 + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                const Real& rho_z_B = Q[0][idx_cell_rho_z_B];
                                const Real& rho_z_F = Q[0][idx_cell_rho_z_F];
                                
                                const Real u_z_B = Q[1][idx_cell_mom_z_B]/rho_z_B;
                                const Real u_z_F = Q[1][idx_cell_mom_z_F]/rho_z_F;
                                
                                const Real v_z_B = Q[2][idx_cell_mom_z_B]/rho_z_B;
                                const Real v_z_F = Q[2][idx_cell_mom_z_F]/rho_z_F;
                                
                                const Real w_z_B = Q[3][idx_cell_mom_z_B]/rho_z_B;
                                const Real w_z_F = Q[3][idx_cell_mom_z_F]/rho_z_F;
                                
                                const Real epsilon_z_B = Q[4][idx_cell_E_z_B]/rho_z_B - half*(u_z_B*u_z_B + v_z_B*v_z_B + w_z_B*w_z_B);
                                const Real epsilon_z_F = Q[4][idx_cell_E_z_F]/rho_z_F - half*(u_z_F*u_z_F + v_z_F*v_z_F + w_z_F*w_z_F);
                                
                                const Real p_z_B = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &rho_z_B,
                                        &epsilon_z_B,
                                        thermo_properties_ptr);
                                
                                const Real p_z_F = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &rho_z_F,
                                        &epsilon_z_F,
                                        thermo_properties_ptr);
                                
                                // Central derivatives.
                                du_dz = (u_z_F - u_z_B)/(Real(2)*Real(dx[2]));
                                // dv_dz = (v_z_F - v_z_B)/(Real(2)*Real(dx[2]));
                                dw_dz = (w_z_F - w_z_B)/(Real(2)*Real(dx[2]));
                                dp_dz = (p_z_F - p_z_B)/(Real(2)*Real(dx[2]));
                            }
                            
                            const Real c_x_R = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getSoundSpeed(
                                    &rho_x_R,
                                    &p_x_R,
                                    thermo_properties_ptr);
                            
                            const Real lambda_5 = u_x_R + c_x_R;
                            
                            // Compute vector Lambda^(-1) * L.
                            
                            Real Lambda_inv_L[5];
                            
                            const Real& p_t         = d_bdry_face_nonreflecting_outflow_p_t[face_loc];
                            const Real& sigma       = d_bdry_face_nonreflecting_outflow_sigma[face_loc];
                            const Real& beta        = d_bdry_face_nonreflecting_outflow_beta[face_loc];
                            const Real& length_char = d_bdry_face_nonreflecting_outflow_length_char[face_loc];
                            
                            const Real T_5 = v_x_R*(dp_dy + rho_x_R*c_x_R*du_dy) + rho_x_R*c_x_R*c_x_R*dv_dy +
                            w_x_R*(dp_dz + rho_x_R*c_x_R*du_dz) + rho_x_R*c_x_R*c_x_R*dw_dz;
                            
                            const Real M_sq = (u_x_R*u_x_R + v_x_R*v_x_R + w_x_R*w_x_R)/(c_x_R*c_x_R);
                            const Real K = sigma*c_x_R*(Real(1) - M_sq)/length_char;
                            
                            Lambda_inv_L[0] = dp_dx - rho_x_R*c_x_R*du_dx;
                            Lambda_inv_L[1] = c_x_R*c_x_R*drho_dx - dp_dx;
                            Lambda_inv_L[2] = dv_dx;
                            Lambda_inv_L[3] = dw_dx;
                            Lambda_inv_L[4] = (Real(1)/lambda_5)*(K*(p_x_R - p_t) - (Real(1) - beta)*T_5);
                            
                            // Compute dV_dx.
                            
                            const Real c_sq_inv  = Real(1)/(c_x_R*c_x_R);
                            const Real rho_c_inv = Real(1)/(rho_x_R*c_x_R);
                            
                            Real dV_dx[5];
                            
                            dV_dx[0] = half*c_sq_inv*(Lambda_inv_L[0] + Lambda_inv_L[4]) + c_sq_inv*Lambda_inv_L[1];
                            dV_dx[1] = half*rho_c_inv*(-Lambda_inv_L[0] + Lambda_inv_L[4]);
                            dV_dx[2] = Lambda_inv_L[2];
                            dV_dx[3] = Lambda_inv_L[3];
                            dV_dx[4] = half*(Lambda_inv_L[0] + Lambda_inv_L[4]);
                            
                            Real V_ghost[5*num_ghosts_to_fill];
                            
                            for (int i = num_ghosts_to_fill - 1; i >= 0; i--)
                            {
                                const int idx_cell_rho = (i + fill_box_lo_idx[0] + num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0]+
                                    (k + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom = (i + fill_box_lo_idx[0] + num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0]+
                                    (k + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E = (i + fill_box_lo_idx[0] + num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0]+
                                    (k + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                if (i == num_ghosts_to_fill - 1)
                                {
                                    V_ghost[i*5 + 0] = rho_x_RR - Real(2)*Real(dx[0])*dV_dx[0];
                                    V_ghost[i*5 + 1] = u_x_RR   - Real(2)*Real(dx[0])*dV_dx[1];
                                    V_ghost[i*5 + 2] = v_x_RR   - Real(2)*Real(dx[0])*dV_dx[2];
                                    V_ghost[i*5 + 3] = w_x_RR   - Real(2)*Real(dx[0])*dV_dx[3];
                                    V_ghost[i*5 + 4] = p_x_RR   - Real(2)*Real(dx[0])*dV_dx[4];
                                }
                                else if (i == num_ghosts_to_fill - 2)
                                {
                                    V_ghost[i*5 + 0] = -Real(2)*rho_x_RR - Real(3)*rho_x_R +
                                        Real(6)*V_ghost[(i + 1)*5 + 0] + Real(6)*Real(dx[0])*dV_dx[0];
                                    
                                    V_ghost[i*5 + 1] = -Real(2)*u_x_RR - Real(3)*u_x_R +
                                        Real(6)*V_ghost[(i + 1)*5 + 1] + Real(6)*Real(dx[0])*dV_dx[1];
                                    
                                    V_ghost[i*5 + 2] = -Real(2)*v_x_RR - Real(3)*v_x_R +
                                        Real(6)*V_ghost[(i + 1)*5 + 2] + Real(6)*Real(dx[0])*dV_dx[2];
                                    
                                    V_ghost[i*5 + 3] = -Real(2)*w_x_RR - Real(3)*w_x_R +
                                        Real(6)*V_ghost[(i + 1)*5 + 3] + Real(6)*Real(dx[0])*dV_dx[3];
                                    
                                    V_ghost[i*5 + 4] = -Real(2)*p_x_RR - Real(3)*p_x_R +
                                        Real(6)*V_ghost[(i + 1)*5 + 4] + Real(6)*Real(dx[0])*dV_dx[4];
                                }
                                else if (i == num_ghosts_to_fill - 3)
                                {
                                    V_ghost[i*5 + 0] = Real(3)*rho_x_RR + Real(10)*rho_x_R -
                                        Real(18)*V_ghost[(i + 2)*5 + 0] + Real(6)*V_ghost[(i + 1)*5 + 0] -
                                        Real(12)*Real(dx[0])*dV_dx[0];
                                    
                                    V_ghost[i*5 + 1] = Real(3)*u_x_RR + Real(10)*u_x_R -
                                        Real(18)*V_ghost[(i + 2)*5 + 1] + Real(6)*V_ghost[(i + 1)*5 + 1] -
                                        Real(12)*Real(dx[0])*dV_dx[1];
                                    
                                    V_ghost[i*5 + 2] = Real(3)*v_x_RR + Real(10)*v_x_R -
                                        Real(18)*V_ghost[(i + 2)*5 + 2] + Real(6)*V_ghost[(i + 1)*5 + 2] -
                                        Real(12)*Real(dx[0])*dV_dx[2];
                                    
                                    V_ghost[i*5 + 3] = Real(3)*w_x_RR + Real(10)*w_x_R -
                                        Real(18)*V_ghost[(i + 2)*5 + 3] + Real(6)*V_ghost[(i + 1)*5 + 3] -
                                        Real(12)*Real(dx[0])*dV_dx[3];
                                    
                                    V_ghost[i*5 + 4] = Real(3)*p_x_RR + Real(10)*p_x_R -
                                        Real(18)*V_ghost[(i + 2)*5 + 4] + Real(6)*V_ghost[(i + 1)*5 + 4] -
                                        Real(12)*Real(dx[0])*dV_dx[4];
                                }
                                else if (i == num_ghosts_to_fill - 4)
                                {
                                    V_ghost[i*5 + 0] = -Real(4)*rho_x_RR - Real(65)/Real(3)*rho_x_R +
                                        Real(40)*V_ghost[(i + 3)*5 + 0] - Real(20)*V_ghost[(i + 2)*5 + 0] +
                                        Real(20)/Real(3)*V_ghost[(i + 1)*5 + 0] + Real(20)*Real(dx[0])*dV_dx[0];
                                    
                                    V_ghost[i*5 + 1] = -Real(4)*u_x_RR - Real(65)/Real(3)*u_x_R +
                                        Real(40)*V_ghost[(i + 3)*5 + 1] - Real(20)*V_ghost[(i + 2)*5 + 1] +
                                        Real(20)/Real(3)*V_ghost[(i + 1)*5 + 1] + Real(20)*Real(dx[0])*dV_dx[1];
                                    
                                    V_ghost[i*5 + 2] = -Real(4)*v_x_RR - Real(65)/Real(3)*v_x_R +
                                        Real(40)*V_ghost[(i + 3)*5 + 2] - Real(20)*V_ghost[(i + 2)*5 + 2] +
                                        Real(20)/Real(3)*V_ghost[(i + 1)*5 + 2] + Real(20)*Real(dx[0])*dV_dx[2];
                                    
                                    V_ghost[i*5 + 3] = -Real(4)*w_x_RR - Real(65)/Real(3)*w_x_R +
                                        Real(40)*V_ghost[(i + 3)*5 + 3] - Real(20)*V_ghost[(i + 2)*5 + 3] +
                                        Real(20)/Real(3)*V_ghost[(i + 1)*5 + 3] + Real(20)*Real(dx[0])*dV_dx[3];
                                    
                                    V_ghost[i*5 + 4] = -Real(4)*p_x_RR - Real(65)/Real(3)*p_x_R +
                                        Real(40)*V_ghost[(i + 3)*5 + 4] - Real(20)*V_ghost[(i + 2)*5 + 4] +
                                        Real(20)/Real(3)*V_ghost[(i + 1)*5 + 4] + Real(20)*Real(dx[0])*dV_dx[4];
                                }
                                else if (i == num_ghosts_to_fill - 5)
                                {
                                    V_ghost[i*5 + 0] = Real(5)*rho_x_RR + Real(77)/Real(2)*rho_x_R -
                                        Real(75)*V_ghost[(i + 4)*5 + 0] + Real(50)*V_ghost[(i + 3)*5 + 0] -
                                        Real(25)*V_ghost[(i + 2)*5 + 0] + Real(15)/Real(2)*V_ghost[(i + 1)*5 + 0] -
                                        Real(30)*Real(dx[0])*dV_dx[0];
                                    
                                    V_ghost[i*5 + 1] = Real(5)*u_x_RR + Real(77)/Real(2)*u_x_R -
                                        Real(75)*V_ghost[(i + 4)*5 + 1] + Real(50)*V_ghost[(i + 3)*5 + 1] -
                                        Real(25)*V_ghost[(i + 2)*5 + 1] + Real(15)/Real(2)*V_ghost[(i + 1)*5 + 1] -
                                        Real(30)*Real(dx[0])*dV_dx[1];
                                    
                                    V_ghost[i*5 + 2] = Real(5)*v_x_RR + Real(77)/Real(2)*v_x_R -
                                        Real(75)*V_ghost[(i + 4)*5 + 2] + Real(50)*V_ghost[(i + 3)*5 + 2] -
                                        Real(25)*V_ghost[(i + 2)*5 + 2] + Real(15)/Real(2)*V_ghost[(i + 1)*5 + 2] -
                                        Real(30)*Real(dx[0])*dV_dx[2];
                                    
                                    V_ghost[i*5 + 3] = Real(5)*w_x_RR + Real(77)/Real(2)*w_x_R -
                                        Real(75)*V_ghost[(i + 4)*5 + 3] + Real(50)*V_ghost[(i + 3)*5 + 3] -
                                        Real(25)*V_ghost[(i + 2)*5 + 3] + Real(15)/Real(2)*V_ghost[(i + 1)*5 + 3] -
                                        Real(30)*Real(dx[0])*dV_dx[3];
                                    
                                    V_ghost[i*5 + 4] = Real(5)*p_x_RR + Real(77)/Real(2)*p_x_R -
                                        Real(75)*V_ghost[(i + 4)*5 + 4] + Real(50)*V_ghost[(i + 3)*5 + 4] -
                                        Real(25)*V_ghost[(i + 2)*5 + 4] + Real(15)/Real(2)*V_ghost[(i + 1)*5 + 4] -
                                        Real(30)*Real(dx[0])*dV_dx[4];
                                }
                                else if (i == num_ghosts_to_fill - 6)
                                {
                                    V_ghost[i*5 + 0] = -Real(6)*rho_x_RR - Real(609)/Real(10)*rho_x_R +
                                        Real(126)*V_ghost[(i + 5)*5 + 0] - Real(105)*V_ghost[(i + 4)*5 + 0] +
                                        Real(70)*V_ghost[(i + 3)*5 + 0] - Real(63)/Real(2)*V_ghost[(i + 2)*5 + 0] +
                                        Real(42)/Real(5)*V_ghost[(i + 1)*5 + 0] + Real(42)*Real(dx[0])*dV_dx[0];
                                    
                                    V_ghost[i*5 + 1] = -Real(6)*u_x_RR - Real(609)/Real(10)*u_x_R +
                                        Real(126)*V_ghost[(i + 5)*5 + 1] - Real(105)*V_ghost[(i + 4)*5 + 1] +
                                        Real(70)*V_ghost[(i + 3)*5 + 1] - Real(63)/Real(2)*V_ghost[(i + 2)*5 + 1] +
                                        Real(42)/Real(5)*V_ghost[(i + 1)*5 + 1] + Real(42)*Real(dx[0])*dV_dx[1];
                                    
                                    V_ghost[i*5 + 2] = -Real(6)*v_x_RR - Real(609)/Real(10)*v_x_R +
                                        Real(126)*V_ghost[(i + 5)*5 + 2] - Real(105)*V_ghost[(i + 4)*5 + 2] +
                                        Real(70)*V_ghost[(i + 3)*5 + 2] - Real(63)/Real(2)*V_ghost[(i + 2)*5 + 2] +
                                        Real(42)/Real(5)*V_ghost[(i + 1)*5 + 2] + Real(42)*Real(dx[0])*dV_dx[2];
                                    
                                    V_ghost[i*5 + 3] = -Real(6)*w_x_RR - Real(609)/Real(10)*w_x_R +
                                        Real(126)*V_ghost[(i + 5)*5 + 3] - Real(105)*V_ghost[(i + 4)*5 + 3] +
                                        Real(70)*V_ghost[(i + 3)*5 + 3] - Real(63)/Real(2)*V_ghost[(i + 2)*5 + 3] +
                                        Real(42)/Real(5)*V_ghost[(i + 1)*5 + 3] + Real(42)*Real(dx[0])*dV_dx[3];
                                    
                                    V_ghost[i*5 + 4] = -Real(6)*p_x_RR - Real(609)/Real(10)*p_x_R +
                                        Real(126)*V_ghost[(i + 5)*5 + 4] - Real(105)*V_ghost[(i + 4)*5 + 4] +
                                        Real(70)*V_ghost[(i + 3)*5 + 4] - Real(63)/Real(2)*V_ghost[(i + 2)*5 + 4] +
                                        Real(42)/Real(5)*V_ghost[(i + 1)*5 + 4] + Real(42)*Real(dx[0])*dV_dx[4];
                                }
                                
                                Q[0][idx_cell_rho] = V_ghost[i*5 + 0];
                                Q[1][idx_cell_mom] = V_ghost[i*5 + 0]*V_ghost[i*5 + 1];
                                Q[2][idx_cell_mom] = V_ghost[i*5 + 0]*V_ghost[i*5 + 2];
                                Q[3][idx_cell_mom] = V_ghost[i*5 + 0]*V_ghost[i*5 + 3];
                                
                                const Real epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getInternalEnergy(
                                        &V_ghost[i*5 + 0],
                                        &V_ghost[i*5 + 4],
                                        thermo_properties_ptr);
                                
                                const Real E = V_ghost[i*5 + 0]*epsilon +
                                    half*(Q[1][idx_cell_mom]*Q[1][idx_cell_mom] + Q[2][idx_cell_mom]*Q[2][idx_cell_mom] +
                                       Q[3][idx_cell_mom]*Q[3][idx_cell_mom])/V_ghost[i*5 + 0];
                                
                                Q[4][idx_cell_E] = E;
                            }
                        }
                    }
                }
                else if (face_loc == BDRY_LOC::XHI)
                {
                    const int num_ghosts_to_fill = fill_box_hi_idx[0] - fill_box_lo_idx[0] + 1;
                    TBOX_ASSERT(fill_box_lo_idx[0] == interior_box_hi_idx[0] + 1);
                    if (num_ghosts_to_fill > 6)
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelBoundaryUtilitiesSingleSpecies::fill3dFaceBoundaryData()\n"
                            << "Non-reflecting outflow BC doesn't support more than six ghost cells yet!");
                    }
                    
                    for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                    {
                        for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                        {
                            // Get the grid spacing.
                            const double* const dx = patch_geom->getDx();
                            
                            const int idx_cell_rho_x_L = (interior_box_hi_idx[0] + num_subghosts_conservative_var[0][0]) +
                                (j + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                (k + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                    subghostcell_dims_conservative_var[0][1];
                            
                            const int idx_cell_rho_x_LL = (interior_box_hi_idx[0] - 1 + num_subghosts_conservative_var[0][0]) +
                                (j + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                (k + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                    subghostcell_dims_conservative_var[0][1];
                            
                            const int idx_cell_rho_x_LLL = (interior_box_hi_idx[0] - 2 + num_subghosts_conservative_var[0][0]) +
                                (j + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                (k + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                    subghostcell_dims_conservative_var[0][1];
                            
                            const int idx_cell_mom_x_L = (interior_box_hi_idx[0] + num_subghosts_conservative_var[1][0]) +
                                (j + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                (k + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                    subghostcell_dims_conservative_var[1][1];
                            
                            const int idx_cell_mom_x_LL = (interior_box_hi_idx[0] - 1 + num_subghosts_conservative_var[1][0]) +
                                (j + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                (k + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                    subghostcell_dims_conservative_var[1][1];
                            
                            const int idx_cell_mom_x_LLL = (interior_box_hi_idx[0] - 2 + num_subghosts_conservative_var[1][0]) +
                                (j + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                (k + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                    subghostcell_dims_conservative_var[1][1];
                            
                            const int idx_cell_E_x_L = (interior_box_hi_idx[0] + num_subghosts_conservative_var[2][0]) +
                                (j + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                (k + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                    subghostcell_dims_conservative_var[2][1];
                            
                            const int idx_cell_E_x_LL = (interior_box_hi_idx[0] - 1 + num_subghosts_conservative_var[2][0]) +
                                (j + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                (k + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                    subghostcell_dims_conservative_var[2][1];
                            
                            const int idx_cell_E_x_LLL = (interior_box_hi_idx[0] - 2 + num_subghosts_conservative_var[2][0]) +
                                (j + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                (k + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                    subghostcell_dims_conservative_var[2][1];
                            
                            const Real& rho_x_L   = Q[0][idx_cell_rho_x_L];
                            const Real& rho_x_LL  = Q[0][idx_cell_rho_x_LL];
                            const Real& rho_x_LLL = Q[0][idx_cell_rho_x_LLL];
                            
                            const Real u_x_L   = Q[1][idx_cell_mom_x_L]/rho_x_L;
                            const Real u_x_LL  = Q[1][idx_cell_mom_x_LL]/rho_x_LL;
                            const Real u_x_LLL = Q[1][idx_cell_mom_x_LLL]/rho_x_LLL;
                            
                            const Real v_x_L   = Q[2][idx_cell_mom_x_L]/rho_x_L;
                            const Real v_x_LL  = Q[2][idx_cell_mom_x_LL]/rho_x_LL;
                            const Real v_x_LLL = Q[2][idx_cell_mom_x_LLL]/rho_x_LLL;
                            
                            const Real w_x_L   = Q[3][idx_cell_mom_x_L]/rho_x_L;
                            const Real w_x_LL  = Q[3][idx_cell_mom_x_LL]/rho_x_LL;
                            const Real w_x_LLL = Q[3][idx_cell_mom_x_LLL]/rho_x_LLL;
                            
                            const Real half = Real(1)/Real(2);
                            const Real epsilon_x_L   = Q[4][idx_cell_E_x_L]/rho_x_L - half*(u_x_L*u_x_L + v_x_L*v_x_L + w_x_L*w_x_L);
                            const Real epsilon_x_LL  = Q[4][idx_cell_E_x_LL]/rho_x_LL - half*(u_x_LL*u_x_LL + v_x_LL*v_x_LL + w_x_LL*w_x_LL);
                            const Real epsilon_x_LLL = Q[4][idx_cell_E_x_LLL]/rho_x_LLL - half*(u_x_LLL*u_x_LLL + v_x_LLL*v_x_LLL + w_x_LLL*w_x_LLL);
                            
                            const Real p_x_L = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getPressure(
                                    &rho_x_L,
                                    &epsilon_x_L,
                                    thermo_properties_ptr);
                            
                            const Real p_x_LL = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getPressure(
                                    &rho_x_LL,
                                    &epsilon_x_LL,
                                    thermo_properties_ptr);
                            
                            const Real p_x_LLL = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getPressure(
                                    &rho_x_LLL,
                                    &epsilon_x_LLL,
                                    thermo_properties_ptr);
                            
                            /*
                             * Compute derivatives at x-direction.
                             */
                            
                            const Real drho_dx = (rho_x_LLL - Real(4)*rho_x_LL + Real(3)*rho_x_L)/(Real(2)*Real(dx[0]));
                            const Real du_dx   = (u_x_LLL - Real(4)*u_x_LL + Real(3)*u_x_L)/(Real(2)*Real(dx[0]));
                            const Real dv_dx   = (v_x_LLL - Real(4)*v_x_LL + Real(3)*v_x_L)/(Real(2)*Real(dx[0]));
                            const Real dw_dx   = (w_x_LLL - Real(4)*w_x_LL + Real(3)*w_x_L)/(Real(2)*Real(dx[0]));
                            const Real dp_dx   = (p_x_LLL - Real(4)*p_x_LL + Real(3)*p_x_L)/(Real(2)*Real(dx[0]));
                            
                            /*
                             * Compute derivatives in y-direction.
                             */
                            
                            Real du_dy = Real(0);
                            Real dv_dy = Real(0);
                            // Real dw_dy = Real(0);
                            Real dp_dy = Real(0);
                            
                            if ((j + num_subghosts_conservative_var[0][1] == 0) ||
                                (j + num_subghosts_conservative_var[1][1] == 0) ||
                                (j + num_subghosts_conservative_var[2][1] == 0))
                            // if (((patch_geom->getTouchesRegularBoundary(1, 0)) && (j == interior_box_lo_idx[1])) ||
                            //     ((j + num_subghosts_conservative_var[0][1] == 0) ||
                            //      (j + num_subghosts_conservative_var[1][1] == 0) ||
                            //      (j + num_subghosts_conservative_var[2][1] == 0)))
                            {
                                // Patch is touching periodic boundary.
                                
                                const int idx_cell_rho_y_T = (interior_box_hi_idx[0] + num_subghosts_conservative_var[0][0]) +
                                    (j + 1 + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (k + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom_y_T = (interior_box_hi_idx[0] + num_subghosts_conservative_var[1][0]) +
                                    (j + 1 + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                    (k + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E_y_T = (interior_box_hi_idx[0] + num_subghosts_conservative_var[2][0]) +
                                    (j + 1 + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                    (k + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                const Real& rho_y_T = Q[0][idx_cell_rho_y_T];
                                const Real u_y_T = Q[1][idx_cell_mom_y_T]/rho_y_T;
                                const Real v_y_T = Q[2][idx_cell_mom_y_T]/rho_y_T;
                                const Real w_y_T = Q[3][idx_cell_mom_y_T]/rho_y_T;
                                const Real epsilon_y_T = Q[4][idx_cell_E_y_T]/rho_y_T - half*(u_y_T*u_y_T + v_y_T*v_y_T + w_y_T*w_y_T);
                                
                                const Real p_y_T = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &rho_y_T,
                                        &epsilon_y_T,
                                        thermo_properties_ptr);
                                
                                // One-sided derivatives.
                                du_dy = (u_y_T - u_x_L)/Real(dx[1]);
                                dv_dy = (v_y_T - v_x_L)/Real(dx[1]);
                                // dw_dy = (w_y_T - w_x_L)/Real(dx[1]);
                                dp_dy = (p_y_T - p_x_L)/Real(dx[1]);
                            }
                            else if ((j + num_subghosts_conservative_var[0][1] + 1 == subghostcell_dims_conservative_var[0][1]) ||
                                     (j + num_subghosts_conservative_var[1][1] + 1 == subghostcell_dims_conservative_var[1][1]) ||
                                     (j + num_subghosts_conservative_var[2][1] + 1 == subghostcell_dims_conservative_var[2][1]))
                            // else if (((patch_geom->getTouchesRegularBoundary(1, 1)) && (j == interior_box_hi_idx[1])) ||
                            //          ((j + num_subghosts_conservative_var[0][1] + 1 == subghostcell_dims_conservative_var[0][1]) ||
                            //           (j + num_subghosts_conservative_var[1][1] + 1 == subghostcell_dims_conservative_var[1][1]) ||
                            //           (j + num_subghosts_conservative_var[2][1] + 1 == subghostcell_dims_conservative_var[2][1])))
                            {
                                // Patch is touching periodic boundary.
                                
                                const int idx_cell_rho_y_B = (interior_box_hi_idx[0] + num_subghosts_conservative_var[0][0]) +
                                    (j - 1 + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (k + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom_y_B = (interior_box_hi_idx[0] + num_subghosts_conservative_var[1][0]) +
                                    (j - 1 + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                    (k + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E_y_B = (interior_box_hi_idx[0] + num_subghosts_conservative_var[2][0]) +
                                    (j - 1 + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                    (k + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                const Real& rho_y_B = Q[0][idx_cell_rho_y_B];
                                const Real u_y_B = Q[1][idx_cell_mom_y_B]/rho_y_B;
                                const Real v_y_B = Q[2][idx_cell_mom_y_B]/rho_y_B;
                                const Real w_y_B = Q[3][idx_cell_mom_y_B]/rho_y_B;
                                const Real epsilon_y_B = Q[4][idx_cell_E_y_B]/rho_y_B - half*(u_y_B*u_y_B + v_y_B*v_y_B + w_y_B*w_y_B);
                                
                                const Real p_y_B = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &rho_y_B,
                                        &epsilon_y_B,
                                        thermo_properties_ptr);
                                
                                // One-sided derivatives.
                                du_dy = (u_x_L - u_y_B)/Real(dx[1]);
                                dv_dy = (v_x_L - v_y_B)/Real(dx[1]);
                                // dw_dy = (w_x_L - w_y_B)/Real(dx[1]);
                                dp_dy = (p_x_L - p_y_B)/Real(dx[1]);
                            }
                            else
                            {
                                const int idx_cell_rho_y_B = (interior_box_hi_idx[0] + num_subghosts_conservative_var[0][0]) +
                                    (j - 1 + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (k + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_rho_y_T = (interior_box_hi_idx[0] + num_subghosts_conservative_var[0][0]) +
                                    (j + 1 + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (k + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom_y_B = (interior_box_hi_idx[0] + num_subghosts_conservative_var[1][0]) +
                                    (j - 1 + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                    (k + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_mom_y_T = (interior_box_hi_idx[0] + num_subghosts_conservative_var[1][0]) +
                                    (j + 1 + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                    (k + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E_y_B = (interior_box_hi_idx[0] + num_subghosts_conservative_var[2][0]) +
                                    (j - 1 + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                    (k + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                const int idx_cell_E_y_T = (interior_box_hi_idx[0] + num_subghosts_conservative_var[2][0]) +
                                    (j + 1 + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                    (k + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                const Real& rho_y_B = Q[0][idx_cell_rho_y_B];
                                const Real& rho_y_T = Q[0][idx_cell_rho_y_T];
                                
                                const Real u_y_B = Q[1][idx_cell_mom_y_B]/rho_y_B;
                                const Real u_y_T = Q[1][idx_cell_mom_y_T]/rho_y_T;
                                
                                const Real v_y_B = Q[2][idx_cell_mom_y_B]/rho_y_B;
                                const Real v_y_T = Q[2][idx_cell_mom_y_T]/rho_y_T;
                                
                                const Real w_y_B = Q[3][idx_cell_mom_y_B]/rho_y_B;
                                const Real w_y_T = Q[3][idx_cell_mom_y_T]/rho_y_T;
                                
                                const Real epsilon_y_B = Q[4][idx_cell_E_y_B]/rho_y_B - half*(u_y_B*u_y_B + v_y_B*v_y_B + w_y_B*w_y_B);
                                const Real epsilon_y_T = Q[4][idx_cell_E_y_T]/rho_y_T - half*(u_y_T*u_y_T + v_y_T*v_y_T + w_y_T*w_y_T);
                                
                                const Real p_y_B = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &rho_y_B,
                                        &epsilon_y_B,
                                        thermo_properties_ptr);
                                
                                const Real p_y_T = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &rho_y_T,
                                        &epsilon_y_T,
                                        thermo_properties_ptr);
                                
                                // Central derivatives.
                                du_dy = (u_y_T - u_y_B)/(Real(2)*Real(dx[1]));
                                dv_dy = (v_y_T - v_y_B)/(Real(2)*Real(dx[1]));
                                // dw_dy = (w_y_T - w_y_B)/(Real(2)*Real(dx[1]));
                                dp_dy = (p_y_T - p_y_B)/(Real(2)*Real(dx[1]));
                            }
                            
                            /*
                             * Compute derivatives in z-direction.
                             */
                            
                            Real du_dz = Real(0);
                            // Real dv_dz = Real(0);
                            Real dw_dz = Real(0);
                            Real dp_dz = Real(0);
                            
                            if ((k + num_subghosts_conservative_var[0][2] == 0) ||
                                (k + num_subghosts_conservative_var[1][2] == 0) ||
                                (k + num_subghosts_conservative_var[2][2] == 0))
                            // if (((patch_geom->getTouchesRegularBoundary(2, 0)) && (k == interior_box_lo_idx[2])) ||
                            //     ((k + num_subghosts_conservative_var[0][2] == 0) ||
                            //      (k + num_subghosts_conservative_var[1][2] == 0) ||
                            //      (k + num_subghosts_conservative_var[2][2] == 0)))
                            {
                                // Patch is touching periodic boundary.
                                
                                const int idx_cell_rho_z_F = (interior_box_hi_idx[0] + num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (k + 1 + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom_z_F = (interior_box_hi_idx[0] + num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                    (k + 1 + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E_z_F = (interior_box_hi_idx[0] + num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                    (k + 1 + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                const Real& rho_z_F = Q[0][idx_cell_rho_z_F];
                                const Real u_z_F = Q[1][idx_cell_mom_z_F]/rho_z_F;
                                const Real v_z_F = Q[2][idx_cell_mom_z_F]/rho_z_F;
                                const Real w_z_F = Q[3][idx_cell_mom_z_F]/rho_z_F;
                                const Real epsilon_z_F = Q[4][idx_cell_E_z_F]/rho_z_F - half*(u_z_F*u_z_F + v_z_F*v_z_F + w_z_F*w_z_F);
                                
                                const Real p_z_F = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &rho_z_F,
                                        &epsilon_z_F,
                                        thermo_properties_ptr);
                                
                                // One-sided derivatives.
                                du_dz = (u_z_F - u_x_L)/Real(dx[2]);
                                // dv_dz = (v_z_F - v_x_L)/Real(dx[2]);
                                dw_dz = (w_z_F - w_x_L)/Real(dx[2]);
                                dp_dz = (p_z_F - p_x_L)/Real(dx[2]);
                            }
                            else if ((k + num_subghosts_conservative_var[0][2] + 1 == subghostcell_dims_conservative_var[0][2]) ||
                                     (k + num_subghosts_conservative_var[1][2] + 1 == subghostcell_dims_conservative_var[1][2]) ||
                                     (k + num_subghosts_conservative_var[2][2] + 1 == subghostcell_dims_conservative_var[2][2]))
                            // else if (((patch_geom->getTouchesRegularBoundary(2, 1)) && (k == interior_box_hi_idx[2])) ||
                            //          ((k + num_subghosts_conservative_var[0][2] + 1 == subghostcell_dims_conservative_var[0][2]) ||
                            //           (k + num_subghosts_conservative_var[1][2] + 1 == subghostcell_dims_conservative_var[1][2]) ||
                            //           (k + num_subghosts_conservative_var[2][2] + 1 == subghostcell_dims_conservative_var[2][2])))
                            {
                                // Patch is touching periodic boundary.
                                
                                const int idx_cell_rho_z_B = (interior_box_hi_idx[0] + num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (k - 1 + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom_z_B = (interior_box_hi_idx[0] + num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                    (k - 1 + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E_z_B = (interior_box_hi_idx[0] + num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                    (k - 1 + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                const Real& rho_z_B = Q[0][idx_cell_rho_z_B];
                                const Real u_z_B = Q[1][idx_cell_mom_z_B]/rho_z_B;
                                const Real v_z_B = Q[2][idx_cell_mom_z_B]/rho_z_B;
                                const Real w_z_B = Q[3][idx_cell_mom_z_B]/rho_z_B;
                                const Real epsilon_z_B = Q[4][idx_cell_E_z_B]/rho_z_B - half*(u_z_B*u_z_B + v_z_B*v_z_B + w_z_B*w_z_B);
                                
                                const Real p_z_B = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &rho_z_B,
                                        &epsilon_z_B,
                                        thermo_properties_ptr);
                                
                                // One-sided derivatives.
                                du_dz = (u_x_L - u_z_B)/Real(dx[2]);
                                // dv_dz = (v_x_L - v_z_B)/Real(dx[2]);
                                dw_dz = (w_x_L - w_z_B)/Real(dx[2]);
                                dp_dz = (p_x_L - p_z_B)/Real(dx[2]);
                            }
                            else
                            {
                                const int idx_cell_rho_z_B = (interior_box_hi_idx[0] + num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (k - 1 + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_rho_z_F = (interior_box_hi_idx[0] + num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (k + 1 + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom_z_B = (interior_box_hi_idx[0] + num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                    (k - 1 + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_mom_z_F = (interior_box_hi_idx[0] + num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                    (k + 1 + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E_z_B = (interior_box_hi_idx[0] + num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                    (k - 1 + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                const int idx_cell_E_z_F = (interior_box_hi_idx[0] + num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                    (k + 1 + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                const Real& rho_z_B = Q[0][idx_cell_rho_z_B];
                                const Real& rho_z_F = Q[0][idx_cell_rho_z_F];
                                
                                const Real u_z_B = Q[1][idx_cell_mom_z_B]/rho_z_B;
                                const Real u_z_F = Q[1][idx_cell_mom_z_F]/rho_z_F;
                                
                                const Real v_z_B = Q[2][idx_cell_mom_z_B]/rho_z_B;
                                const Real v_z_F = Q[2][idx_cell_mom_z_F]/rho_z_F;
                                
                                const Real w_z_B = Q[3][idx_cell_mom_z_B]/rho_z_B;
                                const Real w_z_F = Q[3][idx_cell_mom_z_F]/rho_z_F;
                                
                                const Real epsilon_z_B = Q[4][idx_cell_E_z_B]/rho_z_B - half*(u_z_B*u_z_B + v_z_B*v_z_B + w_z_B*w_z_B);
                                const Real epsilon_z_F = Q[4][idx_cell_E_z_F]/rho_z_F - half*(u_z_F*u_z_F + v_z_F*v_z_F + w_z_F*w_z_F);
                                
                                const Real p_z_B = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &rho_z_B,
                                        &epsilon_z_B,
                                        thermo_properties_ptr);
                                
                                const Real p_z_F = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &rho_z_F,
                                        &epsilon_z_F,
                                        thermo_properties_ptr);
                                
                                // Central derivatives.
                                du_dz = (u_z_F - u_z_B)/(Real(2)*Real(dx[2]));
                                // dv_dz = (v_z_F - v_z_B)/(Real(2)*Real(dx[2]));
                                dw_dz = (w_z_F - w_z_B)/(Real(2)*Real(dx[2]));
                                dp_dz = (p_z_F - p_z_B)/(Real(2)*Real(dx[2]));
                            }
                            
                            const Real c_x_L = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getSoundSpeed(
                                    &rho_x_L,
                                    &p_x_L,
                                    thermo_properties_ptr);
                            
                            const Real lambda_1 = u_x_L - c_x_L;
                            
                            // Compute vector Lambda^(-1) * L.
                            
                            Real Lambda_inv_L[5];
                            
                            const Real& p_t         = d_bdry_face_nonreflecting_outflow_p_t[face_loc];
                            const Real& sigma       = d_bdry_face_nonreflecting_outflow_sigma[face_loc];
                            const Real& beta        = d_bdry_face_nonreflecting_outflow_beta[face_loc];
                            const Real& length_char = d_bdry_face_nonreflecting_outflow_length_char[face_loc];
                            
                            const Real T_1 = v_x_L*(dp_dy - rho_x_L*c_x_L*du_dy) + rho_x_L*c_x_L*c_x_L*dv_dy +
                            w_x_L*(dp_dz - rho_x_L*c_x_L*du_dz) + rho_x_L*c_x_L*c_x_L*dw_dz;
                            
                            const Real M_sq = (u_x_L*u_x_L + v_x_L*v_x_L + w_x_L*w_x_L)/(c_x_L*c_x_L);
                            const Real K = sigma*c_x_L*(Real(1) - M_sq)/length_char;
                            
                            Lambda_inv_L[0] = (Real(1)/lambda_1)*(K*(p_x_L - p_t) - (Real(1) - beta)*T_1);
                            Lambda_inv_L[1] = c_x_L*c_x_L*drho_dx - dp_dx;
                            Lambda_inv_L[2] = dv_dx;
                            Lambda_inv_L[3] = dw_dx;
                            Lambda_inv_L[4] = dp_dx + rho_x_L*c_x_L*du_dx;
                            
                            // Compute dV_dx.
                            
                            const Real c_sq_inv  = Real(1)/(c_x_L*c_x_L);
                            const Real rho_c_inv = Real(1)/(rho_x_L*c_x_L);
                            
                            Real dV_dx[5];
                            
                            dV_dx[0] = half*c_sq_inv*(Lambda_inv_L[0] + Lambda_inv_L[4]) + c_sq_inv*Lambda_inv_L[1];
                            dV_dx[1] = half*rho_c_inv*(-Lambda_inv_L[0] + Lambda_inv_L[4]);
                            dV_dx[2] = Lambda_inv_L[2];
                            dV_dx[3] = Lambda_inv_L[3];
                            dV_dx[4] = half*(Lambda_inv_L[0] + Lambda_inv_L[4]);
                            
                            Real V_ghost[5*num_ghosts_to_fill];
                            
                            for (int i = 0; i < num_ghosts_to_fill; i++)
                            {
                                const int idx_cell_rho = (i + fill_box_lo_idx[0] + num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0]+
                                    (k + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom = (i + fill_box_lo_idx[0] + num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0]+
                                    (k + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E = (i + fill_box_lo_idx[0] + num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0]+
                                    (k + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                if (i == 0)
                                {
                                    V_ghost[i*5 + 0] = rho_x_LL + Real(2)*Real(dx[0])*dV_dx[0];
                                    V_ghost[i*5 + 1] = u_x_LL   + Real(2)*Real(dx[0])*dV_dx[1];
                                    V_ghost[i*5 + 2] = v_x_LL   + Real(2)*Real(dx[0])*dV_dx[2];
                                    V_ghost[i*5 + 3] = w_x_LL   + Real(2)*Real(dx[0])*dV_dx[3];
                                    V_ghost[i*5 + 4] = p_x_LL   + Real(2)*Real(dx[0])*dV_dx[4];
                                }
                                else if (i == 1)
                                {
                                    V_ghost[i*5 + 0] = -Real(2)*rho_x_LL - Real(3)*rho_x_L +
                                        Real(6)*V_ghost[(i - 1)*5 + 0] - Real(6)*Real(dx[0])*dV_dx[0];
                                    
                                    V_ghost[i*5 + 1] = -Real(2)*u_x_LL - Real(3)*u_x_L +
                                        Real(6)*V_ghost[(i - 1)*5 + 1] - Real(6)*Real(dx[0])*dV_dx[1];
                                    
                                    V_ghost[i*5 + 2] = -Real(2)*v_x_LL - Real(3)*v_x_L +
                                        Real(6)*V_ghost[(i - 1)*5 + 2] - Real(6)*Real(dx[0])*dV_dx[2];
                                    
                                    V_ghost[i*5 + 3] = -Real(2)*w_x_LL - Real(3)*w_x_L +
                                        Real(6)*V_ghost[(i - 1)*5 + 3] - Real(6)*Real(dx[0])*dV_dx[3];
                                    
                                    V_ghost[i*5 + 4] = -Real(2)*p_x_LL - Real(3)*p_x_L +
                                        Real(6)*V_ghost[(i - 1)*5 + 4] - Real(6)*Real(dx[0])*dV_dx[4];
                                }
                                else if (i == 2)
                                {
                                    V_ghost[i*5 + 0] = Real(3)*rho_x_LL + Real(10)*rho_x_L -
                                        Real(18)*V_ghost[(i - 2)*5 + 0] + Real(6)*V_ghost[(i - 1)*5 + 0] +
                                        Real(12)*Real(dx[0])*dV_dx[0];
                                    
                                    V_ghost[i*5 + 1] = Real(3)*u_x_LL + Real(10)*u_x_L -
                                        Real(18)*V_ghost[(i - 2)*5 + 1] + Real(6)*V_ghost[(i - 1)*5 + 1] +
                                        Real(12)*Real(dx[0])*dV_dx[1];
                                    
                                    V_ghost[i*5 + 2] = Real(3)*v_x_LL + Real(10)*v_x_L -
                                        Real(18)*V_ghost[(i - 2)*5 + 2] + Real(6)*V_ghost[(i - 1)*5 + 2] +
                                        Real(12)*Real(dx[0])*dV_dx[2];
                                    
                                    V_ghost[i*5 + 3] = Real(3)*w_x_LL + Real(10)*w_x_L -
                                        Real(18)*V_ghost[(i - 2)*5 + 3] + Real(6)*V_ghost[(i - 1)*5 + 3] +
                                        Real(12)*Real(dx[0])*dV_dx[3];
                                    
                                    V_ghost[i*5 + 4] = Real(3)*p_x_LL + Real(10)*p_x_L -
                                        Real(18)*V_ghost[(i - 2)*5 + 4] + Real(6)*V_ghost[(i - 1)*5 + 4] +
                                        Real(12)*Real(dx[0])*dV_dx[4];
                                }
                                else if (i == 3)
                                {
                                    V_ghost[i*5 + 0] = -Real(4)*rho_x_LL - Real(65)/Real(3)*rho_x_L +
                                        Real(40)*V_ghost[(i - 3)*5 + 0] - Real(20)*V_ghost[(i - 2)*5 + 0] +
                                        Real(20)/Real(3)*V_ghost[(i - 1)*5 + 0] - Real(20)*Real(dx[0])*dV_dx[0];
                                    
                                    V_ghost[i*5 + 1] = -Real(4)*u_x_LL - Real(65)/Real(3)*u_x_L +
                                        Real(40)*V_ghost[(i - 3)*5 + 1] - Real(20)*V_ghost[(i - 2)*5 + 1] +
                                        Real(20)/Real(3)*V_ghost[(i - 1)*5 + 1] - Real(20)*Real(dx[0])*dV_dx[1];
                                    
                                    V_ghost[i*5 + 2] = -Real(4)*v_x_LL - Real(65)/Real(3)*v_x_L +
                                        Real(40)*V_ghost[(i - 3)*5 + 2] - Real(20)*V_ghost[(i - 2)*5 + 2] +
                                        Real(20)/Real(3)*V_ghost[(i - 1)*5 + 2] - Real(20)*Real(dx[0])*dV_dx[2];
                                    
                                    V_ghost[i*5 + 3] = -Real(4)*w_x_LL - Real(65)/Real(3)*w_x_L +
                                        Real(40)*V_ghost[(i - 3)*5 + 3] - Real(20)*V_ghost[(i - 2)*5 + 3] +
                                        Real(20)/Real(3)*V_ghost[(i - 1)*5 + 3] - Real(20)*Real(dx[0])*dV_dx[3];
                                    
                                    V_ghost[i*5 + 4] = -Real(4)*p_x_LL - Real(65)/Real(3)*p_x_L +
                                        Real(40)*V_ghost[(i - 3)*5 + 4] - Real(20)*V_ghost[(i - 2)*5 + 4] +
                                        Real(20)/Real(3)*V_ghost[(i - 1)*5 + 4] - Real(20)*Real(dx[0])*dV_dx[4];
                                }
                                else if (i == 4)
                                {
                                    V_ghost[i*5 + 0] = Real(5)*rho_x_LL + Real(77)/Real(2)*rho_x_L -
                                        Real(75)*V_ghost[(i - 4)*5 + 0] + Real(50)*V_ghost[(i - 3)*5 + 0] -
                                        Real(25)*V_ghost[(i - 2)*5 + 0] + Real(15)/Real(2)*V_ghost[(i - 1)*5 + 0] +
                                        Real(30)*Real(dx[0])*dV_dx[0];
                                    
                                    V_ghost[i*5 + 1] = Real(5)*u_x_LL + Real(77)/Real(2)*u_x_L -
                                        Real(75)*V_ghost[(i - 4)*5 + 1] + Real(50)*V_ghost[(i - 3)*5 + 1] -
                                        Real(25)*V_ghost[(i - 2)*5 + 1] + Real(15)/Real(2)*V_ghost[(i - 1)*5 + 1] +
                                        Real(30)*Real(dx[0])*dV_dx[1];
                                    
                                    V_ghost[i*5 + 2] = Real(5)*v_x_LL + Real(77)/Real(2)*v_x_L -
                                        Real(75)*V_ghost[(i - 4)*5 + 2] + Real(50)*V_ghost[(i - 3)*5 + 2] -
                                        Real(25)*V_ghost[(i - 2)*5 + 2] + Real(15)/Real(2)*V_ghost[(i - 1)*5 + 2] +
                                        Real(30)*Real(dx[0])*dV_dx[2];
                                    
                                    V_ghost[i*5 + 3] = Real(5)*w_x_LL + Real(77)/Real(2)*w_x_L -
                                        Real(75)*V_ghost[(i - 4)*5 + 3] + Real(50)*V_ghost[(i - 3)*5 + 3] -
                                        Real(25)*V_ghost[(i - 2)*5 + 3] + Real(15)/Real(2)*V_ghost[(i - 1)*5 + 3] +
                                        Real(30)*Real(dx[0])*dV_dx[3];
                                    
                                    V_ghost[i*5 + 4] = Real(5)*p_x_LL + Real(77)/Real(2)*p_x_L -
                                        Real(75)*V_ghost[(i - 4)*5 + 4] + Real(50)*V_ghost[(i - 3)*5 + 4] -
                                        Real(25)*V_ghost[(i - 2)*5 + 4] + Real(15)/Real(2)*V_ghost[(i - 1)*5 + 4] +
                                        Real(30)*Real(dx[0])*dV_dx[4];
                                }
                                else if (i == 5)
                                {
                                    V_ghost[i*5 + 0] = -Real(6)*rho_x_LL - Real(609)/Real(10)*rho_x_L +
                                        Real(126)*V_ghost[(i - 5)*5 + 0] - Real(105)*V_ghost[(i - 4)*5 + 0] +
                                        Real(70)*V_ghost[(i - 3)*5 + 0] - Real(63)/Real(2)*V_ghost[(i - 2)*5 + 0] +
                                        Real(42)/Real(5)*V_ghost[(i - 1)*5 + 0] - Real(42)*Real(dx[0])*dV_dx[0];
                                    
                                    V_ghost[i*5 + 1] = -Real(6)*u_x_LL - Real(609)/Real(10)*u_x_L +
                                        Real(126)*V_ghost[(i - 5)*5 + 1] - Real(105)*V_ghost[(i - 4)*5 + 1] +
                                        Real(70)*V_ghost[(i - 3)*5 + 1] - Real(63)/Real(2)*V_ghost[(i - 2)*5 + 1] +
                                        Real(42)/Real(5)*V_ghost[(i - 1)*5 + 1] - Real(42)*Real(dx[0])*dV_dx[1];
                                    
                                    V_ghost[i*5 + 2] = -Real(6)*v_x_LL - Real(609)/Real(10)*v_x_L +
                                        Real(126)*V_ghost[(i - 5)*5 + 2] - Real(105)*V_ghost[(i - 4)*5 + 2] +
                                        Real(70)*V_ghost[(i - 3)*5 + 2] - Real(63)/Real(2)*V_ghost[(i - 2)*5 + 2] +
                                        Real(42)/Real(5)*V_ghost[(i - 1)*5 + 2] - Real(42)*Real(dx[0])*dV_dx[2];
                                    
                                    V_ghost[i*5 + 3] = -Real(6)*w_x_LL - Real(609)/Real(10)*w_x_L +
                                        Real(126)*V_ghost[(i - 5)*5 + 3] - Real(105)*V_ghost[(i - 4)*5 + 3] +
                                        Real(70)*V_ghost[(i - 3)*5 + 3] - Real(63)/Real(2)*V_ghost[(i - 2)*5 + 3] +
                                        Real(42)/Real(5)*V_ghost[(i - 1)*5 + 3] - Real(42)*Real(dx[0])*dV_dx[3];
                                    
                                    V_ghost[i*5 + 4] = -Real(6)*p_x_LL - Real(609)/Real(10)*p_x_L +
                                        Real(126)*V_ghost[(i - 5)*5 + 4] - Real(105)*V_ghost[(i - 4)*5 + 4] +
                                        Real(70)*V_ghost[(i - 3)*5 + 4] - Real(63)/Real(2)*V_ghost[(i - 2)*5 + 4] +
                                        Real(42)/Real(5)*V_ghost[(i - 1)*5 + 4] - Real(42)*Real(dx[0])*dV_dx[4];
                                }
                                
                                Q[0][idx_cell_rho] = V_ghost[i*5 + 0];
                                Q[1][idx_cell_mom] = V_ghost[i*5 + 0]*V_ghost[i*5 + 1];
                                Q[2][idx_cell_mom] = V_ghost[i*5 + 0]*V_ghost[i*5 + 2];
                                Q[3][idx_cell_mom] = V_ghost[i*5 + 0]*V_ghost[i*5 + 3];
                                
                                const Real epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getInternalEnergy(
                                        &V_ghost[i*5 + 0],
                                        &V_ghost[i*5 + 4],
                                        thermo_properties_ptr);
                                
                                const Real E = V_ghost[i*5 + 0]*epsilon +
                                    half*(Q[1][idx_cell_mom]*Q[1][idx_cell_mom] + Q[2][idx_cell_mom]*Q[2][idx_cell_mom] +
                                       Q[3][idx_cell_mom]*Q[3][idx_cell_mom])/V_ghost[i*5 + 0];
                                
                                Q[4][idx_cell_E] = E;
                            }
                        }
                    }
                }
                else if (face_loc == BDRY_LOC::YLO)
                {
                    const int num_ghosts_to_fill = fill_box_hi_idx[1] - fill_box_lo_idx[1] + 1;
                    TBOX_ASSERT(fill_box_hi_idx[1] == interior_box_lo_idx[1] - 1);
                    if (num_ghosts_to_fill > 6)
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelBoundaryUtilitiesSingleSpecies::fill3dFaceBoundaryData()\n"
                            << "Non-reflecting outflow BC doesn't support more than six ghost cells yet!");
                    }
                    
                    for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                    {
                        for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                        {
                            // Get the grid spacing.
                            const double* const dx = patch_geom->getDx();
                            
                            const int idx_cell_rho_y_T = (i + num_subghosts_conservative_var[0][0]) +
                                (interior_box_lo_idx[1] + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                (k + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                    subghostcell_dims_conservative_var[0][1];
                            
                            const int idx_cell_rho_y_TT = (i + num_subghosts_conservative_var[0][0]) +
                                (interior_box_lo_idx[1] + 1 + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                (k + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                    subghostcell_dims_conservative_var[0][1];
                            
                            const int idx_cell_rho_y_TTT = (i + num_subghosts_conservative_var[0][0]) +
                                (interior_box_lo_idx[1] + 2 + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                (k + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                    subghostcell_dims_conservative_var[0][1];
                            
                            const int idx_cell_mom_y_T = (i + num_subghosts_conservative_var[1][0]) +
                                ( interior_box_lo_idx[1] + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                (k + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                    subghostcell_dims_conservative_var[1][1];
                            
                            const int idx_cell_mom_y_TT = (i + num_subghosts_conservative_var[1][0]) +
                                (interior_box_lo_idx[1] + 1 + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                (k + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                    subghostcell_dims_conservative_var[1][1];
                            
                            const int idx_cell_mom_y_TTT = (i + num_subghosts_conservative_var[1][0]) +
                                (interior_box_lo_idx[1] + 2 + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                (k + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                    subghostcell_dims_conservative_var[1][1];
                            
                            const int idx_cell_E_y_T = (i + num_subghosts_conservative_var[2][0]) +
                                (interior_box_lo_idx[1] + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                (k + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                    subghostcell_dims_conservative_var[2][1];
                            
                            const int idx_cell_E_y_TT = (i + num_subghosts_conservative_var[2][0]) +
                                (interior_box_lo_idx[1] + 1 + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                (k + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                    subghostcell_dims_conservative_var[2][1];
                            
                            const int idx_cell_E_y_TTT = (i + num_subghosts_conservative_var[2][0]) +
                                (interior_box_lo_idx[1] + 2 + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                (k + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                    subghostcell_dims_conservative_var[2][1];
                            
                            const Real& rho_y_T   = Q[0][idx_cell_rho_y_T];
                            const Real& rho_y_TT  = Q[0][idx_cell_rho_y_TT];
                            const Real& rho_y_TTT = Q[0][idx_cell_rho_y_TTT];
                            
                            const Real u_y_T   = Q[1][idx_cell_mom_y_T]/rho_y_T;
                            const Real u_y_TT  = Q[1][idx_cell_mom_y_TT]/rho_y_TT;
                            const Real u_y_TTT = Q[1][idx_cell_mom_y_TTT]/rho_y_TTT;
                            
                            const Real v_y_T   = Q[2][idx_cell_mom_y_T]/rho_y_T;
                            const Real v_y_TT  = Q[2][idx_cell_mom_y_TT]/rho_y_TT;
                            const Real v_y_TTT = Q[2][idx_cell_mom_y_TTT]/rho_y_TTT;
                            
                            const Real w_y_T   = Q[3][idx_cell_mom_y_T]/rho_y_T;
                            const Real w_y_TT  = Q[3][idx_cell_mom_y_TT]/rho_y_TT;
                            const Real w_y_TTT = Q[3][idx_cell_mom_y_TTT]/rho_y_TTT;
                            
                            const Real half = Real(1)/Real(2);
                            const Real epsilon_y_T   = Q[4][idx_cell_E_y_T]/rho_y_T - half*(u_y_T*u_y_T + v_y_T*v_y_T + w_y_T*w_y_T);
                            const Real epsilon_y_TT  = Q[4][idx_cell_E_y_TT]/rho_y_TT - half*(u_y_TT*u_y_TT + v_y_TT*v_y_TT + w_y_TT*w_y_TT);
                            const Real epsilon_y_TTT = Q[4][idx_cell_E_y_TTT]/rho_y_TTT - half*(u_y_TTT*u_y_TTT + v_y_TTT*v_y_TTT + w_y_TTT*w_y_TTT);
                            
                            const Real p_y_T = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getPressure(
                                    &rho_y_T,
                                    &epsilon_y_T,
                                    thermo_properties_ptr);
                            
                            const Real p_y_TT = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getPressure(
                                    &rho_y_TT,
                                    &epsilon_y_TT,
                                    thermo_properties_ptr);
                            
                            const Real p_y_TTT = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getPressure(
                                    &rho_y_TTT,
                                    &epsilon_y_TTT,
                                    thermo_properties_ptr);
                            
                            /*
                             * Compute derivatives at y-direction.
                             */
                            
                            const Real drho_dy = -(rho_y_TTT - Real(4)*rho_y_TT + Real(3)*rho_y_T)/(Real(2)*Real(dx[1]));
                            const Real du_dy   = -(u_y_TTT - Real(4)*u_y_TT + Real(3)*u_y_T)/(Real(2)*Real(dx[1]));
                            const Real dv_dy   = -(v_y_TTT - Real(4)*v_y_TT + Real(3)*v_y_T)/(Real(2)*Real(dx[1]));
                            const Real dw_dy   = -(w_y_TTT - Real(4)*w_y_TT + Real(3)*w_y_T)/(Real(2)*Real(dx[1]));
                            const Real dp_dy   = -(p_y_TTT - Real(4)*p_y_TT + Real(3)*p_y_T)/(Real(2)*Real(dx[1]));
                            
                            /*
                             * Compute derivatives in x-direction.
                             */
                            
                            Real du_dx = Real(0);
                            Real dv_dx = Real(0);
                            // Real dw_dx = Real(0);
                            Real dp_dx = Real(0);
                            
                            if ((i + num_subghosts_conservative_var[0][0] == 0) ||
                                (i + num_subghosts_conservative_var[1][0] == 0) ||
                                (i + num_subghosts_conservative_var[2][0] == 0))
                            // if (((patch_geom->getTouchesRegularBoundary(0, 0)) && (i == interior_box_lo_idx[0])) ||
                            //     ((i + num_subghosts_conservative_var[0][0] == 0) ||
                            //      (i + num_subghosts_conservative_var[1][0] == 0) ||
                            //      (i + num_subghosts_conservative_var[2][0] == 0)))
                            {
                                // Patch is touching periodic boundary.
                                
                                const int idx_cell_rho_x_R = (i + 1 + num_subghosts_conservative_var[0][0]) +
                                    (interior_box_lo_idx[1] + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (k + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom_x_R = (i + 1 + num_subghosts_conservative_var[1][0]) +
                                    (interior_box_lo_idx[1] + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                    (k + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E_x_R = (i + 1 + num_subghosts_conservative_var[2][0]) +
                                    (interior_box_lo_idx[1] + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                    (k + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                const Real& rho_x_R = Q[0][idx_cell_rho_x_R];
                                const Real u_x_R = Q[1][idx_cell_mom_x_R]/rho_x_R;
                                const Real v_x_R = Q[2][idx_cell_mom_x_R]/rho_x_R;
                                const Real w_x_R = Q[3][idx_cell_mom_x_R]/rho_x_R;
                                const Real epsilon_x_R = Q[4][idx_cell_E_x_R]/rho_x_R - half*(u_x_R*u_x_R + v_x_R*v_x_R + w_x_R*w_x_R);
                                
                                const Real p_x_R = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &rho_x_R,
                                        &epsilon_x_R,
                                        thermo_properties_ptr);
                                
                                // One-sided derivatives.
                                du_dx = (u_x_R - u_y_T)/Real(dx[0]);
                                dv_dx = (v_x_R - v_y_T)/Real(dx[0]);
                                // dw_dx = (w_x_R - w_y_T)/Real(dx[0]);
                                dp_dx = (p_x_R - p_y_T)/Real(dx[0]);
                            }
                            else if ((i + num_subghosts_conservative_var[0][0] + 1 == subghostcell_dims_conservative_var[0][0]) ||
                                     (i + num_subghosts_conservative_var[1][0] + 1 == subghostcell_dims_conservative_var[1][0]) ||
                                     (i + num_subghosts_conservative_var[2][0] + 1 == subghostcell_dims_conservative_var[2][0]))
                            // else if (((patch_geom->getTouchesRegularBoundary(0, 1)) && (i == interior_box_hi_idx[0])) ||
                            //          ((i + num_subghosts_conservative_var[0][0] + 1 == subghostcell_dims_conservative_var[0][0]) ||
                            //           (i + num_subghosts_conservative_var[1][0] + 1 == subghostcell_dims_conservative_var[1][0]) ||
                            //           (i + num_subghosts_conservative_var[2][0] + 1 == subghostcell_dims_conservative_var[2][0])))
                            {
                                // Patch is touching periodic boundary.
                                
                                const int idx_cell_rho_x_L = (i - 1 + num_subghosts_conservative_var[0][0]) +
                                    (interior_box_lo_idx[1] + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (k + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom_x_L = (i - 1 + num_subghosts_conservative_var[1][0]) +
                                    (interior_box_lo_idx[1] + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                    (k + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E_x_L = (i - 1 + num_subghosts_conservative_var[2][0]) +
                                    (interior_box_lo_idx[1] + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                    (k + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                const Real& rho_x_L = Q[0][idx_cell_rho_x_L];
                                const Real u_x_L = Q[1][idx_cell_mom_x_L]/rho_x_L;
                                const Real v_x_L = Q[2][idx_cell_mom_x_L]/rho_x_L;
                                const Real w_x_L = Q[3][idx_cell_mom_x_L]/rho_x_L;
                                const Real epsilon_x_L = Q[4][idx_cell_E_x_L]/rho_x_L - half*(u_x_L*u_x_L + v_x_L*v_x_L + w_x_L*w_x_L);
                                
                                const Real p_x_L = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &rho_x_L,
                                        &epsilon_x_L,
                                        thermo_properties_ptr);
                                
                                // One-sided derivatives.
                                du_dx = (u_y_T - u_x_L)/Real(dx[0]);
                                dv_dx = (v_y_T - v_x_L)/Real(dx[0]);
                                // dw_dx = (w_y_T - w_x_L)/Real(dx[0]);
                                dp_dx = (p_y_T - p_x_L)/Real(dx[0]);
                            }
                            else
                            {
                                const int idx_cell_rho_x_L = (i - 1 + num_subghosts_conservative_var[0][0]) +
                                    (interior_box_lo_idx[1] + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (k + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_rho_x_R = (i + 1 + num_subghosts_conservative_var[0][0]) +
                                    (interior_box_lo_idx[1] + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (k + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom_x_L = (i - 1 + num_subghosts_conservative_var[1][0]) +
                                    (interior_box_lo_idx[1] + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                    (k + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_mom_x_R = (i + 1 + num_subghosts_conservative_var[1][0]) +
                                    (interior_box_lo_idx[1] + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                    (k + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E_x_L = (i - 1 + num_subghosts_conservative_var[2][0]) +
                                    (interior_box_lo_idx[1] + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                    (k + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                const int idx_cell_E_x_R = (i + 1 + num_subghosts_conservative_var[2][0]) +
                                    (interior_box_lo_idx[1] + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                    (k + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                const Real& rho_x_L = Q[0][idx_cell_rho_x_L];
                                const Real& rho_x_R = Q[0][idx_cell_rho_x_R];
                                
                                const Real u_x_L = Q[1][idx_cell_mom_x_L]/rho_x_L;
                                const Real u_x_R = Q[1][idx_cell_mom_x_R]/rho_x_R;
                                
                                const Real v_x_L = Q[2][idx_cell_mom_x_L]/rho_x_L;
                                const Real v_x_R = Q[2][idx_cell_mom_x_R]/rho_x_R;
                                
                                const Real w_x_L = Q[3][idx_cell_mom_x_L]/rho_x_L;
                                const Real w_x_R = Q[3][idx_cell_mom_x_R]/rho_x_R;
                                
                                const Real epsilon_x_L = Q[4][idx_cell_E_x_L]/rho_x_L - half*(u_x_L*u_x_L + v_x_L*v_x_L + w_x_L*w_x_L);
                                const Real epsilon_x_R = Q[4][idx_cell_E_x_R]/rho_x_R - half*(u_x_R*u_x_R + v_x_R*v_x_R + w_x_R*w_x_R);
                                
                                const Real p_x_L = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &rho_x_L,
                                        &epsilon_x_L,
                                        thermo_properties_ptr);
                                
                                const Real p_x_R = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &rho_x_R,
                                        &epsilon_x_R,
                                        thermo_properties_ptr);
                                
                                // Central derivatives.
                                du_dx = (u_x_R - u_x_L)/(Real(2)*Real(dx[0]));
                                dv_dx = (v_x_R - v_x_L)/(Real(2)*Real(dx[0]));
                                // dw_dx = (w_x_R - w_x_L)/(Real(2)*Real(dx[0]));
                                dp_dx = (p_x_R - p_x_L)/(Real(2)*Real(dx[0]));
                            }
                            
                            /*
                             * Compute derivatives in z-direction.
                             */
                            
                            // Real du_dz = Real(0);
                            Real dv_dz = Real(0);
                            Real dw_dz = Real(0);
                            Real dp_dz = Real(0);
                            
                            if ((k + num_subghosts_conservative_var[0][2] == 0) ||
                                (k + num_subghosts_conservative_var[1][2] == 0) ||
                                (k + num_subghosts_conservative_var[2][2] == 0))
                            // if (((patch_geom->getTouchesRegularBoundary(2, 0)) && (k == interior_box_lo_idx[2])) ||
                            //     ((k + num_subghosts_conservative_var[0][2] == 0) ||
                            //      (k + num_subghosts_conservative_var[1][2] == 0) ||
                            //      (k + num_subghosts_conservative_var[2][2] == 0)))
                            {
                                // Patch is touching periodic boundary.
                                
                                const int idx_cell_rho_z_F = (i + num_subghosts_conservative_var[0][0]) +
                                    (interior_box_lo_idx[1] + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (k + 1 + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom_z_F = (i + num_subghosts_conservative_var[1][0]) +
                                    (interior_box_lo_idx[1] + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                    (k + 1 + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E_z_F = (i + num_subghosts_conservative_var[2][0]) +
                                    (interior_box_lo_idx[1] + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                    (k + 1 + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                const Real& rho_z_F = Q[0][idx_cell_rho_z_F];
                                const Real u_z_F = Q[1][idx_cell_mom_z_F]/rho_z_F;
                                const Real v_z_F = Q[2][idx_cell_mom_z_F]/rho_z_F;
                                const Real w_z_F = Q[3][idx_cell_mom_z_F]/rho_z_F;
                                const Real epsilon_z_F = Q[4][idx_cell_E_z_F]/rho_z_F - half*(u_z_F*u_z_F + v_z_F*v_z_F + w_z_F*w_z_F);
                                
                                const Real p_z_F = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &rho_z_F,
                                        &epsilon_z_F,
                                        thermo_properties_ptr);
                                
                                // One-sided derivatives.
                                // du_dz = (u_z_F - u_y_T)/Real(dx[2]);
                                dv_dz = (v_z_F - v_y_T)/Real(dx[2]);
                                dw_dz = (w_z_F - w_y_T)/Real(dx[2]);
                                dp_dz = (p_z_F - p_y_T)/Real(dx[2]);
                            }
                            else if ((k + num_subghosts_conservative_var[0][2] + 1 == subghostcell_dims_conservative_var[0][2]) ||
                                     (k + num_subghosts_conservative_var[1][2] + 1 == subghostcell_dims_conservative_var[1][2]) ||
                                     (k + num_subghosts_conservative_var[2][2] + 1 == subghostcell_dims_conservative_var[2][2]))
                            // else if (((patch_geom->getTouchesRegularBoundary(2, 1)) && (k == interior_box_hi_idx[2])) ||
                            //          ((k + num_subghosts_conservative_var[0][2] + 1 == subghostcell_dims_conservative_var[0][2]) ||
                            //           (k + num_subghosts_conservative_var[1][2] + 1 == subghostcell_dims_conservative_var[1][2]) ||
                            //           (k + num_subghosts_conservative_var[2][2] + 1 == subghostcell_dims_conservative_var[2][2])))
                            {
                                // Patch is touching periodic boundary.
                                
                                const int idx_cell_rho_z_B = (i + num_subghosts_conservative_var[0][0]) +
                                    (interior_box_lo_idx[1] + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (k - 1 + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom_z_B = (i + num_subghosts_conservative_var[1][0]) +
                                    (interior_box_lo_idx[1] + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                    (k - 1 + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E_z_B = (i + num_subghosts_conservative_var[2][0]) +
                                    (interior_box_lo_idx[1] + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                    (k - 1 + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                const Real& rho_z_B = Q[0][idx_cell_rho_z_B];
                                const Real u_z_B = Q[1][idx_cell_mom_z_B]/rho_z_B;
                                const Real v_z_B = Q[2][idx_cell_mom_z_B]/rho_z_B;
                                const Real w_z_B = Q[3][idx_cell_mom_z_B]/rho_z_B;
                                const Real epsilon_z_B = Q[4][idx_cell_E_z_B]/rho_z_B - half*(u_z_B*u_z_B + v_z_B*v_z_B + w_z_B*w_z_B);
                                
                                const Real p_z_B = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &rho_z_B,
                                        &epsilon_z_B,
                                        thermo_properties_ptr);
                                
                                // One-sided derivatives.
                                // du_dz = (u_y_T - u_z_B)/Real(dx[2]);
                                dv_dz = (v_y_T - v_z_B)/Real(dx[2]);
                                dw_dz = (w_y_T - w_z_B)/Real(dx[2]);
                                dp_dz = (p_y_T - p_z_B)/Real(dx[2]);
                            }
                            else
                            {
                                const int idx_cell_rho_z_B = (i + num_subghosts_conservative_var[0][0]) +
                                    (interior_box_lo_idx[1] + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (k - 1 + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_rho_z_F = (i + num_subghosts_conservative_var[0][0]) +
                                    (interior_box_lo_idx[1] + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (k + 1 + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom_z_B = (i + num_subghosts_conservative_var[1][0]) +
                                    (interior_box_lo_idx[1] + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                    (k - 1 + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_mom_z_F = (i + num_subghosts_conservative_var[1][0]) +
                                    (interior_box_lo_idx[1] + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                    (k + 1 + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E_z_B = (i + num_subghosts_conservative_var[2][0]) +
                                    (interior_box_lo_idx[1] + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                    (k - 1 + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                const int idx_cell_E_z_F = (i + num_subghosts_conservative_var[2][0]) +
                                    (interior_box_lo_idx[1] + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                    (k + 1 + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                const Real& rho_z_B = Q[0][idx_cell_rho_z_B];
                                const Real& rho_z_F = Q[0][idx_cell_rho_z_F];
                                
                                const Real u_z_B = Q[1][idx_cell_mom_z_B]/rho_z_B;
                                const Real u_z_F = Q[1][idx_cell_mom_z_F]/rho_z_F;
                                
                                const Real v_z_B = Q[2][idx_cell_mom_z_B]/rho_z_B;
                                const Real v_z_F = Q[2][idx_cell_mom_z_F]/rho_z_F;
                                
                                const Real w_z_B = Q[3][idx_cell_mom_z_B]/rho_z_B;
                                const Real w_z_F = Q[3][idx_cell_mom_z_F]/rho_z_F;
                                
                                const Real epsilon_z_B = Q[4][idx_cell_E_z_B]/rho_z_B - half*(u_z_B*u_z_B + v_z_B*v_z_B + w_z_B*w_z_B);
                                const Real epsilon_z_F = Q[4][idx_cell_E_z_F]/rho_z_F - half*(u_z_F*u_z_F + v_z_F*v_z_F + w_z_F*w_z_F);
                                
                                const Real p_z_B = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &rho_z_B,
                                        &epsilon_z_B,
                                        thermo_properties_ptr);
                                
                                const Real p_z_F = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &rho_z_F,
                                        &epsilon_z_F,
                                        thermo_properties_ptr);
                                
                                // Central derivatives.
                                // du_dz = (u_z_F - u_z_B)/(Real(2)*Real(dx[2]));
                                dv_dz = (v_z_F - v_z_B)/(Real(2)*Real(dx[2]));
                                dw_dz = (w_z_F - w_z_B)/(Real(2)*Real(dx[2]));
                                dp_dz = (p_z_F - p_z_B)/(Real(2)*Real(dx[2]));
                            }
                            
                            const Real c_y_T = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getSoundSpeed(
                                    &rho_y_T,
                                    &p_y_T,
                                    thermo_properties_ptr);
                            
                            const Real lambda_5 = v_y_T + c_y_T;
                            
                            // Compute vector Lambda^(-1) * L.
                            
                            Real Lambda_inv_L[5];
                            
                            const Real& p_t         = d_bdry_face_nonreflecting_outflow_p_t[face_loc];
                            const Real& sigma       = d_bdry_face_nonreflecting_outflow_sigma[face_loc];
                            const Real& beta        = d_bdry_face_nonreflecting_outflow_beta[face_loc];
                            const Real& length_char = d_bdry_face_nonreflecting_outflow_length_char[face_loc];
                            
                            const Real T_5 = u_y_T*(dp_dx + rho_y_T*c_y_T*dv_dx) + rho_y_T*c_y_T*c_y_T*du_dx +
                                w_y_T*(dp_dz + rho_y_T*c_y_T*dv_dz) + rho_y_T*c_y_T*c_y_T*dw_dz;
                            
                            const Real M_sq = (u_y_T*u_y_T + v_y_T*v_y_T + w_y_T*w_y_T)/(c_y_T*c_y_T);
                            const Real K = sigma*c_y_T*(Real(1) - M_sq)/length_char;
                            
                            Lambda_inv_L[0] = dp_dy - rho_y_T*c_y_T*dv_dy;
                            Lambda_inv_L[1] = du_dy;
                            Lambda_inv_L[2] = c_y_T*c_y_T*drho_dy - dp_dy;
                            Lambda_inv_L[3] = dw_dy;
                            Lambda_inv_L[4] = (Real(1)/lambda_5)*(K*(p_y_T - p_t) - (Real(1) - beta)*T_5);
                            
                            // Compute dV_dy.
                            
                            const Real c_sq_inv  = Real(1)/(c_y_T*c_y_T);
                            const Real rho_c_inv = Real(1)/(rho_y_T*c_y_T);
                            
                            Real dV_dy[5];
                            
                            dV_dy[0] = half*c_sq_inv*(Lambda_inv_L[0] + Lambda_inv_L[4]) + c_sq_inv*Lambda_inv_L[2];
                            dV_dy[1] = Lambda_inv_L[1];
                            dV_dy[2] = half*rho_c_inv*(-Lambda_inv_L[0] + Lambda_inv_L[4]);
                            dV_dy[3] = Lambda_inv_L[3];
                            dV_dy[4] = half*(Lambda_inv_L[0] + Lambda_inv_L[4]);
                            
                            Real V_ghost[5*num_ghosts_to_fill];
                            
                            for (int j = num_ghosts_to_fill - 1; j >= 0; j--)
                            {
                                const int idx_cell_rho = (i + num_subghosts_conservative_var[0][0]) +
                                    (j + fill_box_lo_idx[1] + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (k + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom = (i + num_subghosts_conservative_var[1][0]) +
                                    (j + fill_box_lo_idx[1] + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                    (k + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E = (i + num_subghosts_conservative_var[2][0]) +
                                    (j + fill_box_lo_idx[1] + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                    (k + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                if (j == num_ghosts_to_fill - 1)
                                {
                                    V_ghost[j*5 + 0] = rho_y_TT - Real(2)*Real(dx[1])*dV_dy[0];
                                    V_ghost[j*5 + 1] = u_y_TT   - Real(2)*Real(dx[1])*dV_dy[1];
                                    V_ghost[j*5 + 2] = v_y_TT   - Real(2)*Real(dx[1])*dV_dy[2];
                                    V_ghost[j*5 + 3] = w_y_TT   - Real(2)*Real(dx[1])*dV_dy[3];
                                    V_ghost[j*5 + 4] = p_y_TT   - Real(2)*Real(dx[1])*dV_dy[4];
                                }
                                else if (j == num_ghosts_to_fill - 2)
                                {
                                    V_ghost[j*5 + 0] = -Real(2)*rho_y_TT - Real(3)*rho_y_T +
                                        Real(6)*V_ghost[(j + 1)*5 + 0] + Real(6)*Real(dx[1])*dV_dy[0];
                                    
                                    V_ghost[j*5 + 1] = -Real(2)*u_y_TT - Real(3)*u_y_T +
                                        Real(6)*V_ghost[(j + 1)*5 + 1] + Real(6)*Real(dx[1])*dV_dy[1];
                                    
                                    V_ghost[j*5 + 2] = -Real(2)*v_y_TT - Real(3)*v_y_T +
                                        Real(6)*V_ghost[(j + 1)*5 + 2] + Real(6)*Real(dx[1])*dV_dy[2];
                                    
                                    V_ghost[j*5 + 3] = -Real(2)*w_y_TT - Real(3)*w_y_T +
                                        Real(6)*V_ghost[(j + 1)*5 + 3] + Real(6)*Real(dx[1])*dV_dy[3];
                                    
                                    V_ghost[j*5 + 4] = -Real(2)*p_y_TT - Real(3)*p_y_T +
                                        Real(6)*V_ghost[(j + 1)*5 + 4] + Real(6)*Real(dx[1])*dV_dy[4];
                                }
                                else if (j == num_ghosts_to_fill - 3)
                                {
                                    V_ghost[j*5 + 0] = Real(3)*rho_y_TT + Real(10)*rho_y_T -
                                        Real(18)*V_ghost[(j + 2)*5 + 0] + Real(6)*V_ghost[(j + 1)*5 + 0] -
                                        Real(12)*Real(dx[1])*dV_dy[0];
                                    
                                    V_ghost[j*5 + 1] = Real(3)*u_y_TT + Real(10)*u_y_T -
                                        Real(18)*V_ghost[(j + 2)*5 + 1] + Real(6)*V_ghost[(j + 1)*5 + 1] -
                                        Real(12)*Real(dx[1])*dV_dy[1];
                                    
                                    V_ghost[j*5 + 2] = Real(3)*v_y_TT + Real(10)*v_y_T -
                                        Real(18)*V_ghost[(j + 2)*5 + 2] + Real(6)*V_ghost[(j + 1)*5 + 2] -
                                        Real(12)*Real(dx[1])*dV_dy[2];
                                    
                                    V_ghost[j*5 + 3] = Real(3)*w_y_TT + Real(10)*w_y_T -
                                        Real(18)*V_ghost[(j + 2)*5 + 3] + Real(6)*V_ghost[(j + 1)*5 + 3] -
                                        Real(12)*Real(dx[1])*dV_dy[3];
                                    
                                    V_ghost[j*5 + 4] = Real(3)*p_y_TT + Real(10)*p_y_T -
                                        Real(18)*V_ghost[(j + 2)*5 + 4] + Real(6)*V_ghost[(j + 1)*5 + 4] -
                                        Real(12)*Real(dx[1])*dV_dy[4];
                                }
                                else if (j == num_ghosts_to_fill - 4)
                                {
                                    V_ghost[j*5 + 0] = -Real(4)*rho_y_TT - Real(65)/Real(3)*rho_y_T +
                                        Real(40)*V_ghost[(j + 3)*5 + 0] - Real(20)*V_ghost[(j + 2)*5 + 0] +
                                        Real(20)/Real(3)*V_ghost[(j + 1)*5 + 0] + Real(20)*Real(dx[1])*dV_dy[0];
                                    
                                    V_ghost[j*5 + 1] = -Real(4)*u_y_TT - Real(65)/Real(3)*u_y_T +
                                        Real(40)*V_ghost[(j + 3)*5 + 1] - Real(20)*V_ghost[(j + 2)*5 + 1] +
                                        Real(20)/Real(3)*V_ghost[(j + 1)*5 + 1] + Real(20)*Real(dx[1])*dV_dy[1];
                                    
                                    V_ghost[j*5 + 2] = -Real(4)*v_y_TT - Real(65)/Real(3)*v_y_T +
                                        Real(40)*V_ghost[(j + 3)*5 + 2] - Real(20)*V_ghost[(j + 2)*5 + 2] +
                                        Real(20)/Real(3)*V_ghost[(j + 1)*5 + 2] + Real(20)*Real(dx[1])*dV_dy[2];
                                    
                                    V_ghost[j*5 + 3] = -Real(4)*w_y_TT - Real(65)/Real(3)*w_y_T +
                                        Real(40)*V_ghost[(j + 3)*5 + 3] - Real(20)*V_ghost[(j + 2)*5 + 3] +
                                        Real(20)/Real(3)*V_ghost[(j + 1)*5 + 3] + Real(20)*Real(dx[1])*dV_dy[3];
                                    
                                    V_ghost[j*5 + 4] = -Real(4)*p_y_TT - Real(65)/Real(3)*p_y_T +
                                        Real(40)*V_ghost[(j + 3)*5 + 4] - Real(20)*V_ghost[(j + 2)*5 + 4] +
                                        Real(20)/Real(3)*V_ghost[(j + 1)*5 + 4] + Real(20)*Real(dx[1])*dV_dy[4];
                                }
                                else if (j == num_ghosts_to_fill - 5)
                                {
                                    V_ghost[j*5 + 0] = Real(5)*rho_y_TT + Real(77)/Real(2)*rho_y_T -
                                        Real(75)*V_ghost[(j + 4)*5 + 0] + Real(50)*V_ghost[(j + 3)*5 + 0] -
                                        Real(25)*V_ghost[(j + 2)*5 + 0] + Real(15)/Real(2)*V_ghost[(j + 1)*5 + 0] -
                                        Real(30)*Real(dx[1])*dV_dy[0];
                                    
                                    V_ghost[j*5 + 1] = Real(5)*u_y_TT + Real(77)/Real(2)*u_y_T -
                                        Real(75)*V_ghost[(j + 4)*5 + 1] + Real(50)*V_ghost[(j + 3)*5 + 1] -
                                        Real(25)*V_ghost[(j + 2)*5 + 1] + Real(15)/Real(2)*V_ghost[(j + 1)*5 + 1] -
                                        Real(30)*Real(dx[1])*dV_dy[1];
                                    
                                    V_ghost[j*5 + 2] = Real(5)*v_y_TT + Real(77)/Real(2)*v_y_T -
                                        Real(75)*V_ghost[(j + 4)*5 + 2] + Real(50)*V_ghost[(j + 3)*5 + 2] -
                                        Real(25)*V_ghost[(j + 2)*5 + 2] + Real(15)/Real(2)*V_ghost[(j + 1)*5 + 2] -
                                        Real(30)*Real(dx[1])*dV_dy[2];
                                    
                                    V_ghost[j*5 + 3] = Real(5)*w_y_TT + Real(77)/Real(2)*w_y_T -
                                        Real(75)*V_ghost[(j + 4)*5 + 3] + Real(50)*V_ghost[(j + 3)*5 + 3] -
                                        Real(25)*V_ghost[(j + 2)*5 + 3] + Real(15)/Real(2)*V_ghost[(j + 1)*5 + 3] -
                                        Real(30)*Real(dx[1])*dV_dy[3];
                                    
                                    V_ghost[j*5 + 4] = Real(5)*p_y_TT + Real(77)/Real(2)*p_y_T -
                                        Real(75)*V_ghost[(j + 4)*5 + 4] + Real(50)*V_ghost[(j + 3)*5 + 4] -
                                        Real(25)*V_ghost[(j + 2)*5 + 4] + Real(15)/Real(2)*V_ghost[(j + 1)*5 + 4] -
                                        Real(30)*Real(dx[1])*dV_dy[4];
                                }
                                else if (j == num_ghosts_to_fill - 6)
                                {
                                    V_ghost[j*5 + 0] = -Real(6)*rho_y_TT - Real(609)/Real(10)*rho_y_T +
                                        Real(126)*V_ghost[(j + 5)*5 + 0] - Real(105)*V_ghost[(j + 4)*5 + 0] +
                                        Real(70)*V_ghost[(j + 3)*5 + 0] - Real(63)/Real(2)*V_ghost[(j + 2)*5 + 0] +
                                        Real(42)/Real(5)*V_ghost[(j + 1)*5 + 0] + Real(42)*Real(dx[1])*dV_dy[0];
                                    
                                    V_ghost[j*5 + 1] = -Real(6)*u_y_TT - Real(609)/Real(10)*u_y_T +
                                        Real(126)*V_ghost[(j + 5)*5 + 1] - Real(105)*V_ghost[(j + 4)*5 + 1] +
                                        Real(70)*V_ghost[(j + 3)*5 + 1] - Real(63)/Real(2)*V_ghost[(j + 2)*5 + 1] +
                                        Real(42)/Real(5)*V_ghost[(j + 1)*5 + 1] + Real(42)*Real(dx[1])*dV_dy[1];
                                    
                                    V_ghost[j*5 + 2] = -Real(6)*v_y_TT - Real(609)/Real(10)*v_y_T +
                                        Real(126)*V_ghost[(j + 5)*5 + 2] - Real(105)*V_ghost[(j + 4)*5 + 2] +
                                        Real(70)*V_ghost[(j + 3)*5 + 2] - Real(63)/Real(2)*V_ghost[(j + 2)*5 + 2] +
                                        Real(42)/Real(5)*V_ghost[(j + 1)*5 + 2] + Real(42)*Real(dx[1])*dV_dy[2];
                                    
                                    V_ghost[j*5 + 3] = -Real(6)*w_y_TT - Real(609)/Real(10)*w_y_T +
                                        Real(126)*V_ghost[(j + 5)*5 + 3] - Real(105)*V_ghost[(j + 4)*5 + 3] +
                                        Real(70)*V_ghost[(j + 3)*5 + 3] - Real(63)/Real(2)*V_ghost[(j + 2)*5 + 3] +
                                        Real(42)/Real(5)*V_ghost[(j + 1)*5 + 3] + Real(42)*Real(dx[1])*dV_dy[3];
                                    
                                    V_ghost[j*5 + 4] = -Real(6)*p_y_TT - Real(609)/Real(10)*p_y_T +
                                        Real(126)*V_ghost[(j + 5)*5 + 4] - Real(105)*V_ghost[(j + 4)*5 + 4] +
                                        Real(70)*V_ghost[(j + 3)*5 + 4] - Real(63)/Real(2)*V_ghost[(j + 2)*5 + 4] +
                                        Real(42)/Real(5)*V_ghost[(j + 1)*5 + 4] + Real(42)*Real(dx[1])*dV_dy[4];
                                }
                                
                                Q[0][idx_cell_rho] = V_ghost[j*5 + 0];
                                Q[1][idx_cell_mom] = V_ghost[j*5 + 0]*V_ghost[j*5 + 1];
                                Q[2][idx_cell_mom] = V_ghost[j*5 + 0]*V_ghost[j*5 + 2];
                                Q[3][idx_cell_mom] = V_ghost[j*5 + 0]*V_ghost[j*5 + 3];
                                
                                const Real epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getInternalEnergy(
                                        &V_ghost[j*5 + 0],
                                        &V_ghost[j*5 + 4],
                                        thermo_properties_ptr);
                                
                                const Real E = V_ghost[j*5 + 0]*epsilon +
                                    half*(Q[1][idx_cell_mom]*Q[1][idx_cell_mom] + Q[2][idx_cell_mom]*Q[2][idx_cell_mom] +
                                       Q[3][idx_cell_mom]*Q[3][idx_cell_mom])/V_ghost[j*5 + 0];
                                
                                Q[4][idx_cell_E] = E;
                            }
                        }
                    }
                }
                else if (face_loc == BDRY_LOC::YHI)
                {
                    const int num_ghosts_to_fill = fill_box_hi_idx[1] - fill_box_lo_idx[1] + 1;
                    TBOX_ASSERT(fill_box_lo_idx[1] == interior_box_hi_idx[1] + 1);
                    if (num_ghosts_to_fill > 6)
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelBoundaryUtilitiesSingleSpecies::fill3dFaceBoundaryData()\n"
                            << "Non-reflecting outflow BC doesn't support more than six ghost cells yet!");
                    }
                    
                    for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                    {
                        for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                        {
                            // Get the grid spacing.
                            const double* const dx = patch_geom->getDx();
                            
                            const int idx_cell_rho_y_B = (i + num_subghosts_conservative_var[0][0]) +
                                (interior_box_hi_idx[1] + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                (k + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                    subghostcell_dims_conservative_var[0][1];
                            
                            const int idx_cell_rho_y_BB = (i + num_subghosts_conservative_var[0][0]) +
                                (interior_box_hi_idx[1] - 1 + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                (k + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                    subghostcell_dims_conservative_var[0][1];
                            
                            const int idx_cell_rho_y_BBB = (i + num_subghosts_conservative_var[0][0]) +
                                (interior_box_hi_idx[1] - 2 + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                (k + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                    subghostcell_dims_conservative_var[0][1];
                            
                            const int idx_cell_mom_y_B = (i + num_subghosts_conservative_var[1][0]) +
                                ( interior_box_hi_idx[1] + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                (k + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                    subghostcell_dims_conservative_var[1][1];
                            
                            const int idx_cell_mom_y_BB = (i + num_subghosts_conservative_var[1][0]) +
                                (interior_box_hi_idx[1] - 1 + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                (k + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                    subghostcell_dims_conservative_var[1][1];
                            
                            const int idx_cell_mom_y_BBB = (i + num_subghosts_conservative_var[1][0]) +
                                (interior_box_hi_idx[1] - 2 + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                (k + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                    subghostcell_dims_conservative_var[1][1];
                            
                            const int idx_cell_E_y_B = (i + num_subghosts_conservative_var[2][0]) +
                                (interior_box_hi_idx[1] + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                (k + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                    subghostcell_dims_conservative_var[2][1];
                            
                            const int idx_cell_E_y_BB = (i + num_subghosts_conservative_var[2][0]) +
                                (interior_box_hi_idx[1] - 1 + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                (k + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                    subghostcell_dims_conservative_var[2][1];
                            
                            const int idx_cell_E_y_BBB = (i + num_subghosts_conservative_var[2][0]) +
                                (interior_box_hi_idx[1] - 2 + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                (k + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                    subghostcell_dims_conservative_var[2][1];
                            
                            const Real& rho_y_B   = Q[0][idx_cell_rho_y_B];
                            const Real& rho_y_BB  = Q[0][idx_cell_rho_y_BB];
                            const Real& rho_y_BBB = Q[0][idx_cell_rho_y_BBB];
                            
                            const Real u_y_B   = Q[1][idx_cell_mom_y_B]/rho_y_B;
                            const Real u_y_BB  = Q[1][idx_cell_mom_y_BB]/rho_y_BB;
                            const Real u_y_BBB = Q[1][idx_cell_mom_y_BBB]/rho_y_BBB;
                            
                            const Real v_y_B   = Q[2][idx_cell_mom_y_B]/rho_y_B;
                            const Real v_y_BB  = Q[2][idx_cell_mom_y_BB]/rho_y_BB;
                            const Real v_y_BBB = Q[2][idx_cell_mom_y_BBB]/rho_y_BBB;
                            
                            const Real w_y_B   = Q[3][idx_cell_mom_y_B]/rho_y_B;
                            const Real w_y_BB  = Q[3][idx_cell_mom_y_BB]/rho_y_BB;
                            const Real w_y_BBB = Q[3][idx_cell_mom_y_BBB]/rho_y_BBB;
                            
                            const Real half = Real(1)/Real(2);
                            const Real epsilon_y_B   = Q[4][idx_cell_E_y_B]/rho_y_B - half*(u_y_B*u_y_B + v_y_B*v_y_B + w_y_B*w_y_B);
                            const Real epsilon_y_BB  = Q[4][idx_cell_E_y_BB]/rho_y_BB - half*(u_y_BB*u_y_BB + v_y_BB*v_y_BB + w_y_BB*w_y_BB);
                            const Real epsilon_y_BBB = Q[4][idx_cell_E_y_BBB]/rho_y_BBB - half*(u_y_BBB*u_y_BBB + v_y_BBB*v_y_BBB + w_y_BBB*w_y_BBB);
                            
                            const Real p_y_B = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getPressure(
                                    &rho_y_B,
                                    &epsilon_y_B,
                                    thermo_properties_ptr);
                            
                            const Real p_y_BB = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getPressure(
                                    &rho_y_BB,
                                    &epsilon_y_BB,
                                    thermo_properties_ptr);
                            
                            const Real p_y_BBB = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getPressure(
                                    &rho_y_BBB,
                                    &epsilon_y_BBB,
                                    thermo_properties_ptr);
                            
                            /*
                             * Compute derivatives at y-direction.
                             */
                            
                            const Real drho_dy = (rho_y_BBB - Real(4)*rho_y_BB + Real(3)*rho_y_B)/(Real(2)*Real(dx[1]));
                            const Real du_dy   = (u_y_BBB - Real(4)*u_y_BB + Real(3)*u_y_B)/(Real(2)*Real(dx[1]));
                            const Real dv_dy   = (v_y_BBB - Real(4)*v_y_BB + Real(3)*v_y_B)/(Real(2)*Real(dx[1]));
                            const Real dw_dy   = (w_y_BBB - Real(4)*w_y_BB + Real(3)*w_y_B)/(Real(2)*Real(dx[1]));
                            const Real dp_dy   = (p_y_BBB - Real(4)*p_y_BB + Real(3)*p_y_B)/(Real(2)*Real(dx[1]));
                            
                            /*
                             * Compute derivatives in x-direction.
                             */
                            
                            Real du_dx = Real(0);
                            Real dv_dx = Real(0);
                            // Real dw_dx = Real(0);
                            Real dp_dx = Real(0);
                            
                            if ((i + num_subghosts_conservative_var[0][0] == 0) ||
                                (i + num_subghosts_conservative_var[1][0] == 0) ||
                                (i + num_subghosts_conservative_var[2][0] == 0))
                            // if (((patch_geom->getTouchesRegularBoundary(0, 0)) && (i == interior_box_lo_idx[0])) ||
                            //     ((i + num_subghosts_conservative_var[0][0] == 0) ||
                            //      (i + num_subghosts_conservative_var[1][0] == 0) ||
                            //      (i + num_subghosts_conservative_var[2][0] == 0)))
                            {
                                // Patch is touching periodic boundary.
                                
                                const int idx_cell_rho_x_R = (i + 1 + num_subghosts_conservative_var[0][0]) +
                                    (interior_box_hi_idx[1] + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (k + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom_x_R = (i + 1 + num_subghosts_conservative_var[1][0]) +
                                    (interior_box_hi_idx[1] + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                    (k + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E_x_R = (i + 1 + num_subghosts_conservative_var[2][0]) +
                                    (interior_box_hi_idx[1] + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                    (k + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                const Real& rho_x_R = Q[0][idx_cell_rho_x_R];
                                const Real u_x_R = Q[1][idx_cell_mom_x_R]/rho_x_R;
                                const Real v_x_R = Q[2][idx_cell_mom_x_R]/rho_x_R;
                                const Real w_x_R = Q[3][idx_cell_mom_x_R]/rho_x_R;
                                const Real epsilon_x_R = Q[4][idx_cell_E_x_R]/rho_x_R - half*(u_x_R*u_x_R + v_x_R*v_x_R + w_x_R*w_x_R);
                                
                                const Real p_x_R = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &rho_x_R,
                                        &epsilon_x_R,
                                        thermo_properties_ptr);
                                
                                // One-sided derivatives.
                                du_dx = (u_x_R - u_y_B)/Real(dx[0]);
                                dv_dx = (v_x_R - v_y_B)/Real(dx[0]);
                                // dw_dx = (w_x_R - w_y_B)/Real(dx[0]);
                                dp_dx = (p_x_R - p_y_B)/Real(dx[0]);
                            }
                            else if ((i + num_subghosts_conservative_var[0][0] + 1 == subghostcell_dims_conservative_var[0][0]) ||
                                     (i + num_subghosts_conservative_var[1][0] + 1 == subghostcell_dims_conservative_var[1][0]) ||
                                     (i + num_subghosts_conservative_var[2][0] + 1 == subghostcell_dims_conservative_var[2][0]))
                            // else if (((patch_geom->getTouchesRegularBoundary(0, 1)) && (i == interior_box_hi_idx[0])) ||
                            //          ((i + num_subghosts_conservative_var[0][0] + 1 == subghostcell_dims_conservative_var[0][0]) ||
                            //           (i + num_subghosts_conservative_var[1][0] + 1 == subghostcell_dims_conservative_var[1][0]) ||
                            //           (i + num_subghosts_conservative_var[2][0] + 1 == subghostcell_dims_conservative_var[2][0])))
                            {
                                // Patch is touching periodic boundary.
                                
                                const int idx_cell_rho_x_L = (i - 1 + num_subghosts_conservative_var[0][0]) +
                                    (interior_box_hi_idx[1] + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (k + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom_x_L = (i - 1 + num_subghosts_conservative_var[1][0]) +
                                    (interior_box_hi_idx[1] + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                    (k + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E_x_L = (i - 1 + num_subghosts_conservative_var[2][0]) +
                                    (interior_box_hi_idx[1] + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                    (k + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                const Real& rho_x_L = Q[0][idx_cell_rho_x_L];
                                const Real u_x_L = Q[1][idx_cell_mom_x_L]/rho_x_L;
                                const Real v_x_L = Q[2][idx_cell_mom_x_L]/rho_x_L;
                                const Real w_x_L = Q[3][idx_cell_mom_x_L]/rho_x_L;
                                const Real epsilon_x_L = Q[4][idx_cell_E_x_L]/rho_x_L - half*(u_x_L*u_x_L + v_x_L*v_x_L + w_x_L*w_x_L);
                                
                                const Real p_x_L = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &rho_x_L,
                                        &epsilon_x_L,
                                        thermo_properties_ptr);
                                
                                // One-sided derivatives.
                                du_dx = (u_y_B - u_x_L)/Real(dx[0]);
                                dv_dx = (v_y_B - v_x_L)/Real(dx[0]);
                                // dw_dx = (w_y_B - w_x_L)/Real(dx[0]);
                                dp_dx = (p_y_B - p_x_L)/Real(dx[0]);
                            }
                            else
                            {
                                const int idx_cell_rho_x_L = (i - 1 + num_subghosts_conservative_var[0][0]) +
                                    (interior_box_hi_idx[1] + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (k + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_rho_x_R = (i + 1 + num_subghosts_conservative_var[0][0]) +
                                    (interior_box_hi_idx[1] + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (k + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom_x_L = (i - 1 + num_subghosts_conservative_var[1][0]) +
                                    (interior_box_hi_idx[1] + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                    (k + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_mom_x_R = (i + 1 + num_subghosts_conservative_var[1][0]) +
                                    (interior_box_hi_idx[1] + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                    (k + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E_x_L = (i - 1 + num_subghosts_conservative_var[2][0]) +
                                    (interior_box_hi_idx[1] + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                    (k + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                const int idx_cell_E_x_R = (i + 1 + num_subghosts_conservative_var[2][0]) +
                                    (interior_box_hi_idx[1] + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                    (k + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                const Real& rho_x_L = Q[0][idx_cell_rho_x_L];
                                const Real& rho_x_R = Q[0][idx_cell_rho_x_R];
                                
                                const Real u_x_L = Q[1][idx_cell_mom_x_L]/rho_x_L;
                                const Real u_x_R = Q[1][idx_cell_mom_x_R]/rho_x_R;
                                
                                const Real v_x_L = Q[2][idx_cell_mom_x_L]/rho_x_L;
                                const Real v_x_R = Q[2][idx_cell_mom_x_R]/rho_x_R;
                                
                                const Real w_x_L = Q[3][idx_cell_mom_x_L]/rho_x_L;
                                const Real w_x_R = Q[3][idx_cell_mom_x_R]/rho_x_R;
                                
                                const Real epsilon_x_L = Q[4][idx_cell_E_x_L]/rho_x_L - half*(u_x_L*u_x_L + v_x_L*v_x_L + w_x_L*w_x_L);
                                const Real epsilon_x_R = Q[4][idx_cell_E_x_R]/rho_x_R - half*(u_x_R*u_x_R + v_x_R*v_x_R + w_x_R*w_x_R);
                                
                                const Real p_x_L = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &rho_x_L,
                                        &epsilon_x_L,
                                        thermo_properties_ptr);
                                
                                const Real p_x_R = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &rho_x_R,
                                        &epsilon_x_R,
                                        thermo_properties_ptr);
                                
                                // Central derivatives.
                                du_dx = (u_x_R - u_x_L)/(Real(2)*Real(dx[0]));
                                dv_dx = (v_x_R - v_x_L)/(Real(2)*Real(dx[0]));
                                // dw_dx = (w_x_R - w_x_L)/(Real(2)*Real(dx[0]));
                                dp_dx = (p_x_R - p_x_L)/(Real(2)*Real(dx[0]));
                            }
                            
                            /*
                             * Compute derivatives in z-direction.
                             */
                            
                            // Real du_dz = Real(0);
                            Real dv_dz = Real(0);
                            Real dw_dz = Real(0);
                            Real dp_dz = Real(0);
                            
                            if ((k + num_subghosts_conservative_var[0][2] == 0) ||
                                (k + num_subghosts_conservative_var[1][2] == 0) ||
                                (k + num_subghosts_conservative_var[2][2] == 0))
                            // if (((patch_geom->getTouchesRegularBoundary(2, 0)) && (k == interior_box_lo_idx[2])) ||
                            //     ((k + num_subghosts_conservative_var[0][2] == 0) ||
                            //      (k + num_subghosts_conservative_var[1][2] == 0) ||
                            //      (k + num_subghosts_conservative_var[2][2] == 0)))
                            {
                                // Patch is touching periodic boundary.
                                
                                const int idx_cell_rho_z_F = (i + num_subghosts_conservative_var[0][0]) +
                                    (interior_box_hi_idx[1] + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (k + 1 + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom_z_F = (i + num_subghosts_conservative_var[1][0]) +
                                    (interior_box_hi_idx[1] + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                    (k + 1 + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E_z_F = (i + num_subghosts_conservative_var[2][0]) +
                                    (interior_box_hi_idx[1] + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                    (k + 1 + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                const Real& rho_z_F = Q[0][idx_cell_rho_z_F];
                                const Real u_z_F = Q[1][idx_cell_mom_z_F]/rho_z_F;
                                const Real v_z_F = Q[2][idx_cell_mom_z_F]/rho_z_F;
                                const Real w_z_F = Q[3][idx_cell_mom_z_F]/rho_z_F;
                                const Real epsilon_z_F = Q[4][idx_cell_E_z_F]/rho_z_F - half*(u_z_F*u_z_F + v_z_F*v_z_F + w_z_F*w_z_F);
                                
                                const Real p_z_F = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &rho_z_F,
                                        &epsilon_z_F,
                                        thermo_properties_ptr);
                                
                                // One-sided derivatives.
                                // du_dz = (u_z_F - u_y_B)/Real(dx[2]);
                                dv_dz = (v_z_F - v_y_B)/Real(dx[2]);
                                dw_dz = (w_z_F - w_y_B)/Real(dx[2]);
                                dp_dz = (p_z_F - p_y_B)/Real(dx[2]);
                            }
                            else if ((k + num_subghosts_conservative_var[0][2] + 1 == subghostcell_dims_conservative_var[0][2]) ||
                                     (k + num_subghosts_conservative_var[1][2] + 1 == subghostcell_dims_conservative_var[1][2]) ||
                                     (k + num_subghosts_conservative_var[2][2] + 1 == subghostcell_dims_conservative_var[2][2]))
                            // else if (((patch_geom->getTouchesRegularBoundary(2, 1)) && (k == interior_box_hi_idx[2])) ||
                            //          ((k + num_subghosts_conservative_var[0][2] + 1 == subghostcell_dims_conservative_var[0][2]) ||
                            //           (k + num_subghosts_conservative_var[1][2] + 1 == subghostcell_dims_conservative_var[1][2]) ||
                            //           (k + num_subghosts_conservative_var[2][2] + 1 == subghostcell_dims_conservative_var[2][2])))
                            {
                                // Patch is touching periodic boundary.
                                
                                const int idx_cell_rho_z_B = (i + num_subghosts_conservative_var[0][0]) +
                                    (interior_box_hi_idx[1] + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (k - 1 + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom_z_B = (i + num_subghosts_conservative_var[1][0]) +
                                    (interior_box_hi_idx[1] + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                    (k - 1 + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E_z_B = (i + num_subghosts_conservative_var[2][0]) +
                                    (interior_box_hi_idx[1] + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                    (k - 1 + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                const Real& rho_z_B = Q[0][idx_cell_rho_z_B];
                                const Real u_z_B = Q[1][idx_cell_mom_z_B]/rho_z_B;
                                const Real v_z_B = Q[2][idx_cell_mom_z_B]/rho_z_B;
                                const Real w_z_B = Q[3][idx_cell_mom_z_B]/rho_z_B;
                                const Real epsilon_z_B = Q[4][idx_cell_E_z_B]/rho_z_B - half*(u_z_B*u_z_B + v_z_B*v_z_B + w_z_B*w_z_B);
                                
                                const Real p_z_B = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &rho_z_B,
                                        &epsilon_z_B,
                                        thermo_properties_ptr);
                                
                                // One-sided derivatives.
                                // du_dz = (u_y_B - u_z_B)/Real(dx[2]);
                                dv_dz = (v_y_B - v_z_B)/Real(dx[2]);
                                dw_dz = (w_y_B - w_z_B)/Real(dx[2]);
                                dp_dz = (p_y_B - p_z_B)/Real(dx[2]);
                            }
                            else
                            {
                                const int idx_cell_rho_z_B = (i + num_subghosts_conservative_var[0][0]) +
                                    (interior_box_lo_idx[1] + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (k - 1 + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_rho_z_F = (i + num_subghosts_conservative_var[0][0]) +
                                    (interior_box_lo_idx[1] + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (k + 1 + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom_z_B = (i + num_subghosts_conservative_var[1][0]) +
                                    (interior_box_lo_idx[1] + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                    (k - 1 + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_mom_z_F = (i + num_subghosts_conservative_var[1][0]) +
                                    (interior_box_lo_idx[1] + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                    (k + 1 + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E_z_B = (i + num_subghosts_conservative_var[2][0]) +
                                    (interior_box_lo_idx[1] + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                    (k - 1 + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                const int idx_cell_E_z_F = (i + num_subghosts_conservative_var[2][0]) +
                                    (interior_box_lo_idx[1] + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                    (k + 1 + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                const Real& rho_z_B = Q[0][idx_cell_rho_z_B];
                                const Real& rho_z_F = Q[0][idx_cell_rho_z_F];
                                
                                const Real u_z_B = Q[1][idx_cell_mom_z_B]/rho_z_B;
                                const Real u_z_F = Q[1][idx_cell_mom_z_F]/rho_z_F;
                                
                                const Real v_z_B = Q[2][idx_cell_mom_z_B]/rho_z_B;
                                const Real v_z_F = Q[2][idx_cell_mom_z_F]/rho_z_F;
                                
                                const Real w_z_B = Q[3][idx_cell_mom_z_B]/rho_z_B;
                                const Real w_z_F = Q[3][idx_cell_mom_z_F]/rho_z_F;
                                
                                const Real epsilon_z_B = Q[4][idx_cell_E_z_B]/rho_z_B - half*(u_z_B*u_z_B + v_z_B*v_z_B + w_z_B*w_z_B);
                                const Real epsilon_z_F = Q[4][idx_cell_E_z_F]/rho_z_F - half*(u_z_F*u_z_F + v_z_F*v_z_F + w_z_F*w_z_F);
                                
                                const Real p_z_B = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &rho_z_B,
                                        &epsilon_z_B,
                                        thermo_properties_ptr);
                                
                                const Real p_z_F = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &rho_z_F,
                                        &epsilon_z_F,
                                        thermo_properties_ptr);
                                
                                // Central derivatives.
                                // du_dz = (u_z_F - u_z_B)/(Real(2)*Real(dx[2]));
                                dv_dz = (v_z_F - v_z_B)/(Real(2)*Real(dx[2]));
                                dw_dz = (w_z_F - w_z_B)/(Real(2)*Real(dx[2]));
                                dp_dz = (p_z_F - p_z_B)/(Real(2)*Real(dx[2]));
                            }
                            
                            const Real c_y_B = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getSoundSpeed(
                                    &rho_y_B,
                                    &p_y_B,
                                    thermo_properties_ptr);
                            
                            const Real lambda_1 = v_y_B - c_y_B;
                            
                            // Compute vector Lambda^(-1) * L.
                            
                            Real Lambda_inv_L[5];
                            
                            const Real& p_t         = d_bdry_face_nonreflecting_outflow_p_t[face_loc];
                            const Real& sigma       = d_bdry_face_nonreflecting_outflow_sigma[face_loc];
                            const Real& beta        = d_bdry_face_nonreflecting_outflow_beta[face_loc];
                            const Real& length_char = d_bdry_face_nonreflecting_outflow_length_char[face_loc];
                            
                            const Real T_1 = u_y_B*(dp_dx - rho_y_B*c_y_B*dv_dx) + rho_y_B*c_y_B*c_y_B*du_dx +
                                w_y_B*(dp_dz - rho_y_B*c_y_B*dv_dz) + rho_y_B*c_y_B*c_y_B*dw_dz;
                            
                            const Real M_sq = (u_y_B*u_y_B + v_y_B*v_y_B + w_y_B*w_y_B)/(c_y_B*c_y_B);
                            const Real K = sigma*c_y_B*(Real(1) - M_sq)/length_char;
                            
                            Lambda_inv_L[0] = (Real(1)/lambda_1)*(K*(p_y_B - p_t) - (Real(1) - beta)*T_1);
                            Lambda_inv_L[1] = du_dy;
                            Lambda_inv_L[2] = c_y_B*c_y_B*drho_dy - dp_dy;
                            Lambda_inv_L[3] = dw_dy;
                            Lambda_inv_L[4] = dp_dy + rho_y_B*c_y_B*dv_dy;
                            
                            // Compute dV_dy.
                            
                            const Real c_sq_inv  = Real(1)/(c_y_B*c_y_B);
                            const Real rho_c_inv = Real(1)/(rho_y_B*c_y_B);
                            
                            Real dV_dy[5];
                            
                            dV_dy[0] = half*c_sq_inv*(Lambda_inv_L[0] + Lambda_inv_L[4]) + c_sq_inv*Lambda_inv_L[2];
                            dV_dy[1] = Lambda_inv_L[1];
                            dV_dy[2] = half*rho_c_inv*(-Lambda_inv_L[0] + Lambda_inv_L[4]);
                            dV_dy[3] = Lambda_inv_L[3];
                            dV_dy[4] = half*(Lambda_inv_L[0] + Lambda_inv_L[4]);
                            
                            Real V_ghost[5*num_ghosts_to_fill];
                            
                            for (int j = 0; j < num_ghosts_to_fill; j++)
                            {
                                const int idx_cell_rho = (i + num_subghosts_conservative_var[0][0]) +
                                    (j + fill_box_lo_idx[1] + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (k + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom = (i + num_subghosts_conservative_var[1][0]) +
                                    (j + fill_box_lo_idx[1] + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                    (k + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E = (i + num_subghosts_conservative_var[2][0]) +
                                    (j + fill_box_lo_idx[1] + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                    (k + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                if (j == 0)
                                {
                                    V_ghost[j*5 + 0] = rho_y_BB + Real(2)*Real(dx[1])*dV_dy[0];
                                    V_ghost[j*5 + 1] = u_y_BB   + Real(2)*Real(dx[1])*dV_dy[1];
                                    V_ghost[j*5 + 2] = v_y_BB   + Real(2)*Real(dx[1])*dV_dy[2];
                                    V_ghost[j*5 + 3] = w_y_BB   + Real(2)*Real(dx[1])*dV_dy[3];
                                    V_ghost[j*5 + 4] = p_y_BB   + Real(2)*Real(dx[1])*dV_dy[4];
                                }
                                else if (j == 1)
                                {
                                    V_ghost[j*5 + 0] = -Real(2)*rho_y_BB - Real(3)*rho_y_B +
                                        Real(6)*V_ghost[(j - 1)*5 + 0] - Real(6)*Real(dx[1])*dV_dy[0];
                                    
                                    V_ghost[j*5 + 1] = -Real(2)*u_y_BB - Real(3)*u_y_B +
                                        Real(6)*V_ghost[(j - 1)*5 + 1] - Real(6)*Real(dx[1])*dV_dy[1];
                                    
                                    V_ghost[j*5 + 2] = -Real(2)*v_y_BB - Real(3)*v_y_B +
                                        Real(6)*V_ghost[(j - 1)*5 + 2] - Real(6)*Real(dx[1])*dV_dy[2];
                                    
                                    V_ghost[j*5 + 3] = -Real(2)*w_y_BB - Real(3)*w_y_B +
                                        Real(6)*V_ghost[(j - 1)*5 + 3] - Real(6)*Real(dx[1])*dV_dy[3];
                                    
                                    V_ghost[j*5 + 4] = -Real(2)*p_y_BB - Real(3)*p_y_B +
                                        Real(6)*V_ghost[(j - 1)*5 + 4] - Real(6)*Real(dx[1])*dV_dy[4];
                                }
                                else if (j == 2)
                                {
                                    V_ghost[j*5 + 0] = Real(3)*rho_y_BB + Real(10)*rho_y_B -
                                        Real(18)*V_ghost[(j - 2)*5 + 0] + Real(6)*V_ghost[(j - 1)*5 + 0] +
                                        Real(12)*Real(dx[1])*dV_dy[0];
                                    
                                    V_ghost[j*5 + 1] = Real(3)*u_y_BB + Real(10)*u_y_B -
                                        Real(18)*V_ghost[(j - 2)*5 + 1] + Real(6)*V_ghost[(j - 1)*5 + 1] +
                                        Real(12)*Real(dx[1])*dV_dy[1];
                                    
                                    V_ghost[j*5 + 2] = Real(3)*v_y_BB + Real(10)*v_y_B -
                                        Real(18)*V_ghost[(j - 2)*5 + 2] + Real(6)*V_ghost[(j - 1)*5 + 2] +
                                        Real(12)*Real(dx[1])*dV_dy[2];
                                    
                                    V_ghost[j*5 + 3] = Real(3)*w_y_BB + Real(10)*w_y_B -
                                        Real(18)*V_ghost[(j - 2)*5 + 3] + Real(6)*V_ghost[(j - 1)*5 + 3] +
                                        Real(12)*Real(dx[1])*dV_dy[3];
                                    
                                    V_ghost[j*5 + 4] = Real(3)*p_y_BB + Real(10)*p_y_B -
                                        Real(18)*V_ghost[(j - 2)*5 + 4] + Real(6)*V_ghost[(j - 1)*5 + 4] +
                                        Real(12)*Real(dx[1])*dV_dy[4];
                                }
                                else if (j == 3)
                                {
                                    V_ghost[j*5 + 0] = -Real(4)*rho_y_BB - Real(65)/Real(3)*rho_y_B +
                                        Real(40)*V_ghost[(j - 3)*5 + 0] - Real(20)*V_ghost[(j - 2)*5 + 0] +
                                        Real(20)/Real(3)*V_ghost[(j - 1)*5 + 0] - Real(20)*Real(dx[1])*dV_dy[0];
                                    
                                    V_ghost[j*5 + 1] = -Real(4)*u_y_BB - Real(65)/Real(3)*u_y_B +
                                        Real(40)*V_ghost[(j - 3)*5 + 1] - Real(20)*V_ghost[(j - 2)*5 + 1] +
                                        Real(20)/Real(3)*V_ghost[(j - 1)*5 + 1] - Real(20)*Real(dx[1])*dV_dy[1];
                                    
                                    V_ghost[j*5 + 2] = -Real(4)*v_y_BB - Real(65)/Real(3)*v_y_B +
                                        Real(40)*V_ghost[(j - 3)*5 + 2] - Real(20)*V_ghost[(j - 2)*5 + 2] +
                                        Real(20)/Real(3)*V_ghost[(j - 1)*5 + 2] - Real(20)*Real(dx[1])*dV_dy[2];
                                    
                                    V_ghost[j*5 + 3] = -Real(4)*w_y_BB - Real(65)/Real(3)*w_y_B +
                                        Real(40)*V_ghost[(j - 3)*5 + 3] - Real(20)*V_ghost[(j - 2)*5 + 3] +
                                        Real(20)/Real(3)*V_ghost[(j - 1)*5 + 3] - Real(20)*Real(dx[1])*dV_dy[3];
                                    
                                    V_ghost[j*5 + 4] = -Real(4)*p_y_BB - Real(65)/Real(3)*p_y_B +
                                        Real(40)*V_ghost[(j - 3)*5 + 4] - Real(20)*V_ghost[(j - 2)*5 + 4] +
                                        Real(20)/Real(3)*V_ghost[(j - 1)*5 + 4] - Real(20)*Real(dx[1])*dV_dy[4];
                                }
                                else if (j == 4)
                                {
                                    V_ghost[j*5 + 0] = Real(5)*rho_y_BB + Real(77)/Real(2)*rho_y_B -
                                        Real(75)*V_ghost[(j - 4)*5 + 0] + Real(50)*V_ghost[(j - 3)*5 + 0] -
                                        Real(25)*V_ghost[(j - 2)*5 + 0] + Real(15)/Real(2)*V_ghost[(j - 1)*5 + 0] +
                                        Real(30)*Real(dx[1])*dV_dy[0];
                                    
                                    V_ghost[j*5 + 1] = Real(5)*u_y_BB + Real(77)/Real(2)*u_y_B -
                                        Real(75)*V_ghost[(j - 4)*5 + 1] + Real(50)*V_ghost[(j - 3)*5 + 1] -
                                        Real(25)*V_ghost[(j - 2)*5 + 1] + Real(15)/Real(2)*V_ghost[(j - 1)*5 + 1] +
                                        Real(30)*Real(dx[1])*dV_dy[1];
                                    
                                    V_ghost[j*5 + 2] = Real(5)*v_y_BB + Real(77)/Real(2)*v_y_B -
                                        Real(75)*V_ghost[(j - 4)*5 + 2] + Real(50)*V_ghost[(j - 3)*5 + 2] -
                                        Real(25)*V_ghost[(j - 2)*5 + 2] + Real(15)/Real(2)*V_ghost[(j - 1)*5 + 2] +
                                        Real(30)*Real(dx[1])*dV_dy[2];
                                    
                                    V_ghost[j*5 + 3] = Real(5)*w_y_BB + Real(77)/Real(2)*w_y_B -
                                        Real(75)*V_ghost[(j - 4)*5 + 3] + Real(50)*V_ghost[(j - 3)*5 + 3] -
                                        Real(25)*V_ghost[(j - 2)*5 + 3] + Real(15)/Real(2)*V_ghost[(j - 1)*5 + 3] +
                                        Real(30)*Real(dx[1])*dV_dy[3];
                                    
                                    V_ghost[j*5 + 4] = Real(5)*p_y_BB + Real(77)/Real(2)*p_y_B -
                                        Real(75)*V_ghost[(j - 4)*5 + 4] + Real(50)*V_ghost[(j - 3)*5 + 4] -
                                        Real(25)*V_ghost[(j - 2)*5 + 4] + Real(15)/Real(2)*V_ghost[(j - 1)*5 + 4] +
                                        Real(30)*Real(dx[1])*dV_dy[4];
                                }
                                else if (j == 5)
                                {
                                    V_ghost[j*5 + 0] = -Real(6)*rho_y_BB - Real(609)/Real(10)*rho_y_B +
                                        Real(126)*V_ghost[(j - 5)*5 + 0] - Real(105)*V_ghost[(j - 4)*5 + 0] +
                                        Real(70)*V_ghost[(j - 3)*5 + 0] - Real(63)/Real(2)*V_ghost[(j - 2)*5 + 0] +
                                        Real(42)/Real(5)*V_ghost[(j - 1)*5 + 0] - Real(42)*Real(dx[1])*dV_dy[0];
                                    
                                    V_ghost[j*5 + 1] = -Real(6)*u_y_BB - Real(609)/Real(10)*u_y_B +
                                        Real(126)*V_ghost[(j - 5)*5 + 1] - Real(105)*V_ghost[(j - 4)*5 + 1] +
                                        Real(70)*V_ghost[(j - 3)*5 + 1] - Real(63)/Real(2)*V_ghost[(j - 2)*5 + 1] +
                                        Real(42)/Real(5)*V_ghost[(j - 1)*5 + 1] - Real(42)*Real(dx[1])*dV_dy[1];
                                    
                                    V_ghost[j*5 + 2] = -Real(6)*v_y_BB - Real(609)/Real(10)*v_y_B +
                                        Real(126)*V_ghost[(j - 5)*5 + 2] - Real(105)*V_ghost[(j - 4)*5 + 2] +
                                        Real(70)*V_ghost[(j - 3)*5 + 2] - Real(63)/Real(2)*V_ghost[(j - 2)*5 + 2] +
                                        Real(42)/Real(5)*V_ghost[(j - 1)*5 + 2] - Real(42)*Real(dx[1])*dV_dy[2];
                                    
                                    V_ghost[j*5 + 3] = -Real(6)*w_y_BB - Real(609)/Real(10)*w_y_B +
                                        Real(126)*V_ghost[(j - 5)*5 + 3] - Real(105)*V_ghost[(j - 4)*5 + 3] +
                                        Real(70)*V_ghost[(j - 3)*5 + 3] - Real(63)/Real(2)*V_ghost[(j - 2)*5 + 3] +
                                        Real(42)/Real(5)*V_ghost[(j - 1)*5 + 3] - Real(42)*Real(dx[1])*dV_dy[3];
                                    
                                    V_ghost[j*5 + 4] = -Real(6)*p_y_BB - Real(609)/Real(10)*p_y_B +
                                        Real(126)*V_ghost[(j - 5)*5 + 4] - Real(105)*V_ghost[(j - 4)*5 + 4] +
                                        Real(70)*V_ghost[(j - 3)*5 + 4] - Real(63)/Real(2)*V_ghost[(j - 2)*5 + 4] +
                                        Real(42)/Real(5)*V_ghost[(j - 1)*5 + 4] - Real(42)*Real(dx[1])*dV_dy[4];
                                }
                                
                                Q[0][idx_cell_rho] = V_ghost[j*5 + 0];
                                Q[1][idx_cell_mom] = V_ghost[j*5 + 0]*V_ghost[j*5 + 1];
                                Q[2][idx_cell_mom] = V_ghost[j*5 + 0]*V_ghost[j*5 + 2];
                                Q[3][idx_cell_mom] = V_ghost[j*5 + 0]*V_ghost[j*5 + 3];
                                
                                const Real epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getInternalEnergy(
                                        &V_ghost[j*5 + 0],
                                        &V_ghost[j*5 + 4],
                                        thermo_properties_ptr);
                                
                                const Real E = V_ghost[j*5 + 0]*epsilon +
                                    half*(Q[1][idx_cell_mom]*Q[1][idx_cell_mom] + Q[2][idx_cell_mom]*Q[2][idx_cell_mom] +
                                       Q[3][idx_cell_mom]*Q[3][idx_cell_mom])/V_ghost[j*5 + 0];
                                
                                Q[4][idx_cell_E] = E;
                            }
                        }
                    }
                }
                else if (face_loc == BDRY_LOC::ZLO)
                {
                    const int num_ghosts_to_fill = fill_box_hi_idx[2] - fill_box_lo_idx[2] + 1;
                    TBOX_ASSERT(fill_box_hi_idx[2] == interior_box_lo_idx[2] - 1);
                    if (num_ghosts_to_fill > 6)
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelBoundaryUtilitiesSingleSpecies::fill3dFaceBoundaryData()\n"
                            << "Non-reflecting outflow BC doesn't support more than six ghost cells yet!");
                    }
                    
                    for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                    {
                        for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                        {
                            // Get the grid spacing.
                            const double* const dx = patch_geom->getDx();
                            
                            const int idx_cell_rho_z_F = (i + num_subghosts_conservative_var[0][0]) +
                                (j + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                (interior_box_lo_idx[2] + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][1]*
                                    subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_rho_z_FF = (i + num_subghosts_conservative_var[0][0]) +
                                (j + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                (interior_box_lo_idx[2] + 1 + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                    subghostcell_dims_conservative_var[0][1];
                            
                            const int idx_cell_rho_z_FFF = (i + num_subghosts_conservative_var[0][0]) +
                                (j + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                (interior_box_lo_idx[2] + 2 + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                    subghostcell_dims_conservative_var[0][1];
                            
                            const int idx_cell_mom_z_F = (i + num_subghosts_conservative_var[1][0]) +
                                (j + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                (interior_box_lo_idx[2] + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                    subghostcell_dims_conservative_var[1][1];
                            
                            const int idx_cell_mom_z_FF = (i + num_subghosts_conservative_var[1][0]) +
                                (j + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                (interior_box_lo_idx[2] + 1 + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                    subghostcell_dims_conservative_var[1][1];
                            
                            const int idx_cell_mom_z_FFF = (i + num_subghosts_conservative_var[1][0]) +
                                (j + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                (interior_box_lo_idx[2] + 2 + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                    subghostcell_dims_conservative_var[1][1];
                            
                            const int idx_cell_E_z_F = (i + num_subghosts_conservative_var[2][0]) +
                                (j + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                (interior_box_lo_idx[2] + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                    subghostcell_dims_conservative_var[2][1];
                            
                            const int idx_cell_E_z_FF = (i + num_subghosts_conservative_var[2][0]) +
                                (j + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                (interior_box_lo_idx[2] + 1 + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                    subghostcell_dims_conservative_var[2][1];
                            
                            const int idx_cell_E_z_FFF = (i + num_subghosts_conservative_var[2][0]) +
                                (j + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                (interior_box_lo_idx[2] + 2 + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                    subghostcell_dims_conservative_var[2][1];
                            
                            const Real& rho_z_F   = Q[0][idx_cell_rho_z_F];
                            const Real& rho_z_FF  = Q[0][idx_cell_rho_z_FF];
                            const Real& rho_z_FFF = Q[0][idx_cell_rho_z_FFF];
                            
                            const Real u_z_F   = Q[1][idx_cell_mom_z_F]/rho_z_F;
                            const Real u_z_FF  = Q[1][idx_cell_mom_z_FF]/rho_z_FF;
                            const Real u_z_FFF = Q[1][idx_cell_mom_z_FFF]/rho_z_FFF;
                            
                            const Real v_z_F   = Q[2][idx_cell_mom_z_F]/rho_z_F;
                            const Real v_z_FF  = Q[2][idx_cell_mom_z_FF]/rho_z_FF;
                            const Real v_z_FFF = Q[2][idx_cell_mom_z_FFF]/rho_z_FFF;
                            
                            const Real w_z_F   = Q[3][idx_cell_mom_z_F]/rho_z_F;
                            const Real w_z_FF  = Q[3][idx_cell_mom_z_FF]/rho_z_FF;
                            const Real w_z_FFF = Q[3][idx_cell_mom_z_FFF]/rho_z_FFF;
                            
                            const Real half = Real(1)/Real(2);
                            const Real epsilon_z_F   = Q[4][idx_cell_E_z_F]/rho_z_F - half*(u_z_F*u_z_F + v_z_F*v_z_F + w_z_F*w_z_F);
                            const Real epsilon_z_FF  = Q[4][idx_cell_E_z_FF]/rho_z_FF - half*(u_z_FF*u_z_FF + v_z_FF*v_z_FF + w_z_FF*w_z_FF);
                            const Real epsilon_z_FFF = Q[4][idx_cell_E_z_FFF]/rho_z_FFF - half*(u_z_FFF*u_z_FFF + v_z_FFF*v_z_FFF + w_z_FFF*w_z_FFF);
                            
                            const Real p_z_F = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getPressure(
                                    &rho_z_F,
                                    &epsilon_z_F,
                                    thermo_properties_ptr);
                            
                            const Real p_z_FF = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getPressure(
                                    &rho_z_FF,
                                    &epsilon_z_FF,
                                    thermo_properties_ptr);
                            
                            const Real p_z_FFF = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getPressure(
                                    &rho_z_FFF,
                                    &epsilon_z_FFF,
                                    thermo_properties_ptr);
                            
                            /*
                             * Compute derivatives at z-direction.
                             */
                            
                            const Real drho_dz = -(rho_z_FFF - Real(4)*rho_z_FF + Real(3)*rho_z_F)/(Real(2)*Real(dx[2]));
                            const Real du_dz   = -(u_z_FFF - Real(4)*u_z_FF + Real(3)*u_z_F)/(Real(2)*Real(dx[2]));
                            const Real dv_dz   = -(v_z_FFF - Real(4)*v_z_FF + Real(3)*v_z_F)/(Real(2)*Real(dx[2]));
                            const Real dw_dz   = -(w_z_FFF - Real(4)*w_z_FF + Real(3)*w_z_F)/(Real(2)*Real(dx[2]));
                            const Real dp_dz   = -(p_z_FFF - Real(4)*p_z_FF + Real(3)*p_z_F)/(Real(2)*Real(dx[2]));
                            
                            /*
                             * Compute derivatives in x-direction.
                             */
                            
                            Real du_dx = Real(0);
                            // Real dv_dx = Real(0);
                            Real dw_dx = Real(0);
                            Real dp_dx = Real(0);
                            
                            if ((i + num_subghosts_conservative_var[0][0] == 0) ||
                                (i + num_subghosts_conservative_var[1][0] == 0) ||
                                (i + num_subghosts_conservative_var[2][0] == 0))
                            // if (((patch_geom->getTouchesRegularBoundary(0, 0)) && (i == interior_box_lo_idx[0])) ||
                            //     ((i + num_subghosts_conservative_var[0][0] == 0) ||
                            //      (i + num_subghosts_conservative_var[1][0] == 0) ||
                            //      (i + num_subghosts_conservative_var[2][0] == 0)))
                            {
                                // Patch is touching periodic boundary.
                                
                                const int idx_cell_rho_x_R = (i + 1 + num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (interior_box_lo_idx[2] + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom_x_R = (i + 1 + num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                    (interior_box_lo_idx[2] + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E_x_R = (i + 1 + num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                    (interior_box_lo_idx[2] + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                const Real& rho_x_R = Q[0][idx_cell_rho_x_R];
                                const Real u_x_R = Q[1][idx_cell_mom_x_R]/rho_x_R;
                                const Real v_x_R = Q[2][idx_cell_mom_x_R]/rho_x_R;
                                const Real w_x_R = Q[3][idx_cell_mom_x_R]/rho_x_R;
                                const Real epsilon_x_R = Q[4][idx_cell_E_x_R]/rho_x_R - half*(u_x_R*u_x_R + v_x_R*v_x_R + w_x_R*w_x_R);
                                
                                const Real p_x_R = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &rho_x_R,
                                        &epsilon_x_R,
                                        thermo_properties_ptr);
                                
                                // One-sided derivatives.
                                du_dx = (u_x_R - u_z_F)/Real(dx[0]);
                                // dv_dx = (v_x_R - v_z_F)/Real(dx[0]);
                                dw_dx = (w_x_R - w_z_F)/Real(dx[0]);
                                dp_dx = (p_x_R - p_z_F)/Real(dx[0]);
                            }
                            else if ((i + num_subghosts_conservative_var[0][0] + 1 == subghostcell_dims_conservative_var[0][0]) ||
                                     (i + num_subghosts_conservative_var[1][0] + 1 == subghostcell_dims_conservative_var[1][0]) ||
                                     (i + num_subghosts_conservative_var[2][0] + 1 == subghostcell_dims_conservative_var[2][0]))
                            // else if (((patch_geom->getTouchesRegularBoundary(0, 1)) && (i == interior_box_hi_idx[0])) ||
                            //          ((i + num_subghosts_conservative_var[0][0] + 1 == subghostcell_dims_conservative_var[0][0]) ||
                            //           (i + num_subghosts_conservative_var[1][0] + 1 == subghostcell_dims_conservative_var[1][0]) ||
                            //           (i + num_subghosts_conservative_var[2][0] + 1 == subghostcell_dims_conservative_var[2][0])))
                            {
                                // Patch is touching periodic boundary.
                                
                                const int idx_cell_rho_x_L = (i - 1 + num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (interior_box_lo_idx[2] + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom_x_L = (i - 1 + num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                    (interior_box_lo_idx[2] + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E_x_L = (i - 1 + num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                    (interior_box_lo_idx[2] + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                const Real& rho_x_L = Q[0][idx_cell_rho_x_L];
                                const Real u_x_L = Q[1][idx_cell_mom_x_L]/rho_x_L;
                                const Real v_x_L = Q[2][idx_cell_mom_x_L]/rho_x_L;
                                const Real w_x_L = Q[3][idx_cell_mom_x_L]/rho_x_L;
                                const Real epsilon_x_L = Q[4][idx_cell_E_x_L]/rho_x_L - half*(u_x_L*u_x_L + v_x_L*v_x_L + w_x_L*w_x_L);
                                
                                const Real p_x_L = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &rho_x_L,
                                        &epsilon_x_L,
                                        thermo_properties_ptr);
                                
                                // One-sided derivatives.
                                du_dx = (u_z_F - u_x_L)/Real(dx[0]);
                                // dv_dx = (v_z_F - v_x_L)/Real(dx[0]);
                                dw_dx = (w_z_F - w_x_L)/Real(dx[0]);
                                dp_dx = (p_z_F - p_x_L)/Real(dx[0]);
                            }
                            else
                            {
                                const int idx_cell_rho_x_L = (i - 1 + num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (interior_box_lo_idx[2] + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_rho_x_R = (i + 1 + num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (interior_box_lo_idx[2] + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom_x_L = (i - 1 + num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                    (interior_box_lo_idx[2] + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_mom_x_R = (i + 1 + num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                    (interior_box_lo_idx[2] + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E_x_L = (i - 1 + num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                    (interior_box_lo_idx[2] + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                const int idx_cell_E_x_R = (i + 1 + num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                    (interior_box_lo_idx[2] + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                const Real& rho_x_L = Q[0][idx_cell_rho_x_L];
                                const Real& rho_x_R = Q[0][idx_cell_rho_x_R];
                                
                                const Real u_x_L = Q[1][idx_cell_mom_x_L]/rho_x_L;
                                const Real u_x_R = Q[1][idx_cell_mom_x_R]/rho_x_R;
                                
                                const Real v_x_L = Q[2][idx_cell_mom_x_L]/rho_x_L;
                                const Real v_x_R = Q[2][idx_cell_mom_x_R]/rho_x_R;
                                
                                const Real w_x_L = Q[3][idx_cell_mom_x_L]/rho_x_L;
                                const Real w_x_R = Q[3][idx_cell_mom_x_R]/rho_x_R;
                                
                                const Real epsilon_x_L = Q[4][idx_cell_E_x_L]/rho_x_L - half*(u_x_L*u_x_L + v_x_L*v_x_L + w_x_L*w_x_L);
                                const Real epsilon_x_R = Q[4][idx_cell_E_x_R]/rho_x_R - half*(u_x_R*u_x_R + v_x_R*v_x_R + w_x_R*w_x_R);
                                
                                const Real p_x_L = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &rho_x_L,
                                        &epsilon_x_L,
                                        thermo_properties_ptr);
                                
                                const Real p_x_R = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &rho_x_R,
                                        &epsilon_x_R,
                                        thermo_properties_ptr);
                                
                                // Central derivatives.
                                du_dx = (u_x_R - u_x_L)/(Real(2)*Real(dx[0]));
                                // dv_dx = (v_x_R - v_x_L)/(Real(2)*Real(dx[0]));
                                dw_dx = (w_x_R - w_x_L)/(Real(2)*Real(dx[0]));
                                dp_dx = (p_x_R - p_x_L)/(Real(2)*Real(dx[0]));
                            }
                            
                            /*
                             * Compute derivatives in y-direction.
                             */
                            
                            // Real du_dy = Real(0);
                            Real dv_dy = Real(0);
                            Real dw_dy = Real(0);
                            Real dp_dy = Real(0);
                            
                            if ((j + num_subghosts_conservative_var[0][1] == 0) ||
                                (j + num_subghosts_conservative_var[1][1] == 0) ||
                                (j + num_subghosts_conservative_var[2][1] == 0))
                            // if (((patch_geom->getTouchesRegularBoundary(1, 0)) && (j == interior_box_lo_idx[1])) ||
                            //     ((j + num_subghosts_conservative_var[0][1] == 0) ||
                            //      (j + num_subghosts_conservative_var[1][1] == 0) ||
                            //      (j + num_subghosts_conservative_var[2][1] == 0)))
                            {
                                // Patch is touching periodic boundary.
                                
                                const int idx_cell_rho_y_T = (i + num_subghosts_conservative_var[0][0]) +
                                    (j + 1 + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (interior_box_lo_idx[2] + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom_y_T = (i + num_subghosts_conservative_var[1][0]) +
                                    (j + 1 + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                    (interior_box_lo_idx[2] + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E_y_T = (i + num_subghosts_conservative_var[2][0]) +
                                    (j + 1 + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                    (interior_box_lo_idx[2] + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                const Real& rho_y_T = Q[0][idx_cell_rho_y_T];
                                const Real u_y_T = Q[1][idx_cell_mom_y_T]/rho_y_T;
                                const Real v_y_T = Q[2][idx_cell_mom_y_T]/rho_y_T;
                                const Real w_y_T = Q[3][idx_cell_mom_y_T]/rho_y_T;
                                const Real epsilon_y_T = Q[4][idx_cell_E_y_T]/rho_y_T - half*(u_y_T*u_y_T + v_y_T*v_y_T + w_y_T*w_y_T);
                                
                                const Real p_y_T = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &rho_y_T,
                                        &epsilon_y_T,
                                        thermo_properties_ptr);
                                
                                // One-sided derivatives.
                                //du_dy = (u_y_T - u_z_F)/Real(dx[1]);
                                dv_dy = (v_y_T - v_z_F)/Real(dx[1]);
                                dw_dy = (w_y_T - w_z_F)/Real(dx[1]);
                                dp_dy = (p_y_T - p_z_F)/Real(dx[1]);
                            }
                            else if ((j + num_subghosts_conservative_var[0][1] + 1 == subghostcell_dims_conservative_var[0][1]) ||
                                     (j + num_subghosts_conservative_var[1][1] + 1 == subghostcell_dims_conservative_var[1][1]) ||
                                     (j + num_subghosts_conservative_var[2][1] + 1 == subghostcell_dims_conservative_var[2][1]))
                            // else if (((patch_geom->getTouchesRegularBoundary(1, 1)) && (j == interior_box_hi_idx[1])) ||
                            //          ((j + num_subghosts_conservative_var[0][1] + 1 == subghostcell_dims_conservative_var[0][1]) ||
                            //           (j + num_subghosts_conservative_var[1][1] + 1 == subghostcell_dims_conservative_var[1][1]) ||
                            //           (j + num_subghosts_conservative_var[2][1] + 1 == subghostcell_dims_conservative_var[2][1])))
                            {
                                // Patch is touching periodic boundary.
                                
                                const int idx_cell_rho_y_B = (i + num_subghosts_conservative_var[0][0]) +
                                    (j - 1 + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (interior_box_lo_idx[2] + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom_y_B = (i + num_subghosts_conservative_var[1][0]) +
                                    (j - 1 + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                    (interior_box_lo_idx[2] + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E_y_B = (i + num_subghosts_conservative_var[2][0]) +
                                    (j - 1 + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                    (interior_box_lo_idx[2] + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                const Real& rho_y_B = Q[0][idx_cell_rho_y_B];
                                const Real u_y_B = Q[1][idx_cell_mom_y_B]/rho_y_B;
                                const Real v_y_B = Q[2][idx_cell_mom_y_B]/rho_y_B;
                                const Real w_y_B = Q[3][idx_cell_mom_y_B]/rho_y_B;
                                const Real epsilon_y_B = Q[4][idx_cell_E_y_B]/rho_y_B - half*(u_y_B*u_y_B + v_y_B*v_y_B + w_y_B*w_y_B);
                                
                                const Real p_y_B = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &rho_y_B,
                                        &epsilon_y_B,
                                        thermo_properties_ptr);
                                
                                // One-sided derivatives.
                                //du_dy = (u_z_F - u_y_B)/Real(dx[1]);
                                dv_dy = (v_z_F - v_y_B)/Real(dx[1]);
                                dw_dy = (w_z_F - w_y_B)/Real(dx[1]);
                                dp_dy = (p_z_F - p_y_B)/Real(dx[1]);
                            }
                            else
                            {
                                const int idx_cell_rho_y_B = (i + num_subghosts_conservative_var[0][0]) +
                                    (j - 1 + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (interior_box_lo_idx[2] +  num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_rho_y_T = (i + num_subghosts_conservative_var[0][0]) +
                                    (j + 1 + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (interior_box_lo_idx[2] + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom_y_B = (i + num_subghosts_conservative_var[1][0]) +
                                    (j - 1 + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                    (interior_box_lo_idx[2] + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_mom_y_T = (i + num_subghosts_conservative_var[1][0]) +
                                    (j + 1 + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                    (interior_box_lo_idx[2] + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E_y_B = (i + num_subghosts_conservative_var[2][0]) +
                                    (j - 1 + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                    (interior_box_lo_idx[2] + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                const int idx_cell_E_y_T = (i + num_subghosts_conservative_var[2][0]) +
                                    (j + 1 + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                    (interior_box_lo_idx[2] + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                const Real& rho_y_B = Q[0][idx_cell_rho_y_B];
                                const Real& rho_y_T = Q[0][idx_cell_rho_y_T];
                                
                                const Real u_y_B = Q[1][idx_cell_mom_y_B]/rho_y_B;
                                const Real u_y_T = Q[1][idx_cell_mom_y_T]/rho_y_T;
                                
                                const Real v_y_B = Q[2][idx_cell_mom_y_B]/rho_y_B;
                                const Real v_y_T = Q[2][idx_cell_mom_y_T]/rho_y_T;
                                
                                const Real w_y_B = Q[3][idx_cell_mom_y_B]/rho_y_B;
                                const Real w_y_T = Q[3][idx_cell_mom_y_T]/rho_y_T;
                                
                                const Real epsilon_y_B = Q[4][idx_cell_E_y_B]/rho_y_B - half*(u_y_B*u_y_B + v_y_B*v_y_B + w_y_B*w_y_B);
                                const Real epsilon_y_T = Q[4][idx_cell_E_y_T]/rho_y_T - half*(u_y_T*u_y_T + v_y_T*v_y_T + w_y_T*w_y_T);
                                
                                const Real p_y_B = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &rho_y_B,
                                        &epsilon_y_B,
                                        thermo_properties_ptr);
                                
                                const Real p_y_T = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &rho_y_T,
                                        &epsilon_y_T,
                                        thermo_properties_ptr);
                                
                                // Central derivatives.
                                // du_dy = (u_y_T - u_y_B)/(Real(2)*Real(dx[1]));
                                dv_dy = (v_y_T - v_y_B)/(Real(2)*Real(dx[1]));
                                dw_dy = (w_y_T - w_y_B)/(Real(2)*Real(dx[1]));
                                dp_dy = (p_y_T - p_y_B)/(Real(2)*Real(dx[1]));
                            }
                            
                            const Real c_z_F = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getSoundSpeed(
                                    &rho_z_F,
                                    &p_z_F,
                                    thermo_properties_ptr);
                            
                            const Real lambda_5 = w_z_F + c_z_F;
                            
                            // Compute vector Lambda^(-1) * L.
                            
                            Real Lambda_inv_L[5];
                            
                            const Real& p_t         = d_bdry_face_nonreflecting_outflow_p_t[face_loc];
                            const Real& sigma       = d_bdry_face_nonreflecting_outflow_sigma[face_loc];
                            const Real& beta        = d_bdry_face_nonreflecting_outflow_beta[face_loc];
                            const Real& length_char = d_bdry_face_nonreflecting_outflow_length_char[face_loc];
                            
                            const Real T_5 = u_z_F*(dp_dx + rho_z_F*c_z_F*dw_dx) + rho_z_F*c_z_F*c_z_F*du_dx +
                                v_z_F*(dp_dy + rho_z_F*c_z_F*dw_dy) + rho_z_F*c_z_F*c_z_F*dv_dy;
                            
                            const Real M_sq = (u_z_F*u_z_F + v_z_F*v_z_F + w_z_F*w_z_F)/(c_z_F*c_z_F);
                            const Real K = sigma*c_z_F*(Real(1) - M_sq)/length_char;
                            
                            Lambda_inv_L[0] = dp_dz - rho_z_F*c_z_F*dw_dz;
                            Lambda_inv_L[1] = du_dz;
                            Lambda_inv_L[2] = dv_dz;
                            Lambda_inv_L[3] = c_z_F*c_z_F*drho_dz - dp_dz;
                            Lambda_inv_L[4] = (Real(1)/lambda_5)*(K*(p_z_F - p_t) - (Real(1) - beta)*T_5);
                            
                            // Compute dV_dz.
                            
                            const Real c_sq_inv  = Real(1)/(c_z_F*c_z_F);
                            const Real rho_c_inv = Real(1)/(rho_z_F*c_z_F);
                            
                            Real dV_dz[5];
                            
                            dV_dz[0] = half*c_sq_inv*(Lambda_inv_L[0] + Lambda_inv_L[4]) + c_sq_inv*Lambda_inv_L[3];
                            dV_dz[1] = Lambda_inv_L[1];
                            dV_dz[2] = Lambda_inv_L[2];
                            dV_dz[3] = half*rho_c_inv*(-Lambda_inv_L[0] + Lambda_inv_L[4]);
                            dV_dz[4] = half*(Lambda_inv_L[0] + Lambda_inv_L[4]);
                            
                            Real V_ghost[5*num_ghosts_to_fill];
                            
                            for (int k = num_ghosts_to_fill - 1; k >= 0; k--)
                            {
                                const int idx_cell_rho = (i + num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (k + fill_box_lo_idx[2] + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom = (i + num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                    (k + fill_box_lo_idx[2] + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E = (i + num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                    (k + fill_box_lo_idx[2] + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                if (k == num_ghosts_to_fill - 1)
                                {
                                    V_ghost[k*5 + 0] = rho_z_FF - Real(2)*Real(dx[2])*dV_dz[0];
                                    V_ghost[k*5 + 1] = u_z_FF   - Real(2)*Real(dx[2])*dV_dz[1];
                                    V_ghost[k*5 + 2] = v_z_FF   - Real(2)*Real(dx[2])*dV_dz[2];
                                    V_ghost[k*5 + 3] = w_z_FF   - Real(2)*Real(dx[2])*dV_dz[3];
                                    V_ghost[k*5 + 4] = p_z_FF   - Real(2)*Real(dx[2])*dV_dz[4];
                                }
                                else if (k == num_ghosts_to_fill - 2)
                                {
                                    V_ghost[k*5 + 0] = -Real(2)*rho_z_FF - Real(3)*rho_z_F +
                                        Real(6)*V_ghost[(k + 1)*5 + 0] + Real(6)*Real(dx[2])*dV_dz[0];
                                    
                                    V_ghost[k*5 + 1] = -Real(2)*u_z_FF - Real(3)*u_z_F +
                                        Real(6)*V_ghost[(k + 1)*5 + 1] + Real(6)*Real(dx[2])*dV_dz[1];
                                    
                                    V_ghost[k*5 + 2] = -Real(2)*v_z_FF - Real(3)*v_z_F +
                                        Real(6)*V_ghost[(k + 1)*5 + 2] + Real(6)*Real(dx[2])*dV_dz[2];
                                    
                                    V_ghost[k*5 + 3] = -Real(2)*w_z_FF - Real(3)*w_z_F +
                                        Real(6)*V_ghost[(k + 1)*5 + 3] + Real(6)*Real(dx[2])*dV_dz[3];
                                    
                                    V_ghost[k*5 + 4] = -Real(2)*p_z_FF - Real(3)*p_z_F +
                                        Real(6)*V_ghost[(k + 1)*5 + 4] + Real(6)*Real(dx[2])*dV_dz[4];
                                }
                                else if (k == num_ghosts_to_fill - 3)
                                {
                                    V_ghost[k*5 + 0] = Real(3)*rho_z_FF + Real(10)*rho_z_F -
                                        Real(18)*V_ghost[(k + 2)*5 + 0] + Real(6)*V_ghost[(k + 1)*5 + 0] -
                                        Real(12)*Real(dx[2])*dV_dz[0];
                                    
                                    V_ghost[k*5 + 1] = Real(3)*u_z_FF + Real(10)*u_z_F -
                                        Real(18)*V_ghost[(k + 2)*5 + 1] + Real(6)*V_ghost[(k + 1)*5 + 1] -
                                        Real(12)*Real(dx[2])*dV_dz[1];
                                    
                                    V_ghost[k*5 + 2] = Real(3)*v_z_FF + Real(10)*v_z_F -
                                        Real(18)*V_ghost[(k + 2)*5 + 2] + Real(6)*V_ghost[(k + 1)*5 + 2] -
                                        Real(12)*Real(dx[2])*dV_dz[2];
                                    
                                    V_ghost[k*5 + 3] = Real(3)*w_z_FF + Real(10)*w_z_F -
                                        Real(18)*V_ghost[(k + 2)*5 + 3] + Real(6)*V_ghost[(k + 1)*5 + 3] -
                                        Real(12)*Real(dx[2])*dV_dz[3];
                                    
                                    V_ghost[k*5 + 4] = Real(3)*p_z_FF + Real(10)*p_z_F -
                                        Real(18)*V_ghost[(k + 2)*5 + 4] + Real(6)*V_ghost[(k + 1)*5 + 4] -
                                        Real(12)*Real(dx[2])*dV_dz[4];
                                }
                                else if (k == num_ghosts_to_fill - 4)
                                {
                                    V_ghost[k*5 + 0] = -Real(4)*rho_z_FF - Real(65)/Real(3)*rho_z_F +
                                        Real(40)*V_ghost[(k + 3)*5 + 0] - Real(20)*V_ghost[(k + 2)*5 + 0] +
                                        Real(20)/Real(3)*V_ghost[(k + 1)*5 + 0] + Real(20)*Real(dx[2])*dV_dz[0];
                                    
                                    V_ghost[k*5 + 1] = -Real(4)*u_z_FF - Real(65)/Real(3)*u_z_F +
                                        Real(40)*V_ghost[(k + 3)*5 + 1] - Real(20)*V_ghost[(k + 2)*5 + 1] +
                                        Real(20)/Real(3)*V_ghost[(k + 1)*5 + 1] + Real(20)*Real(dx[2])*dV_dz[1];
                                    
                                    V_ghost[k*5 + 2] = -Real(4)*v_z_FF - Real(65)/Real(3)*v_z_F +
                                        Real(40)*V_ghost[(k + 3)*5 + 2] - Real(20)*V_ghost[(k + 2)*5 + 2] +
                                        Real(20)/Real(3)*V_ghost[(k + 1)*5 + 2] + Real(20)*Real(dx[2])*dV_dz[2];
                                    
                                    V_ghost[k*5 + 3] = -Real(4)*w_z_FF - Real(65)/Real(3)*w_z_F +
                                        Real(40)*V_ghost[(k + 3)*5 + 3] - Real(20)*V_ghost[(k + 2)*5 + 3] +
                                        Real(20)/Real(3)*V_ghost[(k + 1)*5 + 3] + Real(20)*Real(dx[2])*dV_dz[3];
                                    
                                    V_ghost[k*5 + 4] = -Real(4)*p_z_FF - Real(65)/Real(3)*p_z_F +
                                        Real(40)*V_ghost[(k + 3)*5 + 4] - Real(20)*V_ghost[(k + 2)*5 + 4] +
                                        Real(20)/Real(3)*V_ghost[(k + 1)*5 + 4] + Real(20)*Real(dx[2])*dV_dz[4];
                                }
                                else if (k == num_ghosts_to_fill - 5)
                                {
                                    V_ghost[k*5 + 0] = Real(5)*rho_z_FF + Real(77)/Real(2)*rho_z_F -
                                        Real(75)*V_ghost[(k + 4)*5 + 0] + Real(50)*V_ghost[(k + 3)*5 + 0] -
                                        Real(25)*V_ghost[(k + 2)*5 + 0] + Real(15)/Real(2)*V_ghost[(k + 1)*5 + 0] -
                                        Real(30)*Real(dx[2])*dV_dz[0];
                                    
                                    V_ghost[k*5 + 1] = Real(5)*u_z_FF + Real(77)/Real(2)*u_z_F -
                                        Real(75)*V_ghost[(k + 4)*5 + 1] + Real(50)*V_ghost[(k + 3)*5 + 1] -
                                        Real(25)*V_ghost[(k + 2)*5 + 1] + Real(15)/Real(2)*V_ghost[(k + 1)*5 + 1] -
                                        Real(30)*Real(dx[2])*dV_dz[1];
                                    
                                    V_ghost[k*5 + 2] = Real(5)*v_z_FF + Real(77)/Real(2)*v_z_F -
                                        Real(75)*V_ghost[(k + 4)*5 + 2] + Real(50)*V_ghost[(k + 3)*5 + 2] -
                                        Real(25)*V_ghost[(k + 2)*5 + 2] + Real(15)/Real(2)*V_ghost[(k + 1)*5 + 2] -
                                        Real(30)*Real(dx[2])*dV_dz[2];
                                    
                                    V_ghost[k*5 + 3] = Real(5)*w_z_FF + Real(77)/Real(2)*w_z_F -
                                        Real(75)*V_ghost[(k + 4)*5 + 3] + Real(50)*V_ghost[(k + 3)*5 + 3] -
                                        Real(25)*V_ghost[(k + 2)*5 + 3] + Real(15)/Real(2)*V_ghost[(k + 1)*5 + 3] -
                                        Real(30)*Real(dx[2])*dV_dz[3];
                                    
                                    V_ghost[k*5 + 4] = Real(5)*p_z_FF + Real(77)/Real(2)*p_z_F -
                                        Real(75)*V_ghost[(k + 4)*5 + 4] + Real(50)*V_ghost[(k + 3)*5 + 4] -
                                        Real(25)*V_ghost[(k + 2)*5 + 4] + Real(15)/Real(2)*V_ghost[(k + 1)*5 + 4] -
                                        Real(30)*Real(dx[2])*dV_dz[4];
                                }
                                else if (k == num_ghosts_to_fill - 6)
                                {
                                    V_ghost[k*5 + 0] = -Real(6)*rho_z_FF - Real(609)/Real(10)*rho_z_F +
                                        Real(126)*V_ghost[(k + 5)*5 + 0] - Real(105)*V_ghost[(k + 4)*5 + 0] +
                                        Real(70)*V_ghost[(k + 3)*5 + 0] - Real(63)/Real(2)*V_ghost[(k + 2)*5 + 0] +
                                        Real(42)/Real(5)*V_ghost[(k + 1)*5 + 0] + Real(42)*Real(dx[2])*dV_dz[0];
                                    
                                    V_ghost[k*5 + 1] = -Real(6)*u_z_FF - Real(609)/Real(10)*u_z_F +
                                        Real(126)*V_ghost[(k + 5)*5 + 1] - Real(105)*V_ghost[(k + 4)*5 + 1] +
                                        Real(70)*V_ghost[(k + 3)*5 + 1] - Real(63)/Real(2)*V_ghost[(k + 2)*5 + 1] +
                                        Real(42)/Real(5)*V_ghost[(k + 1)*5 + 1] + Real(42)*Real(dx[2])*dV_dz[1];
                                    
                                    V_ghost[k*5 + 2] = -Real(6)*v_z_FF - Real(609)/Real(10)*v_z_F +
                                        Real(126)*V_ghost[(k + 5)*5 + 2] - Real(105)*V_ghost[(k + 4)*5 + 2] +
                                        Real(70)*V_ghost[(k + 3)*5 + 2] - Real(63)/Real(2)*V_ghost[(k + 2)*5 + 2] +
                                        Real(42)/Real(5)*V_ghost[(k + 1)*5 + 2] + Real(42)*Real(dx[2])*dV_dz[2];
                                    
                                    V_ghost[k*5 + 3] = -Real(6)*w_z_FF - Real(609)/Real(10)*w_z_F +
                                        Real(126)*V_ghost[(k + 5)*5 + 3] - Real(105)*V_ghost[(k + 4)*5 + 3] +
                                        Real(70)*V_ghost[(k + 3)*5 + 3] - Real(63)/Real(2)*V_ghost[(k + 2)*5 + 3] +
                                        Real(42)/Real(5)*V_ghost[(k + 1)*5 + 3] + Real(42)*Real(dx[2])*dV_dz[3];
                                    
                                    V_ghost[k*5 + 4] = -Real(6)*p_z_FF - Real(609)/Real(10)*p_z_F +
                                        Real(126)*V_ghost[(k + 5)*5 + 4] - Real(105)*V_ghost[(k + 4)*5 + 4] +
                                        Real(70)*V_ghost[(k + 3)*5 + 4] - Real(63)/Real(2)*V_ghost[(k + 2)*5 + 4] +
                                        Real(42)/Real(5)*V_ghost[(k + 1)*5 + 4] + Real(42)*Real(dx[2])*dV_dz[4];
                                }
                                
                                Q[0][idx_cell_rho] = V_ghost[k*5 + 0];
                                Q[1][idx_cell_mom] = V_ghost[k*5 + 0]*V_ghost[k*5 + 1];
                                Q[2][idx_cell_mom] = V_ghost[k*5 + 0]*V_ghost[k*5 + 2];
                                Q[3][idx_cell_mom] = V_ghost[k*5 + 0]*V_ghost[k*5 + 3];
                                
                                const Real epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getInternalEnergy(
                                        &V_ghost[k*5 + 0],
                                        &V_ghost[k*5 + 4],
                                        thermo_properties_ptr);
                                
                                const Real E = V_ghost[k*5 + 0]*epsilon +
                                    half*(Q[1][idx_cell_mom]*Q[1][idx_cell_mom] + Q[2][idx_cell_mom]*Q[2][idx_cell_mom] +
                                       Q[3][idx_cell_mom]*Q[3][idx_cell_mom])/V_ghost[k*5 + 0];
                                
                                Q[4][idx_cell_E] = E;
                            }
                        }
                    }
                }
                else if (face_loc == BDRY_LOC::ZHI)
                {
                    const int num_ghosts_to_fill = fill_box_hi_idx[2] - fill_box_lo_idx[2] + 1;
                    TBOX_ASSERT(fill_box_lo_idx[2] == interior_box_hi_idx[2] + 1);
                    if (num_ghosts_to_fill > 6)
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelBoundaryUtilitiesSingleSpecies::fill3dFaceBoundaryData()\n"
                            << "Non-reflecting outflow BC doesn't support more than six ghost cells yet!");
                    }
                    
                    for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                    {
                        for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                        {
                            // Get the grid spacing.
                            const double* const dx = patch_geom->getDx();
                            
                            const int idx_cell_rho_z_B = (i + num_subghosts_conservative_var[0][0]) +
                                (j + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                (interior_box_hi_idx[2] + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][1]*
                                    subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_rho_z_BB = (i + num_subghosts_conservative_var[0][0]) +
                                (j + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                (interior_box_hi_idx[2] - 1 + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                    subghostcell_dims_conservative_var[0][1];
                            
                            const int idx_cell_rho_z_BBB = (i + num_subghosts_conservative_var[0][0]) +
                                (j + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                (interior_box_hi_idx[2] - 2 + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                    subghostcell_dims_conservative_var[0][1];
                            
                            const int idx_cell_mom_z_B = (i + num_subghosts_conservative_var[1][0]) +
                                (j + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                (interior_box_hi_idx[2] + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                    subghostcell_dims_conservative_var[1][1];
                            
                            const int idx_cell_mom_z_BB = (i + num_subghosts_conservative_var[1][0]) +
                                (j + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                (interior_box_hi_idx[2] - 1 + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                    subghostcell_dims_conservative_var[1][1];
                            
                            const int idx_cell_mom_z_BBB = (i + num_subghosts_conservative_var[1][0]) +
                                (j + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                (interior_box_hi_idx[2] - 2 + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                    subghostcell_dims_conservative_var[1][1];
                            
                            const int idx_cell_E_z_B = (i + num_subghosts_conservative_var[2][0]) +
                                (j + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                (interior_box_hi_idx[2] + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                    subghostcell_dims_conservative_var[2][1];
                            
                            const int idx_cell_E_z_BB = (i + num_subghosts_conservative_var[2][0]) +
                                (j + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                (interior_box_hi_idx[2] - 1 + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                    subghostcell_dims_conservative_var[2][1];
                            
                            const int idx_cell_E_z_BBB = (i + num_subghosts_conservative_var[2][0]) +
                                (j + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                (interior_box_hi_idx[2] - 2 + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                    subghostcell_dims_conservative_var[2][1];
                            
                            const Real& rho_z_B   = Q[0][idx_cell_rho_z_B];
                            const Real& rho_z_BB  = Q[0][idx_cell_rho_z_BB];
                            const Real& rho_z_BBB = Q[0][idx_cell_rho_z_BBB];
                            
                            const Real u_z_B   = Q[1][idx_cell_mom_z_B]/rho_z_B;
                            const Real u_z_BB  = Q[1][idx_cell_mom_z_BB]/rho_z_BB;
                            const Real u_z_BBB = Q[1][idx_cell_mom_z_BBB]/rho_z_BBB;
                            
                            const Real v_z_B   = Q[2][idx_cell_mom_z_B]/rho_z_B;
                            const Real v_z_BB  = Q[2][idx_cell_mom_z_BB]/rho_z_BB;
                            const Real v_z_BBB = Q[2][idx_cell_mom_z_BBB]/rho_z_BBB;
                            
                            const Real w_z_B   = Q[3][idx_cell_mom_z_B]/rho_z_B;
                            const Real w_z_BB  = Q[3][idx_cell_mom_z_BB]/rho_z_BB;
                            const Real w_z_BBB = Q[3][idx_cell_mom_z_BBB]/rho_z_BBB;
                            
                            const Real half = Real(1)/Real(2);
                            const Real epsilon_z_B   = Q[4][idx_cell_E_z_B]/rho_z_B - half*(u_z_B*u_z_B + v_z_B*v_z_B + w_z_B*w_z_B);
                            const Real epsilon_z_BB  = Q[4][idx_cell_E_z_BB]/rho_z_BB - half*(u_z_BB*u_z_BB + v_z_BB*v_z_BB + w_z_BB*w_z_BB);
                            const Real epsilon_z_BBB = Q[4][idx_cell_E_z_BBB]/rho_z_BBB - half*(u_z_BBB*u_z_BBB + v_z_BBB*v_z_BBB + w_z_BBB*w_z_BBB);
                            
                            const Real p_z_B = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getPressure(
                                    &rho_z_B,
                                    &epsilon_z_B,
                                    thermo_properties_ptr);
                            
                            const Real p_z_BB = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getPressure(
                                    &rho_z_BB,
                                    &epsilon_z_BB,
                                    thermo_properties_ptr);
                            
                            const Real p_z_BBB = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getPressure(
                                    &rho_z_BBB,
                                    &epsilon_z_BBB,
                                    thermo_properties_ptr);
                            
                            /*
                             * Compute derivatives at z-direction.
                             */
                            
                            const Real drho_dz = (rho_z_BBB - Real(4)*rho_z_BB + Real(3)*rho_z_B)/(Real(2)*Real(dx[2]));
                            const Real du_dz   = (u_z_BBB - Real(4)*u_z_BB + Real(3)*u_z_B)/(Real(2)*Real(dx[2]));
                            const Real dv_dz   = (v_z_BBB - Real(4)*v_z_BB + Real(3)*v_z_B)/(Real(2)*Real(dx[2]));
                            const Real dw_dz   = (w_z_BBB - Real(4)*w_z_BB + Real(3)*w_z_B)/(Real(2)*Real(dx[2]));
                            const Real dp_dz   = (p_z_BBB - Real(4)*p_z_BB + Real(3)*p_z_B)/(Real(2)*Real(dx[2]));
                            
                            /*
                             * Compute derivatives in x-direction.
                             */
                            
                            Real du_dx = Real(0);
                            // Real dv_dx = Real(0);
                            Real dw_dx = Real(0);
                            Real dp_dx = Real(0);
                            
                            if ((i + num_subghosts_conservative_var[0][0] == 0) ||
                                (i + num_subghosts_conservative_var[1][0] == 0) ||
                                (i + num_subghosts_conservative_var[2][0] == 0))
                            // if (((patch_geom->getTouchesRegularBoundary(0, 0)) && (i == interior_box_lo_idx[0])) ||
                            //     ((i + num_subghosts_conservative_var[0][0] == 0) ||
                            //      (i + num_subghosts_conservative_var[1][0] == 0) ||
                            //      (i + num_subghosts_conservative_var[2][0] == 0)))
                            {
                                // Patch is touching periodic boundary.
                                
                                const int idx_cell_rho_x_R = (i + 1 + num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (interior_box_hi_idx[2] + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom_x_R = (i + 1 + num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                    (interior_box_hi_idx[2] + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E_x_R = (i + 1 + num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                    (interior_box_hi_idx[2] + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                const Real& rho_x_R = Q[0][idx_cell_rho_x_R];
                                const Real u_x_R = Q[1][idx_cell_mom_x_R]/rho_x_R;
                                const Real v_x_R = Q[2][idx_cell_mom_x_R]/rho_x_R;
                                const Real w_x_R = Q[3][idx_cell_mom_x_R]/rho_x_R;
                                const Real epsilon_x_R = Q[4][idx_cell_E_x_R]/rho_x_R - half*(u_x_R*u_x_R + v_x_R*v_x_R + w_x_R*w_x_R);
                                
                                const Real p_x_R = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &rho_x_R,
                                        &epsilon_x_R,
                                        thermo_properties_ptr);
                                
                                // One-sided derivatives.
                                du_dx = (u_x_R - u_z_B)/Real(dx[0]);
                                // dv_dx = (v_x_R - v_z_B)/Real(dx[0]);
                                dw_dx = (w_x_R - w_z_B)/Real(dx[0]);
                                dp_dx = (p_x_R - p_z_B)/Real(dx[0]);
                            }
                            else if ((i + num_subghosts_conservative_var[0][0] + 1 == subghostcell_dims_conservative_var[0][0]) ||
                                     (i + num_subghosts_conservative_var[1][0] + 1 == subghostcell_dims_conservative_var[1][0]) ||
                                     (i + num_subghosts_conservative_var[2][0] + 1 == subghostcell_dims_conservative_var[2][0]))
                            // else if (((patch_geom->getTouchesRegularBoundary(0, 1)) && (i == interior_box_hi_idx[0])) ||
                            //          ((i + num_subghosts_conservative_var[0][0] + 1 == subghostcell_dims_conservative_var[0][0]) ||
                            //           (i + num_subghosts_conservative_var[1][0] + 1 == subghostcell_dims_conservative_var[1][0]) ||
                            //           (i + num_subghosts_conservative_var[2][0] + 1 == subghostcell_dims_conservative_var[2][0])))
                            {
                                // Patch is touching periodic boundary.
                                
                                const int idx_cell_rho_x_L = (i - 1 + num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (interior_box_hi_idx[2] + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom_x_L = (i - 1 + num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                    (interior_box_hi_idx[2] + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E_x_L = (i - 1 + num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                    (interior_box_hi_idx[2] + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                const Real& rho_x_L = Q[0][idx_cell_rho_x_L];
                                const Real u_x_L = Q[1][idx_cell_mom_x_L]/rho_x_L;
                                const Real v_x_L = Q[2][idx_cell_mom_x_L]/rho_x_L;
                                const Real w_x_L = Q[3][idx_cell_mom_x_L]/rho_x_L;
                                const Real epsilon_x_L = Q[4][idx_cell_E_x_L]/rho_x_L - half*(u_x_L*u_x_L + v_x_L*v_x_L + w_x_L*w_x_L);
                                
                                const Real p_x_L = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &rho_x_L,
                                        &epsilon_x_L,
                                        thermo_properties_ptr);
                                
                                // One-sided derivatives.
                                du_dx = (u_z_B - u_x_L)/Real(dx[0]);
                                // dv_dx = (v_z_B - v_x_L)/Real(dx[0]);
                                dw_dx = (w_z_B - w_x_L)/Real(dx[0]);
                                dp_dx = (p_z_B - p_x_L)/Real(dx[0]);
                            }
                            else
                            {
                                const int idx_cell_rho_x_L = (i - 1 + num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (interior_box_hi_idx[2] + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_rho_x_R = (i + 1 + num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (interior_box_hi_idx[2] + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom_x_L = (i - 1 + num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                    (interior_box_hi_idx[2] + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_mom_x_R = (i + 1 + num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                    (interior_box_hi_idx[2] + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E_x_L = (i - 1 + num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                    (interior_box_hi_idx[2] + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                const int idx_cell_E_x_R = (i + 1 + num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                    (interior_box_hi_idx[2] + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                const Real& rho_x_L = Q[0][idx_cell_rho_x_L];
                                const Real& rho_x_R = Q[0][idx_cell_rho_x_R];
                                
                                const Real u_x_L = Q[1][idx_cell_mom_x_L]/rho_x_L;
                                const Real u_x_R = Q[1][idx_cell_mom_x_R]/rho_x_R;
                                
                                const Real v_x_L = Q[2][idx_cell_mom_x_L]/rho_x_L;
                                const Real v_x_R = Q[2][idx_cell_mom_x_R]/rho_x_R;
                                
                                const Real w_x_L = Q[3][idx_cell_mom_x_L]/rho_x_L;
                                const Real w_x_R = Q[3][idx_cell_mom_x_R]/rho_x_R;
                                
                                const Real epsilon_x_L = Q[4][idx_cell_E_x_L]/rho_x_L - half*(u_x_L*u_x_L + v_x_L*v_x_L + w_x_L*w_x_L);
                                const Real epsilon_x_R = Q[4][idx_cell_E_x_R]/rho_x_R - half*(u_x_R*u_x_R + v_x_R*v_x_R + w_x_R*w_x_R);
                                
                                const Real p_x_L = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &rho_x_L,
                                        &epsilon_x_L,
                                        thermo_properties_ptr);
                                
                                const Real p_x_R = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &rho_x_R,
                                        &epsilon_x_R,
                                        thermo_properties_ptr);
                                
                                // Central derivatives.
                                du_dx = (u_x_R - u_x_L)/(Real(2)*Real(dx[0]));
                                // dv_dx = (v_x_R - v_x_L)/(Real(2)*Real(dx[0]));
                                dw_dx = (w_x_R - w_x_L)/(Real(2)*Real(dx[0]));
                                dp_dx = (p_x_R - p_x_L)/(Real(2)*Real(dx[0]));
                            }
                            
                            /*
                             * Compute derivatives in y-direction.
                             */
                            
                            // Real du_dy = Real(0);
                            Real dv_dy = Real(0);
                            Real dw_dy = Real(0);
                            Real dp_dy = Real(0);
                            
                            if ((j + num_subghosts_conservative_var[0][1] == 0) ||
                                (j + num_subghosts_conservative_var[1][1] == 0) ||
                                (j + num_subghosts_conservative_var[2][1] == 0))
                            // if (((patch_geom->getTouchesRegularBoundary(1, 0)) && (j == interior_box_lo_idx[1])) ||
                            //     ((j + num_subghosts_conservative_var[0][1] == 0) ||
                            //      (j + num_subghosts_conservative_var[1][1] == 0) ||
                            //      (j + num_subghosts_conservative_var[2][1] == 0)))
                            {
                                // Patch is touching periodic boundary.
                                
                                const int idx_cell_rho_y_T = (i + num_subghosts_conservative_var[0][0]) +
                                    (j + 1 + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (interior_box_hi_idx[2] + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom_y_T = (i + num_subghosts_conservative_var[1][0]) +
                                    (j + 1 + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                    (interior_box_hi_idx[2] + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E_y_T = (i + num_subghosts_conservative_var[2][0]) +
                                    (j + 1 + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                    (interior_box_hi_idx[2] + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                const Real& rho_y_T = Q[0][idx_cell_rho_y_T];
                                const Real u_y_T = Q[1][idx_cell_mom_y_T]/rho_y_T;
                                const Real v_y_T = Q[2][idx_cell_mom_y_T]/rho_y_T;
                                const Real w_y_T = Q[3][idx_cell_mom_y_T]/rho_y_T;
                                const Real epsilon_y_T = Q[4][idx_cell_E_y_T]/rho_y_T - half*(u_y_T*u_y_T + v_y_T*v_y_T + w_y_T*w_y_T);
                                
                                const Real p_y_T = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &rho_y_T,
                                        &epsilon_y_T,
                                        thermo_properties_ptr);
                                
                                // One-sided derivatives.
                                // du_dy = (u_y_T - u_z_B)/Real(dx[1]);
                                dv_dy = (v_y_T - v_z_B)/Real(dx[1]);
                                dw_dy = (w_y_T - w_z_B)/Real(dx[1]);
                                dp_dy = (p_y_T - p_z_B)/Real(dx[1]);
                            }
                            else if ((j + num_subghosts_conservative_var[0][1] + 1 == subghostcell_dims_conservative_var[0][1]) ||
                                     (j + num_subghosts_conservative_var[1][1] + 1 == subghostcell_dims_conservative_var[1][1]) ||
                                     (j + num_subghosts_conservative_var[2][1] + 1 == subghostcell_dims_conservative_var[2][1]))
                            // else if (((patch_geom->getTouchesRegularBoundary(1, 1)) && (j == interior_box_hi_idx[1])) ||
                            //          ((j + num_subghosts_conservative_var[0][1] + 1 == subghostcell_dims_conservative_var[0][1]) ||
                            //           (j + num_subghosts_conservative_var[1][1] + 1 == subghostcell_dims_conservative_var[1][1]) ||
                            //           (j + num_subghosts_conservative_var[2][1] + 1 == subghostcell_dims_conservative_var[2][1])))
                            {
                                // Patch is touching periodic boundary.
                                
                                const int idx_cell_rho_y_B = (i + num_subghosts_conservative_var[0][0]) +
                                    (j - 1 + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (interior_box_hi_idx[2] + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom_y_B = (i + num_subghosts_conservative_var[1][0]) +
                                    (j - 1 + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                    (interior_box_hi_idx[2] + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E_y_B = (i + num_subghosts_conservative_var[2][0]) +
                                    (j - 1 + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                    (interior_box_hi_idx[2] + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                const Real& rho_y_B = Q[0][idx_cell_rho_y_B];
                                const Real u_y_B = Q[1][idx_cell_mom_y_B]/rho_y_B;
                                const Real v_y_B = Q[2][idx_cell_mom_y_B]/rho_y_B;
                                const Real w_y_B = Q[3][idx_cell_mom_y_B]/rho_y_B;
                                const Real epsilon_y_B = Q[4][idx_cell_E_y_B]/rho_y_B - half*(u_y_B*u_y_B + v_y_B*v_y_B + w_y_B*w_y_B);
                                
                                const Real p_y_B = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &rho_y_B,
                                        &epsilon_y_B,
                                        thermo_properties_ptr);
                                
                                // One-sided derivatives.
                                // du_dy = (u_z_B - u_y_B)/Real(dx[1]);
                                dv_dy = (v_z_B - v_y_B)/Real(dx[1]);
                                dw_dy = (w_z_B - w_y_B)/Real(dx[1]);
                                dp_dy = (p_z_B - p_y_B)/Real(dx[1]);
                            }
                            else
                            {
                                const int idx_cell_rho_y_B = (i + num_subghosts_conservative_var[0][0]) +
                                    (j - 1 + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (interior_box_hi_idx[2] +  num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_rho_y_T = (i + num_subghosts_conservative_var[0][0]) +
                                    (j + 1 + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (interior_box_hi_idx[2] + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom_y_B = (i + num_subghosts_conservative_var[1][0]) +
                                    (j - 1 + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                    (interior_box_hi_idx[2] + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_mom_y_T = (i + num_subghosts_conservative_var[1][0]) +
                                    (j + 1 + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                    (interior_box_hi_idx[2] + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E_y_B = (i + num_subghosts_conservative_var[2][0]) +
                                    (j - 1 + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                    (interior_box_hi_idx[2] + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                const int idx_cell_E_y_T = (i + num_subghosts_conservative_var[2][0]) +
                                    (j + 1 + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                    (interior_box_hi_idx[2] + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                const Real& rho_y_B = Q[0][idx_cell_rho_y_B];
                                const Real& rho_y_T = Q[0][idx_cell_rho_y_T];
                                
                                const Real u_y_B = Q[1][idx_cell_mom_y_B]/rho_y_B;
                                const Real u_y_T = Q[1][idx_cell_mom_y_T]/rho_y_T;
                                
                                const Real v_y_B = Q[2][idx_cell_mom_y_B]/rho_y_B;
                                const Real v_y_T = Q[2][idx_cell_mom_y_T]/rho_y_T;
                                
                                const Real w_y_B = Q[3][idx_cell_mom_y_B]/rho_y_B;
                                const Real w_y_T = Q[3][idx_cell_mom_y_T]/rho_y_T;
                                
                                const Real epsilon_y_B = Q[4][idx_cell_E_y_B]/rho_y_B - half*(u_y_B*u_y_B + v_y_B*v_y_B + w_y_B*w_y_B);
                                const Real epsilon_y_T = Q[4][idx_cell_E_y_T]/rho_y_T - half*(u_y_T*u_y_T + v_y_T*v_y_T + w_y_T*w_y_T);
                                
                                const Real p_y_B = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &rho_y_B,
                                        &epsilon_y_B,
                                        thermo_properties_ptr);
                                
                                const Real p_y_T = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &rho_y_T,
                                        &epsilon_y_T,
                                        thermo_properties_ptr);
                                
                                // Central derivatives.
                                // du_dy = (u_y_T - u_y_B)/(Real(2)*Real(dx[1]));
                                dv_dy = (v_y_T - v_y_B)/(Real(2)*Real(dx[1]));
                                dw_dy = (w_y_T - w_y_B)/(Real(2)*Real(dx[1]));
                                dp_dy = (p_y_T - p_y_B)/(Real(2)*Real(dx[1]));
                            }
                            
                            const Real c_z_B = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getSoundSpeed(
                                    &rho_z_B,
                                    &p_z_B,
                                    thermo_properties_ptr);
                            
                            const Real lambda_1 = w_z_B - c_z_B;
                            
                            // Compute vector Lambda^(-1) * L.
                            
                            Real Lambda_inv_L[5];
                            
                            const Real& p_t         = d_bdry_face_nonreflecting_outflow_p_t[face_loc];
                            const Real& sigma       = d_bdry_face_nonreflecting_outflow_sigma[face_loc];
                            const Real& beta        = d_bdry_face_nonreflecting_outflow_beta[face_loc];
                            const Real& length_char = d_bdry_face_nonreflecting_outflow_length_char[face_loc];
                            
                            const Real T_1 = u_z_B*(dp_dx - rho_z_B*c_z_B*dw_dx) + rho_z_B*c_z_B*c_z_B*du_dx +
                                v_z_B*(dp_dy - rho_z_B*c_z_B*dw_dy) + rho_z_B*c_z_B*c_z_B*dv_dy;
                            
                            const Real M_sq = (u_z_B*u_z_B + v_z_B*v_z_B + w_z_B*w_z_B)/(c_z_B*c_z_B);
                            const Real K = sigma*c_z_B*(Real(1) - M_sq)/length_char;
                            
                            Lambda_inv_L[0] = (Real(1)/lambda_1)*(K*(p_z_B - p_t) - (Real(1) - beta)*T_1);
                            Lambda_inv_L[1] = du_dz;
                            Lambda_inv_L[2] = dv_dz;
                            Lambda_inv_L[3] = c_z_B*c_z_B*drho_dz - dp_dz;
                            Lambda_inv_L[4] = dp_dz + rho_z_B*c_z_B*dw_dz;
                            
                            // Compute dV_dz.
                            
                            const Real c_sq_inv  = Real(1)/(c_z_B*c_z_B);
                            const Real rho_c_inv = Real(1)/(rho_z_B*c_z_B);
                            
                            Real dV_dz[5];
                            
                            dV_dz[0] = half*c_sq_inv*(Lambda_inv_L[0] + Lambda_inv_L[4]) + c_sq_inv*Lambda_inv_L[3];
                            dV_dz[1] = Lambda_inv_L[1];
                            dV_dz[2] = Lambda_inv_L[2];
                            dV_dz[3] = half*rho_c_inv*(-Lambda_inv_L[0] + Lambda_inv_L[4]);
                            dV_dz[4] = half*(Lambda_inv_L[0] + Lambda_inv_L[4]);
                            
                            Real V_ghost[5*num_ghosts_to_fill];
                            
                            for (int k = 0; k < num_ghosts_to_fill; k++)
                            {
                                const int idx_cell_rho = (i + num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (k + fill_box_lo_idx[2] + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom = (i + num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                    (k + fill_box_lo_idx[2] + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E = (i + num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                    (k + fill_box_lo_idx[2] + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                if (k == 0)
                                {
                                    V_ghost[k*5 + 0] = rho_z_BB + Real(2)*Real(dx[2])*dV_dz[0];
                                    V_ghost[k*5 + 1] = u_z_BB   + Real(2)*Real(dx[2])*dV_dz[1];
                                    V_ghost[k*5 + 2] = v_z_BB   + Real(2)*Real(dx[2])*dV_dz[2];
                                    V_ghost[k*5 + 3] = w_z_BB   + Real(2)*Real(dx[2])*dV_dz[3];
                                    V_ghost[k*5 + 4] = p_z_BB   + Real(2)*Real(dx[2])*dV_dz[4];
                                }
                                else if (k == 1)
                                {
                                    V_ghost[k*5 + 0] = -Real(2)*rho_z_BB - Real(3)*rho_z_B +
                                        Real(6)*V_ghost[(k - 1)*5 + 0] - Real(6)*Real(dx[2])*dV_dz[0];
                                    
                                    V_ghost[k*5 + 1] = -Real(2)*u_z_BB - Real(3)*u_z_B +
                                        Real(6)*V_ghost[(k - 1)*5 + 1] - Real(6)*Real(dx[2])*dV_dz[1];
                                    
                                    V_ghost[k*5 + 2] = -Real(2)*v_z_BB - Real(3)*v_z_B +
                                        Real(6)*V_ghost[(k - 1)*5 + 2] - Real(6)*Real(dx[2])*dV_dz[2];
                                    
                                    V_ghost[k*5 + 3] = -Real(2)*w_z_BB - Real(3)*w_z_B +
                                        Real(6)*V_ghost[(k - 1)*5 + 3] - Real(6)*Real(dx[2])*dV_dz[3];
                                    
                                    V_ghost[k*5 + 4] = -Real(2)*p_z_BB - Real(3)*p_z_B +
                                        Real(6)*V_ghost[(k - 1)*5 + 4] - Real(6)*Real(dx[2])*dV_dz[4];
                                }
                                else if (k == 2)
                                {
                                    V_ghost[k*5 + 0] = Real(3)*rho_z_BB + Real(10)*rho_z_B -
                                        Real(18)*V_ghost[(k - 2)*5 + 0] + Real(6)*V_ghost[(k - 1)*5 + 0] +
                                        Real(12)*Real(dx[2])*dV_dz[0];
                                    
                                    V_ghost[k*5 + 1] = Real(3)*u_z_BB + Real(10)*u_z_B -
                                        Real(18)*V_ghost[(k - 2)*5 + 1] + Real(6)*V_ghost[(k - 1)*5 + 1] +
                                        Real(12)*Real(dx[2])*dV_dz[1];
                                    
                                    V_ghost[k*5 + 2] = Real(3)*v_z_BB + Real(10)*v_z_B -
                                        Real(18)*V_ghost[(k - 2)*5 + 2] + Real(6)*V_ghost[(k - 1)*5 + 2] +
                                        Real(12)*Real(dx[2])*dV_dz[2];
                                    
                                    V_ghost[k*5 + 3] = Real(3)*w_z_BB + Real(10)*w_z_B -
                                        Real(18)*V_ghost[(k - 2)*5 + 3] + Real(6)*V_ghost[(k - 1)*5 + 3] +
                                        Real(12)*Real(dx[2])*dV_dz[3];
                                    
                                    V_ghost[k*5 + 4] = Real(3)*p_z_BB + Real(10)*p_z_B -
                                        Real(18)*V_ghost[(k - 2)*5 + 4] + Real(6)*V_ghost[(k - 1)*5 + 4] +
                                        Real(12)*Real(dx[2])*dV_dz[4];
                                }
                                else if (k == 3)
                                {
                                    V_ghost[k*5 + 0] = -Real(4)*rho_z_BB - Real(65)/Real(3)*rho_z_B +
                                        Real(40)*V_ghost[(k - 3)*5 + 0] - Real(20)*V_ghost[(k - 2)*5 + 0] +
                                        Real(20)/Real(3)*V_ghost[(k - 1)*5 + 0] - Real(20)*Real(dx[2])*dV_dz[0];
                                    
                                    V_ghost[k*5 + 1] = -Real(4)*u_z_BB - Real(65)/Real(3)*u_z_B +
                                        Real(40)*V_ghost[(k - 3)*5 + 1] - Real(20)*V_ghost[(k - 2)*5 + 1] +
                                        Real(20)/Real(3)*V_ghost[(k - 1)*5 + 1] - Real(20)*Real(dx[2])*dV_dz[1];
                                    
                                    V_ghost[k*5 + 2] = -Real(4)*v_z_BB - Real(65)/Real(3)*v_z_B +
                                        Real(40)*V_ghost[(k - 3)*5 + 2] - Real(20)*V_ghost[(k - 2)*5 + 2] +
                                        Real(20)/Real(3)*V_ghost[(k - 1)*5 + 2] - Real(20)*Real(dx[2])*dV_dz[2];
                                    
                                    V_ghost[k*5 + 3] = -Real(4)*w_z_BB - Real(65)/Real(3)*w_z_B +
                                        Real(40)*V_ghost[(k - 3)*5 + 3] - Real(20)*V_ghost[(k - 2)*5 + 3] +
                                        Real(20)/Real(3)*V_ghost[(k - 1)*5 + 3] - Real(20)*Real(dx[2])*dV_dz[3];
                                    
                                    V_ghost[k*5 + 4] = -Real(4)*p_z_BB - Real(65)/Real(3)*p_z_B +
                                        Real(40)*V_ghost[(k - 3)*5 + 4] - Real(20)*V_ghost[(k - 2)*5 + 4] +
                                        Real(20)/Real(3)*V_ghost[(k - 1)*5 + 4] - Real(20)*Real(dx[2])*dV_dz[4];
                                }
                                else if (k == 4)
                                {
                                    V_ghost[k*5 + 0] = Real(5)*rho_z_BB + Real(77)/Real(2)*rho_z_B -
                                        Real(75)*V_ghost[(k - 4)*5 + 0] + Real(50)*V_ghost[(k - 3)*5 + 0] -
                                        Real(25)*V_ghost[(k - 2)*5 + 0] + Real(15)/Real(2)*V_ghost[(k - 1)*5 + 0] +
                                        Real(30)*Real(dx[2])*dV_dz[0];
                                    
                                    V_ghost[k*5 + 1] = Real(5)*u_z_BB + Real(77)/Real(2)*u_z_B -
                                        Real(75)*V_ghost[(k - 4)*5 + 1] + Real(50)*V_ghost[(k - 3)*5 + 1] -
                                        Real(25)*V_ghost[(k - 2)*5 + 1] + Real(15)/Real(2)*V_ghost[(k - 1)*5 + 1] +
                                        Real(30)*Real(dx[2])*dV_dz[1];
                                    
                                    V_ghost[k*5 + 2] = Real(5)*v_z_BB + Real(77)/Real(2)*v_z_B -
                                        Real(75)*V_ghost[(k - 4)*5 + 2] + Real(50)*V_ghost[(k - 3)*5 + 2] -
                                        Real(25)*V_ghost[(k - 2)*5 + 2] + Real(15)/Real(2)*V_ghost[(k - 1)*5 + 2] +
                                        Real(30)*Real(dx[2])*dV_dz[2];
                                    
                                    V_ghost[k*5 + 3] = Real(5)*w_z_BB + Real(77)/Real(2)*w_z_B -
                                        Real(75)*V_ghost[(k - 4)*5 + 3] + Real(50)*V_ghost[(k - 3)*5 + 3] -
                                        Real(25)*V_ghost[(k - 2)*5 + 3] + Real(15)/Real(2)*V_ghost[(k - 1)*5 + 3] +
                                        Real(30)*Real(dx[2])*dV_dz[3];
                                    
                                    V_ghost[k*5 + 4] = Real(5)*p_z_BB + Real(77)/Real(2)*p_z_B -
                                        Real(75)*V_ghost[(k - 4)*5 + 4] + Real(50)*V_ghost[(k - 3)*5 + 4] -
                                        Real(25)*V_ghost[(k - 2)*5 + 4] + Real(15)/Real(2)*V_ghost[(k - 1)*5 + 4] +
                                        Real(30)*Real(dx[2])*dV_dz[4];
                                }
                                else if (k == 5)
                                {
                                    V_ghost[k*5 + 0] = -Real(6)*rho_z_BB - Real(609)/Real(10)*rho_z_B +
                                        Real(126)*V_ghost[(k - 5)*5 + 0] - Real(105)*V_ghost[(k - 4)*5 + 0] +
                                        Real(70)*V_ghost[(k - 3)*5 + 0] - Real(63)/Real(2)*V_ghost[(k - 2)*5 + 0] +
                                        Real(42)/Real(5)*V_ghost[(k - 1)*5 + 0] - Real(42)*Real(dx[2])*dV_dz[0];
                                    
                                    V_ghost[k*5 + 1] = -Real(6)*u_z_BB - Real(609)/Real(10)*u_z_B +
                                        Real(126)*V_ghost[(k - 5)*5 + 1] - Real(105)*V_ghost[(k - 4)*5 + 1] +
                                        Real(70)*V_ghost[(k - 3)*5 + 1] - Real(63)/Real(2)*V_ghost[(k - 2)*5 + 1] +
                                        Real(42)/Real(5)*V_ghost[(k - 1)*5 + 1] - Real(42)*Real(dx[2])*dV_dz[1];
                                    
                                    V_ghost[k*5 + 2] = -Real(6)*v_z_BB - Real(609)/Real(10)*v_z_B +
                                        Real(126)*V_ghost[(k - 5)*5 + 2] - Real(105)*V_ghost[(k - 4)*5 + 2] +
                                        Real(70)*V_ghost[(k - 3)*5 + 2] - Real(63)/Real(2)*V_ghost[(k - 2)*5 + 2] +
                                        Real(42)/Real(5)*V_ghost[(k - 1)*5 + 2] - Real(42)*Real(dx[2])*dV_dz[2];
                                    
                                    V_ghost[k*5 + 3] = -Real(6)*w_z_BB - Real(609)/Real(10)*w_z_B +
                                        Real(126)*V_ghost[(k - 5)*5 + 3] - Real(105)*V_ghost[(k - 4)*5 + 3] +
                                        Real(70)*V_ghost[(k - 3)*5 + 3] - Real(63)/Real(2)*V_ghost[(k - 2)*5 + 3] +
                                        Real(42)/Real(5)*V_ghost[(k - 1)*5 + 3] - Real(42)*Real(dx[2])*dV_dz[3];
                                    
                                    V_ghost[k*5 + 4] = -Real(6)*p_z_BB - Real(609)/Real(10)*p_z_B +
                                        Real(126)*V_ghost[(k - 5)*5 + 4] - Real(105)*V_ghost[(k - 4)*5 + 4] +
                                        Real(70)*V_ghost[(k - 3)*5 + 4] - Real(63)/Real(2)*V_ghost[(k - 2)*5 + 4] +
                                        Real(42)/Real(5)*V_ghost[(k - 1)*5 + 4] - Real(42)*Real(dx[2])*dV_dz[4];
                                }
                                
                                Q[0][idx_cell_rho] = V_ghost[k*5 + 0];
                                Q[1][idx_cell_mom] = V_ghost[k*5 + 0]*V_ghost[k*5 + 1];
                                Q[2][idx_cell_mom] = V_ghost[k*5 + 0]*V_ghost[k*5 + 2];
                                Q[3][idx_cell_mom] = V_ghost[k*5 + 0]*V_ghost[k*5 + 3];
                                
                                const Real epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getInternalEnergy(
                                        &V_ghost[k*5 + 0],
                                        &V_ghost[k*5 + 4],
                                        thermo_properties_ptr);
                                
                                const Real E = V_ghost[k*5 + 0]*epsilon +
                                    half*(Q[1][idx_cell_mom]*Q[1][idx_cell_mom] + Q[2][idx_cell_mom]*Q[2][idx_cell_mom] +
                                       Q[3][idx_cell_mom]*Q[3][idx_cell_mom])/V_ghost[k*5 + 0];
                                
                                Q[4][idx_cell_E] = E;
                            }
                        }
                    }
                }
                
                // Remove face locations that have boundary conditions identified.
                bdry_face_locs.erase(std::remove(bdry_face_locs.begin(), bdry_face_locs.end(), face_loc),
                    bdry_face_locs.end());
            }
        }
    }
    
    for (int fi = 0; fi < static_cast<int>(face_bdry.size()); fi++)
    {
        TBOX_ASSERT(face_bdry[fi].getBoundaryType() == BDRY::FACE3D);
        
        int face_loc = face_bdry[fi].getLocationIndex();
        
        if (std::find(bdry_face_locs.begin(), bdry_face_locs.end(), face_loc) !=
            bdry_face_locs.end())
        {
            TBOX_ERROR("FlowModelBoundaryUtilitiesSingleSpecies::fill3dFaceBoundaryData\n"
                << "Invalid face boundary condition!\n"
                << "face_loc = '" << face_loc << "'." << std::endl
                << "bdry_face_conds[face_loc] = '" << bdry_face_conds[face_loc] << "'."
                << std::endl);
        }
    }
}


/*
 * Function to fill 3d edge boundary values for a patch.
 * Edge locations that have boundary conditions identified are removed from the container.
 */
void
FlowModelBoundaryUtilitiesSingleSpecies::fill3dEdgeBoundaryData(
    const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& conservative_var_data,
    const hier::Patch& patch,
    std::vector<int>& bdry_edge_locs,
    const std::vector<int>& bdry_edge_conds,
    const std::vector<std::vector<Real> >& bdry_face_values,
    const hier::IntVector& ghost_width_to_fill)
{
    TBOX_ASSERT(static_cast<int>(conservative_var_data.size()) == 3);
    for (int vi = 0; vi < static_cast<int>(conservative_var_data.size()); vi++)
    {
        TBOX_ASSERT(conservative_var_data[vi]);
    }
    TBOX_ASSERT(static_cast<int>(bdry_edge_locs.size()) <= NUM_3D_EDGES);
    if (static_cast<int>(bdry_edge_locs.size()) > 0)
    {
        TBOX_ASSERT(*min_element(bdry_edge_locs.begin(), bdry_edge_locs.end()) >= 0);
        TBOX_ASSERT(*max_element(bdry_edge_locs.begin(), bdry_edge_locs.end()) < NUM_3D_EDGES);
    }
    TBOX_ASSERT(static_cast<int>(bdry_edge_conds.size()) == NUM_3D_EDGES);
    TBOX_ASSERT(static_cast<int>(bdry_face_values.size()) == 3);
    for (int vi = 0; vi < static_cast<int>(bdry_face_values.size()); vi++)
    {
        TBOX_ASSERT(static_cast<int>(bdry_face_values[vi].size()) ==
                    NUM_3D_FACES*(conservative_var_data[vi]->getDepth()));
    }
    
    TBOX_DIM_ASSERT(ghost_width_to_fill.getDim() == tbox::Dimension(3));
    
    for (int vi = 0; vi < static_cast<int>(conservative_var_data.size()); vi++)
    {
        TBOX_ASSERT_OBJDIM_EQUALITY3(*conservative_var_data[vi], patch, ghost_width_to_fill);
    }
    
    const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
        HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    TBOX_ASSERT(patch_geom);
    
    hier::IntVector num_ghosts(conservative_var_data[0]->getGhostCellWidth());
    for (int vi = 1; vi < static_cast<int>(conservative_var_data.size()); vi++)
    {
        num_ghosts = hier::IntVector::min(
            num_ghosts,
            conservative_var_data[vi]->getGhostCellWidth());
    }
    
    /*
     * Determine the ghost cell width to fill.
     */
    
    hier::IntVector gcw_to_fill(tbox::Dimension(3));
    
    // If the ghost fill width is not used, it is set to the ghost cell width of the data.
    if (ghost_width_to_fill == -hier::IntVector::getOne(tbox::Dimension(3)))
    {
        gcw_to_fill = num_ghosts;
    }
    else
    {
        gcw_to_fill = hier::IntVector::min(
            num_ghosts,
            ghost_width_to_fill);
    }
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box& interior_box(patch.getBox());
    hier::Index interior_box_lo_idx(interior_box.lower());
    hier::Index interior_box_hi_idx(interior_box.upper());
    
    /*
     * Offset the indices.
     */
    
    interior_box_lo_idx = interior_box_lo_idx - interior_box.lower();
    interior_box_hi_idx = interior_box_hi_idx - interior_box.lower();
    
    const std::vector<hier::BoundaryBox>& edge_bdry =
        patch_geom->getCodimensionBoundaries(BDRY::EDGE3D);
    
    for (int ei = 0; ei < static_cast<int>(edge_bdry.size()); ei++)
    {
        TBOX_ASSERT(edge_bdry[ei].getBoundaryType() == BDRY::EDGE3D);
        
        int edge_loc(edge_bdry[ei].getLocationIndex());
        
        if (std::find(bdry_edge_locs.begin(), bdry_edge_locs.end(), edge_loc) !=
            bdry_edge_locs.end())
        {
            hier::Box fill_box(patch_geom->getBoundaryFillBox(
                edge_bdry[ei],
                interior_box,
                gcw_to_fill));
            
            hier::Index fill_box_lo_idx(fill_box.lower());
            hier::Index fill_box_hi_idx(fill_box.upper());
            
            /*
             * Offset the indices.
             */
            
            fill_box_lo_idx = fill_box_lo_idx - interior_box.lower();
            fill_box_hi_idx = fill_box_hi_idx - interior_box.lower();
            
            int face_loc_0 = -1;
            int face_loc_1 = -1;
            int face_loc_2 = -1;
            
            switch (edge_loc)
            {
                case EDGE_BDRY_LOC_3D::XLO_YLO:
                {
                    face_loc_0 = BDRY_LOC::XLO;
                    face_loc_1 = BDRY_LOC::YLO;
                    
                    break;
                }
                case EDGE_BDRY_LOC_3D::XHI_YLO:
                {
                    face_loc_0 = BDRY_LOC::XHI;
                    face_loc_1 = BDRY_LOC::YLO;
                    
                    break;
                }
                case EDGE_BDRY_LOC_3D::XLO_YHI:
                {
                    face_loc_0 = BDRY_LOC::XLO;
                    face_loc_1 = BDRY_LOC::YHI;
                    
                    break;
                }
                case EDGE_BDRY_LOC_3D::XHI_YHI:
                {
                    face_loc_0 = BDRY_LOC::XHI;
                    face_loc_1 = BDRY_LOC::YHI;
                    
                    break;
                }
                case EDGE_BDRY_LOC_3D::XLO_ZLO:
                {
                    face_loc_0 = BDRY_LOC::XLO;
                    face_loc_2 = BDRY_LOC::ZLO;
                    
                    break;
                }
                case EDGE_BDRY_LOC_3D::XHI_ZLO:
                {
                    face_loc_0 = BDRY_LOC::XHI;
                    face_loc_2 = BDRY_LOC::ZLO;
                    
                    break;
                }
                case EDGE_BDRY_LOC_3D::XLO_ZHI:
                {
                    face_loc_0 = BDRY_LOC::XLO;
                    face_loc_2 = BDRY_LOC::ZHI;
                    
                    break;
                }
                case EDGE_BDRY_LOC_3D::XHI_ZHI:
                {
                    face_loc_0 = BDRY_LOC::XHI;
                    face_loc_2 = BDRY_LOC::ZHI;
                    
                    break;
                }
                case EDGE_BDRY_LOC_3D::YLO_ZLO:
                {
                    face_loc_1 = BDRY_LOC::YLO;
                    face_loc_2 = BDRY_LOC::ZLO;
                    
                    break;
                }
                case EDGE_BDRY_LOC_3D::YHI_ZLO:
                {
                    face_loc_1 = BDRY_LOC::YHI;
                    face_loc_2 = BDRY_LOC::ZLO;
                    
                    break;
                }
                case EDGE_BDRY_LOC_3D::YLO_ZHI:
                {
                    face_loc_1 = BDRY_LOC::YLO;
                    face_loc_2 = BDRY_LOC::ZHI;
                    
                    break;
                }
                case EDGE_BDRY_LOC_3D::YHI_ZHI:
                {
                    face_loc_1 = BDRY_LOC::YHI;
                    face_loc_2 = BDRY_LOC::ZHI;
                    
                    break;
                }
            }
            
            if ((bdry_edge_conds[edge_loc] == BDRY_COND::FLOW_MODEL::XADIABATIC_NO_SLIP) ||
                (bdry_edge_conds[edge_loc] == BDRY_COND::FLOW_MODEL::YADIABATIC_NO_SLIP) ||
                (bdry_edge_conds[edge_loc] == BDRY_COND::FLOW_MODEL::ZADIABATIC_NO_SLIP) ||
                (bdry_edge_conds[edge_loc] == BDRY_COND::FLOW_MODEL::XISOTHERMAL_NO_SLIP) ||
                (bdry_edge_conds[edge_loc] == BDRY_COND::FLOW_MODEL::YISOTHERMAL_NO_SLIP) ||
                (bdry_edge_conds[edge_loc] == BDRY_COND::FLOW_MODEL::ZISOTHERMAL_NO_SLIP))
            {
                /*
                 * Get the pointers to the conservative variables.
                 * Get the numbers of ghost cells and the dimensions of the ghost cell boxes of
                 * the conservative variables.
                 */
                
                std::vector<Real*> Q;
                Q.reserve(d_num_eqn);
                
                std::vector<hier::IntVector> num_subghosts_conservative_var;
                std::vector<hier::IntVector> subghostcell_dims_conservative_var;
                
                num_subghosts_conservative_var.reserve(3);
                subghostcell_dims_conservative_var.reserve(3);
                
                int count_eqn = 0;
                
                for (int vi = 0; vi < static_cast<int>(conservative_var_data.size()); vi++)
                {
                    int depth = conservative_var_data[vi]->getDepth();
                    
                    for (int di = 0; di < depth; di++)
                    {
                        // If the last element of the conservative variable vector is not in the
                        // system of equations, ignore it.
                        if (count_eqn >= d_num_eqn)
                            break;
                        
                        Q.push_back(conservative_var_data[vi]->getPointer(di));
                        
                        count_eqn++;
                    }
                    
                    num_subghosts_conservative_var.push_back(
                        conservative_var_data[vi]->getGhostCellWidth());
                    subghostcell_dims_conservative_var.push_back(
                        conservative_var_data[vi]->getGhostBox().numberCells());
                }
                
                // Get the thermodynamic properties of the species.
                std::vector<const Real*> thermo_properties_ptr;
                thermo_properties_ptr.reserve(static_cast<int> (d_thermo_properties.size()));
                for (int ti = 0; ti < static_cast<int> (d_thermo_properties.size()); ti++)
                {
                    thermo_properties_ptr.push_back(&d_thermo_properties[ti]);
                }
                
                if (bdry_edge_conds[edge_loc] == BDRY_COND::FLOW_MODEL::XADIABATIC_NO_SLIP)
                {
                    for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                    {
                        for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                        {
                            for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                            {
                                const int idx_cell_rho = (i + num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0] +
                                    (k + num_subghosts_conservative_var[0][2])*
                                        subghostcell_dims_conservative_var[0][0]*
                                            subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom = (i + num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0] +
                                    (k + num_subghosts_conservative_var[1][2])*
                                        subghostcell_dims_conservative_var[1][0]*
                                            subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E = (i + num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0] +
                                    (k + num_subghosts_conservative_var[2][2])*
                                        subghostcell_dims_conservative_var[2][0]*
                                            subghostcell_dims_conservative_var[2][1];
                                
                                int idx_cell_pivot_rho = idx_cell_rho;
                                int idx_cell_pivot_mom = idx_cell_mom;
                                int idx_cell_pivot_E = idx_cell_E;
                                
                                if (face_loc_0 == BDRY_LOC::XLO)
                                {
                                    idx_cell_pivot_rho = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                            num_subghosts_conservative_var[0][0]) +
                                        (j + num_subghosts_conservative_var[0][1])*
                                            subghostcell_dims_conservative_var[0][0] +
                                        (k + num_subghosts_conservative_var[0][2])*
                                            subghostcell_dims_conservative_var[0][0]*
                                                subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                            num_subghosts_conservative_var[1][0]) +
                                        (j + num_subghosts_conservative_var[1][1])*
                                            subghostcell_dims_conservative_var[1][0] +
                                        (k + num_subghosts_conservative_var[1][2])*
                                            subghostcell_dims_conservative_var[1][0]*
                                                subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                            num_subghosts_conservative_var[2][0]) +
                                        (j + num_subghosts_conservative_var[2][1])*
                                            subghostcell_dims_conservative_var[2][0] +
                                        (k + num_subghosts_conservative_var[2][2])*
                                            subghostcell_dims_conservative_var[2][0]*
                                                subghostcell_dims_conservative_var[2][1];
                                }
                                else if (face_loc_0 == BDRY_LOC::XHI)
                                {
                                    idx_cell_pivot_rho = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                            num_subghosts_conservative_var[0][0]) +
                                        (j + num_subghosts_conservative_var[0][1])*
                                            subghostcell_dims_conservative_var[0][0] +
                                        (k + num_subghosts_conservative_var[0][2])*
                                            subghostcell_dims_conservative_var[0][0]*
                                                subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                            num_subghosts_conservative_var[1][0]) +
                                        (j + num_subghosts_conservative_var[1][1])*
                                            subghostcell_dims_conservative_var[1][0] +
                                        (k + num_subghosts_conservative_var[1][2])*
                                            subghostcell_dims_conservative_var[1][0]*
                                                subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                            num_subghosts_conservative_var[2][0]) +
                                        (j + num_subghosts_conservative_var[2][1])*
                                            subghostcell_dims_conservative_var[2][0] +
                                        (k + num_subghosts_conservative_var[2][2])*
                                            subghostcell_dims_conservative_var[2][0]*
                                                subghostcell_dims_conservative_var[2][1];
                                }
                                
                                /*
                                 * Set the values for density and momentum.
                                 */
                                
                                Q[0][idx_cell_rho] = Q[0][idx_cell_pivot_rho];
                                Q[1][idx_cell_mom] = -Q[1][idx_cell_pivot_mom] +
                                    Real(2)*Q[0][idx_cell_pivot_rho]*d_bdry_face_adiabatic_no_slip_vel[face_loc_0*3];
                                Q[2][idx_cell_mom] = -Q[2][idx_cell_pivot_mom] +
                                    Real(2)*Q[0][idx_cell_pivot_rho]*d_bdry_face_adiabatic_no_slip_vel[face_loc_0*3 + 1];
                                Q[3][idx_cell_mom] = -Q[3][idx_cell_pivot_mom] +
                                    Real(2)*Q[0][idx_cell_pivot_rho]*d_bdry_face_adiabatic_no_slip_vel[face_loc_0*3 + 2];
                                
                                /*
                                 * Set the values for total internal energy.
                                 */
                                
                                Real epsilon_pivot = (Q[4][idx_cell_pivot_E] -
                                    Real(1)/Real(2)*(Q[1][idx_cell_pivot_mom]*Q[1][idx_cell_pivot_mom] +
                                    Q[2][idx_cell_pivot_mom]*Q[2][idx_cell_pivot_mom] +
                                    Q[3][idx_cell_pivot_mom]*Q[3][idx_cell_pivot_mom])/
                                    Q[0][idx_cell_pivot_rho])/Q[0][idx_cell_pivot_rho];
                                
                                Real p_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &Q[0][idx_cell_pivot_rho],
                                        &epsilon_pivot,
                                        thermo_properties_ptr);
                                
                                Real T_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getTemperature(
                                        &Q[0][idx_cell_pivot_rho],
                                        &p_pivot,
                                        thermo_properties_ptr);
                                
                                Real T = T_pivot;
                                
                                Real epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getInternalEnergyFromTemperature(
                                        &Q[0][idx_cell_rho],
                                        &T,
                                        thermo_properties_ptr);
                                
                                Real E = Q[0][idx_cell_rho]*epsilon +
                                    Real(1)/Real(2)*(Q[1][idx_cell_mom]*Q[1][idx_cell_mom] + Q[2][idx_cell_mom]*Q[2][idx_cell_mom] +
                                    Q[3][idx_cell_mom]*Q[3][idx_cell_mom])/
                                    Q[0][idx_cell_rho];
                                
                                Q[4][idx_cell_E] = E;
                            }
                        }
                    }
                    
                    // Remove edge locations that have boundary conditions identified.
                    bdry_edge_locs.erase(std::remove(bdry_edge_locs.begin(), bdry_edge_locs.end(), edge_loc),
                        bdry_edge_locs.end());
                }
                else if (bdry_edge_conds[edge_loc] == BDRY_COND::FLOW_MODEL::YADIABATIC_NO_SLIP)
                {
                    for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                    {
                        for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                        {
                            for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                            {
                                const int idx_cell_rho = (i + num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0] +
                                    (k + num_subghosts_conservative_var[0][2])*
                                        subghostcell_dims_conservative_var[0][0]*
                                            subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom = (i + num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0] +
                                    (k + num_subghosts_conservative_var[1][2])*
                                        subghostcell_dims_conservative_var[1][0]*
                                            subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E = (i + num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0] +
                                    (k + num_subghosts_conservative_var[2][2])*
                                        subghostcell_dims_conservative_var[2][0]*
                                            subghostcell_dims_conservative_var[2][1];
                                
                                int idx_cell_pivot_rho = idx_cell_rho;
                                int idx_cell_pivot_mom = idx_cell_mom;
                                int idx_cell_pivot_E = idx_cell_E;
                                
                                if (face_loc_1 == BDRY_LOC::YLO)
                                {
                                    idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                        (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) +
                                            num_subghosts_conservative_var[0][1])*
                                                subghostcell_dims_conservative_var[0][0] +
                                        (k + num_subghosts_conservative_var[0][2])*
                                            subghostcell_dims_conservative_var[0][0]*
                                                subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                        (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) +
                                            num_subghosts_conservative_var[1][1])*
                                                subghostcell_dims_conservative_var[1][0] +
                                        (k + num_subghosts_conservative_var[1][2])*
                                            subghostcell_dims_conservative_var[1][0]*
                                                subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                        (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) +
                                            num_subghosts_conservative_var[2][1])*
                                                subghostcell_dims_conservative_var[2][0] +
                                        (k + num_subghosts_conservative_var[2][2])*
                                            subghostcell_dims_conservative_var[2][0]*
                                                subghostcell_dims_conservative_var[2][1];
                                }
                                else if (face_loc_1 == BDRY_LOC::YHI)
                                {
                                    idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                        (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) +
                                            num_subghosts_conservative_var[0][1])*
                                                subghostcell_dims_conservative_var[0][0] +
                                        (k + num_subghosts_conservative_var[0][2])*
                                            subghostcell_dims_conservative_var[0][0]*
                                                subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                        (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) +
                                            num_subghosts_conservative_var[1][1])*
                                                subghostcell_dims_conservative_var[1][0] +
                                        (k + num_subghosts_conservative_var[1][2])*
                                            subghostcell_dims_conservative_var[1][0]*
                                                subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                        (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) +
                                            num_subghosts_conservative_var[2][1])*
                                                subghostcell_dims_conservative_var[2][0] +
                                        (k + num_subghosts_conservative_var[2][2])*
                                            subghostcell_dims_conservative_var[2][0]*
                                                subghostcell_dims_conservative_var[2][1];
                                }
                                
                                /*
                                 * Set the values for density and momentum.
                                 */
                                
                                Q[0][idx_cell_rho] = Q[0][idx_cell_pivot_rho];
                                Q[1][idx_cell_mom] = -Q[1][idx_cell_pivot_mom] +
                                    Real(2)*Q[0][idx_cell_pivot_rho]*d_bdry_face_adiabatic_no_slip_vel[face_loc_1*3];
                                Q[2][idx_cell_mom] = -Q[2][idx_cell_pivot_mom] +
                                    Real(2)*Q[0][idx_cell_pivot_rho]*d_bdry_face_adiabatic_no_slip_vel[face_loc_1*3 + 1];
                                Q[3][idx_cell_mom] = -Q[3][idx_cell_pivot_mom] +
                                    Real(2)*Q[0][idx_cell_pivot_rho]*d_bdry_face_adiabatic_no_slip_vel[face_loc_1*3 + 2];
                                
                                /*
                                 * Set the values for total internal energy.
                                 */
                                
                                Real epsilon_pivot = (Q[4][idx_cell_pivot_E] -
                                    Real(1)/Real(2)*(Q[1][idx_cell_pivot_mom]*Q[1][idx_cell_pivot_mom] +
                                    Q[2][idx_cell_pivot_mom]*Q[2][idx_cell_pivot_mom] +
                                    Q[3][idx_cell_pivot_mom]*Q[3][idx_cell_pivot_mom])/
                                    Q[0][idx_cell_pivot_rho])/Q[0][idx_cell_pivot_rho];
                                
                                Real p_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &Q[0][idx_cell_pivot_rho],
                                        &epsilon_pivot,
                                        thermo_properties_ptr);
                                
                                Real T_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getTemperature(
                                        &Q[0][idx_cell_pivot_rho],
                                        &p_pivot,
                                        thermo_properties_ptr);
                                
                                Real T = T_pivot;
                                
                                Real epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getInternalEnergyFromTemperature(
                                        &Q[0][idx_cell_rho],
                                        &T,
                                        thermo_properties_ptr);
                                
                                Real E = Q[0][idx_cell_rho]*epsilon +
                                    Real(1)/Real(2)*(Q[1][idx_cell_mom]*Q[1][idx_cell_mom] + Q[2][idx_cell_mom]*Q[2][idx_cell_mom] +
                                    Q[3][idx_cell_mom]*Q[3][idx_cell_mom])/
                                    Q[0][idx_cell_rho];
                                
                                Q[4][idx_cell_E] = E;
                            }
                        }
                    }
                    
                    // Remove edge locations that have boundary conditions identified.
                    bdry_edge_locs.erase(std::remove(bdry_edge_locs.begin(), bdry_edge_locs.end(), edge_loc),
                        bdry_edge_locs.end());
                }
                else if (bdry_edge_conds[edge_loc] == BDRY_COND::FLOW_MODEL::ZADIABATIC_NO_SLIP)
                {
                    for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                    {
                        for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                        {
                            for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                            {
                                const int idx_cell_rho = (i + num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0] +
                                    (k + num_subghosts_conservative_var[0][2])*
                                        subghostcell_dims_conservative_var[0][0]*
                                            subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom = (i + num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0] +
                                    (k + num_subghosts_conservative_var[1][2])*
                                        subghostcell_dims_conservative_var[1][0]*
                                            subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E = (i + num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0] +
                                    (k + num_subghosts_conservative_var[2][2])*
                                        subghostcell_dims_conservative_var[2][0]*
                                            subghostcell_dims_conservative_var[2][1];
                                
                                int idx_cell_pivot_rho = idx_cell_rho;
                                int idx_cell_pivot_mom = idx_cell_mom;
                                int idx_cell_pivot_E = idx_cell_E;
                                
                                if (face_loc_2 == BDRY_LOC::ZLO)
                                {
                                    idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                        (j + num_subghosts_conservative_var[0][1])*
                                            subghostcell_dims_conservative_var[0][0] +
                                        (interior_box_lo_idx[2] + (fill_box_hi_idx[2] - k) +
                                            num_subghosts_conservative_var[0][2])*
                                                subghostcell_dims_conservative_var[0][0]*
                                                    subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                        (j + num_subghosts_conservative_var[1][1])*
                                            subghostcell_dims_conservative_var[1][0] +
                                        (interior_box_lo_idx[2] + (fill_box_hi_idx[2] - k) +
                                            num_subghosts_conservative_var[1][2])*
                                                subghostcell_dims_conservative_var[1][0]*
                                                    subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                        (j + num_subghosts_conservative_var[2][1])*
                                            subghostcell_dims_conservative_var[2][0] +
                                        (interior_box_lo_idx[2] + (fill_box_hi_idx[2] - k) +
                                            num_subghosts_conservative_var[2][2])*
                                                subghostcell_dims_conservative_var[2][0]*
                                                    subghostcell_dims_conservative_var[2][1];
                                }
                                else if (face_loc_2 == BDRY_LOC::ZHI)
                                {
                                    idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                        (j + num_subghosts_conservative_var[0][1])*
                                            subghostcell_dims_conservative_var[0][0] +
                                        (interior_box_hi_idx[2] - (k - fill_box_lo_idx[2]) +
                                            num_subghosts_conservative_var[0][2])*
                                                subghostcell_dims_conservative_var[0][0]*
                                                    subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                        (j + num_subghosts_conservative_var[1][1])*
                                            subghostcell_dims_conservative_var[1][0] +
                                        (interior_box_hi_idx[2] - (k - fill_box_lo_idx[2]) +
                                            num_subghosts_conservative_var[1][2])*
                                                subghostcell_dims_conservative_var[1][0]*
                                                    subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                        (j + num_subghosts_conservative_var[2][1])*
                                            subghostcell_dims_conservative_var[2][0] +
                                        (interior_box_hi_idx[2] - (k - fill_box_lo_idx[2]) +
                                            num_subghosts_conservative_var[2][2])*
                                                subghostcell_dims_conservative_var[2][0]*
                                                    subghostcell_dims_conservative_var[2][1];
                                }
                                
                                /*
                                 * Set the values for density and momentum.
                                 */
                                
                                Q[0][idx_cell_rho] = Q[0][idx_cell_pivot_rho];
                                Q[1][idx_cell_mom] = -Q[1][idx_cell_pivot_mom] +
                                    Real(2)*Q[0][idx_cell_pivot_rho]*d_bdry_face_adiabatic_no_slip_vel[face_loc_2*3];
                                Q[2][idx_cell_mom] = -Q[2][idx_cell_pivot_mom] +
                                    Real(2)*Q[0][idx_cell_pivot_rho]*d_bdry_face_adiabatic_no_slip_vel[face_loc_2*3 + 1];
                                Q[3][idx_cell_mom] = -Q[3][idx_cell_pivot_mom] +
                                    Real(2)*Q[0][idx_cell_pivot_rho]*d_bdry_face_adiabatic_no_slip_vel[face_loc_2*3 + 2];
                                
                                /*
                                 * Set the values for total internal energy.
                                 */
                                
                                Real epsilon_pivot = (Q[4][idx_cell_pivot_E] -
                                    Real(1)/Real(2)*(Q[1][idx_cell_pivot_mom]*Q[1][idx_cell_pivot_mom] +
                                    Q[2][idx_cell_pivot_mom]*Q[2][idx_cell_pivot_mom] +
                                    Q[3][idx_cell_pivot_mom]*Q[3][idx_cell_pivot_mom])/
                                    Q[0][idx_cell_pivot_rho])/Q[0][idx_cell_pivot_rho];
                                
                                Real p_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &Q[0][idx_cell_pivot_rho],
                                        &epsilon_pivot,
                                        thermo_properties_ptr);
                                
                                Real T_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getTemperature(
                                        &Q[0][idx_cell_pivot_rho],
                                        &p_pivot,
                                        thermo_properties_ptr);
                                
                                Real T = T_pivot;
                                
                                Real epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getInternalEnergyFromTemperature(
                                        &Q[0][idx_cell_rho],
                                        &T,
                                        thermo_properties_ptr);
                                
                                Real E = Q[0][idx_cell_rho]*epsilon +
                                    Real(1)/Real(2)*(Q[1][idx_cell_mom]*Q[1][idx_cell_mom] + Q[2][idx_cell_mom]*Q[2][idx_cell_mom] +
                                    Q[3][idx_cell_mom]*Q[3][idx_cell_mom])/
                                    Q[0][idx_cell_rho];
                                
                                Q[4][idx_cell_E] = E;
                            }
                        }
                    }
                    
                    // Remove edge locations that have boundary conditions identified.
                    bdry_edge_locs.erase(std::remove(bdry_edge_locs.begin(), bdry_edge_locs.end(), edge_loc),
                        bdry_edge_locs.end());
                }
                else if (bdry_edge_conds[edge_loc] == BDRY_COND::FLOW_MODEL::XISOTHERMAL_NO_SLIP)
                {
                    for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                    {
                        for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                        {
                            for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                            {
                                const int idx_cell_rho = (i + num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0] +
                                    (k + num_subghosts_conservative_var[0][2])*
                                        subghostcell_dims_conservative_var[0][0]*
                                            subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom = (i + num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0] +
                                    (k + num_subghosts_conservative_var[1][2])*
                                        subghostcell_dims_conservative_var[1][0]*
                                            subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E = (i + num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0] +
                                    (k + num_subghosts_conservative_var[2][2])*
                                        subghostcell_dims_conservative_var[2][0]*
                                            subghostcell_dims_conservative_var[2][1];
                                
                                int idx_cell_pivot_rho = idx_cell_rho;
                                int idx_cell_pivot_mom = idx_cell_mom;
                                int idx_cell_pivot_E = idx_cell_E;
                                
                                if (face_loc_0 == BDRY_LOC::XLO)
                                {
                                    idx_cell_pivot_rho = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                            num_subghosts_conservative_var[0][0]) +
                                        (j + num_subghosts_conservative_var[0][1])*
                                            subghostcell_dims_conservative_var[0][0] +
                                        (k + num_subghosts_conservative_var[0][2])*
                                            subghostcell_dims_conservative_var[0][0]*
                                                subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                            num_subghosts_conservative_var[1][0]) +
                                        (j + num_subghosts_conservative_var[1][1])*
                                            subghostcell_dims_conservative_var[1][0] +
                                        (k + num_subghosts_conservative_var[1][2])*
                                            subghostcell_dims_conservative_var[1][0]*
                                                subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                            num_subghosts_conservative_var[2][0]) +
                                        (j + num_subghosts_conservative_var[2][1])*
                                            subghostcell_dims_conservative_var[2][0] +
                                        (k + num_subghosts_conservative_var[2][2])*
                                            subghostcell_dims_conservative_var[2][0]*
                                                subghostcell_dims_conservative_var[2][1];
                                }
                                else if (face_loc_0 == BDRY_LOC::XHI)
                                {
                                    idx_cell_pivot_rho = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                            num_subghosts_conservative_var[0][0]) +
                                        (j + num_subghosts_conservative_var[0][1])*
                                            subghostcell_dims_conservative_var[0][0] +
                                        (k + num_subghosts_conservative_var[0][2])*
                                            subghostcell_dims_conservative_var[0][0]*
                                                subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                            num_subghosts_conservative_var[1][0]) +
                                        (j + num_subghosts_conservative_var[1][1])*
                                            subghostcell_dims_conservative_var[1][0] +
                                        (k + num_subghosts_conservative_var[1][2])*
                                            subghostcell_dims_conservative_var[1][0]*
                                                subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                            num_subghosts_conservative_var[2][0]) +
                                        (j + num_subghosts_conservative_var[2][1])*
                                            subghostcell_dims_conservative_var[2][0] +
                                        (k + num_subghosts_conservative_var[2][2])*
                                            subghostcell_dims_conservative_var[2][0]*
                                                subghostcell_dims_conservative_var[2][1];
                                }
                                
                                /*
                                 * Set the values for density, momentum and total internal energy.
                                 */
                                
                                Real epsilon_pivot = (Q[4][idx_cell_pivot_E] -
                                    Real(1)/Real(2)*(Q[1][idx_cell_pivot_mom]*Q[1][idx_cell_pivot_mom] +
                                    Q[2][idx_cell_pivot_mom]*Q[2][idx_cell_pivot_mom] +
                                    Q[3][idx_cell_pivot_mom]*Q[3][idx_cell_pivot_mom])/
                                    Q[0][idx_cell_pivot_rho])/Q[0][idx_cell_pivot_rho];
                                
                                Real p_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &Q[0][idx_cell_pivot_rho],
                                        &epsilon_pivot,
                                        thermo_properties_ptr);
                                
                                Real p = p_pivot;
                                
                                Real T_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getTemperature(
                                        &Q[0][idx_cell_pivot_rho],
                                        &p_pivot,
                                        thermo_properties_ptr);
                                
                                Real T = -T_pivot + Real(2)*d_bdry_face_isothermal_no_slip_T[face_loc_0];
                                
                                Real rho = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getDensity(
                                        &p,
                                        &T,
                                        thermo_properties_ptr);
                                
                                Real u = -Q[1][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                                    Real(2)*d_bdry_edge_isothermal_no_slip_vel[face_loc_0*3];
                                Real v = -Q[2][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                                    Real(2)*d_bdry_edge_isothermal_no_slip_vel[face_loc_0*3 + 1];
                                Real w = -Q[3][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                                    Real(2)*d_bdry_edge_isothermal_no_slip_vel[face_loc_0*3 + 2];
                                
                                Q[0][idx_cell_rho] = rho;
                                Q[1][idx_cell_mom] = rho*u;
                                Q[2][idx_cell_mom] = rho*v;
                                Q[3][idx_cell_mom] = rho*w;
                                
                                Real epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getInternalEnergyFromTemperature(
                                        &rho,
                                        &T,
                                        thermo_properties_ptr);
                                
                                Real E = rho*epsilon +
                                    Real(1)/Real(2)*(Q[1][idx_cell_mom]*Q[1][idx_cell_mom] + Q[2][idx_cell_mom]*Q[2][idx_cell_mom] +
                                    Q[3][idx_cell_mom]*Q[3][idx_cell_mom])/rho;
                                
                                Q[4][idx_cell_E] = E;
                            }
                        }
                    }
                    
                    // Remove edge locations that have boundary conditions identified.
                    bdry_edge_locs.erase(std::remove(bdry_edge_locs.begin(), bdry_edge_locs.end(), edge_loc),
                        bdry_edge_locs.end());
                }
                else if (bdry_edge_conds[edge_loc] == BDRY_COND::FLOW_MODEL::YISOTHERMAL_NO_SLIP)
                {
                    for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                    {
                        for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                        {
                            for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                            {
                                const int idx_cell_rho = (i + num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0] +
                                    (k + num_subghosts_conservative_var[0][2])*
                                        subghostcell_dims_conservative_var[0][0]*
                                            subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom = (i + num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0] +
                                    (k + num_subghosts_conservative_var[1][2])*
                                        subghostcell_dims_conservative_var[1][0]*
                                            subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E = (i + num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0] +
                                    (k + num_subghosts_conservative_var[2][2])*
                                        subghostcell_dims_conservative_var[2][0]*
                                            subghostcell_dims_conservative_var[2][1];
                                
                                int idx_cell_pivot_rho = idx_cell_rho;
                                int idx_cell_pivot_mom = idx_cell_mom;
                                int idx_cell_pivot_E = idx_cell_E;
                                
                                if (face_loc_1 == BDRY_LOC::YLO)
                                {
                                    idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                        (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) +
                                            num_subghosts_conservative_var[0][1])*
                                                subghostcell_dims_conservative_var[0][0] +
                                        (k + num_subghosts_conservative_var[0][2])*
                                            subghostcell_dims_conservative_var[0][0]*
                                                subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                        (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) +
                                            num_subghosts_conservative_var[1][1])*
                                                subghostcell_dims_conservative_var[1][0] +
                                        (k + num_subghosts_conservative_var[1][2])*
                                            subghostcell_dims_conservative_var[1][0]*
                                                subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                        (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) +
                                            num_subghosts_conservative_var[2][1])*
                                                subghostcell_dims_conservative_var[2][0] +
                                        (k + num_subghosts_conservative_var[2][2])*
                                            subghostcell_dims_conservative_var[2][0]*
                                                subghostcell_dims_conservative_var[2][1];
                                }
                                else if (face_loc_1 == BDRY_LOC::YHI)
                                {
                                    idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                        (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) +
                                            num_subghosts_conservative_var[0][1])*
                                                subghostcell_dims_conservative_var[0][0] +
                                        (k + num_subghosts_conservative_var[0][2])*
                                            subghostcell_dims_conservative_var[0][0]*
                                                subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                        (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) +
                                            num_subghosts_conservative_var[1][1])*
                                                subghostcell_dims_conservative_var[1][0] +
                                        (k + num_subghosts_conservative_var[1][2])*
                                            subghostcell_dims_conservative_var[1][0]*
                                                subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                        (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) +
                                            num_subghosts_conservative_var[2][1])*
                                                subghostcell_dims_conservative_var[2][0] +
                                        (k + num_subghosts_conservative_var[2][2])*
                                            subghostcell_dims_conservative_var[2][0]*
                                                subghostcell_dims_conservative_var[2][1];
                                }
                                
                                /*
                                 * Set the values for density, momentum and total internal energy.
                                 */
                                
                                Real epsilon_pivot = (Q[4][idx_cell_pivot_E] -
                                    Real(1)/Real(2)*(Q[1][idx_cell_pivot_mom]*Q[1][idx_cell_pivot_mom] +
                                    Q[2][idx_cell_pivot_mom]*Q[2][idx_cell_pivot_mom] +
                                    Q[3][idx_cell_pivot_mom]*Q[3][idx_cell_pivot_mom])/
                                    Q[0][idx_cell_pivot_rho])/Q[0][idx_cell_pivot_rho];
                                
                                Real p_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &Q[0][idx_cell_pivot_rho],
                                        &epsilon_pivot,
                                        thermo_properties_ptr);
                                
                                Real p = p_pivot;
                                
                                Real T_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getTemperature(
                                        &Q[0][idx_cell_pivot_rho],
                                        &p_pivot,
                                        thermo_properties_ptr);
                                
                                Real T = -T_pivot + Real(2)*d_bdry_face_isothermal_no_slip_T[face_loc_1];
                                
                                Real rho = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getDensity(
                                        &p,
                                        &T,
                                        thermo_properties_ptr);
                                
                                Real u = -Q[1][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                                    Real(2)*d_bdry_edge_isothermal_no_slip_vel[face_loc_1*3];
                                Real v = -Q[2][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                                    Real(2)*d_bdry_edge_isothermal_no_slip_vel[face_loc_1*3 + 1];
                                Real w = -Q[3][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                                    Real(2)*d_bdry_edge_isothermal_no_slip_vel[face_loc_1*3 + 2];
                                
                                Q[0][idx_cell_rho] = rho;
                                Q[1][idx_cell_mom] = rho*u;
                                Q[2][idx_cell_mom] = rho*v;
                                Q[3][idx_cell_mom] = rho*w;
                                
                                Real epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getInternalEnergyFromTemperature(
                                        &rho,
                                        &T,
                                        thermo_properties_ptr);
                                
                                Real E = rho*epsilon +
                                    Real(1)/Real(2)*(Q[1][idx_cell_mom]*Q[1][idx_cell_mom] + Q[2][idx_cell_mom]*Q[2][idx_cell_mom] +
                                    Q[3][idx_cell_mom]*Q[3][idx_cell_mom])/rho;
                                
                                Q[4][idx_cell_E] = E;
                            }
                        }
                    }
                    
                    // Remove edge locations that have boundary conditions identified.
                    bdry_edge_locs.erase(std::remove(bdry_edge_locs.begin(), bdry_edge_locs.end(), edge_loc),
                        bdry_edge_locs.end());
                }
                else if (bdry_edge_conds[edge_loc] == BDRY_COND::FLOW_MODEL::ZISOTHERMAL_NO_SLIP)
                {
                    for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                    {
                        for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                        {
                            for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                            {
                                const int idx_cell_rho = (i + num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0] +
                                    (k + num_subghosts_conservative_var[0][2])*
                                        subghostcell_dims_conservative_var[0][0]*
                                            subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom = (i + num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0] +
                                    (k + num_subghosts_conservative_var[1][2])*
                                        subghostcell_dims_conservative_var[1][0]*
                                            subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E = (i + num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0] +
                                    (k + num_subghosts_conservative_var[2][2])*
                                        subghostcell_dims_conservative_var[2][0]*
                                            subghostcell_dims_conservative_var[2][1];
                                
                                int idx_cell_pivot_rho = idx_cell_rho;
                                int idx_cell_pivot_mom = idx_cell_mom;
                                int idx_cell_pivot_E = idx_cell_E;
                                
                                if (face_loc_2 == BDRY_LOC::ZLO)
                                {
                                    idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                        (j + num_subghosts_conservative_var[0][1])*
                                            subghostcell_dims_conservative_var[0][0] +
                                        (interior_box_lo_idx[2] + (fill_box_hi_idx[2] - k) +
                                            num_subghosts_conservative_var[0][2])*
                                                subghostcell_dims_conservative_var[0][0]*
                                                    subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                        (j + num_subghosts_conservative_var[1][1])*
                                            subghostcell_dims_conservative_var[1][0] +
                                        (interior_box_lo_idx[2] + (fill_box_hi_idx[2] - k) +
                                            num_subghosts_conservative_var[1][2])*
                                                subghostcell_dims_conservative_var[1][0]*
                                                    subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                        (j + num_subghosts_conservative_var[2][1])*
                                            subghostcell_dims_conservative_var[2][0] +
                                        (interior_box_lo_idx[2] + (fill_box_hi_idx[2] - k) +
                                            num_subghosts_conservative_var[2][2])*
                                                subghostcell_dims_conservative_var[2][0]*
                                                    subghostcell_dims_conservative_var[2][1];
                                }
                                else if (face_loc_2 == BDRY_LOC::ZHI)
                                {
                                    idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                        (j + num_subghosts_conservative_var[0][1])*
                                            subghostcell_dims_conservative_var[0][0] +
                                        (interior_box_hi_idx[2] - (k - fill_box_lo_idx[2]) +
                                            num_subghosts_conservative_var[0][2])*
                                                subghostcell_dims_conservative_var[0][0]*
                                                    subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                        (j + num_subghosts_conservative_var[1][1])*
                                            subghostcell_dims_conservative_var[1][0] +
                                        (interior_box_hi_idx[2] - (k - fill_box_lo_idx[2]) +
                                            num_subghosts_conservative_var[1][2])*
                                                subghostcell_dims_conservative_var[1][0]*
                                                    subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                        (j + num_subghosts_conservative_var[2][1])*
                                            subghostcell_dims_conservative_var[2][0] +
                                        (interior_box_hi_idx[2] - (k - fill_box_lo_idx[2]) +
                                            num_subghosts_conservative_var[2][2])*
                                                subghostcell_dims_conservative_var[2][0]*
                                                    subghostcell_dims_conservative_var[2][1];
                                }
                                
                                /*
                                 * Set the values for density, momentum and total internal energy.
                                 */
                                
                                Real epsilon_pivot = (Q[4][idx_cell_pivot_E] -
                                    Real(1)/Real(2)*(Q[1][idx_cell_pivot_mom]*Q[1][idx_cell_pivot_mom] +
                                    Q[2][idx_cell_pivot_mom]*Q[2][idx_cell_pivot_mom] +
                                    Q[3][idx_cell_pivot_mom]*Q[3][idx_cell_pivot_mom])/
                                    Q[0][idx_cell_pivot_rho])/Q[0][idx_cell_pivot_rho];
                                
                                Real p_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &Q[0][idx_cell_pivot_rho],
                                        &epsilon_pivot,
                                        thermo_properties_ptr);
                                
                                Real p = p_pivot;
                                
                                Real T_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getTemperature(
                                        &Q[0][idx_cell_pivot_rho],
                                        &p_pivot,
                                        thermo_properties_ptr);
                                
                                Real T = -T_pivot + Real(2)*d_bdry_face_isothermal_no_slip_T[face_loc_2];
                                
                                Real rho = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getDensity(
                                        &p,
                                        &T,
                                        thermo_properties_ptr);
                                
                                Real u = -Q[1][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                                    Real(2)*d_bdry_edge_isothermal_no_slip_vel[face_loc_2*3];
                                Real v = -Q[2][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                                    Real(2)*d_bdry_edge_isothermal_no_slip_vel[face_loc_2*3 + 1];
                                Real w = -Q[3][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                                    Real(2)*d_bdry_edge_isothermal_no_slip_vel[face_loc_2*3 + 2];
                                
                                Q[0][idx_cell_rho] = rho;
                                Q[1][idx_cell_mom] = rho*u;
                                Q[2][idx_cell_mom] = rho*v;
                                Q[3][idx_cell_mom] = rho*w;
                                
                                Real epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getInternalEnergyFromTemperature(
                                        &rho,
                                        &T,
                                        thermo_properties_ptr);
                                
                                Real E = rho*epsilon +
                                    Real(1)/Real(2)*(Q[1][idx_cell_mom]*Q[1][idx_cell_mom] + Q[2][idx_cell_mom]*Q[2][idx_cell_mom] +
                                    Q[3][idx_cell_mom]*Q[3][idx_cell_mom])/rho;
                                
                                Q[4][idx_cell_E] = E;
                            }
                        }
                    }
                    
                    // Remove edge locations that have boundary conditions identified.
                    bdry_edge_locs.erase(std::remove(bdry_edge_locs.begin(), bdry_edge_locs.end(), edge_loc),
                        bdry_edge_locs.end());
                }
            }
        }
    }
    
    for (int ei = 0; ei < static_cast<int>(edge_bdry.size()); ei++)
    {
        TBOX_ASSERT(edge_bdry[ei].getBoundaryType() == BDRY::EDGE3D);
        
        int edge_loc(edge_bdry[ei].getLocationIndex());
        
        if (std::find(bdry_edge_locs.begin(), bdry_edge_locs.end(), edge_loc) !=
            bdry_edge_locs.end())
        {
            TBOX_ERROR("FlowModelBoundaryUtilitiesSingleSpecies::fill3dEdgeBoundaryData()\n"
                << "Invalid edge boundary condition!\n"
                << "edge_loc = '" << edge_loc << "'." << std::endl
                << "bdry_edge_conds[edge_loc] = '" << bdry_edge_conds[edge_loc] << "'."
                << std::endl);
        }
    }
}


/*
 * Function to fill 3d node boundary values for a patch.
 * Node locations that have boundary conditions identified are removed from the container.
 */
void
FlowModelBoundaryUtilitiesSingleSpecies::fill3dNodeBoundaryData(
    const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& conservative_var_data,
    const hier::Patch& patch,
    std::vector<int>& bdry_node_locs,
    const std::vector<int>& bdry_node_conds,
    const std::vector<std::vector<Real> >& bdry_face_values,
    const hier::IntVector& ghost_width_to_fill)
{
    TBOX_ASSERT(static_cast<int>(conservative_var_data.size()) == 3);
    for (int vi = 0; vi < static_cast<int>(conservative_var_data.size()); vi++)
    {
        TBOX_ASSERT(conservative_var_data[vi]);
    }
    TBOX_ASSERT(static_cast<int>(bdry_node_locs.size()) <= NUM_3D_NODES);
    if (static_cast<int>(bdry_node_locs.size()) > 0)
    {
        TBOX_ASSERT(*min_element(bdry_node_locs.begin(), bdry_node_locs.end()) >= 0);
        TBOX_ASSERT(*max_element(bdry_node_locs.begin(), bdry_node_locs.end()) < NUM_3D_NODES);
    }
    TBOX_ASSERT(static_cast<int>(bdry_node_conds.size()) == NUM_3D_NODES);
    TBOX_ASSERT(static_cast<int>(bdry_face_values.size()) == 3);
    for (int vi = 0; vi < static_cast<int>(bdry_face_values.size()); vi++)
    {
        TBOX_ASSERT(static_cast<int>(bdry_face_values[vi].size()) ==
                    NUM_3D_FACES*(conservative_var_data[vi]->getDepth()));
    }
    
    TBOX_DIM_ASSERT(ghost_width_to_fill.getDim() == tbox::Dimension(3));
    
    for (int vi = 0; vi < static_cast<int>(conservative_var_data.size()); vi++)
    {
        TBOX_ASSERT_OBJDIM_EQUALITY3(*conservative_var_data[vi], patch, ghost_width_to_fill);
    }
    
    const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
        HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    TBOX_ASSERT(patch_geom);
    
    hier::IntVector num_ghosts(conservative_var_data[0]->getGhostCellWidth());
    for (int vi = 1; vi < static_cast<int>(conservative_var_data.size()); vi++)
    {
        num_ghosts = hier::IntVector::min(
            num_ghosts,
            conservative_var_data[vi]->getGhostCellWidth());
    }
    
    /*
     * Determine the ghost cell width to fill.
     */
    
    hier::IntVector gcw_to_fill(tbox::Dimension(3));
    
    // If the ghost fill width is not used, it is set to the ghost cell width of the data.
    if (ghost_width_to_fill == -hier::IntVector::getOne(tbox::Dimension(3)))
    {
        gcw_to_fill = num_ghosts;
    }
    else
    {
        gcw_to_fill = hier::IntVector::min(
            num_ghosts,
            ghost_width_to_fill);
    }
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box& interior_box(patch.getBox());
    hier::Index interior_box_lo_idx(interior_box.lower());
    hier::Index interior_box_hi_idx(interior_box.upper());
    
    /*
     * Offset the indices.
     */
    
    interior_box_lo_idx = interior_box_lo_idx - interior_box.lower();
    interior_box_hi_idx = interior_box_hi_idx - interior_box.lower();
    
    const std::vector<hier::BoundaryBox>& node_bdry =
        patch_geom->getCodimensionBoundaries(BDRY::NODE3D);
    
    for (int ni = 0; ni < static_cast<int>(node_bdry.size()); ni++)
    {
        TBOX_ASSERT(node_bdry[ni].getBoundaryType() == BDRY::NODE3D);
        
        int node_loc(node_bdry[ni].getLocationIndex());
        
        if (std::find(bdry_node_locs.begin(), bdry_node_locs.end(), node_loc) !=
            bdry_node_locs.end())
        {
            hier::Box fill_box(patch_geom->getBoundaryFillBox(
                node_bdry[ni],
                interior_box,
                gcw_to_fill));
            
            hier::Index fill_box_lo_idx(fill_box.lower());
            hier::Index fill_box_hi_idx(fill_box.upper());
            
            /*
             * Offset the indices.
             */
            
            fill_box_lo_idx = fill_box_lo_idx - interior_box.lower();
            fill_box_hi_idx = fill_box_hi_idx - interior_box.lower();
            
            int face_loc_0 = -1;
            int face_loc_1 = -1;
            int face_loc_2 = -1;
            
            switch (node_loc)
            {
                case NODE_BDRY_LOC_3D::XLO_YLO_ZLO:
                {
                    face_loc_0 = BDRY_LOC::XLO;
                    face_loc_1 = BDRY_LOC::YLO;
                    face_loc_2 = BDRY_LOC::ZLO;
                    
                    break;
                }
                case NODE_BDRY_LOC_3D::XHI_YLO_ZLO:
                {
                    face_loc_0 = BDRY_LOC::XHI;
                    face_loc_1 = BDRY_LOC::YLO;
                    face_loc_2 = BDRY_LOC::ZLO;
                    
                    break;
                }
                case NODE_BDRY_LOC_3D::XLO_YHI_ZLO:
                {
                    face_loc_0 = BDRY_LOC::XLO;
                    face_loc_1 = BDRY_LOC::YHI;
                    face_loc_2 = BDRY_LOC::ZLO;
                    
                    break;
                }
                case NODE_BDRY_LOC_3D::XHI_YHI_ZLO:
                {
                    face_loc_0 = BDRY_LOC::XHI;
                    face_loc_1 = BDRY_LOC::YHI;
                    face_loc_2 = BDRY_LOC::ZLO;
                    
                    break;
                }
                case NODE_BDRY_LOC_3D::XLO_YLO_ZHI:
                {
                    face_loc_0 = BDRY_LOC::XLO;
                    face_loc_1 = BDRY_LOC::YLO;
                    face_loc_2 = BDRY_LOC::ZHI;
                    
                    break;
                }
                case NODE_BDRY_LOC_3D::XHI_YLO_ZHI:
                {
                    face_loc_0 = BDRY_LOC::XHI;
                    face_loc_1 = BDRY_LOC::YLO;
                    face_loc_2 = BDRY_LOC::ZHI;
                    
                    break;
                }
                case NODE_BDRY_LOC_3D::XLO_YHI_ZHI:
                {
                    face_loc_0 = BDRY_LOC::XLO;
                    face_loc_1 = BDRY_LOC::YHI;
                    face_loc_2 = BDRY_LOC::ZHI;
                    
                    break;
                }
                case NODE_BDRY_LOC_3D::XHI_YHI_ZHI:
                {
                    face_loc_0 = BDRY_LOC::XHI;
                    face_loc_1 = BDRY_LOC::YHI;
                    face_loc_2 = BDRY_LOC::ZHI;
                    
                    break;
                }
            }
            
            if ((bdry_node_conds[node_loc] == BDRY_COND::FLOW_MODEL::XADIABATIC_NO_SLIP) ||
                (bdry_node_conds[node_loc] == BDRY_COND::FLOW_MODEL::YADIABATIC_NO_SLIP) ||
                (bdry_node_conds[node_loc] == BDRY_COND::FLOW_MODEL::ZADIABATIC_NO_SLIP) ||
                (bdry_node_conds[node_loc] == BDRY_COND::FLOW_MODEL::XISOTHERMAL_NO_SLIP) ||
                (bdry_node_conds[node_loc] == BDRY_COND::FLOW_MODEL::YISOTHERMAL_NO_SLIP) ||
                (bdry_node_conds[node_loc] == BDRY_COND::FLOW_MODEL::ZISOTHERMAL_NO_SLIP))
            {
                /*
                 * Get the pointers to the conservative variables.
                 * Get the numbers of ghost cells and the dimensions of the ghost cell boxes of
                 * the conservative variables.
                 */
                
                std::vector<Real*> Q;
                Q.reserve(d_num_eqn);
                
                std::vector<hier::IntVector> num_subghosts_conservative_var;
                std::vector<hier::IntVector> subghostcell_dims_conservative_var;
                
                num_subghosts_conservative_var.reserve(3);
                subghostcell_dims_conservative_var.reserve(3);
                
                int count_eqn = 0;
                
                for (int vi = 0; vi < static_cast<int>(conservative_var_data.size()); vi++)
                {
                    int depth = conservative_var_data[vi]->getDepth();
                    
                    for (int di = 0; di < depth; di++)
                    {
                        // If the last element of the conservative variable vector is not in the
                        // system of equations, ignore it.
                        if (count_eqn >= d_num_eqn)
                            break;
                        
                        Q.push_back(conservative_var_data[vi]->getPointer(di));
                        
                        count_eqn++;
                    }
                    
                    num_subghosts_conservative_var.push_back(
                        conservative_var_data[vi]->getGhostCellWidth());
                    subghostcell_dims_conservative_var.push_back(
                        conservative_var_data[vi]->getGhostBox().numberCells());
                }
                
                // Get the thermodynamic properties of the species.
                std::vector<const Real*> thermo_properties_ptr;
                thermo_properties_ptr.reserve(static_cast<int> (d_thermo_properties.size()));
                for (int ti = 0; ti < static_cast<int> (d_thermo_properties.size()); ti++)
                {
                    thermo_properties_ptr.push_back(&d_thermo_properties[ti]);
                }
                
                if (bdry_node_conds[node_loc] == BDRY_COND::FLOW_MODEL::XADIABATIC_NO_SLIP)
                {
                    for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                    {
                        for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                        {
                            for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                            {
                                const int idx_cell_rho = (i + num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0] +
                                    (k + num_subghosts_conservative_var[0][2])*
                                        subghostcell_dims_conservative_var[0][0]*
                                            subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom = (i + num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0] +
                                    (k + num_subghosts_conservative_var[1][2])*
                                        subghostcell_dims_conservative_var[1][0]*
                                            subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E = (i + num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0] +
                                    (k + num_subghosts_conservative_var[2][2])*
                                        subghostcell_dims_conservative_var[2][0]*
                                            subghostcell_dims_conservative_var[2][1];
                                
                                int idx_cell_pivot_rho = idx_cell_rho;
                                int idx_cell_pivot_mom = idx_cell_mom;
                                int idx_cell_pivot_E = idx_cell_E;
                                
                                if (face_loc_0 == BDRY_LOC::XLO)
                                {
                                    idx_cell_pivot_rho = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                            num_subghosts_conservative_var[0][0]) +
                                        (j + num_subghosts_conservative_var[0][1])*
                                            subghostcell_dims_conservative_var[0][0] +
                                        (k + num_subghosts_conservative_var[0][2])*
                                            subghostcell_dims_conservative_var[0][0]*
                                                subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                            num_subghosts_conservative_var[1][0]) +
                                        (j + num_subghosts_conservative_var[1][1])*
                                            subghostcell_dims_conservative_var[1][0] +
                                        (k + num_subghosts_conservative_var[1][2])*
                                            subghostcell_dims_conservative_var[1][0]*
                                                subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                            num_subghosts_conservative_var[2][0]) +
                                        (j + num_subghosts_conservative_var[2][1])*
                                            subghostcell_dims_conservative_var[2][0] +
                                        (k + num_subghosts_conservative_var[2][2])*
                                            subghostcell_dims_conservative_var[2][0]*
                                                subghostcell_dims_conservative_var[2][1];
                                }
                                else if (face_loc_0 == BDRY_LOC::XHI)
                                {
                                    idx_cell_pivot_rho = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                            num_subghosts_conservative_var[0][0]) +
                                        (j + num_subghosts_conservative_var[0][1])*
                                            subghostcell_dims_conservative_var[0][0] +
                                        (k + num_subghosts_conservative_var[0][2])*
                                            subghostcell_dims_conservative_var[0][0]*
                                                subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                            num_subghosts_conservative_var[1][0]) +
                                        (j + num_subghosts_conservative_var[1][1])*
                                            subghostcell_dims_conservative_var[1][0] +
                                        (k + num_subghosts_conservative_var[1][2])*
                                            subghostcell_dims_conservative_var[1][0]*
                                                subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                            num_subghosts_conservative_var[2][0]) +
                                        (j + num_subghosts_conservative_var[2][1])*
                                            subghostcell_dims_conservative_var[2][0] +
                                        (k + num_subghosts_conservative_var[2][2])*
                                            subghostcell_dims_conservative_var[2][0]*
                                                subghostcell_dims_conservative_var[2][1];
                                }
                                
                                /*
                                 * Set the values for density and momentum.
                                 */
                                
                                Q[0][idx_cell_rho] = Q[0][idx_cell_pivot_rho];
                                Q[1][idx_cell_mom] = -Q[1][idx_cell_pivot_mom] +
                                    Real(2)*Q[0][idx_cell_pivot_rho]*d_bdry_face_adiabatic_no_slip_vel[face_loc_0*3];
                                Q[2][idx_cell_mom] = -Q[2][idx_cell_pivot_mom] +
                                    Real(2)*Q[0][idx_cell_pivot_rho]*d_bdry_face_adiabatic_no_slip_vel[face_loc_0*3 + 1];
                                Q[3][idx_cell_mom] = -Q[3][idx_cell_pivot_mom] +
                                    Real(2)*Q[0][idx_cell_pivot_rho]*d_bdry_face_adiabatic_no_slip_vel[face_loc_0*3 + 2];
                                
                                /*
                                 * Set the values for total internal energy.
                                 */
                                
                                Real epsilon_pivot = (Q[4][idx_cell_pivot_E] -
                                    Real(1)/Real(2)*(Q[1][idx_cell_pivot_mom]*Q[1][idx_cell_pivot_mom] +
                                    Q[2][idx_cell_pivot_mom]*Q[2][idx_cell_pivot_mom] +
                                    Q[3][idx_cell_pivot_mom]*Q[3][idx_cell_pivot_mom])/
                                    Q[0][idx_cell_pivot_rho])/Q[0][idx_cell_pivot_rho];
                                
                                Real p_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &Q[0][idx_cell_pivot_rho],
                                        &epsilon_pivot,
                                        thermo_properties_ptr);
                                
                                Real T_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getTemperature(
                                        &Q[0][idx_cell_pivot_rho],
                                        &p_pivot,
                                        thermo_properties_ptr);
                                
                                Real T = T_pivot;
                                
                                Real epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getInternalEnergyFromTemperature(
                                        &Q[0][idx_cell_rho],
                                        &T,
                                        thermo_properties_ptr);
                                
                                Real E = Q[0][idx_cell_rho]*epsilon +
                                    Real(1)/Real(2)*(Q[1][idx_cell_mom]*Q[1][idx_cell_mom] + Q[2][idx_cell_mom]*Q[2][idx_cell_mom] +
                                    Q[3][idx_cell_mom]*Q[3][idx_cell_mom])/
                                    Q[0][idx_cell_rho];
                                
                                Q[4][idx_cell_E] = E;
                            }
                        }
                    }
                    
                    // Remove node locations that have boundary conditions identified.
                    bdry_node_locs.erase(std::remove(bdry_node_locs.begin(), bdry_node_locs.end(), node_loc),
                        bdry_node_locs.end());
                }
                else if (bdry_node_conds[node_loc] == BDRY_COND::FLOW_MODEL::YADIABATIC_NO_SLIP)
                {
                    for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                    {
                        for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                        {
                            for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                            {
                                const int idx_cell_rho = (i + num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0] +
                                    (k + num_subghosts_conservative_var[0][2])*
                                        subghostcell_dims_conservative_var[0][0]*
                                            subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom = (i + num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0] +
                                    (k + num_subghosts_conservative_var[1][2])*
                                        subghostcell_dims_conservative_var[1][0]*
                                            subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E = (i + num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0] +
                                    (k + num_subghosts_conservative_var[2][2])*
                                        subghostcell_dims_conservative_var[2][0]*
                                            subghostcell_dims_conservative_var[2][1];
                                
                                int idx_cell_pivot_rho = idx_cell_rho;
                                int idx_cell_pivot_mom = idx_cell_mom;
                                int idx_cell_pivot_E = idx_cell_E;
                                
                                if (face_loc_1 == BDRY_LOC::YLO)
                                {
                                    idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                        (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) +
                                            num_subghosts_conservative_var[0][1])*
                                                subghostcell_dims_conservative_var[0][0] +
                                        (k + num_subghosts_conservative_var[0][2])*
                                            subghostcell_dims_conservative_var[0][0]*
                                                subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                        (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) +
                                            num_subghosts_conservative_var[1][1])*
                                                subghostcell_dims_conservative_var[1][0] +
                                        (k + num_subghosts_conservative_var[1][2])*
                                            subghostcell_dims_conservative_var[1][0]*
                                                subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                        (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) +
                                            num_subghosts_conservative_var[2][1])*
                                                subghostcell_dims_conservative_var[2][0] +
                                        (k + num_subghosts_conservative_var[2][2])*
                                            subghostcell_dims_conservative_var[2][0]*
                                                subghostcell_dims_conservative_var[2][1];
                                }
                                else if (face_loc_1 == BDRY_LOC::YHI)
                                {
                                    idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                        (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) +
                                            num_subghosts_conservative_var[0][1])*
                                                subghostcell_dims_conservative_var[0][0] +
                                        (k + num_subghosts_conservative_var[0][2])*
                                            subghostcell_dims_conservative_var[0][0]*
                                                subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                        (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) +
                                            num_subghosts_conservative_var[1][1])*
                                                subghostcell_dims_conservative_var[1][0] +
                                        (k + num_subghosts_conservative_var[1][2])*
                                            subghostcell_dims_conservative_var[1][0]*
                                                subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                        (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) +
                                            num_subghosts_conservative_var[2][1])*
                                                subghostcell_dims_conservative_var[2][0] +
                                        (k + num_subghosts_conservative_var[2][2])*
                                            subghostcell_dims_conservative_var[2][0]*
                                                subghostcell_dims_conservative_var[2][1];
                                }
                                
                                /*
                                 * Set the values for density and momentum.
                                 */
                                
                                Q[0][idx_cell_rho] = Q[0][idx_cell_pivot_rho];
                                Q[1][idx_cell_mom] = -Q[1][idx_cell_pivot_mom] +
                                    Real(2)*Q[0][idx_cell_pivot_rho]*d_bdry_face_adiabatic_no_slip_vel[face_loc_1*3];
                                Q[2][idx_cell_mom] = -Q[2][idx_cell_pivot_mom] +
                                    Real(2)*Q[0][idx_cell_pivot_rho]*d_bdry_face_adiabatic_no_slip_vel[face_loc_1*3 + 1];
                                Q[3][idx_cell_mom] = -Q[3][idx_cell_pivot_mom] +
                                    Real(2)*Q[0][idx_cell_pivot_rho]*d_bdry_face_adiabatic_no_slip_vel[face_loc_1*3 + 2];
                                
                                /*
                                 * Set the values for total internal energy.
                                 */
                                
                                Real epsilon_pivot = (Q[4][idx_cell_pivot_E] -
                                    Real(1)/Real(2)*(Q[1][idx_cell_pivot_mom]*Q[1][idx_cell_pivot_mom] +
                                    Q[2][idx_cell_pivot_mom]*Q[2][idx_cell_pivot_mom] +
                                    Q[3][idx_cell_pivot_mom]*Q[3][idx_cell_pivot_mom])/
                                    Q[0][idx_cell_pivot_rho])/Q[0][idx_cell_pivot_rho];
                                
                                Real p_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &Q[0][idx_cell_pivot_rho],
                                        &epsilon_pivot,
                                        thermo_properties_ptr);
                                
                                Real T_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getTemperature(
                                        &Q[0][idx_cell_pivot_rho],
                                        &p_pivot,
                                        thermo_properties_ptr);
                                
                                Real T = T_pivot;
                                
                                Real epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getInternalEnergyFromTemperature(
                                        &Q[0][idx_cell_rho],
                                        &T,
                                        thermo_properties_ptr);
                                
                                Real E = Q[0][idx_cell_rho]*epsilon +
                                    Real(1)/Real(2)*(Q[1][idx_cell_mom]*Q[1][idx_cell_mom] + Q[2][idx_cell_mom]*Q[2][idx_cell_mom] +
                                    Q[3][idx_cell_mom]*Q[3][idx_cell_mom])/
                                    Q[0][idx_cell_rho];
                                
                                Q[4][idx_cell_E] = E;
                            }
                        }
                    }
                    
                    // Remove node locations that have boundary conditions identified.
                    bdry_node_locs.erase(std::remove(bdry_node_locs.begin(), bdry_node_locs.end(), node_loc),
                        bdry_node_locs.end());
                }
                else if (bdry_node_conds[node_loc] == BDRY_COND::FLOW_MODEL::ZADIABATIC_NO_SLIP)
                {
                    for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                    {
                        for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                        {
                            for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                            {
                                const int idx_cell_rho = (i + num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0] +
                                    (k + num_subghosts_conservative_var[0][2])*
                                        subghostcell_dims_conservative_var[0][0]*
                                            subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom = (i + num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0] +
                                    (k + num_subghosts_conservative_var[1][2])*
                                        subghostcell_dims_conservative_var[1][0]*
                                            subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E = (i + num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0] +
                                    (k + num_subghosts_conservative_var[2][2])*
                                        subghostcell_dims_conservative_var[2][0]*
                                            subghostcell_dims_conservative_var[2][1];
                                
                                int idx_cell_pivot_rho = idx_cell_rho;
                                int idx_cell_pivot_mom = idx_cell_mom;
                                int idx_cell_pivot_E = idx_cell_E;
                                
                                if (face_loc_2 == BDRY_LOC::ZLO)
                                {
                                    idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                        (j + num_subghosts_conservative_var[0][1])*
                                            subghostcell_dims_conservative_var[0][0] +
                                        (interior_box_lo_idx[2] + (fill_box_hi_idx[2] - k) +
                                            num_subghosts_conservative_var[0][2])*
                                                subghostcell_dims_conservative_var[0][0]*
                                                    subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                        (j + num_subghosts_conservative_var[1][1])*
                                            subghostcell_dims_conservative_var[1][0] +
                                        (interior_box_lo_idx[2] + (fill_box_hi_idx[2] - k) +
                                            num_subghosts_conservative_var[1][2])*
                                                subghostcell_dims_conservative_var[1][0]*
                                                    subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                        (j + num_subghosts_conservative_var[2][1])*
                                            subghostcell_dims_conservative_var[2][0] +
                                        (interior_box_lo_idx[2] + (fill_box_hi_idx[2] - k) +
                                            num_subghosts_conservative_var[2][2])*
                                                subghostcell_dims_conservative_var[2][0]*
                                                    subghostcell_dims_conservative_var[2][1];
                                }
                                else if (face_loc_2 == BDRY_LOC::ZHI)
                                {
                                    idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                        (j + num_subghosts_conservative_var[0][1])*
                                            subghostcell_dims_conservative_var[0][0] +
                                        (interior_box_hi_idx[2] - (k - fill_box_lo_idx[2]) +
                                            num_subghosts_conservative_var[0][2])*
                                                subghostcell_dims_conservative_var[0][0]*
                                                    subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                        (j + num_subghosts_conservative_var[1][1])*
                                            subghostcell_dims_conservative_var[1][0] +
                                        (interior_box_hi_idx[2] - (k - fill_box_lo_idx[2]) +
                                            num_subghosts_conservative_var[1][2])*
                                                subghostcell_dims_conservative_var[1][0]*
                                                    subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                        (j + num_subghosts_conservative_var[2][1])*
                                            subghostcell_dims_conservative_var[2][0] +
                                        (interior_box_hi_idx[2] - (k - fill_box_lo_idx[2]) +
                                            num_subghosts_conservative_var[2][2])*
                                                subghostcell_dims_conservative_var[2][0]*
                                                    subghostcell_dims_conservative_var[2][1];
                                }
                                
                                /*
                                 * Set the values for density and momentum.
                                 */
                                
                                Q[0][idx_cell_rho] = Q[0][idx_cell_pivot_rho];
                                Q[1][idx_cell_mom] = -Q[1][idx_cell_pivot_mom] +
                                    Real(2)*Q[0][idx_cell_pivot_rho]*d_bdry_face_adiabatic_no_slip_vel[face_loc_2*3];
                                Q[2][idx_cell_mom] = -Q[2][idx_cell_pivot_mom] +
                                    Real(2)*Q[0][idx_cell_pivot_rho]*d_bdry_face_adiabatic_no_slip_vel[face_loc_2*3 + 1];
                                Q[3][idx_cell_mom] = -Q[3][idx_cell_pivot_mom] +
                                    Real(2)*Q[0][idx_cell_pivot_rho]*d_bdry_face_adiabatic_no_slip_vel[face_loc_2*3 + 2];
                                
                                /*
                                 * Set the values for total internal energy.
                                 */
                                
                                Real epsilon_pivot = (Q[4][idx_cell_pivot_E] -
                                    Real(1)/Real(2)*(Q[1][idx_cell_pivot_mom]*Q[1][idx_cell_pivot_mom] +
                                    Q[2][idx_cell_pivot_mom]*Q[2][idx_cell_pivot_mom] +
                                    Q[3][idx_cell_pivot_mom]*Q[3][idx_cell_pivot_mom])/
                                    Q[0][idx_cell_pivot_rho])/Q[0][idx_cell_pivot_rho];
                                
                                Real p_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &Q[0][idx_cell_pivot_rho],
                                        &epsilon_pivot,
                                        thermo_properties_ptr);
                                
                                Real T_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getTemperature(
                                        &Q[0][idx_cell_pivot_rho],
                                        &p_pivot,
                                        thermo_properties_ptr);
                                
                                Real T = T_pivot;
                                
                                Real epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getInternalEnergyFromTemperature(
                                        &Q[0][idx_cell_rho],
                                        &T,
                                        thermo_properties_ptr);
                                
                                Real E = Q[0][idx_cell_rho]*epsilon +
                                    Real(1)/Real(2)*(Q[1][idx_cell_mom]*Q[1][idx_cell_mom] + Q[2][idx_cell_mom]*Q[2][idx_cell_mom] +
                                    Q[3][idx_cell_mom]*Q[3][idx_cell_mom])/
                                    Q[0][idx_cell_rho];
                                
                                Q[4][idx_cell_E] = E;
                            }
                        }
                    }
                    
                    // Remove node locations that have boundary conditions identified.
                    bdry_node_locs.erase(std::remove(bdry_node_locs.begin(), bdry_node_locs.end(), node_loc),
                        bdry_node_locs.end());
                }
                else if (bdry_node_conds[node_loc] == BDRY_COND::FLOW_MODEL::XISOTHERMAL_NO_SLIP)
                {
                    for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                    {
                        for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                        {
                            for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                            {
                                const int idx_cell_rho = (i + num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0] +
                                    (k + num_subghosts_conservative_var[0][2])*
                                        subghostcell_dims_conservative_var[0][0]*
                                            subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom = (i + num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0] +
                                    (k + num_subghosts_conservative_var[1][2])*
                                        subghostcell_dims_conservative_var[1][0]*
                                            subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E = (i + num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0] +
                                    (k + num_subghosts_conservative_var[2][2])*
                                        subghostcell_dims_conservative_var[2][0]*
                                            subghostcell_dims_conservative_var[2][1];
                                
                                int idx_cell_pivot_rho = idx_cell_rho;
                                int idx_cell_pivot_mom = idx_cell_mom;
                                int idx_cell_pivot_E = idx_cell_E;
                                
                                if (face_loc_0 == BDRY_LOC::XLO)
                                {
                                    idx_cell_pivot_rho = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                            num_subghosts_conservative_var[0][0]) +
                                        (j + num_subghosts_conservative_var[0][1])*
                                            subghostcell_dims_conservative_var[0][0] +
                                        (k + num_subghosts_conservative_var[0][2])*
                                            subghostcell_dims_conservative_var[0][0]*
                                                subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                            num_subghosts_conservative_var[1][0]) +
                                        (j + num_subghosts_conservative_var[1][1])*
                                            subghostcell_dims_conservative_var[1][0] +
                                        (k + num_subghosts_conservative_var[1][2])*
                                            subghostcell_dims_conservative_var[1][0]*
                                                subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                            num_subghosts_conservative_var[2][0]) +
                                        (j + num_subghosts_conservative_var[2][1])*
                                            subghostcell_dims_conservative_var[2][0] +
                                        (k + num_subghosts_conservative_var[2][2])*
                                            subghostcell_dims_conservative_var[2][0]*
                                                subghostcell_dims_conservative_var[2][1];
                                }
                                else if (face_loc_0 == BDRY_LOC::XHI)
                                {
                                    idx_cell_pivot_rho = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                            num_subghosts_conservative_var[0][0]) +
                                        (j + num_subghosts_conservative_var[0][1])*
                                            subghostcell_dims_conservative_var[0][0] +
                                        (k + num_subghosts_conservative_var[0][2])*
                                            subghostcell_dims_conservative_var[0][0]*
                                                subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                            num_subghosts_conservative_var[1][0]) +
                                        (j + num_subghosts_conservative_var[1][1])*
                                            subghostcell_dims_conservative_var[1][0] +
                                        (k + num_subghosts_conservative_var[1][2])*
                                            subghostcell_dims_conservative_var[1][0]*
                                                subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                            num_subghosts_conservative_var[2][0]) +
                                        (j + num_subghosts_conservative_var[2][1])*
                                            subghostcell_dims_conservative_var[2][0] +
                                        (k + num_subghosts_conservative_var[2][2])*
                                            subghostcell_dims_conservative_var[2][0]*
                                                subghostcell_dims_conservative_var[2][1];
                                }
                                
                                /*
                                 * Set the values for density, momentum and total internal energy.
                                 */
                                
                                Real epsilon_pivot = (Q[4][idx_cell_pivot_E] -
                                    Real(1)/Real(2)*(Q[1][idx_cell_pivot_mom]*Q[1][idx_cell_pivot_mom] +
                                    Q[2][idx_cell_pivot_mom]*Q[2][idx_cell_pivot_mom] +
                                    Q[3][idx_cell_pivot_mom]*Q[3][idx_cell_pivot_mom])/
                                    Q[0][idx_cell_pivot_rho])/Q[0][idx_cell_pivot_rho];
                                
                                Real p_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &Q[0][idx_cell_pivot_rho],
                                        &epsilon_pivot,
                                        thermo_properties_ptr);
                                
                                Real p = p_pivot;
                                
                                Real T_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getTemperature(
                                        &Q[0][idx_cell_pivot_rho],
                                        &p_pivot,
                                        thermo_properties_ptr);
                                
                                Real T = -T_pivot + Real(2)*d_bdry_face_isothermal_no_slip_T[face_loc_0];
                                
                                Real rho = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getDensity(
                                        &p,
                                        &T,
                                        thermo_properties_ptr);
                                
                                Real u = -Q[1][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                                    Real(2)*d_bdry_edge_isothermal_no_slip_vel[face_loc_0*3];
                                Real v = -Q[2][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                                    Real(2)*d_bdry_edge_isothermal_no_slip_vel[face_loc_0*3 + 1];
                                Real w = -Q[3][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                                    Real(2)*d_bdry_edge_isothermal_no_slip_vel[face_loc_0*3 + 2];
                                
                                Q[0][idx_cell_rho] = rho;
                                Q[1][idx_cell_mom] = rho*u;
                                Q[2][idx_cell_mom] = rho*v;
                                Q[3][idx_cell_mom] = rho*w;
                                
                                Real epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getInternalEnergyFromTemperature(
                                        &rho,
                                        &T,
                                        thermo_properties_ptr);
                                
                                Real E = rho*epsilon +
                                    Real(1)/Real(2)*(Q[1][idx_cell_mom]*Q[1][idx_cell_mom] + Q[2][idx_cell_mom]*Q[2][idx_cell_mom] +
                                    Q[3][idx_cell_mom]*Q[3][idx_cell_mom])/rho;
                                
                                Q[4][idx_cell_E] = E;
                            }
                        }
                    }
                    
                    // Remove node locations that have boundary conditions identified.
                    bdry_node_locs.erase(std::remove(bdry_node_locs.begin(), bdry_node_locs.end(), node_loc),
                        bdry_node_locs.end());
                }
                else if (bdry_node_conds[node_loc] == BDRY_COND::FLOW_MODEL::YISOTHERMAL_NO_SLIP)
                {
                    for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                    {
                        for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                        {
                            for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                            {
                                const int idx_cell_rho = (i + num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0] +
                                    (k + num_subghosts_conservative_var[0][2])*
                                        subghostcell_dims_conservative_var[0][0]*
                                            subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom = (i + num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0] +
                                    (k + num_subghosts_conservative_var[1][2])*
                                        subghostcell_dims_conservative_var[1][0]*
                                            subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E = (i + num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0] +
                                    (k + num_subghosts_conservative_var[2][2])*
                                        subghostcell_dims_conservative_var[2][0]*
                                            subghostcell_dims_conservative_var[2][1];
                                
                                int idx_cell_pivot_rho = idx_cell_rho;
                                int idx_cell_pivot_mom = idx_cell_mom;
                                int idx_cell_pivot_E = idx_cell_E;
                                
                                if (face_loc_1 == BDRY_LOC::YLO)
                                {
                                    idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                        (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) +
                                            num_subghosts_conservative_var[0][1])*
                                                subghostcell_dims_conservative_var[0][0] +
                                        (k + num_subghosts_conservative_var[0][2])*
                                            subghostcell_dims_conservative_var[0][0]*
                                                subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                        (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) +
                                            num_subghosts_conservative_var[1][1])*
                                                subghostcell_dims_conservative_var[1][0] +
                                        (k + num_subghosts_conservative_var[1][2])*
                                            subghostcell_dims_conservative_var[1][0]*
                                                subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                        (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) +
                                            num_subghosts_conservative_var[2][1])*
                                                subghostcell_dims_conservative_var[2][0] +
                                        (k + num_subghosts_conservative_var[2][2])*
                                            subghostcell_dims_conservative_var[2][0]*
                                                subghostcell_dims_conservative_var[2][1];
                                }
                                else if (face_loc_1 == BDRY_LOC::YHI)
                                {
                                    idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                        (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) +
                                            num_subghosts_conservative_var[0][1])*
                                                subghostcell_dims_conservative_var[0][0] +
                                        (k + num_subghosts_conservative_var[0][2])*
                                            subghostcell_dims_conservative_var[0][0]*
                                                subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                        (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) +
                                            num_subghosts_conservative_var[1][1])*
                                                subghostcell_dims_conservative_var[1][0] +
                                        (k + num_subghosts_conservative_var[1][2])*
                                            subghostcell_dims_conservative_var[1][0]*
                                                subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                        (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) +
                                            num_subghosts_conservative_var[2][1])*
                                                subghostcell_dims_conservative_var[2][0] +
                                        (k + num_subghosts_conservative_var[2][2])*
                                            subghostcell_dims_conservative_var[2][0]*
                                                subghostcell_dims_conservative_var[2][1];
                                }
                                
                                /*
                                 * Set the values for density, momentum and total internal energy.
                                 */
                                
                                Real epsilon_pivot = (Q[4][idx_cell_pivot_E] -
                                    Real(1)/Real(2)*(Q[1][idx_cell_pivot_mom]*Q[1][idx_cell_pivot_mom] +
                                    Q[2][idx_cell_pivot_mom]*Q[2][idx_cell_pivot_mom] +
                                    Q[3][idx_cell_pivot_mom]*Q[3][idx_cell_pivot_mom])/
                                    Q[0][idx_cell_pivot_rho])/Q[0][idx_cell_pivot_rho];
                                
                                Real p_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &Q[0][idx_cell_pivot_rho],
                                        &epsilon_pivot,
                                        thermo_properties_ptr);
                                
                                Real p = p_pivot;
                                
                                Real T_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getTemperature(
                                        &Q[0][idx_cell_pivot_rho],
                                        &p_pivot,
                                        thermo_properties_ptr);
                                
                                Real T = -T_pivot + Real(2)*d_bdry_face_isothermal_no_slip_T[face_loc_1];
                                
                                Real rho = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getDensity(
                                        &p,
                                        &T,
                                        thermo_properties_ptr);
                                
                                Real u = -Q[1][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                                    Real(2)*d_bdry_edge_isothermal_no_slip_vel[face_loc_1*3];
                                Real v = -Q[2][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                                    Real(2)*d_bdry_edge_isothermal_no_slip_vel[face_loc_1*3 + 1];
                                Real w = -Q[3][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                                    Real(2)*d_bdry_edge_isothermal_no_slip_vel[face_loc_1*3 + 2];
                                
                                Q[0][idx_cell_rho] = rho;
                                Q[1][idx_cell_mom] = rho*u;
                                Q[2][idx_cell_mom] = rho*v;
                                Q[3][idx_cell_mom] = rho*w;
                                
                                Real epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getInternalEnergyFromTemperature(
                                        &rho,
                                        &T,
                                        thermo_properties_ptr);
                                
                                Real E = rho*epsilon +
                                    Real(1)/Real(2)*(Q[1][idx_cell_mom]*Q[1][idx_cell_mom] + Q[2][idx_cell_mom]*Q[2][idx_cell_mom] +
                                    Q[3][idx_cell_mom]*Q[3][idx_cell_mom])/rho;
                                
                                Q[4][idx_cell_E] = E;
                            }
                        }
                    }
                    
                    // Remove node locations that have boundary conditions identified.
                    bdry_node_locs.erase(std::remove(bdry_node_locs.begin(), bdry_node_locs.end(), node_loc),
                        bdry_node_locs.end());
                }
                else if (bdry_node_conds[node_loc] == BDRY_COND::FLOW_MODEL::ZISOTHERMAL_NO_SLIP)
                {
                    for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                    {
                        for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                        {
                            for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                            {
                                const int idx_cell_rho = (i + num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0] +
                                    (k + num_subghosts_conservative_var[0][2])*
                                        subghostcell_dims_conservative_var[0][0]*
                                            subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom = (i + num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0] +
                                    (k + num_subghosts_conservative_var[1][2])*
                                        subghostcell_dims_conservative_var[1][0]*
                                            subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E = (i + num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0] +
                                    (k + num_subghosts_conservative_var[2][2])*
                                        subghostcell_dims_conservative_var[2][0]*
                                            subghostcell_dims_conservative_var[2][1];
                                
                                int idx_cell_pivot_rho = idx_cell_rho;
                                int idx_cell_pivot_mom = idx_cell_mom;
                                int idx_cell_pivot_E = idx_cell_E;
                                
                                if (face_loc_2 == BDRY_LOC::ZLO)
                                {
                                    idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                        (j + num_subghosts_conservative_var[0][1])*
                                            subghostcell_dims_conservative_var[0][0] +
                                        (interior_box_lo_idx[2] + (fill_box_hi_idx[2] - k) +
                                            num_subghosts_conservative_var[0][2])*
                                                subghostcell_dims_conservative_var[0][0]*
                                                    subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                        (j + num_subghosts_conservative_var[1][1])*
                                            subghostcell_dims_conservative_var[1][0] +
                                        (interior_box_lo_idx[2] + (fill_box_hi_idx[2] - k) +
                                            num_subghosts_conservative_var[1][2])*
                                                subghostcell_dims_conservative_var[1][0]*
                                                    subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                        (j + num_subghosts_conservative_var[2][1])*
                                            subghostcell_dims_conservative_var[2][0] +
                                        (interior_box_lo_idx[2] + (fill_box_hi_idx[2] - k) +
                                            num_subghosts_conservative_var[2][2])*
                                                subghostcell_dims_conservative_var[2][0]*
                                                    subghostcell_dims_conservative_var[2][1];
                                }
                                else if (face_loc_2 == BDRY_LOC::ZHI)
                                {
                                    idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                        (j + num_subghosts_conservative_var[0][1])*
                                            subghostcell_dims_conservative_var[0][0] +
                                        (interior_box_hi_idx[2] - (k - fill_box_lo_idx[2]) +
                                            num_subghosts_conservative_var[0][2])*
                                                subghostcell_dims_conservative_var[0][0]*
                                                    subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                        (j + num_subghosts_conservative_var[1][1])*
                                            subghostcell_dims_conservative_var[1][0] +
                                        (interior_box_hi_idx[2] - (k - fill_box_lo_idx[2]) +
                                            num_subghosts_conservative_var[1][2])*
                                                subghostcell_dims_conservative_var[1][0]*
                                                    subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                        (j + num_subghosts_conservative_var[2][1])*
                                            subghostcell_dims_conservative_var[2][0] +
                                        (interior_box_hi_idx[2] - (k - fill_box_lo_idx[2]) +
                                            num_subghosts_conservative_var[2][2])*
                                                subghostcell_dims_conservative_var[2][0]*
                                                    subghostcell_dims_conservative_var[2][1];
                                }
                                
                                /*
                                 * Set the values for density, momentum and total internal energy.
                                 */
                                
                                Real epsilon_pivot = (Q[4][idx_cell_pivot_E] -
                                    Real(1)/Real(2)*(Q[1][idx_cell_pivot_mom]*Q[1][idx_cell_pivot_mom] +
                                    Q[2][idx_cell_pivot_mom]*Q[2][idx_cell_pivot_mom] +
                                    Q[3][idx_cell_pivot_mom]*Q[3][idx_cell_pivot_mom])/
                                    Q[0][idx_cell_pivot_rho])/Q[0][idx_cell_pivot_rho];
                                
                                Real p_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &Q[0][idx_cell_pivot_rho],
                                        &epsilon_pivot,
                                        thermo_properties_ptr);
                                
                                Real p = p_pivot;
                                
                                Real T_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getTemperature(
                                        &Q[0][idx_cell_pivot_rho],
                                        &p_pivot,
                                        thermo_properties_ptr);
                                
                                Real T = -T_pivot + Real(2)*d_bdry_face_isothermal_no_slip_T[face_loc_2];
                                
                                Real rho = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getDensity(
                                        &p,
                                        &T,
                                        thermo_properties_ptr);
                                
                                Real u = -Q[1][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                                    Real(2)*d_bdry_edge_isothermal_no_slip_vel[face_loc_2*3];
                                Real v = -Q[2][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                                    Real(2)*d_bdry_edge_isothermal_no_slip_vel[face_loc_2*3 + 1];
                                Real w = -Q[3][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                                    Real(2)*d_bdry_edge_isothermal_no_slip_vel[face_loc_2*3 + 2];
                                
                                Q[0][idx_cell_rho] = rho;
                                Q[1][idx_cell_mom] = rho*u;
                                Q[2][idx_cell_mom] = rho*v;
                                Q[3][idx_cell_mom] = rho*w;
                                
                                Real epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getInternalEnergyFromTemperature(
                                        &rho,
                                        &T,
                                        thermo_properties_ptr);
                                
                                Real E = rho*epsilon +
                                    Real(1)/Real(2)*(Q[1][idx_cell_mom]*Q[1][idx_cell_mom] + Q[2][idx_cell_mom]*Q[2][idx_cell_mom] +
                                    Q[3][idx_cell_mom]*Q[3][idx_cell_mom])/rho;
                                
                                Q[4][idx_cell_E] = E;
                            }
                        }
                    }
                    
                    // Remove node locations that have boundary conditions identified.
                    bdry_node_locs.erase(std::remove(bdry_node_locs.begin(), bdry_node_locs.end(), node_loc),
                        bdry_node_locs.end());
                }
            }
        }
    }
    
    for (int ni = 0; ni < static_cast<int>(node_bdry.size()); ni++)
    {
        TBOX_ASSERT(node_bdry[ni].getBoundaryType() == BDRY::NODE3D);
        
        int node_loc(node_bdry[ni].getLocationIndex());
        
        if (std::find(bdry_node_locs.begin(), bdry_node_locs.end(), node_loc) !=
            bdry_node_locs.end())
        {
            TBOX_ERROR("FlowModelBoundaryUtilitiesSingleSpecies::fill3dNodeBoundaryData\n"
                << "Invalid node boundary condition!\n"
                << "node_loc = '" << node_loc << "'." << std::endl
                << "bdry_node_conds[node_loc] = '" << bdry_node_conds[node_loc] << "'."
                << std::endl);
        }
    }
}


void
FlowModelBoundaryUtilitiesSingleSpecies::read1dBdryNodes(
    const HAMERS_SHARED_PTR<tbox::Database>& input_db,
    std::vector<int>& node_locs,
    std::vector<int>& node_conds,
    const hier::IntVector& periodic)
{
    TBOX_DIM_ASSERT(periodic.getDim() == tbox::Dimension(1));
    
    TBOX_ASSERT(input_db);
    TBOX_ASSERT(static_cast<int>(node_locs.size()) <= NUM_1D_NODES);
    TBOX_ASSERT(*min_element(node_locs.begin(), node_locs.end()) >= 0);
    TBOX_ASSERT(*max_element(node_locs.begin(), node_locs.end()) < NUM_1D_NODES);
    TBOX_ASSERT(static_cast<int>(node_conds.size()) == NUM_1D_NODES);
    
    int num_per_dirs = 0;
    for (int id = 0; id < 1; id++)
    {
        if (periodic(id))
            ++num_per_dirs;
    }
    
    if (num_per_dirs < 1)
    {
        // Node boundary input required.
        for (int ni = 0; ni < static_cast<int>(node_locs.size()); ni++)
        {
            int s = node_locs[ni];
            
            std::string bdry_loc_str;
            switch (s)
            {
                case BDRY_LOC::XLO:
                {
                    bdry_loc_str = "boundary_node_xlo";
                    break;
                }
                case BDRY_LOC::XHI:
                {
                    bdry_loc_str = "boundary_node_xhi";
                    break;
                }
                default: NULL_STATEMENT;
            }
            
            bool need_data_read = true;
            if (num_per_dirs > 0)
            {
                if (periodic(0) && (s == BDRY_LOC::XLO || s == BDRY_LOC::XHI))
                {
                    need_data_read = false;
                }
            }
            
            if (need_data_read)
            {
                HAMERS_SHARED_PTR<tbox::Database> bdry_loc_db(
                    input_db->getDatabase(bdry_loc_str));
                std::string bdry_cond_str =
                    bdry_loc_db->getString("boundary_condition");
                if (bdry_cond_str == "ADIABATIC_NO_SLIP")
                {
                    node_conds[s] = BDRY_COND::FLOW_MODEL::ADIABATIC_NO_SLIP;
                    
                    readAdiabaticNoSlip(
                        bdry_loc_db,
                        bdry_loc_str,
                        s);
                    
                    node_locs[ni] = BOGUS_BDRY_LOC;
                }
                else if (bdry_cond_str == "ISOTHERMAL_NO_SLIP")
                {
                    node_conds[s] = BDRY_COND::FLOW_MODEL::ISOTHERMAL_NO_SLIP;
                    
                    readIsothermalNoSlip(
                        bdry_loc_db,
                        bdry_loc_str,
                        s);
                    
                    node_locs[ni] = BOGUS_BDRY_LOC;
                }
                else
                {
                    TBOX_ERROR("FlowModelBoundaryUtilitiesSingleSpecies::read1dBdryNodes()\n"
                        << "Unknown node boundary string = '"
                        << bdry_cond_str
                        << "' found in input."
                        << std::endl);
                }
            } // if (need_data_read)
        } // for (int ni = 0 ...
    } // if (num_per_dirs < 1)
    
    // Remove node locations that have boundary conditions identified.
    node_locs.erase(std::remove(node_locs.begin(), node_locs.end(), BOGUS_BDRY_LOC), node_locs.end());
}


void
FlowModelBoundaryUtilitiesSingleSpecies::read2dBdryEdges(
    const HAMERS_SHARED_PTR<tbox::Database>& input_db,
    std::vector<int>& edge_locs,
    std::vector<int>& edge_conds,
    const hier::IntVector& periodic)
{
    TBOX_DIM_ASSERT(periodic.getDim() == tbox::Dimension(2));
    
    TBOX_ASSERT(input_db);
    TBOX_ASSERT(static_cast<int>(edge_locs.size()) <= NUM_2D_EDGES);
    TBOX_ASSERT(*min_element(edge_locs.begin(), edge_locs.end()) >= 0);
    TBOX_ASSERT(*max_element(edge_locs.begin(), edge_locs.end()) < NUM_2D_EDGES);
    TBOX_ASSERT(static_cast<int>(edge_conds.size()) == NUM_2D_EDGES);
    
    int num_per_dirs = 0;
    for (int id = 0; id < 2; id++)
    {
        if (periodic(id))
            ++num_per_dirs;
    }
    
    if (num_per_dirs < 2)
    {
        // Edge boundary input required.
        for (int ei = 0; ei < static_cast<int>(edge_locs.size()); ei++)
        {
            int s = edge_locs[ei];
            
            std::string bdry_loc_str;
            switch (s)
            {
                case BDRY_LOC::XLO:
                {
                    bdry_loc_str = "boundary_edge_xlo";
                    break;
                }
                case BDRY_LOC::XHI:
                {
                    bdry_loc_str = "boundary_edge_xhi";
                    break;
                }
                case BDRY_LOC::YLO:
                {
                    bdry_loc_str = "boundary_edge_ylo";
                    break;
                }
                case BDRY_LOC::YHI:
                {
                    bdry_loc_str = "boundary_edge_yhi";
                    break;
                }
                default: NULL_STATEMENT;
            }
            
            bool need_data_read = true;
            if (num_per_dirs > 0)
            {
                if (periodic(0) && (s == BDRY_LOC::XLO || s == BDRY_LOC::XHI))
                {
                    need_data_read = false;
                }
                else if (periodic(1) && (s == BDRY_LOC::YLO || s == BDRY_LOC::YHI))
                {
                    need_data_read = false;
                }
            }
            
            if (need_data_read)
            {
                HAMERS_SHARED_PTR<tbox::Database> bdry_loc_db(
                    input_db->getDatabase(bdry_loc_str));
                std::string bdry_cond_str =
                    bdry_loc_db->getString("boundary_condition");
                if (bdry_cond_str == "ADIABATIC_NO_SLIP")
                {
                    edge_conds[s] = BDRY_COND::FLOW_MODEL::ADIABATIC_NO_SLIP;
                    
                    readAdiabaticNoSlip(
                        bdry_loc_db,
                        bdry_loc_str,
                        s);
                    
                    edge_locs[ei] = BOGUS_BDRY_LOC;
                }
                else if (bdry_cond_str == "ISOTHERMAL_NO_SLIP")
                {
                    edge_conds[s] = BDRY_COND::FLOW_MODEL::ISOTHERMAL_NO_SLIP;
                    
                    readIsothermalNoSlip(
                        bdry_loc_db,
                        bdry_loc_str,
                        s);
                    
                    edge_locs[ei] = BOGUS_BDRY_LOC;
                }
                else if (bdry_cond_str == "NONREFLECTING_OUTFLOW")
                {
                    edge_conds[s] = BDRY_COND::FLOW_MODEL::NONREFLECTING_OUTFLOW;
                    
                    readNonreflectingOutflow(
                        bdry_loc_db,
                        bdry_loc_str,
                        s);
                    
                    edge_locs[ei] = BOGUS_BDRY_LOC;
                    
                    if (d_bdry_edge_nonreflecting_outflow_beta[s] != Real(1))
                    {
                        d_use_transverse_derivatives_bc |= true;
                        d_num_ghosts_transverse_derivatives_bc = hier::IntVector::max(
                            d_num_ghosts_transverse_derivatives_bc,
                            hier::IntVector::getOne(d_dim));
                    }
                }
                else
                {
                    TBOX_ERROR("FlowModelBoundaryUtilitiesSingleSpecies::read2dBdryEdges()\n"
                        << "Unknown edge boundary string = '"
                        << bdry_cond_str
                        << "' found in input."
                        << std::endl);
                }
            } // if (need_data_read)
       } // for (int ei = 0 ...
    } // if (num_per_dirs < 2)
    
    // Remove edge locations that have boundary conditions identified.
    edge_locs.erase(std::remove(edge_locs.begin(), edge_locs.end(), BOGUS_BDRY_LOC), edge_locs.end());
}


void
FlowModelBoundaryUtilitiesSingleSpecies::read2dBdryNodes(
    const HAMERS_SHARED_PTR<tbox::Database>& input_db,
    std::vector<int>& node_locs,
    const std::vector<int>& edge_conds,
    std::vector<int>& node_conds,
    const hier::IntVector& periodic)
{
    TBOX_DIM_ASSERT(periodic.getDim() == tbox::Dimension(2));
    
    TBOX_ASSERT(input_db);
    TBOX_ASSERT(static_cast<int>(node_locs.size()) <= NUM_2D_NODES);
    TBOX_ASSERT(*min_element(node_locs.begin(), node_locs.end()) >= 0);
    TBOX_ASSERT(*max_element(node_locs.begin(), node_locs.end()) < NUM_2D_NODES);
    TBOX_ASSERT(static_cast<int>(edge_conds.size()) == NUM_2D_EDGES);
    TBOX_ASSERT(static_cast<int>(node_conds.size()) == NUM_2D_NODES);
    
    int num_per_dirs = 0;
    for (int id = 0; id < 2; id++)
    {
        if (periodic(id))
            ++num_per_dirs;
    }
    
    if (num_per_dirs < 1)
    {
        // Node boundary data required.
        for (int ni = 0; ni < static_cast<int>(node_locs.size()); ni++)
        {
            int s = node_locs[ni];
            
            std::string bdry_loc_str;
            switch (s)
            {
                case NODE_BDRY_LOC_2D::XLO_YLO:
                {
                    bdry_loc_str = "boundary_node_xlo_ylo";
                    break;
                }
                case NODE_BDRY_LOC_2D::XHI_YLO:
                {
                    bdry_loc_str = "boundary_node_xhi_ylo";
                    break;
                }
                case NODE_BDRY_LOC_2D::XLO_YHI:
                {
                    bdry_loc_str = "boundary_node_xlo_yhi";
                    break;
                }
                case NODE_BDRY_LOC_2D::XHI_YHI:
                {
                    bdry_loc_str = "boundary_node_xhi_yhi";
                    break;
                }
                default: NULL_STATEMENT;
            }
            
            HAMERS_SHARED_PTR<tbox::Database> bdry_loc_db(
                input_db->getDatabase(bdry_loc_str));
            std::string bdry_cond_str =
                bdry_loc_db->getString("boundary_condition");
            if (bdry_cond_str == "XADIABATIC_NO_SLIP")
            {
                node_conds[s] = BDRY_COND::FLOW_MODEL::XADIABATIC_NO_SLIP;
                node_locs[ni] = BOGUS_BDRY_LOC;
            }
            else if (bdry_cond_str == "YADIABATIC_NO_SLIP")
            {
                node_conds[s] = BDRY_COND::FLOW_MODEL::YADIABATIC_NO_SLIP;
                node_locs[ni] = BOGUS_BDRY_LOC;
            }
            else if (bdry_cond_str == "XISOTHERMAL_NO_SLIP")
            {
                node_conds[s] = BDRY_COND::FLOW_MODEL::XISOTHERMAL_NO_SLIP;
                node_locs[ni] = BOGUS_BDRY_LOC;
            }
            else if (bdry_cond_str == "YISOTHERMAL_NO_SLIP")
            {
                node_conds[s] = BDRY_COND::FLOW_MODEL::YISOTHERMAL_NO_SLIP;
                node_locs[ni] = BOGUS_BDRY_LOC;
            }
            else
            {
                TBOX_ERROR("FlowModelBoundaryUtilitiesSingleSpecies::read2dBdryNodes()\n"
                    << "Unknown node boundary string = '"
                    << bdry_cond_str
                    << "' found in input."
                    << std::endl);
            }
            
            std::string proper_edge;
            std::string proper_edge_data;
            bool no_edge_data_found = false;
            if ((bdry_cond_str == "XADIABATIC_NO_SLIP") ||
                (bdry_cond_str == "XISOTHERMAL_NO_SLIP"))
            {
                if (s == NODE_BDRY_LOC_2D::XLO_YLO ||
                    s == NODE_BDRY_LOC_2D::XLO_YHI)
                {
                    proper_edge = "XLO";
                    if (bdry_cond_str == "XADIABATIC_NO_SLIP" &&
                        edge_conds[BDRY_LOC::XLO] != BDRY_COND::FLOW_MODEL::ADIABATIC_NO_SLIP)
                    {
                        no_edge_data_found = true;
                        proper_edge_data = "ADIABATIC_NO_SLIP";
                    }
                    if (bdry_cond_str == "XISOTHERMAL_NO_SLIP" &&
                        edge_conds[BDRY_LOC::XLO] != BDRY_COND::FLOW_MODEL::ISOTHERMAL_NO_SLIP)
                    {
                        no_edge_data_found = true;
                        proper_edge_data = "ISOTHERMAL_NO_SLIP";
                    }
                }
                else
                {
                    proper_edge = "XHI";
                    if (bdry_cond_str == "XADIABATIC_NO_SLIP" &&
                        edge_conds[BDRY_LOC::XHI] != BDRY_COND::FLOW_MODEL::ADIABATIC_NO_SLIP)
                    {
                        no_edge_data_found = true;
                        proper_edge_data = "ADIABATIC_NO_SLIP";
                    }
                    if (bdry_cond_str == "XISOTHERMAL_NO_SLIP" &&
                        edge_conds[BDRY_LOC::XHI] != BDRY_COND::FLOW_MODEL::ISOTHERMAL_NO_SLIP)
                    {
                        no_edge_data_found = true;
                        proper_edge_data = "ISOTHERMAL_NO_SLIP";
                    }
                }
            }
            else if ((bdry_cond_str == "YADIABATIC_NO_SLIP") ||
                     (bdry_cond_str == "YISOTHERMAL_NO_SLIP"))
            {
                if (s == NODE_BDRY_LOC_2D::XLO_YLO ||
                    s == NODE_BDRY_LOC_2D::XHI_YLO)
                {
                    proper_edge = "YLO";
                    if (bdry_cond_str == "YADIABATIC_NO_SLIP" &&
                        edge_conds[BDRY_LOC::YLO] != BDRY_COND::FLOW_MODEL::ADIABATIC_NO_SLIP)
                    {
                        no_edge_data_found = true;
                        proper_edge_data = "ADIABATIC_NO_SLIP";
                    }
                    if (bdry_cond_str == "YISOTHREMAL_NO_SLIP" &&
                        edge_conds[BDRY_LOC::YLO] != BDRY_COND::FLOW_MODEL::ISOTHERMAL_NO_SLIP)
                    {
                        no_edge_data_found = true;
                        proper_edge_data = "ISOTHERMAL_NO_SLIP";
                    }
                }
                else
                {
                    proper_edge = "YHI";
                    if (bdry_cond_str == "YADIABATIC_NO_SLIP" &&
                        edge_conds[BDRY_LOC::YHI] != BDRY_COND::FLOW_MODEL::ADIABATIC_NO_SLIP)
                    {
                        no_edge_data_found = true;
                        proper_edge_data = "ADIABATIC_NO_SLIP";
                    }
                    if (bdry_cond_str == "YISOTHERMAL_NO_SLIP" &&
                        edge_conds[BDRY_LOC::YHI] != BDRY_COND::FLOW_MODEL::ISOTHERMAL_NO_SLIP)
                    {
                        no_edge_data_found = true;
                        proper_edge_data = "ISOTHERMAL_NO_SLIP";
                    }
                }
            }
            
            if (no_edge_data_found)
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelBoundaryUtilitiesSingleSpecies::read2dBdryNodes()\n"
                    << "Bdry condition '"
                    << bdry_cond_str
                    << "' found for "
                    << bdry_loc_str
                    << "\n but no "
                    << proper_edge_data
                    << " data found for edge "
                    << proper_edge << std::endl);
            }
        } // for (int ni = 0 ...
    } // if (num_per_dirs < 1)
    
    // Remove node locations that have boundary conditions identified.
    node_locs.erase(std::remove(node_locs.begin(), node_locs.end(), BOGUS_BDRY_LOC), node_locs.end());
}


void
FlowModelBoundaryUtilitiesSingleSpecies::read3dBdryFaces(
    const HAMERS_SHARED_PTR<tbox::Database>& input_db,
    std::vector<int>& face_locs,
    std::vector<int>& face_conds,
    const hier::IntVector& periodic)
{
    TBOX_DIM_ASSERT(periodic.getDim() == tbox::Dimension(3));
    
    TBOX_ASSERT(input_db);
    TBOX_ASSERT(static_cast<int>(face_locs.size()) <= NUM_3D_FACES);
    TBOX_ASSERT(*min_element(face_locs.begin(), face_locs.end()) >= 0);
    TBOX_ASSERT(*max_element(face_locs.begin(), face_locs.end()) < NUM_3D_FACES);
    TBOX_ASSERT(static_cast<int>(face_conds.size()) == NUM_3D_FACES);
    
    int num_per_dirs = 0;
    for (int id = 0; id < 3; id++)
    {
        if (periodic(id))
            ++num_per_dirs;
    }
    
    if (num_per_dirs < 3)
    {
        // Face boundary input required.
        for (int fi = 0; fi < static_cast<int>(face_locs.size()); fi++)
        {
            int s = face_locs[fi];
            
            std::string bdry_loc_str;
            switch (s)
            {
                case BDRY_LOC::XLO:
                {
                    bdry_loc_str = "boundary_face_xlo";
                    break;
                }
                case BDRY_LOC::XHI:
                {
                    bdry_loc_str = "boundary_face_xhi";
                    break;
                }
                case BDRY_LOC::YLO:
                {
                    bdry_loc_str = "boundary_face_ylo";
                    break;
                }
                case BDRY_LOC::YHI:
                {
                    bdry_loc_str = "boundary_face_yhi";
                    break;
                }
                case BDRY_LOC::ZLO:
                {
                    bdry_loc_str = "boundary_face_zlo";
                    break;
                }
                case BDRY_LOC::ZHI:
                {
                    bdry_loc_str = "boundary_face_zhi";
                    break;
                }
                default: NULL_STATEMENT;
            }
            
            bool need_data_read = true;
            if (num_per_dirs > 0)
            {
                if (periodic(0) && (s == BDRY_LOC::XLO || s == BDRY_LOC::XHI))
                {
                    need_data_read = false;
                }
                else if (periodic(1) && (s == BDRY_LOC::YLO || s == BDRY_LOC::YHI))
                {
                    need_data_read = false;
                }
                else if (periodic(2) && (s == BDRY_LOC::ZLO || s == BDRY_LOC::ZHI))
                {
                    need_data_read = false;
                }
            }
            
            if (need_data_read)
            {
                HAMERS_SHARED_PTR<tbox::Database> bdry_loc_db(
                    input_db->getDatabase(bdry_loc_str));
                std::string bdry_cond_str =
                    bdry_loc_db->getString("boundary_condition");
                if (bdry_cond_str == "ADIABATIC_NO_SLIP")
                {
                    face_conds[s] = BDRY_COND::FLOW_MODEL::ADIABATIC_NO_SLIP;
                    
                    readAdiabaticNoSlip(
                        bdry_loc_db,
                        bdry_loc_str,
                        s);
                    
                    face_locs[fi] = BOGUS_BDRY_LOC;
                }
                else if (bdry_cond_str == "ISOTHERMAL_NO_SLIP")
                {
                    face_conds[s] = BDRY_COND::FLOW_MODEL::ISOTHERMAL_NO_SLIP;
                    
                    readIsothermalNoSlip(
                        bdry_loc_db,
                        bdry_loc_str,
                        s);
                    
                    face_locs[fi] = BOGUS_BDRY_LOC;
                }
                else if (bdry_cond_str == "NONREFLECTING_OUTFLOW")
                {
                    face_conds[s] = BDRY_COND::FLOW_MODEL::NONREFLECTING_OUTFLOW;
                    
                    readNonreflectingOutflow(
                        bdry_loc_db,
                        bdry_loc_str,
                        s);
                    
                    face_locs[fi] = BOGUS_BDRY_LOC;
                    
                    if (d_bdry_face_nonreflecting_outflow_beta[s] != Real(1))
                    {
                        d_use_transverse_derivatives_bc |= true;
                        d_num_ghosts_transverse_derivatives_bc = hier::IntVector::max(
                            d_num_ghosts_transverse_derivatives_bc,
                            hier::IntVector::getOne(d_dim));
                    }
                }
                else
                {
                    TBOX_ERROR("FlowModelBoundaryUtilitiesSingleSpecies::read3dBdryFaces\n"
                        << "Unknown face boundary string = '"
                        << bdry_cond_str
                        << "' found in input."
                        << std::endl);
                }
            } // if (need_data_read)
        } // for (int fi = 0 ...
    } // if (num_per_dirs < 3)
    
    // Remove face locations that have boundary conditions identified.
    face_locs.erase(std::remove(face_locs.begin(), face_locs.end(), BOGUS_BDRY_LOC), face_locs.end());
}


void
FlowModelBoundaryUtilitiesSingleSpecies::read3dBdryEdges(
    const HAMERS_SHARED_PTR<tbox::Database>& input_db,
    std::vector<int>& edge_locs,
    const std::vector<int>& face_conds,
    std::vector<int>& edge_conds,
    const hier::IntVector& periodic)
{
    TBOX_DIM_ASSERT(periodic.getDim() == tbox::Dimension(3));
    
    TBOX_ASSERT(input_db);
    TBOX_ASSERT(static_cast<int>(edge_locs.size()) <= NUM_3D_EDGES);
    TBOX_ASSERT(*min_element(edge_locs.begin(), edge_locs.end()) >= 0);
    TBOX_ASSERT(*max_element(edge_locs.begin(), edge_locs.end()) < NUM_3D_EDGES);
    TBOX_ASSERT(static_cast<int>(face_conds.size()) == NUM_3D_FACES);
    TBOX_ASSERT(static_cast<int>(edge_conds.size()) == NUM_3D_EDGES);
    
    int num_per_dirs = 0;
    for (int id = 0; id < 3; id++)
    {
        if (periodic(id))
            ++num_per_dirs;
    }
    
    if (num_per_dirs < 2)
    {
        // edge boundary input required
        for (int ei = 0; ei < static_cast<int>(edge_locs.size()); ei++)
        {
            int s = edge_locs[ei];
            
            std::string bdry_loc_str;
            switch (s)
            {
                case EDGE_BDRY_LOC_3D::YLO_ZLO:
                {
                    bdry_loc_str = "boundary_edge_ylo_zlo";
                    break;
                }
                case EDGE_BDRY_LOC_3D::YHI_ZLO:
                {
                    bdry_loc_str = "boundary_edge_yhi_zlo";
                    break;
                }
                case EDGE_BDRY_LOC_3D::YLO_ZHI:
                {
                    bdry_loc_str = "boundary_edge_ylo_zhi";
                    break;
                }
                case EDGE_BDRY_LOC_3D::YHI_ZHI:
                {
                    bdry_loc_str = "boundary_edge_yhi_zhi";
                    break;
                }
                case EDGE_BDRY_LOC_3D::XLO_ZLO:
                {
                    bdry_loc_str = "boundary_edge_xlo_zlo";
                    break;
                }
                case EDGE_BDRY_LOC_3D::XLO_ZHI:
                {
                    bdry_loc_str = "boundary_edge_xlo_zhi";
                    break;
                }
                case EDGE_BDRY_LOC_3D::XHI_ZLO:
                {
                    bdry_loc_str = "boundary_edge_xhi_zlo";
                    break;
                }
                case EDGE_BDRY_LOC_3D::XHI_ZHI:
                {
                    bdry_loc_str = "boundary_edge_xhi_zhi";
                    break;
                }
                case EDGE_BDRY_LOC_3D::XLO_YLO:
                {
                    bdry_loc_str = "boundary_edge_xlo_ylo";
                    break;
                }
                case EDGE_BDRY_LOC_3D::XHI_YLO:
                {
                    bdry_loc_str = "boundary_edge_xhi_ylo";
                    break;
                }
                case EDGE_BDRY_LOC_3D::XLO_YHI:
                {
                    bdry_loc_str = "boundary_edge_xlo_yhi";
                    break;
                }
                case EDGE_BDRY_LOC_3D::XHI_YHI:
                {
                    bdry_loc_str = "boundary_edge_xhi_yhi";
                    break;
                }
                default: NULL_STATEMENT;
            }
            
            bool need_data_read = false;
            if (num_per_dirs == 0)
            {
                need_data_read = true;
            }
            else if (periodic(0) &&
                     (s == EDGE_BDRY_LOC_3D::YLO_ZLO ||
                      s == EDGE_BDRY_LOC_3D::YHI_ZLO ||
                      s == EDGE_BDRY_LOC_3D::YLO_ZHI ||
                      s == EDGE_BDRY_LOC_3D::YHI_ZHI))
            {
                need_data_read = true;
            }
            else if (periodic(1) &&
                     (s == EDGE_BDRY_LOC_3D::XLO_ZLO ||
                      s == EDGE_BDRY_LOC_3D::XLO_ZHI ||
                      s == EDGE_BDRY_LOC_3D::XHI_ZLO ||
                      s == EDGE_BDRY_LOC_3D::XHI_ZHI))
            {
                need_data_read = true;
            }
            else if (periodic(2) &&
                     (s == EDGE_BDRY_LOC_3D::XLO_YLO ||
                      s == EDGE_BDRY_LOC_3D::XHI_YLO ||
                      s == EDGE_BDRY_LOC_3D::XLO_YHI ||
                      s == EDGE_BDRY_LOC_3D::XHI_YHI))
            {
                need_data_read = true;
            }
            
            if (need_data_read)
            {
                HAMERS_SHARED_PTR<tbox::Database> bdry_loc_db(
                   input_db->getDatabase(bdry_loc_str));
                
                std::string bdry_cond_str =
                   bdry_loc_db->getString("boundary_condition");
                if (bdry_cond_str == "XADIABATIC_NO_SLIP")
                {
                    edge_conds[s] = BDRY_COND::FLOW_MODEL::XADIABATIC_NO_SLIP;
                    edge_locs[ei] = BOGUS_BDRY_LOC;
                }
                else if (bdry_cond_str == "YADIABATIC_NO_SLIP")
                {
                    edge_conds[s] = BDRY_COND::FLOW_MODEL::YADIABATIC_NO_SLIP;
                    edge_locs[ei] = BOGUS_BDRY_LOC;
                }
                else if (bdry_cond_str == "ZADIABATIC_NO_SLIP")
                {
                    edge_conds[s] = BDRY_COND::FLOW_MODEL::ZADIABATIC_NO_SLIP;
                    edge_locs[ei] = BOGUS_BDRY_LOC;
                }
                else if (bdry_cond_str == "XISOTHERMAL_NO_SLIP")
                {
                    edge_conds[s] = BDRY_COND::FLOW_MODEL::XISOTHERMAL_NO_SLIP;
                    edge_locs[ei] = BOGUS_BDRY_LOC;
                }
                else if (bdry_cond_str == "YISOTHERMAL_NO_SLIP")
                {
                    edge_conds[s] = BDRY_COND::FLOW_MODEL::YISOTHERMAL_NO_SLIP;
                    edge_locs[ei] = BOGUS_BDRY_LOC;
                }
                else if (bdry_cond_str == "ZISOTHERMAL_NO_SLIP")
                {
                    edge_conds[s] = BDRY_COND::FLOW_MODEL::ZISOTHERMAL_NO_SLIP;
                    edge_locs[ei] = BOGUS_BDRY_LOC;
                }
                else
                {
                    TBOX_ERROR("FlowModelBoundaryUtilitiesSingleSpecies::read3dBdryEdges\n"
                        << "Unknown edge boundary string = '"
                        << bdry_cond_str
                        << "' found in input."
                        << std::endl);
                }
                
                bool ambiguous_type = false;
                if ((bdry_cond_str == "XADIABATIC_NO_SLIP") ||
                    (bdry_cond_str == "XISOTHERMAL_NO_SLIP"))
                {
                    if (s == EDGE_BDRY_LOC_3D::YLO_ZLO ||
                        s == EDGE_BDRY_LOC_3D::YHI_ZLO ||
                        s == EDGE_BDRY_LOC_3D::YLO_ZHI ||
                        s == EDGE_BDRY_LOC_3D::YHI_ZHI)
                    {
                        ambiguous_type = true;
                    }
                }
                else if ((bdry_cond_str == "YADIABATIC_NO_SLIP") ||
                         (bdry_cond_str == "YISOTHERMAL_NO_SLIP"))
                {
                    if (s == EDGE_BDRY_LOC_3D::XLO_ZLO ||
                        s == EDGE_BDRY_LOC_3D::XLO_ZHI ||
                        s == EDGE_BDRY_LOC_3D::XHI_ZLO ||
                        s == EDGE_BDRY_LOC_3D::XHI_ZHI)
                    {
                        ambiguous_type = true;
                    }
                }
                else if ((bdry_cond_str == "ZADIABATIC_NO_SLIP") ||
                         (bdry_cond_str == "ZISOTHERMAL_NO_SLIP"))
                {
                    if (s == EDGE_BDRY_LOC_3D::XLO_YLO ||
                        s == EDGE_BDRY_LOC_3D::XHI_YLO ||
                        s == EDGE_BDRY_LOC_3D::XLO_YHI ||
                        s == EDGE_BDRY_LOC_3D::XHI_YHI)
                    {
                        ambiguous_type = true;
                    }
                }
                if (ambiguous_type)
                {
                    TBOX_ERROR(d_object_name
                        << ": FlowModelBoundaryUtilitiesSingleSpecies::read3dBdryEdges()\n"
                        << "Ambiguous bdry condition '"
                        << bdry_cond_str
                        << "' found for "
                        << bdry_loc_str
                        << std::endl);
                }
                
                std::string proper_face;
                std::string proper_face_data;
                bool no_face_data_found = false;
                if ((bdry_cond_str == "XADIABATIC_NO_SLIP") ||
                    (bdry_cond_str == "XISOTHERMAL_NO_SLIP"))
                {
                    if (s == EDGE_BDRY_LOC_3D::XLO_ZLO ||
                        s == EDGE_BDRY_LOC_3D::XLO_ZHI ||
                        s == EDGE_BDRY_LOC_3D::XLO_YLO ||
                        s == EDGE_BDRY_LOC_3D::XLO_YHI)
                    {
                        proper_face = "XLO";
                        if (bdry_cond_str == "XADIABATIC_NO_SLIP" &&
                            face_conds[BDRY_LOC::XLO] != BDRY_COND::FLOW_MODEL::ADIABATIC_NO_SLIP)
                        {
                            no_face_data_found = true;
                            proper_face_data = "ADIABATIC_NO_SLIP";
                        }
                        if (bdry_cond_str == "XISOTHERMAL_NO_SLIP" &&
                            face_conds[BDRY_LOC::XLO] != BDRY_COND::FLOW_MODEL::ISOTHERMAL_NO_SLIP)
                        {
                            no_face_data_found = true;
                            proper_face_data = "ISOTHERMAL_NO_SLIP";
                        }
                    }
                    else
                    {
                        proper_face = "XHI";
                        if (bdry_cond_str == "XADIABATIC_NO_SLIP" &&
                            face_conds[BDRY_LOC::XHI] != BDRY_COND::FLOW_MODEL::ADIABATIC_NO_SLIP)
                        {
                            no_face_data_found = true;
                            proper_face_data = "ADIABATIC_NO_SLIP";
                        }
                        if (bdry_cond_str == "XISOTHERMAL_NO_SLIP" &&
                            face_conds[BDRY_LOC::XHI] != BDRY_COND::FLOW_MODEL::ISOTHERMAL_NO_SLIP)
                        {
                            no_face_data_found = true;
                            proper_face_data = "ISOTHERMAL_NO_SLIP";
                        }
                    }
                }
                else if ((bdry_cond_str == "YADIABATIC_NO_SLIP") ||
                         (bdry_cond_str == "YISOTHERMAL_NO_SLIP"))
                {
                    if (s == EDGE_BDRY_LOC_3D::YLO_ZLO ||
                        s == EDGE_BDRY_LOC_3D::YLO_ZHI ||
                        s == EDGE_BDRY_LOC_3D::XLO_YLO ||
                        s == EDGE_BDRY_LOC_3D::XHI_YLO)
                    {
                        proper_face = "YLO";
                        if (bdry_cond_str == "YADIABATIC_NO_SLIP" &&
                            face_conds[BDRY_LOC::YLO] != BDRY_COND::FLOW_MODEL::ADIABATIC_NO_SLIP)
                        {
                            no_face_data_found = true;
                            proper_face_data = "ADIABATIC_NO_SLIP";
                        }
                        if (bdry_cond_str == "YISOTHERMAL_NO_SLIP" &&
                            face_conds[BDRY_LOC::YLO] != BDRY_COND::FLOW_MODEL::ISOTHERMAL_NO_SLIP)
                        {
                            no_face_data_found = true;
                            proper_face_data = "ISOTHERMAL_NO_SLIP";
                        }
                    }
                    else
                    {
                        proper_face = "YHI";
                        if (bdry_cond_str == "YADIABATIC_NO_SLIP" &&
                            face_conds[BDRY_LOC::YHI] != BDRY_COND::FLOW_MODEL::ADIABATIC_NO_SLIP)
                        {
                            no_face_data_found = true;
                            proper_face_data = "ADIABATIC_NO_SLIP";
                        }
                        if (bdry_cond_str == "YISOTHERMAL_NO_SLIP" &&
                            face_conds[BDRY_LOC::YHI] != BDRY_COND::FLOW_MODEL::ISOTHERMAL_NO_SLIP)
                        {
                            no_face_data_found = true;
                            proper_face_data = "ISOTHERMAL_NO_SLIP";
                        }
                    }
                }
                else if ((bdry_cond_str == "ZADIABATIC_NO_SLIP") ||
                         (bdry_cond_str == "ZISOTHERMAL_NO_SLIP"))
                {
                    if (s == EDGE_BDRY_LOC_3D::XLO_ZLO ||
                        s == EDGE_BDRY_LOC_3D::YHI_ZLO ||
                        s == EDGE_BDRY_LOC_3D::YLO_ZLO ||
                        s == EDGE_BDRY_LOC_3D::XHI_ZLO)
                    {
                        proper_face = "ZLO";
                        if (bdry_cond_str == "ZADIABATIC_NO_SLIP" &&
                            face_conds[BDRY_LOC::ZLO] != BDRY_COND::FLOW_MODEL::ADIABATIC_NO_SLIP)
                        {
                            no_face_data_found = true;
                            proper_face_data = "ADIABATIC_NO_SLIP";
                        }
                        if (bdry_cond_str == "ZISOTHERMAL_NO_SLIP" &&
                            face_conds[BDRY_LOC::ZLO] != BDRY_COND::FLOW_MODEL::ISOTHERMAL_NO_SLIP)
                        {
                            no_face_data_found = true;
                            proper_face_data = "ISOTHERMAL_NO_SLIP";
                        }
                    }
                    else
                    {
                        proper_face = "ZHI";
                        if (bdry_cond_str == "ZADIABATIC_NO_SLIP" &&
                            face_conds[BDRY_LOC::ZHI] != BDRY_COND::FLOW_MODEL::ADIABATIC_NO_SLIP)
                        {
                            no_face_data_found = true;
                            proper_face_data = "ADIABATIC_NO_SLIP";
                        }
                        if (bdry_cond_str == "ZISOTHERMAL_NO_SLIP" &&
                            face_conds[BDRY_LOC::ZHI] != BDRY_COND::FLOW_MODEL::ISOTHERMAL_NO_SLIP)
                        {
                            no_face_data_found = true;
                            proper_face_data = "ISOTHERMAL_NO_SLIP";
                        }
                    }
                }
                
                if (no_face_data_found)
                {
                    TBOX_ERROR(d_object_name
                        << ": FlowModelBoundaryUtilitiesSingleSpecies::read3dBdryEdges()\n"
                        << "Bdry condition '"
                        << bdry_cond_str
                        << "' found for "
                        << bdry_loc_str
                        << "\n but no "
                        << proper_face_data
                        << " data found for face "
                        << proper_face << std::endl);
                }
            } // if (need_data_read)
        } // for (int ei = 0 ...
    } // if (num_per_dirs < 2)
    
    // Remove edge locations that have boundary conditions identified.
    edge_locs.erase(std::remove(edge_locs.begin(), edge_locs.end(), BOGUS_BDRY_LOC), edge_locs.end());
}


void
FlowModelBoundaryUtilitiesSingleSpecies::read3dBdryNodes(
    const HAMERS_SHARED_PTR<tbox::Database>& input_db,
    std::vector<int>& node_locs,
    const std::vector<int>& face_conds,
    std::vector<int>& node_conds,
    const hier::IntVector& periodic)
{
    TBOX_DIM_ASSERT(periodic.getDim() == tbox::Dimension(3));
    
    TBOX_ASSERT(input_db);
    TBOX_ASSERT(static_cast<int>(node_locs.size()) <= NUM_3D_NODES);
    TBOX_ASSERT(*min_element(node_locs.begin(), node_locs.end()) >= 0);
    TBOX_ASSERT(*max_element(node_locs.begin(), node_locs.end()) < NUM_3D_NODES);
    TBOX_ASSERT(static_cast<int>(face_conds.size()) == NUM_3D_FACES);
    TBOX_ASSERT(static_cast<int>(node_conds.size()) == NUM_3D_NODES);
    
    int num_per_dirs = 0;
    for (int id = 0; id < 3; id++)
    {
        if (periodic(id))
            ++num_per_dirs;
    }
    
    if (num_per_dirs < 1)
    {
        // Node boundary data required.
        for (int ni = 0; ni < static_cast<int>(node_locs.size()); ni++)
        {
            int s = node_locs[ni];
            
            std::string bdry_loc_str;
            switch (s)
            {
                case NODE_BDRY_LOC_3D::XLO_YLO_ZLO:
                {
                    bdry_loc_str = "boundary_node_xlo_ylo_zlo";
                    break;
                }
                case NODE_BDRY_LOC_3D::XHI_YLO_ZLO:
                {
                    bdry_loc_str = "boundary_node_xhi_ylo_zlo";
                    break;
                }
                case NODE_BDRY_LOC_3D::XLO_YHI_ZLO:
                {
                    bdry_loc_str = "boundary_node_xlo_yhi_zlo";
                    break;
                }
                case NODE_BDRY_LOC_3D::XHI_YHI_ZLO:
                {
                    bdry_loc_str = "boundary_node_xhi_yhi_zlo";
                    break;
                }
                case NODE_BDRY_LOC_3D::XLO_YLO_ZHI:
                {
                    bdry_loc_str = "boundary_node_xlo_ylo_zhi";
                    break;
                }
                case NODE_BDRY_LOC_3D::XHI_YLO_ZHI:
                {
                    bdry_loc_str = "boundary_node_xhi_ylo_zhi";
                    break;
                }
                case NODE_BDRY_LOC_3D::XLO_YHI_ZHI:
                {
                    bdry_loc_str = "boundary_node_xlo_yhi_zhi";
                    break;
                }
                case NODE_BDRY_LOC_3D::XHI_YHI_ZHI:
                {
                    bdry_loc_str = "boundary_node_xhi_yhi_zhi";
                    break;
                }
                default: NULL_STATEMENT;
            }
            HAMERS_SHARED_PTR<tbox::Database> bdry_loc_db(
                input_db->getDatabase(bdry_loc_str));
            std::string bdry_cond_str =
                bdry_loc_db->getString("boundary_condition");
            if (bdry_cond_str == "XADIABATIC_NO_SLIP")
            {
                node_conds[s] = BDRY_COND::FLOW_MODEL::XADIABATIC_NO_SLIP;
                node_locs[ni] = BOGUS_BDRY_LOC;
            }
            else if (bdry_cond_str == "YADIABATIC_NO_SLIP")
            {
                node_conds[s] = BDRY_COND::FLOW_MODEL::YADIABATIC_NO_SLIP;
                node_locs[ni] = BOGUS_BDRY_LOC;
            }
            else if (bdry_cond_str == "ZADIABATIC_NO_SLIP")
            {
                node_conds[s] = BDRY_COND::FLOW_MODEL::ZADIABATIC_NO_SLIP;
                node_locs[ni] = BOGUS_BDRY_LOC;
            }
            else if (bdry_cond_str == "XISOTHERMAL_NO_SLIP")
            {
                node_conds[s] = BDRY_COND::FLOW_MODEL::XISOTHERMAL_NO_SLIP;
                node_locs[ni] = BOGUS_BDRY_LOC;
            }
            else if (bdry_cond_str == "YISOTHERMAL_NO_SLIP")
            {
                node_conds[s] = BDRY_COND::FLOW_MODEL::YISOTHERMAL_NO_SLIP;
                node_locs[ni] = BOGUS_BDRY_LOC;
            }
            else if (bdry_cond_str == "ZISOTHERMAL_NO_SLIP")
            {
                node_conds[s] = BDRY_COND::FLOW_MODEL::ZISOTHERMAL_NO_SLIP;
                node_locs[ni] = BOGUS_BDRY_LOC;
            }
            else
            {
                TBOX_ERROR("FlowModelBoundaryUtilitiesSingleSpecies::read3dBdryNodes()\n"
                    << "Unknown node boundary string = '"
                    << bdry_cond_str
                    << "' found in input."
                    << std::endl);
            }
            
            std::string proper_face;
            std::string proper_face_data;
            bool no_face_data_found = false;
            if ((bdry_cond_str == "XADIABATIC_NO_SLIP") ||
                (bdry_cond_str == "XISOTHERMAL_NO_SLIP"))
            {
                if (s == NODE_BDRY_LOC_3D::XLO_YLO_ZLO ||
                    s == NODE_BDRY_LOC_3D::XLO_YHI_ZLO ||
                    s == NODE_BDRY_LOC_3D::XLO_YLO_ZHI ||
                    s == NODE_BDRY_LOC_3D::XLO_YHI_ZHI)
                {
                    proper_face = "XLO";
                    if (bdry_cond_str == "XADIABATIC_NO_SLIP" &&
                        face_conds[BDRY_LOC::XLO] != BDRY_COND::FLOW_MODEL::ADIABATIC_NO_SLIP)
                    {
                        no_face_data_found = true;
                        proper_face_data = "ADIABATIC_NO_SLIP";
                    }
                    if (bdry_cond_str == "XISOTHERMAL_NO_SLIP" &&
                        face_conds[BDRY_LOC::XLO] != BDRY_COND::FLOW_MODEL::ISOTHERMAL_NO_SLIP)
                    {
                        no_face_data_found = true;
                        proper_face_data = "ISOTHERMAL_NO_SLIP";
                    }
                }
                else
                {
                    proper_face = "XHI";
                    if (bdry_cond_str == "XADIABATIC_NO_SLIP" &&
                        face_conds[BDRY_LOC::XHI] != BDRY_COND::FLOW_MODEL::ADIABATIC_NO_SLIP)
                    {
                        no_face_data_found = true;
                        proper_face_data = "ADIABATIC_NO_SLIP";
                    }
                    if (bdry_cond_str == "XISOTHERMAL_NO_SLIP" &&
                        face_conds[BDRY_LOC::XHI] != BDRY_COND::FLOW_MODEL::ISOTHERMAL_NO_SLIP)
                    {
                        no_face_data_found = true;
                        proper_face_data = "ISOTHERMAL_NO_SLIP";
                    }
                }
            }
            else if ((bdry_cond_str == "YADIABATIC_NO_SLIP") ||
                     (bdry_cond_str == "YISOTHERMAL_NO_SLIP"))
            {
                if (s == NODE_BDRY_LOC_3D::XLO_YLO_ZLO ||
                    s == NODE_BDRY_LOC_3D::XHI_YLO_ZLO ||
                    s == NODE_BDRY_LOC_3D::XLO_YLO_ZHI ||
                    s == NODE_BDRY_LOC_3D::XHI_YLO_ZHI)
                {
                    proper_face = "YLO";
                    if (bdry_cond_str == "YADIABATIC_NO_SLIP" &&
                        face_conds[BDRY_LOC::YLO] != BDRY_COND::FLOW_MODEL::ADIABATIC_NO_SLIP)
                    {
                        no_face_data_found = true;
                        proper_face_data = "ADIABATIC_NO_SLIP";
                    }
                    if (bdry_cond_str == "YISOTHERMAL_NO_SLIP" &&
                        face_conds[BDRY_LOC::YLO] != BDRY_COND::FLOW_MODEL::ISOTHERMAL_NO_SLIP)
                    {
                        no_face_data_found = true;
                        proper_face_data = "ISOTHERMAL_NO_SLIP";
                    }
                }
                else
                {
                    proper_face = "YHI";
                    if (bdry_cond_str == "YADIABATIC_NO_SLIP" &&
                        face_conds[BDRY_LOC::YHI] != BDRY_COND::FLOW_MODEL::ADIABATIC_NO_SLIP)
                    {
                        no_face_data_found = true;
                        proper_face_data = "ADIABATIC_NO_SLIP";
                    }
                    if (bdry_cond_str == "YISOTHERMAL_NO_SLIP" &&
                        face_conds[BDRY_LOC::YHI] != BDRY_COND::FLOW_MODEL::ISOTHERMAL_NO_SLIP)
                    {
                        no_face_data_found = true;
                        proper_face_data = "ISOTHERMAL_NO_SLIP";
                    }
                }
            }
            else if ((bdry_cond_str == "ZADIABATIC_NO_SLIP") ||
                     (bdry_cond_str == "ZISOTHERMAL_NO_SLIP"))
            {
                if (s == NODE_BDRY_LOC_3D::XLO_YLO_ZLO ||
                    s == NODE_BDRY_LOC_3D::XHI_YLO_ZLO ||
                    s == NODE_BDRY_LOC_3D::XLO_YHI_ZLO ||
                    s == NODE_BDRY_LOC_3D::XHI_YHI_ZLO)
                {
                    proper_face = "ZLO";
                    if (bdry_cond_str == "ZADIABATIC_NO_SLIP" &&
                        face_conds[BDRY_LOC::ZLO] != BDRY_COND::FLOW_MODEL::ADIABATIC_NO_SLIP)
                    {
                        no_face_data_found = true;
                        proper_face_data = "ADIABATIC_NO_SLIP";
                    }
                    if (bdry_cond_str == "ZISOTHERMAL_NO_SLIP" &&
                        face_conds[BDRY_LOC::ZLO] != BDRY_COND::FLOW_MODEL::ISOTHERMAL_NO_SLIP)
                    {
                        no_face_data_found = true;
                        proper_face_data = "ISOTHERMAL_NO_SLIP";
                    }
                }
                else
                {
                    proper_face = "ZHI";
                    if (bdry_cond_str == "ZADIABATIC_NO_SLIP" &&
                        face_conds[BDRY_LOC::ZHI] != BDRY_COND::FLOW_MODEL::ADIABATIC_NO_SLIP)
                    {
                        no_face_data_found = true;
                        proper_face_data = "ADIABATIC_NO_SLIP";
                    }
                    if (bdry_cond_str == "ZISOTHERMAL_NO_SLIP" &&
                        face_conds[BDRY_LOC::ZHI] != BDRY_COND::FLOW_MODEL::ISOTHERMAL_NO_SLIP)
                    {
                        no_face_data_found = true;
                        proper_face_data = "ISOTHERMAL_NO_SLIP";
                    }
                }
            }
            
            if (no_face_data_found)
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelBoundaryUtilitiesSingleSpecies::read3dBdryNodes()\n"
                    << "Bdry condition '"
                    << bdry_cond_str
                    << "' found for "
                    << bdry_loc_str
                    << "\n but no "
                    << proper_face_data
                    << " data found for face "
                    << proper_face << std::endl);
            }
        } // for (int ni = 0 ...
    } // if (num_per_dirs < 1)
    
    // Remove node locations that have boundary conditions identified.
    node_locs.erase(std::remove(node_locs.begin(), node_locs.end(), BOGUS_BDRY_LOC), node_locs.end());
}


void
FlowModelBoundaryUtilitiesSingleSpecies::readAdiabaticNoSlip(
    const HAMERS_SHARED_PTR<tbox::Database>& db,
    std::string& db_name,
    int bdry_location_index)
{
    TBOX_ASSERT(db);
    TBOX_ASSERT(!db_name.empty());
    
    std::vector<Real> data_vel;
    
    if (db->keyExists("velocity"))
    {
        data_vel = db->getRealVector("velocity");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBoundaryUtilitiesSingleSpecies::readAdiabaticNoSlip()\n"
            << "'velocity' entry missing from '"
            << db_name
            << "' input database."
            << std::endl);
    }
    
    if (d_dim == tbox::Dimension(1))
    {
        if (static_cast<int>(data_vel.size()) != 1)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelBoundaryUtilitiesSingleSpecies::readAdiabaticNoSlip()\n"
                << "'velocity' entry from '"
                << db_name
                << "' input database has incorrect size."
                << std::endl);
        }
        
        d_bdry_node_adiabatic_no_slip_vel[bdry_location_index] = data_vel[0];
    }
    else if (d_dim == tbox::Dimension(2))
    {
        if (static_cast<int>(data_vel.size()) != 2)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelBoundaryUtilitiesSingleSpecies::readAdiabaticNoSlip()\n"
                << "'velocity' entry from '"
                << db_name
                << "' input database has incorrect size."
                << std::endl);
        }
        
        d_bdry_edge_adiabatic_no_slip_vel[bdry_location_index*2] = data_vel[0];
        d_bdry_edge_adiabatic_no_slip_vel[bdry_location_index*2 + 1] = data_vel[1];
    }
    else if (d_dim == tbox::Dimension(3))
    {
        if (static_cast<int>(data_vel.size()) != 3)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelBoundaryUtilitiesSingleSpecies::readAdiabaticNoSlip()\n"
                << "'velocity' entry from '"
                << db_name
                << "' input database has incorrect size."
                << std::endl);
        }
        
        d_bdry_face_adiabatic_no_slip_vel[bdry_location_index*3] = data_vel[0];
        d_bdry_face_adiabatic_no_slip_vel[bdry_location_index*3 + 1] = data_vel[1];
        d_bdry_face_adiabatic_no_slip_vel[bdry_location_index*3 + 2] = data_vel[2];
    }
}


void
FlowModelBoundaryUtilitiesSingleSpecies::readIsothermalNoSlip(
    const HAMERS_SHARED_PTR<tbox::Database>& db,
    std::string& db_name,
    int bdry_location_index)
{
    TBOX_ASSERT(db);
    TBOX_ASSERT(!db_name.empty());
    
    Real data_T = Real(0);
    std::vector<Real> data_vel;
    
    if (db->keyExists("temperature"))
    {
        data_T = db->getFloat("temperature");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBoundaryUtilitiesSingleSpecies::readIsothermalNoSlip()\n"
            << "'temperature' entry missing from '"
            << db_name
            << "' input database."
            << std::endl);
    }
    
    if (db->keyExists("velocity"))
    {
        data_vel = db->getRealVector("velocity");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBoundaryUtilitiesSingleSpecies::readIsothermalNoSlip()\n"
            << "'velocity' entry missing from '"
            << db_name
            << "' input database."
            << std::endl);
    }
    
    if (d_dim == tbox::Dimension(1))
    {
        if (static_cast<int>(data_vel.size()) != 1)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelBoundaryUtilitiesSingleSpecies::readIsothermalNoSlip()\n"
                << "'velocity' entry from '"
                << db_name
                << "' input database has incorrect size."
                << std::endl);
        }
        
        d_bdry_node_isothermal_no_slip_T[bdry_location_index] = data_T;
        d_bdry_node_isothermal_no_slip_vel[bdry_location_index] = data_vel[0];
    }
    else if (d_dim == tbox::Dimension(2))
    {
        if (static_cast<int>(data_vel.size()) != 2)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelBoundaryUtilitiesSingleSpecies::readIsothermalNoSlip()\n"
                << "'velocity' entry from '"
                << db_name
                << "' input database has incorrect size."
                << std::endl);
        }
        
        d_bdry_edge_isothermal_no_slip_T[bdry_location_index] = data_T;
        d_bdry_edge_isothermal_no_slip_vel[bdry_location_index*2] = data_vel[0];
        d_bdry_edge_isothermal_no_slip_vel[bdry_location_index*2 + 1] = data_vel[1];
    }
    else if (d_dim == tbox::Dimension(3))
    {
        if (static_cast<int>(data_vel.size()) != 3)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelBoundaryUtilitiesSingleSpecies::readIsothermalNoSlip()\n"
                << "'velocity' entry from '"
                << db_name
                << "' input database has incorrect size."
                << std::endl);
        }
        
        d_bdry_face_isothermal_no_slip_T[bdry_location_index] = data_T;
        d_bdry_face_isothermal_no_slip_vel[bdry_location_index*3] = data_vel[0];
        d_bdry_face_isothermal_no_slip_vel[bdry_location_index*3 + 1] = data_vel[1];
        d_bdry_face_isothermal_no_slip_vel[bdry_location_index*3 + 2] = data_vel[2];
    }
}

void
FlowModelBoundaryUtilitiesSingleSpecies::readNonreflectingOutflow(
    const HAMERS_SHARED_PTR<tbox::Database>& db,
    std::string& db_name,
    int bdry_location_index)
{
    TBOX_ASSERT(db);
    TBOX_ASSERT(!db_name.empty());
    
    Real p_t = Real(0);
    Real sigma = Real(1)/Real(4); // 0.25
    Real beta = Real(0);
    Real length_char = Real(0);
    
    if (db->keyExists("pressure_target"))
    {
        p_t = db->getFloat("pressure_target");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBoundaryUtilitiesSingleSpecies::readNonreflectingOutflow()\n"
            << "'pressure_target' entry missing from '"
            << db_name
            << "' input database."
            << std::endl);
    }
    
    if (db->keyExists("sigma"))
    {
        sigma = db->getFloat("sigma");
    }
    
    if (db->keyExists("beta"))
    {
        beta = db->getFloat("beta");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBoundaryUtilitiesSingleSpecies::readNonreflectingOutflow()\n"
            << "'beta' entry missing from '"
            << db_name
            << "' input database."
            << std::endl);
    }
    
    if (db->keyExists("length_char"))
    {
        length_char = db->getFloat("length_char");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBoundaryUtilitiesSingleSpecies::readNonreflectingOutflow()\n"
            << "'length_char' entry missing from '"
            << db_name
            << "' input database."
            << std::endl);
    }
    
    if (d_dim == tbox::Dimension(1))
    {
        d_bdry_node_nonreflecting_outflow_p_t[bdry_location_index]         = p_t;
        d_bdry_node_nonreflecting_outflow_sigma[bdry_location_index]       = sigma;
        d_bdry_node_nonreflecting_outflow_beta[bdry_location_index]        = beta;
        d_bdry_node_nonreflecting_outflow_length_char[bdry_location_index] = length_char;
    }
    else if (d_dim == tbox::Dimension(2))
    {
        d_bdry_edge_nonreflecting_outflow_p_t[bdry_location_index]         = p_t;
        d_bdry_edge_nonreflecting_outflow_sigma[bdry_location_index]       = sigma;
        d_bdry_edge_nonreflecting_outflow_beta[bdry_location_index]        = beta;
        d_bdry_edge_nonreflecting_outflow_length_char[bdry_location_index] = length_char;
    }
    else if (d_dim == tbox::Dimension(3))
    {
        d_bdry_face_nonreflecting_outflow_p_t[bdry_location_index]         = p_t;
        d_bdry_face_nonreflecting_outflow_sigma[bdry_location_index]       = sigma;
        d_bdry_face_nonreflecting_outflow_beta[bdry_location_index]        = beta;
        d_bdry_face_nonreflecting_outflow_length_char[bdry_location_index] = length_char;
    }
}