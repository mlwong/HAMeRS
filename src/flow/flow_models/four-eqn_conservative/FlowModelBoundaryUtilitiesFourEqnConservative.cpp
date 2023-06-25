#include "flow/flow_models/four-eqn_conservative/FlowModelBoundaryUtilitiesFourEqnConservative.hpp"

#include "util/basic_boundary_conditions/CartesianBoundaryDefines.hpp"

#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"

#include <algorithm>

// Integer constant for debugging improperly set boundary data.
#define BOGUS_BDRY_LOC (-9999)

FlowModelBoundaryUtilitiesFourEqnConservative::FlowModelBoundaryUtilitiesFourEqnConservative(
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
FlowModelBoundaryUtilitiesFourEqnConservative::getFromInput1d(
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
        TBOX_ERROR(": FlowModelBoundaryUtilitiesFourEqnConservative::getFromInput1d()\n"
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
FlowModelBoundaryUtilitiesFourEqnConservative::getFromInput2d(
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
        TBOX_ERROR(": FlowModelBoundaryUtilitiesFourEqnConservative::getFromInput2d()\n"
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
FlowModelBoundaryUtilitiesFourEqnConservative::getFromInput3d(
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
        TBOX_ERROR(
            ": BasicCartesianBoundaryUtilities3::getFromInput3d()\n"
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
FlowModelBoundaryUtilitiesFourEqnConservative::getEdgeLocationForNodeBdry(
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
FlowModelBoundaryUtilitiesFourEqnConservative::getFaceLocationForEdgeBdry(
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
FlowModelBoundaryUtilitiesFourEqnConservative::getFaceLocationForNodeBdry(
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
FlowModelBoundaryUtilitiesFourEqnConservative::fill1dNodeBoundaryData(
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
                
                if (bdry_node_conds[node_loc] == BDRY_COND::FLOW_MODEL::ADIABATIC_NO_SLIP)
                {
                    for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                    {
                        const int idx_cell_rho_Y = i + num_subghosts_conservative_var[0][0];
                        const int idx_cell_mom = i + num_subghosts_conservative_var[1][0];
                        const int idx_cell_E = i + num_subghosts_conservative_var[2][0];
                        
                        int idx_cell_pivot_rho_Y = idx_cell_rho_Y;
                        int idx_cell_pivot_mom = idx_cell_mom;
                        int idx_cell_pivot_E = idx_cell_E;
                        
                        if (node_loc == BDRY_LOC::XLO)
                        {
                            idx_cell_pivot_rho_Y = interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                num_subghosts_conservative_var[0][0];
                            
                            idx_cell_pivot_mom = interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                num_subghosts_conservative_var[1][0];
                            
                            idx_cell_pivot_E = interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                num_subghosts_conservative_var[2][0];
                            
                        }
                        else if (node_loc == BDRY_LOC::XHI)
                        {
                            idx_cell_pivot_rho_Y = interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                num_subghosts_conservative_var[0][0];
                            
                            idx_cell_pivot_mom = interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                num_subghosts_conservative_var[1][0];
                            
                            idx_cell_pivot_E = interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                num_subghosts_conservative_var[2][0];
                        }
                        
                        /*
                         * Compute the mixture density of the pivot.
                         */
                        
                        Real rho_pivot = Real(0);
                        for (int si = 0; si < d_num_species; si++)
                        {
                            rho_pivot += Q[si][idx_cell_pivot_rho_Y];
                        }
                        
                        /*
                         * Compute the mass fractions of the pivot.
                         */
                        
                        std::vector<Real> Y_pivot;
                        Y_pivot.reserve(d_num_species);
                        for (int si = 0; si < d_num_species; si++)
                        {
                            Y_pivot.push_back(Q[si][idx_cell_pivot_rho_Y]/rho_pivot);
                        }
                        
                        /*
                         * Get the pointers to the mass fractions of the pivot.
                         */
                        
                        std::vector<const Real*> Y_pivot_ptr;
                        Y_pivot_ptr.reserve(d_num_species);
                        for (int si = 0; si < d_num_species; si++)
                        {
                            Y_pivot_ptr.push_back(&Y_pivot[si]);
                        }
                        
                        /*
                         * Set the values for partial densities and momentum.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            Q[si][idx_cell_rho_Y] = Q[si][idx_cell_pivot_rho_Y];
                        }
                        Q[d_num_species][idx_cell_mom] = -Q[d_num_species][idx_cell_pivot_mom] +
                            Real(2)*rho_pivot*d_bdry_node_adiabatic_no_slip_vel[node_loc];
                        
                        /*
                         * Set the values for total internal energy.
                         */
                        
                        Real epsilon_pivot = (Q[d_num_species + 1][idx_cell_pivot_E] -
                            Real(1)/Real(2)*Q[d_num_species][idx_cell_pivot_mom]*Q[d_num_species][idx_cell_pivot_mom]/
                            rho_pivot)/rho_pivot;
                        
                        Real p_pivot = d_equation_of_state_mixing_rules->
                            getPressure(
                                &rho_pivot,
                                &epsilon_pivot,
                                Y_pivot_ptr);
                        
                        Real T_pivot = d_equation_of_state_mixing_rules->
                            getTemperature(
                                &rho_pivot,
                                &p_pivot,
                                Y_pivot_ptr);
                        
                        Real T = T_pivot;
                        
                        Real epsilon = d_equation_of_state_mixing_rules->
                            getInternalEnergyFromTemperature(
                                &rho_pivot,
                                &T,
                                Y_pivot_ptr);
                        
                        Real E = rho_pivot*epsilon +
                            Real(1)/Real(2)*Q[d_num_species][idx_cell_mom]*Q[d_num_species][idx_cell_mom]/rho_pivot;
                        
                        Q[d_num_species + 1][idx_cell_E] = E;
                    }
                    
                    // Remove node locations that have boundary conditions identified.
                    bdry_node_locs.erase(std::remove(bdry_node_locs.begin(), bdry_node_locs.end(), node_loc),
                        bdry_node_locs.end());
                }
                else if (bdry_node_conds[node_loc] == BDRY_COND::FLOW_MODEL::ISOTHERMAL_NO_SLIP)
                {
                    for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                    {
                        const int idx_cell_rho_Y = i + num_subghosts_conservative_var[0][0];
                        const int idx_cell_mom = i + num_subghosts_conservative_var[1][0];
                        const int idx_cell_E = i + num_subghosts_conservative_var[2][0];
                        
                        int idx_cell_pivot_rho_Y = idx_cell_rho_Y;
                        int idx_cell_pivot_mom = idx_cell_mom;
                        int idx_cell_pivot_E = idx_cell_E;
                        
                        if (node_loc == BDRY_LOC::XLO)
                        {
                            idx_cell_pivot_rho_Y = interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                num_subghosts_conservative_var[0][0];
                            
                            idx_cell_pivot_mom = interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                num_subghosts_conservative_var[1][0];
                            
                            idx_cell_pivot_E = interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                num_subghosts_conservative_var[2][0];
                        }
                        else if (node_loc == BDRY_LOC::XHI)
                        {
                            idx_cell_pivot_rho_Y = interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                num_subghosts_conservative_var[0][0];
                            
                            idx_cell_pivot_mom = interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                num_subghosts_conservative_var[1][0];
                            
                            idx_cell_pivot_E = interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                num_subghosts_conservative_var[2][0];
                        }
                        
                        /*
                         * Compute the mixture density of the pivot.
                         */
                        
                        Real rho_pivot = Real(0);
                        for (int si = 0; si < d_num_species; si++)
                        {
                            rho_pivot += Q[si][idx_cell_pivot_rho_Y];
                        }
                        
                        /*
                         * Compute the mass fractions of the pivot.
                         */
                        
                        std::vector<Real> Y_pivot;
                        Y_pivot.reserve(d_num_species);
                        for (int si = 0; si < d_num_species; si++)
                        {
                            Y_pivot.push_back(Q[si][idx_cell_pivot_rho_Y]/rho_pivot);
                        }
                        
                        /*
                         * Get the pointers to the mass fractions of the pivot.
                         */
                        
                        std::vector<const Real*> Y_pivot_ptr;
                        Y_pivot_ptr.reserve(d_num_species);
                        for (int si = 0; si < d_num_species; si++)
                        {
                            Y_pivot_ptr.push_back(&Y_pivot[si]);
                        }
                        
                        /*
                         * Set the values for partial densities, momentum and total internal energy.
                         */
                        
                        Real epsilon_pivot = (Q[d_num_species + 1][idx_cell_pivot_E] -
                            Real(1)/Real(2)*Q[d_num_species][idx_cell_pivot_mom]*Q[d_num_species][idx_cell_pivot_mom]/
                            rho_pivot)/rho_pivot;
                        
                        Real p_pivot = d_equation_of_state_mixing_rules->
                            getPressure(
                                &rho_pivot,
                                &epsilon_pivot,
                                Y_pivot_ptr);
                        
                        Real p = p_pivot;
                        
                        Real T_pivot = d_equation_of_state_mixing_rules->
                            getTemperature(
                                &rho_pivot,
                                &p_pivot,
                                Y_pivot_ptr);
                        
                        Real T = -T_pivot + Real(2)*d_bdry_node_isothermal_no_slip_T[node_loc];
                        
                        Real rho = d_equation_of_state_mixing_rules->
                            getMixtureDensity(
                                &p,
                                &T,
                                Y_pivot_ptr);
                        
                        Real u = -Q[d_num_species][idx_cell_pivot_mom]/rho_pivot +
                            Real(2)*d_bdry_node_isothermal_no_slip_vel[node_loc];
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            Q[si][idx_cell_rho_Y] = rho*Y_pivot[si];
                        }
                        Q[d_num_species][idx_cell_mom] = rho*u;
                        
                        Real epsilon = d_equation_of_state_mixing_rules->
                            getInternalEnergyFromTemperature(
                                &rho,
                                &T,
                                Y_pivot_ptr);
                        
                        Real E = rho*epsilon + Real(1)/Real(2)*Q[d_num_species][idx_cell_mom]*Q[d_num_species][idx_cell_mom]/rho;
                        
                        Q[d_num_species + 1][idx_cell_E] = E;
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
            TBOX_ERROR("FlowModelBoundaryUtilitiesFourEqnConservative::fill1dNodeBoundaryData()\n"
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
FlowModelBoundaryUtilitiesFourEqnConservative::fill2dEdgeBoundaryData(
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
                
                if (bdry_edge_conds[edge_loc] == BDRY_COND::FLOW_MODEL::ADIABATIC_NO_SLIP)
                {
                    for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                    {
                        for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                        {
                            const int idx_cell_rho_Y = (i + num_subghosts_conservative_var[0][0]) +
                                (j + num_subghosts_conservative_var[0][1])*
                                    subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_mom = (i + num_subghosts_conservative_var[1][0]) +
                                (j + num_subghosts_conservative_var[1][1])*
                                    subghostcell_dims_conservative_var[1][0];
                            
                            const int idx_cell_E = (i + num_subghosts_conservative_var[2][0]) +
                                (j + num_subghosts_conservative_var[2][1])*
                                    subghostcell_dims_conservative_var[2][0];
                            
                            int idx_cell_pivot_rho_Y = idx_cell_rho_Y;
                            int idx_cell_pivot_mom = idx_cell_mom;
                            int idx_cell_pivot_E = idx_cell_E;
                            
                            if (edge_loc == BDRY_LOC::XLO)
                            {
                                idx_cell_pivot_rho_Y = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
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
                                idx_cell_pivot_rho_Y = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
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
                                idx_cell_pivot_rho_Y = (i + num_subghosts_conservative_var[0][0]) +
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
                                idx_cell_pivot_rho_Y = (i + num_subghosts_conservative_var[0][0]) +
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
                             * Compute the mixture density of the pivot.
                             */
                            
                            Real rho_pivot = Real(0);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                rho_pivot += Q[si][idx_cell_pivot_rho_Y];
                            }
                            
                            /*
                             * Compute the mass fractions of the pivot.
                             */
                            
                            std::vector<Real> Y_pivot;
                            Y_pivot.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_pivot.push_back(Q[si][idx_cell_pivot_rho_Y]/rho_pivot);
                            }
                            
                            /*
                             * Get the pointers to the mass fractions of the pivot.
                             */
                            
                            std::vector<const Real*> Y_pivot_ptr;
                            Y_pivot_ptr.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_pivot_ptr.push_back(&Y_pivot[si]);
                            }
                            
                            /*
                             * Set the values for partial densities and momentum.
                             */
                            
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Q[si][idx_cell_rho_Y] = Q[si][idx_cell_pivot_rho_Y];
                            }
                            Q[d_num_species][idx_cell_mom] = -Q[d_num_species][idx_cell_pivot_mom] +
                                Real(2)*rho_pivot*d_bdry_edge_adiabatic_no_slip_vel[edge_loc*2];
                            Q[d_num_species + 1][idx_cell_mom] = -Q[d_num_species + 1][idx_cell_pivot_mom] +
                                Real(2)*rho_pivot*d_bdry_edge_adiabatic_no_slip_vel[edge_loc*2 + 1];
                            
                            /*
                             * Set the values for total internal energy.
                             */
                            
                            Real epsilon_pivot = (Q[d_num_species + 2][idx_cell_pivot_E] -
                                Real(1)/Real(2)*(Q[d_num_species][idx_cell_pivot_mom]*Q[d_num_species][idx_cell_pivot_mom] +
                                Q[d_num_species + 1][idx_cell_pivot_mom]*Q[d_num_species + 1][idx_cell_pivot_mom])/
                                rho_pivot)/rho_pivot;
                            
                            Real p_pivot = d_equation_of_state_mixing_rules->
                                getPressure(
                                    &rho_pivot,
                                    &epsilon_pivot,
                                    Y_pivot_ptr);
                            
                            Real T_pivot = d_equation_of_state_mixing_rules->
                                getTemperature(
                                    &rho_pivot,
                                    &p_pivot,
                                    Y_pivot_ptr);
                            
                            Real T = T_pivot;
                            
                            Real epsilon = d_equation_of_state_mixing_rules->
                                getInternalEnergyFromTemperature(
                                    &rho_pivot,
                                    &T,
                                    Y_pivot_ptr);
                            
                            Real E = rho_pivot*epsilon +
                                Real(1)/Real(2)*(Q[d_num_species][idx_cell_mom]*Q[d_num_species][idx_cell_mom] +
                                Q[d_num_species + 1][idx_cell_mom]*Q[d_num_species + 1][idx_cell_mom])/rho_pivot;
                            
                            Q[d_num_species + 2][idx_cell_E] = E;
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
                            const int idx_cell_rho_Y = (i + num_subghosts_conservative_var[0][0]) +
                                (j + num_subghosts_conservative_var[0][1])*
                                    subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_mom = (i + num_subghosts_conservative_var[1][0]) +
                                (j + num_subghosts_conservative_var[1][1])*
                                    subghostcell_dims_conservative_var[1][0];
                            
                            const int idx_cell_E = (i + num_subghosts_conservative_var[2][0]) +
                                (j + num_subghosts_conservative_var[2][1])*
                                    subghostcell_dims_conservative_var[2][0];
                            
                            int idx_cell_pivot_rho_Y = idx_cell_rho_Y;
                            int idx_cell_pivot_mom = idx_cell_mom;
                            int idx_cell_pivot_E = idx_cell_E;
                            
                            if (edge_loc == BDRY_LOC::XLO)
                            {
                                idx_cell_pivot_rho_Y = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
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
                                idx_cell_pivot_rho_Y = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
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
                                idx_cell_pivot_rho_Y = (i + num_subghosts_conservative_var[0][0]) +
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
                                idx_cell_pivot_rho_Y = (i + num_subghosts_conservative_var[0][0]) +
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
                             * Compute the mixture density of the pivot.
                             */
                            
                            Real rho_pivot = Real(0);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                rho_pivot += Q[si][idx_cell_pivot_rho_Y];
                            }
                            
                            /*
                             * Compute the mass fractions of the pivot.
                             */
                            
                            std::vector<Real> Y_pivot;
                            Y_pivot.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_pivot.push_back(Q[si][idx_cell_pivot_rho_Y]/rho_pivot);
                            }
                            
                            /*
                             * Get the pointers to the mass fractions of the pivot.
                             */
                            
                            std::vector<const Real*> Y_pivot_ptr;
                            Y_pivot_ptr.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_pivot_ptr.push_back(&Y_pivot[si]);
                            }
                            
                            /*
                             * Set the values for partial densities, momentum and total internal energy.
                             */
                            
                            Real epsilon_pivot = (Q[d_num_species + 2][idx_cell_pivot_E] -
                                Real(1)/Real(2)*(Q[d_num_species][idx_cell_pivot_mom]*Q[d_num_species][idx_cell_pivot_mom] +
                                Q[d_num_species + 1][idx_cell_pivot_mom]*Q[d_num_species + 1][idx_cell_pivot_mom])/
                                rho_pivot)/rho_pivot;
                            
                            Real p_pivot = d_equation_of_state_mixing_rules->
                                getPressure(
                                    &rho_pivot,
                                    &epsilon_pivot,
                                    Y_pivot_ptr);
                            
                            Real p = p_pivot;
                            
                            Real T_pivot = d_equation_of_state_mixing_rules->
                                getTemperature(
                                    &rho_pivot,
                                    &p_pivot,
                                    Y_pivot_ptr);
                            
                            Real T = -T_pivot + Real(2)*d_bdry_edge_isothermal_no_slip_T[edge_loc];
                            
                            Real rho = d_equation_of_state_mixing_rules->
                                getMixtureDensity(
                                    &p,
                                    &T,
                                    Y_pivot_ptr);
                            
                            Real u = -Q[d_num_species][idx_cell_pivot_mom]/rho_pivot +
                                Real(2)*d_bdry_edge_isothermal_no_slip_vel[edge_loc*2];
                            Real v = -Q[d_num_species + 1][idx_cell_pivot_mom]/rho_pivot +
                                Real(2)*d_bdry_edge_isothermal_no_slip_vel[edge_loc*2 + 1];
                            
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Q[si][idx_cell_rho_Y] = rho*Y_pivot[si];
                            }
                            Q[d_num_species][idx_cell_mom] = rho*u;
                            Q[d_num_species + 1][idx_cell_mom] = rho*v;
                            
                            Real epsilon = d_equation_of_state_mixing_rules->
                                getInternalEnergyFromTemperature(
                                    &rho,
                                    &T,
                                    Y_pivot_ptr);
                            
                            Real E = rho*epsilon +
                                Real(1)/Real(2)*(Q[d_num_species][idx_cell_mom]*Q[d_num_species][idx_cell_mom] +
                                Q[d_num_species + 1][idx_cell_mom]*Q[d_num_species + 1][idx_cell_mom])/rho;
                            
                            Q[d_num_species + 2][idx_cell_E] = E;
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
                
                if (edge_loc == BDRY_LOC::XLO)
                {
                    const int num_ghosts_to_fill = fill_box_hi_idx[0] - fill_box_lo_idx[0] + 1;
                    TBOX_ASSERT(fill_box_hi_idx[0] == interior_box_lo_idx[0] - 1);
                    if (num_ghosts_to_fill > 6)
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelBoundaryUtilitiesFourEqnConservative::fill2dEdgeBoundaryData()\n"
                            << "Non-reflecting outflow BC doesn't support more than six ghost cells yet!");
                    }
                     
                    for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                    {
                        // Get the grid spacing.
                        const double* const dx = patch_geom->getDx();
                        
                        const int idx_cell_rho_Y_x_R = (interior_box_lo_idx[0] + num_subghosts_conservative_var[0][0]) +
                            (j + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0];
                        
                        const int idx_cell_rho_Y_x_RR = (interior_box_lo_idx[0] + 1 + num_subghosts_conservative_var[0][0]) +
                            (j + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0];
                        
                        const int idx_cell_rho_Y_x_RRR = (interior_box_lo_idx[0] + 2 + num_subghosts_conservative_var[0][0]) +
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
                        
                        std::vector<Real> rho_Y_x_R;
                        std::vector<Real> rho_Y_x_RR;
                        std::vector<Real> rho_Y_x_RRR;
                        rho_Y_x_R.reserve(d_num_species);
                        rho_Y_x_RR.reserve(d_num_species);
                        rho_Y_x_RRR.reserve(d_num_species);
                        for (int si = 0; si < d_num_species; si++)
                        {
                            rho_Y_x_R.push_back(Q[si][idx_cell_rho_Y_x_R]);
                            rho_Y_x_RR.push_back(Q[si][idx_cell_rho_Y_x_RR]);
                            rho_Y_x_RRR.push_back(Q[si][idx_cell_rho_Y_x_RRR]);
                        }
                        
                        /*
                         * Compute the mixture density.
                         */
                        
                        Real rho_x_R   = Real(0);
                        Real rho_x_RR  = Real(0);
                        Real rho_x_RRR = Real(0);
                        for (int si = 0; si < d_num_species; si++)
                        {
                            rho_x_R   += Q[si][idx_cell_rho_Y_x_R];
                            rho_x_RR  += Q[si][idx_cell_rho_Y_x_RR];
                            rho_x_RRR += Q[si][idx_cell_rho_Y_x_RRR];
                        }
                        
                        /*
                         * Compute the mass fractions.
                         */
                        
                        std::vector<Real> Y_x_R;
                        std::vector<Real> Y_x_RR;
                        std::vector<Real> Y_x_RRR;
                        Y_x_R.reserve(d_num_species);
                        Y_x_RR.reserve(d_num_species);
                        Y_x_RRR.reserve(d_num_species);
                        for (int si = 0; si < d_num_species; si++)
                        {
                            Y_x_R.push_back(Q[si][idx_cell_rho_Y_x_R]/rho_x_R);
                            Y_x_RR.push_back(Q[si][idx_cell_rho_Y_x_RR]/rho_x_RR);
                            Y_x_RRR.push_back(Q[si][idx_cell_rho_Y_x_RRR]/rho_x_RRR);
                        }
                        
                        /*
                         * Get the pointers to the mass fractions.
                         */
                        
                        std::vector<const Real*> Y_x_R_ptr;
                        std::vector<const Real*> Y_x_RR_ptr;
                        std::vector<const Real*> Y_x_RRR_ptr;
                        Y_x_R_ptr.reserve(d_num_species);
                        Y_x_RR_ptr.reserve(d_num_species);
                        Y_x_RRR_ptr.reserve(d_num_species);
                        for (int si = 0; si < d_num_species; si++)
                        {
                            Y_x_R_ptr.push_back(&Y_x_R[si]);
                            Y_x_RR_ptr.push_back(&Y_x_RR[si]);
                            Y_x_RRR_ptr.push_back(&Y_x_RRR[si]);
                        }
                       
                        const Real u_x_R   = Q[d_num_species][idx_cell_mom_x_R]/rho_x_R;
                        const Real u_x_RR  = Q[d_num_species][idx_cell_mom_x_RR]/rho_x_RR;
                        const Real u_x_RRR = Q[d_num_species][idx_cell_mom_x_RRR]/rho_x_RRR;
                        
                        const Real v_x_R   = Q[d_num_species + 1][idx_cell_mom_x_R]/rho_x_R;
                        const Real v_x_RR  = Q[d_num_species + 1][idx_cell_mom_x_RR]/rho_x_RR;
                        const Real v_x_RRR = Q[d_num_species + 1][idx_cell_mom_x_RRR]/rho_x_RRR;
                        
                        const Real half = Real(1)/Real(2);
                        const Real epsilon_x_R   = Q[d_num_species + 2][idx_cell_E_x_R]/rho_x_R - half*(u_x_R*u_x_R + v_x_R*v_x_R);
                        const Real epsilon_x_RR  = Q[d_num_species + 2][idx_cell_E_x_RR]/rho_x_RR - half*(u_x_RR*u_x_RR + v_x_RR*v_x_RR);
                        const Real epsilon_x_RRR = Q[d_num_species + 2][idx_cell_E_x_RRR]/rho_x_RRR - half*(u_x_RRR*u_x_RRR + v_x_RRR*v_x_RRR);
                        
                        Real p_x_R = d_equation_of_state_mixing_rules->
                            getPressure(
                                &rho_x_R,
                                &epsilon_x_R,
                                Y_x_R_ptr);
                        
                        Real p_x_RR = d_equation_of_state_mixing_rules->
                            getPressure(
                                &rho_x_RR,
                                &epsilon_x_RR,
                                Y_x_RR_ptr);
                        
                        Real p_x_RRR = d_equation_of_state_mixing_rules->
                            getPressure(
                                &rho_x_RRR,
                                &epsilon_x_RRR,
                                Y_x_RRR_ptr);
                        
                        /*
                         * Compute derivatives in x-direction.
                         */
                        
                        std::vector<Real> drho_Y_dx;
                        drho_Y_dx.reserve(d_num_species);
                        for (int si = 0; si < d_num_species; si++)
                        {
                            drho_Y_dx.push_back(-(Q[si][idx_cell_rho_Y_x_RRR] - Real(4)*Q[si][idx_cell_rho_Y_x_RR] +
                                Real(3)*Q[si][idx_cell_rho_Y_x_R])/(Real(2)*Real(dx[0])));
                        }
                        const Real du_dx = -(u_x_RRR - Real(4)*u_x_RR + Real(3)*u_x_R)/(Real(2)*Real(dx[0]));
                        const Real dv_dx = -(v_x_RRR - Real(4)*v_x_RR + Real(3)*v_x_R)/(Real(2)*Real(dx[0]));
                        const Real dp_dx = -(p_x_RRR - Real(4)*p_x_RR + Real(3)*p_x_R)/(Real(2)*Real(dx[0]));
                        
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
                            
                            const int idx_cell_rho_Y_y_T = (interior_box_lo_idx[0] + num_subghosts_conservative_var[0][0]) +
                                (j + 1 + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_mom_y_T = (interior_box_lo_idx[0] + num_subghosts_conservative_var[1][0]) +
                                (j + 1 + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0];
                            
                            const int idx_cell_E_y_T = (interior_box_lo_idx[0] + num_subghosts_conservative_var[2][0]) +
                                (j + 1 + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0];
                            
                            /*
                             * Compute the mixture density.
                             */
                            
                            Real rho_y_T = Real(0);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                rho_y_T += Q[si][idx_cell_rho_Y_y_T];
                            }
                            
                            /*
                             * Compute the mass fractions.
                             */
                            
                            std::vector<Real> Y_y_T;
                            Y_y_T.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_y_T.push_back(Q[si][idx_cell_rho_Y_y_T]/rho_y_T);
                            }
                            
                            /*
                             * Get the pointers to the mass fractions.
                             */
                            
                            std::vector<const Real*> Y_y_T_ptr;
                            Y_y_T_ptr.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_y_T_ptr.push_back(&Y_y_T[si]);
                            }
                            
                            const Real u_y_T = Q[d_num_species][idx_cell_mom_y_T]/rho_y_T;
                            const Real v_y_T = Q[d_num_species + 1][idx_cell_mom_y_T]/rho_y_T;
                            const Real epsilon_y_T = Q[d_num_species + 2][idx_cell_E_y_T]/rho_y_T - half*(u_y_T*u_y_T + v_y_T*v_y_T);
                            
                            Real p_y_T = d_equation_of_state_mixing_rules->
                                getPressure(
                                    &rho_y_T,
                                    &epsilon_y_T,
                                    Y_y_T_ptr);
                            
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
                            
                            const int idx_cell_rho_Y_y_B = (interior_box_lo_idx[0] + num_subghosts_conservative_var[0][0]) +
                                (j - 1 + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_mom_y_B = (interior_box_lo_idx[0] + num_subghosts_conservative_var[1][0]) +
                                (j - 1 + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0];
                            
                            const int idx_cell_E_y_B = (interior_box_lo_idx[0] + num_subghosts_conservative_var[2][0]) +
                                (j - 1 + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0];
                            
                            /*
                             * Compute the mixture density.
                             */
                            
                            Real rho_y_B = Real(0);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                rho_y_B += Q[si][idx_cell_rho_Y_y_B];
                            }
                            
                            /*
                             * Compute the mass fractions.
                             */
                            
                            std::vector<Real> Y_y_B;
                            Y_y_B.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_y_B.push_back(Q[si][idx_cell_rho_Y_y_B]/rho_y_B);
                            }
                            
                            /*
                             * Get the pointers to the mass fractions.
                             */
                            
                            std::vector<const Real*> Y_y_B_ptr;
                            Y_y_B_ptr.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_y_B_ptr.push_back(&Y_y_B[si]);
                            }
                            
                            const Real u_y_B = Q[d_num_species][idx_cell_mom_y_B]/rho_y_B;
                            const Real v_y_B = Q[d_num_species + 1][idx_cell_mom_y_B]/rho_y_B;
                            const Real epsilon_y_B = Q[d_num_species + 2][idx_cell_E_y_B]/rho_y_B - half*(u_y_B*u_y_B + v_y_B*v_y_B);
                            
                            Real p_y_B = d_equation_of_state_mixing_rules->
                                getPressure(
                                    &rho_y_B,
                                    &epsilon_y_B,
                                    Y_y_B_ptr);
                            
                            // One-sided derivatives.
                            du_dy = (u_x_R - u_y_B)/Real(dx[1]);
                            dv_dy = (v_x_R - v_y_B)/Real(dx[1]);
                            dp_dy = (p_x_R - p_y_B)/Real(dx[1]);
                        }
                        else
                        {
                            const int idx_cell_rho_Y_y_B = (interior_box_lo_idx[0] + num_subghosts_conservative_var[0][0]) +
                                (j - 1 + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_rho_Y_y_T = (interior_box_lo_idx[0] + num_subghosts_conservative_var[0][0]) +
                                (j + 1 + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_mom_y_B = (interior_box_lo_idx[0] + num_subghosts_conservative_var[1][0]) +
                                (j - 1 + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0];
                            
                            const int idx_cell_mom_y_T = (interior_box_lo_idx[0] + num_subghosts_conservative_var[1][0]) +
                                (j + 1 + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0];
                            
                            const int idx_cell_E_y_B = (interior_box_lo_idx[0] + num_subghosts_conservative_var[2][0]) +
                                (j - 1 + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0];
                            
                            const int idx_cell_E_y_T = (interior_box_lo_idx[0] + num_subghosts_conservative_var[2][0]) +
                                (j + 1 + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0];
                            
                            /*
                             * Compute the mixture density.
                             */
                            
                            Real rho_y_B = Real(0);
                            Real rho_y_T = Real(0);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                rho_y_B += Q[si][idx_cell_rho_Y_y_B];
                                rho_y_T += Q[si][idx_cell_rho_Y_y_T];
                            }
                            
                            /*
                             * Compute the mass fractions.
                             */
                            
                            std::vector<Real> Y_y_B;
                            std::vector<Real> Y_y_T;
                            Y_y_B.reserve(d_num_species);
                            Y_y_T.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_y_B.push_back(Q[si][idx_cell_rho_Y_y_B]/rho_y_B);
                                Y_y_T.push_back(Q[si][idx_cell_rho_Y_y_T]/rho_y_T);
                            }
                            
                            /*
                             * Get the pointers to the mass fractions.
                             */
                            
                            std::vector<const Real*> Y_y_B_ptr;
                            std::vector<const Real*> Y_y_T_ptr;
                            Y_y_B_ptr.reserve(d_num_species);
                            Y_y_T_ptr.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_y_B_ptr.push_back(&Y_y_B[si]);
                                Y_y_T_ptr.push_back(&Y_y_T[si]);
                            }
                            
                            const Real u_y_B = Q[d_num_species][idx_cell_mom_y_B]/rho_y_B;
                            const Real u_y_T = Q[d_num_species][idx_cell_mom_y_T]/rho_y_T;
                            
                            const Real v_y_B = Q[d_num_species + 1][idx_cell_mom_y_B]/rho_y_B;
                            const Real v_y_T = Q[d_num_species + 1][idx_cell_mom_y_T]/rho_y_T;
                            
                            const Real epsilon_y_B = Q[d_num_species + 2][idx_cell_E_y_B]/rho_y_B - half*(u_y_B*u_y_B + v_y_B*v_y_B);
                            const Real epsilon_y_T = Q[d_num_species + 2][idx_cell_E_y_T]/rho_y_T - half*(u_y_T*u_y_T + v_y_T*v_y_T);
                            
                            Real p_y_B = d_equation_of_state_mixing_rules->
                                getPressure(
                                    &rho_y_B,
                                    &epsilon_y_B,
                                    Y_y_B_ptr);
                            
                            Real p_y_T = d_equation_of_state_mixing_rules->
                                getPressure(
                                    &rho_y_T,
                                    &epsilon_y_T,
                                    Y_y_T_ptr);
                            
                            // Central derivatives.
                            du_dy = (u_y_T - u_y_B)/(Real(2)*Real(dx[1]));
                            dv_dy = (v_y_T - v_y_B)/(Real(2)*Real(dx[1]));
                            dp_dy = (p_y_T - p_y_B)/(Real(2)*Real(dx[1]));
                        }
                        
                        // Compute sound speed.
                        
                        const Real Gamma_x_R = d_equation_of_state_mixing_rules->getGruneisenParameter(
                            &rho_x_R,
                            &p_x_R,
                            Y_x_R_ptr);
                        
                        const std::vector<Real> Psi_x_R = d_equation_of_state_mixing_rules->
                            getPressureDerivativeWithPartialDensities(
                                    &rho_x_R,
                                    &p_x_R,
                                    Y_x_R_ptr);
                        
                        Real c_x_R = Gamma_x_R*p_x_R/rho_x_R;
                        for (int si = 0; si < d_num_species; si++)
                        {
                            c_x_R += Y_x_R[si]*Psi_x_R[si];
                        }
                        c_x_R = std::sqrt(c_x_R);
                        
                        const Real lambda_last = u_x_R + c_x_R;

                        // Compute vector Lambda^(-1) * L.
                        
                        Real Lambda_inv_L[d_num_species + 3];
                        
                        const Real& p_t         = d_bdry_edge_nonreflecting_outflow_p_t[edge_loc];
                        const Real& sigma       = d_bdry_edge_nonreflecting_outflow_sigma[edge_loc];
                        const Real& beta        = d_bdry_edge_nonreflecting_outflow_beta[edge_loc];
                        const Real& length_char = d_bdry_edge_nonreflecting_outflow_length_char[edge_loc];
                        
                        const Real T_last = v_x_R*(dp_dy + rho_x_R*c_x_R*du_dy) + rho_x_R*c_x_R*c_x_R*dv_dy;
                        
                        const Real M_sq = (u_x_R*u_x_R + v_x_R*v_x_R)/(c_x_R*c_x_R);
                        const Real K = sigma*c_x_R*(Real(1) - M_sq)/length_char;
                        
                        Lambda_inv_L[0] = dp_dx - rho_x_R*c_x_R*du_dx;
                        for (int si = 0; si < d_num_species; si++)
                        {
                            Lambda_inv_L[si + 1] = c_x_R*c_x_R*drho_Y_dx[si] - Y_x_R[si]*dp_dx;
                        }
                        Lambda_inv_L[d_num_species + 1] = dv_dx;
                        Lambda_inv_L[d_num_species + 2] = (Real(1)/lambda_last)*(K*(p_x_R - p_t) - (Real(1) - beta)*T_last);
                        
                        // Compute dV_dx.
                        
                        const Real c_sq_inv  = Real(1)/(c_x_R*c_x_R);
                        const Real rho_c_inv = Real(1)/(rho_x_R*c_x_R);
                        
                        Real dV_dx[d_num_species + 3];
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            dV_dx[si] = half*c_sq_inv*Y_x_R[si]*(Lambda_inv_L[0] + Lambda_inv_L[d_num_species + 2]) +
                                c_sq_inv*Lambda_inv_L[si + 1];
                        }
                        dV_dx[d_num_species]     = half*rho_c_inv*(-Lambda_inv_L[0] + Lambda_inv_L[d_num_species + 2]);
                        dV_dx[d_num_species + 1] = Lambda_inv_L[d_num_species + 1];
                        dV_dx[d_num_species + 2] = half*(Lambda_inv_L[0] + Lambda_inv_L[d_num_species + 2]);
                        
                        Real V_ghost[(d_num_species + 3)*num_ghosts_to_fill];
                        
                        for (int i = num_ghosts_to_fill - 1; i >= 0; i--)
                        {
                            const int idx_cell_rho_Y = (i + fill_box_lo_idx[0] + num_subghosts_conservative_var[0][0]) +
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
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    V_ghost[i*(d_num_species + 3) + si] = rho_Y_x_RR[si] - Real(2)*Real(dx[0])*dV_dx[si];
                                }
                                
                                V_ghost[i*(d_num_species + 3) + d_num_species]     = u_x_RR - Real(2)*Real(dx[0])*dV_dx[d_num_species];
                                V_ghost[i*(d_num_species + 3) + d_num_species + 1] = v_x_RR - Real(2)*Real(dx[0])*dV_dx[d_num_species + 1];
                                V_ghost[i*(d_num_species + 3) + d_num_species + 2] = p_x_RR - Real(2)*Real(dx[0])*dV_dx[d_num_species + 2];
                            }
                            else if (i == num_ghosts_to_fill - 2)
                            {
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    V_ghost[i*(d_num_species + 3) + si] = -Real(2)*rho_Y_x_RR[si] - Real(3)*rho_Y_x_R[si] +
                                        Real(6)*V_ghost[(i + 1)*(d_num_species + 3) + si] + Real(6)*Real(dx[0])*dV_dx[si];
                                }
                                
                                V_ghost[i*(d_num_species + 3) + d_num_species] = -Real(2)*u_x_RR - Real(3)*u_x_R +
                                    Real(6)*V_ghost[(i + 1)*(d_num_species + 3) + d_num_species] +
                                    Real(6)*Real(dx[0])*dV_dx[d_num_species];
                                
                                V_ghost[i*(d_num_species + 3) + d_num_species + 1] = -Real(2)*v_x_RR - Real(3)*v_x_R +
                                    Real(6)*V_ghost[(i + 1)*(d_num_species + 3) + d_num_species + 1] +
                                    Real(6)*Real(dx[0])*dV_dx[d_num_species + 1];
                                
                                V_ghost[i*(d_num_species + 3) + d_num_species + 2] = -Real(2)*p_x_RR - Real(3)*p_x_R +
                                    Real(6)*V_ghost[(i + 1)*(d_num_species + 3) + d_num_species + 2] +
                                    Real(6)*Real(dx[0])*dV_dx[d_num_species + 2];
                            }
                            else if (i == num_ghosts_to_fill - 3)
                            {
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    V_ghost[i*(d_num_species + 3) + si] = Real(3)*rho_Y_x_RR[si] + Real(10)*rho_Y_x_R[si] -
                                        Real(18)*V_ghost[(i + 2)*(d_num_species + 3) + si] +
                                        Real(6)*V_ghost[(i + 1)*(d_num_species + 3) + si] -
                                        Real(12)*Real(dx[0])*dV_dx[si];
                                }
                                
                                V_ghost[i*(d_num_species + 3) + d_num_species] = Real(3)*u_x_RR + Real(10)*u_x_R -
                                    Real(18)*V_ghost[(i + 2)*(d_num_species + 3) + d_num_species] +
                                    Real(6)*V_ghost[(i + 1)*(d_num_species + 3) + d_num_species] -
                                    Real(12)*Real(dx[0])*dV_dx[d_num_species];
                                
                                V_ghost[i*(d_num_species + 3) + d_num_species + 1] = Real(3)*v_x_RR + Real(10)*v_x_R -
                                    Real(18)*V_ghost[(i + 2)*(d_num_species + 3) + d_num_species + 1] +
                                    Real(6)*V_ghost[(i + 1)*(d_num_species + 3) + d_num_species + 1] -
                                    Real(12)*Real(dx[0])*dV_dx[d_num_species + 1];
                                
                                V_ghost[i*(d_num_species + 3) + d_num_species + 2] = Real(3)*p_x_RR + Real(10)*p_x_R -
                                    Real(18)*V_ghost[(i + 2)*(d_num_species + 3) + d_num_species + 2] +
                                    Real(6)*V_ghost[(i + 1)*(d_num_species + 3) + d_num_species + 2] -
                                    Real(12)*Real(dx[0])*dV_dx[d_num_species + 2];
                            }
                            else if (i == num_ghosts_to_fill - 4)
                            {
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    V_ghost[i*(d_num_species + 3) + si] = -Real(4)*rho_Y_x_RR[si] -
                                        Real(65)/Real(3)*rho_Y_x_R[si] +
                                        Real(40)*V_ghost[(i + 3)*(d_num_species + 3) + si] -
                                        Real(20)*V_ghost[(i + 2)*(d_num_species + 3) + si] +
                                        Real(20)/Real(3)*V_ghost[(i + 1)*(d_num_species + 3) + si] +
                                        Real(20)*Real(dx[0])*dV_dx[si];
                                }
                                
                                V_ghost[i*(d_num_species + 3) + d_num_species] = -Real(4)*u_x_RR -
                                    Real(65)/Real(3)*u_x_R +
                                    Real(40)*V_ghost[(i + 3)*(d_num_species + 3) + d_num_species] -
                                    Real(20)*V_ghost[(i + 2)*(d_num_species + 3) + d_num_species] +
                                    Real(20)/Real(3)*V_ghost[(i + 1)*(d_num_species + 3) + d_num_species] +
                                    Real(20)*Real(dx[0])*dV_dx[d_num_species];
                                
                                V_ghost[i*(d_num_species + 3) + d_num_species + 1] = -Real(4)*v_x_RR -
                                    Real(65)/Real(3)*v_x_R +
                                    Real(40)*V_ghost[(i + 3)*(d_num_species + 3) + d_num_species + 1] -
                                    Real(20)*V_ghost[(i + 2)*(d_num_species + 3) + d_num_species + 1] +
                                    Real(20)/Real(3)*V_ghost[(i + 1)*(d_num_species + 3) + d_num_species + 1] +
                                    Real(20)*Real(dx[0])*dV_dx[d_num_species + 1];
                                
                                V_ghost[i*(d_num_species + 3) + d_num_species + 2] = -Real(4)*p_x_RR -
                                    Real(65)/Real(3)*p_x_R +
                                    Real(40)*V_ghost[(i + 3)*(d_num_species + 3) + d_num_species + 2] -
                                    Real(20)*V_ghost[(i + 2)*(d_num_species + 3) + d_num_species + 2] +
                                    Real(20)/Real(3)*V_ghost[(i + 1)*(d_num_species + 3) + d_num_species + 2] +
                                    Real(20)*Real(dx[0])*dV_dx[d_num_species + 2];
                            }
                            else if (i == num_ghosts_to_fill - 5)
                            {
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    V_ghost[i*(d_num_species + 3) + si] = Real(5)*rho_Y_x_RR[si] +
                                        Real(77)/Real(2)*rho_Y_x_R[si] -
                                        Real(75)*V_ghost[(i + 4)*(d_num_species + 3) + si] +
                                        Real(50)*V_ghost[(i + 3)*(d_num_species + 3) + si] -
                                        Real(25)*V_ghost[(i + 2)*(d_num_species + 3) + si] +
                                        Real(15)/Real(2)*V_ghost[(i + 1)*(d_num_species + 3) + si] -
                                        Real(30)*Real(dx[0])*dV_dx[si];
                                }
                                
                                V_ghost[i*(d_num_species + 3) + d_num_species] = Real(5)*u_x_RR +
                                    Real(77)/Real(2)*u_x_R -
                                    Real(75)*V_ghost[(i + 4)*(d_num_species + 3) + d_num_species] +
                                    Real(50)*V_ghost[(i + 3)*(d_num_species + 3) + d_num_species] -
                                    Real(25)*V_ghost[(i + 2)*(d_num_species + 3) + d_num_species] +
                                    Real(15)/Real(2)*V_ghost[(i + 1)*(d_num_species + 3) + d_num_species] -
                                    Real(30)*Real(dx[0])*dV_dx[d_num_species];
                                
                                V_ghost[i*(d_num_species + 3) + d_num_species + 1] = Real(5)*v_x_RR +
                                    Real(77)/Real(2)*v_x_R -
                                    Real(75)*V_ghost[(i + 4)*(d_num_species + 3) + d_num_species + 1] +
                                    Real(50)*V_ghost[(i + 3)*(d_num_species + 3) + d_num_species + 1] -
                                    Real(25)*V_ghost[(i + 2)*(d_num_species + 3) + d_num_species + 1] +
                                    Real(15)/Real(2)*V_ghost[(i + 1)*(d_num_species + 3) + d_num_species + 1] -
                                    Real(30)*Real(dx[0])*dV_dx[d_num_species + 1];
                                
                                V_ghost[i*(d_num_species + 3) + d_num_species + 2] = Real(5)*p_x_RR +
                                    Real(77)/Real(2)*p_x_R -
                                    Real(75)*V_ghost[(i + 4)*(d_num_species + 3) + d_num_species + 2] +
                                    Real(50)*V_ghost[(i + 3)*(d_num_species + 3) + d_num_species + 2] -
                                    Real(25)*V_ghost[(i + 2)*(d_num_species + 3) + d_num_species + 2] +
                                    Real(15)/Real(2)*V_ghost[(i + 1)*(d_num_species + 3) + d_num_species + 2] -
                                    Real(30)*Real(dx[0])*dV_dx[d_num_species + 2];
                            }
                            else if (i == num_ghosts_to_fill - 6)
                            {
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    V_ghost[i*(d_num_species + 3) + si] = -Real(6)*rho_Y_x_RR[si] -
                                        Real(609)/Real(10)*rho_Y_x_R[si] +
                                        Real(126)*V_ghost[(i + 5)*(d_num_species + 3) + si] -
                                        Real(105)*V_ghost[(i + 4)*(d_num_species + 3) + si] +
                                        Real(70)*V_ghost[(i + 3)*(d_num_species + 3) + si] -
                                        Real(63)/Real(2)*V_ghost[(i + 2)*(d_num_species + 3) + si] +
                                        Real(42)/Real(5)*V_ghost[(i + 1)*(d_num_species + 3) + si] +
                                        Real(42)*Real(dx[0])*dV_dx[si];
                                }
                                
                                V_ghost[i*(d_num_species + 3) + d_num_species] = -Real(6)*u_x_RR -
                                    Real(609)/Real(10)*u_x_R +
                                    Real(126)*V_ghost[(i + 5)*(d_num_species + 3) + d_num_species] -
                                    Real(105)*V_ghost[(i + 4)*(d_num_species + 3) + d_num_species] +
                                    Real(70)*V_ghost[(i + 3)*(d_num_species + 3) + d_num_species] -
                                    Real(63)/Real(2)*V_ghost[(i + 2)*(d_num_species + 3) + d_num_species] +
                                    Real(42)/Real(5)*V_ghost[(i + 1)*(d_num_species + 3) + d_num_species] +
                                    Real(42)*Real(dx[0])*dV_dx[d_num_species];
                                
                                V_ghost[i*(d_num_species + 3) + d_num_species + 1] = -Real(6)*v_x_RR -
                                    Real(609)/Real(10)*v_x_R +
                                    Real(126)*V_ghost[(i + 5)*(d_num_species + 3) + d_num_species + 1] -
                                    Real(105)*V_ghost[(i + 4)*(d_num_species + 3) + d_num_species + 1] +
                                    Real(70)*V_ghost[(i + 3)*(d_num_species + 3) + d_num_species + 1] -
                                    Real(63)/Real(2)*V_ghost[(i + 2)*(d_num_species + 3) + d_num_species + 1] +
                                    Real(42)/Real(5)*V_ghost[(i + 1)*(d_num_species + 3) + d_num_species + 1] +
                                    Real(42)*Real(dx[0])*dV_dx[d_num_species + 1];
                                
                                V_ghost[i*(d_num_species + 3) + d_num_species + 2] = -Real(6)*p_x_RR -
                                    Real(609)/Real(10)*p_x_R +
                                    Real(126)*V_ghost[(i + 5)*(d_num_species + 3) + d_num_species + 2] -
                                    Real(105)*V_ghost[(i + 4)*(d_num_species + 3) + d_num_species + 2] +
                                    Real(70)*V_ghost[(i + 3)*(d_num_species + 3) + d_num_species + 2] -
                                    Real(63)/Real(2)*V_ghost[(i + 2)*(d_num_species + 3) + d_num_species + 2] +
                                    Real(42)/Real(5)*V_ghost[(i + 1)*(d_num_species + 3) + d_num_species + 2] +
                                    Real(42)*Real(dx[0])*dV_dx[d_num_species + 2];
                            }
                            
                            /*
                             * Compute the mixture density.
                             */
                            
                            Real rho_ghost = Real(0);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                rho_ghost += V_ghost[i*(d_num_species + 3) + si];
                            }
                            
                            /*
                             * Compute the mass fractions.
                             */
                            
                            std::vector<Real> Y_ghost;
                            Y_ghost.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_ghost.push_back(V_ghost[i*(d_num_species + 3) + si]/rho_ghost);
                            }
                            
                            /*
                             * Get the pointers to the mass fractions.
                             */
                            
                            std::vector<const Real*> Y_ghost_ptr;
                            Y_ghost_ptr.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_ghost_ptr.push_back(&Y_ghost[si]);
                            }
                            
                            for(int si=0; si < d_num_species; si++)
                            {
                                Q[si][idx_cell_rho_Y] = V_ghost[i*(d_num_species + 3) + si];
                            }
                            
                            Q[d_num_species][idx_cell_mom]     = rho_ghost*V_ghost[i*(d_num_species + 3) + d_num_species];
                            Q[d_num_species + 1][idx_cell_mom] = rho_ghost*V_ghost[i*(d_num_species + 3) + d_num_species + 1];
                            
                            const Real epsilon = d_equation_of_state_mixing_rules->
                                getInternalEnergy(
                                    &rho_ghost,
                                    &V_ghost[i*(d_num_species + 3) + d_num_species + 2],
                                    Y_ghost_ptr);
                            
                            const Real E = rho_ghost*epsilon +
                                half*(Q[d_num_species][idx_cell_mom]*Q[d_num_species][idx_cell_mom] +
                                    Q[d_num_species + 1][idx_cell_mom]*Q[d_num_species + 1][idx_cell_mom])/
                                    rho_ghost;
                            
                            Q[d_num_species + 2][idx_cell_E] = E;
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
                            << ": FlowModelBoundaryUtilitiesFourEqnConservative::fill2dEdgeBoundaryData()\n"
                            << "Non-reflecting outflow BC doesn't support more than six ghost cells yet!");
                    }
                    
                    for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                    {
                        // Get the grid spacing.
                        const double* const dx = patch_geom->getDx();
                        
                        const int idx_cell_rho_Y_x_L = (interior_box_hi_idx[0] + num_subghosts_conservative_var[0][0]) +
                            (j + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0];
                        
                        const int idx_cell_rho_Y_x_LL = (interior_box_hi_idx[0] - 1 + num_subghosts_conservative_var[0][0]) +
                            (j + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0];
                        
                        const int idx_cell_rho_Y_x_LLL = (interior_box_hi_idx[0] - 2 + num_subghosts_conservative_var[0][0]) +
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

                        std::vector<Real> rho_Y_x_L;
                        std::vector<Real> rho_Y_x_LL;
                        std::vector<Real> rho_Y_x_LLL;
                        rho_Y_x_L.reserve(d_num_species);
                        rho_Y_x_LL.reserve(d_num_species);
                        rho_Y_x_LLL.reserve(d_num_species);
                        for (int si = 0; si < d_num_species; si++)
                        {
                            rho_Y_x_L.push_back(Q[si][idx_cell_rho_Y_x_L]);
                            rho_Y_x_LL.push_back(Q[si][idx_cell_rho_Y_x_LL]);
                            rho_Y_x_LLL.push_back(Q[si][idx_cell_rho_Y_x_LLL]);
                        }
                        
                        /*
                         * Compute the mixture density.
                         */
                        
                        Real rho_x_L   = Real(0);
                        Real rho_x_LL  = Real(0);
                        Real rho_x_LLL = Real(0);
                        for (int si = 0; si < d_num_species; si++)
                        {
                            rho_x_L   += Q[si][idx_cell_rho_Y_x_L];
                            rho_x_LL  += Q[si][idx_cell_rho_Y_x_LL];
                            rho_x_LLL += Q[si][idx_cell_rho_Y_x_LLL];
                        }
                        
                        /*
                         * Compute the mass fractions.
                         */

                        std::vector<Real> Y_x_L;
                        std::vector<Real> Y_x_LL;
                        std::vector<Real> Y_x_LLL;
                        Y_x_L.reserve(d_num_species);
                        Y_x_LL.reserve(d_num_species);
                        Y_x_LLL.reserve(d_num_species);
                        for (int si = 0; si < d_num_species; si++)
                        {
                            Y_x_L.push_back(Q[si][idx_cell_rho_Y_x_L]/rho_x_L);
                            Y_x_LL.push_back(Q[si][idx_cell_rho_Y_x_LL]/rho_x_LL);
                            Y_x_LLL.push_back(Q[si][idx_cell_rho_Y_x_LLL]/rho_x_LLL);
                        }
                        
                        /*
                         * Get the pointers to the mass fractions.
                         */

                        std::vector<const Real*> Y_x_L_ptr;
                        std::vector<const Real*> Y_x_LL_ptr;
                        std::vector<const Real*> Y_x_LLL_ptr;
                        Y_x_L_ptr.reserve(d_num_species);
                        Y_x_LL_ptr.reserve(d_num_species);
                        Y_x_LLL_ptr.reserve(d_num_species);
                        for (int si = 0; si < d_num_species; si++)
                        {
                            Y_x_L_ptr.push_back(&Y_x_L[si]);
                            Y_x_LL_ptr.push_back(&Y_x_LL[si]);
                            Y_x_LLL_ptr.push_back(&Y_x_LLL[si]);
                        }
                       
                        const Real u_x_L   = Q[d_num_species][idx_cell_mom_x_L]/rho_x_L;
                        const Real u_x_LL  = Q[d_num_species][idx_cell_mom_x_LL]/rho_x_LL;
                        const Real u_x_LLL = Q[d_num_species][idx_cell_mom_x_LLL]/rho_x_LLL;
                        
                        const Real v_x_L   = Q[d_num_species + 1][idx_cell_mom_x_L]/rho_x_L;
                        const Real v_x_LL  = Q[d_num_species + 1][idx_cell_mom_x_LL]/rho_x_LL;
                        const Real v_x_LLL = Q[d_num_species + 1][idx_cell_mom_x_LLL]/rho_x_LLL;
                        
                        const Real half = Real(1)/Real(2);
                        const Real epsilon_x_L   = Q[d_num_species + 2][idx_cell_E_x_L]/rho_x_L - half*(u_x_L*u_x_L + v_x_L*v_x_L);
                        const Real epsilon_x_LL  = Q[d_num_species + 2][idx_cell_E_x_LL]/rho_x_LL - half*(u_x_LL*u_x_LL + v_x_LL*v_x_LL);
                        const Real epsilon_x_LLL = Q[d_num_species + 2][idx_cell_E_x_LLL]/rho_x_LLL - half*(u_x_LLL*u_x_LLL + v_x_LLL*v_x_LLL);
                        
                        Real p_x_L = d_equation_of_state_mixing_rules->
                            getPressure(
                                &rho_x_L,
                                &epsilon_x_L,
                                Y_x_L_ptr);
                        
                        Real p_x_LL = d_equation_of_state_mixing_rules->
                            getPressure(
                                &rho_x_LL,
                                &epsilon_x_LL,
                                Y_x_LL_ptr);
                        
                        Real p_x_LLL = d_equation_of_state_mixing_rules->
                            getPressure(
                                &rho_x_LLL,
                                &epsilon_x_LLL,
                                Y_x_LLL_ptr);
                        
                        /*
                         * Compute derivatives in x-direction.
                         */
                        
                        std::vector<Real> drho_Y_dx;
                        drho_Y_dx.reserve(d_num_species);
                        for (int si = 0; si < d_num_species; si++)
                        {
                            drho_Y_dx.push_back((Q[si][idx_cell_rho_Y_x_LLL] - Real(4)*Q[si][idx_cell_rho_Y_x_LL] +
                                Real(3)*Q[si][idx_cell_rho_Y_x_L])/(Real(2)*Real(dx[0])));
                        }
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
                            
                            const int idx_cell_rho_Y_y_T = (interior_box_hi_idx[0] + num_subghosts_conservative_var[0][0]) +
                                (j + 1 + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_mom_y_T = (interior_box_hi_idx[0] + num_subghosts_conservative_var[1][0]) +
                                (j + 1 + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0];
                            
                            const int idx_cell_E_y_T = (interior_box_hi_idx[0] + num_subghosts_conservative_var[2][0]) +
                                (j + 1 + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0];
                            
                            /*
                             * Compute the mixture density.
                             */
                            
                            Real rho_y_T = Real(0);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                rho_y_T += Q[si][idx_cell_rho_Y_y_T];
                            }
                            
                            /*
                             * Compute the mass fractions.
                             */
                            
                            std::vector<Real> Y_y_T;
                            Y_y_T.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_y_T.push_back(Q[si][idx_cell_rho_Y_y_T]/rho_y_T);
                            }
                            
                            /*
                             * Get the pointers to the mass fractions.
                             */
                            
                            std::vector<const Real*> Y_y_T_ptr;
                            Y_y_T_ptr.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_y_T_ptr.push_back(&Y_y_T[si]);
                            }
                            const Real u_y_T = Q[d_num_species][idx_cell_mom_y_T]/rho_y_T;
                            const Real v_y_T = Q[d_num_species + 1][idx_cell_mom_y_T]/rho_y_T;
                            const Real epsilon_y_T = Q[d_num_species + 2][idx_cell_E_y_T]/rho_y_T - half*(u_y_T*u_y_T + v_y_T*v_y_T);
                            
                            Real p_y_T = d_equation_of_state_mixing_rules->
                                getPressure(
                                    &rho_y_T,
                                    &epsilon_y_T,
                                    Y_y_T_ptr);
                            
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
                            
                            const int idx_cell_rho_Y_y_B = (interior_box_hi_idx[0] + num_subghosts_conservative_var[0][0]) +
                                (j - 1 + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_mom_y_B = (interior_box_hi_idx[0] + num_subghosts_conservative_var[1][0]) +
                                (j - 1 + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0];
                            
                            const int idx_cell_E_y_B = (interior_box_hi_idx[0] + num_subghosts_conservative_var[2][0]) +
                                (j - 1 + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0];
                            
                            /*
                             * Compute the mixture density.
                             */
                            
                            Real rho_y_B = Real(0);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                rho_y_B += Q[si][idx_cell_rho_Y_y_B];
                            }
                            
                            /*
                             * Compute the mass fractions.
                             */
                            
                            std::vector<Real> Y_y_B;
                            Y_y_B.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_y_B.push_back(Q[si][idx_cell_rho_Y_y_B]/rho_y_B);
                            }
                            
                            /*
                             * Get the pointers to the mass fractions.
                             */
                            
                            std::vector<const Real*> Y_y_B_ptr;
                            Y_y_B_ptr.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_y_B_ptr.push_back(&Y_y_B[si]);
                            }
                            
                            const Real u_y_B = Q[d_num_species][idx_cell_mom_y_B]/rho_y_B;
                            const Real v_y_B = Q[d_num_species + 1][idx_cell_mom_y_B]/rho_y_B;
                            const Real epsilon_y_B = Q[d_num_species + 2][idx_cell_E_y_B]/rho_y_B - half*(u_y_B*u_y_B + v_y_B*v_y_B);
                            
                            Real p_y_B = d_equation_of_state_mixing_rules->
                                getPressure(
                                    &rho_y_B,
                                    &epsilon_y_B,
                                    Y_y_B_ptr);
                            
                            // One-sided derivatives.
                            du_dy = (u_x_L - u_y_B)/Real(dx[1]);
                            dv_dy = (v_x_L - v_y_B)/Real(dx[1]);
                            dp_dy = (p_x_L - p_y_B)/Real(dx[1]);
                        }
                        else
                        {
                            const int idx_cell_rho_Y_y_B = (interior_box_hi_idx[0] + num_subghosts_conservative_var[0][0]) +
                                (j - 1 + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_rho_Y_y_T = (interior_box_hi_idx[0] + num_subghosts_conservative_var[0][0]) +
                                (j + 1 + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_mom_y_B = (interior_box_hi_idx[0] + num_subghosts_conservative_var[1][0]) +
                                (j - 1 + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0];
                            
                            const int idx_cell_mom_y_T = (interior_box_hi_idx[0] + num_subghosts_conservative_var[1][0]) +
                                (j + 1 + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0];
                            
                            const int idx_cell_E_y_B = (interior_box_hi_idx[0] + num_subghosts_conservative_var[2][0]) +
                                (j - 1 + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0];
                            
                            const int idx_cell_E_y_T = (interior_box_hi_idx[0] + num_subghosts_conservative_var[2][0]) +
                                (j + 1 + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0];
                            
                            /*
                             * Compute the mixture density.
                             */
                            
                            Real rho_y_B = Real(0);
                            Real rho_y_T = Real(0);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                rho_y_B += Q[si][idx_cell_rho_Y_y_B];
                                rho_y_T += Q[si][idx_cell_rho_Y_y_T];
                            }
                            
                            /*
                             * Compute the mass fractions.
                             */
                            
                            std::vector<Real> Y_y_B;
                            std::vector<Real> Y_y_T;
                            Y_y_B.reserve(d_num_species);
                            Y_y_T.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_y_B.push_back(Q[si][idx_cell_rho_Y_y_B]/rho_y_B);
                                Y_y_T.push_back(Q[si][idx_cell_rho_Y_y_T]/rho_y_T);
                            }
                            
                            /*
                             * Get the pointers to the mass fractions.
                             */
                            
                            std::vector<const Real*> Y_y_B_ptr;
                            std::vector<const Real*> Y_y_T_ptr;
                            Y_y_B_ptr.reserve(d_num_species);
                            Y_y_T_ptr.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_y_B_ptr.push_back(&Y_y_B[si]);
                                Y_y_T_ptr.push_back(&Y_y_T[si]);
                            }
                            
                            const Real u_y_B = Q[d_num_species][idx_cell_mom_y_B]/rho_y_B;
                            const Real u_y_T = Q[d_num_species][idx_cell_mom_y_T]/rho_y_T;
                            
                            const Real v_y_B = Q[d_num_species + 1][idx_cell_mom_y_B]/rho_y_B;
                            const Real v_y_T = Q[d_num_species + 1][idx_cell_mom_y_T]/rho_y_T;
                            
                            const Real epsilon_y_B = Q[d_num_species + 2][idx_cell_E_y_B]/rho_y_B - half*(u_y_B*u_y_B + v_y_B*v_y_B);
                            const Real epsilon_y_T = Q[d_num_species + 2][idx_cell_E_y_T]/rho_y_T - half*(u_y_T*u_y_T + v_y_T*v_y_T);
                            
                            Real p_y_B = d_equation_of_state_mixing_rules->
                                getPressure(
                                    &rho_y_B,
                                    &epsilon_y_B,
                                    Y_y_B_ptr);
                            
                            Real p_y_T = d_equation_of_state_mixing_rules->
                                getPressure(
                                    &rho_y_T,
                                    &epsilon_y_T,
                                    Y_y_T_ptr);
                            
                            // Central derivatives.
                            du_dy = (u_y_T - u_y_B)/(Real(2)*Real(dx[1]));
                            dv_dy = (v_y_T - v_y_B)/(Real(2)*Real(dx[1]));
                            dp_dy = (p_y_T - p_y_B)/(Real(2)*Real(dx[1]));
                        }
                        
                        // Compute sound speed.
                        
                        const Real Gamma_x_L = d_equation_of_state_mixing_rules->getGruneisenParameter(
                            &rho_x_L,
                            &p_x_L,
                            Y_x_L_ptr);
                        
                        const std::vector<Real> Psi_x_L = d_equation_of_state_mixing_rules->
                            getPressureDerivativeWithPartialDensities(
                                    &rho_x_L,
                                    &p_x_L,
                                    Y_x_L_ptr);
                        
                        Real c_x_L = Gamma_x_L*p_x_L/rho_x_L;
                        for (int si = 0; si < d_num_species; si++)
                        {
                            c_x_L += Y_x_L[si]*Psi_x_L[si];
                        }
                        c_x_L = std::sqrt(c_x_L);
                        
                        const Real lambda_1 = u_x_L - c_x_L;
                        
                        // Compute vector Lambda^(-1) * L.
                        
                        Real Lambda_inv_L[d_num_species + 3];
                        
                        const Real& p_t         = d_bdry_edge_nonreflecting_outflow_p_t[edge_loc];
                        const Real& sigma       = d_bdry_edge_nonreflecting_outflow_sigma[edge_loc];
                        const Real& beta        = d_bdry_edge_nonreflecting_outflow_beta[edge_loc];
                        const Real& length_char = d_bdry_edge_nonreflecting_outflow_length_char[edge_loc];
                        
                        const Real T_1 = v_x_L*(dp_dy - rho_x_L*c_x_L*du_dy) + rho_x_L*c_x_L*c_x_L*dv_dy;
                        
                        const Real M_sq = (u_x_L*u_x_L + v_x_L*v_x_L)/(c_x_L*c_x_L);
                        const Real K = sigma*c_x_L*(Real(1) - M_sq)/length_char;
                        
                        Lambda_inv_L[0] = (Real(1)/lambda_1)*(K*(p_x_L - p_t) - (Real(1) - beta)*T_1);
                        for (int si = 0; si < d_num_species; si++)
                        {
                            Lambda_inv_L[si + 1] = c_x_L*c_x_L*drho_Y_dx[si] - Y_x_L[si]*dp_dx;
                        }
                        Lambda_inv_L[d_num_species + 1] = dv_dx;
                        Lambda_inv_L[d_num_species + 2] = dp_dx + rho_x_L*c_x_L*du_dx;
                        
                        // Compute dV_dx.
                        
                        const Real c_sq_inv  = Real(1)/(c_x_L*c_x_L);
                        const Real rho_c_inv = Real(1)/(rho_x_L*c_x_L);
                        
                        Real dV_dx[d_num_species + 3];
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            dV_dx[si] = half*c_sq_inv*Y_x_L[si]*(Lambda_inv_L[0] + Lambda_inv_L[d_num_species + 2]) +
                                c_sq_inv*Lambda_inv_L[si + 1];
                        }
                        dV_dx[d_num_species]     = half*rho_c_inv*(-Lambda_inv_L[0] + Lambda_inv_L[d_num_species + 2]);
                        dV_dx[d_num_species + 1] = Lambda_inv_L[d_num_species + 1];
                        dV_dx[d_num_species + 2] = half*(Lambda_inv_L[0] + Lambda_inv_L[d_num_species + 2]);
                        
                        Real V_ghost[(d_num_species + 3)*num_ghosts_to_fill];
                        
                        for (int i = 0; i < num_ghosts_to_fill; i++)
                        {
                            const int idx_cell_rho_Y = (i + fill_box_lo_idx[0] + num_subghosts_conservative_var[0][0]) +
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
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    V_ghost[i*(d_num_species + 3) + si] = rho_Y_x_LL[si] + Real(2)*Real(dx[0])*dV_dx[si];
                                }
                                
                                V_ghost[i*(d_num_species + 3) + d_num_species]     = u_x_LL + Real(2)*Real(dx[0])*dV_dx[d_num_species];
                                V_ghost[i*(d_num_species + 3) + d_num_species + 1] = v_x_LL + Real(2)*Real(dx[0])*dV_dx[d_num_species + 1];
                                V_ghost[i*(d_num_species + 3) + d_num_species + 2] = p_x_LL + Real(2)*Real(dx[0])*dV_dx[d_num_species + 2];
                            }
                            else if (i == 1)
                            {
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    V_ghost[i*(d_num_species + 3) + si] = -Real(2)*rho_Y_x_LL[si] - Real(3)*rho_Y_x_L[si] +
                                        Real(6)*V_ghost[(i - 1)*(d_num_species + 3) + si] - Real(6)*Real(dx[0])*dV_dx[si];
                                }
                                
                                V_ghost[i*(d_num_species + 3) + d_num_species] = -Real(2)*u_x_LL - Real(3)*u_x_L +
                                    Real(6)*V_ghost[(i - 1)*(d_num_species + 3) + d_num_species] -
                                    Real(6)*Real(dx[0])*dV_dx[d_num_species];
                                
                                V_ghost[i*(d_num_species + 3) + d_num_species + 1] = -Real(2)*v_x_LL - Real(3)*v_x_L +
                                    Real(6)*V_ghost[(i - 1)*(d_num_species + 3) + d_num_species + 1] -
                                    Real(6)*Real(dx[0])*dV_dx[d_num_species + 1];
                    
                                V_ghost[i*(d_num_species + 3) + d_num_species + 2] = -Real(2)*p_x_LL - Real(3)*p_x_L +
                                    Real(6)*V_ghost[(i - 1)*(d_num_species + 3) + d_num_species + 2] -
                                    Real(6)*Real(dx[0])*dV_dx[d_num_species + 2];
                            }
                            else if (i == 2)
                            {
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    V_ghost[i*(d_num_species + 3) + si] = Real(3)*rho_Y_x_LL[si] + Real(10)*rho_Y_x_L[si] -
                                        Real(18)*V_ghost[(i - 2)*(d_num_species + 3) + si] +
                                        Real(6)*V_ghost[(i - 1)*(d_num_species + 3) + si] +
                                        Real(12)*Real(dx[0])*dV_dx[si];
                                }
                                
                                V_ghost[i*(d_num_species + 3) + d_num_species] = Real(3)*u_x_LL + Real(10)*u_x_L -
                                    Real(18)*V_ghost[(i - 2)*(d_num_species + 3) + d_num_species] +
                                    Real(6)*V_ghost[(i - 1)*(d_num_species + 3) + d_num_species] +
                                    Real(12)*Real(dx[0])*dV_dx[d_num_species];
                                
                                V_ghost[i*(d_num_species + 3) + d_num_species + 1] = Real(3)*v_x_LL + Real(10)*v_x_L -
                                    Real(18)*V_ghost[(i - 2)*(d_num_species + 3) + d_num_species + 1] +
                                    Real(6)*V_ghost[(i - 1)*(d_num_species + 3) + d_num_species + 1] +
                                    Real(12)*Real(dx[0])*dV_dx[d_num_species + 1];
                                
                                V_ghost[i*(d_num_species + 3) + d_num_species + 2] = Real(3)*p_x_LL + Real(10)*p_x_L -
                                    Real(18)*V_ghost[(i - 2)*(d_num_species + 3) + d_num_species + 2] +
                                    Real(6)*V_ghost[(i - 1)*(d_num_species + 3) + d_num_species + 2] +
                                    Real(12)*Real(dx[0])*dV_dx[d_num_species + 2];
                            }
                            else if (i == 3)
                            {
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    V_ghost[i*(d_num_species + 3) + si] = -Real(4)*rho_Y_x_LL[si] -
                                        Real(65)/Real(3)*rho_Y_x_L[si] +
                                        Real(40)*V_ghost[(i - 3)*(d_num_species + 3) + si] -
                                        Real(20)*V_ghost[(i - 2)*(d_num_species + 3) + si] +
                                        Real(20)/Real(3)*V_ghost[(i - 1)*(d_num_species + 3) + si] -
                                        Real(20)*Real(dx[0])*dV_dx[si];
                                }
                                
                                V_ghost[i*(d_num_species + 3) + d_num_species] = -Real(4)*u_x_LL -
                                    Real(65)/Real(3)*u_x_L +
                                    Real(40)*V_ghost[(i - 3)*(d_num_species + 3) + d_num_species] -
                                    Real(20)*V_ghost[(i - 2)*(d_num_species + 3) + d_num_species] +
                                    Real(20)/Real(3)*V_ghost[(i - 1)*(d_num_species + 3) + d_num_species] -
                                    Real(20)*Real(dx[0])*dV_dx[d_num_species];
                                
                                V_ghost[i*(d_num_species + 3) + d_num_species + 1] = -Real(4)*v_x_LL -
                                    Real(65)/Real(3)*v_x_L +
                                    Real(40)*V_ghost[(i - 3)*(d_num_species + 3) + d_num_species + 1] -
                                    Real(20)*V_ghost[(i - 2)*(d_num_species + 3) + d_num_species + 1] +
                                    Real(20)/Real(3)*V_ghost[(i - 1)*(d_num_species + 3) + d_num_species + 1] -
                                    Real(20)*Real(dx[0])*dV_dx[d_num_species + 1];
                                
                                V_ghost[i*(d_num_species + 3) + d_num_species + 2] = -Real(4)*p_x_LL -
                                    Real(65)/Real(3)*p_x_L +
                                    Real(40)*V_ghost[(i - 3)*(d_num_species + 3) + d_num_species + 2] -
                                    Real(20)*V_ghost[(i - 2)*(d_num_species + 3) + d_num_species + 2] +
                                    Real(20)/Real(3)*V_ghost[(i - 1)*(d_num_species + 3) + d_num_species + 2] -
                                    Real(20)*Real(dx[0])*dV_dx[d_num_species + 2];
                            }
                            else if (i == 4)
                            {
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    V_ghost[i*(d_num_species + 3) + si] = Real(5)*rho_Y_x_LL[si] +
                                        Real(77)/Real(2)*rho_Y_x_L[si] -
                                        Real(75)*V_ghost[(i - 4)*(d_num_species + 3) + si] +
                                        Real(50)*V_ghost[(i - 3)*(d_num_species + 3) + si] -
                                        Real(25)*V_ghost[(i - 2)*(d_num_species + 3) + si] +
                                        Real(15)/Real(2)*V_ghost[(i - 1)*(d_num_species + 3) + si] +
                                        Real(30)*Real(dx[0])*dV_dx[si];
                                }
                                
                                V_ghost[i*(d_num_species + 3) + d_num_species] = Real(5)*u_x_LL +
                                    Real(77)/Real(2)*u_x_L -
                                    Real(75)*V_ghost[(i - 4)*(d_num_species + 3) + d_num_species] +
                                    Real(50)*V_ghost[(i - 3)*(d_num_species + 3) + d_num_species] -
                                    Real(25)*V_ghost[(i - 2)*(d_num_species + 3) + d_num_species] +
                                    Real(15)/Real(2)*V_ghost[(i - 1)*(d_num_species + 3) + d_num_species] +
                                    Real(30)*Real(dx[0])*dV_dx[d_num_species];
                                
                                V_ghost[i*(d_num_species + 3) + d_num_species + 1] = Real(5)*v_x_LL +
                                    Real(77)/Real(2)*v_x_L -
                                    Real(75)*V_ghost[(i - 4)*(d_num_species + 3) + d_num_species + 1] +
                                    Real(50)*V_ghost[(i - 3)*(d_num_species + 3) + d_num_species + 1] -
                                    Real(25)*V_ghost[(i - 2)*(d_num_species + 3) + d_num_species + 1] +
                                    Real(15)/Real(2)*V_ghost[(i - 1)*(d_num_species + 3) + d_num_species + 1] +
                                    Real(30)*Real(dx[0])*dV_dx[d_num_species + 1];
                                
                                V_ghost[i*(d_num_species + 3) + d_num_species + 2] = Real(5)*p_x_LL +
                                    Real(77)/Real(2)*p_x_L -
                                    Real(75)*V_ghost[(i - 4)*(d_num_species + 3) + d_num_species + 2] +
                                    Real(50)*V_ghost[(i - 3)*(d_num_species + 3) + d_num_species + 2] -
                                    Real(25)*V_ghost[(i - 2)*(d_num_species + 3) + d_num_species + 2] +
                                    Real(15)/Real(2)*V_ghost[(i - 1)*(d_num_species + 3) + d_num_species + 2] +
                                    Real(30)*Real(dx[0])*dV_dx[d_num_species + 2];
                            }
                            else if (i == 5)
                            {
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    V_ghost[i*(d_num_species + 3) + si] = -Real(6)*rho_Y_x_LL[si] -
                                        Real(609)/Real(10)*rho_Y_x_L[si] +
                                        Real(126)*V_ghost[(i - 5)*(d_num_species + 3) + si] -
                                        Real(105)*V_ghost[(i - 4)*(d_num_species + 3) + si] +
                                        Real(70)*V_ghost[(i - 3)*(d_num_species + 3) + si] -
                                        Real(63)/Real(2)*V_ghost[(i - 2)*(d_num_species + 3) + si] +
                                        Real(42)/Real(5)*V_ghost[(i - 1)*(d_num_species + 3) + si] -
                                        Real(42)*Real(dx[0])*dV_dx[si];
                                }
                                
                                V_ghost[i*(d_num_species + 3) + d_num_species] = -Real(6)*u_x_LL -
                                    Real(609)/Real(10)*u_x_L +
                                    Real(126)*V_ghost[(i - 5)*(d_num_species + 3) + d_num_species] -
                                    Real(105)*V_ghost[(i - 4)*(d_num_species + 3) + d_num_species] +
                                    Real(70)*V_ghost[(i - 3)*(d_num_species + 3) + d_num_species] -
                                    Real(63)/Real(2)*V_ghost[(i - 2)*(d_num_species + 3) + d_num_species] +
                                    Real(42)/Real(5)*V_ghost[(i - 1)*(d_num_species + 3) + d_num_species] -
                                    Real(42)*Real(dx[0])*dV_dx[d_num_species];
                                
                                V_ghost[i*(d_num_species + 3) + d_num_species + 1] = -Real(6)*v_x_LL -
                                    Real(609)/Real(10)*v_x_L +
                                    Real(126)*V_ghost[(i - 5)*(d_num_species + 3) + d_num_species + 1] -
                                    Real(105)*V_ghost[(i - 4)*(d_num_species + 3) + d_num_species + 1] +
                                    Real(70)*V_ghost[(i - 3)*(d_num_species + 3) + d_num_species + 1] -
                                    Real(63)/Real(2)*V_ghost[(i - 2)*(d_num_species + 3) + d_num_species + 1] +
                                    Real(42)/Real(5)*V_ghost[(i - 1)*(d_num_species + 3) + d_num_species + 1] -
                                    Real(42)*Real(dx[0])*dV_dx[d_num_species + 1];
                                
                                V_ghost[i*(d_num_species + 3) + d_num_species + 2] = -Real(6)*p_x_LL -
                                    Real(609)/Real(10)*p_x_L +
                                    Real(126)*V_ghost[(i - 5)*(d_num_species + 3) + d_num_species + 2] -
                                    Real(105)*V_ghost[(i - 4)*(d_num_species + 3) + d_num_species + 2] +
                                    Real(70)*V_ghost[(i - 3)*(d_num_species + 3) + d_num_species + 2] -
                                    Real(63)/Real(2)*V_ghost[(i - 2)*(d_num_species + 3) + d_num_species + 2] +
                                    Real(42)/Real(5)*V_ghost[(i - 1)*(d_num_species + 3) + d_num_species + 2] -
                                    Real(42)*Real(dx[0])*dV_dx[d_num_species + 2];
                            }
                            
                            /*
                             * Compute the mixture density.
                             */
                            
                            Real rho_ghost = Real(0);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                rho_ghost += V_ghost[i*(d_num_species + 3) + si];
                            }
                            
                            /*
                             * Compute the mass fractions.
                             */
                            
                            std::vector<Real> Y_ghost;
                            Y_ghost.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_ghost.push_back(V_ghost[i*(d_num_species + 3) + si]/rho_ghost);
                            }
                            
                            /*
                             * Get the pointers to the mass fractions.
                             */
                            
                            std::vector<const Real*> Y_ghost_ptr;
                            Y_ghost_ptr.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_ghost_ptr.push_back(&Y_ghost[si]);
                            }
                            
                            for(int si=0; si < d_num_species; si++)
                            {
                                Q[si][idx_cell_rho_Y] = V_ghost[i*(d_num_species + 3) + si];
                            }
                            
                            Q[d_num_species][idx_cell_mom]     = rho_ghost*V_ghost[i*(d_num_species + 3) + d_num_species];
                            Q[d_num_species + 1][idx_cell_mom] = rho_ghost*V_ghost[i*(d_num_species + 3) + d_num_species + 1];
                            
                            const Real epsilon = d_equation_of_state_mixing_rules->
                                getInternalEnergy(
                                    &rho_ghost,
                                    &V_ghost[i*(d_num_species + 3) + d_num_species + 2],
                                    Y_ghost_ptr);
                            
                            const Real E = rho_ghost*epsilon +
                                half*(Q[d_num_species][idx_cell_mom]*Q[d_num_species][idx_cell_mom] +
                                    Q[d_num_species + 1][idx_cell_mom]*Q[d_num_species + 1][idx_cell_mom])/
                                    rho_ghost;
                            
                            Q[d_num_species + 2][idx_cell_E] = E;
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
                            << ": FlowModelBoundaryUtilitiesFourEqnConservative::fill2dEdgeBoundaryData()\n"
                            << "Non-reflecting outflow BC doesn't support more than six ghost cells yet!");
                    }
                    
                    for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                    {
                        // Get the grid spacing.
                        const double* const dx = patch_geom->getDx();
                        
                        const int idx_cell_rho_Y_y_T = (i + num_subghosts_conservative_var[0][0]) +
                            (interior_box_lo_idx[1] + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0];
                        
                        const int idx_cell_rho_Y_y_TT = (i + num_subghosts_conservative_var[0][0]) +
                            (interior_box_lo_idx[1] + 1 + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0];
                        
                        const int idx_cell_rho_Y_y_TTT = (i + num_subghosts_conservative_var[0][0])+
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
                        
                        std::vector<Real> rho_Y_y_T;
                        std::vector<Real> rho_Y_y_TT;
                        std::vector<Real> rho_Y_y_TTT;
                        rho_Y_y_T.reserve(d_num_species);
                        rho_Y_y_TT.reserve(d_num_species);
                        rho_Y_y_TTT.reserve(d_num_species);
                        for (int si = 0; si < d_num_species; si++)
                        {
                            rho_Y_y_T.push_back(Q[si][idx_cell_rho_Y_y_T]);
                            rho_Y_y_TT.push_back(Q[si][idx_cell_rho_Y_y_TT]);
                            rho_Y_y_TTT.push_back(Q[si][idx_cell_rho_Y_y_TTT]);
                        }
                        
                        /*
                         * Compute the mixture density.
                         */
                        
                        Real rho_y_T   = Real(0);
                        Real rho_y_TT  = Real(0);
                        Real rho_y_TTT = Real(0);
                        for (int si = 0; si < d_num_species; si++)
                        {
                            rho_y_T   += Q[si][idx_cell_rho_Y_y_T];
                            rho_y_TT  += Q[si][idx_cell_rho_Y_y_TT];
                            rho_y_TTT += Q[si][idx_cell_rho_Y_y_TTT];
                        }
                        
                        /*
                         * Compute the mass fractions.
                         */
                        
                        std::vector<Real> Y_y_T;
                        std::vector<Real> Y_y_TT;
                        std::vector<Real> Y_y_TTT;
                        Y_y_T.reserve(d_num_species);
                        Y_y_TT.reserve(d_num_species);
                        Y_y_TTT.reserve(d_num_species);
                        for (int si = 0; si < d_num_species; si++)
                        {
                            Y_y_T.push_back(Q[si][idx_cell_rho_Y_y_T]/rho_y_T);
                            Y_y_TT.push_back(Q[si][idx_cell_rho_Y_y_TT]/rho_y_TT);
                            Y_y_TTT.push_back(Q[si][idx_cell_rho_Y_y_TTT]/rho_y_TTT);
                        }
                        
                        /*
                         * Get the pointers to the mass fractions.
                         */
                        
                        std::vector<const Real*> Y_y_T_ptr;
                        std::vector<const Real*> Y_y_TT_ptr;
                        std::vector<const Real*> Y_y_TTT_ptr;
                        Y_y_T_ptr.reserve(d_num_species);
                        Y_y_TT_ptr.reserve(d_num_species);
                        Y_y_TTT_ptr.reserve(d_num_species);
                        for (int si = 0; si < d_num_species; si++)
                        {
                            Y_y_T_ptr.push_back(&Y_y_T[si]);
                            Y_y_TT_ptr.push_back(&Y_y_TT[si]);
                            Y_y_TTT_ptr.push_back(&Y_y_TTT[si]);
                        }
                        
                        const Real u_y_T   = Q[d_num_species][idx_cell_mom_y_T]/rho_y_T;
                        const Real u_y_TT  = Q[d_num_species][idx_cell_mom_y_TT]/rho_y_TT;
                        const Real u_y_TTT = Q[d_num_species][idx_cell_mom_y_TTT]/rho_y_TTT;
                        
                        const Real v_y_T   = Q[d_num_species + 1][idx_cell_mom_y_T]/rho_y_T;
                        const Real v_y_TT  = Q[d_num_species + 1][idx_cell_mom_y_TT]/rho_y_TT;
                        const Real v_y_TTT = Q[d_num_species + 1][idx_cell_mom_y_TTT]/rho_y_TTT;
                        
                        const Real half = Real(1)/Real(2);
                        const Real epsilon_y_T   = Q[d_num_species + 2][idx_cell_E_y_T]/rho_y_T - half*(u_y_T*u_y_T + v_y_T*v_y_T);
                        const Real epsilon_y_TT  = Q[d_num_species + 2][idx_cell_E_y_TT]/rho_y_TT - half*(u_y_TT*u_y_TT + v_y_TT*v_y_TT);
                        const Real epsilon_y_TTT = Q[d_num_species + 2][idx_cell_E_y_TTT]/rho_y_TTT- half*(u_y_TTT*u_y_TTT + v_y_TTT*v_y_TTT);
                        
                        Real p_y_T = d_equation_of_state_mixing_rules->
                            getPressure(
                                &rho_y_T,
                                &epsilon_y_T,
                                Y_y_T_ptr);
                        
                        Real p_y_TT = d_equation_of_state_mixing_rules->
                            getPressure(
                                &rho_y_TT,
                                &epsilon_y_TT,
                                Y_y_TT_ptr);
                        
                        Real p_y_TTT = d_equation_of_state_mixing_rules->
                            getPressure(
                                &rho_y_TTT,
                                &epsilon_y_TTT,
                                Y_y_TTT_ptr);
                        
                        /*
                         * Compute derivatives in y-direction.
                         */
                        
                        std::vector<Real> drho_Y_dy;
                        drho_Y_dy.reserve(d_num_species);
                        for (int si = 0; si < d_num_species; si++)
                        {
                            drho_Y_dy.push_back(-(Q[si][idx_cell_rho_Y_y_TTT] - Real(4)*Q[si][idx_cell_rho_Y_y_TT] +
                                Real(3)*Q[si][idx_cell_rho_Y_y_T])/(Real(2)*Real(dx[1])));
                        }
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
                            
                            const int idx_cell_rho_Y_x_R = (i + 1 + num_subghosts_conservative_var[0][0]) +
                                (interior_box_lo_idx[1] + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_mom_x_R = (i + 1 + num_subghosts_conservative_var[1][0]) +
                                (interior_box_lo_idx[1] + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0];
                            
                            const int idx_cell_E_x_R = (i + 1 + num_subghosts_conservative_var[2][0]) +
                                (interior_box_lo_idx[1] + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0];
                            
                            /*
                             * Compute the mixture density.
                             */
                            
                            Real rho_x_R = Real(0);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                rho_x_R += Q[si][idx_cell_rho_Y_x_R];
                            }
                            
                            /*
                             * Compute the mass fractions.
                             */
                            
                            std::vector<Real> Y_x_R;
                            Y_x_R.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_x_R.push_back(Q[si][idx_cell_rho_Y_x_R]/rho_x_R);
                            }
                            
                            /*
                             * Get the pointers to the mass fractions.
                             */
                            
                            std::vector<const Real*> Y_x_R_ptr;
                            Y_x_R_ptr.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_x_R_ptr.push_back(&Y_x_R[si]);
                            }
                            
                            const Real u_x_R = Q[d_num_species][idx_cell_mom_x_R]/rho_x_R;
                            const Real v_x_R = Q[d_num_species + 1][idx_cell_mom_x_R]/rho_x_R;
                            const Real epsilon_x_R = Q[d_num_species + 2][idx_cell_E_x_R]/rho_x_R - half*(u_x_R*u_x_R + v_x_R*v_x_R);
                            
                            Real p_x_R = d_equation_of_state_mixing_rules->
                                getPressure(
                                    &rho_x_R,
                                    &epsilon_x_R,
                                    Y_x_R_ptr);
                            
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
                            
                            const int idx_cell_rho_Y_x_L = (i - 1 + num_subghosts_conservative_var[0][0]) +
                                (interior_box_lo_idx[1] + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_mom_x_L = (i - 1 + num_subghosts_conservative_var[1][0]) +
                                (interior_box_lo_idx[1] + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0];
                            
                            const int idx_cell_E_x_L = (i - 1 + num_subghosts_conservative_var[2][0]) +
                                (interior_box_lo_idx[1] + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0];
                            
                            /*
                             * Compute the mixture density.
                             */
                            
                            Real rho_x_L = Real(0);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                rho_x_L += Q[si][idx_cell_rho_Y_x_L];
                            }
                            
                            /*
                             * Compute the mass fractions.
                             */
                            
                            std::vector<Real> Y_x_L;
                            Y_x_L.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_x_L.push_back(Q[si][idx_cell_rho_Y_x_L]/rho_x_L);
                            }
                            
                            /*
                             * Get the pointers to the mass fractions.
                             */
                            
                            std::vector<const Real*> Y_x_L_ptr;
                            Y_x_L_ptr.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_x_L_ptr.push_back(&Y_x_L[si]);
                            }
                            
                            const Real u_x_L = Q[d_num_species][idx_cell_mom_x_L]/rho_x_L;
                            const Real v_x_L = Q[d_num_species + 1][idx_cell_mom_x_L]/rho_x_L;
                            const Real epsilon_x_L = Q[d_num_species + 2][idx_cell_E_x_L]/rho_x_L - half*(u_x_L*u_x_L + v_x_L*v_x_L);

                            Real p_x_L = d_equation_of_state_mixing_rules->
                                getPressure(
                                    &rho_x_L,
                                    &epsilon_x_L,
                                    Y_x_L_ptr);
                            
                            // One-sided derivatives.
                            du_dx = (u_y_T - u_x_L)/Real(dx[0]);
                            dv_dx = (v_y_T - v_x_L)/Real(dx[0]);
                            dp_dx = (p_y_T - p_x_L)/Real(dx[0]);
                        }
                        else
                        {
                            const int idx_cell_rho_Y_x_L = (i - 1 + num_subghosts_conservative_var[0][0]) +
                                (interior_box_lo_idx[1] + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_rho_Y_x_R = (i + 1 + num_subghosts_conservative_var[0][0]) +
                                (interior_box_lo_idx[1] + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_mom_x_L = (i - 1 + num_subghosts_conservative_var[1][0]) +
                                (interior_box_lo_idx[1] + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0];
                            
                            const int idx_cell_mom_x_R = (i + 1 + num_subghosts_conservative_var[1][0]) +
                                (interior_box_lo_idx[1] + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0];
                            
                            const int idx_cell_E_x_L = (i - 1 + num_subghosts_conservative_var[2][0]) +
                                (interior_box_lo_idx[1] + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0];
                            
                            const int idx_cell_E_x_R = (i + 1 + num_subghosts_conservative_var[2][0]) +
                                (interior_box_lo_idx[1] + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0];
                            
                            /*
                             * Compute the mixture density.
                             */
                            
                            Real rho_x_L = Real(0);
                            Real rho_x_R = Real(0);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                rho_x_L += Q[si][idx_cell_rho_Y_x_L];
                                rho_x_R += Q[si][idx_cell_rho_Y_x_R];
                            }
                            
                            /*
                             * Compute the mass fractions.
                             */
                            
                            std::vector<Real> Y_x_L;
                            std::vector<Real> Y_x_R;
                            Y_x_L.reserve(d_num_species);
                            Y_x_R.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_x_L.push_back(Q[si][idx_cell_rho_Y_x_L]/rho_x_L);
                                Y_x_R.push_back(Q[si][idx_cell_rho_Y_x_R]/rho_x_R);
                            }
                            
                            /*
                             * Get the pointers to the mass fractions.
                             */
                            
                            std::vector<const Real*> Y_x_L_ptr;
                            std::vector<const Real*> Y_x_R_ptr;
                            Y_x_L_ptr.reserve(d_num_species);
                            Y_x_R_ptr.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_x_L_ptr.push_back(&Y_x_L[si]);
                                Y_x_R_ptr.push_back(&Y_x_R[si]);
                            }
                            
                            const Real u_x_L = Q[d_num_species][idx_cell_mom_x_L]/rho_x_L;
                            const Real u_x_R = Q[d_num_species][idx_cell_mom_x_R]/rho_x_R;
                            
                            const Real v_x_L = Q[d_num_species + 1][idx_cell_mom_x_L]/rho_x_L;
                            const Real v_x_R = Q[d_num_species + 1][idx_cell_mom_x_R]/rho_x_R;
                            
                            const Real epsilon_x_L = Q[d_num_species + 2][idx_cell_E_x_L]/rho_x_L - half*(u_x_L*u_x_L + v_x_L*v_x_L);
                            const Real epsilon_x_R = Q[d_num_species + 2][idx_cell_E_x_R]/rho_x_R - half*(u_x_R*u_x_R + v_x_R*v_x_R);
                            
                            Real p_x_L = d_equation_of_state_mixing_rules->
                                getPressure(
                                    &rho_x_L,
                                    &epsilon_x_L,
                                    Y_x_L_ptr);
                            
                            Real p_x_R = d_equation_of_state_mixing_rules->
                                getPressure(
                                    &rho_x_R,
                                    &epsilon_x_R,
                                    Y_x_R_ptr);
                            
                            // Central derivatives.
                            du_dx = (u_x_R - u_x_L)/(Real(2)*Real(dx[0]));
                            dv_dx = (v_x_R - v_x_L)/(Real(2)*Real(dx[0]));
                            dp_dx = (p_x_R - p_x_L)/(Real(2)*Real(dx[0]));
                        }
                        
                        // Compute sound speed.
                        
                        const Real Gamma_y_T = d_equation_of_state_mixing_rules->getGruneisenParameter(
                            &rho_y_T,
                            &p_y_T,
                            Y_y_T_ptr);
                        
                        const std::vector<Real> Psi_y_T = d_equation_of_state_mixing_rules->
                            getPressureDerivativeWithPartialDensities(
                                    &rho_y_T,
                                    &p_y_T,
                                    Y_y_T_ptr);
                        
                        Real c_y_T = Gamma_y_T*p_y_T/rho_y_T;
                        for (int si = 0; si < d_num_species; si++)
                        {
                            c_y_T += Y_y_T[si]*Psi_y_T[si];
                        }
                        c_y_T = std::sqrt(c_y_T);
                        
                        const Real lambda_last = v_y_T + c_y_T;
                        
                        // Compute vector Lambda^(-1) * L.
                        
                        Real Lambda_inv_L[d_num_species + 3];
                        
                        const Real& p_t         = d_bdry_edge_nonreflecting_outflow_p_t[edge_loc];
                        const Real& sigma       = d_bdry_edge_nonreflecting_outflow_sigma[edge_loc];
                        const Real& beta        = d_bdry_edge_nonreflecting_outflow_beta[edge_loc];
                        const Real& length_char = d_bdry_edge_nonreflecting_outflow_length_char[edge_loc];
                        
                        const Real T_last = u_y_T*(dp_dx + rho_y_T*c_y_T*dv_dx) + rho_y_T*c_y_T*c_y_T*du_dx;
                        
                        const Real M_sq = (v_y_T*v_y_T + u_y_T*u_y_T)/(c_y_T*c_y_T);
                        const Real K = sigma*c_y_T*(Real(1) - M_sq)/length_char;
                        
                        Lambda_inv_L[0] = dp_dy - rho_y_T*c_y_T*dv_dy;
                        Lambda_inv_L[1] = du_dy;
                        for (int si = 0; si < d_num_species; si++)
                        {
                            Lambda_inv_L[si + 2] = c_y_T*c_y_T*drho_Y_dy[si] - Y_y_T[si]*dp_dy;
                        }
                        Lambda_inv_L[d_num_species + 2] = (Real(1)/lambda_last)*(K*(p_y_T - p_t) - (Real(1) - beta)*T_last);
                        
                        // Compute dV_dx.
                        
                        const Real c_sq_inv  = Real(1)/(c_y_T*c_y_T);
                        const Real rho_c_inv = Real(1)/(rho_y_T*c_y_T);
                        
                        Real dV_dy[d_num_species + 3];
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            dV_dy[si] = half*c_sq_inv*Y_y_T[si]*(Lambda_inv_L[0] + Lambda_inv_L[d_num_species + 2]) +
                                c_sq_inv*Lambda_inv_L[si + 2];
                        }
                        dV_dy[d_num_species]     = Lambda_inv_L[1];
                        dV_dy[d_num_species + 1] = half*rho_c_inv*(-Lambda_inv_L[0] + Lambda_inv_L[d_num_species + 2]);
                        dV_dy[d_num_species + 2] = half*(Lambda_inv_L[0] + Lambda_inv_L[d_num_species + 2]);
                        
                        Real V_ghost[(d_num_species + 3)*num_ghosts_to_fill];
                        
                        for (int j = num_ghosts_to_fill - 1; j >= 0; j--)
                        {
                            const int idx_cell_rho_Y = (i + num_subghosts_conservative_var[0][0]) +
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
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    V_ghost[j*(d_num_species + 3) + si] = rho_Y_y_TT[si] - Real(2)*Real(dx[1])*dV_dy[si];
                                }
                                
                                V_ghost[j*(d_num_species + 3) + d_num_species]     = u_y_TT   - Real(2)*Real(dx[1])*dV_dy[d_num_species];
                                V_ghost[j*(d_num_species + 3) + d_num_species + 1] = v_y_TT   - Real(2)*Real(dx[1])*dV_dy[d_num_species + 1];
                                V_ghost[j*(d_num_species + 3) + d_num_species + 2] = p_y_TT   - Real(2)*Real(dx[1])*dV_dy[d_num_species + 2];
                            }
                            else if (j == num_ghosts_to_fill - 2)
                            {
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    V_ghost[j*(d_num_species + 3) + si] = -Real(2)*rho_Y_y_TT[si] - Real(3)*rho_Y_y_T[si] +
                                        Real(6)*V_ghost[(j + 1)*(d_num_species + 3) + si] + Real(6)*Real(dx[1])*dV_dy[si];
                                }
                                
                                V_ghost[j*(d_num_species + 3) + d_num_species] = -Real(2)*u_y_TT - Real(3)*u_y_T +
                                    Real(6)*V_ghost[(j + 1)*(d_num_species + 3) + d_num_species] + Real(6)*Real(dx[1])*dV_dy[d_num_species];
                                
                                V_ghost[j*(d_num_species + 3) + d_num_species + 1] = -Real(2)*v_y_TT - Real(3)*v_y_T +
                                    Real(6)*V_ghost[(j + 1)*(d_num_species + 3) + d_num_species + 1] + Real(6)*Real(dx[1])*dV_dy[d_num_species + 1];
                                
                                V_ghost[j*(d_num_species + 3) +d_num_species + 2] = -Real(2)*p_y_TT - Real(3)*p_y_T +
                                    Real(6)*V_ghost[(j + 1)*(d_num_species + 3) + d_num_species + 2] + Real(6)*Real(dx[1])*dV_dy[d_num_species + 2];
                            }
                            else if (j == num_ghosts_to_fill - 3)
                            {
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    V_ghost[j*(d_num_species + 3) + si] = Real(3)*rho_Y_y_TT[si] + Real(10)*rho_Y_y_T[si] -
                                    Real(18)*V_ghost[(j + 2)*(d_num_species + 3) + si] +
                                    Real(6)*V_ghost[(j + 1)*(d_num_species + 3) + si] -
                                    Real(12)*Real(dx[1])*dV_dy[si];
                                }
                                
                                V_ghost[j*(d_num_species + 3) + d_num_species] = Real(3)*u_y_TT + Real(10)*u_y_T -
                                    Real(18)*V_ghost[(j + 2)*(d_num_species + 3) + d_num_species] +
                                    Real(6)*V_ghost[(j + 1)*(d_num_species + 3) + d_num_species] -
                                    Real(12)*Real(dx[1])*dV_dy[d_num_species];
                                
                                V_ghost[j*(d_num_species + 3) + d_num_species + 1] = Real(3)*v_y_TT + Real(10)*v_y_T -
                                    Real(18)*V_ghost[(j + 2)*(d_num_species + 3) + d_num_species + 1] +
                                    Real(6)*V_ghost[(j + 1)*(d_num_species + 3) + d_num_species + 1] -
                                    Real(12)*Real(dx[1])*dV_dy[d_num_species + 1];
                                
                                V_ghost[j*(d_num_species + 3) + d_num_species + 2] = Real(3)*p_y_TT + Real(10)*p_y_T -
                                    Real(18)*V_ghost[(j + 2)*(d_num_species + 3) + d_num_species + 2] +
                                    Real(6)*V_ghost[(j + 1)*(d_num_species + 3) + d_num_species + 2] -
                                    Real(12)*Real(dx[1])*dV_dy[d_num_species + 2];
                            }
                            else if (j == num_ghosts_to_fill - 4)
                            {
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    V_ghost[j*(d_num_species + 3) + si] = -Real(4)*rho_Y_y_TT[si] -
                                        Real(65)/Real(3)*rho_Y_y_T[si] +
                                        Real(40)*V_ghost[(j + 3)*(d_num_species + 3) + si] -
                                        Real(20)*V_ghost[(j + 2)*(d_num_species + 3) + si] +
                                        Real(20)/Real(3)*V_ghost[(j + 1)*(d_num_species + 3) + si] +
                                        Real(20)*Real(dx[1])*dV_dy[si];
                                } 
                                
                                V_ghost[j*(d_num_species + 3) + d_num_species] = -Real(4)*u_y_TT -
                                    Real(65)/Real(3)*u_y_T +
                                    Real(40)*V_ghost[(j + 3)*(d_num_species + 3) + d_num_species] -
                                    Real(20)*V_ghost[(j + 2)*(d_num_species + 3) + d_num_species] +
                                    Real(20)/Real(3)*V_ghost[(j + 1)*(d_num_species + 3) + d_num_species] +
                                    Real(20)*Real(dx[1])*dV_dy[d_num_species];
                                
                                V_ghost[j*(d_num_species + 3) + d_num_species + 1] = -Real(4)*v_y_TT -
                                    Real(65)/Real(3)*v_y_T +
                                    Real(40)*V_ghost[(j + 3)*(d_num_species + 3) + d_num_species + 1] -
                                    Real(20)*V_ghost[(j + 2)*(d_num_species + 3) + d_num_species + 1] +
                                    Real(20)/Real(3)*V_ghost[(j + 1)*(d_num_species + 3) + d_num_species + 1] +
                                    Real(20)*Real(dx[1])*dV_dy[d_num_species + 1];
                                
                                V_ghost[j*(d_num_species + 3) + d_num_species + 2] = -Real(4)*p_y_TT -
                                    Real(65)/Real(3)*p_y_T +
                                    Real(40)*V_ghost[(j + 3)*(d_num_species + 3) + d_num_species + 2] -
                                    Real(20)*V_ghost[(j + 2)*(d_num_species + 3) + d_num_species + 2] +
                                    Real(20)/Real(3)*V_ghost[(j + 1)*(d_num_species + 3) + d_num_species + 2] +
                                    Real(20)*Real(dx[1])*dV_dy[d_num_species + 2];
                            }
                            else if (j == num_ghosts_to_fill - 5)
                            {
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    V_ghost[j*(d_num_species + 3) + si] = Real(5)*rho_Y_y_TT[si] +
                                        Real(77)/Real(2)*rho_Y_y_T[si] -
                                        Real(75)*V_ghost[(j + 4)*(d_num_species + 3) + si] +
                                        Real(50)*V_ghost[(j + 3)*(d_num_species + 3) + si] -
                                        Real(25)*V_ghost[(j + 2)*(d_num_species + 3) + si] +
                                        Real(15)/Real(2)*V_ghost[(j + 1)*(d_num_species + 3) + si] -
                                        Real(30)*Real(dx[1])*dV_dy[si];
                                }
                                
                                V_ghost[j*(d_num_species + 3) + d_num_species] = Real(5)*u_y_TT +
                                    Real(77)/Real(2)*u_y_T -
                                    Real(75)*V_ghost[(j + 4)*(d_num_species + 3) + d_num_species] +
                                    Real(50)*V_ghost[(j + 3)*(d_num_species + 3) + d_num_species] -
                                    Real(25)*V_ghost[(j + 2)*(d_num_species + 3) + d_num_species] +
                                    Real(15)/Real(2)*V_ghost[(j + 1)*(d_num_species + 3) + d_num_species] -
                                    Real(30)*Real(dx[1])*dV_dy[d_num_species];
                                
                                V_ghost[j*(d_num_species + 3) + d_num_species + 1] = Real(5)*v_y_TT +
                                    Real(77)/Real(2)*v_y_T -
                                    Real(75)*V_ghost[(j + 4)*(d_num_species + 3) + d_num_species + 1] +
                                    Real(50)*V_ghost[(j + 3)*(d_num_species + 3) + d_num_species + 1] -
                                    Real(25)*V_ghost[(j + 2)*(d_num_species + 3) + d_num_species + 1] +
                                    Real(15)/Real(2)*V_ghost[(j + 1)*(d_num_species + 3) + d_num_species + 1] -
                                    Real(30)*Real(dx[1])*dV_dy[d_num_species + 1];
                                
                                V_ghost[j*(d_num_species + 3) + d_num_species + 2] = Real(5)*p_y_TT +
                                    Real(77)/Real(2)*p_y_T -
                                    Real(75)*V_ghost[(j + 4)*(d_num_species + 3) + d_num_species + 2] +
                                    Real(50)*V_ghost[(j + 3)*(d_num_species + 3) + d_num_species + 2] -
                                    Real(25)*V_ghost[(j + 2)*(d_num_species + 3) + d_num_species + 2] +
                                    Real(15)/Real(2)*V_ghost[(j + 1)*(d_num_species + 3) + d_num_species + 2] -
                                    Real(30)*Real(dx[1])*dV_dy[d_num_species + 2];
                            }
                            else if (j == num_ghosts_to_fill - 6)
                            {
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    V_ghost[j*(d_num_species + 3) + si] = -Real(6)*rho_Y_y_TT[si] -
                                        Real(609)/Real(10)*rho_Y_y_T[si] +
                                        Real(126)*V_ghost[(j + 5)*(d_num_species + 3) + si] -
                                        Real(105)*V_ghost[(j + 4)*(d_num_species + 3) + si] +
                                        Real(70)*V_ghost[(j + 3)*(d_num_species + 3) + si] -
                                        Real(63)/Real(2)*V_ghost[(j + 2)*(d_num_species + 3) + si] +
                                        Real(42)/Real(5)*V_ghost[(j + 1)*(d_num_species + 3) + si] +
                                        Real(42)*Real(dx[1])*dV_dy[si];
                                }
                                
                                V_ghost[j*(d_num_species + 3) + d_num_species] = -Real(6)*u_y_TT -
                                    Real(609)/Real(10)*u_y_T +
                                    Real(126)*V_ghost[(j + 5)*(d_num_species + 3) + d_num_species] -
                                    Real(105)*V_ghost[(j + 4)*(d_num_species + 3) + d_num_species] +
                                    Real(70)*V_ghost[(j + 3)*(d_num_species + 3) + d_num_species] -
                                    Real(63)/Real(2)*V_ghost[(j + 2)*(d_num_species + 3) + d_num_species] +
                                    Real(42)/Real(5)*V_ghost[(j + 1)*(d_num_species + 3) + d_num_species] +
                                    Real(42)*Real(dx[1])*dV_dy[d_num_species];
                                
                                V_ghost[j*(d_num_species + 3) + d_num_species + 1] = -Real(6)*v_y_TT -
                                    Real(609)/Real(10)*v_y_T +
                                    Real(126)*V_ghost[(j + 5)*(d_num_species + 3) + d_num_species + 1] -
                                    Real(105)*V_ghost[(j + 4)*(d_num_species + 3) + d_num_species + 1] +
                                    Real(70)*V_ghost[(j + 3)*(d_num_species + 3) + d_num_species + 1] -
                                    Real(63)/Real(2)*V_ghost[(j + 2)*(d_num_species + 3) + d_num_species + 1] +
                                    Real(42)/Real(5)*V_ghost[(j + 1)*(d_num_species + 3) + d_num_species + 1] +
                                    Real(42)*Real(dx[1])*dV_dy[d_num_species + 1];
                                
                                V_ghost[j*(d_num_species + 3) + d_num_species + 2] = -Real(6)*p_y_TT -
                                    Real(609)/Real(10)*p_y_T +
                                    Real(126)*V_ghost[(j + 5)*(d_num_species + 3) + d_num_species + 2] -
                                    Real(105)*V_ghost[(j + 4)*(d_num_species + 3) + d_num_species + 2] +
                                    Real(70)*V_ghost[(j + 3)*(d_num_species + 3) + d_num_species + 2] -
                                    Real(63)/Real(2)*V_ghost[(j + 2)*(d_num_species + 3) + d_num_species + 2] +
                                    Real(42)/Real(5)*V_ghost[(j + 1)*(d_num_species + 3) + d_num_species + 2] +
                                    Real(42)*Real(dx[1])*dV_dy[d_num_species + 2];
                            }
                            
                            /*
                             * Compute the mixture density.
                             */
                            
                            Real rho_ghost = Real(0);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                rho_ghost += V_ghost[j*(d_num_species + 3) + si];
                            }
                            
                            /*
                             * Compute the mass fractions.
                             */
                            
                            std::vector<Real> Y_ghost;
                            Y_ghost.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_ghost.push_back(V_ghost[j*(d_num_species + 3) + si]/rho_ghost);
                            }
                            
                            /*
                             * Get the pointers to the mass fractions.
                             */
                            
                            std::vector<const Real*> Y_ghost_ptr;
                            Y_ghost_ptr.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_ghost_ptr.push_back(&Y_ghost[si]);
                            }
                            
                            for(int si=0; si < d_num_species; si++)
                            {
                                Q[si][idx_cell_rho_Y] = V_ghost[j*(d_num_species + 3) + si];
                            }
                            
                            Q[d_num_species][idx_cell_mom]     = rho_ghost*V_ghost[j*(d_num_species + 3) + d_num_species];
                            Q[d_num_species + 1][idx_cell_mom] = rho_ghost*V_ghost[j*(d_num_species + 3) + d_num_species + 1];
                            
                            const Real epsilon = d_equation_of_state_mixing_rules->
                                getInternalEnergy(
                                    &rho_ghost,
                                    &V_ghost[j*(d_num_species + 3) + d_num_species + 2],
                                    Y_ghost_ptr);
                            
                            const Real E = rho_ghost*epsilon +
                                half*(Q[d_num_species][idx_cell_mom]*Q[d_num_species][idx_cell_mom] +
                                    Q[d_num_species + 1][idx_cell_mom]*Q[d_num_species + 1][idx_cell_mom])/
                                    rho_ghost;
                            
                            Q[d_num_species + 2][idx_cell_E] = E;
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
                            << ": FlowModelBoundaryUtilitiesFourEqnConservative::fill2dEdgeBoundaryData()\n"
                            << "Non-reflecting outflow BC doesn't support more than six ghost cells yet!");
                    }
                    
                    for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                    {
                        // Get the grid spacing.
                        const double* const dx = patch_geom->getDx();
                        
                        const int idx_cell_rho_Y_y_B = (i + num_subghosts_conservative_var[0][0]) +
                            (interior_box_hi_idx[1] + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0];
                        
                        const int idx_cell_rho_Y_y_BB = (i + num_subghosts_conservative_var[0][0]) +
                            (interior_box_hi_idx[1] - 1 + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0];
                        
                        const int idx_cell_rho_Y_y_BBB = (i + num_subghosts_conservative_var[0][0]) +
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
                        
                        std::vector<Real> rho_Y_y_B;
                        std::vector<Real> rho_Y_y_BB;
                        std::vector<Real> rho_Y_y_BBB;
                        rho_Y_y_B.reserve(d_num_species);
                        rho_Y_y_BB.reserve(d_num_species);
                        rho_Y_y_BBB.reserve(d_num_species);
                        for (int si = 0; si < d_num_species; si++)
                        {
                            rho_Y_y_B.push_back(Q[si][idx_cell_rho_Y_y_B]);
                            rho_Y_y_BB.push_back(Q[si][idx_cell_rho_Y_y_BB]);
                            rho_Y_y_BBB.push_back(Q[si][idx_cell_rho_Y_y_BBB]);
                        }
                        
                        /*
                         * Compute the mixture density.
                         */
                        
                        Real rho_y_B   = Real(0);
                        Real rho_y_BB  = Real(0);
                        Real rho_y_BBB = Real(0);
                        for (int si = 0; si < d_num_species; si++)
                        {
                            rho_y_B   += Q[si][idx_cell_rho_Y_y_B];
                            rho_y_BB  += Q[si][idx_cell_rho_Y_y_BB];
                            rho_y_BBB += Q[si][idx_cell_rho_Y_y_BBB];
                        }
                        
                        /*
                         * Compute the mass fractions.
                         */
                        
                        std::vector<Real> Y_y_B;
                        std::vector<Real> Y_y_BB;
                        std::vector<Real> Y_y_BBB;
                        Y_y_B.reserve(d_num_species);
                        Y_y_BB.reserve(d_num_species);
                        Y_y_BBB.reserve(d_num_species);
                        for (int si = 0; si < d_num_species; si++)
                        {
                            Y_y_B.push_back(Q[si][idx_cell_rho_Y_y_B]/rho_y_B);
                            Y_y_BB.push_back(Q[si][idx_cell_rho_Y_y_BB]/rho_y_BB);
                            Y_y_BBB.push_back(Q[si][idx_cell_rho_Y_y_BBB]/rho_y_BBB);
                        }
                        
                        /*
                         * Get the pointers to the mass fractions.
                         */
                        
                        std::vector<const Real*> Y_y_B_ptr;
                        std::vector<const Real*> Y_y_BB_ptr;
                        std::vector<const Real*> Y_y_BBB_ptr;
                        Y_y_B_ptr.reserve(d_num_species);
                        Y_y_BB_ptr.reserve(d_num_species);
                        Y_y_BBB_ptr.reserve(d_num_species);
                        for (int si = 0; si < d_num_species; si++)
                        {
                            Y_y_B_ptr.push_back(&Y_y_B[si]);
                            Y_y_BB_ptr.push_back(&Y_y_BB[si]);
                            Y_y_BBB_ptr.push_back(&Y_y_BBB[si]);
                        }
                       
                        const Real u_y_B   = Q[d_num_species][idx_cell_mom_y_B]/rho_y_B;
                        const Real u_y_BB  = Q[d_num_species][idx_cell_mom_y_BB]/rho_y_BB;
                        const Real u_y_BBB = Q[d_num_species][idx_cell_mom_y_BBB]/rho_y_BBB;
                        
                        const Real v_y_B   = Q[d_num_species + 1][idx_cell_mom_y_B]/rho_y_B;
                        const Real v_y_BB  = Q[d_num_species + 1][idx_cell_mom_y_BB]/rho_y_BB;
                        const Real v_y_BBB = Q[d_num_species + 1][idx_cell_mom_y_BBB]/rho_y_BBB;
                       
                        const Real half = Real(1)/Real(2);
                        const Real epsilon_y_B   = Q[d_num_species + 2][idx_cell_E_y_B]/rho_y_B - half*(u_y_B*u_y_B + v_y_B*v_y_B);
                        const Real epsilon_y_BB  = Q[d_num_species + 2][idx_cell_E_y_BB]/rho_y_BB - half*(u_y_BB*u_y_BB + v_y_BB*v_y_BB);
                        const Real epsilon_y_BBB = Q[d_num_species + 2][idx_cell_E_y_BBB]/rho_y_BBB - half*(u_y_BBB*u_y_BBB + v_y_BBB*v_y_BBB);
                        
                        Real p_y_B = d_equation_of_state_mixing_rules->
                            getPressure(
                                &rho_y_B,
                                &epsilon_y_B,
                                Y_y_B_ptr);
                        
                        Real p_y_BB = d_equation_of_state_mixing_rules->
                            getPressure(
                                &rho_y_BB,
                                &epsilon_y_BB,
                                Y_y_BB_ptr);
                        
                        Real p_y_BBB = d_equation_of_state_mixing_rules->
                            getPressure(
                                &rho_y_BBB,
                                &epsilon_y_BBB,
                                Y_y_BBB_ptr);
                        
                        /*
                         * Compute derivatives in y-direction.
                         */
                        
                        std::vector<Real> drho_Y_dy;
                        drho_Y_dy.reserve(d_num_species);
                        for (int si = 0; si < d_num_species; si++)
                        {
                            drho_Y_dy.push_back((Q[si][idx_cell_rho_Y_y_BBB] - Real(4)*Q[si][idx_cell_rho_Y_y_BB] +
                                Real(3)*Q[si][idx_cell_rho_Y_y_B])/(Real(2)*Real(dx[1])));
                        }
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
                            
                            const int idx_cell_rho_Y_x_R = (i + 1 + num_subghosts_conservative_var[0][0]) +
                                (interior_box_hi_idx[1] + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_mom_x_R = (i + 1 + num_subghosts_conservative_var[1][0]) +
                                (interior_box_hi_idx[1] + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0];
                            
                            const int idx_cell_E_x_R = (i + 1 + num_subghosts_conservative_var[2][0]) +
                                (interior_box_hi_idx[1] + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0];
                            
                            /*
                             * Compute the mixture density.
                             */
                            
                            Real rho_x_R = Real(0);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                rho_x_R += Q[si][idx_cell_rho_Y_x_R];
                            }
                            
                            /*
                             * Compute the mass fractions.
                             */
                            
                            std::vector<Real> Y_x_R;
                            Y_x_R.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_x_R.push_back(Q[si][idx_cell_rho_Y_x_R]/rho_x_R);
                            }
                            
                            /*
                             * Get the pointers to the mass fractions.
                             */
                            
                            std::vector<const Real*> Y_x_R_ptr;
                            Y_x_R_ptr.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_x_R_ptr.push_back(&Y_x_R[si]);
                            }
                            
                            const Real u_x_R = Q[d_num_species][idx_cell_mom_x_R]/rho_x_R;
                            const Real v_x_R = Q[d_num_species + 1][idx_cell_mom_x_R]/rho_x_R;
                            const Real epsilon_x_R = Q[d_num_species + 2][idx_cell_E_x_R]/rho_x_R - half*(u_x_R*u_x_R + v_x_R*v_x_R);
                            
                            Real p_x_R = d_equation_of_state_mixing_rules->
                                getPressure(
                                    &rho_x_R,
                                    &epsilon_x_R,
                                    Y_x_R_ptr);
                            
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
                            
                            const int idx_cell_rho_Y_x_L = (i - 1 + num_subghosts_conservative_var[0][0]) +
                                (interior_box_hi_idx[1] + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_mom_x_L = (i - 1 + num_subghosts_conservative_var[1][0]) +
                                (interior_box_hi_idx[1] + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0];
                            
                            const int idx_cell_E_x_L = (i - 1 + num_subghosts_conservative_var[2][0]) +
                                (interior_box_hi_idx[1] + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0];
                            
                            /*
                             * Compute the mixture density.
                             */
                            
                            Real rho_x_L = Real(0);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                rho_x_L += Q[si][idx_cell_rho_Y_x_L];
                            }
                            
                            /*
                             * Compute the mass fractions.
                             */
                            
                            std::vector<Real> Y_x_L;
                            Y_x_L.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_x_L.push_back(Q[si][idx_cell_rho_Y_x_L]/rho_x_L);
                            }
                            
                            /*
                             * Get the pointers to the mass fractions.
                             */
                            
                            std::vector<const Real*> Y_x_L_ptr;
                            Y_x_L_ptr.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_x_L_ptr.push_back(&Y_x_L[si]);
                            }
                            
                            const Real u_x_L = Q[d_num_species][idx_cell_mom_x_L]/rho_x_L;
                            const Real v_x_L = Q[d_num_species + 1][idx_cell_mom_x_L]/rho_x_L;
                            const Real epsilon_x_L = Q[d_num_species + 2][idx_cell_E_x_L]/rho_x_L - half*(u_x_L*u_x_L + v_x_L*v_x_L);
                            
                            Real p_x_L = d_equation_of_state_mixing_rules->
                                getPressure(
                                    &rho_x_L,
                                    &epsilon_x_L,
                                    Y_x_L_ptr);
                            
                            // One-sided derivatives.
                            du_dx = (u_y_B - u_x_L)/Real(dx[0]);
                            dv_dx = (v_y_B - v_x_L)/Real(dx[0]);
                            dp_dx = (p_y_B - p_x_L)/Real(dx[0]);
                        }
                        else
                        {
                            const int idx_cell_rho_Y_x_L = (i - 1 + num_subghosts_conservative_var[0][0]) +
                                (interior_box_hi_idx[1] + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_rho_Y_x_R = (i + 1 + num_subghosts_conservative_var[0][0]) +
                                (interior_box_hi_idx[1] + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_mom_x_L = (i - 1 + num_subghosts_conservative_var[1][0]) +
                                (interior_box_hi_idx[1] + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0];
                            
                            const int idx_cell_mom_x_R = (i + 1 + num_subghosts_conservative_var[1][0]) +
                                (interior_box_hi_idx[1] + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0];
                            
                            const int idx_cell_E_x_L = (i - 1 + num_subghosts_conservative_var[2][0]) +
                                (interior_box_hi_idx[1] + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0];
                            
                            const int idx_cell_E_x_R = (i + 1 + num_subghosts_conservative_var[2][0]) +
                                (interior_box_hi_idx[1] + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0];
                            
                            /*
                             * Compute the mixture density.
                             */
                            
                            Real rho_x_L = Real(0);
                            Real rho_x_R = Real(0);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                rho_x_L += Q[si][idx_cell_rho_Y_x_L];
                                rho_x_R += Q[si][idx_cell_rho_Y_x_R];
                            }
                            
                            /*
                             * Compute the mass fractions.
                             */
                            
                            std::vector<Real> Y_x_L;
                            std::vector<Real> Y_x_R;
                            Y_x_L.reserve(d_num_species);
                            Y_x_R.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_x_L.push_back(Q[si][idx_cell_rho_Y_x_L]/rho_x_L);
                                Y_x_R.push_back(Q[si][idx_cell_rho_Y_x_R]/rho_x_R);
                            }
                            
                            /*
                             * Get the pointers to the mass fractions.
                             */
                            
                            std::vector<const Real*> Y_x_L_ptr;
                            std::vector<const Real*> Y_x_R_ptr;
                            Y_x_L_ptr.reserve(d_num_species);
                            Y_x_R_ptr.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_x_L_ptr.push_back(&Y_x_L[si]);
                                Y_x_R_ptr.push_back(&Y_x_R[si]);
                            }
                            
                            const Real u_x_L = Q[d_num_species][idx_cell_mom_x_L]/rho_x_L;
                            const Real u_x_R = Q[d_num_species][idx_cell_mom_x_R]/rho_x_R;
                            
                            const Real v_x_L = Q[d_num_species + 1][idx_cell_mom_x_L]/rho_x_L;
                            const Real v_x_R = Q[d_num_species + 1][idx_cell_mom_x_R]/rho_x_R;
                            
                            const Real epsilon_x_L = Q[d_num_species + 2][idx_cell_E_x_L]/rho_x_L - half*(u_x_L*u_x_L + v_x_L*v_x_L);
                            const Real epsilon_x_R = Q[d_num_species + 2][idx_cell_E_x_R]/rho_x_R - half*(u_x_R*u_x_R + v_x_R*v_x_R);
                            
                            Real p_x_L = d_equation_of_state_mixing_rules->
                                getPressure(
                                    &rho_x_L,
                                    &epsilon_x_L,
                                    Y_x_L_ptr);
                            
                            Real p_x_R = d_equation_of_state_mixing_rules->
                                getPressure(
                                    &rho_x_R,
                                    &epsilon_x_R,
                                    Y_x_R_ptr);
                            
                            // Central derivatives.
                            du_dx = (u_x_R - u_x_L)/(Real(2)*Real(dx[0]));
                            dv_dx = (v_x_R - v_x_L)/(Real(2)*Real(dx[0]));
                            dp_dx = (p_x_R - p_x_L)/(Real(2)*Real(dx[0]));
                        }
                        
                        // Compute sound speed.
                        
                        const Real Gamma_y_B = d_equation_of_state_mixing_rules->getGruneisenParameter(
                            &rho_y_B,
                            &p_y_B,
                            Y_y_B_ptr);
                        
                        const std::vector<Real> Psi_y_B = d_equation_of_state_mixing_rules->
                            getPressureDerivativeWithPartialDensities(
                                    &rho_y_B,
                                    &p_y_B,
                                    Y_y_B_ptr);
                        
                        Real c_y_B = Gamma_y_B*p_y_B/rho_y_B;
                        for (int si = 0; si < d_num_species; si++)
                        {
                            c_y_B += Y_y_B[si]*Psi_y_B[si];
                        }
                        c_y_B = std::sqrt(c_y_B);
                        
                        const Real lambda_1 = v_y_B - c_y_B;
                        
                        // Compute vector Lambda^(-1) * L.
                        
                        Real Lambda_inv_L[d_num_species + 3];
                        
                        const Real& p_t         = d_bdry_edge_nonreflecting_outflow_p_t[edge_loc];
                        const Real& sigma       = d_bdry_edge_nonreflecting_outflow_sigma[edge_loc];
                        const Real& beta        = d_bdry_edge_nonreflecting_outflow_beta[edge_loc];
                        const Real& length_char = d_bdry_edge_nonreflecting_outflow_length_char[edge_loc];
                        
                        const Real T_1 = u_y_B*(dp_dx - rho_y_B*c_y_B*dv_dx) + rho_y_B*c_y_B*c_y_B*du_dx;
                        
                        const Real M_sq = (v_y_B*v_y_B + u_y_B*u_y_B)/(c_y_B*c_y_B);
                        const Real K = sigma*c_y_B*(Real(1) - M_sq)/length_char;
                        
                        Lambda_inv_L[0] = (Real(1)/lambda_1)*(K*(p_y_B - p_t) - (Real(1) - beta)*T_1);
                        Lambda_inv_L[1] = du_dy;
                        for (int si = 0; si < d_num_species; si++)
                        {
                            Lambda_inv_L[si + 2] = c_y_B*c_y_B*drho_Y_dy[si] - Y_y_B[si]*dp_dy;
                        }
                        Lambda_inv_L[d_num_species + 2] = dp_dy + rho_y_B*c_y_B*dv_dy;
                        
                       // Compute dV_dy.
                        
                        const Real c_sq_inv  = Real(1)/(c_y_B*c_y_B);
                        const Real rho_c_inv = Real(1)/(rho_y_B*c_y_B);
                        
                        Real dV_dy[d_num_species + 3];
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            dV_dy[si] = half*c_sq_inv*Y_y_B[si]*(Lambda_inv_L[0] + Lambda_inv_L[d_num_species + 2]) +
                                c_sq_inv*Lambda_inv_L[si+2];
                        }
                        dV_dy[d_num_species] = Lambda_inv_L[1];
                        dV_dy[d_num_species + 1] = half*rho_c_inv*(-Lambda_inv_L[0] + Lambda_inv_L[d_num_species + 2]);
                        dV_dy[d_num_species + 2] = half*(Lambda_inv_L[0] + Lambda_inv_L[d_num_species + 2]); 
                        
                        Real V_ghost[(d_num_species + 3)*num_ghosts_to_fill];
                        
                        for (int j = 0; j < num_ghosts_to_fill; j++)
                        {
                            const int idx_cell_rho_Y = (i + num_subghosts_conservative_var[0][0]) +
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
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    V_ghost[j*(d_num_species + 3) + si] = rho_Y_y_BB[si] + Real(2)*Real(dx[1])*dV_dy[si];
                                }
                                
                                V_ghost[j*(d_num_species + 3) + d_num_species]     = u_y_BB   + Real(2)*Real(dx[1])*dV_dy[d_num_species];
                                V_ghost[j*(d_num_species + 3) + d_num_species + 1] = v_y_BB   + Real(2)*Real(dx[1])*dV_dy[d_num_species + 1];
                                V_ghost[j*(d_num_species + 3) + d_num_species + 2] = p_y_BB   + Real(2)*Real(dx[1])*dV_dy[d_num_species + 2];
                            }
                            else if (j == 1)
                            {
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    V_ghost[j*(d_num_species + 3) + si] = -Real(2)*rho_Y_y_BB[si] - Real(3)*rho_Y_y_B[si] +
                                        Real(6)*V_ghost[(j - 1)*(d_num_species + 3) + si] - Real(6)*Real(dx[1])*dV_dy[si];
                                }
                                
                                V_ghost[j*(d_num_species + 3) + d_num_species] = -Real(2)*u_y_BB - Real(3)*u_y_B +
                                    Real(6)*V_ghost[(j - 1)*(d_num_species + 3) + d_num_species] - Real(6)*Real(dx[1])*dV_dy[d_num_species];
                                
                                V_ghost[j*(d_num_species + 3) + d_num_species + 1] = -Real(2)*v_y_BB - Real(3)*v_y_B +
                                    Real(6)*V_ghost[(j - 1)*(d_num_species + 3) + d_num_species + 1] - Real(6)*Real(dx[1])*dV_dy[d_num_species + 1];
                                
                                V_ghost[j*(d_num_species + 3) +d_num_species + 2] = -Real(2)*p_y_BB - Real(3)*p_y_B +
                                    Real(6)*V_ghost[(j - 1)*(d_num_species + 3) + d_num_species + 2] - Real(6)*Real(dx[1])*dV_dy[d_num_species + 2];
                            }
                            else if (j == 2)
                            {
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    V_ghost[j*(d_num_species + 3) + si] = Real(3)*rho_Y_y_BB[si] + Real(10)*rho_Y_y_B[si] -
                                    Real(18)*V_ghost[(j - 2)*(d_num_species + 3) + si] +
                                    Real(6)*V_ghost[(j - 1)*(d_num_species + 3) + si] +
                                    Real(12)*Real(dx[1])*dV_dy[si];
                                }
                                
                                V_ghost[j*(d_num_species + 3) + d_num_species] = Real(3)*u_y_BB + Real(10)*u_y_B -
                                    Real(18)*V_ghost[(j - 2)*(d_num_species + 3) + d_num_species] +
                                    Real(6)*V_ghost[(j - 1)*(d_num_species + 3) + d_num_species] +
                                    Real(12)*Real(dx[1])*dV_dy[d_num_species];
                                
                                V_ghost[j*(d_num_species + 3) + d_num_species + 1] = Real(3)*v_y_BB + Real(10)*v_y_B -
                                    Real(18)*V_ghost[(j - 2)*(d_num_species + 3) + d_num_species + 1] +
                                    Real(6)*V_ghost[(j - 1)*(d_num_species + 3) + d_num_species + 1] +
                                    Real(12)*Real(dx[1])*dV_dy[d_num_species + 1];
                                
                                V_ghost[j*(d_num_species + 3) + d_num_species + 2] = Real(3)*p_y_BB + Real(10)*p_y_B -
                                    Real(18)*V_ghost[(j - 2)*(d_num_species + 3) + d_num_species + 2] +
                                    Real(6)*V_ghost[(j - 1)*(d_num_species + 3) + d_num_species + 2] +
                                    Real(12)*Real(dx[1])*dV_dy[d_num_species + 2];
                            }
                            else if (j == 3)
                            {
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    V_ghost[j*(d_num_species + 3) + si] = -Real(4)*rho_Y_y_BB[si] -
                                        Real(65)/Real(3)*rho_Y_y_B[si] +
                                        Real(40)*V_ghost[(j - 3)*(d_num_species + 3) + si] -
                                        Real(20)*V_ghost[(j - 2)*(d_num_species + 3) + si] +
                                        Real(20)/Real(3)*V_ghost[(j - 1)*(d_num_species + 3) + si] -
                                        Real(20)*Real(dx[1])*dV_dy[si];
                                } 
                                
                                V_ghost[j*(d_num_species + 3) + d_num_species] = -Real(4)*u_y_BB -
                                    Real(65)/Real(3)*u_y_B +
                                    Real(40)*V_ghost[(j - 3)*(d_num_species + 3) + d_num_species] -
                                    Real(20)*V_ghost[(j - 2)*(d_num_species + 3) + d_num_species] +
                                    Real(20)/Real(3)*V_ghost[(j - 1)*(d_num_species + 3) + d_num_species] -
                                    Real(20)*Real(dx[1])*dV_dy[d_num_species];
                                
                                V_ghost[j*(d_num_species + 3) + d_num_species + 1] = -Real(4)*v_y_BB -
                                    Real(65)/Real(3)*v_y_B +
                                    Real(40)*V_ghost[(j - 3)*(d_num_species + 3) + d_num_species + 1] -
                                    Real(20)*V_ghost[(j - 2)*(d_num_species + 3) + d_num_species + 1] +
                                    Real(20)/Real(3)*V_ghost[(j - 1)*(d_num_species + 3) + d_num_species + 1] -
                                    Real(20)*Real(dx[1])*dV_dy[d_num_species + 1];
                                
                                V_ghost[j*(d_num_species + 3) + d_num_species + 2] = -Real(4)*p_y_BB -
                                    Real(65)/Real(3)*p_y_B +
                                    Real(40)*V_ghost[(j - 3)*(d_num_species + 3) + d_num_species + 2] -
                                    Real(20)*V_ghost[(j - 2)*(d_num_species + 3) + d_num_species + 2] +
                                    Real(20)/Real(3)*V_ghost[(j - 1)*(d_num_species + 3) + d_num_species + 2] -
                                    Real(20)*Real(dx[1])*dV_dy[d_num_species + 2];
                            }
                            else if (j == 4)
                            {
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    V_ghost[j*(d_num_species + 3) + si] = Real(5)*rho_Y_y_BB[si] +
                                        Real(77)/Real(2)*rho_Y_y_B[si] -
                                        Real(75)*V_ghost[(j - 4)*(d_num_species + 3) + si] +
                                        Real(50)*V_ghost[(j - 3)*(d_num_species + 3) + si] -
                                        Real(25)*V_ghost[(j - 2)*(d_num_species + 3) + si] +
                                        Real(15)/Real(2)*V_ghost[(j - 1)*(d_num_species + 3) + si] +
                                        Real(30)*Real(dx[1])*dV_dy[si];
                                }
                                
                                V_ghost[j*(d_num_species + 3) + d_num_species] = Real(5)*u_y_BB +
                                    Real(77)/Real(2)*u_y_B -
                                    Real(75)*V_ghost[(j - 4)*(d_num_species + 3) + d_num_species] +
                                    Real(50)*V_ghost[(j - 3)*(d_num_species + 3) + d_num_species] -
                                    Real(25)*V_ghost[(j - 2)*(d_num_species + 3) + d_num_species] +
                                    Real(15)/Real(2)*V_ghost[(j - 1)*(d_num_species + 3) + d_num_species] +
                                    Real(30)*Real(dx[1])*dV_dy[d_num_species];
                                
                                V_ghost[j*(d_num_species + 3) + d_num_species + 1] = Real(5)*v_y_BB +
                                    Real(77)/Real(2)*v_y_B -
                                    Real(75)*V_ghost[(j - 4)*(d_num_species + 3) + d_num_species + 1] +
                                    Real(50)*V_ghost[(j - 3)*(d_num_species + 3) + d_num_species + 1] -
                                    Real(25)*V_ghost[(j - 2)*(d_num_species + 3) + d_num_species + 1] +
                                    Real(15)/Real(2)*V_ghost[(j - 1)*(d_num_species + 3) + d_num_species + 1] +
                                    Real(30)*Real(dx[1])*dV_dy[d_num_species + 1];
                                
                                V_ghost[j*(d_num_species + 3) + d_num_species + 2] = Real(5)*p_y_BB +
                                    Real(77)/Real(2)*p_y_B -
                                    Real(75)*V_ghost[(j - 4)*(d_num_species + 3) + d_num_species + 2] +
                                    Real(50)*V_ghost[(j - 3)*(d_num_species + 3) + d_num_species + 2] -
                                    Real(25)*V_ghost[(j - 2)*(d_num_species + 3) + d_num_species + 2] +
                                    Real(15)/Real(2)*V_ghost[(j - 1)*(d_num_species + 3) + d_num_species + 2] +
                                    Real(30)*Real(dx[1])*dV_dy[d_num_species + 2];
                            }
                            else if (j == 5)
                            {
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    V_ghost[j*(d_num_species + 3) + si] = -Real(6)*rho_Y_y_BB[si] -
                                        Real(609)/Real(10)*rho_Y_y_B[si] +
                                        Real(126)*V_ghost[(j - 5)*(d_num_species + 3) + si] -
                                        Real(105)*V_ghost[(j - 4)*(d_num_species + 3) + si] +
                                        Real(70)*V_ghost[(j - 3)*(d_num_species + 3) + si] -
                                        Real(63)/Real(2)*V_ghost[(j - 2)*(d_num_species + 3) + si] +
                                        Real(42)/Real(5)*V_ghost[(j - 1)*(d_num_species + 3) + si] -
                                        Real(42)*Real(dx[1])*dV_dy[si];
                                }
                                
                                V_ghost[j*(d_num_species + 3) + d_num_species] = -Real(6)*u_y_BB -
                                    Real(609)/Real(10)*u_y_B +
                                    Real(126)*V_ghost[(j - 5)*(d_num_species + 3) + d_num_species] -
                                    Real(105)*V_ghost[(j - 4)*(d_num_species + 3) + d_num_species] +
                                    Real(70)*V_ghost[(j - 3)*(d_num_species + 3) + d_num_species] -
                                    Real(63)/Real(2)*V_ghost[(j - 2)*(d_num_species + 3) + d_num_species] +
                                    Real(42)/Real(5)*V_ghost[(j - 1)*(d_num_species + 3) + d_num_species] -
                                    Real(42)*Real(dx[1])*dV_dy[d_num_species];
                                
                                V_ghost[j*(d_num_species + 3) + d_num_species + 1] = -Real(6)*v_y_BB -
                                    Real(609)/Real(10)*v_y_B +
                                    Real(126)*V_ghost[(j - 5)*(d_num_species + 3) + d_num_species + 1] -
                                    Real(105)*V_ghost[(j - 4)*(d_num_species + 3) + d_num_species + 1] +
                                    Real(70)*V_ghost[(j - 3)*(d_num_species + 3) + d_num_species + 1] -
                                    Real(63)/Real(2)*V_ghost[(j - 2)*(d_num_species + 3) + d_num_species + 1] +
                                    Real(42)/Real(5)*V_ghost[(j - 1)*(d_num_species + 3) + d_num_species + 1] -
                                    Real(42)*Real(dx[1])*dV_dy[d_num_species + 1];
                                
                                V_ghost[j*(d_num_species + 3) + d_num_species + 2] = -Real(6)*p_y_BB -
                                    Real(609)/Real(10)*p_y_B +
                                    Real(126)*V_ghost[(j - 5)*(d_num_species + 3) + d_num_species + 2] -
                                    Real(105)*V_ghost[(j - 4)*(d_num_species + 3) + d_num_species + 2] +
                                    Real(70)*V_ghost[(j - 3)*(d_num_species + 3) + d_num_species + 2] -
                                    Real(63)/Real(2)*V_ghost[(j - 2)*(d_num_species + 3) + d_num_species + 2] +
                                    Real(42)/Real(5)*V_ghost[(j - 1)*(d_num_species + 3) + d_num_species + 2] -
                                    Real(42)*Real(dx[1])*dV_dy[d_num_species + 2];
                            }
                            
                            /*
                             * Compute the mixture density.
                             */
                            
                            Real rho_ghost = Real(0);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                rho_ghost += V_ghost[j*(d_num_species + 3) + si];
                            }
                            
                            /*
                             * Compute the mass fractions.
                             */
                            
                            std::vector<Real> Y_ghost;
                            Y_ghost.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_ghost.push_back(V_ghost[j*(d_num_species + 3) + si]/rho_ghost);
                            }
                            
                            /*
                             * Get the pointers to the mass fractions.
                             */
                            
                            std::vector<const Real*> Y_ghost_ptr;
                            Y_ghost_ptr.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_ghost_ptr.push_back(&Y_ghost[si]);
                            }
                            
                            for(int si=0; si < d_num_species; si++)
                            {
                                Q[si][idx_cell_rho_Y] = V_ghost[j*(d_num_species + 3) + si];
                            }
                            
                            Q[d_num_species][idx_cell_mom]     = rho_ghost*V_ghost[j*(d_num_species + 3) + d_num_species];
                            Q[d_num_species + 1][idx_cell_mom] = rho_ghost*V_ghost[j*(d_num_species + 3) + d_num_species + 1];
                            
                            const Real epsilon = d_equation_of_state_mixing_rules->
                                getInternalEnergy(
                                    &rho_ghost,
                                    &V_ghost[j*(d_num_species + 3) + d_num_species + 2],
                                    Y_ghost_ptr);
                            
                            const Real E = rho_ghost*epsilon +
                                half*(Q[d_num_species][idx_cell_mom]*Q[d_num_species][idx_cell_mom] +
                                    Q[d_num_species + 1][idx_cell_mom]*Q[d_num_species + 1][idx_cell_mom])/
                                    rho_ghost;
                            
                            Q[d_num_species + 2][idx_cell_E] = E;
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
            TBOX_ERROR("FlowModelBoundaryUtilitiesFourEqnConservative::fill2dEdgeBoundaryData()\n"
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
FlowModelBoundaryUtilitiesFourEqnConservative::fill2dNodeBoundaryData(
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
                
                if (bdry_node_conds[node_loc] == BDRY_COND::FLOW_MODEL::XADIABATIC_NO_SLIP)
                {
                    for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                    {
                        for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                        {
                            const int idx_cell_rho_Y = (i + num_subghosts_conservative_var[0][0]) +
                                (j + num_subghosts_conservative_var[0][1])*
                                    subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_mom = (i + num_subghosts_conservative_var[1][0]) +
                                (j + num_subghosts_conservative_var[1][1])*
                                    subghostcell_dims_conservative_var[1][0];
                            
                            const int idx_cell_E = (i + num_subghosts_conservative_var[2][0]) +
                                (j + num_subghosts_conservative_var[2][1])*
                                    subghostcell_dims_conservative_var[2][0];
                            
                            int idx_cell_pivot_rho_Y = idx_cell_rho_Y;
                            int idx_cell_pivot_mom = idx_cell_mom;
                            int idx_cell_pivot_E = idx_cell_E;
                            
                            if (edge_loc_0 == BDRY_LOC::XLO)
                            {
                                idx_cell_pivot_rho_Y = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
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
                                idx_cell_pivot_rho_Y = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
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
                             * Compute the mixture density of the pivot.
                             */
                            
                            Real rho_pivot = Real(0);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                rho_pivot += Q[si][idx_cell_pivot_rho_Y];
                            }
                            
                            /*
                             * Compute the mass fractions of the pivot.
                             */
                            
                            std::vector<Real> Y_pivot;
                            Y_pivot.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_pivot.push_back(Q[si][idx_cell_pivot_rho_Y]/rho_pivot);
                            }
                            
                            /*
                             * Get the pointers to the mass fractions of the pivot.
                             */
                            
                            std::vector<const Real*> Y_pivot_ptr;
                            Y_pivot_ptr.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_pivot_ptr.push_back(&Y_pivot[si]);
                            }
                            
                            /*
                             * Set the values for partial densities and momentum.
                             */
                            
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Q[si][idx_cell_rho_Y] = Q[si][idx_cell_pivot_rho_Y];
                            }
                            Q[d_num_species][idx_cell_mom] = -Q[d_num_species][idx_cell_pivot_mom] +
                                Real(2)*rho_pivot*d_bdry_edge_adiabatic_no_slip_vel[edge_loc_0*2];
                            Q[d_num_species + 1][idx_cell_mom] = -Q[d_num_species + 1][idx_cell_pivot_mom] +
                                Real(2)*rho_pivot*d_bdry_edge_adiabatic_no_slip_vel[edge_loc_0*2 + 1];
                            
                            /*
                             * Set the values for total internal energy.
                             */
                            
                            Real epsilon_pivot = (Q[d_num_species + 2][idx_cell_pivot_E] -
                                Real(1)/Real(2)*(Q[d_num_species][idx_cell_pivot_mom]*Q[d_num_species][idx_cell_pivot_mom] +
                                Q[d_num_species + 1][idx_cell_pivot_mom]*Q[d_num_species + 1][idx_cell_pivot_mom])/
                                rho_pivot)/rho_pivot;
                            
                            Real p_pivot = d_equation_of_state_mixing_rules->
                                getPressure(
                                    &rho_pivot,
                                    &epsilon_pivot,
                                    Y_pivot_ptr);
                            
                            Real T_pivot = d_equation_of_state_mixing_rules->
                                getTemperature(
                                    &rho_pivot,
                                    &p_pivot,
                                    Y_pivot_ptr);
                            
                            Real T = T_pivot;
                            
                            Real epsilon = d_equation_of_state_mixing_rules->
                                getInternalEnergyFromTemperature(
                                    &rho_pivot,
                                    &T,
                                    Y_pivot_ptr);
                            
                            Real E = rho_pivot*epsilon +
                                Real(1)/Real(2)*(Q[d_num_species][idx_cell_mom]*Q[d_num_species][idx_cell_mom] +
                                Q[d_num_species + 1][idx_cell_mom]*Q[d_num_species + 1][idx_cell_mom])/rho_pivot;
                            
                            Q[d_num_species + 2][idx_cell_E] = E;
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
                            const int idx_cell_rho_Y = (i + num_subghosts_conservative_var[0][0]) +
                                (j + num_subghosts_conservative_var[0][1])*
                                    subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_mom = (i + num_subghosts_conservative_var[1][0]) +
                                (j + num_subghosts_conservative_var[1][1])*
                                    subghostcell_dims_conservative_var[1][0];
                            
                            const int idx_cell_E = (i + num_subghosts_conservative_var[2][0]) +
                                (j + num_subghosts_conservative_var[2][1])*
                                    subghostcell_dims_conservative_var[2][0];
                            
                            int idx_cell_pivot_rho_Y = idx_cell_rho_Y;
                            int idx_cell_pivot_mom = idx_cell_mom;
                            int idx_cell_pivot_E = idx_cell_E;
                            
                            if (edge_loc_1 == BDRY_LOC::YLO)
                            {
                                idx_cell_pivot_rho_Y = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
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
                                idx_cell_pivot_rho_Y = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
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
                             * Compute the mixture density of the pivot.
                             */
                            
                            Real rho_pivot = Real(0);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                rho_pivot += Q[si][idx_cell_pivot_rho_Y];
                            }
                            
                            /*
                             * Compute the mass fractions of the pivot.
                             */
                            
                            std::vector<Real> Y_pivot;
                            Y_pivot.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_pivot.push_back(Q[si][idx_cell_pivot_rho_Y]/rho_pivot);
                            }
                            
                            /*
                             * Get the pointers to the mass fractions of the pivot.
                             */
                            
                            std::vector<const Real*> Y_pivot_ptr;
                            Y_pivot_ptr.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_pivot_ptr.push_back(&Y_pivot[si]);
                            }
                            
                            /*
                             * Set the values for partial densities and momentum.
                             */
                            
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Q[si][idx_cell_rho_Y] = Q[si][idx_cell_pivot_rho_Y];
                            }
                            Q[d_num_species][idx_cell_mom] = -Q[d_num_species][idx_cell_pivot_mom] +
                                Real(2)*rho_pivot*d_bdry_edge_adiabatic_no_slip_vel[edge_loc_1*2];
                            Q[d_num_species + 1][idx_cell_mom] = -Q[d_num_species + 1][idx_cell_pivot_mom] +
                                Real(2)*rho_pivot*d_bdry_edge_adiabatic_no_slip_vel[edge_loc_1*2 + 1];
                            
                            /*
                             * Set the values for total internal energy.
                             */
                            
                            Real epsilon_pivot = (Q[d_num_species + 2][idx_cell_pivot_E] -
                                Real(1)/Real(2)*(Q[d_num_species][idx_cell_pivot_mom]*Q[d_num_species][idx_cell_pivot_mom] +
                                Q[d_num_species + 1][idx_cell_pivot_mom]*Q[d_num_species + 1][idx_cell_pivot_mom])/
                                rho_pivot)/rho_pivot;
                            
                            Real p_pivot = d_equation_of_state_mixing_rules->
                                getPressure(
                                    &rho_pivot,
                                    &epsilon_pivot,
                                    Y_pivot_ptr);
                            
                            Real T_pivot = d_equation_of_state_mixing_rules->
                                getTemperature(
                                    &rho_pivot,
                                    &p_pivot,
                                    Y_pivot_ptr);
                            
                            Real T = T_pivot;
                            
                            Real epsilon = d_equation_of_state_mixing_rules->
                                getInternalEnergyFromTemperature(
                                    &rho_pivot,
                                    &T,
                                    Y_pivot_ptr);
                            
                            Real E = rho_pivot*epsilon +
                                Real(1)/Real(2)*(Q[d_num_species][idx_cell_mom]*Q[d_num_species][idx_cell_mom] +
                                Q[d_num_species + 1][idx_cell_mom]*Q[d_num_species + 1][idx_cell_mom])/rho_pivot;
                            
                            Q[d_num_species + 2][idx_cell_E] = E;
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
                            const int idx_cell_rho_Y = (i + num_subghosts_conservative_var[0][0]) +
                                (j + num_subghosts_conservative_var[0][1])*
                                    subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_mom = (i + num_subghosts_conservative_var[1][0]) +
                                (j + num_subghosts_conservative_var[1][1])*
                                    subghostcell_dims_conservative_var[1][0];
                            
                            const int idx_cell_E = (i + num_subghosts_conservative_var[2][0]) +
                                (j + num_subghosts_conservative_var[2][1])*
                                    subghostcell_dims_conservative_var[2][0];
                            
                            int idx_cell_pivot_rho_Y = idx_cell_rho_Y;
                            int idx_cell_pivot_mom = idx_cell_mom;
                            int idx_cell_pivot_E = idx_cell_E;
                            
                            if (edge_loc_0 == BDRY_LOC::XLO)
                            {
                                idx_cell_pivot_rho_Y = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
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
                                idx_cell_pivot_rho_Y = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
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
                             * Compute the mixture density of the pivot.
                             */
                            
                            Real rho_pivot = Real(0);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                rho_pivot += Q[si][idx_cell_pivot_rho_Y];
                            }
                            
                            /*
                             * Compute the mass fractions of the pivot.
                             */
                            
                            std::vector<Real> Y_pivot;
                            Y_pivot.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_pivot.push_back(Q[si][idx_cell_pivot_rho_Y]/rho_pivot);
                            }
                            
                            /*
                             * Get the pointers to the mass fractions of the pivot.
                             */
                            
                            std::vector<const Real*> Y_pivot_ptr;
                            Y_pivot_ptr.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_pivot_ptr.push_back(&Y_pivot[si]);
                            }
                            
                            /*
                             * Set the values for partial densities, momentum and total internal energy.
                             */
                            
                            Real epsilon_pivot = (Q[d_num_species + 2][idx_cell_pivot_E] -
                                Real(1)/Real(2)*(Q[d_num_species][idx_cell_pivot_mom]*Q[d_num_species][idx_cell_pivot_mom] +
                                Q[d_num_species + 1][idx_cell_pivot_mom]*Q[d_num_species + 1][idx_cell_pivot_mom])/
                                rho_pivot)/rho_pivot;
                            
                            Real p_pivot = d_equation_of_state_mixing_rules->
                                getPressure(
                                    &rho_pivot,
                                    &epsilon_pivot,
                                    Y_pivot_ptr);
                            
                            Real p = p_pivot;
                            
                            Real T_pivot = d_equation_of_state_mixing_rules->
                                getTemperature(
                                    &rho_pivot,
                                    &p_pivot,
                                    Y_pivot_ptr);
                            
                            Real T = -T_pivot + Real(2)*d_bdry_edge_isothermal_no_slip_T[edge_loc_0];
                            
                            Real rho = d_equation_of_state_mixing_rules->
                                getMixtureDensity(
                                    &p,
                                    &T,
                                    Y_pivot_ptr);
                            
                            Real u = -Q[d_num_species][idx_cell_pivot_mom]/rho_pivot +
                                Real(2)*d_bdry_edge_isothermal_no_slip_vel[edge_loc_0*2];
                            Real v = -Q[d_num_species + 1][idx_cell_pivot_mom]/rho_pivot +
                                Real(2)*d_bdry_edge_isothermal_no_slip_vel[edge_loc_0*2 + 1];
                            
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Q[si][idx_cell_rho_Y] = rho*Y_pivot[si];
                            }
                            Q[d_num_species][idx_cell_mom] = rho*u;
                            Q[d_num_species + 1][idx_cell_mom] = rho*v;
                            
                            Real epsilon = d_equation_of_state_mixing_rules->
                                getInternalEnergyFromTemperature(
                                    &rho,
                                    &T,
                                    Y_pivot_ptr);
                            
                            Real E = rho*epsilon +
                                Real(1)/Real(2)*(Q[d_num_species][idx_cell_mom]*Q[d_num_species][idx_cell_mom] + 
                                Q[d_num_species + 1][idx_cell_mom]*Q[d_num_species + 1][idx_cell_mom])/rho;
                            
                            Q[d_num_species + 2][idx_cell_E] = E;
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
                            const int idx_cell_rho_Y = (i + num_subghosts_conservative_var[0][0]) +
                                (j + num_subghosts_conservative_var[0][1])*
                                    subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_mom = (i + num_subghosts_conservative_var[1][0]) +
                                (j + num_subghosts_conservative_var[1][1])*
                                    subghostcell_dims_conservative_var[1][0];
                            
                            const int idx_cell_E = (i + num_subghosts_conservative_var[2][0]) +
                                (j + num_subghosts_conservative_var[2][1])*
                                    subghostcell_dims_conservative_var[2][0];
                            
                            int idx_cell_pivot_rho_Y = idx_cell_rho_Y;
                            int idx_cell_pivot_mom = idx_cell_mom;
                            int idx_cell_pivot_E = idx_cell_E;
                            
                            if (edge_loc_1 == BDRY_LOC::YLO)
                            {
                                idx_cell_pivot_rho_Y = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
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
                                idx_cell_pivot_rho_Y = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
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
                             * Compute the mixture density of the pivot.
                             */
                            
                            Real rho_pivot = Real(0);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                rho_pivot += Q[si][idx_cell_pivot_rho_Y];
                            }
                            
                            /*
                             * Compute the mass fractions of the pivot.
                             */
                            
                            std::vector<Real> Y_pivot;
                            Y_pivot.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_pivot.push_back(Q[si][idx_cell_pivot_rho_Y]/rho_pivot);
                            }
                            
                            /*
                             * Get the pointers to the mass fractions of the pivot.
                             */
                            
                            std::vector<const Real*> Y_pivot_ptr;
                            Y_pivot_ptr.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_pivot_ptr.push_back(&Y_pivot[si]);
                            }
                            
                            /*
                             * Set the values for partial densities, momentum and total internal energy.
                             */
                            
                            Real epsilon_pivot = (Q[d_num_species + 2][idx_cell_pivot_E] -
                                Real(1)/Real(2)*(Q[d_num_species][idx_cell_pivot_mom]*Q[d_num_species][idx_cell_pivot_mom] +
                                Q[d_num_species + 1][idx_cell_pivot_mom]*Q[d_num_species + 1][idx_cell_pivot_mom])/
                                rho_pivot)/rho_pivot;
                            
                            Real p_pivot = d_equation_of_state_mixing_rules->
                                getPressure(
                                    &rho_pivot,
                                    &epsilon_pivot,
                                    Y_pivot_ptr);
                            
                            Real p = p_pivot;
                            
                            Real T_pivot = d_equation_of_state_mixing_rules->
                                getTemperature(
                                    &rho_pivot,
                                    &p_pivot,
                                    Y_pivot_ptr);
                            
                            Real T = -T_pivot + Real(2)*d_bdry_edge_isothermal_no_slip_T[edge_loc_1];
                            
                            Real rho = d_equation_of_state_mixing_rules->
                                getMixtureDensity(
                                    &p,
                                    &T,
                                    Y_pivot_ptr);
                            
                            Real u = -Q[d_num_species][idx_cell_pivot_mom]/rho_pivot +
                                Real(2)*d_bdry_edge_isothermal_no_slip_vel[edge_loc_1*2];
                            Real v = -Q[d_num_species + 1][idx_cell_pivot_mom]/rho_pivot +
                                Real(2)*d_bdry_edge_isothermal_no_slip_vel[edge_loc_1*2 + 1];
                            
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Q[si][idx_cell_rho_Y] = rho*Y_pivot[si];
                            }
                            Q[d_num_species][idx_cell_mom] = rho*u;
                            Q[d_num_species + 1][idx_cell_mom] = rho*v;
                            
                            Real epsilon = d_equation_of_state_mixing_rules->
                                getInternalEnergyFromTemperature(
                                    &rho,
                                    &T,
                                    Y_pivot_ptr);
                            
                            Real E = rho*epsilon +
                                Real(1)/Real(2)*(Q[d_num_species][idx_cell_mom]*Q[d_num_species][idx_cell_mom] + 
                                Q[d_num_species + 1][idx_cell_mom]*Q[d_num_species + 1][idx_cell_mom])/rho;
                            
                            Q[d_num_species + 2][idx_cell_E] = E;
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
            TBOX_ERROR("FlowModelBoundaryUtilitiesFourEqnConservative::fill2dNodeBoundaryData\n"
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
FlowModelBoundaryUtilitiesFourEqnConservative::fill3dFaceBoundaryData(
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
                
                if (bdry_face_conds[face_loc] == BDRY_COND::FLOW_MODEL::ADIABATIC_NO_SLIP)
                {
                    for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                    {
                        for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                        {
                            for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                            {
                                const int idx_cell_rho_Y = (i + num_subghosts_conservative_var[0][0]) +
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
                                
                                int idx_cell_pivot_rho_Y = idx_cell_rho_Y;
                                int idx_cell_pivot_mom = idx_cell_mom;
                                int idx_cell_pivot_E = idx_cell_E;
                                
                                if (face_loc == BDRY_LOC::XLO)
                                {
                                    idx_cell_pivot_rho_Y = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
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
                                    idx_cell_pivot_rho_Y = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
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
                                    idx_cell_pivot_rho_Y = (i + num_subghosts_conservative_var[0][0]) +
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
                                    idx_cell_pivot_rho_Y = (i + num_subghosts_conservative_var[0][0]) +
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
                                    idx_cell_pivot_rho_Y = (i + num_subghosts_conservative_var[0][0]) +
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
                                    idx_cell_pivot_rho_Y = (i + num_subghosts_conservative_var[0][0]) +
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
                                 * Compute the mixture density of the pivot.
                                 */
                                
                                Real rho_pivot = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_pivot += Q[si][idx_cell_pivot_rho_Y];
                                }
                                
                                /*
                                 * Compute the mass fractions of the pivot.
                                 */
                                
                                std::vector<Real> Y_pivot;
                                Y_pivot.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_pivot.push_back(Q[si][idx_cell_pivot_rho_Y]/rho_pivot);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions of the pivot.
                                 */
                                
                                std::vector<const Real*> Y_pivot_ptr;
                                Y_pivot_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_pivot_ptr.push_back(&Y_pivot[si]);
                                }
                                
                                /*
                                 * Set the values for partial densities and momentum.
                                 */
                                
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Q[si][idx_cell_rho_Y] = Q[si][idx_cell_pivot_rho_Y];
                                }
                                Q[d_num_species][idx_cell_mom] = -Q[d_num_species][idx_cell_pivot_mom] +
                                    Real(2)*rho_pivot*d_bdry_face_adiabatic_no_slip_vel[face_loc*3];
                                Q[d_num_species + 1][idx_cell_mom] = -Q[d_num_species + 1][idx_cell_pivot_mom] +
                                    Real(2)*rho_pivot*d_bdry_face_adiabatic_no_slip_vel[face_loc*3 + 1];
                                Q[d_num_species + 2][idx_cell_mom] = -Q[d_num_species + 2][idx_cell_pivot_mom] +
                                    Real(2)*rho_pivot*d_bdry_face_adiabatic_no_slip_vel[face_loc*3 + 2];
                                
                                /*
                                 * Set the values for total internal energy.
                                 */
                                
                                Real epsilon_pivot = (Q[d_num_species + 3][idx_cell_pivot_E] -
                                    Real(1)/Real(2)*(Q[d_num_species][idx_cell_pivot_mom]*Q[d_num_species][idx_cell_pivot_mom] +
                                    Q[d_num_species + 1][idx_cell_pivot_mom]*Q[d_num_species + 1][idx_cell_pivot_mom] +
                                    Q[d_num_species + 2][idx_cell_pivot_mom]*Q[d_num_species + 2][idx_cell_pivot_mom])/
                                    rho_pivot)/rho_pivot;
                                
                                Real p_pivot = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_pivot,
                                        &epsilon_pivot,
                                        Y_pivot_ptr);
                                
                                Real T_pivot = d_equation_of_state_mixing_rules->
                                    getTemperature(
                                        &rho_pivot,
                                        &p_pivot,
                                        Y_pivot_ptr);
                                
                                Real T = T_pivot;
                                
                                Real epsilon = d_equation_of_state_mixing_rules->
                                    getInternalEnergyFromTemperature(
                                        &rho_pivot,
                                        &T,
                                        Y_pivot_ptr);
                                
                                Real E = rho_pivot*epsilon +
                                    Real(1)/Real(2)*(Q[d_num_species][idx_cell_mom]*Q[d_num_species][idx_cell_mom] +
                                    Q[d_num_species + 1][idx_cell_mom]*Q[d_num_species + 1][idx_cell_mom] +
                                    Q[d_num_species + 2][idx_cell_mom]*Q[d_num_species + 2][idx_cell_mom])/
                                    rho_pivot;
                                
                                Q[d_num_species + 3][idx_cell_E] = E;
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
                                const int idx_cell_rho_Y = (i + num_subghosts_conservative_var[0][0]) +
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
                                
                                int idx_cell_pivot_rho_Y = idx_cell_rho_Y;
                                int idx_cell_pivot_mom = idx_cell_mom;
                                int idx_cell_pivot_E = idx_cell_E;
                                
                                if (face_loc == BDRY_LOC::XLO)
                                {
                                    idx_cell_pivot_rho_Y = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
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
                                    idx_cell_pivot_rho_Y = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
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
                                    idx_cell_pivot_rho_Y = (i + num_subghosts_conservative_var[0][0]) +
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
                                    idx_cell_pivot_rho_Y = (i + num_subghosts_conservative_var[0][0]) +
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
                                    idx_cell_pivot_rho_Y = (i + num_subghosts_conservative_var[0][0]) +
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
                                    idx_cell_pivot_rho_Y = (i + num_subghosts_conservative_var[0][0]) +
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
                                 * Compute the mixture density of the pivot.
                                 */
                                
                                Real rho_pivot = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_pivot += Q[si][idx_cell_pivot_rho_Y];
                                }
                                
                                /*
                                 * Compute the mass fractions of the pivot.
                                 */
                                
                                std::vector<Real> Y_pivot;
                                Y_pivot.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_pivot.push_back(Q[si][idx_cell_pivot_rho_Y]/rho_pivot);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions of the pivot.
                                 */
                                
                                std::vector<const Real*> Y_pivot_ptr;
                                Y_pivot_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_pivot_ptr.push_back(&Y_pivot[si]);
                                }
                                
                                /*
                                 * Set the values for partial densities, momentum and total internal energy.
                                 */
                                
                                Real epsilon_pivot = (Q[d_num_species + 3][idx_cell_pivot_E] -
                                    Real(1)/Real(2)*(Q[d_num_species][idx_cell_pivot_mom]*Q[d_num_species][idx_cell_pivot_mom] +
                                    Q[d_num_species + 1][idx_cell_pivot_mom]*Q[d_num_species + 1][idx_cell_pivot_mom] +
                                    Q[d_num_species + 2][idx_cell_pivot_mom]*Q[d_num_species + 2][idx_cell_pivot_mom])/
                                    rho_pivot)/rho_pivot;
                                
                                Real p_pivot = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_pivot,
                                        &epsilon_pivot,
                                        Y_pivot_ptr);
                                
                                Real p = p_pivot;
                                
                                Real T_pivot = d_equation_of_state_mixing_rules->
                                    getTemperature(
                                        &rho_pivot,
                                        &p_pivot,
                                        Y_pivot_ptr);
                                
                                Real T = -T_pivot + Real(2)*d_bdry_face_isothermal_no_slip_T[face_loc];
                                
                                Real rho = d_equation_of_state_mixing_rules->
                                    getMixtureDensity(
                                        &p,
                                        &T,
                                        Y_pivot_ptr);
                                
                                Real u = -Q[d_num_species][idx_cell_pivot_mom]/rho_pivot +
                                    Real(2)*d_bdry_face_isothermal_no_slip_vel[face_loc*3];
                                Real v = -Q[d_num_species + 1][idx_cell_pivot_mom]/rho_pivot +
                                    Real(2)*d_bdry_face_isothermal_no_slip_vel[face_loc*3 + 1];
                                Real w = -Q[d_num_species + 2][idx_cell_pivot_mom]/rho_pivot +
                                    Real(2)*d_bdry_face_isothermal_no_slip_vel[face_loc*3 + 2];
                                
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Q[si][idx_cell_rho_Y] = rho*Y_pivot[si];
                                }
                                Q[d_num_species][idx_cell_mom] = rho*u;
                                Q[d_num_species + 1][idx_cell_mom] = rho*v;
                                Q[d_num_species + 2][idx_cell_mom] = rho*w;
                                
                                Real epsilon = d_equation_of_state_mixing_rules->
                                    getInternalEnergyFromTemperature(
                                        &rho,
                                        &T,
                                        Y_pivot_ptr);
                                
                                Real E = rho*epsilon +
                                    Real(1)/Real(2)*(Q[d_num_species][idx_cell_mom]*Q[d_num_species][idx_cell_mom] +
                                    Q[d_num_species + 1][idx_cell_mom]*Q[d_num_species + 1][idx_cell_mom] +
                                    Q[d_num_species + 2][idx_cell_mom]*Q[d_num_species + 2][idx_cell_mom])/rho;
                                
                                Q[d_num_species + 3][idx_cell_E] = E;
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
                
                if (face_loc == BDRY_LOC::XLO)
                {
                    const int num_ghosts_to_fill = fill_box_hi_idx[0] - fill_box_lo_idx[0] + 1;
                    TBOX_ASSERT(fill_box_hi_idx[0] == interior_box_lo_idx[0] - 1);
                    if (num_ghosts_to_fill > 6)
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelBoundaryUtilitiesFourEqnConservative::fill3dFaceBoundaryData()\n"
                            << "Non-reflecting outflow BC doesn't support more than six ghost cells yet!");
                    }
                    
                    for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                    {
                        for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                        {
                            // Get the grid spacing.
                            const double* const dx = patch_geom->getDx();
                            
                            const int idx_cell_rho_Y_x_R = (interior_box_lo_idx[0] + num_subghosts_conservative_var[0][0]) +
                                (j + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                (k + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                    subghostcell_dims_conservative_var[0][1];
                            
                            const int idx_cell_rho_Y_x_RR = (interior_box_lo_idx[0] + 1 + num_subghosts_conservative_var[0][0]) +
                                (j + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                (k + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                    subghostcell_dims_conservative_var[0][1];
                            
                            const int idx_cell_rho_Y_x_RRR = (interior_box_lo_idx[0] + 2 + num_subghosts_conservative_var[0][0]) +
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
                            
                            std::vector<Real> rho_Y_x_R;
                            std::vector<Real> rho_Y_x_RR;
                            std::vector<Real> rho_Y_x_RRR;
                            rho_Y_x_R.reserve(d_num_species);
                            rho_Y_x_RR.reserve(d_num_species);
                            rho_Y_x_RRR.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                rho_Y_x_R.push_back(Q[si][idx_cell_rho_Y_x_R]);
                                rho_Y_x_RR.push_back(Q[si][idx_cell_rho_Y_x_RR]);
                                rho_Y_x_RRR.push_back(Q[si][idx_cell_rho_Y_x_RRR]);
                            }
                            
                            /*
                             * Compute the mixture density.
                             */
                            
                            Real rho_x_R   = Real(0);
                            Real rho_x_RR  = Real(0);
                            Real rho_x_RRR = Real(0);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                rho_x_R   += Q[si][idx_cell_rho_Y_x_R];
                                rho_x_RR  += Q[si][idx_cell_rho_Y_x_RR];
                                rho_x_RRR += Q[si][idx_cell_rho_Y_x_RRR];
                            }
                            
                            /*
                             * Compute the mass fractions.
                             */
                            
                            std::vector<Real> Y_x_R;
                            std::vector<Real> Y_x_RR;
                            std::vector<Real> Y_x_RRR;
                            Y_x_R.reserve(d_num_species);
                            Y_x_RR.reserve(d_num_species);
                            Y_x_RRR.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_x_R.push_back(Q[si][idx_cell_rho_Y_x_R]/rho_x_R);
                                Y_x_RR.push_back(Q[si][idx_cell_rho_Y_x_RR]/rho_x_RR);
                                Y_x_RRR.push_back(Q[si][idx_cell_rho_Y_x_RRR]/rho_x_RRR);
                            }
                            
                            /*
                             * Get the pointers to the mass fractions.
                             */
                            
                            std::vector<const Real*> Y_x_R_ptr;
                            std::vector<const Real*> Y_x_RR_ptr;
                            std::vector<const Real*> Y_x_RRR_ptr;
                            Y_x_R_ptr.reserve(d_num_species);
                            Y_x_RR_ptr.reserve(d_num_species);
                            Y_x_RRR_ptr.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_x_R_ptr.push_back(&Y_x_R[si]);
                                Y_x_RR_ptr.push_back(&Y_x_RR[si]);
                                Y_x_RRR_ptr.push_back(&Y_x_RRR[si]);
                            }
                            
                            const Real u_x_R   = Q[d_num_species][idx_cell_mom_x_R]/rho_x_R;
                            const Real u_x_RR  = Q[d_num_species][idx_cell_mom_x_RR]/rho_x_RR;
                            const Real u_x_RRR = Q[d_num_species][idx_cell_mom_x_RRR]/rho_x_RRR;
                            
                            const Real v_x_R   = Q[d_num_species + 1][idx_cell_mom_x_R]/rho_x_R;
                            const Real v_x_RR  = Q[d_num_species + 1][idx_cell_mom_x_RR]/rho_x_RR;
                            const Real v_x_RRR = Q[d_num_species + 1][idx_cell_mom_x_RRR]/rho_x_RRR;
                            
                            const Real w_x_R   = Q[d_num_species + 2][idx_cell_mom_x_R]/rho_x_R;
                            const Real w_x_RR  = Q[d_num_species + 2][idx_cell_mom_x_RR]/rho_x_RR;
                            const Real w_x_RRR = Q[d_num_species + 2][idx_cell_mom_x_RRR]/rho_x_RRR;
                            
                            const Real half = Real(1)/Real(2);
                            const Real epsilon_x_R   = Q[d_num_species + 3][idx_cell_E_x_R]/rho_x_R - half*(u_x_R*u_x_R + v_x_R*v_x_R + w_x_R*w_x_R);
                            const Real epsilon_x_RR  = Q[d_num_species + 3][idx_cell_E_x_RR]/rho_x_RR - half*(u_x_RR*u_x_RR + v_x_RR*v_x_RR + w_x_RR*w_x_RR);
                            const Real epsilon_x_RRR = Q[d_num_species + 3][idx_cell_E_x_RRR]/rho_x_RRR - half*(u_x_RRR*u_x_RRR + v_x_RRR*v_x_RRR + w_x_RRR*w_x_RRR);
                            
                            Real p_x_R = d_equation_of_state_mixing_rules->
                                getPressure(
                                    &rho_x_R,
                                    &epsilon_x_R,
                                    Y_x_R_ptr);
                            
                            Real p_x_RR = d_equation_of_state_mixing_rules->
                                getPressure(
                                    &rho_x_RR,
                                    &epsilon_x_RR,
                                    Y_x_RR_ptr);
                            
                            Real p_x_RRR = d_equation_of_state_mixing_rules->
                                getPressure(
                                    &rho_x_RRR,
                                    &epsilon_x_RRR,
                                    Y_x_RRR_ptr);
                            
                            /*
                             * Compute derivatives in x-direction.
                             */
                            
                            std::vector<Real> drho_Y_dx;
                            drho_Y_dx.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                drho_Y_dx.push_back(-(Q[si][idx_cell_rho_Y_x_RRR] - Real(4)*Q[si][idx_cell_rho_Y_x_RR] +
                                    Real(3)*Q[si][idx_cell_rho_Y_x_R])/(Real(2)*Real(dx[0])));
                            }
                            const Real du_dx = -(u_x_RRR - Real(4)*u_x_RR + Real(3)*u_x_R)/(Real(2)*Real(dx[0]));
                            const Real dv_dx = -(v_x_RRR - Real(4)*v_x_RR + Real(3)*v_x_R)/(Real(2)*Real(dx[0]));
                            const Real dw_dx = -(w_x_RRR - Real(4)*w_x_RR + Real(3)*w_x_R)/(Real(2)*Real(dx[0]));
                            const Real dp_dx = -(p_x_RRR - Real(4)*p_x_RR + Real(3)*p_x_R)/(Real(2)*Real(dx[0]));
                            
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
                                
                                const int idx_cell_rho_Y_y_T = (interior_box_lo_idx[0] + num_subghosts_conservative_var[0][0]) +
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
                                
                                /*
                                 * Compute the mixture density.
                                 */
                                
                                Real rho_y_T = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_y_T += Q[si][idx_cell_rho_Y_y_T];
                                }
                                
                                /*
                                 * Compute the mass fractions.
                                 */
                                
                                std::vector<Real> Y_y_T;
                                Y_y_T.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_y_T.push_back(Q[si][idx_cell_rho_Y_y_T]/rho_y_T);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions.
                                 */
                            
                                std::vector<const Real*> Y_y_T_ptr;
                                Y_y_T_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_y_T_ptr.push_back(&Y_y_T[si]);
                                }
                                
                                const Real u_y_T = Q[d_num_species][idx_cell_mom_y_T]/rho_y_T;
                                const Real v_y_T = Q[d_num_species + 1][idx_cell_mom_y_T]/rho_y_T;
                                const Real w_y_T = Q[d_num_species + 2][idx_cell_mom_y_T]/rho_y_T;
                                const Real epsilon_y_T = Q[d_num_species + 3][idx_cell_E_y_T]/rho_y_T - half*(u_y_T*u_y_T + v_y_T*v_y_T + w_y_T*w_y_T);
                                
                                Real p_y_T = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_y_T,
                                        &epsilon_y_T,
                                        Y_y_T_ptr);
                                
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
                                const int idx_cell_rho_Y_y_B = (interior_box_lo_idx[0] + num_subghosts_conservative_var[0][0]) +
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
                                
                                /*
                                 * Compute the mixture density.
                                 */
                                
                                Real rho_y_B = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_y_B += Q[si][idx_cell_rho_Y_y_B];
                                }
                                
                                /*
                                 * Compute the mass fractions.
                                 */
                                
                                std::vector<Real> Y_y_B;
                                Y_y_B.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_y_B.push_back(Q[si][idx_cell_rho_Y_y_B]/rho_y_B);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions.
                                 */
                                
                                std::vector<const Real*> Y_y_B_ptr;
                                Y_y_B_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_y_B_ptr.push_back(&Y_y_B[si]);
                                }
                                
                                const Real u_y_B = Q[d_num_species][idx_cell_mom_y_B]/rho_y_B;
                                const Real v_y_B = Q[d_num_species + 1][idx_cell_mom_y_B]/rho_y_B;
                                const Real w_y_B = Q[d_num_species + 2][idx_cell_mom_y_B]/rho_y_B;
                                const Real epsilon_y_B = Q[d_num_species + 3][idx_cell_E_y_B]/rho_y_B - half*(u_y_B*u_y_B + v_y_B*v_y_B + w_y_B*w_y_B);
                                
                                Real p_y_B = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_y_B,
                                        &epsilon_y_B,
                                        Y_y_B_ptr);
                                
                                // One-sided derivatives.
                                du_dy = (u_x_R - u_y_B)/Real(dx[1]);
                                dv_dy = (v_x_R - v_y_B)/Real(dx[1]);
                                // dw_dy = (w_x_R - w_y_B)/Real(dx[1]);
                                dp_dy = (p_x_R - p_y_B)/Real(dx[1]);
                            }
                            else
                            {
                                const int idx_cell_rho_Y_y_B = (interior_box_lo_idx[0] + num_subghosts_conservative_var[0][0]) +
                                    (j - 1 + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (k + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_rho_Y_y_T = (interior_box_lo_idx[0] + num_subghosts_conservative_var[0][0]) +
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
                                
                                /*
                                 * Compute the mixture density.
                                 */
                                
                                Real rho_y_B = Real(0);
                                Real rho_y_T = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_y_B += Q[si][idx_cell_rho_Y_y_B];
                                    rho_y_T += Q[si][idx_cell_rho_Y_y_T];
                                }
                                
                                /*
                                 * Compute the mass fractions.
                                 */
                                
                                std::vector<Real> Y_y_B;
                                std::vector<Real> Y_y_T;
                                Y_y_B.reserve(d_num_species);
                                Y_y_T.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_y_B.push_back(Q[si][idx_cell_rho_Y_y_B]/rho_y_B);
                                    Y_y_T.push_back(Q[si][idx_cell_rho_Y_y_T]/rho_y_T);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions.
                                 */
                                
                                std::vector<const Real*> Y_y_B_ptr;
                                std::vector<const Real*> Y_y_T_ptr;
                                Y_y_B_ptr.reserve(d_num_species);
                                Y_y_T_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_y_B_ptr.push_back(&Y_y_B[si]);
                                    Y_y_T_ptr.push_back(&Y_y_T[si]);
                                }
                                
                                const Real u_y_B = Q[d_num_species][idx_cell_mom_y_B]/rho_y_B;
                                const Real u_y_T = Q[d_num_species][idx_cell_mom_y_T]/rho_y_T;
                                
                                const Real v_y_B = Q[d_num_species + 1][idx_cell_mom_y_B]/rho_y_B;
                                const Real v_y_T = Q[d_num_species + 1][idx_cell_mom_y_T]/rho_y_T;
                                
                                const Real w_y_B = Q[d_num_species + 2][idx_cell_mom_y_B]/rho_y_B;
                                const Real w_y_T = Q[d_num_species + 2][idx_cell_mom_y_T]/rho_y_T;
                                
                                const Real epsilon_y_B = Q[d_num_species + 3][idx_cell_E_y_B]/rho_y_B - half*(u_y_B*u_y_B + v_y_B*v_y_B + w_y_B*w_y_B);
                                const Real epsilon_y_T = Q[d_num_species + 3][idx_cell_E_y_T]/rho_y_T - half*(u_y_T*u_y_T + v_y_T*v_y_T + w_y_T*w_y_T);
                                
                                Real p_y_B = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_y_B,
                                        &epsilon_y_B,
                                        Y_y_B_ptr);
                            
                                Real p_y_T = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_y_T,
                                        &epsilon_y_T,
                                        Y_y_T_ptr);
                                
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
                                
                                const int idx_cell_rho_Y_z_F = (interior_box_lo_idx[0] + num_subghosts_conservative_var[0][0]) +
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
                                
                                /*
                                 * Compute the mixture density.
                                 */
                                
                                Real rho_z_F = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_z_F += Q[si][idx_cell_rho_Y_z_F];
                                }
                                
                                /*
                                 * Compute the mass fractions.
                                 */
                                
                                std::vector<Real> Y_z_F;
                                Y_z_F.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_z_F.push_back(Q[si][idx_cell_rho_Y_z_F]/rho_z_F);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions.
                                 */
                                
                                std::vector<const Real*> Y_z_F_ptr;
                                Y_z_F_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_z_F_ptr.push_back(&Y_z_F[si]);
                                }
                                
                                const Real u_z_F = Q[d_num_species][idx_cell_mom_z_F]/rho_z_F;
                                const Real v_z_F = Q[d_num_species + 1][idx_cell_mom_z_F]/rho_z_F;
                                const Real w_z_F = Q[d_num_species + 2][idx_cell_mom_z_F]/rho_z_F;
                                const Real epsilon_z_F = Q[d_num_species + 3][idx_cell_E_z_F]/rho_z_F - half*(u_z_F*u_z_F + v_z_F*v_z_F + w_z_F*w_z_F);
                                
                                Real p_z_F = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_z_F,
                                        &epsilon_z_F,
                                        Y_z_F_ptr);
                                
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
                                const int idx_cell_rho_Y_z_B = (interior_box_lo_idx[0] + num_subghosts_conservative_var[0][0]) +
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
                                
                                /*
                                 * Compute the mixture density.
                                 */
                                
                                Real rho_z_B = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_z_B += Q[si][idx_cell_rho_Y_z_B];
                                }
                                
                                /*
                                 * Compute the mass fractions.
                                 */
                                
                                std::vector<Real> Y_z_B;
                                Y_z_B.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_z_B.push_back(Q[si][idx_cell_rho_Y_z_B]/rho_z_B);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions.
                                 */
                                
                                std::vector<const Real*> Y_z_B_ptr;
                                Y_z_B_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_z_B_ptr.push_back(&Y_z_B[si]);
                                }
                                
                                const Real u_z_B = Q[d_num_species][idx_cell_mom_z_B]/rho_z_B;
                                const Real v_z_B = Q[d_num_species + 1][idx_cell_mom_z_B]/rho_z_B;
                                const Real w_z_B = Q[d_num_species + 2][idx_cell_mom_z_B]/rho_z_B;
                                const Real epsilon_z_B = Q[d_num_species + 3][idx_cell_E_z_B]/rho_z_B - half*(u_z_B*u_z_B + v_z_B*v_z_B + w_z_B*w_z_B);
                                
                                Real p_z_B = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_z_B,
                                        &epsilon_z_B,
                                        Y_z_B_ptr);
                                
                                // One-sided derivatives.
                                du_dz = (u_x_R - u_z_B)/Real(dx[2]);
                                // dv_dz = (v_x_R - v_z_B)/Real(dx[2]);
                                dw_dz = (w_x_R - w_z_B)/Real(dx[2]);
                                dp_dz = (p_x_R - p_z_B)/Real(dx[2]);
                            }
                            else
                            {
                                const int idx_cell_rho_Y_z_B = (interior_box_lo_idx[0] + num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (k - 1 + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_rho_Y_z_F = (interior_box_lo_idx[0] + num_subghosts_conservative_var[0][0]) +
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
                                
                                /*
                                 * Compute the mixture density.
                                 */
                                
                                Real rho_z_B = Real(0);
                                Real rho_z_F = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_z_B += Q[si][idx_cell_rho_Y_z_B];
                                    rho_z_F += Q[si][idx_cell_rho_Y_z_F];
                                }
                                
                                /*
                                 * Compute the mass fractions.
                                 */
                                
                                std::vector<Real> Y_z_B;
                                std::vector<Real> Y_z_F;
                                Y_z_B.reserve(d_num_species);
                                Y_z_F.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_z_B.push_back(Q[si][idx_cell_rho_Y_z_B]/rho_z_B);
                                    Y_z_F.push_back(Q[si][idx_cell_rho_Y_z_F]/rho_z_F);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions.
                                 */
                                
                                std::vector<const Real*> Y_z_B_ptr;
                                std::vector<const Real*> Y_z_F_ptr;
                                Y_z_B_ptr.reserve(d_num_species);
                                Y_z_F_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_z_B_ptr.push_back(&Y_z_B[si]);
                                    Y_z_F_ptr.push_back(&Y_z_F[si]);
                                }
                                
                                const Real u_z_B = Q[d_num_species][idx_cell_mom_z_B]/rho_z_B;
                                const Real u_z_F = Q[d_num_species][idx_cell_mom_z_F]/rho_z_F;
                                
                                const Real v_z_B = Q[d_num_species + 1][idx_cell_mom_z_B]/rho_z_B;
                                const Real v_z_F = Q[d_num_species + 1][idx_cell_mom_z_F]/rho_z_F;
                                
                                const Real w_z_B = Q[d_num_species + 2][idx_cell_mom_z_B]/rho_z_B;
                                const Real w_z_F = Q[d_num_species + 2][idx_cell_mom_z_F]/rho_z_F;
                                
                                const Real epsilon_z_B = Q[d_num_species + 3][idx_cell_E_z_B]/rho_z_B - half*(u_z_B*u_z_B + v_z_B*v_z_B + w_z_B*w_z_B);
                                const Real epsilon_z_F = Q[d_num_species + 3][idx_cell_E_z_F]/rho_z_F - half*(u_z_F*u_z_F + v_z_F*v_z_F + w_z_F*w_z_F);
                                
                                Real p_z_B = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_z_B,
                                        &epsilon_z_B,
                                        Y_z_B_ptr);
                                
                                Real p_z_F = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_z_F,
                                        &epsilon_z_F,
                                        Y_z_F_ptr);
                                
                                // Central derivatives.
                                du_dz = (u_z_F - u_z_B)/(Real(2)*Real(dx[2]));
                                // dv_dz = (v_z_F - v_z_B)/(Real(2)*Real(dx[2]));
                                dw_dz = (w_z_F - w_z_B)/(Real(2)*Real(dx[2]));
                                dp_dz = (p_z_F - p_z_B)/(Real(2)*Real(dx[2]));
                            }
                            
                            // Compute sound speed.
                            
                            const Real Gamma_x_R = d_equation_of_state_mixing_rules->getGruneisenParameter(
                                &rho_x_R,
                                &p_x_R,
                                Y_x_R_ptr);
                            
                            const std::vector<Real> Psi_x_R = d_equation_of_state_mixing_rules->
                                getPressureDerivativeWithPartialDensities(
                                        &rho_x_R,
                                        &p_x_R,
                                        Y_x_R_ptr);
                            
                            Real c_x_R = Gamma_x_R*p_x_R/rho_x_R;
                            for (int si = 0; si < d_num_species; si++)
                            {
                                c_x_R += Y_x_R[si]*Psi_x_R[si];
                            }
                            c_x_R = std::sqrt(c_x_R);
                            
                            const Real lambda_last = u_x_R + c_x_R;
                            
                            // Compute vector Lambda^(-1) * L.
                            
                            Real Lambda_inv_L[d_num_species + 4];
                            
                            const Real& p_t         = d_bdry_face_nonreflecting_outflow_p_t[face_loc];
                            const Real& sigma       = d_bdry_face_nonreflecting_outflow_sigma[face_loc];
                            const Real& beta        = d_bdry_face_nonreflecting_outflow_beta[face_loc];
                            const Real& length_char = d_bdry_face_nonreflecting_outflow_length_char[face_loc];
                            
                            const Real T_last = v_x_R*(dp_dy + rho_x_R*c_x_R*du_dy) + rho_x_R*c_x_R*c_x_R*dv_dy + 
                            w_x_R*(dp_dz + rho_x_R*c_x_R*du_dz) + rho_x_R*c_x_R*c_x_R*dw_dz;
                            
                            const Real M_sq = (u_x_R*u_x_R + v_x_R*v_x_R + w_x_R*w_x_R)/(c_x_R*c_x_R);
                            const Real K = sigma*c_x_R*(Real(1) - M_sq)/length_char;
                            
                            Lambda_inv_L[0] = dp_dx - rho_x_R*c_x_R*du_dx;
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Lambda_inv_L[si + 1] = c_x_R*c_x_R*drho_Y_dx[si] - Y_x_R[si]*dp_dx;
                            }
                            Lambda_inv_L[d_num_species + 1] = dv_dx;
                            Lambda_inv_L[d_num_species + 2] = dw_dx;
                            Lambda_inv_L[d_num_species + 3] = (Real(1)/lambda_last)*(K*(p_x_R - p_t) - (Real(1) - beta)*T_last);
                            
                            // Compute dV_dx.
                            
                            const Real c_sq_inv  = Real(1)/(c_x_R*c_x_R);
                            const Real rho_c_inv = Real(1)/(rho_x_R*c_x_R);
                            
                            Real dV_dx[d_num_species + 4];
                            
                            for (int si = 0; si < d_num_species; si++)
                            {
                                dV_dx[si] = half*c_sq_inv*Y_x_R[si]*(Lambda_inv_L[0] + Lambda_inv_L[d_num_species + 3]) +
                                    c_sq_inv*Lambda_inv_L[si + 1];
                            }
                            dV_dx[d_num_species]     = half*rho_c_inv*(-Lambda_inv_L[0] + Lambda_inv_L[d_num_species + 3]);
                            dV_dx[d_num_species + 1] = Lambda_inv_L[d_num_species + 1];
                            dV_dx[d_num_species + 2] = Lambda_inv_L[d_num_species + 2];
                            dV_dx[d_num_species + 3] = half*(Lambda_inv_L[0] + Lambda_inv_L[d_num_species + 3]);
                            
                            Real V_ghost[(d_num_species + 4)*num_ghosts_to_fill];
                            
                            for (int i = num_ghosts_to_fill - 1; i >= 0; i--)
                            {
                                const int idx_cell_rho_Y = (i + fill_box_lo_idx[0] + num_subghosts_conservative_var[0][0]) +
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
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        V_ghost[i*(d_num_species + 4) + si] = rho_Y_x_RR[si] - Real(2)*Real(dx[0])*dV_dx[si];
                                    }
                                    
                                    V_ghost[i*(d_num_species + 4) + d_num_species]     = u_x_RR - Real(2)*Real(dx[0])*dV_dx[d_num_species];
                                    V_ghost[i*(d_num_species + 4) + d_num_species + 1] = v_x_RR - Real(2)*Real(dx[0])*dV_dx[d_num_species + 1];
                                    V_ghost[i*(d_num_species + 4) + d_num_species + 2] = w_x_RR - Real(2)*Real(dx[0])*dV_dx[d_num_species + 2];
                                    V_ghost[i*(d_num_species + 4) + d_num_species + 3] = p_x_RR - Real(2)*Real(dx[0])*dV_dx[d_num_species + 3];
                                }
                                else if (i == num_ghosts_to_fill - 2)
                                {
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        V_ghost[i*(d_num_species + 4) + si] = -Real(2)*rho_Y_x_RR[si] - Real(3)*rho_Y_x_R[si] +
                                            Real(6)*V_ghost[(i + 1)*(d_num_species + 4) + si] + Real(6)*Real(dx[0])*dV_dx[si];
                                    }
                                    
                                    V_ghost[i*(d_num_species + 4) + d_num_species] = -Real(2)*u_x_RR - Real(3)*u_x_R +
                                        Real(6)*V_ghost[(i + 1)*(d_num_species + 4) + d_num_species] +
                                        Real(6)*Real(dx[0])*dV_dx[d_num_species];
                                    
                                    V_ghost[i*(d_num_species + 4) + d_num_species + 1] = -Real(2)*v_x_RR - Real(3)*v_x_R +
                                        Real(6)*V_ghost[(i + 1)*(d_num_species + 4) + d_num_species + 1] +
                                        Real(6)*Real(dx[0])*dV_dx[d_num_species + 1];
                                    
                                    V_ghost[i*(d_num_species + 4) + d_num_species + 2] = -Real(2)*w_x_RR - Real(3)*w_x_R +
                                        Real(6)*V_ghost[(i + 1)*(d_num_species + 4) + d_num_species + 2] +
                                        Real(6)*Real(dx[0])*dV_dx[d_num_species + 2];
                                    
                                    V_ghost[i*(d_num_species + 4) + d_num_species + 3] = -Real(2)*p_x_RR - Real(3)*p_x_R +
                                        Real(6)*V_ghost[(i + 1)*(d_num_species + 4) + d_num_species + 3] +
                                        Real(6)*Real(dx[0])*dV_dx[d_num_species + 3];
                                }
                                else if (i == num_ghosts_to_fill - 3)
                                {
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        V_ghost[i*(d_num_species + 4) + si] = Real(3)*rho_Y_x_RR[si] + Real(10)*rho_Y_x_R[si] -
                                            Real(18)*V_ghost[(i + 2)*(d_num_species + 4) + si] +
                                            Real(6)*V_ghost[(i + 1)*(d_num_species + 4) + si] -
                                            Real(12)*Real(dx[0])*dV_dx[si];
                                    }
                                    
                                    V_ghost[i*(d_num_species + 4) + d_num_species] = Real(3)*u_x_RR + Real(10)*u_x_R -
                                        Real(18)*V_ghost[(i + 2)*(d_num_species + 4) + d_num_species] +
                                        Real(6)*V_ghost[(i + 1)*(d_num_species + 4) + d_num_species] -
                                        Real(12)*Real(dx[0])*dV_dx[d_num_species];
                                    
                                    V_ghost[i*(d_num_species + 4) + d_num_species + 1] = Real(3)*v_x_RR + Real(10)*v_x_R -
                                        Real(18)*V_ghost[(i + 2)*(d_num_species + 4) + d_num_species + 1] +
                                        Real(6)*V_ghost[(i + 1)*(d_num_species + 4) + d_num_species + 1] -
                                        Real(12)*Real(dx[0])*dV_dx[d_num_species + 1];
                                    
                                    V_ghost[i*(d_num_species + 4) + d_num_species + 2] = Real(3)*w_x_RR + Real(10)*w_x_R -
                                        Real(18)*V_ghost[(i + 2)*(d_num_species + 4) + d_num_species + 2] +
                                        Real(6)*V_ghost[(i + 1)*(d_num_species + 4) + d_num_species + 2] -
                                        Real(12)*Real(dx[0])*dV_dx[d_num_species + 2];
                                    
                                    V_ghost[i*(d_num_species + 4) + d_num_species + 3] = Real(3)*p_x_RR + Real(10)*p_x_R -
                                        Real(18)*V_ghost[(i + 2)*(d_num_species + 4) + d_num_species + 3] +
                                        Real(6)*V_ghost[(i + 1)*(d_num_species + 4) + d_num_species + 3] -
                                        Real(12)*Real(dx[0])*dV_dx[d_num_species + 3];
                                }
                                else if (i == num_ghosts_to_fill - 4)
                                {
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        V_ghost[i*(d_num_species + 4) + si] = -Real(4)*rho_Y_x_RR[si] -
                                            Real(65)/Real(3)*rho_Y_x_R[si] +
                                            Real(40)*V_ghost[(i + 3)*(d_num_species + 4) + si] -
                                            Real(20)*V_ghost[(i + 2)*(d_num_species + 4) + si] +
                                            Real(20)/Real(3)*V_ghost[(i + 1)*(d_num_species + 4) + si] +
                                            Real(20)*Real(dx[0])*dV_dx[si];
                                    }
                                    
                                    V_ghost[i*(d_num_species + 4) + d_num_species] = -Real(4)*u_x_RR -
                                        Real(65)/Real(3)*u_x_R +
                                        Real(40)*V_ghost[(i + 3)*(d_num_species + 4) + d_num_species] -
                                        Real(20)*V_ghost[(i + 2)*(d_num_species + 4) + d_num_species] +
                                        Real(20)/Real(3)*V_ghost[(i + 1)*(d_num_species + 4) + d_num_species] +
                                        Real(20)*Real(dx[0])*dV_dx[d_num_species];
                                    
                                    V_ghost[i*(d_num_species + 4) + d_num_species + 1] = -Real(4)*v_x_RR -
                                        Real(65)/Real(3)*v_x_R +
                                        Real(40)*V_ghost[(i + 3)*(d_num_species + 4) + d_num_species + 1] -
                                        Real(20)*V_ghost[(i + 2)*(d_num_species + 4) + d_num_species + 1] +
                                        Real(20)/Real(3)*V_ghost[(i + 1)*(d_num_species + 4) + d_num_species + 1] +
                                        Real(20)*Real(dx[0])*dV_dx[d_num_species + 1];
                                    
                                    V_ghost[i*(d_num_species + 4) + d_num_species + 2] = -Real(4)*w_x_RR -
                                        Real(65)/Real(3)*w_x_R +
                                        Real(40)*V_ghost[(i + 3)*(d_num_species + 4) + d_num_species + 2] -
                                        Real(20)*V_ghost[(i + 2)*(d_num_species + 4) + d_num_species + 2] +
                                        Real(20)/Real(3)*V_ghost[(i + 1)*(d_num_species + 4) + d_num_species + 2] +
                                        Real(20)*Real(dx[0])*dV_dx[d_num_species + 2];
                                    
                                    V_ghost[i*(d_num_species + 4) + d_num_species + 3] = -Real(4)*p_x_RR -
                                        Real(65)/Real(3)*p_x_R +
                                        Real(40)*V_ghost[(i + 3)*(d_num_species + 4) + d_num_species + 3] -
                                        Real(20)*V_ghost[(i + 2)*(d_num_species + 4) + d_num_species + 3] +
                                        Real(20)/Real(3)*V_ghost[(i + 1)*(d_num_species + 4) + d_num_species + 3] +
                                        Real(20)*Real(dx[0])*dV_dx[d_num_species + 3];
                                }
                                else if (i == num_ghosts_to_fill - 5)
                                {
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        V_ghost[i*(d_num_species + 4) + si] = Real(5)*rho_Y_x_RR[si] +
                                            Real(77)/Real(2)*rho_Y_x_R[si] -
                                            Real(75)*V_ghost[(i + 4)*(d_num_species + 4) + si] +
                                            Real(50)*V_ghost[(i + 3)*(d_num_species + 4) + si] -
                                            Real(25)*V_ghost[(i + 2)*(d_num_species + 4) + si] +
                                            Real(15)/Real(2)*V_ghost[(i + 1)*(d_num_species + 4) + si] -
                                            Real(30)*Real(dx[0])*dV_dx[si];
                                    }
                                    
                                    V_ghost[i*(d_num_species + 4) + d_num_species] = Real(5)*u_x_RR +
                                        Real(77)/Real(2)*u_x_R -
                                        Real(75)*V_ghost[(i + 4)*(d_num_species + 4) + d_num_species] +
                                        Real(50)*V_ghost[(i + 3)*(d_num_species + 4) + d_num_species] -
                                        Real(25)*V_ghost[(i + 2)*(d_num_species + 4) + d_num_species] +
                                        Real(15)/Real(2)*V_ghost[(i + 1)*(d_num_species + 4) + d_num_species] -
                                        Real(30)*Real(dx[0])*dV_dx[d_num_species];
                                    
                                    V_ghost[i*(d_num_species + 4) + d_num_species + 1] = Real(5)*v_x_RR +
                                        Real(77)/Real(2)*v_x_R -
                                        Real(75)*V_ghost[(i + 4)*(d_num_species + 4) + d_num_species + 1] +
                                        Real(50)*V_ghost[(i + 3)*(d_num_species + 4) + d_num_species + 1] -
                                        Real(25)*V_ghost[(i + 2)*(d_num_species + 4) + d_num_species + 1] +
                                        Real(15)/Real(2)*V_ghost[(i + 1)*(d_num_species + 4) + d_num_species + 1] -
                                        Real(30)*Real(dx[0])*dV_dx[d_num_species + 1];
                                    
                                    V_ghost[i*(d_num_species + 4) + d_num_species + 2] = Real(5)*w_x_RR +
                                        Real(77)/Real(2)*w_x_R -
                                        Real(75)*V_ghost[(i + 4)*(d_num_species + 4) + d_num_species + 2] +
                                        Real(50)*V_ghost[(i + 3)*(d_num_species + 4) + d_num_species + 2] -
                                        Real(25)*V_ghost[(i + 2)*(d_num_species + 4) + d_num_species + 2] +
                                        Real(15)/Real(2)*V_ghost[(i + 1)*(d_num_species + 4) + d_num_species + 2] -
                                        Real(30)*Real(dx[0])*dV_dx[d_num_species + 2];
                                    
                                    V_ghost[i*(d_num_species + 4) + d_num_species + 3] = Real(5)*p_x_RR +
                                        Real(77)/Real(2)*p_x_R -
                                        Real(75)*V_ghost[(i + 4)*(d_num_species + 4) + d_num_species + 3] +
                                        Real(50)*V_ghost[(i + 3)*(d_num_species + 4) + d_num_species + 3] -
                                        Real(25)*V_ghost[(i + 2)*(d_num_species + 4) + d_num_species + 3] +
                                        Real(15)/Real(2)*V_ghost[(i + 1)*(d_num_species + 4) + d_num_species + 3] -
                                        Real(30)*Real(dx[0])*dV_dx[d_num_species + 3];
                                }
                                else if (i == num_ghosts_to_fill - 6)
                                {
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        V_ghost[i*(d_num_species + 4) + si] = -Real(6)*rho_Y_x_RR[si] -
                                            Real(609)/Real(10)*rho_Y_x_R[si] +
                                            Real(126)*V_ghost[(i + 5)*(d_num_species + 4) + si] -
                                            Real(105)*V_ghost[(i + 4)*(d_num_species + 4) + si] +
                                            Real(70)*V_ghost[(i + 3)*(d_num_species + 4) + si] -
                                            Real(63)/Real(2)*V_ghost[(i + 2)*(d_num_species + 4) + si] +
                                            Real(42)/Real(5)*V_ghost[(i + 1)*(d_num_species + 4) + si] +
                                            Real(42)*Real(dx[0])*dV_dx[si];
                                    }
                                    
                                    V_ghost[i*(d_num_species + 4) + d_num_species] = -Real(6)*u_x_RR -
                                        Real(609)/Real(10)*u_x_R +
                                        Real(126)*V_ghost[(i + 5)*(d_num_species + 4) + d_num_species] -
                                        Real(105)*V_ghost[(i + 4)*(d_num_species + 4) + d_num_species] +
                                        Real(70)*V_ghost[(i + 3)*(d_num_species + 4) + d_num_species] -
                                        Real(63)/Real(2)*V_ghost[(i + 2)*(d_num_species + 4) + d_num_species] +
                                        Real(42)/Real(5)*V_ghost[(i + 1)*(d_num_species + 4) + d_num_species] +
                                        Real(42)*Real(dx[0])*dV_dx[d_num_species];
                                    
                                    V_ghost[i*(d_num_species + 4) + d_num_species + 1] = -Real(6)*v_x_RR -
                                        Real(609)/Real(10)*v_x_R +
                                        Real(126)*V_ghost[(i + 5)*(d_num_species + 4) + d_num_species + 1] -
                                        Real(105)*V_ghost[(i + 4)*(d_num_species + 4) + d_num_species + 1] +
                                        Real(70)*V_ghost[(i + 3)*(d_num_species + 4) + d_num_species + 1] -
                                        Real(63)/Real(2)*V_ghost[(i + 2)*(d_num_species + 4) + d_num_species + 1] +
                                        Real(42)/Real(5)*V_ghost[(i + 1)*(d_num_species + 4) + d_num_species + 1] +
                                        Real(42)*Real(dx[0])*dV_dx[d_num_species + 1];
                                    
                                    V_ghost[i*(d_num_species + 4) + d_num_species + 2] = -Real(6)*w_x_RR -
                                        Real(609)/Real(10)*w_x_R +
                                        Real(126)*V_ghost[(i + 5)*(d_num_species + 4) + d_num_species + 2] -
                                        Real(105)*V_ghost[(i + 4)*(d_num_species + 4) + d_num_species + 2] +
                                        Real(70)*V_ghost[(i + 3)*(d_num_species + 4) + d_num_species + 2] -
                                        Real(63)/Real(2)*V_ghost[(i + 2)*(d_num_species + 4) + d_num_species + 2] +
                                        Real(42)/Real(5)*V_ghost[(i + 1)*(d_num_species + 4) + d_num_species + 2] +
                                        Real(42)*Real(dx[0])*dV_dx[d_num_species + 2];
                                    
                                    V_ghost[i*(d_num_species + 4) + d_num_species + 3] = -Real(6)*p_x_RR -
                                        Real(609)/Real(10)*p_x_R +
                                        Real(126)*V_ghost[(i + 5)*(d_num_species + 4) + d_num_species + 3] -
                                        Real(105)*V_ghost[(i + 4)*(d_num_species + 4) + d_num_species + 3] +
                                        Real(70)*V_ghost[(i + 3)*(d_num_species + 4) + d_num_species + 3] -
                                        Real(63)/Real(2)*V_ghost[(i + 2)*(d_num_species + 4) + d_num_species + 3] +
                                        Real(42)/Real(5)*V_ghost[(i + 1)*(d_num_species + 4) + d_num_species + 3] +
                                        Real(42)*Real(dx[0])*dV_dx[d_num_species + 3];
                                }
                                
                                /*
                                 * Compute the mixture density.
                                 */
                                
                                Real rho_ghost = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_ghost += V_ghost[i*(d_num_species + 4) + si];
                                }
                                
                                /*
                                 * Compute the mass fractions.
                                 */
                                
                                std::vector<Real> Y_ghost;
                                Y_ghost.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_ghost.push_back(V_ghost[i*(d_num_species + 4) + si]/rho_ghost);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions.
                                 */
                                
                                std::vector<const Real*> Y_ghost_ptr;
                                Y_ghost_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_ghost_ptr.push_back(&Y_ghost[si]);
                                }
                                
                                for(int si=0; si < d_num_species; si++)
                                {
                                    Q[si][idx_cell_rho_Y] = V_ghost[i*(d_num_species + 4) + si];
                                }
                                
                                Q[d_num_species][idx_cell_mom]     = rho_ghost*V_ghost[i*(d_num_species + 4) + d_num_species];
                                Q[d_num_species + 1][idx_cell_mom] = rho_ghost*V_ghost[i*(d_num_species + 4) + d_num_species + 1];
                                Q[d_num_species + 2][idx_cell_mom] = rho_ghost*V_ghost[i*(d_num_species + 4) + d_num_species + 2];
                                
                                const Real epsilon = d_equation_of_state_mixing_rules->
                                    getInternalEnergy(
                                        &rho_ghost,
                                        &V_ghost[i*(d_num_species + 4) + d_num_species + 3],
                                        Y_ghost_ptr);
                                
                                const Real E = rho_ghost*epsilon +
                                    half*(Q[d_num_species][idx_cell_mom]*Q[d_num_species][idx_cell_mom] +
                                        Q[d_num_species + 1][idx_cell_mom]*Q[d_num_species + 1][idx_cell_mom] +
                                        Q[d_num_species + 2][idx_cell_mom]*Q[d_num_species + 2][idx_cell_mom])/
                                        rho_ghost;
                                
                                Q[d_num_species + 3][idx_cell_E] = E;
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
                            << ": FlowModelBoundaryUtilitiesFourEqnConservative::fill3dFaceBoundaryData()\n"
                            << "Non-reflecting outflow BC doesn't support more than six ghost cells yet!");
                    }
                    
                    for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                    {
                        for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                        {
                            // Get the grid spacing.
                            const double* const dx = patch_geom->getDx();
                            
                            const int idx_cell_rho_Y_x_L = (interior_box_hi_idx[0] + num_subghosts_conservative_var[0][0]) +
                                (j + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                (k + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                    subghostcell_dims_conservative_var[0][1];
                            
                            const int idx_cell_rho_Y_x_LL = (interior_box_hi_idx[0] - 1 + num_subghosts_conservative_var[0][0]) +
                                (j + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                (k + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                    subghostcell_dims_conservative_var[0][1];
                            
                            const int idx_cell_rho_Y_x_LLL = (interior_box_hi_idx[0] - 2 + num_subghosts_conservative_var[0][0]) +
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
                            
                            std::vector<Real> rho_Y_x_L;
                            std::vector<Real> rho_Y_x_LL;
                            std::vector<Real> rho_Y_x_LLL;
                            rho_Y_x_L.reserve(d_num_species);
                            rho_Y_x_LL.reserve(d_num_species);
                            rho_Y_x_LLL.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                rho_Y_x_L.push_back(Q[si][idx_cell_rho_Y_x_L]);
                                rho_Y_x_LL.push_back(Q[si][idx_cell_rho_Y_x_LL]);
                                rho_Y_x_LLL.push_back(Q[si][idx_cell_rho_Y_x_LLL]);
                            }
                            
                            /*
                             * Compute the mixture density.
                             */
                            
                            Real rho_x_L   = Real(0);
                            Real rho_x_LL  = Real(0);
                            Real rho_x_LLL = Real(0);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                rho_x_L   += Q[si][idx_cell_rho_Y_x_L];
                                rho_x_LL  += Q[si][idx_cell_rho_Y_x_LL];
                                rho_x_LLL += Q[si][idx_cell_rho_Y_x_LLL];
                            }
                            
                            /*
                             * Compute the mass fractions.
                             */
                            
                            std::vector<Real> Y_x_L;
                            std::vector<Real> Y_x_LL;
                            std::vector<Real> Y_x_LLL;
                            Y_x_L.reserve(d_num_species);
                            Y_x_LL.reserve(d_num_species);
                            Y_x_LLL.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_x_L.push_back(Q[si][idx_cell_rho_Y_x_L]/rho_x_L);
                                Y_x_LL.push_back(Q[si][idx_cell_rho_Y_x_LL]/rho_x_LL);
                                Y_x_LLL.push_back(Q[si][idx_cell_rho_Y_x_LLL]/rho_x_LLL);
                            }
                            
                            /*
                             * Get the pointers to the mass fractions.
                             */
                            
                            std::vector<const Real*> Y_x_L_ptr;
                            std::vector<const Real*> Y_x_LL_ptr;
                            std::vector<const Real*> Y_x_LLL_ptr;
                            Y_x_L_ptr.reserve(d_num_species);
                            Y_x_LL_ptr.reserve(d_num_species);
                            Y_x_LLL_ptr.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_x_L_ptr.push_back(&Y_x_L[si]);
                                Y_x_LL_ptr.push_back(&Y_x_LL[si]);
                                Y_x_LLL_ptr.push_back(&Y_x_LLL[si]);
                            }
                            
                            const Real u_x_L   = Q[d_num_species][idx_cell_mom_x_L]/rho_x_L;
                            const Real u_x_LL  = Q[d_num_species][idx_cell_mom_x_LL]/rho_x_LL;
                            const Real u_x_LLL = Q[d_num_species][idx_cell_mom_x_LLL]/rho_x_LLL;
                            
                            const Real v_x_L   = Q[d_num_species + 1][idx_cell_mom_x_L]/rho_x_L;
                            const Real v_x_LL  = Q[d_num_species + 1][idx_cell_mom_x_LL]/rho_x_LL;
                            const Real v_x_LLL = Q[d_num_species + 1][idx_cell_mom_x_LLL]/rho_x_LLL;
                            
                            const Real w_x_L   = Q[d_num_species + 2][idx_cell_mom_x_L]/rho_x_L;
                            const Real w_x_LL  = Q[d_num_species + 2][idx_cell_mom_x_LL]/rho_x_LL;
                            const Real w_x_LLL = Q[d_num_species + 2][idx_cell_mom_x_LLL]/rho_x_LLL;
                            
                            const Real half = Real(1)/Real(2);
                            const Real epsilon_x_L   = Q[d_num_species + 3][idx_cell_E_x_L]/rho_x_L - half*(u_x_L*u_x_L + v_x_L*v_x_L + w_x_L*w_x_L);
                            const Real epsilon_x_LL  = Q[d_num_species + 3][idx_cell_E_x_LL]/rho_x_LL - half*(u_x_LL*u_x_LL + v_x_LL*v_x_LL + w_x_LL*w_x_LL);
                            const Real epsilon_x_LLL = Q[d_num_species + 3][idx_cell_E_x_LLL]/rho_x_LLL - half*(u_x_LLL*u_x_LLL + v_x_LLL*v_x_LLL + w_x_LLL*w_x_LLL);

                            Real p_x_L = d_equation_of_state_mixing_rules->
                                getPressure(
                                    &rho_x_L,
                                    &epsilon_x_L,
                                    Y_x_L_ptr);
                            
                            Real p_x_LL = d_equation_of_state_mixing_rules->
                                getPressure(
                                    &rho_x_LL,
                                    &epsilon_x_LL,
                                    Y_x_LL_ptr);
                            
                            Real p_x_LLL = d_equation_of_state_mixing_rules->
                                getPressure(
                                    &rho_x_LLL,
                                    &epsilon_x_LLL,
                                    Y_x_LLL_ptr);
                            
                            /*
                             * Compute derivatives in x-direction.
                             */
                            
                            std::vector<Real> drho_Y_dx;
                            drho_Y_dx.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                drho_Y_dx.push_back((Q[si][idx_cell_rho_Y_x_LLL] - Real(4)*Q[si][idx_cell_rho_Y_x_LL] +
                                    Real(3)*Q[si][idx_cell_rho_Y_x_L])/(Real(2)*Real(dx[0])));
                            }
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
                                
                                const int idx_cell_rho_Y_y_T = (interior_box_hi_idx[0] + num_subghosts_conservative_var[0][0]) +
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
                                
                                /*
                                 * Compute the mixture density.
                                 */
                                
                                Real rho_y_T = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_y_T += Q[si][idx_cell_rho_Y_y_T];
                                }
                                
                                /*
                                 * Compute the mass fractions.
                                 */
                                
                                std::vector<Real> Y_y_T;
                                Y_y_T.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_y_T.push_back(Q[si][idx_cell_rho_Y_y_T]/rho_y_T);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions.
                                 */
                            
                                std::vector<const Real*> Y_y_T_ptr;
                                Y_y_T_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_y_T_ptr.push_back(&Y_y_T[si]);
                                }
                                
                                const Real u_y_T = Q[d_num_species][idx_cell_mom_y_T]/rho_y_T;
                                const Real v_y_T = Q[d_num_species + 1][idx_cell_mom_y_T]/rho_y_T;
                                const Real w_y_T = Q[d_num_species + 2][idx_cell_mom_y_T]/rho_y_T;
                                const Real epsilon_y_T = Q[d_num_species + 3][idx_cell_E_y_T]/rho_y_T - half*(u_y_T*u_y_T + v_y_T*v_y_T + w_y_T*w_y_T);
                                
                                Real p_y_T = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_y_T,
                                        &epsilon_y_T,
                                        Y_y_T_ptr);
                                
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
                                // Patch is touching periodIic boundary.
                                
                                const int idx_cell_rho_Y_y_B = (interior_box_hi_idx[0] + num_subghosts_conservative_var[0][0]) +
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
                                
                                /*
                                 * Compute the mixture density.
                                 */
                                
                                Real rho_y_B = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_y_B += Q[si][idx_cell_rho_Y_y_B];
                                }
                                
                                /*
                                 * Compute the mass fractions.
                                 */
                                
                                std::vector<Real> Y_y_B;
                                Y_y_B.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_y_B.push_back(Q[si][idx_cell_rho_Y_y_B]/rho_y_B);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions.
                                 */
                                
                                std::vector<const Real*> Y_y_B_ptr;
                                Y_y_B_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_y_B_ptr.push_back(&Y_y_B[si]);
                                }
                                
                                const Real u_y_B = Q[d_num_species][idx_cell_mom_y_B]/rho_y_B;
                                const Real v_y_B = Q[d_num_species + 1][idx_cell_mom_y_B]/rho_y_B;
                                const Real w_y_B = Q[d_num_species + 2][idx_cell_mom_y_B]/rho_y_B;
                                const Real epsilon_y_B = Q[d_num_species + 3][idx_cell_E_y_B]/rho_y_B - half*(u_y_B*u_y_B + v_y_B*v_y_B + w_y_B*w_y_B);
                                
                                Real p_y_B = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_y_B,
                                        &epsilon_y_B,
                                        Y_y_B_ptr);
                                
                                // One-sided derivatives.
                                du_dy = (u_x_L - u_y_B)/Real(dx[1]);
                                dv_dy = (v_x_L - v_y_B)/Real(dx[1]);
                                // dw_dy = (w_x_L - w_y_B)/Real(dx[1]);
                                dp_dy = (p_x_L - p_y_B)/Real(dx[1]);
                            }
                            else
                            {
                                const int idx_cell_rho_Y_y_B = (interior_box_hi_idx[0] + num_subghosts_conservative_var[0][0]) +
                                    (j - 1 + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (k + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_rho_Y_y_T = (interior_box_hi_idx[0] + num_subghosts_conservative_var[0][0]) +
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
                                
                                /*
                                 * Compute the mixture density.
                                 */
                                
                                Real rho_y_B = Real(0);
                                Real rho_y_T = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_y_B += Q[si][idx_cell_rho_Y_y_B];
                                    rho_y_T += Q[si][idx_cell_rho_Y_y_T];
                                }
                                
                                /*
                                 * Compute the mass fractions.
                                 */
                                
                                std::vector<Real> Y_y_B;
                                std::vector<Real> Y_y_T;
                                Y_y_B.reserve(d_num_species);
                                Y_y_T.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_y_B.push_back(Q[si][idx_cell_rho_Y_y_B]/rho_y_B);
                                    Y_y_T.push_back(Q[si][idx_cell_rho_Y_y_T]/rho_y_T);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions.
                                 */
                                
                                std::vector<const Real*> Y_y_B_ptr;
                                std::vector<const Real*> Y_y_T_ptr;
                                Y_y_B_ptr.reserve(d_num_species);
                                Y_y_T_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_y_B_ptr.push_back(&Y_y_B[si]);
                                    Y_y_T_ptr.push_back(&Y_y_T[si]);
                                }
                                
                                const Real u_y_B = Q[d_num_species][idx_cell_mom_y_B]/rho_y_B;
                                const Real u_y_T = Q[d_num_species][idx_cell_mom_y_T]/rho_y_T;
                                
                                const Real v_y_B = Q[d_num_species + 1][idx_cell_mom_y_B]/rho_y_B;
                                const Real v_y_T = Q[d_num_species + 1][idx_cell_mom_y_T]/rho_y_T;
                                
                                const Real w_y_B = Q[d_num_species + 2][idx_cell_mom_y_B]/rho_y_B;
                                const Real w_y_T = Q[d_num_species + 2][idx_cell_mom_y_T]/rho_y_T;
                                
                                const Real epsilon_y_B = Q[d_num_species + 3][idx_cell_E_y_B]/rho_y_B - half*(u_y_B*u_y_B + v_y_B*v_y_B + w_y_B*w_y_B);
                                const Real epsilon_y_T = Q[d_num_species + 3][idx_cell_E_y_T]/rho_y_T - half*(u_y_T*u_y_T + v_y_T*v_y_T + w_y_T*w_y_T);
                                
                                Real p_y_B = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_y_B,
                                        &epsilon_y_B,
                                        Y_y_B_ptr);
                            
                                Real p_y_T = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_y_T,
                                        &epsilon_y_T,
                                        Y_y_T_ptr);
                                
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
                                
                                const int idx_cell_rho_Y_z_F = (interior_box_hi_idx[0] + num_subghosts_conservative_var[0][0]) +
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
                                
                                /*
                                 * Compute the mixture density.
                                 */
                                
                                Real rho_z_F = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_z_F += Q[si][idx_cell_rho_Y_z_F];
                                }
                                
                                /*
                                 * Compute the mass fractions.
                                 */
                                
                                std::vector<Real> Y_z_F;
                                Y_z_F.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_z_F.push_back(Q[si][idx_cell_rho_Y_z_F]/rho_z_F);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions.
                                 */
                                
                                std::vector<const Real*> Y_z_F_ptr;
                                Y_z_F_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_z_F_ptr.push_back(&Y_z_F[si]);
                                }
                                
                                const Real u_z_F = Q[d_num_species][idx_cell_mom_z_F]/rho_z_F;
                                const Real v_z_F = Q[d_num_species + 1][idx_cell_mom_z_F]/rho_z_F;
                                const Real w_z_F = Q[d_num_species + 2][idx_cell_mom_z_F]/rho_z_F;
                                const Real epsilon_z_F = Q[d_num_species + 3][idx_cell_E_z_F]/rho_z_F - half*(u_z_F*u_z_F + v_z_F*v_z_F + w_z_F*w_z_F);
                                
                                Real p_z_F = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_z_F,
                                        &epsilon_z_F,
                                        Y_z_F_ptr);
                                
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
                                
                                const int idx_cell_rho_Y_z_B = (interior_box_hi_idx[0] + num_subghosts_conservative_var[0][0]) +
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
                                
                                /*
                                 * Compute the mixture density.
                                 */
                                
                                Real rho_z_B = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_z_B += Q[si][idx_cell_rho_Y_z_B];
                                }
                                
                                /*
                                 * Compute the mass fractions.
                                 */
                                
                                std::vector<Real> Y_z_B;
                                Y_z_B.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_z_B.push_back(Q[si][idx_cell_rho_Y_z_B]/rho_z_B);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions.
                                 */
                                
                                std::vector<const Real*> Y_z_B_ptr;
                                Y_z_B_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_z_B_ptr.push_back(&Y_z_B[si]);
                                }
                                
                                const Real u_z_B = Q[d_num_species][idx_cell_mom_z_B]/rho_z_B;
                                const Real v_z_B = Q[d_num_species + 1][idx_cell_mom_z_B]/rho_z_B;
                                const Real w_z_B = Q[d_num_species + 2][idx_cell_mom_z_B]/rho_z_B;
                                const Real epsilon_z_B = Q[d_num_species + 3][idx_cell_E_z_B]/rho_z_B - half*(u_z_B*u_z_B + v_z_B*v_z_B + w_z_B*w_z_B);
                                
                                Real p_z_B = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_z_B,
                                        &epsilon_z_B,
                                        Y_z_B_ptr);
                                
                                // One-sided derivatives.
                                du_dz = (u_x_L - u_z_B)/Real(dx[2]);
                                // dv_dz = (v_x_L - v_z_B)/Real(dx[2]);
                                dw_dz = (w_x_L - w_z_B)/Real(dx[2]);
                                dp_dz = (p_x_L - p_z_B)/Real(dx[2]);
                            }
                            else
                            {
                                const int idx_cell_rho_Y_z_B = (interior_box_hi_idx[0] + num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (k - 1 + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_rho_Y_z_F = (interior_box_hi_idx[0] + num_subghosts_conservative_var[0][0]) +
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
                                
                                /*
                                 * Compute the mixture density.
                                 */
                                
                                Real rho_z_B = Real(0);
                                Real rho_z_F = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_z_B += Q[si][idx_cell_rho_Y_z_B];
                                    rho_z_F += Q[si][idx_cell_rho_Y_z_F];
                                }
                                
                                /*
                                 * Compute the mass fractions.
                                 */
                                
                                std::vector<Real> Y_z_B;
                                std::vector<Real> Y_z_F;
                                Y_z_B.reserve(d_num_species);
                                Y_z_F.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_z_B.push_back(Q[si][idx_cell_rho_Y_z_B]/rho_z_B);
                                    Y_z_F.push_back(Q[si][idx_cell_rho_Y_z_F]/rho_z_F);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions.
                                 */
                                
                                std::vector<const Real*> Y_z_B_ptr;
                                std::vector<const Real*> Y_z_F_ptr;
                                Y_z_B_ptr.reserve(d_num_species);
                                Y_z_F_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_z_B_ptr.push_back(&Y_z_B[si]);
                                    Y_z_F_ptr.push_back(&Y_z_F[si]);
                                }
                                
                                const Real u_z_B = Q[d_num_species][idx_cell_mom_z_B]/rho_z_B;
                                const Real u_z_F = Q[d_num_species][idx_cell_mom_z_F]/rho_z_F;
                                
                                const Real v_z_B = Q[d_num_species + 1][idx_cell_mom_z_B]/rho_z_B;
                                const Real v_z_F = Q[d_num_species + 1][idx_cell_mom_z_F]/rho_z_F;
                                
                                const Real w_z_B = Q[d_num_species + 2][idx_cell_mom_z_B]/rho_z_B;
                                const Real w_z_F = Q[d_num_species + 2][idx_cell_mom_z_F]/rho_z_F;
                                
                                const Real epsilon_z_B = Q[d_num_species + 3][idx_cell_E_z_B]/rho_z_B - half*(u_z_B*u_z_B + v_z_B*v_z_B + w_z_B*w_z_B);
                                const Real epsilon_z_F = Q[d_num_species + 3][idx_cell_E_z_F]/rho_z_F - half*(u_z_F*u_z_F + v_z_F*v_z_F + w_z_F*w_z_F);
                                
                                Real p_z_B = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_z_B,
                                        &epsilon_z_B,
                                        Y_z_B_ptr);
                                
                                Real p_z_F = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_z_F,
                                        &epsilon_z_F,
                                        Y_z_F_ptr);
                                
                                // Central derivatives.
                                du_dz = (u_z_F - u_z_B)/(Real(2)*Real(dx[2]));
                                // dv_dz = (v_z_F - v_z_B)/(Real(2)*Real(dx[2]));
                                dw_dz = (w_z_F - w_z_B)/(Real(2)*Real(dx[2]));
                                dp_dz = (p_z_F - p_z_B)/(Real(2)*Real(dx[2]));
                            }
                            
                            // Compute sound speed.
                            
                            const Real Gamma_x_L = d_equation_of_state_mixing_rules->getGruneisenParameter(
                                &rho_x_L,
                                &p_x_L,
                                Y_x_L_ptr);
                            
                            const std::vector<Real> Psi_x_L = d_equation_of_state_mixing_rules->
                                getPressureDerivativeWithPartialDensities(
                                        &rho_x_L,
                                        &p_x_L,
                                        Y_x_L_ptr);
                            
                            Real c_x_L = Gamma_x_L*p_x_L/rho_x_L;
                            for (int si = 0; si < d_num_species; si++)
                            {
                                c_x_L += Y_x_L[si]*Psi_x_L[si];
                            }
                            c_x_L = std::sqrt(c_x_L);
                            
                            const Real lambda_1 = u_x_L - c_x_L;
                            
                            // Compute vector Lambda^(-1) * L.
                            
                            Real Lambda_inv_L[d_num_species + 4];
                            
                            const Real& p_t         = d_bdry_face_nonreflecting_outflow_p_t[face_loc];
                            const Real& sigma       = d_bdry_face_nonreflecting_outflow_sigma[face_loc];
                            const Real& beta        = d_bdry_face_nonreflecting_outflow_beta[face_loc];
                            const Real& length_char = d_bdry_face_nonreflecting_outflow_length_char[face_loc];
                            
                            const Real T_1 = v_x_L*(dp_dy - rho_x_L*c_x_L*du_dy) + rho_x_L*c_x_L*c_x_L*dv_dy +
                                w_x_L*(dp_dz - rho_x_L*c_x_L*du_dz) + rho_x_L*c_x_L*c_x_L*dw_dz;
                            
                            const Real M_sq = (u_x_L*u_x_L + v_x_L*v_x_L + w_x_L*w_x_L)/(c_x_L*c_x_L);
                            const Real K = sigma*c_x_L*(Real(1) - M_sq)/length_char;
                            
                            Lambda_inv_L[0] = (Real(1)/lambda_1)*(K*(p_x_L - p_t) - (Real(1) - beta)*T_1);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Lambda_inv_L[si + 1] = c_x_L*c_x_L*drho_Y_dx[si] - Y_x_L[si]*dp_dx;
                            }
                            Lambda_inv_L[d_num_species + 1] = dv_dx;
                            Lambda_inv_L[d_num_species + 2] = dw_dx;
                            Lambda_inv_L[d_num_species + 3] = dp_dx + rho_x_L*c_x_L*du_dx;
                            
                            // Compute dV_dx.
                            
                            const Real c_sq_inv  = Real(1)/(c_x_L*c_x_L);
                            const Real rho_c_inv = Real(1)/(rho_x_L*c_x_L);
                            
                            Real dV_dx[d_num_species + 4];
                            
                            for (int si = 0; si < d_num_species; si++)
                            {
                                dV_dx[si] = half*c_sq_inv*Y_x_L[si]*(Lambda_inv_L[0] + Lambda_inv_L[d_num_species + 3]) + 
                                    c_sq_inv*Lambda_inv_L[si + 1];
                            }
                            dV_dx[d_num_species] = half*rho_c_inv*(-Lambda_inv_L[0] + Lambda_inv_L[d_num_species + 3]);
                            dV_dx[d_num_species + 1] = Lambda_inv_L[d_num_species + 1];
                            dV_dx[d_num_species + 2] = Lambda_inv_L[d_num_species + 2];
                            dV_dx[d_num_species + 3] = half*(Lambda_inv_L[0] + Lambda_inv_L[d_num_species + 3]);
                            
                            Real V_ghost[(d_num_species + 4)*num_ghosts_to_fill];
                            
                            for (int i = 0; i < num_ghosts_to_fill; i++)
                            {
                                const int idx_cell_rho_Y = (i + fill_box_lo_idx[0] + num_subghosts_conservative_var[0][0]) +
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
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        V_ghost[i*(d_num_species + 4) + si] = rho_Y_x_LL[si] + Real(2)*Real(dx[0])*dV_dx[si];
                                    }
                                    
                                    V_ghost[i*(d_num_species + 4) + d_num_species]     = u_x_LL + Real(2)*Real(dx[0])*dV_dx[d_num_species];
                                    V_ghost[i*(d_num_species + 4) + d_num_species + 1] = v_x_LL + Real(2)*Real(dx[0])*dV_dx[d_num_species + 1];
                                    V_ghost[i*(d_num_species + 4) + d_num_species + 2] = w_x_LL + Real(2)*Real(dx[0])*dV_dx[d_num_species + 2];
                                    V_ghost[i*(d_num_species + 4) + d_num_species + 3] = p_x_LL + Real(2)*Real(dx[0])*dV_dx[d_num_species + 3];
                                }
                                else if (i == 1)
                                {
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        V_ghost[i*(d_num_species + 4) + si] = -Real(2)*rho_Y_x_LL[si] - Real(3)*rho_Y_x_L[si] +
                                            Real(6)*V_ghost[(i - 1)*(d_num_species + 4) + si] - Real(6)*Real(dx[0])*dV_dx[si];
                                    }
                                    
                                    V_ghost[i*(d_num_species + 4) + d_num_species] = -Real(2)*u_x_LL - Real(3)*u_x_L +
                                        Real(6)*V_ghost[(i - 1)*(d_num_species + 4) + d_num_species] -
                                        Real(6)*Real(dx[0])*dV_dx[d_num_species];
                                    
                                    V_ghost[i*(d_num_species + 4) + d_num_species + 1] = -Real(2)*v_x_LL - Real(3)*v_x_L +
                                        Real(6)*V_ghost[(i - 1)*(d_num_species + 4) + d_num_species + 1] -
                                        Real(6)*Real(dx[0])*dV_dx[d_num_species + 1];
                                    
                                    V_ghost[i*(d_num_species + 4) + d_num_species + 2] = -Real(2)*w_x_LL - Real(3)*w_x_L +
                                        Real(6)*V_ghost[(i - 1)*(d_num_species + 4) + d_num_species + 2] -
                                        Real(6)*Real(dx[0])*dV_dx[d_num_species + 2];
                                    
                                    V_ghost[i*(d_num_species + 4) + d_num_species + 3] = -Real(2)*p_x_LL - Real(3)*p_x_L +
                                        Real(6)*V_ghost[(i - 1)*(d_num_species + 4) + d_num_species + 3] -
                                        Real(6)*Real(dx[0])*dV_dx[d_num_species + 3];
                                }
                                else if (i == 2)
                                {
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        V_ghost[i*(d_num_species + 4) + si] = Real(3)*rho_Y_x_LL[si] + Real(10)*rho_Y_x_L[si] -
                                            Real(18)*V_ghost[(i - 2)*(d_num_species + 4) + si] +
                                            Real(6)*V_ghost[(i - 1)*(d_num_species + 4) + si] +
                                            Real(12)*Real(dx[0])*dV_dx[si];
                                    }
                                    
                                    V_ghost[i*(d_num_species + 4) + d_num_species] = Real(3)*u_x_LL + Real(10)*u_x_L -
                                        Real(18)*V_ghost[(i - 2)*(d_num_species + 4) + d_num_species] +
                                        Real(6)*V_ghost[(i - 1)*(d_num_species + 4) + d_num_species] +
                                        Real(12)*Real(dx[0])*dV_dx[d_num_species];
                                    
                                    V_ghost[i*(d_num_species + 4) + d_num_species + 1] = Real(3)*v_x_LL + Real(10)*v_x_L -
                                        Real(18)*V_ghost[(i - 2)*(d_num_species + 4) + d_num_species + 1] +
                                        Real(6)*V_ghost[(i - 1)*(d_num_species + 4) + d_num_species + 1] +
                                        Real(12)*Real(dx[0])*dV_dx[d_num_species + 1];
                                    
                                    V_ghost[i*(d_num_species + 4) + d_num_species + 2] = Real(3)*w_x_LL + Real(10)*w_x_L -
                                        Real(18)*V_ghost[(i - 2)*(d_num_species + 4) + d_num_species + 2] +
                                        Real(6)*V_ghost[(i - 1)*(d_num_species + 4) + d_num_species + 2] +
                                        Real(12)*Real(dx[0])*dV_dx[d_num_species + 2];
                                    
                                    V_ghost[i*(d_num_species + 4) + d_num_species + 3] = Real(3)*p_x_LL + Real(10)*p_x_L -
                                        Real(18)*V_ghost[(i - 2)*(d_num_species + 4) + d_num_species + 3] +
                                        Real(6)*V_ghost[(i - 1)*(d_num_species + 4) + d_num_species + 3] +
                                        Real(12)*Real(dx[0])*dV_dx[d_num_species + 3];
                                }
                                else if (i == 3)
                                {
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        V_ghost[i*(d_num_species + 4) + si] = -Real(4)*rho_Y_x_LL[si] -
                                            Real(65)/Real(3)*rho_Y_x_L[si] +
                                            Real(40)*V_ghost[(i - 3)*(d_num_species + 4) + si] -
                                            Real(20)*V_ghost[(i - 2)*(d_num_species + 4) + si] +
                                            Real(20)/Real(3)*V_ghost[(i - 1)*(d_num_species + 4) + si] -
                                            Real(20)*Real(dx[0])*dV_dx[si];
                                    }
                                    
                                    V_ghost[i*(d_num_species + 4) + d_num_species] = -Real(4)*u_x_LL -
                                        Real(65)/Real(3)*u_x_L +
                                        Real(40)*V_ghost[(i - 3)*(d_num_species + 4) + d_num_species] -
                                        Real(20)*V_ghost[(i - 2)*(d_num_species + 4) + d_num_species] +
                                        Real(20)/Real(3)*V_ghost[(i - 1)*(d_num_species + 4) + d_num_species] -
                                        Real(20)*Real(dx[0])*dV_dx[d_num_species];
                                    
                                    V_ghost[i*(d_num_species + 4) + d_num_species + 1] = -Real(4)*v_x_LL -
                                        Real(65)/Real(3)*v_x_L +
                                        Real(40)*V_ghost[(i - 3)*(d_num_species + 4) + d_num_species + 1] -
                                        Real(20)*V_ghost[(i - 2)*(d_num_species + 4) + d_num_species + 1] +
                                        Real(20)/Real(3)*V_ghost[(i - 1)*(d_num_species + 4) + d_num_species + 1] -
                                        Real(20)*Real(dx[0])*dV_dx[d_num_species + 1];
                                    
                                    V_ghost[i*(d_num_species + 4) + d_num_species + 2] = -Real(4)*v_x_LL -
                                        Real(65)/Real(3)*v_x_L +
                                        Real(40)*V_ghost[(i - 3)*(d_num_species + 4) + d_num_species + 2] -
                                        Real(20)*V_ghost[(i - 2)*(d_num_species + 4) + d_num_species + 2] +
                                        Real(20)/Real(3)*V_ghost[(i - 1)*(d_num_species + 4) + d_num_species + 2] -
                                        Real(20)*Real(dx[0])*dV_dx[d_num_species + 2];
                                    
                                    V_ghost[i*(d_num_species + 4) + d_num_species + 3] = -Real(4)*p_x_LL -
                                        Real(65)/Real(3)*p_x_L +
                                        Real(40)*V_ghost[(i - 3)*(d_num_species + 4) + d_num_species + 3] -
                                        Real(20)*V_ghost[(i - 2)*(d_num_species + 4) + d_num_species + 3] +
                                        Real(20)/Real(3)*V_ghost[(i - 1)*(d_num_species + 4) + d_num_species + 3] -
                                        Real(20)*Real(dx[0])*dV_dx[d_num_species + 3];
                                }
                                else if (i == 4)
                                {
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        V_ghost[i*(d_num_species + 4) + si] = Real(5)*rho_Y_x_LL[si] +
                                            Real(77)/Real(2)*rho_Y_x_L[si] -
                                            Real(75)*V_ghost[(i - 4)*(d_num_species + 4) + si] +
                                            Real(50)*V_ghost[(i - 3)*(d_num_species + 4) + si] -
                                            Real(25)*V_ghost[(i - 2)*(d_num_species + 4) + si] +
                                            Real(15)/Real(2)*V_ghost[(i - 1)*(d_num_species + 4) + si] +
                                            Real(30)*Real(dx[0])*dV_dx[si];
                                    }
                                    
                                    V_ghost[i*(d_num_species + 4) + d_num_species] = Real(5)*u_x_LL +
                                        Real(77)/Real(2)*u_x_L -
                                        Real(75)*V_ghost[(i - 4)*(d_num_species + 4) + d_num_species] +
                                        Real(50)*V_ghost[(i - 3)*(d_num_species + 4) + d_num_species] -
                                        Real(25)*V_ghost[(i - 2)*(d_num_species + 4) + d_num_species] +
                                        Real(15)/Real(2)*V_ghost[(i - 1)*(d_num_species + 4) + d_num_species] +
                                        Real(30)*Real(dx[0])*dV_dx[d_num_species];
                                    
                                    V_ghost[i*(d_num_species + 4) + d_num_species + 1] = Real(5)*v_x_LL +
                                        Real(77)/Real(2)*v_x_L -
                                        Real(75)*V_ghost[(i - 4)*(d_num_species + 4) + d_num_species + 1] +
                                        Real(50)*V_ghost[(i - 3)*(d_num_species + 4) + d_num_species + 1] -
                                        Real(25)*V_ghost[(i - 2)*(d_num_species + 4) + d_num_species + 1] +
                                        Real(15)/Real(2)*V_ghost[(i - 1)*(d_num_species + 4) + d_num_species + 1] +
                                        Real(30)*Real(dx[0])*dV_dx[d_num_species + 1];
                                    
                                    V_ghost[i*(d_num_species + 4) + d_num_species + 2] = Real(5)*w_x_LL +
                                        Real(77)/Real(2)*w_x_L -
                                        Real(75)*V_ghost[(i - 4)*(d_num_species + 4) + d_num_species + 2] +
                                        Real(50)*V_ghost[(i - 3)*(d_num_species + 4) + d_num_species + 2] -
                                        Real(25)*V_ghost[(i - 2)*(d_num_species + 4) + d_num_species + 2] +
                                        Real(15)/Real(2)*V_ghost[(i - 1)*(d_num_species + 4) + d_num_species + 2] +
                                        Real(30)*Real(dx[0])*dV_dx[d_num_species + 2];
                                    
                                    V_ghost[i*(d_num_species + 4) + d_num_species + 3] = Real(5)*p_x_LL +
                                        Real(77)/Real(2)*p_x_L -
                                        Real(75)*V_ghost[(i - 4)*(d_num_species + 4) + d_num_species + 3] +
                                        Real(50)*V_ghost[(i - 3)*(d_num_species + 4) + d_num_species + 3] -
                                        Real(25)*V_ghost[(i - 2)*(d_num_species + 4) + d_num_species + 3] +
                                        Real(15)/Real(2)*V_ghost[(i - 1)*(d_num_species + 4) + d_num_species + 3] +
                                        Real(30)*Real(dx[0])*dV_dx[d_num_species + 3];
                                }
                                else if (i == 5)
                                {
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        V_ghost[i*(d_num_species + 4) + si] = -Real(6)*rho_Y_x_LL[si] -
                                            Real(609)/Real(10)*rho_Y_x_L[si] +
                                            Real(126)*V_ghost[(i - 5)*(d_num_species + 4) + si] -
                                            Real(105)*V_ghost[(i - 4)*(d_num_species + 4) + si] +
                                            Real(70)*V_ghost[(i - 3)*(d_num_species + 4) + si] -
                                            Real(63)/Real(2)*V_ghost[(i - 2)*(d_num_species + 4) + si] +
                                            Real(42)/Real(5)*V_ghost[(i - 1)*(d_num_species + 4) + si] -
                                            Real(42)*Real(dx[0])*dV_dx[si];
                                    }
                                    
                                    V_ghost[i*(d_num_species + 4) + d_num_species] = -Real(6)*u_x_LL -
                                        Real(609)/Real(10)*u_x_L +
                                        Real(126)*V_ghost[(i - 5)*(d_num_species + 4) + d_num_species] -
                                        Real(105)*V_ghost[(i - 4)*(d_num_species + 4) + d_num_species] +
                                        Real(70)*V_ghost[(i - 3)*(d_num_species + 4) + d_num_species] -
                                        Real(63)/Real(2)*V_ghost[(i - 2)*(d_num_species + 4) + d_num_species] +
                                        Real(42)/Real(5)*V_ghost[(i - 1)*(d_num_species + 4) + d_num_species] -
                                        Real(42)*Real(dx[0])*dV_dx[d_num_species];
                                    
                                    V_ghost[i*(d_num_species + 4) + d_num_species + 1] = -Real(6)*v_x_LL -
                                        Real(609)/Real(10)*v_x_L +
                                        Real(126)*V_ghost[(i - 5)*(d_num_species + 4) + d_num_species + 1] -
                                        Real(105)*V_ghost[(i - 4)*(d_num_species + 4) + d_num_species + 1] +
                                        Real(70)*V_ghost[(i - 3)*(d_num_species + 4) + d_num_species + 1] -
                                        Real(63)/Real(2)*V_ghost[(i - 2)*(d_num_species + 4) + d_num_species + 1] +
                                        Real(42)/Real(5)*V_ghost[(i - 1)*(d_num_species + 4) + d_num_species + 1] -
                                        Real(42)*Real(dx[0])*dV_dx[d_num_species + 1];
                                    
                                    V_ghost[i*(d_num_species + 4) + d_num_species + 2] = -Real(6)*w_x_LL -
                                        Real(609)/Real(10)*w_x_L +
                                        Real(126)*V_ghost[(i - 5)*(d_num_species + 4) + d_num_species + 2] -
                                        Real(105)*V_ghost[(i - 4)*(d_num_species + 4) + d_num_species + 2] +
                                        Real(70)*V_ghost[(i - 3)*(d_num_species + 4) + d_num_species + 2] -
                                        Real(63)/Real(2)*V_ghost[(i - 2)*(d_num_species + 4) + d_num_species + 2] +
                                        Real(42)/Real(5)*V_ghost[(i - 1)*(d_num_species + 4) + d_num_species + 2] -
                                        Real(42)*Real(dx[0])*dV_dx[d_num_species + 2];
                                    
                                    V_ghost[i*(d_num_species + 4) + d_num_species + 3] = -Real(6)*p_x_LL -
                                        Real(609)/Real(10)*p_x_L +
                                        Real(126)*V_ghost[(i - 5)*(d_num_species + 4) + d_num_species + 3] -
                                        Real(105)*V_ghost[(i - 4)*(d_num_species + 4) + d_num_species + 3] +
                                        Real(70)*V_ghost[(i - 3)*(d_num_species + 4) + d_num_species + 3] -
                                        Real(63)/Real(2)*V_ghost[(i - 2)*(d_num_species + 4) + d_num_species + 3] +
                                        Real(42)/Real(5)*V_ghost[(i - 1)*(d_num_species + 4) + d_num_species + 3] -
                                        Real(42)*Real(dx[0])*dV_dx[d_num_species + 3];
                                }
                                
                                /*
                                 * Compute the mixture density.
                                 */
                                
                                Real rho_ghost = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_ghost += V_ghost[i*(d_num_species + 4) + si];
                                }
                                
                                /*
                                 * Compute the mass fractions.
                                 */
                                
                                std::vector<Real> Y_ghost;
                                Y_ghost.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_ghost.push_back(V_ghost[i*(d_num_species + 4) + si]/rho_ghost);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions.
                                 */
                                
                                std::vector<const Real*> Y_ghost_ptr;
                                Y_ghost_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_ghost_ptr.push_back(&Y_ghost[si]);
                                }
                                
                                for(int si=0; si < d_num_species; si++)
                                {
                                    Q[si][idx_cell_rho_Y] = V_ghost[i*(d_num_species + 4) + si];
                                }
                                
                                Q[d_num_species][idx_cell_mom]     = rho_ghost*V_ghost[i*(d_num_species + 4) + d_num_species];
                                Q[d_num_species + 1][idx_cell_mom] = rho_ghost*V_ghost[i*(d_num_species + 4) + d_num_species + 1];
                                Q[d_num_species + 2][idx_cell_mom] = rho_ghost*V_ghost[i*(d_num_species + 4) + d_num_species + 2];
                                
                                const Real epsilon = d_equation_of_state_mixing_rules->
                                    getInternalEnergy(
                                        &rho_ghost,
                                        &V_ghost[i*(d_num_species + 4) + d_num_species + 3],
                                        Y_ghost_ptr);
                                
                                const Real E = rho_ghost*epsilon +
                                    half*(Q[d_num_species][idx_cell_mom]*Q[d_num_species][idx_cell_mom] +
                                        Q[d_num_species + 1][idx_cell_mom]*Q[d_num_species + 1][idx_cell_mom] +
                                        Q[d_num_species + 2][idx_cell_mom]*Q[d_num_species + 2][idx_cell_mom])/
                                        rho_ghost;
                                
                                Q[d_num_species + 3][idx_cell_E] = E;
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
                            << ": FlowModelBoundaryUtilitiesFourEqnConservative::fill3dFaceBoundaryData()\n"
                            << "Non-reflecting outflow BC doesn't support more than six ghost cells yet!");
                    }
                    
                    for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                    {
                        for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                        {
                            // Get the grid spacing.
                            const double* const dx = patch_geom->getDx();
                            
                            const int idx_cell_rho_Y_y_T = (i + num_subghosts_conservative_var[0][0]) +
                                (interior_box_lo_idx[1] + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                (k + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                    subghostcell_dims_conservative_var[0][1];
                            
                            const int idx_cell_rho_Y_y_TT = (i + num_subghosts_conservative_var[0][0]) +
                                (interior_box_lo_idx[1] + 1 + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                (k + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                    subghostcell_dims_conservative_var[0][1];
                            
                            const int idx_cell_rho_Y_y_TTT = (i + num_subghosts_conservative_var[0][0]) +
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
                            
                            std::vector<Real> rho_Y_y_T;
                            std::vector<Real> rho_Y_y_TT;
                            std::vector<Real> rho_Y_y_TTT;
                            rho_Y_y_T.reserve(d_num_species);
                            rho_Y_y_TT.reserve(d_num_species);
                            rho_Y_y_TTT.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                rho_Y_y_T.push_back(Q[si][idx_cell_rho_Y_y_T]);
                                rho_Y_y_TT.push_back(Q[si][idx_cell_rho_Y_y_TT]);
                                rho_Y_y_TTT.push_back(Q[si][idx_cell_rho_Y_y_TTT]);
                            }
                            
                            /*
                             * Compute the mixture density.
                             */
                            
                            Real rho_y_T   = Real(0);
                            Real rho_y_TT  = Real(0);
                            Real rho_y_TTT = Real(0);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                rho_y_T   += Q[si][idx_cell_rho_Y_y_T];
                                rho_y_TT  += Q[si][idx_cell_rho_Y_y_TT];
                                rho_y_TTT += Q[si][idx_cell_rho_Y_y_TTT];
                            }
                            
                            /*
                             * Compute the mass fractions.
                             */
                            
                            std::vector<Real> Y_y_T;
                            std::vector<Real> Y_y_TT;
                            std::vector<Real> Y_y_TTT;
                            Y_y_T.reserve(d_num_species);
                            Y_y_TT.reserve(d_num_species);
                            Y_y_TTT.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_y_T.push_back(Q[si][idx_cell_rho_Y_y_T]/rho_y_T);
                                Y_y_TT.push_back(Q[si][idx_cell_rho_Y_y_TT]/rho_y_TT);
                                Y_y_TTT.push_back(Q[si][idx_cell_rho_Y_y_TTT]/rho_y_TTT);
                            }
                            
                            /*
                             * Get the pointers to the mass fractions.
                             */
                            
                            std::vector<const Real*> Y_y_T_ptr;
                            std::vector<const Real*> Y_y_TT_ptr;
                            std::vector<const Real*> Y_y_TTT_ptr;
                            Y_y_T_ptr.reserve(d_num_species);
                            Y_y_TT_ptr.reserve(d_num_species);
                            Y_y_TTT_ptr.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_y_T_ptr.push_back(&Y_y_T[si]);
                                Y_y_TT_ptr.push_back(&Y_y_TT[si]);
                                Y_y_TTT_ptr.push_back(&Y_y_TTT[si]);
                            }
                            
                            const Real u_y_T   = Q[d_num_species][idx_cell_mom_y_T]/rho_y_T;
                            const Real u_y_TT  = Q[d_num_species][idx_cell_mom_y_TT]/rho_y_TT;
                            const Real u_y_TTT = Q[d_num_species][idx_cell_mom_y_TTT]/rho_y_TTT;
                            
                            const Real v_y_T   = Q[d_num_species + 1][idx_cell_mom_y_T]/rho_y_T;
                            const Real v_y_TT  = Q[d_num_species + 1][idx_cell_mom_y_TT]/rho_y_TT;
                            const Real v_y_TTT = Q[d_num_species + 1][idx_cell_mom_y_TTT]/rho_y_TTT;
                            
                            const Real w_y_T   = Q[d_num_species + 2][idx_cell_mom_y_T]/rho_y_T;
                            const Real w_y_TT  = Q[d_num_species + 2][idx_cell_mom_y_TT]/rho_y_TT;
                            const Real w_y_TTT = Q[d_num_species + 2][idx_cell_mom_y_TTT]/rho_y_TTT;
                            
                            const Real half = Real(1)/Real(2);
                            const Real epsilon_y_T   = Q[d_num_species + 3][idx_cell_E_y_T]/rho_y_T - half*(u_y_T*u_y_T + v_y_T*v_y_T + w_y_T*w_y_T);
                            const Real epsilon_y_TT  = Q[d_num_species + 3][idx_cell_E_y_TT]/rho_y_TT - half*(u_y_TT*u_y_TT + v_y_TT*v_y_TT + w_y_TT*w_y_TT);
                            const Real epsilon_y_TTT = Q[d_num_species + 3][idx_cell_E_y_TTT]/rho_y_TTT - half*(u_y_TTT*u_y_TTT + v_y_TTT*v_y_TTT + w_y_TTT*w_y_TTT);

                            Real p_y_T = d_equation_of_state_mixing_rules->
                                getPressure(
                                    &rho_y_T,
                                    &epsilon_y_T,
                                    Y_y_T_ptr);
                            
                            Real p_y_TT = d_equation_of_state_mixing_rules->
                                getPressure(
                                    &rho_y_TT,
                                    &epsilon_y_TT,
                                    Y_y_TT_ptr);
                            
                            Real p_y_TTT = d_equation_of_state_mixing_rules->
                                getPressure(
                                    &rho_y_TTT,
                                    &epsilon_y_TTT,
                                    Y_y_TTT_ptr);
                            
                            /*
                             * Compute derivatives in y-direction.
                             */
                            
                            std::vector<Real> drho_Y_dy;
                            drho_Y_dy.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                drho_Y_dy.push_back(-(Q[si][idx_cell_rho_Y_y_TTT] - Real(4)*Q[si][idx_cell_rho_Y_y_TT] +
                                    Real(3)*Q[si][idx_cell_rho_Y_y_T])/(Real(2)*Real(dx[1])));
                            }
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
                                
                                const int idx_cell_rho_Y_x_R = (i + 1 + num_subghosts_conservative_var[0][0]) +
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
                                
                                /*
                                 * Compute the mixture density.
                                 */
                                
                                Real rho_x_R = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_x_R += Q[si][idx_cell_rho_Y_x_R];
                                }
                                
                                /*
                                 * Compute the mass fractions.
                                 */
                                
                                std::vector<Real> Y_x_R;
                                Y_x_R.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_x_R.push_back(Q[si][idx_cell_rho_Y_x_R]/rho_x_R);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions.
                                 */
                            
                                std::vector<const Real*> Y_x_R_ptr;
                                Y_x_R_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_x_R_ptr.push_back(&Y_x_R[si]);
                                }
                                
                                const Real u_x_R = Q[d_num_species][idx_cell_mom_x_R]/rho_x_R;
                                const Real v_x_R = Q[d_num_species + 1][idx_cell_mom_x_R]/rho_x_R;
                                const Real w_x_R = Q[d_num_species + 2][idx_cell_mom_x_R]/rho_x_R;
                                const Real epsilon_x_R = Q[d_num_species + 3][idx_cell_E_x_R]/rho_x_R - half*(u_x_R*u_x_R + v_x_R*v_x_R + w_x_R*w_x_R);
                                
                                Real p_x_R = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_x_R,
                                        &epsilon_x_R,
                                        Y_x_R_ptr);
                                
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
                                
                                const int idx_cell_rho_Y_x_L = (i - 1 + num_subghosts_conservative_var[0][0]) +
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
                                
                                /*
                                 * Compute the mixture density.
                                 */
                                
                                Real rho_x_L = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_x_L += Q[si][idx_cell_rho_Y_x_L];
                                }
                                
                                /*
                                 * Compute the mass fractions.
                                 */
                                
                                std::vector<Real> Y_x_L;
                                Y_x_L.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_x_L.push_back(Q[si][idx_cell_rho_Y_x_L]/rho_x_L);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions.
                                 */
                                
                                std::vector<const Real*> Y_x_L_ptr;
                                Y_x_L_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_x_L_ptr.push_back(&Y_x_L[si]);
                                }
                                
                                const Real u_x_L = Q[d_num_species][idx_cell_mom_x_L]/rho_x_L;
                                const Real v_x_L = Q[d_num_species + 1][idx_cell_mom_x_L]/rho_x_L;
                                const Real w_x_L = Q[d_num_species + 2][idx_cell_mom_x_L]/rho_x_L;
                                const Real epsilon_x_L = Q[d_num_species + 3][idx_cell_E_x_L]/rho_x_L - half*(u_x_L*u_x_L + v_x_L*v_x_L + w_x_L*w_x_L);
                                
                                Real p_x_L = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_x_L,
                                        &epsilon_x_L,
                                        Y_x_L_ptr);
                                
                                // One-sided derivatives.
                                du_dx = (u_y_T - u_x_L)/Real(dx[0]);
                                dv_dx = (v_y_T - v_x_L)/Real(dx[0]);
                                // dw_dx = (w_y_T - w_x_L)/Real(dx[0]);
                                dp_dx = (p_y_T - p_x_L)/Real(dx[0]);
                            }
                            else
                            {
                                const int idx_cell_rho_Y_x_L = (i - 1 + num_subghosts_conservative_var[0][0]) +
                                    (interior_box_lo_idx[1] + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (k + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_rho_Y_x_R = (i + 1 + num_subghosts_conservative_var[0][0]) +
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
                                
                                /*
                                 * Compute the mixture density.
                                 */
                                
                                Real rho_x_L = Real(0);
                                Real rho_x_R = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_x_L += Q[si][idx_cell_rho_Y_x_L];
                                    rho_x_R += Q[si][idx_cell_rho_Y_x_R];
                                }
                                
                                /*
                                 * Compute the mass fractions.
                                 */
                                
                                std::vector<Real> Y_x_L;
                                std::vector<Real> Y_x_R;
                                Y_x_L.reserve(d_num_species);
                                Y_x_R.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_x_L.push_back(Q[si][idx_cell_rho_Y_x_L]/rho_x_L);
                                    Y_x_R.push_back(Q[si][idx_cell_rho_Y_x_R]/rho_x_R);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions.
                                 */
                                
                                std::vector<const Real*> Y_x_L_ptr;
                                std::vector<const Real*> Y_x_R_ptr;
                                
                                Y_x_L_ptr.reserve(d_num_species);
                                Y_x_R_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_x_L_ptr.push_back(&Y_x_L[si]);
                                    Y_x_R_ptr.push_back(&Y_x_R[si]);
                                }
                                
                                const Real u_x_L = Q[d_num_species][idx_cell_mom_x_L]/rho_x_L;
                                const Real u_x_R = Q[d_num_species][idx_cell_mom_x_R]/rho_x_R;
                                
                                const Real v_x_L = Q[d_num_species + 1][idx_cell_mom_x_L]/rho_x_L;
                                const Real v_x_R = Q[d_num_species + 1][idx_cell_mom_x_R]/rho_x_R;
                                
                                const Real w_x_L = Q[d_num_species + 2][idx_cell_mom_x_L]/rho_x_L;
                                const Real w_x_R = Q[d_num_species + 2][idx_cell_mom_x_R]/rho_x_R;
                                
                                const Real epsilon_x_L = Q[d_num_species + 3][idx_cell_E_x_L]/rho_x_L - half*(u_x_L*u_x_L + v_x_L*v_x_L + w_x_L*w_x_L);
                                const Real epsilon_x_R = Q[d_num_species + 3][idx_cell_E_x_R]/rho_x_R - half*(u_x_R*u_x_R + v_x_R*v_x_R + w_x_R*w_x_R);
                                
                                Real p_x_L = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_x_L,
                                        &epsilon_x_L,
                                        Y_x_L_ptr);
                                
                                Real p_x_R = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_x_R,
                                        &epsilon_x_R,
                                        Y_x_R_ptr);
                                
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
                                
                                const int idx_cell_rho_Y_z_F = (i + num_subghosts_conservative_var[0][0]) +
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
                                
                                /*
                                 * Compute the mixture density.
                                 */
                                
                                Real rho_z_F = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_z_F += Q[si][idx_cell_rho_Y_z_F];
                                }
                                
                                /*
                                 * Compute the mass fractions.
                                 */
                                
                                std::vector<Real> Y_z_F;
                                Y_z_F.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_z_F.push_back(Q[si][idx_cell_rho_Y_z_F]/rho_z_F);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions.
                                 */
                                
                                std::vector<const Real*> Y_z_F_ptr;
                                Y_z_F_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_z_F_ptr.push_back(&Y_z_F[si]);
                                }
                                
                                const Real u_z_F = Q[d_num_species][idx_cell_mom_z_F]/rho_z_F;
                                const Real v_z_F = Q[d_num_species + 1][idx_cell_mom_z_F]/rho_z_F;
                                const Real w_z_F = Q[d_num_species + 2][idx_cell_mom_z_F]/rho_z_F;
                                const Real epsilon_z_F = Q[d_num_species + 3][idx_cell_E_z_F]/rho_z_F - half*(u_z_F*u_z_F + v_z_F*v_z_F + w_z_F*w_z_F);
                                
                                Real p_z_F = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_z_F,
                                        &epsilon_z_F,
                                        Y_z_F_ptr);
                                
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
                                
                                const int idx_cell_rho_Y_z_B = (i + num_subghosts_conservative_var[0][0]) +
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
                                
                                /*
                                 * Compute the mixture density.
                                 */
                                
                                Real rho_z_B = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_z_B += Q[si][idx_cell_rho_Y_z_B];
                                }
                                
                                /*
                                 * Compute the mass fractions.
                                 */
                                
                                std::vector<Real> Y_z_B;
                                Y_z_B.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_z_B.push_back(Q[si][idx_cell_rho_Y_z_B]/rho_z_B);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions.
                                 */
                                
                                std::vector<const Real*> Y_z_B_ptr;
                                Y_z_B_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_z_B_ptr.push_back(&Y_z_B[si]);
                                }
                                
                                const Real u_z_B = Q[d_num_species][idx_cell_mom_z_B]/rho_z_B;
                                const Real v_z_B = Q[d_num_species + 1][idx_cell_mom_z_B]/rho_z_B;
                                const Real w_z_B = Q[d_num_species + 2][idx_cell_mom_z_B]/rho_z_B;
                                const Real epsilon_z_B = Q[d_num_species + 3][idx_cell_E_z_B]/rho_z_B - half*(u_z_B*u_z_B + v_z_B*v_z_B + w_z_B*w_z_B);
                                
                                Real p_z_B = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_z_B,
                                        &epsilon_z_B,
                                        Y_z_B_ptr);
                                
                                // One-sided derivatives.
                                // du_dz = (u_y_T - u_z_B)/Real(dx[2]);
                                dv_dz = (v_y_T - v_z_B)/Real(dx[2]);
                                dw_dz = (w_y_T - w_z_B)/Real(dx[2]);
                                dp_dz = (p_y_T - p_z_B)/Real(dx[2]);
                            }
                            else
                            {
                                const int idx_cell_rho_Y_z_B = (i + num_subghosts_conservative_var[0][0]) +
                                    (interior_box_lo_idx[1] + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (k - 1 + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_rho_Y_z_F = (i + num_subghosts_conservative_var[0][0]) +
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
                                
                                /*
                                 * Compute the mixture density.
                                 */
                                
                                Real rho_z_B = Real(0);
                                Real rho_z_F = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_z_B += Q[si][idx_cell_rho_Y_z_B];
                                    rho_z_F += Q[si][idx_cell_rho_Y_z_F];
                                }
                                
                                /*
                                 * Compute the mass fractions.
                                 */
                                
                                std::vector<Real> Y_z_B;
                                std::vector<Real> Y_z_F;
                                Y_z_B.reserve(d_num_species);
                                Y_z_F.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_z_B.push_back(Q[si][idx_cell_rho_Y_z_B]/rho_z_B);
                                    Y_z_F.push_back(Q[si][idx_cell_rho_Y_z_F]/rho_z_F);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions.
                                 */
                                
                                std::vector<const Real*> Y_z_B_ptr;
                                std::vector<const Real*> Y_z_F_ptr;
                                Y_z_B_ptr.reserve(d_num_species);
                                Y_z_F_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_z_B_ptr.push_back(&Y_z_B[si]);
                                    Y_z_F_ptr.push_back(&Y_z_F[si]);
                                }
                                
                                const Real u_z_B = Q[d_num_species][idx_cell_mom_z_B]/rho_z_B;
                                const Real u_z_F = Q[d_num_species][idx_cell_mom_z_F]/rho_z_F;
                                
                                const Real v_z_B = Q[d_num_species + 1][idx_cell_mom_z_B]/rho_z_B;
                                const Real v_z_F = Q[d_num_species + 1][idx_cell_mom_z_F]/rho_z_F;
                                
                                const Real w_z_B = Q[d_num_species + 2][idx_cell_mom_z_B]/rho_z_B;
                                const Real w_z_F = Q[d_num_species + 2][idx_cell_mom_z_F]/rho_z_F;
                                
                                const Real epsilon_z_B = Q[d_num_species + 3][idx_cell_E_z_B]/rho_z_B - half*(u_z_B*u_z_B + v_z_B*v_z_B + w_z_B*w_z_B);
                                const Real epsilon_z_F = Q[d_num_species + 3][idx_cell_E_z_F]/rho_z_F - half*(u_z_F*u_z_F + v_z_F*v_z_F + w_z_F*w_z_F);
                                
                                Real p_z_B = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_z_B,
                                        &epsilon_z_B,
                                        Y_z_B_ptr);
                                
                                Real p_z_F = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_z_F,
                                        &epsilon_z_F,
                                        Y_z_F_ptr);
                                
                                // Central derivatives.
                                // du_dz = (u_z_F - u_z_B)/(Real(2)*Real(dx[2]));
                                dv_dz = (v_z_F - v_z_B)/(Real(2)*Real(dx[2]));
                                dw_dz = (w_z_F - w_z_B)/(Real(2)*Real(dx[2]));
                                dp_dz = (p_z_F - p_z_B)/(Real(2)*Real(dx[2]));
                            }
                            
                            // Compute sound speed.
                            
                            const Real Gamma_y_T = d_equation_of_state_mixing_rules->getGruneisenParameter(
                                &rho_y_T,
                                &p_y_T,
                                Y_y_T_ptr);
                            
                            const std::vector<Real> Psi_y_T = d_equation_of_state_mixing_rules->
                                getPressureDerivativeWithPartialDensities(
                                        &rho_y_T,
                                        &p_y_T,
                                        Y_y_T_ptr);
                            
                            Real c_y_T = Gamma_y_T*p_y_T/rho_y_T;
                            for (int si = 0; si < d_num_species; si++)
                            {
                                c_y_T += Y_y_T[si]*Psi_y_T[si];
                            }
                            c_y_T = std::sqrt(c_y_T);
                            
                            const Real lambda_last = v_y_T + c_y_T;
                            
                            // Compute vector Lambda^(-1) * L.
                            
                            Real Lambda_inv_L[d_num_species + 4];
                            
                            const Real& p_t         = d_bdry_face_nonreflecting_outflow_p_t[face_loc];
                            const Real& sigma       = d_bdry_face_nonreflecting_outflow_sigma[face_loc];
                            const Real& beta        = d_bdry_face_nonreflecting_outflow_beta[face_loc];
                            const Real& length_char = d_bdry_face_nonreflecting_outflow_length_char[face_loc];
                            
                            const Real T_last = u_y_T*(dp_dx + rho_y_T*c_y_T*dv_dx) + rho_y_T*c_y_T*c_y_T*du_dx + 
                                w_y_T*(dp_dz + rho_y_T*c_y_T*dv_dz) + rho_y_T*c_y_T*c_y_T*dw_dz;
                            
                            const Real M_sq = (u_y_T*u_y_T + v_y_T*v_y_T + w_y_T*w_y_T)/(c_y_T*c_y_T);
                            const Real K = sigma*c_y_T*(Real(1) - M_sq)/length_char;
                            
                            Lambda_inv_L[0] = dp_dy - rho_y_T*c_y_T*dv_dy;
                            Lambda_inv_L[1] = du_dy;
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Lambda_inv_L[si + 2] = c_y_T*c_y_T*drho_Y_dy[si] - Y_y_T[si]*dp_dy;
                            }
                            Lambda_inv_L[d_num_species + 2] = dw_dy;
                            Lambda_inv_L[d_num_species + 3] = (Real(1)/lambda_last)*(K*(p_y_T - p_t) - (Real(1) - beta)*T_last);
                            
                            // Compute dV_dy.
                            
                            const Real c_sq_inv  = Real(1)/(c_y_T*c_y_T);
                            const Real rho_c_inv = Real(1)/(rho_y_T*c_y_T);
                            
                            Real dV_dy[d_num_species + 4];

                            for (int si = 0; si < d_num_species; si++)
                            {
                                dV_dy[si] = half*c_sq_inv*Y_y_T[si]*(Lambda_inv_L[0] + Lambda_inv_L[d_num_species + 3]) +
                                    c_sq_inv*Lambda_inv_L[si + 2];
                            }
                            dV_dy[d_num_species]     = Lambda_inv_L[1];
                            dV_dy[d_num_species + 1] = half*rho_c_inv*(-Lambda_inv_L[0] + Lambda_inv_L[d_num_species + 3]);
                            dV_dy[d_num_species + 2] = Lambda_inv_L[d_num_species + 2];
                            dV_dy[d_num_species + 3] = half*(Lambda_inv_L[0] + Lambda_inv_L[d_num_species + 3]);
                            
                            Real V_ghost[(d_num_species + 4)*num_ghosts_to_fill];
                            
                            for (int j = num_ghosts_to_fill - 1; j >= 0; j--)
                            {
                                const int idx_cell_rho_Y = (i + num_subghosts_conservative_var[0][0]) +
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
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        V_ghost[j*(d_num_species + 4) + si] = rho_Y_y_TT[si] - Real(2)*Real(dx[1])*dV_dy[si];
                                    }
                                    
                                    V_ghost[j*(d_num_species + 4) + d_num_species]     = u_y_TT - Real(2)*Real(dx[1])*dV_dy[d_num_species];
                                    V_ghost[j*(d_num_species + 4) + d_num_species + 1] = v_y_TT - Real(2)*Real(dx[1])*dV_dy[d_num_species + 1];
                                    V_ghost[j*(d_num_species + 4) + d_num_species + 2] = w_y_TT - Real(2)*Real(dx[1])*dV_dy[d_num_species + 2];
                                    V_ghost[j*(d_num_species + 4) + d_num_species + 3] = p_y_TT - Real(2)*Real(dx[1])*dV_dy[d_num_species + 3];
                                }
                                else if (j == num_ghosts_to_fill - 2)
                                {
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        V_ghost[j*(d_num_species + 4) + si] = -Real(2)*rho_Y_y_TT[si] - Real(3)*rho_Y_y_T[si] +
                                            Real(6)*V_ghost[(j + 1)*(d_num_species + 4) + si] + Real(6)*Real(dx[1])*dV_dy[si];
                                    }
                                    
                                    V_ghost[j*(d_num_species + 4) + d_num_species] = -Real(2)*u_y_TT - Real(3)*u_y_T +
                                        Real(6)*V_ghost[(j + 1)*(d_num_species + 4) + d_num_species] +
                                        Real(6)*Real(dx[1])*dV_dy[d_num_species];
                                    
                                    V_ghost[j*(d_num_species + 4) + d_num_species + 1] = -Real(2)*v_y_TT - Real(3)*v_y_T +
                                        Real(6)*V_ghost[(j + 1)*(d_num_species + 4) + d_num_species + 1] +
                                        Real(6)*Real(dx[1])*dV_dy[d_num_species + 1];
                                    
                                    V_ghost[j*(d_num_species + 4) + d_num_species + 2] = -Real(2)*w_y_TT - Real(3)*w_y_T +
                                        Real(6)*V_ghost[(j + 1)*(d_num_species + 4) + d_num_species + 2] +
                                        Real(6)*Real(dx[1])*dV_dy[d_num_species + 2];
                                    
                                    V_ghost[j*(d_num_species + 4) + d_num_species + 3] = -Real(2)*p_y_TT - Real(3)*p_y_T +
                                        Real(6)*V_ghost[(j + 1)*(d_num_species + 4) + d_num_species + 3] +
                                        Real(6)*Real(dx[1])*dV_dy[d_num_species + 3];
                                }
                                else if (j == num_ghosts_to_fill - 3)
                                {
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        V_ghost[j*(d_num_species + 4) + si] = Real(3)*rho_Y_y_TT[si] + Real(10)*rho_Y_y_T[si] -
                                            Real(18)*V_ghost[(j + 2)*(d_num_species + 4) + si] +
                                            Real(6)*V_ghost[(j + 1)*(d_num_species + 4) + si] -
                                            Real(12)*Real(dx[1])*dV_dy[si];
                                    }
                                    
                                    V_ghost[j*(d_num_species + 4) + d_num_species] = Real(3)*u_y_TT + Real(10)*u_y_T -
                                        Real(18)*V_ghost[(j + 2)*(d_num_species + 4) + d_num_species] +
                                        Real(6)*V_ghost[(j + 1)*(d_num_species + 4) + d_num_species] -
                                        Real(12)*Real(dx[1])*dV_dy[d_num_species];
                                    
                                    V_ghost[j*(d_num_species + 4) + d_num_species + 1] = Real(3)*v_y_TT + Real(10)*v_y_T -
                                        Real(18)*V_ghost[(j + 2)*(d_num_species + 4) + d_num_species + 1] +
                                        Real(6)*V_ghost[(j + 1)*(d_num_species + 4) + d_num_species + 1] -
                                        Real(12)*Real(dx[1])*dV_dy[d_num_species + 1];
                                    
                                    V_ghost[j*(d_num_species + 4) + d_num_species + 2] = Real(3)*w_y_TT + Real(10)*w_y_T -
                                        Real(18)*V_ghost[(j + 2)*(d_num_species + 4) + d_num_species + 2] +
                                        Real(6)*V_ghost[(j + 1)*(d_num_species + 4) + d_num_species + 2] -
                                        Real(12)*Real(dx[1])*dV_dy[d_num_species + 2];
                                    
                                    V_ghost[j*(d_num_species + 4) + d_num_species + 3] = Real(3)*p_y_TT + Real(10)*p_y_T -
                                        Real(18)*V_ghost[(j + 2)*(d_num_species + 4) + d_num_species + 3] +
                                        Real(6)*V_ghost[(j + 1)*(d_num_species + 4) + d_num_species + 3] -
                                        Real(12)*Real(dx[1])*dV_dy[d_num_species + 3];
                                }
                                else if (j == num_ghosts_to_fill - 4)
                                {
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        V_ghost[j*(d_num_species + 4) + si] = -Real(4)*rho_Y_y_TT[si] -
                                            Real(65)/Real(3)*rho_Y_y_T[si] +
                                            Real(40)*V_ghost[(j + 3)*(d_num_species + 4) + si] -
                                            Real(20)*V_ghost[(j + 2)*(d_num_species + 4) + si] +
                                            Real(20)/Real(3)*V_ghost[(j + 1)*(d_num_species + 4) + si] +
                                            Real(20)*Real(dx[1])*dV_dy[si];
                                    }
                                    
                                    V_ghost[j*(d_num_species + 4) + d_num_species] = -Real(4)*u_y_TT -
                                        Real(65)/Real(3)*u_y_T +
                                        Real(40)*V_ghost[(j + 3)*(d_num_species + 4) + d_num_species] -
                                        Real(20)*V_ghost[(j + 2)*(d_num_species + 4) + d_num_species] +
                                        Real(20)/Real(3)*V_ghost[(j + 1)*(d_num_species + 4) + d_num_species] +
                                        Real(20)*Real(dx[1])*dV_dy[d_num_species];
                                    
                                    V_ghost[j*(d_num_species + 4) + d_num_species + 1] = -Real(4)*v_y_TT -
                                        Real(65)/Real(3)*v_y_T +
                                        Real(40)*V_ghost[(j + 3)*(d_num_species + 4) + d_num_species + 1] -
                                        Real(20)*V_ghost[(j + 2)*(d_num_species + 4) + d_num_species + 1] +
                                        Real(20)/Real(3)*V_ghost[(j + 1)*(d_num_species + 4) + d_num_species + 1] +
                                        Real(20)*Real(dx[1])*dV_dy[d_num_species + 1];
                                    
                                    V_ghost[j*(d_num_species + 4) + d_num_species + 2] = -Real(4)*w_y_TT -
                                        Real(65)/Real(3)*w_y_T +
                                        Real(40)*V_ghost[(j + 3)*(d_num_species + 4) + d_num_species + 2] -
                                        Real(20)*V_ghost[(j + 2)*(d_num_species + 4) + d_num_species + 2] +
                                        Real(20)/Real(3)*V_ghost[(j + 1)*(d_num_species + 4) + d_num_species + 2] +
                                        Real(20)*Real(dx[1])*dV_dy[d_num_species + 2];
                                    
                                    V_ghost[j*(d_num_species + 4) + d_num_species + 3] = -Real(4)*p_y_TT -
                                        Real(65)/Real(3)*p_y_T +
                                        Real(40)*V_ghost[(j + 3)*(d_num_species + 4) + d_num_species + 3] -
                                        Real(20)*V_ghost[(j + 2)*(d_num_species + 4) + d_num_species + 3] +
                                        Real(20)/Real(3)*V_ghost[(j + 1)*(d_num_species + 4) + d_num_species + 3] +
                                        Real(20)*Real(dx[1])*dV_dy[d_num_species + 3];
                                }
                                else if (j == num_ghosts_to_fill - 5)
                                {
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        V_ghost[j*(d_num_species + 4) + si] = Real(5)*rho_Y_y_TT[si] +
                                            Real(77)/Real(2)*rho_Y_y_T[si] -
                                            Real(75)*V_ghost[(j + 4)*(d_num_species + 4) + si] +
                                            Real(50)*V_ghost[(j + 3)*(d_num_species + 4) + si] -
                                            Real(25)*V_ghost[(j + 2)*(d_num_species + 4) + si] +
                                            Real(15)/Real(2)*V_ghost[(j + 1)*(d_num_species + 4) + si] -
                                            Real(30)*Real(dx[1])*dV_dy[si];
                                    }
                                    
                                    V_ghost[j*(d_num_species + 4) + d_num_species] = Real(5)*u_y_TT +
                                        Real(77)/Real(2)*u_y_T -
                                        Real(75)*V_ghost[(j + 4)*(d_num_species + 4) + d_num_species] +
                                        Real(50)*V_ghost[(j + 3)*(d_num_species + 4) + d_num_species] -
                                        Real(25)*V_ghost[(j + 2)*(d_num_species + 4) + d_num_species] +
                                        Real(15)/Real(2)*V_ghost[(j + 1)*(d_num_species + 4) + d_num_species] -
                                        Real(30)*Real(dx[1])*dV_dy[d_num_species];
                                    
                                    V_ghost[j*(d_num_species + 4) + d_num_species + 1] = Real(5)*v_y_TT +
                                        Real(77)/Real(2)*v_y_T -
                                        Real(75)*V_ghost[(j + 4)*(d_num_species + 4) + d_num_species + 1] +
                                        Real(50)*V_ghost[(j + 3)*(d_num_species + 4) + d_num_species + 1] -
                                        Real(25)*V_ghost[(j + 2)*(d_num_species + 4) + d_num_species + 1] +
                                        Real(15)/Real(2)*V_ghost[(j + 1)*(d_num_species + 4) + d_num_species + 1] -
                                        Real(30)*Real(dx[1])*dV_dy[d_num_species + 1];
                                    
                                    V_ghost[j*(d_num_species + 4) + d_num_species + 2] = Real(5)*w_y_TT +
                                        Real(77)/Real(2)*w_y_T -
                                        Real(75)*V_ghost[(j + 4)*(d_num_species + 4) + d_num_species + 2] +
                                        Real(50)*V_ghost[(j + 3)*(d_num_species + 4) + d_num_species + 2] -
                                        Real(25)*V_ghost[(j + 2)*(d_num_species + 4) + d_num_species + 2] +
                                        Real(15)/Real(2)*V_ghost[(j + 1)*(d_num_species + 4) + d_num_species + 2] -
                                        Real(30)*Real(dx[1])*dV_dy[d_num_species + 2];
                                    
                                    V_ghost[j*(d_num_species + 4) + d_num_species + 3] = Real(5)*p_y_TT +
                                        Real(77)/Real(2)*p_y_T -
                                        Real(75)*V_ghost[(j + 4)*(d_num_species + 4) + d_num_species + 3] +
                                        Real(50)*V_ghost[(j + 3)*(d_num_species + 4) + d_num_species + 3] -
                                        Real(25)*V_ghost[(j + 2)*(d_num_species + 4) + d_num_species + 3] +
                                        Real(15)/Real(2)*V_ghost[(j + 1)*(d_num_species + 4) + d_num_species + 3] -
                                        Real(30)*Real(dx[1])*dV_dy[d_num_species + 3];
                                }
                                else if (j == num_ghosts_to_fill - 6)
                                {
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        V_ghost[j*(d_num_species + 4) + si] = -Real(6)*rho_Y_y_TT[si] -
                                            Real(609)/Real(10)*rho_Y_y_T[si] +
                                            Real(126)*V_ghost[(j + 5)*(d_num_species + 4) + si] -
                                            Real(105)*V_ghost[(j + 4)*(d_num_species + 4) + si] +
                                            Real(70)*V_ghost[(j + 3)*(d_num_species + 4) + si] -
                                            Real(63)/Real(2)*V_ghost[(j + 2)*(d_num_species + 4) + si] +
                                            Real(42)/Real(5)*V_ghost[(j + 1)*(d_num_species + 4) + si] +
                                            Real(42)*Real(dx[1])*dV_dy[si];
                                    }
                                    
                                    V_ghost[j*(d_num_species + 4) + d_num_species] = -Real(6)*u_y_TT -
                                        Real(609)/Real(10)*u_y_T +
                                        Real(126)*V_ghost[(j + 5)*(d_num_species + 4) + d_num_species] -
                                        Real(105)*V_ghost[(j + 4)*(d_num_species + 4) + d_num_species] +
                                        Real(70)*V_ghost[(j + 3)*(d_num_species + 4) + d_num_species] -
                                        Real(63)/Real(2)*V_ghost[(j + 2)*(d_num_species + 4) + d_num_species] +
                                        Real(42)/Real(5)*V_ghost[(j + 1)*(d_num_species + 4) + d_num_species] +
                                        Real(42)*Real(dx[1])*dV_dy[d_num_species];
                                    
                                    V_ghost[j*(d_num_species + 4) + d_num_species + 1] = -Real(6)*v_y_TT -
                                        Real(609)/Real(10)*v_y_T +
                                        Real(126)*V_ghost[(j + 5)*(d_num_species + 4) + d_num_species + 1] -
                                        Real(105)*V_ghost[(j + 4)*(d_num_species + 4) + d_num_species + 1] +
                                        Real(70)*V_ghost[(j + 3)*(d_num_species + 4) + d_num_species + 1] -
                                        Real(63)/Real(2)*V_ghost[(j + 2)*(d_num_species + 4) + d_num_species + 1] +
                                        Real(42)/Real(5)*V_ghost[(j + 1)*(d_num_species + 4) + d_num_species + 1] +
                                        Real(42)*Real(dx[1])*dV_dy[d_num_species + 1];
                                    
                                    V_ghost[j*(d_num_species + 4) + d_num_species + 2] = -Real(6)*w_y_TT -
                                        Real(609)/Real(10)*w_y_T +
                                        Real(126)*V_ghost[(j + 5)*(d_num_species + 4) + d_num_species + 2] -
                                        Real(105)*V_ghost[(j + 4)*(d_num_species + 4) + d_num_species + 2] +
                                        Real(70)*V_ghost[(j + 3)*(d_num_species + 4) + d_num_species + 2] -
                                        Real(63)/Real(2)*V_ghost[(j + 2)*(d_num_species + 4) + d_num_species + 2] +
                                        Real(42)/Real(5)*V_ghost[(j + 1)*(d_num_species + 4) + d_num_species + 2] +
                                        Real(42)*Real(dx[1])*dV_dy[d_num_species + 2];
                                    
                                    V_ghost[j*(d_num_species + 4) + d_num_species + 3] = -Real(6)*p_y_TT -
                                        Real(609)/Real(10)*p_y_T +
                                        Real(126)*V_ghost[(j + 5)*(d_num_species + 4) + d_num_species + 3] -
                                        Real(105)*V_ghost[(j + 4)*(d_num_species + 4) + d_num_species + 3] +
                                        Real(70)*V_ghost[(j + 3)*(d_num_species + 4) + d_num_species + 3] -
                                        Real(63)/Real(2)*V_ghost[(j + 2)*(d_num_species + 4) + d_num_species + 3] +
                                        Real(42)/Real(5)*V_ghost[(j + 1)*(d_num_species + 4) + d_num_species + 3] +
                                        Real(42)*Real(dx[1])*dV_dy[d_num_species + 3];
                                }
                                
                                /*
                                 * Compute the mixture density.
                                 */
                                
                                Real rho_ghost = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_ghost += V_ghost[j*(d_num_species + 4) + si];
                                }
                                
                                /*
                                 * Compute the mass fractions.
                                 */
                                
                                std::vector<Real> Y_ghost;
                                Y_ghost.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_ghost.push_back(V_ghost[j*(d_num_species + 4) + si]/rho_ghost);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions.
                                 */
                                
                                std::vector<const Real*> Y_ghost_ptr;
                                Y_ghost_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_ghost_ptr.push_back(&Y_ghost[si]);
                                }
                                
                                for(int si=0; si < d_num_species; si++)
                                {
                                    Q[si][idx_cell_rho_Y] = V_ghost[j*(d_num_species + 4) + si];
                                }
                                
                                Q[d_num_species][idx_cell_mom]     = rho_ghost*V_ghost[j*(d_num_species + 4) + d_num_species];
                                Q[d_num_species + 1][idx_cell_mom] = rho_ghost*V_ghost[j*(d_num_species + 4) + d_num_species + 1];
                                Q[d_num_species + 2][idx_cell_mom] = rho_ghost*V_ghost[j*(d_num_species + 4) + d_num_species + 2];
                                
                                const Real epsilon = d_equation_of_state_mixing_rules->
                                    getInternalEnergy(
                                        &rho_ghost,
                                        &V_ghost[j*(d_num_species + 4) + d_num_species + 3],
                                        Y_ghost_ptr);
                                
                                const Real E = rho_ghost*epsilon +
                                    half*(Q[d_num_species][idx_cell_mom]*Q[d_num_species][idx_cell_mom] +
                                        Q[d_num_species + 1][idx_cell_mom]*Q[d_num_species + 1][idx_cell_mom] +
                                        Q[d_num_species + 2][idx_cell_mom]*Q[d_num_species + 2][idx_cell_mom])/
                                        rho_ghost;
                                
                                Q[d_num_species + 3][idx_cell_E] = E;
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
                            << ": FlowModelBoundaryUtilitiesFourEqnConservative::fill3dFaceBoundaryData()\n"
                            << "Non-reflecting outflow BC doesn't support more than six ghost cells yet!");
                    }
                    
                    for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                    {
                        for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                        {
                            // Get the grid spacing.
                            const double* const dx = patch_geom->getDx();
                            
                            const int idx_cell_rho_Y_y_B = (i + num_subghosts_conservative_var[0][0]) +
                                (interior_box_hi_idx[1] + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                (k + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                    subghostcell_dims_conservative_var[0][1];
                            
                            const int idx_cell_rho_Y_y_BB = (i + num_subghosts_conservative_var[0][0]) +
                                (interior_box_hi_idx[1] - 1 + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                (k + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                    subghostcell_dims_conservative_var[0][1];
                            
                            const int idx_cell_rho_Y_y_BBB = (i + num_subghosts_conservative_var[0][0]) +
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
                            
                            std::vector<Real> rho_Y_y_B;
                            std::vector<Real> rho_Y_y_BB;
                            std::vector<Real> rho_Y_y_BBB;
                            rho_Y_y_B.reserve(d_num_species);
                            rho_Y_y_BB.reserve(d_num_species);
                            rho_Y_y_BBB.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                rho_Y_y_B.push_back(Q[si][idx_cell_rho_Y_y_B]);
                                rho_Y_y_BB.push_back(Q[si][idx_cell_rho_Y_y_BB]);
                                rho_Y_y_BBB.push_back(Q[si][idx_cell_rho_Y_y_BBB]);
                            }
                            
                            /*
                             * Compute the mixture density.
                             */
                            
                            Real rho_y_B   = Real(0);
                            Real rho_y_BB  = Real(0);
                            Real rho_y_BBB = Real(0);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                rho_y_B   += Q[si][idx_cell_rho_Y_y_B];
                                rho_y_BB  += Q[si][idx_cell_rho_Y_y_BB];
                                rho_y_BBB += Q[si][idx_cell_rho_Y_y_BBB];
                            }
                            
                            /*
                             * Compute the mass fractions.
                             */
                            
                            std::vector<Real> Y_y_B;
                            std::vector<Real> Y_y_BB;
                            std::vector<Real> Y_y_BBB;
                            Y_y_B.reserve(d_num_species);
                            Y_y_BB.reserve(d_num_species);
                            Y_y_BBB.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_y_B.push_back(Q[si][idx_cell_rho_Y_y_B]/rho_y_B);
                                Y_y_BB.push_back(Q[si][idx_cell_rho_Y_y_BB]/rho_y_BB);
                                Y_y_BBB.push_back(Q[si][idx_cell_rho_Y_y_BBB]/rho_y_BBB);
                            }
                            
                            /*
                             * Get the pointers to the mass fractions.
                             */
                            
                            std::vector<const Real*> Y_y_B_ptr;
                            std::vector<const Real*> Y_y_BB_ptr;
                            std::vector<const Real*> Y_y_BBB_ptr;
                            Y_y_B_ptr.reserve(d_num_species);
                            Y_y_BB_ptr.reserve(d_num_species);
                            Y_y_BBB_ptr.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_y_B_ptr.push_back(&Y_y_B[si]);
                                Y_y_BB_ptr.push_back(&Y_y_BB[si]);
                                Y_y_BBB_ptr.push_back(&Y_y_BBB[si]);
                            }
                            
                            const Real u_y_B   = Q[d_num_species][idx_cell_mom_y_B]/rho_y_B;
                            const Real u_y_BB  = Q[d_num_species][idx_cell_mom_y_BB]/rho_y_BB;
                            const Real u_y_BBB = Q[d_num_species][idx_cell_mom_y_BBB]/rho_y_BBB;
                            
                            const Real v_y_B   = Q[d_num_species + 1][idx_cell_mom_y_B]/rho_y_B;
                            const Real v_y_BB  = Q[d_num_species + 1][idx_cell_mom_y_BB]/rho_y_BB;
                            const Real v_y_BBB = Q[d_num_species + 1][idx_cell_mom_y_BBB]/rho_y_BBB;
                            
                            const Real w_y_B   = Q[d_num_species + 2][idx_cell_mom_y_B]/rho_y_B;
                            const Real w_y_BB  = Q[d_num_species + 2][idx_cell_mom_y_BB]/rho_y_BB;
                            const Real w_y_BBB = Q[d_num_species + 2][idx_cell_mom_y_BBB]/rho_y_BBB;
                            
                            const Real half = Real(1)/Real(2);
                            const Real epsilon_y_B   = Q[d_num_species + 3][idx_cell_E_y_B]/rho_y_B - half*(u_y_B*u_y_B + v_y_B*v_y_B + w_y_B*w_y_B);
                            const Real epsilon_y_BB  = Q[d_num_species + 3][idx_cell_E_y_BB]/rho_y_BB - half*(u_y_BB*u_y_BB + v_y_BB*v_y_BB + w_y_BB*w_y_BB);
                            const Real epsilon_y_BBB = Q[d_num_species + 3][idx_cell_E_y_BBB]/rho_y_BBB - half*(u_y_BBB*u_y_BBB + v_y_BBB*v_y_BBB + w_y_BBB*w_y_BBB);

                            Real p_y_B = d_equation_of_state_mixing_rules->
                                getPressure(
                                    &rho_y_B,
                                    &epsilon_y_B,
                                    Y_y_B_ptr);
                            
                            Real p_y_BB = d_equation_of_state_mixing_rules->
                                getPressure(
                                    &rho_y_BB,
                                    &epsilon_y_BB,
                                    Y_y_BB_ptr);
                            
                            Real p_y_BBB = d_equation_of_state_mixing_rules->
                                getPressure(
                                    &rho_y_BBB,
                                    &epsilon_y_BBB,
                                    Y_y_BBB_ptr);
                            
                            /*
                             * Compute derivatives in y-direction.
                             */
                            
                            std::vector<Real> drho_Y_dy;
                            drho_Y_dy.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                drho_Y_dy.push_back((Q[si][idx_cell_rho_Y_y_BBB] - Real(4)*Q[si][idx_cell_rho_Y_y_BB] +
                                    Real(3)*Q[si][idx_cell_rho_Y_y_B])/(Real(2)*Real(dx[1])));
                            }
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
                                
                                const int idx_cell_rho_Y_x_R = (i + 1 + num_subghosts_conservative_var[0][0]) +
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
                                
                                /*
                                 * Compute the mixture density.
                                 */
                                
                                Real rho_x_R = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_x_R += Q[si][idx_cell_rho_Y_x_R];
                                }
                                
                                /*
                                 * Compute the mass fractions.
                                 */
                                
                                std::vector<Real> Y_x_R;
                                Y_x_R.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_x_R.push_back(Q[si][idx_cell_rho_Y_x_R]/rho_x_R);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions.
                                 */
                                
                                std::vector<const Real*> Y_x_R_ptr;
                                Y_x_R_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_x_R_ptr.push_back(&Y_x_R[si]);
                                }
                                
                                const Real u_x_R = Q[d_num_species][idx_cell_mom_x_R]/rho_x_R;
                                const Real v_x_R = Q[d_num_species + 1][idx_cell_mom_x_R]/rho_x_R;
                                const Real w_x_R = Q[d_num_species + 2][idx_cell_mom_x_R]/rho_x_R;
                                const Real epsilon_x_R = Q[d_num_species + 3][idx_cell_E_x_R]/rho_x_R - half*(u_x_R*u_x_R + v_x_R*v_x_R + w_x_R*w_x_R);
                                
                                Real p_x_R = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_x_R,
                                        &epsilon_x_R,
                                        Y_x_R_ptr);
                                
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
                                
                                const int idx_cell_rho_Y_x_L = (i - 1 + num_subghosts_conservative_var[0][0]) +
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
                                
                                /*
                                 * Compute the mixture density.
                                 */
                                
                                Real rho_x_L = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_x_L += Q[si][idx_cell_rho_Y_x_L];
                                }
                                
                                /*
                                 * Compute the mass fractions.
                                 */
                                
                                std::vector<Real> Y_x_L;
                                Y_x_L.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_x_L.push_back(Q[si][idx_cell_rho_Y_x_L]/rho_x_L);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions.
                                 */
                            
                                std::vector<const Real*> Y_x_L_ptr;
                                Y_x_L_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_x_L_ptr.push_back(&Y_x_L[si]);
                                }
                                
                                const Real u_x_L = Q[d_num_species][idx_cell_mom_x_L]/rho_x_L;
                                const Real v_x_L = Q[d_num_species + 1][idx_cell_mom_x_L]/rho_x_L;
                                const Real w_x_L = Q[d_num_species + 2][idx_cell_mom_x_L]/rho_x_L;
                                const Real epsilon_x_L = Q[d_num_species + 3][idx_cell_E_x_L]/rho_x_L - half*(u_x_L*u_x_L + v_x_L*v_x_L + w_x_L*w_x_L);
                                
                                Real p_x_L = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_x_L,
                                        &epsilon_x_L,
                                        Y_x_L_ptr);
                                
                                // One-sided derivatives.
                                du_dx = (u_y_B - u_x_L)/Real(dx[0]);
                                dv_dx = (v_y_B - v_x_L)/Real(dx[0]);
                                // dw_dx = (w_y_B - w_x_L)/Real(dx[0]);
                                dp_dx = (p_y_B - p_x_L)/Real(dx[0]);
                            }
                            else
                            {
                                const int idx_cell_rho_Y_x_L = (i - 1 + num_subghosts_conservative_var[0][0]) +
                                    (interior_box_hi_idx[1] + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (k + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_rho_Y_x_R = (i + 1 + num_subghosts_conservative_var[0][0]) +
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
                                
                                /*
                                 * Compute the mixture density.
                                 */
                                
                                Real rho_x_L = Real(0);
                                Real rho_x_R = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_x_L += Q[si][idx_cell_rho_Y_x_L];
                                    rho_x_R += Q[si][idx_cell_rho_Y_x_R];
                                }
                                
                                /*
                                 * Compute the mass fractions.
                                 */
                                
                                std::vector<Real> Y_x_L;
                                std::vector<Real> Y_x_R;
                                Y_x_L.reserve(d_num_species);
                                Y_x_R.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_x_L.push_back(Q[si][idx_cell_rho_Y_x_L]/rho_x_L);
                                    Y_x_R.push_back(Q[si][idx_cell_rho_Y_x_R]/rho_x_R);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions.
                                 */
                                
                                std::vector<const Real*> Y_x_L_ptr;
                                std::vector<const Real*> Y_x_R_ptr;
                                                                    
                                Y_x_L_ptr.reserve(d_num_species);
                                Y_x_R_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_x_L_ptr.push_back(&Y_x_L[si]);
                                    Y_x_R_ptr.push_back(&Y_x_R[si]);
                                }
                                
                                const Real u_x_L = Q[d_num_species][idx_cell_mom_x_L]/rho_x_L;
                                const Real u_x_R = Q[d_num_species][idx_cell_mom_x_R]/rho_x_R;

                                const Real v_x_L = Q[d_num_species + 1][idx_cell_mom_x_L]/rho_x_L;
                                const Real v_x_R = Q[d_num_species + 1][idx_cell_mom_x_R]/rho_x_R;

                                const Real w_x_L = Q[d_num_species + 2][idx_cell_mom_x_L]/rho_x_L;
                                const Real w_x_R = Q[d_num_species + 2][idx_cell_mom_x_R]/rho_x_R;
                                
                                const Real epsilon_x_L = Q[d_num_species + 3][idx_cell_E_x_L]/rho_x_L - half*(u_x_L*u_x_L + v_x_L*v_x_L + w_x_L*w_x_L);
                                const Real epsilon_x_R = Q[d_num_species + 3][idx_cell_E_x_R]/rho_x_R - half*(u_x_R*u_x_R + v_x_R*v_x_R + w_x_R*w_x_R);
                                
                                Real p_x_L = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_x_L,
                                        &epsilon_x_L,
                                        Y_x_L_ptr);
                            
                                Real p_x_R = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_x_R,
                                        &epsilon_x_R,
                                        Y_x_R_ptr);

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
                            
                            if (((k + num_subghosts_conservative_var[0][2] == 0) ||
                                (k + num_subghosts_conservative_var[1][2] == 0) ||
                                (k + num_subghosts_conservative_var[2][2] == 0)))
                            // if (((patch_geom->getTouchesRegularBoundary(2, 0)) && (k == interior_box_lo_idx[2])) ||
                            //     ((k + num_subghosts_conservative_var[0][2] == 0) ||
                            //      (k + num_subghosts_conservative_var[1][2] == 0) ||
                            //      (k + num_subghosts_conservative_var[2][2] == 0)))
                            {
                                // Patch is touching periodic boundary.
                                
                                const int idx_cell_rho_Y_z_F = (i + num_subghosts_conservative_var[0][0]) +
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
                                
                                /*
                                 * Compute the mixture density.
                                 */
                                
                                Real rho_z_F = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_z_F += Q[si][idx_cell_rho_Y_z_F];
                                }
                                
                                /*
                                 * Compute the mass fractions.
                                 */
                                
                                std::vector<Real> Y_z_F;
                                Y_z_F.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_z_F.push_back(Q[si][idx_cell_rho_Y_z_F]/rho_z_F);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions.
                                 */
                                
                                std::vector<const Real*> Y_z_F_ptr;
                                Y_z_F_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_z_F_ptr.push_back(&Y_z_F[si]);
                                }
                                
                                const Real u_z_F = Q[d_num_species][idx_cell_mom_z_F]/rho_z_F;
                                const Real v_z_F = Q[d_num_species + 1][idx_cell_mom_z_F]/rho_z_F;
                                const Real w_z_F = Q[d_num_species + 2][idx_cell_mom_z_F]/rho_z_F;
                                const Real epsilon_z_F = Q[d_num_species + 3][idx_cell_E_z_F]/rho_z_F - half*(u_z_F*u_z_F + v_z_F*v_z_F + w_z_F*w_z_F);
                                
                                Real p_z_F = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_z_F,
                                        &epsilon_z_F,
                                        Y_z_F_ptr);
                                
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
                                
                                const int idx_cell_rho_Y_z_B = (i + num_subghosts_conservative_var[0][0]) +
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
                                
                                /*
                                 * Compute the mixture density.
                                 */
                                
                                Real rho_z_B = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_z_B += Q[si][idx_cell_rho_Y_z_B];
                                }
                                
                                /*
                                 * Compute the mass fractions.
                                 */
                                
                                std::vector<Real> Y_z_B;
                                Y_z_B.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_z_B.push_back(Q[si][idx_cell_rho_Y_z_B]/rho_z_B);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions.
                                 */
                                
                                std::vector<const Real*> Y_z_B_ptr;
                                Y_z_B_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_z_B_ptr.push_back(&Y_z_B[si]);
                                }
                                
                                const Real u_z_B = Q[d_num_species][idx_cell_mom_z_B]/rho_z_B;
                                const Real v_z_B = Q[d_num_species + 1][idx_cell_mom_z_B]/rho_z_B;
                                const Real w_z_B = Q[d_num_species + 2][idx_cell_mom_z_B]/rho_z_B;
                                const Real epsilon_z_B = Q[d_num_species + 3][idx_cell_E_z_B]/rho_z_B - half*(u_z_B*u_z_B + v_z_B*v_z_B + w_z_B*w_z_B);
                                
                                Real p_z_B = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_z_B,
                                        &epsilon_z_B,
                                        Y_z_B_ptr);
                                
                                // One-sided derivatives.
                                // du_dz = (u_y_B - u_z_B)/Real(dx[2]);
                                dv_dz = (v_y_B - v_z_B)/Real(dx[2]);
                                dw_dz = (w_y_B - w_z_B)/Real(dx[2]);
                                dp_dz = (p_y_B - p_z_B)/Real(dx[2]);
                            }
                            else
                            {
                                const int idx_cell_rho_Y_z_B = (i + num_subghosts_conservative_var[0][0]) +
                                    (interior_box_hi_idx[1] + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (k - 1 + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_rho_Y_z_F = (i + num_subghosts_conservative_var[0][0]) +
                                    (interior_box_hi_idx[1] + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (k + 1 + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom_z_B = (i + num_subghosts_conservative_var[1][0]) +
                                    (interior_box_hi_idx[1] + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                    (k - 1 + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_mom_z_F = (i + num_subghosts_conservative_var[1][0]) +
                                    (interior_box_hi_idx[1] + num_subghosts_conservative_var[1][1])*subghostcell_dims_conservative_var[1][0] +
                                    (k + 1 + num_subghosts_conservative_var[1][2])*subghostcell_dims_conservative_var[1][0]*
                                        subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E_z_B = (i + num_subghosts_conservative_var[2][0]) +
                                    (interior_box_hi_idx[1] + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                    (k - 1 + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                const int idx_cell_E_z_F = (i + num_subghosts_conservative_var[2][0]) +
                                    (interior_box_hi_idx[1] + num_subghosts_conservative_var[2][1])*subghostcell_dims_conservative_var[2][0] +
                                    (k + 1 + num_subghosts_conservative_var[2][2])*subghostcell_dims_conservative_var[2][0]*
                                        subghostcell_dims_conservative_var[2][1];
                                
                                /*
                                 * Compute the mixture density.
                                 */
                                
                                Real rho_z_B = Real(0);
                                Real rho_z_F = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_z_B += Q[si][idx_cell_rho_Y_z_B];
                                    rho_z_F += Q[si][idx_cell_rho_Y_z_F];
                                }
                                
                                /*
                                 * Compute the mass fractions.
                                 */
                                
                                std::vector<Real> Y_z_B;
                                std::vector<Real> Y_z_F;
                                Y_z_B.reserve(d_num_species);
                                Y_z_F.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_z_B.push_back(Q[si][idx_cell_rho_Y_z_B]/rho_z_B);
                                    Y_z_F.push_back(Q[si][idx_cell_rho_Y_z_F]/rho_z_F);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions.
                                 */
                                
                                std::vector<const Real*> Y_z_B_ptr;
                                std::vector<const Real*> Y_z_F_ptr;
                                Y_z_B_ptr.reserve(d_num_species);
                                Y_z_F_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_z_B_ptr.push_back(&Y_z_B[si]);
                                    Y_z_F_ptr.push_back(&Y_z_F[si]);
                                }
                                
                                const Real u_z_B = Q[d_num_species][idx_cell_mom_z_B]/rho_z_B;
                                const Real u_z_F = Q[d_num_species][idx_cell_mom_z_F]/rho_z_F;
                                
                                const Real v_z_B = Q[d_num_species + 1][idx_cell_mom_z_B]/rho_z_B;
                                const Real v_z_F = Q[d_num_species + 1][idx_cell_mom_z_F]/rho_z_F;
                                
                                const Real w_z_B = Q[d_num_species + 2][idx_cell_mom_z_B]/rho_z_B;
                                const Real w_z_F = Q[d_num_species + 2][idx_cell_mom_z_F]/rho_z_F;
                                
                                const Real epsilon_z_B = Q[d_num_species + 3][idx_cell_E_z_B]/rho_z_B - half*(u_z_B*u_z_B + v_z_B*v_z_B + w_z_B*w_z_B);
                                const Real epsilon_z_F = Q[d_num_species + 3][idx_cell_E_z_F]/rho_z_F - half*(u_z_F*u_z_F + v_z_F*v_z_F + w_z_F*w_z_F);
                                
                                Real p_z_B = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_z_B,
                                        &epsilon_z_B,
                                        Y_z_B_ptr);
                                
                                Real p_z_F = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_z_F,
                                        &epsilon_z_F,
                                        Y_z_F_ptr);
                                
                                // Central derivatives.
                                // du_dz = (u_z_F - u_z_B)/(Real(2)*Real(dx[2]));
                                dv_dz = (v_z_F - v_z_B)/(Real(2)*Real(dx[2]));
                                dw_dz = (w_z_F - w_z_B)/(Real(2)*Real(dx[2]));
                                dp_dz = (p_z_F - p_z_B)/(Real(2)*Real(dx[2]));
                            }
                            
                            // Compute sound speed.
                            
                            const Real Gamma_y_B = d_equation_of_state_mixing_rules->getGruneisenParameter(
                                &rho_y_B,
                                &p_y_B,
                                Y_y_B_ptr);
                            
                            const std::vector<Real> Psi_y_B = d_equation_of_state_mixing_rules->
                                getPressureDerivativeWithPartialDensities(
                                        &rho_y_B,
                                        &p_y_B,
                                        Y_y_B_ptr);
                            
                            Real c_y_B = Gamma_y_B*p_y_B/rho_y_B;
                            for (int si = 0; si < d_num_species; si++)
                            {
                                c_y_B += Y_y_B[si]*Psi_y_B[si];
                            }
                            c_y_B = std::sqrt(c_y_B);
                            
                            const Real lambda_1 = v_y_B - c_y_B;

                            // Compute vector Lambda^(-1) * L.
                            
                            Real Lambda_inv_L[d_num_species + 4];
                            
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
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Lambda_inv_L[si + 2] = c_y_B*c_y_B*drho_Y_dy[si] - Y_y_B[si]*dp_dy;
                            }
                            Lambda_inv_L[d_num_species + 2] = dw_dy;
                            Lambda_inv_L[d_num_species + 3] = dp_dy + rho_y_B*c_y_B*dv_dy;
                            
                            // Compute dV_dy.
                            
                            const Real c_sq_inv  = Real(1)/(c_y_B*c_y_B);
                            const Real rho_c_inv = Real(1)/(rho_y_B*c_y_B);
                            
                            Real dV_dy[d_num_species + 4];
                            
                            for (int si = 0; si < d_num_species; si++)
                            {
                                dV_dy[si] = half*c_sq_inv*Y_y_B[si]*(Lambda_inv_L[0] + Lambda_inv_L[d_num_species + 3]) +
                                    c_sq_inv*Lambda_inv_L[si + 2];
                            }
                            dV_dy[d_num_species]     = Lambda_inv_L[1];
                            dV_dy[d_num_species + 1] = half*rho_c_inv*(-Lambda_inv_L[0] + Lambda_inv_L[d_num_species + 3]);
                            dV_dy[d_num_species + 2] = Lambda_inv_L[d_num_species + 2];
                            dV_dy[d_num_species + 3] = half*(Lambda_inv_L[0] + Lambda_inv_L[d_num_species + 3]);
                            
                            Real V_ghost[(d_num_species + 4)*num_ghosts_to_fill];
                            
                            for (int j = 0; j < num_ghosts_to_fill; j++)
                            {
                                const int idx_cell_rho_Y = (i + num_subghosts_conservative_var[0][0]) +
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
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        V_ghost[j*(d_num_species + 4) + si] = rho_Y_y_BB[si] + Real(2)*Real(dx[1])*dV_dy[si];
                                    }
                                    
                                    V_ghost[j*(d_num_species + 4) + d_num_species]     = u_y_BB + Real(2)*Real(dx[1])*dV_dy[d_num_species];
                                    V_ghost[j*(d_num_species + 4) + d_num_species + 1] = v_y_BB + Real(2)*Real(dx[1])*dV_dy[d_num_species + 1];
                                    V_ghost[j*(d_num_species + 4) + d_num_species + 2] = w_y_BB + Real(2)*Real(dx[1])*dV_dy[d_num_species + 2];
                                    V_ghost[j*(d_num_species + 4) + d_num_species + 3] = p_y_BB + Real(2)*Real(dx[1])*dV_dy[d_num_species + 3];
                                }
                                else if (j == 1)
                                {
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        V_ghost[j*(d_num_species + 4) + si] = -Real(2)*rho_Y_y_BB[si] - Real(3)*rho_Y_y_B[si] +
                                            Real(6)*V_ghost[(j - 1)*(d_num_species + 4) + si] - Real(6)*Real(dx[1])*dV_dy[si];
                                    }
                                    
                                    V_ghost[j*(d_num_species + 4) + d_num_species] = -Real(2)*u_y_BB - Real(3)*u_y_B +
                                        Real(6)*V_ghost[(j - 1)*(d_num_species + 4) + d_num_species] -
                                        Real(6)*Real(dx[1])*dV_dy[d_num_species];
                                    
                                    V_ghost[j*(d_num_species + 4) + d_num_species + 1] = -Real(2)*v_y_BB - Real(3)*v_y_B +
                                        Real(6)*V_ghost[(j - 1)*(d_num_species + 4) + d_num_species + 1] -
                                        Real(6)*Real(dx[1])*dV_dy[d_num_species + 1];
                                    
                                    V_ghost[j*(d_num_species + 4) + d_num_species + 2] = -Real(2)*w_y_BB - Real(3)*w_y_B +
                                        Real(6)*V_ghost[(j - 1)*(d_num_species + 4) + d_num_species + 2] -
                                        Real(6)*Real(dx[1])*dV_dy[d_num_species + 2];
                                    
                                    V_ghost[j*(d_num_species + 4) + d_num_species + 3] = -Real(2)*p_y_BB - Real(3)*p_y_B +
                                        Real(6)*V_ghost[(j - 1)*(d_num_species + 4) + d_num_species + 3] -
                                        Real(6)*Real(dx[1])*dV_dy[d_num_species + 3];
                                }
                                else if (j == 2)
                                {
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        V_ghost[j*(d_num_species + 4) + si] = Real(3)*rho_Y_y_BB[si] + Real(10)*rho_Y_y_B[si] -
                                            Real(18)*V_ghost[(j - 2)*(d_num_species + 4) + si] +
                                            Real(6)*V_ghost[(j - 1)*(d_num_species + 4) + si] +
                                            Real(12)*Real(dx[1])*dV_dy[si];
                                    }
                                    
                                    V_ghost[j*(d_num_species + 4) + d_num_species] = Real(3)*u_y_BB + Real(10)*u_y_B -
                                        Real(18)*V_ghost[(j - 2)*(d_num_species + 4) + d_num_species] +
                                        Real(6)*V_ghost[(j - 1)*(d_num_species + 4) + d_num_species] +
                                        Real(12)*Real(dx[1])*dV_dy[d_num_species];
                                    
                                    V_ghost[j*(d_num_species + 4) + d_num_species + 1] = Real(3)*v_y_BB + Real(10)*v_y_B -
                                        Real(18)*V_ghost[(j - 2)*(d_num_species + 4) + d_num_species + 1] +
                                        Real(6)*V_ghost[(j - 1)*(d_num_species + 4) + d_num_species + 1] +
                                        Real(12)*Real(dx[1])*dV_dy[d_num_species + 1];
                                    
                                    V_ghost[j*(d_num_species + 4) + d_num_species + 2] = Real(3)*w_y_BB + Real(10)*w_y_B -
                                        Real(18)*V_ghost[(j - 2)*(d_num_species + 4) + d_num_species + 2] +
                                        Real(6)*V_ghost[(j - 1)*(d_num_species + 4) + d_num_species + 2] +
                                        Real(12)*Real(dx[1])*dV_dy[d_num_species + 2];
                                    
                                    V_ghost[j*(d_num_species + 4) + d_num_species + 3] = Real(3)*p_y_BB + Real(10)*p_y_B -
                                        Real(18)*V_ghost[(j - 2)*(d_num_species + 4) + d_num_species + 3] +
                                        Real(6)*V_ghost[(j - 1)*(d_num_species + 4) + d_num_species + 3] +
                                        Real(12)*Real(dx[1])*dV_dy[d_num_species + 3];
                                }
                                else if (j == 3)
                                {
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        V_ghost[j*(d_num_species + 4) + si] = -Real(4)*rho_Y_y_BB[si] -
                                            Real(65)/Real(3)*rho_Y_y_B[si] +
                                            Real(40)*V_ghost[(j - 3)*(d_num_species + 4) + si] -
                                            Real(20)*V_ghost[(j - 2)*(d_num_species + 4) + si] +
                                            Real(20)/Real(3)*V_ghost[(j - 1)*(d_num_species + 4) + si] -
                                            Real(20)*Real(dx[1])*dV_dy[si];
                                    }
                                    
                                    V_ghost[j*(d_num_species + 4) + d_num_species] = -Real(4)*u_y_BB -
                                        Real(65)/Real(3)*u_y_B +
                                        Real(40)*V_ghost[(j - 3)*(d_num_species + 4) + d_num_species] -
                                        Real(20)*V_ghost[(j - 2)*(d_num_species + 4) + d_num_species] +
                                        Real(20)/Real(3)*V_ghost[(j - 1)*(d_num_species + 4) + d_num_species] -
                                        Real(20)*Real(dx[1])*dV_dy[d_num_species];
                                    
                                    V_ghost[j*(d_num_species + 4) + d_num_species + 1] = -Real(4)*v_y_BB -
                                        Real(65)/Real(3)*v_y_B +
                                        Real(40)*V_ghost[(j - 3)*(d_num_species + 4) + d_num_species + 1] -
                                        Real(20)*V_ghost[(j - 2)*(d_num_species + 4) + d_num_species + 1] +
                                        Real(20)/Real(3)*V_ghost[(j - 1)*(d_num_species + 4) + d_num_species + 1] -
                                        Real(20)*Real(dx[1])*dV_dy[d_num_species + 1];
                                    
                                    V_ghost[j*(d_num_species + 4) + d_num_species + 2] = -Real(4)*w_y_BB -
                                        Real(65)/Real(3)*w_y_B +
                                        Real(40)*V_ghost[(j - 3)*(d_num_species + 4) + d_num_species + 2] -
                                        Real(20)*V_ghost[(j - 2)*(d_num_species + 4) + d_num_species + 2] +
                                        Real(20)/Real(3)*V_ghost[(j - 1)*(d_num_species + 4) + d_num_species + 2] -
                                        Real(20)*Real(dx[1])*dV_dy[d_num_species + 2];
                                    
                                    V_ghost[j*(d_num_species + 4) + d_num_species + 3] = -Real(4)*p_y_BB -
                                        Real(65)/Real(3)*p_y_B +
                                        Real(40)*V_ghost[(j - 3)*(d_num_species + 4) + d_num_species + 3] -
                                        Real(20)*V_ghost[(j - 2)*(d_num_species + 4) + d_num_species + 3] +
                                        Real(20)/Real(3)*V_ghost[(j - 1)*(d_num_species + 4) + d_num_species + 3] -
                                        Real(20)*Real(dx[1])*dV_dy[d_num_species + 3];
                                }
                                else if (j == 4)
                                {
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        V_ghost[j*(d_num_species + 4) + si] = Real(5)*rho_Y_y_BB[si] +
                                            Real(77)/Real(2)*rho_Y_y_B[si] -
                                            Real(75)*V_ghost[(j - 4)*(d_num_species + 4) + si] +
                                            Real(50)*V_ghost[(j - 3)*(d_num_species + 4) + si] -
                                            Real(25)*V_ghost[(j - 2)*(d_num_species + 4) + si] +
                                            Real(15)/Real(2)*V_ghost[(j - 1)*(d_num_species + 4) + si] +
                                            Real(30)*Real(dx[1])*dV_dy[si];
                                    }
                                    
                                    V_ghost[j*(d_num_species + 4) + d_num_species] = Real(5)*u_y_BB +
                                        Real(77)/Real(2)*u_y_B -
                                        Real(75)*V_ghost[(j - 4)*(d_num_species + 4) + d_num_species] +
                                        Real(50)*V_ghost[(j - 3)*(d_num_species + 4) + d_num_species] -
                                        Real(25)*V_ghost[(j - 2)*(d_num_species + 4) + d_num_species] +
                                        Real(15)/Real(2)*V_ghost[(j - 1)*(d_num_species + 4) + d_num_species] +
                                        Real(30)*Real(dx[1])*dV_dy[d_num_species];
                                    
                                    V_ghost[j*(d_num_species + 4) + d_num_species + 1] = Real(5)*v_y_BB +
                                        Real(77)/Real(2)*v_y_B -
                                        Real(75)*V_ghost[(j - 4)*(d_num_species + 4) + d_num_species + 1] +
                                        Real(50)*V_ghost[(j - 3)*(d_num_species + 4) + d_num_species + 1] -
                                        Real(25)*V_ghost[(j - 2)*(d_num_species + 4) + d_num_species + 1] +
                                        Real(15)/Real(2)*V_ghost[(j - 1)*(d_num_species + 4) + d_num_species + 1] +
                                        Real(30)*Real(dx[1])*dV_dy[d_num_species + 1];
                                    
                                    V_ghost[j*(d_num_species + 4) + d_num_species + 2] = Real(5)*w_y_BB +
                                        Real(77)/Real(2)*w_y_B -
                                        Real(75)*V_ghost[(j - 4)*(d_num_species + 4) + d_num_species + 2] +
                                        Real(50)*V_ghost[(j - 3)*(d_num_species + 4) + d_num_species + 2] -
                                        Real(25)*V_ghost[(j - 2)*(d_num_species + 4) + d_num_species + 2] +
                                        Real(15)/Real(2)*V_ghost[(j - 1)*(d_num_species + 4) + d_num_species + 2] +
                                        Real(30)*Real(dx[1])*dV_dy[d_num_species + 2];
                                    
                                    V_ghost[j*(d_num_species + 4) + d_num_species + 3] = Real(5)*p_y_BB +
                                        Real(77)/Real(2)*p_y_B -
                                        Real(75)*V_ghost[(j - 4)*(d_num_species + 4) + d_num_species + 3] +
                                        Real(50)*V_ghost[(j - 3)*(d_num_species + 4) + d_num_species + 3] -
                                        Real(25)*V_ghost[(j - 2)*(d_num_species + 4) + d_num_species + 3] +
                                        Real(15)/Real(2)*V_ghost[(j - 1)*(d_num_species + 4) + d_num_species + 3] +
                                        Real(30)*Real(dx[1])*dV_dy[d_num_species + 3];
                                }
                                else if (j == 5)
                                {
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        V_ghost[j*(d_num_species + 4) + si] = -Real(6)*rho_Y_y_BB[si] -
                                            Real(609)/Real(10)*rho_Y_y_B[si] +
                                            Real(126)*V_ghost[(j - 5)*(d_num_species + 4) + si] -
                                            Real(105)*V_ghost[(j - 4)*(d_num_species + 4) + si] +
                                            Real(70)*V_ghost[(j - 3)*(d_num_species + 4) + si] -
                                            Real(63)/Real(2)*V_ghost[(j - 2)*(d_num_species + 4) + si] +
                                            Real(42)/Real(5)*V_ghost[(j - 1)*(d_num_species + 4) + si] -
                                            Real(42)*Real(dx[1])*dV_dy[si];
                                    }
                                    
                                    V_ghost[j*(d_num_species + 4) + d_num_species] = -Real(6)*u_y_BB -
                                        Real(609)/Real(10)*u_y_B +
                                        Real(126)*V_ghost[(j - 5)*(d_num_species + 4) + d_num_species] -
                                        Real(105)*V_ghost[(j - 4)*(d_num_species + 4) + d_num_species] +
                                        Real(70)*V_ghost[(j - 3)*(d_num_species + 4) + d_num_species] -
                                        Real(63)/Real(2)*V_ghost[(j - 2)*(d_num_species + 4) + d_num_species] +
                                        Real(42)/Real(5)*V_ghost[(j - 1)*(d_num_species + 4) + d_num_species] -
                                        Real(42)*Real(dx[1])*dV_dy[d_num_species];
                                    
                                    V_ghost[j*(d_num_species + 4) + d_num_species + 1] = -Real(6)*v_y_BB -
                                        Real(609)/Real(10)*v_y_B +
                                        Real(126)*V_ghost[(j - 5)*(d_num_species + 4) + d_num_species + 1] -
                                        Real(105)*V_ghost[(j - 4)*(d_num_species + 4) + d_num_species + 1] +
                                        Real(70)*V_ghost[(j - 3)*(d_num_species + 4) + d_num_species + 1] -
                                        Real(63)/Real(2)*V_ghost[(j - 2)*(d_num_species + 4) + d_num_species + 1] +
                                        Real(42)/Real(5)*V_ghost[(j - 1)*(d_num_species + 4) + d_num_species + 1] -
                                        Real(42)*Real(dx[1])*dV_dy[d_num_species + 1];
                                    
                                    V_ghost[j*(d_num_species + 4) + d_num_species + 2] = -Real(6)*w_y_BB -
                                        Real(609)/Real(10)*w_y_B +
                                        Real(126)*V_ghost[(j - 5)*(d_num_species + 4) + d_num_species + 2] -
                                        Real(105)*V_ghost[(j - 4)*(d_num_species + 4) + d_num_species + 2] +
                                        Real(70)*V_ghost[(j - 3)*(d_num_species + 4) + d_num_species + 2] -
                                        Real(63)/Real(2)*V_ghost[(j - 2)*(d_num_species + 4) + d_num_species + 2] +
                                        Real(42)/Real(5)*V_ghost[(j - 1)*(d_num_species + 4) + d_num_species + 2] -
                                        Real(42)*Real(dx[1])*dV_dy[d_num_species + 2];
                                    
                                    V_ghost[j*(d_num_species + 4) + d_num_species + 3] = -Real(6)*p_y_BB -
                                        Real(609)/Real(10)*p_y_B +
                                        Real(126)*V_ghost[(j - 5)*(d_num_species + 4) + d_num_species + 3] -
                                        Real(105)*V_ghost[(j - 4)*(d_num_species + 4) + d_num_species + 3] +
                                        Real(70)*V_ghost[(j - 3)*(d_num_species + 4) + d_num_species + 3] -
                                        Real(63)/Real(2)*V_ghost[(j - 2)*(d_num_species + 4) + d_num_species + 3] +
                                        Real(42)/Real(5)*V_ghost[(j - 1)*(d_num_species + 4) + d_num_species + 3] -
                                        Real(42)*Real(dx[1])*dV_dy[d_num_species + 3];
                                }
                                
                                /*
                                 * Compute the mixture density.
                                 */
                                
                                Real rho_ghost = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_ghost += V_ghost[j*(d_num_species + 4) + si];
                                }
                                
                                /*
                                 * Compute the mass fractions.
                                 */
                                
                                std::vector<Real> Y_ghost;
                                Y_ghost.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_ghost.push_back(V_ghost[j*(d_num_species + 4) + si]/rho_ghost);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions.
                                 */
                                
                                std::vector<const Real*> Y_ghost_ptr;
                                Y_ghost_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_ghost_ptr.push_back(&Y_ghost[si]);
                                }
                                
                                for(int si=0; si < d_num_species; si++)
                                {
                                    Q[si][idx_cell_rho_Y] = V_ghost[j*(d_num_species + 4) + si];
                                }
                                
                                Q[d_num_species][idx_cell_mom]     = rho_ghost*V_ghost[j*(d_num_species + 4) + d_num_species];
                                Q[d_num_species + 1][idx_cell_mom] = rho_ghost*V_ghost[j*(d_num_species + 4) + d_num_species + 1];
                                Q[d_num_species + 2][idx_cell_mom] = rho_ghost*V_ghost[j*(d_num_species + 4) + d_num_species + 2];
                                
                                const Real epsilon = d_equation_of_state_mixing_rules->
                                    getInternalEnergy(
                                        &rho_ghost,
                                        &V_ghost[j*(d_num_species + 4) + d_num_species + 3],
                                        Y_ghost_ptr);
                                
                                const Real E = rho_ghost*epsilon +
                                    half*(Q[d_num_species][idx_cell_mom]*Q[d_num_species][idx_cell_mom] +
                                        Q[d_num_species + 1][idx_cell_mom]*Q[d_num_species + 1][idx_cell_mom] +
                                        Q[d_num_species + 2][idx_cell_mom]*Q[d_num_species + 2][idx_cell_mom])/
                                        rho_ghost;
                                
                                Q[d_num_species + 3][idx_cell_E] = E;
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
                            << ": FlowModelBoundaryUtilitiesFourEqnConservative::fill3dFaceBoundaryData()\n"
                            << "Non-reflecting outflow BC doesn't support more than six ghost cells yet!");
                    }
                    
                    for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                    {
                        for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                        {
                            // Get the grid spacing.
                            const double* const dx = patch_geom->getDx();
                            
                            const int idx_cell_rho_Y_z_F = (i + num_subghosts_conservative_var[0][0]) +
                                (j + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                (interior_box_lo_idx[2] + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][1]*
                                    subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_rho_Y_z_FF = (i + num_subghosts_conservative_var[0][0]) +
                                (j + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                (interior_box_lo_idx[2] + 1 + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                    subghostcell_dims_conservative_var[0][1];
                            
                            const int idx_cell_rho_Y_z_FFF = (i + num_subghosts_conservative_var[0][0]) +
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
                            
                            std::vector<Real> rho_Y_z_F;
                            std::vector<Real> rho_Y_z_FF;
                            std::vector<Real> rho_Y_z_FFF;
                            rho_Y_z_F.reserve(d_num_species);
                            rho_Y_z_FF.reserve(d_num_species);
                            rho_Y_z_FFF.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                rho_Y_z_F.push_back(Q[si][idx_cell_rho_Y_z_F]);
                                rho_Y_z_FF.push_back(Q[si][idx_cell_rho_Y_z_FF]);
                                rho_Y_z_FFF.push_back(Q[si][idx_cell_rho_Y_z_FFF]);
                            }
                            
                            /*
                             * Compute the mixture density.
                             */
                            
                            Real rho_z_F   = Real(0);
                            Real rho_z_FF  = Real(0);
                            Real rho_z_FFF = Real(0);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                rho_z_F   += Q[si][idx_cell_rho_Y_z_F];
                                rho_z_FF  += Q[si][idx_cell_rho_Y_z_FF];
                                rho_z_FFF += Q[si][idx_cell_rho_Y_z_FFF];
                            }
                            
                            /*
                             * Compute the mass fractions.
                             */
                            
                            std::vector<Real> Y_z_F;
                            std::vector<Real> Y_z_FF;
                            std::vector<Real> Y_z_FFF;
                            Y_z_F.reserve(d_num_species);
                            Y_z_FF.reserve(d_num_species);
                            Y_z_FFF.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_z_F.push_back(Q[si][idx_cell_rho_Y_z_F]/rho_z_F);
                                Y_z_FF.push_back(Q[si][idx_cell_rho_Y_z_FF]/rho_z_FF);
                                Y_z_FFF.push_back(Q[si][idx_cell_rho_Y_z_FFF]/rho_z_FFF);
                            }
                            
                            /*
                             * Get the pointers to the mass fractions.
                             */
                            
                            std::vector<const Real*> Y_z_F_ptr;
                            std::vector<const Real*> Y_z_FF_ptr;
                            std::vector<const Real*> Y_z_FFF_ptr;
                            Y_z_F_ptr.reserve(d_num_species);
                            Y_z_FF_ptr.reserve(d_num_species);
                            Y_z_FFF_ptr.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_z_F_ptr.push_back(&Y_z_F[si]);
                                Y_z_FF_ptr.push_back(&Y_z_FF[si]);
                                Y_z_FFF_ptr.push_back(&Y_z_FFF[si]);
                            }
                            
                            const Real u_z_F   = Q[d_num_species][idx_cell_mom_z_F]/rho_z_F;
                            const Real u_z_FF  = Q[d_num_species][idx_cell_mom_z_FF]/rho_z_FF;
                            const Real u_z_FFF = Q[d_num_species][idx_cell_mom_z_FFF]/rho_z_FFF;
                            
                            const Real v_z_F   = Q[d_num_species + 1][idx_cell_mom_z_F]/rho_z_F;
                            const Real v_z_FF  = Q[d_num_species + 1][idx_cell_mom_z_FF]/rho_z_FF;
                            const Real v_z_FFF = Q[d_num_species + 1][idx_cell_mom_z_FFF]/rho_z_FFF;
                            
                            const Real w_z_F   = Q[d_num_species + 2][idx_cell_mom_z_F]/rho_z_F;
                            const Real w_z_FF  = Q[d_num_species + 2][idx_cell_mom_z_FF]/rho_z_FF;
                            const Real w_z_FFF = Q[d_num_species + 2][idx_cell_mom_z_FFF]/rho_z_FFF;
                            
                            const Real half = Real(1)/Real(2);
                            const Real epsilon_z_F   = Q[d_num_species + 3][idx_cell_E_z_F]/rho_z_F - half*(u_z_F*u_z_F + v_z_F*v_z_F + w_z_F*w_z_F);
                            const Real epsilon_z_FF  = Q[d_num_species + 3][idx_cell_E_z_FF]/rho_z_FF - half*(u_z_FF*u_z_FF + v_z_FF*v_z_FF + w_z_FF*w_z_FF);
                            const Real epsilon_z_FFF = Q[d_num_species + 3][idx_cell_E_z_FFF]/rho_z_FFF - half*(u_z_FFF*u_z_FFF + v_z_FFF*v_z_FFF + w_z_FFF*w_z_FFF);
                            
                            Real p_z_F = d_equation_of_state_mixing_rules->
                                getPressure(
                                    &rho_z_F,
                                    &epsilon_z_F,
                                    Y_z_F_ptr);
                            
                            Real p_z_FF = d_equation_of_state_mixing_rules->
                                getPressure(
                                    &rho_z_FF,
                                    &epsilon_z_FF,
                                    Y_z_FF_ptr);
                            
                            Real p_z_FFF = d_equation_of_state_mixing_rules->
                                getPressure(
                                    &rho_z_FFF,
                                    &epsilon_z_FFF,
                                    Y_z_FFF_ptr);
                            
                            /*
                             * Compute derivatives in z-direction.
                             */
                            
                            std::vector<Real> drho_Y_dz;
                            drho_Y_dz.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                drho_Y_dz.push_back(-(Q[si][idx_cell_rho_Y_z_FFF] - Real(4)*Q[si][idx_cell_rho_Y_z_FF] +
                                    Real(3)*Q[si][idx_cell_rho_Y_z_F])/(Real(2)*Real(dx[2])));
                            }
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
                                
                                const int idx_cell_rho_Y_x_R = (i + 1 + num_subghosts_conservative_var[0][0]) +
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
                                
                                /*
                                 * Compute the mixture density.
                                 */
                                
                                Real rho_x_R = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_x_R += Q[si][idx_cell_rho_Y_x_R];
                                }
                                
                                /*
                                 * Compute the mass fractions.
                                 */
                                
                                std::vector<Real> Y_x_R;
                                Y_x_R.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_x_R.push_back(Q[si][idx_cell_rho_Y_x_R]/rho_x_R);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions.
                                 */
                            
                                std::vector<const Real*> Y_x_R_ptr;
                                Y_x_R_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_x_R_ptr.push_back(&Y_x_R[si]);
                                }
                                
                                const Real u_x_R = Q[d_num_species][idx_cell_mom_x_R]/rho_x_R;
                                const Real v_x_R = Q[d_num_species + 1][idx_cell_mom_x_R]/rho_x_R;
                                const Real w_x_R = Q[d_num_species + 2][idx_cell_mom_x_R]/rho_x_R;
                                const Real epsilon_x_R = Q[d_num_species + 3][idx_cell_E_x_R]/rho_x_R - half*(u_x_R*u_x_R + v_x_R*v_x_R + w_x_R*w_x_R);
                                
                                Real p_x_R = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_x_R,
                                        &epsilon_x_R,
                                        Y_x_R_ptr);
                                
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
                                
                                const int idx_cell_rho_Y_x_L = (i - 1 + num_subghosts_conservative_var[0][0]) +
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
                                
                                /*
                                 * Compute the mixture density.
                                 */
                                
                                Real rho_x_L = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_x_L += Q[si][idx_cell_rho_Y_x_L];
                                }
                                
                                /*
                                 * Compute the mass fractions.
                                 */
                                
                                std::vector<Real> Y_x_L;
                                Y_x_L.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_x_L.push_back(Q[si][idx_cell_rho_Y_x_L]/rho_x_L);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions.
                                 */
                            
                                std::vector<const Real*> Y_x_L_ptr;
                                Y_x_L_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_x_L_ptr.push_back(&Y_x_L[si]);
                                }
                                
                                const Real u_x_L = Q[d_num_species][idx_cell_mom_x_L]/rho_x_L;
                                const Real v_x_L = Q[d_num_species + 1][idx_cell_mom_x_L]/rho_x_L;
                                const Real w_x_L = Q[d_num_species + 2][idx_cell_mom_x_L]/rho_x_L;
                                const Real epsilon_x_L = Q[d_num_species + 3][idx_cell_E_x_L]/rho_x_L - half*(u_x_L*u_x_L + v_x_L*v_x_L + w_x_L*w_x_L);
                                
                                Real p_x_L = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_x_L,
                                        &epsilon_x_L,
                                        Y_x_L_ptr);
                                
                                // One-sided derivatives.
                                du_dx = (u_z_F - u_x_L)/Real(dx[0]);
                                // dv_dx = (v_z_F - v_x_L)/Real(dx[0]);
                                dw_dx = (w_z_F - w_x_L)/Real(dx[0]);
                                dp_dx = (p_z_F - p_x_L)/Real(dx[0]);
                            }
                            else
                            {
                                const int idx_cell_rho_Y_x_L = (i - 1 + num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (interior_box_lo_idx[2] + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_rho_Y_x_R = (i + 1 + num_subghosts_conservative_var[0][0]) +
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
                                
                                /*
                                 * Compute the mixture density.
                                 */
                                
                                Real rho_x_L = Real(0);
                                Real rho_x_R = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_x_L += Q[si][idx_cell_rho_Y_x_L];
                                    rho_x_R += Q[si][idx_cell_rho_Y_x_R];
                                }
                                
                                /*
                                 * Compute the mass fractions.
                                 */
                                
                                std::vector<Real> Y_x_L;
                                std::vector<Real> Y_x_R;
                                Y_x_L.reserve(d_num_species);
                                Y_x_R.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_x_L.push_back(Q[si][idx_cell_rho_Y_x_L]/rho_x_L);
                                    Y_x_R.push_back(Q[si][idx_cell_rho_Y_x_R]/rho_x_R);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions.
                                 */
                                
                                std::vector<const Real*> Y_x_L_ptr;
                                std::vector<const Real*> Y_x_R_ptr;
                                                                    
                                Y_x_L_ptr.reserve(d_num_species);
                                Y_x_R_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_x_L_ptr.push_back(&Y_x_L[si]);
                                    Y_x_R_ptr.push_back(&Y_x_R[si]);
                                }
                                
                                const Real u_x_L = Q[d_num_species][idx_cell_mom_x_L]/rho_x_L;
                                const Real u_x_R = Q[d_num_species][idx_cell_mom_x_R]/rho_x_R;
                                
                                const Real v_x_L = Q[d_num_species + 1][idx_cell_mom_x_L]/rho_x_L;
                                const Real v_x_R = Q[d_num_species + 1][idx_cell_mom_x_R]/rho_x_R;
                                
                                const Real w_x_L = Q[d_num_species + 2][idx_cell_mom_x_L]/rho_x_L;
                                const Real w_x_R = Q[d_num_species + 2][idx_cell_mom_x_R]/rho_x_R;
                                
                                const Real epsilon_x_L = Q[d_num_species + 3][idx_cell_E_x_L]/rho_x_L - half*(u_x_L*u_x_L + v_x_L*v_x_L + w_x_L*w_x_L);
                                const Real epsilon_x_R = Q[d_num_species + 3][idx_cell_E_x_R]/rho_x_R - half*(u_x_R*u_x_R + v_x_R*v_x_R + w_x_R*w_x_R);
                                
                                Real p_x_L = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_x_L,
                                        &epsilon_x_L,
                                        Y_x_L_ptr);
                            
                                Real p_x_R = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_x_R,
                                        &epsilon_x_R,
                                        Y_x_R_ptr);
                                
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
                                
                                const int idx_cell_rho_Y_y_T = (i + num_subghosts_conservative_var[0][0]) +
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
                                
                                /*
                                 * Compute the mixture density.
                                 */
                                
                                Real rho_y_T = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_y_T += Q[si][idx_cell_rho_Y_y_T];
                                }
                                
                                /*
                                 * Compute the mass fractions.
                                 */
                                
                                std::vector<Real> Y_y_T;
                                Y_y_T.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_y_T.push_back(Q[si][idx_cell_rho_Y_y_T]/rho_y_T);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions.
                                 */
                                
                                std::vector<const Real*> Y_y_T_ptr;
                                Y_y_T_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_y_T_ptr.push_back(&Y_y_T[si]);
                                }
                                
                                const Real u_y_T = Q[d_num_species][idx_cell_mom_y_T]/rho_y_T;
                                const Real v_y_T = Q[d_num_species + 1][idx_cell_mom_y_T]/rho_y_T;
                                const Real w_y_T = Q[d_num_species + 2][idx_cell_mom_y_T]/rho_y_T;
                                const Real epsilon_y_T = Q[d_num_species + 3][idx_cell_E_y_T]/rho_y_T - half*(u_y_T*u_y_T + v_y_T*v_y_T + w_y_T*w_y_T);
                                
                                Real p_y_T = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_y_T,
                                        &epsilon_y_T,
                                        Y_y_T_ptr);
                                
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
                                
                                const int idx_cell_rho_Y_y_B = (i + num_subghosts_conservative_var[0][0]) +
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
                                
                                /*
                                 * Compute the mixture density.
                                 */
                                
                                Real rho_y_B = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_y_B += Q[si][idx_cell_rho_Y_y_B];
                                }
                                
                                /*
                                 * Compute the mass fractions.
                                 */
                                
                                std::vector<Real> Y_y_B;
                                Y_y_B.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_y_B.push_back(Q[si][idx_cell_rho_Y_y_B]/rho_y_B);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions.
                                 */
                                
                                std::vector<const Real*> Y_y_B_ptr;
                                Y_y_B_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_y_B_ptr.push_back(&Y_y_B[si]);
                                }
                                
                                const Real u_y_B = Q[d_num_species][idx_cell_mom_y_B]/rho_y_B;
                                const Real v_y_B = Q[d_num_species + 1][idx_cell_mom_y_B]/rho_y_B;
                                const Real w_y_B = Q[d_num_species + 2][idx_cell_mom_y_B]/rho_y_B;
                                const Real epsilon_y_B = Q[d_num_species + 3][idx_cell_E_y_B]/rho_y_B - half*(u_y_B*u_y_B + v_y_B*v_y_B + w_y_B*w_y_B);
                                
                                Real p_y_B = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_y_B,
                                        &epsilon_y_B,
                                        Y_y_B_ptr);
                                
                                // One-sided derivatives.
                                //du_dy = (u_z_F - u_y_B)/Real(dx[1]);
                                dv_dy = (v_z_F - v_y_B)/Real(dx[1]);
                                dw_dy = (w_z_F - w_y_B)/Real(dx[1]);
                                dp_dy = (p_z_F - p_y_B)/Real(dx[1]);
                            }
                            else
                            {
                                const int idx_cell_rho_Y_y_B = (i + num_subghosts_conservative_var[0][0]) +
                                    (j - 1 + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (interior_box_lo_idx[2] +  num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_rho_Y_y_T = (i + num_subghosts_conservative_var[0][0]) +
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
                                
                                /*
                                 * Compute the mixture density.
                                 */
                                
                                Real rho_y_B = Real(0);
                                Real rho_y_T = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_y_B += Q[si][idx_cell_rho_Y_y_B];
                                    rho_y_T += Q[si][idx_cell_rho_Y_y_T];
                                }
                                
                                /*
                                 * Compute the mass fractions.
                                 */
                                
                                std::vector<Real> Y_y_B;
                                std::vector<Real> Y_y_T;
                                Y_y_B.reserve(d_num_species);
                                Y_y_T.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_y_B.push_back(Q[si][idx_cell_rho_Y_y_B]/rho_y_B);
                                    Y_y_T.push_back(Q[si][idx_cell_rho_Y_y_T]/rho_y_T);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions.
                                 */
                                
                                std::vector<const Real*> Y_y_B_ptr;
                                std::vector<const Real*> Y_y_T_ptr;
                                Y_y_B_ptr.reserve(d_num_species);
                                Y_y_T_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_y_B_ptr.push_back(&Y_y_B[si]);
                                    Y_y_T_ptr.push_back(&Y_y_T[si]);
                                }
                                
                                const Real u_y_B = Q[d_num_species][idx_cell_mom_y_B]/rho_y_B;
                                const Real u_y_T = Q[d_num_species][idx_cell_mom_y_T]/rho_y_T;
                                
                                const Real v_y_B = Q[d_num_species + 1][idx_cell_mom_y_B]/rho_y_B;
                                const Real v_y_T = Q[d_num_species + 1][idx_cell_mom_y_T]/rho_y_T;
                                
                                const Real w_y_B = Q[d_num_species + 2][idx_cell_mom_y_B]/rho_y_B;
                                const Real w_y_T = Q[d_num_species + 2][idx_cell_mom_y_T]/rho_y_T;
                                
                                const Real epsilon_y_B = Q[d_num_species + 3][idx_cell_E_y_B]/rho_y_B - half*(u_y_B*u_y_B + v_y_B*v_y_B + w_y_B*w_y_B);
                                const Real epsilon_y_T = Q[d_num_species + 3][idx_cell_E_y_T]/rho_y_T - half*(u_y_T*u_y_T + v_y_T*v_y_T + w_y_T*w_y_T);
                                
                                Real p_y_B = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_y_B,
                                        &epsilon_y_B,
                                        Y_y_B_ptr);
                                
                                Real p_y_T = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_y_T,
                                        &epsilon_y_T,
                                        Y_y_T_ptr);
                                
                                // Central derivatives.
                                // du_dy = (u_y_T - u_y_B)/(Real(2)*Real(dx[1]));
                                dv_dy = (v_y_T - v_y_B)/(Real(2)*Real(dx[1]));
                                dw_dy = (w_y_T - w_y_B)/(Real(2)*Real(dx[1]));
                                dp_dy = (p_y_T - p_y_B)/(Real(2)*Real(dx[1]));
                            }
                            
                            // Compute sound speed.
                            
                            const Real Gamma_z_F = d_equation_of_state_mixing_rules->getGruneisenParameter(
                                &rho_z_F,
                                &p_z_F,
                                Y_z_F_ptr);

                            const std::vector<Real> Psi_z_F = d_equation_of_state_mixing_rules->
                                getPressureDerivativeWithPartialDensities(
                                        &rho_z_F,
                                        &p_z_F,
                                        Y_z_F_ptr);
                            
                            Real c_z_F = Gamma_z_F*p_z_F/rho_z_F;
                            for (int si = 0; si < d_num_species; si++)
                            {
                                c_z_F += Y_z_F[si]*Psi_z_F[si];
                            }
                            c_z_F = std::sqrt(c_z_F);
                            
                            const Real lambda_last = w_z_F + c_z_F;
                            
                            // Compute vector Lambda^(-1) * L.
                            
                            Real Lambda_inv_L[d_num_species + 4];
                            
                            const Real& p_t         = d_bdry_face_nonreflecting_outflow_p_t[face_loc];
                            const Real& sigma       = d_bdry_face_nonreflecting_outflow_sigma[face_loc];
                            const Real& beta        = d_bdry_face_nonreflecting_outflow_beta[face_loc];
                            const Real& length_char = d_bdry_face_nonreflecting_outflow_length_char[face_loc];
                            
                            const Real T_last = u_z_F*(dp_dx + rho_z_F*c_z_F*dw_dx) + rho_z_F*c_z_F*c_z_F*du_dx +
                                v_z_F*(dp_dy + rho_z_F*c_z_F*dw_dy) + rho_z_F*c_z_F*c_z_F*dv_dy;
                            
                            const Real M_sq = (u_z_F*u_z_F + v_z_F*v_z_F + w_z_F*w_z_F)/(c_z_F*c_z_F);
                            const Real K = sigma*c_z_F*(Real(1) - M_sq)/length_char;
                            
                            Lambda_inv_L[0] = dp_dz - rho_z_F*c_z_F*dw_dz;
                            Lambda_inv_L[1] = du_dz;
                            Lambda_inv_L[2] = dv_dz;
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Lambda_inv_L[si + 3] = c_z_F*c_z_F*drho_Y_dz[si] - Y_z_F[si]*dp_dz;
                            }
                            Lambda_inv_L[d_num_species + 3] = (Real(1)/lambda_last)*(K*(p_z_F - p_t) - (Real(1) - beta)*T_last);
                            
                            // Compute dV_dz.
                            
                            const Real c_sq_inv  = Real(1)/(c_z_F*c_z_F);
                            const Real rho_c_inv = Real(1)/(rho_z_F*c_z_F);
                            
                            Real dV_dz[d_num_species + 4];
                            
                            for (int si = 0; si < d_num_species; si++)
                            {
                                dV_dz[si] = half*c_sq_inv*Y_z_F[si]*(Lambda_inv_L[0] + Lambda_inv_L[d_num_species + 3]) +
                                    c_sq_inv*Lambda_inv_L[si + 3];
                            }
                            dV_dz[d_num_species]     = Lambda_inv_L[1];
                            dV_dz[d_num_species + 1] = Lambda_inv_L[2];
                            dV_dz[d_num_species + 2] = half*rho_c_inv*(-Lambda_inv_L[0] + Lambda_inv_L[d_num_species + 3]);
                            dV_dz[d_num_species + 3] = half*(Lambda_inv_L[0] + Lambda_inv_L[d_num_species + 3]);
                            
                            Real V_ghost[(d_num_species + 4)*num_ghosts_to_fill];
                            
                            for (int k = num_ghosts_to_fill - 1; k >= 0; k--)
                            {
                                const int idx_cell_rho_Y = (i + num_subghosts_conservative_var[0][0]) +
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
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        V_ghost[k*(d_num_species + 4) + si] = rho_Y_z_FF[si] - Real(2)*Real(dx[2])*dV_dz[si];
                                    }
                                    
                                    V_ghost[k*(d_num_species + 4) + d_num_species]     = u_z_FF - Real(2)*Real(dx[2])*dV_dz[d_num_species];
                                    V_ghost[k*(d_num_species + 4) + d_num_species + 1] = v_z_FF - Real(2)*Real(dx[2])*dV_dz[d_num_species + 1];
                                    V_ghost[k*(d_num_species + 4) + d_num_species + 2] = w_z_FF - Real(2)*Real(dx[2])*dV_dz[d_num_species + 2];
                                    V_ghost[k*(d_num_species + 4) + d_num_species + 3] = p_z_FF - Real(2)*Real(dx[2])*dV_dz[d_num_species + 3];
                                }
                                else if (k == num_ghosts_to_fill - 2)
                                {
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        V_ghost[k*(d_num_species + 4) + si] = -Real(2)*rho_Y_z_FF[si] - Real(3)*rho_Y_z_F[si] +
                                            Real(6)*V_ghost[(k + 1)*(d_num_species + 4) + si] + Real(6)*Real(dx[2])*dV_dz[si];
                                    }
                                    
                                    V_ghost[k*(d_num_species + 4) + d_num_species] = -Real(2)*u_z_FF - Real(3)*u_z_F +
                                        Real(6)*V_ghost[(k + 1)*(d_num_species + 4) + d_num_species] +
                                        Real(6)*Real(dx[2])*dV_dz[d_num_species];
                                    
                                    V_ghost[k*(d_num_species + 4) + d_num_species + 1] = -Real(2)*v_z_FF - Real(3)*v_z_F +
                                        Real(6)*V_ghost[(k + 1)*(d_num_species + 4) + d_num_species + 1] +
                                        Real(6)*Real(dx[2])*dV_dz[d_num_species + 1];
                                    
                                    V_ghost[k*(d_num_species + 4) + d_num_species + 2] = -Real(2)*w_z_FF - Real(3)*w_z_F +
                                        Real(6)*V_ghost[(k + 1)*(d_num_species + 4) + d_num_species + 2] +
                                        Real(6)*Real(dx[2])*dV_dz[d_num_species + 2];
                                    
                                    V_ghost[k*(d_num_species + 4) + d_num_species + 3] = -Real(2)*p_z_FF - Real(3)*p_z_F +
                                        Real(6)*V_ghost[(k + 1)*(d_num_species + 4) + d_num_species + 3] +
                                        Real(6)*Real(dx[2])*dV_dz[d_num_species + 3];
                                }
                                else if (k == num_ghosts_to_fill - 3)
                                {
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        V_ghost[k*(d_num_species + 4) + si] = Real(3)*rho_Y_z_FF[si] + Real(10)*rho_Y_z_F[si] -
                                            Real(18)*V_ghost[(k + 2)*(d_num_species + 4) + si] +
                                            Real(6)*V_ghost[(k + 1)*(d_num_species + 4) + si] -
                                            Real(12)*Real(dx[2])*dV_dz[si];
                                    }
                                    
                                    V_ghost[k*(d_num_species + 4) + d_num_species] = Real(3)*u_z_FF + Real(10)*u_z_F -
                                        Real(18)*V_ghost[(k + 2)*(d_num_species + 4) + d_num_species] +
                                        Real(6)*V_ghost[(k + 1)*(d_num_species + 4) + d_num_species] -
                                        Real(12)*Real(dx[2])*dV_dz[d_num_species];
                                    
                                    V_ghost[k*(d_num_species + 4) + d_num_species + 1] = Real(3)*v_z_FF + Real(10)*v_z_F -
                                        Real(18)*V_ghost[(k + 2)*(d_num_species + 4) + d_num_species + 1] +
                                        Real(6)*V_ghost[(k + 1)*(d_num_species + 4) + d_num_species + 1] -
                                        Real(12)*Real(dx[2])*dV_dz[d_num_species + 1];
                                    
                                    V_ghost[k*(d_num_species + 4) + d_num_species + 2] = Real(3)*w_z_FF + Real(10)*w_z_F -
                                        Real(18)*V_ghost[(k + 2)*(d_num_species + 4) + d_num_species + 2] +
                                        Real(6)*V_ghost[(k + 1)*(d_num_species + 4) + d_num_species + 2] -
                                        Real(12)*Real(dx[2])*dV_dz[d_num_species + 2];
                                    
                                    V_ghost[k*(d_num_species + 4) + d_num_species + 3] = Real(3)*p_z_FF + Real(10)*p_z_F -
                                        Real(18)*V_ghost[(k + 2)*(d_num_species + 4) + d_num_species + 3] +
                                        Real(6)*V_ghost[(k + 1)*(d_num_species + 4) + d_num_species + 3] -
                                        Real(12)*Real(dx[2])*dV_dz[d_num_species + 3];
                                }
                                else if (k == num_ghosts_to_fill - 4)
                                {
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        V_ghost[k*(d_num_species + 4) + si] = -Real(4)*rho_Y_z_FF[si] -
                                            Real(65)/Real(3)*rho_Y_z_F[si] +
                                            Real(40)*V_ghost[(k + 3)*(d_num_species + 4) + si] -
                                            Real(20)*V_ghost[(k + 2)*(d_num_species + 4) + si] +
                                            Real(20)/Real(3)*V_ghost[(k + 1)*(d_num_species + 4) + si] +
                                            Real(20)*Real(dx[2])*dV_dz[si];
                                    }
                                    
                                    V_ghost[k*(d_num_species + 4) + d_num_species] = -Real(4)*u_z_FF -
                                        Real(65)/Real(3)*u_z_F +
                                        Real(40)*V_ghost[(k + 3)*(d_num_species + 4) + d_num_species] -
                                        Real(20)*V_ghost[(k + 2)*(d_num_species + 4) + d_num_species] +
                                        Real(20)/Real(3)*V_ghost[(k + 1)*(d_num_species + 4) + d_num_species] +
                                        Real(20)*Real(dx[2])*dV_dz[d_num_species];
                                    
                                    V_ghost[k*(d_num_species + 4) + d_num_species + 1] = -Real(4)*v_z_FF -
                                        Real(65)/Real(3)*v_z_F +
                                        Real(40)*V_ghost[(k + 3)*(d_num_species + 4) + d_num_species + 1] -
                                        Real(20)*V_ghost[(k + 2)*(d_num_species + 4) + d_num_species + 1] +
                                        Real(20)/Real(3)*V_ghost[(k + 1)*(d_num_species + 4) + d_num_species + 1] +
                                        Real(20)*Real(dx[2])*dV_dz[d_num_species + 1];
                                    
                                    V_ghost[k*(d_num_species + 4) + d_num_species + 2] = -Real(4)*w_z_FF -
                                        Real(65)/Real(3)*w_z_F +
                                        Real(40)*V_ghost[(k + 3)*(d_num_species + 4) + d_num_species + 2] -
                                        Real(20)*V_ghost[(k + 2)*(d_num_species + 4) + d_num_species + 2] +
                                        Real(20)/Real(3)*V_ghost[(k + 1)*(d_num_species + 4) + d_num_species + 2] +
                                        Real(20)*Real(dx[2])*dV_dz[d_num_species + 2];
                                    
                                    V_ghost[k*(d_num_species + 4) + d_num_species + 3] = -Real(4)*p_z_FF -
                                        Real(65)/Real(3)*p_z_F +
                                        Real(40)*V_ghost[(k + 3)*(d_num_species + 4) + d_num_species + 3] -
                                        Real(20)*V_ghost[(k + 2)*(d_num_species + 4) + d_num_species + 3] +
                                        Real(20)/Real(3)*V_ghost[(k + 1)*(d_num_species + 4) + d_num_species + 3] +
                                        Real(20)*Real(dx[2])*dV_dz[d_num_species + 3];
                                }
                                else if (k == num_ghosts_to_fill - 5)
                                {
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        V_ghost[k*(d_num_species + 4) + si] = Real(5)*rho_Y_z_FF[si] +
                                            Real(77)/Real(2)*rho_Y_z_F[si] -
                                            Real(75)*V_ghost[(k + 4)*(d_num_species + 4) + si] +
                                            Real(50)*V_ghost[(k + 3)*(d_num_species + 4) + si] -
                                            Real(25)*V_ghost[(k + 2)*(d_num_species + 4) + si] +
                                            Real(15)/Real(2)*V_ghost[(k + 1)*(d_num_species + 4) + si] -
                                            Real(30)*Real(dx[2])*dV_dz[si];
                                    }
                                    
                                    V_ghost[k*(d_num_species + 4) + d_num_species] = Real(5)*u_z_FF +
                                        Real(77)/Real(2)*u_z_F -
                                        Real(75)*V_ghost[(k + 4)*(d_num_species + 4) + d_num_species] +
                                        Real(50)*V_ghost[(k + 3)*(d_num_species + 4) + d_num_species] -
                                        Real(25)*V_ghost[(k + 2)*(d_num_species + 4) + d_num_species] +
                                        Real(15)/Real(2)*V_ghost[(k + 1)*(d_num_species + 4) + d_num_species] -
                                        Real(30)*Real(dx[2])*dV_dz[d_num_species];
                                    
                                    V_ghost[k*(d_num_species + 4) + d_num_species + 1] = Real(5)*v_z_FF +
                                        Real(77)/Real(2)*v_z_F -
                                        Real(75)*V_ghost[(k + 4)*(d_num_species + 4) + d_num_species + 1] +
                                        Real(50)*V_ghost[(k + 3)*(d_num_species + 4) + d_num_species + 1] -
                                        Real(25)*V_ghost[(k + 2)*(d_num_species + 4) + d_num_species + 1] +
                                        Real(15)/Real(2)*V_ghost[(k + 1)*(d_num_species + 4) + d_num_species + 1] -
                                        Real(30)*Real(dx[2])*dV_dz[d_num_species + 1];
                                    
                                    V_ghost[k*(d_num_species + 4) + d_num_species + 2] = Real(5)*w_z_FF +
                                        Real(77)/Real(2)*w_z_F -
                                        Real(75)*V_ghost[(k + 4)*(d_num_species + 4) + d_num_species + 2] +
                                        Real(50)*V_ghost[(k + 3)*(d_num_species + 4) + d_num_species + 2] -
                                        Real(25)*V_ghost[(k + 2)*(d_num_species + 4) + d_num_species + 2] +
                                        Real(15)/Real(2)*V_ghost[(k + 1)*(d_num_species + 4) + d_num_species + 2] -
                                        Real(30)*Real(dx[2])*dV_dz[d_num_species + 2];
                                    
                                    V_ghost[k*(d_num_species + 4) + d_num_species + 3] = Real(5)*p_z_FF +
                                        Real(77)/Real(2)*p_z_F -
                                        Real(75)*V_ghost[(k + 4)*(d_num_species + 4) + d_num_species + 3] +
                                        Real(50)*V_ghost[(k + 3)*(d_num_species + 4) + d_num_species + 3] -
                                        Real(25)*V_ghost[(k + 2)*(d_num_species + 4) + d_num_species + 3] +
                                        Real(15)/Real(2)*V_ghost[(k + 1)*(d_num_species + 4) + d_num_species + 3] -
                                        Real(30)*Real(dx[2])*dV_dz[d_num_species + 3];
                                }
                                else if (k == num_ghosts_to_fill - 6)
                                {
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        V_ghost[k*(d_num_species + 4) + si] = -Real(6)*rho_Y_z_FF[si] -
                                            Real(609)/Real(10)*rho_Y_z_F[si] +
                                            Real(126)*V_ghost[(k + 5)*(d_num_species + 4) + si] -
                                            Real(105)*V_ghost[(k + 4)*(d_num_species + 4) + si] +
                                            Real(70)*V_ghost[(k + 3)*(d_num_species + 4) + si] -
                                            Real(63)/Real(2)*V_ghost[(k + 2)*(d_num_species + 4) + si] +
                                            Real(42)/Real(5)*V_ghost[(k + 1)*(d_num_species + 4) + si] +
                                            Real(42)*Real(dx[2])*dV_dz[si];
                                    }
                                    
                                    V_ghost[k*(d_num_species + 4) + d_num_species] = -Real(6)*u_z_FF -
                                        Real(609)/Real(10)*u_z_F +
                                        Real(126)*V_ghost[(k + 5)*(d_num_species + 4) + d_num_species] -
                                        Real(105)*V_ghost[(k + 4)*(d_num_species + 4) + d_num_species] +
                                        Real(70)*V_ghost[(k + 3)*(d_num_species + 4) + d_num_species] -
                                        Real(63)/Real(2)*V_ghost[(k + 2)*(d_num_species + 4) + d_num_species] +
                                        Real(42)/Real(5)*V_ghost[(k + 1)*(d_num_species + 4) + d_num_species] +
                                        Real(42)*Real(dx[2])*dV_dz[d_num_species];
                                    
                                    V_ghost[k*(d_num_species + 4) + d_num_species + 1] = -Real(6)*v_z_FF -
                                        Real(609)/Real(10)*v_z_F +
                                        Real(126)*V_ghost[(k + 5)*(d_num_species + 4) + d_num_species + 1] -
                                        Real(105)*V_ghost[(k + 4)*(d_num_species + 4) + d_num_species + 1] +
                                        Real(70)*V_ghost[(k + 3)*(d_num_species + 4) + d_num_species + 1] -
                                        Real(63)/Real(2)*V_ghost[(k + 2)*(d_num_species + 4) + d_num_species + 1] +
                                        Real(42)/Real(5)*V_ghost[(k + 1)*(d_num_species + 4) + d_num_species + 1] +
                                        Real(42)*Real(dx[2])*dV_dz[d_num_species + 1];
                                    
                                    V_ghost[k*(d_num_species + 4) + d_num_species + 2] = -Real(6)*w_z_FF -
                                        Real(609)/Real(10)*w_z_F +
                                        Real(126)*V_ghost[(k + 5)*(d_num_species + 4) + d_num_species + 2] -
                                        Real(105)*V_ghost[(k + 4)*(d_num_species + 4) + d_num_species + 2] +
                                        Real(70)*V_ghost[(k + 3)*(d_num_species + 4) + d_num_species + 2] -
                                        Real(63)/Real(2)*V_ghost[(k + 2)*(d_num_species + 4) + d_num_species + 2] +
                                        Real(42)/Real(5)*V_ghost[(k + 1)*(d_num_species + 4) + d_num_species + 2] +
                                        Real(42)*Real(dx[2])*dV_dz[d_num_species + 2];
                                    
                                    V_ghost[k*(d_num_species + 4) + d_num_species + 3] = -Real(6)*p_z_FF -
                                        Real(609)/Real(10)*p_z_F +
                                        Real(126)*V_ghost[(k + 5)*(d_num_species + 4) + d_num_species + 3] -
                                        Real(105)*V_ghost[(k + 4)*(d_num_species + 4) + d_num_species + 3] +
                                        Real(70)*V_ghost[(k + 3)*(d_num_species + 4) + d_num_species + 3] -
                                        Real(63)/Real(2)*V_ghost[(k + 2)*(d_num_species + 4) + d_num_species + 3] +
                                        Real(42)/Real(5)*V_ghost[(k + 1)*(d_num_species + 4) + d_num_species + 3] +
                                        Real(42)*Real(dx[2])*dV_dz[d_num_species + 3];
                                }
                                
                                /*
                                 * Compute the mixture density.
                                 */
                                
                                Real rho_ghost = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_ghost += V_ghost[k*(d_num_species + 4) + si];
                                }
                                
                                /*
                                 * Compute the mass fractions.
                                 */
                                
                                std::vector<Real> Y_ghost;
                                Y_ghost.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_ghost.push_back(V_ghost[k*(d_num_species + 4) + si]/rho_ghost);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions.
                                 */
                                
                                std::vector<const Real*> Y_ghost_ptr;
                                Y_ghost_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_ghost_ptr.push_back(&Y_ghost[si]);
                                }
                                
                                for(int si=0; si < d_num_species; si++)
                                {
                                    Q[si][idx_cell_rho_Y] = V_ghost[k*(d_num_species + 4) + si];
                                }
                                
                                Q[d_num_species][idx_cell_mom]     = rho_ghost*V_ghost[k*(d_num_species + 4) + d_num_species];
                                Q[d_num_species + 1][idx_cell_mom] = rho_ghost*V_ghost[k*(d_num_species + 4) + d_num_species + 1];
                                Q[d_num_species + 2][idx_cell_mom] = rho_ghost*V_ghost[k*(d_num_species + 4) + d_num_species + 2];
                                
                                const Real epsilon = d_equation_of_state_mixing_rules->
                                    getInternalEnergy(
                                        &rho_ghost,
                                        &V_ghost[k*(d_num_species + 4) + d_num_species + 3],
                                        Y_ghost_ptr);
                                
                                const Real E = rho_ghost*epsilon +
                                    half*(Q[d_num_species][idx_cell_mom]*Q[d_num_species][idx_cell_mom] +
                                        Q[d_num_species + 1][idx_cell_mom]*Q[d_num_species + 1][idx_cell_mom] +
                                        Q[d_num_species + 2][idx_cell_mom]*Q[d_num_species + 2][idx_cell_mom])/
                                        rho_ghost;
                                
                                Q[d_num_species + 3][idx_cell_E] = E;
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
                            << ": FlowModelBoundaryUtilitiesFourEqnConservative::fill3dFaceBoundaryData()\n"
                            << "Non-reflecting outflow BC doesn't support more than six ghost cells yet!");
                    }
                    
                    for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                    {
                        for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                        {
                            // Get the grid spacing.
                            const double* const dx = patch_geom->getDx();
                            
                            const int idx_cell_rho_Y_z_B = (i + num_subghosts_conservative_var[0][0]) +
                                (j + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                (interior_box_hi_idx[2] + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][1]*
                                    subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_rho_Y_z_BB = (i + num_subghosts_conservative_var[0][0]) +
                                (j + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                (interior_box_hi_idx[2] - 1 + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                    subghostcell_dims_conservative_var[0][1];
                            
                            const int idx_cell_rho_Y_z_BBB = (i + num_subghosts_conservative_var[0][0]) +
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
                            
                            std::vector<Real> rho_Y_z_B;
                            std::vector<Real> rho_Y_z_BB;
                            std::vector<Real> rho_Y_z_BBB;
                            rho_Y_z_B.reserve(d_num_species);
                            rho_Y_z_BB.reserve(d_num_species);
                            rho_Y_z_BBB.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                rho_Y_z_B.push_back(Q[si][idx_cell_rho_Y_z_B]);
                                rho_Y_z_BB.push_back(Q[si][idx_cell_rho_Y_z_BB]);
                                rho_Y_z_BBB.push_back(Q[si][idx_cell_rho_Y_z_BBB]);
                            }
                            
                            /*
                             * Compute the mixture density.
                             */
                            
                            Real rho_z_B   = Real(0);
                            Real rho_z_BB  = Real(0);
                            Real rho_z_BBB = Real(0);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                rho_z_B   += Q[si][idx_cell_rho_Y_z_B];
                                rho_z_BB  += Q[si][idx_cell_rho_Y_z_BB];
                                rho_z_BBB += Q[si][idx_cell_rho_Y_z_BBB];
                            }
                            
                            /*
                             * Compute the mass fractions.
                             */
                            
                            std::vector<Real> Y_z_B;
                            std::vector<Real> Y_z_BB;
                            std::vector<Real> Y_z_BBB;
                            Y_z_B.reserve(d_num_species);
                            Y_z_BB.reserve(d_num_species);
                            Y_z_BBB.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_z_B.push_back(Q[si][idx_cell_rho_Y_z_B]/rho_z_B);
                                Y_z_BB.push_back(Q[si][idx_cell_rho_Y_z_BB]/rho_z_BB);
                                Y_z_BBB.push_back(Q[si][idx_cell_rho_Y_z_BBB]/rho_z_BBB);
                            }
                            
                            /*
                             * Get the pointers to the mass fractions.
                             */
                            
                            std::vector<const Real*> Y_z_B_ptr;
                            std::vector<const Real*> Y_z_BB_ptr;
                            std::vector<const Real*> Y_z_BBB_ptr;
                            Y_z_B_ptr.reserve(d_num_species);
                            Y_z_BB_ptr.reserve(d_num_species);
                            Y_z_BBB_ptr.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_z_B_ptr.push_back(&Y_z_B[si]);
                                Y_z_BB_ptr.push_back(&Y_z_BB[si]);
                                Y_z_BBB_ptr.push_back(&Y_z_BBB[si]);
                            }
                            
                            const Real u_z_B   = Q[d_num_species][idx_cell_mom_z_B]/rho_z_B;
                            const Real u_z_BB  = Q[d_num_species][idx_cell_mom_z_BB]/rho_z_BB;
                            const Real u_z_BBB = Q[d_num_species][idx_cell_mom_z_BBB]/rho_z_BBB;
                            
                            const Real v_z_B   = Q[d_num_species + 1][idx_cell_mom_z_B]/rho_z_B;
                            const Real v_z_BB  = Q[d_num_species + 1][idx_cell_mom_z_BB]/rho_z_BB;
                            const Real v_z_BBB = Q[d_num_species + 1][idx_cell_mom_z_BBB]/rho_z_BBB;
                            
                            const Real w_z_B   = Q[d_num_species + 2][idx_cell_mom_z_B]/rho_z_B;
                            const Real w_z_BB  = Q[d_num_species + 2][idx_cell_mom_z_BB]/rho_z_BB;
                            const Real w_z_BBB = Q[d_num_species + 2][idx_cell_mom_z_BBB]/rho_z_BBB;
                            
                            const Real half = Real(1)/Real(2);
                            const Real epsilon_z_B   = Q[d_num_species + 3][idx_cell_E_z_B]/rho_z_B - half*(u_z_B*u_z_B + v_z_B*v_z_B + w_z_B*w_z_B);
                            const Real epsilon_z_BB  = Q[d_num_species + 3][idx_cell_E_z_BB]/rho_z_BB - half*(u_z_BB*u_z_BB + v_z_BB*v_z_BB + w_z_BB*w_z_BB);
                            const Real epsilon_z_BBB = Q[d_num_species + 3][idx_cell_E_z_BBB]/rho_z_BBB - half*(u_z_BBB*u_z_BBB + v_z_BBB*v_z_BBB + w_z_BBB*w_z_BBB);
                            
                            Real p_z_B = d_equation_of_state_mixing_rules->
                                getPressure(
                                    &rho_z_B,
                                    &epsilon_z_B,
                                    Y_z_B_ptr);
                            
                            Real p_z_BB = d_equation_of_state_mixing_rules->
                                getPressure(
                                    &rho_z_BB,
                                    &epsilon_z_BB,
                                    Y_z_BB_ptr);
                            
                            Real p_z_BBB = d_equation_of_state_mixing_rules->
                                getPressure(
                                    &rho_z_BBB,
                                    &epsilon_z_BBB,
                                    Y_z_BBB_ptr);
                            
                            /*
                             * Compute derivatives in z-direction.
                             */
                            
                            std::vector<Real> drho_Y_dz;
                            drho_Y_dz.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                drho_Y_dz.push_back((Q[si][idx_cell_rho_Y_z_BBB] - Real(4)*Q[si][idx_cell_rho_Y_z_BB] +
                                    Real(3)*Q[si][idx_cell_rho_Y_z_B])/(Real(2)*Real(dx[2])));
                            }
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
                                
                                const int idx_cell_rho_Y_x_R = (i + 1 + num_subghosts_conservative_var[0][0]) +
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
                                
                                /*
                                 * Compute the mixture density.
                                 */
                                
                                Real rho_x_R = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_x_R += Q[si][idx_cell_rho_Y_x_R];
                                }
                                
                                /*
                                 * Compute the mass fractions.
                                 */
                                
                                std::vector<Real> Y_x_R;
                                Y_x_R.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_x_R.push_back(Q[si][idx_cell_rho_Y_x_R]/rho_x_R);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions.
                                 */
                                
                                std::vector<const Real*> Y_x_R_ptr;
                                Y_x_R_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_x_R_ptr.push_back(&Y_x_R[si]);
                                }
                                
                                const Real u_x_R = Q[d_num_species][idx_cell_mom_x_R]/rho_x_R;
                                const Real v_x_R = Q[d_num_species + 1][idx_cell_mom_x_R]/rho_x_R;
                                const Real w_x_R = Q[d_num_species + 2][idx_cell_mom_x_R]/rho_x_R;
                                const Real epsilon_x_R = Q[d_num_species + 3][idx_cell_E_x_R]/rho_x_R - half*(u_x_R*u_x_R + v_x_R*v_x_R + w_x_R*w_x_R);
                                
                                Real p_x_R = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_x_R,
                                        &epsilon_x_R,
                                        Y_x_R_ptr);
                                
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
                                
                                const int idx_cell_rho_Y_x_L = (i - 1 + num_subghosts_conservative_var[0][0]) +
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
                                
                                /*
                                 * Compute the mixture density.
                                 */
                                
                                Real rho_x_L = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_x_L += Q[si][idx_cell_rho_Y_x_L];
                                }
                                
                                /*
                                 * Compute the mass fractions.
                                 */
                                
                                std::vector<Real> Y_x_L;
                                Y_x_L.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_x_L.push_back(Q[si][idx_cell_rho_Y_x_L]/rho_x_L);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions.
                                 */
                            
                                std::vector<const Real*> Y_x_L_ptr;
                                Y_x_L_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_x_L_ptr.push_back(&Y_x_L[si]);
                                }
                                
                                const Real u_x_L = Q[d_num_species][idx_cell_mom_x_L]/rho_x_L;
                                const Real v_x_L = Q[d_num_species + 1][idx_cell_mom_x_L]/rho_x_L;
                                const Real w_x_L = Q[d_num_species + 2][idx_cell_mom_x_L]/rho_x_L;
                                const Real epsilon_x_L = Q[d_num_species + 3][idx_cell_E_x_L]/rho_x_L - half*(u_x_L*u_x_L + v_x_L*v_x_L + w_x_L*w_x_L);
                                
                                Real p_x_L = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_x_L,
                                        &epsilon_x_L,
                                        Y_x_L_ptr);
                                
                                // One-sided derivatives.
                                du_dx = (u_z_B - u_x_L)/Real(dx[0]);
                                // dv_dx = (v_z_B - v_x_L)/Real(dx[0]);
                                dw_dx = (w_z_B - w_x_L)/Real(dx[0]);
                                dp_dx = (p_z_B - p_x_L)/Real(dx[0]);
                            }
                            else
                            {
                                const int idx_cell_rho_Y_x_L = (i - 1 + num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (interior_box_hi_idx[2] + num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_rho_Y_x_R = (i + 1 + num_subghosts_conservative_var[0][0]) +
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
                                
                                /*
                                 * Compute the mixture density.
                                 */
                                
                                Real rho_x_L = Real(0);
                                Real rho_x_R = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_x_L += Q[si][idx_cell_rho_Y_x_L];
                                    rho_x_R += Q[si][idx_cell_rho_Y_x_R];
                                }
                                
                                /*
                                 * Compute the mass fractions.
                                 */
                                
                                std::vector<Real> Y_x_L;
                                std::vector<Real> Y_x_R;
                                Y_x_L.reserve(d_num_species);
                                Y_x_R.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_x_L.push_back(Q[si][idx_cell_rho_Y_x_L]/rho_x_L);
                                    Y_x_R.push_back(Q[si][idx_cell_rho_Y_x_R]/rho_x_R);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions.
                                 */
                                
                                std::vector<const Real*> Y_x_L_ptr;
                                std::vector<const Real*> Y_x_R_ptr;
                                                                    
                                Y_x_L_ptr.reserve(d_num_species);
                                Y_x_R_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_x_L_ptr.push_back(&Y_x_L[si]);
                                    Y_x_R_ptr.push_back(&Y_x_R[si]);
                                }
                                
                                const Real u_x_L = Q[d_num_species][idx_cell_mom_x_L]/rho_x_L;
                                const Real u_x_R = Q[d_num_species][idx_cell_mom_x_R]/rho_x_R;

                                const Real v_x_L = Q[d_num_species + 1][idx_cell_mom_x_L]/rho_x_L;
                                const Real v_x_R = Q[d_num_species + 1][idx_cell_mom_x_R]/rho_x_R;

                                const Real w_x_L = Q[d_num_species + 2][idx_cell_mom_x_L]/rho_x_L;
                                const Real w_x_R = Q[d_num_species + 2][idx_cell_mom_x_R]/rho_x_R;
                                
                                const Real epsilon_x_L = Q[d_num_species + 3][idx_cell_E_x_L]/rho_x_L - half*(u_x_L*u_x_L + v_x_L*v_x_L + w_x_L*w_x_L);
                                const Real epsilon_x_R = Q[d_num_species + 3][idx_cell_E_x_R]/rho_x_R - half*(u_x_R*u_x_R + v_x_R*v_x_R + w_x_R*w_x_R);
                                
                                Real p_x_L = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_x_L,
                                        &epsilon_x_L,
                                        Y_x_L_ptr);
                            
                                Real p_x_R = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_x_R,
                                        &epsilon_x_R,
                                        Y_x_R_ptr);
                                
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
                                
                                const int idx_cell_rho_Y_y_T = (i + num_subghosts_conservative_var[0][0]) +
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
                                
                                /*
                                 * Compute the mixture density.
                                 */
                                
                                Real rho_y_T = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_y_T += Q[si][idx_cell_rho_Y_y_T];
                                }
                                
                                /*
                                 * Compute the mass fractions.
                                 */
                                
                                std::vector<Real> Y_y_T;
                                Y_y_T.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_y_T.push_back(Q[si][idx_cell_rho_Y_y_T]/rho_y_T);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions.
                                 */
                                
                                std::vector<const Real*> Y_y_T_ptr;
                                Y_y_T_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_y_T_ptr.push_back(&Y_y_T[si]);
                                }
                                
                                const Real u_y_T = Q[d_num_species][idx_cell_mom_y_T]/rho_y_T;
                                const Real v_y_T = Q[d_num_species + 1][idx_cell_mom_y_T]/rho_y_T;
                                const Real w_y_T = Q[d_num_species + 2][idx_cell_mom_y_T]/rho_y_T;
                                const Real epsilon_y_T = Q[d_num_species + 3][idx_cell_E_y_T]/rho_y_T - half*(u_y_T*u_y_T + v_y_T*v_y_T + w_y_T*w_y_T);
                                
                                Real p_y_T = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_y_T,
                                        &epsilon_y_T,
                                        Y_y_T_ptr);
                                
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
                                
                                const int idx_cell_rho_Y_y_B = (i + num_subghosts_conservative_var[0][0]) +
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
                                
                                /*
                                 * Compute the mixture density.
                                 */
                                
                                Real rho_y_B = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_y_B += Q[si][idx_cell_rho_Y_y_B];
                                }
                                
                                /*
                                 * Compute the mass fractions.
                                 */
                                
                                std::vector<Real> Y_y_B;
                                Y_y_B.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_y_B.push_back(Q[si][idx_cell_rho_Y_y_B]/rho_y_B);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions.
                                 */
                                
                                std::vector<const Real*> Y_y_B_ptr;
                                Y_y_B_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_y_B_ptr.push_back(&Y_y_B[si]);
                                }
                                
                                const Real u_y_B = Q[d_num_species][idx_cell_mom_y_B]/rho_y_B;
                                const Real v_y_B = Q[d_num_species + 1][idx_cell_mom_y_B]/rho_y_B;
                                const Real w_y_B = Q[d_num_species + 2][idx_cell_mom_y_B]/rho_y_B;
                                const Real epsilon_y_B = Q[d_num_species + 3][idx_cell_E_y_B]/rho_y_B - half*(u_y_B*u_y_B + v_y_B*v_y_B + w_y_B*w_y_B);
                                
                                Real p_y_B = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_y_B,
                                        &epsilon_y_B,
                                        Y_y_B_ptr);
                                
                                // One-sided derivatives.
                                // du_dy = (u_z_B - u_y_B)/Real(dx[1]);
                                dv_dy = (v_z_B - v_y_B)/Real(dx[1]);
                                dw_dy = (w_z_B - w_y_B)/Real(dx[1]);
                                dp_dy = (p_z_B - p_y_B)/Real(dx[1]);
                            }
                            else
                            {
                                const int idx_cell_rho_Y_y_B = (i + num_subghosts_conservative_var[0][0]) +
                                    (j - 1 + num_subghosts_conservative_var[0][1])*subghostcell_dims_conservative_var[0][0] +
                                    (interior_box_hi_idx[2] +  num_subghosts_conservative_var[0][2])*subghostcell_dims_conservative_var[0][0]*
                                        subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_rho_Y_y_T = (i + num_subghosts_conservative_var[0][0]) +
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
                                
                                /*
                                 * Compute the mixture density.
                                 */
                                
                                Real rho_y_B = Real(0);
                                Real rho_y_T = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_y_B += Q[si][idx_cell_rho_Y_y_B];
                                    rho_y_T += Q[si][idx_cell_rho_Y_y_T];
                                }
                                
                                /*
                                 * Compute the mass fractions.
                                 */
                                
                                std::vector<Real> Y_y_B;
                                std::vector<Real> Y_y_T;
                                Y_y_B.reserve(d_num_species);
                                Y_y_T.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_y_B.push_back(Q[si][idx_cell_rho_Y_y_B]/rho_y_B);
                                    Y_y_T.push_back(Q[si][idx_cell_rho_Y_y_T]/rho_y_T);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions.
                                 */
                                
                                std::vector<const Real*> Y_y_B_ptr;
                                std::vector<const Real*> Y_y_T_ptr;
                                Y_y_B_ptr.reserve(d_num_species);
                                Y_y_T_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_y_B_ptr.push_back(&Y_y_B[si]);
                                    Y_y_T_ptr.push_back(&Y_y_T[si]);
                                }
                                
                                const Real u_y_B = Q[d_num_species][idx_cell_mom_y_B]/rho_y_B;
                                const Real u_y_T = Q[d_num_species][idx_cell_mom_y_T]/rho_y_T;
                                
                                const Real v_y_B = Q[d_num_species + 1][idx_cell_mom_y_B]/rho_y_B;
                                const Real v_y_T = Q[d_num_species + 1][idx_cell_mom_y_T]/rho_y_T;
                                
                                const Real w_y_B = Q[d_num_species + 2][idx_cell_mom_y_B]/rho_y_B;
                                const Real w_y_T = Q[d_num_species + 2][idx_cell_mom_y_T]/rho_y_T;
                                
                                const Real epsilon_y_B = Q[d_num_species + 3][idx_cell_E_y_B]/rho_y_B - half*(u_y_B*u_y_B + v_y_B*v_y_B + w_y_B*w_y_B);
                                const Real epsilon_y_T = Q[d_num_species + 3][idx_cell_E_y_T]/rho_y_T - half*(u_y_T*u_y_T + v_y_T*v_y_T + w_y_T*w_y_T);
                                
                                Real p_y_B = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_y_B,
                                        &epsilon_y_B,
                                        Y_y_B_ptr);
                            
                                Real p_y_T = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_y_T,
                                        &epsilon_y_T,
                                        Y_y_T_ptr);
                                
                                // Central derivatives.
                                // du_dy = (u_y_T - u_y_B)/(Real(2)*Real(dx[1]));
                                dv_dy = (v_y_T - v_y_B)/(Real(2)*Real(dx[1]));
                                dw_dy = (w_y_T - w_y_B)/(Real(2)*Real(dx[1]));
                                dp_dy = (p_y_T - p_y_B)/(Real(2)*Real(dx[1]));
                            }
                            
                            // Compute sound speed.
                            
                            const Real Gamma_z_B = d_equation_of_state_mixing_rules->getGruneisenParameter(
                                &rho_z_B,
                                &p_z_B,
                                Y_z_B_ptr);
                            
                            const std::vector<Real> Psi_z_B = d_equation_of_state_mixing_rules->
                                getPressureDerivativeWithPartialDensities(
                                        &rho_z_B,
                                        &p_z_B,
                                        Y_z_B_ptr);
                            
                            Real c_z_B = Gamma_z_B*p_z_B/rho_z_B;
                            for (int si = 0; si < d_num_species; si++)
                            {
                                c_z_B += Y_z_B[si]*Psi_z_B[si];
                            }
                            c_z_B = std::sqrt(c_z_B);
                            
                            const Real lambda_1 = w_z_B - c_z_B;
                            
                            // Compute vector Lambda^(-1) * L.
                            
                            Real Lambda_inv_L[d_num_species + 4];
                            
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
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Lambda_inv_L[si + 3] = c_z_B*c_z_B*drho_Y_dz[si] - Y_z_B[si]*dp_dz;
                            }
                            Lambda_inv_L[d_num_species + 3] = dp_dz + rho_z_B*c_z_B*dw_dz;
                            
                            // Compute dV_dz.
                            
                            const Real c_sq_inv  = Real(1)/(c_z_B*c_z_B);
                            const Real rho_c_inv = Real(1)/(rho_z_B*c_z_B);
                            
                            Real dV_dz[d_num_species + 4];
                            
                            for (int si = 0; si < d_num_species; si++)
                            {
                                dV_dz[si] = half*c_sq_inv*Y_z_B[si]*(Lambda_inv_L[0] + Lambda_inv_L[d_num_species + 3]) +
                                    c_sq_inv*Lambda_inv_L[si + 3];
                            }
                            dV_dz[d_num_species]     = Lambda_inv_L[1];
                            dV_dz[d_num_species + 1] = Lambda_inv_L[2];
                            dV_dz[d_num_species + 2] = half*rho_c_inv*(-Lambda_inv_L[0] + Lambda_inv_L[d_num_species + 3]);
                            dV_dz[d_num_species + 3] = half*(Lambda_inv_L[0] + Lambda_inv_L[d_num_species + 3]);
                            
                            Real V_ghost[(d_num_species + 4)*num_ghosts_to_fill];
                            
                            for (int k = 0; k < num_ghosts_to_fill; k++)
                            {
                                const int idx_cell_rho_Y = (i + num_subghosts_conservative_var[0][0]) +
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
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        V_ghost[k*(d_num_species + 4) + si] = rho_Y_z_BB[si] + Real(2)*Real(dx[2])*dV_dz[si];
                                    }
                                    
                                    V_ghost[k*(d_num_species + 4) + d_num_species]     = u_z_BB + Real(2)*Real(dx[2])*dV_dz[d_num_species];
                                    V_ghost[k*(d_num_species + 4) + d_num_species + 1] = v_z_BB + Real(2)*Real(dx[2])*dV_dz[d_num_species + 1];
                                    V_ghost[k*(d_num_species + 4) + d_num_species + 2] = w_z_BB + Real(2)*Real(dx[2])*dV_dz[d_num_species + 2];
                                    V_ghost[k*(d_num_species + 4) + d_num_species + 3] = p_z_BB + Real(2)*Real(dx[2])*dV_dz[d_num_species + 3];
                                }
                                else if (k == 1)
                                {
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        V_ghost[k*(d_num_species + 4) + si] = -Real(2)*rho_Y_z_BB[si] - Real(3)*rho_Y_z_B[si] +
                                            Real(6)*V_ghost[(k - 1)*(d_num_species + 4) + si] - Real(6)*Real(dx[2])*dV_dz[si];
                                    }
                                    
                                    V_ghost[k*(d_num_species + 4) + d_num_species] = -Real(2)*u_z_BB - Real(3)*u_z_B +
                                        Real(6)*V_ghost[(k - 1)*(d_num_species + 4) + d_num_species] -
                                        Real(6)*Real(dx[2])*dV_dz[d_num_species];
                                    
                                    V_ghost[k*(d_num_species + 4) + d_num_species + 1] = -Real(2)*v_z_BB - Real(3)*v_z_B +
                                        Real(6)*V_ghost[(k - 1)*(d_num_species + 4) + d_num_species + 1] -
                                        Real(6)*Real(dx[2])*dV_dz[d_num_species + 1];
                                    
                                    V_ghost[k*(d_num_species + 4) + d_num_species + 2] = -Real(2)*w_z_BB - Real(3)*w_z_B +
                                        Real(6)*V_ghost[(k - 1)*(d_num_species + 4) + d_num_species + 2] -
                                        Real(6)*Real(dx[2])*dV_dz[d_num_species + 2];
                                    
                                    V_ghost[k*(d_num_species + 4) + d_num_species + 3] = -Real(2)*p_z_BB - Real(3)*p_z_B +
                                        Real(6)*V_ghost[(k - 1)*(d_num_species + 4) + d_num_species + 3] -
                                        Real(6)*Real(dx[2])*dV_dz[d_num_species + 3];
                                }
                                else if (k == 2)
                                {
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        V_ghost[k*(d_num_species + 4) + si] = Real(3)*rho_Y_z_BB[si] + Real(10)*rho_Y_z_B[si] -
                                            Real(18)*V_ghost[(k - 2)*(d_num_species + 4) + si] +
                                            Real(6)*V_ghost[(k - 1)*(d_num_species + 4) + si] +
                                            Real(12)*Real(dx[2])*dV_dz[si];
                                    }
                                    
                                    V_ghost[k*(d_num_species + 4) + d_num_species] = Real(3)*u_z_BB + Real(10)*u_z_B -
                                        Real(18)*V_ghost[(k - 2)*(d_num_species + 4) + d_num_species] +
                                        Real(6)*V_ghost[(k - 1)*(d_num_species + 4) + d_num_species] +
                                        Real(12)*Real(dx[2])*dV_dz[d_num_species];
                                    
                                    V_ghost[k*(d_num_species + 4) + d_num_species + 1] = Real(3)*v_z_BB + Real(10)*v_z_B -
                                        Real(18)*V_ghost[(k - 2)*(d_num_species + 4) + d_num_species + 1] +
                                        Real(6)*V_ghost[(k - 1)*(d_num_species + 4) + d_num_species + 1] +
                                        Real(12)*Real(dx[2])*dV_dz[d_num_species + 1];
                                    
                                    V_ghost[k*(d_num_species + 4) + d_num_species + 2] = Real(3)*w_z_BB + Real(10)*w_z_B -
                                        Real(18)*V_ghost[(k - 2)*(d_num_species + 4) + d_num_species + 2] +
                                        Real(6)*V_ghost[(k - 1)*(d_num_species + 4) + d_num_species + 2] +
                                        Real(12)*Real(dx[2])*dV_dz[d_num_species + 2];
                                    
                                    V_ghost[k*(d_num_species + 4) + d_num_species + 3] = Real(3)*p_z_BB + Real(10)*p_z_B -
                                        Real(18)*V_ghost[(k - 2)*(d_num_species + 4) + d_num_species + 3] +
                                        Real(6)*V_ghost[(k - 1)*(d_num_species + 4) + d_num_species + 3] +
                                        Real(12)*Real(dx[2])*dV_dz[d_num_species + 3];
                                }
                                else if (k == 3)
                                {
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        V_ghost[k*(d_num_species + 4) + si] = -Real(4)*rho_Y_z_BB[si] -
                                            Real(65)/Real(3)*rho_Y_z_B[si] +
                                            Real(40)*V_ghost[(k - 3)*(d_num_species + 4) + si] -
                                            Real(20)*V_ghost[(k - 2)*(d_num_species + 4) + si] +
                                            Real(20)/Real(3)*V_ghost[(k - 1)*(d_num_species + 4) + si] -
                                            Real(20)*Real(dx[2])*dV_dz[si];
                                    }
                                    
                                    V_ghost[k*(d_num_species + 4) + d_num_species] = -Real(4)*u_z_BB -
                                        Real(65)/Real(3)*u_z_B +
                                        Real(40)*V_ghost[(k - 3)*(d_num_species + 4) + d_num_species] -
                                        Real(20)*V_ghost[(k - 2)*(d_num_species + 4) + d_num_species] +
                                        Real(20)/Real(3)*V_ghost[(k - 1)*(d_num_species + 4) + d_num_species] -
                                        Real(20)*Real(dx[2])*dV_dz[d_num_species];
                                    
                                    V_ghost[k*(d_num_species + 4) + d_num_species + 1] = -Real(4)*v_z_BB -
                                        Real(65)/Real(3)*v_z_B +
                                        Real(40)*V_ghost[(k - 3)*(d_num_species + 4) + d_num_species + 1] -
                                        Real(20)*V_ghost[(k - 2)*(d_num_species + 4) + d_num_species + 1] +
                                        Real(20)/Real(3)*V_ghost[(k - 1)*(d_num_species + 4) + d_num_species + 1] -
                                        Real(20)*Real(dx[2])*dV_dz[d_num_species + 1];
                                    
                                    V_ghost[k*(d_num_species + 4) + d_num_species + 2] = -Real(4)*w_z_BB -
                                        Real(65)/Real(3)*w_z_B +
                                        Real(40)*V_ghost[(k - 3)*(d_num_species + 4) + d_num_species + 2] -
                                        Real(20)*V_ghost[(k - 2)*(d_num_species + 4) + d_num_species + 2] +
                                        Real(20)/Real(3)*V_ghost[(k - 1)*(d_num_species + 4) + d_num_species + 2] -
                                        Real(20)*Real(dx[2])*dV_dz[d_num_species + 2];
                                    
                                    V_ghost[k*(d_num_species + 4) + d_num_species + 3] = -Real(4)*p_z_BB -
                                        Real(65)/Real(3)*p_z_B +
                                        Real(40)*V_ghost[(k - 3)*(d_num_species + 4) + d_num_species + 3] -
                                        Real(20)*V_ghost[(k - 2)*(d_num_species + 4) + d_num_species + 3] +
                                        Real(20)/Real(3)*V_ghost[(k - 1)*(d_num_species + 4) + d_num_species + 3] -
                                        Real(20)*Real(dx[2])*dV_dz[d_num_species + 3];
                                }
                                else if (k == 4)
                                {
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        V_ghost[k*(d_num_species + 4) + si] = Real(5)*rho_Y_z_BB[si] +
                                            Real(77)/Real(2)*rho_Y_z_B[si] -
                                            Real(75)*V_ghost[(k - 4)*(d_num_species + 4) + si] +
                                            Real(50)*V_ghost[(k - 3)*(d_num_species + 4) + si] -
                                            Real(25)*V_ghost[(k - 2)*(d_num_species + 4) + si] +
                                            Real(15)/Real(2)*V_ghost[(k - 1)*(d_num_species + 4) + si] +
                                            Real(30)*Real(dx[2])*dV_dz[si];
                                    }
                                    
                                    V_ghost[k*(d_num_species + 4) + d_num_species] = Real(5)*u_z_BB +
                                        Real(77)/Real(2)*u_z_B -
                                        Real(75)*V_ghost[(k - 4)*(d_num_species + 4) + d_num_species] +
                                        Real(50)*V_ghost[(k - 3)*(d_num_species + 4) + d_num_species] -
                                        Real(25)*V_ghost[(k - 2)*(d_num_species + 4) + d_num_species] +
                                        Real(15)/Real(2)*V_ghost[(k - 1)*(d_num_species + 4) + d_num_species] +
                                        Real(30)*Real(dx[2])*dV_dz[d_num_species];
                                    
                                    V_ghost[k*(d_num_species + 4) + d_num_species + 1] = Real(5)*v_z_BB +
                                        Real(77)/Real(2)*v_z_B -
                                        Real(75)*V_ghost[(k - 4)*(d_num_species + 4) + d_num_species + 1] +
                                        Real(50)*V_ghost[(k - 3)*(d_num_species + 4) + d_num_species + 1] -
                                        Real(25)*V_ghost[(k - 2)*(d_num_species + 4) + d_num_species + 1] +
                                        Real(15)/Real(2)*V_ghost[(k - 1)*(d_num_species + 4) + d_num_species + 1] +
                                        Real(30)*Real(dx[2])*dV_dz[d_num_species + 1];
                                    
                                    V_ghost[k*(d_num_species + 4) + d_num_species + 2] = Real(5)*w_z_BB +
                                        Real(77)/Real(2)*w_z_B -
                                        Real(75)*V_ghost[(k - 4)*(d_num_species + 4) + d_num_species + 2] +
                                        Real(50)*V_ghost[(k - 3)*(d_num_species + 4) + d_num_species + 2] -
                                        Real(25)*V_ghost[(k - 2)*(d_num_species + 4) + d_num_species + 2] +
                                        Real(15)/Real(2)*V_ghost[(k - 1)*(d_num_species + 4) + d_num_species + 2] +
                                        Real(30)*Real(dx[2])*dV_dz[d_num_species + 2];
                                    
                                    V_ghost[k*(d_num_species + 4) + d_num_species + 3] = Real(5)*p_z_BB +
                                        Real(77)/Real(2)*p_z_B -
                                        Real(75)*V_ghost[(k - 4)*(d_num_species + 4) + d_num_species + 3] +
                                        Real(50)*V_ghost[(k - 3)*(d_num_species + 4) + d_num_species + 3] -
                                        Real(25)*V_ghost[(k - 2)*(d_num_species + 4) + d_num_species + 3] +
                                        Real(15)/Real(2)*V_ghost[(k - 1)*(d_num_species + 4) + d_num_species + 3] +
                                        Real(30)*Real(dx[2])*dV_dz[d_num_species + 3];
                                }
                                else if (k == 5)
                                {
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        V_ghost[k*(d_num_species + 4) + si] = -Real(6)*rho_Y_z_BB[si] -
                                            Real(609)/Real(10)*rho_Y_z_B[si] +
                                            Real(126)*V_ghost[(k - 5)*(d_num_species + 4) + si] -
                                            Real(105)*V_ghost[(k - 4)*(d_num_species + 4) + si] +
                                            Real(70)*V_ghost[(k - 3)*(d_num_species + 4) + si] -
                                            Real(63)/Real(2)*V_ghost[(k - 2)*(d_num_species + 4) + si] +
                                            Real(42)/Real(5)*V_ghost[(k - 1)*(d_num_species + 4) + si] -
                                            Real(42)*Real(dx[2])*dV_dz[si];
                                    }
                                    
                                    V_ghost[k*(d_num_species + 4) + d_num_species] = -Real(6)*u_z_BB -
                                        Real(609)/Real(10)*u_z_B +
                                        Real(126)*V_ghost[(k - 5)*(d_num_species + 4) + d_num_species] -
                                        Real(105)*V_ghost[(k - 4)*(d_num_species + 4) + d_num_species] +
                                        Real(70)*V_ghost[(k - 3)*(d_num_species + 4) + d_num_species] -
                                        Real(63)/Real(2)*V_ghost[(k - 2)*(d_num_species + 4) + d_num_species] +
                                        Real(42)/Real(5)*V_ghost[(k - 1)*(d_num_species + 4) + d_num_species] -
                                        Real(42)*Real(dx[2])*dV_dz[d_num_species];
                                    
                                    V_ghost[k*(d_num_species + 4) + d_num_species + 1] = -Real(6)*v_z_BB -
                                        Real(609)/Real(10)*v_z_B +
                                        Real(126)*V_ghost[(k - 5)*(d_num_species + 4) + d_num_species + 1] -
                                        Real(105)*V_ghost[(k - 4)*(d_num_species + 4) + d_num_species + 1] +
                                        Real(70)*V_ghost[(k - 3)*(d_num_species + 4) + d_num_species + 1] -
                                        Real(63)/Real(2)*V_ghost[(k - 2)*(d_num_species + 4) + d_num_species + 1] +
                                        Real(42)/Real(5)*V_ghost[(k - 1)*(d_num_species + 4) + d_num_species + 1] -
                                        Real(42)*Real(dx[2])*dV_dz[d_num_species + 1];
                                    
                                    V_ghost[k*(d_num_species + 4) + d_num_species + 2] = -Real(6)*w_z_BB -
                                        Real(609)/Real(10)*w_z_B +
                                        Real(126)*V_ghost[(k - 5)*(d_num_species + 4) + d_num_species + 2] -
                                        Real(105)*V_ghost[(k - 4)*(d_num_species + 4) + d_num_species + 2] +
                                        Real(70)*V_ghost[(k - 3)*(d_num_species + 4) + d_num_species + 2] -
                                        Real(63)/Real(2)*V_ghost[(k - 2)*(d_num_species + 4) + d_num_species + 2] +
                                        Real(42)/Real(5)*V_ghost[(k - 1)*(d_num_species + 4) + d_num_species + 2] -
                                        Real(42)*Real(dx[2])*dV_dz[d_num_species + 2];
                                    
                                    V_ghost[k*(d_num_species + 4) + d_num_species + 3] = -Real(6)*p_z_BB -
                                        Real(609)/Real(10)*p_z_B +
                                        Real(126)*V_ghost[(k - 5)*(d_num_species + 4) + d_num_species + 3] -
                                        Real(105)*V_ghost[(k - 4)*(d_num_species + 4) + d_num_species + 3] +
                                        Real(70)*V_ghost[(k - 3)*(d_num_species + 4) + d_num_species + 3] -
                                        Real(63)/Real(2)*V_ghost[(k - 2)*(d_num_species + 4) + d_num_species + 3] +
                                        Real(42)/Real(5)*V_ghost[(k - 1)*(d_num_species + 4) + d_num_species + 3] -
                                        Real(42)*Real(dx[2])*dV_dz[d_num_species + 3];
                                }
                                
                                /*
                                 * Compute the mixture density.
                                 */
                                
                                Real rho_ghost = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_ghost += V_ghost[k*(d_num_species + 4) + si];
                                }
                                
                                /*
                                 * Compute the mass fractions.
                                 */
                                
                                std::vector<Real> Y_ghost;
                                Y_ghost.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_ghost.push_back(V_ghost[k*(d_num_species + 4) + si]/rho_ghost);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions.
                                 */
                                
                                std::vector<const Real*> Y_ghost_ptr;
                                Y_ghost_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_ghost_ptr.push_back(&Y_ghost[si]);
                                }
                                
                                for(int si=0; si < d_num_species; si++)
                                {
                                    Q[si][idx_cell_rho_Y] = V_ghost[k*(d_num_species + 4) + si];
                                }
                                
                                Q[d_num_species][idx_cell_mom]     = rho_ghost*V_ghost[k*(d_num_species + 4) + d_num_species];
                                Q[d_num_species + 1][idx_cell_mom] = rho_ghost*V_ghost[k*(d_num_species + 4) + d_num_species + 1];
                                Q[d_num_species + 2][idx_cell_mom] = rho_ghost*V_ghost[k*(d_num_species + 4) + d_num_species + 2];
                                
                                const Real epsilon = d_equation_of_state_mixing_rules->
                                    getInternalEnergy(
                                        &rho_ghost,
                                        &V_ghost[k*(d_num_species + 4) + d_num_species + 3],
                                        Y_ghost_ptr);
                                
                                const Real E = rho_ghost*epsilon +
                                    half*(Q[d_num_species][idx_cell_mom]*Q[d_num_species][idx_cell_mom] +
                                        Q[d_num_species + 1][idx_cell_mom]*Q[d_num_species + 1][idx_cell_mom] +
                                        Q[d_num_species + 2][idx_cell_mom]*Q[d_num_species + 2][idx_cell_mom])/
                                        rho_ghost;
                                
                                Q[d_num_species + 3][idx_cell_E] = E;
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
            TBOX_ERROR("FlowModelBoundaryUtilitiesFourEqnConservative::fill3dFaceBoundaryData\n"
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
FlowModelBoundaryUtilitiesFourEqnConservative::fill3dEdgeBoundaryData(
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
                
                if (bdry_edge_conds[edge_loc] == BDRY_COND::FLOW_MODEL::XADIABATIC_NO_SLIP)
                {
                    for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                    {
                        for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                        {
                            for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                            {
                                const int idx_cell_rho_Y = (i + num_subghosts_conservative_var[0][0]) +
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
                                
                                int idx_cell_pivot_rho_Y = idx_cell_rho_Y;
                                int idx_cell_pivot_mom = idx_cell_mom;
                                int idx_cell_pivot_E = idx_cell_E;
                                
                                if (face_loc_0 == BDRY_LOC::XLO)
                                {
                                    idx_cell_pivot_rho_Y = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
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
                                    idx_cell_pivot_rho_Y = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
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
                                 * Compute the mixture density of the pivot.
                                 */
                                
                                Real rho_pivot = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_pivot += Q[si][idx_cell_pivot_rho_Y];
                                }
                                
                                /*
                                 * Compute the mass fractions of the pivot.
                                 */
                                
                                std::vector<Real> Y_pivot;
                                Y_pivot.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_pivot.push_back(Q[si][idx_cell_pivot_rho_Y]/rho_pivot);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions of the pivot.
                                 */
                                
                                std::vector<const Real*> Y_pivot_ptr;
                                Y_pivot_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_pivot_ptr.push_back(&Y_pivot[si]);
                                }
                                
                                /*
                                 * Set the values for partial densities and momentum.
                                 */
                                
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Q[si][idx_cell_rho_Y] = Q[si][idx_cell_pivot_rho_Y];
                                }
                                Q[d_num_species][idx_cell_mom] = -Q[d_num_species][idx_cell_pivot_mom] +
                                    Real(2)*rho_pivot*d_bdry_face_adiabatic_no_slip_vel[face_loc_0*3];
                                Q[d_num_species + 1][idx_cell_mom] = -Q[d_num_species + 1][idx_cell_pivot_mom] +
                                    Real(2)*rho_pivot*d_bdry_face_adiabatic_no_slip_vel[face_loc_0*3 + 1];
                                Q[d_num_species + 2][idx_cell_mom] = -Q[d_num_species + 2][idx_cell_pivot_mom] +
                                    Real(2)*rho_pivot*d_bdry_face_adiabatic_no_slip_vel[face_loc_0*3 + 2];
                                
                                /*
                                 * Set the values for total internal energy.
                                 */
                                
                                Real epsilon_pivot = (Q[d_num_species + 3][idx_cell_pivot_E] -
                                    Real(1)/Real(2)*(Q[d_num_species][idx_cell_pivot_mom]*Q[d_num_species][idx_cell_pivot_mom] +
                                    Q[d_num_species + 1][idx_cell_pivot_mom]*Q[d_num_species + 1][idx_cell_pivot_mom] +
                                    Q[d_num_species + 2][idx_cell_pivot_mom]*Q[d_num_species + 2][idx_cell_pivot_mom])/
                                    rho_pivot)/rho_pivot;
                                
                                Real p_pivot = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_pivot,
                                        &epsilon_pivot,
                                        Y_pivot_ptr);
                                
                                Real T_pivot = d_equation_of_state_mixing_rules->
                                    getTemperature(
                                        &rho_pivot,
                                        &p_pivot,
                                        Y_pivot_ptr);
                                
                                Real T = T_pivot;
                                
                                Real epsilon = d_equation_of_state_mixing_rules->
                                    getInternalEnergyFromTemperature(
                                        &rho_pivot,
                                        &T,
                                        Y_pivot_ptr);
                                
                                Real E = rho_pivot*epsilon +
                                    Real(1)/Real(2)*(Q[d_num_species][idx_cell_mom]*Q[d_num_species][idx_cell_mom] + 
                                    Q[d_num_species + 1][idx_cell_mom]*Q[d_num_species + 1][idx_cell_mom] +
                                    Q[d_num_species + 2][idx_cell_mom]*Q[d_num_species + 2][idx_cell_mom])/
                                    rho_pivot;
                                
                                Q[d_num_species + 3][idx_cell_E] = E;
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
                                const int idx_cell_rho_Y = (i + num_subghosts_conservative_var[0][0]) +
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
                                
                                int idx_cell_pivot_rho_Y = idx_cell_rho_Y;
                                int idx_cell_pivot_mom = idx_cell_mom;
                                int idx_cell_pivot_E = idx_cell_E;
                                
                                if (face_loc_1 == BDRY_LOC::YLO)
                                {
                                    idx_cell_pivot_rho_Y = (i + num_subghosts_conservative_var[0][0]) +
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
                                    idx_cell_pivot_rho_Y = (i + num_subghosts_conservative_var[0][0]) +
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
                                 * Compute the mixture density of the pivot.
                                 */
                                
                                Real rho_pivot = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_pivot += Q[si][idx_cell_pivot_rho_Y];
                                }
                                
                                /*
                                 * Compute the mass fractions of the pivot.
                                 */
                                
                                std::vector<Real> Y_pivot;
                                Y_pivot.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_pivot.push_back(Q[si][idx_cell_pivot_rho_Y]/rho_pivot);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions of the pivot.
                                 */
                                
                                std::vector<const Real*> Y_pivot_ptr;
                                Y_pivot_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_pivot_ptr.push_back(&Y_pivot[si]);
                                }
                                
                                /*
                                 * Set the values for partial densities and momentum.
                                 */
                                
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Q[si][idx_cell_rho_Y] = Q[si][idx_cell_pivot_rho_Y];
                                }
                                Q[d_num_species][idx_cell_mom] = -Q[d_num_species][idx_cell_pivot_mom] +
                                    Real(2)*rho_pivot*d_bdry_face_adiabatic_no_slip_vel[face_loc_1*3];
                                Q[d_num_species + 1][idx_cell_mom] = -Q[d_num_species + 1][idx_cell_pivot_mom] +
                                    Real(2)*rho_pivot*d_bdry_face_adiabatic_no_slip_vel[face_loc_1*3 + 1];
                                Q[d_num_species + 2][idx_cell_mom] = -Q[d_num_species + 2][idx_cell_pivot_mom] +
                                    Real(2)*rho_pivot*d_bdry_face_adiabatic_no_slip_vel[face_loc_1*3 + 2];
                                
                                /*
                                 * Set the values for total internal energy.
                                 */
                                
                                Real epsilon_pivot = (Q[d_num_species + 3][idx_cell_pivot_E] -
                                    Real(1)/Real(2)*(Q[d_num_species][idx_cell_pivot_mom]*Q[d_num_species][idx_cell_pivot_mom] +
                                    Q[d_num_species + 1][idx_cell_pivot_mom]*Q[d_num_species + 1][idx_cell_pivot_mom] +
                                    Q[d_num_species + 2][idx_cell_pivot_mom]*Q[d_num_species + 2][idx_cell_pivot_mom])/
                                    rho_pivot)/rho_pivot;
                                
                                Real p_pivot = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_pivot,
                                        &epsilon_pivot,
                                        Y_pivot_ptr);
                                
                                Real T_pivot = d_equation_of_state_mixing_rules->
                                    getTemperature(
                                        &rho_pivot,
                                        &p_pivot,
                                        Y_pivot_ptr);
                                
                                Real T = T_pivot;
                                
                                Real epsilon = d_equation_of_state_mixing_rules->
                                    getInternalEnergyFromTemperature(
                                        &rho_pivot,
                                        &T,
                                        Y_pivot_ptr);
                                
                                Real E = rho_pivot*epsilon +
                                    Real(1)/Real(2)*(Q[d_num_species][idx_cell_mom]*Q[d_num_species][idx_cell_mom] + 
                                    Q[d_num_species + 1][idx_cell_mom]*Q[d_num_species + 1][idx_cell_mom] +
                                    Q[d_num_species + 2][idx_cell_mom]*Q[d_num_species + 2][idx_cell_mom])/
                                    rho_pivot;
                                
                                Q[d_num_species + 3][idx_cell_E] = E;
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
                                const int idx_cell_rho_Y = (i + num_subghosts_conservative_var[0][0]) +
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
                                
                                int idx_cell_pivot_rho_Y = idx_cell_rho_Y;
                                int idx_cell_pivot_mom = idx_cell_mom;
                                int idx_cell_pivot_E = idx_cell_E;
                                
                                if (face_loc_2 == BDRY_LOC::ZLO)
                                {
                                    idx_cell_pivot_rho_Y = (i + num_subghosts_conservative_var[0][0]) +
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
                                    idx_cell_pivot_rho_Y = (i + num_subghosts_conservative_var[0][0]) +
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
                                 * Compute the mixture density of the pivot.
                                 */
                                
                                Real rho_pivot = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_pivot += Q[si][idx_cell_pivot_rho_Y];
                                }
                                
                                /*
                                 * Compute the mass fractions of the pivot.
                                 */
                                
                                std::vector<Real> Y_pivot;
                                Y_pivot.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_pivot.push_back(Q[si][idx_cell_pivot_rho_Y]/rho_pivot);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions of the pivot.
                                 */
                                
                                std::vector<const Real*> Y_pivot_ptr;
                                Y_pivot_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_pivot_ptr.push_back(&Y_pivot[si]);
                                }
                                
                                /*
                                 * Set the values for partial densities and momentum.
                                 */
                                
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Q[si][idx_cell_rho_Y] = Q[si][idx_cell_pivot_rho_Y];
                                }
                                Q[d_num_species][idx_cell_mom] = -Q[d_num_species][idx_cell_pivot_mom] +
                                    Real(2)*rho_pivot*d_bdry_face_adiabatic_no_slip_vel[face_loc_2*3];
                                Q[d_num_species + 1][idx_cell_mom] = -Q[d_num_species + 1][idx_cell_pivot_mom] +
                                    Real(2)*rho_pivot*d_bdry_face_adiabatic_no_slip_vel[face_loc_2*3 + 1];
                                Q[d_num_species + 2][idx_cell_mom] = -Q[d_num_species + 2][idx_cell_pivot_mom] +
                                    Real(2)*rho_pivot*d_bdry_face_adiabatic_no_slip_vel[face_loc_2*3 + 2];
                                
                                /*
                                 * Set the values for total internal energy.
                                 */
                                
                                Real epsilon_pivot = (Q[d_num_species + 3][idx_cell_pivot_E] -
                                    Real(1)/Real(2)*(Q[d_num_species][idx_cell_pivot_mom]*Q[d_num_species][idx_cell_pivot_mom] +
                                    Q[d_num_species + 1][idx_cell_pivot_mom]*Q[d_num_species + 1][idx_cell_pivot_mom] +
                                    Q[d_num_species + 2][idx_cell_pivot_mom]*Q[d_num_species + 2][idx_cell_pivot_mom])/
                                    rho_pivot)/rho_pivot;
                                
                                Real p_pivot = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_pivot,
                                        &epsilon_pivot,
                                        Y_pivot_ptr);
                                
                                Real T_pivot = d_equation_of_state_mixing_rules->
                                    getTemperature(
                                        &rho_pivot,
                                        &p_pivot,
                                        Y_pivot_ptr);
                                
                                Real T = T_pivot;
                                
                                Real epsilon = d_equation_of_state_mixing_rules->
                                    getInternalEnergyFromTemperature(
                                        &rho_pivot,
                                        &T,
                                        Y_pivot_ptr);
                                
                                Real E = rho_pivot*epsilon +
                                    Real(1)/Real(2)*(Q[d_num_species][idx_cell_mom]*Q[d_num_species][idx_cell_mom] + 
                                    Q[d_num_species + 1][idx_cell_mom]*Q[d_num_species + 1][idx_cell_mom] +
                                    Q[d_num_species + 2][idx_cell_mom]*Q[d_num_species + 2][idx_cell_mom])/
                                    rho_pivot;
                                
                                Q[d_num_species + 3][idx_cell_E] = E;
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
                                const int idx_cell_rho_Y = (i + num_subghosts_conservative_var[0][0]) +
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
                                
                                int idx_cell_pivot_rho_Y = idx_cell_rho_Y;
                                int idx_cell_pivot_mom = idx_cell_mom;
                                int idx_cell_pivot_E = idx_cell_E;
                                
                                if (face_loc_0 == BDRY_LOC::XLO)
                                {
                                    idx_cell_pivot_rho_Y = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
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
                                    idx_cell_pivot_rho_Y = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
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
                                 * Compute the mixture density of the pivot.
                                 */
                                
                                Real rho_pivot = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_pivot += Q[si][idx_cell_pivot_rho_Y];
                                }
                                
                                /*
                                 * Compute the mass fractions of the pivot.
                                 */
                                
                                std::vector<Real> Y_pivot;
                                Y_pivot.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_pivot.push_back(Q[si][idx_cell_pivot_rho_Y]/rho_pivot);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions of the pivot.
                                 */
                                
                                std::vector<const Real*> Y_pivot_ptr;
                                Y_pivot_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_pivot_ptr.push_back(&Y_pivot[si]);
                                }
                                
                                /*
                                 * Set the values for partial densities, momentum and total internal energy.
                                 */
                                
                                Real epsilon_pivot = (Q[d_num_species + 3][idx_cell_pivot_E] -
                                    Real(1)/Real(2)*(Q[d_num_species][idx_cell_pivot_mom]*Q[d_num_species][idx_cell_pivot_mom] +
                                    Q[d_num_species + 1][idx_cell_pivot_mom]*Q[d_num_species + 1][idx_cell_pivot_mom] +
                                    Q[d_num_species + 2][idx_cell_pivot_mom]*Q[d_num_species + 2][idx_cell_pivot_mom])/
                                    rho_pivot)/rho_pivot;
                                
                                Real p_pivot = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_pivot,
                                        &epsilon_pivot,
                                        Y_pivot_ptr);
                                
                                Real p = p_pivot;
                                
                                Real T_pivot = d_equation_of_state_mixing_rules->
                                    getTemperature(
                                        &rho_pivot,
                                        &p_pivot,
                                        Y_pivot_ptr);
                                
                                Real T = -T_pivot + Real(2)*d_bdry_face_isothermal_no_slip_T[face_loc_0];
                                
                                Real rho = d_equation_of_state_mixing_rules->
                                    getMixtureDensity(
                                        &p,
                                        &T,
                                        Y_pivot_ptr);
                                
                                Real u = -Q[d_num_species][idx_cell_pivot_mom]/rho_pivot +
                                    Real(2)*d_bdry_edge_isothermal_no_slip_vel[face_loc_0*3];
                                Real v = -Q[d_num_species + 1][idx_cell_pivot_mom]/rho_pivot +
                                    Real(2)*d_bdry_edge_isothermal_no_slip_vel[face_loc_0*3 + 1];
                                Real w = -Q[d_num_species + 2][idx_cell_pivot_mom]/rho_pivot +
                                    Real(2)*d_bdry_edge_isothermal_no_slip_vel[face_loc_0*3 + 2];
                                
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Q[si][idx_cell_rho_Y] = rho*Y_pivot[si];
                                }
                                Q[d_num_species][idx_cell_mom] = rho*u;
                                Q[d_num_species + 1][idx_cell_mom] = rho*v;
                                Q[d_num_species + 2][idx_cell_mom] = rho*w;
                                
                                Real epsilon = d_equation_of_state_mixing_rules->
                                    getInternalEnergyFromTemperature(
                                        &rho,
                                        &T,
                                        Y_pivot_ptr);
                                
                                Real E = rho*epsilon +
                                    Real(1)/Real(2)*(Q[d_num_species][idx_cell_mom]*Q[d_num_species][idx_cell_mom] +
                                    Q[d_num_species + 1][idx_cell_mom]*Q[d_num_species + 1][idx_cell_mom] +
                                    Q[d_num_species + 2][idx_cell_mom]*Q[d_num_species + 2][idx_cell_mom])/rho;
                                
                                Q[d_num_species + 3][idx_cell_E] = E;
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
                                const int idx_cell_rho_Y = (i + num_subghosts_conservative_var[0][0]) +
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
                                
                                int idx_cell_pivot_rho_Y = idx_cell_rho_Y;
                                int idx_cell_pivot_mom = idx_cell_mom;
                                int idx_cell_pivot_E = idx_cell_E;
                                
                                if (face_loc_1 == BDRY_LOC::YLO)
                                {
                                    idx_cell_pivot_rho_Y = (i + num_subghosts_conservative_var[0][0]) +
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
                                    idx_cell_pivot_rho_Y = (i + num_subghosts_conservative_var[0][0]) +
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
                                 * Compute the mixture density of the pivot.
                                 */
                                
                                Real rho_pivot = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_pivot += Q[si][idx_cell_pivot_rho_Y];
                                }
                                
                                /*
                                 * Compute the mass fractions of the pivot.
                                 */
                                
                                std::vector<Real> Y_pivot;
                                Y_pivot.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_pivot.push_back(Q[si][idx_cell_pivot_rho_Y]/rho_pivot);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions of the pivot.
                                 */
                                
                                std::vector<const Real*> Y_pivot_ptr;
                                Y_pivot_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_pivot_ptr.push_back(&Y_pivot[si]);
                                }
                                
                                /*
                                 * Set the values for partial densities, momentum and total internal energy.
                                 */
                                
                                Real epsilon_pivot = (Q[d_num_species + 3][idx_cell_pivot_E] -
                                    Real(1)/Real(2)*(Q[d_num_species][idx_cell_pivot_mom]*Q[d_num_species][idx_cell_pivot_mom] +
                                    Q[d_num_species + 1][idx_cell_pivot_mom]*Q[d_num_species + 1][idx_cell_pivot_mom] +
                                    Q[d_num_species + 2][idx_cell_pivot_mom]*Q[d_num_species + 2][idx_cell_pivot_mom])/
                                    rho_pivot)/rho_pivot;
                                
                                Real p_pivot = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_pivot,
                                        &epsilon_pivot,
                                        Y_pivot_ptr);
                                
                                Real p = p_pivot;
                                
                                Real T_pivot = d_equation_of_state_mixing_rules->
                                    getTemperature(
                                        &rho_pivot,
                                        &p_pivot,
                                        Y_pivot_ptr);
                                
                                Real T = -T_pivot + Real(2)*d_bdry_face_isothermal_no_slip_T[face_loc_1];
                                
                                Real rho = d_equation_of_state_mixing_rules->
                                    getMixtureDensity(
                                        &p,
                                        &T,
                                        Y_pivot_ptr);
                                
                                Real u = -Q[d_num_species][idx_cell_pivot_mom]/rho_pivot +
                                    Real(2)*d_bdry_edge_isothermal_no_slip_vel[face_loc_1*3];
                                Real v = -Q[d_num_species + 1][idx_cell_pivot_mom]/rho_pivot +
                                    Real(2)*d_bdry_edge_isothermal_no_slip_vel[face_loc_1*3 + 1];
                                Real w = -Q[d_num_species + 2][idx_cell_pivot_mom]/rho_pivot +
                                    Real(2)*d_bdry_edge_isothermal_no_slip_vel[face_loc_1*3 + 2];
                                
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Q[si][idx_cell_rho_Y] = rho*Y_pivot[si];
                                }
                                Q[d_num_species][idx_cell_mom] = rho*u;
                                Q[d_num_species + 1][idx_cell_mom] = rho*v;
                                Q[d_num_species + 2][idx_cell_mom] = rho*w;
                                
                                Real epsilon = d_equation_of_state_mixing_rules->
                                    getInternalEnergyFromTemperature(
                                        &rho,
                                        &T,
                                        Y_pivot_ptr);
                                
                                Real E = rho*epsilon +
                                    Real(1)/Real(2)*(Q[d_num_species][idx_cell_mom]*Q[d_num_species][idx_cell_mom] +
                                    Q[d_num_species + 1][idx_cell_mom]*Q[d_num_species + 1][idx_cell_mom] +
                                    Q[d_num_species + 2][idx_cell_mom]*Q[d_num_species + 2][idx_cell_mom])/rho;
                                
                                Q[d_num_species + 3][idx_cell_E] = E;
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
                                const int idx_cell_rho_Y = (i + num_subghosts_conservative_var[0][0]) +
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
                                
                                int idx_cell_pivot_rho_Y = idx_cell_rho_Y;
                                int idx_cell_pivot_mom = idx_cell_mom;
                                int idx_cell_pivot_E = idx_cell_E;
                                
                                if (face_loc_2 == BDRY_LOC::ZLO)
                                {
                                    idx_cell_pivot_rho_Y = (i + num_subghosts_conservative_var[0][0]) +
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
                                    idx_cell_pivot_rho_Y = (i + num_subghosts_conservative_var[0][0]) +
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
                                 * Compute the mixture density of the pivot.
                                 */
                                
                                Real rho_pivot = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_pivot += Q[si][idx_cell_pivot_rho_Y];
                                }
                                
                                /*
                                 * Compute the mass fractions of the pivot.
                                 */
                                
                                std::vector<Real> Y_pivot;
                                Y_pivot.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_pivot.push_back(Q[si][idx_cell_pivot_rho_Y]/rho_pivot);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions of the pivot.
                                 */
                                
                                std::vector<const Real*> Y_pivot_ptr;
                                Y_pivot_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_pivot_ptr.push_back(&Y_pivot[si]);
                                }
                                
                                /*
                                 * Set the values for partial densities, momentum and total internal energy.
                                 */
                                
                                Real epsilon_pivot = (Q[d_num_species + 3][idx_cell_pivot_E] -
                                    Real(1)/Real(2)*(Q[d_num_species][idx_cell_pivot_mom]*Q[d_num_species][idx_cell_pivot_mom] +
                                    Q[d_num_species + 1][idx_cell_pivot_mom]*Q[d_num_species + 1][idx_cell_pivot_mom] +
                                    Q[d_num_species + 2][idx_cell_pivot_mom]*Q[d_num_species + 2][idx_cell_pivot_mom])/
                                    rho_pivot)/rho_pivot;
                                
                                Real p_pivot = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_pivot,
                                        &epsilon_pivot,
                                        Y_pivot_ptr);
                                
                                Real p = p_pivot;
                                
                                Real T_pivot = d_equation_of_state_mixing_rules->
                                    getTemperature(
                                        &rho_pivot,
                                        &p_pivot,
                                        Y_pivot_ptr);
                                
                                Real T = -T_pivot + Real(2)*d_bdry_face_isothermal_no_slip_T[face_loc_2];
                                
                                Real rho = d_equation_of_state_mixing_rules->
                                    getMixtureDensity(
                                        &p,
                                        &T,
                                        Y_pivot_ptr);
                                
                                Real u = -Q[d_num_species][idx_cell_pivot_mom]/rho_pivot +
                                    Real(2)*d_bdry_edge_isothermal_no_slip_vel[face_loc_2*3];
                                Real v = -Q[d_num_species + 1][idx_cell_pivot_mom]/rho_pivot +
                                    Real(2)*d_bdry_edge_isothermal_no_slip_vel[face_loc_2*3 + 1];
                                Real w = -Q[d_num_species + 2][idx_cell_pivot_mom]/rho_pivot +
                                    Real(2)*d_bdry_edge_isothermal_no_slip_vel[face_loc_2*3 + 2];
                                
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Q[si][idx_cell_rho_Y] = rho*Y_pivot[si];
                                }
                                Q[d_num_species][idx_cell_mom] = rho*u;
                                Q[d_num_species + 1][idx_cell_mom] = rho*v;
                                Q[d_num_species + 2][idx_cell_mom] = rho*w;
                                
                                Real epsilon = d_equation_of_state_mixing_rules->
                                    getInternalEnergyFromTemperature(
                                        &rho,
                                        &T,
                                        Y_pivot_ptr);
                                
                                Real E = rho*epsilon +
                                    Real(1)/Real(2)*(Q[d_num_species][idx_cell_mom]*Q[d_num_species][idx_cell_mom] +
                                    Q[d_num_species + 1][idx_cell_mom]*Q[d_num_species + 1][idx_cell_mom] +
                                    Q[d_num_species + 2][idx_cell_mom]*Q[d_num_species + 2][idx_cell_mom])/rho;
                                
                                Q[d_num_species + 3][idx_cell_E] = E;
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
            TBOX_ERROR("FlowModelBoundaryUtilitiesFourEqnConservative::fill3dEdgeBoundaryData()\n"
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
FlowModelBoundaryUtilitiesFourEqnConservative::fill3dNodeBoundaryData(
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
                
                if (bdry_node_conds[node_loc] == BDRY_COND::FLOW_MODEL::XADIABATIC_NO_SLIP)
                {
                    for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                    {
                        for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                        {
                            for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                            {
                                const int idx_cell_rho_Y = (i + num_subghosts_conservative_var[0][0]) +
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
                                
                                int idx_cell_pivot_rho_Y = idx_cell_rho_Y;
                                int idx_cell_pivot_mom = idx_cell_mom;
                                int idx_cell_pivot_E = idx_cell_E;
                                
                                if (face_loc_0 == BDRY_LOC::XLO)
                                {
                                    idx_cell_pivot_rho_Y = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
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
                                    idx_cell_pivot_rho_Y = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
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
                                 * Compute the mixture density of the pivot.
                                 */
                                
                                Real rho_pivot = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_pivot += Q[si][idx_cell_pivot_rho_Y];
                                }
                                
                                /*
                                 * Compute the mass fractions of the pivot.
                                 */
                                
                                std::vector<Real> Y_pivot;
                                Y_pivot.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_pivot.push_back(Q[si][idx_cell_pivot_rho_Y]/rho_pivot);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions of the pivot.
                                 */
                                
                                std::vector<const Real*> Y_pivot_ptr;
                                Y_pivot_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_pivot_ptr.push_back(&Y_pivot[si]);
                                }
                                
                                /*
                                 * Set the values for partial densities and momentum.
                                 */
                                
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Q[si][idx_cell_rho_Y] = Q[si][idx_cell_pivot_rho_Y];
                                }
                                Q[d_num_species][idx_cell_mom] = -Q[d_num_species][idx_cell_pivot_mom] +
                                    Real(2)*rho_pivot*d_bdry_face_adiabatic_no_slip_vel[face_loc_0*3];
                                Q[d_num_species + 1][idx_cell_mom] = -Q[d_num_species + 1][idx_cell_pivot_mom] +
                                    Real(2)*rho_pivot*d_bdry_face_adiabatic_no_slip_vel[face_loc_0*3 + 1];
                                Q[d_num_species + 2][idx_cell_mom] = -Q[d_num_species + 2][idx_cell_pivot_mom] +
                                    Real(2)*rho_pivot*d_bdry_face_adiabatic_no_slip_vel[face_loc_0*3 + 2];
                                
                                /*
                                 * Set the values for total internal energy.
                                 */
                                
                                Real epsilon_pivot = (Q[d_num_species + 3][idx_cell_pivot_E] -
                                    Real(1)/Real(2)*(Q[d_num_species][idx_cell_pivot_mom]*Q[d_num_species][idx_cell_pivot_mom] +
                                    Q[d_num_species + 1][idx_cell_pivot_mom]*Q[d_num_species + 1][idx_cell_pivot_mom] +
                                    Q[d_num_species + 2][idx_cell_pivot_mom]*Q[d_num_species + 2][idx_cell_pivot_mom])/
                                    rho_pivot)/rho_pivot;
                                
                                Real p_pivot = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_pivot,
                                        &epsilon_pivot,
                                        Y_pivot_ptr);
                                
                                Real T_pivot = d_equation_of_state_mixing_rules->
                                    getTemperature(
                                        &rho_pivot,
                                        &p_pivot,
                                        Y_pivot_ptr);
                                
                                Real T = T_pivot;
                                
                                Real epsilon = d_equation_of_state_mixing_rules->
                                    getInternalEnergyFromTemperature(
                                        &rho_pivot,
                                        &T,
                                        Y_pivot_ptr);
                                
                                Real E = rho_pivot*epsilon +
                                    Real(1)/Real(2)*(Q[d_num_species][idx_cell_mom]*Q[d_num_species][idx_cell_mom] +
                                    Q[d_num_species + 1][idx_cell_mom]*Q[d_num_species + 1][idx_cell_mom] +
                                    Q[d_num_species + 2][idx_cell_mom]*Q[d_num_species + 2][idx_cell_mom])/
                                    rho_pivot;
                                
                                Q[d_num_species + 3][idx_cell_E] = E;
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
                                const int idx_cell_rho_Y = (i + num_subghosts_conservative_var[0][0]) +
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
                                
                                int idx_cell_pivot_rho_Y = idx_cell_rho_Y;
                                int idx_cell_pivot_mom = idx_cell_mom;
                                int idx_cell_pivot_E = idx_cell_E;
                                
                                if (face_loc_1 == BDRY_LOC::YLO)
                                {
                                    idx_cell_pivot_rho_Y = (i + num_subghosts_conservative_var[0][0]) +
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
                                    idx_cell_pivot_rho_Y = (i + num_subghosts_conservative_var[0][0]) +
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
                                 * Compute the mixture density of the pivot.
                                 */
                                
                                Real rho_pivot = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_pivot += Q[si][idx_cell_pivot_rho_Y];
                                }
                                
                                /*
                                 * Compute the mass fractions of the pivot.
                                 */
                                
                                std::vector<Real> Y_pivot;
                                Y_pivot.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_pivot.push_back(Q[si][idx_cell_pivot_rho_Y]/rho_pivot);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions of the pivot.
                                 */
                                
                                std::vector<const Real*> Y_pivot_ptr;
                                Y_pivot_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_pivot_ptr.push_back(&Y_pivot[si]);
                                }
                                
                                /*
                                 * Set the values for partial densities and momentum.
                                 */
                                
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Q[si][idx_cell_rho_Y] = Q[si][idx_cell_pivot_rho_Y];
                                }
                                Q[d_num_species][idx_cell_mom] = -Q[d_num_species][idx_cell_pivot_mom] +
                                    Real(2)*rho_pivot*d_bdry_face_adiabatic_no_slip_vel[face_loc_1*3];
                                Q[d_num_species + 1][idx_cell_mom] = -Q[d_num_species + 1][idx_cell_pivot_mom] +
                                    Real(2)*rho_pivot*d_bdry_face_adiabatic_no_slip_vel[face_loc_1*3 + 1];
                                Q[d_num_species + 2][idx_cell_mom] = -Q[d_num_species + 2][idx_cell_pivot_mom] +
                                    Real(2)*rho_pivot*d_bdry_face_adiabatic_no_slip_vel[face_loc_1*3 + 2];
                                
                                /*
                                 * Set the values for total internal energy.
                                 */
                                
                                Real epsilon_pivot = (Q[d_num_species + 3][idx_cell_pivot_E] -
                                    Real(1)/Real(2)*(Q[d_num_species][idx_cell_pivot_mom]*Q[d_num_species][idx_cell_pivot_mom] +
                                    Q[d_num_species + 1][idx_cell_pivot_mom]*Q[d_num_species + 1][idx_cell_pivot_mom] +
                                    Q[d_num_species + 2][idx_cell_pivot_mom]*Q[d_num_species + 2][idx_cell_pivot_mom])/
                                    rho_pivot)/rho_pivot;
                                
                                Real p_pivot = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_pivot,
                                        &epsilon_pivot,
                                        Y_pivot_ptr);
                                
                                Real T_pivot = d_equation_of_state_mixing_rules->
                                    getTemperature(
                                        &rho_pivot,
                                        &p_pivot,
                                        Y_pivot_ptr);
                                
                                Real T = T_pivot;
                                
                                Real epsilon = d_equation_of_state_mixing_rules->
                                    getInternalEnergyFromTemperature(
                                        &rho_pivot,
                                        &T,
                                        Y_pivot_ptr);
                                
                                Real E = rho_pivot*epsilon +
                                    Real(1)/Real(2)*(Q[d_num_species][idx_cell_mom]*Q[d_num_species][idx_cell_mom] +
                                    Q[d_num_species + 1][idx_cell_mom]*Q[d_num_species + 1][idx_cell_mom] +
                                    Q[d_num_species + 2][idx_cell_mom]*Q[d_num_species + 2][idx_cell_mom])/
                                    rho_pivot;
                                
                                Q[d_num_species + 3][idx_cell_E] = E;
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
                                const int idx_cell_rho_Y = (i + num_subghosts_conservative_var[0][0]) +
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
                                
                                int idx_cell_pivot_rho_Y = idx_cell_rho_Y;
                                int idx_cell_pivot_mom = idx_cell_mom;
                                int idx_cell_pivot_E = idx_cell_E;
                                
                                if (face_loc_2 == BDRY_LOC::ZLO)
                                {
                                    idx_cell_pivot_rho_Y = (i + num_subghosts_conservative_var[0][0]) +
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
                                    idx_cell_pivot_rho_Y = (i + num_subghosts_conservative_var[0][0]) +
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
                                 * Compute the mixture density of the pivot.
                                 */
                                
                                Real rho_pivot = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_pivot += Q[si][idx_cell_pivot_rho_Y];
                                }
                                
                                /*
                                 * Compute the mass fractions of the pivot.
                                 */
                                
                                std::vector<Real> Y_pivot;
                                Y_pivot.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_pivot.push_back(Q[si][idx_cell_pivot_rho_Y]/rho_pivot);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions of the pivot.
                                 */
                                
                                std::vector<const Real*> Y_pivot_ptr;
                                Y_pivot_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_pivot_ptr.push_back(&Y_pivot[si]);
                                }
                                
                                /*
                                 * Set the values for partial densities and momentum.
                                 */
                                
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Q[si][idx_cell_rho_Y] = Q[si][idx_cell_pivot_rho_Y];
                                }
                                Q[d_num_species][idx_cell_mom] = -Q[d_num_species][idx_cell_pivot_mom] +
                                    Real(2)*rho_pivot*d_bdry_face_adiabatic_no_slip_vel[face_loc_2*3];
                                Q[d_num_species + 1][idx_cell_mom] = -Q[d_num_species + 1][idx_cell_pivot_mom] +
                                    Real(2)*rho_pivot*d_bdry_face_adiabatic_no_slip_vel[face_loc_2*3 + 1];
                                Q[d_num_species + 2][idx_cell_mom] = -Q[d_num_species + 2][idx_cell_pivot_mom] +
                                    Real(2)*rho_pivot*d_bdry_face_adiabatic_no_slip_vel[face_loc_2*3 + 2];
                                
                                /*
                                 * Set the values for total internal energy.
                                 */
                                
                                Real epsilon_pivot = (Q[d_num_species + 3][idx_cell_pivot_E] -
                                    Real(1)/Real(2)*(Q[d_num_species][idx_cell_pivot_mom]*Q[d_num_species][idx_cell_pivot_mom] +
                                    Q[d_num_species + 1][idx_cell_pivot_mom]*Q[d_num_species + 1][idx_cell_pivot_mom] +
                                    Q[d_num_species + 2][idx_cell_pivot_mom]*Q[d_num_species + 2][idx_cell_pivot_mom])/
                                    rho_pivot)/rho_pivot;
                                
                                Real p_pivot = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_pivot,
                                        &epsilon_pivot,
                                        Y_pivot_ptr);
                                
                                Real T_pivot = d_equation_of_state_mixing_rules->
                                    getTemperature(
                                        &rho_pivot,
                                        &p_pivot,
                                        Y_pivot_ptr);
                                
                                Real T = T_pivot;
                                
                                Real epsilon = d_equation_of_state_mixing_rules->
                                    getInternalEnergyFromTemperature(
                                        &rho_pivot,
                                        &T,
                                        Y_pivot_ptr);
                                
                                Real E = rho_pivot*epsilon +
                                    Real(1)/Real(2)*(Q[d_num_species][idx_cell_mom]*Q[d_num_species][idx_cell_mom] +
                                    Q[d_num_species + 1][idx_cell_mom]*Q[d_num_species + 1][idx_cell_mom] +
                                    Q[d_num_species + 2][idx_cell_mom]*Q[d_num_species + 2][idx_cell_mom])/
                                    rho_pivot;
                                
                                Q[d_num_species + 3][idx_cell_E] = E;
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
                                const int idx_cell_rho_Y = (i + num_subghosts_conservative_var[0][0]) +
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
                                
                                int idx_cell_pivot_rho_Y = idx_cell_rho_Y;
                                int idx_cell_pivot_mom = idx_cell_mom;
                                int idx_cell_pivot_E = idx_cell_E;
                                
                                if (face_loc_0 == BDRY_LOC::XLO)
                                {
                                    idx_cell_pivot_rho_Y = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
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
                                    idx_cell_pivot_rho_Y = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
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
                                 * Compute the mixture density of the pivot.
                                 */
                                
                                Real rho_pivot = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_pivot += Q[si][idx_cell_pivot_rho_Y];
                                }
                                
                                /*
                                 * Compute the mass fractions of the pivot.
                                 */
                                
                                std::vector<Real> Y_pivot;
                                Y_pivot.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_pivot.push_back(Q[si][idx_cell_pivot_rho_Y]/rho_pivot);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions of the pivot.
                                 */
                                
                                std::vector<const Real*> Y_pivot_ptr;
                                Y_pivot_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_pivot_ptr.push_back(&Y_pivot[si]);
                                }
                                
                                /*
                                 * Set the values for partial densities, momentum and total internal energy.
                                 */
                                
                                Real epsilon_pivot = (Q[d_num_species + 3][idx_cell_pivot_E] -
                                    Real(1)/Real(2)*(Q[d_num_species][idx_cell_pivot_mom]*Q[d_num_species][idx_cell_pivot_mom] +
                                    Q[d_num_species + 1][idx_cell_pivot_mom]*Q[d_num_species + 1][idx_cell_pivot_mom] +
                                    Q[d_num_species + 2][idx_cell_pivot_mom]*Q[d_num_species + 2][idx_cell_pivot_mom])/
                                    rho_pivot)/rho_pivot;
                                
                                Real p_pivot = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_pivot,
                                        &epsilon_pivot,
                                        Y_pivot_ptr);
                                
                                Real p = p_pivot;
                                
                                Real T_pivot = d_equation_of_state_mixing_rules->
                                    getTemperature(
                                        &rho_pivot,
                                        &p_pivot,
                                        Y_pivot_ptr);
                                
                                Real T = -T_pivot + Real(2)*d_bdry_face_isothermal_no_slip_T[face_loc_0];
                                
                                Real rho = d_equation_of_state_mixing_rules->
                                    getMixtureDensity(
                                        &p,
                                        &T,
                                        Y_pivot_ptr);
                                
                                Real u = -Q[d_num_species][idx_cell_pivot_mom]/rho_pivot +
                                    Real(2)*d_bdry_edge_isothermal_no_slip_vel[face_loc_0*3];
                                Real v = -Q[d_num_species + 1][idx_cell_pivot_mom]/rho_pivot +
                                    Real(2)*d_bdry_edge_isothermal_no_slip_vel[face_loc_0*3 + 1];
                                Real w = -Q[d_num_species + 2][idx_cell_pivot_mom]/rho_pivot +
                                    Real(2)*d_bdry_edge_isothermal_no_slip_vel[face_loc_0*3 + 2];
                                
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Q[si][idx_cell_rho_Y] = rho*Y_pivot[si];
                                }
                                Q[d_num_species][idx_cell_mom] = rho*u;
                                Q[d_num_species + 1][idx_cell_mom] = rho*v;
                                Q[d_num_species + 2][idx_cell_mom] = rho*w;
                                
                                Real epsilon = d_equation_of_state_mixing_rules->
                                    getInternalEnergyFromTemperature(
                                        &rho,
                                        &T,
                                        Y_pivot_ptr);
                                
                                Real E = rho*epsilon +
                                    Real(1)/Real(2)*(Q[d_num_species][idx_cell_mom]*Q[d_num_species][idx_cell_mom] +
                                    Q[d_num_species + 1][idx_cell_mom]*Q[d_num_species + 1][idx_cell_mom] +
                                    Q[d_num_species + 2][idx_cell_mom]*Q[d_num_species + 2][idx_cell_mom])/rho;
                                
                                Q[d_num_species + 3][idx_cell_E] = E;
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
                                const int idx_cell_rho_Y = (i + num_subghosts_conservative_var[0][0]) +
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
                                
                                int idx_cell_pivot_rho_Y = idx_cell_rho_Y;
                                int idx_cell_pivot_mom = idx_cell_mom;
                                int idx_cell_pivot_E = idx_cell_E;
                                
                                if (face_loc_1 == BDRY_LOC::YLO)
                                {
                                    idx_cell_pivot_rho_Y = (i + num_subghosts_conservative_var[0][0]) +
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
                                    idx_cell_pivot_rho_Y = (i + num_subghosts_conservative_var[0][0]) +
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
                                 * Compute the mixture density of the pivot.
                                 */
                                
                                Real rho_pivot = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_pivot += Q[si][idx_cell_pivot_rho_Y];
                                }
                                
                                /*
                                 * Compute the mass fractions of the pivot.
                                 */
                                
                                std::vector<Real> Y_pivot;
                                Y_pivot.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_pivot.push_back(Q[si][idx_cell_pivot_rho_Y]/rho_pivot);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions of the pivot.
                                 */
                                
                                std::vector<const Real*> Y_pivot_ptr;
                                Y_pivot_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_pivot_ptr.push_back(&Y_pivot[si]);
                                }
                                
                                /*
                                 * Set the values for partial densities, momentum and total internal energy.
                                 */
                                
                                Real epsilon_pivot = (Q[d_num_species + 3][idx_cell_pivot_E] -
                                    Real(1)/Real(2)*(Q[d_num_species][idx_cell_pivot_mom]*Q[d_num_species][idx_cell_pivot_mom] +
                                    Q[d_num_species + 1][idx_cell_pivot_mom]*Q[d_num_species + 1][idx_cell_pivot_mom] +
                                    Q[d_num_species + 2][idx_cell_pivot_mom]*Q[d_num_species + 2][idx_cell_pivot_mom])/
                                    rho_pivot)/rho_pivot;
                                
                                Real p_pivot = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_pivot,
                                        &epsilon_pivot,
                                        Y_pivot_ptr);
                                
                                Real p = p_pivot;
                                
                                Real T_pivot = d_equation_of_state_mixing_rules->
                                    getTemperature(
                                        &rho_pivot,
                                        &p_pivot,
                                        Y_pivot_ptr);
                                
                                Real T = -T_pivot + Real(2)*d_bdry_face_isothermal_no_slip_T[face_loc_1];
                                
                                Real rho = d_equation_of_state_mixing_rules->
                                    getMixtureDensity(
                                        &p,
                                        &T,
                                        Y_pivot_ptr);
                                
                                Real u = -Q[d_num_species][idx_cell_pivot_mom]/rho_pivot +
                                    Real(2)*d_bdry_edge_isothermal_no_slip_vel[face_loc_1*3];
                                Real v = -Q[d_num_species + 1][idx_cell_pivot_mom]/rho_pivot +
                                    Real(2)*d_bdry_edge_isothermal_no_slip_vel[face_loc_1*3 + 1];
                                Real w = -Q[d_num_species + 2][idx_cell_pivot_mom]/rho_pivot +
                                    Real(2)*d_bdry_edge_isothermal_no_slip_vel[face_loc_1*3 + 2];
                                
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Q[si][idx_cell_rho_Y] = rho*Y_pivot[si];
                                }
                                Q[d_num_species][idx_cell_mom] = rho*u;
                                Q[d_num_species + 1][idx_cell_mom] = rho*v;
                                Q[d_num_species + 2][idx_cell_mom] = rho*w;
                                
                                Real epsilon = d_equation_of_state_mixing_rules->
                                    getInternalEnergyFromTemperature(
                                        &rho,
                                        &T,
                                        Y_pivot_ptr);
                                
                                Real E = rho*epsilon +
                                    Real(1)/Real(2)*(Q[d_num_species][idx_cell_mom]*Q[d_num_species][idx_cell_mom] +
                                    Q[d_num_species + 1][idx_cell_mom]*Q[d_num_species + 1][idx_cell_mom] +
                                    Q[d_num_species + 2][idx_cell_mom]*Q[d_num_species + 2][idx_cell_mom])/rho;
                                
                                Q[d_num_species + 3][idx_cell_E] = E;
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
                                const int idx_cell_rho_Y = (i + num_subghosts_conservative_var[0][0]) +
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
                                
                                int idx_cell_pivot_rho_Y = idx_cell_rho_Y;
                                int idx_cell_pivot_mom = idx_cell_mom;
                                int idx_cell_pivot_E = idx_cell_E;
                                
                                if (face_loc_2 == BDRY_LOC::ZLO)
                                {
                                    idx_cell_pivot_rho_Y = (i + num_subghosts_conservative_var[0][0]) +
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
                                    idx_cell_pivot_rho_Y = (i + num_subghosts_conservative_var[0][0]) +
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
                                 * Compute the mixture density of the pivot.
                                 */
                                
                                Real rho_pivot = Real(0);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    rho_pivot += Q[si][idx_cell_pivot_rho_Y];
                                }
                                
                                /*
                                 * Compute the mass fractions of the pivot.
                                 */
                                
                                std::vector<Real> Y_pivot;
                                Y_pivot.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_pivot.push_back(Q[si][idx_cell_pivot_rho_Y]/rho_pivot);
                                }
                                
                                /*
                                 * Get the pointers to the mass fractions of the pivot.
                                 */
                                
                                std::vector<const Real*> Y_pivot_ptr;
                                Y_pivot_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_pivot_ptr.push_back(&Y_pivot[si]);
                                }
                                
                                /*
                                 * Set the values for partial densities, momentum and total internal energy.
                                 */
                                
                                Real epsilon_pivot = (Q[d_num_species + 3][idx_cell_pivot_E] -
                                    Real(1)/Real(2)*(Q[d_num_species][idx_cell_pivot_mom]*Q[d_num_species][idx_cell_pivot_mom] +
                                    Q[d_num_species + 1][idx_cell_pivot_mom]*Q[d_num_species + 1][idx_cell_pivot_mom] +
                                    Q[d_num_species + 2][idx_cell_pivot_mom]*Q[d_num_species + 2][idx_cell_pivot_mom])/
                                    rho_pivot)/rho_pivot;
                                
                                Real p_pivot = d_equation_of_state_mixing_rules->
                                    getPressure(
                                        &rho_pivot,
                                        &epsilon_pivot,
                                        Y_pivot_ptr);
                                
                                Real p = p_pivot;
                                
                                Real T_pivot = d_equation_of_state_mixing_rules->
                                    getTemperature(
                                        &rho_pivot,
                                        &p_pivot,
                                        Y_pivot_ptr);
                                
                                Real T = -T_pivot + Real(2)*d_bdry_face_isothermal_no_slip_T[face_loc_2];
                                
                                Real rho = d_equation_of_state_mixing_rules->
                                    getMixtureDensity(
                                        &p,
                                        &T,
                                        Y_pivot_ptr);
                                
                                Real u = -Q[d_num_species][idx_cell_pivot_mom]/rho_pivot +
                                    Real(2)*d_bdry_edge_isothermal_no_slip_vel[face_loc_2*3];
                                Real v = -Q[d_num_species + 1][idx_cell_pivot_mom]/rho_pivot +
                                    Real(2)*d_bdry_edge_isothermal_no_slip_vel[face_loc_2*3 + 1];
                                Real w = -Q[d_num_species + 2][idx_cell_pivot_mom]/rho_pivot +
                                    Real(2)*d_bdry_edge_isothermal_no_slip_vel[face_loc_2*3 + 2];
                                
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Q[si][idx_cell_rho_Y] = rho*Y_pivot[si];
                                }
                                Q[d_num_species][idx_cell_mom] = rho*u;
                                Q[d_num_species + 1][idx_cell_mom] = rho*v;
                                Q[d_num_species + 2][idx_cell_mom] = rho*w;
                                
                                Real epsilon = d_equation_of_state_mixing_rules->
                                    getInternalEnergyFromTemperature(
                                        &rho,
                                        &T,
                                        Y_pivot_ptr);
                                
                                Real E = rho*epsilon +
                                    Real(1)/Real(2)*(Q[d_num_species][idx_cell_mom]*Q[d_num_species][idx_cell_mom] +
                                    Q[d_num_species + 1][idx_cell_mom]*Q[d_num_species + 1][idx_cell_mom] +
                                    Q[d_num_species + 2][idx_cell_mom]*Q[d_num_species + 2][idx_cell_mom])/rho;
                                
                                Q[d_num_species + 3][idx_cell_E] = E;
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
            TBOX_ERROR("FlowModelBoundaryUtilitiesFourEqnConservative::fill3dNodeBoundaryData\n"
                << "Invalid node boundary condition!\n"
                << "node_loc = '" << node_loc << "'." << std::endl
                << "bdry_node_conds[node_loc] = '" << bdry_node_conds[node_loc] << "'."
                << std::endl);
        }
    }
}


void
FlowModelBoundaryUtilitiesFourEqnConservative::read1dBdryNodes(
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
                    TBOX_ERROR("FlowModelBoundaryUtilitiesFourEqnConservative::read1dBdryNodes()\n"
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
FlowModelBoundaryUtilitiesFourEqnConservative::read2dBdryEdges(
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
                    TBOX_ERROR("FlowModelBoundaryUtilitiesFourEqnConservative::read2dBdryEdges()\n"
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
FlowModelBoundaryUtilitiesFourEqnConservative::read2dBdryNodes(
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
                TBOX_ERROR("FlowModelBoundaryUtilitiesFourEqnConservative::read2dBdryNodes()\n"
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
                    << ": FlowModelBoundaryUtilitiesFourEqnConservative::read2dBdryNodes()\n"
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
FlowModelBoundaryUtilitiesFourEqnConservative::read3dBdryFaces(
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
                    TBOX_ERROR("FlowModelBoundaryUtilitiesFourEqnConservative::read3dBdryFaces\n"
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
FlowModelBoundaryUtilitiesFourEqnConservative::read3dBdryEdges(
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
                    TBOX_ERROR("FlowModelBoundaryUtilitiesFourEqnConservative::read3dBdryEdges\n"
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
                        << ": FlowModelBoundaryUtilitiesFourEqnConservative::read3dBdryEdges()\n"
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
                        << ": FlowModelBoundaryUtilitiesFourEqnConservative::read3dBdryEdges()\n"
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
FlowModelBoundaryUtilitiesFourEqnConservative::read3dBdryNodes(
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
                TBOX_ERROR("FlowModelBoundaryUtilitiesFourEqnConservative::read3dBdryNodes()\n"
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
                    << ": FlowModelBoundaryUtilitiesFourEqnConservative::read3dBdryNodes()\n"
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
FlowModelBoundaryUtilitiesFourEqnConservative::readAdiabaticNoSlip(
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
            << ": FlowModelBoundaryUtilitiesFourEqnConservative::readAdiabaticNoSlip()\n"
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
                << ": FlowModelBoundaryUtilitiesFourEqnConservative::readAdiabaticNoSlip()\n"
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
                << ": FlowModelBoundaryUtilitiesFourEqnConservative::readAdiabaticNoSlip()\n"
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
                << ": FlowModelBoundaryUtilitiesFourEqnConservative::readAdiabaticNoSlip()\n"
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
FlowModelBoundaryUtilitiesFourEqnConservative::readIsothermalNoSlip(
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
            << ": FlowModelBoundaryUtilitiesFourEqnConservative::readIsothermalNoSlip()\n"
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
            << ": FlowModelBoundaryUtilitiesFourEqnConservative::readIsothermalNoSlip()\n"
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
                << ": FlowModelBoundaryUtilitiesFourEqnConservative::readIsothermalNoSlip()\n"
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
                << ": FlowModelBoundaryUtilitiesFourEqnConservative::readIsothermalNoSlip()\n"
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
                << ": FlowModelBoundaryUtilitiesFourEqnConservative::readIsothermalNoSlip()\n"
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
FlowModelBoundaryUtilitiesFourEqnConservative::readNonreflectingOutflow(
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
            << ": FlowModelBoundaryUtilitiesFourEqnConservative::readNonreflectingOutflow()\n"
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
            << ": FlowModelBoundaryUtilitiesFourEqnConservative::readNonreflectingOutflow()\n"
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
            << ": FlowModelBoundaryUtilitiesFourEqnConservative::readNonreflectingOutflow()\n"
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
