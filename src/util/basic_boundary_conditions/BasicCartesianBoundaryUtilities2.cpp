/*************************************************************************
 *
 * This file is modified from CartesianBoundaryUtilities2.c of the SAMRAI
 * distribution. For full copyright information, see COPYRIGHT and
 * COPYING.LESSER of SAMRAI distribution.
 *
 * Copyright:     (c) 1997-2014 Lawrence Livermore National Security, LLC
 * Description:   Utility routines for manipulating 2D Cartesian boundary data
 *
 ************************************************************************/

#include "util/basic_boundary_conditions/BasicCartesianBoundaryUtilities2.hpp"

#include "util/basic_boundary_conditions/BasicBoundaryConditions.hpp"
#include "util/basic_boundary_conditions/CartesianBoundaryDefines.hpp"

#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"

#include <algorithm>

/*
 * This function reads 2D boundary data from given input database.
 * The integer boundary condition types are placed in the integer
 * arrays supplied by the caller (typically, the concrete BoundaryStrategy
 * provided).  When DIRICHLET or NEUMANN conditions are specified, control
 * is passed to the BoundaryStrategy to read the boundary state data
 * specific to the problem.
 *
 * Errors will be reported and the program will abort whenever necessary
 * boundary condition information is missing in the input database, or
 * when the data read in is either unknown or inconsistent.  The periodic
 * domain information is used to determine which boundary edges or
 * node entries are not required from input.  Error checking requires
 * that node boundary conditions are consistent with those
 * specified for the edges.
 *
 * Arguments are:
 *    bdry_strategy .... object that reads DIRICHLET or NEUMANN data
 *    input_db ......... input database containing all boundary data
 *    edge_locs ........ array of locations of edges for applying
 *                       boundary conditions.
 *    node_locs ........ array of locations of nodes for applying
 *                       boundary conditions.
 *    edge_conds ....... array into which integer boundary conditions
 *                       for edges are read
 *    node_conds ....... array into which integer boundary conditions
 *                       for nodes are read
 *    periodic ......... integer vector specifying which coordinate
 *                       directions are periodic (value returned from
 *                       GridGeometry2::getPeriodicShift())
 */
void
BasicCartesianBoundaryUtilities2::getFromInput(
    BoundaryUtilityStrategy* bdry_strategy,
    const HAMERS_SHARED_PTR<tbox::Database>& input_db,
    const std::vector<int>& edge_locs,
    const std::vector<int>& node_locs,
    std::vector<int>& edge_conds,
    std::vector<int>& node_conds,
    const hier::IntVector& periodic)
{
    TBOX_DIM_ASSERT(periodic.getDim() == tbox::Dimension(2));
    TBOX_ASSERT(bdry_strategy != 0);
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
        TBOX_ERROR("BasicCartesianBoundaryUtilities2::getFromInput()\n"
            << "No input database supplied."
            << std::endl);
    }
    
    if (static_cast<int>(edge_locs.size()) > 0)
    {
        read2dBdryEdges(
            bdry_strategy,
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
 * Function to fill edge boundary values.
 *
 * Arguments are:
 *    var_name ............. name of variable (for error reporting)
 *    var_data ............. cell-centered patch data object to check
 *    patch ................ patch on which data object lives
 *    bdry_edge_locs ....... array of locations of edges for applying
 *                           boundary conditions.
 *    bdry_edge_conds ...... array of boundary conditions for patch edges
 *    bdry_edge_values ..... array of boundary values for edges
 *                           (this must be consistent with boundary
 *                           condition types)
 *    ghost_width_to_fill .. width of ghost region to fill
 */
void
BasicCartesianBoundaryUtilities2::fillEdgeBoundaryData(
    const std::string& var_name,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& var_data,
    const hier::Patch& patch,
    const std::vector<int>& bdry_edge_locs,
    const std::vector<int>& bdry_edge_conds,
    const std::vector<double>& bdry_edge_values,
    const hier::IntVector& ghost_width_to_fill)
{
    TBOX_ASSERT(!var_name.empty());
    TBOX_ASSERT(var_data);
    TBOX_ASSERT(static_cast<int>(bdry_edge_locs.size()) <= NUM_2D_EDGES);
    TBOX_ASSERT(*min_element(bdry_edge_locs.begin(), bdry_edge_locs.end()) >= 0);
    TBOX_ASSERT(*max_element(bdry_edge_locs.begin(), bdry_edge_locs.end()) < NUM_2D_EDGES);
    TBOX_ASSERT(static_cast<int>(bdry_edge_conds.size()) == NUM_2D_EDGES);
    TBOX_ASSERT(static_cast<int>(bdry_edge_values.size()) == NUM_2D_EDGES*(var_data->getDepth()));
    
    TBOX_DIM_ASSERT(ghost_width_to_fill.getDim() == tbox::Dimension(2));
    TBOX_ASSERT_OBJDIM_EQUALITY3(*var_data, patch, ghost_width_to_fill);
    
    NULL_USE(var_name);
    
    const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
        HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    TBOX_ASSERT(patch_geom);
    
    const hier::IntVector& num_ghosts(var_data->getGhostCellWidth());
    
    /*
     * Determine the ghost cell width to fill.
     */
    hier::IntVector gcw_to_fill(tbox::Dimension(2));
    
    // If the ghost fill width is not used, it is set to the ghost cell width of the data.
    if (ghost_width_to_fill == -hier::IntVector::getOne(tbox::Dimension(2)))
    {
        gcw_to_fill = var_data->getGhostCellWidth();
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
    
    // Get the dimensions of box that covers interior of patch plus
    // ghost cells.
    const hier::IntVector ghostcell_dims = var_data->getGhostBox().numberCells();
    
    const std::vector<hier::BoundaryBox>& edge_bdry =
        patch_geom->getCodimensionBoundaries(BDRY::EDGE2D);
    
    const int var_depth = var_data->getDepth();
    
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
            
            if (bdry_edge_conds[edge_loc] == BDRY_COND::BASIC::DIRICHLET)
            {
                for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                {
                    for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                    {
                        const int idx_cell = (i + num_ghosts[0]) +
                            (j + num_ghosts[1])*ghostcell_dims[0];
                        
                        for (int di = 0; di < var_depth; di++)
                        {
                            var_data->getPointer(di)[idx_cell] = bdry_edge_values[edge_loc*var_depth + di];
                        }
                    }
                }
            }
            else if (bdry_edge_conds[edge_loc] == BDRY_COND::BASIC::NEUMANN)
            {
                // NOT YET IMPLEMENTED
            }
            else if (bdry_edge_conds[edge_loc] == BDRY_COND::BASIC::FLOW)
            {
                for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                {
                    for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                    {
                        const int idx_cell = (i + num_ghosts[0]) +
                            (j + num_ghosts[1])*ghostcell_dims[0];
                        
                        int idx_cell_pivot = idx_cell;
                        
                        if (edge_loc == BDRY_LOC::XLO)
                        {
                            idx_cell_pivot = (interior_box_lo_idx[0] + num_ghosts[0]) +
                                (j + num_ghosts[1])*ghostcell_dims[0];
                        }
                        else if (edge_loc == BDRY_LOC::XHI)
                        {
                            idx_cell_pivot = (interior_box_hi_idx[0] + num_ghosts[0]) +
                                (j + num_ghosts[1])*ghostcell_dims[0];
                        }
                        else if (edge_loc == BDRY_LOC::YLO)
                        {
                            idx_cell_pivot = (i + num_ghosts[0]) +
                                (interior_box_lo_idx[1] + num_ghosts[1])*ghostcell_dims[0];
                        }
                        else if (edge_loc == BDRY_LOC::YHI)
                        {
                            idx_cell_pivot = (i + num_ghosts[0]) +
                                (interior_box_hi_idx[1] + num_ghosts[1])*ghostcell_dims[0];
                        }
                        
                        for (int di = 0; di < var_depth; di++)
                        {
                            var_data->getPointer(di)[idx_cell] = var_data->getPointer(di)[idx_cell_pivot];
                        }
                    }
                }
            }
            else if (bdry_edge_conds[edge_loc] == BDRY_COND::BASIC::REFLECT)
            {
                for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                {
                    for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                    {
                        const int idx_cell = (i + num_ghosts[0]) +
                            (j + num_ghosts[1])*ghostcell_dims[0];
                        
                        int idx_cell_pivot = idx_cell;
                        
                        if (edge_loc == BDRY_LOC::XLO)
                        {
                            idx_cell_pivot = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) + num_ghosts[0]) +
                                (j + num_ghosts[1])*ghostcell_dims[0];
                        }
                        else if (edge_loc == BDRY_LOC::XHI)
                        {
                            idx_cell_pivot = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) + num_ghosts[0]) +
                                (j + num_ghosts[1])*ghostcell_dims[0];
                        }
                        else if (edge_loc == BDRY_LOC::YLO)
                        {
                            idx_cell_pivot = (i + num_ghosts[0]) +
                                (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) + num_ghosts[1])*ghostcell_dims[0];
                        }
                        else if (edge_loc == BDRY_LOC::YHI)
                        {
                            idx_cell_pivot = (i + num_ghosts[0]) +
                                (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) + num_ghosts[1])*ghostcell_dims[0];
                        }
                        
                        for (int di = 0; di < var_depth; di++)
                        {
                            var_data->getPointer(di)[idx_cell] = var_data->getPointer(di)[idx_cell_pivot];
                        }
                        
                        if (edge_loc == BDRY_LOC::XLO || edge_loc == BDRY_LOC::XHI)
                        {
                            var_data->getPointer(0)[idx_cell] = -var_data->getPointer(0)[idx_cell_pivot];
                        }
                        else if (edge_loc == BDRY_LOC::YLO || edge_loc == BDRY_LOC::YHI)
                        {
                            var_data->getPointer(1)[idx_cell] = -var_data->getPointer(1)[idx_cell_pivot];
                        }
                    }
                }
            }
            else if (bdry_edge_conds[edge_loc] == BDRY_COND::BASIC::SYMMETRY)
            {
                for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                {
                    for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                    {
                        const int idx_cell = (i + num_ghosts[0]) +
                            (j + num_ghosts[1])*ghostcell_dims[0];
                        
                        int idx_cell_pivot = idx_cell;
                        
                        if (edge_loc == BDRY_LOC::XLO)
                        {
                            idx_cell_pivot = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) + num_ghosts[0]) +
                                (j + num_ghosts[1])*ghostcell_dims[0];
                        }
                        else if (edge_loc == BDRY_LOC::XHI)
                        {
                            idx_cell_pivot = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) + num_ghosts[0]) +
                                (j + num_ghosts[1])*ghostcell_dims[0];
                        }
                        else if (edge_loc == BDRY_LOC::YLO)
                        {
                            idx_cell_pivot = (i + num_ghosts[0]) +
                                (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) + num_ghosts[1])*ghostcell_dims[0];
                        }
                        else if (edge_loc == BDRY_LOC::YHI)
                        {
                            idx_cell_pivot = (i + num_ghosts[0]) +
                                (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) + num_ghosts[1])*ghostcell_dims[0];
                        }
                        
                        for (int di = 0; di < var_depth; di++)
                        {
                            var_data->getPointer(di)[idx_cell] = var_data->getPointer(di)[idx_cell_pivot];
                        }
                    }
                }
            }
            else
            {
                TBOX_ERROR("BasicCartesianBoundaryUtilities2::fillEdgeBoundaryData()\n"
                    << "Invalid edge boundary condition!\n"
                    << "edge_loc = '" << edge_loc << "'." << std::endl
                    << "bdry_edge_conds[edge_loc] = '" << bdry_edge_conds[edge_loc] << "'."
                    << std::endl);
            }
        }
    }
}


/*
 * Function to fill node boundary values.
 *
 * Arguments are:
 *    var_name ............. name of variable (for error reporting)
 *    var_data ............. cell-centered patch data object to check
 *    patch ................ patch on which data object lives
 *    bdry_node_locs ....... array of locations of nodes for applying
 *                           boundary conditions.
 *    bdry_node_conds ...... array of boundary conditions for patch nodes
 *    bdry_edge_values ..... array of boundary values for edges
 *                           (this must be consistent with boundary
 *                           condition types)
 *    ghost_width_to_fill .. width of ghost region to fill
 */
void
BasicCartesianBoundaryUtilities2::fillNodeBoundaryData(
    const std::string& var_name,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& var_data,
    const hier::Patch& patch,
    const std::vector<int>& bdry_node_locs,
    const std::vector<int>& bdry_node_conds,
    const std::vector<double>& bdry_edge_values,
    const hier::IntVector& ghost_width_to_fill)
{
    TBOX_ASSERT(!var_name.empty());
    TBOX_ASSERT(var_data);
    TBOX_ASSERT(static_cast<int>(bdry_node_locs.size()) <= NUM_2D_NODES);
    TBOX_ASSERT(*min_element(bdry_node_locs.begin(), bdry_node_locs.end()) >= 0);
    TBOX_ASSERT(*max_element(bdry_node_locs.begin(), bdry_node_locs.end()) < NUM_2D_NODES);
    TBOX_ASSERT(static_cast<int>(bdry_node_conds.size()) == NUM_2D_NODES);
    TBOX_ASSERT(static_cast<int>(bdry_edge_values.size()) == NUM_2D_EDGES*(var_data->getDepth()));
    
    TBOX_DIM_ASSERT(ghost_width_to_fill.getDim() == tbox::Dimension(2));
    TBOX_ASSERT_OBJDIM_EQUALITY3(*var_data, patch, ghost_width_to_fill);
    
    NULL_USE(var_name);
    
    const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
        HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    TBOX_ASSERT(patch_geom);
    
    const hier::IntVector& num_ghosts(var_data->getGhostCellWidth());
    
    /*
     * Determine the ghost cell width to fill.
     */
    hier::IntVector gcw_to_fill(tbox::Dimension(2));
    
    // If the ghost fill width is not used, it is set to the ghost cell width of the data.
    if (ghost_width_to_fill == -hier::IntVector::getOne(tbox::Dimension(2)))
    {
        gcw_to_fill = var_data->getGhostCellWidth();
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
    
    // Get the dimensions of box that covers interior of patch plus
    // ghost cells.
    const hier::IntVector ghostcell_dims = var_data->getGhostBox().numberCells();
    
    const std::vector<hier::BoundaryBox>& node_bdry =
        patch_geom->getCodimensionBoundaries(BDRY::NODE2D);
    
    const int var_depth = var_data->getDepth();
    
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
            
            if (bdry_node_conds[node_loc] == BDRY_COND::BASIC::XDIRICHLET)
            {
                for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                {
                    for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                    {
                        const int idx_cell = (i + num_ghosts[0]) +
                            (j + num_ghosts[1])*ghostcell_dims[0];
                        
                        for (int di = 0; di < var_depth; di++)
                        {
                            var_data->getPointer(di)[idx_cell] = bdry_edge_values[edge_loc_0*var_depth + di];
                        }
                        
                    }
                }
            }
            else if (bdry_node_conds[node_loc] == BDRY_COND::BASIC::YDIRICHLET)
            {
                for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                {
                    for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                    {
                        const int idx_cell = (i + num_ghosts[0]) +
                            (j + num_ghosts[1])*ghostcell_dims[0];
                        
                        for (int di = 0; di < var_depth; di++)
                        {
                            var_data->getPointer(di)[idx_cell] = bdry_edge_values[edge_loc_1*var_depth + di];
                        }
                        
                    }
                }
            }
            else if (bdry_node_conds[node_loc] == BDRY_COND::BASIC::XNEUMANN)
            {
                // NOT YET IMPLEMENTED
            }
            else if (bdry_node_conds[node_loc] == BDRY_COND::BASIC::YNEUMANN)
            {
                // NOT YET IMPLEMENTED
            }
            else if (bdry_node_conds[node_loc] == BDRY_COND::BASIC::XFLOW)
            {
                for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                {
                    for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                    {
                        const int idx_cell = (i + num_ghosts[0]) +
                            (j + num_ghosts[1])*ghostcell_dims[0];
                        
                        int idx_cell_pivot = idx_cell;
                        
                        if (edge_loc_0 == BDRY_LOC::XLO)
                        {
                            idx_cell_pivot = (interior_box_lo_idx[0] + num_ghosts[0]) +
                                (j + num_ghosts[1])*ghostcell_dims[0];
                        }
                        else if (edge_loc_0 == BDRY_LOC::XHI)
                        {
                            idx_cell_pivot = (interior_box_hi_idx[0] + num_ghosts[0]) +
                                (j + num_ghosts[1])*ghostcell_dims[0];
                        }
                        
                        for (int di = 0; di < var_depth; di++)
                        {
                            var_data->getPointer(di)[idx_cell] = var_data->getPointer(di)[idx_cell_pivot];
                        }
                    }
                }
            }
            else if (bdry_node_conds[node_loc] == BDRY_COND::BASIC::YFLOW)
            {
                for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                {
                    for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                    {
                        const int idx_cell = (i + num_ghosts[0]) +
                            (j + num_ghosts[1])*ghostcell_dims[0];
                        
                        int idx_cell_pivot = idx_cell;
                        
                        if (edge_loc_1 == BDRY_LOC::YLO)
                        {
                            idx_cell_pivot = (i + num_ghosts[0]) +
                                (interior_box_lo_idx[1] + num_ghosts[1])*ghostcell_dims[0];
                        }
                        else if (edge_loc_1 == BDRY_LOC::YHI)
                        {
                            idx_cell_pivot = (i + num_ghosts[0]) +
                                (interior_box_hi_idx[1] + num_ghosts[1])*ghostcell_dims[0];
                        }
                        
                        for (int di = 0; di < var_depth; di++)
                        {
                            var_data->getPointer(di)[idx_cell] = var_data->getPointer(di)[idx_cell_pivot];
                        }
                    }
                }
            }
            else if (bdry_node_conds[node_loc] == BDRY_COND::BASIC::XREFLECT)
            {
                for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                {
                    for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                    {
                        const int idx_cell = (i + num_ghosts[0]) +
                            (j + num_ghosts[1])*ghostcell_dims[0];
                        
                        int idx_cell_pivot = idx_cell;
                        
                        if (edge_loc_0 == BDRY_LOC::XLO)
                        {
                            idx_cell_pivot = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) + num_ghosts[0]) +
                                (j + num_ghosts[1])*ghostcell_dims[0];
                        }
                        else if (edge_loc_0 == BDRY_LOC::XHI)
                        {
                            idx_cell_pivot = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) + num_ghosts[0]) +
                                (j + num_ghosts[1])*ghostcell_dims[0];
                        }
                        
                        for (int di = 0; di < var_depth; di++)
                        {
                            var_data->getPointer(di)[idx_cell] = var_data->getPointer(di)[idx_cell_pivot];
                        }
                        
                        var_data->getPointer(0)[idx_cell] = -var_data->getPointer(0)[idx_cell_pivot];
                    }
                }
            }
            else if (bdry_node_conds[node_loc] == BDRY_COND::BASIC::YREFLECT)
            {
                for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                {
                    for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                    {
                        const int idx_cell = (i + num_ghosts[0]) +
                            (j + num_ghosts[1])*ghostcell_dims[0];
                        
                        int idx_cell_pivot = idx_cell;
                        
                        if (edge_loc_1 == BDRY_LOC::YLO)
                        {
                            idx_cell_pivot = (i + num_ghosts[0]) +
                                (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) + num_ghosts[1])*ghostcell_dims[0];
                        }
                        else if (edge_loc_1 == BDRY_LOC::YHI)
                        {
                            idx_cell_pivot = (i + num_ghosts[0]) +
                                (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) + num_ghosts[1])*ghostcell_dims[0];
                        }
                        
                        for (int di = 0; di < var_depth; di++)
                        {
                            var_data->getPointer(di)[idx_cell] = var_data->getPointer(di)[idx_cell_pivot];
                        }
                        
                        var_data->getPointer(1)[idx_cell] = -var_data->getPointer(1)[idx_cell_pivot];
                    }
                }
            }
            else if (bdry_node_conds[node_loc] == BDRY_COND::BASIC::XSYMMETRY)
            {
                for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                {
                    for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                    {
                        const int idx_cell = (i + num_ghosts[0]) +
                            (j + num_ghosts[1])*ghostcell_dims[0];
                        
                        int idx_cell_pivot = idx_cell;
                        
                        if (edge_loc_0 == BDRY_LOC::XLO)
                        {
                            idx_cell_pivot = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) + num_ghosts[0]) +
                                (j + num_ghosts[1])*ghostcell_dims[0];
                        }
                        else if (edge_loc_0 == BDRY_LOC::XHI)
                        {
                            idx_cell_pivot = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) + num_ghosts[0]) +
                                (j + num_ghosts[1])*ghostcell_dims[0];
                        }
                        
                        for (int di = 0; di < var_depth; di++)
                        {
                            var_data->getPointer(di)[idx_cell] = var_data->getPointer(di)[idx_cell_pivot];
                        }
                    }
                }
            }
            else if (bdry_node_conds[node_loc] == BDRY_COND::BASIC::YSYMMETRY)
            {
                for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                {
                    for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                    {
                        const int idx_cell = (i + num_ghosts[0]) +
                            (j + num_ghosts[1])*ghostcell_dims[0];
                        
                        int idx_cell_pivot = idx_cell;
                        
                        if (edge_loc_1 == BDRY_LOC::YLO)
                        {
                            idx_cell_pivot = (i + num_ghosts[0]) +
                                (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) + num_ghosts[1])*ghostcell_dims[0];
                        }
                        else if (edge_loc_1 == BDRY_LOC::YHI)
                        {
                            idx_cell_pivot = (i + num_ghosts[0]) +
                                (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) + num_ghosts[1])*ghostcell_dims[0];
                        }
                        
                        for (int di = 0; di < var_depth; di++)
                        {
                            var_data->getPointer(di)[idx_cell] = var_data->getPointer(di)[idx_cell_pivot];
                        }
                    }
                }
            }
            else
            {
                TBOX_ERROR("BasicCartesianBoundaryUtilities2::fillNodeBoundaryData()\n"
                    << "Invalid node boundary condition!\n"
                    << "node_loc = '" << node_loc << "'." << std::endl
                    << "bdry_node_conds[node_loc] = '" << bdry_node_conds[node_loc] << "'."
                    << std::endl);
            }
        }
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
BasicCartesianBoundaryUtilities2::getEdgeLocationForNodeBdry(
    int node_loc,
    int node_btype)
{
    int ret_edge = -1;
    
    switch (node_btype)
    {
        case BDRY_COND::BASIC::XFLOW:
        case BDRY_COND::BASIC::XREFLECT:
        case BDRY_COND::BASIC::XSYMMETRY:
        case BDRY_COND::BASIC::XDIRICHLET:
        case BDRY_COND::BASIC::XNEUMANN:
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
        case BDRY_COND::BASIC::YFLOW:
        case BDRY_COND::BASIC::YREFLECT:
        case BDRY_COND::BASIC::YSYMMETRY:
        case BDRY_COND::BASIC::YDIRICHLET:
        case BDRY_COND::BASIC::YNEUMANN:
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
            // TBOX_ERROR("BasicCartesianBoundaryUtilities2::getEdgeLocationForNodeBdry()\n"
            //     << "Unknown node boundary condition type = '"
            //     << node_btype
            //     << "' passed."
            //     << std::endl);
        }
    }
    
    // if (ret_edge == -1)
    // {
    //     TBOX_ERROR("BasicCartesianBoundaryUtilities2::getEdgeLocationForNodeBdry()\n"
    //         << "Node boundary condition type = '"
    //         << node_btype << "' and \n"
    //         << "node location = '" << node_loc
    //         << "' passed are inconsistent."
    //         << std::endl);
    // }
    
    return ret_edge;
}


/*
 * Private function to read 2D edge boundary data from input database.
 */
void
BasicCartesianBoundaryUtilities2::read2dBdryEdges(
    BoundaryUtilityStrategy* bdry_strategy,
    const HAMERS_SHARED_PTR<tbox::Database>& input_db,
    const std::vector<int>& edge_locs,
    std::vector<int>& edge_conds,
    const hier::IntVector& periodic)
{
    TBOX_DIM_ASSERT(periodic.getDim() == tbox::Dimension(2));
    
    TBOX_ASSERT(bdry_strategy != 0);
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
        // edge boundary input required
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
                if (bdry_cond_str == "FLOW")
                {
                    edge_conds[s] = BDRY_COND::BASIC::FLOW;
                }
                else if (bdry_cond_str == "REFLECT")
                {
                    edge_conds[s] = BDRY_COND::BASIC::REFLECT;
                }
                else if (bdry_cond_str == "SYMMETRY")
                {
                    edge_conds[s] = BDRY_COND::BASIC::SYMMETRY;
                }
                else if (bdry_cond_str == "DIRICHLET")
                {
                    edge_conds[s] = BDRY_COND::BASIC::DIRICHLET;
                    bdry_strategy->readDirichletBoundaryDataEntry(
                        bdry_loc_db,
                        bdry_loc_str,
                        s);
                }
                else if (bdry_cond_str == "NEUMANN")
                {
                    edge_conds[s] = BDRY_COND::BASIC::NEUMANN;
                    bdry_strategy->readNeumannBoundaryDataEntry(
                        bdry_loc_db,
                        bdry_loc_str,
                        s);
                }
                else
                {
                    TBOX_ERROR("BasicCartesianBoundaryUtilities2::read2dBdryEdges()\n"
                        << "Unknown edge boundary string = '"
                        << bdry_cond_str
                        << "' found in input."
                        << std::endl);
                }
            } // if (need_data_read)
       } // for (int ei = 0 ...
    } // if (num_per_dirs < 2)
}


/*
 * Private function to read 2D node boundary data from input database.
 */
void
BasicCartesianBoundaryUtilities2::read2dBdryNodes(
    const HAMERS_SHARED_PTR<tbox::Database>& input_db,
    const std::vector<int>& node_locs,
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
        // node boundary data required
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
            if (bdry_cond_str == "XFLOW")
            {
                node_conds[s] = BDRY_COND::BASIC::XFLOW;
            }
            else if (bdry_cond_str == "YFLOW")
            {
                node_conds[s] = BDRY_COND::BASIC::YFLOW;
            }
            else if (bdry_cond_str == "XREFLECT")
            {
                node_conds[s] = BDRY_COND::BASIC::XREFLECT;
            }
            else if (bdry_cond_str == "YREFLECT")
            {
                node_conds[s] = BDRY_COND::BASIC::YREFLECT;
            }
            else if (bdry_cond_str == "XSYMMETRY")
            {
                node_conds[s] = BDRY_COND::BASIC::XSYMMETRY;
            }
            else if (bdry_cond_str == "YSYMMETRY")
            {
                node_conds[s] = BDRY_COND::BASIC::YSYMMETRY;
            }
            else if (bdry_cond_str == "XDIRICHLET")
            {
                node_conds[s] = BDRY_COND::BASIC::XDIRICHLET;
            }
            else if (bdry_cond_str == "YDIRICHLET")
            {
                node_conds[s] = BDRY_COND::BASIC::YDIRICHLET;
            }
            else if (bdry_cond_str == "XNEUMANN")
            {
                node_conds[s] = BDRY_COND::BASIC::XNEUMANN;
            }
            else if (bdry_cond_str == "YNEUMANN")
            {
                node_conds[s] = BDRY_COND::BASIC::YNEUMANN;
            }
            else
            {
                TBOX_ERROR("BasicCartesianBoundaryUtilities2::read2dBdryNodes()\n"
                    << "Unknown node boundary string = '"
                    << bdry_cond_str
                    << "' found in input."
                    << std::endl);
            }
            
            std::string proper_edge;
            std::string proper_edge_data;
            bool no_edge_data_found = false;
            if (bdry_cond_str == "XFLOW" ||
                bdry_cond_str == "XDIRICHLET" ||
                bdry_cond_str == "XNEUMANN" ||
                bdry_cond_str == "XREFLECT" ||
                bdry_cond_str == "XSYMMETRY")
            {
                if (s == NODE_BDRY_LOC_2D::XLO_YLO ||
                    s == NODE_BDRY_LOC_2D::XLO_YHI)
                {
                    proper_edge = "XLO";
                    if (bdry_cond_str == "XFLOW" &&
                        edge_conds[BDRY_LOC::XLO] != BDRY_COND::BASIC::FLOW)
                    {
                        no_edge_data_found = true;
                        proper_edge_data = "FLOW";
                    }
                    if (bdry_cond_str == "XDIRICHLET" &&
                        edge_conds[BDRY_LOC::XLO] != BDRY_COND::BASIC::DIRICHLET)
                    {
                        no_edge_data_found = true;
                        proper_edge_data = "DIRICHLET";
                    }
                    if (bdry_cond_str == "XNEUMANN" &&
                        edge_conds[BDRY_LOC::XLO] != BDRY_COND::BASIC::NEUMANN)
                    {
                        no_edge_data_found = true;
                        proper_edge_data = "NEUMANN";
                    }
                    if (bdry_cond_str == "XREFLECT" &&
                        edge_conds[BDRY_LOC::XLO] != BDRY_COND::BASIC::REFLECT)
                    {
                        no_edge_data_found = true;
                        proper_edge_data = "REFLECT";
                    }
                    if (bdry_cond_str == "XSYMMETRY" &&
                        edge_conds[BDRY_LOC::XLO] != BDRY_COND::BASIC::SYMMETRY)
                    {
                        no_edge_data_found = true;
                        proper_edge_data = "SYMMETRY";
                    }
                }
                else
                {
                    proper_edge = "XHI";
                    if (bdry_cond_str == "XFLOW" &&
                        edge_conds[BDRY_LOC::XHI] != BDRY_COND::BASIC::FLOW)
                    {
                        no_edge_data_found = true;
                        proper_edge_data = "FLOW";
                    }
                    if (bdry_cond_str == "XDIRICHLET" &&
                        edge_conds[BDRY_LOC::XHI] != BDRY_COND::BASIC::DIRICHLET)
                    {
                        no_edge_data_found = true;
                        proper_edge_data = "DIRICHLET";
                    }
                    if (bdry_cond_str == "XNEUMANN" &&
                        edge_conds[BDRY_LOC::XHI] != BDRY_COND::BASIC::NEUMANN)
                    {
                        no_edge_data_found = true;
                        proper_edge_data = "NEUMANN";
                    }
                    if (bdry_cond_str == "XREFLECT" &&
                        edge_conds[BDRY_LOC::XHI] != BDRY_COND::BASIC::REFLECT)
                    {
                        no_edge_data_found = true;
                        proper_edge_data = "REFLECT";
                    }
                    if (bdry_cond_str == "XSYMMETRY" &&
                        edge_conds[BDRY_LOC::XHI] != BDRY_COND::BASIC::SYMMETRY)
                    {
                        no_edge_data_found = true;
                        proper_edge_data = "SYMMETRY";
                    }
                }
            }
            else if (bdry_cond_str == "YFLOW" ||
                     bdry_cond_str == "YDIRICHLET" ||
                     bdry_cond_str == "YNEUMANN" ||
                     bdry_cond_str == "YREFLECT" ||
                     bdry_cond_str == "YSYMMETRY")
            {
                if (s == NODE_BDRY_LOC_2D::XLO_YLO ||
                    s == NODE_BDRY_LOC_2D::XHI_YLO)
                {
                    proper_edge = "YLO";
                    if (bdry_cond_str == "YFLOW" &&
                        edge_conds[BDRY_LOC::YLO] != BDRY_COND::BASIC::FLOW)
                    {
                        no_edge_data_found = true;
                        proper_edge_data = "FLOW";
                    }
                    if (bdry_cond_str == "YDIRICHLET" &&
                        edge_conds[BDRY_LOC::YLO] != BDRY_COND::BASIC::DIRICHLET)
                    {
                        no_edge_data_found = true;
                        proper_edge_data = "DIRICHLET";
                    }
                    if (bdry_cond_str == "YNEUMANN" &&
                        edge_conds[BDRY_LOC::YLO] != BDRY_COND::BASIC::NEUMANN)
                    {
                        no_edge_data_found = true;
                        proper_edge_data = "NEUMANN";
                    }
                    if (bdry_cond_str == "YREFLECT" &&
                        edge_conds[BDRY_LOC::YLO] != BDRY_COND::BASIC::REFLECT)
                    {
                        no_edge_data_found = true;
                        proper_edge_data = "REFLECT";
                    }
                    if (bdry_cond_str == "YSYMMETRY" &&
                        edge_conds[BDRY_LOC::YLO] != BDRY_COND::BASIC::SYMMETRY)
                    {
                        no_edge_data_found = true;
                        proper_edge_data = "SYMMETRY";
                    }
                }
                else
                {
                    proper_edge = "YHI";
                    if (bdry_cond_str == "YFLOW" &&
                        edge_conds[BDRY_LOC::YHI] != BDRY_COND::BASIC::FLOW)
                    {
                        no_edge_data_found = true;
                        proper_edge_data = "FLOW";
                    }
                    if (bdry_cond_str == "YDIRICHLET" &&
                        edge_conds[BDRY_LOC::YHI] != BDRY_COND::BASIC::DIRICHLET)
                    {
                        no_edge_data_found = true;
                        proper_edge_data = "DIRICHLET";
                    }
                    if (bdry_cond_str == "YNEUMANN" &&
                        edge_conds[BDRY_LOC::YHI] != BDRY_COND::BASIC::NEUMANN)
                    {
                        no_edge_data_found = true;
                        proper_edge_data = "NEUMANN";
                    }
                    if (bdry_cond_str == "YREFLECT" &&
                        edge_conds[BDRY_LOC::YHI] != BDRY_COND::BASIC::REFLECT)
                    {
                        no_edge_data_found = true;
                        proper_edge_data = "REFLECT";
                    }
                    if (bdry_cond_str == "YSYMMETRY" &&
                        edge_conds[BDRY_LOC::YHI] != BDRY_COND::BASIC::SYMMETRY)
                    {
                        no_edge_data_found = true;
                        proper_edge_data = "SYMMETRY";
                    }
                }
            }
            if (no_edge_data_found)
            {
                TBOX_ERROR("BasicCartesianBoundaryUtilities2::read2dBdryNodes()\n"
                    << "Bdry condition '"
                    << bdry_cond_str
                    << "' found for '"
                    << bdry_loc_str
                    << "' but no '"
                    << proper_edge_data
                    << "' data found for edge '"
                    << proper_edge
                    << "'."
                    << std::endl);
            }
        } // for (int ni = 0 ...
    } // if (num_per_dirs < 1)
}
