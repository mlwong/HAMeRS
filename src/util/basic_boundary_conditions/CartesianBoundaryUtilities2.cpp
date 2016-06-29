/*************************************************************************
 *
 * This file is modified from CartesianBoundaryUtilities2.c of the SAMRAI distribution.
 * For full copyright information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2014 Lawrence Livermore National Security, LLC
 * Description:   Utility routines for manipulating 2D Cartesian boundary data
 *
 ************************************************************************/

#include "util/basic_boundary_conditions/CartesianBoundaryUtilities2.hpp"

#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/tbox/Utilities.h"

#include "util/basic_boundary_conditions/CartesianBoundaryDefines.hpp"

/*
 *************************************************************************
 *
 * External declarations for FORTRAN 77 routines used in
 * boundary condition implementation.
 *
 *************************************************************************
 */


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
 *    bdry_db .......... input database containing all boundary data
 *    edge_conds ....... array into which integer boundary conditions
 *                       for edges are read
 *    node_conds ....... array into which integer boundary conditions
 *                       for nodes are read
 *    periodic ......... integer vector specifying which coordinate
 *                       directions are periodic (value returned from
 *                       GridGeometry2::getPeriodicShift())
 */
void
CartesianBoundaryUtilities2::getFromInput(
    BoundaryUtilityStrategy* bdry_strategy,
    const boost::shared_ptr<tbox::Database>& input_db,
    std::vector<int>& edge_conds,
    std::vector<int>& node_conds,
    const hier::IntVector& periodic)
{
    TBOX_DIM_ASSERT(periodic.getDim() == tbox::Dimension(2));
    TBOX_ASSERT(bdry_strategy != 0);
    TBOX_ASSERT(static_cast<int>(edge_conds.size()) == NUM_2D_EDGES);
    TBOX_ASSERT(static_cast<int>(node_conds.size()) == NUM_2D_NODES);
    
    if (!input_db)
    {
        TBOX_ERROR(": CartesianBoundaryUtility2::getFromInput()\n"
            << "no input database supplied"
            << std::endl);
    }
    
    read2dBdryEdges(
        bdry_strategy,
        input_db,
        edge_conds,
        periodic);
    
    read2dBdryNodes(
        input_db,
        edge_conds,
        node_conds,
        periodic);
}


/*
 * Function to fill edge boundary values.
 *
 * Arguments are:
 *    var_name ............. name of variable (for error reporting)
 *    var_data ............. cell-centered patch data object to check
 *    patch ................ patch on which data object lives
 *    bdry_edge_conds ...... array of boundary conditions for patch edges
 *    bdry_edge_values ..... array of boundary values for edges
 *                           (this must be consistent with boundary
 *                           condition types)
 *    ghost_width_to_fill .. width of ghost region to fill
 */
void
CartesianBoundaryUtilities2::fillEdgeBoundaryData(
    const std::string& var_name,
    const boost::shared_ptr<pdat::CellData<double> >& var_data,
    const hier::Patch& patch,
    const std::vector<int>& bdry_edge_conds,
    const std::vector<double>& bdry_edge_values,
    const hier::IntVector& ghost_fill_width)
{
    TBOX_ASSERT(!var_name.empty());
    TBOX_ASSERT(var_data);
    TBOX_ASSERT(static_cast<int>(bdry_edge_conds.size()) == NUM_2D_EDGES);
    TBOX_ASSERT(static_cast<int>(bdry_edge_values.size()) == NUM_2D_EDGES*(var_data->getDepth()));
    
    TBOX_DIM_ASSERT(ghost_fill_width.getDim() == tbox::Dimension(2));
    TBOX_ASSERT_OBJDIM_EQUALITY3(*var_data, patch, ghost_fill_width);
    
    NULL_USE(var_name);
    
    const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    TBOX_ASSERT(patch_geom);
    
    const hier::IntVector& num_ghosts(var_data->getGhostCellWidth());
    
    /*
     * Determine the ghost cell width to fill.
     */
    hier::IntVector gcw_to_fill(tbox::Dimension(2));
    
    // If the ghost fill width is not used, it is set to the ghost cell width of the data.
    if (ghost_fill_width == -hier::IntVector::getOne(tbox::Dimension(2)))
    {
        gcw_to_fill = var_data->getGhostCellWidth();
    }
    else
    {
        gcw_to_fill = hier::IntVector::min(
            num_ghosts,
            ghost_fill_width);
    }
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box& interior_box(patch.getBox());
    hier::Index interior_box_lo_idx(interior_box.lower());
    hier::Index interior_box_hi_idx(interior_box.upper());
    
    // Get the dimensions of box that covers interior of patch plus
    // ghost cells.
    const hier::IntVector ghostcell_dims = var_data->getGhostBox().numberCells();
    
    const std::vector<hier::BoundaryBox>& edge_bdry =
        patch_geom->getCodimensionBoundaries(Bdry::EDGE2D);
    
    const int var_depth = var_data->getDepth();
    
    for (int ei = 0; ei < static_cast<int>(edge_bdry.size()); ei++)
    {
        TBOX_ASSERT(edge_bdry[ei].getBoundaryType() == Bdry::EDGE2D);
        
        int edge_loc = edge_bdry[ei].getLocationIndex();
        
        hier::Box fill_box(patch_geom->getBoundaryFillBox(
            edge_bdry[ei],
            interior_box,
            gcw_to_fill));
        
        hier::Index fill_box_lo_idx(fill_box.lower());
        hier::Index fill_box_hi_idx(fill_box.upper());
        
if (patch_geom->getTouchesRegularBoundary(0, 1) &&
    patch_geom->getTouchesRegularBoundary(1, 1))
{
    tbox::plog << "interior_box_lo_idx: " << interior_box_lo_idx << std::endl;
    tbox::plog << "interior_box_hi_idx: " << interior_box_hi_idx << std::endl;
    tbox::plog << "fill_box_lo_idx: " << fill_box_lo_idx << std::endl;
    tbox::plog << "fill_box_hi_idx: " << fill_box_hi_idx << std::endl;
}
        
        /*
         * Offset the indices.
         */
        interior_box_lo_idx = interior_box_lo_idx - interior_box.lower();
        interior_box_hi_idx = interior_box_hi_idx - interior_box.lower();
        fill_box_lo_idx = fill_box_lo_idx - interior_box.lower();
        fill_box_hi_idx = fill_box_hi_idx - interior_box.lower();
        
        if (bdry_edge_conds[edge_loc] == BdryCond::Basic::DIRICHLET)
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
        else if (bdry_edge_conds[edge_loc] == BdryCond::Basic::NEUMANN)
        {
            // NOT YET IMPLEMENTED
        }
        else if (bdry_edge_conds[edge_loc] == BdryCond::Basic::FLOW)
        {
            for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
            {
                for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                {
                    const int idx_cell = (i + num_ghosts[0]) +
                        (j + num_ghosts[1])*ghostcell_dims[0];
                    
                    int idx_cell_pivot = idx_cell;
                    
                    if (edge_loc == BdryLoc::XLO)
                    {
                        idx_cell_pivot = (interior_box_lo_idx[0] + num_ghosts[0]) +
                            (j + num_ghosts[1])*ghostcell_dims[0];
                    }
                    else if (edge_loc == BdryLoc::XHI)
                    {
                        idx_cell_pivot = (interior_box_hi_idx[0] + num_ghosts[0]) +
                            (j + num_ghosts[1])*ghostcell_dims[0];
                    }
                    else if (edge_loc == BdryLoc::YLO)
                    {
                        idx_cell_pivot = (i + num_ghosts[0]) +
                            (interior_box_lo_idx[1] + num_ghosts[1])*ghostcell_dims[0];
                    }
                    else if (edge_loc == BdryLoc::YHI)
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
        else if (bdry_edge_conds[edge_loc] == BdryCond::Basic::REFLECT)
        {
            for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
            {
                for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                {
                    const int idx_cell = (i + num_ghosts[0]) +
                        (j + num_ghosts[1])*ghostcell_dims[0];
                    
                    int idx_cell_pivot = idx_cell;
                    
                    if (edge_loc == BdryLoc::XLO)
                    {
                        idx_cell_pivot = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) + num_ghosts[0]) +
                            (j + num_ghosts[1])*ghostcell_dims[0];
                    }
                    else if (edge_loc == BdryLoc::XHI)
                    {
                        idx_cell_pivot = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) + num_ghosts[0]) +
                            (j + num_ghosts[1])*ghostcell_dims[0];
                    }
                    else if (edge_loc == BdryLoc::YLO)
                    {
                        idx_cell_pivot = (i + num_ghosts[0]) +
                            (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) + num_ghosts[1])*ghostcell_dims[0];
                    }
                    else if (edge_loc == BdryLoc::YHI)
                    {
                        idx_cell_pivot = (i + num_ghosts[0]) +
                            (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) + num_ghosts[1])*ghostcell_dims[0];
                    }
                    
                    for (int di = 0; di < var_depth; di++)
                    {
                        var_data->getPointer(di)[idx_cell] = var_data->getPointer(di)[idx_cell_pivot];
                    }
                    
                    if (edge_loc == BdryLoc::XLO || edge_loc == BdryLoc::XHI)
                    {
                        var_data->getPointer(0)[idx_cell] = -var_data->getPointer(0)[idx_cell_pivot];
                    }
                    else if (edge_loc == BdryLoc::YLO || edge_loc == BdryLoc::YHI)
                    {
                        var_data->getPointer(1)[idx_cell] = -var_data->getPointer(1)[idx_cell_pivot];
                    }
                }
            }
        }
        else if (bdry_edge_conds[edge_loc] == BdryCond::Basic::SYMMETRY)
        {
            for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
            {
                for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                {
                    const int idx_cell = (i + num_ghosts[0]) +
                        (j + num_ghosts[1])*ghostcell_dims[0];
                    
                    int idx_cell_pivot = idx_cell;
                    
                    if (edge_loc == BdryLoc::XLO)
                    {
                        idx_cell_pivot = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) + num_ghosts[0]) +
                            (j + num_ghosts[1])*ghostcell_dims[0];
                    }
                    else if (edge_loc == BdryLoc::XHI)
                    {
                        idx_cell_pivot = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) + num_ghosts[0]) +
                            (j + num_ghosts[1])*ghostcell_dims[0];
                    }
                    else if (edge_loc == BdryLoc::YLO)
                    {
                        idx_cell_pivot = (i + num_ghosts[0]) +
                            (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) + num_ghosts[1])*ghostcell_dims[0];
                    }
                    else if (edge_loc == BdryLoc::YHI)
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
            TBOX_ERROR("CartesianBoundaryUtilities2: fillEdgeBoundaryData2D()\n"
                << "Invalid edge boundary condition!\n"
                << "edge_loc = " << edge_loc << std::endl
                << "bdry_edge_conds[edge_loc] = " << bdry_edge_conds[edge_loc]
                << std::endl);
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
 *    bdry_node_conds ...... array of boundary conditions for patch nodes
 *    bdry_edge_values ..... array of boundary values for edges
 *                           (this must be consistent with boundary
 *                           condition types)
 *    ghost_width_to_fill .. width of ghost region to fill
 */
void
CartesianBoundaryUtilities2::fillNodeBoundaryData(
    const std::string& var_name,
    const boost::shared_ptr<pdat::CellData<double> >& var_data,
    const hier::Patch& patch,
    const std::vector<int>& bdry_node_conds,
    const std::vector<double>& bdry_edge_values,
    const hier::IntVector& ghost_fill_width)
{
    TBOX_ASSERT(!var_name.empty());
    TBOX_ASSERT(var_data);
    TBOX_ASSERT(static_cast<int>(bdry_node_conds.size()) == NUM_2D_NODES);
    TBOX_ASSERT(static_cast<int>(bdry_edge_values.size()) == NUM_2D_EDGES*(var_data->getDepth()));
    
    TBOX_DIM_ASSERT(ghost_fill_width.getDim() == tbox::Dimension(2));
    TBOX_ASSERT_OBJDIM_EQUALITY3(*var_data, patch, ghost_fill_width);
    
    NULL_USE(var_name);
    
    const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    TBOX_ASSERT(patch_geom);
    
    const hier::IntVector& num_ghosts(var_data->getGhostCellWidth());
    
    /*
     * Determine the ghost cell width to fill.
     */
    hier::IntVector gcw_to_fill(tbox::Dimension(2));
    
    // If the ghost fill width is not used, it is set to the ghost cell width of the data.
    if (ghost_fill_width == -hier::IntVector::getOne(tbox::Dimension(2)))
    {
        gcw_to_fill = var_data->getGhostCellWidth();
    }
    else
    {
        gcw_to_fill = hier::IntVector::min(
            num_ghosts,
            ghost_fill_width);
    }
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box& interior_box(patch.getBox());
    hier::Index interior_box_lo_idx(interior_box.lower());
    hier::Index interior_box_hi_idx(interior_box.upper());
    
    // Get the dimensions of box that covers interior of patch plus
    // ghost cells.
    const hier::IntVector ghostcell_dims = var_data->getGhostBox().numberCells();
    
    const std::vector<hier::BoundaryBox>& node_bdry =
        patch_geom->getCodimensionBoundaries(Bdry::NODE2D);
    
    const int var_depth = var_data->getDepth();
    
    for (int ei = 0; ei < static_cast<int>(node_bdry.size()); ei++)
    {
        TBOX_ASSERT(node_bdry[ei].getBoundaryType() == Bdry::NODE2D);
        
        int node_loc = node_bdry[ei].getLocationIndex();
        
        hier::Box fill_box(patch_geom->getBoundaryFillBox(
            node_bdry[ei],
            interior_box,
            gcw_to_fill));
        
        hier::Index fill_box_lo_idx(fill_box.lower());
        hier::Index fill_box_hi_idx(fill_box.upper());
        
        /*
         * Offset the indices.
         */
        interior_box_lo_idx = interior_box_lo_idx - interior_box.lower();
        interior_box_hi_idx = interior_box_hi_idx - interior_box.lower();
        fill_box_lo_idx = fill_box_lo_idx - interior_box.lower();
        fill_box_hi_idx = fill_box_hi_idx - interior_box.lower();
        
        int edge_loc_0 = 0;
        int edge_loc_1 = 0;
        
        switch (node_loc)
        {
            case NodeBdyLoc2D::XLO_YLO:
            {
                edge_loc_0 = BdryLoc::XLO;
                edge_loc_1 = BdryLoc::YLO;
                
                break;
            }
            case NodeBdyLoc2D::XHI_YLO:
            {
                edge_loc_0 = BdryLoc::XHI;
                edge_loc_1 = BdryLoc::YLO;
                
                break;
            }
            case NodeBdyLoc2D::XLO_YHI:
            {
                edge_loc_0 = BdryLoc::XLO;
                edge_loc_1 = BdryLoc::YHI;
                
                break;
            }
            case NodeBdyLoc2D::XHI_YHI:
            {
                edge_loc_0 = BdryLoc::XHI;
                edge_loc_1 = BdryLoc::YHI;
                
                break;
            }
        }
        
        if (bdry_node_conds[node_loc] == BdryCond::Basic::XDIRICHLET)
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
        else if (bdry_node_conds[node_loc] == BdryCond::Basic::YDIRICHLET)
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
        else if (bdry_node_conds[node_loc] == BdryCond::Basic::XNEUMANN)
        {
            // NOT YET IMPLEMENTED
        }
        else if (bdry_node_conds[node_loc] == BdryCond::Basic::YNEUMANN)
        {
            // NOT YET IMPLEMENTED
        }
        else if (bdry_node_conds[node_loc] == BdryCond::Basic::XFLOW)
        {
            for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
            {
                for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                {
                    const int idx_cell = (i + num_ghosts[0]) +
                        (j + num_ghosts[1])*ghostcell_dims[0];
                    
                    int idx_cell_pivot = idx_cell;
                    
                    if (edge_loc_0 == BdryLoc::XLO)
                    {
                        idx_cell_pivot = (interior_box_lo_idx[0] + num_ghosts[0]) +
                            (j + num_ghosts[1])*ghostcell_dims[0];
                    }
                    else if (edge_loc_0 == BdryLoc::XHI)
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
        else if (bdry_node_conds[node_loc] == BdryCond::Basic::YFLOW)
        {
            for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
            {
                for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                {
                    const int idx_cell = (i + num_ghosts[0]) +
                        (j + num_ghosts[1])*ghostcell_dims[0];
                    
                    int idx_cell_pivot = idx_cell;
                    
                    if (edge_loc_1 == BdryLoc::YLO)
                    {
                        idx_cell_pivot = (i + num_ghosts[0]) +
                            (interior_box_lo_idx[1] + num_ghosts[1])*ghostcell_dims[0];
                    }
                    else if (edge_loc_1 == BdryLoc::YHI)
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
        else if (bdry_node_conds[node_loc] == BdryCond::Basic::XREFLECT)
        {
            for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
            {
                for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                {
                    const int idx_cell = (i + num_ghosts[0]) +
                        (j + num_ghosts[1])*ghostcell_dims[0];
                    
                    int idx_cell_pivot = idx_cell;
                    
                    if (edge_loc_0 == BdryLoc::XLO)
                    {
                        idx_cell_pivot = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) + num_ghosts[0]) +
                            (j + num_ghosts[1])*ghostcell_dims[0];
                    }
                    else if (edge_loc_0 == BdryLoc::XHI)
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
        else if (bdry_node_conds[node_loc] == BdryCond::Basic::YREFLECT)
        {
            for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
            {
                for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                {
                    const int idx_cell = (i + num_ghosts[0]) +
                        (j + num_ghosts[1])*ghostcell_dims[0];
                    
                    int idx_cell_pivot = idx_cell;
                    
                    if (edge_loc_1 == BdryLoc::YLO)
                    {
                        idx_cell_pivot = (i + num_ghosts[0]) +
                            (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) + num_ghosts[1])*ghostcell_dims[0];
                    }
                    else if (edge_loc_1 == BdryLoc::YHI)
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
        else if (bdry_node_conds[node_loc] == BdryCond::Basic::XSYMMETRY)
        {
            for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
            {
                for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                {
                    const int idx_cell = (i + num_ghosts[0]) +
                        (j + num_ghosts[1])*ghostcell_dims[0];
                    
                    int idx_cell_pivot = idx_cell;
                    
                    if (edge_loc_0 == BdryLoc::XLO)
                    {
                        idx_cell_pivot = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) + num_ghosts[0]) +
                            (j + num_ghosts[1])*ghostcell_dims[0];
                    }
                    else if (edge_loc_0 == BdryLoc::XHI)
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
        else if (bdry_node_conds[node_loc] == BdryCond::Basic::YSYMMETRY)
        {
            for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
            {
                for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                {
                    const int idx_cell = (i + num_ghosts[0]) +
                        (j + num_ghosts[1])*ghostcell_dims[0];
                    
                    int idx_cell_pivot = idx_cell;
                    
                    if (edge_loc_1 == BdryLoc::YLO)
                    {
                        idx_cell_pivot = (i + num_ghosts[0]) +
                            (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) + num_ghosts[1])*ghostcell_dims[0];
                    }
                    else if (edge_loc_1 == BdryLoc::YHI)
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
            TBOX_ERROR("CartesianBoundaryUtilities2: fillNodeBoundaryData2D()\n"
                << "Invalid node boundary condition!\n"
                << "node_loc = " << node_loc << std::endl
                << "bdry_node_conds[node_loc] = " << bdry_node_conds[node_loc]
                << std::endl);
        }
    }
}


/*
 * Function that returns the integer edge boundary location
 * corresponding to the given node location and node boundary
 * condition.
 *
 * If the node boundary condition type or node location are unknown,
 * or the boundary condition type is inconsistant with the node location
 * an error results.
 */
int
CartesianBoundaryUtilities2::getEdgeLocationForNodeBdry(
    int node_loc,
    int node_btype)
{
    int ret_edge = -1;
    
    switch (node_btype)
    {
        case BdryCond::Basic::XFLOW:
        case BdryCond::Basic::XREFLECT:
        case BdryCond::Basic::XSYMMETRY:
        case BdryCond::Basic::XDIRICHLET:
        case BdryCond::Basic::XNEUMANN:
        {
            if (node_loc == NodeBdyLoc2D::XLO_YLO ||
                node_loc == NodeBdyLoc2D::XLO_YHI)
            {
                ret_edge = BdryLoc::XLO;
            }
            else
            {
                ret_edge = BdryLoc::XHI;
            }
            break;
        }
        case BdryCond::Basic::YFLOW:
        case BdryCond::Basic::YREFLECT:
        case BdryCond::Basic::YSYMMETRY:
        case BdryCond::Basic::YDIRICHLET:
        case BdryCond::Basic::YNEUMANN:
        {
            if (node_loc == NodeBdyLoc2D::XLO_YLO ||
                node_loc == NodeBdyLoc2D::XHI_YLO)
            {
                ret_edge = BdryLoc::YLO;
            }
            else
            {
                ret_edge = BdryLoc::YHI;
            }
            break;
        }
        default:
        {
            TBOX_ERROR("Unknown node boundary condition type = "
                << node_btype
                << " passed to \n"
                << "CartesianBoundaryUtilities2::getEdgeLocationForNodeBdry"
                << std::endl);
        }
    }
    
    if (ret_edge == -1)
    {
        TBOX_ERROR("Node boundary condition type = "
            << node_btype << " and node location = " << node_loc
            << "\n passed to "
            << "CartesianBoundaryUtilities2::getEdgeLocationForNodeBdry"
            << " are inconsistant." << std::endl);
    }
    
    return ret_edge;
}


/*
 * Private function to read 2D edge boundary data from input database.
 */
void
CartesianBoundaryUtilities2::read2dBdryEdges(
    BoundaryUtilityStrategy* bdry_strategy,
    const boost::shared_ptr<tbox::Database>& input_db,
    std::vector<int>& edge_conds,
    const hier::IntVector& periodic)
{
    TBOX_DIM_ASSERT(periodic.getDim() == tbox::Dimension(2));
    
    TBOX_ASSERT(bdry_strategy != 0);
    TBOX_ASSERT(input_db);
    TBOX_ASSERT(static_cast<int>(edge_conds.size()) == NUM_2D_EDGES);
    
    int num_per_dirs = 0;
    for (int id = 0; id < 2; id++)
    {
        if (periodic(id))
            ++num_per_dirs;
    }
    
    if (num_per_dirs < 2)
    {
        // face boundary input required
        for (int s = 0; s < NUM_2D_EDGES; s++)
        {
            std::string bdry_loc_str;
            switch (s)
            {
                case BdryLoc::XLO:
                {
                    bdry_loc_str = "boundary_edge_xlo";
                    break;
                }
                case BdryLoc::XHI:
                {
                    bdry_loc_str = "boundary_edge_xhi";
                    break;
                }
                case BdryLoc::YLO:
                {
                    bdry_loc_str = "boundary_edge_ylo";
                    break;
                }
                case BdryLoc::YHI:
                {
                    bdry_loc_str = "boundary_edge_yhi";
                    break;
                }
                default: NULL_STATEMENT;
            }
            
            bool need_data_read = true;
            if (num_per_dirs > 0)
            {
                if (periodic(0) && (s == BdryLoc::XLO || s == BdryLoc::XHI))
                {
                    need_data_read = false;
                }
                else if (periodic(1) && (s == BdryLoc::YLO || s == BdryLoc::YHI))
                {
                    need_data_read = false;
                }
            }
            
            if (need_data_read)
            {
                boost::shared_ptr<tbox::Database> bdry_loc_db(
                    input_db->getDatabase(bdry_loc_str));
                std::string bdry_cond_str =
                    bdry_loc_db->getString("boundary_condition");
                if (bdry_cond_str == "FLOW")
                {
                    edge_conds[s] = BdryCond::Basic::FLOW;
                }
                else if (bdry_cond_str == "REFLECT")
                {
                    edge_conds[s] = BdryCond::Basic::REFLECT;
                }
                else if (bdry_cond_str == "SYMMETRY")
                {
                    edge_conds[s] = BdryCond::Basic::SYMMETRY;
                }
                else if (bdry_cond_str == "DIRICHLET")
                {
                    edge_conds[s] = BdryCond::Basic::DIRICHLET;
                    bdry_strategy->readDirichletBoundaryDataEntry(
                        bdry_loc_db,
                        bdry_loc_str,
                        s);
                }
                else if (bdry_cond_str == "NEUMANN")
                {
                    edge_conds[s] = BdryCond::Basic::NEUMANN;
                    bdry_strategy->readNeumannBoundaryDataEntry(
                        bdry_loc_db,
                        bdry_loc_str,
                        s);
                }
                else
                {
                    TBOX_ERROR("Unknown edge boundary string = "
                        << bdry_cond_str << " found in input." << std::endl);
                }
            } // if (need_data_read)
       } // for (int s = 0 ...
    } // if (num_per_dirs < 2)
}


/*
 * Private function to read 2D node boundary data from input database.
 */
void
CartesianBoundaryUtilities2::read2dBdryNodes(
    const boost::shared_ptr<tbox::Database>& input_db,
    const std::vector<int>& edge_conds,
    std::vector<int>& node_conds,
    const hier::IntVector& periodic)
{
    TBOX_DIM_ASSERT(periodic.getDim() == tbox::Dimension(2));
    
    TBOX_ASSERT(input_db);
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
        for (int s = 0; s < NUM_2D_NODES; ++s)
        {
            std::string bdry_loc_str;
            switch (s)
            {
                case NodeBdyLoc2D::XLO_YLO:
                {
                    bdry_loc_str = "boundary_node_xlo_ylo";
                    break;
                }
                case NodeBdyLoc2D::XHI_YLO:
                {
                    bdry_loc_str = "boundary_node_xhi_ylo";
                    break;
                }
                case NodeBdyLoc2D::XLO_YHI:
                {
                    bdry_loc_str = "boundary_node_xlo_yhi";
                    break;
                }
                case NodeBdyLoc2D::XHI_YHI:
                {
                    bdry_loc_str = "boundary_node_xhi_yhi";
                    break;
                }
                default: NULL_STATEMENT;
            }
            
            boost::shared_ptr<tbox::Database> bdry_loc_db(
                input_db->getDatabase(bdry_loc_str));
            std::string bdry_cond_str =
                bdry_loc_db->getString("boundary_condition");
            if (bdry_cond_str == "XFLOW")
            {
                node_conds[s] = BdryCond::Basic::XFLOW;
            }
            else if (bdry_cond_str == "YFLOW")
            {
                node_conds[s] = BdryCond::Basic::YFLOW;
            }
            else if (bdry_cond_str == "XREFLECT")
            {
                node_conds[s] = BdryCond::Basic::XREFLECT;
            }
            else if (bdry_cond_str == "YREFLECT")
            {
                node_conds[s] = BdryCond::Basic::YREFLECT;
            }
            else if (bdry_cond_str == "XSYMMETRY")
            {
                node_conds[s] = BdryCond::Basic::XSYMMETRY;
            }
            else if (bdry_cond_str == "YSYMMETRY")
            {
                node_conds[s] = BdryCond::Basic::YSYMMETRY;
            }
            else if (bdry_cond_str == "XDIRICHLET")
            {
                node_conds[s] = BdryCond::Basic::XDIRICHLET;
            }
            else if (bdry_cond_str == "YDIRICHLET")
            {
                node_conds[s] = BdryCond::Basic::YDIRICHLET;
            }
            else if (bdry_cond_str == "XNEUMANN")
            {
                node_conds[s] = BdryCond::Basic::XNEUMANN;
            }
            else if (bdry_cond_str == "YNEUMANN")
            {
                node_conds[s] = BdryCond::Basic::YNEUMANN;
            }
            else
            {
                TBOX_ERROR("Unknown node boundary string = "
                    << bdry_cond_str
                    << " found in input."
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
                if (s == NodeBdyLoc2D::XLO_YLO ||
                    s == NodeBdyLoc2D::XLO_YHI)
                {
                    proper_edge = "XLO";
                    if (bdry_cond_str == "XFLOW" &&
                        edge_conds[BdryLoc::XLO] != BdryCond::Basic::FLOW)
                    {
                        no_edge_data_found = true;
                        proper_edge_data = "FLOW";
                    }
                    if (bdry_cond_str == "XDIRICHLET" &&
                        edge_conds[BdryLoc::XLO] != BdryCond::Basic::DIRICHLET)
                    {
                        no_edge_data_found = true;
                        proper_edge_data = "DIRICHLET";
                    }
                    if (bdry_cond_str == "XNEUMANN" &&
                        edge_conds[BdryLoc::XLO] != BdryCond::Basic::NEUMANN)
                    {
                        no_edge_data_found = true;
                        proper_edge_data = "NEUMANN";
                    }
                    if (bdry_cond_str == "XREFLECT" &&
                        edge_conds[BdryLoc::XLO] != BdryCond::Basic::REFLECT)
                    {
                        no_edge_data_found = true;
                        proper_edge_data = "REFLECT";
                    }
                    if (bdry_cond_str == "XSYMMETRY" &&
                        edge_conds[BdryLoc::XLO] != BdryCond::Basic::SYMMETRY)
                    {
                        no_edge_data_found = true;
                        proper_edge_data = "SYMMETRY";
                    }
                }
                else
                {
                    proper_edge = "XHI";
                    if (bdry_cond_str == "XFLOW" &&
                        edge_conds[BdryLoc::XHI] != BdryCond::Basic::FLOW)
                    {
                        no_edge_data_found = true;
                        proper_edge_data = "FLOW";
                    }
                    if (bdry_cond_str == "XDIRICHLET" &&
                        edge_conds[BdryLoc::XHI] != BdryCond::Basic::DIRICHLET)
                    {
                        no_edge_data_found = true;
                        proper_edge_data = "DIRICHLET";
                    }
                    if (bdry_cond_str == "XNEUMANN" &&
                        edge_conds[BdryLoc::XHI] != BdryCond::Basic::NEUMANN)
                    {
                        no_edge_data_found = true;
                        proper_edge_data = "NEUMANN";
                    }
                    if (bdry_cond_str == "XREFLECT" &&
                        edge_conds[BdryLoc::XHI] != BdryCond::Basic::REFLECT)
                    {
                        no_edge_data_found = true;
                        proper_edge_data = "REFLECT";
                    }
                    if (bdry_cond_str == "XSYMMETRY" &&
                        edge_conds[BdryLoc::XHI] != BdryCond::Basic::SYMMETRY)
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
                if (s == NodeBdyLoc2D::XLO_YLO ||
                    s == NodeBdyLoc2D::XHI_YLO)
                {
                    proper_edge = "YLO";
                    if (bdry_cond_str == "YFLOW" &&
                        edge_conds[BdryLoc::YLO] != BdryCond::Basic::FLOW)
                    {
                        no_edge_data_found = true;
                        proper_edge_data = "FLOW";
                    }
                    if (bdry_cond_str == "YDIRICHLET" &&
                        edge_conds[BdryLoc::YLO] != BdryCond::Basic::DIRICHLET)
                    {
                        no_edge_data_found = true;
                        proper_edge_data = "DIRICHLET";
                    }
                    if (bdry_cond_str == "YNEUMANN" &&
                        edge_conds[BdryLoc::YLO] != BdryCond::Basic::NEUMANN)
                    {
                        no_edge_data_found = true;
                        proper_edge_data = "NEUMANN";
                    }
                    if (bdry_cond_str == "YREFLECT" &&
                        edge_conds[BdryLoc::YLO] != BdryCond::Basic::REFLECT)
                    {
                        no_edge_data_found = true;
                        proper_edge_data = "REFLECT";
                    }
                    if (bdry_cond_str == "YSYMMETRY" &&
                        edge_conds[BdryLoc::YLO] != BdryCond::Basic::SYMMETRY)
                    {
                        no_edge_data_found = true;
                        proper_edge_data = "SYMMETRY";
                    }
                }
                else
                {
                    proper_edge = "YHI";
                    if (bdry_cond_str == "YFLOW" &&
                        edge_conds[BdryLoc::YHI] != BdryCond::Basic::FLOW)
                    {
                        no_edge_data_found = true;
                        proper_edge_data = "FLOW";
                    }
                    if (bdry_cond_str == "YDIRICHLET" &&
                        edge_conds[BdryLoc::YHI] != BdryCond::Basic::DIRICHLET)
                    {
                        no_edge_data_found = true;
                        proper_edge_data = "DIRICHLET";
                    }
                    if (bdry_cond_str == "YNEUMANN" &&
                        edge_conds[BdryLoc::YHI] != BdryCond::Basic::NEUMANN)
                    {
                        no_edge_data_found = true;
                        proper_edge_data = "NEUMANN";
                    }
                    if (bdry_cond_str == "YREFLECT" &&
                        edge_conds[BdryLoc::YHI] != BdryCond::Basic::REFLECT)
                    {
                        no_edge_data_found = true;
                        proper_edge_data = "REFLECT";
                    }
                    if (bdry_cond_str == "YSYMMETRY" &&
                        edge_conds[BdryLoc::YHI] != BdryCond::Basic::SYMMETRY)
                    {
                        no_edge_data_found = true;
                        proper_edge_data = "SYMMETRY";
                    }
                }
            }
            if (no_edge_data_found)
            {
                TBOX_ERROR(
                    "Bdry condition "
                    << bdry_cond_str
                    << " found for "
                    << bdry_loc_str
                    << "\n but no "
                    << proper_edge_data
                    << " data found for edge "
                    << proper_edge << std::endl);
            }
        } // for (int s = 0 ...
    } // if (num_per_dirs < 1)
}
