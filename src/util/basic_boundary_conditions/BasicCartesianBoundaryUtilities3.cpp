/*************************************************************************
 *
 * This file is modified from CartesianBoundaryUtilities3.c of the SAMRAI
 * distribution. For full copyright information, see COPYRIGHT and
 * COPYING.LESSER of SAMRAI distribution.
 *
 * Copyright:     (c) 1997-2014 Lawrence Livermore National Security, LLC
 * Description:   Utility routines for manipulating 3D Cartesian boundary data
 *
 ************************************************************************/

#include "util/basic_boundary_conditions/BasicCartesianBoundaryUtilities3.hpp"

#include "util/basic_boundary_conditions/BasicBoundaryConditions.hpp"
#include "util/basic_boundary_conditions/CartesianBoundaryDefines.hpp"

#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"

#include <algorithm>

// Integer constant for debugging improperly set boundary data.
#define BOGUS_BDRY_LOC (-9999)

/*
 * This function reads 3D boundary data from given input database.
 * The integer boundary condition types are placed in the integer
 * arrays supplied by the caller (typically, the concrete BoundaryStrategy
 * provided).  When DIRICHLET or NEUMANN conditions are specified, control
 * is passed to the BoundaryStrategy to read the boundary state data
 * specific to the problem.
 *
 * Errors will be reported and the program will abort whenever necessary
 * boundary condition information is missing in the input database, or
 * when the data read in is either unknown or inconsistent.  The periodic
 * domain information is used to determine which boundary face, edge, or
 * node entries are not required from input.  Error checking requires
 * that node and edge boundary conditions are consistent with those
 * specified for the faces.
 *
 * Arguments are:
 *    bdry_strategy .... object that reads DIRICHLET or NEUMANN conditions
 *    input_db ......... input database containing all boundary data
 *    face_locs ........ array of locations of faces for applying
 *                       boundary conditions.
 *    edge_locs ........ array of locations of edges for applying
 *                       boundary conditions.
 *    node_locs ........ array of locations of nodes for applying
 *                       boundary conditions.
 *    face_conds ....... array into which integer boundary conditions
 *                       for faces are read
 *    edge_conds ....... array into which integer boundary conditions
 *                       for edges are read
 *    node_conds ....... array into which integer boundary conditions
 *                       for nodes are read
 *    periodic ......... integer vector specifying which coordinate
 *                       directions are periodic (value returned from
 *                       GridGeometry3::getPeriodicShift())
 */
void
BasicCartesianBoundaryUtilities3::getFromInput(
    BoundaryUtilityStrategy* bdry_strategy,
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
    
    TBOX_ASSERT(bdry_strategy != 0);
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
        TBOX_ERROR("BasicCartesianBoundaryUtilities3::getFromInput()\n"
            << "No input database supplied."
            << std::endl);
    }
    
    if (static_cast<int>(face_locs.size()) > 0)
    {
        read3dBdryFaces(
            bdry_strategy,
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
 * Function to remove 3d boundary faces with boundary conditions filled by this class for a patch.
 *
 * Arguments are:
 *    bdry_face_locs ....... array of locations of faces for applying
 *                           boundary conditions.
 *    patch ................ patch on which data object lives
 *    bdry_face_conds ...... array of boundary conditions for patch faces
 */
void
BasicCartesianBoundaryUtilities3::removeBoundaryFaceLocations(
    std::vector<int>& bdry_face_locs,
    const hier::Patch& patch,
    const std::vector<int>& bdry_face_conds)
{
    TBOX_ASSERT(static_cast<int>(bdry_face_locs.size()) <= NUM_3D_FACES);
    TBOX_ASSERT(*min_element(bdry_face_locs.begin(), bdry_face_locs.end()) >= 0);
    TBOX_ASSERT(*max_element(bdry_face_locs.begin(), bdry_face_locs.end()) < NUM_3D_FACES);
    TBOX_ASSERT(static_cast<int>(bdry_face_conds.size()) == NUM_3D_FACES);
    
    const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
        HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    TBOX_ASSERT(patch_geom);
    
    const std::vector<hier::BoundaryBox>& face_bdry =
        patch_geom->getCodimensionBoundaries(BDRY::FACE3D);
    
    for (int fi = 0; fi < static_cast<int>(face_bdry.size()); fi++)
    {
        TBOX_ASSERT(face_bdry[fi].getBoundaryType() == BDRY::FACE3D);
        
        int face_loc = face_bdry[fi].getLocationIndex();
        
        if (std::find(bdry_face_locs.begin(), bdry_face_locs.end(), face_loc) !=
            bdry_face_locs.end())
        {
            if (bdry_face_conds[face_loc] == BDRY_COND::BASIC::DIRICHLET)
            {
                // Remove face locations that have boundary conditions identified.
                bdry_face_locs.erase(std::remove(bdry_face_locs.begin(), bdry_face_locs.end(), face_loc),
                    bdry_face_locs.end());
            }
            else if (bdry_face_conds[face_loc] == BDRY_COND::BASIC::NEUMANN)
            {
                // NOT YET IMPLEMENTED
            }
            else if (bdry_face_conds[face_loc] == BDRY_COND::BASIC::FLOW)
            {
                // Remove face locations that have boundary conditions identified.
                bdry_face_locs.erase(std::remove(bdry_face_locs.begin(), bdry_face_locs.end(), face_loc),
                    bdry_face_locs.end());
            }
            else if (bdry_face_conds[face_loc] == BDRY_COND::BASIC::REFLECT)
            {
                // Remove face locations that have boundary conditions identified.
                bdry_face_locs.erase(std::remove(bdry_face_locs.begin(), bdry_face_locs.end(), face_loc),
                    bdry_face_locs.end());
            }
            else if (bdry_face_conds[face_loc] == BDRY_COND::BASIC::SYMMETRY)
            {
                // Remove face locations that have boundary conditions identified.
                bdry_face_locs.erase(std::remove(bdry_face_locs.begin(), bdry_face_locs.end(), face_loc),
                    bdry_face_locs.end());
            }
        }
    }
}


/*
 * Function to fill face boundary values.
 *
 * Arguments are:
 *    var_name ............. name of variable (for error reporting)
 *    var_data ............. cell-centered patch data object to check
 *    patch ................ patch on which data object lives
 *    bdry_face_locs ....... array of locations of faces for applying
 *                           boundary conditions.
 *    bdry_face_conds ...... array of boundary conditions for patch faces
 *    bdry_face_values ..... array of boundary values for faces
 *                           (this must be consistent with boundary
 *                           condition types)
 *    ghost_width_to_fill .. width of ghost region to fill
 */
void
BasicCartesianBoundaryUtilities3::fillFaceBoundaryData(
    const std::string& var_name,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& var_data,
    const hier::Patch& patch,
    const std::vector<int>& bdry_face_locs,
    const std::vector<int>& bdry_face_conds,
    const std::vector<double>& bdry_face_values,
    const hier::IntVector& ghost_width_to_fill)
{
    TBOX_ASSERT(!var_name.empty());
    TBOX_ASSERT(var_data);
    TBOX_ASSERT(static_cast<int>(bdry_face_locs.size()) <= NUM_3D_FACES);
    TBOX_ASSERT(*min_element(bdry_face_locs.begin(), bdry_face_locs.end()) >= 0);
    TBOX_ASSERT(*max_element(bdry_face_locs.begin(), bdry_face_locs.end()) < NUM_3D_FACES);
    TBOX_ASSERT(static_cast<int>(bdry_face_conds.size()) == NUM_3D_FACES);
    TBOX_ASSERT(static_cast<int>(bdry_face_values.size()) == NUM_3D_FACES*(var_data->getDepth()));
    
    TBOX_DIM_ASSERT(ghost_width_to_fill.getDim() == tbox::Dimension(3));
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
    hier::IntVector gcw_to_fill(tbox::Dimension(3));
    
    // If the ghost fill width is not used, it is set to the ghost cell width of the data.
    if (ghost_width_to_fill == -hier::IntVector::getOne(tbox::Dimension(3)))
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
    
    const std::vector<hier::BoundaryBox>& face_bdry =
        patch_geom->getCodimensionBoundaries(BDRY::FACE3D);
    
    const int var_depth = var_data->getDepth();
    
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
            
            if (bdry_face_conds[face_loc] == BDRY_COND::BASIC::DIRICHLET)
            {
                for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                {
                    for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                    {
                        for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                        {
                            const int idx_cell = (i + num_ghosts[0]) +
                                (j + num_ghosts[1])*ghostcell_dims[0] +
                                (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            
                            for (int di = 0; di < var_depth; di++)
                            {
                                var_data->getPointer(di)[idx_cell] = bdry_face_values[face_loc*var_depth + di];
                            }
                        }
                    }
                }
            }
            else if (bdry_face_conds[face_loc] == BDRY_COND::BASIC::NEUMANN)
            {
                // NOT YET IMPLEMENTED
            }
            else if (bdry_face_conds[face_loc] == BDRY_COND::BASIC::FLOW)
            {
                for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                {
                    for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                    {
                        for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                        {
                            const int idx_cell = (i + num_ghosts[0]) +
                                (j + num_ghosts[1])*ghostcell_dims[0] +
                                (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            
                            int idx_cell_pivot = idx_cell;
                            
                            if (face_loc == BDRY_LOC::XLO)
                            {
                                idx_cell_pivot = (interior_box_lo_idx[0] + num_ghosts[0]) +
                                    (j + num_ghosts[1])*ghostcell_dims[0] +
                                    (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            }
                            else if (face_loc == BDRY_LOC::XHI)
                            {
                                idx_cell_pivot = (interior_box_hi_idx[0] + num_ghosts[0]) +
                                    (j + num_ghosts[1])*ghostcell_dims[0] +
                                    (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            }
                            else if (face_loc == BDRY_LOC::YLO)
                            {
                                idx_cell_pivot = (i + num_ghosts[0]) +
                                    (interior_box_lo_idx[1] + num_ghosts[1])*ghostcell_dims[0] +
                                    (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            }
                            else if (face_loc == BDRY_LOC::YHI)
                            {
                                idx_cell_pivot = (i + num_ghosts[0]) +
                                    (interior_box_hi_idx[1] + num_ghosts[1])*ghostcell_dims[0] +
                                    (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            }
                            else if (face_loc == BDRY_LOC::ZLO)
                            {
                                idx_cell_pivot = (i + num_ghosts[0]) +
                                    (j + num_ghosts[1])*ghostcell_dims[0] +
                                    (interior_box_lo_idx[2] + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            }
                            else if (face_loc == BDRY_LOC::ZHI)
                            {
                                idx_cell_pivot = (i + num_ghosts[0]) +
                                    (j + num_ghosts[1])*ghostcell_dims[0] +
                                    (interior_box_hi_idx[2] + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            }
                            
                            for (int di = 0; di < var_depth; di++)
                            {
                                var_data->getPointer(di)[idx_cell] = var_data->getPointer(di)[idx_cell_pivot];
                            }
                        }
                    }
                }
            }
            else if (bdry_face_conds[face_loc] == BDRY_COND::BASIC::REFLECT)
            {
                for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                {
                    for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                    {
                        for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                        {
                            const int idx_cell = (i + num_ghosts[0]) +
                                (j + num_ghosts[1])*ghostcell_dims[0] +
                                (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            
                            int idx_cell_pivot = idx_cell;
                            
                            if (face_loc == BDRY_LOC::XLO)
                            {
                                idx_cell_pivot = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) + num_ghosts[0]) +
                                    (j + num_ghosts[1])*ghostcell_dims[0] +
                                    (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            }
                            else if (face_loc == BDRY_LOC::XHI)
                            {
                                idx_cell_pivot = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) + num_ghosts[0]) +
                                    (j + num_ghosts[1])*ghostcell_dims[0] +
                                    (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            }
                            else if (face_loc == BDRY_LOC::YLO)
                            {
                                idx_cell_pivot = (i + num_ghosts[0]) +
                                    (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) + num_ghosts[1])*ghostcell_dims[0] +
                                    (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            }
                            else if (face_loc == BDRY_LOC::YHI)
                            {
                                idx_cell_pivot = (i + num_ghosts[0]) +
                                    (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) + num_ghosts[1])*ghostcell_dims[0] +
                                    (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            }
                            else if (face_loc == BDRY_LOC::ZLO)
                            {
                                idx_cell_pivot = (i + num_ghosts[0]) +
                                    (j + num_ghosts[1])*ghostcell_dims[0] +
                                    (interior_box_lo_idx[2] + (fill_box_hi_idx[2] - k) + num_ghosts[2])*ghostcell_dims[0]*
                                        ghostcell_dims[1];
                            }
                            else if (face_loc == BDRY_LOC::ZHI)
                            {
                                idx_cell_pivot = (i + num_ghosts[0]) +
                                    (j + num_ghosts[1])*ghostcell_dims[0] +
                                    (interior_box_hi_idx[2] - (k - fill_box_lo_idx[2]) + num_ghosts[2])*ghostcell_dims[0]*
                                        ghostcell_dims[1];
                            }
                            
                            for (int di = 0; di < var_depth; di++)
                            {
                                var_data->getPointer(di)[idx_cell] = var_data->getPointer(di)[idx_cell_pivot];
                            }
                            
                            if (face_loc == BDRY_LOC::XLO || face_loc == BDRY_LOC::XHI)
                            {
                                var_data->getPointer(0)[idx_cell] = -var_data->getPointer(0)[idx_cell_pivot];
                            }
                            else if (face_loc == BDRY_LOC::YLO || face_loc == BDRY_LOC::YHI)
                            {
                                var_data->getPointer(1)[idx_cell] = -var_data->getPointer(1)[idx_cell_pivot];
                            }
                            else if (face_loc == BDRY_LOC::ZLO || face_loc == BDRY_LOC::ZHI)
                            {
                                var_data->getPointer(2)[idx_cell] = -var_data->getPointer(2)[idx_cell_pivot];
                            }
                        }
                    }
                }
            }
            else if (bdry_face_conds[face_loc] == BDRY_COND::BASIC::SYMMETRY)
            {
                for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                {
                    for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                    {
                        for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                        {
                            const int idx_cell = (i + num_ghosts[0]) +
                                (j + num_ghosts[1])*ghostcell_dims[0] +
                                (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            
                            int idx_cell_pivot = idx_cell;
                            
                            if (face_loc == BDRY_LOC::XLO)
                            {
                                idx_cell_pivot = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) + num_ghosts[0]) +
                                    (j + num_ghosts[1])*ghostcell_dims[0] +
                                    (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            }
                            else if (face_loc == BDRY_LOC::XHI)
                            {
                                idx_cell_pivot = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) + num_ghosts[0]) +
                                    (j + num_ghosts[1])*ghostcell_dims[0] +
                                    (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            }
                            else if (face_loc == BDRY_LOC::YLO)
                            {
                                idx_cell_pivot = (i + num_ghosts[0]) +
                                    (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) + num_ghosts[1])*ghostcell_dims[0] +
                                    (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            }
                            else if (face_loc == BDRY_LOC::YHI)
                            {
                                idx_cell_pivot = (i + num_ghosts[0]) +
                                    (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) + num_ghosts[1])*ghostcell_dims[0] +
                                    (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            }
                            else if (face_loc == BDRY_LOC::ZLO)
                            {
                                idx_cell_pivot = (i + num_ghosts[0]) +
                                    (j + num_ghosts[1])*ghostcell_dims[0] +
                                    (interior_box_lo_idx[2] + (fill_box_hi_idx[2] - k) + num_ghosts[2])*ghostcell_dims[0]*
                                        ghostcell_dims[1];
                            }
                            else if (face_loc == BDRY_LOC::ZHI)
                            {
                                idx_cell_pivot = (i + num_ghosts[0]) +
                                    (j + num_ghosts[1])*ghostcell_dims[0] +
                                    (interior_box_hi_idx[2] - (k - fill_box_lo_idx[2]) + num_ghosts[2])*ghostcell_dims[0]*
                                        ghostcell_dims[1];
                            }
                            
                            for (int di = 0; di < var_depth; di++)
                            {
                                var_data->getPointer(di)[idx_cell] = var_data->getPointer(di)[idx_cell_pivot];
                            }
                        }
                    }
                }
            }
        }
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
 *    bdry_face_values ..... array of boundary values for faces
 *                           (this must be consistent with boundary
 *                           condition types)
 *    ghost_width_to_fill .. width of ghost region to fill
 */
void
BasicCartesianBoundaryUtilities3::removeBoundaryEdgeLocations(
            std::vector<int>& bdry_edge_locs,
            const hier::Patch& patch,
            const std::vector<int>& bdry_edge_conds)
{
    TBOX_ASSERT(static_cast<int>(bdry_edge_locs.size()) <= NUM_3D_EDGES);
    TBOX_ASSERT(*min_element(bdry_edge_locs.begin(), bdry_edge_locs.end()) >= 0);
    TBOX_ASSERT(*max_element(bdry_edge_locs.begin(), bdry_edge_locs.end()) < NUM_3D_EDGES);
    TBOX_ASSERT(static_cast<int>(bdry_edge_conds.size()) == NUM_3D_EDGES);
    
    const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
        HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    TBOX_ASSERT(patch_geom);
    
    const std::vector<hier::BoundaryBox>& edge_bdry =
        patch_geom->getCodimensionBoundaries(BDRY::EDGE3D);
    
    for (int ei = 0; ei < static_cast<int>(edge_bdry.size()); ei++)
    {
        TBOX_ASSERT(edge_bdry[ei].getBoundaryType() == BDRY::EDGE3D);
        
        int edge_loc(edge_bdry[ei].getLocationIndex());
        
        if (std::find(bdry_edge_locs.begin(), bdry_edge_locs.end(), edge_loc) !=
            bdry_edge_locs.end())
        {
            if (bdry_edge_conds[edge_loc] == BDRY_COND::BASIC::XDIRICHLET)
            {
                // Remove edge locations that have boundary conditions identified.
                bdry_edge_locs.erase(std::remove(bdry_edge_locs.begin(), bdry_edge_locs.end(), edge_loc),
                    bdry_edge_locs.end());
            }
            else if (bdry_edge_conds[edge_loc] == BDRY_COND::BASIC::YDIRICHLET)
            {
                // Remove edge locations that have boundary conditions identified.
                bdry_edge_locs.erase(std::remove(bdry_edge_locs.begin(), bdry_edge_locs.end(), edge_loc),
                    bdry_edge_locs.end());
            }
            else if (bdry_edge_conds[edge_loc] == BDRY_COND::BASIC::ZDIRICHLET)
            {
                // Remove edge locations that have boundary conditions identified.
                bdry_edge_locs.erase(std::remove(bdry_edge_locs.begin(), bdry_edge_locs.end(), edge_loc),
                    bdry_edge_locs.end());
            }
            else if (bdry_edge_conds[edge_loc] == BDRY_COND::BASIC::XNEUMANN)
            {
                // NOT YET IMPLEMENTED
            }
            else if (bdry_edge_conds[edge_loc] == BDRY_COND::BASIC::YNEUMANN)
            {
                // NOT YET IMPLEMENTED
            }
            else if (bdry_edge_conds[edge_loc] == BDRY_COND::BASIC::ZNEUMANN)
            {
                // NOT YET IMPLEMENTED
            }
            else if (bdry_edge_conds[edge_loc] == BDRY_COND::BASIC::XFLOW)
            {
                // Remove edge locations that have boundary conditions identified.
                bdry_edge_locs.erase(std::remove(bdry_edge_locs.begin(), bdry_edge_locs.end(), edge_loc),
                    bdry_edge_locs.end());
            }
            else if (bdry_edge_conds[edge_loc] == BDRY_COND::BASIC::YFLOW)
            {
                // Remove edge locations that have boundary conditions identified.
                bdry_edge_locs.erase(std::remove(bdry_edge_locs.begin(), bdry_edge_locs.end(), edge_loc),
                    bdry_edge_locs.end());
            }
            else if (bdry_edge_conds[edge_loc] == BDRY_COND::BASIC::ZFLOW)
            {
                // Remove edge locations that have boundary conditions identified.
                bdry_edge_locs.erase(std::remove(bdry_edge_locs.begin(), bdry_edge_locs.end(), edge_loc),
                    bdry_edge_locs.end());
            }
            else if (bdry_edge_conds[edge_loc] == BDRY_COND::BASIC::XREFLECT)
            {
                // Remove edge locations that have boundary conditions identified.
                bdry_edge_locs.erase(std::remove(bdry_edge_locs.begin(), bdry_edge_locs.end(), edge_loc),
                    bdry_edge_locs.end());
            }
            else if (bdry_edge_conds[edge_loc] == BDRY_COND::BASIC::YREFLECT)
            {
                // Remove edge locations that have boundary conditions identified.
                bdry_edge_locs.erase(std::remove(bdry_edge_locs.begin(), bdry_edge_locs.end(), edge_loc),
                    bdry_edge_locs.end());
            }
            else if (bdry_edge_conds[edge_loc] == BDRY_COND::BASIC::ZREFLECT)
            {
                // Remove edge locations that have boundary conditions identified.
                bdry_edge_locs.erase(std::remove(bdry_edge_locs.begin(), bdry_edge_locs.end(), edge_loc),
                    bdry_edge_locs.end());
            }
            else if (bdry_edge_conds[edge_loc] == BDRY_COND::BASIC::XSYMMETRY)
            {
                // Remove edge locations that have boundary conditions identified.
                bdry_edge_locs.erase(std::remove(bdry_edge_locs.begin(), bdry_edge_locs.end(), edge_loc),
                    bdry_edge_locs.end());
            }
            else if (bdry_edge_conds[edge_loc] == BDRY_COND::BASIC::YSYMMETRY)
            {
                // Remove edge locations that have boundary conditions identified.
                bdry_edge_locs.erase(std::remove(bdry_edge_locs.begin(), bdry_edge_locs.end(), edge_loc),
                    bdry_edge_locs.end());
            }
            else if (bdry_edge_conds[edge_loc] == BDRY_COND::BASIC::ZSYMMETRY)
            {
                // Remove edge locations that have boundary conditions identified.
                bdry_edge_locs.erase(std::remove(bdry_edge_locs.begin(), bdry_edge_locs.end(), edge_loc),
                    bdry_edge_locs.end());
            }
        }
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
 *    bdry_face_values ..... array of boundary values for faces
 *                           (this must be consistent with boundary
 *                           condition types)
 *    ghost_width_to_fill .. width of ghost region to fill
 */
void
BasicCartesianBoundaryUtilities3::fillEdgeBoundaryData(
    const std::string& var_name,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& var_data,
    const hier::Patch& patch,
    const std::vector<int>& bdry_edge_locs,
    const std::vector<int>& bdry_edge_conds,
    const std::vector<double>& bdry_face_values,
    const hier::IntVector& ghost_width_to_fill)
{
    TBOX_ASSERT(!var_name.empty());
    TBOX_ASSERT(var_data);
    TBOX_ASSERT(static_cast<int>(bdry_edge_locs.size()) <= NUM_3D_EDGES);
    TBOX_ASSERT(*min_element(bdry_edge_locs.begin(), bdry_edge_locs.end()) >= 0);
    TBOX_ASSERT(*max_element(bdry_edge_locs.begin(), bdry_edge_locs.end()) < NUM_3D_EDGES);
    TBOX_ASSERT(static_cast<int>(bdry_edge_conds.size()) == NUM_3D_EDGES);
    TBOX_ASSERT(static_cast<int>(bdry_face_values.size()) == NUM_3D_FACES*(var_data->getDepth()));
    
    TBOX_DIM_ASSERT(ghost_width_to_fill.getDim() == tbox::Dimension(3));
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
    hier::IntVector gcw_to_fill(tbox::Dimension(3));
    
    // If the ghost fill width is not used, it is set to the ghost cell width of the data.
    if (ghost_width_to_fill == -hier::IntVector::getOne(tbox::Dimension(3)))
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
        patch_geom->getCodimensionBoundaries(BDRY::EDGE3D);
    
    const int var_depth = var_data->getDepth();
    
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
            
            if (bdry_edge_conds[edge_loc] == BDRY_COND::BASIC::XDIRICHLET)
            {
                for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                {
                    for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                    {
                        for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                        {
                            const int idx_cell = (i + num_ghosts[0]) +
                                (j + num_ghosts[1])*ghostcell_dims[0] +
                                (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                
                            for (int di = 0; di < var_depth; di++)
                            {
                                var_data->getPointer(di)[idx_cell] = bdry_face_values[face_loc_0*var_depth + di];
                            }
                        }
                    }
                }
            }
            else if (bdry_edge_conds[edge_loc] == BDRY_COND::BASIC::YDIRICHLET)
            {
                for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                {
                    for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                    {
                        for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                        {
                            const int idx_cell = (i + num_ghosts[0]) +
                                (j + num_ghosts[1])*ghostcell_dims[0] +
                                (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                
                            for (int di = 0; di < var_depth; di++)
                            {
                                var_data->getPointer(di)[idx_cell] = bdry_face_values[face_loc_1*var_depth + di];
                            }
                        }
                    }
                }
            }
            else if (bdry_edge_conds[edge_loc] == BDRY_COND::BASIC::ZDIRICHLET)
            {
                for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                {
                    for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                    {
                        for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                        {
                            const int idx_cell = (i + num_ghosts[0]) +
                                (j + num_ghosts[1])*ghostcell_dims[0] +
                                (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                
                            for (int di = 0; di < var_depth; di++)
                            {
                                var_data->getPointer(di)[idx_cell] = bdry_face_values[face_loc_2*var_depth + di];
                            }
                        }
                    }
                }
            }
            else if (bdry_edge_conds[edge_loc] == BDRY_COND::BASIC::XNEUMANN)
            {
                // NOT YET IMPLEMENTED
            }
            else if (bdry_edge_conds[edge_loc] == BDRY_COND::BASIC::YNEUMANN)
            {
                // NOT YET IMPLEMENTED
            }
            else if (bdry_edge_conds[edge_loc] == BDRY_COND::BASIC::ZNEUMANN)
            {
                // NOT YET IMPLEMENTED
            }
            else if (bdry_edge_conds[edge_loc] == BDRY_COND::BASIC::XFLOW)
            {
                for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                {
                    for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                    {
                        for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                        {
                            const int idx_cell = (i + num_ghosts[0]) +
                                (j + num_ghosts[1])*ghostcell_dims[0] +
                                (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            
                            int idx_cell_pivot = idx_cell;
                            
                            if (face_loc_0 == BDRY_LOC::XLO)
                            {
                                idx_cell_pivot = (interior_box_lo_idx[0] + num_ghosts[0]) +
                                    (j + num_ghosts[1])*ghostcell_dims[0] +
                                    (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            }
                            else if (face_loc_0 == BDRY_LOC::XHI)
                            {
                                idx_cell_pivot = (interior_box_hi_idx[0] + num_ghosts[0]) +
                                    (j + num_ghosts[1])*ghostcell_dims[0] +
                                    (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            }
                            
                            for (int di = 0; di < var_depth; di++)
                            {
                                var_data->getPointer(di)[idx_cell] = var_data->getPointer(di)[idx_cell_pivot];
                            }
                        }
                    }
                }
            }
            else if (bdry_edge_conds[edge_loc] == BDRY_COND::BASIC::YFLOW)
            {
                for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                {
                    for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                    {
                        for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                        {
                            const int idx_cell = (i + num_ghosts[0]) +
                                (j + num_ghosts[1])*ghostcell_dims[0] +
                                (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            
                            int idx_cell_pivot = idx_cell;
                            
                            if (face_loc_1 == BDRY_LOC::YLO)
                            {
                                idx_cell_pivot = (i + num_ghosts[0]) +
                                    (interior_box_lo_idx[1] + num_ghosts[1])*ghostcell_dims[0] +
                                    (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            }
                            else if (face_loc_1 == BDRY_LOC::YHI)
                            {
                                idx_cell_pivot = (i + num_ghosts[0]) +
                                    (interior_box_hi_idx[1] + num_ghosts[1])*ghostcell_dims[0] +
                                    (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            }
                            
                            for (int di = 0; di < var_depth; di++)
                            {
                                var_data->getPointer(di)[idx_cell] = var_data->getPointer(di)[idx_cell_pivot];
                            }
                        }
                    }
                }
            }
            else if (bdry_edge_conds[edge_loc] == BDRY_COND::BASIC::ZFLOW)
            {
                for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                {
                    for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                    {
                        for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                        {
                            const int idx_cell = (i + num_ghosts[0]) +
                                (j + num_ghosts[1])*ghostcell_dims[0] +
                                (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            
                            int idx_cell_pivot = idx_cell;
                            
                            if (face_loc_2 == BDRY_LOC::ZLO)
                            {
                                idx_cell_pivot = (i + num_ghosts[0]) +
                                    (j + num_ghosts[1])*ghostcell_dims[0] +
                                    (interior_box_lo_idx[2] + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            }
                            else if (face_loc_2 == BDRY_LOC::ZHI)
                            {
                                idx_cell_pivot = (i + num_ghosts[0]) +
                                    (j + num_ghosts[1])*ghostcell_dims[0] +
                                    (interior_box_hi_idx[2] + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            }
                            
                            for (int di = 0; di < var_depth; di++)
                            {
                                var_data->getPointer(di)[idx_cell] = var_data->getPointer(di)[idx_cell_pivot];
                            }
                        }
                    }
                }
            }
            else if (bdry_edge_conds[edge_loc] == BDRY_COND::BASIC::XREFLECT)
            {
                for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                {
                    for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                    {
                        for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                        {
                            const int idx_cell = (i + num_ghosts[0]) +
                                (j + num_ghosts[1])*ghostcell_dims[0] +
                                (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            
                            int idx_cell_pivot = idx_cell;
                            
                            if (face_loc_0 == BDRY_LOC::XLO)
                            {
                                idx_cell_pivot = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) + num_ghosts[0]) +
                                    (j + num_ghosts[1])*ghostcell_dims[0] +
                                    (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            }
                            else if (face_loc_0 == BDRY_LOC::XHI)
                            {
                                idx_cell_pivot = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) + num_ghosts[0]) +
                                    (j + num_ghosts[1])*ghostcell_dims[0] +
                                    (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            }
                            
                            for (int di = 0; di < var_depth; di++)
                            {
                                var_data->getPointer(di)[idx_cell] = var_data->getPointer(di)[idx_cell_pivot];
                            }
                            
                            var_data->getPointer(0)[idx_cell] = -var_data->getPointer(0)[idx_cell_pivot];
                        }
                    }
                }
            }
            else if (bdry_edge_conds[edge_loc] == BDRY_COND::BASIC::YREFLECT)
            {
                for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                {
                    for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                    {
                        for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                        {
                            const int idx_cell = (i + num_ghosts[0]) +
                                (j + num_ghosts[1])*ghostcell_dims[0] +
                                (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            
                            int idx_cell_pivot = idx_cell;
                            
                            if (face_loc_1 == BDRY_LOC::YLO)
                            {
                                idx_cell_pivot = (i + num_ghosts[0]) +
                                    (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) + num_ghosts[1])*ghostcell_dims[0] +
                                    (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            }
                            else if (face_loc_1 == BDRY_LOC::YHI)
                            {
                                idx_cell_pivot = (i + num_ghosts[0]) +
                                    (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) + num_ghosts[1])*ghostcell_dims[0] +
                                    (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            }
                            
                            for (int di = 0; di < var_depth; di++)
                            {
                                var_data->getPointer(di)[idx_cell] = var_data->getPointer(di)[idx_cell_pivot];
                            }
                            
                            var_data->getPointer(1)[idx_cell] = -var_data->getPointer(1)[idx_cell_pivot];
                        }
                    }
                }
            }
            else if (bdry_edge_conds[edge_loc] == BDRY_COND::BASIC::ZREFLECT)
            {
                for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                {
                    for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                    {
                        for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                        {
                            const int idx_cell = (i + num_ghosts[0]) +
                                (j + num_ghosts[1])*ghostcell_dims[0] +
                                (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            
                            int idx_cell_pivot = idx_cell;
                            
                            if (face_loc_2 == BDRY_LOC::ZLO)
                            {
                                idx_cell_pivot = (i + num_ghosts[0]) +
                                    (j + num_ghosts[1])*ghostcell_dims[0] +
                                    (interior_box_lo_idx[2] + (fill_box_hi_idx[2] - k) + num_ghosts[2])
                                        *ghostcell_dims[0]*ghostcell_dims[1];
                            }
                            else if (face_loc_2 == BDRY_LOC::ZHI)
                            {
                                idx_cell_pivot = (i + num_ghosts[0]) +
                                    (j + num_ghosts[1])*ghostcell_dims[0] +
                                    (interior_box_hi_idx[2] - (k - fill_box_lo_idx[2]) + num_ghosts[2])
                                        *ghostcell_dims[0]*ghostcell_dims[1];
                            }
                            
                            for (int di = 0; di < var_depth; di++)
                            {
                                var_data->getPointer(di)[idx_cell] = var_data->getPointer(di)[idx_cell_pivot];
                            }
                            
                            var_data->getPointer(2)[idx_cell] = -var_data->getPointer(2)[idx_cell_pivot];
                        }
                    }
                }
            }
            else if (bdry_edge_conds[edge_loc] == BDRY_COND::BASIC::XSYMMETRY)
            {
                for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                {
                    for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                    {
                        for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                        {
                            const int idx_cell = (i + num_ghosts[0]) +
                                (j + num_ghosts[1])*ghostcell_dims[0] +
                                (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            
                            int idx_cell_pivot = idx_cell;
                            
                            if (face_loc_0 == BDRY_LOC::XLO)
                            {
                                idx_cell_pivot = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) + num_ghosts[0]) +
                                    (j + num_ghosts[1])*ghostcell_dims[0] +
                                    (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            }
                            else if (face_loc_0 == BDRY_LOC::XHI)
                            {
                                idx_cell_pivot = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) + num_ghosts[0]) +
                                    (j + num_ghosts[1])*ghostcell_dims[0] +
                                    (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            }
                            
                            for (int di = 0; di < var_depth; di++)
                            {
                                var_data->getPointer(di)[idx_cell] = var_data->getPointer(di)[idx_cell_pivot];
                            }
                        }
                    }
                }
            }
            else if (bdry_edge_conds[edge_loc] == BDRY_COND::BASIC::YSYMMETRY)
            {
                for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                {
                    for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                    {
                        for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                        {
                            const int idx_cell = (i + num_ghosts[0]) +
                                (j + num_ghosts[1])*ghostcell_dims[0] +
                                (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            
                            int idx_cell_pivot = idx_cell;
                            
                            if (face_loc_1 == BDRY_LOC::YLO)
                            {
                                idx_cell_pivot = (i + num_ghosts[0]) +
                                    (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) + num_ghosts[1])*ghostcell_dims[0] +
                                    (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            }
                            else if (face_loc_1 == BDRY_LOC::YHI)
                            {
                                idx_cell_pivot = (i + num_ghosts[0]) +
                                    (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) + num_ghosts[1])*ghostcell_dims[0] +
                                    (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            }
                            
                            for (int di = 0; di < var_depth; di++)
                            {
                                var_data->getPointer(di)[idx_cell] = var_data->getPointer(di)[idx_cell_pivot];
                            }
                        }
                    }
                }
            }
            else if (bdry_edge_conds[edge_loc] == BDRY_COND::BASIC::ZSYMMETRY)
            {
                for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                {
                    for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                    {
                        for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                        {
                            const int idx_cell = (i + num_ghosts[0]) +
                                (j + num_ghosts[1])*ghostcell_dims[0] +
                                (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            
                            int idx_cell_pivot = idx_cell;
                            
                            if (face_loc_2 == BDRY_LOC::ZLO)
                            {
                                idx_cell_pivot = (i + num_ghosts[0]) +
                                    (j + num_ghosts[1])*ghostcell_dims[0] +
                                    (interior_box_lo_idx[2] + (fill_box_hi_idx[2] - k) + num_ghosts[2])
                                        *ghostcell_dims[0]*ghostcell_dims[1];
                            }
                            else if (face_loc_2 == BDRY_LOC::ZHI)
                            {
                                idx_cell_pivot = (i + num_ghosts[0]) +
                                    (j + num_ghosts[1])*ghostcell_dims[0] +
                                    (interior_box_hi_idx[2] - (k - fill_box_lo_idx[2]) + num_ghosts[2])
                                        *ghostcell_dims[0]*ghostcell_dims[1];
                            }
                            
                            for (int di = 0; di < var_depth; di++)
                            {
                                var_data->getPointer(di)[idx_cell] = var_data->getPointer(di)[idx_cell_pivot];
                            }
                        }
                    }
                }
            }
        }
    }
}


/*
 * Function to remove 3d boundary nodes with boundary conditions filled by this class for a patch.
 *
 * Arguments are:
 *    bdry_node_locs ....... array of locations of nodes for applying
 *                           boundary conditions.
 *    patch ................ patch on which data object lives
 *    bdry_node_conds ...... array of boundary conditions for patch nodes
 */
void
BasicCartesianBoundaryUtilities3::removeBoundaryNodeLocations(
    std::vector<int>& bdry_node_locs,
    const hier::Patch& patch,
    const std::vector<int>& bdry_node_conds)
{
    TBOX_ASSERT(static_cast<int>(bdry_node_locs.size()) <= NUM_3D_NODES);
    TBOX_ASSERT(*min_element(bdry_node_locs.begin(), bdry_node_locs.end()) >= 0);
    TBOX_ASSERT(*max_element(bdry_node_locs.begin(), bdry_node_locs.end()) < NUM_3D_NODES);
    TBOX_ASSERT(static_cast<int>(bdry_node_conds.size()) == NUM_3D_NODES);
    
    const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
        HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    TBOX_ASSERT(patch_geom);
    
    const std::vector<hier::BoundaryBox>& node_bdry =
        patch_geom->getCodimensionBoundaries(BDRY::NODE3D);
    
    for (int ni = 0; ni < static_cast<int>(node_bdry.size()); ni++)
    {
        TBOX_ASSERT(node_bdry[ni].getBoundaryType() == BDRY::NODE3D);
        
        int node_loc(node_bdry[ni].getLocationIndex());
        
        if (std::find(bdry_node_locs.begin(), bdry_node_locs.end(), node_loc) !=
            bdry_node_locs.end())
        {
            if (bdry_node_conds[node_loc] == BDRY_COND::BASIC::XDIRICHLET)
            {
                // Remove node locations that have boundary conditions identified.
                bdry_node_locs.erase(std::remove(bdry_node_locs.begin(), bdry_node_locs.end(), node_loc),
                    bdry_node_locs.end());
            }
            else if (bdry_node_conds[node_loc] == BDRY_COND::BASIC::YDIRICHLET)
            {
                // Remove node locations that have boundary conditions identified.
                bdry_node_locs.erase(std::remove(bdry_node_locs.begin(), bdry_node_locs.end(), node_loc),
                    bdry_node_locs.end());
            }
            else if (bdry_node_conds[node_loc] == BDRY_COND::BASIC::ZDIRICHLET)
            {
                // Remove node locations that have boundary conditions identified.
                bdry_node_locs.erase(std::remove(bdry_node_locs.begin(), bdry_node_locs.end(), node_loc),
                    bdry_node_locs.end());
            }
            else if (bdry_node_conds[node_loc] == BDRY_COND::BASIC::XNEUMANN)
            {
                // NOT YET IMPLEMENTED
            }
            else if (bdry_node_conds[node_loc] == BDRY_COND::BASIC::YNEUMANN)
            {
                // NOT YET IMPLEMENTED
            }
            else if (bdry_node_conds[node_loc] == BDRY_COND::BASIC::ZNEUMANN)
            {
                // NOT YET IMPLEMENTED
            }
            else if (bdry_node_conds[node_loc] == BDRY_COND::BASIC::XFLOW)
            {
                // Remove node locations that have boundary conditions identified.
                bdry_node_locs.erase(std::remove(bdry_node_locs.begin(), bdry_node_locs.end(), node_loc),
                    bdry_node_locs.end());
            }
            else if (bdry_node_conds[node_loc] == BDRY_COND::BASIC::YFLOW)
            {
                // Remove node locations that have boundary conditions identified.
                bdry_node_locs.erase(std::remove(bdry_node_locs.begin(), bdry_node_locs.end(), node_loc),
                    bdry_node_locs.end());
            }
            else if (bdry_node_conds[node_loc] == BDRY_COND::BASIC::ZFLOW)
            {
                // Remove node locations that have boundary conditions identified.
                bdry_node_locs.erase(std::remove(bdry_node_locs.begin(), bdry_node_locs.end(), node_loc),
                    bdry_node_locs.end());
            }
            else if (bdry_node_conds[node_loc] == BDRY_COND::BASIC::XREFLECT)
            {
                // Remove node locations that have boundary conditions identified.
                bdry_node_locs.erase(std::remove(bdry_node_locs.begin(), bdry_node_locs.end(), node_loc),
                    bdry_node_locs.end());
            }
            else if (bdry_node_conds[node_loc] == BDRY_COND::BASIC::YREFLECT)
            {
                // Remove node locations that have boundary conditions identified.
                bdry_node_locs.erase(std::remove(bdry_node_locs.begin(), bdry_node_locs.end(), node_loc),
                    bdry_node_locs.end());
            }
            else if (bdry_node_conds[node_loc] == BDRY_COND::BASIC::ZREFLECT)
            {
                // Remove node locations that have boundary conditions identified.
                bdry_node_locs.erase(std::remove(bdry_node_locs.begin(), bdry_node_locs.end(), node_loc),
                    bdry_node_locs.end());
            }
            else if (bdry_node_conds[node_loc] == BDRY_COND::BASIC::XSYMMETRY)
            {
                // Remove node locations that have boundary conditions identified.
                bdry_node_locs.erase(std::remove(bdry_node_locs.begin(), bdry_node_locs.end(), node_loc),
                    bdry_node_locs.end());
            }
            else if (bdry_node_conds[node_loc] == BDRY_COND::BASIC::YSYMMETRY)
            {
                // Remove node locations that have boundary conditions identified.
                bdry_node_locs.erase(std::remove(bdry_node_locs.begin(), bdry_node_locs.end(), node_loc),
                    bdry_node_locs.end());
            }
            else if (bdry_node_conds[node_loc] == BDRY_COND::BASIC::ZSYMMETRY)
            {
                // Remove node locations that have boundary conditions identified.
                bdry_node_locs.erase(std::remove(bdry_node_locs.begin(), bdry_node_locs.end(), node_loc),
                    bdry_node_locs.end());
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
 *    bdry_face_values ..... array of boundary values for faces
 *                           (this must be consistent with boundary
 *                           condition types)
 *    ghost_width_to_fill .. width of ghost region to fill
 */
void
BasicCartesianBoundaryUtilities3::fillNodeBoundaryData(
    const std::string& var_name,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& var_data,
    const hier::Patch& patch,
    const std::vector<int>& bdry_node_locs,
    const std::vector<int>& bdry_node_conds,
    const std::vector<double>& bdry_face_values,
    const hier::IntVector& ghost_width_to_fill)
{
    TBOX_ASSERT(!var_name.empty());
    TBOX_ASSERT(var_data);
    TBOX_ASSERT(static_cast<int>(bdry_node_locs.size()) <= NUM_3D_NODES);
    TBOX_ASSERT(*min_element(bdry_node_locs.begin(), bdry_node_locs.end()) >= 0);
    TBOX_ASSERT(*max_element(bdry_node_locs.begin(), bdry_node_locs.end()) < NUM_3D_NODES);
    TBOX_ASSERT(static_cast<int>(bdry_node_conds.size()) == NUM_3D_NODES);
    TBOX_ASSERT(static_cast<int>(bdry_face_values.size()) == NUM_3D_FACES*(var_data->getDepth()));
    
    TBOX_DIM_ASSERT(ghost_width_to_fill.getDim() == tbox::Dimension(3));
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
    hier::IntVector gcw_to_fill(tbox::Dimension(3));
    
    // If the ghost fill width is not used, it is set to the ghost cell width of the data.
    if (ghost_width_to_fill == -hier::IntVector::getOne(tbox::Dimension(3)))
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
        patch_geom->getCodimensionBoundaries(BDRY::NODE3D);
    
    const int var_depth = var_data->getDepth();
    
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
            
            if (bdry_node_conds[node_loc] == BDRY_COND::BASIC::XDIRICHLET)
            {
                for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                {
                    for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                    {
                        for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                        {
                            const int idx_cell = (i + num_ghosts[0]) +
                                (j + num_ghosts[1])*ghostcell_dims[0] +
                                (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                
                            for (int di = 0; di < var_depth; di++)
                            {
                                var_data->getPointer(di)[idx_cell] = bdry_face_values[face_loc_0*var_depth + di];
                            }
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
                        for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                        {
                            const int idx_cell = (i + num_ghosts[0]) +
                                (j + num_ghosts[1])*ghostcell_dims[0] +
                                (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                
                            for (int di = 0; di < var_depth; di++)
                            {
                                var_data->getPointer(di)[idx_cell] = bdry_face_values[face_loc_1*var_depth + di];
                            }
                        }
                    }
                }
            }
            else if (bdry_node_conds[node_loc] == BDRY_COND::BASIC::ZDIRICHLET)
            {
                for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                {
                    for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                    {
                        for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                        {
                            const int idx_cell = (i + num_ghosts[0]) +
                                (j + num_ghosts[1])*ghostcell_dims[0] +
                                (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                
                            for (int di = 0; di < var_depth; di++)
                            {
                                var_data->getPointer(di)[idx_cell] = bdry_face_values[face_loc_2*var_depth + di];
                            }
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
            else if (bdry_node_conds[node_loc] == BDRY_COND::BASIC::ZNEUMANN)
            {
                // NOT YET IMPLEMENTED
            }
            else if (bdry_node_conds[node_loc] == BDRY_COND::BASIC::XFLOW)
            {
                for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                {
                    for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                    {
                        for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                        {
                            const int idx_cell = (i + num_ghosts[0]) +
                                (j + num_ghosts[1])*ghostcell_dims[0] +
                                (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            
                            int idx_cell_pivot = idx_cell;
                            
                            if (face_loc_0 == BDRY_LOC::XLO)
                            {
                                idx_cell_pivot = (interior_box_lo_idx[0] + num_ghosts[0]) +
                                    (j + num_ghosts[1])*ghostcell_dims[0] +
                                    (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            }
                            else if (face_loc_0 == BDRY_LOC::XHI)
                            {
                                idx_cell_pivot = (interior_box_hi_idx[0] + num_ghosts[0]) +
                                    (j + num_ghosts[1])*ghostcell_dims[0] +
                                    (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            }
                            
                            for (int di = 0; di < var_depth; di++)
                            {
                                var_data->getPointer(di)[idx_cell] = var_data->getPointer(di)[idx_cell_pivot];
                            }
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
                        for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                        {
                            const int idx_cell = (i + num_ghosts[0]) +
                                (j + num_ghosts[1])*ghostcell_dims[0] +
                                (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            
                            int idx_cell_pivot = idx_cell;
                            
                            if (face_loc_1 == BDRY_LOC::YLO)
                            {
                                idx_cell_pivot = (i + num_ghosts[0]) +
                                    (interior_box_lo_idx[1] + num_ghosts[1])*ghostcell_dims[0] +
                                    (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            }
                            else if (face_loc_1 == BDRY_LOC::YHI)
                            {
                                idx_cell_pivot = (i + num_ghosts[0]) +
                                    (interior_box_hi_idx[1] + num_ghosts[1])*ghostcell_dims[0] +
                                    (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            }
                            
                            for (int di = 0; di < var_depth; di++)
                            {
                                var_data->getPointer(di)[idx_cell] = var_data->getPointer(di)[idx_cell_pivot];
                            }
                        }
                    }
                }
            }
            else if (bdry_node_conds[node_loc] == BDRY_COND::BASIC::ZFLOW)
            {
                for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                {
                    for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                    {
                        for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                        {
                            const int idx_cell = (i + num_ghosts[0]) +
                                (j + num_ghosts[1])*ghostcell_dims[0] +
                                (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            
                            int idx_cell_pivot = idx_cell;
                            
                            if (face_loc_2 == BDRY_LOC::ZLO)
                            {
                                idx_cell_pivot = (i + num_ghosts[0]) +
                                    (j + num_ghosts[1])*ghostcell_dims[0] +
                                    (interior_box_lo_idx[2] + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            }
                            else if (face_loc_2 == BDRY_LOC::ZHI)
                            {
                                idx_cell_pivot = (i + num_ghosts[0]) +
                                    (j + num_ghosts[1])*ghostcell_dims[0] +
                                    (interior_box_hi_idx[2] + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            }
                            
                            for (int di = 0; di < var_depth; di++)
                            {
                                var_data->getPointer(di)[idx_cell] = var_data->getPointer(di)[idx_cell_pivot];
                            }
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
                        for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                        {
                            const int idx_cell = (i + num_ghosts[0]) +
                                (j + num_ghosts[1])*ghostcell_dims[0] +
                                (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            
                            int idx_cell_pivot = idx_cell;
                            
                            if (face_loc_0 == BDRY_LOC::XLO)
                            {
                                idx_cell_pivot = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) + num_ghosts[0]) +
                                    (j + num_ghosts[1])*ghostcell_dims[0] +
                                    (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            }
                            else if (face_loc_0 == BDRY_LOC::XHI)
                            {
                                idx_cell_pivot = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) + num_ghosts[0]) +
                                    (j + num_ghosts[1])*ghostcell_dims[0] +
                                    (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            }
                            
                            for (int di = 0; di < var_depth; di++)
                            {
                                var_data->getPointer(di)[idx_cell] = var_data->getPointer(di)[idx_cell_pivot];
                            }
                            
                            var_data->getPointer(0)[idx_cell] = -var_data->getPointer(0)[idx_cell_pivot];
                        }
                    }
                }
            }
            else if (bdry_node_conds[node_loc] == BDRY_COND::BASIC::YREFLECT)
            {
                for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                {
                    for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                    {
                        for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                        {
                            const int idx_cell = (i + num_ghosts[0]) +
                                (j + num_ghosts[1])*ghostcell_dims[0] +
                                (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            
                            int idx_cell_pivot = idx_cell;
                            
                            if (face_loc_1 == BDRY_LOC::YLO)
                            {
                                idx_cell_pivot = (i + num_ghosts[0]) +
                                    (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) + num_ghosts[1])*ghostcell_dims[0] +
                                    (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            }
                            else if (face_loc_1 == BDRY_LOC::YHI)
                            {
                                idx_cell_pivot = (i + num_ghosts[0]) +
                                    (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) + num_ghosts[1])*ghostcell_dims[0] +
                                    (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            }
                            
                            for (int di = 0; di < var_depth; di++)
                            {
                                var_data->getPointer(di)[idx_cell] = var_data->getPointer(di)[idx_cell_pivot];
                            }
                            
                            var_data->getPointer(1)[idx_cell] = -var_data->getPointer(1)[idx_cell_pivot];
                        }
                    }
                }
            }
            else if (bdry_node_conds[node_loc] == BDRY_COND::BASIC::ZREFLECT)
            {
                for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                {
                    for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                    {
                        for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                        {
                            const int idx_cell = (i + num_ghosts[0]) +
                                (j + num_ghosts[1])*ghostcell_dims[0] +
                                (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            
                            int idx_cell_pivot = idx_cell;
                            
                            if (face_loc_2 == BDRY_LOC::ZLO)
                            {
                                idx_cell_pivot = (i + num_ghosts[0]) +
                                    (j + num_ghosts[1])*ghostcell_dims[0] +
                                    (interior_box_lo_idx[2] + (fill_box_hi_idx[2] - k) + num_ghosts[2])
                                        *ghostcell_dims[0]*ghostcell_dims[1];
                            }
                            else if (face_loc_2 == BDRY_LOC::ZHI)
                            {
                                idx_cell_pivot = (i + num_ghosts[0]) +
                                    (j + num_ghosts[1])*ghostcell_dims[0] +
                                    (interior_box_hi_idx[2] - (k - fill_box_lo_idx[2]) + num_ghosts[2])
                                        *ghostcell_dims[0]*ghostcell_dims[1];
                            }
                            
                            for (int di = 0; di < var_depth; di++)
                            {
                                var_data->getPointer(di)[idx_cell] = var_data->getPointer(di)[idx_cell_pivot];
                            }
                            
                            var_data->getPointer(2)[idx_cell] = -var_data->getPointer(2)[idx_cell_pivot];
                        }
                    }
                }
            }
            else if (bdry_node_conds[node_loc] == BDRY_COND::BASIC::XSYMMETRY)
            {
                for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                {
                    for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                    {
                        for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                        {
                            const int idx_cell = (i + num_ghosts[0]) +
                                (j + num_ghosts[1])*ghostcell_dims[0] +
                                (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            
                            int idx_cell_pivot = idx_cell;
                            
                            if (face_loc_0 == BDRY_LOC::XLO)
                            {
                                idx_cell_pivot = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) + num_ghosts[0]) +
                                    (j + num_ghosts[1])*ghostcell_dims[0] +
                                    (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            }
                            else if (face_loc_0 == BDRY_LOC::XHI)
                            {
                                idx_cell_pivot = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) + num_ghosts[0]) +
                                    (j + num_ghosts[1])*ghostcell_dims[0] +
                                    (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            }
                            
                            for (int di = 0; di < var_depth; di++)
                            {
                                var_data->getPointer(di)[idx_cell] = var_data->getPointer(di)[idx_cell_pivot];
                            }
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
                        for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                        {
                            const int idx_cell = (i + num_ghosts[0]) +
                                (j + num_ghosts[1])*ghostcell_dims[0] +
                                (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            
                            int idx_cell_pivot = idx_cell;
                            
                            if (face_loc_1 == BDRY_LOC::YLO)
                            {
                                idx_cell_pivot = (i + num_ghosts[0]) +
                                    (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) + num_ghosts[1])*ghostcell_dims[0] +
                                    (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            }
                            else if (face_loc_1 == BDRY_LOC::YHI)
                            {
                                idx_cell_pivot = (i + num_ghosts[0]) +
                                    (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) + num_ghosts[1])*ghostcell_dims[0] +
                                    (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            }
                            
                            for (int di = 0; di < var_depth; di++)
                            {
                                var_data->getPointer(di)[idx_cell] = var_data->getPointer(di)[idx_cell_pivot];
                            }
                        }
                    }
                }
            }
            else if (bdry_node_conds[node_loc] == BDRY_COND::BASIC::ZSYMMETRY)
            {
                for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                {
                    for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                    {
                        for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                        {
                            const int idx_cell = (i + num_ghosts[0]) +
                                (j + num_ghosts[1])*ghostcell_dims[0] +
                                (k + num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            
                            int idx_cell_pivot = idx_cell;
                            
                            if (face_loc_2 == BDRY_LOC::ZLO)
                            {
                                idx_cell_pivot = (i + num_ghosts[0]) +
                                    (j + num_ghosts[1])*ghostcell_dims[0] +
                                    (interior_box_lo_idx[2] + (fill_box_hi_idx[2] - k) + num_ghosts[2])
                                        *ghostcell_dims[0]*ghostcell_dims[1];
                            }
                            else if (face_loc_2 == BDRY_LOC::ZHI)
                            {
                                idx_cell_pivot = (i + num_ghosts[0]) +
                                    (j + num_ghosts[1])*ghostcell_dims[0] +
                                    (interior_box_hi_idx[2] - (k - fill_box_lo_idx[2]) + num_ghosts[2])
                                        *ghostcell_dims[0]*ghostcell_dims[1];
                            }
                            
                            for (int di = 0; di < var_depth; di++)
                            {
                                var_data->getPointer(di)[idx_cell] = var_data->getPointer(di)[idx_cell_pivot];
                            }
                        }
                    }
                }
            }
            // else
            // {
            //     TBOX_ERROR("BasicCartesianBoundaryUtilities3::fillNodeBoundaryData()\n"
            //         << "Invalid node boundary condition!\n"
            //         << "node_loc = '" << node_loc << "'." << std::endl
            //         << "bdry_node_conds[node_loc] = '" << bdry_node_conds[node_loc] << "'."
            //         << std::endl);
            // }
        }
    }
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
BasicCartesianBoundaryUtilities3::getFaceLocationForEdgeBdry(
    int edge_loc,
    int edge_btype)
{
    int ret_face = -1;
    
    switch (edge_btype)
    {
        case BDRY_COND::BASIC::XFLOW:
        case BDRY_COND::BASIC::XREFLECT:
        case BDRY_COND::BASIC::XSYMMETRY:
        case BDRY_COND::BASIC::XDIRICHLET:
        case BDRY_COND::BASIC::XNEUMANN:
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
        case BDRY_COND::BASIC::YFLOW:
        case BDRY_COND::BASIC::YREFLECT:
        case BDRY_COND::BASIC::YSYMMETRY:
        case BDRY_COND::BASIC::YDIRICHLET:
        case BDRY_COND::BASIC::YNEUMANN:
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
        case BDRY_COND::BASIC::ZFLOW:
        case BDRY_COND::BASIC::ZREFLECT:
        case BDRY_COND::BASIC::ZSYMMETRY:
        case BDRY_COND::BASIC::ZDIRICHLET:
        case BDRY_COND::BASIC::ZNEUMANN:
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
            // TBOX_ERROR("BasicCartesianBoundaryUtilities3::getFaceLocationForEdgeBdry()\n"
            //     << "Unknown edge boundary condition type = '"
            //     << edge_btype
            //     << "' passed."
            //     << std::endl);
        }
    }
    
    // if (ret_face == -1)
    // {
    //     TBOX_ERROR("BasicCartesianBoundaryUtilities3::getFaceLocationForEdgeBdry()\n"
    //         << "Edge boundary condition type = '"
    //         << edge_btype << "' and \n"
    //         << "edge location = '" << edge_loc
    //         << "' passed are inconsistent."
    //         << std::endl);
    // }
    
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
BasicCartesianBoundaryUtilities3::getFaceLocationForNodeBdry(
    int node_loc,
    int node_btype)
{
    int ret_face = -1;
    
    switch (node_btype)
    {
        case BDRY_COND::BASIC::XFLOW:
        case BDRY_COND::BASIC::XREFLECT:
        case BDRY_COND::BASIC::XSYMMETRY:
        case BDRY_COND::BASIC::XDIRICHLET:
        case BDRY_COND::BASIC::XNEUMANN:
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
        case BDRY_COND::BASIC::YFLOW:
        case BDRY_COND::BASIC::YREFLECT:
        case BDRY_COND::BASIC::YSYMMETRY:
        case BDRY_COND::BASIC::YDIRICHLET:
        case BDRY_COND::BASIC::YNEUMANN:
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
        case BDRY_COND::BASIC::ZFLOW:
        case BDRY_COND::BASIC::ZREFLECT:
        case BDRY_COND::BASIC::ZSYMMETRY:
        case BDRY_COND::BASIC::ZDIRICHLET:
        case BDRY_COND::BASIC::ZNEUMANN:
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
            // TBOX_ERROR("BasicCartesianBoundaryUtilities3::getFaceLocationForNodeBdry()\n"
            //     << "Unknown node boundary condition type = '"
            //     << node_btype << "' passed."
            //     << std::endl);
        }
    }
    
    // if (ret_face == -1)
    // {
    //     TBOX_ERROR("BasicCartesianBoundaryUtilities3::getFaceLocationForNodeBdry()\n"
    //         << "Node boundary condition type = '"
    //         << node_btype << "' and \n"
    //         << "node location = '" << node_loc
    //         << "' passed are inconsistent." << std::endl);
    // }
    
    return ret_face;
}


/*
 * Private function to read 3D face boundary data from input database.
 */
void
BasicCartesianBoundaryUtilities3::read3dBdryFaces(
    BoundaryUtilityStrategy* bdry_strategy,
    const HAMERS_SHARED_PTR<tbox::Database>& input_db,
    std::vector<int>& face_locs,
    std::vector<int>& face_conds,
    const hier::IntVector& periodic)
{
    TBOX_DIM_ASSERT(periodic.getDim() == tbox::Dimension(3));
    
    TBOX_ASSERT(bdry_strategy != 0);
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
        // face boundary input required
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
                if (bdry_cond_str == "FLOW")
                {
                    face_conds[s] = BDRY_COND::BASIC::FLOW;
                    
                    face_locs[fi] = BOGUS_BDRY_LOC;
                }
                else if (bdry_cond_str == "REFLECT")
                {
                    face_conds[s] = BDRY_COND::BASIC::REFLECT;
                    
                    face_locs[fi] = BOGUS_BDRY_LOC;
                }
                else if (bdry_cond_str == "SYMMETRY")
                {
                    face_conds[s] = BDRY_COND::BASIC::SYMMETRY;
                    
                    face_locs[fi] = BOGUS_BDRY_LOC;
                }
                else if (bdry_cond_str == "DIRICHLET")
                {
                    face_conds[s] = BDRY_COND::BASIC::DIRICHLET;
                    bdry_strategy->readDirichletBoundaryDataEntry(
                        bdry_loc_db,
                        bdry_loc_str,
                        s);
                    
                    face_locs[fi] = BOGUS_BDRY_LOC;
                }
                else if (bdry_cond_str == "NEUMANN")
                {
                    face_conds[s] = BDRY_COND::BASIC::NEUMANN;
                    bdry_strategy->readNeumannBoundaryDataEntry(
                        bdry_loc_db,
                        bdry_loc_str,
                        s);
                    
                    face_locs[fi] = BOGUS_BDRY_LOC;
                }
            } // if (need_data_read)
        } // for (int fi = 0 ...
    } // if (num_per_dirs < 3)
    
    // Remove face locations that have boundary conditions identified.
    face_locs.erase(std::remove(face_locs.begin(), face_locs.end(), BOGUS_BDRY_LOC), face_locs.end());
}


/*
 * Private function to read 3D edge boundary data from input database.
 */
void
BasicCartesianBoundaryUtilities3::read3dBdryEdges(
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
                if (bdry_cond_str == "XFLOW")
                {
                    edge_conds[s] = BDRY_COND::BASIC::XFLOW;
                    edge_locs[ei] = BOGUS_BDRY_LOC;
                }
                else if (bdry_cond_str == "YFLOW")
                {
                    edge_conds[s] = BDRY_COND::BASIC::YFLOW;
                    edge_locs[ei] = BOGUS_BDRY_LOC;
                }
                else if (bdry_cond_str == "ZFLOW")
                {
                    edge_conds[s] = BDRY_COND::BASIC::ZFLOW;
                    edge_locs[ei] = BOGUS_BDRY_LOC;
                }
                else if (bdry_cond_str == "XREFLECT")
                {
                    edge_conds[s] = BDRY_COND::BASIC::XREFLECT;
                    edge_locs[ei] = BOGUS_BDRY_LOC;
                }
                else if (bdry_cond_str == "YREFLECT")
                {
                    edge_conds[s] = BDRY_COND::BASIC::YREFLECT;
                    edge_locs[ei] = BOGUS_BDRY_LOC;
                }
                else if (bdry_cond_str == "ZREFLECT")
                {
                    edge_conds[s] = BDRY_COND::BASIC::ZREFLECT;
                    edge_locs[ei] = BOGUS_BDRY_LOC;
                }
                else if (bdry_cond_str == "XSYMMETRY")
                {
                    edge_conds[s] = BDRY_COND::BASIC::XSYMMETRY;
                    edge_locs[ei] = BOGUS_BDRY_LOC;
                }
                else if (bdry_cond_str == "YSYMMETRY")
                {
                    edge_conds[s] = BDRY_COND::BASIC::YSYMMETRY;
                    edge_locs[ei] = BOGUS_BDRY_LOC;
                }
                else if (bdry_cond_str == "ZSYMMETRY")
                {
                    edge_conds[s] = BDRY_COND::BASIC::ZSYMMETRY;
                    edge_locs[ei] = BOGUS_BDRY_LOC;
                }
                else if (bdry_cond_str == "XDIRICHLET")
                {
                    edge_conds[s] = BDRY_COND::BASIC::XDIRICHLET;
                    edge_locs[ei] = BOGUS_BDRY_LOC;
                }
                else if (bdry_cond_str == "YDIRICHLET")
                {
                    edge_conds[s] = BDRY_COND::BASIC::YDIRICHLET;
                    edge_locs[ei] = BOGUS_BDRY_LOC;
                }
                else if (bdry_cond_str == "ZDIRICHLET")
                {
                    edge_conds[s] = BDRY_COND::BASIC::ZDIRICHLET;
                    edge_locs[ei] = BOGUS_BDRY_LOC;
                }
                else if (bdry_cond_str == "XNEUMANN")
                {
                    edge_conds[s] = BDRY_COND::BASIC::XNEUMANN;
                    edge_locs[ei] = BOGUS_BDRY_LOC;
                }
                else if (bdry_cond_str == "YNEUMANN")
                {
                    edge_conds[s] = BDRY_COND::BASIC::YNEUMANN;
                    edge_locs[ei] = BOGUS_BDRY_LOC;
                }
                else if (bdry_cond_str == "ZNEUMANN")
                {
                    edge_conds[s] = BDRY_COND::BASIC::ZNEUMANN;
                    edge_locs[ei] = BOGUS_BDRY_LOC;
                }
                
                bool ambiguous_type = false;
                if (bdry_cond_str == "XFLOW" ||
                    bdry_cond_str == "XREFLECT" ||
                    bdry_cond_str == "XSYMMETRY" ||
                    bdry_cond_str == "XDIRICHLET" ||
                    bdry_cond_str == "XNEUMANN")
                {
                    if (s == EDGE_BDRY_LOC_3D::YLO_ZLO ||
                        s == EDGE_BDRY_LOC_3D::YHI_ZLO ||
                        s == EDGE_BDRY_LOC_3D::YLO_ZHI ||
                        s == EDGE_BDRY_LOC_3D::YHI_ZHI)
                    {
                        ambiguous_type = true;
                    }
                }
                else if (bdry_cond_str == "YFLOW" ||
                         bdry_cond_str == "YREFLECT" ||
                         bdry_cond_str == "YSYMMETRY" ||
                         bdry_cond_str == "YDIRICHLET" ||
                         bdry_cond_str == "YNEUMANN")
                {
                    if (s == EDGE_BDRY_LOC_3D::XLO_ZLO ||
                        s == EDGE_BDRY_LOC_3D::XLO_ZHI ||
                        s == EDGE_BDRY_LOC_3D::XHI_ZLO ||
                        s == EDGE_BDRY_LOC_3D::XHI_ZHI)
                    {
                        ambiguous_type = true;
                    }
                }
                else if (bdry_cond_str == "ZFLOW" ||
                         bdry_cond_str == "ZREFLECT" ||
                         bdry_cond_str == "ZSYMMETRY" ||
                         bdry_cond_str == "ZDIRICHLET" ||
                         bdry_cond_str == "ZNEUMANN")
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
                    TBOX_ERROR("BasicCartesianBoundaryUtilities3::read3dBdryEdges()\n"
                        << "Ambiguous bdry condition '"
                        << bdry_cond_str
                        << "' found for '"
                        << bdry_loc_str
                        << "'."
                        << std::endl);
                }
                
                std::string proper_face;
                std::string proper_face_data;
                bool no_face_data_found = false;
                if (bdry_cond_str == "XFLOW" ||
                    bdry_cond_str == "XDIRICHLET" ||
                    bdry_cond_str == "XNEUMANN" ||
                    bdry_cond_str == "XREFLECT" ||
                    bdry_cond_str == "XSYMMETRY")
                {
                    if (s == EDGE_BDRY_LOC_3D::XLO_ZLO ||
                        s == EDGE_BDRY_LOC_3D::XLO_ZHI ||
                        s == EDGE_BDRY_LOC_3D::XLO_YLO ||
                        s == EDGE_BDRY_LOC_3D::XLO_YHI)
                    {
                        proper_face = "XLO";
                        if (bdry_cond_str == "XFLOW" &&
                            face_conds[BDRY_LOC::XLO] != BDRY_COND::BASIC::FLOW)
                        {
                            no_face_data_found = true;
                            proper_face_data = "FLOW";
                        }
                        if (bdry_cond_str == "XDIRICHLET" &&
                            face_conds[BDRY_LOC::XLO] != BDRY_COND::BASIC::DIRICHLET)
                        {
                            no_face_data_found = true;
                            proper_face_data = "DIRICHLET";
                        }
                        if (bdry_cond_str == "XNEUMANN" &&
                            face_conds[BDRY_LOC::XLO] != BDRY_COND::BASIC::NEUMANN)
                        {
                            no_face_data_found = true;
                            proper_face_data = "NEUMANN";
                        }
                        if (bdry_cond_str == "XREFLECT" &&
                            face_conds[BDRY_LOC::XLO] != BDRY_COND::BASIC::REFLECT)
                        {
                            no_face_data_found = true;
                            proper_face_data = "REFLECT";
                        }
                        if (bdry_cond_str == "XSYMMETRY" &&
                            face_conds[BDRY_LOC::XLO] != BDRY_COND::BASIC::SYMMETRY)
                        {
                            no_face_data_found = true;
                            proper_face_data = "SYMMETRY";
                        }
                    }
                    else
                    {
                        proper_face = "XHI";
                        if (bdry_cond_str == "XFLOW" &&
                            face_conds[BDRY_LOC::XHI] != BDRY_COND::BASIC::FLOW)
                        {
                            no_face_data_found = true;
                            proper_face_data = "FLOW";
                        }
                        if (bdry_cond_str == "XDIRICHLET" &&
                            face_conds[BDRY_LOC::XHI] != BDRY_COND::BASIC::DIRICHLET)
                        {
                            no_face_data_found = true;
                            proper_face_data = "DIRICHLET";
                        }
                        if (bdry_cond_str == "XNEUMANN" &&
                            face_conds[BDRY_LOC::XHI] != BDRY_COND::BASIC::NEUMANN)
                        {
                            no_face_data_found = true;
                            proper_face_data = "NEUMANN";
                        }
                        if (bdry_cond_str == "XREFLECT" &&
                            face_conds[BDRY_LOC::XHI] != BDRY_COND::BASIC::REFLECT)
                        {
                            no_face_data_found = true;
                            proper_face_data = "REFLECT";
                        }
                        if (bdry_cond_str == "XSYMMETRY" &&
                            face_conds[BDRY_LOC::XHI] != BDRY_COND::BASIC::SYMMETRY)
                        {
                            no_face_data_found = true;
                            proper_face_data = "SYMMETRY";
                        }
                    }
                }
                else if (bdry_cond_str == "YFLOW" ||
                         bdry_cond_str == "YDIRICHLET" ||
                         bdry_cond_str == "YNEUMANN" ||
                         bdry_cond_str == "YREFLECT" ||
                         bdry_cond_str == "YSYMMETRY")
                {
                    if (s == EDGE_BDRY_LOC_3D::YLO_ZLO ||
                        s == EDGE_BDRY_LOC_3D::YLO_ZHI ||
                        s == EDGE_BDRY_LOC_3D::XLO_YLO ||
                        s == EDGE_BDRY_LOC_3D::XHI_YLO)
                    {
                        proper_face = "YLO";
                        if (bdry_cond_str == "YFLOW" &&
                            face_conds[BDRY_LOC::YLO] != BDRY_COND::BASIC::FLOW)
                        {
                            no_face_data_found = true;
                            proper_face_data = "FLOW";
                        }
                        if (bdry_cond_str == "YDIRICHLET" &&
                            face_conds[BDRY_LOC::YLO] != BDRY_COND::BASIC::DIRICHLET)
                        {
                            no_face_data_found = true;
                            proper_face_data = "DIRICHLET";
                        }
                        if (bdry_cond_str == "YNEUMANN" &&
                            face_conds[BDRY_LOC::YLO] != BDRY_COND::BASIC::NEUMANN)
                        {
                            no_face_data_found = true;
                            proper_face_data = "NEUMANN";
                        }
                        if (bdry_cond_str == "YREFLECT" &&
                            face_conds[BDRY_LOC::YLO] != BDRY_COND::BASIC::REFLECT)
                        {
                            no_face_data_found = true;
                            proper_face_data = "REFLECT";
                        }
                        if (bdry_cond_str == "YSYMMETRY" &&
                            face_conds[BDRY_LOC::YLO] != BDRY_COND::BASIC::SYMMETRY)
                        {
                            no_face_data_found = true;
                            proper_face_data = "SYMMETRY";
                        }
                    }
                    else
                    {
                        proper_face = "YHI";
                        if (bdry_cond_str == "YFLOW" &&
                            face_conds[BDRY_LOC::YHI] != BDRY_COND::BASIC::FLOW)
                        {
                            no_face_data_found = true;
                            proper_face_data = "FLOW";
                        }
                        if (bdry_cond_str == "YDIRICHLET" &&
                            face_conds[BDRY_LOC::YHI] != BDRY_COND::BASIC::DIRICHLET)
                        {
                            no_face_data_found = true;
                            proper_face_data = "DIRICHLET";
                        }
                        if (bdry_cond_str == "YNEUMANN" &&
                            face_conds[BDRY_LOC::YHI] != BDRY_COND::BASIC::NEUMANN)
                        {
                            no_face_data_found = true;
                            proper_face_data = "NEUMANN";
                        }
                        if (bdry_cond_str == "YREFLECT" &&
                            face_conds[BDRY_LOC::YHI] != BDRY_COND::BASIC::REFLECT)
                        {
                            no_face_data_found = true;
                            proper_face_data = "REFLECT";
                        }
                        if (bdry_cond_str == "YSYMMETRY" &&
                            face_conds[BDRY_LOC::YHI] != BDRY_COND::BASIC::SYMMETRY)
                        {
                            no_face_data_found = true;
                            proper_face_data = "SYMMETRY";
                        }
                    }
                }
                else if (bdry_cond_str == "ZFLOW" ||
                         bdry_cond_str == "ZDIRICHLET" ||
                         bdry_cond_str == "ZNEUMANN" ||
                         bdry_cond_str == "ZREFLECT" ||
                         bdry_cond_str == "ZSYMMETRY")
                {
                    if (s == EDGE_BDRY_LOC_3D::XLO_ZLO ||
                        s == EDGE_BDRY_LOC_3D::YHI_ZLO ||
                        s == EDGE_BDRY_LOC_3D::YLO_ZLO ||
                        s == EDGE_BDRY_LOC_3D::XHI_ZLO)
                    {
                        proper_face = "ZLO";
                        if (bdry_cond_str == "ZFLOW" &&
                            face_conds[BDRY_LOC::ZLO] != BDRY_COND::BASIC::FLOW)
                        {
                            no_face_data_found = true;
                            proper_face_data = "FLOW";
                        }
                        if (bdry_cond_str == "ZDIRICHLET" &&
                            face_conds[BDRY_LOC::ZLO] != BDRY_COND::BASIC::DIRICHLET)
                        {
                            no_face_data_found = true;
                            proper_face_data = "DIRICHLET";
                        }
                        if (bdry_cond_str == "ZNEUMANN" &&
                            face_conds[BDRY_LOC::ZLO] != BDRY_COND::BASIC::NEUMANN)
                        {
                            no_face_data_found = true;
                            proper_face_data = "NEUMANN";
                        }
                        if (bdry_cond_str == "ZREFLECT" &&
                            face_conds[BDRY_LOC::ZLO] != BDRY_COND::BASIC::REFLECT)
                        {
                            no_face_data_found = true;
                            proper_face_data = "REFLECT";
                        }
                        if (bdry_cond_str == "ZSYMMETRY" &&
                            face_conds[BDRY_LOC::ZLO] != BDRY_COND::BASIC::SYMMETRY)
                        {
                            no_face_data_found = true;
                            proper_face_data = "SYMMETRY";
                        }
                    }
                    else
                    {
                        proper_face = "ZHI";
                        if (bdry_cond_str == "ZFLOW" &&
                            face_conds[BDRY_LOC::ZHI] != BDRY_COND::BASIC::FLOW)
                        {
                            no_face_data_found = true;
                            proper_face_data = "FLOW";
                        }
                        if (bdry_cond_str == "ZDIRICHLET" &&
                            face_conds[BDRY_LOC::ZHI] != BDRY_COND::BASIC::DIRICHLET)
                        {
                            no_face_data_found = true;
                            proper_face_data = "DIRICHLET";
                        }
                        if (bdry_cond_str == "ZNEUMANN" &&
                            face_conds[BDRY_LOC::ZHI] != BDRY_COND::BASIC::NEUMANN)
                        {
                            no_face_data_found = true;
                            proper_face_data = "NEUMANN";
                        }
                        if (bdry_cond_str == "ZREFLECT" &&
                            face_conds[BDRY_LOC::ZHI] != BDRY_COND::BASIC::REFLECT)
                        {
                            no_face_data_found = true;
                            proper_face_data = "REFLECT";
                        }
                        if (bdry_cond_str == "ZSYMMETRY" &&
                            face_conds[BDRY_LOC::ZHI] != BDRY_COND::BASIC::SYMMETRY)
                        {
                            no_face_data_found = true;
                            proper_face_data = "SYMMETRY";
                        }
                    }
                }
                if (no_face_data_found)
                {
                    TBOX_ERROR("BasicCartesianBoundaryUtilities3::read3dBdryEdges()\n"
                        << "Bdry condition '"
                        << bdry_cond_str
                        << "' found for '"
                        << bdry_loc_str
                        << "' but no '"
                        << proper_face_data
                        << "' data found for face '"
                        << proper_face
                        << "'."
                        << std::endl);
                }
            } // if (need_data_read)
        } // for (int ei = 0 ...
    } // if (num_per_dirs < 2)
    
    // Remove edge locations that have boundary conditions identified.
    edge_locs.erase(std::remove(edge_locs.begin(), edge_locs.end(), BOGUS_BDRY_LOC), edge_locs.end());
}


/*
 * Private function to read 3D node boundary data from input database.
 */
void
BasicCartesianBoundaryUtilities3::read3dBdryNodes(
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
        // node boundary data required
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
            if (bdry_cond_str == "XFLOW")
            {
                node_conds[s] = BDRY_COND::BASIC::XFLOW;
                node_locs[ni] = BOGUS_BDRY_LOC;
            }
            else if (bdry_cond_str == "YFLOW")
            {
                node_conds[s] = BDRY_COND::BASIC::YFLOW;
                node_locs[ni] = BOGUS_BDRY_LOC;
            }
            else if (bdry_cond_str == "ZFLOW")
            {
                node_conds[s] = BDRY_COND::BASIC::ZFLOW;
                node_locs[ni] = BOGUS_BDRY_LOC;
            }
            else if (bdry_cond_str == "XREFLECT")
            {
                node_conds[s] = BDRY_COND::BASIC::XREFLECT;
                node_locs[ni] = BOGUS_BDRY_LOC;
            }
            else if (bdry_cond_str == "YREFLECT")
            {
                node_conds[s] = BDRY_COND::BASIC::YREFLECT;
                node_locs[ni] = BOGUS_BDRY_LOC;
            }
            else if (bdry_cond_str == "ZREFLECT")
            {
                node_conds[s] = BDRY_COND::BASIC::ZREFLECT;
                node_locs[ni] = BOGUS_BDRY_LOC;
            }
            else if (bdry_cond_str == "XSYMMETRY")
            {
                node_conds[s] = BDRY_COND::BASIC::XSYMMETRY;
                node_locs[ni] = BOGUS_BDRY_LOC;
            }
            else if (bdry_cond_str == "YSYMMETRY")
            {
                node_conds[s] = BDRY_COND::BASIC::YSYMMETRY;
                node_locs[ni] = BOGUS_BDRY_LOC;
            }
            else if (bdry_cond_str == "ZSYMMETRY")
            {
                node_conds[s] = BDRY_COND::BASIC::ZSYMMETRY;
                node_locs[ni] = BOGUS_BDRY_LOC;
            }
            else if (bdry_cond_str == "XDIRICHLET")
            {
                node_conds[s] = BDRY_COND::BASIC::XDIRICHLET;
                node_locs[ni] = BOGUS_BDRY_LOC;
            }
            else if (bdry_cond_str == "YDIRICHLET")
            {
                node_conds[s] = BDRY_COND::BASIC::YDIRICHLET;
                node_locs[ni] = BOGUS_BDRY_LOC;
            }
            else if (bdry_cond_str == "ZDIRICHLET")
            {
                node_conds[s] = BDRY_COND::BASIC::ZDIRICHLET;
                node_locs[ni] = BOGUS_BDRY_LOC;
            }
            else if (bdry_cond_str == "XNEUMANN")
            {
                node_conds[s] = BDRY_COND::BASIC::XNEUMANN;
                node_locs[ni] = BOGUS_BDRY_LOC;
            }
            else if (bdry_cond_str == "YNEUMANN")
            {
                node_conds[s] = BDRY_COND::BASIC::YNEUMANN;
                node_locs[ni] = BOGUS_BDRY_LOC;
            }
            else if (bdry_cond_str == "ZNEUMANN")
            {
                node_conds[s] = BDRY_COND::BASIC::ZNEUMANN;
                node_locs[ni] = BOGUS_BDRY_LOC;
            }
            
            std::string proper_face;
            std::string proper_face_data;
            bool no_face_data_found = false;
            if (bdry_cond_str == "XFLOW" ||
                bdry_cond_str == "XDIRICHLET" ||
                bdry_cond_str == "XNEUMANN" ||
                bdry_cond_str == "XREFLECT" ||
                bdry_cond_str == "XSYMMETRY")
            {
                if (s == NODE_BDRY_LOC_3D::XLO_YLO_ZLO ||
                    s == NODE_BDRY_LOC_3D::XLO_YHI_ZLO ||
                    s == NODE_BDRY_LOC_3D::XLO_YLO_ZHI ||
                    s == NODE_BDRY_LOC_3D::XLO_YHI_ZHI)
                {
                    proper_face = "XLO";
                    if (bdry_cond_str == "XFLOW" &&
                        face_conds[BDRY_LOC::XLO] != BDRY_COND::BASIC::FLOW)
                    {
                        no_face_data_found = true;
                        proper_face_data = "FLOW";
                    }
                    if (bdry_cond_str == "XDIRICHLET" &&
                        face_conds[BDRY_LOC::XLO] != BDRY_COND::BASIC::DIRICHLET)
                    {
                        no_face_data_found = true;
                        proper_face_data = "DIRICHLET";
                    }
                    if (bdry_cond_str == "XNEUMANN" &&
                        face_conds[BDRY_LOC::XLO] != BDRY_COND::BASIC::NEUMANN)
                    {
                        no_face_data_found = true;
                        proper_face_data = "NEUMANN";
                    }
                    if (bdry_cond_str == "XREFLECT" &&
                        face_conds[BDRY_LOC::XLO] != BDRY_COND::BASIC::REFLECT)
                    {
                        no_face_data_found = true;
                        proper_face_data = "REFLECT";
                    }
                    if (bdry_cond_str == "XSYMMETRY" &&
                        face_conds[BDRY_LOC::XLO] != BDRY_COND::BASIC::SYMMETRY)
                    {
                        no_face_data_found = true;
                        proper_face_data = "SYMMETRY";
                    }
                }
                else
                {
                    proper_face = "XHI";
                    if (bdry_cond_str == "XFLOW" &&
                        face_conds[BDRY_LOC::XHI] != BDRY_COND::BASIC::FLOW)
                    {
                        no_face_data_found = true;
                        proper_face_data = "FLOW";
                    }
                    if (bdry_cond_str == "XDIRICHLET" &&
                        face_conds[BDRY_LOC::XHI] != BDRY_COND::BASIC::DIRICHLET)
                    {
                        no_face_data_found = true;
                        proper_face_data = "DIRICHLET";
                    }
                    if (bdry_cond_str == "XNEUMANN" &&
                        face_conds[BDRY_LOC::XHI] != BDRY_COND::BASIC::NEUMANN)
                    {
                        no_face_data_found = true;
                        proper_face_data = "NEUMANN";
                    }
                    if (bdry_cond_str == "XREFLECT" &&
                        face_conds[BDRY_LOC::XHI] != BDRY_COND::BASIC::REFLECT)
                    {
                       no_face_data_found = true;
                       proper_face_data = "REFLECT";
                    }
                    if (bdry_cond_str == "XSYMMETRY" &&
                        face_conds[BDRY_LOC::XHI] != BDRY_COND::BASIC::SYMMETRY)
                    {
                       no_face_data_found = true;
                       proper_face_data = "SYMMETRY";
                    }
                }
            }
            else if (bdry_cond_str == "YFLOW" ||
                     bdry_cond_str == "YDIRICHLET" ||
                     bdry_cond_str == "YNEUMANN" ||
                     bdry_cond_str == "YREFLECT" ||
                     bdry_cond_str == "YSYMMETRY")
            {
                if (s == NODE_BDRY_LOC_3D::XLO_YLO_ZLO ||
                    s == NODE_BDRY_LOC_3D::XHI_YLO_ZLO ||
                    s == NODE_BDRY_LOC_3D::XLO_YLO_ZHI ||
                    s == NODE_BDRY_LOC_3D::XHI_YLO_ZHI)
                {
                    proper_face = "YLO";
                    if (bdry_cond_str == "YFLOW" &&
                        face_conds[BDRY_LOC::YLO] != BDRY_COND::BASIC::FLOW)
                    {
                        no_face_data_found = true;
                        proper_face_data = "FLOW";
                    }
                    if (bdry_cond_str == "YDIRICHLET" &&
                        face_conds[BDRY_LOC::YLO] != BDRY_COND::BASIC::DIRICHLET)
                    {
                        no_face_data_found = true;
                        proper_face_data = "DIRICHLET";
                    }
                    if (bdry_cond_str == "YNEUMANN" &&
                        face_conds[BDRY_LOC::YLO] != BDRY_COND::BASIC::NEUMANN)
                    {
                        no_face_data_found = true;
                        proper_face_data = "NEUMANN";
                    }
                    if (bdry_cond_str == "YREFLECT" &&
                        face_conds[BDRY_LOC::YLO] != BDRY_COND::BASIC::REFLECT)
                    {
                        no_face_data_found = true;
                        proper_face_data = "REFLECT";
                    }
                    if (bdry_cond_str == "YSYMMETRY" &&
                        face_conds[BDRY_LOC::YLO] != BDRY_COND::BASIC::SYMMETRY)
                    {
                        no_face_data_found = true;
                        proper_face_data = "SYMMETRY";
                    }
                }
                else
                {
                    proper_face = "YHI";
                    if (bdry_cond_str == "YFLOW" &&
                        face_conds[BDRY_LOC::YHI] != BDRY_COND::BASIC::FLOW)
                    {
                        no_face_data_found = true;
                        proper_face_data = "FLOW";
                    }
                    if (bdry_cond_str == "YDIRICHLET" &&
                        face_conds[BDRY_LOC::YHI] != BDRY_COND::BASIC::DIRICHLET)
                    {
                        no_face_data_found = true;
                        proper_face_data = "DIRICHLET";
                    }
                    if (bdry_cond_str == "YNEUMANN" &&
                        face_conds[BDRY_LOC::YHI] != BDRY_COND::BASIC::NEUMANN)
                    {
                        no_face_data_found = true;
                        proper_face_data = "NEUMANN";
                    }
                    if (bdry_cond_str == "YREFLECT" &&
                        face_conds[BDRY_LOC::YHI] != BDRY_COND::BASIC::REFLECT)
                    {
                        no_face_data_found = true;
                        proper_face_data = "REFLECT";
                    }
                    if (bdry_cond_str == "YSYMMETRY" &&
                        face_conds[BDRY_LOC::YHI] != BDRY_COND::BASIC::SYMMETRY)
                    {
                        no_face_data_found = true;
                        proper_face_data = "SYMMETRY";
                    }
                }
            }
            else if (bdry_cond_str == "ZFLOW" ||
                     bdry_cond_str == "ZDIRICHLET" ||
                     bdry_cond_str == "ZNEUMANN" ||
                     bdry_cond_str == "ZREFLECT" ||
                     bdry_cond_str == "ZSYMMETRY")
            {
                if (s == NODE_BDRY_LOC_3D::XLO_YLO_ZLO ||
                    s == NODE_BDRY_LOC_3D::XHI_YLO_ZLO ||
                    s == NODE_BDRY_LOC_3D::XLO_YHI_ZLO ||
                    s == NODE_BDRY_LOC_3D::XHI_YHI_ZLO)
                {
                    proper_face = "ZLO";
                    if (bdry_cond_str == "ZFLOW" &&
                        face_conds[BDRY_LOC::ZLO] != BDRY_COND::BASIC::FLOW)
                    {
                        no_face_data_found = true;
                        proper_face_data = "FLOW";
                    }
                    if (bdry_cond_str == "ZDIRICHLET" &&
                        face_conds[BDRY_LOC::ZLO] != BDRY_COND::BASIC::DIRICHLET)
                    {
                        no_face_data_found = true;
                        proper_face_data = "DIRICHLET";
                    }
                    if (bdry_cond_str == "ZNEUMANN" &&
                        face_conds[BDRY_LOC::ZLO] != BDRY_COND::BASIC::NEUMANN)
                    {
                        no_face_data_found = true;
                        proper_face_data = "NEUMANN";
                    }
                    if (bdry_cond_str == "ZREFLECT" &&
                        face_conds[BDRY_LOC::ZLO] != BDRY_COND::BASIC::REFLECT)
                    {
                        no_face_data_found = true;
                        proper_face_data = "REFLECT";
                    }
                    if (bdry_cond_str == "ZSYMMETRY" &&
                        face_conds[BDRY_LOC::ZLO] != BDRY_COND::BASIC::SYMMETRY)
                    {
                        no_face_data_found = true;
                        proper_face_data = "SYMMETRY";
                    }
                }
                else
                {
                    proper_face = "ZHI";
                    if (bdry_cond_str == "ZFLOW" &&
                        face_conds[BDRY_LOC::ZHI] != BDRY_COND::BASIC::FLOW)
                    {
                        no_face_data_found = true;
                        proper_face_data = "FLOW";
                    }
                    if (bdry_cond_str == "ZDIRICHLET" &&
                        face_conds[BDRY_LOC::ZHI] != BDRY_COND::BASIC::DIRICHLET)
                    {
                        no_face_data_found = true;
                        proper_face_data = "DIRICHLET";
                    }
                    if (bdry_cond_str == "ZNEUMANN" &&
                        face_conds[BDRY_LOC::ZHI] != BDRY_COND::BASIC::NEUMANN)
                    {
                        no_face_data_found = true;
                        proper_face_data = "NEUMANN";
                    }
                    if (bdry_cond_str == "ZREFLECT" &&
                        face_conds[BDRY_LOC::ZHI] != BDRY_COND::BASIC::REFLECT)
                    {
                        no_face_data_found = true;
                        proper_face_data = "REFLECT";
                    }
                    if (bdry_cond_str == "ZSYMMETRY" &&
                        face_conds[BDRY_LOC::ZHI] != BDRY_COND::BASIC::SYMMETRY)
                    {
                        no_face_data_found = true;
                        proper_face_data = "SYMMETRY";
                    }
                }
            }
            if (no_face_data_found)
            {
                TBOX_ERROR("BasicCartesianBoundaryUtilities3::read3dBdryNodes()\n"
                    << "Bdry condition '"
                    << bdry_cond_str
                    << "' found for '"
                    << bdry_loc_str
                    << "' but no '"
                    << proper_face_data
                    << "' data found for face '"
                    << proper_face
                    << "'."
                    << std::endl);
            }
        } // for (int ni = 0 ...
    } // if (num_per_dirs < 1)
    
    // Remove node locations that have boundary conditions identified.
    node_locs.erase(std::remove(node_locs.begin(), node_locs.end(), BOGUS_BDRY_LOC), node_locs.end());
}
