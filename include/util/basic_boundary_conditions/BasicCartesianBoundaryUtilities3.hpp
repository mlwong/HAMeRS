/*************************************************************************
 *
 * This file is modified from CartesianBoundaryUtilities3.h of the SAMRAI
 * distribution. For full copyright information, see COPYRIGHT and
 * COPYING.LESSER of SAMRAI distribution.
 *
 * Copyright:     (c) 1997-2014 Lawrence Livermore National Security, LLC
 * Description:   Utility routines for manipulating Cartesian 3d boundary data
 *
 ************************************************************************/

#ifndef CARTESIAN_BOUNDARY_UTILITIES3_HPP
#define CARTESIAN_BOUNDARY_UTILITIES3_HPP

#include "HAMeRS_config.hpp"

#include "HAMeRS_memory.hpp"

#include "util/basic_boundary_conditions/BoundaryUtilityStrategy.hpp"

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/tbox/Database.h"

#include <string>
#include <vector>

using namespace SAMRAI;

/*!
 * @brief Class BasicCartesianBoundaryUtilities3 is a utility class that
 * simplifies the implementation of simple physical boundary data in
 * 3 spatial dimensions.  It contains routines for reading boundary data
 * information from input files, applying those boundary conditions,
 * and error checking boundary data.  These routines apply to the
 * case of cell-centered double data only.  One may use all of these
 * capabilities, or use the input reading, boundary setting, and error
 * checking routines independently.
 *
 * <b> Input Parameters </b>
 *
 * To use the boundary condition input reading capabilities, the format
 * of the input file section containing the boundary information must
 * be as described next.  Boundary face, node, and edge entries are only
 * required for those that are not filled automatically when periodic
 * conditions apply.
 *
 * The boundary condition for face "*" is provided in a section as follows:
 *
 * @code
 *    boundary_face_* {
 *       boundary_condition  = ...  // boundary condition string identifier
 *    }
 * @endcode
 *
 * Allowable face identifiers (i.e., values for "*") are: <br>
 *       xlo, xhi, ylo, yhi, zlo, zhi <br>
 * Supported face boundary_condition string values are: <br>
 *       "FLOW", "REFLECT", "DIRICHLET", "NEUMANN"
 *
 * The boundary condition for edge "*" is provided in a section as follows:
 *
 * @code
 *    boundary_edge_* {
 *       boundary_condition  = ...  // boundary condition string identifier
 *    }
 * @endcode
 *
 * Allowable edge identifiers (i.e., values for "*") are: <br>
 *       ylo_zlo, yhi_zlo, ylo_zhi, yhi_zhi,
 *       xlo_zlo, xlo_zhi, xhi_zlo, xhi_zhi,
 *       xlo_ylo, xhi_ylo, xlo_yhi, xhi_yhi <br>
 * Supported edge boundary_condition string values are: <br>
 *       "XFLOW", "YFLOW", "ZFLOW",
 *       "XREFLECT", "YREFLECT", "ZREFLECT",
 *       "XDIRICHLET", "YDIRICHLET", "ZDIRICHLET"
 *       "XNEUMANN", "YNEUMANN", "ZNEUMANN"
 *
 * Note that edge conditions must be consistent with adjacent face conditions.
 *
 * The boundary condition for node "*" is provided in a section as follows:
 *
 * @code
 *    boundary_node_* {
 *       boundary_condition  = ...  // boundary condition string identifier
 *    }
 * @endcode
 *
 * Allowable node identifiers (i.e., values for "*") are: <br>
 *       xlo_ylo_zlo, xhi_ylo_zlo, xlo_yhi_zlo, xhi_yhi_zlo,
 *       xlo_ylo_zhi, xhi_ylo_zhi, xlo_yhi_zhi, xhi_yhi_zhi <br>
 * Supported node boundary_condition values are: <br>
 *       "XFLOW", "YFLOW", "ZFLOW",
 *       "XREFLECT", "YREFLECT", "ZREFLECT",
 *       "XDIRICHLET", "YDIRICHLET", "ZDIRICHLET"
 *       "XNEUMANN", "YNEUMANN", "ZNEUMANN"
 *
 * Note that node conditions must be consistent with adjacent face conditions.
 *
 * See the include file CartesianBoundaryDefines.h for integer constant
 * definitions that apply for the various boundary types, locations,
 * and boundary conditions.  If you choose to use the input reading
 * capabilities only and write your own boundary condition routines in
 * FORTRAN, you should note that the integer constants for the various
 * boundary condition types and locations are automatically "stuffed" into
 * FORTRAN common blocks.  This avoids potential problems with
 * inconsistencies between C++ and FORTRAN usage.  Please see the
 * FORTRAN include file cartbdryparams3d.i for details.
 *
 * @see appu::BoundaryUtilityStrategy3
 */

struct BasicCartesianBoundaryUtilities3
{
    public:
        /*!
         * Function to read 3d boundary data from input database.
         * The integer boundary condition types are placed in the integer
         * arrays supplied by the caller (typically, the concrete
         * BoundaryUtilityStrategy object provided).  When DIRICHLET or
         * NEUMANN conditions are specified, control is passed to the
         * BoundaryUtilityStrategy to read the boundary state data specific to
         * the problem.
         *
         * Errors will be reported and the program will abort whenever necessary
         * boundary condition information is missing in the input database, or
         * when the data read in is either unknown or inconsistent.  The periodic
         * domain information is used to determine which boundary face, edge, or
         * node entries are not required from input.  Error checking
         * requires that node and edge boundary conditions are consistent
         * with those specified for the faces.
         *
         *
         * When assertion checking is active, assertions will result when any
         * of the pointer arguments is null, or an array is passed in with the
         * the wrong size.
         *
         * @param bdry_strategy user-defined object that reads DIRICHLET or NEUMANN
         *                      conditions
         * @param input_db      input database containing all boundary data
         * @param face_locs     array of locations of faces for applying
         *                      boundary conditions.
         * @param edge_locs     array of locations of edges for applying
         *                      boundary conditions.
         * @param node_locs     array of locations of nodes for applying
         *                      boundary conditions.
         * @param face_conds    array into which integer face boundary condition
         *                      types are read
         * @param edge_conds    array into which integer edge boundary condition
         *                      types are read
         * @param node_conds    array into which integer node boundary condition
         *                      types are read
         * @param periodic      integer vector specifying which coordinate
         *                      directions are periodic (e.g., value returned from
         *                      GridGeometry2::getPeriodicShift())
         *
         * @pre input_db
         * @pre periodic.getDim() == tbox::Dimension(3)
         * @pre bdry_strategy != 0
         * @pre edge_locs.size() > 0 and <= NUM_3D_FACES
         * @pre edge_locs.size() > 0 and <= NUM_3D_EDGES
         * @pre node_locs.size() > 0 and <= NUM_3D_NODES
         * @pre face_conds.size() == NUM_3D_FACES
         * @pre edge_conds.size() == NUM_3D_EDGES
         * @pre node_conds.size() == NUM_3D_NODES
         */
        static void
        getFromInput(
            BoundaryUtilityStrategy* bdry_strategy,
            const HAMERS_SHARED_PTR<tbox::Database>& input_db,
            const std::vector<int>& face_locs,
            const std::vector<int>& edge_locs,
            const std::vector<int>& node_locs,
            std::vector<int>& face_conds,
            std::vector<int>& edge_conds,
            std::vector<int>& node_conds,
            const hier::IntVector& periodic);
        
        /*!
         * Function to remove 3d boundary faces with boundary conditions filled by this class for a patch.
         *
         * @param bdry_face_locs      Array of locations of faces for applying
         *                            boundary conditions.
         * @param patch               hier::Patch on which data object lives.
         * @param bdry_face_conds     Array of boundary condition types for
         *                            patch faces.
         *
         * @pre bdry_face_locs.size() <= NUM_3D_FACES
         * @pre bdry_face_conds.size() == NUM_3D_FACES
         */
        static void
        removeBoundaryFaceLocations(
            std::vector<int>& bdry_face_locs,
            const hier::Patch& patch,
            const std::vector<int>& bdry_face_conds);
        
        /*!
         * Function to fill 3d face boundary values for a patch.
         *
         * When assertion checking is active, assertions will result when any
         * of the pointer arguments is null, or an array is passed in with the
         * the wrong size.
         *
         * @param var_name            String name of variable (for error reporting).
         * @param var_data            Cell-centered patch data object to fill.
         * @param patch               hier::Patch on which data object lives.
         * @param bdry_face_locs      Array of locations of faces for applying
         *                            boundary conditions.
         * @param bdry_face_conds     Array of boundary condition types for
         *                            patch faces.
         * @param bdry_face_values    Array of boundary values for patch
         *                            faces.
         * @param ghost_width_to_fill Width of ghost region to fill.
         *
         * @pre !var_name.empty()
         * @pre var_data
         * @pre bdry_face_locs.size() <= NUM_3D_FACES
         * @pre bdry_face_conds.size() == NUM_3D_FACES
         * @pre bdry_face_values.size() == NUM_3D_FACES * (var_data->getDepth())
         * @pre ghost_fill_width.getDim() == tbox::Dimension(3)
         * @pre (var_data->getDim() == patch.getDim()) &&
         *      (var_data->getDim() == ghost_fill_width.getDim())
         */
        static void
        fillFaceBoundaryData(
            const std::string& var_name,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& var_data,
            const hier::Patch& patch,
            const std::vector<int>& bdry_face_locs,
            const std::vector<int>& bdry_face_conds,
            const std::vector<double>& bdry_face_values,
            const hier::IntVector& ghost_width_to_fill = -hier::IntVector::getOne(tbox::Dimension(3)));
        
        /*!
         * Function to remove 3d boundary edges with boundary conditions filled by this class for a patch.
         *
         * @param bdry_edge_locs      Array of locations of edges for applying
         *                            boundary conditions.
         * @param patch               hier::Patch on which data object lives.
         * @param bdry_edge_conds     Array of boundary condition types for
         *                            patch edges.
         *
         * @pre bdry_edge_locs.size() <= NUM_3D_EDGES
         * @pre bdry_edge_conds.size() == NUM_3D_EDGES
         */
        static void
        removeBoundaryEdgeLocations(
            std::vector<int>& bdry_edge_locs,
            const hier::Patch& patch,
            const std::vector<int>& bdry_edge_conds);
        
        /*!
         * Function to fill 3d edge boundary values for a patch.
         *
         * When assertion checking is active, assertions will result when any
         * of the pointer arguments is null, or an array is passed in with the
         * the wrong size.
         *
         * @param var_name            String name of variable (for error reporting).
         * @param var_data            Cell-centered patch data object to fill.
         * @param patch               hier::Patch on which data object lives.
         * @param bdry_edge_locs      Array of locations of edges for applying
         *                            boundary conditions.
         * @param bdry_edge_conds     Array of boundary condition types for
         *                            patch edges.
         * @param bdry_face_values    Array of boundary values for patch
         *                            faces.
         * @param ghost_width_to_fill Width of ghost region to fill.
         *
         * @pre !var_name.empty()
         * @pre var_data
         * @pre bdry_edge_locs.size() <= NUM_3D_EDGES
         * @pre bdry_edge_conds.size() == NUM_3D_EDGES
         * @pre bdry_face_values.size() == NUM_3D_FACES * (var_data->getDepth())
         * @pre ghost_fill_width.getDim() == tbox::Dimension(3)
         * @pre (var_data->getDim() == patch.getDim()) &&
         *      (var_data->getDim() == ghost_fill_width.getDim())
         */
        static void
        fillEdgeBoundaryData(
            const std::string& var_name,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& var_data,
            const hier::Patch& patch,
            const std::vector<int>& bdry_edge_locs,
            const std::vector<int>& bdry_edge_conds,
            const std::vector<double>& bdry_face_values,
            const hier::IntVector& ghost_width_to_fill = -hier::IntVector::getOne(tbox::Dimension(3)));
        
        /*!
         * Function to remove 3d boundary nodes with boundary conditions filled by this class for a patch.
         *
         * @param bdry_node_locs      Array of locations of nodes for applying
         *                            boundary conditions.
         * @param patch               hier::Patch on which data object lives.
         * @param bdry_node_conds     Array of boundary condition types for
         *                            patch nodes.
         *
         * @pre bdry_node_locs.size() <= NUM_3D_NODES
         * @pre bdry_node_conds.size() == NUM_3D_NODES
         */
        static void
        removeBoundaryNodeLocations(
            std::vector<int>& bdry_node_locs,
            const hier::Patch& patch,
            const std::vector<int>& bdry_node_conds);
        
        /*!
         * Function to fill 3d node boundary values for a patch.
         *
         * When assertion checking is active, assertions will result when any
         * of the pointer arguments is null, or an array is passed in with the
         * the wrong size.
         *
         * @param var_name            String name of variable (for error reporting).
         * @param var_data            Cell-centered patch data object to fill.
         * @param patch               hier::Patch on which data object lives.
         * @param bdry_node_locs      Array of locations of nodes for applying
         *                            boundary conditions.
         * @param bdry_node_conds     Array of boundary condition types for
         *                            patch nodes.
         * @param bdry_face_values    Array of boundary values for patch
         *                            faces.
         * @param ghost_width_to_fill Width of ghost region to fill.
         *
         * @pre !var_name.empty()
         * @pre var_data
         * @pre bdry_node_locs.size() <= NUM_3D_NODES
         * @pre bdry_node_conds.size() == NUM_3D_NODES
         * @pre bdry_face_values.size() == NUM_3D_FACES * (var_data->getDepth())
         * @pre ghost_fill_width.getDim() == tbox::Dimension(3)
         * @pre (var_data->getDim() == patch.getDim()) &&
         *      (var_data->getDim() == ghost_fill_width.getDim())
         */
        static void
        fillNodeBoundaryData(
            const std::string& var_name,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& var_data,
            const hier::Patch& patch,
            const std::vector<int>& bdry_node_locs,
            const std::vector<int>& bdry_node_conds,
            const std::vector<double>& bdry_face_values,
            const hier::IntVector& ghost_width_to_fill = -hier::IntVector::getOne(tbox::Dimension(3)));
        
        /*!
         * Function that returns the integer face boundary location
         * corresponding to the given edge location and edge boundary
         * condition.
         *
         * If the edge boundary condition type or edge location are unknown,
         * or the boundary condition type is inconsistent with the edge location,
         * an error code (-1) is returned.
         *
         * @return Integer face location for edge location and boundary condition
         *         type.
         *
         * @param edge_loc   Integer location for edge.
         * @param edge_btype Integer boundary condition type for edge.
         *
         * @pre (edge_btype == BdryCond::Basic::XFLOW) ||
         *      (edge_btype == BdryCond::Basic::XREFLECT) ||
         *      (edge_btype == BdryCond::Basic::XSYMMETRY) ||
         *      (edge_btype == BdryCond::Basic::XDIRICHLET) ||
         *      (edge_btype == BdryCond::Basic::XNEUMANN) ||
         *      (edge_btype == BdryCond::Basic::YFLOW) ||
         *      (edge_btype == BdryCond::Basic::YREFLECT) ||
         *      (edge_btype == BdryCond::Basic::YSYMMETRY) ||
         *      (edge_btype == BdryCond::Basic::YDIRICHLET) ||
         *      (edge_btype == BdryCond::Basic::YNEUMANN) ||
         *      (edge_btype == BdryCond::Basic::ZFLOW) ||
         *      (edge_btype == BdryCond::Basic::ZREFLECT) ||
         *      (edge_btype == BdryCond::Basic::ZSYMMETRY) ||
         *      (edge_btype == BdryCond::Basic::ZDIRICHLET) ||
         *      (edge_btype == BdryCond::Basic::ZNEUMANN)
         */
        static int
        getFaceLocationForEdgeBdry(
            int edge_loc,
            int edge_btype);
        
        /*!
         * Function that returns the integer face boundary location
         * corresponding to the given node location and node boundary
         * condition.
         *
         * If the node boundary condition type or node location are unknown,
         * or the boundary condition type is inconsistent with the node location,
         * an error code (-1) is returned.
         *
         * @return Integer face location for node location and boundary condition
         *         type.
         *
         * @param node_loc   Integer location for node.
         * @param node_btype Integer boundary condition type for node.
         *
         * @pre (edge_btype == BdryCond::Basic::XFLOW) ||
         *      (edge_btype == BdryCond::Basic::XREFLECT) ||
         *      (edge_btype == BdryCond::Basic::XSYMMETRY) ||
         *      (edge_btype == BdryCond::Basic::XDIRICHLET) ||
         *      (edge_btype == BdryCond::Basic::XNEUMANN) ||
         *      (edge_btype == BdryCond::Basic::YFLOW) ||
         *      (edge_btype == BdryCond::Basic::YREFLECT) ||
         *      (edge_btype == BdryCond::Basic::YSYMMETRY) ||
         *      (edge_btype == BdryCond::Basic::YDIRICHLET) ||
         *      (edge_btype == BdryCond::Basic::YNEUMANN) ||
         *      (edge_btype == BdryCond::Basic::ZFLOW) ||
         *      (edge_btype == BdryCond::Basic::ZREFLECT) ||
         *      (edge_btype == BdryCond::Basic::ZSYMMETRY) ||
         *      (edge_btype == BdryCond::Basic::ZDIRICHLET) ||
         *      (edge_btype == BdryCond::Basic::ZNEUMANN)
         */
        static int
        getFaceLocationForNodeBdry(
            int node_loc,
            int node_btype);
        
    private:
        static void
        read3dBdryFaces(
            BoundaryUtilityStrategy* bdry_strategy,
            const HAMERS_SHARED_PTR<tbox::Database>& input_db,
            const std::vector<int>& face_locs,
            std::vector<int>& face_conds,
            const hier::IntVector& periodic);
        
        static void
        read3dBdryEdges(
            const HAMERS_SHARED_PTR<tbox::Database>& input_db,
            const std::vector<int>& edge_locs,
            const std::vector<int>& face_conds,
            std::vector<int>& edge_conds,
            const hier::IntVector& periodic);
        
        static void
        read3dBdryNodes(
            const HAMERS_SHARED_PTR<tbox::Database>& input_db,
            const std::vector<int>& node_locs,
            const std::vector<int>& face_conds,
            std::vector<int>& node_conds,
            const hier::IntVector& periodic);
        
};

#endif /* CARTESIAN_BOUNDARY_UTILITIES3_HPP */
