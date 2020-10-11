/*************************************************************************
 *
 * This file is modified from CartesianBoundaryUtilities2.h of the
 * SAMRAI distribution. For full copyright information, see COPYRIGHT and
 * COPYING.LESSER of SAMRAI distribution.
 *
 * Copyright:     (c) 1997-2014 Lawrence Livermore National Security, LLC
 * Description:   Utility routines for manipulating Cartesian 2d boundary
 *                data
 *
 ************************************************************************/

#ifndef BASIC_CARTESIAN_BOUNDARY_UTILITIES2_HPP
#define BASIC_CARTESIAN_BOUNDARY_UTILITIES2_HPP

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

/*
 * Definitions for Face, Edge, and Node boundary conditions:
 *
 * Note that FLOW, REFLECT, DIRICHLET and NEUMANN are used only for:
 * - Node boundary conditions in 1d, or
 * - Edge boundary conditions in 2d, or
 * - Face boundary conditions in 3d.
 *
 * Note that [X, Y, Z]FLOW, [X, Y, Z]REFLECT, [X, Y, Z]DIRICHLET, and
 * [X, Y, Z]NEUMANN are used only for:
 * - Node boundary conditions in 2d (X and Y cases only), or
 * - Edge and Node boundary conditions in 3d.
 */


/*!
 * @brief Class BasicCartesianBoundaryUtilities2 is a utility class that
 * simplifies the implementation of simple physical boundary data in
 * 2 spatial dimensions.  It contains routines for reading boundary data
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
 * be as described next.  Boundary node and edge entries are only
 * required for those that are not filled automatically when periodic
 * conditions apply.
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
 *       xlo, xhi, ylo, yhi <br>
 * Supported edge boundary_condition string values are: <br>
 *       "FLOW", "REFLECT", "DIRICHLET", "NEUMANN"
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
 *       xlo_ylo, xhi_ylo, xlo_yhi, xhi_yhi <br>
 * Supported node boundary_condition string values are: <br>
 *       "XFLOW", "YFLOW",
 *       "XREFLECT", "YREFLECT",
 *       "XDIRICHLET", "YDIRICHLET",
 *       "XNEUMANN", "YNEUMANN"
 *
 * Note that node conditions must be consistent with adjacent edge conditions.
 *
 * See the include file CartesianBoundaryDefines.h for integer constant
 * definitions that apply for the various boundary types, locations,
 * and boundary conditions.  If you choose to use the input reading
 * capabilities only and write your own boundary condition routines in
 * FORTRAN, you should note that the integer constants for the various
 * boundary condition types and locations are automatically "stuffed" into
 * FORTRAN common blocks.  This avoids potential problems with
 * inconsistencies between C++ and FORTRAN usage.  Please see the
 * FORTRAN include file cartbdryparams2d.i for details.
 *
 * @see BoundaryUtilityStrategy2
 */

struct BasicCartesianBoundaryUtilities2
{
    public:
        /*!
         * Function to read 2d boundary data from input database.
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
         * domain information is used to determine which boundary edge or
         * node entries are not required from input.  Error checking
         * requires that node boundary conditions are consistent
         * with those specified for the edges.
         *
         * When assertion checking is active, assertions will result when any
         * of the pointer arguments is null, or an array is passed in with the
         * the wrong size.
         *
         * @param bdry_strategy user-defined object that reads DIRICHLET or NEUMANN
         *                      conditions
         * @param input_db      input database containing all boundary data
         * @param edge_locs     array of locations of edges for applying
         *                      boundary conditions.
         * @param node_locs     array of locations of nodes for applying
         *                      boundary conditions.
         * @param edge_conds    array into which integer edge boundary condition
         *                      types are read
         * @param node_conds    array into which integer node boundary condition
         *                      types are read
         * @param periodic      integer vector specifying which coordinate
         *                      directions are periodic (e.g., value returned from
         *                      GridGeometry2::getPeriodicShift())
         *
         * @pre input_db
         * @pre periodic.getDim() == tbox::Dimension(2)
         * @pre bdry_strategy != 0
         * @pre edge_locs.size() > 0 and <= NUM_2D_EDGES
         * @pre node_locs.size() > 0 and <= NUM_2D_NODES
         * @pre edge_conds.size() == NUM_2D_EDGES
         * @pre node_conds.size() == NUM_2D_NODES
         */
        static void
        getFromInput(
            BoundaryUtilityStrategy* bdry_strategy,
            const boost::shared_ptr<tbox::Database>& input_db,
            const std::vector<int>& edge_locs,
            const std::vector<int>& node_locs,
            std::vector<int>& edge_conds,
            std::vector<int>& node_conds,
            const hier::IntVector& periodic);
        
        /*!
         * Function to fill 2d edge boundary values for a patch.
         *
         * When assertion checking is active, assertions will result when any
         * of the pointer arguments is null, or an array is passed in with the
         * the wrong size.
         *
         * @param var_name            String name of variable (for error reporting).
         * @param var_data            Cell-centered patch data object to fill.
         * @param patch               hier::Patch on which data object lives.
         * @param bdry_edge_locs      tbox::Array of locations of edges for applying
         *                            boundary conditions.
         * @param bdry_edge_conds     tbox::Array of boundary condition types for
         *                            patch edges.
         * @param bdry_edge_values    tbox::Array of boundary values for patch
         *                            edges.
         * @param ghost_width_to_fill Width of ghost region to fill.
         *
         * @pre !var_name.empty());
         * @pre var_data
         * @pre bdry_edge_locs.size() > 0 and <= NUM_2D_EDGES
         * @pre bdry_edge_conds.size() == NUM_2D_EDGES
         * @pre bdry_edge_values.size() == NUM_2D_EDGES * (var_data->getDepth())
         * @pre (var_data->getDim() == patch.getDim()) &&
         *      (var_data->getDim() == ghost_fill_width.getDim())
         * @pre ghost_fill_width.getDim() == tbox::Dimension(2)
         */
        static void
        fillEdgeBoundaryData(
            const std::string& var_name,
            const boost::shared_ptr<pdat::CellData<double> >& var_data,
            const hier::Patch& patch,
            const std::vector<int>& bdry_edge_locs,
            const std::vector<int>& bdry_edge_conds,
            const std::vector<double>& bdry_edge_values,
            const hier::IntVector& ghost_width_to_fill = -hier::IntVector::getOne(tbox::Dimension(2)));
        
        /*!
         * Function to fill 2d node boundary values for a patch.
         *
         * When assertion checking is active, assertions will result when any
         * of the pointer arguments is null, or an array is passed in with the
         * the wrong size.
         *
         * @param var_name            String name of variable (for error reporting).
         * @param var_data            Cell-centered patch data object to fill.
         * @param patch               hier::Patch on which data object lives.
         * @param bdry_node_locs      tbox::Array of locations of nodes for applying
         *                            boundary conditions.
         * @param bdry_node_conds     tbox::Array of boundary condition types for
         *                            patch nodes.
         * @param bdry_edge_values    tbox::Array of boundary values for patch
         *                            edges.
         * @param ghost_width_to_fill Width of ghost region to fill.
         *
         * @pre !var_name.empty()
         * @pre var_data
         * @pre bdry_node_conds.size() > 0 and <= NUM_2D_EDGES
         * @pre bdry_node_conds.size() == NUM_2D_NODES
         * @pre bdry_edge_values.size() == NUM_2D_EDGES * (var_data->getDepth())
         * @pre (var_data->getDim() == patch.getDim()) &&
         *      (var_data->getDim() == ghost_fill_width.getDim())
         * @pre ghost_fill_width.getDim() == tbox::Dimension(2)
         */
        static void
        fillNodeBoundaryData(
            const std::string& var_name,
            const boost::shared_ptr<pdat::CellData<double> >& var_data,
            const hier::Patch& patch,
            const std::vector<int>& bdry_node_locs,
            const std::vector<int>& bdry_node_conds,
            const std::vector<double>& bdry_edge_values,
            const hier::IntVector& ghost_width_to_fill = -hier::IntVector::getOne(tbox::Dimension(2)));
        
        /*!
         * Function that returns the integer edge boundary location
         * corresponding to the given node location and node boundary
         * condition.
         *
         * If the node boundary condition type or node location are unknown,
         * or the boundary condition type is inconsistent with the node location,
         * an error code (-1) is returned.
         *
         * @return Integer edge location for node location and boundary condition
         *         type.
         *
         * @param node_loc   Integer location for node.
         * @param node_btype Integer boundary condition type for node.
         *
         * @pre (node_btype == BdryCond::XFLOW) ||
         *      (node_btype == BdryCond::XREFLECT) ||
         *      (node_btype == BdryCond::XDIRICHLET) ||
         *      (node_btype == BdryCond::XNEUMANN) ||
         *      (node_btype == BdryCond::YFLOW) ||
         *      (node_btype == BdryCond::YREFLECT) ||
         *      (node_btype == BdryCond::YDIRICHLET) ||
         *      (node_btype == BdryCond::YNEUMANN)
         */
        static int
        getEdgeLocationForNodeBdry(
            int node_loc,
            int node_btype);
        
    private:
        static void
        read2dBdryEdges(
            BoundaryUtilityStrategy* bdry_strategy,
            const boost::shared_ptr<tbox::Database>& input_db,
            const std::vector<int>& edge_locs,
            std::vector<int>& edge_conds,
            const hier::IntVector& periodic);
        
        static void
        read2dBdryNodes(
            const boost::shared_ptr<tbox::Database>& input_db,
            const std::vector<int>& node_locs,
            const std::vector<int>& edge_conds,
            std::vector<int>& node_conds,
            const hier::IntVector& periodic);
        
};

#endif /* BASIC_CARTESIAN_BOUNDARY_UTILITIES2_HPP */
