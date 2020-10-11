/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER of SAMRAI distribution.
 *
 * Copyright:     (c) 1997-2014 Lawrence Livermore National Security, LLC
 * Description:   Interface for processing user-defined boundary data in
 *                CartesianBoundaryUtilities classes
 *
 ************************************************************************/

#ifndef BOUNDARY_UTILITY_STRATEGY_HPP
#define BOUNDARY_UTILITY_STRATEGY_HPP

#include "HAMeRS_config.hpp"

#include "HAMeRS_memory.hpp"

#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Utilities.h"

#include <string>

using namespace SAMRAI;

/*!
 * @brief Class BoundaryUtilityStrategy is an abstract base
 * class that declares an interface that allows application code to
 * read problem-specific boundary data when using the boundary
 * utilities.  Currently, there are two virtual member functions defined.
 * One allows users to read problem-specific DIRICHLET boundary values
 * from an input database; the other does the same for NEUMANN boundary
 * values.  More virtual functions may be added in the future
 * as additional boundary conditions are supported.
 *
 * See the include file CartesianBoundaryDefines.hp for integer constant
 * definitions that apply for the various boundary types, locations,
 * and boundary conditions.
 *
 * @see BasicBoundaryUtilities2
 * @see BasicBoundaryUtilities3
 */

class BoundaryUtilityStrategy
{
    public:
        /*!
         * The default constructor for BoundaryUtilityStrategy does nothing
         * interesting.
         */
        BoundaryUtilityStrategy();
        
        /*!
         * The destructor for BoundaryUtilityStrategy does nothing
         * interesting.
         */
        virtual ~BoundaryUtilityStrategy();
        
        /*!
         * Read DIRICHLET boundary state values for an edge (in 2d) or a face
         * (in 3d) from a given database.
         *
         * @param db      Input database from which to read boundary values.
         * @param db_name Name of input database (e.g., for error reporting).
         * @param bdry_location_index Integer index for location of edge (in 2d)
         *                            or face (in 3d) boundary.
         */
        virtual void
        readDirichletBoundaryDataEntry(
            const boost::shared_ptr<tbox::Database>& db,
            std::string& db_name,
            int bdry_location_index) = 0;
        
        /*!
         * Read NEUMANN boundary state values for an edge (in 2d) or a face
         * (in 3d) from a given database.
         *
         * @param db      Input database from which to read boundary values.
         * @param db_name Name of input database (e.g., for error reporting).
         * @param bdry_location_index Integer index for location of edge (in 2d)
         *                            or face (in 3d) boundary.
         */
        virtual void
        readNeumannBoundaryDataEntry(
            const boost::shared_ptr<tbox::Database>& db,
            std::string& db_name,
            int bdry_location_index) = 0;
        
};

#endif /* BOUNDARY_UTILITY_STRATEGY_HPP */
