#ifndef MPI_HELPER_HPP
#define MPI_HELPER_HPP

#include "HAMeRS_config.hpp"

#include "HAMeRS_memory.hpp"

#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/tbox/Dimension.h"

#include <string>

using namespace SAMRAI;

class MPIHelper
{
    public:
        MPIHelper(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy):
                d_object_name(object_name),
                d_dim(dim),
                d_grid_geometry(grid_geometry),
                d_patch_hierarchy(patch_hierarchy),
                d_mpi(tbox::SAMRAI_MPI::getSAMRAIWorld())
        {}
        
        /*
         * Get number of points in the x-direction of the refined domain.
         */
        int
        getRefinedDomainNumberOfPointsX() const;
        
        /*
         * Get grid spacing in the x-direction of the refined domain.
         */
        double
        getRefinedDomainGridSpacingX() const;
        
    protected:
        /*
         * The object name is used for error/warning reporting.
         */
        const std::string d_object_name;
        
        /*
         * Problem dimension.
         */
        const tbox::Dimension d_dim;
        
        /*
         * HAMERS_SHARED_PTR to the grid geometry.
         */
        const HAMERS_SHARED_PTR<geom::CartesianGridGeometry> d_grid_geometry;
        
        /*
         * HAMERS_SHARED_PTR to patch hierarchy.
         */
        HAMERS_SHARED_PTR<hier::PatchHierarchy> d_patch_hierarchy;
        
        /*
         * MPI communication object.
         */
        const tbox::SAMRAI_MPI& d_mpi;
        
};

#endif /* MPI_HELPER_HPP */
