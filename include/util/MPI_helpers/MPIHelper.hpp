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
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy);
        
        /*
         * Get refinement ratio from the finest level to the coarsest level.
         */
        const hier::IntVector&
        getRatioFinestLevelToCoarsestLevel() const
        {
            return d_ratio_finest_level_to_coarsest_level;
        }
        
        /*
         * Get number of points in the coarsest domain.
         */
        const hier::IntVector&
        getCoarsestDomainNumberOfPoints() const
        {
            return d_coarsest_level_dims;
        }
        
        /*
         * Get grid spacing of the coarsest domain.
         */
        const std::vector<Real>&
        getCoarsestDomainGridSpacing() const
        {
            return d_dx_coarsest_level_dims;
        }
        
        /*
         * Get number of points in the finest refined domain.
         */
        const hier::IntVector&
        getFinestRefinedDomainNumberOfPoints() const
        {
            return d_finest_level_dims;
        }
        
        /*
         * Get grid spacing of the finest refined domain.
         */
        const std::vector<Real>&
        getFinestRefinedDomainGridSpacing() const
        {
            return d_dx_finest_level_dims;
        }
        
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
        
        /*
         * Refinement ratio from the finest level to the coarsest level.
         */
        hier::IntVector d_ratio_finest_level_to_coarsest_level;
        
        /*
         * Number of points in the coarsest domain.
         */
        hier::IntVector d_coarsest_level_dims;
        
        /*
         * Grid spacing of the coarsest domain.
         */
        std::vector<Real> d_dx_coarsest_level_dims;
        
        /*
         * Number of points in the finest refined domain.
         */
        hier::IntVector d_finest_level_dims;
        
        /*
         * Grid spacing of the finest refined domain.
         */
        std::vector<Real> d_dx_finest_level_dims;
        
};

#endif /* MPI_HELPER_HPP */
