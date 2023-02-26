#ifndef MPI_HELPER_NUMBER_OF_CELLS_HPP
#define MPI_HELPER_NUMBER_OF_CELLS_HPP

#include "util/MPI_helpers/MPIHelper.hpp"

#include "SAMRAI/pdat/CellVariable.h"

class MPIHelperGrid: public MPIHelper
{
    public:
        MPIHelperGrid(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy):
                MPIHelper(
                    object_name,
                    dim,
                    grid_geometry,
                    patch_hierarchy)
        {}
        
        /*
         * Compute number of cells.
         */
        Real getNumberOfCells() const;
        
        /*
         * Compute weighted number of cells.
         */
        Real getWeightedNumberOfCells() const;
        
        /*
         * Compute averaged grid level number with only x-direction as inhomogeneous direction.
         */
        std::vector<Real>
        getAveragedGridLevelNumberWithInhomogeneousXDirection() const;
        
        /*
         * Compute averaged grid level number with only y-direction as inhomogeneous direction.
         */
        std::vector<Real>
        getAveragedGridLevelNumberWithInhomogeneousYDirection() const;
        
        /*
         * Compute averaged grid level number with only z-direction as inhomogeneous direction.
         */
        std::vector<Real>
        getAveragedGridLevelNumberWithInhomogeneousZDirection() const;
        
    private:
        
};

#endif /* MPI_HELPER_NUMBER_OF_CELLS_HPP */
