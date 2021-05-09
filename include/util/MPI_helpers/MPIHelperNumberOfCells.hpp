#ifndef MPI_HELPER_NUMBER_OF_CELLS_HPP
#define MPI_HELPER_NUMBER_OF_CELLS_HPP

#include "util/MPI_helpers/MPIHelper.hpp"

#include "SAMRAI/pdat/CellVariable.h"

class MPIHelperNumberOfCells: public MPIHelper
{
    public:
        MPIHelperNumberOfCells(
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
        double getNumberOfCells() const;
        
        /*
         * Compute weighted number of cells.
         */
        double getWeightedNumberOfCells() const;
        
    private:
        
};

#endif /* MPI_HELPER_NUMBER_OF_CELLS_HPP */
