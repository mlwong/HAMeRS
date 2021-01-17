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
        double getNumberOfCells(
            HAMERS_SHARED_PTR<pdat::CellVariable<double> >& variable_quantity,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Compute weighted number of cells.
         */
        double getWeightedNumberOfCells(
            HAMERS_SHARED_PTR<pdat::CellVariable<double> >& variable_quantity,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
    private:
        
};

#endif /* MPI_HELPER_NUMBER_OF_CELLS_HPP */
