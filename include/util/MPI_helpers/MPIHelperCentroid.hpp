#ifndef MPI_HELPER_CENTROID_HPP
#define MPI_HELPER_CENTROID_HPP

#include "util/MPI_helpers/MPIHelper.hpp"

#include "SAMRAI/pdat/CellVariable.h"

class MPIHelperCentroid: public MPIHelper
{
    public:
        MPIHelperCentroid(
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
         * Compute centroid in x-direction.
         */
        Real getCentroidInXDirection(
            HAMERS_SHARED_PTR<pdat::CellVariable<Real> >& variable_quantity,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Compute centroid in y-direction.
         */
        Real getCentroidInYDirection(
            HAMERS_SHARED_PTR<pdat::CellVariable<Real> >& variable_quantity,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Compute centroid in z-direction.
         */
        Real getCentroidInZDirection(
            HAMERS_SHARED_PTR<pdat::CellVariable<Real> >& variable_quantity,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
    private:
        
};

#endif /* MPI_HELPER_CENTROID_HPP */
