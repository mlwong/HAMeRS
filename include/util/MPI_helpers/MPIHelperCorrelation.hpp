#ifndef MPI_HELPER_CORRELATION_HPP
#define MPI_HELPER_CORRELATION_HPP

#include "util/MPI_helpers/MPIHelper.hpp"

#include "SAMRAI/pdat/CellVariable.h"

class MPIHelperCorrelation: public MPIHelper
{
    public:
        MPIHelperCorrelation(
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
         * Compute correlation with only x-direction as inhomogeneous direction.
         */
        std::vector<double> getQuantityCorrelationWithInhomogeneousXDirection(
            std::vector<HAMERS_SHARED_PTR<pdat::CellVariable<double> > >& variable_quantities,
            const std::vector<int>& component_indices,
            const std::vector<std::vector<double> >& averaged_quantities,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Compute correlation with only x-direction as inhomogeneous direction.
         */
        std::vector<double> getQuantityCorrelationWithInhomogeneousXDirection(
            std::vector<HAMERS_SHARED_PTR<pdat::CellVariable<double> > >& variable_quantities,
            const std::vector<int>& component_indices,
            const std::vector<bool>& use_reciprocal,
            const std::vector<std::vector<double> >& averaged_quantities,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
    private:
        
};

#endif /* MPI_HELPER_CORRELATION_HPP */
