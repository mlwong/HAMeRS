#ifndef MPI_HELPER_MAX_MIN_HPP
#define MPI_HELPER_MAX_MIN_HPP

#include "util/MPI_helpers/MPIHelper.hpp"

#include "SAMRAI/pdat/CellVariable.h"

class MPIHelperMaxMin: public MPIHelper
{
    public:
        MPIHelperMaxMin(
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
         * Compute maximum value with only x-direction as inhomogeneous direction.
         */
        std::vector<double> getMaxQuantityWithInhomogeneousXDirection(
            HAMERS_SHARED_PTR<pdat::CellVariable<double> >& variable_quantity,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Compute minimum value with only x-direction as inhomogeneous direction.
         */
        std::vector<double> getMinQuantityWithInhomogeneousXDirection(
            HAMERS_SHARED_PTR<pdat::CellVariable<double> >& variable_quantity,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Compute maximum location within quantity bounds in x-direction.
         */
        double getMaxLocationWithinQuantityBoundsInXDirection(
            HAMERS_SHARED_PTR<pdat::CellVariable<double> >& variable_quantity,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double bound_lo,
            const double bound_hi) const;
        
        /*
         * Compute minimum location within quantity bounds in x-direction..
         */
        double getMinLocationWithinQuantityBoundsInXDirection(
            HAMERS_SHARED_PTR<pdat::CellVariable<double> >& variable_quantity,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double bound_lo,
            const double bound_hi) const;
        
        /*
         * Compute maximum location within quantity bounds in y-direction.
         */
        double getMaxLocationWithinQuantityBoundsInYDirection(
            HAMERS_SHARED_PTR<pdat::CellVariable<double> >& variable_quantity,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double bound_lo,
            const double bound_hi) const;
        
        /*
         * Compute minimum location within quantity bounds in y-direction.
         */
        double getMinLocationWithinQuantityBoundsInYDirection(
            HAMERS_SHARED_PTR<pdat::CellVariable<double> >& variable_quantity,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double bound_lo,
            const double bound_hi) const;
        
        /*
         * Compute maximum location within quantity bounds in z-direction.
         */
        double getMaxLocationWithinQuantityBoundsInZDirection(
            HAMERS_SHARED_PTR<pdat::CellVariable<double> >& variable_quantity,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double bound_lo,
            const double bound_hi) const;
        
        /*
         * Compute minimum location within quantity bounds in z-direction.
         */
        double getMinLocationWithinQuantityBoundsInZDirection(
            HAMERS_SHARED_PTR<pdat::CellVariable<double> >& variable_quantity,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double bound_lo,
            const double bound_hi) const;
        
    private:
        
};

#endif /* MPI_HELPER_MAX_MIN_HPP */
