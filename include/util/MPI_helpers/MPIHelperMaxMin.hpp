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
         * Compute maximum value over the entire domain.
         */
        Real getMaxQuantity(
            HAMERS_SHARED_PTR<pdat::CellVariable<Real> >& variable_quantity,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Compute minimum value over the entire domain.
         */
        Real getMinQuantity(
            HAMERS_SHARED_PTR<pdat::CellVariable<Real> >& variable_quantity,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Compute maximum value with only x-direction as inhomogeneous direction.
         */
        std::vector<Real> getMaxQuantityWithInhomogeneousXDirection(
            HAMERS_SHARED_PTR<pdat::CellVariable<Real> >& variable_quantity,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Compute minimum value with only x-direction as inhomogeneous direction.
         */
        std::vector<Real> getMinQuantityWithInhomogeneousXDirection(
            HAMERS_SHARED_PTR<pdat::CellVariable<Real> >& variable_quantity,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Compute maximum value with only y-direction as inhomogeneous direction.
         */
        std::vector<Real> getMaxQuantityWithInhomogeneousYDirection(
            HAMERS_SHARED_PTR<pdat::CellVariable<Real> >& variable_quantity,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Compute minimum value with only y-direction as inhomogeneous direction.
         */
        std::vector<Real> getMinQuantityWithInhomogeneousYDirection(
            HAMERS_SHARED_PTR<pdat::CellVariable<Real> >& variable_quantity,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Compute maximum value with only z-direction as inhomogeneous direction.
         */
        std::vector<Real> getMaxQuantityWithInhomogeneousZDirection(
            HAMERS_SHARED_PTR<pdat::CellVariable<Real> >& variable_quantity,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Compute minimum value with only z-direction as inhomogeneous direction.
         */
        std::vector<Real> getMinQuantityWithInhomogeneousZDirection(
            HAMERS_SHARED_PTR<pdat::CellVariable<Real> >& variable_quantity,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Compute maximum location within quantity bounds in x-direction.
         */
        Real getMaxLocationWithinQuantityBoundsInXDirection(
            HAMERS_SHARED_PTR<pdat::CellVariable<Real> >& variable_quantity,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const Real bound_lo,
            const Real bound_hi) const;
        
        /*
         * Compute minimum location within quantity bounds in x-direction..
         */
        Real getMinLocationWithinQuantityBoundsInXDirection(
            HAMERS_SHARED_PTR<pdat::CellVariable<Real> >& variable_quantity,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const Real bound_lo,
            const Real bound_hi) const;
        
        /*
         * Compute maximum location within quantity bounds in y-direction.
         */
        Real getMaxLocationWithinQuantityBoundsInYDirection(
            HAMERS_SHARED_PTR<pdat::CellVariable<Real> >& variable_quantity,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const Real bound_lo,
            const Real bound_hi) const;
        
        /*
         * Compute minimum location within quantity bounds in y-direction.
         */
        Real getMinLocationWithinQuantityBoundsInYDirection(
            HAMERS_SHARED_PTR<pdat::CellVariable<Real> >& variable_quantity,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const Real bound_lo,
            const Real bound_hi) const;
        
        /*
         * Compute maximum location within quantity bounds in z-direction.
         */
        Real getMaxLocationWithinQuantityBoundsInZDirection(
            HAMERS_SHARED_PTR<pdat::CellVariable<Real> >& variable_quantity,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const Real bound_lo,
            const Real bound_hi) const;
        
        /*
         * Compute minimum location within quantity bounds in z-direction.
         */
        Real getMinLocationWithinQuantityBoundsInZDirection(
            HAMERS_SHARED_PTR<pdat::CellVariable<Real> >& variable_quantity,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const Real bound_lo,
            const Real bound_hi) const;
        
    private:
        
};

#endif /* MPI_HELPER_MAX_MIN_HPP */
