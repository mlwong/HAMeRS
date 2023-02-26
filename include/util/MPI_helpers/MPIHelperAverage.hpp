#ifndef MPI_HELPER_AVERAGE_HPP
#define MPI_HELPER_AVERAGE_HPP

#include "util/MPI_helpers/MPIHelper.hpp"

#include "SAMRAI/pdat/CellVariable.h"

class MPIHelperAverage: public MPIHelper
{
    public:
        MPIHelperAverage(
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
         * Compute averaged value over the entire domain.
         */
        Real getAveragedQuantity(
            HAMERS_SHARED_PTR<pdat::CellVariable<Real> >& variable_quantity,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Compute averaged value (on product of variables) over the entire domain.
         */
        Real getAveragedQuantity(
            std::vector<HAMERS_SHARED_PTR<pdat::CellVariable<Real> > >& variable_quantities,
            const std::vector<int>& component_indices,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Compute averaged value with only x-direction as inhomogeneous direction.
         */
        std::vector<Real> getAveragedQuantityWithInhomogeneousXDirection(
            HAMERS_SHARED_PTR<pdat::CellVariable<Real> >& variable_quantity,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Compute averaged value (on product of variables) with only x-direction as inhomogeneous direction.
         */
        std::vector<Real> getAveragedQuantityWithInhomogeneousXDirection(
            std::vector<HAMERS_SHARED_PTR<pdat::CellVariable<Real> > >& variable_quantities,
            const std::vector<int>& component_indices,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Compute averaged value with only y-direction as inhomogeneous direction.
         */
        std::vector<Real> getAveragedQuantityWithInhomogeneousYDirection(
            HAMERS_SHARED_PTR<pdat::CellVariable<Real> >& variable_quantity,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Compute averaged value (on product of variables) with only y-direction as inhomogeneous direction.
         */
        std::vector<Real> getAveragedQuantityWithInhomogeneousYDirection(
            std::vector<HAMERS_SHARED_PTR<pdat::CellVariable<Real> > >& variable_quantities,
            const std::vector<int>& component_indices,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Compute averaged value with only z-direction as inhomogeneous direction.
         */
        std::vector<Real> getAveragedQuantityWithInhomogeneousZDirection(
            HAMERS_SHARED_PTR<pdat::CellVariable<Real> >& variable_quantity,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Compute averaged value (on product of variables) with only z-direction as inhomogeneous direction.
         */
        std::vector<Real> getAveragedQuantityWithInhomogeneousZDirection(
            std::vector<HAMERS_SHARED_PTR<pdat::CellVariable<Real> > >& variable_quantities,
            const std::vector<int>& component_indices,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
    private:
        
};

#endif /* MPI_HELPER_AVERAGE_HPP */
