#ifndef DIFFERENCE_HPP
#define DIFFERENCE_HPP

#include "HAMeRS_config.hpp"

#include "HAMeRS_memory.hpp"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/tbox/Dimension.h"

#include <string>

using namespace SAMRAI;

class Difference
{
    public:
        Difference(
            const std::string& object_name,
            const tbox::Dimension& dim):
                d_object_name(object_name),
                d_dim(dim),
                d_num_difference_ghosts(hier::IntVector::getZero(d_dim))
        {}
        
        virtual ~Difference() {}
        
        /*
         * Get the number of ghost cells needed by the difference operation.
         */
        hier::IntVector
        getDifferenceNumberOfGhostCells() const
        {
            return d_num_difference_ghosts;
        }
        
        /*
         * Compute the difference with the given cell data.
         */
        virtual void
        computeDifference(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& difference,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& cell_data,
            const int depth = 0) = 0;
        
        /*
         * Compute the difference with the given cell data.
         */
        virtual void
        computeDifference(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& difference,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& cell_data,
            const hier::Box& domain,
            const int depth = 0) = 0;
        
        /*
         * Compute the difference and the local mean of the given cell data.
         */
        virtual void
        computeDifferenceWithVariableLocalMean(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& difference,
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& variable_local_mean,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& cell_data,
            const int depth = 0) = 0;
        
        /*
         * Compute the difference and the local mean of the given cell data.
         */
        virtual void
        computeDifferenceWithVariableLocalMean(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& difference,
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& variable_local_mean,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& cell_data,
            const hier::Box& domain,
            const int depth = 0) = 0;
        
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
         * Number of ghost cells needed by the gradient sensor.
         */
        hier::IntVector d_num_difference_ghosts;
        
};

#endif /* DIFFERENCE_HPP */
