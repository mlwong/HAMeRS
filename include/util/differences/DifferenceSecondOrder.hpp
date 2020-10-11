#ifndef DIFFERENCE_SECOND_ORDER_HPP
#define DIFFERENCE_SECOND_ORDER_HPP

#include "util/differences/Difference.hpp"

class DifferenceSecondOrder: public Difference
{
    public:
        DifferenceSecondOrder(
            const std::string& object_name,
            const tbox::Dimension& dim);
        
        ~DifferenceSecondOrder() {}
        
        /*
         * Compute the difference with the given cell data.
         */
        void
        computeDifference(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& difference,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& cell_data,
            const int depth = 0)
        {
            const hier::Box empty_box(d_dim);
            computeDifference(
                difference,
                cell_data,
                empty_box,
                depth);
        }
        
        /*
         * Compute the difference with the given cell data.
         */
        void
        computeDifference(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& difference,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& cell_data,
            const hier::Box& domain,
            const int depth = 0);
        
        /*
         * Compute the difference and the local mean of the given cell data.
         */
        void
        computeDifferenceWithVariableLocalMean(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& difference,
            HAMERS_SHARED_PTR<pdat::CellData<double> >& variable_local_mean,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& cell_data,
            const int depth = 0)
        {
            const hier::Box empty_box(d_dim);
            computeDifferenceWithVariableLocalMean(
                difference,
                variable_local_mean,
                cell_data,
                empty_box,
                depth);
        }
        
        /*
         * Compute the difference and the local mean of the given cell data.
         */
        void
        computeDifferenceWithVariableLocalMean(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& difference,
            HAMERS_SHARED_PTR<pdat::CellData<double> >& variable_local_mean,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& cell_data,
            const hier::Box& domain,
            const int depth = 0);
        
};

#endif /* DIFFERENCE_SECOND_ORDER_HPP */
