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
            boost::shared_ptr<pdat::CellData<double> >& difference,
            const boost::shared_ptr<pdat::CellData<double> >& cell_data,
            const hier::Box& domain,
            const int depth = 0);
        
        /*
         * Compute the difference and the local mean of the given cell data.
         */
        void
        computeDifferenceWithVariableLocalMean(
            boost::shared_ptr<pdat::CellData<double> >& difference,
            boost::shared_ptr<pdat::CellData<double> >& variable_local_mean,
            const boost::shared_ptr<pdat::CellData<double> >& cell_data,
            const hier::Box& domain,
            const int depth = 0);
        
};

#endif /* DIFFERENCE_SECOND_ORDER_HPP */
