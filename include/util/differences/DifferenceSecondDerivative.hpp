#ifndef DIFFERENCE_SECOND_DERIVATIVE_HPP
#define DIFFERENCE_SECOND_DERIVATIVE_HPP

#include "util/differences/Difference.hpp"

class DifferenceSecondDerivative: public Difference
{
    public:
        DifferenceSecondDerivative(
            const std::string& object_name,
            const tbox::Dimension& dim);
        
        ~DifferenceSecondDerivative() {}
        
        /*
         * Compute the difference with the given cell data.
         */
        void
        computeDifference(
            boost::shared_ptr<pdat::CellData<double> >& difference,
            const boost::shared_ptr<pdat::CellData<double> >& cell_data,
            hier::Patch& patch,
            const int depth = 0);
        
        /*
         * Compute the difference and the local mean of the given cell data.
         */
        void
        computeDifferenceWithVariableLocalMean(
            boost::shared_ptr<pdat::CellData<double> >& difference,
            boost::shared_ptr<pdat::CellData<double> >& variable_local_mean,
            const boost::shared_ptr<pdat::CellData<double> >& cell_data,
            hier::Patch& patch,
            const int depth = 0);
        
};

#endif /* DIFFERENCE_SECOND_DERIVATIVE_HPP */
