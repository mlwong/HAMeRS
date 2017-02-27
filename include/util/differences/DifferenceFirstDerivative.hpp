#ifndef DIFFERENCE_FIRST_DERIVATIVE_HPP
#define DIFFERENCE_FIRST_DERIVATIVE_HPP

#include "util/differences/Difference.hpp"

class DifferenceFirstDerivative: public Difference
{
    public:
        DifferenceFirstDerivative(
            const std::string& object_name,
            const tbox::Dimension& dim);
        
        ~DifferenceFirstDerivative() {}
        
        /*
         * Compute the difference with the given cell data.
         */
        void
        computeDifference(
            hier::Patch& patch,
            boost::shared_ptr<pdat::CellData<double> > cell_data,
            boost::shared_ptr<pdat::CellData<double> > difference,
            int depth = 0);
        
        /*
         * Compute the difference and the local mean of the given cell data.
         */
        void
        computeDifferenceWithVariableLocalMean(
            hier::Patch& patch,
            boost::shared_ptr<pdat::CellData<double> > cell_data,
            boost::shared_ptr<pdat::CellData<double> > difference,
            boost::shared_ptr<pdat::CellData<double> > variable_local_mean,
            int depth = 0);
        
};

#endif /* DIFFERENCE_FIRST_DERIVATIVE_HPP */
