#ifndef DERIVATIVE_FIRST_ORDER_HPP
#define DERIVATIVE_FIRST_ORDER_HPP

#include "util/derivatives/Derivative.hpp"

class DerivativeFirstOrder: public Derivative
{
    public:
        DerivativeFirstOrder(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const DIRECTION::TYPE& direction,
            const int num_derivative_ghosts);
        
        ~DerivativeFirstOrder() {}
        
        /*
         * Compute the derivative with the given cell data.
         */
        void
        computeDerivative(
            boost::shared_ptr<pdat::CellData<double> >& derivative,
            const boost::shared_ptr<pdat::CellData<double> >& cell_data,
            const double dx,
            const int depth = 0)
        {
            const hier::Box empty_box(d_dim);
            computeDerivative(
                derivative,
                cell_data,
                dx,
                empty_box,
                depth);
        }
        
        /*
         * Compute the derivative with the given cell data.
         */
        void
        computeDerivative(
            boost::shared_ptr<pdat::CellData<double> >& derivative,
            const boost::shared_ptr<pdat::CellData<double> >& cell_data,
            const double dx,
            const hier::Box& domain,
            const int depth = 0);
        
};

#endif /* DERIVATIVE_FIRST_ORDER_HPP */
