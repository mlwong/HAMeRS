#ifndef DERIVATIVE_SECOND_ORDER_HPP
#define DERIVATIVE_SECOND_ORDER_HPP

#include "util/derivatives/Derivative.hpp"

class DerivativeSecondOrder: public Derivative
{
    public:
        DerivativeSecondOrder(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const DIRECTION::TYPE& direction,
            const int num_derivative_ghosts);
        
        ~DerivativeSecondOrder() {}
        
        /*
         * Compute the derivative with the given cell data.
         */
        void
        computeDerivative(
            boost::shared_ptr<pdat::CellData<double> >& derivative,
            const boost::shared_ptr<pdat::CellData<double> >& data,
            const double dx,
            const int depth_derivative = 0,
            const int depth_data = 0)
        {
            const hier::Box empty_box(d_dim);
            computeDerivative(
                derivative,
                data,
                dx,
                empty_box,
                depth_derivative,
                depth_data);
        }
        
        /*
         * Compute the derivative with the given cell data.
         */
        void
        computeDerivative(
            boost::shared_ptr<pdat::CellData<double> >& derivative,
            const boost::shared_ptr<pdat::CellData<double> >& data,
            const double dx,
            const hier::Box& domain,
            const int depth_derivative = 0,
            const int depth_data = 0);
        
};

#endif /* DERIVATIVE_SECOND_ORDER_HPP */
