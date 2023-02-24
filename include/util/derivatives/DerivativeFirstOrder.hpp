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
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& derivative,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data,
            const Real dx,
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
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& derivative,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data,
            const Real dx,
            const hier::Box& domain,
            const int depth_derivative = 0,
            const int depth_data = 0);
        
};

#endif /* DERIVATIVE_FIRST_ORDER_HPP */
