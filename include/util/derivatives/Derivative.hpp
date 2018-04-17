#ifndef DERIVATIVE_HPP
#define DERIVATIVE_HPP

#include "HAMeRS_config.hpp"

#include "util/Directions.hpp"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/tbox/Dimension.h"

#include "boost/shared_ptr.hpp"
#include <string>
#include <vector>

using namespace SAMRAI;

class Derivative
{
    public:
        Derivative(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const DIRECTION::TYPE& direction,
            const int num_derivative_ghosts);
        
        virtual ~Derivative() {}
        
        /*
         * Compute the derivative with the given cell data.
         */
        virtual void
        computeDerivative(
            boost::shared_ptr<pdat::CellData<double> >& derivative,
            const boost::shared_ptr<pdat::CellData<double> >& data,
            const double dx,
            const int depth_derivative = 0,
            const int depth_data = 0) = 0;
        
        /*
         * Compute the derivative with the given cell data.
         */
        virtual void
        computeDerivative(
            boost::shared_ptr<pdat::CellData<double> >& derivative,
            const boost::shared_ptr<pdat::CellData<double> >& data,
            const double dx,
            const hier::Box& domain,
            const int depth_derivative = 0,
            const int depth_data = 0) = 0;
        
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
         * The direction to take derivative.
         */
        const DIRECTION::TYPE d_direction;
        
        /*
         * Number of ghost cells needed to take derivative.
         */
        hier::IntVector d_num_derivative_ghosts;
        
};

#endif /* DERIVATIVE_HPP */
