#ifndef DIFFERENCE_HPP
#define DIFFERENCE_HPP

#include "HAMeRS_config.hpp"

#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/tbox/Dimension.h"

#include "boost/shared_ptr.hpp"
#include <string>
#include <vector>

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
        getDifferenceNumberOfGhostCells(void) const
        {
            return d_num_difference_ghosts;
        }
        
        /*
         * Compute the difference with the given cell data.
         */
        virtual void
        computeDifference(
            boost::shared_ptr<pdat::CellData<double> >& difference,
            const boost::shared_ptr<pdat::CellData<double> >& cell_data,
            hier::Patch& patch,
            const int depth = 0) = 0;
        
        /*
         * Compute the difference and the local mean of the given cell data.
         */
        virtual void
        computeDifferenceWithVariableLocalMean(
            boost::shared_ptr<pdat::CellData<double> >& difference,
            boost::shared_ptr<pdat::CellData<double> >& variable_local_mean,
            const boost::shared_ptr<pdat::CellData<double> >& cell_data,
            hier::Patch& patch,
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
