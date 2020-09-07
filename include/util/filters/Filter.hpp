#ifndef FILTER_HPP
#define FILTER_HPP

#include "HAMeRS_config.hpp"

#include "util/Directions.hpp"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/tbox/Dimension.h"

#include "boost/shared_ptr.hpp"
#include <string>

using namespace SAMRAI;

class Filter
{
    public:
        Filter(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const DIRECTION::TYPE& direction);
        
        virtual ~Filter() {}
        
        /*
         * Get the number of ghost cells needed by the filtering operation.
         */
        hier::IntVector
        getFilterNumberOfGhostCells(void) const
        {
            return d_num_filter_ghosts;
        }
        
        /*
         * Apply filter to the given cell data.
         */
        virtual void
        applyFilter(
            boost::shared_ptr<pdat::CellData<double> >& filtered_cell_data,
            const boost::shared_ptr<pdat::CellData<double> >& cell_data,
            const int depth_filtered_cell_data,
            const int depth_cell_data) = 0;
        
        /*
         * Apply filter to the given cell data.
         */
        virtual void
        applyFilter(
            boost::shared_ptr<pdat::CellData<double> >& filtered_cell_data,
            const boost::shared_ptr<pdat::CellData<double> >& cell_data,
            const int depth_filtered_cell_data,
            const int depth_cell_data,
            const hier::Box& domain) = 0;
        
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
         * The direction to apply derivative.
         */
        const DIRECTION::TYPE d_direction;
        
        /*
         * Number of ghost cells needed by the filter.
         */
        hier::IntVector d_num_filter_ghosts;
        
};

#endif /* FILTER_HPP */
