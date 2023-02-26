#ifndef FILTER_NONE_HPP
#define FILTER_NONE_HPP

#include "util/filters/Filter.hpp"

/*
 * This class is the implementation of a filter that has no filtering effect.
 */
class FilterNone: public Filter
{
    public:
        FilterNone(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const DIRECTION::TYPE& direction):
                Filter(
                object_name,
                dim,
                direction)
        {}
        
        ~FilterNone() {}
        
        /*
         * Apply filter to the given cell data.
         */
        void
        applyFilter(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& filtered_cell_data,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& cell_data,
            const int depth_filtered_cell_data,
            const int depth_cell_data)
        {
            const hier::Box empty_box(d_dim);
            applyFilter(
                filtered_cell_data,
                cell_data,
                depth_filtered_cell_data,
                depth_cell_data,
                empty_box);
        }
        
        /*
         * Apply filter to the given cell data.
         */
        void
        applyFilter(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& filtered_cell_data,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& cell_data,
            const int depth_filtered_cell_data,
            const int depth_cell_data,
            const hier::Box& domain);
        
};

#endif /* FILTER_NONE_HPP */
