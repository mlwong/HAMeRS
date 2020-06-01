#ifndef FILTER_TRUNCATED_GAUSSIAN_HPP
#define FILTER_TRUNCATED_GAUSSIAN_HPP

#include "util/filters/Filter.hpp"

class FilterTruncatedGaussian: public Filter
{
    public:
        FilterTruncatedGaussian(
            const std::string& object_name,
            const tbox::Dimension& dim);
        
        ~FilterTruncatedGaussian() {}
        
        /*
         * Apply filter to the given cell data.
         */
        void
        applyFilter(
            boost::shared_ptr<pdat::CellData<double> >& filtered_cell_data,
            const boost::shared_ptr<pdat::CellData<double> >& cell_data,
            const int depth = 0)
        {
            const hier::Box empty_box(d_dim);
            applyFilter(
                filtered_cell_data,
                cell_data,
                empty_box,
                depth);
        }
        
        /*
         * Apply filter to the given cell data.
         */
        void
        applyFilter(
            boost::shared_ptr<pdat::CellData<double> >& filtered_cell_data,
            const boost::shared_ptr<pdat::CellData<double> >& cell_data,
            const hier::Box& domain,
            const int depth = 0);
        
};

#endif /* FILTER_TRUNCATED_GAUSSIAN_HPP */
