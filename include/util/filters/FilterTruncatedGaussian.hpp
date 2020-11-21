#ifndef FILTER_TRUNCATED_GAUSSIAN_HPP
#define FILTER_TRUNCATED_GAUSSIAN_HPP

#include "util/filters/Filter.hpp"

/*
 * This class is the implementation of the truncated Gaussian filter by Cook, Andrew W.,
 * "Artificial fluid properties for large-eddy simulation of compressible turbulent mixing.",
 * Physics of fluids 19.5 (2007): 055103.
 * The cutoff of the low-pass filter is around 4 Delta, which is the grid spaceing.
 */
class FilterTruncatedGaussian: public Filter
{
    public:
        FilterTruncatedGaussian(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const DIRECTION::TYPE& direction);
        
        ~FilterTruncatedGaussian() {}
        
        /*
         * Apply filter to the given cell data.
         */
        void
        applyFilter(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& filtered_cell_data,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& cell_data,
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
            HAMERS_SHARED_PTR<pdat::CellData<double> >& filtered_cell_data,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& cell_data,
            const int depth_filtered_cell_data,
            const int depth_cell_data,
            const hier::Box& domain);
        
    private:
        /*
         * Coefficients for truncated Gaussian filter.
         */
        const double a_G, b_G, c_G, d_G, e_G;
        
};

#endif /* FILTER_TRUNCATED_GAUSSIAN_HPP */
