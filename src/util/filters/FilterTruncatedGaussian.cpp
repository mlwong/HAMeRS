#include "util/filters/FilterTruncatedGaussian.hpp"

FilterTruncatedGaussian::FilterTruncatedGaussian(
    const std::string& object_name,
    const tbox::Dimension& dim):
        Filter(
            object_name,
            dim)
{
    d_num_filter_ghosts = hier::IntVector::getOne(d_dim)*4;
}

/*
 * Apply filter to the given cell data.
 */
void
FilterTruncatedGaussian::applyFilter(
    boost::shared_ptr<pdat::CellData<double> >& filtered_cell_data,
    const boost::shared_ptr<pdat::CellData<double> >& cell_data,
    const hier::Box& domain,
    const int depth)
{
}
