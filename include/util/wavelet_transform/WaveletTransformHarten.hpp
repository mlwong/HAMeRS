#ifndef WAVELET_TRANSFORM_HARTEN_HPP
#define WAVELET_TRANSFORM_HARTEN_HPP

#include "util/wavelet_transform/WaveletTransform.hpp"

class WaveletTransformHarten: public WaveletTransform
{
    public:
        WaveletTransformHarten(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const int num_level,
            const int num_vanishing_moments);
        
        ~WaveletTransformHarten() {}
        
        /*
         * Perform the wavelet transformation on the given cell data.
         */
        void
        computeWaveletCoefficients(
            std::vector<boost::shared_ptr<pdat::CellData<double> > >& wavelet_coeffs,
            const boost::shared_ptr<pdat::CellData<double> >& cell_data,
            hier::Patch& patch,
            const int depth = 0,
            const bool smooth_cell_data = false);
        
        /*
         * Perform the wavelet transformation and compute the local mean of the given cell data.
         */
        void
        computeWaveletCoefficientsWithVariableLocalMeans(
            std::vector<boost::shared_ptr<pdat::CellData<double> > >& wavelet_coeffs,
            std::vector<boost::shared_ptr<pdat::CellData<double> > >& variable_local_means,
            const boost::shared_ptr<pdat::CellData<double> >& cell_data,
            hier::Patch& patch,
            const int depth = 0,
            const bool smooth_cell_data = false);
        
    private:
        /*
         * The number of vanishing moments.
         */
        const int d_k;
        
        /*
         * Smooth the given cell data in different directions.
         */
        boost::shared_ptr<pdat::CellData<double> >
        smoothCellData(
            const boost::shared_ptr<pdat::CellData<double> >& cell_data,
            hier::Patch& patch,
            const int depth = 0);
};

#endif /* WAVELET_TRANSFORM_HARTEN_HPP */
