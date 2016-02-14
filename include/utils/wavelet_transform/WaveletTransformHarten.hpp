#ifndef WAVELET_TRANSFORM_HARTEN_HPP
#define WAVELET_TRANSFORM_HARTEN_HPP

#include "utils/wavelet_transform/WaveletTransform.hpp"

class WaveletTransformHarten: public WaveletTransform
{
    public:
        WaveletTransformHarten(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const int num_level,
            const int num_vanishing_moments);
        
        /*
         * Perform the wavelet transformation on the given cell data.
         */
        void
        computeWaveletCoefficients(
            hier::Patch& patch,
            boost::shared_ptr<pdat::CellData<double> > cell_data,
            std::vector<boost::shared_ptr<pdat::CellData<double> > > wavelet_coeffs,
            int depth = 0);
    
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
            hier::Patch& patch,
            boost::shared_ptr<pdat::CellData<double> > cell_data,
            int depth = 0);
};

#endif /* WAVELET_TRANSFORM_HARTEN_HPP */
