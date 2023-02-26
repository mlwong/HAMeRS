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
            std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& wavelet_coeffs,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& cell_data,
            hier::Patch& patch,
            const int depth = 0,
            const bool smooth_cell_data = false);
        
        /*
         * Perform the wavelet transformation and compute the local mean of the given cell data.
         */
        void
        computeWaveletCoefficientsWithVariableLocalMeans(
            std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& wavelet_coeffs,
            std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& variable_local_means,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& cell_data,
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
        HAMERS_SHARED_PTR<pdat::CellData<Real> >
        smoothCellData(
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& cell_data,
            hier::Patch& patch,
            const int depth = 0);
};

#endif /* WAVELET_TRANSFORM_HARTEN_HPP */
