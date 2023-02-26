#ifndef WAVELET_TRANSFORM_HPP
#define WAVELET_TRANSFORM_HPP

#include "HAMeRS_config.hpp"

#include "HAMeRS_memory.hpp"

#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/tbox/Dimension.h"

#include <string>
#include <vector>

using namespace SAMRAI;

class WaveletTransform
{
    public:
        WaveletTransform(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const int num_level):
                d_object_name(object_name),
                d_dim(dim),
                d_num_level(num_level),
                d_num_wavelet_ghosts(hier::IntVector::getZero(d_dim))
        {}
        
        virtual ~WaveletTransform() {}
        
        /*
         * Get the number of ghost cells needed by the wavelet transformation.
         */
        hier::IntVector
        getWaveletTransformNumberOfGhostCells(void) const
        {
            return d_num_wavelet_ghosts;
        }
        
        /*
         * Get the left stencil width.
         */
        int
        getLeftStencilWidth() const
        {
            return d_p;
        }
        
        /*
         * Get the right stencil width.
         */
        int
        getRightStencilWidth() const
        {
            return d_q;
        }
        
        /*
         * Perform the wavelet transformation on the cell data.
         */
        virtual void
        computeWaveletCoefficients(
            std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& wavelet_coeffs,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& cell_data,
            hier::Patch& patch,
            const int depth = 0,
            const bool smooth_cell_data = false) = 0;
        
        /*
         * Perform the wavelet transformation and compute the local mean of the given cell data.
         */
        virtual void
        computeWaveletCoefficientsWithVariableLocalMeans(
            std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& wavelet_coeffs,
            std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& variable_local_means,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& cell_data,
            hier::Patch& patch,
            const int depth = 0,
            const bool smooth_cell_data = false) = 0;
        
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
         * Number of wavelet levels.
         */
        const int d_num_level;
        
        /*
         * Number of ghost cells needed by the wavelet scheme.
         */
        hier::IntVector d_num_wavelet_ghosts;
        
        /*
         * Stencil width of the wavelet function.
         * p: number of points on the left.
         * q: number of points on the right.
         */
        int d_p;
        int d_q;
        
};

#endif /* WAVELET_TRANSFORM_HPP */
