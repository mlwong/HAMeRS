#ifndef WAVELET_TRANSFORM_HPP
#define WAVELET_TRANSFORM_HPP

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/tbox/Dimension.h"

#include "boost/shared_ptr.hpp"
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
                d_num_wavelet_ghosts(hier::IntVector::getZero(d_dim)),
                d_num_ghosts(hier::IntVector::getZero(d_dim))
        {}
        
        /*
         * Get the number of ghost cells needed by the wavelet transformation.
         */
        hier::IntVector
        getWaveletTransformNumberOfGhostCells(void) const
        {
            return d_num_wavelet_ghosts;
        }
        
        /*
         * Set the number of ghost cells needed.
         */
        void
        setNumberOfGhostCells(const hier::IntVector& num_ghosts)
        {
            d_num_ghosts = num_ghosts;
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
            hier::Patch& patch,
            boost::shared_ptr<pdat::CellData<double> > cell_data,
            std::vector<boost::shared_ptr<pdat::CellData<double> > > wavelet_coeffs,
            int depth = 0) = 0;
    
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
         * Number of ghost cells for time-independent variables.
         */
        hier::IntVector d_num_ghosts;
        
        /*
         * Stencil width of the wavelet function.
         * p: number of points on the left.
         * q: number of points on the right.
         */
        int d_p;
        int d_q;
        
};

#endif /* WAVELET_TRANSFORM_HPP */
