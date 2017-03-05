#include "util/wavelet_transform/WaveletTransformHarten.hpp"

#include "SAMRAI/geom/CartesianPatchGeometry.h"

#include <cfloat>

WaveletTransformHarten::WaveletTransformHarten(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const int num_level,
    const int num_vanishing_moments):
        WaveletTransform(
            object_name,
            dim,
            num_level),
        d_k(num_vanishing_moments)
{
    switch (num_vanishing_moments)
    {
        case 2:
        {
            d_p = 1;
            d_q = 1;
            
            break;
        }
        case 4:
        {
            d_p = 2;
            d_q = 2;
            
            break;
        }
        default:
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "num_vanishing_moments = "
                << d_k
                << " not supported. \n"
                << "Only 2 or 4 vanishing moments are allowed."
                << std::endl);
        }
    }
    
    if (num_level < 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Only number of wavelet levels larger than 1 is allowed. \n"
            << "num_level = "
            << num_level
            << " is provided."
            << std::endl);
    }
    
    int num_wavelet_ghosts = 0;
    
    for (int li = 0; li < num_level; li++)
    {
        num_wavelet_ghosts += fmax(d_p, d_q)*pow(2, li);
    }
    
    num_wavelet_ghosts += fmax(d_p, d_q)*pow(2, num_level);
    
    d_num_wavelet_ghosts = hier::IntVector::getOne(d_dim)*num_wavelet_ghosts;
}


/*
 * Perform the wavelet transformation on the given cell data.
 */
void
WaveletTransformHarten::computeWaveletCoefficients(
    hier::Patch& patch,
    boost::shared_ptr<pdat::CellData<double> > cell_data,
    std::vector<boost::shared_ptr<pdat::CellData<double> > >& wavelet_coeffs,
    int depth,
    bool smooth_cell_data)
{
    // Declare an empty vector.
    std::vector<boost::shared_ptr<pdat::CellData<double> > > variable_local_means;
    
    computeWaveletCoefficientsWithVariableLocalMeans(
        patch,
        cell_data,
        wavelet_coeffs,
        variable_local_means,
        depth,
        smooth_cell_data);
}


/*
 * Perform the wavelet transformation and compute the local mean of the given cell data.
 */
void
WaveletTransformHarten::computeWaveletCoefficientsWithVariableLocalMeans(
    hier::Patch& patch,
    boost::shared_ptr<pdat::CellData<double> > cell_data,
    std::vector<boost::shared_ptr<pdat::CellData<double> > >& wavelet_coeffs,
    std::vector<boost::shared_ptr<pdat::CellData<double> > >& variable_local_means,
    int depth,
    bool smooth_cell_data)
{
    TBOX_ASSERT(wavelet_coeffs.size() == variable_local_means.size());
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the number of ghost cells of the cell data, wavelet coefficients and local means.
    const hier::IntVector num_ghosts_cell_data = cell_data->getGhostCellWidth();
    const hier::IntVector num_ghosts_wavelet_coeffs = wavelet_coeffs[0]->getGhostCellWidth();
    
    TBOX_ASSERT(num_ghosts_cell_data >= num_ghosts_wavelet_coeffs);
    TBOX_ASSERT(num_ghosts_wavelet_coeffs >= d_num_wavelet_ghosts);
    
    // Get the dimensions of boxes that cover interior of patch plus ghost cells.
    const hier::Box ghost_box_cell_data = cell_data->getGhostBox();
    const hier::IntVector ghostcell_dims_cell_data = ghost_box_cell_data.numberCells();
    
    const hier::Box ghost_box_wavelet_coeffs = wavelet_coeffs[0]->getGhostBox();
    const hier::IntVector ghostcell_dims_wavelet_coeffs = ghost_box_wavelet_coeffs.numberCells();
    
    for (int li = 0; li < static_cast<int>(wavelet_coeffs.size()); li++)
    {
        TBOX_ASSERT(wavelet_coeffs[li]->getGhostBox().numberCells() == ghostcell_dims_wavelet_coeffs);
    }
    
    // Determine whether local means at different levels are required to be computed.
    bool compute_variable_local_means = false;
    if (!variable_local_means.empty())
    {
        compute_variable_local_means = true;
        
        for (int li = 0; li < static_cast<int>(wavelet_coeffs.size()); li++)
        {
            TBOX_ASSERT(variable_local_means[li]->getGhostBox().numberCells() == ghostcell_dims_wavelet_coeffs);
        }
    }
    
    // Get the pointer to the desired depth component of the given cell data.
    double* f = cell_data->getPointer(depth);
    
    // Get the pointer to the wavelet coefficients at different levels.
    std::vector<double*> w;
    for (int li = 0; li < d_num_level; li++)
    {
        w.push_back(wavelet_coeffs[li]->getPointer(0));
        wavelet_coeffs[li]->fillAll(0.0);
    }
    
    // Get the pointer to the local means at different levels if it is required to compute them.
    std::vector<double*> f_mean;
    if (compute_variable_local_means)
    {
        for (int li = 0; li < d_num_level; li++)
        {
            f_mean.push_back(variable_local_means[li]->getPointer(0));
        }
    }
    
    /*
     * Smooth the depth component of the given cell data.
     */
    boost::shared_ptr<pdat::CellData<double> > smoothed_cell_data;
    
    if (smooth_cell_data)
    {
        smoothed_cell_data = smoothCellData(
            patch,
            cell_data,
            depth);
    }
    
    if (d_dim == tbox::Dimension(1))
    {
        const int interior_dim_0 = interior_dims[0];
        
        const int num_ghosts_0_cell_data = num_ghosts_cell_data[0];
        const int num_ghosts_0_wavelet_coeffs = num_ghosts_wavelet_coeffs[0];
        
        // Allocate scaling function coefficients.
        std::vector<boost::shared_ptr<pdat::CellData<double> > > scaling_coeffs_x;
        for (int li = 0; li < d_num_level; li++)
        {
            scaling_coeffs_x.push_back(boost::make_shared<pdat::CellData<double> >(
                interior_box, 1, num_ghosts_wavelet_coeffs));
        }
        
        // Get the pointers to the scaling function coefficients and initialize them with zero values.
        std::vector<double*> f_x;
        for (int li = 0; li < d_num_level; li++)
        {
            f_x.push_back(scaling_coeffs_x[li]->getPointer(0));
            scaling_coeffs_x[li]->fillAll(0.0);
        }
        
        /*
         * Compute scaling and wavelet coefficients at first level.
         */
        
        switch (d_k)
        {
            case 2:
            {
                /*
                 * Compute scaling and wavelet coefficients in the x-direction.
                 */
                
                // Get the pointer to the smoothed cell data in the x-direction.
                if (smooth_cell_data)
                {
                    f = smoothed_cell_data->getPointer(0);
                }
                else
                {
                    f = cell_data->getPointer(depth);
                }
                
                // Compute the starting and ending indices.
                const int start_index_x_i = -d_num_wavelet_ghosts[0] + 1;
                const int end_index_x_i = interior_dim_0 + d_num_wavelet_ghosts[0] - 1;
                
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = start_index_x_i; i < end_index_x_i; i++)
                {
                    // Compute the linear indices.
                    const int idx = i + num_ghosts_0_wavelet_coeffs;
                    
                    const int idx_cell_data     = i + num_ghosts_0_cell_data;
                    const int idx_cell_data_x_L = i - 1 + num_ghosts_0_cell_data;
                    const int idx_cell_data_x_R = i + 1 + num_ghosts_0_cell_data;
                    
                    f_x[0][idx] = 0.5*(f[idx_cell_data_x_L] + f[idx_cell_data_x_R]);
                    
                    w[0][idx] = fabs(-0.5*(f[idx_cell_data_x_L] - 2.0*f[idx_cell_data] +
                        f[idx_cell_data_x_R]));
                }
                
                if (compute_variable_local_means)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx = i + num_ghosts_0_wavelet_coeffs;
                        
                        const int idx_cell_data     = i + num_ghosts_0_cell_data;
                        const int idx_cell_data_x_L = i - 1 + num_ghosts_0_cell_data;
                        const int idx_cell_data_x_R = i + 1 + num_ghosts_0_cell_data;
                        
                        f_mean[0][idx] = 0.5*(f[idx_cell_data_x_L] + 2.0*f[idx_cell_data] +
                            f[idx_cell_data_x_R]);
                    }
                }
                
                break;
            }
            case 4:
            {
                /*
                 * Compute scaling and wavelet coefficients in the x-direction.
                 */
                
                // Get the pointer to the smoothed cell data in the x-direction.
                if (smooth_cell_data)
                {
                    f = smoothed_cell_data->getPointer(0);
                }
                else
                {
                    f = cell_data->getPointer(depth);
                }
                
                // Compute the starting and ending indices.
                const int start_index_x_i = -d_num_wavelet_ghosts[0] + 2;
                const int end_index_x_i = interior_dim_0 + d_num_wavelet_ghosts[0] - 2;
                
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = start_index_x_i; i < end_index_x_i; i++)
                {
                    // Compute the linear indices.
                    const int idx = i + num_ghosts_0_wavelet_coeffs;
                    
                    const int idx_cell_data      = i + num_ghosts_0_cell_data;
                    const int idx_cell_data_x_LL = i - 2 + num_ghosts_0_cell_data;
                    const int idx_cell_data_x_L  = i - 1 + num_ghosts_0_cell_data;
                    const int idx_cell_data_x_R  = i + 1 + num_ghosts_0_cell_data;
                    const int idx_cell_data_x_RR = i + 2 + num_ghosts_0_cell_data;
                    
                    f_x[0][idx] = 1.0/6.0*(-f[idx_cell_data_x_LL] + 4.0*f[idx_cell_data_x_L] +
                        4.0*f[idx_cell_data_x_R] - f[idx_cell_data_x_RR]);
                    
                    w[0][idx] = fabs(1.0/6.0*(f[idx_cell_data_x_LL] - 4.0*f[idx_cell_data_x_L] +
                        6.0*f[idx_cell_data] - 4.0*f[idx_cell_data_x_R] + f[idx_cell_data_x_RR]));
                }
                
                if (compute_variable_local_means)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx = i + num_ghosts_0_wavelet_coeffs;
                        
                        const int idx_cell_data      = i + num_ghosts_0_cell_data;
                        const int idx_cell_data_x_LL = i - 2 + num_ghosts_0_cell_data;
                        const int idx_cell_data_x_L  = i - 1 + num_ghosts_0_cell_data;
                        const int idx_cell_data_x_R  = i + 1 + num_ghosts_0_cell_data;
                        const int idx_cell_data_x_RR = i + 2 + num_ghosts_0_cell_data;
                        
                        f_mean[0][idx] = 1.0/6.0*(f[idx_cell_data_x_LL] + 4.0*f[idx_cell_data_x_L] +
                            6.0*f[idx_cell_data] + 4.0*f[idx_cell_data_x_R] + f[idx_cell_data_x_RR]);
                    }
                }
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "number of vanishing moments = "
                    << d_k
                    << " not supported."
                    << std::endl);
            }
        }
        
        /*
         * Compute scaling and wavelet coefficients at higher levels.
         */
        
        for (int li = 1; li < d_num_level; li++)
        {
            const int offset = pow(2, li);
            
            switch (d_k)
            {
                case 2:
                {
                    /*
                     * Compute scaling and wavelet coefficients in the x-direction.
                     */
                    
                    // Compute the starting and ending indices.
                    const int start_index_x_i = -d_num_wavelet_ghosts[0] + offset;
                    const int end_index_x_i = interior_dim_0 + d_num_wavelet_ghosts[0] - offset;
                    
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = start_index_x_i; i < end_index_x_i; i++)
                    {
                        // Compute the linear indices.
                        const int idx     = i + num_ghosts_0_wavelet_coeffs;
                        const int idx_x_L = i - offset + num_ghosts_0_wavelet_coeffs;
                        const int idx_x_R = i + offset + num_ghosts_0_wavelet_coeffs;
                        
                        f_x[li][idx] = 0.5*(f_x[li-1][idx_x_L] + f_x[li-1][idx_x_R]);
                        
                        w[li][idx]   = fabs(-0.5*(f_x[li-1][idx_x_L] - 2.0*f_x[li-1][idx] +
                            f_x[li-1][idx_x_R]));
                    }
                    
                    if (compute_variable_local_means)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx     = i + num_ghosts_0_wavelet_coeffs;
                            const int idx_x_L = i - offset + num_ghosts_0_wavelet_coeffs;
                            const int idx_x_R = i + offset + num_ghosts_0_wavelet_coeffs;
                            
                            f_mean[li][idx] = 0.5*(f_x[li-1][idx_x_L] + 2.0*f_x[li-1][idx] +
                                f_x[li-1][idx_x_R]);
                        }
                    }
                    
                    break;
                }
                case 4:
                {
                    /*
                     * Compute scaling and wavelet coefficients in the x-direction.
                     */
                    
                    // Compute the starting and ending indices.
                    const int start_index_x_i = -d_num_wavelet_ghosts[0] + 2*offset;
                    const int end_index_x_i = interior_dim_0 + d_num_wavelet_ghosts[0] - 2*offset;
                    
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = start_index_x_i; i < end_index_x_i; i++)
                    {
                        // Compute indices.
                        const int idx      = i + num_ghosts_0_wavelet_coeffs;
                        const int idx_x_LL = i - 2*offset + num_ghosts_0_wavelet_coeffs;
                        const int idx_x_L  = i - offset + num_ghosts_0_wavelet_coeffs;
                        const int idx_x_R  = i + offset + num_ghosts_0_wavelet_coeffs;
                        const int idx_x_RR = i + 2*offset + num_ghosts_0_wavelet_coeffs;
                        
                        f_x[li][idx] = 1.0/6.0*(-f_x[li-1][idx_x_LL] + 4.0*f_x[li-1][idx_x_L] +
                            4.0*f_x[li-1][idx_x_R] - f_x[li-1][idx_x_RR]);
                        
                        w[li][idx] = fabs(1.0/6.0*(f_x[li-1][idx_x_LL] - 4.0*f_x[li-1][idx_x_L] +
                            6.0*f_x[li-1][idx] - 4.0*f_x[li-1][idx_x_R] + f_x[li-1][idx_x_RR]));
                    }
                    
                    if (compute_variable_local_means)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute indices.
                            const int idx      = i + num_ghosts_0_wavelet_coeffs;
                            const int idx_x_LL = i - 2*offset + num_ghosts_0_wavelet_coeffs;
                            const int idx_x_L  = i - offset + num_ghosts_0_wavelet_coeffs;
                            const int idx_x_R  = i + offset + num_ghosts_0_wavelet_coeffs;
                            const int idx_x_RR = i + 2*offset + num_ghosts_0_wavelet_coeffs;
                            
                            f_mean[li][idx] = 1.0/6.0*(f_x[li-1][idx_x_LL] + 4.0*f_x[li-1][idx_x_L] +
                                6.0*f_x[li-1][idx] + 4.0*f_x[li-1][idx_x_R] + f_x[li-1][idx_x_RR]);
                        }
                    }
                    
                    break;
                }
                default:
                {
                    TBOX_ERROR(d_object_name
                        << ": "
                        << "number of vanishing moments = "
                        << d_k
                        << " not supported."
                        << std::endl);
                }
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        
        const int num_ghosts_0_cell_data = num_ghosts_cell_data[0];
        const int num_ghosts_1_cell_data = num_ghosts_cell_data[1];
        const int ghostcell_dim_0_cell_data = ghostcell_dims_cell_data[0];
        
        const int num_ghosts_0_wavelet_coeffs = num_ghosts_wavelet_coeffs[0];
        const int num_ghosts_1_wavelet_coeffs = num_ghosts_wavelet_coeffs[1];
        const int ghostcell_dim_0_wavelet_coeffs = ghostcell_dims_wavelet_coeffs[0];
        
        // Allocate wavelet coefficients in different dimensions.
        std::vector<boost::shared_ptr<pdat::CellData<double> > > wavelet_coeffs_x;
        std::vector<boost::shared_ptr<pdat::CellData<double> > > wavelet_coeffs_y;
        std::vector<boost::shared_ptr<pdat::CellData<double> > > scaling_coeffs_x;
        std::vector<boost::shared_ptr<pdat::CellData<double> > > scaling_coeffs_y;
        for (int li = 0; li < d_num_level; li++)
        {
            wavelet_coeffs_x.push_back(boost::make_shared<pdat::CellData<double> >(
                interior_box, 1, num_ghosts_wavelet_coeffs));
            wavelet_coeffs_y.push_back(boost::make_shared<pdat::CellData<double> >(
                interior_box, 1, num_ghosts_wavelet_coeffs));
            scaling_coeffs_x.push_back(boost::make_shared<pdat::CellData<double> >(
                interior_box, 1, num_ghosts_wavelet_coeffs));
            scaling_coeffs_y.push_back(boost::make_shared<pdat::CellData<double> >(
                interior_box, 1, num_ghosts_wavelet_coeffs));
        }
        
        // Get the pointers to the one-dimensional wavelet coefficients and initialize them with zero values.
        std::vector<double*> w_x;
        std::vector<double*> w_y;
        std::vector<double*> f_x;
        std::vector<double*> f_y;
        for (int li = 0; li < d_num_level; li++)
        {
            w_x.push_back(wavelet_coeffs_x[li]->getPointer(0));
            w_y.push_back(wavelet_coeffs_y[li]->getPointer(0));
            f_x.push_back(scaling_coeffs_x[li]->getPointer(0));
            f_y.push_back(scaling_coeffs_y[li]->getPointer(0));
            
            wavelet_coeffs_x[li]->fillAll(0.0);
            wavelet_coeffs_y[li]->fillAll(0.0);
            scaling_coeffs_x[li]->fillAll(0.0);
            scaling_coeffs_y[li]->fillAll(0.0);
        }
        
        /*
         * Compute scaling and wavelet coefficients at first level.
         */
        
        switch (d_k)
        {
            case 2:
            {
                /*
                 * Compute scaling and wavelet coefficients in the x-direction.
                 */
                
                // Get the pointer to the smoothed cell data in the x-direction.
                if (smooth_cell_data)
                {
                    f = smoothed_cell_data->getPointer(0);
                }
                else
                {
                    f = cell_data->getPointer(depth);
                }
                
                // Compute the starting and ending indices.
                const int start_index_x_i = -d_num_wavelet_ghosts[0] + 1;
                const int end_index_x_i = interior_dim_0 + d_num_wavelet_ghosts[0] - 1;
                const int start_index_x_j = -d_num_wavelet_ghosts[1];
                const int end_index_x_j = interior_dim_1 + d_num_wavelet_ghosts[1];
                
                for (int j = start_index_x_j; j < end_index_x_j; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = start_index_x_i; i < end_index_x_i; i++)
                    {
                        // Compute the linear indices.
                        const int idx = (i + num_ghosts_0_wavelet_coeffs) +
                            (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs;
                        
                        const int idx_cell_data = (i + num_ghosts_0_cell_data) +
                            (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        const int idx_cell_data_x_L = (i - 1 + num_ghosts_0_cell_data) +
                            (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        const int idx_cell_data_x_R = (i + 1 + num_ghosts_0_cell_data) +
                            (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        f_x[0][idx] = 0.5*(f[idx_cell_data_x_L] + f[idx_cell_data_x_R]);
                        
                        w_x[0][idx] = -0.5*(f[idx_cell_data_x_L] - 2.0*f[idx_cell_data] +
                            f[idx_cell_data_x_R]);
                    }
                }
                
                /*
                 * Compute scaling and wavelet coefficients in the y-direction.
                 */
                
                // Get the pointer to the smoothed cell data in the y-direction.
                if (smooth_cell_data)
                {
                    f = smoothed_cell_data->getPointer(1);
                }
                else
                {
                    f = cell_data->getPointer(depth);
                }
                
                // Compute the starting and ending indices.
                const int start_index_y_i = -d_num_wavelet_ghosts[0];
                const int end_index_y_i = interior_dim_0 + d_num_wavelet_ghosts[0];
                const int start_index_y_j = -d_num_wavelet_ghosts[1] + 1;
                const int end_index_y_j = interior_dim_1 + d_num_wavelet_ghosts[1] - 1;
                
                for (int j = start_index_y_j; j < end_index_y_j; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = start_index_y_i; i < end_index_y_i; i++)
                    {
                        // Compute the linear indices.
                        const int idx = (i + num_ghosts_0_wavelet_coeffs) +
                            (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs;
                        
                        const int idx_cell_data = (i + num_ghosts_0_cell_data) +
                            (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        const int idx_cell_data_y_B = (i + num_ghosts_0_cell_data) +
                            (j - 1 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        const int idx_cell_data_y_T = (i + num_ghosts_0_cell_data) +
                            (j + 1 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        f_y[0][idx] = 0.5*(f[idx_cell_data_y_B] + f[idx_cell_data_y_T]);
                        w_y[0][idx] = -0.5*(f[idx_cell_data_y_B] - 2.0*f[idx_cell_data] +
                            f[idx_cell_data_y_T]);
                    }
                }
                
                if (compute_variable_local_means)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx = (i + num_ghosts_0_wavelet_coeffs) +
                                (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs;
                            
                            const int idx_cell_data = (i + num_ghosts_0_cell_data) +
                                (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                            
                            const int idx_cell_data_x_L = (i - 1 + num_ghosts_0_cell_data) +
                                (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                            
                            const int idx_cell_data_x_R = (i + 1 + num_ghosts_0_cell_data) +
                                (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                            
                            const int idx_cell_data_y_B = (i + num_ghosts_0_cell_data) +
                                (j - 1 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                            
                            const int idx_cell_data_y_T = (i + num_ghosts_0_cell_data) +
                                (j + 1 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                            
                            double f_mean_x = f[idx_cell_data_x_L] + 2.0*f[idx_cell_data] +
                                              f[idx_cell_data_x_R];
                            
                            double f_mean_y = f[idx_cell_data_y_B] + 2.0*f[idx_cell_data] +
                                              f[idx_cell_data_y_T];
                            
                            f_mean[0][idx] = 0.5*sqrt(f_mean_x*f_mean_x + f_mean_y*f_mean_y);
                        }
                    }
                }
                
                break;
            }
            case 4:
            {
                /*
                 * Compute scaling and wavelet coefficients in the x-direction.
                 */
                
                // Get the pointer to the smoothed cell data in the x-direction.
                if (smooth_cell_data)
                {
                    f = smoothed_cell_data->getPointer(0);
                }
                else
                {
                    f = cell_data->getPointer(depth);
                }
                
                // Compute the starting and ending indices.
                const int start_index_x_i = -d_num_wavelet_ghosts[0] + 2;
                const int end_index_x_i = interior_dim_0 + d_num_wavelet_ghosts[0] - 2;
                const int start_index_x_j = -d_num_wavelet_ghosts[1];
                const int end_index_x_j = interior_dim_1 + d_num_wavelet_ghosts[1];
                
                for (int j = start_index_x_j; j < end_index_x_j; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = start_index_x_i; i < end_index_x_i; i++)
                    {
                        // Compute the linear indices.
                        const int idx = (i + num_ghosts_0_wavelet_coeffs) +
                            (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs;
                        
                        const int idx_cell_data = (i + num_ghosts_0_cell_data) +
                            (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        const int idx_cell_data_x_LL = (i - 2 + num_ghosts_0_cell_data) +
                            (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        const int idx_cell_data_x_L = (i - 1 + num_ghosts_0_cell_data) +
                            (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        const int idx_cell_data_x_R = (i + 1 + num_ghosts_0_cell_data) +
                            (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        const int idx_cell_data_x_RR = (i + 2 + num_ghosts_0_cell_data) +
                            (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        f_x[0][idx] = 1.0/6.0*(-f[idx_cell_data_x_LL] + 4.0*f[idx_cell_data_x_L] +
                            4.0*f[idx_cell_data_x_R] - f[idx_cell_data_x_RR]);
                        
                        w_x[0][idx] = 1.0/6.0*(f[idx_cell_data_x_LL] - 4.0*f[idx_cell_data_x_L] +
                            6.0*f[idx_cell_data] - 4.0*f[idx_cell_data_x_R] + f[idx_cell_data_x_RR]);
                    }
                }
                
                /*
                 * Compute scaling and wavelet coefficients in the y-direction.
                 */
                
                // Get the pointer to the smoothed cell data in the y-direction.
                if (smooth_cell_data)
                {
                    f = smoothed_cell_data->getPointer(1);
                }
                else
                {
                    f = cell_data->getPointer(depth);
                }
                
                // Compute the starting and ending indices.
                const int start_index_y_i = -d_num_wavelet_ghosts[0];
                const int end_index_y_i = interior_dim_0 + d_num_wavelet_ghosts[0];
                const int start_index_y_j = -d_num_wavelet_ghosts[1] + 2;
                const int end_index_y_j = interior_dim_1 + d_num_wavelet_ghosts[1] - 2;
                
                for (int j = start_index_y_j; j < end_index_y_j; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = start_index_y_i; i < end_index_y_i; i++)
                    {
                        // Compute indices.
                        const int idx = (i + num_ghosts_0_wavelet_coeffs) +
                            (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs;
                        
                        const int idx_cell_data = (i + num_ghosts_0_cell_data) +
                            (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        const int idx_cell_data_y_BB = (i + num_ghosts_0_cell_data) +
                            (j - 2 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        const int idx_cell_data_y_B = (i + num_ghosts_0_cell_data) +
                            (j - 1 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        const int idx_cell_data_y_T = (i + num_ghosts_0_cell_data) +
                            (j + 1 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        const int idx_cell_data_y_TT = (i + num_ghosts_0_cell_data) +
                            (j + 2 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        f_y[0][idx] = 1.0/6.0*(-f[idx_cell_data_y_BB] + 4.0*f[idx_cell_data_y_B] +
                            4.0*f[idx_cell_data_y_T] - f[idx_cell_data_y_TT]);
                        
                        w_y[0][idx] = 1.0/6.0*(f[idx_cell_data_y_BB] - 4.0*f[idx_cell_data_y_B] +
                            6.0*f[idx_cell_data] - 4.0*f[idx_cell_data_y_T] + f[idx_cell_data_y_TT]);
                    }
                }
                
                if (compute_variable_local_means)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx = (i + num_ghosts_0_wavelet_coeffs) +
                                (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs;
                            
                            const int idx_cell_data = (i + num_ghosts_0_cell_data) +
                                (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                            
                            const int idx_cell_data_x_LL = (i - 2 + num_ghosts_0_cell_data) +
                                (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                            
                            const int idx_cell_data_x_L = (i - 1 + num_ghosts_0_cell_data) +
                                (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                            
                            const int idx_cell_data_x_R = (i + 1 + num_ghosts_0_cell_data) +
                                (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                            
                            const int idx_cell_data_x_RR = (i + 2 + num_ghosts_0_cell_data) +
                                (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                            
                            const int idx_cell_data_y_BB = (i + num_ghosts_0_cell_data) +
                                (j - 2 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                            
                            const int idx_cell_data_y_B = (i + num_ghosts_0_cell_data) +
                                (j - 1 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                            
                            const int idx_cell_data_y_T = (i + num_ghosts_0_cell_data) +
                                (j + 1 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                            
                            const int idx_cell_data_y_TT = (i + num_ghosts_0_cell_data) +
                                (j + 2 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                            
                            double f_mean_x = f[idx_cell_data_x_LL] + 4.0*f[idx_cell_data_x_L] +
                                6.0*f[idx_cell_data] + 4.0*f[idx_cell_data_x_R] + f[idx_cell_data_x_RR];
                            
                            double f_mean_y = f[idx_cell_data_y_BB] + 4.0*f[idx_cell_data_y_B] +
                                6.0*f[idx_cell_data] + 4.0*f[idx_cell_data_y_T] + f[idx_cell_data_y_TT];
                            
                            f_mean[0][idx] = 1.0/6.0*sqrt(f_mean_x*f_mean_x + f_mean_y*f_mean_y);
                        }
                    }
                }
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "number of vanishing moments = "
                    << d_k
                    << " not supported."
                    << std::endl);
            }
        }
        
        /*
         * Compute scaling and wavelet coefficients at higher levels.
         */
        
        for (int li = 1; li < d_num_level; li++)
        {
            const int offset = pow(2, li);
            
            switch (d_k)
            {
                case 2:
                {
                    /*
                     * Compute scaling and wavelet coefficients in the x-direction.
                     */
                    
                    // Compute the starting and ending indices.
                    const int start_index_x_i = -d_num_wavelet_ghosts[0] + offset;
                    const int end_index_x_i = interior_dim_0 + d_num_wavelet_ghosts[0] - offset;
                    const int start_index_x_j = -d_num_wavelet_ghosts[1];
                    const int end_index_x_j = interior_dim_1 + d_num_wavelet_ghosts[1];
                    
                    for (int j = start_index_x_j; j < end_index_x_j; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = start_index_x_i; i < end_index_x_i; i++)
                        {
                            // Compute the linear indices.
                            const int idx = (i + num_ghosts_0_wavelet_coeffs) +
                                (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs;
                            
                            const int idx_x_L = (i - offset + num_ghosts_0_wavelet_coeffs) +
                                (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs;
                            
                            const int idx_x_R = (i + offset + num_ghosts_0_wavelet_coeffs) +
                                (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs;
                            
                            f_x[li][idx] = 0.5*(f_x[li-1][idx_x_L] + f_x[li-1][idx_x_R]);
                            
                            w_x[li][idx] = -0.5*(f_x[li-1][idx_x_L] - 2.0*f_x[li-1][idx] +
                                f_x[li-1][idx_x_R]);
                        }
                    }
                    
                    /*
                     * Compute scaling and wavelet coefficients in the y-direction.
                     */
                    
                    // Compute the starting and ending indices.
                    const int start_index_y_i = -d_num_wavelet_ghosts[0];
                    const int end_index_y_i = interior_dim_0 + d_num_wavelet_ghosts[0];
                    const int start_index_y_j = -d_num_wavelet_ghosts[1] + offset;
                    const int end_index_y_j = interior_dim_1 + d_num_wavelet_ghosts[1] - offset;
                    
                    for (int j = start_index_y_j; j < end_index_y_j; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = start_index_y_i; i < end_index_y_i; i++)
                        {
                            // Compute the linear indices.
                            const int idx = (i + num_ghosts_0_wavelet_coeffs) +
                                (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs;
                            
                            const int idx_y_B = (i + num_ghosts_0_wavelet_coeffs) +
                                (j - offset + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs;
                            
                            const int idx_y_T = (i + num_ghosts_0_wavelet_coeffs) +
                                (j + offset + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs;
                            
                            f_y[li][idx] = 0.5*(f_y[li-1][idx_y_B] + f_y[li-1][idx_y_T]);
                            
                            w_y[li][idx] = -0.5*(f_y[li-1][idx_y_B] - 2.0*f_y[li-1][idx] +
                                f_y[li-1][idx_y_T]);
                        }
                    }
                    
                    if (compute_variable_local_means)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
#ifdef HAMERS_ENABLE_SIMD
                            #pragma omp simd
#endif
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                // Compute the linear indices.
                                const int idx = (i + num_ghosts_0_wavelet_coeffs) +
                                    (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs;
                                
                                const int idx_x_L = (i - offset + num_ghosts_0_wavelet_coeffs) +
                                    (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs;
                                
                                const int idx_x_R = (i + offset + num_ghosts_0_wavelet_coeffs) +
                                    (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs;
                                
                                const int idx_y_B = (i + num_ghosts_0_wavelet_coeffs) +
                                    (j - offset + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs;
                                
                                const int idx_y_T = (i + num_ghosts_0_wavelet_coeffs) +
                                    (j + offset + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs;
                                
                                double f_mean_x = f_x[li-1][idx_x_L] + 2.0*f_x[li-1][idx] +
                                    f_x[li-1][idx_x_R];
                                
                                double f_mean_y = f_y[li-1][idx_y_B] + 2.0*f_y[li-1][idx] +
                                    f_y[li-1][idx_y_T];
                                
                                f_mean[li][idx] = 0.5*sqrt(f_mean_x*f_mean_x + f_mean_y*f_mean_y);
                            }
                        }
                    }
                    
                    break;
                }
                case 4:
                {
                    /*
                     * Compute scaling and wavelet coefficients in the x-direction.
                     */
                    
                    // Compute the starting and ending indices.
                    const int start_index_x_i = -d_num_wavelet_ghosts[0] + 2*offset;
                    const int end_index_x_i = interior_dim_0 + d_num_wavelet_ghosts[0] - 2*offset;
                    const int start_index_x_j = -d_num_wavelet_ghosts[1];
                    const int end_index_x_j = interior_dim_1 + d_num_wavelet_ghosts[1];
                    
                    for (int j = start_index_x_j; j < end_index_x_j; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = start_index_x_i; i < end_index_x_i; i++)
                        {
                            // Compute the linear indices.
                            const int idx = (i + num_ghosts_0_wavelet_coeffs) +
                                (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs;
                            
                            const int idx_x_LL = (i - 2*offset + num_ghosts_0_wavelet_coeffs) +
                                (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs;
                            
                            const int idx_x_L = (i - offset + num_ghosts_0_wavelet_coeffs) +
                                (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs;
                            
                            const int idx_x_R = (i + offset + num_ghosts_0_wavelet_coeffs) +
                                (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs;
                            
                            const int idx_x_RR = (i + 2*offset + num_ghosts_0_wavelet_coeffs) +
                                (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs;
                            
                            f_x[li][idx] = 1.0/6.0*(-f_x[li-1][idx_x_LL] + 4.0*f_x[li-1][idx_x_L] +
                                4.0*f_x[li-1][idx_x_R] - f_x[li-1][idx_x_RR]);
                            
                            w_x[li][idx] = 1.0/6.0*(f_x[li-1][idx_x_LL] - 4.0*f_x[li-1][idx_x_L] +
                                6.0*f_x[li-1][idx] - 4.0*f_x[li-1][idx_x_R] + f_x[li-1][idx_x_RR]);
                        }
                    }
                    
                    /*
                     * Compute scaling and wavelet coefficients in the y-direction.
                     */
                    
                    // Compute the starting and ending indices.
                    const int start_index_y_i = -d_num_wavelet_ghosts[0];
                    const int end_index_y_i = interior_dim_0 + d_num_wavelet_ghosts[0];
                    const int start_index_y_j = -d_num_wavelet_ghosts[1] + 2*offset;
                    const int end_index_y_j = interior_dim_1 + d_num_wavelet_ghosts[1] - 2*offset;
                    
                    for (int j = start_index_y_j; j < end_index_y_j; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = start_index_y_i; i < end_index_y_i; i++)
                        {
                            // Compute the linear indices.
                            const int idx = (i + num_ghosts_0_wavelet_coeffs) +
                                (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs;
                            
                            const int idx_y_BB = (i + num_ghosts_0_wavelet_coeffs) +
                                (j - 2*offset + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs;
                            
                            const int idx_y_B = (i + num_ghosts_0_wavelet_coeffs) +
                                (j - offset + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs;
                            
                            const int idx_y_T = (i + num_ghosts_0_wavelet_coeffs) +
                                (j + offset + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs;
                            
                            const int idx_y_TT = (i + num_ghosts_0_wavelet_coeffs) +
                                (j + 2*offset + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs;
                            
                            f_y[li][idx] = 1.0/6.0*(-f_y[li-1][idx_y_BB] + 4.0*f_y[li-1][idx_y_B] +
                                4.0*f_y[li-1][idx_y_T] - f_y[li-1][idx_y_TT]);
                            
                            w_y[li][idx] = 1.0/6.0*(f_y[li-1][idx_y_BB] - 4.0*f_y[li-1][idx_y_B] +
                                6.0*f_y[li-1][idx] - 4.0*f_y[li-1][idx_y_T] + f_y[li-1][idx_y_TT]);
                        }
                    }
                    
                    if (compute_variable_local_means)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
#ifdef HAMERS_ENABLE_SIMD
                            #pragma omp simd
#endif
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                // Compute the linear indices.
                                const int idx = (i + num_ghosts_0_wavelet_coeffs) +
                                    (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs;
                                
                                const int idx_x_LL = (i - 2*offset + num_ghosts_0_wavelet_coeffs) +
                                    (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs;
                                
                                const int idx_x_L = (i - offset + num_ghosts_0_wavelet_coeffs) +
                                    (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs;
                                
                                const int idx_x_R = (i + offset + num_ghosts_0_wavelet_coeffs) +
                                    (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs;
                                
                                const int idx_x_RR = (i + 2*offset + num_ghosts_0_wavelet_coeffs) +
                                    (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs;
                                
                                const int idx_y_BB = (i + num_ghosts_0_wavelet_coeffs) +
                                    (j - 2*offset + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs;
                                
                                const int idx_y_B = (i + num_ghosts_0_wavelet_coeffs) +
                                    (j - offset + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs;
                                
                                const int idx_y_T = (i + num_ghosts_0_wavelet_coeffs) +
                                    (j + offset + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs;
                                
                                const int idx_y_TT = (i + num_ghosts_0_wavelet_coeffs) +
                                    (j + 2*offset + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs;
                                
                                double f_mean_x = f_x[li-1][idx_x_LL] + 4.0*f_x[li-1][idx_x_L] +
                                    6.0*f_x[li-1][idx] + 4.0*f_x[li-1][idx_x_R] + f_x[li-1][idx_x_RR];
                                
                                double f_mean_y = f_y[li-1][idx_y_BB] + 4.0*f_y[li-1][idx_y_B] +
                                    6.0*f_y[li-1][idx] + 4.0*f_y[li-1][idx_y_T] + f_y[li-1][idx_y_TT];
                                
                                f_mean[li][idx] = 1.0/6.0*sqrt(f_mean_x*f_mean_x + f_mean_y*f_mean_y);
                            }
                        }
                    }
                    
                    break;
                }
                default:
                {
                    TBOX_ERROR(d_object_name
                        << ": "
                        << "number of vanishing moments = "
                        << d_k
                        << " not supported."
                        << std::endl);
                }
            }
        }
        
        /*
         * Compute the two-dimensional wavelet coefficients from the one-dimensional wavelet coefficients.
         */
        
        for (int li = 0; li < d_num_level; li++)
        {
            const int offset = pow(2, li + 1);
            
            // Compute the starting and ending indices.
            const int start_index_i = -d_p*offset;
            const int end_index_i = interior_dim_0 + d_q*offset;
            const int start_index_j = -d_p*offset;
            const int end_index_j = interior_dim_1 + d_q*offset;
            
            for (int j = start_index_j; j < end_index_j; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = start_index_i; i < end_index_i; i++)
                {
                    const int idx = (i + num_ghosts_0_wavelet_coeffs) +
                        (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs;
                    
                    w[li][idx] = sqrt(w_x[li][idx]*w_x[li][idx] + w_y[li][idx]*w_y[li][idx]);
                }
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        const int interior_dim_2 = interior_dims[2];
        
        const int num_ghosts_0_cell_data = num_ghosts_cell_data[0];
        const int num_ghosts_1_cell_data = num_ghosts_cell_data[1];
        const int num_ghosts_2_cell_data = num_ghosts_cell_data[2];
        const int ghostcell_dim_0_cell_data = ghostcell_dims_cell_data[0];
        const int ghostcell_dim_1_cell_data = ghostcell_dims_cell_data[1];
        
        const int num_ghosts_0_wavelet_coeffs = num_ghosts_wavelet_coeffs[0];
        const int num_ghosts_1_wavelet_coeffs = num_ghosts_wavelet_coeffs[1];
        const int num_ghosts_2_wavelet_coeffs = num_ghosts_wavelet_coeffs[2];
        const int ghostcell_dim_0_wavelet_coeffs = ghostcell_dims_wavelet_coeffs[0];
        const int ghostcell_dim_1_wavelet_coeffs = ghostcell_dims_wavelet_coeffs[1];
        
        // Allocate wavelet coefficients in different dimensions.
        std::vector<boost::shared_ptr<pdat::CellData<double> > > wavelet_coeffs_x;
        std::vector<boost::shared_ptr<pdat::CellData<double> > > wavelet_coeffs_y;
        std::vector<boost::shared_ptr<pdat::CellData<double> > > wavelet_coeffs_z;
        std::vector<boost::shared_ptr<pdat::CellData<double> > > scaling_coeffs_x;
        std::vector<boost::shared_ptr<pdat::CellData<double> > > scaling_coeffs_y;
        std::vector<boost::shared_ptr<pdat::CellData<double> > > scaling_coeffs_z;
        for (int li = 0; li < d_num_level; li++)
        {
            wavelet_coeffs_x.push_back(boost::make_shared<pdat::CellData<double> >(
                interior_box, 1, num_ghosts_wavelet_coeffs));
            wavelet_coeffs_y.push_back(boost::make_shared<pdat::CellData<double> >(
                interior_box, 1, num_ghosts_wavelet_coeffs));
            wavelet_coeffs_z.push_back(boost::make_shared<pdat::CellData<double> >(
                interior_box, 1, num_ghosts_wavelet_coeffs));
            scaling_coeffs_x.push_back(boost::make_shared<pdat::CellData<double> >(
                interior_box, 1, num_ghosts_wavelet_coeffs));
            scaling_coeffs_y.push_back(boost::make_shared<pdat::CellData<double> >(
                interior_box, 1, num_ghosts_wavelet_coeffs));
            scaling_coeffs_z.push_back(boost::make_shared<pdat::CellData<double> >(
                interior_box, 1, num_ghosts_wavelet_coeffs));
        }
        
        // Get the pointers to the one-dimensional wavelet coefficients and initialize them with zero values.
        std::vector<double*> w_x;
        std::vector<double*> w_y;
        std::vector<double*> w_z;
        std::vector<double*> f_x;
        std::vector<double*> f_y;
        std::vector<double*> f_z;
        for (int li = 0; li < d_num_level; li++)
        {
            w_x.push_back(wavelet_coeffs_x[li]->getPointer(0));
            w_y.push_back(wavelet_coeffs_y[li]->getPointer(0));
            w_z.push_back(wavelet_coeffs_z[li]->getPointer(0));
            f_x.push_back(scaling_coeffs_x[li]->getPointer(0));
            f_y.push_back(scaling_coeffs_y[li]->getPointer(0));
            f_z.push_back(scaling_coeffs_z[li]->getPointer(0));
            
            wavelet_coeffs_x[li]->fillAll(0.0);
            wavelet_coeffs_y[li]->fillAll(0.0);
            wavelet_coeffs_z[li]->fillAll(0.0);
            scaling_coeffs_x[li]->fillAll(0.0);
            scaling_coeffs_y[li]->fillAll(0.0);
            scaling_coeffs_z[li]->fillAll(0.0);
        }
        
        /*
         * Compute scaling and wavelet coefficients at first level.
         */
        switch (d_k)
        {
            case 2:
            {
                /*
                 * Compute scaling and wavelet coefficients in the x-direction.
                 */
                
                // Get the pointer to the smoothed cell data in the x-direction.
                if (smooth_cell_data)
                {
                    f = smoothed_cell_data->getPointer(0);
                }
                else
                {
                    f = cell_data->getPointer(depth);
                }
                
                // Compute the starting and ending indices.
                const int start_index_x_i = -d_num_wavelet_ghosts[0] + 1;
                const int end_index_x_i = interior_dim_0 + d_num_wavelet_ghosts[0] - 1;
                const int start_index_x_j = -d_num_wavelet_ghosts[1];
                const int end_index_x_j = interior_dim_1 + d_num_wavelet_ghosts[1];
                const int start_index_x_k = -d_num_wavelet_ghosts[2];
                const int end_index_x_k = interior_dim_2 + d_num_wavelet_ghosts[2];
                
                for (int k = start_index_x_k; k < end_index_x_k; k++)
                {
                    for (int j = start_index_x_j; j < end_index_x_j; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = start_index_x_i; i < end_index_x_i; i++)
                        {
                            // Compute the linear indices.
                            const int idx = (i + num_ghosts_0_wavelet_coeffs) +
                                (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                (k + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                    ghostcell_dim_1_wavelet_coeffs;
                            
                            const int idx_cell_data = (i + num_ghosts_0_cell_data) +
                                (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            const int idx_cell_data_x_L = (i - 1 + num_ghosts_0_cell_data) +
                                (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            const int idx_cell_data_x_R = (i + 1 + num_ghosts_0_cell_data) +
                                (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            f_x[0][idx] = 0.5*(f[idx_cell_data_x_L] + f[idx_cell_data_x_R]);
                            
                            w_x[0][idx] = -0.5*(f[idx_cell_data_x_L] - 2.0*f[idx_cell_data] +
                                f[idx_cell_data_x_R]);
                        }
                    }
                }
                
                /*
                 * Compute scaling and wavelet coefficients in the y-direction.
                 */
                
                // Get the pointer to the smoothed cell data in the y-direction.
                if (smooth_cell_data)
                {
                    f = smoothed_cell_data->getPointer(1);
                }
                else
                {
                    f = cell_data->getPointer(depth);
                }
                
                // Compute the starting and ending indices.
                const int start_index_y_i = -d_num_wavelet_ghosts[0];
                const int end_index_y_i = interior_dim_0 + d_num_wavelet_ghosts[0];
                const int start_index_y_j = -d_num_wavelet_ghosts[1] + 1;
                const int end_index_y_j = interior_dim_1 + d_num_wavelet_ghosts[1] - 1;
                const int start_index_y_k = -d_num_wavelet_ghosts[2];
                const int end_index_y_k = interior_dim_2 + d_num_wavelet_ghosts[2];
                
                for (int k = start_index_y_k; k < end_index_y_k; k++)
                {
                    for (int j = start_index_y_j; j < end_index_y_j; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = start_index_y_i; i < end_index_y_i; i++)
                        {
                            // Compute the linear indices.
                            const int idx = (i + num_ghosts_0_wavelet_coeffs) +
                                (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                (k + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                    ghostcell_dim_1_wavelet_coeffs;
                            
                            const int idx_cell_data = (i + num_ghosts_0_cell_data) +
                                (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            const int idx_cell_data_y_B = (i + num_ghosts_0_cell_data) +
                                (j - 1 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            const int idx_cell_data_y_T = (i + num_ghosts_0_cell_data) +
                                (j + 1 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            f_y[0][idx] = 0.5*(f[idx_cell_data_y_B] + f[idx_cell_data_y_T]);
                            
                            w_y[0][idx] = -0.5*(f[idx_cell_data_y_B] - 2.0*f[idx_cell_data] +
                                f[idx_cell_data_y_T]);
                        }
                    }
                }
                
                /*
                 * Compute scaling and wavelet coefficients in the z-direction.
                 */
                
                // Get the pointer to the smoothed cell data in the z-direction.
                if (smooth_cell_data)
                {
                    f = smoothed_cell_data->getPointer(2);
                }
                else
                {
                    f = cell_data->getPointer(depth);
                }
                
                // Compute the starting and ending indices.
                const int start_index_z_i = -d_num_wavelet_ghosts[0];
                const int end_index_z_i = interior_dim_0 + d_num_wavelet_ghosts[0];
                const int start_index_z_j = -d_num_wavelet_ghosts[1];
                const int end_index_z_j = interior_dim_1 + d_num_wavelet_ghosts[1];
                const int start_index_z_k = -d_num_wavelet_ghosts[2] + 1;
                const int end_index_z_k = interior_dim_2 + d_num_wavelet_ghosts[1] - 1;
                
                for (int k = start_index_z_k; k < end_index_z_k; k++)
                {
                    for (int j = start_index_z_j; j < end_index_z_j; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = start_index_z_i; i < end_index_z_i; i++)
                        {
                            // Compute the linear indices.
                            const int idx = (i + num_ghosts_0_wavelet_coeffs) +
                                (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                (k + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                    ghostcell_dim_1_wavelet_coeffs;
                            
                            const int idx_cell_data = (i + num_ghosts_0_cell_data) +
                                (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            const int idx_cell_data_z_B = (i + num_ghosts_0_cell_data) +
                                (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k - 1 + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            const int idx_cell_data_z_F = (i + num_ghosts_0_cell_data) +
                                (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + 1 + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            f_z[0][idx] = 0.5*(f[idx_cell_data_z_B] + f[idx_cell_data_z_F]);
                            
                            w_z[0][idx] = -0.5*(f[idx_cell_data_z_B] - 2.0*f[idx_cell_data] +
                                f[idx_cell_data_z_F]);
                        }
                    }
                }
                
                if (compute_variable_local_means)
                {
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
#ifdef HAMERS_ENABLE_SIMD
                            #pragma omp simd
#endif
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                // Compute the indices.
                                const int idx = (i + num_ghosts_0_wavelet_coeffs) +
                                    (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                    (k + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                        ghostcell_dim_1_wavelet_coeffs;
                                
                                const int idx_cell_data = (i + num_ghosts_0_cell_data) +
                                    (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                    (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                        ghostcell_dim_1_cell_data;
                                
                                const int idx_cell_data_x_L = (i - 1 + num_ghosts_0_cell_data) +
                                    (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                    (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                        ghostcell_dim_1_cell_data;
                                
                                const int idx_cell_data_x_R = (i + 1 + num_ghosts_0_cell_data) +
                                    (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                    (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                        ghostcell_dim_1_cell_data;
                                
                                const int idx_cell_data_y_B = (i + num_ghosts_0_cell_data) +
                                    (j - 1 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                    (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                        ghostcell_dim_1_cell_data;
                                
                                const int idx_cell_data_y_T = (i + num_ghosts_0_cell_data) +
                                    (j + 1 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                    (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                        ghostcell_dim_1_cell_data;
                                
                                const int idx_cell_data_z_B = (i + num_ghosts_0_cell_data) +
                                    (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                    (k - 1 + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                        ghostcell_dim_1_cell_data;
                                
                                const int idx_cell_data_z_F = (i + num_ghosts_0_cell_data) +
                                    (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                    (k + 1 + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                        ghostcell_dim_1_cell_data;
                                
                                double f_mean_x = f[idx_cell_data_x_L] + 2.0*f[idx_cell_data] +
                                    f[idx_cell_data_x_R];
                                
                                double f_mean_y = f[idx_cell_data_y_B] + 2.0*f[idx_cell_data] +
                                    f[idx_cell_data_y_T];
                                
                                double f_mean_z = f[idx_cell_data_z_B] + 2.0*f[idx_cell_data] +
                                    f[idx_cell_data_z_F];
                                
                                f_mean[0][idx] = 0.5*sqrt(f_mean_x*f_mean_x + f_mean_y*f_mean_y +
                                    f_mean_z*f_mean_z);
                            }
                        }
                    }
                }
                
                break;
            }
            case 4:
            {
                /*
                 * Compute scaling and wavelet coefficients in the x-direction.
                 */
                
                // Get the pointer to the smoothed cell data in the x-direction.
                if (smooth_cell_data)
                {
                    f = smoothed_cell_data->getPointer(0);
                }
                else
                {
                    f = cell_data->getPointer(depth);
                }
                
                // Compute the starting and ending indices.
                const int start_index_x_i = -d_num_wavelet_ghosts[0] + 2;
                const int end_index_x_i = interior_dim_0 + d_num_wavelet_ghosts[0] - 2;
                const int start_index_x_j = -d_num_wavelet_ghosts[1];
                const int end_index_x_j = interior_dim_1 + d_num_wavelet_ghosts[1];
                const int start_index_x_k = -d_num_wavelet_ghosts[2];
                const int end_index_x_k = interior_dim_2 + d_num_wavelet_ghosts[2];
                
                for (int k = start_index_x_k; k < end_index_x_k; k++)
                {
                    for (int j = start_index_x_j; j < end_index_x_j; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = start_index_x_i; i < end_index_x_i; i++)
                        {
                            // Compute the indices.
                            const int idx = (i + num_ghosts_0_wavelet_coeffs) +
                                (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                (k + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                    ghostcell_dim_1_wavelet_coeffs;
                            
                            const int idx_cell_data = (i + num_ghosts_0_cell_data) +
                                (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            const int idx_cell_data_x_LL = (i - 2 + num_ghosts_0_cell_data) +
                                (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            const int idx_cell_data_x_L = (i - 1 + num_ghosts_0_cell_data) +
                                (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            const int idx_cell_data_x_R = (i + 1 + num_ghosts_0_cell_data) +
                                (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            const int idx_cell_data_x_RR = (i + 2 + num_ghosts_0_cell_data) +
                                (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            f_x[0][idx] = 1.0/6.0*(-f[idx_cell_data_x_LL] + 4.0*f[idx_cell_data_x_L] +
                                4.0*f[idx_cell_data_x_R] - f[idx_cell_data_x_RR]);
                            
                            w_x[0][idx] = 1.0/6.0*(f[idx_cell_data_x_LL] - 4.0*f[idx_cell_data_x_L] +
                                6.0*f[idx_cell_data] - 4.0*f[idx_cell_data_x_R] + f[idx_cell_data_x_RR]);
                        }
                    }
                }
                
                /*
                 * Compute scaling and wavelet coefficients in the y-direction.
                 */
                
                // Get the pointer to the smoothed cell data in the y-direction.
                if (smooth_cell_data)
                {
                    f = smoothed_cell_data->getPointer(1);
                }
                else
                {
                    f = cell_data->getPointer(depth);
                }
                
                // Compute the starting and ending indices.
                const int start_index_y_i = -d_num_wavelet_ghosts[0];
                const int end_index_y_i = interior_dim_0 + d_num_wavelet_ghosts[0];
                const int start_index_y_j = -d_num_wavelet_ghosts[1] + 2;
                const int end_index_y_j = interior_dim_1 + d_num_wavelet_ghosts[1] - 2;
                const int start_index_y_k = -d_num_wavelet_ghosts[2];
                const int end_index_y_k = interior_dim_2 + d_num_wavelet_ghosts[2];
                
                for (int k = start_index_y_k; k < end_index_y_k; k++)
                {
                    for (int j = start_index_y_j; j < end_index_y_j; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = start_index_y_i; i < end_index_y_i; i++)
                        {
                            // Compute the linear indices.
                            const int idx = (i + num_ghosts_0_wavelet_coeffs) +
                                (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                (k + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                    ghostcell_dim_1_wavelet_coeffs;
                            
                            const int idx_cell_data = (i + num_ghosts_0_cell_data) +
                                (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            const int idx_cell_data_y_BB = (i + num_ghosts_0_cell_data) +
                                (j - 2 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            const int idx_cell_data_y_B = (i + num_ghosts_0_cell_data) +
                                (j - 1 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            const int idx_cell_data_y_T = (i + num_ghosts_0_cell_data) +
                                (j + 1 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            const int idx_cell_data_y_TT = (i + num_ghosts_0_cell_data) +
                                (j + 2 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            f_y[0][idx] = 1.0/6.0*(-f[idx_cell_data_y_BB] + 4.0*f[idx_cell_data_y_B] +
                                4.0*f[idx_cell_data_y_T] - f[idx_cell_data_y_TT]);
                            
                            w_y[0][idx] = 1.0/6.0*(f[idx_cell_data_y_BB] - 4.0*f[idx_cell_data_y_B] +
                                6.0*f[idx_cell_data] - 4.0*f[idx_cell_data_y_T] + f[idx_cell_data_y_TT]);
                        }
                    }
                }
                
                /*
                 * Compute scaling and wavelet coefficients in the z-direction.
                 */
                
                // Get the pointer to the smoothed cell data in the z-direction.
                if (smooth_cell_data)
                {
                    f = smoothed_cell_data->getPointer(2);
                }
                else
                {
                    f = cell_data->getPointer(depth);
                }
                
                // Compute the starting and ending indices.
                const int start_index_z_i = -d_num_wavelet_ghosts[0];
                const int end_index_z_i = interior_dim_0 + d_num_wavelet_ghosts[0];
                const int start_index_z_j = -d_num_wavelet_ghosts[1];
                const int end_index_z_j = interior_dim_1 + d_num_wavelet_ghosts[1];
                const int start_index_z_k = -d_num_wavelet_ghosts[2] + 2;
                const int end_index_z_k = interior_dim_2 + d_num_wavelet_ghosts[1] - 2;
                
                for (int k = start_index_z_k; k < end_index_z_k; k++)
                {
                    for (int j = start_index_z_j; j < end_index_z_j; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = start_index_z_i; i < end_index_z_i; i++)
                        {
                            // Compute the linear indices.
                            const int idx = (i + num_ghosts_0_wavelet_coeffs) +
                                (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                (k + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                    ghostcell_dim_1_wavelet_coeffs;
                            
                            const int idx_cell_data = (i + num_ghosts_0_cell_data) +
                                (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            const int idx_cell_data_z_BB = (i + num_ghosts_0_cell_data) +
                                (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k - 2 + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            const int idx_cell_data_z_B = (i + num_ghosts_0_cell_data) +
                                (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k - 1 + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            const int idx_cell_data_z_F = (i + num_ghosts_0_cell_data) +
                                (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + 1 + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            const int idx_cell_data_z_FF = (i + num_ghosts_0_cell_data) +
                                (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + 2 + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            f_z[0][idx] = 1.0/6.0*(-f[idx_cell_data_z_BB] + 4.0*f[idx_cell_data_z_B] +
                                4.0*f[idx_cell_data_z_F] - f[idx_cell_data_z_FF]);
                            
                            w_z[0][idx] = 1.0/6.0*(f[idx_cell_data_z_BB] - 4.0*f[idx_cell_data_z_B] +
                                6.0*f[idx_cell_data] - 4.0*f[idx_cell_data_z_F] + f[idx_cell_data_z_FF]);
                        }
                    }
                }
                
                if (compute_variable_local_means)
                {
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
#ifdef HAMERS_ENABLE_SIMD
                            #pragma omp simd
#endif
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                // Compute the linear indices.
                                const int idx = (i + num_ghosts_0_wavelet_coeffs) +
                                    (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                    (k + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                        ghostcell_dim_1_wavelet_coeffs;
                                
                                const int idx_cell_data = (i + num_ghosts_0_cell_data) +
                                    (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                    (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                        ghostcell_dim_1_cell_data;
                                
                                const int idx_cell_data_x_LL = (i - 2 + num_ghosts_0_cell_data) +
                                    (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                    (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                        ghostcell_dim_1_cell_data;
                                
                                const int idx_cell_data_x_L = (i - 1 + num_ghosts_0_cell_data) +
                                    (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                    (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                        ghostcell_dim_1_cell_data;
                                
                                const int idx_cell_data_x_R = (i + 1 + num_ghosts_0_cell_data) +
                                    (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                    (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                        ghostcell_dim_1_cell_data;
                                
                                const int idx_cell_data_x_RR = (i + 2 + num_ghosts_0_cell_data) +
                                    (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                    (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                        ghostcell_dim_1_cell_data;
                                
                                const int idx_cell_data_y_BB = (i + num_ghosts_0_cell_data) +
                                    (j - 2 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                    (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                        ghostcell_dim_1_cell_data;
                                
                                const int idx_cell_data_y_B = (i + num_ghosts_0_cell_data) +
                                    (j - 1 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                    (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                        ghostcell_dim_1_cell_data;
                                
                                const int idx_cell_data_y_T = (i + num_ghosts_0_cell_data) +
                                    (j + 1 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                    (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                        ghostcell_dim_1_cell_data;
                                
                                const int idx_cell_data_y_TT = (i + num_ghosts_0_cell_data) +
                                    (j + 2 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                    (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                        ghostcell_dim_1_cell_data;
                                
                                const int idx_cell_data_z_BB = (i + num_ghosts_0_cell_data) +
                                    (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                    (k - 2 + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                        ghostcell_dim_1_cell_data;
                                
                                const int idx_cell_data_z_B = (i + num_ghosts_0_cell_data) +
                                    (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                    (k - 1 + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                        ghostcell_dim_1_cell_data;
                                
                                const int idx_cell_data_z_F = (i + num_ghosts_0_cell_data) +
                                    (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                    (k + 1 + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                        ghostcell_dim_1_cell_data;
                                
                                const int idx_cell_data_z_FF = (i + num_ghosts_0_cell_data) +
                                    (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                    (k + 2 + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                        ghostcell_dim_1_cell_data;
                                
                                double f_mean_x = f[idx_cell_data_x_LL] + 4.0*f[idx_cell_data_x_L] +
                                    6.0*f[idx_cell_data] + 4.0*f[idx_cell_data_x_R] + f[idx_cell_data_x_RR];
                                
                                double f_mean_y = f[idx_cell_data_y_BB] + 4.0*f[idx_cell_data_y_B] +
                                    6.0*f[idx_cell_data] + 4.0*f[idx_cell_data_y_T] + f[idx_cell_data_y_TT];
                                
                                double f_mean_z = f[idx_cell_data_z_BB] + 4.0*f[idx_cell_data_z_B] +
                                    6.0*f[idx_cell_data] + 4.0*f[idx_cell_data_z_F] + f[idx_cell_data_z_FF];
                                
                                f_mean[0][idx] = 1.0/6.0*sqrt(f_mean_x*f_mean_x + f_mean_y*f_mean_y +
                                    f_mean_z*f_mean_z);
                            }
                        }
                    }
                }
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "number of vanishing moments = "
                    << d_k
                    << " not supported."
                    << std::endl);
            }
        }
        
        /*
         * Compute scaling and wavelet coefficients at higher levels.
         */
        
        for (int li = 1; li < d_num_level; li++)
        {
            const int offset = pow(2, li);
            
            switch (d_k)
            {
                case 2:
                {
                    /*
                     * Compute scaling and wavelet coefficients in the x-direction.
                     */
                    
                    // Compute the starting and ending indices.
                    const int start_index_x_i = -d_num_wavelet_ghosts[0] + offset;
                    const int end_index_x_i = interior_dim_0 + d_num_wavelet_ghosts[0] - offset;
                    const int start_index_x_j = -d_num_wavelet_ghosts[1];
                    const int end_index_x_j = interior_dim_1 + d_num_wavelet_ghosts[1];
                    const int start_index_x_k = -d_num_wavelet_ghosts[2];
                    const int end_index_x_k = interior_dim_2 + d_num_wavelet_ghosts[2];
                    
                    for (int k = start_index_x_k; k < end_index_x_k; k++)
                    {
                        for (int j = start_index_x_j; j < end_index_x_j; j++)
                        {
#ifdef HAMERS_ENABLE_SIMD
                            #pragma omp simd
#endif
                            for (int i = start_index_x_i; i < end_index_x_i; i++)
                            {
                                // Compute the linear indices.
                                const int idx = (i + num_ghosts_0_wavelet_coeffs) +
                                    (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                    (k + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                        ghostcell_dim_1_wavelet_coeffs;
                                
                                const int idx_x_L = (i - offset + num_ghosts_0_wavelet_coeffs) +
                                    (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                    (k + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                        ghostcell_dim_1_wavelet_coeffs;
                                
                                const int idx_x_R = (i + offset + num_ghosts_0_wavelet_coeffs) +
                                    (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                    (k + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                        ghostcell_dim_1_wavelet_coeffs;
                                
                                f_x[li][idx] = 0.5*(f_x[li-1][idx_x_L] + f_x[li-1][idx_x_R]);
                                
                                w_x[li][idx] = -0.5*(f_x[li-1][idx_x_L] - 2.0*f_x[li-1][idx] +
                                    f_x[li-1][idx_x_R]);
                            }
                        }
                    }
                    
                    /*
                     * Compute scaling and wavelet coefficients in the y-direction.
                     */
                    
                    // Compute the starting and ending indices.
                    const int start_index_y_i = -d_num_wavelet_ghosts[0];
                    const int end_index_y_i = interior_dim_0 + d_num_wavelet_ghosts[0];
                    const int start_index_y_j = -d_num_wavelet_ghosts[1] + offset ;
                    const int end_index_y_j = interior_dim_1 + d_num_wavelet_ghosts[1] - offset;
                    const int start_index_y_k = -d_num_wavelet_ghosts[2];
                    const int end_index_y_k = interior_dim_2 + d_num_wavelet_ghosts[2];
                    
                    for (int k = start_index_y_k; k < end_index_y_k; k++)
                    {
                        for (int j = start_index_y_j; j < end_index_y_j; j++)
                        {
#ifdef HAMERS_ENABLE_SIMD
                            #pragma omp simd
#endif
                            for (int i = start_index_y_i; i < end_index_y_i; i++)
                            {
                                // Compute the linear indices.
                                const int idx = (i + num_ghosts_0_wavelet_coeffs) +
                                    (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                    (k + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                        ghostcell_dim_1_wavelet_coeffs;
                                
                                const int idx_y_B = (i + num_ghosts_0_wavelet_coeffs) +
                                    (j - offset + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                    (k + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                        ghostcell_dim_1_wavelet_coeffs;
                                
                                const int idx_y_T = (i + num_ghosts_0_wavelet_coeffs) +
                                    (j + offset + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                    (k + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                        ghostcell_dim_1_wavelet_coeffs;
                                
                                f_y[li][idx] = 0.5*(f_y[li-1][idx_y_B] + f_y[li-1][idx_y_T]);
                                
                                w_y[li][idx] = -0.5*(f_y[li-1][idx_y_B] - 2.0*f_y[li-1][idx] +
                                    f_y[li-1][idx_y_T]);
                            }
                        }
                    }
                    
                    /*
                     * Compute scaling and wavelet coefficients in the z-direction.
                     */
                    
                    // Compute the starting and ending indices.
                    const int start_index_z_i = -d_num_wavelet_ghosts[0];
                    const int end_index_z_i = interior_dim_0 + d_num_wavelet_ghosts[0];
                    const int start_index_z_j = -d_num_wavelet_ghosts[1];
                    const int end_index_z_j = interior_dim_1 + d_num_wavelet_ghosts[1];
                    const int start_index_z_k = -d_num_wavelet_ghosts[2] + offset ;
                    const int end_index_z_k = interior_dim_2 + d_num_wavelet_ghosts[2] - offset;
                    
                    for (int k = start_index_z_k; k < end_index_z_k; k++)
                    {
                        for (int j = start_index_z_j; j < end_index_z_j; j++)
                        {
#ifdef HAMERS_ENABLE_SIMD
                            #pragma omp simd
#endif
                            for (int i = start_index_z_i; i < end_index_z_i; i++)
                            {
                                // Compute the linear indices.
                                const int idx = (i + num_ghosts_0_wavelet_coeffs) +
                                    (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                    (k + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                        ghostcell_dim_1_wavelet_coeffs;
                                
                                const int idx_z_B = (i + num_ghosts_0_wavelet_coeffs) +
                                    (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                    (k - offset + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                        ghostcell_dim_1_wavelet_coeffs;
                                
                                const int idx_z_F = (i + num_ghosts_0_wavelet_coeffs) +
                                    (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                    (k + offset + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                        ghostcell_dim_1_wavelet_coeffs;
                                
                                f_z[li][idx] = 0.5*(f_z[li-1][idx_z_B] + f_z[li-1][idx_z_F]);
                                
                                w_z[li][idx] = -0.5*(f_z[li-1][idx_z_B] - 2.0*f_z[li-1][idx] +
                                    f_z[li-1][idx_z_F]);
                            }
                        }
                    }
                    
                    if (compute_variable_local_means)
                    {
                        for (int k = 0; k < interior_dim_2; k++)
                        {
                            for (int j = 0; j < interior_dim_1; j++)
                            {
#ifdef HAMERS_ENABLE_SIMD
                                #pragma omp simd
#endif
                                for (int i = 0; i < interior_dim_0; i++)
                                {
                                    // Compute the linear indices.
                                    const int idx = (i + num_ghosts_0_wavelet_coeffs) +
                                        (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                        (k + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                            ghostcell_dim_1_wavelet_coeffs;
                                    
                                    const int idx_x_L = (i - offset + num_ghosts_0_wavelet_coeffs) +
                                        (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                        (k + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                            ghostcell_dim_1_wavelet_coeffs;
                                    
                                    const int idx_x_R = (i + offset + num_ghosts_0_wavelet_coeffs) +
                                        (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                        (k + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                            ghostcell_dim_1_wavelet_coeffs;
                                    
                                    const int idx_y_B = (i + num_ghosts_0_wavelet_coeffs) +
                                        (j - offset + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                        (k + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                            ghostcell_dim_1_wavelet_coeffs;
                                    
                                    const int idx_y_T = (i + num_ghosts_0_wavelet_coeffs) +
                                        (j + offset + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                        (k + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                            ghostcell_dim_1_wavelet_coeffs;
                                    
                                    const int idx_z_B = (i + num_ghosts_0_wavelet_coeffs) +
                                        (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                        (k - offset + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                            ghostcell_dim_1_wavelet_coeffs;
                                    
                                    const int idx_z_F = (i + num_ghosts_0_wavelet_coeffs) +
                                        (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                        (k + offset + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                            ghostcell_dim_1_wavelet_coeffs;
                                    
                                    double f_mean_x = f_x[li-1][idx_x_L] + 2.0*f_x[li-1][idx] +
                                        f_x[li-1][idx_x_R];
                                    
                                    double f_mean_y = f_y[li-1][idx_y_B] + 2.0*f_y[li-1][idx] +
                                        f_y[li-1][idx_y_T];
                                    
                                    double f_mean_z = f_z[li-1][idx_z_B] + 2.0*f_z[li-1][idx] +
                                        f_z[li-1][idx_z_F];
                                    
                                    f_mean[li][idx] = 0.5*sqrt(f_mean_x*f_mean_x + f_mean_y*f_mean_y +
                                        f_mean_z*f_mean_z);
                                }
                            }
                        }
                    }
                    
                    break;
                }
                case 4:
                {
                    /*
                     * Compute scaling and wavelet coefficients in the x-direction.
                     */
                    
                    // Compute the starting and ending indices.
                    const int start_index_x_i = -d_num_wavelet_ghosts[0] + 2*offset;
                    const int end_index_x_i = interior_dim_0 + d_num_wavelet_ghosts[0] - 2*offset;
                    const int start_index_x_j = -d_num_wavelet_ghosts[1];
                    const int end_index_x_j = interior_dim_1 + d_num_wavelet_ghosts[1];
                    const int start_index_x_k = -d_num_wavelet_ghosts[2];
                    const int end_index_x_k = interior_dim_2 + d_num_wavelet_ghosts[2];
                    
                    for (int k = start_index_x_k; k < end_index_x_k; k++)
                    {
                        for (int j = start_index_x_j; j < end_index_x_j; j++)
                        {
#ifdef HAMERS_ENABLE_SIMD
                            #pragma omp simd
#endif
                            for (int i = start_index_x_i; i < end_index_x_i; i++)
                            {
                                // Compute the linear indices.
                                const int idx = (i + num_ghosts_0_wavelet_coeffs) +
                                    (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                    (k + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                        ghostcell_dim_1_wavelet_coeffs;
                                
                                const int idx_x_LL = (i - 2*offset + num_ghosts_0_wavelet_coeffs) +
                                    (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                    (k + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                        ghostcell_dim_1_wavelet_coeffs;
                                
                                const int idx_x_L = (i - offset + num_ghosts_0_wavelet_coeffs) +
                                    (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                    (k + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                        ghostcell_dim_1_wavelet_coeffs;
                                
                                const int idx_x_R = (i + offset + num_ghosts_0_wavelet_coeffs) +
                                    (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                    (k + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                        ghostcell_dim_1_wavelet_coeffs;
                                
                                const int idx_x_RR = (i + 2*offset + num_ghosts_0_wavelet_coeffs) +
                                    (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                    (k + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                        ghostcell_dim_1_wavelet_coeffs;
                                
                                f_x[li][idx] = 1.0/6.0*(-f_x[li-1][idx_x_LL] + 4.0*f_x[li-1][idx_x_L] +
                                    4.0*f_x[li-1][idx_x_R] - f_x[li-1][idx_x_RR]);
                                
                                w_x[li][idx] = 1.0/6.0*(f_x[li-1][idx_x_LL] - 4.0*f_x[li-1][idx_x_L] +
                                    6.0*f_x[li-1][idx] - 4.0*f_x[li-1][idx_x_R] + f_x[li-1][idx_x_RR]);
                            }
                        }
                    }
                    
                    /*
                     * Compute scaling and wavelet coefficients in the y-direction.
                     */
                    
                    // Compute the starting and ending indices.
                    const int start_index_y_i = -d_num_wavelet_ghosts[0];
                    const int end_index_y_i = interior_dim_0 + d_num_wavelet_ghosts[0];
                    const int start_index_y_j = -d_num_wavelet_ghosts[1] + 2*offset;
                    const int end_index_y_j = interior_dim_1 + d_num_wavelet_ghosts[1] - 2*offset;
                    const int start_index_y_k = -d_num_wavelet_ghosts[2];
                    const int end_index_y_k = interior_dim_2 + d_num_wavelet_ghosts[2];
                    
                    for (int k = start_index_y_k; k < end_index_y_k; k++)
                    {
                        for (int j = start_index_y_j; j < end_index_y_j; j++)
                        {
#ifdef HAMERS_ENABLE_SIMD
                            #pragma omp simd
#endif
                            for (int i = start_index_y_i; i < end_index_y_i; i++)
                            {
                                // Compute the linear indices.
                                const int idx = (i + num_ghosts_0_wavelet_coeffs) +
                                    (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                    (k + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                        ghostcell_dim_1_wavelet_coeffs;
                                
                                const int idx_y_BB = (i + num_ghosts_0_wavelet_coeffs) +
                                    (j - 2*offset + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                    (k + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                        ghostcell_dim_1_wavelet_coeffs;
                                
                                const int idx_y_B = (i + num_ghosts_0_wavelet_coeffs) +
                                    (j - offset + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                    (k + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                        ghostcell_dim_1_wavelet_coeffs;
                                
                                const int idx_y_T = (i + num_ghosts_0_wavelet_coeffs) +
                                    (j + offset + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                    (k + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                        ghostcell_dim_1_wavelet_coeffs;
                                
                                const int idx_y_TT = (i + num_ghosts_0_wavelet_coeffs) +
                                    (j + 2*offset + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                    (k + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                        ghostcell_dim_1_wavelet_coeffs;
                                
                                f_y[li][idx] = 1.0/6.0*(-f_y[li-1][idx_y_BB] + 4.0*f_y[li-1][idx_y_B] +
                                    4.0*f_y[li-1][idx_y_T] - f_y[li-1][idx_y_TT]);
                                
                                w_y[li][idx] = 1.0/6.0*(f_y[li-1][idx_y_BB] - 4.0*f_y[li-1][idx_y_B] +
                                    6.0*f_y[li-1][idx] - 4.0*f_y[li-1][idx_y_T] + f_y[li-1][idx_y_TT]);
                            }
                        }
                    }
                    
                    /*
                     * Compute scaling and wavelet coefficients in the z-direction.
                     */
                    
                    // Compute the starting and ending indices.
                    const int start_index_z_i = -d_num_wavelet_ghosts[0];
                    const int end_index_z_i = interior_dim_0 + d_num_wavelet_ghosts[0];
                    const int start_index_z_j = -d_num_wavelet_ghosts[1];
                    const int end_index_z_j = interior_dim_1 + d_num_wavelet_ghosts[1];
                    const int start_index_z_k = -d_num_wavelet_ghosts[2] + 2*offset;
                    const int end_index_z_k = interior_dim_2 + d_num_wavelet_ghosts[2] - 2*offset;
                    
                    for (int k = start_index_z_k; k < end_index_z_k; k++)
                    {
                        for (int j = start_index_z_j; j < end_index_z_j; j++)
                        {
#ifdef HAMERS_ENABLE_SIMD
                            #pragma omp simd
#endif
                            for (int i = start_index_z_i; i < end_index_z_i; i++)
                            {
                                // Compute indices.
                                const int idx = (i + num_ghosts_0_wavelet_coeffs) +
                                    (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                    (k + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                        ghostcell_dim_1_wavelet_coeffs;
                                
                                const int idx_z_BB = (i + num_ghosts_0_wavelet_coeffs) +
                                    (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                    (k - 2*offset + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                        ghostcell_dim_1_wavelet_coeffs;
                                
                                const int idx_z_B = (i + num_ghosts_0_wavelet_coeffs) +
                                    (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                    (k - offset + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                        ghostcell_dim_1_wavelet_coeffs;
                                
                                const int idx_z_F = (i + num_ghosts_0_wavelet_coeffs) +
                                    (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                    (k + offset + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                        ghostcell_dim_1_wavelet_coeffs;
                                
                                const int idx_z_FF = (i + num_ghosts_0_wavelet_coeffs) +
                                    (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                    (k + 2*offset + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                        ghostcell_dim_1_wavelet_coeffs;
                                
                                f_z[li][idx] = 1.0/6.0*(-f_z[li-1][idx_z_BB] + 4.0*f_z[li-1][idx_z_B] +
                                    4.0*f_z[li-1][idx_z_F] - f_z[li-1][idx_z_FF]);
                                
                                w_z[li][idx] = 1.0/6.0*(f_z[li-1][idx_z_BB] - 4.0*f_z[li-1][idx_z_B] +
                                    6.0*f_z[li-1][idx] - 4.0*f_z[li-1][idx_z_F] + f_z[li-1][idx_z_FF]);
                            }
                        }
                    }
                    
                    if (compute_variable_local_means)
                    {
                        for (int k = 0; k < interior_dim_2; k++)
                        {
                            for (int j = 0; j < interior_dim_1; j++)
                            {
#ifdef HAMERS_ENABLE_SIMD
                                #pragma omp simd
#endif
                                for (int i = 0; i < interior_dim_0; i++)
                                {
                                    // Compute the linear indices.
                                    const int idx = (i + num_ghosts_0_wavelet_coeffs) +
                                        (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                        (k + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                            ghostcell_dim_1_wavelet_coeffs;
                                    
                                    const int idx_x_LL = (i - 2*offset + num_ghosts_0_wavelet_coeffs) +
                                        (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                        (k + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                            ghostcell_dim_1_wavelet_coeffs;
                                    
                                    const int idx_x_L = (i - offset + num_ghosts_0_wavelet_coeffs) +
                                        (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                        (k + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                            ghostcell_dim_1_wavelet_coeffs;
                                    
                                    const int idx_x_R = (i + offset + num_ghosts_0_wavelet_coeffs) +
                                        (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                        (k + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                            ghostcell_dim_1_wavelet_coeffs;
                                    
                                    const int idx_x_RR = (i + 2*offset + num_ghosts_0_wavelet_coeffs) +
                                        (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                        (k + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                            ghostcell_dim_1_wavelet_coeffs;
                                    
                                    const int idx_y_BB = (i + num_ghosts_0_wavelet_coeffs) +
                                        (j - 2*offset + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                        (k + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                            ghostcell_dim_1_wavelet_coeffs;
                                    
                                    const int idx_y_B = (i + num_ghosts_0_wavelet_coeffs) +
                                        (j - offset + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                        (k + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                            ghostcell_dim_1_wavelet_coeffs;
                                    
                                    const int idx_y_T = (i + num_ghosts_0_wavelet_coeffs) +
                                        (j + offset + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                        (k + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                            ghostcell_dim_1_wavelet_coeffs;
                                    
                                    const int idx_y_TT = (i + num_ghosts_0_wavelet_coeffs) +
                                        (j + 2*offset + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                        (k + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                            ghostcell_dim_1_wavelet_coeffs;
                                    
                                    const int idx_z_BB = (i + num_ghosts_0_wavelet_coeffs) +
                                        (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                        (k - 2*offset + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                            ghostcell_dim_1_wavelet_coeffs;
                                    
                                    const int idx_z_B = (i + num_ghosts_0_wavelet_coeffs) +
                                        (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                        (k - offset + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                            ghostcell_dim_1_wavelet_coeffs;
                                    
                                    const int idx_z_F = (i + num_ghosts_0_wavelet_coeffs) +
                                        (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                        (k + offset + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                            ghostcell_dim_1_wavelet_coeffs;
                                    
                                    const int idx_z_FF = (i + num_ghosts_0_wavelet_coeffs) +
                                        (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                                        (k + 2*offset + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                            ghostcell_dim_1_wavelet_coeffs;
                                    
                                    double f_mean_x = f_x[li-1][idx_x_LL] + 4.0*f_x[li-1][idx_x_L] +
                                        6.0*f_x[li-1][idx] + 4.0*f_x[li-1][idx_x_R] + f_x[li-1][idx_x_RR];
                                    
                                    double f_mean_y = f_y[li-1][idx_y_BB] + 4.0*f_y[li-1][idx_y_B] +
                                        6.0*f_y[li-1][idx] + 4.0*f_y[li-1][idx_y_T] + f_y[li-1][idx_y_TT];
                                    
                                    double f_mean_z = f_z[li-1][idx_z_BB] + 4.0*f_z[li-1][idx_z_B] +
                                        6.0*f_z[li-1][idx] + 4.0*f_z[li-1][idx_z_F] + f_z[li-1][idx_z_FF];
                                    
                                    f_mean[li][idx] = 1.0/6.0*sqrt(f_mean_x*f_mean_x + f_mean_y*f_mean_y +
                                        f_mean_z*f_mean_z);
                                }
                            }
                        }
                    }
                    
                    break;
                }
                default:
                {
                    TBOX_ERROR(d_object_name
                        << ": "
                        << "number of vanishing moments = "
                        << d_k
                        << " not supported."
                        << std::endl);
                }
            }
        }
        
        
        /*
         * Compute the three-dimensional wavelet coefficients from the one-dimensional wavelet coefficients.
         */
        
        for (int li = 0; li < d_num_level; li++)
        {
            const int offset = pow(2, li + 1);
            
            // Compute the starting and ending indices.
            const int start_index_i = -d_p*offset;
            const int end_index_i = interior_dim_0 + d_q*offset;
            const int start_index_j = -d_p*offset;
            const int end_index_j = interior_dim_1 + d_q*offset;
            const int start_index_k = -d_p*offset;
            const int end_index_k = interior_dim_2 + d_q*offset;
            
            for (int k = start_index_k; k < end_index_k; k++)
            {
                for (int j = start_index_j; j < end_index_j; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = start_index_i; i < end_index_i; i++)
                    {
                        const int idx = (i + num_ghosts_0_wavelet_coeffs) +
                            (j + num_ghosts_1_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs +
                            (k + num_ghosts_2_wavelet_coeffs)*ghostcell_dim_0_wavelet_coeffs*
                                ghostcell_dim_1_wavelet_coeffs;
                        
                        w[li][idx] = sqrt(w_x[li][idx]*w_x[li][idx] + w_y[li][idx]*w_y[li][idx] +
                            w_z[li][idx]*w_z[li][idx]);
                    }
                }
            }
        }
    }
}


/*
 * Smooth the given cell data in different directions.
 */
boost::shared_ptr<pdat::CellData<double> >
WaveletTransformHarten::smoothCellData(
    hier::Patch& patch,
    boost::shared_ptr<pdat::CellData<double> > cell_data,
    int depth)
{
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the number of ghost cells of the cell data, wavelet coefficients and local means.
    const hier::IntVector num_ghosts_cell_data = cell_data->getGhostCellWidth();
    
    // Get the dimensions of boxes that cover interior of patch plus ghost cells.
    const hier::Box ghost_box_cell_data = cell_data->getGhostBox();
    const hier::IntVector ghostcell_dims_cell_data = ghost_box_cell_data.numberCells();
    
    // Get the pointer to the desired depth component of the given cell data.
    double* f = cell_data->getPointer(depth);
    
    // Allocate memory for the smoothed cell data.
    boost::shared_ptr<pdat::CellData<double> > smoothed_cell_data(
        new pdat::CellData<double>(interior_box, d_dim.getValue(), num_ghosts_cell_data));
    
    if (d_dim == tbox::Dimension(1))
    {
        const int num_ghosts_0_cell_data = num_ghosts_cell_data[0];
        
        // Get the pointer to the cell data smoothed in the x-direction.
        double* f_smoothed = smoothed_cell_data->getPointer(0);
        
        // Compute the starting and ending indices.
        const int start_index_x_i = -num_ghosts_cell_data[0];
        const int end_index_x_i = interior_dims[0] + num_ghosts_cell_data[0];
        
        for (int i = start_index_x_i; i < end_index_x_i; i++)
        {
            // Compute the linear index.
            const int idx = i + num_ghosts_0_cell_data;
            
            f_smoothed[idx] = f[idx];
            
            int count = 1;
            
            // Sum over the neighboring cells.
            for (int ii = fmax(i - 1, start_index_x_i); ii < fmin(i + 2, end_index_x_i); ii++)
            {
                if (ii != i)
                {
                    const int idx_ii = ii + num_ghosts_cell_data[0];
                        
                    f_smoothed[idx] += f[idx_ii];
                    count++;
                }
            }
            
            // Compute the average value.
            f_smoothed[idx] = f_smoothed[idx]/count;
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        int start_index_i, end_index_i,
            start_index_j, end_index_j;
        
        // Get the pointer to the cell data smoothed in the x-direction.
        double* f_smoothed = smoothed_cell_data->getPointer(0);
        
        // Compute the starting and ending indices.
        start_index_i = -num_ghosts_cell_data[0];
        end_index_i   = interior_dims[0] + num_ghosts_cell_data[0];
        start_index_j = -num_ghosts_cell_data[1];
        end_index_j   = interior_dims[1] + num_ghosts_cell_data[1];
        
        for (int j = start_index_j; j < end_index_j; j++)
        {
            for (int i = start_index_i; i < end_index_i; i++)
            {
                // Compute the linear index.
                const int idx = (i + num_ghosts_cell_data[0]) +
                    (j + num_ghosts_cell_data[1])*ghostcell_dims_cell_data[0];
                
                f_smoothed[idx] = f[idx];
                
                int count = 1;
                
                // Sum over the neighboring cells.
                for (int ii = fmax(i - 1, start_index_i); ii < fmin(i + 2, end_index_i); ii++)
                {
                    if (ii != i)
                    {
                        const int idx_ii = (ii + num_ghosts_cell_data[0]) +
                            (j + num_ghosts_cell_data[1])*ghostcell_dims_cell_data[0];
                            
                        f_smoothed[idx] += f[idx_ii];
                        count++;
                    }
                }
                
                // Compute the average value.
                f_smoothed[idx] = f_smoothed[idx]/count;
            }
        }
        
        // Get the pointer to the cell data smoothed in the y-direction.
        f_smoothed = smoothed_cell_data->getPointer(1);
        
        // Compute the starting and ending indices.
        start_index_i = -num_ghosts_cell_data[0];
        end_index_i   = interior_dims[0] + num_ghosts_cell_data[0];
        start_index_j = -num_ghosts_cell_data[1];
        end_index_j   = interior_dims[1] + num_ghosts_cell_data[1];
        
        for (int i = start_index_i; i < end_index_i; i++)
        {
            for (int j = start_index_j; j < end_index_j; j++)
            {
                // Compute the linear index.
                const int idx = (i + num_ghosts_cell_data[0]) +
                    (j + num_ghosts_cell_data[1])*ghostcell_dims_cell_data[0];
                
                f_smoothed[idx] = f[idx];
                
                int count = 1;
                
                // Sum over the neighboring cells.
                for (int jj = fmax(j - 1, start_index_j); jj < fmin(j + 2, end_index_j); jj++)
                {
                    if (jj != j)
                    {
                        const int idx_jj = (i + num_ghosts_cell_data[0]) +
                            (jj + num_ghosts_cell_data[1])*ghostcell_dims_cell_data[0];
                            
                        f_smoothed[idx] += f[idx_jj];
                        count++;
                    }
                }
                
                // Compute the average value.
                f_smoothed[idx] = f_smoothed[idx]/count;
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        int start_index_i, end_index_i,
            start_index_j, end_index_j,
            start_index_k, end_index_k;
        
        // Get the pointer to the cell data smoothed in the x-direction.
        double* f_smoothed = smoothed_cell_data->getPointer(0);
        
        // Compute the starting and ending indices.
        start_index_i = -num_ghosts_cell_data[0];
        end_index_i   = interior_dims[0] + num_ghosts_cell_data[0];
        start_index_j = -num_ghosts_cell_data[1];
        end_index_j   = interior_dims[1] + num_ghosts_cell_data[1];
        start_index_k = -num_ghosts_cell_data[2];
        end_index_k   = interior_dims[2] + num_ghosts_cell_data[2];
        
        for (int k = start_index_k; k < end_index_k; k++)
        {
            for (int j = start_index_j; j < end_index_j; j++)
            {
                for (int i = start_index_i; i < end_index_i; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_cell_data[0]) +
                        (j + num_ghosts_cell_data[1])*ghostcell_dims_cell_data[0] +
                        (k + num_ghosts_cell_data[2])*ghostcell_dims_cell_data[0]*
                            ghostcell_dims_cell_data[1];
                    
                    f_smoothed[idx] = f[idx];
                    
                    int count = 1;
                    
                    // Sum over the neighboring cells.
                    for (int ii = fmax(i - 1, start_index_i); ii < fmin(i + 2, end_index_i); ii++)
                    {
                        if (ii != i)
                        {
                            const int idx_ii = (ii + num_ghosts_cell_data[0]) +
                                (j + num_ghosts_cell_data[1])*ghostcell_dims_cell_data[0] +
                                (k + num_ghosts_cell_data[2])*ghostcell_dims_cell_data[0]*
                                    ghostcell_dims_cell_data[1];
                                
                            f_smoothed[idx] += f[idx_ii];
                            count++;
                        }
                    }
                    
                    // Compute the average value.
                    f_smoothed[idx] = f_smoothed[idx]/count;
                }
            }
        }
        
        // Get the pointer to the cell data smoothed in the y-direction.
        f_smoothed = smoothed_cell_data->getPointer(1);
        
        // Compute the starting and ending indices.
        start_index_i = -num_ghosts_cell_data[0];
        end_index_i   = interior_dims[0] + num_ghosts_cell_data[0];
        start_index_j = -num_ghosts_cell_data[1];
        end_index_j   = interior_dims[1] + num_ghosts_cell_data[1];
        start_index_k = -num_ghosts_cell_data[2];
        end_index_k   = interior_dims[2] + num_ghosts_cell_data[2];
        
        for (int i = start_index_i; i < end_index_i; i++)
        {
            for (int k = start_index_k; k < end_index_k; k++)
            {
                for (int j = start_index_j; j < end_index_j; j++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_cell_data[0]) +
                        (j + num_ghosts_cell_data[1])*ghostcell_dims_cell_data[0] +
                        (k + num_ghosts_cell_data[2])*ghostcell_dims_cell_data[0]*
                            ghostcell_dims_cell_data[1];
                    
                    f_smoothed[idx] = f[idx];
                    
                    int count = 1;
                    
                    for (int jj = fmax(j - 1, start_index_j); jj < fmin(j + 2, end_index_j); jj++)
                    {
                        if (jj != j)
                        {
                            const int idx_jj = (i + num_ghosts_cell_data[0]) +
                                (jj + num_ghosts_cell_data[1])*ghostcell_dims_cell_data[0] +
                                (k + num_ghosts_cell_data[2])*ghostcell_dims_cell_data[0]*
                                    ghostcell_dims_cell_data[1];
                                
                            f_smoothed[idx] += f[idx_jj];
                            count++;
                        }
                    }
                    
                    // Compute the average value.
                    f_smoothed[idx] = f_smoothed[idx]/count;
                }
            }
        }
        
        // Get the pointer to the cell data smoothed in the z-direction.
        f_smoothed = smoothed_cell_data->getPointer(2);
        
        // Compute the starting and ending indices.
        start_index_i = -num_ghosts_cell_data[0];
        end_index_i   = interior_dims[0] + num_ghosts_cell_data[0];
        start_index_j = -num_ghosts_cell_data[1];
        end_index_j   = interior_dims[1] + num_ghosts_cell_data[1];
        start_index_k = -num_ghosts_cell_data[2];
        end_index_k   = interior_dims[2] + num_ghosts_cell_data[2];
        
        for (int j = start_index_j; j < end_index_j; j++)
        {
            for (int i = start_index_i; i < end_index_i; i++)
            {
                for (int k = start_index_k; k < end_index_k; k++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_cell_data[0]) +
                        (j + num_ghosts_cell_data[1])*ghostcell_dims_cell_data[0] +
                        (k + num_ghosts_cell_data[2])*ghostcell_dims_cell_data[0]*
                            ghostcell_dims_cell_data[1];
                    
                    f_smoothed[idx] = f[idx];
                    
                    int count = 1;
                    
                    for (int kk = fmax(k - 1, start_index_k); kk < fmin(k + 2, end_index_k); kk++)
                    {
                        if (kk != k)
                        {
                            const int idx_kk = (i + num_ghosts_cell_data[0]) +
                                (j + num_ghosts_cell_data[1])*ghostcell_dims_cell_data[0] +
                                (kk + num_ghosts_cell_data[2])*ghostcell_dims_cell_data[0]*
                                    ghostcell_dims_cell_data[1];
                                
                            f_smoothed[idx] += f[idx_kk];
                            count++;
                        }
                    }
                    
                    // Compute the average value.
                    f_smoothed[idx] = f_smoothed[idx]/count;
                }
            }
        }
    }
    
    return smoothed_cell_data;
}
