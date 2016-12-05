#include "flow/convective_flux_reconstructors/WCNS56/ConvectiveFluxReconstructorWCNS5-JS-HLLC-HLL.hpp"

#define EPSILON 1e-40

ConvectiveFluxReconstructorWCNS5_JS_HLLC_HLL::ConvectiveFluxReconstructorWCNS5_JS_HLLC_HLL(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
    const int& num_eqn,
    const int& num_species,
    const boost::shared_ptr<FlowModel>& flow_model,
    const boost::shared_ptr<tbox::Database>& convective_flux_reconstructor_db):
        ConvectiveFluxReconstructorWCNS56(
            object_name,
            dim,
            grid_geometry,
            num_eqn,
            num_species,
            flow_model,
            convective_flux_reconstructor_db),
        W_array(
            boost::extents[6][d_num_eqn],
            boost::fortran_storage_order()),
        W_minus(d_num_eqn),
        W_plus(d_num_eqn),
        beta(3),
        beta_tilde(3)
{
    /*
     * Set the constants that are used in the scheme.
     */
    
    d_constant_q = d_convective_flux_reconstructor_db->getIntegerWithDefault("constant_q", 2);
    d_constant_q = d_convective_flux_reconstructor_db->getIntegerWithDefault("d_constant_q", d_constant_q);
    
    d_weights_d.reserve(3);
    d_weights_d.push_back(1.0/16);
    d_weights_d.push_back(10.0/16);
    d_weights_d.push_back(5.0/16);
    
    d_weights_c.resize(boost::extents[4][3]);
    d_weights_c[0][0] = 3.0/8;
    d_weights_c[0][1] = -5.0/4;
    d_weights_c[0][2] = 15.0/8;
    d_weights_c[1][0] = -1.0/8;
    d_weights_c[1][1] = 3.0/4;
    d_weights_c[1][2] = 3.0/8;
    d_weights_c[2][0] = 3.0/8;
    d_weights_c[2][1] = 3.0/4;
    d_weights_c[2][2] = -1.0/8;
    d_weights_c[3][0] = 15.0/8;
    d_weights_c[3][1] = -5.0/4;
    d_weights_c[3][2] = 3.0/8;
}


/*
 * Print all characteristics of the convective flux reconstruction class.
 */
void
ConvectiveFluxReconstructorWCNS5_JS_HLLC_HLL::printClassData(
    std::ostream& os) const
{
    os << "\nPrint ConvectiveFluxReconstructorWCNS5_JS_HLLC_HLL object..."
       << std::endl;
    
    os << std::endl;
    
    os << "ConvectiveFluxReconstructorWCNS5_JS_HLLC_HLL: this = "
       << (ConvectiveFluxReconstructorWCNS5_JS_HLLC_HLL *)this
       << std::endl;
    os << "d_object_name = "
       << d_object_name
       << std::endl;
    os << "d_constant_q = "
       << d_constant_q
       << std::endl;
}


/*
 * Put the characteristics of the convective flux reconstruction class
 * into the restart database.
 */
void
ConvectiveFluxReconstructorWCNS5_JS_HLLC_HLL::putToRestart(
   const boost::shared_ptr<tbox::Database>& restart_db) const
{
    restart_db->putInteger("d_constant_q", d_constant_q);
}


/*
 * Compute beta's.
 */
void
ConvectiveFluxReconstructorWCNS5_JS_HLLC_HLL::computeBeta(
    std::vector<double>& beta,
    const boost::multi_array_ref<double, 2>::const_array_view<1>::type& W_array)
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(beta.size()) == 3);
    TBOX_ASSERT(static_cast<int>(W_array.shape()[0]) == 6);
#endif
    
    beta[0] = 1.0/3*(W_array[0]*(4*W_array[0] - 19*W_array[1] + 11*W_array[2]) +
        W_array[1]*(25*W_array[1] - 31*W_array[2]) + 10*W_array[2]*W_array[2]);
    
    beta[1] = 1.0/3*(W_array[1]*(4*W_array[1] - 13*W_array[2] + 5*W_array[3]) +
        13*W_array[2]*(W_array[2] - W_array[3]) + 4*W_array[3]*W_array[3]);
    
    beta[2] = 1.0/3*(W_array[2]*(10*W_array[2] - 31*W_array[3] + 11*W_array[4]) +
        W_array[3]*(25*W_array[3] - 19*W_array[4]) + 4*W_array[4]*W_array[4]);
}


/*
 * Compute beta_tilde's.
 */
void
ConvectiveFluxReconstructorWCNS5_JS_HLLC_HLL::computeBetaTilde(
    std::vector<double>& beta_tilde,
    const boost::multi_array_ref<double, 2>::const_array_view<1>::type& W_array)
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(beta_tilde.size()) == 3);
    TBOX_ASSERT(static_cast<int>(W_array.shape()[0]) == 6);
#endif
    
    beta_tilde[0] = 1.0/3*(W_array[5]*(4*W_array[5] - 19*W_array[4] + 11*W_array[3]) +
        W_array[4]*(25*W_array[4] - 31*W_array[3]) + 10*W_array[3]*W_array[3]);
    
    beta_tilde[1] = 1.0/3*(W_array[4]*(4*W_array[4] - 13*W_array[3] + 5*W_array[2]) +
        13*W_array[3]*(W_array[3] - W_array[2]) + 4*W_array[2]*W_array[2]);
    
    beta_tilde[2] = 1.0/3*(W_array[3]*(10*W_array[3] - 31*W_array[2] + 11*W_array[1]) +
        W_array[2]*(25*W_array[2] - 19*W_array[1]) + 4*W_array[1]*W_array[1]);
}


/*
 * Perform WENO interpolation.
 */
void
ConvectiveFluxReconstructorWCNS5_JS_HLLC_HLL::performWENOInterpolation(
    std::vector<double>& U_minus,
    std::vector<double>& U_plus,
    const boost::multi_array<const double*, 2>& U_array,
    const hier::Index& cell_index_minus,
    const hier::Index& cell_index_plus,
    const DIRECTION::TYPE& direction)
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(U_array.shape()[0]) == 6);
    TBOX_ASSERT(static_cast<int>(U_array.shape()[1]) == d_num_eqn);
    TBOX_ASSERT(static_cast<int>(U_minus.size()) == d_num_eqn);
    TBOX_ASSERT(static_cast<int>(U_plus.size()) == d_num_eqn);
#endif
    
    /*
     * Compute the projection matrix.
     * Transform the physical variables into the characteristic variables.
     */
    
    d_flow_model->computeLocalFaceProjectionMatrixOfPrimitiveVariables(
        R_inv_intercell,
        cell_index_minus,
        cell_index_plus,
        direction);
    
    projectPhysicalVariablesToCharacteristicFields(W_array, U_array, R_inv_intercell);
    
    /*
     * Perform the WENO interpolation.
     */
    
    const int& q = d_constant_q;
    
    const std::vector<double>& d = d_weights_d;
    const boost::multi_array<double, 2>& c = d_weights_c;
    
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        boost::multi_array_ref<double, 2>::const_array_view<1>::type W_array_ei =
            W_array[boost::indices[boost::multi_array_ref<double, 2>::index_range()][ei]];
        
        // Compute beta's.
        computeBeta(beta, W_array_ei);
        computeBetaTilde(beta_tilde, W_array_ei);
    
        /*
         * Compute W_minus of the current characteristic variable.
         */
        
        // Compute the weights alpha.
        double alpha[3];
        double alpha_sum = 0.0;
            
        for (int r = 0; r < 3; r++)
        {
            // Compute the weights alpha.
            alpha[r] = d[r]/pow((beta[r] + EPSILON), q);
            
            // Sum up the weights alpha.
            alpha_sum += alpha[r];
        }
    
        W_minus[ei] = 0.0;
        // Compute the W_minus.
        for (int r = 0; r < 3; r++)
        {
            // Compute the linear interpolated value.
            double W_minus_r = 0.0;
            for (int m = r; m < 3 + r; m++)
            {
                W_minus_r += c[r][m - r]*W_array[m][ei];
            }
            
            // Compute omega.
            const double omega = alpha[r]/alpha_sum;
            
            // Compute the nonlinear interpolated value.
            W_minus[ei] += omega*W_minus_r;
        }
    
        /*
         * Compute W_plus of the current characteristic variable.
         */
        
        // Compute the weights alpha_tilde.
        double alpha_tilde[3];
        double alpha_tilde_sum = 0.0;
        
        for (int r = 0; r < 3; r++)
        {
            // Compute the weights alpha_tilde.
            alpha_tilde[r] = d[r]/pow((beta_tilde[r] + EPSILON), q);
            
            // Sum up the weights alpha.
            alpha_tilde_sum += alpha_tilde[r];
        }
        
        // Compute the W_plus.
        W_plus[ei] = 0.0;
        for (int r = 0; r < 3; r++)
        {
            // Compute the linear interpolated value.
            double W_plus_r = 0.0;
            for (int m = r; m < 3 + r; m++)
            {
                W_plus_r += c[r][m - r]*W_array[6 - m - 1][ei];
            }
            
            // Compute omega_tilde;
            const double omega_tilde = alpha_tilde[r]/alpha_tilde_sum;
            
            // Compute the nonlinear interpolated value.
            W_plus[ei] += omega_tilde*W_plus_r;
        }
    }
    
    /*
     * Compute the inverse of projection matrix.
     * Transform the characteristic variables back to physcial variables.
     */
    
    d_flow_model->computeLocalFaceProjectionMatrixInverseOfPrimitiveVariables(
        R_intercell,
        cell_index_minus,
        cell_index_plus,
        direction);
    
    projectCharacteristicVariablesToPhysicalFields(U_minus, W_minus, R_intercell);
    projectCharacteristicVariablesToPhysicalFields(U_plus, W_plus, R_intercell);
}
