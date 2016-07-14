#include "flow/convective_flux_reconstructors/ConvectiveFluxReconstructorWCNS6-HW-LD-HLLC-HLL.hpp"

#define EPSILON 1e-40

ConvectiveFluxReconstructorWCNS6_HW_LD_HLLC_HLL::ConvectiveFluxReconstructorWCNS6_HW_LD_HLLC_HLL(
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
        TV(4),
        beta(4),
        beta_tilde(4)
{
    /*
     * Set the constants that are used in the scheme.
     */
    
    d_constant_C          = d_convective_flux_reconstructor_db->getDoubleWithDefault("constant_C", 1.0e10);
    d_constant_C          = d_convective_flux_reconstructor_db->getDoubleWithDefault("d_constant_C", d_constant_C);
    
    d_constant_p          = d_convective_flux_reconstructor_db->getIntegerWithDefault("constant_p", 2);
    d_constant_p          = d_convective_flux_reconstructor_db->getIntegerWithDefault("d_constant_p", d_constant_p);
    
    d_constant_q          = d_convective_flux_reconstructor_db->getIntegerWithDefault("constant_q", 4);
    d_constant_q          = d_convective_flux_reconstructor_db->getIntegerWithDefault("d_constant_q", d_constant_q);
    
    d_constant_alpha_TV   = d_convective_flux_reconstructor_db->getDoubleWithDefault("constant_alpha_TV", 2.0);
    d_constant_alpha_TV   = d_convective_flux_reconstructor_db->getDoubleWithDefault("d_constant_alpha_TV", d_constant_alpha_TV);
    
    d_constant_alpha_beta = d_convective_flux_reconstructor_db->getDoubleWithDefault("constant_alpha_beta", 50.0);
    d_constant_alpha_beta = d_convective_flux_reconstructor_db->getDoubleWithDefault("d_constant_alpha_beta", d_constant_alpha_beta);
    
    d_weights_d.reserve(4);
    d_weights_d.push_back(1.0/32.0);
    d_weights_d.push_back(15.0/32.0);
    d_weights_d.push_back(15.0/32.0);
    d_weights_d.push_back(1.0/32.0);
    
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
ConvectiveFluxReconstructorWCNS6_HW_LD_HLLC_HLL::printClassData(
    std::ostream& os) const
{
    os << "\nPrint ConvectiveFluxReconstructorWCNS6_HW_LD_HLLC_HLL object..."
       << std::endl;
    
    os << std::endl;
    
    os << "ConvectiveFluxReconstructorWCNS6_HW_LD_HLLC_HLL: this = "
       << (ConvectiveFluxReconstructorWCNS6_HW_LD_HLLC_HLL *)this
       << std::endl;
    os << "d_object_name = "
       << d_object_name
       << std::endl;
    os << "d_constant_C = "
       << d_constant_C
       << std::endl;
    os << "d_constant_p = "
       << d_constant_p
       << std::endl;
    os << "d_constant_q = "
       << d_constant_q
       << std::endl;
    os << "d_constant_alpha_TV = "
       << d_constant_alpha_TV
       << std::endl;
    os << "d_constant_alpha_beta = "
       << d_constant_alpha_beta
       << std::endl;
}


/*
 * Put the characteristics of the convective flux reconstruction class
 * into the restart database.
 */
void
ConvectiveFluxReconstructorWCNS6_HW_LD_HLLC_HLL::putToRestart(
   const boost::shared_ptr<tbox::Database>& restart_db) const
{
    restart_db->putDouble("d_constant_C", d_constant_C);
    restart_db->putInteger("d_constant_p", d_constant_p);
    restart_db->putInteger("d_constant_q", d_constant_q);
    restart_db->putDouble("d_constant_alpha_TV", d_constant_alpha_TV);
    restart_db->putDouble("d_constant_alpha_beta", d_constant_alpha_beta);
}


/*
 * Compute both sigma's and TV's.
 */
void
ConvectiveFluxReconstructorWCNS6_HW_LD_HLLC_HLL::computeSigmaAndTV(
    double& sigma,
    std::vector<double>& TV,
    const boost::multi_array_ref<double, 2>::const_array_view<1>::type& W_array)
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(TV.size()) == 4);
    TBOX_ASSERT(static_cast<int>(W_array.shape()[0]) == 6);
#endif
    
    double W_array_diff[5];
    for (int r = 0; r < 5; r++)
    {
        W_array_diff[r] = fabs(W_array[r + 1] - W_array[r]);
    }
    
    /*
     * Compute the sigma.
     */
    
    const double r_1 = 2*W_array_diff[1]*W_array_diff[0]/
        (W_array_diff[1]*W_array_diff[1] + W_array_diff[0]*W_array_diff[0] + EPSILON);
    
    const double r_2 = 2*W_array_diff[2]*W_array_diff[1]/
        (W_array_diff[2]*W_array_diff[2] + W_array_diff[1]*W_array_diff[1] + EPSILON);
    
    const double r_3 = 2*W_array_diff[3]*W_array_diff[2]/
        (W_array_diff[3]*W_array_diff[3] + W_array_diff[2]*W_array_diff[2] + EPSILON);
        
    const double r_4 = 2*W_array_diff[4]*W_array_diff[3]/
        (W_array_diff[4]*W_array_diff[4] + W_array_diff[3]*W_array_diff[3] + EPSILON);
    
    sigma = fmin(fmin(r_1, r_2), fmin(r_3, r_4));
    
    /*
     * Compute the TV.
     */
    
    TV[0] = W_array_diff[0] + W_array_diff[1];
    TV[1] = W_array_diff[1] + W_array_diff[2];
    TV[2] = W_array_diff[2] + W_array_diff[3];
    TV[3] = W_array_diff[3] + W_array_diff[4];
}


/*
 * Compute beta's.
 */
void
ConvectiveFluxReconstructorWCNS6_HW_LD_HLLC_HLL::computeBeta(
    std::vector<double>& beta,
    const boost::multi_array_ref<double, 2>::const_array_view<1>::type& W_array)
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(beta.size()) == 4);
    TBOX_ASSERT(static_cast<int>(W_array.shape()[0]) == 6);
#endif
    
    beta[0] = 1.0/3*(W_array[0]*(4*W_array[0] - 19*W_array[1] + 11*W_array[2]) +
        W_array[1]*(25*W_array[1] - 31*W_array[2]) + 10*W_array[2]*W_array[2]);
    
    beta[1] = 1.0/3*(W_array[1]*(4*W_array[1] - 13*W_array[2] + 5*W_array[3]) +
        13*W_array[2]*(W_array[2] - W_array[3]) + 4*W_array[3]*W_array[3]);
    
    beta[2] = 1.0/3*(W_array[2]*(10*W_array[2] - 31*W_array[3] + 11*W_array[4]) +
        W_array[3]*(25*W_array[3] - 19*W_array[4]) + 4*W_array[4]*W_array[4]);
    
    beta[3] = 1.0/232243200*(W_array[0]*(525910327*W_array[0] - 4562164630*W_array[1] +
        7799501420*W_array[2] - 6610694540*W_array[3] + 2794296070*W_array[4] -
        472758974*W_array[5]) + 5*W_array[1]*(2146987907*W_array[1] - 7722406988*W_array[2] +
        6763559276*W_array[3] - 2926461814*W_array[4] + 503766638*W_array[5]) +
        20*W_array[2]*(1833221603*W_array[2] - 3358664662*W_array[3] + 1495974539*W_array[4] -
        263126407*W_array[5]) + 20*W_array[3]*(1607794163*W_array[3] - 1486026707*W_array[4] +
        268747951*W_array[5]) +  5*W_array[4]*(1432381427*W_array[4] - 536951582*W_array[5]) +
        263126407*W_array[5]*W_array[5]);
}


/*
 * Compute beta_tilde's.
 */
void
ConvectiveFluxReconstructorWCNS6_HW_LD_HLLC_HLL::computeBetaTilde(
    std::vector<double>& beta_tilde,
    const boost::multi_array_ref<double, 2>::const_array_view<1>::type& W_array)
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(beta_tilde.size()) == 4);
    TBOX_ASSERT(static_cast<int>(W_array.shape()[0]) == 6);
#endif
    
    beta_tilde[0] = 1.0/3*(W_array[5]*(4*W_array[5] - 19*W_array[4] + 11*W_array[3]) +
        W_array[4]*(25*W_array[4] - 31*W_array[3]) + 10*W_array[3]*W_array[3]);
    
    beta_tilde[1] = 1.0/3*(W_array[4]*(4*W_array[4] - 13*W_array[3] + 5*W_array[2]) +
        13*W_array[3]*(W_array[3] - W_array[2]) + 4*W_array[2]*W_array[2]);
    
    beta_tilde[2] = 1.0/3*(W_array[3]*(10*W_array[3] - 31*W_array[2] + 11*W_array[1]) +
        W_array[2]*(25*W_array[2] - 19*W_array[1]) + 4*W_array[1]*W_array[1]);
    
    beta_tilde[3] = 1.0/232243200*(W_array[5]*(525910327*W_array[5] - 4562164630*W_array[4] +
        7799501420*W_array[3] - 6610694540*W_array[2] + 2794296070*W_array[1] -
        472758974*W_array[0]) + 5*W_array[4]*(2146987907*W_array[4] - 7722406988*W_array[3] +
        6763559276*W_array[2] - 2926461814*W_array[1] + 503766638*W_array[0]) +
        20*W_array[3]*(1833221603*W_array[3] - 3358664662*W_array[2] + 1495974539*W_array[1] -
        263126407*W_array[0]) + 20*W_array[2]*(1607794163*W_array[2] -
        1486026707*W_array[1] + 268747951*W_array[0]) + 5*W_array[1]*(1432381427*W_array[1] -
        536951582*W_array[0])+263126407*W_array[0]*W_array[0]);
}


/*
 * Perform WENO interpolation.
 */
void
ConvectiveFluxReconstructorWCNS6_HW_LD_HLLC_HLL::performWENOInterpolation(
    std::vector<double>& U_minus,
    std::vector<double>& U_plus,
    const boost::multi_array<const double*, 2>& U_array,
    const hier::Index& cell_index_minus,
    const hier::Index& cell_index_plus,
    const DIRECTION& direction)
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
    
    const double& C = d_constant_C;
    const int& p = d_constant_p;
    const int& q = d_constant_q;
    const double& alpha_TV = d_constant_alpha_TV;
    const double& alpha_beta = d_constant_alpha_beta;
    
    const std::vector<double>& d = d_weights_d;
    const boost::multi_array<double, 2>& c = d_weights_c;
    
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        boost::multi_array_ref<double, 2>::const_array_view<1>::type W_array_ei =
            W_array[boost::indices[boost::multi_array_ref<double, 2>::index_range()][ei]];
        
        // Compute sigma and TV's.
        double sigma;
        computeSigmaAndTV(sigma, TV, W_array_ei);
        
        // Compute beta's.
        computeBeta(beta, W_array_ei);
        computeBetaTilde(beta_tilde, W_array_ei);
        
        // Compute ratio of TV's.
        
        const double TV_ratio = fmax(fmax(fmax(TV[0], TV[1]), TV[2]), TV[3])/
            (fmin(fmin(fmin(TV[0], TV[1]), TV[2]), TV[3]) + EPSILON);
        
        // Compute ratio of beta's.
        
        const double beta_ratio = fmax(fmax(fmax(beta[0], beta[1]), beta[2]), beta[3])/
            (fmin(fmin(fmin(beta[0], beta[1]), beta[2]), beta[3]) + EPSILON);
        
        const double beta_tilde_ratio = fmax(fmax(fmax(beta_tilde[0], beta_tilde[1]), beta_tilde[2]), beta_tilde[3])/
            (fmin(fmin(fmin(beta_tilde[0], beta_tilde[1]), beta_tilde[2]), beta_tilde[3]) + EPSILON);
        
        /*
         * Compute W_minus of the current characteristic variable.
         */
        
        // Compute the reference smoothness indicators tau_6.
        const double tau_6 = fabs(beta[3] - 1.0/8*(beta[0] + beta[2] + 6*beta[1]));
        
        // Compute the weights alpha.
        double alpha[4];
        double alpha_sum = 0.0;
        
        if ((TV_ratio > alpha_TV) &&
            (beta_ratio > alpha_beta))
        {
            // Compute hybrid linear weights d_h.
            double d_h[4];
            const double sigma_to_p = pow(sigma, p);
            d_h[0] = (2.0 - sigma_to_p)/32.0;
            d_h[1] = 5.0*(4.0 - sigma_to_p)/32.0;
            d_h[2] = 5.0*(2.0 + sigma_to_p)/32.0;
            d_h[3] = sigma_to_p/32.0;
            
            for (int r = 0; r < 4; r++)
            {
                // Compute the weights alpha.
                alpha[r] = d_h[r]*(EPSILON + sigma_to_p*C + pow(tau_6/(beta[r] + EPSILON), q));
                
                // Sum up the weights alpha.
                alpha_sum += alpha[r];
            }
        }
        else
        {
            for (int r = 0; r < 4; r++)
            {
                // Compute the weights alpha.
                alpha[r] = d[r]*(C + pow(tau_6/(beta[r] + EPSILON), q));
                
                // Sum up the weights alpha.
                alpha_sum += alpha[r];
            }
        }
        
        W_minus[ei] = 0.0;
        // Compute the W_minus.
        for (int r = 0; r < 4; r++)
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
        
        // Compute the reference smoothness indicators tau_6_tilde.
        const double tau_6_tilde = fabs(beta_tilde[3] - 1.0/8*(beta_tilde[0] + beta_tilde[2] + 6*beta_tilde[1]));
        
        // Compute the weights alpha_tilde.
        double alpha_tilde[4];
        double alpha_tilde_sum = 0.0;
        
        if ((TV_ratio > alpha_TV)  &&
            (beta_tilde_ratio > alpha_beta))
        {
            // Compute hybrid linear weights d_h.
            const double sigma_to_p = pow(sigma, p);
            double d_h[4];
            d_h[0] = (2.0 - sigma_to_p)/32.0;
            d_h[1] = 5.0*(4.0 - sigma_to_p)/32.0;
            d_h[2] = 5.0*(2.0 + sigma_to_p)/32.0;
            d_h[3] = sigma_to_p/32.0;
            
            for (int r = 0; r < 4; r++)
            {
                // Compute the weights alpha_tilde.
                alpha_tilde[r] = d_h[r]*(EPSILON + sigma_to_p*C + pow(tau_6_tilde/(beta_tilde[r] + EPSILON), q));
                
                // Sum up the weights alpha.
                alpha_tilde_sum += alpha_tilde[r];
            }
        }
        else
        {
            for (int r = 0; r < 4; r++)
            {
                // Compute the weights alpha_tilde.
                alpha_tilde[r] = d[r]*(C + pow(tau_6_tilde/(beta_tilde[r] + EPSILON), q));
                
                // Sum up the weights alpha.
                alpha_tilde_sum += alpha_tilde[r];
            }
        }
        
        // Compute the W_plus.
        W_plus[ei] = 0.0;
        for (int r = 0; r < 4; r++)
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
