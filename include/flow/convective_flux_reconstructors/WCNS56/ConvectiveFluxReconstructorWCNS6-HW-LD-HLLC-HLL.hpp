#ifndef CONVECTIVE_FLUX_RECONSTRUCTOR_WCNS6_HW_LD_HLLC_HLL_HPP
#define CONVECTIVE_FLUX_RECONSTRUCTOR_WCNS6_HW_LD_HLLC_HLL_HPP

#include "flow/convective_flux_reconstructors/WCNS56/ConvectiveFluxReconstructorWCNS56-HLLC-HLL.hpp"

class ConvectiveFluxReconstructorWCNS6_HW_LD_HLLC_HLL: public ConvectiveFluxReconstructorWCNS56
{
    public:
        ConvectiveFluxReconstructorWCNS6_HW_LD_HLLC_HLL(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_eqn,
            const int& num_species,
            const boost::shared_ptr<FlowModel>& flow_model,
            const boost::shared_ptr<tbox::Database>& convective_flux_reconstructor_db);
        
        ~ConvectiveFluxReconstructorWCNS6_HW_LD_HLLC_HLL() {}
        
        /*
         * Print all characteristics of the convective flux reconstruction class.
         */
        void
        printClassData(std::ostream& os) const;
        
        /*
         * Put the characteristics of the convective flux reconstruction class
         * into the restart database.
         */
        void
        putToRestart(
            const boost::shared_ptr<tbox::Database>& restart_db) const;
        
    private:
        /*
         * Compute both sigma's and TV's.
         */
        void
        computeSigmaAndTV(
            double& sigma,
            std::vector<double>& TV,
            const boost::multi_array_ref<double, 2>::const_array_view<1>::type& W_array);
        
        /*
         * Compute beta's.
         */
        void
        computeBeta(
            std::vector<double>& beta,
            const boost::multi_array_ref<double, 2>::const_array_view<1>::type& W_array);
        
        /*
         * Compute beta_tilde's.
         */
        void
        computeBetaTilde(
            std::vector<double>& beta_tilde,
            const boost::multi_array_ref<double, 2>::const_array_view<1>::type& W_array);
        
        /*
         * Perform WENO interpolation.
         */
        void
        performWENOInterpolation(
            std::vector<double>& U_minus,
            std::vector<double>& U_plus,
            const boost::multi_array<const double*, 2>& U_array,
            const hier::Index& cell_index_minus,
            const hier::Index& cell_index_plus,
            const DIRECTION& direction);
        
        /*
         * Constants used by the scheme.
         */
        double d_constant_C;
        int d_constant_p;
        int d_constant_q;
        double d_constant_alpha_TV;
        double d_constant_alpha_beta;
        
        /*
         * Weights used in WENO interpolations.
         */
        std::vector<double> d_weights_d;
        boost::multi_array<double, 2> d_weights_c;
        
        /*
         * Projection matrix and its inversed used in WENO interpolation.
         */
        boost::multi_array<double, 2> R_inv_intercell;
        boost::multi_array<double, 2> R_intercell;
        
        /*
         * Container to store the characteristic variables in WENO interpolation.
         */
        boost::multi_array<double, 2> W_array;
        
        /*
         * Containers to store the interpolated characteristic variables in WENO interpolation.
         */
        std::vector<double> W_minus;
        std::vector<double> W_plus;
        
        /*
         * Containers used in WENO interpolation.
         */
        std::vector<double> TV;
        std::vector<double> beta;
        std::vector<double> beta_tilde;
};

#endif /* CONVECTIVE_FLUX_RECONSTRUCTOR_WCNS6_HW_LD_HLLC_HLL_HPP */
