#ifndef CONVECTIVE_FLUX_RECONSTRUCTOR_WCNS_TEST_HPP
#define CONVECTIVE_FLUX_RECONSTRUCTOR_WCNS_TEST_HPP

#include "flow_model/convective_flux_reconstructor/ConvectiveFluxReconstructor.hpp"

#include "Directions.hpp"

#include "flow_model/Riemann_solver/RiemannSolverHLLC.hpp"
#include "flow_model/Riemann_solver/RiemannSolverHLLC_HLL.hpp"

#include "boost/multi_array.hpp"

class ConvectiveFluxReconstructorWCNS_Test: public ConvectiveFluxReconstructor
{
    public:
        ConvectiveFluxReconstructorWCNS_Test(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
            const FLOW_MODEL& flow_model,
            const int& num_eqn,
            const int& num_species,
            const boost::shared_ptr<EquationOfState>& equation_of_state,
            const boost::shared_ptr<tbox::Database>& shock_capturing_scheme_db);
        
        ~ConvectiveFluxReconstructorWCNS_Test();
        
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
        
        /*
         * Compute the convective fluxes and sources due to hyperbolization
         * of the equations.
         */
        void
        computeConvectiveFluxesAndSources(
            hier::Patch& patch,
            const double time,
            const double dt,
            const int RK_step_number,
            const boost::shared_ptr<hier::VariableContext>& data_context);
    
    private:
        /*
         * Transform physcial variables into characteristic variables.
         * W_array: Characteristic variables.
         * U_array: Physical variables.
         * R_inv_intercell: Projection matrix.
         */
        void
        projectPhysicalVariablesToCharacteristicFields(
            boost::multi_array<double, 2>& W_array,
            const boost::multi_array<const double*, 2>& U_array,
            const boost::multi_array<double, 2>& R_inv_intercell);
        
        /*
         * Transform characteristic variables into primitive variables.
         * U: Physcial variables.
         * W: Characteristic variables.
         * R_intercell: Inverse of projection matrix.
         */
        void
        projectCharacteristicVariablesToPhysicalFields(
            std::vector<double>& U,
            const std::vector<double>& W,
            const boost::multi_array<double, 2>& R_intercell);
        
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
         * Constants used to bound the mass fraction.
         */
        double d_Y_bnd_lo;
        double d_Y_bnd_up;
        
        /*
         * Constants used to bound the volume fraction.
         */
        double d_Z_bnd_lo;
        double d_Z_bnd_up;
        
        /*
         * Weights used in WENO interpolations.
         */
        boost::multi_array<double, 2> d_weights_c;
        
        /*
         * Riemann solvers used for computing mid-point fluxes.
         */
        RiemannSolverHLLC d_Riemann_solver_HLLC;
        RiemannSolverHLLC_HLL d_Riemann_solver_HLLC_HLL;
        
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
        
        /*
         * Timers interspersed throughout the class.
         */
        static boost::shared_ptr<tbox::Timer> t_compute_v_p_f;
        static boost::shared_ptr<tbox::Timer> t_compute_theta_omega;
        static boost::shared_ptr<tbox::Timer> t_compute_projection_matrices;
        static boost::shared_ptr<tbox::Timer> t_compute_midpoint_fluxes;
        static boost::shared_ptr<tbox::Timer> t_compute_fluxes;
        static boost::shared_ptr<tbox::Timer> t_miscellaneous;
        static boost::shared_ptr<tbox::Timer> t_allocate_memory;
        
        static boost::shared_ptr<tbox::Timer> t_projection;
        static boost::shared_ptr<tbox::Timer> t_WENO_interpolation;
        static boost::shared_ptr<tbox::Timer> t_Riemann_solver;
        
        static boost::shared_ptr<tbox::Timer> t_WENO_interpolation_compute_sigma_TV;
        static boost::shared_ptr<tbox::Timer> t_WENO_interpolation_compute_beta;
        static boost::shared_ptr<tbox::Timer> t_WENO_interpolation_compute_interpolated_values;
        static boost::shared_ptr<tbox::Timer> t_WENO_interpolation_compute_alpha;
        static boost::shared_ptr<tbox::Timer> t_WENO_interpolation_compute_omega;
        
};

#endif /* CONVECTIVE_FLUX_RECONSTRUCTOR_TEST_HPP */
