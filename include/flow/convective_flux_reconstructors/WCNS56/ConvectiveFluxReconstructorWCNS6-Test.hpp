#ifndef CONVECTIVE_FLUX_RECONSTRUCTOR_WCNS6_TEST_HPP
#define CONVECTIVE_FLUX_RECONSTRUCTOR_WCNS6_TEST_HPP

#include "SAMRAI/pdat/SideVariable.h"

#include "flow/convective_flux_reconstructors/ConvectiveFluxReconstructor.hpp"
#include "util/Directions.hpp"

#include "boost/multi_array.hpp"

class ConvectiveFluxReconstructorWCNS6_Test: public ConvectiveFluxReconstructor
{
    public:
        ConvectiveFluxReconstructorWCNS6_Test(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_eqn,
            const int& num_species,
            const boost::shared_ptr<FlowModel>& flow_model,
            const boost::shared_ptr<tbox::Database>& convective_flux_reconstructor_db);
        
        ~ConvectiveFluxReconstructorWCNS6_Test();
        
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
            const boost::shared_ptr<pdat::FaceVariable<double> >& variable_convective_flux,
            const boost::shared_ptr<pdat::CellVariable<double> >& variable_source,
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
         * Compute sigma's.
         */
        void
        computeSigma(
            double& sigma,
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
         * Compute sigma's.
         */
        void
        computeSigma(
            double& sigma,
            const boost::multi_array_ref<const double*, 2>::const_array_view<1>::type& U_array);
        
        /*
         * Compute beta's.
         */
        void
        computeBeta(
            std::vector<double>& beta,
            const boost::multi_array_ref<const double*, 2>::const_array_view<1>::type& U_array);
        
        /*
         * Compute beta_tilde's.
         */
        void
        computeBetaTilde(   
            std::vector<double>& beta_tilde,
            const boost::multi_array_ref<const double*, 2>::const_array_view<1>::type& U_array);
        
        /*
         * Compute sigma's.
         */
        void
        computeSigma(
            double& sigma,
            const std::vector<double*>& U_array,
            const int& idx_face);
        
        /*
         * Compute beta's.
         */
        void
        computeBeta(
            std::vector<double>& beta,
            const std::vector<double*>& U_array,
            const int& idx_face);
        
        /*
         * Compute beta_tilde's.
         */
        void
        computeBetaTilde(
            std::vector<double>& beta_tilde,
            const std::vector<double*>& U_array,
            const int& idx_face);
        
        /*
         * Perform WENO interpolation.
         */
        void
        performWENOInterpolation_new(
            double* U_minus,
            double* U_plus,
            const std::vector<double*>& U_array,
            const int& idx_face);
        
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
            const DIRECTION::TYPE& direction);
        
        void
        computeSigma(
            std::vector<boost::shared_ptr<pdat::SideData<double> > >& variables_sigma,
            const std::vector<std::vector<boost::shared_ptr<pdat::SideData<double> > > >& variables_array);

        /*
         * Constants used by the scheme.
         */
        double d_constant_C;
        int d_constant_p;
        int d_constant_q;
        double d_constant_alpha_beta;
        
        /*
         * Weights used in WENO interpolations.
         */
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
        std::vector<double> beta;
        std::vector<double> beta_tilde;
        
        /*
         * Timers interspersed throughout the class.
         */
        
        static boost::shared_ptr<tbox::Timer> t_characteristic_decomposition;
        static boost::shared_ptr<tbox::Timer> t_WENO_interpolation;
        static boost::shared_ptr<tbox::Timer> t_Riemann_solver;
        static boost::shared_ptr<tbox::Timer> t_reconstruct_flux;
        
};

#endif /* CONVECTIVE_FLUX_RECONSTRUCTOR_WCNS6_TEST_HPP */
