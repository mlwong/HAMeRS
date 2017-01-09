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
         * Compute local sigma's.
         */
        void
        computeLocalSigma(
            double* sigma,
            const std::vector<double*>& U_array,
            const int& idx_side);
        
        /*
         * Compute local beta's.
         */
        void
        computeLocalBeta(
            double* beta_0,
            double* beta_1,
            double* beta_2,
            double* beta_3,
            const std::vector<double*>& U_array,
            const int& idx_side);
        
        /*
         * Compute local beta_tilde's.
         */
        void
        computeLocalBetaTilde(
            double* beta_tilde_0,
            double* beta_tilde_1,
            double* beta_tilde_2,
            double* beta_tilde_3,
            const std::vector<double*>& U_array,
            const int& idx_side);
        
        /*
         * Perform local WENO interpolation.
         */
        void
        performLocalWENOInterpolation(
            double* U_minus,
            double* U_plus,
            const std::vector<double*>& U_array,
            const int& idx_side);
        
        /*
         * Perform WENO interpolation.
         */
        void
        performWENOInterpolation(
            std::vector<boost::shared_ptr<pdat::SideData<double> > >& variables_minus,
            std::vector<boost::shared_ptr<pdat::SideData<double> > >& variables_plus,
            const std::vector<std::vector<boost::shared_ptr<pdat::SideData<double> > > >& variables);
        
        /*
         * Constants used by the scheme.
         */
        double d_constant_C;
        int    d_constant_p;
        int    d_constant_q;
        double d_constant_alpha_tau;
        
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
        static boost::shared_ptr<tbox::Timer> t_compute_source;
        
};

#endif /* CONVECTIVE_FLUX_RECONSTRUCTOR_WCNS6_TEST_HPP */
