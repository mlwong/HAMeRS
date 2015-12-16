#ifndef CONVECTIVE_FLUX_RECONSTRUCTOR_WENO_JS5_LLF_HPP
#define CONVECTIVE_FLUX_RECONSTRUCTOR_WENO_JS5_LLF_HPP

#include "flow_model/convective_flux_reconstructor/ConvectiveFluxReconstructor.hpp"

#include "Directions.hpp"

#include "boost/multi_array.hpp"

#define EPSILON 1e-40

class ConvectiveFluxReconstructorWENO_JS5_LLF: public ConvectiveFluxReconstructor
{
    public:
        ConvectiveFluxReconstructorWENO_JS5_LLF(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geom,
            const FLOW_MODEL& flow_model,
            const int& num_eqn,
            const int& num_species,
            const boost::shared_ptr<EquationOfState>& equation_of_state,
            const boost::shared_ptr<tbox::Database>& shock_capturing_scheme_db);
        
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
         * Weights used in WENO reconstruction.
         */
        std::vector<double> d_weights_d;
        boost::multi_array<double, 2> d_weights_c;
        
        /*
         * Convert fluxes/conservative variables into characteristic variables.
         */
        void
        projectVariablesToCharacteristicFields(
            std::vector<double*> characteristic_variables,
            const std::vector<const double*> physical_variables,
            const boost::multi_array<const double*, 2> projection_matrix);
        
        /*
         * Convert characteristic variables into fluxes/conservative variables.
         */
        void
        projectCharacteristicVariablesToPhysicalFields(
            std::vector<double*> physical_variables,
            const std::vector<const double*> characteristic_variables,
            const boost::multi_array<const double*, 2> projection_matrix_inv);
        
        /*
         * Compute beta_plus's.
         */
        boost::multi_array<double, 2>
        computeBetaPlus(
            const boost::multi_array<double, 2>& G_plus_array);
        
        /*
         * Compute beta_minus's.
         */
        boost::multi_array<double, 2>
        computeBetaMinus(
            const boost::multi_array<double, 2>& G_minus_array);
        
        /*
         * Perform WENO reconstruction.
         */
        void
        performWENOReconstruction(
            std::vector<double>& G_plus,
            std::vector<double>& G_minus,
            const boost::multi_array<double, 2>& G_plus_array,
            const boost::multi_array<double, 2>& G_minus_array);
};

#endif /* CONVECTIVE_FLUX_RECONSTRUCTOR_WCNS_JS5_HLLC_HLL_HPP */
