#ifndef CONVECTIVE_FLUX_RECONSTRUCTOR_LLF_HPP
#define CONVECTIVE_FLUX_RECONSTRUCTOR_LLF_HPP

#include "ConvectiveFluxReconstructor.hpp"

class ConvectiveFluxReconstructorLLF: public ConvectiveFluxReconstructor
{
    public:
        ConvectiveFluxReconstructorLLF(
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
         * of the equtions.
         */
        void computeConvectiveFluxAndSource(
            hier::Patch& patch,
            const double time,
            const double dt,
            const boost::shared_ptr<hier::VariableContext> data_context);
    
};

#endif /* CONVECTIVE_FLUX_RECONSTRUCTOR_HPP */
