#ifndef CONVECTIVE_FLUX_RECONSTRUCTOR_FIRST_ORDER_HLLC_HPP
#define CONVECTIVE_FLUX_RECONSTRUCTOR_FIRST_ORDER_HLLC_HPP

#include "flow_model/convective_flux_reconstructor/ConvectiveFluxReconstructor.hpp"

#include "flow_model/Riemann_solver/RiemannSolverHLLC.hpp"

class ConvectiveFluxReconstructorFirstOrderHLLC: public ConvectiveFluxReconstructor
{
    public:
        ConvectiveFluxReconstructorFirstOrderHLLC(
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
        void computeConvectiveFluxesAndSources(
            hier::Patch& patch,
            const double time,
            const double dt,
            const int RK_step_number,
            const boost::shared_ptr<hier::VariableContext>& data_context);
    
    private:
        RiemannSolverHLLC d_Riemann_solver;
};

#endif /* CONVECTIVE_FLUX_RECONSTRUCTOR_FIRST_ORDER_HLLC_HPP */
