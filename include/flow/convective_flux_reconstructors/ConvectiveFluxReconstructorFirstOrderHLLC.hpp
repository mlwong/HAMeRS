#ifndef CONVECTIVE_FLUX_RECONSTRUCTOR_FIRST_ORDER_HLLC_HPP
#define CONVECTIVE_FLUX_RECONSTRUCTOR_FIRST_ORDER_HLLC_HPP

#include "flow/convective_flux_reconstructors/ConvectiveFluxReconstructor.hpp"

class ConvectiveFluxReconstructorFirstOrderHLLC: public ConvectiveFluxReconstructor
{
    public:
        ConvectiveFluxReconstructorFirstOrderHLLC(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_eqn,
            const int& num_species,
            const boost::shared_ptr<FlowModel>& flow_model,
            const boost::shared_ptr<tbox::Database>& convective_flux_reconstructor_db);
        
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
            const boost::shared_ptr<pdat::FaceVariable<double> >& variable_convective_flux,
            const boost::shared_ptr<pdat::CellVariable<double> >& variable_source,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
    private:
        
};

#endif /* CONVECTIVE_FLUX_RECONSTRUCTOR_FIRST_ORDER_HLLC_HPP */
