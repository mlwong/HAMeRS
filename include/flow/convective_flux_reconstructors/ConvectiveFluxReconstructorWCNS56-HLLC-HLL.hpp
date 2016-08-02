#ifndef CONVECTIVE_FLUX_RECONSTRUCTOR_56_HLLC_HLL_HPP
#define CONVECTIVE_FLUX_RECONSTRUCTOR_56_HLLC_HLL_HPP

#include "flow/convective_flux_reconstructors/ConvectiveFluxReconstructor.hpp"

#include "Directions.hpp"

#include "boost/multi_array.hpp"

class ConvectiveFluxReconstructorWCNS56: public ConvectiveFluxReconstructor
{
    public:
        ConvectiveFluxReconstructorWCNS56(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_eqn,
            const int& num_species,
            const boost::shared_ptr<FlowModel>& flow_model,
            const boost::shared_ptr<tbox::Database>& convective_flux_reconstructor_db);
        
        virtual ~ConvectiveFluxReconstructorWCNS56() {}
        
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
        
        /*
         * Print all characteristics of the convective flux reconstruction class.
         */
        virtual void
        printClassData(std::ostream& os) const = 0;
        
        /*
         * Put the characteristics of the convective flux reconstruction class
         * into the restart database.
         */
        virtual void
        putToRestart(
            const boost::shared_ptr<tbox::Database>& restart_db) const = 0;
    
    protected:
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
         * Perform WENO interpolation.
         */
        virtual void
        performWENOInterpolation(
            std::vector<double>& U_minus,
            std::vector<double>& U_plus,
            const boost::multi_array<const double*, 2>& U_array,
            const hier::Index& cell_index_minus,
            const hier::Index& cell_index_plus,
            const DIRECTION& direction) = 0;
        
};

#endif /* CONVECTIVE_FLUX_RECONSTRUCTOR_56_HLLC_HLL_HPP */
