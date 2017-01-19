#ifndef DIFFUSIVE_FLUX_RECONSTRUCTOR_SIXTH_ORDER_HPP
#define DIFFUSIVE_FLUX_RECONSTRUCTOR_SIXTH_ORDER_HPP

#include "flow/diffusive_flux_reconstructors/DiffusiveFluxReconstructor.hpp"

class DiffusiveFluxReconstructorSixthOrder: public DiffusiveFluxReconstructor
{
    public:
        DiffusiveFluxReconstructorSixthOrder(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_eqn,
            const int& num_species,
            const boost::shared_ptr<FlowModel>& flow_model,
            const boost::shared_ptr<tbox::Database>& diffusive_flux_reconstructor_db);
        
        ~DiffusiveFluxReconstructorSixthOrder() {}
        
        /*
         * Print all characteristics of the diffusive flux reconstruction class.
         */
        void
        printClassData(std::ostream& os) const;
        
        /*
         * Put the characteristics of the diffusive flux reconstruction class
         * into the restart database.
         */
        void
        putToRestart(
            const boost::shared_ptr<tbox::Database>& restart_db) const;
        
        /*
         * Compute the diffusive fluxes.
         */
        void computeDiffusiveFluxes(
            hier::Patch& patch,
            const double time,
            const double dt,
            const int RK_step_number,
            const boost::shared_ptr<pdat::FaceVariable<double> >& variable_diffusive_flux,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
    private:
        /*
         * Compute the derivatives in the x-direction for diffusive flux.
         */
        void computeDerivativesInXForDiffusiveFlux(
            hier::Patch& patch,
            std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > >& derivative_x,
            std::map<double*, boost::shared_ptr<pdat::CellData<double> > >& derivative_x_computed,
            const std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > >& derivative_var_data_x,
            const std::vector<std::vector<int> >& derivative_var_component_idx_x);
        
        /*
         * Compute the derivatives in the y-direction for diffusive flux.
         */
        void computeDerivativesInYForDiffusiveFlux(
            hier::Patch& patch,
            std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > >& derivative_y,
            std::map<double*, boost::shared_ptr<pdat::CellData<double> > >& derivative_y_computed,
            const std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > >& derivative_var_data_y,
            const std::vector<std::vector<int> >& derivative_var_component_idx_y);
        
        /*
         * Compute the derivatives in the z-direction for diffusive flux.
         */
        void computeDerivativesInZForDiffusiveFlux(
            hier::Patch& patch,
            std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > >& derivative_z,
            std::map<double*, boost::shared_ptr<pdat::CellData<double> > >& derivative_z_computed,
            const std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > >& derivative_var_data_z,
            const std::vector<std::vector<int> >& derivative_var_component_idx_z);
        
};

#endif /* DIFFUSIVE_FLUX_RECONSTRUCTOR_SIXTH_ORDER_HPP */
