#ifndef NONCONSERVATIVE_DIFFUSIVE_FLUX_DIVERGENCE_OPERATOR_EIGHTH_ORDER_HPP
#define NONCONSERVATIVE_DIFFUSIVE_FLUX_DIVERGENCE_OPERATOR_EIGHTH_ORDER_HPP

#include "flow/nonconservative_diffusive_flux_divergence_operators/NonconservativeDiffusiveFluxDivergenceOperator.hpp"

class NonconservativeDiffusiveFluxDivergenceOperatorEighthOrder: public NonconservativeDiffusiveFluxDivergenceOperator
{
    public:
        NonconservativeDiffusiveFluxDivergenceOperatorEighthOrder(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_eqn,
            const HAMERS_SHARED_PTR<FlowModel>& flow_model,
            const HAMERS_SHARED_PTR<tbox::Database>& nonconservative_diffusive_flux_divergence_operator_db);
        
        ~NonconservativeDiffusiveFluxDivergenceOperatorEighthOrder() {}
        
        /*
         * Print all characteristics of the non-conservative diffusive flux divergence operator class.
         */
        void
        printClassData(std::ostream& os) const;
        
        /*
         * Put the characteristics of the non-conservative diffusive flux divergence operator class
         * into the restart database.
         */
        void
        putToRestart(
            const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const;
        
    private:
        /*
         * Compute the first derivatives in the x-direction.
         */
        void computeFirstDerivativesInX(
            std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > >& derivative_x,
            std::map<Real*, HAMERS_SHARED_PTR<pdat::CellData<Real> > >& derivative_x_computed,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > >& data_x,
            const std::vector<std::vector<int> >& data_component_idx_x,
            const hier::Patch& patch);
        
        /*
         * Compute the first derivatives in the y-direction.
         */
        void computeFirstDerivativesInY(
            std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > >& derivative_y,
            std::map<Real*, HAMERS_SHARED_PTR<pdat::CellData<Real> > >& derivative_y_computed,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > >& data_y,
            const std::vector<std::vector<int> >& data_component_idx_y,
            const hier::Patch& patch);
        
        /*
         * Compute the first derivatives in the z-direction.
         */
        void computeFirstDerivativesInZ(
            std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > >& derivative_z,
            std::map<Real*, HAMERS_SHARED_PTR<pdat::CellData<Real> > >& derivative_z_computed,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > >& data_z,
            const std::vector<std::vector<int> >& data_component_idx_z,
            const hier::Patch& patch);
        
        /*
         * Compute the second derivatives in the x-direction.
         */
        void computeSecondDerivativesInX(
            std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > >& derivative_x,
            std::map<Real*, HAMERS_SHARED_PTR<pdat::CellData<Real> > >& derivative_x_computed,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > >& data_x,
            const std::vector<std::vector<int> >& data_component_idx_x,
            const hier::Patch& patch);
        
        /*
         * Compute the second derivatives in the y-direction.
         */
        void computeSecondDerivativesInY(
            std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > >& derivative_y,
            std::map<Real*, HAMERS_SHARED_PTR<pdat::CellData<Real> > >& derivative_y_computed,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > >& data_y,
            const std::vector<std::vector<int> >& data_component_idx_y,
            const hier::Patch& patch);
        
        /*
         * Compute the second derivatives in the z-direction.
         */
        void computeSecondDerivativesInZ(
            std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > >& derivative_z,
            std::map<Real*, HAMERS_SHARED_PTR<pdat::CellData<Real> > >& derivative_z_computed,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > >& data_z,
            const std::vector<std::vector<int> >& data_component_idx_z,
            const hier::Patch& patch);
        
};

#endif /* NONCONSERVATIVE_DIFFUSIVE_FLUX_DIVERGENCE_OPERATOR_EIGHTH_ORDER_HPP */
