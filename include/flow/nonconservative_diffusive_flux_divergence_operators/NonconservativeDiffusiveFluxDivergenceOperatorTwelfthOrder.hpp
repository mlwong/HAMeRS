#ifndef NONCONSERVATIVE_DIFFUSIVE_FLUX_DIVERGENCE_OPERATOR_TWELFTH_ORDER_HPP
#define NONCONSERVATIVE_DIFFUSIVE_FLUX_DIVERGENCE_OPERATOR_TWELFTH_ORDER_HPP

#include "flow/nonconservative_diffusive_flux_divergence_operators/NonconservativeDiffusiveFluxDivergenceOperator.hpp"

class NonconservativeDiffusiveFluxDivergenceOperatorTwelfthOrder: public NonconservativeDiffusiveFluxDivergenceOperator
{
    public:
        NonconservativeDiffusiveFluxDivergenceOperatorTwelfthOrder(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_eqn,
            const HAMERS_SHARED_PTR<FlowModel>& flow_model,
            const HAMERS_SHARED_PTR<tbox::Database>& nonconservative_diffusive_flux_divergence_operator_db);
        
        ~NonconservativeDiffusiveFluxDivergenceOperatorTwelfthOrder() {}
        
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
            hier::Patch& patch,
            std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& derivative_x,
            std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > >& derivative_x_computed,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& data_x,
            const std::vector<std::vector<int> >& data_component_idx_x);
        
        /*
         * Compute the first derivatives in the y-direction.
         */
        void computeFirstDerivativesInY(
            hier::Patch& patch,
            std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& derivative_y,
            std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > >& derivative_y_computed,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& data_y,
            const std::vector<std::vector<int> >& data_component_idx_y);
        
        /*
         * Compute the first derivatives in the z-direction.
         */
        void computeFirstDerivativesInZ(
            hier::Patch& patch,
            std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& derivative_z,
            std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > >& derivative_z_computed,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& data_z,
            const std::vector<std::vector<int> >& data_component_idx_z);
        
        /*
         * Compute the second derivatives in the x-direction.
         */
        void computeSecondDerivativesInX(
            hier::Patch& patch,
            std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& derivative_x,
            std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > >& derivative_x_computed,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& data_x,
            const std::vector<std::vector<int> >& data_component_idx_x);
        
        /*
         * Compute the second derivatives in the y-direction.
         */
        void computeSecondDerivativesInY(
            hier::Patch& patch,
            std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& derivative_y,
            std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > >& derivative_y_computed,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& data_y,
            const std::vector<std::vector<int> >& data_component_idx_y);
        
        /*
         * Compute the second derivatives in the z-direction.
         */
        void computeSecondDerivativesInZ(
            hier::Patch& patch,
            std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& derivative_z,
            std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > >& derivative_z_computed,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& data_z,
            const std::vector<std::vector<int> >& data_component_idx_z);
        
};

#endif /* NONCONSERVATIVE_DIFFUSIVE_FLUX_DIVERGENCE_OPERATOR_TWELFTH_ORDER_HPP */
