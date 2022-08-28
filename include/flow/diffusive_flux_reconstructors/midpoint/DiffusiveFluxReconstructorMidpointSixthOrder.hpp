#ifndef DIFFUSIVE_FLUX_RECONSTRUCTOR_MIDPOINT_SIXTH_ORDER_HPP
#define DIFFUSIVE_FLUX_RECONSTRUCTOR_MIDPOINT_SIXTH_ORDER_HPP

#include "flow/diffusive_flux_reconstructors/DiffusiveFluxReconstructor.hpp"

class DiffusiveFluxReconstructorMidpointSixthOrder: public DiffusiveFluxReconstructor
{
    public:
        DiffusiveFluxReconstructorMidpointSixthOrder(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_eqn,
            const HAMERS_SHARED_PTR<FlowModel>& flow_model,
            const HAMERS_SHARED_PTR<tbox::Database>& diffusive_flux_reconstructor_db);
        
        ~DiffusiveFluxReconstructorMidpointSixthOrder() {}
        
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
            const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const;
        
        /*
         * Compute the diffusive flux on a patch.
         */
        void computeDiffusiveFluxOnPatch(
            hier::Patch& patch,
            const HAMERS_SHARED_PTR<pdat::SideVariable<double> >& variable_diffusive_flux,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double time,
            const double dt,
            const int RK_step_number);
        
    private:
        /*
         * Compute the derivatives in x-direction at midpoints.
         */
        void computeFirstDerivativesInXAtMidpoint(
            hier::Patch& patch,
            std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > >& derivative_x_midpoint,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& data_x,
            const std::vector<std::vector<int> >& data_component_idx_x);
        
        /*
         * Compute the derivatives in y-direction at midpoints.
         */
        void computeFirstDerivativesInYAtMidpoint(
            hier::Patch& patch,
            std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > >& derivative_y_midpoint,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& data_y,
            const std::vector<std::vector<int> >& data_component_idx_y);
        
        /*
         * Compute the derivatives in z-direction at midpoints.
         */
        void computeFirstDerivativesInZAtMidpoint(
            hier::Patch& patch,
            std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > >& derivative_z_midpoint,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& data_z,
            const std::vector<std::vector<int> >& data_component_idx_z);
        
        /*
         * Compute the derivatives in x-direction at nodes.
         */
        void computeFirstDerivativesInXAtNode(
            hier::Patch& patch,
            std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& derivative_x_node,
            std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > >& derivative_x_node_computed,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& var_data_x,
            const std::vector<std::vector<int> >& var_component_idx_x);
        
        /*
         * Compute the derivatives in y-direction at nodes.
         */
        void computeFirstDerivativesInYAtNode(
            hier::Patch& patch,
            std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& derivative_y_node,
            std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > >& derivative_y_node_computed,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& var_data_y,
            const std::vector<std::vector<int> >& var_component_idx_y);
        
        /*
         * Compute the derivatives in z-direction at nodes.
         */
        void computeFirstDerivativesInZAtNode(
            hier::Patch& patch,
            std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& derivative_z_node,
            std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > >& derivative_z_node_computed,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& var_data_z,
            const std::vector<std::vector<int> >& var_component_idx_z);
        
        /*
         * Interpolate the variables from nodes to midpoints in x-direction.
         */
        void interpolateFromNodeToMidpointX(
            hier::Patch& patch,
            std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > >& var_midpoint,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& var_node);
        
        /*
         * Interpolate the variables from nodes to midpoints in y-direction.
         */
        void interpolateFromNodeToMidpointY(
            hier::Patch& patch,
            std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > >& var_midpoint,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& var_node);
        
        /*
         * Interpolate the variables from nodes to midpoints in z-direction.
         */
        void interpolateFromNodeToMidpointZ(
            hier::Patch& patch,
            std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > >& var_midpoint,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& var_node);
        
};

#endif /* DIFFUSIVE_FLUX_RECONSTRUCTOR_MIDPOINT_SIXTH_ORDER_HPP */
