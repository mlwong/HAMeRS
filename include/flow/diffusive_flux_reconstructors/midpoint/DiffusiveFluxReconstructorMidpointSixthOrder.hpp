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
        void computeFirstDerivativesInXAtMidpointX(
            std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& derivatives_x_midpoint_x,
            std::map<double*, HAMERS_SHARED_PTR<pdat::SideData<double> > >& derivatives_x_midpoint_x_computed,
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& data_x,
            const std::vector<int>& data_component_idx_x,
            const hier::Patch& patch,
            const bool allocate_scratch_derivatives_midpoint);
        
        /*
         * Compute the derivatives in x-direction at midpoints.
         */
        void computeFirstDerivativesInXAtMidpointX(
            std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > >& derivatives_x_midpoint_x,
            std::map<double*, HAMERS_SHARED_PTR<pdat::SideData<double> > >& derivatives_x_midpoint_x_computed,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& data_x,
            const std::vector<std::vector<int> >& data_component_idx_x,
            const hier::Patch& patch,
            const bool allocate_scratch_derivatives_midpoint);
        
        /*
         * Compute the derivatives in y-direction at midpoints.
         */
        void computeFirstDerivativesInYAtMidpointY(
            std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& derivatives_y_midpoint_y,
            std::map<double*, HAMERS_SHARED_PTR<pdat::SideData<double> > >& derivatives_y_midpoint_y_computed,
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& data_y,
            const std::vector<int>& data_component_idx_y,
            const hier::Patch& patch,
            const bool allocate_scratch_derivatives_midpoint);
        
        /*
         * Compute the derivatives in y-direction at midpoints.
         */
        void computeFirstDerivativesInYAtMidpointY(
            std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > >& derivatives_y_midpoint_y,
            std::map<double*, HAMERS_SHARED_PTR<pdat::SideData<double> > >& derivatives_y_midpoint_y_computed,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& data_y,
            const std::vector<std::vector<int> >& data_component_idx_y,
            const hier::Patch& patch,
            const bool allocate_scratch_derivatives_midpoint);
        
        /*
         * Compute the derivatives in z-direction at midpoints.
         */
        void computeFirstDerivativesInZAtMidpointZ(
            std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& derivatives_z_midpoint_z,
            std::map<double*, HAMERS_SHARED_PTR<pdat::SideData<double> > >& derivatives_z_midpoint_z_computed,
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& data_z,
            const std::vector<int>& data_component_idx_z,
            const hier::Patch& patch,
            const bool allocate_scratch_derivatives_midpoint);
        
        /*
         * Compute the derivatives in z-direction at midpoints.
         */
        void computeFirstDerivativesInZAtMidpointZ(
            std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > >& derivatives_z_midpoint_z,
            std::map<double*, HAMERS_SHARED_PTR<pdat::SideData<double> > >& derivatives_z_midpoint_z_computed,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& data_z,
            const std::vector<std::vector<int> >& data_component_idx_z,
            const hier::Patch& patch,
            const bool allocate_scratch_derivatives_midpoint);
        
        /*
         * Compute the derivatives in x-direction at nodes.
         */
        void computeFirstDerivativesInXAtNode(
            std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& derivatives_x_node,
            std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > >& derivatives_x_node_computed,
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& data_x,
            const std::vector<int>& data_component_idx_x,
            const hier::Patch& patch,
            const bool allocate_scratch_derivatives_node);
        
        /*
         * Compute the derivatives in x-direction at nodes.
         */
        void computeFirstDerivativesInXAtNode(
            std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& derivatives_x_node,
            std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > >& derivatives_x_node_computed,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& data_x,
            const std::vector<std::vector<int> >& data_component_idx_x,
            const hier::Patch& patch,
            const bool allocate_scratch_derivatives_node);
        
        /*
         * Compute the derivatives in y-direction at nodes.
         */
        void computeFirstDerivativesInYAtNode(
            std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& derivatives_y_node,
            std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > >& derivatives_y_node_computed,
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& data_y,
            const std::vector<int>& data_component_idx_y,
            const hier::Patch& patch,
            const bool allocate_scratch_derivatives_node);
        
        /*
         * Compute the derivatives in y-direction at nodes.
         */
        void computeFirstDerivativesInYAtNode(
            std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& derivatives_y_node,
            std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > >& derivatives_y_node_computed,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& data_y,
            const std::vector<std::vector<int> >& data_component_idx_y,
            const hier::Patch& patch,
            const bool allocate_scratch_derivatives_node);
        
        /*
         * Compute the derivatives in z-direction at nodes.
         */
        void computeFirstDerivativesInZAtNode(
            std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& derivatives_z_node,
            std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > >& derivatives_z_node_computed,
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& data_z,
            const std::vector<int>& data_component_idx_z,
            const hier::Patch& patch,
            const bool allocate_scratch_derivatives_node);
        
        /*
         * Compute the derivatives in z-direction at nodes.
         */
        void computeFirstDerivativesInZAtNode(
            std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& derivatives_z_node,
            std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > >& derivatives_z_node_computed,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& data_z,
            const std::vector<std::vector<int> >& data_component_idx_z,
            const hier::Patch& patch,
            const bool allocate_scratch_derivatives_node);
        
        /*
         * Interpolate the diffusivities from nodes to midpoints.
         */
        void interpolateDiffusivitiesFromNodeToMidpoint(
            std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& var_side_data_for_diffusivities,
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& var_cell_data_for_diffusivities,
            const std::vector<int>& var_cell_data_for_diffusivities_component_idx,
            const hier::Patch& patch,
            const bool allocate_scratch_diffusivities_midpoint);
        
        /*
         * Interpolate the derivatives from nodes to midpoints in x-direction.
         */
        void interpolateDerivativesFromNodeToMidpointX(
            std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& derivatives_midpoint_x,
            std::map<double*, HAMERS_SHARED_PTR<pdat::SideData<double> > >& derivatives_midpoint_x_computed,
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& derivative_node,
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& data,
            const std::vector<int>& data_component_idx,
            const hier::Patch& patch,
            const bool allocate_scratch_derivatives_midpoint);
        
        /*
         * Interpolate the derivatives from nodes to midpoints in x-direction.
         */
        void interpolateDerivativesFromNodeToMidpointX(
            std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > >& derivatives_midpoint_x,
            std::map<double*, HAMERS_SHARED_PTR<pdat::SideData<double> > >& derivatives_midpoint_x_computed,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& derivative_node,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& data,
            const std::vector<std::vector<int> >& data_component_idx,
            const hier::Patch& patch,
            const bool allocate_scratch_derivatives_midpoint);
        
        /*
         * Interpolate the derivatives from nodes to midpoints in y-direction.
         */
        void interpolateDerivativesFromNodeToMidpointY(
            std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& derivatives_midpoint_y,
            std::map<double*, HAMERS_SHARED_PTR<pdat::SideData<double> > >& derivatives_midpoint_y_computed,
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& derivative_node,
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& data,
            const std::vector<int>& data_component_idx,
            const hier::Patch& patch,
            const bool allocate_scratch_derivatives_midpoint);
        
        /*
         * Interpolate the derivatives from nodes to midpoints in y-direction.
         */
        void interpolateDerivativesFromNodeToMidpointY(
            std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > >& derivatives_midpoint_y,
            std::map<double*, HAMERS_SHARED_PTR<pdat::SideData<double> > >& derivatives_midpoint_y_computed,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& derivative_node,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& data,
            const std::vector<std::vector<int> >& data_component_idx,
            const hier::Patch& patch,
            const bool allocate_scratch_derivatives_midpoint);
        
        /*
         * Interpolate the derivatives from nodes to midpoints in z-direction.
         */
        void interpolateDerivativesFromNodeToMidpointZ(
            std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& derivatives_midpoint_z,
            std::map<double*, HAMERS_SHARED_PTR<pdat::SideData<double> > >& derivatives_midpoint_z_computed,
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& derivative_node,
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& data,
            const std::vector<int>& data_component_idx,
            const hier::Patch& patch,
            const bool allocate_scratch_derivatives_midpoint);
        
        /*
         * Interpolate the derivatives from nodes to midpoints in z-direction.
         */
        void interpolateDerivativesFromNodeToMidpointZ(
            std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > >& derivatives_midpoint_z,
            std::map<double*, HAMERS_SHARED_PTR<pdat::SideData<double> > >& derivatives_midpoint_z_computed,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& derivative_node,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& data,
            const std::vector<std::vector<int> >& data_component_idx,
            const hier::Patch& patch,
            const bool allocate_scratch_derivatives_midpoint);
        
        /*
         * Reconstruct the flux using flux at midpoints in x-direction.
         */
        void reconstructFluxX(
            HAMERS_SHARED_PTR<pdat::SideData<double> >& diffusive_flux,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& diffusive_flux_midpoint,
            const hier::Patch& patch,
            const double dt) const;
        
        /*
         * Reconstruct the flux using flux at midpoints in y-direction.
         */
        void reconstructFluxY(
            HAMERS_SHARED_PTR<pdat::SideData<double> >& diffusive_flux,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& diffusive_flux_midpoint,
            const hier::Patch& patch,
            const double dt) const;
        
        /*
         * Reconstruct the flux using flux at midpoints in z-direction.
         */
        void reconstructFluxZ(
            HAMERS_SHARED_PTR<pdat::SideData<double> >& diffusive_flux,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& diffusive_flux_midpoint,
            const hier::Patch& patch,
            const double dt) const;
        
        ///
        
        /*
         * Kernel to compute the derivatives in x-direction at midpoints.
         */
        void computeFirstDerivativesInXAtMidpointX(
            double* dudx,
            const double* const u,
            const hier::IntVector& num_ghosts_derivative_midpoint,
            const hier::IntVector& num_ghosts_data_node,
            const hier::IntVector& ghostcell_dims_derivative_midpoint,
            const hier::IntVector& ghostcell_dims_data_node,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims,
            const double& dx_0_inv) const;
        
        /*
         * Kernel to compute the derivatives in y-direction at midpoints.
         */
        void computeFirstDerivativesInYAtMidpointY(
            double* dudy,
            const double* const u,
            const hier::IntVector& num_ghosts_derivative_midpoint,
            const hier::IntVector& num_ghosts_data_node,
            const hier::IntVector& ghostcell_dims_derivative_midpoint,
            const hier::IntVector& ghostcell_dims_data_node,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims,
            const double& dx_1_inv) const;
        
        /*
         * Kernel to compute the derivatives in z-direction at midpoints.
         */
        void computeFirstDerivativesInZAtMidpointZ(
            double* dudz,
            const double* const u,
            const hier::IntVector& num_ghosts_derivative_midpoint,
            const hier::IntVector& num_ghosts_data_node,
            const hier::IntVector& ghostcell_dims_derivative_midpoint,
            const hier::IntVector& ghostcell_dims_data_node,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims,
            const double& dx_2_inv) const;
        
        /*
         * Kernel to compute the derivatives in x-direction at nodes.
         */
        void computeFirstDerivativesInXAtNode(
            double* dudx,
            const double* const u,
            const hier::IntVector& num_ghosts_derivative_node,
            const hier::IntVector& num_ghosts_data_node,
            const hier::IntVector& ghostcell_dims_derivative_node,
            const hier::IntVector& ghostcell_dims_data_node,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims,
            const double& dx_0_inv) const;
        
        /*
         * Kernel to compute the derivatives in y-direction at nodes.
         */
        void computeFirstDerivativesInYAtNode(
            double* dudy,
            const double* const u,
            const hier::IntVector& num_ghosts_derivative_node,
            const hier::IntVector& num_ghosts_data_node,
            const hier::IntVector& ghostcell_dims_derivative_node,
            const hier::IntVector& ghostcell_dims_data_node,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims,
            const double& dx_1_inv) const;
        
        /*
         * Kernel to compute the derivatives in z-direction at nodes.
         */
        void computeFirstDerivativesInZAtNode(
            double* dudz,
            const double* const u,
            const hier::IntVector& num_ghosts_derivative_node,
            const hier::IntVector& num_ghosts_data_node,
            const hier::IntVector& ghostcell_dims_derivative_node,
            const hier::IntVector& ghostcell_dims_data_node,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims,
            const double& dx_2_inv) const;
        
        /*
         * Kernel to interpolate the data from nodes to midpoints in x-direction.
         */
        void interpolateDataFromNodeToMidpointX(
            double* u_midpoint_x,
            const double* const u_node,
            const hier::IntVector& num_ghosts_data_midpoint,
            const hier::IntVector& num_ghosts_data_node,
            const hier::IntVector& ghostcell_dims_data_midpoint,
            const hier::IntVector& ghostcell_dims_data_node,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims) const;
        
        /*
         * Kernel to interpolate the data from nodes to midpoints in y-direction.
         */
        void interpolateDataFromNodeToMidpointY(
            double* u_midpoint_y,
            const double* const u_node,
            const hier::IntVector& num_ghosts_data_midpoint,
            const hier::IntVector& num_ghosts_data_node,
            const hier::IntVector& ghostcell_dims_data_midpoint,
            const hier::IntVector& ghostcell_dims_data_node,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims) const;
        
        /*
         * Kernel to interpolate the data from nodes to midpoints in z-direction.
         */
        void interpolateDataFromNodeToMidpointZ(
            double* u_midpoint_z,
            const double* const u_node,
            const hier::IntVector& num_ghosts_derivative_midpoint,
            const hier::IntVector& num_ghosts_derivative_node,
            const hier::IntVector& ghostcell_dims_data_midpoint,
            const hier::IntVector& ghostcell_dims_data_node,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims) const;
        
        /*
         * Kernel to reconstruct the flux using flux at midpoints in x-direction.
         */
        void
        reconstructFluxX(
            double* F_face_x,
            const double* const F_midpoint_x,
            const hier::IntVector& num_ghosts_flux_midpoint,
            const hier::IntVector& ghostcell_dims_flux_midpoint,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims,
            const hier::IntVector& interior_dims,
            const double& dt) const;
        
        /*
         * Kernel to reconstruct the flux using flux at midpoints in y-direction.
         */
        void
        reconstructFluxY(
            double* F_face_y,
            const double* const F_midpoint_y,
            const hier::IntVector& num_ghosts_flux_midpoint,
            const hier::IntVector& ghostcell_dims_flux_midpoint,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims,
            const hier::IntVector& interior_dims,
            const double& dt) const;
        
        /*
         * Kernel to reconstruct the flux using flux at midpoints in z-direction.
         */
        void
        reconstructFluxZ(
            double* F_face_z,
            const double* const F_midpoint_z,
            const hier::IntVector& num_ghosts_flux_midpoint,
            const hier::IntVector& ghostcell_dims_flux_midpoint,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims,
            const hier::IntVector& interior_dims,
            const double& dt) const;
        
        /*
         * Numbers of ghost cells needed for the stencil operations.
         */
        hier::IntVector d_num_der_midpoint_ghosts;
        hier::IntVector d_num_der_node_ghosts;
        hier::IntVector d_num_interp_midpoint_ghosts;
        hier::IntVector d_num_flux_reconstruct_ghosts;
        
        /*
         * Numbers of scratch data containers used.
         */
        int d_num_scratch_derivatives_node_used;
        int d_num_scratch_derivatives_midpoint_x_used;
        int d_num_scratch_derivatives_midpoint_y_used;
        int d_num_scratch_derivatives_midpoint_z_used;
        int d_num_scratch_diffusivities_midpoint_used;
        
        /*
         * Scratch data containers.
         */
        std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > d_scratch_derivatives_node;
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > d_scratch_derivatives_midpoint_x;
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > d_scratch_derivatives_midpoint_y;
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > d_scratch_derivatives_midpoint_z;
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > d_scratch_diffusivities_midpoint;
        
};

#endif /* DIFFUSIVE_FLUX_RECONSTRUCTOR_MIDPOINT_SIXTH_ORDER_HPP */
