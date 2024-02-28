#ifndef DIFFUSIVE_FLUX_RECONSTRUCTOR_MIDPOINT_HPP
#define DIFFUSIVE_FLUX_RECONSTRUCTOR_MIDPOINT_HPP

#include "flow/diffusive_flux_reconstructors/DiffusiveFluxReconstructor.hpp"

class DiffusiveFluxReconstructorMidpoint: public DiffusiveFluxReconstructor
{
    public:
        DiffusiveFluxReconstructorMidpoint(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_eqn,
            const HAMERS_SHARED_PTR<FlowModel>& flow_model,
            const HAMERS_SHARED_PTR<tbox::Database>& diffusive_flux_reconstructor_db);
        
        virtual ~DiffusiveFluxReconstructorMidpoint() {}
        
        /*
         * Print all characteristics of the diffusive flux reconstruction class.
         */
        virtual void
        printClassData(std::ostream& os) const = 0;
        
        /*
         * Put the characteristics of the diffusive flux reconstruction class
         * into the restart database.
         */
        virtual void
        putToRestart(
            const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const = 0;
        
        /*
         * Compute the diffusive flux on a patch.
         */
        void computeDiffusiveFluxOnPatch(
            hier::Patch& patch,
            const HAMERS_SHARED_PTR<hier::CoarseFineBoundary> coarse_fine_bdry,
            const HAMERS_SHARED_PTR<pdat::SideVariable<Real> >& variable_diffusive_flux,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double time,
            const double dt,
            const int RK_step_number);
        
    private:
        /*
         * Compute the derivatives in x-direction at midpoints.
         */
        void computeFirstDerivativesInXAtMidpointX(
            std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& derivatives_x_midpoint_x,
            std::map<Real*, HAMERS_SHARED_PTR<pdat::SideData<Real> > >& derivatives_x_midpoint_x_computed,
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& data_x,
            const std::vector<int>& data_component_idx_x,
            const hier::Patch& patch,
            const bool allocate_scratch_derivatives_midpoint);
        
        /*
         * Compute the derivatives in x-direction at midpoints.
         */
        void computeFirstDerivativesInXAtMidpointX(
            std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > >& derivatives_x_midpoint_x,
            std::map<Real*, HAMERS_SHARED_PTR<pdat::SideData<Real> > >& derivatives_x_midpoint_x_computed,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > >& data_x,
            const std::vector<std::vector<int> >& data_component_idx_x,
            const hier::Patch& patch,
            const bool allocate_scratch_derivatives_midpoint);
        
        /*
         * Compute the derivatives in y-direction at midpoints.
         */
        void computeFirstDerivativesInYAtMidpointY(
            std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& derivatives_y_midpoint_y,
            std::map<Real*, HAMERS_SHARED_PTR<pdat::SideData<Real> > >& derivatives_y_midpoint_y_computed,
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& data_y,
            const std::vector<int>& data_component_idx_y,
            const hier::Patch& patch,
            const bool allocate_scratch_derivatives_midpoint);
        
        /*
         * Compute the derivatives in y-direction at midpoints.
         */
        void computeFirstDerivativesInYAtMidpointY(
            std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > >& derivatives_y_midpoint_y,
            std::map<Real*, HAMERS_SHARED_PTR<pdat::SideData<Real> > >& derivatives_y_midpoint_y_computed,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > >& data_y,
            const std::vector<std::vector<int> >& data_component_idx_y,
            const hier::Patch& patch,
            const bool allocate_scratch_derivatives_midpoint);
        
        /*
         * Compute the derivatives in z-direction at midpoints.
         */
        void computeFirstDerivativesInZAtMidpointZ(
            std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& derivatives_z_midpoint_z,
            std::map<Real*, HAMERS_SHARED_PTR<pdat::SideData<Real> > >& derivatives_z_midpoint_z_computed,
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& data_z,
            const std::vector<int>& data_component_idx_z,
            const hier::Patch& patch,
            const bool allocate_scratch_derivatives_midpoint);
        
        /*
         * Compute the derivatives in z-direction at midpoints.
         */
        void computeFirstDerivativesInZAtMidpointZ(
            std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > >& derivatives_z_midpoint_z,
            std::map<Real*, HAMERS_SHARED_PTR<pdat::SideData<Real> > >& derivatives_z_midpoint_z_computed,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > >& data_z,
            const std::vector<std::vector<int> >& data_component_idx_z,
            const hier::Patch& patch,
            const bool allocate_scratch_derivatives_midpoint);
        
        /*
         * Compute the derivatives in x-direction at nodes.
         */
        void computeFirstDerivativesInXAtNode(
            std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& derivatives_x_node,
            std::map<Real*, HAMERS_SHARED_PTR<pdat::CellData<Real> > >& derivatives_x_node_computed,
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& data_x,
            const std::vector<int>& data_component_idx_x,
            const hier::Patch& patch,
            const bool allocate_scratch_derivatives_node);
        
        /*
         * Compute the derivatives in x-direction at nodes.
         */
        void computeFirstDerivativesInXAtNode(
            std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > >& derivatives_x_node,
            std::map<Real*, HAMERS_SHARED_PTR<pdat::CellData<Real> > >& derivatives_x_node_computed,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > >& data_x,
            const std::vector<std::vector<int> >& data_component_idx_x,
            const hier::Patch& patch,
            const bool allocate_scratch_derivatives_node);
        
        /*
         * Compute the derivatives in y-direction at nodes.
         */
        void computeFirstDerivativesInYAtNode(
            std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& derivatives_y_node,
            std::map<Real*, HAMERS_SHARED_PTR<pdat::CellData<Real> > >& derivatives_y_node_computed,
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& data_y,
            const std::vector<int>& data_component_idx_y,
            const hier::Patch& patch,
            const bool allocate_scratch_derivatives_node);
        
        /*
         * Compute the derivatives in y-direction at nodes.
         */
        void computeFirstDerivativesInYAtNode(
            std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > >& derivatives_y_node,
            std::map<Real*, HAMERS_SHARED_PTR<pdat::CellData<Real> > >& derivatives_y_node_computed,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > >& data_y,
            const std::vector<std::vector<int> >& data_component_idx_y,
            const hier::Patch& patch,
            const bool allocate_scratch_derivatives_node);
        
        /*
         * Compute the derivatives in z-direction at nodes.
         */
        void computeFirstDerivativesInZAtNode(
            std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& derivatives_z_node,
            std::map<Real*, HAMERS_SHARED_PTR<pdat::CellData<Real> > >& derivatives_z_node_computed,
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& data_z,
            const std::vector<int>& data_component_idx_z,
            const hier::Patch& patch,
            const bool allocate_scratch_derivatives_node);
        
        /*
         * Compute the derivatives in z-direction at nodes.
         */
        void computeFirstDerivativesInZAtNode(
            std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > >& derivatives_z_node,
            std::map<Real*, HAMERS_SHARED_PTR<pdat::CellData<Real> > >& derivatives_z_node_computed,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > >& data_z,
            const std::vector<std::vector<int> >& data_component_idx_z,
            const hier::Patch& patch,
            const bool allocate_scratch_derivatives_node);
        
        /*
         * Interpolate the diffusivities from nodes to midpoints.
         */
        void interpolateDiffusivitiesFromNodeToMidpoint(
            std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& var_side_data_for_diffusivities,
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& var_cell_data_for_diffusivities,
            const std::vector<int>& var_cell_data_for_diffusivities_component_idx,
            const hier::Patch& patch,
            const bool allocate_scratch_diffusivities_midpoint);
        
        /*
         * Interpolate the derivatives from nodes to midpoints in x-direction.
         */
        void interpolateDerivativesFromNodeToMidpointX(
            std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& derivatives_midpoint_x,
            std::map<Real*, HAMERS_SHARED_PTR<pdat::SideData<Real> > >& derivatives_midpoint_x_computed,
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& derivative_node,
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& data,
            const std::vector<int>& data_component_idx,
            const hier::Patch& patch,
            const bool allocate_scratch_derivatives_midpoint);
        
        /*
         * Interpolate the derivatives from nodes to midpoints in x-direction.
         */
        void interpolateDerivativesFromNodeToMidpointX(
            std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > >& derivatives_midpoint_x,
            std::map<Real*, HAMERS_SHARED_PTR<pdat::SideData<Real> > >& derivatives_midpoint_x_computed,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > >& derivative_node,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > >& data,
            const std::vector<std::vector<int> >& data_component_idx,
            const hier::Patch& patch,
            const bool allocate_scratch_derivatives_midpoint);
        
        /*
         * Interpolate the derivatives from nodes to midpoints in y-direction.
         */
        void interpolateDerivativesFromNodeToMidpointY(
            std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& derivatives_midpoint_y,
            std::map<Real*, HAMERS_SHARED_PTR<pdat::SideData<Real> > >& derivatives_midpoint_y_computed,
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& derivative_node,
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& data,
            const std::vector<int>& data_component_idx,
            const hier::Patch& patch,
            const bool allocate_scratch_derivatives_midpoint);
        
        /*
         * Interpolate the derivatives from nodes to midpoints in y-direction.
         */
        void interpolateDerivativesFromNodeToMidpointY(
            std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > >& derivatives_midpoint_y,
            std::map<Real*, HAMERS_SHARED_PTR<pdat::SideData<Real> > >& derivatives_midpoint_y_computed,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > >& derivative_node,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > >& data,
            const std::vector<std::vector<int> >& data_component_idx,
            const hier::Patch& patch,
            const bool allocate_scratch_derivatives_midpoint);
        
        /*
         * Interpolate the derivatives from nodes to midpoints in z-direction.
         */
        void interpolateDerivativesFromNodeToMidpointZ(
            std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& derivatives_midpoint_z,
            std::map<Real*, HAMERS_SHARED_PTR<pdat::SideData<Real> > >& derivatives_midpoint_z_computed,
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& derivative_node,
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& data,
            const std::vector<int>& data_component_idx,
            const hier::Patch& patch,
            const bool allocate_scratch_derivatives_midpoint);
        
        /*
         * Interpolate the derivatives from nodes to midpoints in z-direction.
         */
        void interpolateDerivativesFromNodeToMidpointZ(
            std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > >& derivatives_midpoint_z,
            std::map<Real*, HAMERS_SHARED_PTR<pdat::SideData<Real> > >& derivatives_midpoint_z_computed,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > >& derivative_node,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > >& data,
            const std::vector<std::vector<int> >& data_component_idx,
            const hier::Patch& patch,
            const bool allocate_scratch_derivatives_midpoint);
        
        /*
         * Reconstruct the flux using flux at midpoints in x-direction.
         */
        void reconstructFluxX(
            HAMERS_SHARED_PTR<pdat::SideData<Real> >& diffusive_flux,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& diffusive_flux_midpoint,
            const hier::Patch& patch,
            const double dt) const;
        
        /*
         * Reconstruct the flux using flux at midpoints in y-direction.
         */
        void reconstructFluxY(
            HAMERS_SHARED_PTR<pdat::SideData<Real> >& diffusive_flux,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& diffusive_flux_midpoint,
            const hier::Patch& patch,
            const double dt) const;
        
        /*
         * Reconstruct the flux using flux at midpoints in z-direction.
         */
        void reconstructFluxZ(
            HAMERS_SHARED_PTR<pdat::SideData<Real> >& diffusive_flux,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& diffusive_flux_midpoint,
            const hier::Patch& patch,
            const double dt) const;
        
    protected:
        /*
         * Kernel to compute the derivatives in x-direction at midpoints.
         */
        virtual void computeFirstDerivativesInXAtMidpointX(
            Real* dudx,
            const Real* const u,
            const hier::IntVector& num_ghosts_derivative_midpoint,
            const hier::IntVector& num_ghosts_data_node,
            const hier::IntVector& ghostcell_dims_derivative_midpoint,
            const hier::IntVector& ghostcell_dims_data_node,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims,
            const Real& dx_0_inv) const = 0;
        
        /*
         * Kernel to compute the derivatives in y-direction at midpoints.
         */
        virtual void computeFirstDerivativesInYAtMidpointY(
            Real* dudy,
            const Real* const u,
            const hier::IntVector& num_ghosts_derivative_midpoint,
            const hier::IntVector& num_ghosts_data_node,
            const hier::IntVector& ghostcell_dims_derivative_midpoint,
            const hier::IntVector& ghostcell_dims_data_node,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims,
            const Real& dx_1_inv) const = 0;
        
        /*
         * Kernel to compute the derivatives in z-direction at midpoints.
         */
        virtual void computeFirstDerivativesInZAtMidpointZ(
            Real* dudz,
            const Real* const u,
            const hier::IntVector& num_ghosts_derivative_midpoint,
            const hier::IntVector& num_ghosts_data_node,
            const hier::IntVector& ghostcell_dims_derivative_midpoint,
            const hier::IntVector& ghostcell_dims_data_node,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims,
            const Real& dx_2_inv) const = 0;
        
        /*
         * Kernel to compute the derivatives in x-direction at nodes.
         */
        virtual void computeFirstDerivativesInXAtNode(
            Real* dudx,
            const Real* const u,
            const hier::IntVector& num_ghosts_derivative_node,
            const hier::IntVector& num_ghosts_data_node,
            const hier::IntVector& ghostcell_dims_derivative_node,
            const hier::IntVector& ghostcell_dims_data_node,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims,
            const Real& dx_0_inv) const = 0;
        
        /*
         * Kernel to compute the derivatives in y-direction at nodes.
         */
        virtual void computeFirstDerivativesInYAtNode(
            Real* dudy,
            const Real* const u,
            const hier::IntVector& num_ghosts_derivative_node,
            const hier::IntVector& num_ghosts_data_node,
            const hier::IntVector& ghostcell_dims_derivative_node,
            const hier::IntVector& ghostcell_dims_data_node,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims,
            const Real& dx_1_inv) const = 0;
        
        /*
         * Kernel to compute the derivatives in z-direction at nodes.
         */
        virtual void computeFirstDerivativesInZAtNode(
            Real* dudz,
            const Real* const u,
            const hier::IntVector& num_ghosts_derivative_node,
            const hier::IntVector& num_ghosts_data_node,
            const hier::IntVector& ghostcell_dims_derivative_node,
            const hier::IntVector& ghostcell_dims_data_node,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims,
            const Real& dx_2_inv) const = 0;
        
        /*
         * Kernel to interpolate the data from nodes to midpoints in x-direction.
         */
        virtual void interpolateDataFromNodeToMidpointX(
            Real* u_midpoint_x,
            const Real* const u_node,
            const hier::IntVector& num_ghosts_data_midpoint,
            const hier::IntVector& num_ghosts_data_node,
            const hier::IntVector& ghostcell_dims_data_midpoint,
            const hier::IntVector& ghostcell_dims_data_node,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims) const = 0;
        
        /*
         * Kernel to interpolate the data from nodes to midpoints in y-direction.
         */
        virtual void interpolateDataFromNodeToMidpointY(
            Real* u_midpoint_y,
            const Real* const u_node,
            const hier::IntVector& num_ghosts_data_midpoint,
            const hier::IntVector& num_ghosts_data_node,
            const hier::IntVector& ghostcell_dims_data_midpoint,
            const hier::IntVector& ghostcell_dims_data_node,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims) const = 0;
        
        /*
         * Kernel to interpolate the data from nodes to midpoints in z-direction.
         */
        virtual void interpolateDataFromNodeToMidpointZ(
            Real* u_midpoint_z,
            const Real* const u_node,
            const hier::IntVector& num_ghosts_derivative_midpoint,
            const hier::IntVector& num_ghosts_derivative_node,
            const hier::IntVector& ghostcell_dims_data_midpoint,
            const hier::IntVector& ghostcell_dims_data_node,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims) const = 0;
        
        /*
         * Kernel to reconstruct the flux using flux at midpoints in x-direction.
         */
        virtual void
        reconstructFluxX(
            Real* F_face_x,
            const Real* const F_midpoint_x,
            const hier::IntVector& num_ghosts_flux_midpoint,
            const hier::IntVector& ghostcell_dims_flux_midpoint,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims,
            const hier::IntVector& interior_dims,
            const double& dt) const = 0;
        
        /*
         * Kernel to reconstruct the flux using flux at midpoints in y-direction.
         */
        virtual void
        reconstructFluxY(
            Real* F_face_y,
            const Real* const F_midpoint_y,
            const hier::IntVector& num_ghosts_flux_midpoint,
            const hier::IntVector& ghostcell_dims_flux_midpoint,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims,
            const hier::IntVector& interior_dims,
            const double& dt) const = 0;
        
        /*
         * Kernel to reconstruct the flux using flux at midpoints in z-direction.
         */
        virtual void
        reconstructFluxZ(
            Real* F_face_z,
            const Real* const F_midpoint_z,
            const hier::IntVector& num_ghosts_flux_midpoint,
            const hier::IntVector& ghostcell_dims_flux_midpoint,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims,
            const hier::IntVector& interior_dims,
            const double& dt) const = 0;
        
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
        std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > d_scratch_derivatives_node;
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > d_scratch_derivatives_midpoint_x;
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > d_scratch_derivatives_midpoint_y;
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > d_scratch_derivatives_midpoint_z;
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > d_scratch_diffusivities_midpoint;
        
};

#endif /* DIFFUSIVE_FLUX_RECONSTRUCTOR_MIDPOINT_HPP */
