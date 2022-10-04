#ifndef DIFFUSIVE_FLUX_RECONSTRUCTOR_SIXTH_ORDER_HPP
#define DIFFUSIVE_FLUX_RECONSTRUCTOR_SIXTH_ORDER_HPP

#include "flow/diffusive_flux_reconstructors/DiffusiveFluxReconstructor.hpp"

class DiffusiveFluxReconstructorSixthOrder: public DiffusiveFluxReconstructor
{
    public:
        DiffusiveFluxReconstructorSixthOrder(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_eqn,
            const HAMERS_SHARED_PTR<FlowModel>& flow_model,
            const HAMERS_SHARED_PTR<tbox::Database>& diffusive_flux_reconstructor_db);
        
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
         * Compute the first derivatives in the x-direction.
         */
        void computeFirstDerivativesInX(
            std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& derivatives_x,
            std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > >& derivatives_x_computed,
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& data_x,
            const std::vector<int>& data_component_idx_x,
            const hier::Patch& patch,
            const bool allocate_scratch_derivatives_node);
            
        /*
         * Compute the first derivatives in the x-direction.
         */
        void computeFirstDerivativesInX(
            std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& derivatives_x,
            std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > >& derivatives_x_computed,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& data_x,
            const std::vector<std::vector<int> >& data_component_idx_x,
            const hier::Patch& patch,
            const bool allocate_scratch_derivatives_node);
        
        /*
         * Compute the first derivatives in the y-direction.
         */
        void computeFirstDerivativesInY(
            std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& derivatives_y,
            std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > >& derivatives_y_computed,
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& data_y,
            const std::vector<int>& data_component_idx_y,
            const hier::Patch& patch,
            const bool allocate_scratch_derivatives_node);
        
        /*
         * Compute the first derivatives in the y-direction.
         */
        void computeFirstDerivativesInY(
            std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& derivatives_y,
            std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > >& derivatives_y_computed,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& data_y,
            const std::vector<std::vector<int> >& data_component_idx_y,
            const hier::Patch& patch,
            const bool allocate_scratch_derivatives_node);
        
        /*
         * Compute the first derivatives in the z-direction.
         */
        void computeFirstDerivativesInZ(
            std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& derivatives_z,
            std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > >& derivatives_z_computed,
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& data_z,
            const std::vector<int>& data_component_idx_z,
            const hier::Patch& patch,
            const bool allocate_scratch_derivatives_node);
        
        /*
         * Compute the first derivatives in the z-direction.
         */
        void computeFirstDerivativesInZ(
            std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& derivatives_z,
            std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > >& derivatives_z_computed,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& data_z,
            const std::vector<std::vector<int> >& data_component_idx_z,
            const hier::Patch& patch,
            const bool allocate_scratch_derivatives_node);
        
        /*
         * Reconstruct the flux using flux at nodes in x-direction.
         */
        void reconstructFluxX(
            HAMERS_SHARED_PTR<pdat::SideData<double> >& diffusive_flux,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& diffusive_flux_node,
            const hier::Patch& patch,
            const double dt) const;
        
        /*
         * Reconstruct the flux using flux at nodes in y-direction.
         */
        void reconstructFluxY(
            HAMERS_SHARED_PTR<pdat::SideData<double> >& diffusive_flux,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& diffusive_flux_node,
            const hier::Patch& patch,
            const double dt) const;
        
        /*
         * Reconstruct the flux using flux at nodes in z-direction.
         */
        void reconstructFluxZ(
            HAMERS_SHARED_PTR<pdat::SideData<double> >& diffusive_flux,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& diffusive_flux_node,
            const hier::Patch& patch,
            const double dt) const;
        
        ///
        
        /*
         * Kernel to compute the first derivatives in the x-direction.
         */
        void computeFirstDerivativesInX(
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
         * Kernel to compute the first derivatives in the y-direction.
         */
        void computeFirstDerivativesInY(
            double* dudy,
            const double* const u,
            const hier::IntVector& num_ghosts_derivative_node,
            const hier::IntVector& num_ghosts_data_node,
            const hier::IntVector& ghostcell_dims_derivative_node,
            const hier::IntVector& ghostcell_dims_data_node,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims,
            const double& dx_1_in) const;
        
        /*
         * Kernel to compute the first derivatives in the z-direction.
         */
        void computeFirstDerivativesInZ(
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
         * Kernel to reconstruct the flux using flux at nodes in x-direction.
         */
        void reconstructFluxX(
            double* F_face_x,
            const double* const F_node_x,
            const hier::IntVector& num_ghosts_flux_node,
            const hier::IntVector& ghostcell_dims_flux_node,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims,
            const hier::IntVector& interior_dims,
            const double& dt) const;
        
        /*
         * Kernel to reconstruct the flux using flux at nodes in y-direction.
         */
        void reconstructFluxY(
            double* F_face_y,
            const double* const F_node_y,
            const hier::IntVector& num_ghosts_flux_node,
            const hier::IntVector& ghostcell_dims_flux_node,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims,
            const hier::IntVector& interior_dims,
            const double& dt) const;
        
        /*
         * Kernel to reconstruct the flux using flux at nodes in z-direction.
         */
        void reconstructFluxZ(
            double* F_face_z,
            const double* const F_node_z,
            const hier::IntVector& num_ghosts_flux_node,
            const hier::IntVector& ghostcell_dims_flux_node,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims,
            const hier::IntVector& interior_dims,
            const double& dt) const;
        
        /*
         * Numbers of ghost cells needed for the stencil operations.
         */
        hier::IntVector d_num_der_node_ghosts;
        hier::IntVector d_num_flux_reconstruct_ghosts;
        
        /*
         * Number of scratch data containers used.
         */
        int d_num_scratch_derivatives_node_used;
        
        /*
         * Scratch data containers.
         */
        std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > d_scratch_derivatives_node;
        
};

#endif /* DIFFUSIVE_FLUX_RECONSTRUCTOR_SIXTH_ORDER_HPP */
