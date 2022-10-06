#ifndef DIFFUSIVE_FLUX_RECONSTRUCTOR_NODE_SIXTH_ORDER_HPP
#define DIFFUSIVE_FLUX_RECONSTRUCTOR_NODE_SIXTH_ORDER_HPP

#include "flow/diffusive_flux_reconstructors/node/DiffusiveFluxReconstructorNode.hpp"

class DiffusiveFluxReconstructorNodeSixthOrder: public DiffusiveFluxReconstructorNode
{
    public:
        DiffusiveFluxReconstructorNodeSixthOrder(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_eqn,
            const HAMERS_SHARED_PTR<FlowModel>& flow_model,
            const HAMERS_SHARED_PTR<tbox::Database>& diffusive_flux_reconstructor_db);
        
        ~DiffusiveFluxReconstructorNodeSixthOrder() {}
        
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
        
    private:
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
        
};

#endif /* DIFFUSIVE_FLUX_RECONSTRUCTOR_NODE_SIXTH_ORDER_HPP */
