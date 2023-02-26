#ifndef DIFFUSIVE_FLUX_RECONSTRUCTOR_MIDPOINT_SECOND_ORDER_HPP
#define DIFFUSIVE_FLUX_RECONSTRUCTOR_MIDPOINT_SECOND_ORDER_HPP

#include "flow/diffusive_flux_reconstructors/midpoint/DiffusiveFluxReconstructorMidpoint.hpp"

class DiffusiveFluxReconstructorMidpointSecondOrder: public DiffusiveFluxReconstructorMidpoint
{
    public:
        DiffusiveFluxReconstructorMidpointSecondOrder(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_eqn,
            const HAMERS_SHARED_PTR<FlowModel>& flow_model,
            const HAMERS_SHARED_PTR<tbox::Database>& diffusive_flux_reconstructor_db);
        
        ~DiffusiveFluxReconstructorMidpointSecondOrder() {}
        
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
         * Kernel to compute the derivatives in x-direction at midpoints.
         */
        void computeFirstDerivativesInXAtMidpointX(
            Real* dudx,
            const Real* const u,
            const hier::IntVector& num_ghosts_derivative_midpoint,
            const hier::IntVector& num_ghosts_data_node,
            const hier::IntVector& ghostcell_dims_derivative_midpoint,
            const hier::IntVector& ghostcell_dims_data_node,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims,
            const Real& dx_0_inv) const;
        
        /*
         * Kernel to compute the derivatives in y-direction at midpoints.
         */
        void computeFirstDerivativesInYAtMidpointY(
            Real* dudy,
            const Real* const u,
            const hier::IntVector& num_ghosts_derivative_midpoint,
            const hier::IntVector& num_ghosts_data_node,
            const hier::IntVector& ghostcell_dims_derivative_midpoint,
            const hier::IntVector& ghostcell_dims_data_node,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims,
            const Real& dx_1_inv) const;
        
        /*
         * Kernel to compute the derivatives in z-direction at midpoints.
         */
        void computeFirstDerivativesInZAtMidpointZ(
            Real* dudz,
            const Real* const u,
            const hier::IntVector& num_ghosts_derivative_midpoint,
            const hier::IntVector& num_ghosts_data_node,
            const hier::IntVector& ghostcell_dims_derivative_midpoint,
            const hier::IntVector& ghostcell_dims_data_node,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims,
            const Real& dx_2_inv) const;
        
        /*
         * Kernel to compute the derivatives in x-direction at nodes.
         */
        void computeFirstDerivativesInXAtNode(
            Real* dudx,
            const Real* const u,
            const hier::IntVector& num_ghosts_derivative_node,
            const hier::IntVector& num_ghosts_data_node,
            const hier::IntVector& ghostcell_dims_derivative_node,
            const hier::IntVector& ghostcell_dims_data_node,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims,
            const Real& dx_0_inv) const;
        
        /*
         * Kernel to compute the derivatives in y-direction at nodes.
         */
        void computeFirstDerivativesInYAtNode(
            Real* dudy,
            const Real* const u,
            const hier::IntVector& num_ghosts_derivative_node,
            const hier::IntVector& num_ghosts_data_node,
            const hier::IntVector& ghostcell_dims_derivative_node,
            const hier::IntVector& ghostcell_dims_data_node,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims,
            const Real& dx_1_inv) const;
        
        /*
         * Kernel to compute the derivatives in z-direction at nodes.
         */
        void computeFirstDerivativesInZAtNode(
            Real* dudz,
            const Real* const u,
            const hier::IntVector& num_ghosts_derivative_node,
            const hier::IntVector& num_ghosts_data_node,
            const hier::IntVector& ghostcell_dims_derivative_node,
            const hier::IntVector& ghostcell_dims_data_node,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims,
            const Real& dx_2_inv) const;
        
        /*
         * Kernel to interpolate the data from nodes to midpoints in x-direction.
         */
        void interpolateDataFromNodeToMidpointX(
            Real* u_midpoint_x,
            const Real* const u_node,
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
            Real* u_midpoint_y,
            const Real* const u_node,
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
            Real* u_midpoint_z,
            const Real* const u_node,
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
            Real* F_face_x,
            const Real* const F_midpoint_x,
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
            Real* F_face_y,
            const Real* const F_midpoint_y,
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
            Real* F_face_z,
            const Real* const F_midpoint_z,
            const hier::IntVector& num_ghosts_flux_midpoint,
            const hier::IntVector& ghostcell_dims_flux_midpoint,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims,
            const hier::IntVector& interior_dims,
            const double& dt) const;
        
};

#endif /* DIFFUSIVE_FLUX_RECONSTRUCTOR_MIDPOINT_SECOND_ORDER_HPP */
