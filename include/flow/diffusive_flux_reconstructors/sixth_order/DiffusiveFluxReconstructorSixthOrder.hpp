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
            hier::Patch& patch,
            std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& derivatives_x,
            std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > >& derivatives_x_computed,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& data_x,
            const std::vector<std::vector<int> >& data_component_idx_x,
            const bool allocate_scratch_derivatives_node = true);
        
        /*
         * Compute the first derivatives in the y-direction.
         */
        void computeFirstDerivativesInY(
            hier::Patch& patch,
            std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& derivatives_y,
            std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > >& derivatives_y_computed,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& data_y,
            const std::vector<std::vector<int> >& data_component_idx_y,
            const bool allocate_scratch_derivatives_node = true);
        
        /*
         * Compute the first derivatives in the z-direction.
         */
        void computeFirstDerivativesInZ(
            hier::Patch& patch,
            std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& derivatives_z,
            std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > >& derivatives_z_computed,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& data_z,
            const std::vector<std::vector<int> >& data_component_idx_z,
            const bool allocate_scratch_derivatives_node = true);
        
        /*
         * Reconstruct the flux using flux at nodes in x-direction.
         */
        void reconstructFluxX(
            hier::Patch& patch,
            HAMERS_SHARED_PTR<pdat::SideData<double> >& diffusive_flux,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& diffusive_flux_node,
            const double dt);
        
        /*
         * Reconstruct the flux using flux at nodes in y-direction.
         */
        void reconstructFluxY(
            hier::Patch& patch,
            HAMERS_SHARED_PTR<pdat::SideData<double> >& diffusive_flux,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& diffusive_flux_node,
            const double dt);
        
        /*
         * Reconstruct the flux using flux at nodes in z-direction.
         */
        void reconstructFluxZ(
            hier::Patch& patch,
            HAMERS_SHARED_PTR<pdat::SideData<double> >& diffusive_flux,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& diffusive_flux_node,
            const double dt);
        
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
