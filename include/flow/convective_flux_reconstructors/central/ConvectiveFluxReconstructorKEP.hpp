#ifndef CONVECTIVE_FLUX_RECONSTRUCTOR_KEP_HPP
#define CONVECTIVE_FLUX_RECONSTRUCTOR_KEP_HPP

#include "flow/convective_flux_reconstructors/ConvectiveFluxReconstructor.hpp"
#include "util/derivatives/DerivativeFirstOrder.hpp"
#include "util/Directions.hpp"

#include "SAMRAI/pdat/SideVariable.h"

// Follow the kinetic energy preserving schemes in
// Pirozzoli, Sergio.
// "Generalized conservative approximations of split convective derivative operators."
// Journal of Computational Physics 229.19 (2010): 7180-7190.

class ConvectiveFluxReconstructorKEP: public ConvectiveFluxReconstructor
{
    public:
        ConvectiveFluxReconstructorKEP(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_eqn,
            const FLOW_MODEL::TYPE& flow_model_type,
            const HAMERS_SHARED_PTR<FlowModel>& flow_model,
            const HAMERS_SHARED_PTR<tbox::Database>& convective_flux_reconstructor_db);
        
        ~ConvectiveFluxReconstructorKEP();
        
        /*
         * Print all characteristics of the convective flux reconstruction class.
         */
        void
        printClassData(std::ostream& os) const;
        
        /*
         * Put the characteristics of the convective flux reconstruction class
         * into the restart database.
         */
        void
        putToRestart(
            const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const;
        
        /*
         * Compute the convective flux and source due to splitting of convective term on a patch.
         */
        void
        computeConvectiveFluxAndSourceOnPatch(
            hier::Patch& patch,
            const HAMERS_SHARED_PTR<pdat::SideVariable<Real> >& variable_convective_flux,
            const HAMERS_SHARED_PTR<pdat::CellVariable<Real> >& variable_source,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double time,
            const double dt,
            const int RK_step_number);
        
    private:
        /*
         * Add linear term to convective flux in x-direction.
         */
        void
        addLinearTermToConvectiveFluxX(
            HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_convective_flux,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_f,
            const int component_idx_flux,
            const int component_idx_f,
            const double dt) const;
        
        /*
         * Add quadratic term to convective flux in x-direction.
         */
        void
        addQuadraticTermToConvectiveFluxX(
            HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_convective_flux,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_f,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_g,
            const int component_idx_flux,
            const int component_idx_f,
            const int component_idx_g,
            const double dt) const;
        
        /*
         * Add cubic term to convective flux in x-direction.
         */
        void
        addCubicTermToConvectiveFluxX(
            HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_convective_flux,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_f,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_g,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_h,
            const int component_idx_flux,
            const int component_idx_f,
            const int component_idx_g,
            const int component_idx_h,
            const double dt) const;
        
        /*
         * Add linear term to convective flux in y-direction.
         */
        void
        addLinearTermToConvectiveFluxY(
            HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_convective_flux,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_f,
            const int component_idx_flux,
            const int component_idx_f,
            const double dt) const;
        
        /*
         * Add quadratic term to convective flux in y-direction.
         */
        void
        addQuadraticTermToConvectiveFluxY(
            HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_convective_flux,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_f,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_g,
            const int component_idx_flux,
            const int component_idx_f,
            const int component_idx_g,
            const double dt) const;
        
        /*
         * Add cubic term to convective flux in y-direction.
         */
        void
        addCubicTermToConvectiveFluxY(
            HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_convective_flux,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_f,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_g,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_h,
            const int component_idx_flux,
            const int component_idx_f,
            const int component_idx_g,
            const int component_idx_h,
            const double dt) const;
        
        /*
         * Add linear term to convective flux in z-direction.
         */
        void
        addLinearTermToConvectiveFluxZ(
            HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_convective_flux,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_f,
            const int component_idx_flux,
            const int component_idx_f,
            const double dt) const;
        
        /*
         * Add quadratic term to convective flux in z-direction.
         */
        void
        addQuadraticTermToConvectiveFluxZ(
            HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_convective_flux,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_f,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_g,
            const int component_idx_flux,
            const int component_idx_f,
            const int component_idx_g,
            const double dt) const;
        
        /*
         * Add cubic term to convective flux in z-direction.
         */
        void
        addCubicTermToConvectiveFluxZ(
            HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_convective_flux,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_f,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_g,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_h,
            const int component_idx_flux,
            const int component_idx_f,
            const int component_idx_g,
            const int component_idx_h,
            const double dt) const;
        
        /*
         * Add source terms to the advection equations of volume fractions.
         * (for five-equation model by Allaire et al.)
         */
        void
        addSourceTermsToVolumeFractionEquations(
            HAMERS_SHARED_PTR<pdat::CellData<Real> > data_source,
            HAMERS_SHARED_PTR<pdat::CellData<Real> > data_velocity,
            HAMERS_SHARED_PTR<pdat::CellData<Real> > data_volume_fractions,
            const double* const dx,
            const double dt) const;
        
        /*
         * Options of the scheme.
         */
        bool d_use_DRP4;
        int d_stencil_width;
        int d_order;
        
        /*
         * Forms of equations.
         */
        std::vector<EQN_FORM::TYPE> d_eqn_form;
        bool d_has_advective_eqn_form;
        
        /*
         * Timers interspersed throughout the class.
         */
        static HAMERS_SHARED_PTR<tbox::Timer> t_reconstruct_flux;
        static HAMERS_SHARED_PTR<tbox::Timer> t_compute_source;
        
        /*
         * Coefficients for finite differencing.
         */
         
        Real d_coef_a;
        Real d_coef_b;
        Real d_coef_c;
        Real d_coef_d;
        Real d_coef_e;
        Real d_coef_f;
};

#endif /* CONVECTIVE_FLUX_RECONSTRUCTOR_KEP_HPP */
