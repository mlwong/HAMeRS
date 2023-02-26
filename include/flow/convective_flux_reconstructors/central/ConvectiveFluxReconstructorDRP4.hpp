#ifndef CONVECTIVE_FLUX_RECONSTRUCTOR_DRP4_HPP
#define CONVECTIVE_FLUX_RECONSTRUCTOR_DRP4_HPP

#include "flow/convective_flux_reconstructors/ConvectiveFluxReconstructor.hpp"
#include "util/derivatives/DerivativeFirstOrder.hpp"
#include "util/Directions.hpp"

#include "SAMRAI/pdat/SideVariable.h"

// Follow the DRP schemes in
// Bogey, Christophe, and Christophe Bailly.
// "A family of low dispersive and low dissipative explicit schemes for flow and noise computations."
// Journal of Computational physics 194.1 (2004): 194-214.

class ConvectiveFluxReconstructorDRP4: public ConvectiveFluxReconstructor
{
    public:
        ConvectiveFluxReconstructorDRP4(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_eqn,
            const FLOW_MODEL::TYPE& flow_model_type,
            const HAMERS_SHARED_PTR<FlowModel>& flow_model,
            const HAMERS_SHARED_PTR<tbox::Database>& convective_flux_reconstructor_db);
        
        ~ConvectiveFluxReconstructorDRP4();
        
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
         * Options of the scheme.
         */
        int d_stencil_width;
        
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
        
};

#endif /* CONVECTIVE_FLUX_RECONSTRUCTOR_DRP4_HPP */
