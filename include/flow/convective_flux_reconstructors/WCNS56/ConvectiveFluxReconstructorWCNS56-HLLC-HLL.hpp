#ifndef CONVECTIVE_FLUX_RECONSTRUCTOR_56_HLLC_HLL_HPP
#define CONVECTIVE_FLUX_RECONSTRUCTOR_56_HLLC_HLL_HPP

#include "flow/convective_flux_reconstructors/ConvectiveFluxReconstructor.hpp"
#include "util/derivatives/DerivativeFirstOrder.hpp"
#include "util/Directions.hpp"

class ConvectiveFluxReconstructorWCNS56: public ConvectiveFluxReconstructor
{
    public:
        ConvectiveFluxReconstructorWCNS56(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_eqn,
            const FLOW_MODEL::TYPE& flow_model_type,
            const HAMERS_SHARED_PTR<FlowModel>& flow_model,
            const HAMERS_SHARED_PTR<tbox::Database>& convective_flux_reconstructor_db);
        
        virtual ~ConvectiveFluxReconstructorWCNS56() {}
        
        /*
         * Print all characteristics of the convective flux reconstruction class.
         */
        virtual void
        printClassData(std::ostream& os) const = 0;
        
        /*
         * Put the characteristics of the convective flux reconstruction class
         * into the restart database.
         */
        virtual void
        putToRestart(
            const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const = 0;
        
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
        
    protected:
        /*
         * Perform WENO interpolation.
         */
        virtual void
        performWENOInterpolation(
            std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& variables_minus,
            std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& variables_plus,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > >& variables) = 0;
        
        /*
         * Forms of equations.
         */
        std::vector<EQN_FORM::TYPE> d_eqn_form;
        bool d_has_advective_eqn_form;
        
};

#endif /* CONVECTIVE_FLUX_RECONSTRUCTOR_56_HLLC_HLL_HPP */
