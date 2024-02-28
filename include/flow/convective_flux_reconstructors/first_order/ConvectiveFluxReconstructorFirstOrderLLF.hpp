#ifndef CONVECTIVE_FLUX_RECONSTRUCTOR_LLF_HPP
#define CONVECTIVE_FLUX_RECONSTRUCTOR_LLF_HPP

#include "flow/convective_flux_reconstructors/ConvectiveFluxReconstructor.hpp"

class ConvectiveFluxReconstructorFirstOrderLLF: public ConvectiveFluxReconstructor
{
    public:
        ConvectiveFluxReconstructorFirstOrderLLF(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_eqn,
            const FLOW_MODEL::TYPE& flow_model_type,
            const HAMERS_SHARED_PTR<FlowModel>& flow_model,
            const HAMERS_SHARED_PTR<tbox::Database>& convective_flux_reconstructor_db);
        
        ~ConvectiveFluxReconstructorFirstOrderLLF() {}
        
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
        void computeConvectiveFluxAndSourceOnPatch(
            hier::Patch& patch,
            const HAMERS_SHARED_PTR<hier::CoarseFineBoundary> coarse_fine_bdry,
            const HAMERS_SHARED_PTR<pdat::SideVariable<Real> >& variable_convective_flux,
            const HAMERS_SHARED_PTR<pdat::CellVariable<Real> >& variable_source,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double time,
            const double dt,
            const int RK_step_number);
    
};

#endif /* CONVECTIVE_FLUX_RECONSTRUCTOR_HPP */
