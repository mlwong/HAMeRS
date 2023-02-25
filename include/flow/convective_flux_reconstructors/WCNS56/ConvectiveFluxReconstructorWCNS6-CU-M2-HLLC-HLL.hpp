#ifndef CONVECTIVE_FLUX_RECONSTRUCTOR_WCNS6_CU_M2_HLLC_HLL_HPP
#define CONVECTIVE_FLUX_RECONSTRUCTOR_WCNS6_CU_M2_HLLC_HLL_HPP

#include "flow/convective_flux_reconstructors/WCNS56/ConvectiveFluxReconstructorWCNS56-HLLC-HLL.hpp"

class ConvectiveFluxReconstructorWCNS6_CU_M2_HLLC_HLL: public ConvectiveFluxReconstructorWCNS56
{
    public:
        ConvectiveFluxReconstructorWCNS6_CU_M2_HLLC_HLL(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_eqn,
            const FLOW_MODEL::TYPE& flow_model_type,
            const HAMERS_SHARED_PTR<FlowModel>& flow_model,
            const HAMERS_SHARED_PTR<tbox::Database>& convective_flux_reconstructor_db);
        
        ~ConvectiveFluxReconstructorWCNS6_CU_M2_HLLC_HLL() {}
        
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
        
    private:
        /*
         * Perform WENO interpolation.
         */
        void
        performWENOInterpolation(
            std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& variables_minus,
            std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& variables_plus,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > >& variables);
        
        /*
         * Constants used by the scheme.
         */
        int d_constant_q;
        Real d_constant_C;
        Real d_constant_Chi;
        Real d_constant_epsilon;
        
};

#endif /* CONVECTIVE_FLUX_RECONSTRUCTOR_WCNS6_CU_M2_HLLC_HLL_HPP */
