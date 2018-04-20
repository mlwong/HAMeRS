#ifndef CONVECTIVE_FLUX_RECONSTRUCTOR_WCNS6_LD_HLLC_HLL_HPP
#define CONVECTIVE_FLUX_RECONSTRUCTOR_WCNS6_LD_HLLC_HLL_HPP

#include "flow/convective_flux_reconstructors/WCNS56/ConvectiveFluxReconstructorWCNS56-HLLC-HLL.hpp"

class ConvectiveFluxReconstructorWCNS6_LD_HLLC_HLL: public ConvectiveFluxReconstructorWCNS56
{
    public:
        ConvectiveFluxReconstructorWCNS6_LD_HLLC_HLL(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_eqn,
            const boost::shared_ptr<FlowModel>& flow_model,
            const boost::shared_ptr<tbox::Database>& convective_flux_reconstructor_db);
        
        ~ConvectiveFluxReconstructorWCNS6_LD_HLLC_HLL() {}
        
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
            const boost::shared_ptr<tbox::Database>& restart_db) const;
        
    private:
        /*
         * Perform WENO interpolation.
         */
        void
        performWENOInterpolation(
            std::vector<boost::shared_ptr<pdat::SideData<double> > >& variables_minus,
            std::vector<boost::shared_ptr<pdat::SideData<double> > >& variables_plus,
            const std::vector<std::vector<boost::shared_ptr<pdat::SideData<double> > > >& variables);
        
        /*
         * Constants used by the scheme.
         */
        int d_constant_p;
        int d_constant_q;
        double d_constant_C;
        double d_constant_alpha_tau;
        
};

#endif /* CONVECTIVE_FLUX_RECONSTRUCTOR_WCNS6_LD_HLLC_HLL_HPP */
