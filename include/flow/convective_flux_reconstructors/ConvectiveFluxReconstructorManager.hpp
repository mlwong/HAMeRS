#ifndef CONVECTIVE_FLUX_RECONSTRUCTOR_MANAGER_HPP
#define CONVECTIVE_FLUX_RECONSTRUCTOR_MANAGER_HPP

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/tbox/Dimension.h"

#include "flow/convective_flux_reconstructors/ConvectiveFluxReconstructors.hpp"

#include "boost/shared_ptr.hpp"
#include <string>

using namespace SAMRAI;

class ConvectiveFluxReconstructorManager
{
    public:
        ConvectiveFluxReconstructorManager(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_eqn,
            const int& num_species,
            const boost::shared_ptr<FlowModel>& flow_model,
            const boost::shared_ptr<tbox::Database>& convective_flux_reconstructor_db,
            const std::string& convective_flux_reconstructor_str);
        
        /*
         * Get the label of equation of state.
         */
        const CONVECTIVE_FLUX_RECONSTRUCTOR_LABEL&
        getConvectiveFluxReconstructorLabel() const
        {
            return d_convective_flux_reconstructor_label;
        }
        
        /*
         * Get the convective flux reconstructor.
         */
        boost::shared_ptr<ConvectiveFluxReconstructor>
        getConvectiveFluxReconstructor() const
        {
            return d_conv_flux_reconstructor;
        }
        
        /*
         * Print all characteristics of convective flux reconstructor manager.
         */
        void
        printClassData(std::ostream& os) const;
        
    private:
        /*
         * The object name is used for error/warning reporting.
         */
        const std::string d_object_name;
        
        /*
         * Label of convective flux reconstructor.
         */
        CONVECTIVE_FLUX_RECONSTRUCTOR_LABEL d_convective_flux_reconstructor_label;
        
        /*
         * boost::shared_ptr to the convective flux reconstructor.
         */
        boost::shared_ptr<ConvectiveFluxReconstructor> d_conv_flux_reconstructor;
        
};

#endif /* CONVECTIVE_FLUX_RECONSTRUCTOR_MANAGER_HPP */
