#ifndef CONVECTIVE_FLUX_RECONSTRUCTOR_MANAGER_HPP
#define CONVECTIVE_FLUX_RECONSTRUCTOR_MANAGER_HPP

#include "HAMeRS_config.hpp"

#include "HAMeRS_memory.hpp"

#include "flow/convective_flux_reconstructors/ConvectiveFluxReconstructors.hpp"

#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/tbox/Dimension.h"

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
            const boost::shared_ptr<FlowModel>& flow_model,
            const boost::shared_ptr<tbox::Database>& convective_flux_reconstructor_db,
            const std::string& convective_flux_reconstructor_str);
        
        /*
         * Get the type of convective flux reconstuctor.
         */
        const CONVECTIVE_FLUX_RECONSTRUCTOR::TYPE&
        getConvectiveFluxReconstructorType() const
        {
            return d_convective_flux_reconstructor_type;
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
         * Type of convective flux reconstructor.
         */
        CONVECTIVE_FLUX_RECONSTRUCTOR::TYPE d_convective_flux_reconstructor_type;
        
        /*
         * boost::shared_ptr to the convective flux reconstructor.
         */
        boost::shared_ptr<ConvectiveFluxReconstructor> d_conv_flux_reconstructor;
        
};

#endif /* CONVECTIVE_FLUX_RECONSTRUCTOR_MANAGER_HPP */
