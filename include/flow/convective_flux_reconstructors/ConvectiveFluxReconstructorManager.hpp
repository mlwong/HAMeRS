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
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_eqn,
            const HAMERS_SHARED_PTR<FlowModel>& flow_model,
            const HAMERS_SHARED_PTR<tbox::Database>& convective_flux_reconstructor_db,
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
        HAMERS_SHARED_PTR<ConvectiveFluxReconstructor>
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
         * HAMERS_SHARED_PTR to the convective flux reconstructor.
         */
        HAMERS_SHARED_PTR<ConvectiveFluxReconstructor> d_conv_flux_reconstructor;
        
};

#endif /* CONVECTIVE_FLUX_RECONSTRUCTOR_MANAGER_HPP */
