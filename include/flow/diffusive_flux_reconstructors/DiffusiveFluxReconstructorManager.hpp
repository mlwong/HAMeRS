#ifndef DIFFUSIVE_FLUX_RECONSTRUCTOR_MANAGER_HPP
#define DIFFUSIVE_FLUX_RECONSTRUCTOR_MANAGER_HPP

#include "HAMeRS_config.hpp"

#include "HAMeRS_memory.hpp"

#include "flow/diffusive_flux_reconstructors/DiffusiveFluxReconstructors.hpp"

#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/tbox/Dimension.h"

#include <string>

using namespace SAMRAI;

class DiffusiveFluxReconstructorManager
{
    public:
        DiffusiveFluxReconstructorManager(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_eqn,
            const boost::shared_ptr<FlowModel>& flow_model,
            const boost::shared_ptr<tbox::Database>& diffusive_flux_reconstructor_db,
            const std::string& diffusive_flux_reconstructor_str);
        
        /*
         * Get the type of diffusive flux reconstructor.
         */
        const DIFFUSIVE_FLUX_RECONSTRUCTOR::TYPE&
        getDiffusiveFluxReconstructorType() const
        {
            return d_diffusive_flux_reconstructor_type;
        }
        
        /*
         * Get the diffusive flux reconstructor.
         */
        boost::shared_ptr<DiffusiveFluxReconstructor>
        getDiffusiveFluxReconstructor() const
        {
            return d_diff_flux_reconstructor;
        }
        
        /*
         * Print all characteristics of diffusive flux reconstructor manager.
         */
        void
        printClassData(std::ostream& os) const;
        
    private:
        /*
         * The object name is used for error/warning reporting.
         */
        const std::string d_object_name;
        
        /*
         * Type of diffusive flux reconstructor.
         */
        DIFFUSIVE_FLUX_RECONSTRUCTOR::TYPE d_diffusive_flux_reconstructor_type;
        
        /*
         * boost::shared_ptr to the diffusive flux reconstructor.
         */
        boost::shared_ptr<DiffusiveFluxReconstructor> d_diff_flux_reconstructor;
        
};

#endif /* DIFFUSIVE_FLUX_RECONSTRUCTOR_MANAGER_HPP */
