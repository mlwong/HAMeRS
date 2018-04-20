#ifndef NONCONSERVATIVE_DIFFUSIVE_FLUX_DIVERGENCE_OPERATOR_MANAGER_HPP
#define NONCONSERVATIVE_DIFFUSIVE_FLUX_DIVERGENCE_OPERATOR_MANAGER_HPP

#include "HAMeRS_config.hpp"

#include "flow/nonconservative_diffusive_flux_divergence_operators/NonconservativeDiffusiveFluxDivergenceOperators.hpp"

#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/tbox/Dimension.h"

#include "boost/shared_ptr.hpp"
#include <string>

using namespace SAMRAI;

class NonconservativeDiffusiveFluxDivergenceOperatorManager
{
    public:
        NonconservativeDiffusiveFluxDivergenceOperatorManager(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_eqn,
            const boost::shared_ptr<FlowModel>& flow_model,
            const boost::shared_ptr<tbox::Database>& nonconservative_diffusive_flux_divergence_operator_db,
            const std::string& nonconservative_diffusive_flux_divergence_operator_str);
        
        /*
         * Get the type of non-conservative diffusive flux divergence operator.
         */
        const NONCONSERVATIVE_DIFFUSIVE_FLUX_DIVERGENCE_OPERATOR::TYPE&
        getNonconservativeDiffusiveFluxDivergenceOperatorType() const
        {
            return d_nonconservative_diffusive_flux_divergence_operator_type;
        }
        
        /*
         * Get the non-conservative diffusive flux divergence operator.
         */
        boost::shared_ptr<NonconservativeDiffusiveFluxDivergenceOperator>
        getNonconservativeDiffusiveFluxDivergenceOperator() const
        {
            return d_noncons_diff_flux_div_op;
        }
        
        /*
         * Print all characteristics of diffusive flux divergence operator manager.
         */
        void
        printClassData(std::ostream& os) const;
        
    private:
        /*
         * The object name is used for error/warning reporting.
         */
        const std::string d_object_name;
        
        /*
         * Type of diffusive flux divergence operator.
         */
        NONCONSERVATIVE_DIFFUSIVE_FLUX_DIVERGENCE_OPERATOR::TYPE d_nonconservative_diffusive_flux_divergence_operator_type;
        
        /*
         * boost::shared_ptr to the diffusive flux divergence operator.
         */
        boost::shared_ptr<NonconservativeDiffusiveFluxDivergenceOperator> d_noncons_diff_flux_div_op;
        
};

#endif /* NONCONSERVATIVE_DIFFUSIVE_FLUX_DIVERGENCE_OPERATOR_MANAGER_HPP */
