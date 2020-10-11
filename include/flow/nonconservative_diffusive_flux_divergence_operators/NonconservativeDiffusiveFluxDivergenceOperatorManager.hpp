#ifndef NONCONSERVATIVE_DIFFUSIVE_FLUX_DIVERGENCE_OPERATOR_MANAGER_HPP
#define NONCONSERVATIVE_DIFFUSIVE_FLUX_DIVERGENCE_OPERATOR_MANAGER_HPP

#include "HAMeRS_config.hpp"

#include "HAMeRS_memory.hpp"

#include "flow/nonconservative_diffusive_flux_divergence_operators/NonconservativeDiffusiveFluxDivergenceOperators.hpp"

#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/tbox/Dimension.h"

#include <string>

using namespace SAMRAI;

class NonconservativeDiffusiveFluxDivergenceOperatorManager
{
    public:
        NonconservativeDiffusiveFluxDivergenceOperatorManager(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_eqn,
            const HAMERS_SHARED_PTR<FlowModel>& flow_model,
            const HAMERS_SHARED_PTR<tbox::Database>& nonconservative_diffusive_flux_divergence_operator_db,
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
        HAMERS_SHARED_PTR<NonconservativeDiffusiveFluxDivergenceOperator>
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
         * HAMERS_SHARED_PTR to the diffusive flux divergence operator.
         */
        HAMERS_SHARED_PTR<NonconservativeDiffusiveFluxDivergenceOperator> d_noncons_diff_flux_div_op;
        
};

#endif /* NONCONSERVATIVE_DIFFUSIVE_FLUX_DIVERGENCE_OPERATOR_MANAGER_HPP */
