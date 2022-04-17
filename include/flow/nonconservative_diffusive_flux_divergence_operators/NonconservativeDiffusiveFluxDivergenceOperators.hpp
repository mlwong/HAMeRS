#ifndef NONCONSERVATIVE_DIFFUSIVE_FLUX_DIVERGENCE_OPERATORS_HPP
#define NONCONSERVATIVE_DIFFUSIVE_FLUX_DIVERGENCE_OPERATORS_HPP

#include "flow/nonconservative_diffusive_flux_divergence_operators/NonconservativeDiffusiveFluxDivergenceOperatorSixthOrder.hpp"
#include "flow/nonconservative_diffusive_flux_divergence_operators/NonconservativeDiffusiveFluxDivergenceOperatorEighthOrder.hpp"
#include "flow/nonconservative_diffusive_flux_divergence_operators/NonconservativeDiffusiveFluxDivergenceOperatorTenthOrder.hpp"
#include "flow/nonconservative_diffusive_flux_divergence_operators/NonconservativeDiffusiveFluxDivergenceOperatorTwelfthOrder.hpp"

#include <map>
#include <string>

namespace NONCONSERVATIVE_DIFFUSIVE_FLUX_DIVERGENCE_OPERATOR
{
    enum TYPE { SIXTH_ORDER,
                EIGHTH_ORDER,
                TENTH_ORDER,
                TWELFTH_ORDER };
}

/*
 * Function to print out enum NONCONSERVATIVE_DIFFUSIVE_FLUX_DIVERGENCE_OPERATOR::TYPE value as text.
 */
inline std::ostream& operator<<(
    std::ostream& os,
    const NONCONSERVATIVE_DIFFUSIVE_FLUX_DIVERGENCE_OPERATOR::TYPE& value)
{
    static std::map<NONCONSERVATIVE_DIFFUSIVE_FLUX_DIVERGENCE_OPERATOR::TYPE, std::string> strings;
    
    if (strings.size() == 0)
    {
#define INSERT_ELEMENT(p) strings[p] = #p
        INSERT_ELEMENT(NONCONSERVATIVE_DIFFUSIVE_FLUX_DIVERGENCE_OPERATOR::SIXTH_ORDER);
        INSERT_ELEMENT(NONCONSERVATIVE_DIFFUSIVE_FLUX_DIVERGENCE_OPERATOR::EIGHTH_ORDER);
        INSERT_ELEMENT(NONCONSERVATIVE_DIFFUSIVE_FLUX_DIVERGENCE_OPERATOR::TENTH_ORDER);
        INSERT_ELEMENT(NONCONSERVATIVE_DIFFUSIVE_FLUX_DIVERGENCE_OPERATOR::TWELFTH_ORDER);
#undef INSERT_ELEMENT
    }
    
    return os << strings[value];
}

#endif /* NONCONSERVATIVE_DIFFUSIVE_FLUX_DIVERGENCE_OPERATORS_HPP */
