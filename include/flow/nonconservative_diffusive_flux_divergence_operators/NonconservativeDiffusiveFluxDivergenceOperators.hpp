#ifndef NONCONSERVATIVE_DIFFUSIVE_FLUX_DIVERGENCE_OPERATORS_HPP
#define NONCONSERVATIVE_DIFFUSIVE_FLUX_DIVERGENCE_OPERATORS_HPP

#include "flow/nonconservative_diffusive_flux_divergence_operators/sixth_order/NonconservativeDiffusiveFluxDivergenceOperatorSixthOrder.hpp"

#include <map>
#include <string>

namespace NONCONSERVATIVE_DIFFUSIVE_FLUX_DIVERGENCE_OPERATOR
{
    enum TYPE { SIXTH_ORDER };
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
#undef INSERT_ELEMENT
    }
    
    return os << strings[value];
}

#endif /* NONCONSERVATIVE_DIFFUSIVE_FLUX_DIVERGENCE_OPERATORS_HPP */
