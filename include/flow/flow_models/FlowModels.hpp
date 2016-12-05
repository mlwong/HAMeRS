#ifndef FLOW_MODELS_HPP
#define FLOW_MODELS_HPP

#include "flow/flow_models/FlowModel.hpp"
#include "flow/flow_models/single-species/FlowModelSingleSpecies.hpp"
#include "flow/flow_models/four-eqn_conservative/FlowModelFourEqnConservative.hpp"
#include "flow/flow_models/five-eqn_Allaire/FlowModelFiveEqnAllaire.hpp"

#include <map>
#include <string>

namespace FLOW_MODEL
{
    enum TYPE { SINGLE_SPECIES,
                FOUR_EQN_CONSERVATIVE,
                FIVE_EQN_ALLAIRE };
}

/*
 * Function to print out enum FLOW_MODEL::TYPE value as text.
 */
inline std::ostream& operator<<(std::ostream& os, const FLOW_MODEL::TYPE& value)
{
    static std::map<FLOW_MODEL::TYPE, std::string> strings;
    
    if (strings.size() == 0)
    {
#define INSERT_ELEMENT(p) strings[p] = #p
        INSERT_ELEMENT(FLOW_MODEL::SINGLE_SPECIES);
        INSERT_ELEMENT(FLOW_MODEL::FOUR_EQN_CONSERVATIVE);
        INSERT_ELEMENT(FLOW_MODEL::FIVE_EQN_ALLAIRE);
#undef INSERT_ELEMENT
    }
    
    return os << strings[value];
}

#endif /* FLOW_MODELS_HPP */
