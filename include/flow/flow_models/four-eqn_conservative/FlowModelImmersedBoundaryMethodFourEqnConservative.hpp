#ifndef FLOW_MODEL_IMMERSED_BOUNDARY_METHOD_FOUR_EQN_CONSERVATIVE_HPP
#define FLOW_MODEL_IMMERSED_BOUNDARY_METHOD_FOUR_EQN_CONSERVATIVE_HPP

#include "flow/flow_models/FlowModelImmersedBoundaryMethod.hpp"

class FlowModelImmersedBoundaryMethodFourEqnConservative: public FlowModelImmersedBoundaryMethod
{
    public:
        FlowModelImmersedBoundaryMethodFourEqnConservative(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_species,
            const int& num_eqn,
            const HAMERS_SHARED_PTR<ImmersedBoundaries>& immersed_boundaries,
            const HAMERS_SHARED_PTR<tbox::Database>& immersed_boundary_method_db);
        
        ~FlowModelImmersedBoundaryMethodFourEqnConservative() {}
        
};

#endif /* FLOW_MODEL_IMMERSED_BOUNDARY_METHOD_FOUR_EQN_CONSERVATIVE_HPP */
