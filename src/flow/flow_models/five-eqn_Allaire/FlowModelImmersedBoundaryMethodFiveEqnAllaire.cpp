#include "flow/flow_models/five-eqn_Allaire/FlowModelImmersedBoundaryMethodFiveEqnAllaire.hpp"

FlowModelImmersedBoundaryMethodFiveEqnAllaire::FlowModelImmersedBoundaryMethodFiveEqnAllaire(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
    const int& num_species,
    const int& num_eqn,
    const HAMERS_SHARED_PTR<ImmersedBoundaries>& immersed_boundaries,
    const HAMERS_SHARED_PTR<tbox::Database>& flow_model_db):
        FlowModelImmersedBoundaryMethod(
            object_name,
            dim,
            grid_geometry,
            num_species,
            num_eqn,
            immersed_boundaries,
            flow_model_db)
{
}
