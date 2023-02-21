#ifndef FLOW_MODEL_IMMERSED_BOUNDARY_METHOD_SINGLE_SPECIES_HPP
#define FLOW_MODEL_IMMERSED_BOUNDARY_METHOD_SINGLE_SPECIES_HPP

#include "flow/flow_models/FlowModelImmersedBoundaryMethod.hpp"

class FlowModelImmersedBoundaryMethodSingleSpecies: public FlowModelImmersedBoundaryMethod
{
    public:
        FlowModelImmersedBoundaryMethodSingleSpecies(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_species,
            const int& num_eqn,
            const HAMERS_SHARED_PTR<ImmersedBoundaries>& immersed_boundaries,
            const HAMERS_SHARED_PTR<tbox::Database>& immersed_boundary_method_db);
        
        ~FlowModelImmersedBoundaryMethodSingleSpecies() {}
        
};

#endif /* FLOW_MODEL_IMMERSED_BOUNDARY_METHOD_SINGLE_SPECIES_HPP */
