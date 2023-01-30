#ifndef FLOW_MODEL_IMMERSED_BOUNDARY_METHOD_HPP
#define FLOW_MODEL_IMMERSED_BOUNDARY_METHOD_HPP

#include "HAMeRS_config.hpp"

#include "HAMeRS_memory.hpp"

#include "flow/flow_models/FlowModel.hpp"
#include "util/immersed_boundaries/ImmersedBoundaries.hpp"

#include "SAMRAI/pdat/CellData.h"

class FlowModel;

class FlowModelImmersedBoundaryMethod
{
    public:
        FlowModelImmersedBoundaryMethod(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_species,
            const int& num_eqn,
            const HAMERS_SHARED_PTR<ImmersedBoundaries>& immersed_boundaries,
            const HAMERS_SHARED_PTR<tbox::Database>& immersed_boundary_method_db):
                d_object_name(object_name),
                d_dim(dim),
                d_grid_geometry(grid_geometry),
                d_num_species(num_species),
                d_num_eqn(num_eqn),
                d_immersed_boundaries(immersed_boundaries)
        {}
        
        virtual ~FlowModelImmersedBoundaryMethod() {}
        
        /*
         * Set the weak pointer to the flow model from the parent FlowModel class.
         */
        void setFlowModel(const HAMERS_WEAK_PTR<FlowModel>& flow_model)
        {
            d_flow_model = flow_model;
        }
        
        /*
         * Get the pointer to the immersed boundaries.
         */
        HAMERS_SHARED_PTR<ImmersedBoundaries> getImmersedBoundaries() const
        {
            return d_immersed_boundaries;
        }
        
        void putToRestart(const HAMERS_SHARED_PTR<tbox::Database>& immersed_boundary_method_db) const
        {
        }
        
    protected:
        /*
         * The object name is used for error/warning reporting.
         */
        const std::string d_object_name;
        
        /*
         * Problem dimension.
         */
        const tbox::Dimension d_dim;
        
        /*
         * HAMERS_SHARED_PTR to the grid geometry.
         */
        const HAMERS_SHARED_PTR<geom::CartesianGridGeometry> d_grid_geometry;
        
        /*
         * Number of species.
         */
        const int d_num_species;
        
        /*
         * Number of equations.
         */
        const int d_num_eqn;
        
        /*
         * Pointer to immersed boundaries.
         */
        HAMERS_SHARED_PTR<ImmersedBoundaries> d_immersed_boundaries;
        
        /*
         * HAMERS_WEAK_PTR to FlowModel.
         */
        HAMERS_WEAK_PTR<FlowModel> d_flow_model;
        
};

#endif /* FLOW_MODEL_BASIC_UTILITIES_HPP */
