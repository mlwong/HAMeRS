#ifndef FLOW_MODEL_SUBGRID_SCALE_MODEL_HPP
#define FLOW_MODEL_SUBGRID_SCALE_MODEL_HPP

#include "HAMeRS_config.hpp"

#include "HAMeRS_memory.hpp"

#include "util/Directions.hpp"

#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/SideData.h"

#include <string>
#include <vector>

using namespace SAMRAI;

namespace SUBGRID_SCALE_MODEL
{
    enum TYPE { VREMAN };
}

class FlowModelSubgridScaleModel
{
    public:
        FlowModelSubgridScaleModel(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_species,
            const int& num_eqn,
            const HAMERS_SHARED_PTR<tbox::Database>& subgrid_scale_model_db):
                d_object_name(object_name),
                d_dim(dim),
                d_grid_geometry(grid_geometry),
                d_num_species(num_species),
                d_num_eqn(num_eqn)
        {
            NULL_USE(subgrid_scale_model_db);
        }
        
        virtual ~FlowModelSubgridScaleModel() {}
        
        /*
         * Return names of different derived variables related to this class in the registered patch.
         */
        virtual std::vector<std::string>
        getDerivedVariablesToRegister() const
        {
            std::vector<std::string> var_to_register;
            return var_to_register;
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
        
};


#endif /* FLOW_MODEL_SUBGRID_SCALE_MODEL_HPP */