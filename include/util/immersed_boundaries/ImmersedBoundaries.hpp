#ifndef IMMERSED_BOUNDARIES_HPP
#define IMMERSED_BOUNDARIES_HPP

#include "HAMeRS_config.hpp"

#include "HAMeRS_memory.hpp"

#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/pdat/CellData.h"

#include <string>

using namespace SAMRAI;

namespace IB_MASK
{
    enum TYPE { FLUID = 0,
                IB_GHOST,
                BODY };
}

class ImmersedBoundaries
{
    public:
        ImmersedBoundaries(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry):
                d_object_name(object_name),
                d_dim(dim),
                d_grid_geometry(grid_geometry)
        {}
        
        virtual ~ImmersedBoundaries() {}
        
        void setImmersedBoundaryVariablesOnPatch(
            const hier::Patch& patch,
            const HAMERS_SHARED_PTR<pdat::CellData<int> >& data_mask,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_wall_distance,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_surface_normal,
            const double data_time,
            const bool initial_time);
        
    private:
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
        
};

#endif /* IMMERSED_BOUNDARIES_HPP */
