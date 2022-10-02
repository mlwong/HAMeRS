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
            const HAMERS_SHARED_PTR<tbox::Database>& subgrid_scale_model_db);
        
        virtual ~FlowModelSubgridScaleModel() {}
        
        /*
         * Return names of different derived variables required to register.
         */
        virtual std::vector<std::string> getDerivedVariablesToRegister() const;
        
        /*
         * Return different derived variables required for interpolation.
         */
        void getDerivedVariablesForInterpolationToSideData(
            std::vector<std::string>& var_to_interpolate,
            std::vector<int>& var_to_interpolate_component_idx) const;
        
        /*
         * Get the variables for the derivatives used at computing subgrid scale diffusivity/viscosity at sides.
         */
        virtual void
        getCellDataOfVariablesForSideDerivativeForSubgridScaleViscosity(
            std::vector<std::string>& derivative_var_data_str,
            std::vector<int>& derivative_var_component_idx,
            const DIRECTION::TYPE& side_direction,
            const DIRECTION::TYPE& derivative_direction);
        
        /*
         * Modify the side data of the diffusivities/viscosities at sides with subgrid scale diffusivity/viscosity.
         */
        virtual void
        updateSideDataOfDiffusiveFluxDiffusivities(
            std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& var_data_for_diffusivities,
            const std::map<DIRECTION::TYPE, std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > >& derivatives,
            const DIRECTION::TYPE& side_direction);
        
        /*
         * Put the characteristics of this class into the restart database.
         */
        virtual void
        putToRestart(
            const HAMERS_SHARED_PTR<tbox::Database>& restart_subgrid_scale_model_db) const;
        
    protected:
        /*
         * Put the characteristics of base class into the restart database.
         */
        void
        putToRestartBase(
            const HAMERS_SHARED_PTR<tbox::Database>& restart_subgrid_scale_model_db) const;
        
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
         * Whether to use subgrid-scale model.
         */
        SUBGRID_SCALE_MODEL::TYPE d_subgrid_scale_model_type;
        
};

#endif /* FLOW_MODEL_SUBGRID_SCALE_MODEL_HPP */