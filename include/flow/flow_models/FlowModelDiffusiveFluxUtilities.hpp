#ifndef FLOW_MODEL_DIFFUSIVE_FLUX_UTILITIES_HPP
#define FLOW_MODEL_DIFFUSIVE_FLUX_UTILITIES_HPP

#include "HAMeRS_config.hpp"

#include "HAMeRS_memory.hpp"

#include "flow/flow_models/FlowModel.hpp"
#include "flow/flow_models/FlowModelSubgridScaleModel.hpp"

#include "util/Directions.hpp"

#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/SideData.h"

#include <string>
#include <unordered_map>

class FlowModel;

class FlowModelDiffusiveFluxUtilities
{
    public:
        FlowModelDiffusiveFluxUtilities(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_species,
            const int& num_eqn,
            const HAMERS_SHARED_PTR<tbox::Database>& flow_model_db);
        
        virtual ~FlowModelDiffusiveFluxUtilities() {}
        
        /*
         * Set the weak pointer to the flow model from the parent FlowModel class.
         */
        void setFlowModel(const HAMERS_WEAK_PTR<FlowModel>& flow_model)
        {
            d_flow_model = flow_model;
        }
        
        /*
         * Register different derived variables related to this class in the registered patch. The
         * derived variables to be registered are given as entries in a map of the variable name to
         * the number of sub-ghost cells required.
         */
        virtual void
        registerDerivedVariables(
            const std::unordered_map<std::string, hier::IntVector>& num_subghosts_of_data);
        
        /*
         * Register the required variables for the computation of diffusive fluxes in the registered patch.
         */
        virtual void
        registerDerivedVariablesForDiffusiveFluxes(
            const hier::IntVector& num_subghosts,
            const bool need_side_diffusivities = false);
        
        /*
         * Allocate memory for cell data of different registered derived variables related to this
         * class in the registered patch.
         */
        virtual void allocateMemoryForDerivedCellData();
        
        /*
         * Allocate memory for side data of the diffusivities.
         */
        virtual void allocateMemoryForSideDataOfDiffusiveFluxDiffusivities();
        
        /*
         * Clear cell and side data of different derived variables related to this class in the registered patch.
         */
        virtual void clearCellAndSideData();
        
        /*
         * Compute cell data of different registered derived variables related to this class.
         */
        virtual void computeDerivedCellData();
        
        /*
         * Get the cell data of one cell variable related to this class in the registered patch.
         */
        virtual HAMERS_SHARED_PTR<pdat::CellData<double> >
        getCellData(const std::string& variable_key);
        
        /*
         * Get the cell data of different cell variables related to this class in the registered patch.
         */
        virtual std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >
        getCellData(
            const std::vector<std::string>& variable_keys);
        
        /*
         * Get the variables for the derivatives in the diffusive fluxes.
         */
        virtual void
        getCellDataOfDiffusiveFluxVariablesForDerivative(
            std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& derivative_var_data,
            std::vector<std::vector<int> >& derivative_var_component_idx,
            const DIRECTION::TYPE& flux_direction,
            const DIRECTION::TYPE& derivative_direction);
        
        /*
         * Get the diffusivities in the diffusive flux.
         */
        virtual void
        getCellDataOfDiffusiveFluxDiffusivities(
            std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& diffusivities_data,
            std::vector<std::vector<int> >& diffusivities_component_idx,
            const DIRECTION::TYPE& flux_direction,
            const DIRECTION::TYPE& derivative_direction);
        
        /*
         * Get the cell data that needs interpolation to sides for computing side data of diffusivities in the
         * diffusive flux.
         */
        virtual void
        getCellDataForInterpolationToSideDataForDiffusiveFluxDiffusivities(
            std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& var_data_for_diffusivities,
            std::vector<int>& var_data_for_diffusivities_component_idx);
        
        /*
         * Compute the side data of the diffusivities in the diffusive flux with the interpolated side data.
         */
        void
        virtual computeSideDataOfDiffusiveFluxDiffusivities(
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& var_data_for_diffusivities);
        
        /*
         * Get the side data of the diffusivities in the diffusive fluxa.
         */
        void
        virtual getSideDataOfDiffusiveFluxDiffusivities(
            std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > >& diffusivities_data,
            std::vector<std::vector<int> >& diffusivities_component_idx,
            const DIRECTION::TYPE& flux_direction,
            const DIRECTION::TYPE& derivative_direction);
        
        bool useSubgridScaleModel() const
        {
            return d_use_subgrid_scale_model;
        }
        
        /*
         * Get the variables for the derivatives used at computing subgrid scale diffusivity/viscosity at sides.
         */
        void
        getCellDataOfVariablesForSideDerivativeForSubgridScaleViscosity(
            std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& derivative_var_data,
            std::vector<int>& derivative_var_component_idx,
            const DIRECTION::TYPE& side_direction,
            const DIRECTION::TYPE& derivative_direction);
        
        /*
         * Modify the side data of the diffusivities/viscosities at sides with subgrid scale diffusivity/viscosity.
         */
        void
        updateSideDataOfDiffusiveFluxDiffusivitiesWithSubgridScaleModel(
            std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& var_data_for_diffusivities,
            const std::map<DIRECTION::TYPE, std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > >& derivatives,
            const DIRECTION::TYPE& side_direction);
        
        /*
         * Put the characteristics of this class into the restart database.
         */
        virtual void
        putToRestart(
            const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const;
        
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
         * Whether all derived cell data related to this class is computed in full domain or sub-domain.
         */
        bool d_derived_cell_data_computed;
        
        /*
         * Number of sub-ghost cells of diffusivities.
         */
        hier::IntVector d_num_subghosts_diffusivities;
        
        /*
         * Box with sub-ghost cells of diffusivities.
         */
        hier::Box d_subghost_box_diffusivities;
        
        /*
         * Dimensions of box with sub-ghost cells of diffusivities.
         */
        hier::IntVector d_subghostcell_dims_diffusivities;
        
        /*
         * HAMERS_SHARED_PTR to cell data of diffusivities.
         */
        HAMERS_SHARED_PTR<pdat::CellData<double> > d_data_diffusivities;
        
        /*
         * Whether cell data of diffusivities is computed.
         */
        bool d_cell_data_computed_diffusivities;
        
        /*
         * HAMERS_SHARED_PTR to side data of diffusivities.
         */
        HAMERS_SHARED_PTR<pdat::SideData<double> > d_side_data_diffusivities;
        
        /*
         * Whether side data of diffusivities is computed.
         */
        bool d_side_data_diffusivities_computed;
        
        /*
         * HAMERS_WEAK_PTR to FlowModel.
         */
        HAMERS_WEAK_PTR<FlowModel> d_flow_model;
        
        /*
         * Whether side data of diffusivities is needed.
         */
        bool d_need_side_diffusivities;
        
        /*
         * Whether to use subgrid-scale model.
         */
        bool d_use_subgrid_scale_model;
        
        /*
         * HAMERS_SHARED_PTR to subgrid-scale model object.
         */
        HAMERS_SHARED_PTR<FlowModelSubgridScaleModel> d_flow_model_subgrid_scale_model;
};

#endif /* FLOW_MODEL_DIFFUSIVE_FLUX_UTILITIES_HPP */
