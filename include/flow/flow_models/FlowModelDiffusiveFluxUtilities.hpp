#ifndef FLOW_MODEL_DIFFUSIVE_FLUX_UTILITIES_HPP
#define FLOW_MODEL_DIFFUSIVE_FLUX_UTILITIES_HPP

#include "HAMeRS_config.hpp"

#include "flow/flow_models/FlowModel.hpp"

#include "util/Directions.hpp"

#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/SideData.h"

#include "boost/weak_ptr.hpp"
#include <string>

class FlowModel;

class FlowModelDiffusiveFluxUtilities
{
    public:
        FlowModelDiffusiveFluxUtilities(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_species,
            const int& num_eqn);
        
        virtual ~FlowModelDiffusiveFluxUtilities() {}
        
        /*
         * Set the weak pointer to the flow model from the parent FlowModel class.
         */
        void setFlowModel(const boost::weak_ptr<FlowModel>& flow_model)
        {
            d_flow_model = flow_model;
        }
        
        /*
         * Register the required variables for the computation of diffusive fluxes in the registered patch.
         */
        virtual void
        registerDerivedVariablesForDiffusiveFluxes(
            const hier::IntVector& num_subghosts);
        
        /*
         * The cell data of all derived variables in the patch for this class are dumped.
         */
        virtual void clearCellData()
        {}
        
        /*
         * Get the variables for the derivatives in the diffusive fluxes.
         */
        virtual void
        getCellDataOfDiffusiveFluxVariablesForDerivative(
            std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > >& derivative_var_data,
            std::vector<std::vector<int> >& derivative_var_component_idx,
            const DIRECTION::TYPE& flux_direction,
            const DIRECTION::TYPE& derivative_direction);
        
        /*
         * Get the diffusivities in the diffusive flux.
         */
        virtual void
        getCellDataOfDiffusiveFluxDiffusivities(
            std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > >& diffusivities_data,
            std::vector<std::vector<int> >& diffusivities_component_idx,
            const DIRECTION::TYPE& flux_direction,
            const DIRECTION::TYPE& derivative_direction);
        
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
         * boost::shared_ptr to the grid geometry.
         */
        const boost::shared_ptr<geom::CartesianGridGeometry> d_grid_geometry;
        
        /*
         * Number of species.
         */
        const int d_num_species;
        
        /*
         * Number of equations.
         */
        const int d_num_eqn;
        
        /*
         * Whether all or part of global derived cell data related to this class is computed.
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
         * boost::shared_ptr to cell data of diffusivities.
         */
        boost::shared_ptr<pdat::CellData<double> > d_data_diffusivities;
        
        /*
         * boost::weak_ptr to FlowModel.
         */
        boost::weak_ptr<FlowModel> d_flow_model;
        
};

#endif /* FLOW_MODEL_DIFFUSIVE_FLUX_UTILITIES_HPP */
