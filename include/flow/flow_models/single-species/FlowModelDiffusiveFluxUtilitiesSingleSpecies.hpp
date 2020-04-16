#ifndef FLOW_MODEL_DIFFUSIVE_FLUX_UTILITIES_SINGLE_SPECIES_HPP
#define FLOW_MODEL_DIFFUSIVE_FLUX_UTILITIES_SINGLE_SPECIES_HPP

#include "flow/flow_models/FlowModelDiffusiveFluxUtilities.hpp"

class FlowModelDiffusiveFluxUtilitiesSingleSpecies: public FlowModelDiffusiveFluxUtilities
{
    public:
        FlowModelDiffusiveFluxUtilitiesSingleSpecies(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_species):
                FlowModelDiffusiveFluxUtilities(
                    object_name,
                    dim,
                    grid_geometry,
                    num_species,
                    2 + dim.getValue())
        {}
        
        ~FlowModelDiffusiveFluxUtilitiesSingleSpecies() {}
        
        /*
         * Register the required variables for the computation of diffusive fluxes in the registered patch.
         */
        void
        registerDiffusiveFluxes(
            const hier::IntVector& num_subghosts);
        
        /*
         * Get the variables for the derivatives in the diffusive fluxes.
         */
        void
        getDiffusiveFluxVariablesForDerivative(
            std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > >& derivative_var_data,
            std::vector<std::vector<int> >& derivative_var_component_idx,
            const DIRECTION::TYPE& flux_direction,
            const DIRECTION::TYPE& derivative_direction);
        
        /*
         * Get the diffusivities in the diffusive flux.
         */
        void
        getDiffusiveFluxDiffusivities(
            std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > >& diffusivities_data,
            std::vector<std::vector<int> >& diffusivities_component_idx,
            const DIRECTION::TYPE& flux_direction,
            const DIRECTION::TYPE& derivative_direction);
        
    private:
        
};

#endif /* FLOW_MODEL_DIFFUSIVE_FLUX_UTILITIES_SINGLE_SPECIES_HPP */
