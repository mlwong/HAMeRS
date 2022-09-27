#ifndef FLOW_MODEL_SUBGRID_SCALE_MODEL_SINGLE_SPECIES_HPP
#define FLOW_MODEL_SUBGRID_SCALE_MODEL_SINGLE_SPECIES_HPP

#include "flow/flow_models/FlowModelSubgridScaleModel.hpp"

class FlowModelSubgridScaleModelSingleSpecies: public FlowModelSubgridScaleModel
{
    public:
        FlowModelSubgridScaleModelSingleSpecies(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_species,
            const HAMERS_SHARED_PTR<tbox::Database>& subgrid_scale_model_db):
                FlowModelSubgridScaleModel(
                    object_name,
                    dim,
                    grid_geometry,
                    num_species,
                    2 + dim.getValue(),
                    subgrid_scale_model_db)
        {}
        
        ~FlowModelSubgridScaleModelSingleSpecies() {}
        
        /*
         * Return names of different derived variables related to this class in the registered patch.
         */
        std::vector<std::string>
        getDerivedVariablesToRegister() const;
        
        /*
         * Get the variables for the derivatives used at computing subgrid scale diffusivity/viscosity at midpoints.
         */
        void
        getCellDataOfVariablesForDerivativeForSubgridScaleViscosity(
            std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& derivative_var_data,
            std::vector<std::vector<int> >& derivative_var_component_idx,
            const DIRECTION::TYPE& midpoint_direction,
            const DIRECTION::TYPE& derivative_direction);
        
        /*
         * Modify the side data of the diffusivities/viscosities with subgrid scale diffusivity/viscosity.
         */
        void
        updateSideDataOfDiffusiveFluxDiffusivities(
            std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& var_data_for_diffusivities,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > >& derivatives_x,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > >& derivatives_y,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > >& derivatives_z);
        
    private:
    
};


#endif /* FLOW_MODEL_SUBGRID_SCALE_MODEL_SINGLE_SPECIES_HPP */