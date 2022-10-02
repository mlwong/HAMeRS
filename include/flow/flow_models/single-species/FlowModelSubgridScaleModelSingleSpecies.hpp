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
            const HAMERS_SHARED_PTR<tbox::Database>& subgrid_scale_model_db);
        
        ~FlowModelSubgridScaleModelSingleSpecies() {}
        
        /*
         * Return names of different derived variables required to register.
         */
        std::vector<std::string> getDerivedVariablesToRegister() const;
        
        /*
         * Return different derived variables required for interpolation.
         */
        void getDerivedVariablesForInterpolationToSideData(
            std::vector<std::string>& var_to_interpolate,
            std::vector<int>& var_to_interpolate_component_idx) const;
        
        /*
         * Get the variables for the derivatives used at computing subgrid scale diffusivity/viscosity at sides.
         */
        void
        getCellDataOfVariablesForSideDerivativeForSubgridScaleViscosity(
            std::vector<std::string>& derivative_var_data_str,
            std::vector<int>& derivative_var_component_idx,
            const DIRECTION::TYPE& side_direction,
            const DIRECTION::TYPE& derivative_direction);
        
        /*
         * Modify the side data of the diffusivities/viscosities at sides with subgrid scale diffusivity/viscosity.
         */
        void
        updateSideDataOfDiffusiveFluxDiffusivities(
            std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& var_data_for_diffusivities,
            const std::map<DIRECTION::TYPE, std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > >& derivatives,
            const DIRECTION::TYPE& side_direction);
        
        /*
         * Put the characteristics of this class into the restart database.
         */
        void
        putToRestart(
            const HAMERS_SHARED_PTR<tbox::Database>& restart_subgrid_scale_model_db) const;
        
    private:
        /*
         * Constants used to compute subgrid scale viscosity.
         */
        double d_constant_sgs;
        double d_species_Pr_t;
        double d_species_c_p;
};


#endif /* FLOW_MODEL_SUBGRID_SCALE_MODEL_SINGLE_SPECIES_HPP */