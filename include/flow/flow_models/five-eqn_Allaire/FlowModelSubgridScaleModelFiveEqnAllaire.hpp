#ifndef FLOW_MODEL_SUBGRID_SCALE_MODEL_FIVE_EQN_ALLAIRE_HPP
#define FLOW_MODEL_SUBGRID_SCALE_MODEL_FIVE_EQN_ALLAIRE_HPP

#include "flow/flow_models/FlowModelSubgridScaleModel.hpp"

class FlowModelSubgridScaleModelFiveEqnAllaire: public FlowModelSubgridScaleModel
{
    public:
        FlowModelSubgridScaleModelFiveEqnAllaire(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_species,
            const HAMERS_SHARED_PTR<tbox::Database>& subgrid_scale_model_db);
        
        ~FlowModelSubgridScaleModelFiveEqnAllaire() {}
        
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
            const DIRECTION::TYPE& side_direction,
            const hier::Patch& patch);
        
        /*
         * Put the characteristics of this class into the restart database.
         */
        void
        putToRestart(
            const HAMERS_SHARED_PTR<tbox::Database>& restart_subgrid_scale_model_db) const;
        
    private:
        /*
         * Kernal to modify the side data of the diffusivities/viscosities at sides with Vreman's subgrid scale
         * diffusivity/viscosity.
         */
        void updateSideDataOfDiffusiveFluxDiffusivitiesVreman(
            double* mu,
            const double* const rho,
            const double* const ddx_u,
            const double* const ddx_v,
            const double* const ddx_w,
            const double* const ddy_u,
            const double* const ddy_v,
            const double* const ddy_w,
            const double* const ddz_u,
            const double* const ddz_v,
            const double* const ddz_w,
            const hier::IntVector& num_ghosts_diffus,
            const hier::IntVector& num_ghosts_der,
            const hier::IntVector& ghostcell_dims_diffus,
            const hier::IntVector& ghostcell_dims_der,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims,
            const double& delta) const;
        
        /*
         * Constant used to compute subgrid scale viscosity.
         */
        double d_constant_sgs;
};

#endif /* FLOW_MODEL_SUBGRID_SCALE_MODEL_FIVE_EQN_ALLAIRE_HPP */