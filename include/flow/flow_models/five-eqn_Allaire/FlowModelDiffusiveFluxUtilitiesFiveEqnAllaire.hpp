#ifndef FLOW_MODEL_DIFFUSIVE_FLUX_UTILITIES_FIVE_EQN_ALLAIRE_HPP
#define FLOW_MODEL_DIFFUSIVE_FLUX_UTILITIES_FIVE_EQN_ALLAIRE_HPP

#include "flow/flow_models/FlowModelDiffusiveFluxUtilities.hpp"
#include "util/mixing_rules/equations_of_bulk_viscosity/EquationOfBulkViscosityMixingRulesManager.hpp"
#include "util/mixing_rules/equations_of_shear_viscosity/EquationOfShearViscosityMixingRulesManager.hpp"

class FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire: public FlowModelDiffusiveFluxUtilities
{
    public:
        FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_species,
            const boost::shared_ptr<EquationOfShearViscosityMixingRules> equation_of_shear_viscosity_mixing_rules,
            const boost::shared_ptr<EquationOfBulkViscosityMixingRules> equation_of_bulk_viscosity_mixing_rules):
                FlowModelDiffusiveFluxUtilities(
                    object_name,
                    dim,
                    grid_geometry,
                    num_species,
                    dim.getValue() + 2*num_species),
                d_equation_of_shear_viscosity_mixing_rules(equation_of_shear_viscosity_mixing_rules),
                d_equation_of_bulk_viscosity_mixing_rules(equation_of_bulk_viscosity_mixing_rules)
        {}
        
        ~FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire() {}
        
        /*
         * Register the required variables for the computation of diffusive fluxes in the registered patch.
         */
        void
        registerDiffusiveFluxes(
            const hier::IntVector& num_subghosts);
        
        /*
         * The cell data of all derived variables in the patch for this class are dumped.
         */
        void clearData();
        
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
        /*
         * boost::shared_ptr to EquationOfShearViscosityMixingRules.
         */
        const boost::shared_ptr<EquationOfShearViscosityMixingRules>
            d_equation_of_shear_viscosity_mixing_rules;
        
        /*
         * boost::shared_ptr to EquationOfBulkViscosityMixingRules.
         */
        const boost::shared_ptr<EquationOfBulkViscosityMixingRules>
            d_equation_of_bulk_viscosity_mixing_rules;
        
};

#endif /* FLOW_MODEL_DIFFUSIVE_FLUX_UTILITIES_FIVE_EQN_ALLAIRE_HPP */
