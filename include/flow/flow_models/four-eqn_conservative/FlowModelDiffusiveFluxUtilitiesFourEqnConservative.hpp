#ifndef FLOW_MODEL_DIFFUSIVE_FLUX_UTILITIES_FOUR_EQN_CONSERVATIVE_HPP
#define FLOW_MODEL_DIFFUSIVE_FLUX_UTILITIES_FOUR_EQN_CONSERVATIVE_HPP

#include "flow/flow_models/FlowModelDiffusiveFluxUtilities.hpp"
#include "util/mixing_rules/equations_of_bulk_viscosity/EquationOfBulkViscosityMixingRulesManager.hpp"
#include "util/mixing_rules/equations_of_mass_diffusivity/EquationOfMassDiffusivityMixingRulesManager.hpp"
#include "util/mixing_rules/equations_of_shear_viscosity/EquationOfShearViscosityMixingRulesManager.hpp"
#include "util/mixing_rules/equations_of_thermal_conductivity/EquationOfThermalConductivityMixingRulesManager.hpp"

class FlowModelDiffusiveFluxUtilitiesFourEqnConservative: public FlowModelDiffusiveFluxUtilities
{
    public:
        FlowModelDiffusiveFluxUtilitiesFourEqnConservative(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_species,
            const boost::shared_ptr<EquationOfStateMixingRules> equation_of_state_mixing_rules,
            const boost::shared_ptr<EquationOfMassDiffusivityMixingRules> equation_of_mass_diffusivity_mixing_rules,
            const boost::shared_ptr<EquationOfShearViscosityMixingRules> equation_of_shear_viscosity_mixing_rules,
            const boost::shared_ptr<EquationOfBulkViscosityMixingRules> equation_of_bulk_viscosity_mixing_rules,
            const boost::shared_ptr<EquationOfThermalConductivityMixingRules> equation_of_thermal_conductivity_mixing_rules):
                FlowModelDiffusiveFluxUtilities(
                    object_name,
                    dim,
                    grid_geometry,
                    num_species,
                    num_species + dim.getValue() + 1),
                d_num_subghosts_mass_diffusivities(-hier::IntVector::getOne(d_dim)),
                d_num_subghosts_shear_viscosity(-hier::IntVector::getOne(d_dim)),
                d_num_subghosts_bulk_viscosity(-hier::IntVector::getOne(d_dim)),
                d_num_subghosts_thermal_conductivity(-hier::IntVector::getOne(d_dim)),
                d_subghost_box_mass_diffusivities(hier::Box::getEmptyBox(dim)),
                d_subghost_box_shear_viscosity(hier::Box::getEmptyBox(dim)),
                d_subghost_box_bulk_viscosity(hier::Box::getEmptyBox(dim)),
                d_subghost_box_thermal_conductivity(hier::Box::getEmptyBox(dim)),
                d_subghostcell_dims_mass_diffusivities(hier::IntVector::getZero(d_dim)),
                d_subghostcell_dims_shear_viscosity(hier::IntVector::getZero(d_dim)),
                d_subghostcell_dims_bulk_viscosity(hier::IntVector::getZero(d_dim)),
                d_subghostcell_dims_thermal_conductivity(hier::IntVector::getZero(d_dim)),
                d_equation_of_state_mixing_rules(equation_of_state_mixing_rules),
                d_equation_of_mass_diffusivity_mixing_rules(equation_of_mass_diffusivity_mixing_rules),
                d_equation_of_shear_viscosity_mixing_rules(equation_of_shear_viscosity_mixing_rules),
                d_equation_of_bulk_viscosity_mixing_rules(equation_of_bulk_viscosity_mixing_rules),
                d_equation_of_thermal_conductivity_mixing_rules(equation_of_thermal_conductivity_mixing_rules)
        {}
        
        ~FlowModelDiffusiveFluxUtilitiesFourEqnConservative() {}
        
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
         * Number of sub-ghost cells of derived cell data for this class.
         */
        hier::IntVector d_num_subghosts_mass_diffusivities;
        hier::IntVector d_num_subghosts_shear_viscosity;
        hier::IntVector d_num_subghosts_bulk_viscosity;
        hier::IntVector d_num_subghosts_thermal_conductivity;
        
        /*
         * Boxes with sub-ghost cells of derived cell data for this class.
         */
        hier::Box d_subghost_box_mass_diffusivities;
        hier::Box d_subghost_box_shear_viscosity;
        hier::Box d_subghost_box_bulk_viscosity;
        hier::Box d_subghost_box_thermal_conductivity;
        
        /*
         * Dimensions of boxes with sub-ghost cells of derived cell data for this class.
         */
        hier::IntVector d_subghostcell_dims_mass_diffusivities;
        hier::IntVector d_subghostcell_dims_shear_viscosity;
        hier::IntVector d_subghostcell_dims_bulk_viscosity;
        hier::IntVector d_subghostcell_dims_thermal_conductivity;
        
        /*
         * boost::shared_ptr to derived cell data for this class.
         */
        boost::shared_ptr<pdat::CellData<double> > d_data_mass_diffusivities;
        boost::shared_ptr<pdat::CellData<double> > d_data_shear_viscosity;
        boost::shared_ptr<pdat::CellData<double> > d_data_bulk_viscosity;
        boost::shared_ptr<pdat::CellData<double> > d_data_thermal_conductivity;
        
        /*
         * boost::shared_ptr to EquationOfStateMixingRules.
         */
        const boost::shared_ptr<EquationOfStateMixingRules>
            d_equation_of_state_mixing_rules;
        
        /*
         * boost::shared_ptr to EquationOfMassDiffusivityMixingRules.
         */
        const boost::shared_ptr<EquationOfMassDiffusivityMixingRules>
            d_equation_of_mass_diffusivity_mixing_rules;
        
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
        
        /*
         * boost::shared_ptr to EquationOfThermalConductivityMixingRules.
         */
        const boost::shared_ptr<EquationOfThermalConductivityMixingRules>
            d_equation_of_thermal_conductivity_mixing_rules;
        
};

#endif /* FLOW_MODEL_DIFFUSIVE_FLUX_UTILITIES_FOUR_EQN_CONSERVATIVE_HPP */
