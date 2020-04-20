#ifndef FLOW_MODEL_DIFFUSIVE_FLUX_UTILITIES_SINGLE_SPECIES_HPP
#define FLOW_MODEL_DIFFUSIVE_FLUX_UTILITIES_SINGLE_SPECIES_HPP

#include "flow/flow_models/FlowModelDiffusiveFluxUtilities.hpp"
#include "util/mixing_rules/equations_of_bulk_viscosity/EquationOfBulkViscosityMixingRulesManager.hpp"
#include "util/mixing_rules/equations_of_shear_viscosity/EquationOfShearViscosityMixingRulesManager.hpp"
#include "util/mixing_rules/equations_of_thermal_conductivity/EquationOfThermalConductivityMixingRulesManager.hpp"

class FlowModelDiffusiveFluxUtilitiesSingleSpecies: public FlowModelDiffusiveFluxUtilities
{
    public:
        FlowModelDiffusiveFluxUtilitiesSingleSpecies(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_species,
            const boost::shared_ptr<EquationOfShearViscosityMixingRules> equation_of_shear_viscosity_mixing_rules,
            const boost::shared_ptr<EquationOfBulkViscosityMixingRules> equation_of_bulk_viscosity_mixing_rules,
            const boost::shared_ptr<EquationOfThermalConductivityMixingRules> equation_of_thermal_conductivity_mixing_rules);
        
        ~FlowModelDiffusiveFluxUtilitiesSingleSpecies() {}
        
        /*
         * Register the required variables for the computation of diffusive fluxes in the registered patch.
         */
        void
        registerDerivedVariablesForDiffusiveFluxes(
            const hier::IntVector& num_subghosts);
        
        /*
         * The cell data of all derived variables in the patch for this class are dumped.
         */
        void clearCellData();
        
        /*
         * Get the variables for the derivatives in the diffusive fluxes.
         */
        void
        getCellDataOfDiffusiveFluxVariablesForDerivative(
            std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > >& derivative_var_data,
            std::vector<std::vector<int> >& derivative_var_component_idx,
            const DIRECTION::TYPE& flux_direction,
            const DIRECTION::TYPE& derivative_direction);
        
        /*
         * Get the diffusivities in the diffusive flux.
         */
        void
        getCellDataOfDiffusiveFluxDiffusivities(
            std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > >& diffusivities_data,
            std::vector<std::vector<int> >& diffusivities_component_idx,
            const DIRECTION::TYPE& flux_direction,
            const DIRECTION::TYPE& derivative_direction);
        
    private:
        /*
         * Number of sub-ghost cells of derived cell data for this class.
         */
        hier::IntVector d_num_subghosts_shear_viscosity;
        hier::IntVector d_num_subghosts_bulk_viscosity;
        hier::IntVector d_num_subghosts_thermal_conductivity;
        
        /*
         * Boxes with sub-ghost cells of derived cell data for this class.
         */
        hier::Box d_subghost_box_shear_viscosity;
        hier::Box d_subghost_box_bulk_viscosity;
        hier::Box d_subghost_box_thermal_conductivity;
        
        /*
         * Dimensions of boxes with sub-ghost cells of derived cell data for this class.
         */
        hier::IntVector d_subghostcell_dims_shear_viscosity;
        hier::IntVector d_subghostcell_dims_bulk_viscosity;
        hier::IntVector d_subghostcell_dims_thermal_conductivity;
        
        /*
         * boost::shared_ptr to derived cell data for this class.
         */
        boost::shared_ptr<pdat::CellData<double> > d_data_shear_viscosity;
        boost::shared_ptr<pdat::CellData<double> > d_data_bulk_viscosity;
        boost::shared_ptr<pdat::CellData<double> > d_data_thermal_conductivity;
        
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

#endif /* FLOW_MODEL_DIFFUSIVE_FLUX_UTILITIES_SINGLE_SPECIES_HPP */
