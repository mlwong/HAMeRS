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
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_species,
            const HAMERS_SHARED_PTR<tbox::Database>& flow_model_db,
            const HAMERS_SHARED_PTR<EquationOfMassDiffusivityMixingRules> equation_of_mass_diffusivity_mixing_rules,
            const HAMERS_SHARED_PTR<EquationOfShearViscosityMixingRules> equation_of_shear_viscosity_mixing_rules,
            const HAMERS_SHARED_PTR<EquationOfBulkViscosityMixingRules> equation_of_bulk_viscosity_mixing_rules,
            const HAMERS_SHARED_PTR<EquationOfThermalConductivityMixingRules> equation_of_thermal_conductivity_mixing_rules);
        
        ~FlowModelDiffusiveFluxUtilitiesFourEqnConservative() {}
        
        /*
         * Register different derived variables related to this class in the registered patch. The
         * derived variables to be registered are given as entries in a map of the variable name to
         * the number of sub-ghost cells required.
         */
        void
        registerDerivedVariables(
            const std::unordered_map<std::string, hier::IntVector>& num_subghosts_of_data);
        
        /*
         * Register the required variables for the computation of diffusive fluxes in the registered patch.
         */
        void
        registerDerivedVariablesForDiffusiveFluxes(
            const hier::IntVector& num_subghosts,
            const bool need_side_diffusivities = false);
        
        /*
         * Allocate memory for cell data of different registered derived variables related to this
         * class in the registered patch.
         */
        void allocateMemoryForDerivedCellData();
        
        /*
         * Allocate memory for side data of the diffusivities.
         */
        void allocateMemoryForSideDataOfDiffusiveFluxDiffusivities();
        
        /*
         * Clear cell and side data of different derived variables related to this class in the registered patch.
         */
        void clearCellAndSideData();
        
        /*
         * Compute cell data of different registered derived variables related to this class.
         */
        void computeDerivedCellData();
        
        /*
         * Get the cell data of one cell variable related to this class in the registered patch.
         */
        HAMERS_SHARED_PTR<pdat::CellData<Real> >
        getCellData(const std::string& variable_key);
        
        /*
         * Get the cell data of different cell variables related to this class in the registered patch.
         */
        std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >
        getCellData(
            const std::vector<std::string>& variable_keys);
        
        /*
         * Get the variables for the derivatives in the diffusive fluxes.
         */
        void
        getCellDataOfDiffusiveFluxVariablesForDerivative(
            std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > >& derivative_var_data,
            std::vector<std::vector<int> >& derivative_var_component_idx,
            const DIRECTION::TYPE& flux_direction,
            const DIRECTION::TYPE& derivative_direction);
        
        /*
         * Get the diffusivities in the diffusive flux.
         */
        void
        getCellDataOfDiffusiveFluxDiffusivities(
            std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > >& diffusivities_data,
            std::vector<std::vector<int> >& diffusivities_component_idx,
            const DIRECTION::TYPE& flux_direction,
            const DIRECTION::TYPE& derivative_direction);
        
        /*
         * Get the cell data that needs interpolation to sides for computing side data of diffusivities in the
         * diffusive flux.
         */
        void
        getCellDataForInterpolationToSideDataForDiffusiveFluxDiffusivities(
            std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& var_data_for_diffusivities,
            std::vector<int>& var_data_for_diffusivities_component_idx);
        
        /*
         * Compute the side data of the diffusivities in the diffusive flux with the interpolated side data.
         */
        void
        computeSideDataOfDiffusiveFluxDiffusivities(
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& var_data_for_diffusivities);
        
        /*
         * Get the side data of the diffusivities in the diffusive fluxa.
         */
        void
        getSideDataOfDiffusiveFluxDiffusivities(
            std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > >& diffusivities_data,
            std::vector<std::vector<int> >& diffusivities_component_idx,
            const DIRECTION::TYPE& flux_direction,
            const DIRECTION::TYPE& derivative_direction);
        
    private:
        /*
         * Set the number of sub-ghost cells of a variable.
         * This function can be called recursively if the variables are computed recursively.
         */
        void
        setNumberOfSubGhosts(
            const hier::IntVector& num_subghosts,
            const std::string& variable_name,
            const std::string& parent_variable_name);
        
        /*
         * Set the ghost boxes of derived cell variables.
         */
        void
        setDerivedCellVariableGhostBoxes();
        
        /*
         * Compute the cell data of mass diffusivities in the registered patch.
         */
        void computeCellDataOfMassDiffusivities();
        
        /*
         * Compute the cell data of shear viscosity in the registered patch.
         */
        void computeCellDataOfShearViscosity();
        
        /*
         * Compute the cell data of bulk viscosity in the registered patch.
         */
        void computeCellDataOfBulkViscosity();
        
        /*
         * Compute the cell data of thermal conductivity in the registered patch.
         */
        void computeCellDataOfThermalConductivity();
        
        /*
         * Compute the cell data of diffusivities in the registered patch.
         */
        void
        computeCellDataOfDiffusivities();
        
        /*
         * Number of sub-ghost cells of derived cell data related to this class.
         */
        hier::IntVector d_num_subghosts_mass_diffusivities;
        hier::IntVector d_num_subghosts_shear_viscosity;
        hier::IntVector d_num_subghosts_bulk_viscosity;
        hier::IntVector d_num_subghosts_thermal_conductivity;
        
        /*
         * Boxes with sub-ghost cells of derived cell data related to this class.
         */
        hier::Box d_subghost_box_mass_diffusivities;
        hier::Box d_subghost_box_shear_viscosity;
        hier::Box d_subghost_box_bulk_viscosity;
        hier::Box d_subghost_box_thermal_conductivity;
        
        /*
         * Dimensions of boxes with sub-ghost cells of derived cell data related to this class.
         */
        hier::IntVector d_subghostcell_dims_mass_diffusivities;
        hier::IntVector d_subghostcell_dims_shear_viscosity;
        hier::IntVector d_subghostcell_dims_bulk_viscosity;
        hier::IntVector d_subghostcell_dims_thermal_conductivity;
        
        /*
         * HAMERS_SHARED_PTR to derived cell data related to this class.
         */
        HAMERS_SHARED_PTR<pdat::CellData<Real> > d_data_mass_diffusivities;
        HAMERS_SHARED_PTR<pdat::CellData<Real> > d_data_shear_viscosity;
        HAMERS_SHARED_PTR<pdat::CellData<Real> > d_data_bulk_viscosity;
        HAMERS_SHARED_PTR<pdat::CellData<Real> > d_data_thermal_conductivity;
        
        /*
         * Whether derived cell data related to this class is computed.
         */
        bool d_cell_data_computed_mass_diffusivities;
        bool d_cell_data_computed_shear_viscosity;
        bool d_cell_data_computed_bulk_viscosity;
        bool d_cell_data_computed_thermal_conductivity;
        
        /*
         * HAMERS_SHARED_PTR to EquationOfMassDiffusivityMixingRules.
         */
        const HAMERS_SHARED_PTR<EquationOfMassDiffusivityMixingRules>
            d_equation_of_mass_diffusivity_mixing_rules;
        
        /*
         * HAMERS_SHARED_PTR to EquationOfShearViscosityMixingRules.
         */
        const HAMERS_SHARED_PTR<EquationOfShearViscosityMixingRules>
            d_equation_of_shear_viscosity_mixing_rules;
        
        /*
         * HAMERS_SHARED_PTR to EquationOfBulkViscosityMixingRules.
         */
        const HAMERS_SHARED_PTR<EquationOfBulkViscosityMixingRules>
            d_equation_of_bulk_viscosity_mixing_rules;
        
        /*
         * HAMERS_SHARED_PTR to EquationOfThermalConductivityMixingRules.
         */
        const HAMERS_SHARED_PTR<EquationOfThermalConductivityMixingRules>
            d_equation_of_thermal_conductivity_mixing_rules;
        
};

#endif /* FLOW_MODEL_DIFFUSIVE_FLUX_UTILITIES_FOUR_EQN_CONSERVATIVE_HPP */
