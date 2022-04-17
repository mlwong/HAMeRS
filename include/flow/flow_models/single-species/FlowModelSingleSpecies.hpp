#ifndef FLOW_MODEL_SINGLE_SPECIES_HPP
#define FLOW_MODEL_SINGLE_SPECIES_HPP

#include "flow/flow_models/FlowModel.hpp"
#include "util/mixing_rules/equations_of_shear_viscosity/EquationOfShearViscosityMixingRulesManager.hpp"
#include "util/mixing_rules/equations_of_bulk_viscosity/EquationOfBulkViscosityMixingRulesManager.hpp"
#include "util/mixing_rules/equations_of_thermal_conductivity/EquationOfThermalConductivityMixingRulesManager.hpp"

class FlowModelSingleSpecies: public FlowModel
{
    public:
        FlowModelSingleSpecies(
            const std::string& object_name,
            const std::string& project_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_species,
            const HAMERS_SHARED_PTR<tbox::Database>& flow_model_db);
        
        ~FlowModelSingleSpecies() {}
        
        /*
         * Print all characteristics of the flow model class.
         */
        void printClassData(std::ostream& os) const;
        
        /*
         * Put the characteristics of the flow model class into the restart database.
         */
        void putToRestart(
            const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const;
        
        /*
         * Register the conservative variables.
         */
        void
        registerConservativeVariables(
            RungeKuttaLevelIntegrator* integrator,
            const hier::IntVector& num_ghosts,
            const hier::IntVector& num_ghosts_intermediate);
        
        /*
         * Get the names of conservative variables.
         */
        std::vector<std::string> getNamesOfConservativeVariables(bool have_underscores = false);
        
        /*
         * Get the names of primitive variables.
         */
        std::vector<std::string> getNamesOfPrimitiveVariables(bool have_underscores = false);
        
        /*
         * Get the variable types of conservative variables.
         */
        std::vector<std::string> getVariableTypesOfConservativeVariables();
        
        /*
         * Get the variable types of primitive variables.
         */
        std::vector<std::string> getVariableTypesOfPrimitiveVariables();
        
        /*
         * Get the conservative variables.
         */
        std::vector<HAMERS_SHARED_PTR<pdat::CellVariable<double> > >
        getConservativeVariables();
        
        /*
         * Register a patch with a data context.
         */
        void
        registerPatchWithDataContext(
            const hier::Patch& patch,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Register different derived variables in the registered patch. The derived variables to be registered
         * are given as entires in a map of the variable name to the number of sub-ghost cells required.
         * If the variable to be registered is one of the conservative variable, the corresponding entry
         * in the map is ignored.
         */
        void
        registerDerivedVariables(
            const std::unordered_map<std::string, hier::IntVector>& num_subghosts_of_data);
        
        /*
         * Unregister the registered patch. The registered data context and the cell data of all derived variables in
         * the patch are cleared.
         */
        void unregisterPatch();
        
        /*
         * Allocate memory for cell data of different registered derived variables.
         */
        void allocateMemoryForDerivedCellData();
        
        /*
         * Compute the cell data of different registered derived variables with the registered data context.
         */
        void
        computeDerivedCellData();
        
        /*
         * Get the cell data of one cell variable in the registered patch.
         */
        HAMERS_SHARED_PTR<pdat::CellData<double> >
        getCellData(const std::string& variable_key);
        
        /*
         * Get the cell data of different cell variables in the registered patch.
         */
        std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >
        getCellData(const std::vector<std::string>& variable_keys);
        
        /*
         * Get the cell data of species cell variables in the registered patch.
         */
        std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >
        getSpeciesCellData(const std::string& variable_key);
        
        /*
         * Fill the cell data of conservative variables in the interior box with value zero.
         */
        void
        fillCellDataOfConservativeVariablesWithZero();
        
        /*
         * Update the cell data of conservative variables in the interior box after time advancement.
         */
        void
        updateCellDataOfConservativeVariables();
        
        /*
         * Get the cell data of the conservative variables in the registered patch.
         */
        std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >
        getCellDataOfConservativeVariables();
        
        /*
         * Get the cell data of the primitive variables in the registered patch.
         */
        std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >
        getCellDataOfPrimitiveVariables();
        
        /*
         * Compute derived plot quantities registered with the VisIt data writers from data that
         * is maintained on each patch in the hierarchy.
         */
        bool
        packDerivedDataIntoDoubleBuffer(
            double* buffer,
            const hier::Patch& patch,
            const hier::Box& region,
            const std::string& variable_name,
            int depth_id,
            double simulation_time) const;
        
        /*
         * Register the plotting quantities.
         */
#ifdef HAVE_HDF5
        void
        registerPlotQuantities(
            const HAMERS_SHARED_PTR<ExtendedVisItDataWriter>& visit_writer);
#endif
        
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
         * Get the cell data of density in the registered patch.
         */
        HAMERS_SHARED_PTR<pdat::CellData<double> >
        getCellDataOfDensity();
        
        /*
         * Get the cell data of momentum in the registered patch.
         */
        HAMERS_SHARED_PTR<pdat::CellData<double> >
        getCellDataOfMomentum();
        
        /*
         * Get the cell data of total energy in the registered patch.
         */
        HAMERS_SHARED_PTR<pdat::CellData<double> >
        getCellDataOfTotalEnergy();
        
        /*
         * Compute the cell data of velocity in the registered patch.
         */
        void computeCellDataOfVelocity(
            const hier::Box& domain);
        
        /*
         * Compute the cell data of internal energy with velocity in the registered patch.
         */
        void computeCellDataOfInternalEnergyWithVelocity(
            const hier::Box& domain);
        
        /*
         * Compute the cell data of pressure with internal energy in the registered patch.
         */
        void computeCellDataOfPressureWithInternalEnergy(
            const hier::Box& domain);
        
        /*
         * Compute the cell data of sound speed with pressure in the registered patch.
         */
        void computeCellDataOfSoundSpeedWithPressure(
            const hier::Box& domain);
        
        /*
         * Compute the cell data of temperature with pressure in the registered patch.
         */
        void computeCellDataOfTemperatureWithPressure(
            const hier::Box& domain);
        
        /*
         * Compute the cell data of convective flux with velocity and pressure in the registered patch.
         */
        void computeCellDataOfConvectiveFluxWithVelocityAndPressure(
            const DIRECTION::TYPE& direction,
            const hier::Box& domain);
        
        /*
         * Compute the cell data of maximum wave speed with velocity and sound speed in the registered patch.
         */
        void computeCellDataOfMaxWaveSpeedWithVelocityAndSoundSpeed(
            const DIRECTION::TYPE& direction,
            const hier::Box& domain);
        
        /*
         * Compute the cell data of maximum diffusivity with pressure and temperature in the registered patch.
         */
        void computeCellDataOfMaxDiffusivityWithPressureAndTemperature(
            const hier::Box& domain);
        
        /*
         * HAMERS_SHARED_PTR to registered conservative variables.
         */
        static HAMERS_SHARED_PTR<pdat::CellVariable<double> > s_variable_density;
        static HAMERS_SHARED_PTR<pdat::CellVariable<double> > s_variable_momentum;
        static HAMERS_SHARED_PTR<pdat::CellVariable<double> > s_variable_total_energy;
        
        /*
         * Number of sub-ghost cells of derived cell data.
         */
        hier::IntVector d_num_subghosts_velocity;
        hier::IntVector d_num_subghosts_internal_energy;
        hier::IntVector d_num_subghosts_pressure;
        hier::IntVector d_num_subghosts_sound_speed;
        hier::IntVector d_num_subghosts_temperature;
        hier::IntVector d_num_subghosts_convective_flux_x;
        hier::IntVector d_num_subghosts_convective_flux_y;
        hier::IntVector d_num_subghosts_convective_flux_z;
        hier::IntVector d_num_subghosts_max_wave_speed_x;
        hier::IntVector d_num_subghosts_max_wave_speed_y;
        hier::IntVector d_num_subghosts_max_wave_speed_z;
        hier::IntVector d_num_subghosts_max_diffusivity;
        
        /*
         * Boxes with sub-ghost cells of derived cell data.
         */
        hier::Box d_subghost_box_velocity;
        hier::Box d_subghost_box_internal_energy;
        hier::Box d_subghost_box_pressure;
        hier::Box d_subghost_box_sound_speed;
        hier::Box d_subghost_box_temperature;
        hier::Box d_subghost_box_convective_flux_x;
        hier::Box d_subghost_box_convective_flux_y;
        hier::Box d_subghost_box_convective_flux_z;
        hier::Box d_subghost_box_max_wave_speed_x;
        hier::Box d_subghost_box_max_wave_speed_y;
        hier::Box d_subghost_box_max_wave_speed_z;
        hier::Box d_subghost_box_max_diffusivity;
        
        /*
         * Dimensions of boxes with sub-ghost cells of derived cell data.
         */
        hier::IntVector d_subghostcell_dims_velocity;
        hier::IntVector d_subghostcell_dims_internal_energy;
        hier::IntVector d_subghostcell_dims_pressure;
        hier::IntVector d_subghostcell_dims_sound_speed;
        hier::IntVector d_subghostcell_dims_temperature;
        hier::IntVector d_subghostcell_dims_convective_flux_x;
        hier::IntVector d_subghostcell_dims_convective_flux_y;
        hier::IntVector d_subghostcell_dims_convective_flux_z;
        hier::IntVector d_subghostcell_dims_max_wave_speed_x;
        hier::IntVector d_subghostcell_dims_max_wave_speed_y;
        hier::IntVector d_subghostcell_dims_max_wave_speed_z;
        hier::IntVector d_subghostcell_dims_max_diffusivity;
        
        /*
         * HAMERS_SHARED_PTR to derived cell data.
         */
        HAMERS_SHARED_PTR<pdat::CellData<double> > d_data_velocity;
        HAMERS_SHARED_PTR<pdat::CellData<double> > d_data_internal_energy;
        HAMERS_SHARED_PTR<pdat::CellData<double> > d_data_pressure;
        HAMERS_SHARED_PTR<pdat::CellData<double> > d_data_sound_speed;
        HAMERS_SHARED_PTR<pdat::CellData<double> > d_data_temperature;
        HAMERS_SHARED_PTR<pdat::CellData<double> > d_data_convective_flux_x;
        HAMERS_SHARED_PTR<pdat::CellData<double> > d_data_convective_flux_y;
        HAMERS_SHARED_PTR<pdat::CellData<double> > d_data_convective_flux_z;
        HAMERS_SHARED_PTR<pdat::CellData<double> > d_data_max_wave_speed_x;
        HAMERS_SHARED_PTR<pdat::CellData<double> > d_data_max_wave_speed_y;
        HAMERS_SHARED_PTR<pdat::CellData<double> > d_data_max_wave_speed_z;
        HAMERS_SHARED_PTR<pdat::CellData<double> > d_data_max_diffusivity;
        
        /*
         * Whether derived cell data is computed.
         */
        bool d_cell_data_computed_velocity;
        bool d_cell_data_computed_internal_energy;
        bool d_cell_data_computed_pressure;
        bool d_cell_data_computed_sound_speed;
        bool d_cell_data_computed_temperature;
        bool d_cell_data_computed_convective_flux_x;
        bool d_cell_data_computed_convective_flux_y;
        bool d_cell_data_computed_convective_flux_z;
        bool d_cell_data_computed_max_wave_speed_x;
        bool d_cell_data_computed_max_wave_speed_y;
        bool d_cell_data_computed_max_wave_speed_z;
        bool d_cell_data_computed_max_diffusivity;
        
        /*
         * A string variable to describe the equation of shear viscosity used.
         */
        std::string d_equation_of_shear_viscosity_str;
        
        /*
         * A string variable to describe the equation of bulk viscosity used.
         */
        std::string d_equation_of_bulk_viscosity_str;
        
        /*
         * A string variable to describe the equation of thermal conductivity used.
         */
        std::string d_equation_of_thermal_conductivity_str;
        
        /*
         * HAMERS_SHARED_PTR to EquationOfShearViscosityMixingRules.
         */
        HAMERS_SHARED_PTR<EquationOfShearViscosityMixingRules>
            d_equation_of_shear_viscosity_mixing_rules;
        
        /*
         * HAMERS_SHARED_PTR to EquationOfShearViscosityMixingRulesManager.
         */
        HAMERS_SHARED_PTR<EquationOfShearViscosityMixingRulesManager>
            d_equation_of_shear_viscosity_mixing_rules_manager;
        
        /*
         * HAMERS_SHARED_PTR to EquationOfBulkViscosityMixingRules.
         */
        HAMERS_SHARED_PTR<EquationOfBulkViscosityMixingRules>
            d_equation_of_bulk_viscosity_mixing_rules;
        
        /*
         * HAMERS_SHARED_PTR to EquationOfBulkViscosityMixingRulesManager.
         */
        HAMERS_SHARED_PTR<EquationOfBulkViscosityMixingRulesManager>
            d_equation_of_bulk_viscosity_mixing_rules_manager;
        
        /*
         * HAMERS_SHARED_PTR to EquationOfThermalConductivityMixingRules.
         */
        HAMERS_SHARED_PTR<EquationOfThermalConductivityMixingRules>
            d_equation_of_thermal_conductivity_mixing_rules;
        
        /*
         * HAMERS_SHARED_PTR to EquationOfThermalConductivityMixingRulesManager.
         */
        HAMERS_SHARED_PTR<EquationOfThermalConductivityMixingRulesManager>
            d_equation_of_thermal_conductivity_mixing_rules_manager;
            
        /*
         * Thermodynamic properties of the species.
         */
        std::vector<double> d_thermo_properties;
        
        /*
         * Molecular properties of the species for shear viscosity.
         */
        std::vector<double> d_molecular_properties_shear_viscosity;
        
        /*
         * Molecular properties of the species for bulk viscosity.
         */
        std::vector<double> d_molecular_properties_bulk_viscosity;
        
        /*
         * Molecular properties of the species for thermal conductivity.
         */
        std::vector<double> d_molecular_properties_thermal_conductivity;
        
};

#endif /* FLOW_MODEL_SINGLE_SPECIES_HPP */
