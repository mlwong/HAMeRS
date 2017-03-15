#ifndef FLOW_MODEL_FOUR_EQN_CONSERVATIVE_HPP
#define FLOW_MODEL_FOUR_EQN_CONSERVATIVE_HPP

#include "flow/flow_models/FlowModel.hpp"
#include "flow/flow_models/four-eqn_conservative/Riemann_solvers/RiemannSolverFourEqnConservativeHLLC.hpp"
#include "flow/flow_models/four-eqn_conservative/Riemann_solvers/RiemannSolverFourEqnConservativeHLLC_HLL.hpp"
#include "util/mixing_rules/equations_of_mass_diffusivity/EquationOfMassDiffusivityMixingRulesManager.hpp"
#include "util/mixing_rules/equations_of_shear_viscosity/EquationOfShearViscosityMixingRulesManager.hpp"
#include "util/mixing_rules/equations_of_bulk_viscosity/EquationOfBulkViscosityMixingRulesManager.hpp"
#include "util/mixing_rules/equations_of_thermal_conductivity/EquationOfThermalConductivityMixingRulesManager.hpp"

class FlowModelFourEqnConservative: public FlowModel
{
    public:
        FlowModelFourEqnConservative(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_species,
            const boost::shared_ptr<tbox::Database>& flow_model_db);
        
        ~FlowModelFourEqnConservative() {}
        
        /*
         * Print all characteristics of the flow model class.
         */
        void printClassData(std::ostream& os) const;
        
        /*
         * Put the characteristics of the flow model class into the restart database.
         */
        void putToRestart(
            const boost::shared_ptr<tbox::Database>& restart_db) const;
        
        /*
         * Register the conservative variables.
         */
        void
        registerConservativeVariables(
            RungeKuttaLevelIntegrator* integrator,
            const hier::IntVector& num_ghosts);
        
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
        std::vector<boost::shared_ptr<pdat::CellVariable<double> > >
        getConservativeVariables();
        
        /*
         * Register a patch with a data context.
         */
        void
        registerPatchWithDataContext(
            const hier::Patch& patch,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Register different derived variables in the registered patch. The derived variables to be registered
         * are given as entires in a map of the variable name to the number of sub-ghost cells required.
         * If the variable to be registered is one of the conservative variable, the corresponding entry
         * in the map is ignored.
         */
        void
        registerDerivedCellVariable(
            const std::unordered_map<std::string, hier::IntVector>& num_subghosts_of_data);
        
        /*
         * Register the required derived variables for transformation between conservative
         * variables and characteristic variables.
         */
        void
        registerDerivedVariablesForCharacteristicProjectionOfConservativeVariables(
            const hier::IntVector& num_subghosts,
            const AVERAGING::TYPE& averaging);
        
        /*
         * Register the required derived variables for transformation between primitive variables
         * and characteristic variables.
         */
        void
        registerDerivedVariablesForCharacteristicProjectionOfPrimitiveVariables(
            const hier::IntVector& num_subghosts,
            const AVERAGING::TYPE& averaging);
        
        /*
         * Register the required variables for the computation of diffusive flux in the
         * registered patch.
         */
        void
        registerDiffusiveFlux(
            const hier::IntVector& num_subghosts);
        
        /*
         * Unregister the registered patch. The registered data context and all global derived
         * cell data in the patch are dumped.
         */
        void unregisterPatch();
        
        /*
         * Compute global cell data of different registered derived variables with the registered data context.
         */
        void
        computeGlobalDerivedCellData();
        
        /*
         * Get the global cell data of one cell variable in the registered patch.
         */
        boost::shared_ptr<pdat::CellData<double> >
        getGlobalCellData(const std::string& variable_key);
        
        /*
         * Get the global cell data of different cell variables in the registered patch.
         */
        std::vector<boost::shared_ptr<pdat::CellData<double> > >
        getGlobalCellData(const std::vector<std::string>& variable_keys);
        
        /*
         * Fill the interior global cell data of conservative variables with zeros.
         */
        void
        fillZeroGlobalCellDataConservativeVariables();
        
        /*
         * Update the interior global cell data of conservative variables.
         */
        void
        updateGlobalCellDataConservativeVariables();
        
        /*
         * Get the global cell data of the conservative variables in the registered patch.
         */
        std::vector<boost::shared_ptr<pdat::CellData<double> > >
        getGlobalCellDataConservativeVariables();
        
        /*
         * Get the global cell data of the primitive variables in the registered patch.
         */
        std::vector<boost::shared_ptr<pdat::CellData<double> > >
        getGlobalCellDataPrimitiveVariables();
        
        /*
         * Get the number of projection variables for transformation between conservative
         * variables and characteristic variables.
         */
        int
        getNumberOfProjectionVariablesForConservativeVariables() const;
        
        /*
         * Get the number of projection variables for transformation between primitive variables
         * and characteristic variables.
         */
        int
        getNumberOfProjectionVariablesForPrimitiveVariables() const;
        
        /*
         * Compute global side data of the projection variables for transformation between
         * conservative variables and characteristic variables.
         */
        void
        computeGlobalSideDataProjectionVariablesForConservativeVariables(
            std::vector<boost::shared_ptr<pdat::SideData<double> > >& projection_variables);
        
        /*
         * Compute global side data of the projection variables for transformation between
         * primitive variables and characteristic variables.
         */
        void
        computeGlobalSideDataProjectionVariablesForPrimitiveVariables(
            std::vector<boost::shared_ptr<pdat::SideData<double> > >& projection_variables);
        
        /*
         * Compute global side data of characteristic variables from conservative variables.
         */
        void
        computeGlobalSideDataCharacteristicVariablesFromConservativeVariables(
            std::vector<boost::shared_ptr<pdat::SideData<double> > >& characteristic_variables,
            const std::vector<boost::shared_ptr<pdat::CellData<double> > >& conservative_variables,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& projection_variables,
            const int& idx_offset);
        
        /*
         * Compute global side data of characteristic variables from primitive variables.
         */
        void
        computeGlobalSideDataCharacteristicVariablesFromPrimitiveVariables(
            std::vector<boost::shared_ptr<pdat::SideData<double> > >& characteristic_variables,
            const std::vector<boost::shared_ptr<pdat::CellData<double> > >& primitive_variables,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& projection_variables,
            const int& idx_offset);
        
        /*
         * Compute global side data of conservative variables from characteristic variables.
         */
        void
        computeGlobalSideDataConservativeVariablesFromCharacteristicVariables(
            std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& characteristic_variables,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& projection_variables);
        
        /*
         * Compute global side data of primitive variables from characteristic variables.
         */
        void
        computeGlobalSideDataPrimitiveVariablesFromCharacteristicVariables(
            std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& characteristic_variables,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& projection_variables);
        
        /*
         * Compute the local intercell quantities with conservative variables on each side of the face
         * from Riemann solver at face.
         * flux_face: Convective flux at face.
         * velocity_face: Velocity at face.
         */
        void
        computeLocalFaceFluxAndVelocityFromRiemannSolverWithConservativeVariables(
            std::vector<boost::reference_wrapper<double> >& flux_face,
            std::vector<boost::reference_wrapper<double> >& velocity_face,
            const std::vector<boost::reference_wrapper<double> >& conservative_variables_minus,
            const std::vector<boost::reference_wrapper<double> >& conservative_variables_plus,
            const DIRECTION::TYPE& direction,
            const RIEMANN_SOLVER::TYPE& Riemann_solver);
        
        /*
         * Compute the local intercell quantities with primitive variables on each side of the face
         * from Riemann solver at face.
         * flux_face: Convective flux at face.
         * velocity_face: Velocity at face.
         */
        void
        computeLocalFaceFluxAndVelocityFromRiemannSolverWithPrimitiveVariables(
            std::vector<boost::reference_wrapper<double> >& flux_face,
            std::vector<boost::reference_wrapper<double> >& velocity_face,
            const std::vector<boost::reference_wrapper<double> >& primitive_variables_minus,
            const std::vector<boost::reference_wrapper<double> >& primitive_variables_plus,
            const DIRECTION::TYPE& direction,
            const RIEMANN_SOLVER::TYPE& Riemann_solver);
        
        /*
         * Check whether the given side conservative variables are within the bounds.
         */
        void
        checkGlobalSideDataConservativeVariablesBounded(
            boost::shared_ptr<pdat::SideData<int> >& bounded_flag,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables);
        
        /*
         * Check whether the given side primitive variables are within the bounds.
         */
        void
        checkGlobalSideDataPrimitiveVariablesBounded(
            boost::shared_ptr<pdat::SideData<int> >& bounded_flag,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables);
        
        /*
         * Convert vector of pointers of conservative cell data to vectors of pointers of primitive cell data.
         */
        void
        convertLocalCellDataPointersConservativeVariablesToPrimitiveVariables(
            const std::vector<const double*>& conservative_variables,
            const std::vector<double*>& primitive_variables);
        
        /*
         * Convert vector of pointers of primitive cell data to vectors of pointers of conservative cell data.
         */
        void
        convertLocalCellDataPointersPrimitiveVariablesToConservativeVariables(
            const std::vector<const double*>& primitive_variables,
            const std::vector<double*>& conservative_variables);
        
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
            const boost::shared_ptr<appu::VisItDataWriter>& visit_writer);
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
         * Set the ghost boxes and their dimensions of derived cell variables.
         */
        void
        setGhostBoxesAndDimensionsDerivedCellVariables();
        
        /*
         * Get the global cell data of partial density in the registered patch.
         */
        boost::shared_ptr<pdat::CellData<double> >
        getGlobalCellDataPartialDensity();
        
        /*
         * Get the global cell data of momentum in the registered patch.
         */
        boost::shared_ptr<pdat::CellData<double> >
        getGlobalCellDataMomentum();
        
        /*
         * Get the global cell data of total energy in the registered patch.
         */
        boost::shared_ptr<pdat::CellData<double> >
        getGlobalCellDataTotalEnergy();
        
        /*
         * Compute the global cell data of density in the registered patch.
         */
        void computeGlobalCellDataDensity();
        
        /*
         * Compute the global cell data of mass fraction with density in the registered patch.
         */
        void computeGlobalCellDataMassFractionWithDensity();
        
        /*
         * Compute the global cell data of velocity with density in the registered patch.
         */
        void computeGlobalCellDataVelocityWithDensity();
        
        /*
         * Compute the global cell data of internal energy with density and velocity in the registered
         * patch.
         */
        void computeGlobalCellDataInternalEnergyWithDensityAndVelocity();
        
        /*
         * Compute the global cell data of pressure with density, mass fraction and internal energy in
         * the registered patch.
         */
        void computeGlobalCellDataPressureWithDensityMassFractionAndInternalEnergy();
        
        /*
         * Compute the global cell data of sound speed with density, mass fraction and pressure in the
         * registered patch.
         */
        void computeGlobalCellDataSoundSpeedWithDensityMassFractionAndPressure();
        
        /*
         * Compute the global cell data of temperature with density, mass fraction and pressure in the
         * registered patch.
         */
        void computeGlobalCellDataTemperatureWithDensityMassFractionAndPressure();
        
        /*
         * Compute the global cell data of dilatation with density and velocity in the registered patch.
         */
        void computeGlobalCellDataDilatationWithDensityAndVelocity();
        
        /*
         * Compute the global cell data of vorticity with density and velocity in the registered patch.
         */
        void computeGlobalCellDataVorticityWithDensityAndVelocity();
        
        /*
         * Compute the global cell data of enstrophy with vorticity in the registered patch.
         */
        void computeGlobalCellDataEnstrophyWithVorticity();
        
        /*
         * Compute the global cell data of convective flux with velocity and pressure in the registered
         * patch.
         */
        void computeGlobalCellDataConvectiveFluxWithVelocityAndPressure(
            DIRECTION::TYPE direction);
        
        /*
         * Compute the global cell data of maximum wave speed with velocity and sound speed in the
         * registered patch.
         */
        void computeGlobalCellDataMaxWaveSpeedWithVelocityAndSoundSpeed(
            DIRECTION::TYPE direction);
        
        /*
         * boost::shared_ptr to registered conservative variables.
         */
        boost::shared_ptr<pdat::CellVariable<double> > d_variable_partial_density;
        boost::shared_ptr<pdat::CellVariable<double> > d_variable_momentum;
        boost::shared_ptr<pdat::CellVariable<double> > d_variable_total_energy;
        
        /*
         * Number of sub-ghost cells of derived cell data.
         */
        hier::IntVector d_num_subghosts_density;
        hier::IntVector d_num_subghosts_mass_fraction;
        hier::IntVector d_num_subghosts_velocity;
        hier::IntVector d_num_subghosts_internal_energy;
        hier::IntVector d_num_subghosts_pressure;
        hier::IntVector d_num_subghosts_sound_speed;
        hier::IntVector d_num_subghosts_temperature;
        hier::IntVector d_num_subghosts_dilatation;
        hier::IntVector d_num_subghosts_vorticity;
        hier::IntVector d_num_subghosts_enstrophy;
        hier::IntVector d_num_subghosts_convective_flux_x;
        hier::IntVector d_num_subghosts_convective_flux_y;
        hier::IntVector d_num_subghosts_convective_flux_z;
        hier::IntVector d_num_subghosts_max_wave_speed_x;
        hier::IntVector d_num_subghosts_max_wave_speed_y;
        hier::IntVector d_num_subghosts_max_wave_speed_z;
        hier::IntVector d_num_subghosts_diffusivities;
        
        /*
         * Boxes with sub-ghost cells of derived cell data.
         */
        hier::Box d_subghost_box_density;
        hier::Box d_subghost_box_mass_fraction;
        hier::Box d_subghost_box_velocity;
        hier::Box d_subghost_box_internal_energy;
        hier::Box d_subghost_box_pressure;
        hier::Box d_subghost_box_sound_speed;
        hier::Box d_subghost_box_temperature;
        hier::Box d_subghost_box_dilatation;
        hier::Box d_subghost_box_vorticity;
        hier::Box d_subghost_box_enstrophy;
        hier::Box d_subghost_box_convective_flux_x;
        hier::Box d_subghost_box_convective_flux_y;
        hier::Box d_subghost_box_convective_flux_z;
        hier::Box d_subghost_box_max_wave_speed_x;
        hier::Box d_subghost_box_max_wave_speed_y;
        hier::Box d_subghost_box_max_wave_speed_z;
        hier::Box d_subghost_box_diffusivities;
        
        /*
         * Dimensions of boxes with sub-ghost cells of derived cell data.
         */
        hier::IntVector d_subghostcell_dims_density;
        hier::IntVector d_subghostcell_dims_mass_fraction;
        hier::IntVector d_subghostcell_dims_velocity;
        hier::IntVector d_subghostcell_dims_internal_energy;
        hier::IntVector d_subghostcell_dims_pressure;
        hier::IntVector d_subghostcell_dims_sound_speed;
        hier::IntVector d_subghostcell_dims_temperature;
        hier::IntVector d_subghostcell_dims_dilatation;
        hier::IntVector d_subghostcell_dims_vorticity;
        hier::IntVector d_subghostcell_dims_enstrophy;
        hier::IntVector d_subghostcell_dims_convective_flux_x;
        hier::IntVector d_subghostcell_dims_convective_flux_y;
        hier::IntVector d_subghostcell_dims_convective_flux_z;
        hier::IntVector d_subghostcell_dims_max_wave_speed_x;
        hier::IntVector d_subghostcell_dims_max_wave_speed_y;
        hier::IntVector d_subghostcell_dims_max_wave_speed_z;
        hier::IntVector d_subghostcell_dims_diffusivities;
        
        /*
         * boost::shared_ptr to derived cell data.
         */
        boost::shared_ptr<pdat::CellData<double> > d_data_density;
        boost::shared_ptr<pdat::CellData<double> > d_data_mass_fraction;
        boost::shared_ptr<pdat::CellData<double> > d_data_velocity;
        boost::shared_ptr<pdat::CellData<double> > d_data_internal_energy;
        boost::shared_ptr<pdat::CellData<double> > d_data_pressure;
        boost::shared_ptr<pdat::CellData<double> > d_data_sound_speed;
        boost::shared_ptr<pdat::CellData<double> > d_data_temperature;
        boost::shared_ptr<pdat::CellData<double> > d_data_dilatation;
        boost::shared_ptr<pdat::CellData<double> > d_data_vorticity;
        boost::shared_ptr<pdat::CellData<double> > d_data_enstrophy;
        boost::shared_ptr<pdat::CellData<double> > d_data_convective_flux_x;
        boost::shared_ptr<pdat::CellData<double> > d_data_convective_flux_y;
        boost::shared_ptr<pdat::CellData<double> > d_data_convective_flux_z;
        boost::shared_ptr<pdat::CellData<double> > d_data_max_wave_speed_x;
        boost::shared_ptr<pdat::CellData<double> > d_data_max_wave_speed_y;
        boost::shared_ptr<pdat::CellData<double> > d_data_max_wave_speed_z;
        boost::shared_ptr<pdat::CellData<double> > d_data_diffusivities;
        
        /*
         * boost::shared_ptr to Riemann solvers.
         */
        boost::shared_ptr<RiemannSolverFourEqnConservativeHLLC>     d_Riemann_solver_HLLC;
        boost::shared_ptr<RiemannSolverFourEqnConservativeHLLC_HLL> d_Riemann_solver_HLLC_HLL;
        
        /*
         * Upper and lower bounds on variables.
         */
        double d_Y_bound_lo;
        double d_Y_bound_up;
        
        /*
         * A string variable to describe the equation of mass diffusivity used.
         */
        std::string d_equation_of_mass_diffusivity_str;
        
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
         * boost::shared_ptr to EquationOfMassDiffusivityMixingRules.
         */
        boost::shared_ptr<EquationOfMassDiffusivityMixingRules>
            d_equation_of_mass_diffusivity_mixing_rules;
        
        /*
         * boost::shared_ptr to EquationOfMassDiffusivityMixingRulesManager.
         */
        boost::shared_ptr<EquationOfMassDiffusivityMixingRulesManager>
            d_equation_of_mass_diffusivity_mixing_rules_manager;
        
        /*
         * boost::shared_ptr to EquationOfShearViscosityMixingRules.
         */
        boost::shared_ptr<EquationOfShearViscosityMixingRules>
            d_equation_of_shear_viscosity_mixing_rules;
        
        /*
         * boost::shared_ptr to EquationOfShearViscosityMixingRulesManager.
         */
        boost::shared_ptr<EquationOfShearViscosityMixingRulesManager>
            d_equation_of_shear_viscosity_mixing_rules_manager;
        
        /*
         * boost::shared_ptr to EquationOfBulkViscosityMixingRules.
         */
        boost::shared_ptr<EquationOfBulkViscosityMixingRules>
            d_equation_of_bulk_viscosity_mixing_rules;
        
        /*
         * boost::shared_ptr to EquationOfBulkViscosityMixingRulesManager.
         */
        boost::shared_ptr<EquationOfBulkViscosityMixingRulesManager>
            d_equation_of_bulk_viscosity_mixing_rules_manager;
        
        /*
         * boost::shared_ptr to EquationOfThermalConductivityMixingRules.
         */
        boost::shared_ptr<EquationOfThermalConductivityMixingRules>
            d_equation_of_thermal_conductivity_mixing_rules;
        
        /*
         * boost::shared_ptr to EquationOfThermalConductivityMixingRulesManager.
         */
        boost::shared_ptr<EquationOfThermalConductivityMixingRulesManager>
            d_equation_of_thermal_conductivity_mixing_rules_manager;
        
};

#endif /* FLOW_MODEL_FOUR_EQN_CONSERVATIVE_HPP */
