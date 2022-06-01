#include "flow/flow_models/four-eqn_conservative/FlowModelStatisticsUtilitiesFourEqnConservative.hpp"

#include "extn/patch_hierarchies/ExtendedFlattenedHierarchy.hpp"
#include "flow/flow_models/MPI_helpers/FlowModelMPIHelperAverage.hpp"
#include "flow/flow_models/MPI_helpers/FlowModelMPIHelperCorrelation.hpp"
#include "flow/flow_models/MPI_helpers/FlowModelMPIHelperMaxMin.hpp"
#include "util/MPI_helpers/MPIHelperNumberOfCells.hpp"

#include <fstream>

class EnsembleStatisticsRTIRMI: public EnsembleStatistics
{
    public:
        EnsembleStatisticsRTIRMI(const std::string& object_name):
            EnsembleStatistics(
                object_name)
        {
            setVariablesNotComputed();
        }
        
        void setVariablesNotComputed()
        {
            Y_0_avg_computed      = false;
            Y_0_sq_avg_computed   = false;
            Y_0_Y_1_avg_computed  = false;
            rho_avg_computed      = false;
            rho_sq_avg_computed   = false;
            rho_inv_avg_computed  = false;
            u_avg_computed        = false;
            v_avg_computed        = false;
            w_avg_computed        = false;
            u_sq_avg_computed     = false;
            v_sq_avg_computed     = false;
            w_sq_avg_computed     = false;
            rho_u_avg_computed    = false;
            rho_v_avg_computed    = false;
            rho_w_avg_computed    = false;
            rho_u_u_avg_computed  = false;
            rho_v_v_avg_computed  = false;
            rho_w_w_avg_computed  = false;
            c_avg_computed        = false;
            Omega_avg_computed    = false;
            chi_avg_computed      = false;
            mu_avg_computed       = false;
        }
        
        void clearAllData()
        {
            Y_0_avg_realizations.clear();
            Y_0_sq_avg_realizations.clear();
            Y_0_Y_1_avg_realizations.clear();
            rho_avg_realizations.clear();
            rho_sq_avg_realizations.clear();
            rho_inv_avg_realizations.clear();
            u_avg_realizations.clear();
            v_avg_realizations.clear();
            w_avg_realizations.clear();
            u_sq_avg_realizations.clear();
            v_sq_avg_realizations.clear();
            w_sq_avg_realizations.clear();
            rho_u_avg_realizations.clear();
            rho_v_avg_realizations.clear();
            rho_w_avg_realizations.clear();
            rho_u_u_avg_realizations.clear();
            rho_v_v_avg_realizations.clear();
            rho_w_w_avg_realizations.clear();
            c_avg_realizations.clear();
            Omega_avg_realizations.clear();
            chi_avg_realizations.clear();
            mu_avg_realizations.clear();
            
            setVariablesNotComputed();
        }
        
        // Scratch arrays.
        // Number of realizalizations; number of cells.
        std::vector<std::vector<double> > Y_0_avg_realizations;
        std::vector<std::vector<double> > Y_0_sq_avg_realizations;
        std::vector<std::vector<double> > Y_0_Y_1_avg_realizations;
        std::vector<std::vector<double> > rho_avg_realizations;
        std::vector<std::vector<double> > rho_sq_avg_realizations;
        std::vector<std::vector<double> > rho_inv_avg_realizations;
        std::vector<std::vector<double> > u_avg_realizations;
        std::vector<std::vector<double> > v_avg_realizations;
        std::vector<std::vector<double> > w_avg_realizations;
        std::vector<std::vector<double> > u_sq_avg_realizations;
        std::vector<std::vector<double> > v_sq_avg_realizations;
        std::vector<std::vector<double> > w_sq_avg_realizations;
        std::vector<std::vector<double> > rho_u_avg_realizations;
        std::vector<std::vector<double> > rho_v_avg_realizations;
        std::vector<std::vector<double> > rho_w_avg_realizations;
        std::vector<std::vector<double> > rho_u_u_avg_realizations;
        std::vector<std::vector<double> > rho_v_v_avg_realizations;
        std::vector<std::vector<double> > rho_w_w_avg_realizations;
        std::vector<std::vector<double> > c_avg_realizations;
        std::vector<std::vector<double> > Omega_avg_realizations;
        std::vector<std::vector<double> > chi_avg_realizations;
        std::vector<std::vector<double> > mu_avg_realizations;
        
        // Whether the scratch arrays are filled.
        bool Y_0_avg_computed;
        bool Y_0_sq_avg_computed;
        bool Y_0_Y_1_avg_computed;
        bool rho_avg_computed;
        bool rho_sq_avg_computed;
        bool rho_inv_avg_computed;
        bool u_avg_computed;
        bool v_avg_computed;
        bool w_avg_computed;
        bool u_sq_avg_computed;
        bool v_sq_avg_computed;
        bool w_sq_avg_computed;
        bool rho_u_avg_computed;
        bool rho_v_avg_computed;
        bool rho_w_avg_computed;
        bool rho_u_u_avg_computed;
        bool rho_v_v_avg_computed;
        bool rho_w_w_avg_computed;
        bool c_avg_computed;
        bool Omega_avg_computed;
        bool chi_avg_computed;
        bool mu_avg_computed;
        
    private:
        
};


class RTIRMIStatisticsUtilities
{
    public:
        RTIRMIStatisticsUtilities(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_species,
            const HAMERS_WEAK_PTR<FlowModel> flow_model,
            const HAMERS_SHARED_PTR<EquationOfStateMixingRules> equation_of_state_mixing_rules,
            const HAMERS_SHARED_PTR<EquationOfMassDiffusivityMixingRules> equation_of_mass_diffusivity_mixing_rules,
            const HAMERS_SHARED_PTR<EquationOfShearViscosityMixingRules> equation_of_shear_viscosity_mixing_rules,
            const HAMERS_SHARED_PTR<EquationOfBulkViscosityMixingRules> equation_of_bulk_viscosity_mixing_rules,
            const HAMERS_SHARED_PTR<EquationOfThermalConductivityMixingRules> equation_of_thermal_conductivity_mixing_rules,
            const HAMERS_SHARED_PTR<EnsembleStatisticsRTIRMI> ensemble_statistics):
                d_object_name(object_name),
                d_dim(dim),
                d_grid_geometry(grid_geometry),
                d_num_species(num_species),
                d_flow_model(flow_model),
                d_equation_of_state_mixing_rules(equation_of_state_mixing_rules),
                d_equation_of_mass_diffusivity_mixing_rules(equation_of_mass_diffusivity_mixing_rules),
                d_equation_of_shear_viscosity_mixing_rules(equation_of_shear_viscosity_mixing_rules),
                d_equation_of_bulk_viscosity_mixing_rules(equation_of_bulk_viscosity_mixing_rules),
                d_equation_of_thermal_conductivity_mixing_rules(equation_of_thermal_conductivity_mixing_rules),
                d_num_ghosts_derivative(3),
                d_ensemble_statistics(ensemble_statistics)
        {}
        
        /*
         * Compute averaged mass fraction with assumed homogeneity in y-direction (2D) or yz-plane (3D).
         */
        void
        computeAveragedMassFractionWithHomogeneityInYDirectionOrInYZPlane(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute mass fraction variance with assumed homogeneity in y-direction (2D) or yz-plane (3D).
         */
        void
        computeMassFractionVarianceWithHomogeneityInYDirectionOrInYZPlane(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute mass fraction product with assumed homogeneity in y-direction (2D) or yz-plane (3D).
         */
        void
        computeMassFractionProductWithHomogeneityInYDirectionOrInYZPlane(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute averaged density with assumed homogeneity in y-direction (2D) or yz-plane (3D).
         */
        void
        computeAveragedDensityWithHomogeneityInYDirectionOrInYZPlane(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute density variance with assumed homogeneity in y-direction (2D) or yz-plane (3D).
         */
        void
        computeDensityVarianceWithHomogeneityInYDirectionOrInYZPlane(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute averaged specific volume with assumed homogeneity in y-direction (2D) or yz-plane (3D).
         */
        void
        computeAveragedSpecificVolumeWithHomogeneityInYDirectionOrInYZPlane(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute averaged velocity x-component with assumed homogeneity in y-direction (2D) or yz-plane (3D).
         */
        void
        computeAveragedVelocityXWithHomogeneityInYDirectionOrInYZPlane(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute averaged velocity y-component with assumed homogeneity in y-direction (2D) or yz-plane (3D).
         */
        void
        computeAveragedVelocityYWithHomogeneityInYDirectionOrInYZPlane(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute averaged velocity z-component with assumed homogeneity in yz-plane (3D).
         */
        void
        computeAveragedVelocityZWithHomogeneityInYZPlane(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute velocity x-component variance with assumed homogeneity in y-direction (2D) or yz-plane (3D).
         */
        void
        computeVelocityXVarianceWithHomogeneityInYDirectionOrInYZPlane(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute velocity y-component variance with assumed homogeneity in y-direction (2D) or yz-plane (3D).
         */
        void
        computeVelocityYVarianceWithHomogeneityInYDirectionOrInYZPlane(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute velocity z-component variance with assumed homogeneity in yz-plane (3D).
         */
        void
        computeVelocityZVarianceWithHomogeneityInYZPlane(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute averaged momentum x-component with assumed homogeneity in y-direction (2D) or yz-plane (3D).
         */
        void
        computeAveragedMomentumXWithHomogeneityInYDirectionOrInYZPlane(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute averaged momentum y-component with assumed homogeneity in y-direction (2D) or yz-plane (3D).
         */
        void
        computeAveragedMomentumYWithHomogeneityInYDirectionOrInYZPlane(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute averaged momentum z-component with assumed homogeneity in yz-plane (3D).
         */
        void
        computeAveragedMomentumZWithHomogeneityInYZPlane(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute Reynolds normal stress x-component with assumed homogeneity in y-direction (2D) or yz-plane (3D)
         * to a file.
         */
        void
        computeReynoldsNormalStressXWithHomogeneityInYDirectionOrInYZPlane(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute Reynolds normal stress y-component with assumed homogeneity in y-direction (2D) or yz-plane (3D)
         * to a file.
         */
        void
        computeReynoldsNormalStressYWithHomogeneityInYDirectionOrInYZPlane(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute Reynolds normal stress z-component with assumed homogeneity in yz-plane (3D) to a file.
         */
        void
        computeReynoldsNormalStressZWithHomogeneityInYZPlane(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute averaged sound speed with assumed homogeneity in y-direction (2D) or yz-plane (3D).
         */
        void
        computeAveragedSoundSpeedWithHomogeneityInYDirectionOrInYZPlane(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute enstrophy with assumed homogeneity in y-direction (2D) or yz-plane (3D) to a file.
         */
        void
        computeAveragedEnstrophyWithHomogeneityInYDirectionOrInYZPlane(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute scalar dissipation rate with assumed homogeneity in y-direction (2D) or yz-plane (3D) to a file.
         */
        void
        computeAveragedScalarDissipationRateWithHomogeneityInYDirectionOrInYZPlane(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute averaged dynamic shear viscosity with assumed homogeneity in y-direction (2D) or yz-plane (3D) to a file.
         */
        void
        computeAveragedShearViscosityWithHomogeneityInYDirectionOrInYZPlane(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Output spatial profile of ensemble averaged mass fraction with assumed homogeneity in y-direction (2D) or
         * yz-plane (3D) to a file.
         */
        void
        outputSpatialProfileEnsembleAveragedMassFractionWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double output_time) const;
        
        /*
         * Output spatial profile of ensemble variance of mass fraction with assumed homogeneity in y-direction (2D) or
         * yz-plane (3D) to a file.
         */
        void
        outputSpatialProfileEnsembleMassFractionVarianceWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double output_time) const;
        
        /*
         * Output spatial profile of ensemble averaged density with assumed homogeneity in y-direction (2D) or
         * yz-plane (3D) to a file.
         */
        void
        outputSpatialProfileEnsembleAveragedDensityWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double output_time) const;
        
        /*
         * Output spatial profile of ensemble variance of density with assumed homogeneity in y-direction (2D) or
         * yz-plane (3D) to a file.
         */
        void
        outputSpatialProfileEnsembleDensityVarianceWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double output_time) const;
        
        /*
         * Output spatial profile of ensemble averaged specific volume with assumed homogeneity in y-direction (2D) or
         * yz-plane (3D) to a file.
         */
        void
        outputSpatialProfileEnsembleAveragedSpecificVolumeWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double output_time) const;
        
        /*
         * Output spatial profile of ensemble averaged density-specific-volume covariance with assumed homogeneity in
         * y-direction (2D) or yz-plane (3D) to a file.
         */
        void
        outputSpatialProfileEnsembleAveragedDensitySpecificVolumeCovarianceWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double output_time) const;
        
        /*
         * Output spatial profile of ensemble averaged velocity x-component with assumed homogeneity in y-direction (2D)
         * or yz-plane (3D) to a file.
         */
        void
        outputSpatialProfileEnsembleAveragedVelocityXWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double output_time) const;
        
        /*
         * Output spatial profile of ensemble averaged velocity y-component with assumed homogeneity in y-direction (2D)
         * or yz-plane (3D) to a file.
         */
        void
        outputSpatialProfileEnsembleAveragedVelocityYWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double output_time) const;
        
        /*
         * Output spatial profile of ensemble averaged velocity z-component with assumed homogeneity in yz-plane (3D)
         * to a file.
         */
        void
        outputSpatialProfileEnsembleAveragedVelocityZWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double output_time) const;
        
        /*
         * Output spatial profile of ensemble averaged momentum x-component with assumed homogeneity in y-direction (2D)
         * or yz-plane (3D) to a file.
         */
        void
        outputSpatialProfileEnsembleAveragedMomentumXWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double output_time) const;
        
        /*
         * Output spatial profile of ensemble averaged momentum y-component with assumed homogeneity in y-direction (2D)
         * or yz-plane (3D) to a file.
         */
        void
        outputSpatialProfileEnsembleAveragedMomentumYWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double output_time) const;
        
        /*
         * Output spatial profile of ensemble averaged momentum z-component with assumed homogeneity in yz-plane (3D)
         * to a file.
         */
        void
        outputSpatialProfileEnsembleAveragedMomentumZWithHomogeneitynYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double output_time) const;
        
        /*
         * Output spatial profile of ensemble turbulent mass flux in x-direction with assumed homogeneity in
         * y-direction (2D) or yz-plane (3D) to a file.
         */
        void
        outputSpatialProfileEnsembleTurbulentMassFluxXWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double output_time) const;
        
        /*
         * Output spatial profile of ensemble Reynolds normal stress x-component with assumed homogeneity in
         * y-direction (2D) or yz-plane (3D) to a file.
         */
        void
        outputSpatialProfileEnsembleReynoldsNormalStressXWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double output_time) const;
        
        /*
         * Output spatial profile of ensemble Reynolds normal stress y-component with assumed homogeneity in
         * y-direction (2D) or yz-plane (3D) to a file.
         */
        void
        outputSpatialProfileEnsembleReynoldsNormalStressYWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double output_time) const;
        
        /*
         * Output spatial profile of ensemble Reynolds normal stress z-component with assumed homogeneity in
         * yz-plane (3D) to a file.
         */
        void
        outputSpatialProfileEnsembleReynoldsNormalStressZWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double output_time) const;
        
        /*
         * Output spatial profile of ensemble averaged enstrophy with assumed homogeneity in y-direction (2D) or
         * yz-plane (3D) to a file.
         */
        void
        outputSpatialProfileEnsembleAveragedEnstrophyWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double output_time) const;
        
        /*
         * Output spatial profile of ensemble averaged scalar dissipation rate with assumed homogeneity in
         * y-direction (2D) or yz-plane (3D) to a file.
         */
        void
        outputSpatialProfileEnsembleAveragedScalarDissipationRateWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double output_time) const;
        
        /// Non-spatial profiles.
        
        /*
         * Output ensemble mixing width in x-direction to a file.
         */
        void
        outputEnsembleMixingWidthInXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output ensemble minimum interface location in x-direction to a file.
         */
        void
        outputEnsembleInterfaceMinInXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output ensemble maximum interface location in x-direction to a file.
         */
        void
        outputEnsembleInterfaceMaxInXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output ensemble mixedness in x-direction to a file.
         */
        void
        outputEnsembleMixednessInXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output ensemble mean velocity associated with turbulent mass flux component in x-direction with
         * assumed homogeneity in y-direction (2D) or yz-plane (3D) to a file.
         */
        void
        outputEnsembleTurbulentMassFluxVelocityXWithInhomogeneousXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output ensemble mean density-specific-volume covariance with assumed homogeneity in y-direction (2D)
         * or yz-plane (3D) to a file.
         */
        void
        outputEnsembleDensitySpecificVolumeCovarianceWithInhomogeneousXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output ensemble TKE integrated with assumed homogeneity in y-direction (2D) or yz-plane (3D) to a file.
         */
        void
        outputEnsembleTKEIntegrateWithInhomogeneousXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output ensemble Reynolds normal stress component in x-direction with assumed homogeneity in y-direction (2D)
         * or yz-plane (3D) to a file.
         */
        void
        outputEnsembleReynoldsNormalStressXWithInhomogeneousXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output ensemble Reynolds normal stress component in y-direction with assumed homogeneity in y-direction (2D)
         * or yz-plane (3D) to a file.
         */
        void
        outputEnsembleReynoldsNormalStressYWithInhomogeneousXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output ensemble Reynolds normal stress component in z-direction with assumed homogeneity in yz-plane (3D)
         * to a file.
         */
        void
        outputEnsembleReynoldsNormalStressZWithInhomogeneousXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output ensemble enstrophy integrated to a file.
         */
        void
        outputEnsembleEnstrophyIntegrated(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output ensemble scalar dissipation rate of first species integrated to a file.
         */
        void
        outputEnsembleScalarDissipationRateIntegrated(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output turbulent Reynolds number based on mixing width with assumed homogeneity in y-direction (2D)
         * or yz-plane (3D) to a file.to a file.
         */
        void
        outputEnsembleMixingWidthReynoldsNumberWithInhomogeneousXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output turbulent Mach number with assumed homogeneity in y-direction (2D) or yz-plane (3D) to a file.to a file.
         */
        void
        outputEnsembleTurbulentMachNumberWithInhomogeneousXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Store ensemble statistics.
         */
        HAMERS_SHARED_PTR<EnsembleStatisticsRTIRMI> d_ensemble_statistics;
        
    private:
        /*
         * The object name is used for error/warning reporting.
         */
        const std::string d_object_name;
        
        /*
         * Problem dimension.
         */
        const tbox::Dimension d_dim;
        
        /*
         * HAMERS_SHARED_PTR to the grid geometry.
         */
        const HAMERS_SHARED_PTR<geom::CartesianGridGeometry> d_grid_geometry;
        
        /*
         * Number of species.
         */
        const int d_num_species;
        
        /*
         * HAMERS_WEAK_PTR to FlowModel.
         */
        const HAMERS_WEAK_PTR<FlowModel> d_flow_model;
        
        /*
         * HAMERS_SHARED_PTR to EquationOfStateMixingRules.
         */
        const HAMERS_SHARED_PTR<EquationOfStateMixingRules>
            d_equation_of_state_mixing_rules;
        
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
        
        /*
         * Number of ghost cells to use in taking derivatives.
         */
        const int d_num_ghosts_derivative;
        
};


/*
 * Compute averaged mass fraction with assumed homogeneity in y-direction (2D) or yz-plane (3D).
 */
void
RTIRMIStatisticsUtilities::computeAveragedMassFractionWithHomogeneityInYDirectionOrInYZPlane(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'MASS_FRACTION_AVG_SP' can be computed with two species only."
            << std::endl);
    }
    
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
    FlowModelMPIHelperAverage MPI_helper_average = FlowModelMPIHelperAverage(
        "MPI_helper_average",
        d_dim,
        d_grid_geometry,
        patch_hierarchy,
        flow_model_tmp);
    
    std::vector<double> Y_0_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        "MASS_FRACTIONS",
        0,
        data_context);
    
    std::vector<std::vector<double> >& Y_0_avg_realizations = d_ensemble_statistics->Y_0_avg_realizations;
    Y_0_avg_realizations.push_back(Y_0_avg_global);
    
    d_ensemble_statistics->Y_0_avg_computed = true;
}


/*
 * Compute mass fraction variance with assumed homogeneity in y-direction (2D) or yz-plane (3D).
 */
void
RTIRMIStatisticsUtilities::computeMassFractionVarianceWithHomogeneityInYDirectionOrInYZPlane(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'MASS_FRACTION_VAR_SP' can be computed with two species only."
            << std::endl);
    }
    
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
    FlowModelMPIHelperAverage MPI_helper_average = FlowModelMPIHelperAverage(
        "MPI_helper_average",
        d_dim,
        d_grid_geometry,
        patch_hierarchy,
        flow_model_tmp);
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    
    quantity_names.push_back("MASS_FRACTIONS");
    component_indices.push_back(0);
    
    quantity_names.push_back("MASS_FRACTIONS");
    component_indices.push_back(0);
    
    std::vector<double> Y_0_sq_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    
    std::vector<std::vector<double> >& Y_0_sq_avg_realizations = d_ensemble_statistics->Y_0_sq_avg_realizations;
    Y_0_sq_avg_realizations.push_back(Y_0_sq_avg_global);
    
    d_ensemble_statistics->Y_0_sq_avg_computed = true;
}


/*
 * Compute mass fraction product with assumed homogeneity in y-direction (2D) or yz-plane (3D).
 */
void
RTIRMIStatisticsUtilities::computeMassFractionProductWithHomogeneityInYDirectionOrInYZPlane(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'MASS_FRACTION_PROD_SP' can be computed with two species only."
            << std::endl);
    }
    
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
    FlowModelMPIHelperAverage MPI_helper_average = FlowModelMPIHelperAverage(
        "MPI_helper_average",
        d_dim,
        d_grid_geometry,
        patch_hierarchy,
        flow_model_tmp);
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    
    quantity_names.push_back("MASS_FRACTIONS");
    component_indices.push_back(0);
    
    quantity_names.push_back("MASS_FRACTIONS");
    component_indices.push_back(1);
    
    std::vector<double> Y_0_Y_1_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    
    std::vector<std::vector<double> >& Y_0_Y_1_avg_realizations = d_ensemble_statistics->Y_0_Y_1_avg_realizations;
    Y_0_Y_1_avg_realizations.push_back(Y_0_Y_1_avg_global);
    
    d_ensemble_statistics->Y_0_Y_1_avg_computed = true;
}


/*
 * Compute averaged density with assumed homogeneity in y-direction (2D) or yz-plane (3D).
 */
void
RTIRMIStatisticsUtilities::computeAveragedDensityWithHomogeneityInYDirectionOrInYZPlane(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
    FlowModelMPIHelperAverage MPI_helper_average = FlowModelMPIHelperAverage(
        "MPI_helper_average",
        d_dim,
        d_grid_geometry,
        patch_hierarchy,
        flow_model_tmp);
    
    std::vector<double> rho_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        data_context);
    
    std::vector<std::vector<double> >& rho_avg_realizations = d_ensemble_statistics->rho_avg_realizations;
    rho_avg_realizations.push_back(rho_avg_global);
    
    d_ensemble_statistics->rho_avg_computed = true;
}


/*
 * Compute density variance with assumed homogeneity in y-direction (2D) or yz-plane (3D).
 */
void
RTIRMIStatisticsUtilities::computeDensityVarianceWithHomogeneityInYDirectionOrInYZPlane(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
    FlowModelMPIHelperAverage MPI_helper_average = FlowModelMPIHelperAverage(
        "MPI_helper_average",
        d_dim,
        d_grid_geometry,
        patch_hierarchy,
        flow_model_tmp);
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    
    std::vector<double> rho_sq_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    
    std::vector<std::vector<double> >& rho_sq_avg_realizations = d_ensemble_statistics->rho_sq_avg_realizations;
    rho_sq_avg_realizations.push_back(rho_sq_avg_global);
    
    d_ensemble_statistics->rho_sq_avg_computed = true;
}


/*
 * Compute averaged specific volume with assumed homogeneity in y-direction (2D) or yz-plane (3D).
 */
void
RTIRMIStatisticsUtilities::computeAveragedSpecificVolumeWithHomogeneityInYDirectionOrInYZPlane(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
    FlowModelMPIHelperAverage MPI_helper_average = FlowModelMPIHelperAverage(
        "MPI_helper_average",
        d_dim,
        d_grid_geometry,
        patch_hierarchy,
        flow_model_tmp);
    
    std::vector<double> rho_inv_avg_global = MPI_helper_average.getAveragedReciprocalOfQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        data_context);
    
    std::vector<std::vector<double> >& rho_inv_avg_realizations = d_ensemble_statistics->rho_inv_avg_realizations;
    rho_inv_avg_realizations.push_back(rho_inv_avg_global);
    
    d_ensemble_statistics->rho_inv_avg_computed = true;
}


/*
 * Compute averaged velocity x-component with assumed homogeneity in y-direction (2D) or yz-plane (3D).
 */
void
RTIRMIStatisticsUtilities::computeAveragedVelocityXWithHomogeneityInYDirectionOrInYZPlane(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
    FlowModelMPIHelperAverage MPI_helper_average = FlowModelMPIHelperAverage(
        "MPI_helper_average",
        d_dim,
        d_grid_geometry,
        patch_hierarchy,
        flow_model_tmp);
    
    std::vector<double> u_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        "VELOCITY",
        0,
        data_context);
    
    std::vector<std::vector<double> >& u_avg_realizations = d_ensemble_statistics->u_avg_realizations;
    u_avg_realizations.push_back(u_avg_global);
    
    d_ensemble_statistics->u_avg_computed = true;
}


/*
 * Compute averaged velocity y-component with assumed homogeneity in y-direction (2D) or yz-plane (3D).
 */
void
RTIRMIStatisticsUtilities::computeAveragedVelocityYWithHomogeneityInYDirectionOrInYZPlane(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
    FlowModelMPIHelperAverage MPI_helper_average = FlowModelMPIHelperAverage(
        "MPI_helper_average",
        d_dim,
        d_grid_geometry,
        patch_hierarchy,
        flow_model_tmp);
    
    std::vector<double> v_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        "VELOCITY",
        1,
        data_context);
    
    std::vector<std::vector<double> >& v_avg_realizations = d_ensemble_statistics->v_avg_realizations;
    v_avg_realizations.push_back(v_avg_global);
    
    d_ensemble_statistics->v_avg_computed = true;
}


/*
 * Compute averaged velocity z-component with assumed homogeneity in yz-plane (3D).
 */
void
RTIRMIStatisticsUtilities::computeAveragedVelocityZWithHomogeneityInYZPlane(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
    FlowModelMPIHelperAverage MPI_helper_average = FlowModelMPIHelperAverage(
        "MPI_helper_average",
        d_dim,
        d_grid_geometry,
        patch_hierarchy,
        flow_model_tmp);
    
    std::vector<double> w_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        "VELOCITY",
        2,
        data_context);
    
    std::vector<std::vector<double> >& w_avg_realizations = d_ensemble_statistics->w_avg_realizations;
    w_avg_realizations.push_back(w_avg_global);
    
    d_ensemble_statistics->w_avg_computed = true;
}


/*
 * Compute velocity x-component variance with assumed homogeneity in y-direction (2D) or yz-plane (3D).
 */
void
RTIRMIStatisticsUtilities::computeVelocityXVarianceWithHomogeneityInYDirectionOrInYZPlane(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
    FlowModelMPIHelperAverage MPI_helper_average = FlowModelMPIHelperAverage(
        "MPI_helper_average",
        d_dim,
        d_grid_geometry,
        patch_hierarchy,
        flow_model_tmp);
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    
    std::vector<double> u_sq_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    
    std::vector<std::vector<double> >& u_sq_avg_realizations = d_ensemble_statistics->u_sq_avg_realizations;
    u_sq_avg_realizations.push_back(u_sq_avg_global);
    
    d_ensemble_statistics->u_sq_avg_computed = true;
}


/*
 * Compute velocity y-component variance with assumed homogeneity in y-direction (2D) or yz-plane (3D).
 */
void
RTIRMIStatisticsUtilities::computeVelocityYVarianceWithHomogeneityInYDirectionOrInYZPlane(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
    FlowModelMPIHelperAverage MPI_helper_average = FlowModelMPIHelperAverage(
        "MPI_helper_average",
        d_dim,
        d_grid_geometry,
        patch_hierarchy,
        flow_model_tmp);
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(1);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(1);
    
    std::vector<double> v_sq_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    
    std::vector<std::vector<double> >& v_sq_avg_realizations = d_ensemble_statistics->v_sq_avg_realizations;
    v_sq_avg_realizations.push_back(v_sq_avg_global);
    
    d_ensemble_statistics->v_sq_avg_computed = true;
}


/*
 * Compute velocity z-component variance with assumed homogeneity in yz-plane (3D).
 */
void
RTIRMIStatisticsUtilities::computeVelocityZVarianceWithHomogeneityInYZPlane(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
    FlowModelMPIHelperAverage MPI_helper_average = FlowModelMPIHelperAverage(
        "MPI_helper_average",
        d_dim,
        d_grid_geometry,
        patch_hierarchy,
        flow_model_tmp);
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(2);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(2);
    
    std::vector<double> w_sq_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    
    std::vector<std::vector<double> >& w_sq_avg_realizations = d_ensemble_statistics->w_sq_avg_realizations;
    w_sq_avg_realizations.push_back(w_sq_avg_global);
    
    d_ensemble_statistics->w_sq_avg_computed = true;
}


/*
 * Compute averaged momentum x-component with assumed homogeneity in y-direction (2D) or yz-plane (3D).
 */
void
RTIRMIStatisticsUtilities::computeAveragedMomentumXWithHomogeneityInYDirectionOrInYZPlane(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
    FlowModelMPIHelperAverage MPI_helper_average = FlowModelMPIHelperAverage(
        "MPI_helper_average",
        d_dim,
        d_grid_geometry,
        patch_hierarchy,
        flow_model_tmp);
    
    std::vector<double> rho_u_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        "MOMENTUM",
        0,
        data_context);
    
    std::vector<std::vector<double> >& rho_u_avg_realizations = d_ensemble_statistics->rho_u_avg_realizations;
    rho_u_avg_realizations.push_back(rho_u_avg_global);
    
    d_ensemble_statistics->rho_u_avg_computed = true;
}


/*
 * Compute averaged momentum y-component with assumed homogeneity in y-direction (2D) or yz-plane (3D).
 */
void
RTIRMIStatisticsUtilities::computeAveragedMomentumYWithHomogeneityInYDirectionOrInYZPlane(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
    FlowModelMPIHelperAverage MPI_helper_average = FlowModelMPIHelperAverage(
        "MPI_helper_average",
        d_dim,
        d_grid_geometry,
        patch_hierarchy,
        flow_model_tmp);
    
    std::vector<double> rho_v_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        "MOMENTUM",
        1,
        data_context);
    
    std::vector<std::vector<double> >& rho_v_avg_realizations = d_ensemble_statistics->rho_v_avg_realizations;
    rho_v_avg_realizations.push_back(rho_v_avg_global);
    
    d_ensemble_statistics->rho_v_avg_computed = true;
}


/*
 * Compute averaged momentum z-component with assumed homogeneity in yz-plane (3D).
 */
void
RTIRMIStatisticsUtilities::computeAveragedMomentumZWithHomogeneityInYZPlane(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
    FlowModelMPIHelperAverage MPI_helper_average = FlowModelMPIHelperAverage(
        "MPI_helper_average",
        d_dim,
        d_grid_geometry,
        patch_hierarchy,
        flow_model_tmp);
    
    std::vector<double> rho_w_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        "MOMENTUM",
        2,
        data_context);
    
    std::vector<std::vector<double> >& rho_w_avg_realizations = d_ensemble_statistics->rho_w_avg_realizations;
    rho_w_avg_realizations.push_back(rho_w_avg_global);
    
    d_ensemble_statistics->rho_w_avg_computed = true;
}


/*
 * Compute Reynolds normal stress x-component with assumed homogeneity in y-direction (2D) or yz-plane (3D)
 * to a file.
 */
void
RTIRMIStatisticsUtilities::computeReynoldsNormalStressXWithHomogeneityInYDirectionOrInYZPlane(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
    FlowModelMPIHelperAverage MPI_helper_average = FlowModelMPIHelperAverage(
        "MPI_helper_average",
        d_dim,
        d_grid_geometry,
        patch_hierarchy,
        flow_model_tmp);
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    
    std::vector<double> rho_u_u_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    
    std::vector<std::vector<double> >& rho_u_u_avg_realizations = d_ensemble_statistics->rho_u_u_avg_realizations;
    rho_u_u_avg_realizations.push_back(rho_u_u_avg_global);
    
    d_ensemble_statistics->rho_u_u_avg_computed = true;
}


/*
 * Compute Reynolds normal stress y-component with assumed homogeneity in y-direction (2D) or yz-plane (3D)
 * to a file.
 */
void
RTIRMIStatisticsUtilities::computeReynoldsNormalStressYWithHomogeneityInYDirectionOrInYZPlane(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
    FlowModelMPIHelperAverage MPI_helper_average = FlowModelMPIHelperAverage(
        "MPI_helper_average",
        d_dim,
        d_grid_geometry,
        patch_hierarchy,
        flow_model_tmp);
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(1);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(1);
    
    std::vector<double> rho_v_v_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    
    std::vector<std::vector<double> >& rho_v_v_avg_realizations = d_ensemble_statistics->rho_v_v_avg_realizations;
    rho_v_v_avg_realizations.push_back(rho_v_v_avg_global);
    
    d_ensemble_statistics->rho_v_v_avg_computed = true;
}


/*
 * Compute Reynolds normal stress z-component with assumed homogeneity in yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::computeReynoldsNormalStressZWithHomogeneityInYZPlane(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
    FlowModelMPIHelperAverage MPI_helper_average = FlowModelMPIHelperAverage(
        "MPI_helper_average",
        d_dim,
        d_grid_geometry,
        patch_hierarchy,
        flow_model_tmp);
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(2);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(2);
    
    std::vector<double> rho_w_w_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    
    std::vector<std::vector<double> >& rho_w_w_avg_realizations = d_ensemble_statistics->rho_w_w_avg_realizations;
    rho_w_w_avg_realizations.push_back(rho_w_w_avg_global);
    
    d_ensemble_statistics->rho_w_w_avg_computed = true;
}


/*
 * Compute averaged sound speed with assumed homogeneity in y-direction (2D) or yz-plane (3D).
 */
void
RTIRMIStatisticsUtilities::computeAveragedSoundSpeedWithHomogeneityInYDirectionOrInYZPlane(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
    FlowModelMPIHelperAverage MPI_helper_average = FlowModelMPIHelperAverage(
        "MPI_helper_average",
        d_dim,
        d_grid_geometry,
        patch_hierarchy,
        flow_model_tmp);
    
    std::vector<double> c_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        "SOUND_SPEED",
        0,
        data_context);
    
    std::vector<std::vector<double> >& c_avg_realizations = d_ensemble_statistics->c_avg_realizations;
    c_avg_realizations.push_back(c_avg_global);
    
    d_ensemble_statistics->c_avg_computed = true;
}


/*
 * Compute enstrophy with assumed homogeneity in y-direction (2D) or yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::computeAveragedEnstrophyWithHomogeneityInYDirectionOrInYZPlane(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
    FlowModelMPIHelperAverage MPI_helper_average = FlowModelMPIHelperAverage(
        "MPI_helper_average",
        d_dim,
        d_grid_geometry,
        patch_hierarchy,
        flow_model_tmp);
    
    const hier::IntVector& finest_level_dims = MPI_helper_average.getFinestRefinedDomainNumberOfPoints();
    
    std::vector<std::vector<double> >& Omega_avg_realizations = d_ensemble_statistics->Omega_avg_realizations;
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no 'ENSTROPHY_AVG_SP' for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        std::vector<std::string> quantity_names;
        std::vector<int> component_indices;
        std::vector<bool> use_derivative;
        std::vector<int> derivative_directions;
        
        quantity_names.push_back("DENSITY");
        component_indices.push_back(0);
        use_derivative.push_back(false);
        derivative_directions.push_back(-1);
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(1);
        use_derivative.push_back(true);
        derivative_directions.push_back(0);
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(1);
        use_derivative.push_back(true);
        derivative_directions.push_back(0);
        
        std::vector<double> enstrophy_part_1 = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
            quantity_names,
            component_indices,
            use_derivative,
            derivative_directions,
            d_num_ghosts_derivative,
            data_context);
        
        quantity_names.clear();
        component_indices.clear();
        use_derivative.clear();
        derivative_directions.clear();
        
        quantity_names.push_back("DENSITY");
        component_indices.push_back(0);
        use_derivative.push_back(false);
        derivative_directions.push_back(-1);
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(0);
        use_derivative.push_back(true);
        derivative_directions.push_back(1);
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(0);
        use_derivative.push_back(true);
        derivative_directions.push_back(1);
        
        std::vector<double> enstrophy_part_2 = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
            quantity_names,
            component_indices,
            use_derivative,
            derivative_directions,
            d_num_ghosts_derivative,
            data_context);
        
        quantity_names.clear();
        component_indices.clear();
        use_derivative.clear();
        derivative_directions.clear();
        
        quantity_names.push_back("DENSITY");
        component_indices.push_back(0);
        use_derivative.push_back(false);
        derivative_directions.push_back(-1);
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(1);
        use_derivative.push_back(true);
        derivative_directions.push_back(0);
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(0);
        use_derivative.push_back(true);
        derivative_directions.push_back(1);
        
        std::vector<double> enstrophy_part_3 = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
            quantity_names,
            component_indices,
            use_derivative,
            derivative_directions,
            d_num_ghosts_derivative,
            data_context);
        
        quantity_names.clear();
        component_indices.clear();
        use_derivative.clear();
        derivative_directions.clear();
        
        std::vector<double> enstrophy_full(finest_level_dims[0], double(0));
        for (int i = 0; i < finest_level_dims[0]; i++)
        {
            enstrophy_full[i] = enstrophy_part_1[i] + enstrophy_part_2[i] - double(2)*enstrophy_part_3[i];
        }
        
        Omega_avg_realizations.push_back(enstrophy_full);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        std::vector<std::string> quantity_names;
        std::vector<int> component_indices;
        std::vector<bool> use_derivative;
        std::vector<int> derivative_directions;
        
        quantity_names.push_back("DENSITY");
        component_indices.push_back(0);
        use_derivative.push_back(false);
        derivative_directions.push_back(-1);
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(2);
        use_derivative.push_back(true);
        derivative_directions.push_back(1);
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(2);
        use_derivative.push_back(true);
        derivative_directions.push_back(1);
        
        std::vector<double> enstrophy_part_1 = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
            quantity_names,
            component_indices,
            use_derivative,
            derivative_directions,
            d_num_ghosts_derivative,
            data_context);
        
        quantity_names.clear();
        component_indices.clear();
        use_derivative.clear();
        derivative_directions.clear();
        
        quantity_names.push_back("DENSITY");
        component_indices.push_back(0);
        use_derivative.push_back(false);
        derivative_directions.push_back(-1);
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(1);
        use_derivative.push_back(true);
        derivative_directions.push_back(2);
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(1);
        use_derivative.push_back(true);
        derivative_directions.push_back(2);
        
        std::vector<double> enstrophy_part_2 = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
            quantity_names,
            component_indices,
            use_derivative,
            derivative_directions,
            d_num_ghosts_derivative,
            data_context);
        
        quantity_names.clear();
        component_indices.clear();
        use_derivative.clear();
        derivative_directions.clear();
        
        quantity_names.push_back("DENSITY");
        component_indices.push_back(0);
        use_derivative.push_back(false);
        derivative_directions.push_back(-1);
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(2);
        use_derivative.push_back(true);
        derivative_directions.push_back(1);
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(1);
        use_derivative.push_back(true);
        derivative_directions.push_back(2);
        
        std::vector<double> enstrophy_part_3 = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
            quantity_names,
            component_indices,
            use_derivative,
            derivative_directions,
            d_num_ghosts_derivative,
            data_context);
        
        quantity_names.clear();
        component_indices.clear();
        use_derivative.clear();
        derivative_directions.clear();
        
        quantity_names.push_back("DENSITY");
        component_indices.push_back(0);
        use_derivative.push_back(false);
        derivative_directions.push_back(-1);
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(0);
        use_derivative.push_back(true);
        derivative_directions.push_back(2);
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(0);
        use_derivative.push_back(true);
        derivative_directions.push_back(2);
        
        std::vector<double> enstrophy_part_4 = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
            quantity_names,
            component_indices,
            use_derivative,
            derivative_directions,
            d_num_ghosts_derivative,
            data_context);
        
        quantity_names.clear();
        component_indices.clear();
        use_derivative.clear();
        derivative_directions.clear();
        
        quantity_names.push_back("DENSITY");
        component_indices.push_back(0);
        use_derivative.push_back(false);
        derivative_directions.push_back(-1);
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(2);
        use_derivative.push_back(true);
        derivative_directions.push_back(0);
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(2);
        use_derivative.push_back(true);
        derivative_directions.push_back(0);
        
        std::vector<double> enstrophy_part_5 = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
            quantity_names,
            component_indices,
            use_derivative,
            derivative_directions,
            d_num_ghosts_derivative,
            data_context);
        
        quantity_names.clear();
        component_indices.clear();
        use_derivative.clear();
        derivative_directions.clear();
        
        quantity_names.push_back("DENSITY");
        component_indices.push_back(0);
        use_derivative.push_back(false);
        derivative_directions.push_back(-1);
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(0);
        use_derivative.push_back(true);
        derivative_directions.push_back(2);
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(2);
        use_derivative.push_back(true);
        derivative_directions.push_back(0);
        
        std::vector<double> enstrophy_part_6 = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
            quantity_names,
            component_indices,
            use_derivative,
            derivative_directions,
            d_num_ghosts_derivative,
            data_context);
        
        quantity_names.clear();
        component_indices.clear();
        use_derivative.clear();
        derivative_directions.clear();
        
        quantity_names.push_back("DENSITY");
        component_indices.push_back(0);
        use_derivative.push_back(false);
        derivative_directions.push_back(-1);
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(1);
        use_derivative.push_back(true);
        derivative_directions.push_back(0);
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(1);
        use_derivative.push_back(true);
        derivative_directions.push_back(0);
        
        std::vector<double> enstrophy_part_7 = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
            quantity_names,
            component_indices,
            use_derivative,
            derivative_directions,
            d_num_ghosts_derivative,
            data_context);
        
        quantity_names.clear();
        component_indices.clear();
        use_derivative.clear();
        derivative_directions.clear();
        
        quantity_names.push_back("DENSITY");
        component_indices.push_back(0);
        use_derivative.push_back(false);
        derivative_directions.push_back(-1);
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(0);
        use_derivative.push_back(true);
        derivative_directions.push_back(1);
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(0);
        use_derivative.push_back(true);
        derivative_directions.push_back(1);
        
        std::vector<double> enstrophy_part_8 = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
            quantity_names,
            component_indices,
            use_derivative,
            derivative_directions,
            d_num_ghosts_derivative,
            data_context);
        
        quantity_names.clear();
        component_indices.clear();
        use_derivative.clear();
        derivative_directions.clear();
        
        quantity_names.push_back("DENSITY");
        component_indices.push_back(0);
        use_derivative.push_back(false);
        derivative_directions.push_back(-1);
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(1);
        use_derivative.push_back(true);
        derivative_directions.push_back(0);
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(0);
        use_derivative.push_back(true);
        derivative_directions.push_back(1);
        
        std::vector<double> enstrophy_part_9 = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
            quantity_names,
            component_indices,
            use_derivative,
            derivative_directions,
            d_num_ghosts_derivative,
            data_context);
        
        quantity_names.clear();
        component_indices.clear();
        use_derivative.clear();
        derivative_directions.clear();
        
        std::vector<double> enstrophy_full(finest_level_dims[0], double(0));
        for (int i = 0; i < finest_level_dims[0]; i++)
        {
            enstrophy_full[i] = enstrophy_part_1[i] + enstrophy_part_2[i] - double(2)*enstrophy_part_3[i] +
                enstrophy_part_4[i] + enstrophy_part_5[i] - double(2)*enstrophy_part_6[i] +
                enstrophy_part_7[i] + enstrophy_part_8[i] - double(2)*enstrophy_part_9[i];
        }
        
        Omega_avg_realizations.push_back(enstrophy_full);
    }
    
    d_ensemble_statistics->Omega_avg_computed = true;
}


/*
 * Compute scalar dissipation rate with assumed homogeneity in y-direction (2D) or yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::computeAveragedScalarDissipationRateWithHomogeneityInYDirectionOrInYZPlane(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
    FlowModelMPIHelperAverage MPI_helper_average = FlowModelMPIHelperAverage(
        "MPI_helper_average",
        d_dim,
        d_grid_geometry,
        patch_hierarchy,
        flow_model_tmp,
        true);
    
    const hier::IntVector& finest_level_dims = MPI_helper_average.getFinestRefinedDomainNumberOfPoints();
    
    std::vector<std::vector<double> >& chi_avg_realizations = d_ensemble_statistics->chi_avg_realizations;
    
    if (d_dim == tbox::Dimension(1))
    {
        std::vector<std::string> quantity_names;
        std::vector<int> component_indices;
        std::vector<bool> use_derivative;
        std::vector<int> derivative_directions;
        
        quantity_names.push_back("MASS_DIFFUSIVITIES");
        component_indices.push_back(0);
        use_derivative.push_back(false);
        derivative_directions.push_back(-1);
        quantity_names.push_back("MASS_FRACTIONS");
        component_indices.push_back(0);
        use_derivative.push_back(true);
        derivative_directions.push_back(0);
        quantity_names.push_back("MASS_FRACTIONS");
        component_indices.push_back(0);
        use_derivative.push_back(true);
        derivative_directions.push_back(0);
        
        std::vector<double> scalar_dissipation_rate = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
            quantity_names,
            component_indices,
            use_derivative,
            derivative_directions,
            d_num_ghosts_derivative,
            data_context);
        
        quantity_names.clear();
        component_indices.clear();
        use_derivative.clear();
        derivative_directions.clear();
        
        chi_avg_realizations.push_back(scalar_dissipation_rate);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        std::vector<std::string> quantity_names;
        std::vector<int> component_indices;
        std::vector<bool> use_derivative;
        std::vector<int> derivative_directions;
        
        quantity_names.push_back("MASS_DIFFUSIVITIES");
        component_indices.push_back(0);
        use_derivative.push_back(false);
        derivative_directions.push_back(-1);
        quantity_names.push_back("MASS_FRACTIONS");
        component_indices.push_back(0);
        use_derivative.push_back(true);
        derivative_directions.push_back(0);
        quantity_names.push_back("MASS_FRACTIONS");
        component_indices.push_back(0);
        use_derivative.push_back(true);
        derivative_directions.push_back(0);
        
        std::vector<double> scalar_dissipation_rate_part_1 = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
            quantity_names,
            component_indices,
            use_derivative,
            derivative_directions,
            d_num_ghosts_derivative,
            data_context);
        
        quantity_names.clear();
        component_indices.clear();
        use_derivative.clear();
        derivative_directions.clear();
        
        quantity_names.push_back("MASS_DIFFUSIVITIES");
        component_indices.push_back(0);
        use_derivative.push_back(false);
        derivative_directions.push_back(-1);
        quantity_names.push_back("MASS_FRACTIONS");
        component_indices.push_back(0);
        use_derivative.push_back(true);
        derivative_directions.push_back(1);
        quantity_names.push_back("MASS_FRACTIONS");
        component_indices.push_back(0);
        use_derivative.push_back(true);
        derivative_directions.push_back(1);
        
        std::vector<double> scalar_dissipation_rate_part_2 = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
            quantity_names,
            component_indices,
            use_derivative,
            derivative_directions,
            d_num_ghosts_derivative,
            data_context);
        
        quantity_names.clear();
        component_indices.clear();
        use_derivative.clear();
        derivative_directions.clear();
        
        std::vector<double> scalar_dissipation_rate_full(finest_level_dims[0], double(0));
        for (int i = 0; i < finest_level_dims[0]; i++)
        {
            scalar_dissipation_rate_full[i] = scalar_dissipation_rate_part_1[i] + scalar_dissipation_rate_part_2[i];
        }
        
        chi_avg_realizations.push_back(scalar_dissipation_rate_full);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        std::vector<std::string> quantity_names;
        std::vector<int> component_indices;
        std::vector<bool> use_derivative;
        std::vector<int> derivative_directions;
        
        quantity_names.push_back("MASS_DIFFUSIVITIES");
        component_indices.push_back(0);
        use_derivative.push_back(false);
        derivative_directions.push_back(-1);
        quantity_names.push_back("MASS_FRACTIONS");
        component_indices.push_back(0);
        use_derivative.push_back(true);
        derivative_directions.push_back(0);
        quantity_names.push_back("MASS_FRACTIONS");
        component_indices.push_back(0);
        use_derivative.push_back(true);
        derivative_directions.push_back(0);
        
        std::vector<double> scalar_dissipation_rate_part_1 = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
            quantity_names,
            component_indices,
            use_derivative,
            derivative_directions,
            d_num_ghosts_derivative,
            data_context);
        
        quantity_names.clear();
        component_indices.clear();
        use_derivative.clear();
        derivative_directions.clear();
        
        quantity_names.push_back("MASS_DIFFUSIVITIES");
        component_indices.push_back(0);
        use_derivative.push_back(false);
        derivative_directions.push_back(-1);
        quantity_names.push_back("MASS_FRACTIONS");
        component_indices.push_back(0);
        use_derivative.push_back(true);
        derivative_directions.push_back(1);
        quantity_names.push_back("MASS_FRACTIONS");
        component_indices.push_back(0);
        use_derivative.push_back(true);
        derivative_directions.push_back(1);
        
        std::vector<double> scalar_dissipation_rate_part_2 = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
            quantity_names,
            component_indices,
            use_derivative,
            derivative_directions,
            d_num_ghosts_derivative,
            data_context);
        
        quantity_names.clear();
        component_indices.clear();
        use_derivative.clear();
        derivative_directions.clear();
        
        quantity_names.push_back("MASS_DIFFUSIVITIES");
        component_indices.push_back(0);
        use_derivative.push_back(false);
        derivative_directions.push_back(-1);
        quantity_names.push_back("MASS_FRACTIONS");
        component_indices.push_back(0);
        use_derivative.push_back(true);
        derivative_directions.push_back(2);
        quantity_names.push_back("MASS_FRACTIONS");
        component_indices.push_back(0);
        use_derivative.push_back(true);
        derivative_directions.push_back(2);
        
        std::vector<double> scalar_dissipation_rate_part_3 = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
            quantity_names,
            component_indices,
            use_derivative,
            derivative_directions,
            d_num_ghosts_derivative,
            data_context);
        
        quantity_names.clear();
        component_indices.clear();
        use_derivative.clear();
        derivative_directions.clear();
        
        std::vector<double> scalar_dissipation_rate_full(finest_level_dims[0], double(0));
        for (int i = 0; i < finest_level_dims[0]; i++)
        {
            scalar_dissipation_rate_full[i] =
                scalar_dissipation_rate_part_1[i] + scalar_dissipation_rate_part_2[i] + scalar_dissipation_rate_part_3[i];
        }
        
        chi_avg_realizations.push_back(scalar_dissipation_rate_full);
    }
    
    d_ensemble_statistics->chi_avg_computed = true;
}


/*
 * Compute averaged dynamic shear viscosity with assumed homogeneity in y-direction (2D) or yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::computeAveragedShearViscosityWithHomogeneityInYDirectionOrInYZPlane(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
    FlowModelMPIHelperAverage MPI_helper_average = FlowModelMPIHelperAverage(
        "MPI_helper_average",
        d_dim,
        d_grid_geometry,
        patch_hierarchy,
        flow_model_tmp,
        true);
    
    const hier::IntVector& finest_level_dims = MPI_helper_average.getFinestRefinedDomainNumberOfPoints();
    
    std::vector<std::vector<double> >& mu_avg_realizations = d_ensemble_statistics->mu_avg_realizations;
    
    std::vector<double> shear_viscosity = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        "SHEAR_VISCOSITY",
        0,
        data_context);
    
    mu_avg_realizations.push_back(shear_viscosity);
    
    d_ensemble_statistics->mu_avg_computed = true;
}


/*
 * Output spatial profile of ensemble averaged mass fraction with assumed homogeneity in y-direction (2D) or
 * yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputSpatialProfileEnsembleAveragedMassFractionWithHomogeneityInYDirectionOrInYZPlane(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    std::ofstream f_out;
    
    /*
     * Output the spatial profile (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        f_out.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        const std::vector<std::vector<double> >& Y_0_avg_realizations =
            d_ensemble_statistics->Y_0_avg_realizations;
        
        const int num_realizations = static_cast<int>(Y_0_avg_realizations.size());
        
        TBOX_ASSERT(d_ensemble_statistics->getNumberOfEnsembles() == num_realizations);
        TBOX_ASSERT(num_realizations > 0);
        
        const int num_cells = static_cast<int>(Y_0_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> Y_0_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                Y_0_avg_global[i] += weight*Y_0_avg_realizations[ri][i];
            }
        }
        
        f_out.write((char*)&output_time, sizeof(double));
        f_out.write((char*)&Y_0_avg_global[0], sizeof(double)*Y_0_avg_global.size());
        
        f_out.close();
    }
}


/*
 * Output spatial profile of ensemble variance of mass fraction with assumed homogeneity in y-direction (2D) or
 * yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputSpatialProfileEnsembleMassFractionVarianceWithHomogeneityInYDirectionOrInYZPlane(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    std::ofstream f_out;
    
    /*
     * Output the spatial profile (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        f_out.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        const std::vector<std::vector<double> >& Y_0_avg_realizations =
            d_ensemble_statistics->Y_0_avg_realizations;
        
        const std::vector<std::vector<double> >& Y_0_sq_avg_realizations = 
            d_ensemble_statistics->Y_0_sq_avg_realizations;
        
        const int num_realizations = static_cast<int>(Y_0_avg_realizations.size());
        
        TBOX_ASSERT(d_ensemble_statistics->getNumberOfEnsembles() == num_realizations);
        TBOX_ASSERT(num_realizations > 0);
        TBOX_ASSERT(num_realizations == static_cast<int>(Y_0_sq_avg_realizations.size()));
        
        const int num_cells = static_cast<int>(Y_0_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> Y_0_avg_global(num_cells, double(0));
        std::vector<double> Y_0_sq_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                Y_0_avg_global[i]    += weight*Y_0_avg_realizations[ri][i];
                Y_0_sq_avg_global[i] += weight*Y_0_sq_avg_realizations[ri][i];
            }
        }
        
        std::vector<double> Y_0_var_global(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            Y_0_var_global[i] = Y_0_sq_avg_global[i] - Y_0_avg_global[i]*Y_0_avg_global[i];
        }
        
        f_out.write((char*)&output_time, sizeof(double));
        f_out.write((char*)&Y_0_var_global[0], sizeof(double)*Y_0_var_global.size());
        
        f_out.close();
    }
}


/*
 * Output spatial profile of ensemble averaged density with assumed homogeneity in y-direction (2D) or
 * yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputSpatialProfileEnsembleAveragedDensityWithHomogeneityInYDirectionOrInYZPlane(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    std::ofstream f_out;
    
    /*
     * Output the spatial profile (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        f_out.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        const std::vector<std::vector<double> >& rho_avg_realizations =
            d_ensemble_statistics->rho_avg_realizations;
        
        const int num_realizations = static_cast<int>(rho_avg_realizations.size());
        
        TBOX_ASSERT(d_ensemble_statistics->getNumberOfEnsembles() == num_realizations);
        TBOX_ASSERT(num_realizations > 0);
        
        const int num_cells = static_cast<int>(rho_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> rho_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                rho_avg_global[i] += weight*rho_avg_realizations[ri][i];
            }
        }
        
        f_out.write((char*)&output_time, sizeof(double));
        f_out.write((char*)&rho_avg_global[0], sizeof(double)*rho_avg_global.size());
        
        f_out.close();
    }
}


/*
 * Output spatial profile of ensemble variance of density with assumed homogeneity in y-direction (2D) or
 * yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputSpatialProfileEnsembleDensityVarianceWithHomogeneityInYDirectionOrInYZPlane(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    std::ofstream f_out;
    
    /*
     * Output the spatial profile (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        f_out.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        const std::vector<std::vector<double> >& rho_avg_realizations =
            d_ensemble_statistics->rho_avg_realizations;
        
        const std::vector<std::vector<double> >& rho_sq_avg_realizations = 
            d_ensemble_statistics->rho_sq_avg_realizations;
        
        const int num_realizations = static_cast<int>(rho_avg_realizations.size());
        
        TBOX_ASSERT(d_ensemble_statistics->getNumberOfEnsembles() == num_realizations);
        TBOX_ASSERT(num_realizations > 0);
        TBOX_ASSERT(num_realizations == static_cast<int>(rho_sq_avg_realizations.size()));
        
        const int num_cells = static_cast<int>(rho_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> rho_avg_global(num_cells, double(0));
        std::vector<double> rho_sq_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                rho_avg_global[i]    += weight*rho_avg_realizations[ri][i];
                rho_sq_avg_global[i] += weight*rho_sq_avg_realizations[ri][i];
            }
        }
        
        std::vector<double> rho_var_global(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            rho_var_global[i] = rho_sq_avg_global[i] - rho_avg_global[i]*rho_avg_global[i];
        }
        
        f_out.write((char*)&output_time, sizeof(double));
        f_out.write((char*)&rho_var_global[0], sizeof(double)*rho_var_global.size());
        
        f_out.close();
    }
}


/*
 * Output spatial profile of ensemble averaged specific volume with assumed homogeneity in y-direction (2D) or
 * yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputSpatialProfileEnsembleAveragedSpecificVolumeWithHomogeneityInYDirectionOrInYZPlane(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    std::ofstream f_out;
    
    /*
     * Output the spatial profile (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        f_out.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        const std::vector<std::vector<double> >& rho_inv_avg_realizations =
            d_ensemble_statistics->rho_inv_avg_realizations;
        
        const int num_realizations = static_cast<int>(rho_inv_avg_realizations.size());
        
        TBOX_ASSERT(d_ensemble_statistics->getNumberOfEnsembles() == num_realizations);
        TBOX_ASSERT(num_realizations > 0);
        
        const int num_cells = static_cast<int>(rho_inv_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> rho_inv_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                rho_inv_avg_global[i] += weight*rho_inv_avg_realizations[ri][i];
            }
        }
        
        f_out.write((char*)&output_time, sizeof(double));
        f_out.write((char*)&rho_inv_avg_global[0], sizeof(double)*rho_inv_avg_global.size());
        
        f_out.close();
    }
}


/*
 * Output spatial profile of ensemble averaged density-specific-volume covariance with assumed homogeneity in
 * y-direction (2D) or yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputSpatialProfileEnsembleAveragedDensitySpecificVolumeCovarianceWithHomogeneityInYDirectionOrInYZPlane(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    std::ofstream f_out;
    
    /*
     * Output the spatial profile (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        f_out.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        const std::vector<std::vector<double> >& rho_avg_realizations =
            d_ensemble_statistics->rho_avg_realizations;
        
        const std::vector<std::vector<double> >& rho_inv_avg_realizations =
            d_ensemble_statistics->rho_inv_avg_realizations;
        
        const int num_realizations = static_cast<int>(rho_avg_realizations.size());
        
        TBOX_ASSERT(d_ensemble_statistics->getNumberOfEnsembles() == num_realizations);
        TBOX_ASSERT(num_realizations > 0);
        TBOX_ASSERT(num_realizations == static_cast<int>(rho_inv_avg_realizations.size()));
        
        const int num_cells = static_cast<int>(rho_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> rho_avg_global(num_cells, double(0));
        std::vector<double> rho_inv_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                rho_avg_global[i]     += weight*rho_avg_realizations[ri][i];
                rho_inv_avg_global[i] += weight*rho_inv_avg_realizations[ri][i];
            }
        }
        
        std::vector<double> rho_p_rho_inv_p(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            rho_p_rho_inv_p[i] = double(1) - rho_avg_global[i]*rho_inv_avg_global[i];
        }
        
        f_out.write((char*)&output_time, sizeof(double));
        f_out.write((char*)&rho_p_rho_inv_p[0], sizeof(double)*rho_p_rho_inv_p.size());
        
        f_out.close();
    }
}


/*
 * Output spatial profile of ensemble averaged velocity x-component with assumed homogeneity in y-direction (2D)
 * or yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputSpatialProfileEnsembleAveragedVelocityXWithHomogeneityInYDirectionOrInYZPlane(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    std::ofstream f_out;
    
    /*
     * Output the spatial profile (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        f_out.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        const std::vector<std::vector<double> >& u_avg_realizations =
            d_ensemble_statistics->u_avg_realizations;
        
        const int num_realizations = static_cast<int>(u_avg_realizations.size());
        
        TBOX_ASSERT(d_ensemble_statistics->getNumberOfEnsembles() == num_realizations);
        TBOX_ASSERT(num_realizations > 0);
        
        const int num_cells = static_cast<int>(u_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> u_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                u_avg_global[i] += weight*u_avg_realizations[ri][i];
            }
        }
        
        f_out.write((char*)&output_time, sizeof(double));
        f_out.write((char*)&u_avg_global[0], sizeof(double)*u_avg_global.size());
        
        f_out.close();
    }
}


/*
 * Output spatial profile of ensemble averaged velocity y-component with assumed homogeneity in y-direction (2D)
 * or yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputSpatialProfileEnsembleAveragedVelocityYWithHomogeneityInYDirectionOrInYZPlane(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    std::ofstream f_out;
    
    /*
     * Output the spatial profile (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        f_out.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        const std::vector<std::vector<double> >& v_avg_realizations =
            d_ensemble_statistics->v_avg_realizations;
        
        const int num_realizations = static_cast<int>(v_avg_realizations.size());
        
        TBOX_ASSERT(d_ensemble_statistics->getNumberOfEnsembles() == num_realizations);
        TBOX_ASSERT(num_realizations > 0);
        
        const int num_cells = static_cast<int>(v_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> v_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                v_avg_global[i] += weight*v_avg_realizations[ri][i];
            }
        }
        
        f_out.write((char*)&output_time, sizeof(double));
        f_out.write((char*)&v_avg_global[0], sizeof(double)*v_avg_global.size());
        
        f_out.close();
    }
}


/*
 * Output spatial profile of ensemble averaged velocity z-component with assumed homogeneity in yz-plane (3D)
 * to a file.
 */
void
RTIRMIStatisticsUtilities::outputSpatialProfileEnsembleAveragedVelocityZWithHomogeneityInYZPlane(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    std::ofstream f_out;
    
    /*
     * Output the spatial profile (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        f_out.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        const std::vector<std::vector<double> >& w_avg_realizations =
            d_ensemble_statistics->w_avg_realizations;
        
        const int num_realizations = static_cast<int>(w_avg_realizations.size());
        
        TBOX_ASSERT(d_ensemble_statistics->getNumberOfEnsembles() == num_realizations);
        TBOX_ASSERT(num_realizations > 0);
        
        const int num_cells = static_cast<int>(w_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> w_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                w_avg_global[i] += weight*w_avg_realizations[ri][i];
            }
        }
        
        f_out.write((char*)&output_time, sizeof(double));
        f_out.write((char*)&w_avg_global[0], sizeof(double)*w_avg_global.size());
        
        f_out.close();
    }
}


/*
 * Output spatial profile of ensemble averaged momentum x-component with assumed homogeneity in y-direction (2D)
 * or yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputSpatialProfileEnsembleAveragedMomentumXWithHomogeneityInYDirectionOrInYZPlane(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    std::ofstream f_out;
    
    /*
     * Output the spatial profile (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        f_out.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        const std::vector<std::vector<double> >& rho_u_avg_realizations =
            d_ensemble_statistics->rho_u_avg_realizations;
        
        const int num_realizations = static_cast<int>(rho_u_avg_realizations.size());
        
        TBOX_ASSERT(d_ensemble_statistics->getNumberOfEnsembles() == num_realizations);
        TBOX_ASSERT(num_realizations > 0);
        
        const int num_cells = static_cast<int>(rho_u_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> rho_u_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                rho_u_avg_global[i] += weight*rho_u_avg_realizations[ri][i];
            }
        }
        
        f_out.write((char*)&output_time, sizeof(double));
        f_out.write((char*)&rho_u_avg_global[0], sizeof(double)*rho_u_avg_global.size());
        
        f_out.close();
    }
}


/*
 * Output spatial profile of ensemble averaged momentum y-component with assumed homogeneity in y-direction (2D)
 * or yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputSpatialProfileEnsembleAveragedMomentumYWithHomogeneityInYDirectionOrInYZPlane(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    std::ofstream f_out;
    
    /*
     * Output the spatial profile (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        f_out.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        const std::vector<std::vector<double> >& rho_v_avg_realizations =
            d_ensemble_statistics->rho_v_avg_realizations;
        
        const int num_realizations = static_cast<int>(rho_v_avg_realizations.size());
        
        TBOX_ASSERT(d_ensemble_statistics->getNumberOfEnsembles() == num_realizations);
        TBOX_ASSERT(num_realizations > 0);
        
        const int num_cells = static_cast<int>(rho_v_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> rho_v_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                rho_v_avg_global[i] += weight*rho_v_avg_realizations[ri][i];
            }
        }
        
        f_out.write((char*)&output_time, sizeof(double));
        f_out.write((char*)&rho_v_avg_global[0], sizeof(double)*rho_v_avg_global.size());
        
        f_out.close();
    }
}


/*
 * Output spatial profile of ensemble averaged momentum z-component with assumed homogeneity in yz-plane (3D)
 * to a file.
 */
void
RTIRMIStatisticsUtilities::outputSpatialProfileEnsembleAveragedMomentumZWithHomogeneitynYZPlane(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    std::ofstream f_out;
    
    /*
     * Output the spatial profile (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        f_out.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        const std::vector<std::vector<double> >& rho_w_avg_realizations =
            d_ensemble_statistics->rho_w_avg_realizations;
        
        const int num_realizations = static_cast<int>(rho_w_avg_realizations.size());
        
        TBOX_ASSERT(d_ensemble_statistics->getNumberOfEnsembles() == num_realizations);
        TBOX_ASSERT(num_realizations > 0);
        
        const int num_cells = static_cast<int>(rho_w_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> rho_w_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                rho_w_avg_global[i] += weight*rho_w_avg_realizations[ri][i];
            }
        }
        
        f_out.write((char*)&output_time, sizeof(double));
        f_out.write((char*)&rho_w_avg_global[0], sizeof(double)*rho_w_avg_global.size());
        
        f_out.close();
    }
}


/*
 * Output spatial profile of ensemble turbulent mass flux in x-direction with assumed homogeneity in
 * y-direction (2D) or yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputSpatialProfileEnsembleTurbulentMassFluxXWithHomogeneityInYDirectionOrInYZPlane(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    std::ofstream f_out;
    
    /*
     * Output the spatial profile (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        f_out.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        const std::vector<std::vector<double> >& rho_avg_realizations =
            d_ensemble_statistics->rho_avg_realizations;
        
        const std::vector<std::vector<double> >& u_avg_realizations =
            d_ensemble_statistics->u_avg_realizations;
        
        const std::vector<std::vector<double> >& rho_u_avg_realizations =
            d_ensemble_statistics->rho_u_avg_realizations;
        
        const int num_realizations = static_cast<int>(rho_avg_realizations.size());
        
        TBOX_ASSERT(d_ensemble_statistics->getNumberOfEnsembles() == num_realizations);
        TBOX_ASSERT(num_realizations > 0);
        TBOX_ASSERT(num_realizations == static_cast<int>(u_avg_realizations.size()));
        TBOX_ASSERT(num_realizations == static_cast<int>(rho_u_avg_realizations.size()));
        
        const int num_cells = static_cast<int>(rho_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> rho_avg_global(num_cells, double(0));
        std::vector<double> u_avg_global(num_cells, double(0));
        std::vector<double> rho_u_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                rho_avg_global[i]   += weight*rho_avg_realizations[ri][i];
                u_avg_global[i]     += weight*u_avg_realizations[ri][i];
                rho_u_avg_global[i] += weight*rho_u_avg_realizations[ri][i];
            }
        }
        
        std::vector<double> rho_p_u_p(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            rho_p_u_p[i] = rho_u_avg_global[i] - rho_avg_global[i]*u_avg_global[i];
        }
        
        f_out.write((char*)&output_time, sizeof(double));
        f_out.write((char*)&rho_p_u_p[0], sizeof(double)*rho_p_u_p.size());
        
        f_out.close();
    }
}


/*
 * Output spatial profile of ensemble Reynolds normal stress x-component with assumed homogeneity in
 * y-direction (2D) or yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputSpatialProfileEnsembleReynoldsNormalStressXWithHomogeneityInYDirectionOrInYZPlane(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    std::ofstream f_out;
    
    /*
     * Output the spatial profile (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        f_out.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        const std::vector<std::vector<double> >& rho_avg_realizations =
            d_ensemble_statistics->rho_avg_realizations;
        
        const std::vector<std::vector<double> >& rho_u_avg_realizations =
            d_ensemble_statistics->rho_u_avg_realizations;
        
        const std::vector<std::vector<double> >& rho_u_u_avg_realizations =
            d_ensemble_statistics->rho_u_u_avg_realizations;
        
        const int num_realizations = static_cast<int>(rho_avg_realizations.size());
        
        TBOX_ASSERT(d_ensemble_statistics->getNumberOfEnsembles() == num_realizations);
        TBOX_ASSERT(num_realizations > 0);
        TBOX_ASSERT(num_realizations == static_cast<int>(rho_avg_realizations.size()));
        TBOX_ASSERT(num_realizations == static_cast<int>(rho_u_avg_realizations.size()));
        TBOX_ASSERT(num_realizations == static_cast<int>(rho_u_u_avg_realizations.size()));
        
        const int num_cells = static_cast<int>(rho_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> rho_avg_global(num_cells, double(0));
        std::vector<double> rho_u_avg_global(num_cells, double(0));
        std::vector<double> rho_u_u_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                rho_avg_global[i]     += weight*rho_avg_realizations[ri][i];
                rho_u_avg_global[i]   += weight*rho_u_avg_realizations[ri][i];
                rho_u_u_avg_global[i] += weight*rho_u_u_avg_realizations[ri][i];
            }
        }
        
        std::vector<double> R_11(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            const double u_tilde = rho_u_avg_global[i]/rho_avg_global[i];
            const double rho_u_pp_u_pp = rho_u_u_avg_global[i] - rho_u_avg_global[i]*u_tilde;
            R_11[i] = rho_u_pp_u_pp/rho_avg_global[i];
        }
        
        f_out.write((char*)&output_time, sizeof(double));
        f_out.write((char*)&R_11[0], sizeof(double)*R_11.size());
        
        f_out.close();
    }
}


/*
 * Output spatial profile of ensemble Reynolds normal stress y-component with assumed homogeneity in
 * y-direction (2D) or yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputSpatialProfileEnsembleReynoldsNormalStressYWithHomogeneityInYDirectionOrInYZPlane(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    std::ofstream f_out;
    
    /*
     * Output the spatial profile (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        f_out.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        const std::vector<std::vector<double> >& rho_avg_realizations =
            d_ensemble_statistics->rho_avg_realizations;
        
        const std::vector<std::vector<double> >& rho_v_avg_realizations =
            d_ensemble_statistics->rho_v_avg_realizations;
        
        const std::vector<std::vector<double> >& rho_v_v_avg_realizations =
            d_ensemble_statistics->rho_v_v_avg_realizations;
        
        const int num_realizations = static_cast<int>(rho_avg_realizations.size());
        
        TBOX_ASSERT(d_ensemble_statistics->getNumberOfEnsembles() == num_realizations);
        TBOX_ASSERT(num_realizations > 0);
        TBOX_ASSERT(num_realizations == static_cast<int>(rho_avg_realizations.size()));
        TBOX_ASSERT(num_realizations == static_cast<int>(rho_v_avg_realizations.size()));
        TBOX_ASSERT(num_realizations == static_cast<int>(rho_v_v_avg_realizations.size()));
        
        const int num_cells = static_cast<int>(rho_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> rho_avg_global(num_cells, double(0));
        std::vector<double> rho_v_avg_global(num_cells, double(0));
        std::vector<double> rho_v_v_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                rho_avg_global[i]     += weight*rho_avg_realizations[ri][i];
                rho_v_avg_global[i]   += weight*rho_v_avg_realizations[ri][i];
                rho_v_v_avg_global[i] += weight*rho_v_v_avg_realizations[ri][i];
            }
        }
        
        std::vector<double> R_22(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            const double v_tilde = rho_v_avg_global[i]/rho_avg_global[i];
            const double rho_v_pp_v_pp = rho_v_v_avg_global[i] - rho_v_avg_global[i]*v_tilde;
            R_22[i] = rho_v_pp_v_pp/rho_avg_global[i];
        }
        
        f_out.write((char*)&output_time, sizeof(double));
        f_out.write((char*)&R_22[0], sizeof(double)*R_22.size());
        
        f_out.close();
    }
}


/*
 * Output spatial profile of ensemble Reynolds normal stress z-component with assumed homogeneity in
 * yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputSpatialProfileEnsembleReynoldsNormalStressZWithHomogeneityInYZPlane(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    std::ofstream f_out;
    
    /*
     * Output the spatial profile (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        f_out.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        const std::vector<std::vector<double> >& rho_avg_realizations =
            d_ensemble_statistics->rho_avg_realizations;
        
        const std::vector<std::vector<double> >& rho_w_avg_realizations =
            d_ensemble_statistics->rho_w_avg_realizations;
        
        const std::vector<std::vector<double> >& rho_w_w_avg_realizations =
            d_ensemble_statistics->rho_w_w_avg_realizations;
        
        const int num_realizations = static_cast<int>(rho_avg_realizations.size());
        
        TBOX_ASSERT(d_ensemble_statistics->getNumberOfEnsembles() == num_realizations);
        TBOX_ASSERT(num_realizations > 0);
        TBOX_ASSERT(num_realizations == static_cast<int>(rho_avg_realizations.size()));
        TBOX_ASSERT(num_realizations == static_cast<int>(rho_w_avg_realizations.size()));
        TBOX_ASSERT(num_realizations == static_cast<int>(rho_w_w_avg_realizations.size()));
        
        const int num_cells = static_cast<int>(rho_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> rho_avg_global(num_cells, double(0));
        std::vector<double> rho_w_avg_global(num_cells, double(0));
        std::vector<double> rho_w_w_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                rho_avg_global[i]     += weight*rho_avg_realizations[ri][i];
                rho_w_avg_global[i]   += weight*rho_w_avg_realizations[ri][i];
                rho_w_w_avg_global[i] += weight*rho_w_w_avg_realizations[ri][i];
            }
        }
        
        std::vector<double> R_33(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            const double w_tilde = rho_w_avg_global[i]/rho_avg_global[i];
            const double rho_w_pp_w_pp = rho_w_w_avg_global[i] - rho_w_avg_global[i]*w_tilde;
            R_33[i] = rho_w_pp_w_pp/rho_avg_global[i];
        }
        
        f_out.write((char*)&output_time, sizeof(double));
        f_out.write((char*)&R_33[0], sizeof(double)*R_33.size());
        
        f_out.close();
    }
}


/*
 * Output spatial profile of ensemble averaged enstrophy with assumed homogeneity in y-direction (2D) or
 * yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputSpatialProfileEnsembleAveragedEnstrophyWithHomogeneityInYDirectionOrInYZPlane(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    std::ofstream f_out;
    
    /*
     * Output the spatial profile (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        f_out.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        const std::vector<std::vector<double> >& Omega_avg_realizations =
            d_ensemble_statistics->Omega_avg_realizations;
        
        const int num_realizations = static_cast<int>(Omega_avg_realizations.size());
        
        TBOX_ASSERT(d_ensemble_statistics->getNumberOfEnsembles() == num_realizations);
        TBOX_ASSERT(num_realizations > 0);
        
        const int num_cells = static_cast<int>(Omega_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> Omega_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                Omega_avg_global[i] += weight*Omega_avg_realizations[ri][i];
            }
        }
        
        f_out.write((char*)&output_time, sizeof(double));
        f_out.write((char*)&Omega_avg_global[0], sizeof(double)*Omega_avg_global.size());
        
        f_out.close();
    }
}


/*
 * Output spatial profile of ensemble averaged scalar dissipation rate with assumed homogeneity in
 * y-direction (2D) or yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputSpatialProfileEnsembleAveragedScalarDissipationRateWithHomogeneityInYDirectionOrInYZPlane(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    std::ofstream f_out;
    
    /*
     * Output the spatial profile (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        f_out.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        const std::vector<std::vector<double> >& chi_avg_realizations =
            d_ensemble_statistics->chi_avg_realizations;
        
        const int num_realizations = static_cast<int>(chi_avg_realizations.size());
        
        TBOX_ASSERT(d_ensemble_statistics->getNumberOfEnsembles() == num_realizations);
        TBOX_ASSERT(num_realizations > 0);
        
        const int num_cells = static_cast<int>(chi_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> chi_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                chi_avg_global[i] += weight*chi_avg_realizations[ri][i];
            }
        }
        
        f_out.write((char*)&output_time, sizeof(double));
        f_out.write((char*)&chi_avg_global[0], sizeof(double)*chi_avg_global.size());
        
        f_out.close();
    }
}


/*
 * Output ensemble mixing width in x-direction to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleMixingWidthInXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'MIXING_WIDTH_X' can be computed with two species only."
            << std::endl);
    }
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    MPIHelper MPI_helper = MPIHelper(
        "MPI_helper",
        d_dim,
        d_grid_geometry,
        patch_hierarchy);
    
    const std::vector<double>& dx_finest = MPI_helper.getFinestRefinedDomainGridSpacing();
    
    /*
     * Compute and output the quantity (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_out;
        
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        const std::vector<std::vector<double> >& Y_0_avg_realizations =
            d_ensemble_statistics->Y_0_avg_realizations;
        
        const int num_realizations = static_cast<int>(Y_0_avg_realizations.size());
        
        TBOX_ASSERT(d_ensemble_statistics->getNumberOfEnsembles() == num_realizations);
        TBOX_ASSERT(num_realizations > 0);
        
        const int num_cells = static_cast<int>(Y_0_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> Y_0_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                Y_0_avg_global[i] += weight*Y_0_avg_realizations[ri][i];
            }
        }
        
        double W = double(0);
        
        for (int i = 0; i < num_cells; i++)
        {
            W += Y_0_avg_global[i]*(double(1) - Y_0_avg_global[i]);
        }
        
        W = double(4)*W*dx_finest[0];
        
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << W;
        
        f_out.close();
    }
}


/*
 * Output ensemble minimum interface location in x-direction to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleInterfaceMinInXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'INTERFACE_MIN_X' can be computed with two species only."
            << std::endl);
    }
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    MPIHelper MPI_helper = MPIHelper(
        "MPI_helper",
        d_dim,
        d_grid_geometry,
        patch_hierarchy);
    
    const std::vector<double>& dx_finest = MPI_helper.getFinestRefinedDomainGridSpacing();
    
    /*
     * Compute and output the quantity (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_out;
        
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        const std::vector<std::vector<double> >& Y_0_avg_realizations =
            d_ensemble_statistics->Y_0_avg_realizations;
        
        const int num_realizations = static_cast<int>(Y_0_avg_realizations.size());
        
        TBOX_ASSERT(d_ensemble_statistics->getNumberOfEnsembles() == num_realizations);
        TBOX_ASSERT(num_realizations > 0);
        
        const int num_cells = static_cast<int>(Y_0_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> Y_0_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                Y_0_avg_global[i] += weight*Y_0_avg_realizations[ri][i];
            }
        }
        
        const double* x_lo = d_grid_geometry->getXLower();
        const double* x_hi = d_grid_geometry->getXUpper();
        double interface_min = x_hi[0];
        
        for (int i = num_cells - 1; i >= 0;  i--)
        {
            if (Y_0_avg_global[i] > 0.01 && Y_0_avg_global[i] < 0.99)
            {
                const double x_loc = x_lo[0] + 0.5*dx_finest[0] + i*dx_finest[0];
                if (x_loc < interface_min)
                {
                   interface_min = x_loc;
                }
            }
        }
        
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << interface_min;
        
        f_out.close();
    }
}


/*
 * Output ensemble maximum interface location in x-direction to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleInterfaceMaxInXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'INTERFACE_MIN_X' can be computed with two species only."
            << std::endl);
    }
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    MPIHelper MPI_helper = MPIHelper(
        "MPI_helper",
        d_dim,
        d_grid_geometry,
        patch_hierarchy);
    
    const std::vector<double>& dx_finest = MPI_helper.getFinestRefinedDomainGridSpacing();
    
    /*
     * Compute and output the quantity (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_out;
        
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        const std::vector<std::vector<double> >& Y_0_avg_realizations =
            d_ensemble_statistics->Y_0_avg_realizations;
        
        const int num_realizations = static_cast<int>(Y_0_avg_realizations.size());
        
        TBOX_ASSERT(d_ensemble_statistics->getNumberOfEnsembles() == num_realizations);
        TBOX_ASSERT(num_realizations > 0);
        
        const int num_cells = static_cast<int>(Y_0_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> Y_0_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                Y_0_avg_global[i] += weight*Y_0_avg_realizations[ri][i];
            }
        }
        
        const double* x_lo = d_grid_geometry->getXLower();
        // const double* x_hi = d_grid_geometry->getXUpper();
        double interface_max = x_lo[0];
        
        for (int i = 0; i < num_cells; i++)
        {
            if (Y_0_avg_global[i] > 0.01 && Y_0_avg_global[i] < 0.99)
            {
                const double x_loc = x_lo[0] + 0.5*dx_finest[0] + i*dx_finest[0];
                if (x_loc > interface_max)
                {
                   interface_max = x_loc;
                }
            }
        }
        
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << interface_max;
        
        f_out.close();
    }
}


/*
 * Output ensemble mixedness in x-direction to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleMixednessInXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'MIXEDNESS_X' can be computed with two species only."
            << std::endl);
    }
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Compute and output the quantity (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_out;
        
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        const std::vector<std::vector<double> >& Y_0_avg_realizations =
            d_ensemble_statistics->Y_0_avg_realizations;
        
        const std::vector<std::vector<double> >& Y_0_Y_1_avg_realizations =
            d_ensemble_statistics->Y_0_Y_1_avg_realizations;
        
        const int num_realizations = static_cast<int>(Y_0_avg_realizations.size());
        
        TBOX_ASSERT(d_ensemble_statistics->getNumberOfEnsembles() == num_realizations);
        TBOX_ASSERT(num_realizations > 0);
        TBOX_ASSERT(num_realizations == static_cast<int>(Y_0_Y_1_avg_realizations.size()));
        
        const int num_cells = static_cast<int>(Y_0_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> Y_0_avg_global(num_cells, double(0));
        std::vector<double> Y_0_Y_1_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                Y_0_avg_global[i]     += weight*Y_0_avg_realizations[ri][i];
                Y_0_Y_1_avg_global[i] += weight*Y_0_Y_1_avg_realizations[ri][i];
            }
        }
        
        double num = double(0);
        double den = double(0);
        
        for (int i = 0; i < num_cells; i++)
        {
            num += Y_0_Y_1_avg_global[i];
            den += Y_0_avg_global[i]*(double(1) - Y_0_avg_global[i]);
        }
        
        const double Theta = num/den;
        
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << Theta;
        
        f_out.close();
    }
}


/*
 * Output ensemble mean velocity associated with turbulent mass flux component in x-direction with
 * assumed homogeneity in y-direction (2D) or yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleTurbulentMassFluxVelocityXWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Compute and output the quantity (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_out;
        
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        const std::vector<std::vector<double> >& rho_avg_realizations =
            d_ensemble_statistics->rho_avg_realizations;
        
        const std::vector<std::vector<double> >& u_avg_realizations =
            d_ensemble_statistics->u_avg_realizations;
        
        const std::vector<std::vector<double> >& rho_u_avg_realizations =
            d_ensemble_statistics->rho_u_avg_realizations;
        
        const std::vector<std::vector<double> >& Y_0_avg_realizations =
            d_ensemble_statistics->Y_0_avg_realizations;
        
        const int num_realizations = static_cast<int>(rho_avg_realizations.size());
        
        TBOX_ASSERT(d_ensemble_statistics->getNumberOfEnsembles() == num_realizations);
        TBOX_ASSERT(num_realizations > 0);
        TBOX_ASSERT(num_realizations == static_cast<int>(u_avg_realizations.size()));
        TBOX_ASSERT(num_realizations == static_cast<int>(rho_u_avg_realizations.size()));
        TBOX_ASSERT(num_realizations == static_cast<int>(Y_0_avg_realizations.size()));
        
        const int num_cells = static_cast<int>(rho_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> rho_avg_global(num_cells, double(0));
        std::vector<double> u_avg_global(num_cells, double(0));
        std::vector<double> rho_u_avg_global(num_cells, double(0));
        std::vector<double> Y_0_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                rho_avg_global[i]   += weight*rho_avg_realizations[ri][i];
                u_avg_global[i]     += weight*u_avg_realizations[ri][i];
                rho_u_avg_global[i] += weight*rho_u_avg_realizations[ri][i];
                Y_0_avg_global[i]   += weight*Y_0_avg_realizations[ri][i];
            }
        }
        
        std::vector<double> rho_p_u_p(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            rho_p_u_p[i] = rho_u_avg_global[i] - rho_avg_global[i]*u_avg_global[i];
        }
        
        double a_sum = double(0);
        int count = 0;
        
        for (int i = 0; i < num_cells; i++)
        {
            const double mixing_metric = double(4)*Y_0_avg_global[i]*(double(1) - Y_0_avg_global[i]);
            if (mixing_metric > double(9)/double(10))
            {
                a_sum += rho_p_u_p[i]/rho_avg_global[i];
                count++;
            }
        }
        
        const double a_mean = a_sum/count;
        
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << a_mean;
        
        f_out.close();
    }
}


/*
 * Output ensemble mean density-specific-volume covariance with assumed homogeneity in y-direction (2D)
 * or yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleDensitySpecificVolumeCovarianceWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Compute and output the quantity (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_out;
        
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        const std::vector<std::vector<double> >& rho_avg_realizations =
            d_ensemble_statistics->rho_avg_realizations;
        
        const std::vector<std::vector<double> >& rho_inv_avg_realizations =
            d_ensemble_statistics->rho_inv_avg_realizations;
        
        const std::vector<std::vector<double> >& Y_0_avg_realizations =
            d_ensemble_statistics->Y_0_avg_realizations;
        
        const int num_realizations = static_cast<int>(rho_avg_realizations.size());
        
        TBOX_ASSERT(d_ensemble_statistics->getNumberOfEnsembles() == num_realizations);
        TBOX_ASSERT(num_realizations > 0);
        TBOX_ASSERT(num_realizations == static_cast<int>(rho_inv_avg_realizations.size()));
        TBOX_ASSERT(num_realizations == static_cast<int>(Y_0_avg_realizations.size()));
        
        const int num_cells = static_cast<int>(rho_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> rho_avg_global(num_cells, double(0));
        std::vector<double> rho_inv_avg_global(num_cells, double(0));
        std::vector<double> Y_0_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                rho_avg_global[i]     += weight*rho_avg_realizations[ri][i];
                rho_inv_avg_global[i] += weight*rho_inv_avg_realizations[ri][i];
                Y_0_avg_global[i]     += weight*Y_0_avg_realizations[ri][i];
            }
        }
        
        std::vector<double> rho_p_rho_inv_p(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            rho_p_rho_inv_p[i] = double(1) - rho_avg_global[i]*rho_inv_avg_global[i];
        }
        
        double b_sum = double(0);
        int count = 0;
        
        for (int i = 0; i < num_cells; i++)
        {
            const double mixing_metric = double(4)*Y_0_avg_global[i]*(double(1) - Y_0_avg_global[i]);
            if (mixing_metric > double(9)/double(10))
            {
                b_sum -= rho_p_rho_inv_p[i];
                count++;
            }
        }
        
        const double b_mean = b_sum/count;
        
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << b_mean;
        
        f_out.close();
    }
}


/*
 * Output ensemble TKE integrated with assumed homogeneity in y-direction (2D) or yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleTKEIntegrateWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    MPIHelper MPI_helper = MPIHelper(
        "MPI_helper",
        d_dim,
        d_grid_geometry,
        patch_hierarchy);
    
    const std::vector<double>& dx_finest = MPI_helper.getFinestRefinedDomainGridSpacing();
    
    /*
     * Compute and output the quantity (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_out;
        
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        double TKE_integrated_global = double(0);
        const double half = double(1)/double(2);
        
        if (d_dim == tbox::Dimension(1))
        {
            const std::vector<std::vector<double> >& rho_avg_realizations =
                d_ensemble_statistics->rho_avg_realizations;
            
            const std::vector<std::vector<double> >& rho_u_avg_realizations =
                d_ensemble_statistics->rho_u_avg_realizations;
            
            const std::vector<std::vector<double> >& rho_u_u_avg_realizations =
                d_ensemble_statistics->rho_u_u_avg_realizations;
            
            const int num_realizations = static_cast<int>(rho_avg_realizations.size());
            
            TBOX_ASSERT(d_ensemble_statistics->getNumberOfEnsembles() == num_realizations);
            TBOX_ASSERT(num_realizations > 0);
            TBOX_ASSERT(num_realizations == static_cast<int>(rho_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(rho_u_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(rho_u_u_avg_realizations.size()));
            
            const int num_cells = static_cast<int>(rho_avg_realizations[0].size());
            const double weight = double(1)/double(num_realizations);
            
            std::vector<double> rho_avg_global(num_cells, double(0));
            std::vector<double> rho_u_avg_global(num_cells, double(0));
            std::vector<double> rho_u_u_avg_global(num_cells, double(0));
            
            for (int ri = 0; ri < num_realizations; ri++)
            {
                for (int i = 0; i < num_cells; i++)
                {
                    rho_avg_global[i]     += weight*rho_avg_realizations[ri][i];
                    rho_u_avg_global[i]   += weight*rho_u_avg_realizations[ri][i];
                    rho_u_u_avg_global[i] += weight*rho_u_u_avg_realizations[ri][i];
                }
            }
            
            for (int i = 0; i < num_cells; i++)
            {
                const double u_tilde = rho_u_avg_global[i]/rho_avg_global[i];
                const double rho_u_pp_u_pp = rho_u_u_avg_global[i] - rho_u_avg_global[i]*u_tilde;
                
                TKE_integrated_global += (half*rho_u_pp_u_pp*dx_finest[0]);
            }
        }
        else if (d_dim == tbox::Dimension(2))
        {
            const std::vector<std::vector<double> >& rho_avg_realizations =
                d_ensemble_statistics->rho_avg_realizations;
            
            const std::vector<std::vector<double> >& rho_u_avg_realizations =
                d_ensemble_statistics->rho_u_avg_realizations;
            
            const std::vector<std::vector<double> >& rho_v_avg_realizations =
                d_ensemble_statistics->rho_u_avg_realizations;
            
            const std::vector<std::vector<double> >& rho_u_u_avg_realizations =
                d_ensemble_statistics->rho_u_u_avg_realizations;
            
            const std::vector<std::vector<double> >& rho_v_v_avg_realizations =
                d_ensemble_statistics->rho_u_u_avg_realizations;
            
            const int num_realizations = static_cast<int>(rho_avg_realizations.size());
            
            TBOX_ASSERT(d_ensemble_statistics->getNumberOfEnsembles() == num_realizations);
            TBOX_ASSERT(num_realizations > 0);
            TBOX_ASSERT(num_realizations == static_cast<int>(rho_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(rho_u_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(rho_v_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(rho_u_u_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(rho_v_v_avg_realizations.size()));
            
            const int num_cells = static_cast<int>(rho_avg_realizations[0].size());
            const double weight = double(1)/double(num_realizations);
            
            std::vector<double> rho_avg_global(num_cells, double(0));
            std::vector<double> rho_u_avg_global(num_cells, double(0));
            std::vector<double> rho_v_avg_global(num_cells, double(0));
            std::vector<double> rho_u_u_avg_global(num_cells, double(0));
            std::vector<double> rho_v_v_avg_global(num_cells, double(0));
            
            for (int ri = 0; ri < num_realizations; ri++)
            {
                for (int i = 0; i < num_cells; i++)
                {
                    rho_avg_global[i]     += weight*rho_avg_realizations[ri][i];
                    rho_u_avg_global[i]   += weight*rho_u_avg_realizations[ri][i];
                    rho_v_avg_global[i]   += weight*rho_v_avg_realizations[ri][i];
                    rho_u_u_avg_global[i] += weight*rho_u_u_avg_realizations[ri][i];
                    rho_v_v_avg_global[i] += weight*rho_v_v_avg_realizations[ri][i];
                }
            }
            
            for (int i = 0; i < num_cells; i++)
            {
                const double u_tilde = rho_u_avg_global[i]/rho_avg_global[i];
                const double v_tilde = rho_v_avg_global[i]/rho_avg_global[i];
                const double rho_u_pp_u_pp = rho_u_u_avg_global[i] - rho_u_avg_global[i]*u_tilde;
                const double rho_v_pp_v_pp = rho_v_v_avg_global[i] - rho_v_avg_global[i]*v_tilde;
                
                TKE_integrated_global += (half*(rho_u_pp_u_pp + rho_v_pp_v_pp)*dx_finest[0]);
            }
            
            const double* x_lo = d_grid_geometry->getXLower();
            const double* x_hi = d_grid_geometry->getXUpper();
            const double L_y = x_hi[1] - x_lo[1];
            TKE_integrated_global *= L_y;
        }
        else if (d_dim == tbox::Dimension(3))
        {
            const std::vector<std::vector<double> >& rho_avg_realizations =
                d_ensemble_statistics->rho_avg_realizations;
            
            const std::vector<std::vector<double> >& rho_u_avg_realizations =
                d_ensemble_statistics->rho_u_avg_realizations;
            
            const std::vector<std::vector<double> >& rho_v_avg_realizations =
                d_ensemble_statistics->rho_u_avg_realizations;
            
            const std::vector<std::vector<double> >& rho_w_avg_realizations =
                d_ensemble_statistics->rho_u_avg_realizations;
            
            const std::vector<std::vector<double> >& rho_u_u_avg_realizations =
                d_ensemble_statistics->rho_u_u_avg_realizations;
            
            const std::vector<std::vector<double> >& rho_v_v_avg_realizations =
                d_ensemble_statistics->rho_u_u_avg_realizations;
            
            const std::vector<std::vector<double> >& rho_w_w_avg_realizations =
                d_ensemble_statistics->rho_u_u_avg_realizations;
            
            const int num_realizations = static_cast<int>(rho_avg_realizations.size());
            
            TBOX_ASSERT(d_ensemble_statistics->getNumberOfEnsembles() == num_realizations);
            TBOX_ASSERT(num_realizations > 0);
            TBOX_ASSERT(num_realizations == static_cast<int>(rho_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(rho_u_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(rho_v_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(rho_w_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(rho_u_u_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(rho_v_v_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(rho_w_w_avg_realizations.size()));
            
            const int num_cells = static_cast<int>(rho_avg_realizations[0].size());
            const double weight = double(1)/double(num_realizations);
            
            std::vector<double> rho_avg_global(num_cells, double(0));
            std::vector<double> rho_u_avg_global(num_cells, double(0));
            std::vector<double> rho_v_avg_global(num_cells, double(0));
            std::vector<double> rho_w_avg_global(num_cells, double(0));
            std::vector<double> rho_u_u_avg_global(num_cells, double(0));
            std::vector<double> rho_v_v_avg_global(num_cells, double(0));
            std::vector<double> rho_w_w_avg_global(num_cells, double(0));
            
            for (int ri = 0; ri < num_realizations; ri++)
            {
                for (int i = 0; i < num_cells; i++)
                {
                    rho_avg_global[i]     += weight*rho_avg_realizations[ri][i];
                    rho_u_avg_global[i]   += weight*rho_u_avg_realizations[ri][i];
                    rho_v_avg_global[i]   += weight*rho_v_avg_realizations[ri][i];
                    rho_w_avg_global[i]   += weight*rho_w_avg_realizations[ri][i];
                    rho_u_u_avg_global[i] += weight*rho_u_u_avg_realizations[ri][i];
                    rho_v_v_avg_global[i] += weight*rho_v_v_avg_realizations[ri][i];
                    rho_w_w_avg_global[i] += weight*rho_w_w_avg_realizations[ri][i];
                }
            }
            
            for (int i = 0; i < num_cells; i++)
            {
                const double u_tilde = rho_u_avg_global[i]/rho_avg_global[i];
                const double v_tilde = rho_v_avg_global[i]/rho_avg_global[i];
                const double w_tilde = rho_w_avg_global[i]/rho_avg_global[i];
                const double rho_u_pp_u_pp = rho_u_u_avg_global[i] - rho_u_avg_global[i]*u_tilde;
                const double rho_v_pp_v_pp = rho_v_v_avg_global[i] - rho_v_avg_global[i]*v_tilde;
                const double rho_w_pp_w_pp = rho_w_w_avg_global[i] - rho_w_avg_global[i]*w_tilde;
                
                TKE_integrated_global += (half*(rho_u_pp_u_pp + rho_v_pp_v_pp + rho_w_pp_w_pp)*dx_finest[0]);
            }
            
            const double* x_lo = d_grid_geometry->getXLower();
            const double* x_hi = d_grid_geometry->getXUpper();
            const double L_y = x_hi[1] - x_lo[1];
            const double L_z = x_hi[2] - x_lo[2];
            TKE_integrated_global *= (L_y*L_z);
        }
        
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << TKE_integrated_global;
        
        f_out.close();
    }
}


/*
 * Output ensemble Reynolds normal stress component in x-direction with assumed homogeneity in y-direction (2D)
 * or yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleReynoldsNormalStressXWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Compute and output the quantity (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_out;
        
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        const std::vector<std::vector<double> >& rho_avg_realizations =
            d_ensemble_statistics->rho_avg_realizations;
        
        const std::vector<std::vector<double> >& rho_u_avg_realizations =
            d_ensemble_statistics->rho_u_avg_realizations;
        
        const std::vector<std::vector<double> >& rho_u_u_avg_realizations =
            d_ensemble_statistics->rho_u_u_avg_realizations;
        
        const std::vector<std::vector<double> >& Y_0_avg_realizations =
            d_ensemble_statistics->Y_0_avg_realizations;
        
        const int num_realizations = static_cast<int>(rho_avg_realizations.size());
        
        TBOX_ASSERT(d_ensemble_statistics->getNumberOfEnsembles() == num_realizations);
        TBOX_ASSERT(num_realizations > 0);
        TBOX_ASSERT(num_realizations == static_cast<int>(rho_avg_realizations.size()));
        TBOX_ASSERT(num_realizations == static_cast<int>(rho_u_avg_realizations.size()));
        TBOX_ASSERT(num_realizations == static_cast<int>(rho_u_u_avg_realizations.size()));
        TBOX_ASSERT(num_realizations == static_cast<int>(Y_0_avg_realizations.size()));
        
        const int num_cells = static_cast<int>(rho_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> rho_avg_global(num_cells, double(0));
        std::vector<double> rho_u_avg_global(num_cells, double(0));
        std::vector<double> rho_u_u_avg_global(num_cells, double(0));
        std::vector<double> Y_0_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                rho_avg_global[i]     += weight*rho_avg_realizations[ri][i];
                rho_u_avg_global[i]   += weight*rho_u_avg_realizations[ri][i];
                rho_u_u_avg_global[i] += weight*rho_u_u_avg_realizations[ri][i];
                Y_0_avg_global[i]     += weight*Y_0_avg_realizations[ri][i];
            }
        }
        
        std::vector<double> R_11(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            const double u_tilde = rho_u_avg_global[i]/rho_avg_global[i];
            const double rho_u_pp_u_pp = rho_u_u_avg_global[i] - rho_u_avg_global[i]*u_tilde;
            R_11[i] = rho_u_pp_u_pp/rho_avg_global[i];
        }
        
        double R_11_sum = double(0);
        int count = 0;
        
        for (int i = 0; i < num_cells; i++)
        {
            const double mixing_metric = double(4)*Y_0_avg_global[i]*(double(1) - Y_0_avg_global[i]);
            if (mixing_metric > double(9)/double(10))
            {
                R_11_sum += R_11[i];
                count++;
            }
        }
        
        const double R_11_mean = R_11_sum/count;
        
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << R_11_mean;
        
        f_out.close();
    }
}


/*
 * Output ensemble Reynolds normal stress component in y-direction with assumed homogeneity in y-direction (2D)
 * or yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleReynoldsNormalStressYWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Compute and output the quantity (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_out;
        
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        const std::vector<std::vector<double> >& rho_avg_realizations =
            d_ensemble_statistics->rho_avg_realizations;
        
        const std::vector<std::vector<double> >& rho_v_avg_realizations =
            d_ensemble_statistics->rho_v_avg_realizations;
        
        const std::vector<std::vector<double> >& rho_v_v_avg_realizations =
            d_ensemble_statistics->rho_v_v_avg_realizations;
        
        const std::vector<std::vector<double> >& Y_0_avg_realizations =
            d_ensemble_statistics->Y_0_avg_realizations;
        
        const int num_realizations = static_cast<int>(rho_avg_realizations.size());
        
        TBOX_ASSERT(d_ensemble_statistics->getNumberOfEnsembles() == num_realizations);
        TBOX_ASSERT(num_realizations > 0);
        TBOX_ASSERT(num_realizations == static_cast<int>(rho_avg_realizations.size()));
        TBOX_ASSERT(num_realizations == static_cast<int>(rho_v_avg_realizations.size()));
        TBOX_ASSERT(num_realizations == static_cast<int>(rho_v_v_avg_realizations.size()));
        TBOX_ASSERT(num_realizations == static_cast<int>(Y_0_avg_realizations.size()));
        
        const int num_cells = static_cast<int>(rho_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> rho_avg_global(num_cells, double(0));
        std::vector<double> rho_v_avg_global(num_cells, double(0));
        std::vector<double> rho_v_v_avg_global(num_cells, double(0));
        std::vector<double> Y_0_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                rho_avg_global[i]     += weight*rho_avg_realizations[ri][i];
                rho_v_avg_global[i]   += weight*rho_v_avg_realizations[ri][i];
                rho_v_v_avg_global[i] += weight*rho_v_v_avg_realizations[ri][i];
                Y_0_avg_global[i]     += weight*Y_0_avg_realizations[ri][i];
            }
        }
        
        std::vector<double> R_22(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            const double v_tilde = rho_v_avg_global[i]/rho_avg_global[i];
            const double rho_v_pp_v_pp = rho_v_v_avg_global[i] - rho_v_avg_global[i]*v_tilde;
            R_22[i] = rho_v_pp_v_pp/rho_avg_global[i];
        }
        
        double R_22_sum = double(0);
        int count = 0;
        
        for (int i = 0; i < num_cells; i++)
        {
            const double mixing_metric = double(4)*Y_0_avg_global[i]*(double(1) - Y_0_avg_global[i]);
            if (mixing_metric > double(9)/double(10))
            {
                R_22_sum += R_22[i];
                count++;
            }
        }
        
        const double R_22_mean = R_22_sum/count;
        
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << R_22_mean;
        
        f_out.close();
    }
}


/*
 * Output ensemble Reynolds normal stress component in z-direction with assumed homogeneity in yz-plane (3D)
 * to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleReynoldsNormalStressZWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Compute and output the quantity (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_out;
        
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        const std::vector<std::vector<double> >& rho_avg_realizations =
            d_ensemble_statistics->rho_avg_realizations;
        
        const std::vector<std::vector<double> >& rho_w_avg_realizations =
            d_ensemble_statistics->rho_w_avg_realizations;
        
        const std::vector<std::vector<double> >& rho_w_w_avg_realizations =
            d_ensemble_statistics->rho_w_w_avg_realizations;
        
        const std::vector<std::vector<double> >& Y_0_avg_realizations =
            d_ensemble_statistics->Y_0_avg_realizations;
        
        const int num_realizations = static_cast<int>(rho_avg_realizations.size());
        
        TBOX_ASSERT(d_ensemble_statistics->getNumberOfEnsembles() == num_realizations);
        TBOX_ASSERT(num_realizations > 0);
        TBOX_ASSERT(num_realizations == static_cast<int>(rho_avg_realizations.size()));
        TBOX_ASSERT(num_realizations == static_cast<int>(rho_w_avg_realizations.size()));
        TBOX_ASSERT(num_realizations == static_cast<int>(rho_w_w_avg_realizations.size()));
        TBOX_ASSERT(num_realizations == static_cast<int>(Y_0_avg_realizations.size()));
        
        const int num_cells = static_cast<int>(rho_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> rho_avg_global(num_cells, double(0));
        std::vector<double> rho_w_avg_global(num_cells, double(0));
        std::vector<double> rho_w_w_avg_global(num_cells, double(0));
        std::vector<double> Y_0_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                rho_avg_global[i]     += weight*rho_avg_realizations[ri][i];
                rho_w_avg_global[i]   += weight*rho_w_avg_realizations[ri][i];
                rho_w_w_avg_global[i] += weight*rho_w_w_avg_realizations[ri][i];
                Y_0_avg_global[i]     += weight*Y_0_avg_realizations[ri][i];
            }
        }
        
        std::vector<double> R_33(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            const double w_tilde = rho_w_avg_global[i]/rho_avg_global[i];
            const double rho_w_pp_w_pp = rho_w_w_avg_global[i] - rho_w_avg_global[i]*w_tilde;
            R_33[i] = rho_w_pp_w_pp/rho_avg_global[i];
        }
        
        double R_33_sum = double(0);
        int count = 0;
        
        for (int i = 0; i < num_cells; i++)
        {
            const double mixing_metric = double(4)*Y_0_avg_global[i]*(double(1) - Y_0_avg_global[i]);
            if (mixing_metric > double(9)/double(10))
            {
                R_33_sum += R_33[i];
                count++;
            }
        }
        
        const double R_33_mean = R_33_sum/count;
        
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << R_33_mean;
        
        f_out.close();
    }
}


/*
 * Output ensemble enstrophy integrated to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleEnstrophyIntegrated(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    MPIHelper MPI_helper = MPIHelper(
        "MPI_helper",
        d_dim,
        d_grid_geometry,
        patch_hierarchy);
    
    const std::vector<double>& dx_finest = MPI_helper.getFinestRefinedDomainGridSpacing();
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_out;
        
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        const std::vector<std::vector<double> >& Omega_avg_realizations =
            d_ensemble_statistics->Omega_avg_realizations;
        
        const int num_realizations = static_cast<int>(Omega_avg_realizations.size());
        
        TBOX_ASSERT(d_ensemble_statistics->getNumberOfEnsembles() == num_realizations);
        TBOX_ASSERT(num_realizations > 0);
        
        const int num_cells = static_cast<int>(Omega_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> Omega_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                Omega_avg_global[i] += weight*Omega_avg_realizations[ri][i];
            }
        }
        
        double Omega_integrated_global = double(0);
        for (int i = 0; i < num_cells; i++)
        {
            Omega_integrated_global += Omega_avg_global[i]*dx_finest[0];
        }
        
        const double* x_lo = d_grid_geometry->getXLower();
        const double* x_hi = d_grid_geometry->getXUpper();
        if (d_dim == tbox::Dimension(2))
        {
            const double L_y = x_hi[1] - x_lo[1];
            Omega_integrated_global *= L_y;
        }
        else if (d_dim == tbox::Dimension(3))
        {
            const double L_y = x_hi[1] - x_lo[1];
            const double L_z = x_hi[2] - x_lo[2];
            Omega_integrated_global *= (L_y*L_z);
        }
        
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << Omega_integrated_global;
    }
}


/*
 * Output ensemble scalar dissipation rate of first species integrated to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleScalarDissipationRateIntegrated(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    MPIHelper MPI_helper = MPIHelper(
        "MPI_helper",
        d_dim,
        d_grid_geometry,
        patch_hierarchy);
    
    const std::vector<double>& dx_finest = MPI_helper.getFinestRefinedDomainGridSpacing();
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_out;
        
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        const std::vector<std::vector<double> >& chi_avg_realizations =
            d_ensemble_statistics->chi_avg_realizations;
        
        const int num_realizations = static_cast<int>(chi_avg_realizations.size());
        
        TBOX_ASSERT(d_ensemble_statistics->getNumberOfEnsembles() == num_realizations);
        TBOX_ASSERT(num_realizations > 0);
        
        const int num_cells = static_cast<int>(chi_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> chi_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                chi_avg_global[i] += weight*chi_avg_realizations[ri][i];
            }
        }
        
        double chi_integrated_global = double(0);
        for (int i = 0; i < num_cells; i++)
        {
            chi_integrated_global += chi_avg_global[i]*dx_finest[0];
        }
        
        const double* x_lo = d_grid_geometry->getXLower();
        const double* x_hi = d_grid_geometry->getXUpper();
        if (d_dim == tbox::Dimension(2))
        {
            const double L_y = x_hi[1] - x_lo[1];
            chi_integrated_global *= L_y;
        }
        else if (d_dim == tbox::Dimension(3))
        {
            const double L_y = x_hi[1] - x_lo[1];
            const double L_z = x_hi[2] - x_lo[2];
            chi_integrated_global *= (L_y*L_z);
        }
        
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << chi_integrated_global;
    }
}


/*
 * Output turbulent Reynolds number based on mixing width with assumed homogeneity in y-direction (2D)
 * or yz-plane (3D) to a file.to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleMixingWidthReynoldsNumberWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'RE_W_INHOMO_X' can be computed with two species only."
            << std::endl);
    }
    
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    MPIHelper MPI_helper = MPIHelper(
        "MPI_helper",
        d_dim,
        d_grid_geometry,
        patch_hierarchy);
    
    const std::vector<double>& dx_finest = MPI_helper.getFinestRefinedDomainGridSpacing();
    
    /*
     * Compute and output the quantity (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_out;
        
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        const std::vector<std::vector<double> >& Y_0_avg_realizations =
            d_ensemble_statistics->Y_0_avg_realizations;
        
        const std::vector<std::vector<double> >& rho_avg_realizations =
            d_ensemble_statistics->rho_avg_realizations;
        
        const std::vector<std::vector<double> >& mu_avg_realizations =
            d_ensemble_statistics->mu_avg_realizations;
        
        const std::vector<std::vector<double> >& u_avg_realizations =
            d_ensemble_statistics->u_sq_avg_realizations;
        
        const std::vector<std::vector<double> >& v_avg_realizations =
            d_ensemble_statistics->v_sq_avg_realizations;
        
        const std::vector<std::vector<double> >& w_avg_realizations =
            d_ensemble_statistics->w_sq_avg_realizations;
        
        const std::vector<std::vector<double> >& u_sq_avg_realizations =
            d_ensemble_statistics->u_sq_avg_realizations;
        
        const std::vector<std::vector<double> >& v_sq_avg_realizations =
            d_ensemble_statistics->v_sq_avg_realizations;
        
        const std::vector<std::vector<double> >& w_sq_avg_realizations =
            d_ensemble_statistics->w_sq_avg_realizations;
        
        const std::vector<std::vector<double> >& rho_u_avg_realizations =
            d_ensemble_statistics->rho_u_avg_realizations;
        
        const std::vector<std::vector<double> >& rho_v_avg_realizations =
            d_ensemble_statistics->rho_v_avg_realizations;
        
        const std::vector<std::vector<double> >& rho_w_avg_realizations =
            d_ensemble_statistics->rho_w_avg_realizations;
        
        const int num_realizations = static_cast<int>(Y_0_avg_realizations.size());
        
        TBOX_ASSERT(d_ensemble_statistics->getNumberOfEnsembles() == num_realizations);
        TBOX_ASSERT(num_realizations > 0);
        TBOX_ASSERT(num_realizations == static_cast<int>(rho_avg_realizations.size()));
        TBOX_ASSERT(num_realizations == static_cast<int>(mu_avg_realizations.size()));
        TBOX_ASSERT(num_realizations == static_cast<int>(u_avg_realizations.size()));
        TBOX_ASSERT(num_realizations == static_cast<int>(u_sq_avg_realizations.size()));
        TBOX_ASSERT(num_realizations == static_cast<int>(rho_u_avg_realizations.size()));
        
        const int num_cells = static_cast<int>(Y_0_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> Y_0_avg_global(num_cells, double(0));
        std::vector<double> rho_avg_global(num_cells, double(0));
        std::vector<double> mu_avg_global(num_cells, double(0));
        std::vector<double> u_avg_global(num_cells, double(0));
        std::vector<double> v_avg_global(num_cells, double(0));
        std::vector<double> w_avg_global(num_cells, double(0));
        std::vector<double> u_sq_avg_global(num_cells, double(0));
        std::vector<double> v_sq_avg_global(num_cells, double(0));
        std::vector<double> w_sq_avg_global(num_cells, double(0));
        std::vector<double> rho_u_avg_global(num_cells, double(0));
        std::vector<double> rho_v_avg_global(num_cells, double(0));
        std::vector<double> rho_w_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                Y_0_avg_global[i]   += weight*Y_0_avg_realizations[ri][i];
                rho_avg_global[i]   += weight*rho_avg_realizations[ri][i];
                mu_avg_global[i]    += weight*mu_avg_realizations[ri][i];
                u_avg_global[i]     += weight*u_avg_realizations[ri][i];
                u_sq_avg_global[i]  += weight*u_sq_avg_realizations[ri][i];
                rho_u_avg_global[i] += weight*rho_u_avg_realizations[ri][i];
            }
        }
        
        if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
        {
            TBOX_ASSERT(num_realizations == static_cast<int>(v_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(v_sq_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(rho_v_avg_realizations.size()));
            
            for (int ri = 0; ri < num_realizations; ri++)
            {
                for (int i = 0; i < num_cells; i++)
                {
                    v_avg_global[i]     += weight*v_avg_realizations[ri][i];
                    v_sq_avg_global[i]  += weight*v_sq_avg_realizations[ri][i];
                    rho_v_avg_global[i] += weight*rho_v_avg_realizations[ri][i];
                }
            }
        }
        
        if (d_dim == tbox::Dimension(3))
        {
            TBOX_ASSERT(num_realizations == static_cast<int>(w_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(w_sq_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(rho_w_avg_realizations.size()));
            
            for (int ri = 0; ri < num_realizations; ri++)
            {
                for (int i = 0; i < num_cells; i++)
                {
                    w_avg_global[i]     += weight*w_avg_realizations[ri][i];
                    w_sq_avg_global[i]  += weight*w_sq_avg_realizations[ri][i];
                    rho_w_avg_global[i] += weight*rho_w_avg_realizations[ri][i];
                }
            }
        }
        
        const double two = double(2);
        std::vector<double> u_i_pp_u_i_pp(num_cells, double(0));
        
        for (int i = 0; i < num_cells; i++)
        {
            const double u_tilde = rho_u_avg_global[i]/rho_avg_global[i];
            const double u_pp_sq = u_sq_avg_global[i] + u_tilde*u_tilde - two*u_tilde*u_avg_global[i];
            u_i_pp_u_i_pp[i] += u_pp_sq;
        }
        
        if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
        {
            for (int i = 0; i < num_cells; i++)
            {
                const double v_tilde = rho_v_avg_global[i]/rho_avg_global[i];
                const double v_pp_sq = v_sq_avg_global[i] + v_tilde*v_tilde - two*v_tilde*v_avg_global[i];
                u_i_pp_u_i_pp[i] += v_pp_sq;
            }
        }
        
        if (d_dim == tbox::Dimension(3))
        {
            for (int i = 0; i < num_cells; i++)
            {
                const double w_tilde = rho_w_avg_global[i]/rho_avg_global[i];
                const double w_pp_sq = w_sq_avg_global[i] + w_tilde*w_tilde - two*w_tilde*w_avg_global[i];
                u_i_pp_u_i_pp[i] += w_pp_sq;
            }
        }
        
        double W = double(0);
        
        for (int i = 0; i < num_cells; i++)
        {
            W += Y_0_avg_global[i]*(double(1) - Y_0_avg_global[i]);
        }
        
        W = double(4)*W*dx_finest[0];
        
        double rho_sum           = double(0);
        double mu_sum            = double(0);
        double u_i_pp_u_i_pp_sum = double(0);
        int count = 0;
        
        for (int i = 0; i < num_cells; i++)
        {
            const double mixing_metric = double(4)*Y_0_avg_global[i]*(double(1) - Y_0_avg_global[i]);
            if (mixing_metric > double(9)/double(10))
            {
                rho_sum           += rho_avg_global[i];
                mu_sum            += mu_avg_global[i];
                u_i_pp_u_i_pp_sum += u_i_pp_u_i_pp[i];
                count++;
            }
        }
        
        const double rho_mean = rho_sum/count;
        const double mu_mean  = mu_sum/count;
        
        double u_rms = double(0);
        if (d_dim == tbox::Dimension(1))
        {
            u_rms = sqrt(u_i_pp_u_i_pp_sum/count);
        }
        else if (d_dim == tbox::Dimension(2))
        {
            const double half = double(1)/double(2);
            u_rms = sqrt(half*u_i_pp_u_i_pp_sum/count);
        }
        else if (d_dim == tbox::Dimension(3))
        {
            const double one_third = double(1)/double(3);
            u_rms = sqrt(one_third*u_i_pp_u_i_pp_sum/count);
        }
        
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << rho_mean*u_rms*W/mu_mean;
        
        f_out.close();
    }
}


/*
 * Output turbulent Mach number with assumed homogeneity in y-direction (2D) or yz-plane (3D) to a file.to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleTurbulentMachNumberWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'MA_T_INHOMO_X' can be computed with two species only."
            << std::endl);
    }
    
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Compute and output the quantity (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_out;
        
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        const std::vector<std::vector<double> >& Y_0_avg_realizations =
            d_ensemble_statistics->Y_0_avg_realizations;
        
        const std::vector<std::vector<double> >& rho_avg_realizations =
            d_ensemble_statistics->rho_avg_realizations;
        
        const std::vector<std::vector<double> >& c_avg_realizations =
            d_ensemble_statistics->c_avg_realizations;
        
        const std::vector<std::vector<double> >& u_avg_realizations =
            d_ensemble_statistics->u_sq_avg_realizations;
        
        const std::vector<std::vector<double> >& v_avg_realizations =
            d_ensemble_statistics->v_sq_avg_realizations;
        
        const std::vector<std::vector<double> >& w_avg_realizations =
            d_ensemble_statistics->w_sq_avg_realizations;
        
        const std::vector<std::vector<double> >& u_sq_avg_realizations =
            d_ensemble_statistics->u_sq_avg_realizations;
        
        const std::vector<std::vector<double> >& v_sq_avg_realizations =
            d_ensemble_statistics->v_sq_avg_realizations;
        
        const std::vector<std::vector<double> >& w_sq_avg_realizations =
            d_ensemble_statistics->w_sq_avg_realizations;
        
        const std::vector<std::vector<double> >& rho_u_avg_realizations =
            d_ensemble_statistics->rho_u_avg_realizations;
        
        const std::vector<std::vector<double> >& rho_v_avg_realizations =
            d_ensemble_statistics->rho_v_avg_realizations;
        
        const std::vector<std::vector<double> >& rho_w_avg_realizations =
            d_ensemble_statistics->rho_w_avg_realizations;
        
        const int num_realizations = static_cast<int>(Y_0_avg_realizations.size());
        
        TBOX_ASSERT(d_ensemble_statistics->getNumberOfEnsembles() == num_realizations);
        TBOX_ASSERT(num_realizations > 0);
        TBOX_ASSERT(num_realizations == static_cast<int>(rho_avg_realizations.size()));
        TBOX_ASSERT(num_realizations == static_cast<int>(c_avg_realizations.size()));
        TBOX_ASSERT(num_realizations == static_cast<int>(u_avg_realizations.size()));
        TBOX_ASSERT(num_realizations == static_cast<int>(u_sq_avg_realizations.size()));
        TBOX_ASSERT(num_realizations == static_cast<int>(rho_u_avg_realizations.size()));
        
        const int num_cells = static_cast<int>(Y_0_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> Y_0_avg_global(num_cells, double(0));
        std::vector<double> rho_avg_global(num_cells, double(0));
        std::vector<double> c_avg_global(num_cells, double(0));
        std::vector<double> u_avg_global(num_cells, double(0));
        std::vector<double> v_avg_global(num_cells, double(0));
        std::vector<double> w_avg_global(num_cells, double(0));
        std::vector<double> u_sq_avg_global(num_cells, double(0));
        std::vector<double> v_sq_avg_global(num_cells, double(0));
        std::vector<double> w_sq_avg_global(num_cells, double(0));
        std::vector<double> rho_u_avg_global(num_cells, double(0));
        std::vector<double> rho_v_avg_global(num_cells, double(0));
        std::vector<double> rho_w_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                Y_0_avg_global[i]   += weight*Y_0_avg_realizations[ri][i];
                rho_avg_global[i]   += weight*rho_avg_realizations[ri][i];
                c_avg_global[i]     += weight*c_avg_realizations[ri][i];
                u_avg_global[i]     += weight*u_avg_realizations[ri][i];
                u_sq_avg_global[i]  += weight*u_sq_avg_realizations[ri][i];
                rho_u_avg_global[i] += weight*rho_u_avg_realizations[ri][i];
            }
        }
        
        if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
        {
            TBOX_ASSERT(num_realizations == static_cast<int>(v_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(v_sq_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(rho_v_avg_realizations.size()));
            
            for (int ri = 0; ri < num_realizations; ri++)
            {
                for (int i = 0; i < num_cells; i++)
                {
                    v_avg_global[i]     += weight*v_avg_realizations[ri][i];
                    v_sq_avg_global[i]  += weight*v_sq_avg_realizations[ri][i];
                    rho_v_avg_global[i] += weight*rho_v_avg_realizations[ri][i];
                }
            }
        }
        
        if (d_dim == tbox::Dimension(3))
        {
            TBOX_ASSERT(num_realizations == static_cast<int>(w_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(w_sq_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(rho_w_avg_realizations.size()));
            
            for (int ri = 0; ri < num_realizations; ri++)
            {
                for (int i = 0; i < num_cells; i++)
                {
                    w_avg_global[i]     += weight*w_avg_realizations[ri][i];
                    w_sq_avg_global[i]  += weight*w_sq_avg_realizations[ri][i];
                    rho_w_avg_global[i] += weight*rho_w_avg_realizations[ri][i];
                }
            }
        }
        
        const double two = double(2);
        std::vector<double> u_i_pp_u_i_pp(num_cells, double(0));
        
        for (int i = 0; i < num_cells; i++)
        {
            const double u_tilde = rho_u_avg_global[i]/rho_avg_global[i];
            const double u_pp_sq = u_sq_avg_global[i] + u_tilde*u_tilde - two*u_tilde*u_avg_global[i];
            u_i_pp_u_i_pp[i] += u_pp_sq;
        }
        
        if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
        {
            for (int i = 0; i < num_cells; i++)
            {
                const double v_tilde = rho_v_avg_global[i]/rho_avg_global[i];
                const double v_pp_sq = v_sq_avg_global[i] + v_tilde*v_tilde - two*v_tilde*v_avg_global[i];
                u_i_pp_u_i_pp[i] += v_pp_sq;
            }
        }
        
        if (d_dim == tbox::Dimension(3))
        {
            for (int i = 0; i < num_cells; i++)
            {
                const double w_tilde = rho_w_avg_global[i]/rho_avg_global[i];
                const double w_pp_sq = w_sq_avg_global[i] + w_tilde*w_tilde - two*w_tilde*w_avg_global[i];
                u_i_pp_u_i_pp[i] += w_pp_sq;
            }
        }
        
        double c_sum             = double(0);
        double u_i_pp_u_i_pp_sum = double(0);
        int count = 0;
        
        for (int i = 0; i < num_cells; i++)
        {
            const double mixing_metric = double(4)*Y_0_avg_global[i]*(double(1) - Y_0_avg_global[i]);
            if (mixing_metric > double(9)/double(10))
            {
                c_sum             += c_avg_global[i];
                u_i_pp_u_i_pp_sum += u_i_pp_u_i_pp[i];
                count++;
            }
        }
        
        const double u_rms  = sqrt(u_i_pp_u_i_pp_sum/count);
        const double c_mean = c_sum/count;
        
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << u_rms/c_mean;
        
        f_out.close();
    }
}


/*
 * Output names of statistical quantities to output to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputStatisticalQuantitiesNames(
    const std::string& stat_dump_filename)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    if (mpi.getRank() == 0)
    {
        std::ofstream f_out;
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
        
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        // Loop over statistical quantities.
        for (int qi = 0; qi < static_cast<int>(d_statistical_quantities.size()); qi++)
        {
            // Get the key of the current variable.
            std::string statistical_quantity_key = d_statistical_quantities[qi];
            
            if (statistical_quantity_key == "MIXING_WIDTH_X")
            {
                f_out << "\t" << "MIXING_WIDTH_X       ";
            }
            else if (statistical_quantity_key == "INTERFACE_MIN_X")
            {
                f_out << "\t" << "INTERFACE_MIN_X      ";
            }
            else if (statistical_quantity_key == "INTERFACE_MAX_X")
            {
                f_out << "\t" << "INTERFACE_MAX_X      ";
            }
            else if (statistical_quantity_key == "MIXEDNESS_X")
            {
                f_out << "\t" << "MIXEDNESS_X          ";
            }
            else if (statistical_quantity_key == "TKE_INT")
            {
                f_out << "\t" << "TKE_INT              ";
            }
            else if (statistical_quantity_key == "a1_MEAN_INHOMO_X")
            {
                f_out << "\t" << "a1_MEAN_INHOMO_X     ";
            }
            else if (statistical_quantity_key == "b_MEAN_INHOMO_X")
            {
                f_out << "\t" << "b_MEAN_INHOMO_X      ";
            }
            else if (statistical_quantity_key == "R11_MEAN_INHOMO_X")
            {
                f_out << "\t" << "R11_MEAN_INHOMO_X    ";
            }
            else if (statistical_quantity_key == "R22_MEAN_INHOMO_X")
            {
                f_out << "\t" << "R22_MEAN_INHOMO_X    ";
            }
            else if (statistical_quantity_key == "R33_MEAN_INHOMO_X")
            {
                f_out << "\t" << "R33_MEAN_INHOMO_X    ";
            }
            else if (statistical_quantity_key == "ENSTROPHY_INT")
            {
                f_out << "\t" << "ENSTROPHY_INT        ";
            }
            else if (statistical_quantity_key == "SCAL_DISS_RAT_INT")
            {
                f_out << "\t" << "SCAL_DISS_RAT_INT    ";
            }
            else if (statistical_quantity_key == "RE_W_INHOMO_X")
            {
                f_out << "\t" << "RE_W_INHOMO_X        ";
            }
            else if (statistical_quantity_key == "MA_T_INHOMO_X")
            {
                f_out << "\t" << "MA_T_INHOMO_X        ";
            }
        }
        
        f_out.close();
    }
}


/*
 * Compute statisitcal quantities.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::computeStatisticalQuantities(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double statistics_data_time)
{
    NULL_USE(statistics_data_time);
    
    if (d_is_ensemble_statistics_initialized == false)
    {
        d_ensemble_statistics = HAMERS_MAKE_SHARED<EnsembleStatisticsRTIRMI>("d_ensemble_statistics");
        d_is_ensemble_statistics_initialized = true;
    }
    
    HAMERS_SHARED_PTR<RTIRMIStatisticsUtilities> rti_rmi_statistics_utilities(
        new RTIRMIStatisticsUtilities(
            "RTI RMI spatial profiles utilities",
            d_dim,
            d_grid_geometry,
            d_num_species,
            d_flow_model,
            d_equation_of_state_mixing_rules,
            d_equation_of_mass_diffusivity_mixing_rules,
            d_equation_of_shear_viscosity_mixing_rules,
            d_equation_of_bulk_viscosity_mixing_rules,
            d_equation_of_thermal_conductivity_mixing_rules,
            HAMERS_DYNAMIC_POINTER_CAST<EnsembleStatisticsRTIRMI>(d_ensemble_statistics)));
    
    // Statistics are not computed for this realization yet.
    rti_rmi_statistics_utilities->d_ensemble_statistics->setVariablesNotComputed();
    
    // Loop over statistical quantities.
    for (int qi = 0; qi < static_cast<int>(d_statistical_quantities.size()); qi++)
    {
        // Get the key of the current variable.
        std::string statistical_quantity_key = d_statistical_quantities[qi];
        
        // Spatial profiles.
        if (statistical_quantity_key == "MASS_FRACTION_AVG_SP")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->Y_0_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedMassFractionWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
        }
        else if (statistical_quantity_key == "MASS_FRACTION_VAR_SP")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->Y_0_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedMassFractionWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
            
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->Y_0_sq_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeMassFractionVarianceWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
        }
        else if (statistical_quantity_key == "DENSITY_AVG_SP")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->rho_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedDensityWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
        }
        else if (statistical_quantity_key == "DENSITY_VAR_SP")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->rho_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedDensityWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
            
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->rho_sq_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeDensityVarianceWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
        }
        else if (statistical_quantity_key == "SPECIFIC_VOLUME_AVG_SP")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->rho_inv_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedSpecificVolumeWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
        }
        else if (statistical_quantity_key == "DENSITY_SPEC_VOL_COV_SP")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->rho_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedDensityWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
            
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->rho_inv_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedSpecificVolumeWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
        }
        else if (statistical_quantity_key == "VELOCITY_X_AVG_SP")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->u_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedVelocityXWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
        }
        else if (statistical_quantity_key == "VELOCITY_Y_AVG_SP")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->v_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedVelocityYWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
        }
        else if (statistical_quantity_key == "VELOCITY_Z_AVG_SP")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->w_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedVelocityZWithHomogeneityInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
        }
        else if (statistical_quantity_key == "MOMENTUM_X_AVG_SP")
        {
            rti_rmi_statistics_utilities->
                computeAveragedMomentumXWithHomogeneityInYDirectionOrInYZPlane(
                    patch_hierarchy,
                    data_context);
        }
        else if (statistical_quantity_key == "MOMENTUM_Y_AVG_SP")
        {
            rti_rmi_statistics_utilities->
                computeAveragedMomentumYWithHomogeneityInYDirectionOrInYZPlane(
                    patch_hierarchy,
                    data_context);
        }
        else if (statistical_quantity_key == "MOMENTUM_Z_AVG_SP")
        {
            rti_rmi_statistics_utilities->
                computeAveragedMomentumZWithHomogeneityInYZPlane(
                    patch_hierarchy,
                    data_context);
        }
        else if (statistical_quantity_key == "TURB_MASS_FLUX_X_SP")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->rho_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedDensityWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
            
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->u_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedVelocityXWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
            
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->rho_u_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedMomentumXWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
        }
        else if (statistical_quantity_key == "R11_SP")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->rho_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedDensityWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
            
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->rho_u_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedMomentumXWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
            
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->rho_u_u_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeReynoldsNormalStressXWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
        }
        else if (statistical_quantity_key == "R22_SP")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->rho_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedDensityWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
            
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->rho_v_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedMomentumYWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
            
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->rho_v_v_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeReynoldsNormalStressYWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
        }
        else if (statistical_quantity_key == "R33_SP")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->rho_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedDensityWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
            
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->rho_w_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedMomentumZWithHomogeneityInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
            
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->rho_w_w_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeReynoldsNormalStressZWithHomogeneityInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
        }
        else if (statistical_quantity_key == "ENSTROPHY_AVG_SP")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->Omega_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedEnstrophyWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
        }
        else if (statistical_quantity_key == "SCAL_DISS_RAT_AVG_SP")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->chi_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedScalarDissipationRateWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
        }
        // Non-spatial profiles.
        else if (statistical_quantity_key == "MIXING_WIDTH_X")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->Y_0_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedMassFractionWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
        }
        else if (statistical_quantity_key == "INTERFACE_MIN_X")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->Y_0_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedMassFractionWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
        }
        else if (statistical_quantity_key == "INTERFACE_MAX_X")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->Y_0_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedMassFractionWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
        }
        else if (statistical_quantity_key == "MIXEDNESS_X")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->Y_0_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedMassFractionWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
            
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->Y_0_Y_1_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeMassFractionProductWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
        }
        else if (statistical_quantity_key == "TKE_INT")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->rho_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedDensityWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
            
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->rho_u_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedMomentumXWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
            
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->rho_u_u_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeReynoldsNormalStressXWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
            
            if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
            {
                if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->rho_v_avg_computed))
                {
                    rti_rmi_statistics_utilities->
                        computeAveragedMomentumYWithHomogeneityInYDirectionOrInYZPlane(
                            patch_hierarchy,
                            data_context);
                }
                
                if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->rho_v_v_avg_computed))
                {
                    rti_rmi_statistics_utilities->
                        computeReynoldsNormalStressYWithHomogeneityInYDirectionOrInYZPlane(
                            patch_hierarchy,
                            data_context);
                }
            }
            
            if (d_dim == tbox::Dimension(3))
            {
                if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->rho_w_avg_computed))
                {
                    rti_rmi_statistics_utilities->
                        computeAveragedMomentumZWithHomogeneityInYZPlane(
                            patch_hierarchy,
                            data_context);
                }
                
                if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->rho_w_w_avg_computed))
                {
                    rti_rmi_statistics_utilities->
                        computeReynoldsNormalStressZWithHomogeneityInYZPlane(
                            patch_hierarchy,
                            data_context);
                }
            }
        }
        else if (statistical_quantity_key == "a1_MEAN_INHOMO_X")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->rho_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedDensityWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
            
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->u_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedVelocityXWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
            
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->rho_u_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedMomentumXWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
            
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->Y_0_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedMassFractionWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
        }
        else if (statistical_quantity_key == "b_MEAN_INHOMO_X")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->rho_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedDensityWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
            
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->rho_inv_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedSpecificVolumeWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
        }
        else if (statistical_quantity_key == "R11_MEAN_INHOMO_X")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->rho_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedDensityWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
            
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->rho_u_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedMomentumXWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
            
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->rho_u_u_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeReynoldsNormalStressXWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
            
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->Y_0_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedMassFractionWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
        }
        else if (statistical_quantity_key == "R22_MEAN_INHOMO_X")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->rho_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedDensityWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
            
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->rho_v_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedMomentumYWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
            
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->rho_v_v_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeReynoldsNormalStressYWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
            
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->Y_0_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedMassFractionWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
        }
        else if (statistical_quantity_key == "R33_MEAN_INHOMO_X")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->rho_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedDensityWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
            
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->rho_w_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedMomentumZWithHomogeneityInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
            
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->rho_w_w_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeReynoldsNormalStressZWithHomogeneityInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
            
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->Y_0_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedMassFractionWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
        }
        else if (statistical_quantity_key == "ENSTROPHY_INT")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->Omega_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedEnstrophyWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
        }
        else if (statistical_quantity_key == "SCAL_DISS_RAT_INT")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->chi_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedScalarDissipationRateWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
        }
        else if (statistical_quantity_key == "RE_W_INHOMO_X")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->Y_0_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedMassFractionWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
            
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->rho_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedDensityWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
            
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->mu_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedShearViscosityWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
            
            // The following are for computing rms of velocity.
            
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->u_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedVelocityXWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
            
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->u_sq_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeVelocityXVarianceWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
            
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->rho_u_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedMomentumXWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
            
            if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
            {
                if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->v_avg_computed))
                {
                    rti_rmi_statistics_utilities->
                        computeAveragedVelocityYWithHomogeneityInYDirectionOrInYZPlane(
                            patch_hierarchy,
                            data_context);
                }
                
                if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->v_sq_avg_computed))
                {
                    rti_rmi_statistics_utilities->
                        computeVelocityYVarianceWithHomogeneityInYDirectionOrInYZPlane(
                            patch_hierarchy,
                            data_context);
                }
                
                if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->rho_v_avg_computed))
                {
                    rti_rmi_statistics_utilities->
                        computeAveragedMomentumYWithHomogeneityInYDirectionOrInYZPlane(
                            patch_hierarchy,
                            data_context);
                }
            }
            
            if (d_dim == tbox::Dimension(3))
            {
                if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->w_avg_computed))
                {
                    rti_rmi_statistics_utilities->
                        computeAveragedVelocityZWithHomogeneityInYZPlane(
                            patch_hierarchy,
                            data_context);
                }
                
                if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->w_sq_avg_computed))
                {
                    rti_rmi_statistics_utilities->
                        computeVelocityZVarianceWithHomogeneityInYZPlane(
                            patch_hierarchy,
                            data_context);
                }
                
                if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->rho_w_avg_computed))
                {
                    rti_rmi_statistics_utilities->
                        computeAveragedMomentumZWithHomogeneityInYZPlane(
                            patch_hierarchy,
                            data_context);
                }
            }
        }
        else if (statistical_quantity_key == "MA_T_INHOMO_X")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->Y_0_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedMassFractionWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
            
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->rho_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedDensityWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
            
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->c_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedSoundSpeedWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
            
            // The following are for computing rms of velocity.
            
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->u_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedVelocityXWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
            
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->u_sq_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeVelocityXVarianceWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
            
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->rho_u_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedMomentumXWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
            
            if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
            {
                if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->v_avg_computed))
                {
                    rti_rmi_statistics_utilities->
                        computeAveragedVelocityYWithHomogeneityInYDirectionOrInYZPlane(
                            patch_hierarchy,
                            data_context);
                }
                
                if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->v_sq_avg_computed))
                {
                    rti_rmi_statistics_utilities->
                        computeVelocityYVarianceWithHomogeneityInYDirectionOrInYZPlane(
                            patch_hierarchy,
                            data_context);
                }
                
                if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->rho_v_avg_computed))
                {
                    rti_rmi_statistics_utilities->
                        computeAveragedMomentumYWithHomogeneityInYDirectionOrInYZPlane(
                            patch_hierarchy,
                            data_context);
                }
            }
            
            if (d_dim == tbox::Dimension(3))
            {
                if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->w_avg_computed))
                {
                    rti_rmi_statistics_utilities->
                        computeAveragedVelocityZWithHomogeneityInYZPlane(
                            patch_hierarchy,
                            data_context);
                }
                
                if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->w_sq_avg_computed))
                {
                    rti_rmi_statistics_utilities->
                        computeVelocityZVarianceWithHomogeneityInYZPlane(
                            patch_hierarchy,
                            data_context);
                }
                
                if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->rho_w_avg_computed))
                {
                    rti_rmi_statistics_utilities->
                        computeAveragedMomentumZWithHomogeneityInYZPlane(
                            patch_hierarchy,
                            data_context);
                }
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Unknown statistical quantity key = '"
                << statistical_quantity_key
                << "' found."
                << std::endl);
        }
    }
    
    d_ensemble_statistics->incrementNumberOfEnsembles();
}


/*
 * Output statisitcal quantities to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputStatisticalQuantities(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double output_time)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_is_ensemble_statistics_initialized);
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    HAMERS_SHARED_PTR<RTIRMIStatisticsUtilities> rti_rmi_statistics_utilities(
        new RTIRMIStatisticsUtilities(
            "RTI RMI spatial profiles utilities",
            d_dim,
            d_grid_geometry,
            d_num_species,
            d_flow_model,
            d_equation_of_state_mixing_rules,
            d_equation_of_mass_diffusivity_mixing_rules,
            d_equation_of_shear_viscosity_mixing_rules,
            d_equation_of_bulk_viscosity_mixing_rules,
            d_equation_of_thermal_conductivity_mixing_rules,
            HAMERS_DYNAMIC_POINTER_CAST<EnsembleStatisticsRTIRMI>(d_ensemble_statistics)));
    
    // Loop over statistical quantities.
    for (int qi = 0; qi < static_cast<int>(d_statistical_quantities.size()); qi++)
    {
        // Get the key of the current variable.
        std::string statistical_quantity_key = d_statistical_quantities[qi];
        
        // Spatial profiles.
        if (statistical_quantity_key == "MASS_FRACTION_AVG_SP")
        {
            rti_rmi_statistics_utilities->
                outputSpatialProfileEnsembleAveragedMassFractionWithHomogeneityInYDirectionOrInYZPlane(
                    "Y_avg.dat",
                    patch_hierarchy,
                    data_context,
                    output_time);
        }
        else if (statistical_quantity_key == "MASS_FRACTION_VAR_SP")
        {
            rti_rmi_statistics_utilities->
                outputSpatialProfileEnsembleMassFractionVarianceWithHomogeneityInYDirectionOrInYZPlane(
                    "Y_var.dat",
                    patch_hierarchy,
                    data_context,
                    output_time);
        }
        else if (statistical_quantity_key == "DENSITY_AVG_SP")
        {
            rti_rmi_statistics_utilities->
                outputSpatialProfileEnsembleAveragedDensityWithHomogeneityInYDirectionOrInYZPlane(
                    "rho_avg.dat",
                    patch_hierarchy,
                    data_context,
                    output_time);
        }
        else if (statistical_quantity_key == "DENSITY_VAR_SP")
        {
            rti_rmi_statistics_utilities->
                outputSpatialProfileEnsembleDensityVarianceWithHomogeneityInYDirectionOrInYZPlane(
                    "rho_var.dat",
                    patch_hierarchy,
                    data_context,
                    output_time);
        }
        else if (statistical_quantity_key == "SPECIFIC_VOLUME_AVG_SP")
        {
            rti_rmi_statistics_utilities->
                outputSpatialProfileEnsembleAveragedSpecificVolumeWithHomogeneityInYDirectionOrInYZPlane(
                    "rho_inv_avg.dat",
                    patch_hierarchy,
                    data_context,
                    output_time);
        }
        else if (statistical_quantity_key == "DENSITY_SPEC_VOL_COV_SP")
        {
            rti_rmi_statistics_utilities->
                outputSpatialProfileEnsembleAveragedDensitySpecificVolumeCovarianceWithHomogeneityInYDirectionOrInYZPlane(
                    "b.dat",
                    patch_hierarchy,
                    data_context,
                    output_time);
        }
        else if (statistical_quantity_key == "VELOCITY_X_AVG_SP")
        {
            rti_rmi_statistics_utilities->
                outputSpatialProfileEnsembleAveragedVelocityXWithHomogeneityInYDirectionOrInYZPlane(
                    "u_avg.dat",
                    patch_hierarchy,
                    data_context,
                    output_time);
        }
        else if (statistical_quantity_key == "VELOCITY_Y_AVG_SP")
        {
            rti_rmi_statistics_utilities->
                outputSpatialProfileEnsembleAveragedVelocityYWithHomogeneityInYDirectionOrInYZPlane(
                    "v_avg.dat",
                    patch_hierarchy,
                    data_context,
                    output_time);
        }
        else if (statistical_quantity_key == "VELOCITY_Z_AVG_SP")
        {
            rti_rmi_statistics_utilities->
                outputSpatialProfileEnsembleAveragedVelocityZWithHomogeneityInYZPlane(
                    "w_avg.dat",
                    patch_hierarchy,
                    data_context,
                    output_time);
        }
        else if (statistical_quantity_key == "MOMENTUM_X_AVG_SP")
        {
            rti_rmi_statistics_utilities->
                outputSpatialProfileEnsembleAveragedMomentumXWithHomogeneityInYDirectionOrInYZPlane(
                    "rho_u_avg.dat",
                    patch_hierarchy,
                    data_context,
                    output_time);
        }
        else if (statistical_quantity_key == "MOMENTUM_Y_AVG_SP")
        {
            rti_rmi_statistics_utilities->
                outputSpatialProfileEnsembleAveragedMomentumYWithHomogeneityInYDirectionOrInYZPlane(
                    "rho_v_avg.dat",
                    patch_hierarchy,
                    data_context,
                    output_time);
        }
        else if (statistical_quantity_key == "MOMENTUM_Z_AVG_SP")
        {
            rti_rmi_statistics_utilities->
                outputSpatialProfileEnsembleAveragedMomentumZWithHomogeneitynYZPlane(
                    "rho_w_avg.dat",
                    patch_hierarchy,
                    data_context,
                    output_time);
        }
        else if (statistical_quantity_key == "TURB_MASS_FLUX_X_SP")
        {
            rti_rmi_statistics_utilities->
                outputSpatialProfileEnsembleTurbulentMassFluxXWithHomogeneityInYDirectionOrInYZPlane(
                    "rho_a1.dat",
                    patch_hierarchy,
                    data_context,
                    output_time);
        }
        else if (statistical_quantity_key == "R11_SP")
        {
            rti_rmi_statistics_utilities->
                outputSpatialProfileEnsembleReynoldsNormalStressXWithHomogeneityInYDirectionOrInYZPlane(
                    "R11.dat",
                    patch_hierarchy,
                    data_context,
                    output_time);
        }
        else if (statistical_quantity_key == "R22_SP")
        {
            rti_rmi_statistics_utilities->
                outputSpatialProfileEnsembleReynoldsNormalStressYWithHomogeneityInYDirectionOrInYZPlane(
                    "R22.dat",
                    patch_hierarchy,
                    data_context,
                    output_time);
        }
        else if (statistical_quantity_key == "R33_SP")
        {
            rti_rmi_statistics_utilities->
                outputSpatialProfileEnsembleReynoldsNormalStressZWithHomogeneityInYZPlane(
                    "R33.dat",
                    patch_hierarchy,
                    data_context,
                    output_time);
        }
        else if (statistical_quantity_key == "ENSTROPHY_AVG_SP")
        {
            rti_rmi_statistics_utilities->
                outputSpatialProfileEnsembleAveragedEnstrophyWithHomogeneityInYDirectionOrInYZPlane(
                    "Omega_avg.dat",
                    patch_hierarchy,
                    data_context,
                    output_time);
        }
        else if (statistical_quantity_key == "SCAL_DISS_RAT_AVG_SP")
        {
            rti_rmi_statistics_utilities->
                outputSpatialProfileEnsembleAveragedScalarDissipationRateWithHomogeneityInYDirectionOrInYZPlane(
                    "chi_avg.dat",
                    patch_hierarchy,
                    data_context,
                    output_time);
        }
        // Non-spatial profiles.
        else if (statistical_quantity_key == "MIXING_WIDTH_X")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleMixingWidthInXDirection(
                    stat_dump_filename,
                    patch_hierarchy,
                    data_context);
        }
        else if (statistical_quantity_key == "INTERFACE_MIN_X")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleInterfaceMinInXDirection(
                    stat_dump_filename,
                    patch_hierarchy,
                    data_context);
        }
        else if (statistical_quantity_key == "INTERFACE_MAX_X")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleInterfaceMaxInXDirection(
                    stat_dump_filename,
                    patch_hierarchy,
                    data_context);
        }
        else if (statistical_quantity_key == "MIXEDNESS_X")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleMixednessInXDirection(
                    stat_dump_filename,
                    patch_hierarchy,
                    data_context);
        }
        else if (statistical_quantity_key == "TKE_INT")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleTKEIntegrateWithInhomogeneousXDirection(
                    stat_dump_filename,
                    patch_hierarchy,
                    data_context);
        }
        else if (statistical_quantity_key == "a1_MEAN_INHOMO_X")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleTurbulentMassFluxVelocityXWithInhomogeneousXDirection(
                    stat_dump_filename,
                    patch_hierarchy,
                    data_context);
        }
        else if (statistical_quantity_key == "b_MEAN_INHOMO_X")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleDensitySpecificVolumeCovarianceWithInhomogeneousXDirection(
                    stat_dump_filename,
                    patch_hierarchy,
                    data_context);
        }
        else if (statistical_quantity_key == "R11_MEAN_INHOMO_X")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleReynoldsNormalStressXWithInhomogeneousXDirection(
                    stat_dump_filename,
                    patch_hierarchy,
                    data_context);
        }
        else if (statistical_quantity_key == "R22_MEAN_INHOMO_X")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleReynoldsNormalStressYWithInhomogeneousXDirection(
                    stat_dump_filename,
                    patch_hierarchy,
                    data_context);
        }
        else if (statistical_quantity_key == "R33_MEAN_INHOMO_X")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleReynoldsNormalStressZWithInhomogeneousXDirection(
                    stat_dump_filename,
                    patch_hierarchy,
                    data_context);
        }
        else if (statistical_quantity_key == "ENSTROPHY_INT")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleEnstrophyIntegrated(
                    stat_dump_filename,
                    patch_hierarchy,
                    data_context);
        }
        else if (statistical_quantity_key == "SCAL_DISS_RAT_INT")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleScalarDissipationRateIntegrated(
                    stat_dump_filename,
                    patch_hierarchy,
                    data_context);
        }
        else if (statistical_quantity_key == "RE_W_INHOMO_X")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleMixingWidthReynoldsNumberWithInhomogeneousXDirection(
                    stat_dump_filename,
                    patch_hierarchy,
                    data_context);
        }
        else if (statistical_quantity_key == "MA_T_INHOMO_X")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleTurbulentMachNumberWithInhomogeneousXDirection(
                    stat_dump_filename,
                    patch_hierarchy,
                    data_context);
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Unknown statistical quantity key = '"
                << statistical_quantity_key
                << "' found."
                << std::endl);
        }
    }
}
