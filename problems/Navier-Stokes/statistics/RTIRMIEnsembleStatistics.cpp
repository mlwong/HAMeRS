#include "flow/flow_models/four-eqn_conservative/FlowModelStatisticsUtilitiesFourEqnConservative.hpp"

#include "extn/patch_hierarchies/ExtendedFlattenedHierarchy.hpp"
#include "flow/flow_models/MPI_helpers/FlowModelMPIHelperAverage.hpp"
#include "flow/flow_models/MPI_helpers/FlowModelMPIHelperCorrelation.hpp"
#include "flow/flow_models/MPI_helpers/FlowModelMPIHelperMaxMin.hpp"
#include "util/MPI_helpers/MPIHelperGrid.hpp"

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
            X_0_avg_computed      = false;
            X_0_sq_avg_computed   = false;
            X_0_X_1_avg_computed  = false;
            Y_0_avg_computed      = false;
            Y_0_sq_avg_computed   = false;
            Y_0_Y_1_avg_computed  = false;
            Z_0_avg_computed      = false;
            Z_0_sq_avg_computed   = false;
            Z_0_Z_1_avg_computed  = false;
            rho_avg_computed      = false;
            rho_sq_avg_computed   = false;
            rho_inv_avg_computed  = false;
            p_avg_computed        = false;
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
            D_avg_computed        = false;
            
            ddx_tau11_avg_computed = false;
            
            // For computing baroclinic torque.
            
            B1_sq_avg_computed = false;
            B2_sq_avg_computed = false;
            B3_sq_avg_computed = false;
            
            // For computing TKE dissipation.
            
            TKE_diss_computed = false;
            
            ddx_p_avg_computed = false;
            ddy_p_avg_computed = false;
            ddz_p_avg_computed = false;
            
            ddx_u_avg_computed = false;
            ddy_u_avg_computed = false;
            ddz_u_avg_computed = false;
            
            ddx_v_avg_computed = false;
            ddy_v_avg_computed = false;
            ddz_v_avg_computed = false;
            
            ddx_w_avg_computed = false;
            ddy_w_avg_computed = false;
            ddz_w_avg_computed = false;
            
            tau11_avg_computed = false;
            tau12_avg_computed = false;
            tau13_avg_computed = false;
            tau22_avg_computed = false;
            tau23_avg_computed = false;
            tau33_avg_computed = false;
            
            tau11_ddx_u_avg_computed = false;
            tau12_ddy_u_avg_computed = false;
            tau13_ddz_u_avg_computed = false;
            
            tau12_ddx_v_avg_computed = false;
            tau22_ddy_v_avg_computed = false;
            tau23_ddz_v_avg_computed = false;
            
            tau13_ddx_w_avg_computed = false;
            tau23_ddy_w_avg_computed = false;
            tau33_ddz_w_avg_computed = false;
        }
        
        void clearAllData()
        {
            X_0_avg_realizations.clear();
            X_0_sq_avg_realizations.clear();
            X_0_X_1_avg_realizations.clear();
            Y_0_avg_realizations.clear();
            Y_0_sq_avg_realizations.clear();
            Y_0_Y_1_avg_realizations.clear();
            Z_0_avg_realizations.clear();
            Z_0_sq_avg_realizations.clear();
            Z_0_Z_1_avg_realizations.clear();
            rho_avg_realizations.clear();
            rho_sq_avg_realizations.clear();
            rho_inv_avg_realizations.clear();
            p_avg_realizations.clear();
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
            D_avg_realizations.clear();
            
            ddx_p_avg_realizations.clear();
            ddy_p_avg_realizations.clear();
            ddz_p_avg_realizations.clear();
            
            ddx_tau11_avg_realizations.clear();
            
            // For computing baroclinic torque.
            
            B1_sq_avg_realizations.clear();
            B2_sq_avg_realizations.clear();
            B3_sq_avg_realizations.clear();
            
            // For computing TKE dissipation.
            
            ddx_u_avg_realizations.clear();
            ddy_u_avg_realizations.clear();
            ddz_u_avg_realizations.clear();
            
            ddx_v_avg_realizations.clear();
            ddy_v_avg_realizations.clear();
            ddz_v_avg_realizations.clear();
            
            ddx_w_avg_realizations.clear();
            ddy_w_avg_realizations.clear();
            ddz_w_avg_realizations.clear();
            
            tau11_avg_realizations.clear();
            tau12_avg_realizations.clear();
            tau13_avg_realizations.clear();
            tau22_avg_realizations.clear();
            tau23_avg_realizations.clear();
            tau33_avg_realizations.clear();
            
            tau11_ddx_u_avg_realizations.clear();
            tau12_ddy_u_avg_realizations.clear();
            tau13_ddz_u_avg_realizations.clear();
            
            tau12_ddx_v_avg_realizations.clear();
            tau22_ddy_v_avg_realizations.clear();
            tau23_ddz_v_avg_realizations.clear();
            
            tau13_ddx_w_avg_realizations.clear();
            tau23_ddy_w_avg_realizations.clear();
            tau33_ddz_w_avg_realizations.clear();
            
            setVariablesNotComputed();
        }
        
        // Scratch arrays.
        // Number of realizalizations; number of cells.
        std::vector<std::vector<double> > X_0_avg_realizations;
        std::vector<std::vector<double> > X_0_sq_avg_realizations;
        std::vector<std::vector<double> > X_0_X_1_avg_realizations;
        std::vector<std::vector<double> > Y_0_avg_realizations;
        std::vector<std::vector<double> > Y_0_sq_avg_realizations;
        std::vector<std::vector<double> > Y_0_Y_1_avg_realizations;
        std::vector<std::vector<double> > Z_0_avg_realizations;
        std::vector<std::vector<double> > Z_0_sq_avg_realizations;
        std::vector<std::vector<double> > Z_0_Z_1_avg_realizations;
        std::vector<std::vector<double> > rho_avg_realizations;
        std::vector<std::vector<double> > rho_sq_avg_realizations;
        std::vector<std::vector<double> > rho_inv_avg_realizations;
        std::vector<std::vector<double> > p_avg_realizations;
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
        std::vector<std::vector<double> > D_avg_realizations;
        
        std::vector<std::vector<double> > ddx_p_avg_realizations;
        std::vector<std::vector<double> > ddy_p_avg_realizations;
        std::vector<std::vector<double> > ddz_p_avg_realizations;
        
        std::vector<std::vector<double> > ddx_tau11_avg_realizations;
        
        // For computing baroclinic torque.
        
        std::vector<std::vector<double> > B1_sq_avg_realizations;
        std::vector<std::vector<double> > B2_sq_avg_realizations;
        std::vector<std::vector<double> > B3_sq_avg_realizations;
        
        // For computing TKE dissipation.
        
        std::vector<std::vector<double> > ddx_u_avg_realizations;
        std::vector<std::vector<double> > ddy_u_avg_realizations;
        std::vector<std::vector<double> > ddz_u_avg_realizations;
        
        std::vector<std::vector<double> > ddx_v_avg_realizations;
        std::vector<std::vector<double> > ddy_v_avg_realizations;
        std::vector<std::vector<double> > ddz_v_avg_realizations;
        
        std::vector<std::vector<double> > ddx_w_avg_realizations;
        std::vector<std::vector<double> > ddy_w_avg_realizations;
        std::vector<std::vector<double> > ddz_w_avg_realizations;
        
        std::vector<std::vector<double> > tau11_avg_realizations;
        std::vector<std::vector<double> > tau12_avg_realizations;
        std::vector<std::vector<double> > tau13_avg_realizations;
        std::vector<std::vector<double> > tau22_avg_realizations;
        std::vector<std::vector<double> > tau23_avg_realizations;
        std::vector<std::vector<double> > tau33_avg_realizations;
        
        std::vector<std::vector<double> > tau11_ddx_u_avg_realizations;
        std::vector<std::vector<double> > tau12_ddy_u_avg_realizations;
        std::vector<std::vector<double> > tau13_ddz_u_avg_realizations;
        
        std::vector<std::vector<double> > tau12_ddx_v_avg_realizations;
        std::vector<std::vector<double> > tau22_ddy_v_avg_realizations;
        std::vector<std::vector<double> > tau23_ddz_v_avg_realizations;
        
        std::vector<std::vector<double> > tau13_ddx_w_avg_realizations;
        std::vector<std::vector<double> > tau23_ddy_w_avg_realizations;
        std::vector<std::vector<double> > tau33_ddz_w_avg_realizations;
        
        // Whether the scratch arrays are filled.
        
        bool X_0_avg_computed;
        bool X_0_sq_avg_computed;
        bool X_0_X_1_avg_computed;
        bool Y_0_avg_computed;
        bool Y_0_sq_avg_computed;
        bool Y_0_Y_1_avg_computed;
        bool Z_0_avg_computed;
        bool Z_0_sq_avg_computed;
        bool Z_0_Z_1_avg_computed;
        bool rho_avg_computed;
        bool rho_sq_avg_computed;
        bool rho_inv_avg_computed;
        bool p_avg_computed;
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
        bool D_avg_computed;
        
        bool ddx_p_avg_computed;
        bool ddy_p_avg_computed;
        bool ddz_p_avg_computed;
        
        bool ddx_tau11_avg_computed;
        
        // For computing baroclinic torque.
        
        bool B1_sq_avg_computed;
        bool B2_sq_avg_computed;
        bool B3_sq_avg_computed;
        
        // For computing TKE dissipation.
        
        bool TKE_diss_computed;
        
        bool ddx_u_avg_computed;
        bool ddy_u_avg_computed;
        bool ddz_u_avg_computed;
        
        bool ddx_v_avg_computed;
        bool ddy_v_avg_computed;
        bool ddz_v_avg_computed;
        
        bool ddx_w_avg_computed;
        bool ddy_w_avg_computed;
        bool ddz_w_avg_computed;
        
        bool tau11_avg_computed;
        bool tau12_avg_computed;
        bool tau13_avg_computed;
        bool tau22_avg_computed;
        bool tau23_avg_computed;
        bool tau33_avg_computed;
        
        bool tau11_ddx_u_avg_computed;
        bool tau12_ddy_u_avg_computed;
        bool tau13_ddz_u_avg_computed;
        
        bool tau12_ddx_v_avg_computed;
        bool tau22_ddy_v_avg_computed;
        bool tau23_ddz_v_avg_computed;
        
        bool tau13_ddx_w_avg_computed;
        bool tau23_ddy_w_avg_computed;
        bool tau33_ddz_w_avg_computed;
        
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
                d_ensemble_statistics(ensemble_statistics),
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
                d_num_ghosts_derivative(3)
        {}
        
        /*
         * Compute averaged mass fraction with assumed homogeneity in y-direction (2D) or yz-plane (3D).
         */
        void
        computeAveragedMassFractionWithHomogeneityInYDirectionOrInYZPlane(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute averaged mole fraction with assumed homogeneity in y-direction (2D) or yz-plane (3D).
         */
        void
        computeAveragedMoleFractionWithHomogeneityInYDirectionOrInYZPlane(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute averaged volume fraction with assumed homogeneity in y-direction (2D) or yz-plane (3D).
         */
        void
        computeAveragedVolumeFractionWithHomogeneityInYDirectionOrInYZPlane(
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
         * Compute mole fraction variance with assumed homogeneity in y-direction (2D) or yz-plane (3D).
         */
        void
        computeMoleFractionVarianceWithHomogeneityInYDirectionOrInYZPlane(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute volume fraction variance with assumed homogeneity in y-direction (2D) or yz-plane (3D).
         */
        void
        computeVolumeFractionVarianceWithHomogeneityInYDirectionOrInYZPlane(
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
         * Compute mole fraction product with assumed homogeneity in y-direction (2D) or yz-plane (3D).
         */
        void
        computeMoleFractionProductWithHomogeneityInYDirectionOrInYZPlane(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute volume fraction product with assumed homogeneity in y-direction (2D) or yz-plane (3D).
         */
        void
        computeVolumeFractionProductWithHomogeneityInYDirectionOrInYZPlane(
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
         * Compute averaged pressure with assumed homogeneity in y-direction (2D) or yz-plane (3D).
         */
        void
        computeAveragedPressureWithHomogeneityInYDirectionOrInYZPlane(
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
         * Compute Reynolds normal stress x-component with assumed homogeneity in y-direction (2D) or yz-plane (3D).
         */
        void
        computeReynoldsNormalStressXWithHomogeneityInYDirectionOrInYZPlane(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute Reynolds normal stress y-component with assumed homogeneity in y-direction (2D) or yz-plane (3D).
         */
        void
        computeReynoldsNormalStressYWithHomogeneityInYDirectionOrInYZPlane(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute Reynolds normal stress z-component with assumed homogeneity in yz-plane (3D).
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
         * Compute enstrophy with assumed homogeneity in y-direction (2D) or yz-plane (3D).
         */
        void
        computeAveragedEnstrophyWithHomogeneityInYDirectionOrInYZPlane(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute scalar dissipation rate with assumed homogeneity in y-direction (2D) or yz-plane (3D).
         */
        void
        computeAveragedScalarDissipationRateWithHomogeneityInYDirectionOrInYZPlane(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute averaged dynamic shear viscosity with assumed homogeneity in y-direction (2D) or yz-plane (3D).
         */
        void
        computeAveragedShearViscosityWithHomogeneityInYDirectionOrInYZPlane(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute averaged mass diffusivity with assumed homogeneity in y-direction (2D) or yz-plane (3D).
         */
        void
        computeAveragedMassDiffusivityWithHomogeneityInYDirectionOrInYZPlane(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Computed averaged x-derivative of pressure with assumed homogeneity in y-direction (2D) or yz-plane (3D).
         */
        void
        computeAveragedXDerivativeOfPressureWithHomogeneityInYDirectionOrInYZPlane(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Computed averaged y-derivative of pressure with assumed homogeneity in y-direction (2D) or yz-plane (3D).
         */
        void
        computeAveragedYDerivativeOfPressureWithHomogeneityInYDirectionOrInYZPlane(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Computed averaged z-derivative of pressure with assumed homogeneity in yz-plane (3D).
         */
        void
        computeAveragedZDerivativeOfPressureWithHomogeneityInYZPlane(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute averaged x-derivative of normal shear stress in x-drection with assumed homogeneity in yz-plane (3D) to a file.
         */
        void
        computeAveragedXDerivativeNormalShearStressXWithHomogeneityInYDirectionOrInYZPlane(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Computed averaged squared baroclinic torque in z-direction with assumed homogeneity in y-direction (2D)
         * or yz-plane (3D).
         */
        void
        computeAveragedSquaredBaroclinicTorqueZWithHomogeneityInYDirectionOrInYZPlane(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Computed averaged squared baroclinic torque in x-direction with assumed homogeneity in yz-plane (3D).
         */
        void
        computeAveragedSquaredBaroclinicTorqueXWithHomogeneityInYZPlane(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Computed averaged squared baroclinic torque in y-direction with assumed homogeneity in yz-plane (3D).
         */
        void
        computeAveragedSquaredBaroclinicTorqueYWithHomogeneityInYZPlane(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Computed averaged quantities for TKE dissipation with assumed homogeneity in y-direction (2D) or yz-plane (3D) to
         * a file.
         */
        void
        computeAveragedQuantiitesForTKEDissipationWithHomogeneityInYDirectionOrInYZPlane(
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
            const double output_time) const;
        
        /*
         * Output spatial profile of ensemble averaged mole fraction with assumed homogeneity in y-direction (2D) or
         * yz-plane (3D) to a file.
         */
        void
        outputSpatialProfileEnsembleAveragedMoleFractionWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const double output_time) const;
        
        /*
         * Output spatial profile of ensemble averaged volume fraction with assumed homogeneity in y-direction (2D) or
         * yz-plane (3D) to a file.
         */
        void
        outputSpatialProfileEnsembleAveragedVolumeFractionWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const double output_time) const;
        
        /*
         * Output spatial profile of ensemble variance of mass fraction with assumed homogeneity in y-direction (2D) or
         * yz-plane (3D) to a file.
         */
        void
        outputSpatialProfileEnsembleMassFractionVarianceWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const double output_time) const;
        
        /*
         * Output spatial profile of ensemble variance of mole fraction with assumed homogeneity in y-direction (2D) or
         * yz-plane (3D) to a file.
         */
        void
        outputSpatialProfileEnsembleMoleFractionVarianceWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const double output_time) const;
        
        /*
         * Output spatial profile of ensemble variance of volume fraction with assumed homogeneity in y-direction (2D) or
         * yz-plane (3D) to a file.
         */
        void
        outputSpatialProfileEnsembleVolumeFractionVarianceWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const double output_time) const;
        
        /*
         * Output spatial profile of ensemble averaged density with assumed homogeneity in y-direction (2D) or
         * yz-plane (3D) to a file.
         */
        void
        outputSpatialProfileEnsembleAveragedDensityWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const double output_time) const;
        
        /*
         * Output spatial profile of ensemble variance of density with assumed homogeneity in y-direction (2D) or
         * yz-plane (3D) to a file.
         */
        void
        outputSpatialProfileEnsembleDensityVarianceWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const double output_time) const;
        
        /*
         * Output spatial profile of ensemble averaged specific volume with assumed homogeneity in y-direction (2D) or
         * yz-plane (3D) to a file.
         */
        void
        outputSpatialProfileEnsembleAveragedSpecificVolumeWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const double output_time) const;
        
        /*
         * Output spatial profile of ensemble averaged density-specific-volume covariance with assumed homogeneity in
         * y-direction (2D) or yz-plane (3D) to a file.
         */
        void
        outputSpatialProfileEnsembleAveragedDensitySpecificVolumeCovarianceWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const double output_time) const;
        
        /*
         * Output spatial profile of ensemble averaged pressure with assumed homogeneity in y-direction (2D) or
         * yz-plane (3D) to a file.
         */
        void
        outputSpatialProfileEnsembleAveragedPressureWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const double output_time) const;
        
        /*
         * Output spatial profile of ensemble averaged velocity x-component with assumed homogeneity in y-direction (2D)
         * or yz-plane (3D) to a file.
         */
        void
        outputSpatialProfileEnsembleAveragedVelocityXWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const double output_time) const;
        
        /*
         * Output spatial profile of ensemble averaged velocity y-component with assumed homogeneity in y-direction (2D)
         * or yz-plane (3D) to a file.
         */
        void
        outputSpatialProfileEnsembleAveragedVelocityYWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const double output_time) const;
        
        /*
         * Output spatial profile of ensemble averaged velocity z-component with assumed homogeneity in yz-plane (3D)
         * to a file.
         */
        void
        outputSpatialProfileEnsembleAveragedVelocityZWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const double output_time) const;
        
        /*
         * Output spatial profile of ensemble averaged momentum x-component with assumed homogeneity in y-direction (2D)
         * or yz-plane (3D) to a file.
         */
        void
        outputSpatialProfileEnsembleAveragedMomentumXWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const double output_time) const;
        
        /*
         * Output spatial profile of ensemble averaged momentum y-component with assumed homogeneity in y-direction (2D)
         * or yz-plane (3D) to a file.
         */
        void
        outputSpatialProfileEnsembleAveragedMomentumYWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const double output_time) const;
        
        /*
         * Output spatial profile of ensemble averaged momentum z-component with assumed homogeneity in yz-plane (3D)
         * to a file.
         */
        void
        outputSpatialProfileEnsembleAveragedMomentumZWithHomogeneitynYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const double output_time) const;
        
        /*
         * Output spatial profile of ensemble turbulent mass flux in x-direction with assumed homogeneity in
         * y-direction (2D) or yz-plane (3D) to a file.
         */
        void
        outputSpatialProfileEnsembleTurbulentMassFluxXWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const double output_time) const;
        
        /*
         * Output spatial profile of ensemble Reynolds normal stress x-component with assumed homogeneity in
         * y-direction (2D) or yz-plane (3D) to a file.
         */
        void
        outputSpatialProfileEnsembleReynoldsNormalStressXWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const double output_time) const;
        
        /*
         * Output spatial profile of ensemble Reynolds normal stress y-component with assumed homogeneity in
         * y-direction (2D) or yz-plane (3D) to a file.
         */
        void
        outputSpatialProfileEnsembleReynoldsNormalStressYWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const double output_time) const;
        
        /*
         * Output spatial profile of ensemble Reynolds normal stress z-component with assumed homogeneity in
         * yz-plane (3D) to a file.
         */
        void
        outputSpatialProfileEnsembleReynoldsNormalStressZWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const double output_time) const;
        
        /*
         * Output spatial profile of ensemble averaged enstrophy with assumed homogeneity in y-direction (2D) or
         * yz-plane (3D) to a file.
         */
        void
        outputSpatialProfileEnsembleAveragedEnstrophyWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const double output_time) const;
        
        /*
         * Output spatial profile of ensemble averaged scalar dissipation rate with assumed homogeneity in
         * y-direction (2D) or yz-plane (3D) to a file.
         */
        void
        outputSpatialProfileEnsembleAveragedScalarDissipationRateWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const double output_time) const;
        
        /*
         * Output spatial profile of ensemble RMS of baroclinic torque in z-direction with assumed homogeneity in
         * y-direction (2D) or yz-plane (3D) to a file.
         */
        void
        outputSpatialProfileEnsembleRMSBaroclinicTorqueZWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const double output_time) const;
        
        /*
         * Output spatial profile of ensemble RMS of baroclinic torque in x-direction with assumed homogeneity in
         * yz-plane (3D) to a file.
         */
        void
        outputSpatialProfileEnsembleRMSBaroclinicTorqueXWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const double output_time) const;
        
        /*
         * Output spatial profile of ensemble RMS of baroclinic torque in y-direction with assumed homogeneity in
         * yz-plane (3D) to a file.
         */
        void
        outputSpatialProfileEnsembleRMSBaroclinicTorqueYWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const double output_time) const;
        
        /// Non-spatial profiles.
        
        /*
         * Output ensemble mixing width in x-direction to a file.
         */
        void
        outputEnsembleMixingWidthInXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const;
        
        /*
         * Output ensemble mixing width in x-direction using mole fractions to a file.
         */
        void
        outputEnsembleMixingWidthInXDirectionWithMoleFractions(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const;
        
        /*
         * Output ensemble mixing width in x-direction using volume fractions to a file.
         */
        void
        outputEnsembleMixingWidthInXDirectionWithVolumeFractions(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const;
        
        /*
         * Output ensemble minimum interface location in x-direction to a file.
         */
        void
        outputEnsembleInterfaceMinInXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const;
        
        /*
         * Output ensemble minimum interface location in x-direction using mole fractions to a file.
         */
        void
        outputEnsembleInterfaceMinInXDirectionWithMoleFractions(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const;
        
        /*
         * Output ensemble minimum interface location in x-direction using volume fractions to a file.
         */
        void
        outputEnsembleInterfaceMinInXDirectionWithVolumeFractions(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const;
        
        /*
         * Output ensemble maximum interface location in x-direction to a file.
         */
        void
        outputEnsembleInterfaceMaxInXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const;
        
        /*
         * Output ensemble maximum interface location in x-direction using mole fractions to a file.
         */
        void
        outputEnsembleInterfaceMaxInXDirectionWithMoleFractions(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const;
        
        /*
         * Output ensemble maximum interface location in x-direction using volume fractions to a file.
         */
        void
        outputEnsembleInterfaceMaxInXDirectionWithVolumeFractions(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const;
        
        /*
         * Output ensemble mixedness in x-direction to a file.
         */
        void
        outputEnsembleMixednessInXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const;
        
        /*
         * Output ensemble mixedness in x-direction using mole fractions to a file.
         */
        void
        outputEnsembleMixednessInXDirectionWithMoleFractions(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const;
        
        /*
         * Output ensemble mixedness in x-direction using volume fractions to a file.
         */
        void
        outputEnsembleMixednessInXDirectionWithVolumeFractions(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const;
        
        /*
         * Output mean of ensemble density mean in mixing layer with with assumed homogeneity in y-direction (2D) or
         * yz-plane (3D) to a file.
         */
        void
        outputEnsembleDensityMeanWithInhomogeneousXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const;
        
        /*
         * Output mean of ensemble specific volume mean in mixing layer with with assumed homogeneity in
         * y-direction (2D) or yz-plane (3D) to a file.
         */
        void
        outputEnsembleSpecificVolumeMeanWithInhomogeneousXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const;
        
        /*
         * Output mean of ensemble pressure mean in mixing layer with with assumed homogeneity in y-direction (2D) or
         * yz-plane (3D) to a file.
         */
        void
        outputEnsemblePressureMeanWithInhomogeneousXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const;
        
        /*
         * Output mean of ensemble velocity associated with turbulent mass flux component in x-direction in
         * mixing layer with assumed homogeneity in y-direction (2D) or yz-plane (3D) to a file.
         */
        void
        outputEnsembleTurbulentMassFluxVelocityXMeanWithInhomogeneousXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const;
        
        /*
         * Output mean of ensemble density-specific-volume covariance in mixing layer with assumed homogeneity in
         * y-direction (2D) or yz-plane (3D) to a file.
         */
        void
        outputEnsembleDensitySpecificVolumeCovarianceMeanWithInhomogeneousXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const;
        
        /*
         * Output ensemble TKE integrated with assumed homogeneity in y-direction (2D) or yz-plane (3D) to a file.
         */
        void
        outputEnsembleTKEIntegrateWithInhomogeneousXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const;
        
        /*
         * Output mean of ensemble Reynolds normal stress component in x-direction in mixing layer with assumed
         * homogeneity in y-direction (2D) or yz-plane (3D) to a file.
         */
        void
        outputEnsembleReynoldsNormalStressXMeanWithInhomogeneousXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const;
        
        /*
         * Output mean of ensemble Reynolds normal stress component in y-direction in mixing layer with assumed
         * homogeneity in y-direction (2D) or yz-plane (3D) to a file.
         */
        void
        outputEnsembleReynoldsNormalStressYMeanWithInhomogeneousXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const;
        
        /*
         * Output mean of ensemble Reynolds normal stress component in z-direction in mixing layer with assumed
         * homogeneity in yz-plane (3D) to a file.
         */
        void
        outputEnsembleReynoldsNormalStressZMeanWithInhomogeneousXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const;
        
        /*
         * Output ensemble enstrophy integrated to a file.
         */
        void
        outputEnsembleEnstrophyIntegrated(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const;
        
        /*
         * Output mean of ensemble enstrophy in mixing layer with assumed homogeneity in y-direction (2D) or yz-plane (3D)
         * to a file.
         */
        void
        outputEnsembleEnstrophyMeanWithInhomogeneousXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const;
        
        /*
         * Output ensemble scalar dissipation rate of first species integrated to a file.
         */
        void
        outputEnsembleScalarDissipationRateIntegrated(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const;
        
        /*
         * Output mean of ensemble scalar dissipation rate of first species in mixing layer with assumed homogeneity in
         * y-direction (2D) or yz-plane (3D) to a file.
         */
        void
        outputEnsembleScalarDissipationRateMeanWithInhomogeneousXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const;
        
        /*
         * Output turbulent Reynolds number based on mixing width with assumed homogeneity in y-direction (2D)
         * or yz-plane (3D) to a file.
         */
        void
        outputEnsembleMixingWidthReynoldsNumberWithInhomogeneousXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const;
        
        /*
         * Output turbulent Mach number with assumed homogeneity in y-direction (2D) or yz-plane (3D) to a file.
         */
        void
        outputEnsembleTurbulentMachNumberWithInhomogeneousXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const;
        
        /*
         * Output effective Atwood number with assumed homogeneity in y-direction (2D) or yz-plane (3D) to a file.
         */
        void
        outputEnsembleEffectiveAtwoodNumberWithInhomogeneousXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const;
        
        /*
         * Output ensemble RMS of baroclinic torque in z-direction integrated to a file.
         */
        void
        outputEnsembleRMSBaroclinicTorqueZIntegrated(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const;
        
        /*
         * Output mean of ensemble RMS of baroclinic torque in z-direction in mixing layer with assumed homogeneity in
         * y-direction (2D) or yz-plane (3D) to a file.
         */
        void
        outputEnsembleRMSBaroclinicTorqueZMeanWithInhomogeneousXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const;
        
        /*
         * Output ensemble RMS of baroclinic torque in x-direction integrated to a file.
         */
        void
        outputEnsembleRMSBaroclinicTorqueXIntegrated(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const;
        
        /*
         * Output mean of ensemble RMS of baroclinic torque in x-direction in mixing layer with assumed homogeneity in
         * yz-plane (3D) to a file.
         */
        void
        outputEnsembleRMSBaroclinicTorqueXMeanWithInhomogeneousXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const;
        
        /*
         * Output ensemble RMS of baroclinic torque in y-direction integrated to a file.
         */
        void
        outputEnsembleRMSBaroclinicTorqueYIntegrated(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const;
        
        /*
         * Output mean of ensemble RMS of baroclinic torque in y-direction in mixing layer with assumed homogeneity in
         * yz-plane (3D) to a file.
         */
        void
        outputEnsembleRMSBaroclinicTorqueYMeanWithInhomogeneousXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const;
        
        /*
         * Output ensemble TKE production term 1 integrated to a file.
         */
        void
        outputEnsembleTKEProduction1Integrated(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const;
        
        /*
         * Output mean of ensemble TKE production term 1 in mixing layer with assumed homogeneity in y-direction (2D) or
         * yz-plane (3D) to a file.
         */
        void
        outputEnsembleTKEProduction1MeanWithInhomogeneousXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const;
        
        /*
         * Output ensemble TKE production term 2 integrated to a file.
         */
        void
        outputEnsembleTKEProduction2Integrated(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const;
        
        /*
         * Output mean of ensemble TKE production term 2 in mixing layer with assumed homogeneity in y-direction (2D) or
         * yz-plane (3D) to a file.
         */
        void
        outputEnsembleTKEProduction2MeanWithInhomogeneousXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const;
        
        /*
         * Output ensemble TKE production term 3 integrated to a file.
         */
        void
        outputEnsembleTKEProduction3Integrated(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const;
        
        /*
         * Output mean of ensemble TKE production term 3 in mixing layer with assumed homogeneity in y-direction (2D) or
         * yz-plane (3D) to a file.
         */
        void
        outputEnsembleTKEProduction3MeanWithInhomogeneousXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const;
        
        /*
         * Output ensemble TKE dissipation integrated to a file.
         */
        void
        outputEnsembleTKEDissipationIntegrated(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const;
        
        /*
         * Output mean of ensemble TKE dissipation in mixing layer with assumed homogeneity in y-direction (2D) or
         * yz-plane (3D) to a file.
         */
        void
        outputEnsembleTKEDissipationMeanWithInhomogeneousXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const;
        
        /*
         * Output mean of ensemble dynamic shear viscosity in mixing layer with assumed homogeneity in y-direction (2D)
         * or yz-plane (3D) to a file.
         */
        void
        outputEnsembleShearViscosityMeanWithInhomogeneousXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const;
        
        /*
         * Output mean of ensemble mass diffusivity in mixing layer with assumed homogeneity in y-direction (2D) or
         * yz-plane (3D) to a file.
         */
        void
        outputEnsembleMassDiffusivityMeanWithInhomogeneousXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const;
        
        /*
         * Compute averaged shear stress component with only x direction as inhomogeneous direction.
         * component_idx:
         * 0: tau11
         * 1: tau12
         * 2: tau13
         * 3: tau22
         * 4: tau23
         * 5: tau33
         */
        std::vector<double>
        getAveragedShearStressComponentWithInhomogeneousXDirection(
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Compute averaged derivative of shear stress component with only x direction as inhomogeneous direction.
         * component_idx:
         * 0: tau11
         * 1: tau12
         * 2: tau13
         * 3: tau22
         * 4: tau23
         * 5: tau33
         */
        std::vector<double>
        getAveragedDerivativeOfShearStressComponentWithInhomogeneousXDirection(
            const int component_idx,
            const int derivative_direction,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Compute averaged value (on product of variable derivatives and shear stress component) with only x direction
         * as inhomogeneous direction.
         * component_idx:
         * 0: tau11
         * 1: tau12
         * 2: tau13
         * 3: tau22
         * 4: tau23
         * 5: tau33
         */
        std::vector<double>
        getAveragedQuantityWithShearStressComponentWithInhomogeneousXDirection(
            const std::vector<std::string>& quantity_names,
            const std::vector<int>& component_indices,
            const std::vector<bool>& use_derivative,
            const std::vector<int>& derivative_directions,
            const std::vector<bool>& use_reciprocal,
            const int shear_stress_component_idx,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Get gravity vector from the flow model database.
         */
        std::vector<double> getGravityVector() const;
        
        /*
         * Compute the one-dimensional derivative given a vector.
         */
        std::vector<double> computeDerivativeOfVector1D(
            const std::vector<double> quantity_vector,
            const double dx) const;
        
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
 * Compute averaged mole fraction with assumed homogeneity in y-direction (2D) or yz-plane (3D).
 */
void
RTIRMIStatisticsUtilities::computeAveragedMoleFractionWithHomogeneityInYDirectionOrInYZPlane(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'MOLE_FRACTION_AVG_SP' can be computed with two species only."
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
    
    std::vector<double> X_0_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        "MOLE_FRACTIONS",
        0,
        data_context);
    
    std::vector<std::vector<double> >& X_0_avg_realizations = d_ensemble_statistics->X_0_avg_realizations;
    X_0_avg_realizations.push_back(X_0_avg_global);
    
    d_ensemble_statistics->X_0_avg_computed = true;
}


/*
 * Compute averaged volume fraction with assumed homogeneity in y-direction (2D) or yz-plane (3D).
 */
void
RTIRMIStatisticsUtilities::computeAveragedVolumeFractionWithHomogeneityInYDirectionOrInYZPlane(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'VOLUME_FRACTION_AVG_SP' can be computed with two species only."
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
    
    std::vector<double> Z_0_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        "VOLUME_FRACTIONS",
        0,
        data_context);
    
    std::vector<std::vector<double> >& Z_0_avg_realizations = d_ensemble_statistics->Z_0_avg_realizations;
    Z_0_avg_realizations.push_back(Z_0_avg_global);
    
    d_ensemble_statistics->Z_0_avg_computed = true;
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
 * Compute mole fraction variance with assumed homogeneity in y-direction (2D) or yz-plane (3D).
 */
void
RTIRMIStatisticsUtilities::computeMoleFractionVarianceWithHomogeneityInYDirectionOrInYZPlane(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'MOLE_FRACTION_VAR_SP' can be computed with two species only."
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
    
    quantity_names.push_back("MOLE_FRACTIONS");
    component_indices.push_back(0);
    
    quantity_names.push_back("MOLE_FRACTIONS");
    component_indices.push_back(0);
    
    std::vector<double> X_0_sq_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    
    std::vector<std::vector<double> >& X_0_sq_avg_realizations = d_ensemble_statistics->X_0_sq_avg_realizations;
    X_0_sq_avg_realizations.push_back(X_0_sq_avg_global);
    
    d_ensemble_statistics->X_0_sq_avg_computed = true;
}


/*
 * Compute volume fraction variance with assumed homogeneity in y-direction (2D) or yz-plane (3D).
 */
void
RTIRMIStatisticsUtilities::computeVolumeFractionVarianceWithHomogeneityInYDirectionOrInYZPlane(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'VOLUME_FRACTION_VAR_SP' can be computed with two species only."
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
    
    quantity_names.push_back("VOLUME_FRACTIONS");
    component_indices.push_back(0);
    
    quantity_names.push_back("VOLUME_FRACTIONS");
    component_indices.push_back(0);
    
    std::vector<double> Z_0_sq_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    
    std::vector<std::vector<double> >& Z_0_sq_avg_realizations = d_ensemble_statistics->Z_0_sq_avg_realizations;
    Z_0_sq_avg_realizations.push_back(Z_0_sq_avg_global);
    
    d_ensemble_statistics->Z_0_sq_avg_computed = true;
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
 * Compute mole fraction product with assumed homogeneity in y-direction (2D) or yz-plane (3D).
 */
void
RTIRMIStatisticsUtilities::computeMoleFractionProductWithHomogeneityInYDirectionOrInYZPlane(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'MOLE_FRACTION_PROD_SP' can be computed with two species only."
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
    
    quantity_names.push_back("MOLE_FRACTIONS");
    component_indices.push_back(0);
    
    quantity_names.push_back("MOLE_FRACTIONS");
    component_indices.push_back(1);
    
    std::vector<double> X_0_X_1_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    
    std::vector<std::vector<double> >& X_0_X_1_avg_realizations = d_ensemble_statistics->X_0_X_1_avg_realizations;
    X_0_X_1_avg_realizations.push_back(X_0_X_1_avg_global);
    
    d_ensemble_statistics->X_0_X_1_avg_computed = true;
}


/*
 * Compute volume fraction product with assumed homogeneity in y-direction (2D) or yz-plane (3D).
 */
void
RTIRMIStatisticsUtilities::computeVolumeFractionProductWithHomogeneityInYDirectionOrInYZPlane(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'VOLUME_FRACTION_PROD_SP' can be computed with two species only."
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
    
    quantity_names.push_back("VOLUME_FRACTIONS");
    component_indices.push_back(0);
    
    quantity_names.push_back("VOLUME_FRACTIONS");
    component_indices.push_back(1);
    
    std::vector<double> Z_0_Z_1_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    
    std::vector<std::vector<double> >& Z_0_Z_1_avg_realizations = d_ensemble_statistics->Z_0_Z_1_avg_realizations;
    Z_0_Z_1_avg_realizations.push_back(Z_0_Z_1_avg_global);
    
    d_ensemble_statistics->Z_0_Z_1_avg_computed = true;
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
 * Compute averaged pressure with assumed homogeneity in y-direction (2D) or yz-plane (3D).
 */
void
RTIRMIStatisticsUtilities::computeAveragedPressureWithHomogeneityInYDirectionOrInYZPlane(
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
    
    std::vector<double> p_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        "PRESSURE",
        0,
        data_context);
    
    std::vector<std::vector<double> >& p_avg_realizations = d_ensemble_statistics->p_avg_realizations;
    p_avg_realizations.push_back(p_avg_global);
    
    d_ensemble_statistics->p_avg_computed = true;
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
 * Compute Reynolds normal stress x-component with assumed homogeneity in y-direction (2D) or yz-plane (3D).
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
 * Compute Reynolds normal stress y-component with assumed homogeneity in y-direction (2D) or yz-plane (3D).
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
 * Compute Reynolds normal stress z-component with assumed homogeneity in yz-plane (3D).
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
 * Compute enstrophy with assumed homogeneity in y-direction (2D) or yz-plane (3D).
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
 * Compute scalar dissipation rate with assumed homogeneity in y-direction (2D) or yz-plane (3D).
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
 * Compute averaged dynamic shear viscosity with assumed homogeneity in y-direction (2D) or yz-plane (3D).
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
    
    std::vector<std::vector<double> >& mu_avg_realizations = d_ensemble_statistics->mu_avg_realizations;
    
    std::vector<double> shear_viscosity = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        "SHEAR_VISCOSITY",
        0,
        data_context);
    
    mu_avg_realizations.push_back(shear_viscosity);
    
    d_ensemble_statistics->mu_avg_computed = true;
}


/*
 * Compute averaged mass diffusivity with assumed homogeneity in y-direction (2D) or yz-plane (3D).
 */
void
RTIRMIStatisticsUtilities::computeAveragedMassDiffusivityWithHomogeneityInYDirectionOrInYZPlane(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'MASS_DIFFUSIVITY_AVG_SP' can be computed with two species only."
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
        flow_model_tmp,
        true);
    
    std::vector<std::vector<double> >& D_avg_realizations = d_ensemble_statistics->D_avg_realizations;
    
    std::vector<double> mass_diffusivity = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        "MASS_DIFFUSIVITIES",
        0,
        data_context);
    
    D_avg_realizations.push_back(mass_diffusivity);
    
    d_ensemble_statistics->D_avg_computed = true;
}


/*
 * Computed averaged x-derivative of pressure with assumed homogeneity in y-direction (2D) or yz-plane (3D).
 */
void
RTIRMIStatisticsUtilities::computeAveragedXDerivativeOfPressureWithHomogeneityInYDirectionOrInYZPlane(
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
    
    std::vector<std::vector<double> >& ddx_p_avg_realizations = d_ensemble_statistics->ddx_p_avg_realizations;
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    std::vector<bool> use_derivative;
    std::vector<int> derivative_directions;
    
    quantity_names.push_back("PRESSURE");
    component_indices.push_back(0);
    use_derivative.push_back(true);
    derivative_directions.push_back(0);
    
    std::vector<double> ddx_p_avg = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
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
    
    ddx_p_avg_realizations.push_back(ddx_p_avg);
    
    d_ensemble_statistics->ddx_p_avg_computed = true;
}


/*
 * Compute averaged x-derivative of normal shear stress in x-drection with assumed homogeneity in yz-plane (3D).
 */
void
RTIRMIStatisticsUtilities::computeAveragedXDerivativeNormalShearStressXWithHomogeneityInYDirectionOrInYZPlane(
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
    
    std::vector<std::vector<double> >& ddx_tau11_avg_realizations = d_ensemble_statistics->ddx_tau11_avg_realizations;
    
    std::vector<double> ddx_tau11_avg = getAveragedDerivativeOfShearStressComponentWithInhomogeneousXDirection(
        0,
        0,
        patch_hierarchy,
        data_context);
    
    ddx_tau11_avg_realizations.push_back(ddx_tau11_avg);
    
    d_ensemble_statistics->ddx_tau11_avg_computed = true;
}


/*
 * Computed averaged y-derivative of pressure with assumed homogeneity in y-direction (2D) or yz-plane (3D).
 */
void
RTIRMIStatisticsUtilities::computeAveragedYDerivativeOfPressureWithHomogeneityInYDirectionOrInYZPlane(
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
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no 'DDY_PRESSURE_AVG_SP' for one-dimensional problem."
            << std::endl);
    }
    
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
    FlowModelMPIHelperAverage MPI_helper_average = FlowModelMPIHelperAverage(
        "MPI_helper_average",
        d_dim,
        d_grid_geometry,
        patch_hierarchy,
        flow_model_tmp);
    
    std::vector<std::vector<double> >& ddy_p_avg_realizations = d_ensemble_statistics->ddy_p_avg_realizations;
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    std::vector<bool> use_derivative;
    std::vector<int> derivative_directions;
    
    quantity_names.push_back("PRESSURE");
    component_indices.push_back(0);
    use_derivative.push_back(true);
    derivative_directions.push_back(1);
    
    std::vector<double> ddy_p_avg = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
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
    
    ddy_p_avg_realizations.push_back(ddy_p_avg);
    
    d_ensemble_statistics->ddy_p_avg_computed = true;
}


/*
 * Computed averaged z-derivative of pressure with assumed homogeneity in yz-plane (3D).
 */
void
RTIRMIStatisticsUtilities::computeAveragedZDerivativeOfPressureWithHomogeneityInYZPlane(
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
    
    if (d_dim == tbox::Dimension(1) || d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no 'DDZ_PRESSURE_AVG_SP' for one-dimensional or two-dimensional problem."
            << std::endl);
    }
    
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
    FlowModelMPIHelperAverage MPI_helper_average = FlowModelMPIHelperAverage(
        "MPI_helper_average",
        d_dim,
        d_grid_geometry,
        patch_hierarchy,
        flow_model_tmp);
    
    std::vector<std::vector<double> >& ddz_p_avg_realizations = d_ensemble_statistics->ddz_p_avg_realizations;
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    std::vector<bool> use_derivative;
    std::vector<int> derivative_directions;
    
    quantity_names.push_back("PRESSURE");
    component_indices.push_back(0);
    use_derivative.push_back(true);
    derivative_directions.push_back(2);
    
    std::vector<double> ddz_p_avg = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
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
    
    ddz_p_avg_realizations.push_back(ddz_p_avg);
    
    d_ensemble_statistics->ddz_p_avg_computed = true;
}


/*
 * Computed averaged squared baroclinic torque in z-direction with assumed homogeneity in y-direction (2D)
 * or yz-plane (3D).
 */
void
RTIRMIStatisticsUtilities::computeAveragedSquaredBaroclinicTorqueZWithHomogeneityInYDirectionOrInYZPlane(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no 'BARO_TQ_Z_RMS_SP' for one-dimensional problem."
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
    
    const hier::IntVector& finest_level_dims = MPI_helper_average.getFinestRefinedDomainNumberOfPoints();
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    std::vector<bool> use_derivative;
    std::vector<int> derivative_directions;
    std::vector<bool> use_reciprocal;
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(false);
    derivative_directions.push_back(-1);
    use_reciprocal.push_back(true);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(false);
    derivative_directions.push_back(-1);
    use_reciprocal.push_back(true);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(false);
    derivative_directions.push_back(-1);
    use_reciprocal.push_back(true);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(false);
    derivative_directions.push_back(-1);
    use_reciprocal.push_back(true);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(true);
    derivative_directions.push_back(0);
    use_reciprocal.push_back(false);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(true);
    derivative_directions.push_back(0);
    use_reciprocal.push_back(false);
    
    quantity_names.push_back("PRESSURE");
    component_indices.push_back(0);
    use_derivative.push_back(true);
    derivative_directions.push_back(1);
    use_reciprocal.push_back(false);
    
    quantity_names.push_back("PRESSURE");
    component_indices.push_back(0);
    use_derivative.push_back(true);
    derivative_directions.push_back(1);
    use_reciprocal.push_back(false);
    
    std::vector<double> B3_part_1 = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        use_derivative,
        derivative_directions,
        use_reciprocal,
        d_num_ghosts_derivative,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    use_derivative.clear();
    derivative_directions.clear();
    use_reciprocal.clear();
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(false);
    derivative_directions.push_back(-1);
    use_reciprocal.push_back(true);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(false);
    derivative_directions.push_back(-1);
    use_reciprocal.push_back(true);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(false);
    derivative_directions.push_back(-1);
    use_reciprocal.push_back(true);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(false);
    derivative_directions.push_back(-1);
    use_reciprocal.push_back(true);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(true);
    derivative_directions.push_back(1);
    use_reciprocal.push_back(false);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(true);
    derivative_directions.push_back(1);
    use_reciprocal.push_back(false);
    
    quantity_names.push_back("PRESSURE");
    component_indices.push_back(0);
    use_derivative.push_back(true);
    derivative_directions.push_back(0);
    use_reciprocal.push_back(false);
    
    quantity_names.push_back("PRESSURE");
    component_indices.push_back(0);
    use_derivative.push_back(true);
    derivative_directions.push_back(0);
    use_reciprocal.push_back(false);
    
    std::vector<double> B3_part_2 = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        use_derivative,
        derivative_directions,
        use_reciprocal,
        d_num_ghosts_derivative,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    use_derivative.clear();
    derivative_directions.clear();
    use_reciprocal.clear();
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(false);
    derivative_directions.push_back(-1);
    use_reciprocal.push_back(true);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(false);
    derivative_directions.push_back(-1);
    use_reciprocal.push_back(true);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(false);
    derivative_directions.push_back(-1);
    use_reciprocal.push_back(true);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(false);
    derivative_directions.push_back(-1);
    use_reciprocal.push_back(true);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(true);
    derivative_directions.push_back(0);
    use_reciprocal.push_back(false);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(true);
    derivative_directions.push_back(1);
    use_reciprocal.push_back(false);
    
    quantity_names.push_back("PRESSURE");
    component_indices.push_back(0);
    use_derivative.push_back(true);
    derivative_directions.push_back(0);
    use_reciprocal.push_back(false);
    
    quantity_names.push_back("PRESSURE");
    component_indices.push_back(0);
    use_derivative.push_back(true);
    derivative_directions.push_back(1);
    use_reciprocal.push_back(false);
    
    std::vector<double> B3_part_3 = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        use_derivative,
        derivative_directions,
        use_reciprocal,
        d_num_ghosts_derivative,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    use_derivative.clear();
    derivative_directions.clear();
    use_reciprocal.clear();
    
    std::vector<double> B3_sq_avg(finest_level_dims[0], double(0));
    for (int i = 0; i < finest_level_dims[0]; i++)
    {
        B3_sq_avg[i] = B3_part_1[i] + B3_part_2[i] - double(2)*B3_part_3[i];
    }
    
    std::vector<std::vector<double> >& B3_sq_avg_realizations = d_ensemble_statistics->B3_sq_avg_realizations;
    B3_sq_avg_realizations.push_back(B3_sq_avg);
    
    d_ensemble_statistics->B3_sq_avg_computed = true;
}


/*
 * Computed averaged squared baroclinic torque in x-direction with assumed homogeneity in yz-plane (3D).
 */
void
RTIRMIStatisticsUtilities::computeAveragedSquaredBaroclinicTorqueXWithHomogeneityInYZPlane(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    if (d_dim == tbox::Dimension(1) || d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no 'BARO_TQ_X_RMS_SP' for one-dimensional or two-dimensional problem."
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
    
    const hier::IntVector& finest_level_dims = MPI_helper_average.getFinestRefinedDomainNumberOfPoints();
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    std::vector<bool> use_derivative;
    std::vector<int> derivative_directions;
    std::vector<bool> use_reciprocal;
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(false);
    derivative_directions.push_back(-1);
    use_reciprocal.push_back(true);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(false);
    derivative_directions.push_back(-1);
    use_reciprocal.push_back(true);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(false);
    derivative_directions.push_back(-1);
    use_reciprocal.push_back(true);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(false);
    derivative_directions.push_back(-1);
    use_reciprocal.push_back(true);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(true);
    derivative_directions.push_back(1);
    use_reciprocal.push_back(false);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(true);
    derivative_directions.push_back(1);
    use_reciprocal.push_back(false);
    
    quantity_names.push_back("PRESSURE");
    component_indices.push_back(0);
    use_derivative.push_back(true);
    derivative_directions.push_back(2);
    use_reciprocal.push_back(false);
    
    quantity_names.push_back("PRESSURE");
    component_indices.push_back(0);
    use_derivative.push_back(true);
    derivative_directions.push_back(2);
    use_reciprocal.push_back(false);
    
    std::vector<double> B1_part_1 = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        use_derivative,
        derivative_directions,
        use_reciprocal,
        d_num_ghosts_derivative,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    use_derivative.clear();
    derivative_directions.clear();
    use_reciprocal.clear();
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(false);
    derivative_directions.push_back(-1);
    use_reciprocal.push_back(true);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(false);
    derivative_directions.push_back(-1);
    use_reciprocal.push_back(true);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(false);
    derivative_directions.push_back(-1);
    use_reciprocal.push_back(true);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(false);
    derivative_directions.push_back(-1);
    use_reciprocal.push_back(true);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(true);
    derivative_directions.push_back(2);
    use_reciprocal.push_back(false);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(true);
    derivative_directions.push_back(2);
    use_reciprocal.push_back(false);
    
    quantity_names.push_back("PRESSURE");
    component_indices.push_back(0);
    use_derivative.push_back(true);
    derivative_directions.push_back(1);
    use_reciprocal.push_back(false);
    
    quantity_names.push_back("PRESSURE");
    component_indices.push_back(0);
    use_derivative.push_back(true);
    derivative_directions.push_back(1);
    use_reciprocal.push_back(false);
    
    std::vector<double> B1_part_2 = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        use_derivative,
        derivative_directions,
        use_reciprocal,
        d_num_ghosts_derivative,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    use_derivative.clear();
    derivative_directions.clear();
    use_reciprocal.clear();
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(false);
    derivative_directions.push_back(-1);
    use_reciprocal.push_back(true);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(false);
    derivative_directions.push_back(-1);
    use_reciprocal.push_back(true);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(false);
    derivative_directions.push_back(-1);
    use_reciprocal.push_back(true);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(false);
    derivative_directions.push_back(-1);
    use_reciprocal.push_back(true);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(true);
    derivative_directions.push_back(1);
    use_reciprocal.push_back(false);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(true);
    derivative_directions.push_back(2);
    use_reciprocal.push_back(false);
    
    quantity_names.push_back("PRESSURE");
    component_indices.push_back(0);
    use_derivative.push_back(true);
    derivative_directions.push_back(1);
    use_reciprocal.push_back(false);
    
    quantity_names.push_back("PRESSURE");
    component_indices.push_back(0);
    use_derivative.push_back(true);
    derivative_directions.push_back(2);
    use_reciprocal.push_back(false);
    
    std::vector<double> B1_part_3 = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        use_derivative,
        derivative_directions,
        use_reciprocal,
        d_num_ghosts_derivative,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    use_derivative.clear();
    derivative_directions.clear();
    use_reciprocal.clear();
    
    
    std::vector<double> B1_sq_avg(finest_level_dims[0], double(0));
    for (int i = 0; i < finest_level_dims[0]; i++)
    {
        B1_sq_avg[i] = B1_part_1[i] + B1_part_2[i] - double(2)*B1_part_3[i];
    }
    
    std::vector<std::vector<double> >& B1_sq_avg_realizations = d_ensemble_statistics->B1_sq_avg_realizations;
    B1_sq_avg_realizations.push_back(B1_sq_avg);
    
    d_ensemble_statistics->B1_sq_avg_computed = true;
}


/*
 * Computed averaged squared baroclinic torque in y-direction with assumed homogeneity in yz-plane (3D).
 */
void
RTIRMIStatisticsUtilities::computeAveragedSquaredBaroclinicTorqueYWithHomogeneityInYZPlane(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    if (d_dim == tbox::Dimension(1) || d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no 'BARO_TQ_Y_RMS_SP' for one-dimensional or two-dimensional problem."
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
    
    const hier::IntVector& finest_level_dims = MPI_helper_average.getFinestRefinedDomainNumberOfPoints();
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    std::vector<bool> use_derivative;
    std::vector<int> derivative_directions;
    std::vector<bool> use_reciprocal;
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(false);
    derivative_directions.push_back(-1);
    use_reciprocal.push_back(true);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(false);
    derivative_directions.push_back(-1);
    use_reciprocal.push_back(true);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(false);
    derivative_directions.push_back(-1);
    use_reciprocal.push_back(true);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(false);
    derivative_directions.push_back(-1);
    use_reciprocal.push_back(true);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(true);
    derivative_directions.push_back(2);
    use_reciprocal.push_back(false);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(true);
    derivative_directions.push_back(2);
    use_reciprocal.push_back(false);
    
    quantity_names.push_back("PRESSURE");
    component_indices.push_back(0);
    use_derivative.push_back(true);
    derivative_directions.push_back(0);
    use_reciprocal.push_back(false);
    
    quantity_names.push_back("PRESSURE");
    component_indices.push_back(0);
    use_derivative.push_back(true);
    derivative_directions.push_back(0);
    use_reciprocal.push_back(false);
    
    std::vector<double> B2_part_1 = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        use_derivative,
        derivative_directions,
        use_reciprocal,
        d_num_ghosts_derivative,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    use_derivative.clear();
    derivative_directions.clear();
    use_reciprocal.clear();
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(false);
    derivative_directions.push_back(-1);
    use_reciprocal.push_back(true);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(false);
    derivative_directions.push_back(-1);
    use_reciprocal.push_back(true);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(false);
    derivative_directions.push_back(-1);
    use_reciprocal.push_back(true);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(false);
    derivative_directions.push_back(-1);
    use_reciprocal.push_back(true);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(true);
    derivative_directions.push_back(0);
    use_reciprocal.push_back(false);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(true);
    derivative_directions.push_back(0);
    use_reciprocal.push_back(false);
    
    quantity_names.push_back("PRESSURE");
    component_indices.push_back(0);
    use_derivative.push_back(true);
    derivative_directions.push_back(2);
    use_reciprocal.push_back(false);
    
    quantity_names.push_back("PRESSURE");
    component_indices.push_back(0);
    use_derivative.push_back(true);
    derivative_directions.push_back(2);
    use_reciprocal.push_back(false);
    
    std::vector<double> B2_part_2 = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        use_derivative,
        derivative_directions,
        use_reciprocal,
        d_num_ghosts_derivative,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    use_derivative.clear();
    derivative_directions.clear();
    use_reciprocal.clear();
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(false);
    derivative_directions.push_back(-1);
    use_reciprocal.push_back(true);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(false);
    derivative_directions.push_back(-1);
    use_reciprocal.push_back(true);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(false);
    derivative_directions.push_back(-1);
    use_reciprocal.push_back(true);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(false);
    derivative_directions.push_back(-1);
    use_reciprocal.push_back(true);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(true);
    derivative_directions.push_back(0);
    use_reciprocal.push_back(false);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(true);
    derivative_directions.push_back(2);
    use_reciprocal.push_back(false);
    
    quantity_names.push_back("PRESSURE");
    component_indices.push_back(0);
    use_derivative.push_back(true);
    derivative_directions.push_back(0);
    use_reciprocal.push_back(false);
    
    quantity_names.push_back("PRESSURE");
    component_indices.push_back(0);
    use_derivative.push_back(true);
    derivative_directions.push_back(2);
    use_reciprocal.push_back(false);
    
    std::vector<double> B2_part_3 = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        use_derivative,
        derivative_directions,
        use_reciprocal,
        d_num_ghosts_derivative,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    use_derivative.clear();
    derivative_directions.clear();
    use_reciprocal.clear();
    
    std::vector<double> B2_sq_avg(finest_level_dims[0], double(0));
    for (int i = 0; i < finest_level_dims[0]; i++)
    {
        B2_sq_avg[i] = B2_part_1[i] + B2_part_2[i] - double(2)*B2_part_3[i];
    }
    
    std::vector<std::vector<double> >& B2_sq_avg_realizations = d_ensemble_statistics->B2_sq_avg_realizations;
    B2_sq_avg_realizations.push_back(B2_sq_avg);
    
    d_ensemble_statistics->B2_sq_avg_computed = true;
}


/*
 * Computed averaged quantities for TKE dissipation with assumed homogeneity in y-direction (2D) or yz-plane (3D).
 */
void
RTIRMIStatisticsUtilities::computeAveragedQuantiitesForTKEDissipationWithHomogeneityInYDirectionOrInYZPlane(
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
    std::vector<bool> use_derivative;
    std::vector<int> derivative_directions;
    std::vector<bool> use_reciprocal;
    
    if (!d_ensemble_statistics->ddx_u_avg_computed)
    {
        std::vector<std::vector<double> >& ddx_u_avg_realizations = d_ensemble_statistics->ddx_u_avg_realizations;
        
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(0);
        use_derivative.push_back(true);
        derivative_directions.push_back(0);
        
        std::vector<double> ddx_u_avg = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
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
        
        ddx_u_avg_realizations.push_back(ddx_u_avg);
        
        d_ensemble_statistics->ddx_u_avg_computed = true;
    }
    
    if (!d_ensemble_statistics->tau11_avg_computed)
    {
        std::vector<std::vector<double> >& tau11_avg_realizations = d_ensemble_statistics->tau11_avg_realizations;
        
        std::vector<double> tau11_avg = getAveragedShearStressComponentWithInhomogeneousXDirection(
            0,
            patch_hierarchy,
            data_context);
        
        tau11_avg_realizations.push_back(tau11_avg);
        
        d_ensemble_statistics->tau11_avg_computed = true;
    }
    
    if (!d_ensemble_statistics->tau11_ddx_u_avg_computed)
    {
        std::vector<std::vector<double> >& tau11_ddx_u_avg_realizations = d_ensemble_statistics->tau11_ddx_u_avg_realizations;
        
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(0);
        use_derivative.push_back(true);
        derivative_directions.push_back(0);
        use_reciprocal.push_back(false);
        
        std::vector<double> tau11_ddx_u_avg = getAveragedQuantityWithShearStressComponentWithInhomogeneousXDirection(
            quantity_names,
            component_indices,
            use_derivative,
            derivative_directions,
            use_reciprocal,
            0,
            patch_hierarchy,
            data_context);
        
        quantity_names.clear();
        component_indices.clear();
        use_derivative.clear();
        derivative_directions.clear();
        use_reciprocal.clear();
        
        tau11_ddx_u_avg_realizations.push_back(tau11_ddx_u_avg);
        
        d_ensemble_statistics->tau11_ddx_u_avg_computed = true;
    }
    
    if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
    {
        if (!d_ensemble_statistics->ddy_u_avg_computed)
        {
            std::vector<std::vector<double> >& ddy_u_avg_realizations = d_ensemble_statistics->ddy_u_avg_realizations;
            
            quantity_names.push_back("VELOCITY");
            component_indices.push_back(0);
            use_derivative.push_back(true);
            derivative_directions.push_back(1);
            
            std::vector<double> ddy_u_avg = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
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
            
            ddy_u_avg_realizations.push_back(ddy_u_avg);
            
            d_ensemble_statistics->ddy_u_avg_computed = true;
        }
        
        if (!d_ensemble_statistics->ddx_v_avg_computed)
        {
            std::vector<std::vector<double> >& ddx_v_avg_realizations = d_ensemble_statistics->ddx_v_avg_realizations;
            
            quantity_names.push_back("VELOCITY");
            component_indices.push_back(1);
            use_derivative.push_back(true);
            derivative_directions.push_back(0);
            
            std::vector<double> ddx_v_avg = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
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
            
            ddx_v_avg_realizations.push_back(ddx_v_avg);
            
            d_ensemble_statistics->ddx_v_avg_computed = true;
        }
        
        if (!d_ensemble_statistics->ddy_v_avg_computed)
        {
            std::vector<std::vector<double> >& ddy_v_avg_realizations = d_ensemble_statistics->ddy_v_avg_realizations;
            
            quantity_names.push_back("VELOCITY");
            component_indices.push_back(1);
            use_derivative.push_back(true);
            derivative_directions.push_back(1);
            
            std::vector<double> ddy_v_avg = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
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
            
            ddy_v_avg_realizations.push_back(ddy_v_avg);
            
            d_ensemble_statistics->ddy_v_avg_computed = true;
        }
        
        if (!d_ensemble_statistics->tau12_avg_computed)
        {
            std::vector<std::vector<double> >& tau12_avg_realizations = d_ensemble_statistics->tau12_avg_realizations;
            
            std::vector<double> tau12_avg = getAveragedShearStressComponentWithInhomogeneousXDirection(
                1,
                patch_hierarchy,
                data_context);
            
            tau12_avg_realizations.push_back(tau12_avg);
            
            d_ensemble_statistics->tau12_avg_computed = true;
        }
        
        if (!d_ensemble_statistics->tau22_avg_computed)
        {
            std::vector<std::vector<double> >& tau22_avg_realizations = d_ensemble_statistics->tau22_avg_realizations;
            
            std::vector<double> tau22_avg = getAveragedShearStressComponentWithInhomogeneousXDirection(
                3,
                patch_hierarchy,
                data_context);
            
            tau22_avg_realizations.push_back(tau22_avg);
            
            d_ensemble_statistics->tau22_avg_computed = true;
        }
        
        if (!d_ensemble_statistics->tau12_ddy_u_avg_computed)
        {
            std::vector<std::vector<double> >& tau12_ddy_u_avg_realizations = d_ensemble_statistics->tau12_ddy_u_avg_realizations;
            
            quantity_names.push_back("VELOCITY");
            component_indices.push_back(0);
            use_derivative.push_back(true);
            derivative_directions.push_back(1);
            use_reciprocal.push_back(false);
            
            std::vector<double> tau12_ddy_u_avg = getAveragedQuantityWithShearStressComponentWithInhomogeneousXDirection(
                quantity_names,
                component_indices,
                use_derivative,
                derivative_directions,
                use_reciprocal,
                1,
                patch_hierarchy,
                data_context);
            
            quantity_names.clear();
            component_indices.clear();
            use_derivative.clear();
            derivative_directions.clear();
            use_reciprocal.clear();
            
            tau12_ddy_u_avg_realizations.push_back(tau12_ddy_u_avg);
            
            d_ensemble_statistics->tau12_ddy_u_avg_computed = true;
        }
        
        if (!d_ensemble_statistics->tau12_ddx_v_avg_computed)
        {
            std::vector<std::vector<double> >& tau12_ddx_v_avg_realizations = d_ensemble_statistics->tau12_ddx_v_avg_realizations;
            
            quantity_names.push_back("VELOCITY");
            component_indices.push_back(1);
            use_derivative.push_back(true);
            derivative_directions.push_back(0);
            use_reciprocal.push_back(false);
            
            std::vector<double> tau12_ddx_v_avg = getAveragedQuantityWithShearStressComponentWithInhomogeneousXDirection(
                quantity_names,
                component_indices,
                use_derivative,
                derivative_directions,
                use_reciprocal,
                1,
                patch_hierarchy,
                data_context);
            
            quantity_names.clear();
            component_indices.clear();
            use_derivative.clear();
            derivative_directions.clear();
            use_reciprocal.clear();
            
            tau12_ddx_v_avg_realizations.push_back(tau12_ddx_v_avg);
            
            d_ensemble_statistics->tau12_ddx_v_avg_computed = true;
        }
        
        if (!d_ensemble_statistics->tau22_ddy_v_avg_computed)
        {
            std::vector<std::vector<double> >& tau22_ddy_v_avg_realizations = d_ensemble_statistics->tau22_ddy_v_avg_realizations;
            
            quantity_names.push_back("VELOCITY");
            component_indices.push_back(1);
            use_derivative.push_back(true);
            derivative_directions.push_back(1);
            use_reciprocal.push_back(false);
            
            std::vector<double> tau22_ddy_v_avg = getAveragedQuantityWithShearStressComponentWithInhomogeneousXDirection(
                quantity_names,
                component_indices,
                use_derivative,
                derivative_directions,
                use_reciprocal,
                3,
                patch_hierarchy,
                data_context);
            
            quantity_names.clear();
            component_indices.clear();
            use_derivative.clear();
            derivative_directions.clear();
            use_reciprocal.clear();
            
            tau22_ddy_v_avg_realizations.push_back(tau22_ddy_v_avg);
            
            d_ensemble_statistics->tau22_ddy_v_avg_computed = true;
        }
    } // if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
    
    if (d_dim == tbox::Dimension(3))
    {
        if (!d_ensemble_statistics->ddz_u_avg_computed)
        {
            std::vector<std::vector<double> >& ddz_u_avg_realizations = d_ensemble_statistics->ddz_u_avg_realizations;
            
            quantity_names.push_back("VELOCITY");
            component_indices.push_back(0);
            use_derivative.push_back(true);
            derivative_directions.push_back(2);
            
            std::vector<double> ddz_u_avg = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
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
            
            ddz_u_avg_realizations.push_back(ddz_u_avg);
            
            d_ensemble_statistics->ddz_u_avg_computed = true;
        }
        
        if (!d_ensemble_statistics->ddz_v_avg_computed)
        {
            std::vector<std::vector<double> >& ddz_v_avg_realizations = d_ensemble_statistics->ddz_v_avg_realizations;
            
            quantity_names.push_back("VELOCITY");
            component_indices.push_back(1);
            use_derivative.push_back(true);
            derivative_directions.push_back(2);
            
            std::vector<double> ddz_v_avg = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
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
            
            ddz_v_avg_realizations.push_back(ddz_v_avg);
            
            d_ensemble_statistics->ddz_v_avg_computed = true;
        }
        
        if (!d_ensemble_statistics->ddx_w_avg_computed)
        {
            std::vector<std::vector<double> >& ddx_w_avg_realizations = d_ensemble_statistics->ddx_w_avg_realizations;
            
            quantity_names.push_back("VELOCITY");
            component_indices.push_back(2);
            use_derivative.push_back(true);
            derivative_directions.push_back(0);
            
            std::vector<double> ddx_w_avg = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
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
            
            ddx_w_avg_realizations.push_back(ddx_w_avg);
            
            d_ensemble_statistics->ddx_w_avg_computed = true;
        }
        
        if (!d_ensemble_statistics->ddy_w_avg_computed)
        {
            std::vector<std::vector<double> >& ddy_w_avg_realizations = d_ensemble_statistics->ddy_w_avg_realizations;
            
            quantity_names.push_back("VELOCITY");
            component_indices.push_back(2);
            use_derivative.push_back(true);
            derivative_directions.push_back(1);
            
            std::vector<double> ddy_w_avg = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
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
            
            ddy_w_avg_realizations.push_back(ddy_w_avg);
            
            d_ensemble_statistics->ddy_w_avg_computed = true;
        }
        
        if (!d_ensemble_statistics->ddz_w_avg_computed)
        {
            std::vector<std::vector<double> >& ddz_w_avg_realizations = d_ensemble_statistics->ddz_w_avg_realizations;
            
            quantity_names.push_back("VELOCITY");
            component_indices.push_back(2);
            use_derivative.push_back(true);
            derivative_directions.push_back(2);
            
            std::vector<double> ddz_w_avg = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
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
            
            ddz_w_avg_realizations.push_back(ddz_w_avg);
            
            d_ensemble_statistics->ddz_w_avg_computed = true;
        }
        
        if (!d_ensemble_statistics->tau13_avg_computed)
        {
            std::vector<std::vector<double> >& tau13_avg_realizations = d_ensemble_statistics->tau13_avg_realizations;
            
            std::vector<double> tau13_avg = getAveragedShearStressComponentWithInhomogeneousXDirection(
                2,
                patch_hierarchy,
                data_context);
            
            tau13_avg_realizations.push_back(tau13_avg);
            
            d_ensemble_statistics->tau13_avg_computed = true;
        }
        
        if (!d_ensemble_statistics->tau23_avg_computed)
        {
            std::vector<std::vector<double> >& tau23_avg_realizations = d_ensemble_statistics->tau23_avg_realizations;
            
            std::vector<double> tau23_avg = getAveragedShearStressComponentWithInhomogeneousXDirection(
                4,
                patch_hierarchy,
                data_context);
            
            tau23_avg_realizations.push_back(tau23_avg);
            
            d_ensemble_statistics->tau23_avg_computed = true;
        }
        
        if (!d_ensemble_statistics->tau33_avg_computed)
        {
            std::vector<std::vector<double> >& tau33_avg_realizations = d_ensemble_statistics->tau33_avg_realizations;
            
            std::vector<double> tau33_avg = getAveragedShearStressComponentWithInhomogeneousXDirection(
                5,
                patch_hierarchy,
                data_context);
            
            tau33_avg_realizations.push_back(tau33_avg);
            
            d_ensemble_statistics->tau33_avg_computed = true;
        }
        
        if (!d_ensemble_statistics->tau13_ddz_u_avg_computed)
        {
            std::vector<std::vector<double> >& tau13_ddz_u_avg_realizations = d_ensemble_statistics->tau13_ddz_u_avg_realizations;
            
            quantity_names.push_back("VELOCITY");
            component_indices.push_back(0);
            use_derivative.push_back(true);
            derivative_directions.push_back(2);
            use_reciprocal.push_back(false);
            
            std::vector<double> tau13_ddz_u_avg = getAveragedQuantityWithShearStressComponentWithInhomogeneousXDirection(
                quantity_names,
                component_indices,
                use_derivative,
                derivative_directions,
                use_reciprocal,
                2,
                patch_hierarchy,
                data_context);
            
            quantity_names.clear();
            component_indices.clear();
            use_derivative.clear();
            derivative_directions.clear();
            use_reciprocal.clear();
            
            tau13_ddz_u_avg_realizations.push_back(tau13_ddz_u_avg);
            
            d_ensemble_statistics->tau13_ddz_u_avg_computed = true;
        }
        
        if (!d_ensemble_statistics->tau23_ddz_v_avg_computed)
        {
            std::vector<std::vector<double> >& tau23_ddz_v_avg_realizations = d_ensemble_statistics->tau23_ddz_v_avg_realizations;
            
            quantity_names.push_back("VELOCITY");
            component_indices.push_back(1);
            use_derivative.push_back(true);
            derivative_directions.push_back(2);
            use_reciprocal.push_back(false);
            
            std::vector<double> tau23_ddz_v_avg = getAveragedQuantityWithShearStressComponentWithInhomogeneousXDirection(
                quantity_names,
                component_indices,
                use_derivative,
                derivative_directions,
                use_reciprocal,
                4,
                patch_hierarchy,
                data_context);
            
            quantity_names.clear();
            component_indices.clear();
            use_derivative.clear();
            derivative_directions.clear();
            use_reciprocal.clear();
            
            tau23_ddz_v_avg_realizations.push_back(tau23_ddz_v_avg);
            
            d_ensemble_statistics->tau23_ddz_v_avg_computed = true;
        }
        
        if (!d_ensemble_statistics->tau13_ddx_w_avg_computed)
        {
            std::vector<std::vector<double> >& tau13_ddx_w_avg_realizations = d_ensemble_statistics->tau13_ddx_w_avg_realizations;
            
            quantity_names.push_back("VELOCITY");
            component_indices.push_back(2);
            use_derivative.push_back(true);
            derivative_directions.push_back(0);
            use_reciprocal.push_back(false);
            
            std::vector<double> tau13_ddx_w_avg = getAveragedQuantityWithShearStressComponentWithInhomogeneousXDirection(
                quantity_names,
                component_indices,
                use_derivative,
                derivative_directions,
                use_reciprocal,
                2,
                patch_hierarchy,
                data_context);
            
            quantity_names.clear();
            component_indices.clear();
            use_derivative.clear();
            derivative_directions.clear();
            use_reciprocal.clear();
            
            tau13_ddx_w_avg_realizations.push_back(tau13_ddx_w_avg);
            
            d_ensemble_statistics->tau13_ddx_w_avg_computed = true;
        }
        
        if (!d_ensemble_statistics->tau23_ddy_w_avg_computed)
        {
            std::vector<std::vector<double> >& tau23_ddy_w_avg_realizations = d_ensemble_statistics->tau23_ddy_w_avg_realizations;
            
            quantity_names.push_back("VELOCITY");
            component_indices.push_back(2);
            use_derivative.push_back(true);
            derivative_directions.push_back(1);
            use_reciprocal.push_back(false);
            
            std::vector<double> tau23_ddy_w_avg = getAveragedQuantityWithShearStressComponentWithInhomogeneousXDirection(
                quantity_names,
                component_indices,
                use_derivative,
                derivative_directions,
                use_reciprocal,
                4,
                patch_hierarchy,
                data_context);
            
            quantity_names.clear();
            component_indices.clear();
            use_derivative.clear();
            derivative_directions.clear();
            use_reciprocal.clear();
            
            tau23_ddy_w_avg_realizations.push_back(tau23_ddy_w_avg);
            
            d_ensemble_statistics->tau23_ddy_w_avg_computed = true;
        }
        
        if (!d_ensemble_statistics->tau33_ddz_w_avg_computed)
        {
            std::vector<std::vector<double> >& tau33_ddz_w_avg_realizations = d_ensemble_statistics->tau33_ddz_w_avg_realizations;
            
            quantity_names.push_back("VELOCITY");
            component_indices.push_back(2);
            use_derivative.push_back(true);
            derivative_directions.push_back(2);
            use_reciprocal.push_back(false);
            
            std::vector<double> tau33_ddz_w_avg = getAveragedQuantityWithShearStressComponentWithInhomogeneousXDirection(
                quantity_names,
                component_indices,
                use_derivative,
                derivative_directions,
                use_reciprocal,
                5,
                patch_hierarchy,
                data_context);
            
            quantity_names.clear();
            component_indices.clear();
            use_derivative.clear();
            derivative_directions.clear();
            use_reciprocal.clear();
            
            tau33_ddz_w_avg_realizations.push_back(tau33_ddz_w_avg);
            
            d_ensemble_statistics->tau33_ddz_w_avg_computed = true;
        }
    } // if (d_dim == tbox::Dimension(3))
    
    d_ensemble_statistics->TKE_diss_computed = true;
}


/*
 * Output spatial profile of ensemble averaged mass fraction with assumed homogeneity in y-direction (2D) or
 * yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputSpatialProfileEnsembleAveragedMassFractionWithHomogeneityInYDirectionOrInYZPlane(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Output the spatial profile (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_out;
        
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
 * Output spatial profile of ensemble averaged mole fraction with assumed homogeneity in y-direction (2D) or
 * yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputSpatialProfileEnsembleAveragedMoleFractionWithHomogeneityInYDirectionOrInYZPlane(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Output the spatial profile (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_out;
        
        f_out.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        const std::vector<std::vector<double> >& X_0_avg_realizations =
            d_ensemble_statistics->X_0_avg_realizations;
        
        const int num_realizations = static_cast<int>(X_0_avg_realizations.size());
        
        TBOX_ASSERT(d_ensemble_statistics->getNumberOfEnsembles() == num_realizations);
        TBOX_ASSERT(num_realizations > 0);
        
        const int num_cells = static_cast<int>(X_0_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> X_0_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                X_0_avg_global[i] += weight*X_0_avg_realizations[ri][i];
            }
        }
        
        f_out.write((char*)&output_time, sizeof(double));
        f_out.write((char*)&X_0_avg_global[0], sizeof(double)*X_0_avg_global.size());
        
        f_out.close();
    }
}


/*
 * Output spatial profile of ensemble averaged volume fraction with assumed homogeneity in y-direction (2D) or
 * yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputSpatialProfileEnsembleAveragedVolumeFractionWithHomogeneityInYDirectionOrInYZPlane(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Output the spatial profile (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_out;
        
        f_out.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        const std::vector<std::vector<double> >& Z_0_avg_realizations =
            d_ensemble_statistics->Z_0_avg_realizations;
        
        const int num_realizations = static_cast<int>(Z_0_avg_realizations.size());
        
        TBOX_ASSERT(d_ensemble_statistics->getNumberOfEnsembles() == num_realizations);
        TBOX_ASSERT(num_realizations > 0);
        
        const int num_cells = static_cast<int>(Z_0_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> Z_0_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                Z_0_avg_global[i] += weight*Z_0_avg_realizations[ri][i];
            }
        }
        
        f_out.write((char*)&output_time, sizeof(double));
        f_out.write((char*)&Z_0_avg_global[0], sizeof(double)*Z_0_avg_global.size());
        
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
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Output the spatial profile (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_out;
        
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
 * Output spatial profile of ensemble variance of mole fraction with assumed homogeneity in y-direction (2D) or
 * yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputSpatialProfileEnsembleMoleFractionVarianceWithHomogeneityInYDirectionOrInYZPlane(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Output the spatial profile (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_out;
        
        f_out.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        const std::vector<std::vector<double> >& X_0_avg_realizations =
            d_ensemble_statistics->X_0_avg_realizations;
        
        const std::vector<std::vector<double> >& X_0_sq_avg_realizations = 
            d_ensemble_statistics->X_0_sq_avg_realizations;
        
        const int num_realizations = static_cast<int>(X_0_avg_realizations.size());
        
        TBOX_ASSERT(d_ensemble_statistics->getNumberOfEnsembles() == num_realizations);
        TBOX_ASSERT(num_realizations > 0);
        TBOX_ASSERT(num_realizations == static_cast<int>(X_0_sq_avg_realizations.size()));
        
        const int num_cells = static_cast<int>(X_0_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> X_0_avg_global(num_cells, double(0));
        std::vector<double> X_0_sq_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                X_0_avg_global[i]    += weight*X_0_avg_realizations[ri][i];
                X_0_sq_avg_global[i] += weight*X_0_sq_avg_realizations[ri][i];
            }
        }
        
        std::vector<double> X_0_var_global(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            X_0_var_global[i] = X_0_sq_avg_global[i] - X_0_avg_global[i]*X_0_avg_global[i];
        }
        
        f_out.write((char*)&output_time, sizeof(double));
        f_out.write((char*)&X_0_var_global[0], sizeof(double)*X_0_var_global.size());
        
        f_out.close();
    }
}


/*
 * Output spatial profile of ensemble variance of volume fraction with assumed homogeneity in y-direction (2D) or
 * yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputSpatialProfileEnsembleVolumeFractionVarianceWithHomogeneityInYDirectionOrInYZPlane(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Output the spatial profile (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_out;
        
        f_out.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        const std::vector<std::vector<double> >& Z_0_avg_realizations =
            d_ensemble_statistics->Z_0_avg_realizations;
        
        const std::vector<std::vector<double> >& Z_0_sq_avg_realizations = 
            d_ensemble_statistics->Z_0_sq_avg_realizations;
        
        const int num_realizations = static_cast<int>(Z_0_avg_realizations.size());
        
        TBOX_ASSERT(d_ensemble_statistics->getNumberOfEnsembles() == num_realizations);
        TBOX_ASSERT(num_realizations > 0);
        TBOX_ASSERT(num_realizations == static_cast<int>(Z_0_sq_avg_realizations.size()));
        
        const int num_cells = static_cast<int>(Z_0_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> Z_0_avg_global(num_cells, double(0));
        std::vector<double> Z_0_sq_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                Z_0_avg_global[i]    += weight*Z_0_avg_realizations[ri][i];
                Z_0_sq_avg_global[i] += weight*Z_0_sq_avg_realizations[ri][i];
            }
        }
        
        std::vector<double> Z_0_var_global(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            Z_0_var_global[i] = Z_0_sq_avg_global[i] - Z_0_avg_global[i]*Z_0_avg_global[i];
        }
        
        f_out.write((char*)&output_time, sizeof(double));
        f_out.write((char*)&Z_0_var_global[0], sizeof(double)*Z_0_var_global.size());
        
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
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Output the spatial profile (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_out;
        
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
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Output the spatial profile (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_out;
        
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
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Output the spatial profile (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_out;
        
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
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Output the spatial profile (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_out;
        
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
 * Output spatial profile of ensemble averaged pressure with assumed homogeneity in y-direction (2D) or
 * yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputSpatialProfileEnsembleAveragedPressureWithHomogeneityInYDirectionOrInYZPlane(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Output the spatial profile (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_out;
        
        f_out.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        const std::vector<std::vector<double> >& p_avg_realizations =
            d_ensemble_statistics->p_avg_realizations;
        
        const int num_realizations = static_cast<int>(p_avg_realizations.size());
        
        TBOX_ASSERT(d_ensemble_statistics->getNumberOfEnsembles() == num_realizations);
        TBOX_ASSERT(num_realizations > 0);
        
        const int num_cells = static_cast<int>(p_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> p_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                p_avg_global[i] += weight*p_avg_realizations[ri][i];
            }
        }
        
        f_out.write((char*)&output_time, sizeof(double));
        f_out.write((char*)&p_avg_global[0], sizeof(double)*p_avg_global.size());
        
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
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Output the spatial profile (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_out;
        
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
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Output the spatial profile (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_out;
        
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
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Output the spatial profile (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_out;
        
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
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Output the spatial profile (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_out;
        
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
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Output the spatial profile (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_out;
        
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
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Output the spatial profile (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_out;
        
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
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Output the spatial profile (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_out;
        
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
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Output the spatial profile (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_out;
        
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
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Output the spatial profile (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_out;
        
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
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Output the spatial profile (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_out;
        
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
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Output the spatial profile (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_out;
        
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
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Output the spatial profile (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_out;
        
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
 * Output spatial profile of ensemble RMS of baroclinic torque in z-direction with assumed homogeneity in
 * y-direction (2D) or yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputSpatialProfileEnsembleRMSBaroclinicTorqueZWithHomogeneityInYDirectionOrInYZPlane(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Output the spatial profile (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_out;
        
        f_out.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        const std::vector<std::vector<double> >& B3_sq_avg_realizations =
            d_ensemble_statistics->B3_sq_avg_realizations;
        
        const int num_realizations = static_cast<int>(B3_sq_avg_realizations.size());
        
        TBOX_ASSERT(d_ensemble_statistics->getNumberOfEnsembles() == num_realizations);
        TBOX_ASSERT(num_realizations > 0);
        
        const int num_cells = static_cast<int>(B3_sq_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> B3_rms_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                B3_rms_global[i] += weight*B3_sq_avg_realizations[ri][i];
            }
        }
        
        for (int i = 0; i < num_cells; i++)
        {
            B3_rms_global[i] = sqrt(B3_rms_global[i]);
        }
        
        f_out.write((char*)&output_time, sizeof(double));
        f_out.write((char*)&B3_rms_global[0], sizeof(double)*B3_rms_global.size());
        
        f_out.close();
    }
}


/*
 * Output spatial profile of ensemble RMS of baroclinic torque in x-direction with assumed homogeneity in
 * yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputSpatialProfileEnsembleRMSBaroclinicTorqueXWithHomogeneityInYZPlane(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Output the spatial profile (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_out;
        
        f_out.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        const std::vector<std::vector<double> >& B1_sq_avg_realizations =
            d_ensemble_statistics->B1_sq_avg_realizations;
        
        const int num_realizations = static_cast<int>(B1_sq_avg_realizations.size());
        
        TBOX_ASSERT(d_ensemble_statistics->getNumberOfEnsembles() == num_realizations);
        TBOX_ASSERT(num_realizations > 0);
        
        const int num_cells = static_cast<int>(B1_sq_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> B1_rms_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                B1_rms_global[i] += weight*B1_sq_avg_realizations[ri][i];
            }
        }
        
        for (int i = 0; i < num_cells; i++)
        {
            B1_rms_global[i] = sqrt(B1_rms_global[i]);
        }
        
        f_out.write((char*)&output_time, sizeof(double));
        f_out.write((char*)&B1_rms_global[0], sizeof(double)*B1_rms_global.size());
        
        f_out.close();
    }
}


/*
 * Output spatial profile of ensemble RMS of baroclinic torque in y-direction with assumed homogeneity in
 * yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputSpatialProfileEnsembleRMSBaroclinicTorqueYWithHomogeneityInYZPlane(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Output the spatial profile (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_out;
        
        f_out.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        const std::vector<std::vector<double> >& B2_sq_avg_realizations =
            d_ensemble_statistics->B2_sq_avg_realizations;
        
        const int num_realizations = static_cast<int>(B2_sq_avg_realizations.size());
        
        TBOX_ASSERT(d_ensemble_statistics->getNumberOfEnsembles() == num_realizations);
        TBOX_ASSERT(num_realizations > 0);
        
        const int num_cells = static_cast<int>(B2_sq_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> B2_rms_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                B2_rms_global[i] += weight*B2_sq_avg_realizations[ri][i];
            }
        }
        
        for (int i = 0; i < num_cells; i++)
        {
            B2_rms_global[i] = sqrt(B2_rms_global[i]);
        }
        
        f_out.write((char*)&output_time, sizeof(double));
        f_out.write((char*)&B2_rms_global[0], sizeof(double)*B2_rms_global.size());
        
        f_out.close();
    }
}


/*
 * Output ensemble mixing width in x-direction to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleMixingWidthInXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const
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
 * Output ensemble mixing width in x-direction using mole fractions to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleMixingWidthInXDirectionWithMoleFractions(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'MIXING_WIDTH_X_MOL_F' can be computed with two species only."
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
        
        const std::vector<std::vector<double> >& X_0_avg_realizations =
            d_ensemble_statistics->X_0_avg_realizations;
        
        const int num_realizations = static_cast<int>(X_0_avg_realizations.size());
        
        TBOX_ASSERT(d_ensemble_statistics->getNumberOfEnsembles() == num_realizations);
        TBOX_ASSERT(num_realizations > 0);
        
        const int num_cells = static_cast<int>(X_0_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> X_0_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                X_0_avg_global[i] += weight*X_0_avg_realizations[ri][i];
            }
        }
        
        double W = double(0);
        
        for (int i = 0; i < num_cells; i++)
        {
            W += X_0_avg_global[i]*(double(1) - X_0_avg_global[i]);
        }
        
        W = double(4)*W*dx_finest[0];
        
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << W;
        
        f_out.close();
    }
}


/*
 * Output ensemble mixing width in x-direction using volume fractions to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleMixingWidthInXDirectionWithVolumeFractions(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'MIXING_WIDTH_X_VOL_F' can be computed with two species only."
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
        
        const std::vector<std::vector<double> >& Z_0_avg_realizations =
            d_ensemble_statistics->Z_0_avg_realizations;
        
        const int num_realizations = static_cast<int>(Z_0_avg_realizations.size());
        
        TBOX_ASSERT(d_ensemble_statistics->getNumberOfEnsembles() == num_realizations);
        TBOX_ASSERT(num_realizations > 0);
        
        const int num_cells = static_cast<int>(Z_0_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> Z_0_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                Z_0_avg_global[i] += weight*Z_0_avg_realizations[ri][i];
            }
        }
        
        double W = double(0);
        
        for (int i = 0; i < num_cells; i++)
        {
            W += Z_0_avg_global[i]*(double(1) - Z_0_avg_global[i]);
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
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const
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
 * Output ensemble minimum interface location in x-direction using mole fractions to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleInterfaceMinInXDirectionWithMoleFractions(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'INTERFACE_MIN_X_MOL_F' can be computed with two species only."
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
        
        const std::vector<std::vector<double> >& X_0_avg_realizations =
            d_ensemble_statistics->X_0_avg_realizations;
        
        const int num_realizations = static_cast<int>(X_0_avg_realizations.size());
        
        TBOX_ASSERT(d_ensemble_statistics->getNumberOfEnsembles() == num_realizations);
        TBOX_ASSERT(num_realizations > 0);
        
        const int num_cells = static_cast<int>(X_0_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> X_0_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                X_0_avg_global[i] += weight*X_0_avg_realizations[ri][i];
            }
        }
        
        const double* x_lo = d_grid_geometry->getXLower();
        const double* x_hi = d_grid_geometry->getXUpper();
        double interface_min = x_hi[0];
        
        for (int i = num_cells - 1; i >= 0;  i--)
        {
            if (X_0_avg_global[i] > 0.01 && X_0_avg_global[i] < 0.99)
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
 * Output ensemble minimum interface location in x-direction using volume fractions to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleInterfaceMinInXDirectionWithVolumeFractions(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'INTERFACE_MIN_X_VOL_F' can be computed with two species only."
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
        
        const std::vector<std::vector<double> >& Z_0_avg_realizations =
            d_ensemble_statistics->Z_0_avg_realizations;
        
        const int num_realizations = static_cast<int>(Z_0_avg_realizations.size());
        
        TBOX_ASSERT(d_ensemble_statistics->getNumberOfEnsembles() == num_realizations);
        TBOX_ASSERT(num_realizations > 0);
        
        const int num_cells = static_cast<int>(Z_0_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> Z_0_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                Z_0_avg_global[i] += weight*Z_0_avg_realizations[ri][i];
            }
        }
        
        const double* x_lo = d_grid_geometry->getXLower();
        const double* x_hi = d_grid_geometry->getXUpper();
        double interface_min = x_hi[0];
        
        for (int i = num_cells - 1; i >= 0;  i--)
        {
            if (Z_0_avg_global[i] > 0.01 && Z_0_avg_global[i] < 0.99)
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
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'INTERFACE_MAX_X' can be computed with two species only."
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
 * Output ensemble maximum interface location in x-direction using mole fractions to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleInterfaceMaxInXDirectionWithMoleFractions(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'INTERFACE_MAX_X_MOL_F' can be computed with two species only."
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
        
        const std::vector<std::vector<double> >& X_0_avg_realizations =
            d_ensemble_statistics->X_0_avg_realizations;
        
        const int num_realizations = static_cast<int>(X_0_avg_realizations.size());
        
        TBOX_ASSERT(d_ensemble_statistics->getNumberOfEnsembles() == num_realizations);
        TBOX_ASSERT(num_realizations > 0);
        
        const int num_cells = static_cast<int>(X_0_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> X_0_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                X_0_avg_global[i] += weight*X_0_avg_realizations[ri][i];
            }
        }
        
        const double* x_lo = d_grid_geometry->getXLower();
        // const double* x_hi = d_grid_geometry->getXUpper();
        double interface_max = x_lo[0];
        
        for (int i = 0; i < num_cells; i++)
        {
            if (X_0_avg_global[i] > 0.01 && X_0_avg_global[i] < 0.99)
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
 * Output ensemble maximum interface location in x-direction using volume fractions to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleInterfaceMaxInXDirectionWithVolumeFractions(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'INTERFACE_MAX_X_VOL_F' can be computed with two species only."
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
        
        const std::vector<std::vector<double> >& Z_0_avg_realizations =
            d_ensemble_statistics->Z_0_avg_realizations;
        
        const int num_realizations = static_cast<int>(Z_0_avg_realizations.size());
        
        TBOX_ASSERT(d_ensemble_statistics->getNumberOfEnsembles() == num_realizations);
        TBOX_ASSERT(num_realizations > 0);
        
        const int num_cells = static_cast<int>(Z_0_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> Z_0_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                Z_0_avg_global[i] += weight*Z_0_avg_realizations[ri][i];
            }
        }
        
        const double* x_lo = d_grid_geometry->getXLower();
        // const double* x_hi = d_grid_geometry->getXUpper();
        double interface_max = x_lo[0];
        
        for (int i = 0; i < num_cells; i++)
        {
            if (Z_0_avg_global[i] > 0.01 && Z_0_avg_global[i] < 0.99)
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
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const
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
 * Output ensemble mixedness in x-direction using mole fractions to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleMixednessInXDirectionWithMoleFractions(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'MIXEDNESS_X_MOL_F' can be computed with two species only."
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
        
        const std::vector<std::vector<double> >& X_0_avg_realizations =
            d_ensemble_statistics->X_0_avg_realizations;
        
        const std::vector<std::vector<double> >& X_0_X_1_avg_realizations =
            d_ensemble_statistics->X_0_X_1_avg_realizations;
        
        const int num_realizations = static_cast<int>(X_0_avg_realizations.size());
        
        TBOX_ASSERT(d_ensemble_statistics->getNumberOfEnsembles() == num_realizations);
        TBOX_ASSERT(num_realizations > 0);
        TBOX_ASSERT(num_realizations == static_cast<int>(X_0_X_1_avg_realizations.size()));
        
        const int num_cells = static_cast<int>(X_0_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> X_0_avg_global(num_cells, double(0));
        std::vector<double> X_0_X_1_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                X_0_avg_global[i]     += weight*X_0_avg_realizations[ri][i];
                X_0_X_1_avg_global[i] += weight*X_0_X_1_avg_realizations[ri][i];
            }
        }
        
        double num = double(0);
        double den = double(0);
        
        for (int i = 0; i < num_cells; i++)
        {
            num += X_0_X_1_avg_global[i];
            den += X_0_avg_global[i]*(double(1) - X_0_avg_global[i]);
        }
        
        const double Theta = num/den;
        
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << Theta;
        
        f_out.close();
    }
}


/*
 * Output ensemble mixedness in x-direction using volume fractions to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleMixednessInXDirectionWithVolumeFractions(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'MIXEDNESS_X_VOL_F' can be computed with two species only."
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
        
        const std::vector<std::vector<double> >& Z_0_avg_realizations =
            d_ensemble_statistics->Z_0_avg_realizations;
        
        const std::vector<std::vector<double> >& Z_0_Z_1_avg_realizations =
            d_ensemble_statistics->Z_0_Z_1_avg_realizations;
        
        const int num_realizations = static_cast<int>(Z_0_avg_realizations.size());
        
        TBOX_ASSERT(d_ensemble_statistics->getNumberOfEnsembles() == num_realizations);
        TBOX_ASSERT(num_realizations > 0);
        TBOX_ASSERT(num_realizations == static_cast<int>(Z_0_Z_1_avg_realizations.size()));
        
        const int num_cells = static_cast<int>(Z_0_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> Z_0_avg_global(num_cells, double(0));
        std::vector<double> Z_0_Z_1_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                Z_0_avg_global[i]     += weight*Z_0_avg_realizations[ri][i];
                Z_0_Z_1_avg_global[i] += weight*Z_0_Z_1_avg_realizations[ri][i];
            }
        }
        
        double num = double(0);
        double den = double(0);
        
        for (int i = 0; i < num_cells; i++)
        {
            num += Z_0_Z_1_avg_global[i];
            den += Z_0_avg_global[i]*(double(1) - Z_0_avg_global[i]);
        }
        
        const double Theta = num/den;
        
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << Theta;
        
        f_out.close();
    }
}


/*
 * Output mean of ensemble density mean in mixing layer with with assumed homogeneity in y-direction (2D) or
 * yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleDensityMeanWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const
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
        
        const std::vector<std::vector<double> >& Y_0_avg_realizations =
            d_ensemble_statistics->Y_0_avg_realizations;
        
        const int num_realizations = d_ensemble_statistics->getNumberOfEnsembles();
        
        TBOX_ASSERT(num_realizations > 0);
        TBOX_ASSERT(num_realizations == static_cast<int>(rho_avg_realizations.size()));
        TBOX_ASSERT(num_realizations == static_cast<int>(Y_0_avg_realizations.size()));
        
        const int num_cells = static_cast<int>(rho_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> rho_avg_global(num_cells, double(0));
        std::vector<double> Y_0_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                rho_avg_global[i] += weight*rho_avg_realizations[ri][i];
                Y_0_avg_global[i] += weight*Y_0_avg_realizations[ri][i];
            }
        }
        
        double rho_sum = double(0);
        int count = 0;
        
        for (int i = 0; i < num_cells; i++)
        {
            const double mixing_metric = double(4)*Y_0_avg_global[i]*(double(1) - Y_0_avg_global[i]);
            if (mixing_metric > double(9)/double(10))
            {
                rho_sum += rho_avg_global[i];
                count++;
            }
        }
        
        const double rho_mean = rho_sum/count;
        
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << rho_mean;
        
        f_out.close();
    }
}


/*
 * Output mean of ensemble specific volume mean in mixing layer with with assumed homogeneity in
 * y-direction (2D) or yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleSpecificVolumeMeanWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const
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
        
        const std::vector<std::vector<double> >& rho_inv_avg_realizations =
            d_ensemble_statistics->rho_inv_avg_realizations;
        
        const std::vector<std::vector<double> >& Y_0_avg_realizations =
            d_ensemble_statistics->Y_0_avg_realizations;
        
        const int num_realizations = d_ensemble_statistics->getNumberOfEnsembles();
        
        TBOX_ASSERT(num_realizations > 0);
        TBOX_ASSERT(num_realizations == static_cast<int>(rho_inv_avg_realizations.size()));
        TBOX_ASSERT(num_realizations == static_cast<int>(Y_0_avg_realizations.size()));
        
        const int num_cells = static_cast<int>(rho_inv_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> rho_inv_avg_global(num_cells, double(0));
        std::vector<double> Y_0_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                rho_inv_avg_global[i] += weight*rho_inv_avg_realizations[ri][i];
                Y_0_avg_global[i]     += weight*Y_0_avg_realizations[ri][i];
            }
        }
        
        double rho_inv_sum = double(0);
        int count = 0;
        
        for (int i = 0; i < num_cells; i++)
        {
            const double mixing_metric = double(4)*Y_0_avg_global[i]*(double(1) - Y_0_avg_global[i]);
            if (mixing_metric > double(9)/double(10))
            {
                rho_inv_sum += rho_inv_avg_global[i];
                count++;
            }
        }
        
        const double rho_inv_mean = rho_inv_sum/count;
        
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << rho_inv_mean;
        
        f_out.close();
    }
}


/*
 * Output mean of ensemble pressure mean in mixing layer with with assumed homogeneity in y-direction (2D) or
 * yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsemblePressureMeanWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const
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
        
        const std::vector<std::vector<double> >& p_avg_realizations =
            d_ensemble_statistics->p_avg_realizations;
        
        const std::vector<std::vector<double> >& Y_0_avg_realizations =
            d_ensemble_statistics->Y_0_avg_realizations;
        
        const int num_realizations = d_ensemble_statistics->getNumberOfEnsembles();
        
        TBOX_ASSERT(num_realizations > 0);
        TBOX_ASSERT(num_realizations == static_cast<int>(p_avg_realizations.size()));
        TBOX_ASSERT(num_realizations == static_cast<int>(Y_0_avg_realizations.size()));
        
        const int num_cells = static_cast<int>(p_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> p_avg_global(num_cells, double(0));
        std::vector<double> Y_0_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                p_avg_global[i]   += weight*p_avg_realizations[ri][i];
                Y_0_avg_global[i] += weight*Y_0_avg_realizations[ri][i];
            }
        }
        
        double p_sum = double(0);
        int count = 0;
        
        for (int i = 0; i < num_cells; i++)
        {
            const double mixing_metric = double(4)*Y_0_avg_global[i]*(double(1) - Y_0_avg_global[i]);
            if (mixing_metric > double(9)/double(10))
            {
                p_sum += p_avg_global[i];
                count++;
            }
        }
        
        const double p_mean = p_sum/count;
        
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << p_mean;
        
        f_out.close();
    }
}


/*
 * Output mean of ensemble velocity associated with turbulent mass flux component in x-direction in
 * mixing layer with assumed homogeneity in y-direction (2D) or yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleTurbulentMassFluxVelocityXMeanWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const
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
 * Output mean of ensemble density-specific-volume covariance in mixing layer with assumed homogeneity in
 * y-direction (2D) or yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleDensitySpecificVolumeCovarianceMeanWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const
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
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const
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
 * Output mean of ensemble Reynolds normal stress component in x-direction in mixing layer with assumed
 * homogeneity in y-direction (2D) or yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleReynoldsNormalStressXMeanWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const
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
 * Output mean of ensemble Reynolds normal stress component in y-direction in mixing layer with assumed
 * homogeneity in y-direction (2D) or yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleReynoldsNormalStressYMeanWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const
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
 * Output mean of ensemble Reynolds normal stress component in z-direction in mixing layer with assumed
 * homogeneity in yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleReynoldsNormalStressZMeanWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const
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
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const
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
        
        f_out.close();
    }
}


/*
 * Output mean of ensemble enstrophy in mixing layer with assumed homogeneity in y-direction (2D) or yz-plane (3D)
 * to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleEnstrophyMeanWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const
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
        
        const std::vector<std::vector<double> >& Omega_avg_realizations =
            d_ensemble_statistics->Omega_avg_realizations;
        
        const std::vector<std::vector<double> >& Y_0_avg_realizations =
            d_ensemble_statistics->Y_0_avg_realizations;
        
        const int num_realizations = static_cast<int>(Omega_avg_realizations.size());
        
        TBOX_ASSERT(d_ensemble_statistics->getNumberOfEnsembles() == num_realizations);
        TBOX_ASSERT(num_realizations > 0);
        TBOX_ASSERT(num_realizations == static_cast<int>(Y_0_avg_realizations.size()));
        
        const int num_cells = static_cast<int>(Omega_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> Omega_avg_global(num_cells, double(0));
        std::vector<double> Y_0_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                Omega_avg_global[i] += weight*Omega_avg_realizations[ri][i];
                Y_0_avg_global[i]   += weight*Y_0_avg_realizations[ri][i];
            }
        }
        
        double Omega_sum = double(0);
        int count = 0;
        
        for (int i = 0; i < num_cells; i++)
        {
            const double mixing_metric = double(4)*Y_0_avg_global[i]*(double(1) - Y_0_avg_global[i]);
            if (mixing_metric > double(9)/double(10))
            {
                Omega_sum += Omega_avg_global[i];
                count++;
            }
        }
        
        const double Omega_mean = Omega_sum/count;
        
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << Omega_mean;
        
        f_out.close();
    }
}


/*
 * Output ensemble scalar dissipation rate of first species integrated to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleScalarDissipationRateIntegrated(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const
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
        
        f_out.close();
    }
}


/*
 * Output mean of ensemble scalar dissipation rate of first species in mixing layer with assumed homogeneity in
 * y-direction (2D) or yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleScalarDissipationRateMeanWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const
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
        
        const std::vector<std::vector<double> >& chi_avg_realizations =
            d_ensemble_statistics->chi_avg_realizations;
        
        const std::vector<std::vector<double> >& Y_0_avg_realizations =
            d_ensemble_statistics->Y_0_avg_realizations;
        
        const int num_realizations = static_cast<int>(chi_avg_realizations.size());
        
        TBOX_ASSERT(d_ensemble_statistics->getNumberOfEnsembles() == num_realizations);
        TBOX_ASSERT(num_realizations > 0);
        TBOX_ASSERT(num_realizations == static_cast<int>(Y_0_avg_realizations.size()));
        
        const int num_cells = static_cast<int>(chi_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> chi_avg_global(num_cells, double(0));
        std::vector<double> Y_0_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                chi_avg_global[i] += weight*chi_avg_realizations[ri][i];
                Y_0_avg_global[i] += weight*Y_0_avg_realizations[ri][i];
            }
        }
        
        double chi_sum = double(0);
        int count = 0;
        
        for (int i = 0; i < num_cells; i++)
        {
            const double mixing_metric = double(4)*Y_0_avg_global[i]*(double(1) - Y_0_avg_global[i]);
            if (mixing_metric > double(9)/double(10))
            {
                chi_sum += chi_avg_global[i];
                count++;
            }
        }
        
        const double chi_mean = chi_sum/count;
        
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << chi_mean;
        
        f_out.close();
    }
}


/*
 * Output turbulent Reynolds number based on mixing width with assumed homogeneity in y-direction (2D)
 * or yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleMixingWidthReynoldsNumberWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const
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
 * Output turbulent Mach number with assumed homogeneity in y-direction (2D) or yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleTurbulentMachNumberWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const
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
 * Output effective Atwood number with assumed homogeneity in y-direction (2D) or yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleEffectiveAtwoodNumberWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const
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
        
        const std::vector<std::vector<double> >& rho_sq_avg_realizations =
            d_ensemble_statistics->rho_sq_avg_realizations;
        
        const std::vector<std::vector<double> >& Y_0_avg_realizations =
            d_ensemble_statistics->Y_0_avg_realizations;
        
        const int num_realizations = d_ensemble_statistics->getNumberOfEnsembles();
        
        TBOX_ASSERT(num_realizations > 0);
        TBOX_ASSERT(num_realizations == static_cast<int>(rho_avg_realizations.size()));
        TBOX_ASSERT(num_realizations == static_cast<int>(rho_sq_avg_realizations.size()));
        TBOX_ASSERT(num_realizations == static_cast<int>(Y_0_avg_realizations.size()));
        
        const int num_cells = static_cast<int>(rho_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> rho_avg_global(num_cells, double(0));
        std::vector<double> rho_sq_avg_global(num_cells, double(0));
        std::vector<double> Y_0_avg_global(num_cells, double(0));
        
        std::vector<double> rho_p_rho_p(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                rho_avg_global[i]    += weight*rho_avg_realizations[ri][i];
                rho_sq_avg_global[i] += weight*rho_sq_avg_realizations[ri][i];
                Y_0_avg_global[i]    += weight*Y_0_avg_realizations[ri][i];
            }
        }
        
        for (int i = 0; i < num_cells; i++)
        {
            rho_p_rho_p[i] = rho_sq_avg_global[i] - rho_avg_global[i]*rho_avg_global[i];
        }
        
        double At_sum = double(0);
        int count = 0;
        
        for (int i = 0; i < num_cells; i++)
        {
            const double mixing_metric = double(4)*Y_0_avg_global[i]*(double(1) - Y_0_avg_global[i]);
            if (mixing_metric > double(9)/double(10))
            {
                At_sum += sqrt(rho_p_rho_p[i])/rho_avg_global[i];
                count++;
            }
        }
        
        const double At_mean = At_sum/count;
        
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << At_mean;
        
        f_out.close();
    }
}


/*
 * Output ensemble RMS of baroclinic torque in z-direction integrated to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleRMSBaroclinicTorqueZIntegrated(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
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
        
        f_out.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        const std::vector<std::vector<double> >& B3_sq_avg_realizations =
            d_ensemble_statistics->B3_sq_avg_realizations;
        
        const int num_realizations = d_ensemble_statistics->getNumberOfEnsembles();
        
        TBOX_ASSERT(num_realizations > 0);
        TBOX_ASSERT(num_realizations == static_cast<int>(B3_sq_avg_realizations.size()));
        
        const int num_cells = static_cast<int>(B3_sq_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> B3_sq_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                B3_sq_avg_global[i] += weight*B3_sq_avg_realizations[ri][i];
            }
        }
        
        double B3_rms_integrated_global = double(0);
        for (int i = 0; i < num_cells; i++)
        {
            B3_rms_integrated_global += sqrt(B3_sq_avg_global[i])*dx_finest[0];
        }
        
        const double* x_lo = d_grid_geometry->getXLower();
        const double* x_hi = d_grid_geometry->getXUpper();
        if (d_dim == tbox::Dimension(2))
        {
            const double L_y = x_hi[1] - x_lo[1];
            B3_rms_integrated_global *= L_y;
        }
        else if (d_dim == tbox::Dimension(3))
        {
            const double L_y = x_hi[1] - x_lo[1];
            const double L_z = x_hi[2] - x_lo[2];
            B3_rms_integrated_global *= (L_y*L_z);
        }
        
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << B3_rms_integrated_global;
        
        f_out.close();
    }
}


/*
 * Output mean of ensemble RMS of baroclinic torque in z-direction in mixing layer with assumed homogeneity in
 * y-direction (2D) or yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleRMSBaroclinicTorqueZMeanWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const
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
        
        f_out.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        const std::vector<std::vector<double> >& B3_sq_avg_realizations =
            d_ensemble_statistics->B3_sq_avg_realizations;
        
        const std::vector<std::vector<double> >& Y_0_avg_realizations =
            d_ensemble_statistics->Y_0_avg_realizations;
        
        const int num_realizations = d_ensemble_statistics->getNumberOfEnsembles();
        
        TBOX_ASSERT(num_realizations > 0);
        TBOX_ASSERT(num_realizations == static_cast<int>(B3_sq_avg_realizations.size()));
        TBOX_ASSERT(num_realizations == static_cast<int>(Y_0_avg_realizations.size()));
        
        const int num_cells = static_cast<int>(B3_sq_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> B3_sq_avg_global(num_cells, double(0));
        std::vector<double> Y_0_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                B3_sq_avg_global[i] += weight*B3_sq_avg_realizations[ri][i];
                Y_0_avg_global[i]   += weight*Y_0_avg_realizations[ri][i];
            }
        }
        
        double B3_rms_sum = double(0);
        int count = 0;
        for (int i = 0; i < num_cells; i++)
        {
            const double mixing_metric = double(4)*Y_0_avg_global[i]*(double(1) - Y_0_avg_global[i]);
            if (mixing_metric > double(9)/double(10))
            {
                B3_rms_sum += sqrt(B3_sq_avg_global[i]);
                count++;
            }
        }
        
        const double B3_rms_mean = B3_rms_sum/count;
        
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << B3_rms_mean;
        
        f_out.close();
    }
}


/*
 * Output ensemble RMS of baroclinic torque in x-direction integrated to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleRMSBaroclinicTorqueXIntegrated(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
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
        
        f_out.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        const std::vector<std::vector<double> >& B1_sq_avg_realizations =
            d_ensemble_statistics->B1_sq_avg_realizations;
        
        const int num_realizations = d_ensemble_statistics->getNumberOfEnsembles();
        
        TBOX_ASSERT(num_realizations > 0);
        TBOX_ASSERT(num_realizations == static_cast<int>(B1_sq_avg_realizations.size()));
        
        const int num_cells = static_cast<int>(B1_sq_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> B1_sq_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                B1_sq_avg_global[i] += weight*B1_sq_avg_realizations[ri][i];
            }
        }
        
        double B1_rms_integrated_global = double(0);
        for (int i = 0; i < num_cells; i++)
        {
            B1_rms_integrated_global += sqrt(B1_sq_avg_global[i])*dx_finest[0];
        }
        
        const double* x_lo = d_grid_geometry->getXLower();
        const double* x_hi = d_grid_geometry->getXUpper();
        if (d_dim == tbox::Dimension(2))
        {
            const double L_y = x_hi[1] - x_lo[1];
            B1_rms_integrated_global *= L_y;
        }
        else if (d_dim == tbox::Dimension(3))
        {
            const double L_y = x_hi[1] - x_lo[1];
            const double L_z = x_hi[2] - x_lo[2];
            B1_rms_integrated_global *= (L_y*L_z);
        }
        
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << B1_rms_integrated_global;
        
        f_out.close();
    }
}


/*
 * Output mean of ensemble RMS of baroclinic torque in x-direction in mixing layer with assumed homogeneity in
 * yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleRMSBaroclinicTorqueXMeanWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const
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
        
        f_out.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        const std::vector<std::vector<double> >& B1_sq_avg_realizations =
            d_ensemble_statistics->B1_sq_avg_realizations;
        
        const std::vector<std::vector<double> >& Y_0_avg_realizations =
            d_ensemble_statistics->Y_0_avg_realizations;
        
        const int num_realizations = d_ensemble_statistics->getNumberOfEnsembles();
        
        TBOX_ASSERT(num_realizations > 0);
        TBOX_ASSERT(num_realizations == static_cast<int>(B1_sq_avg_realizations.size()));
        TBOX_ASSERT(num_realizations == static_cast<int>(Y_0_avg_realizations.size()));
        
        const int num_cells = static_cast<int>(B1_sq_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> B1_sq_avg_global(num_cells, double(0));
        std::vector<double> Y_0_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                B1_sq_avg_global[i] += weight*B1_sq_avg_realizations[ri][i];
                Y_0_avg_global[i]   += weight*Y_0_avg_realizations[ri][i];
            }
        }
        
        double B1_rms_sum = double(0);
        int count = 0;
        for (int i = 0; i < num_cells; i++)
        {
            const double mixing_metric = double(4)*Y_0_avg_global[i]*(double(1) - Y_0_avg_global[i]);
            if (mixing_metric > double(9)/double(10))
            {
                B1_rms_sum += sqrt(B1_sq_avg_global[i]);
                count++;
            }
        }
        
        const double B1_rms_mean = B1_rms_sum/count;
        
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << B1_rms_mean;
        
        f_out.close();
    }
}


/*
 * Output ensemble RMS of baroclinic torque in y-direction integrated to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleRMSBaroclinicTorqueYIntegrated(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
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
        
        f_out.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        const std::vector<std::vector<double> >& B2_sq_avg_realizations =
            d_ensemble_statistics->B2_sq_avg_realizations;
        
        const int num_realizations = d_ensemble_statistics->getNumberOfEnsembles();
        
        TBOX_ASSERT(num_realizations > 0);
        TBOX_ASSERT(num_realizations == static_cast<int>(B2_sq_avg_realizations.size()));
        
        const int num_cells = static_cast<int>(B2_sq_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> B2_sq_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                B2_sq_avg_global[i] += weight*B2_sq_avg_realizations[ri][i];
            }
        }
        
        double B2_rms_integrated_global = double(0);
        for (int i = 0; i < num_cells; i++)
        {
            B2_rms_integrated_global += sqrt(B2_sq_avg_global[i])*dx_finest[0];
        }
        
        const double* x_lo = d_grid_geometry->getXLower();
        const double* x_hi = d_grid_geometry->getXUpper();
        if (d_dim == tbox::Dimension(2))
        {
            const double L_y = x_hi[1] - x_lo[1];
            B2_rms_integrated_global *= L_y;
        }
        else if (d_dim == tbox::Dimension(3))
        {
            const double L_y = x_hi[1] - x_lo[1];
            const double L_z = x_hi[2] - x_lo[2];
            B2_rms_integrated_global *= (L_y*L_z);
        }
        
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << B2_rms_integrated_global;
        
        f_out.close();
    }
}


/*
 * Output mean of ensemble RMS of baroclinic torque in y-direction in mixing layer with assumed homogeneity in
 * yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleRMSBaroclinicTorqueYMeanWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const
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
        
        f_out.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        const std::vector<std::vector<double> >& B2_sq_avg_realizations =
            d_ensemble_statistics->B2_sq_avg_realizations;
        
        const std::vector<std::vector<double> >& Y_0_avg_realizations =
            d_ensemble_statistics->Y_0_avg_realizations;
        
        const int num_realizations = d_ensemble_statistics->getNumberOfEnsembles();
        
        TBOX_ASSERT(num_realizations > 0);
        TBOX_ASSERT(num_realizations == static_cast<int>(B2_sq_avg_realizations.size()));
        TBOX_ASSERT(num_realizations == static_cast<int>(Y_0_avg_realizations.size()));
        
        const int num_cells = static_cast<int>(B2_sq_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> B2_sq_avg_global(num_cells, double(0));
        std::vector<double> Y_0_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                B2_sq_avg_global[i] += weight*B2_sq_avg_realizations[ri][i];
                Y_0_avg_global[i]   += weight*Y_0_avg_realizations[ri][i];
            }
        }
        
        double B2_rms_sum = double(0);
        int count = 0;
        for (int i = 0; i < num_cells; i++)
        {
            const double mixing_metric = double(4)*Y_0_avg_global[i]*(double(1) - Y_0_avg_global[i]);
            if (mixing_metric > double(9)/double(10))
            {
                B2_rms_sum += sqrt(B2_sq_avg_global[i]);
                count++;
            }
        }
        
        const double B2_rms_mean = B2_rms_sum/count;
        
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << B2_rms_mean;
        
        f_out.close();
    }
}


/*
 * Output ensemble TKE production term 1 integrated to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleTKEProduction1Integrated(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
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
        
        const std::vector<std::vector<double> >& v_avg_realizations =
            d_ensemble_statistics->v_avg_realizations;
        
        const std::vector<std::vector<double> >& w_avg_realizations =
            d_ensemble_statistics->w_avg_realizations;
        
        const std::vector<std::vector<double> >& rho_u_avg_realizations =
            d_ensemble_statistics->rho_u_avg_realizations;
        
        const std::vector<std::vector<double> >& rho_v_avg_realizations =
            d_ensemble_statistics->rho_v_avg_realizations;
        
        const std::vector<std::vector<double> >& rho_w_avg_realizations =
            d_ensemble_statistics->rho_w_avg_realizations;
        
        const std::vector<std::vector<double> >& ddx_p_avg_realizations =
            d_ensemble_statistics->ddx_p_avg_realizations;
        
        const std::vector<std::vector<double> >& ddy_p_avg_realizations =
            d_ensemble_statistics->ddy_p_avg_realizations;
        
        const std::vector<std::vector<double> >& ddz_p_avg_realizations =
            d_ensemble_statistics->ddz_p_avg_realizations;
        
        const int num_realizations = static_cast<int>(rho_avg_realizations.size());
        
        TBOX_ASSERT(d_ensemble_statistics->getNumberOfEnsembles() == num_realizations);
        TBOX_ASSERT(num_realizations > 0);
        TBOX_ASSERT(num_realizations == static_cast<int>(u_avg_realizations.size()));
        TBOX_ASSERT(num_realizations == static_cast<int>(rho_u_avg_realizations.size()));
        TBOX_ASSERT(num_realizations == static_cast<int>(ddx_p_avg_realizations.size()));
        if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
        {
            TBOX_ASSERT(num_realizations == static_cast<int>(v_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(rho_v_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(ddy_p_avg_realizations.size()));
        }
        if (d_dim == tbox::Dimension(3))
        {
            TBOX_ASSERT(num_realizations == static_cast<int>(w_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(rho_w_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(ddz_p_avg_realizations.size()));
        }
        
        const int num_cells = static_cast<int>(rho_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> rho_avg_global(num_cells, double(0));
        std::vector<double> u_avg_global(num_cells, double(0));
        std::vector<double> v_avg_global(num_cells, double(0));
        std::vector<double> w_avg_global(num_cells, double(0));
        std::vector<double> rho_u_avg_global(num_cells, double(0));
        std::vector<double> rho_v_avg_global(num_cells, double(0));
        std::vector<double> rho_w_avg_global(num_cells, double(0));
        std::vector<double> ddx_p_avg_global(num_cells, double(0));
        std::vector<double> ddy_p_avg_global(num_cells, double(0));
        std::vector<double> ddz_p_avg_global(num_cells, double(0));
        
        std::vector<double> a1(num_cells, double(0));
        std::vector<double> a2(num_cells, double(0));
        std::vector<double> a3(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                rho_avg_global[i]   += weight*rho_avg_realizations[ri][i];
                u_avg_global[i]     += weight*u_avg_realizations[ri][i];
                rho_u_avg_global[i] += weight*rho_u_avg_realizations[ri][i];
                ddx_p_avg_global[i] += weight*ddx_p_avg_realizations[ri][i];
            }
        }
        
        for (int i = 0; i < num_cells; i++)
        {
            a1[i] = (rho_u_avg_global[i] - rho_avg_global[i]*u_avg_global[i])/rho_avg_global[i];
        }
        
        if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
        {
            for (int ri = 0; ri < num_realizations; ri++)
            {
                for (int i = 0; i < num_cells; i++)
                {
                    v_avg_global[i]     += weight*v_avg_realizations[ri][i];
                    rho_v_avg_global[i] += weight*rho_v_avg_realizations[ri][i];
                    ddy_p_avg_global[i] += weight*ddy_p_avg_realizations[ri][i];
                }
            }
            
            for (int i = 0; i < num_cells; i++)
            {
                a2[i] = (rho_v_avg_global[i] - rho_avg_global[i]*v_avg_global[i])/rho_avg_global[i];
            }
        }
        
        if (d_dim == tbox::Dimension(3))
        {
            for (int ri = 0; ri < num_realizations; ri++)
            {
                for (int i = 0; i < num_cells; i++)
                {
                    w_avg_global[i]     += weight*w_avg_realizations[ri][i];
                    rho_w_avg_global[i] += weight*rho_w_avg_realizations[ri][i];
                    ddz_p_avg_global[i] += weight*ddz_p_avg_realizations[ri][i];
                }
            }
            
            for (int i = 0; i < num_cells; i++)
            {
                a3[i] = (rho_w_avg_global[i] - rho_avg_global[i]*w_avg_global[i])/rho_avg_global[i];
            }
        }
        
        double TKE_prod_integrated_global = double(0);
        if (d_dim == tbox::Dimension(1))
        {
            for (int i = 0; i < num_cells; i++)
            {
                TKE_prod_integrated_global += (a1[i]*ddx_p_avg_global[i])*dx_finest[0];
            }
        }
        else if (d_dim == tbox::Dimension(2))
        {
            for (int i = 0; i < num_cells; i++)
            {
                TKE_prod_integrated_global += (a1[i]*ddx_p_avg_global[i] + a2[i]*ddy_p_avg_global[i])*dx_finest[0];
            }
        }
        else if (d_dim == tbox::Dimension(3))
        {
            for (int i = 0; i < num_cells; i++)
            {
                TKE_prod_integrated_global += (a1[i]*ddx_p_avg_global[i] + a2[i]*ddy_p_avg_global[i] + a3[i]*ddz_p_avg_global[i])*dx_finest[0];
            }
        }
        
        const double* x_lo = d_grid_geometry->getXLower();
        const double* x_hi = d_grid_geometry->getXUpper();
        if (d_dim == tbox::Dimension(2))
        {
            const double L_y = x_hi[1] - x_lo[1];
            TKE_prod_integrated_global *= L_y;
        }
        else if (d_dim == tbox::Dimension(3))
        {
            const double L_y = x_hi[1] - x_lo[1];
            const double L_z = x_hi[2] - x_lo[2];
            TKE_prod_integrated_global *= (L_y*L_z);
        }
        
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << TKE_prod_integrated_global;
        
        f_out.close();
    }
}


/*
 * Output mean of ensemble TKE production term 1 in mixing layer with assumed homogeneity in y-direction (2D) or
 * yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleTKEProduction1MeanWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const
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
        
        const std::vector<std::vector<double> >& v_avg_realizations =
            d_ensemble_statistics->v_avg_realizations;
        
        const std::vector<std::vector<double> >& w_avg_realizations =
            d_ensemble_statistics->w_avg_realizations;
        
        const std::vector<std::vector<double> >& rho_u_avg_realizations =
            d_ensemble_statistics->rho_u_avg_realizations;
        
        const std::vector<std::vector<double> >& rho_v_avg_realizations =
            d_ensemble_statistics->rho_v_avg_realizations;
        
        const std::vector<std::vector<double> >& rho_w_avg_realizations =
            d_ensemble_statistics->rho_w_avg_realizations;
        
        const std::vector<std::vector<double> >& ddx_p_avg_realizations =
            d_ensemble_statistics->ddx_p_avg_realizations;
        
        const std::vector<std::vector<double> >& ddy_p_avg_realizations =
            d_ensemble_statistics->ddy_p_avg_realizations;
        
        const std::vector<std::vector<double> >& ddz_p_avg_realizations =
            d_ensemble_statistics->ddz_p_avg_realizations;
        
        const std::vector<std::vector<double> >& Y_0_avg_realizations =
            d_ensemble_statistics->Y_0_avg_realizations;
        
        const int num_realizations = static_cast<int>(rho_avg_realizations.size());
        
        TBOX_ASSERT(d_ensemble_statistics->getNumberOfEnsembles() == num_realizations);
        TBOX_ASSERT(num_realizations > 0);
        TBOX_ASSERT(num_realizations == static_cast<int>(u_avg_realizations.size()));
        TBOX_ASSERT(num_realizations == static_cast<int>(rho_u_avg_realizations.size()));
        TBOX_ASSERT(num_realizations == static_cast<int>(ddx_p_avg_realizations.size()));
        if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
        {
            TBOX_ASSERT(num_realizations == static_cast<int>(v_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(rho_v_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(ddy_p_avg_realizations.size()));
        }
        if (d_dim == tbox::Dimension(3))
        {
            TBOX_ASSERT(num_realizations == static_cast<int>(w_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(rho_w_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(ddz_p_avg_realizations.size()));
        }
        TBOX_ASSERT(num_realizations == static_cast<int>(Y_0_avg_realizations.size()));
        
        const int num_cells = static_cast<int>(rho_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> rho_avg_global(num_cells, double(0));
        std::vector<double> u_avg_global(num_cells, double(0));
        std::vector<double> v_avg_global(num_cells, double(0));
        std::vector<double> w_avg_global(num_cells, double(0));
        std::vector<double> rho_u_avg_global(num_cells, double(0));
        std::vector<double> rho_v_avg_global(num_cells, double(0));
        std::vector<double> rho_w_avg_global(num_cells, double(0));
        std::vector<double> ddx_p_avg_global(num_cells, double(0));
        std::vector<double> ddy_p_avg_global(num_cells, double(0));
        std::vector<double> ddz_p_avg_global(num_cells, double(0));
        std::vector<double> Y_0_avg_global(num_cells, double(0));
        
        std::vector<double> a1(num_cells, double(0));
        std::vector<double> a2(num_cells, double(0));
        std::vector<double> a3(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                rho_avg_global[i]   += weight*rho_avg_realizations[ri][i];
                u_avg_global[i]     += weight*u_avg_realizations[ri][i];
                rho_u_avg_global[i] += weight*rho_u_avg_realizations[ri][i];
                ddx_p_avg_global[i] += weight*ddx_p_avg_realizations[ri][i];
                Y_0_avg_global[i]   += weight*Y_0_avg_realizations[ri][i];
            }
        }
        
        for (int i = 0; i < num_cells; i++)
        {
            a1[i] = (rho_u_avg_global[i] - rho_avg_global[i]*u_avg_global[i])/rho_avg_global[i];
        }
        
        if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
        {
            for (int ri = 0; ri < num_realizations; ri++)
            {
                for (int i = 0; i < num_cells; i++)
                {
                    v_avg_global[i]     += weight*v_avg_realizations[ri][i];
                    rho_v_avg_global[i] += weight*rho_v_avg_realizations[ri][i];
                    ddy_p_avg_global[i] += weight*ddy_p_avg_realizations[ri][i];
                }
            }
            
            for (int i = 0; i < num_cells; i++)
            {
                a2[i] = (rho_v_avg_global[i] - rho_avg_global[i]*v_avg_global[i])/rho_avg_global[i];
            }
        }
        
        if (d_dim == tbox::Dimension(3))
        {
            for (int ri = 0; ri < num_realizations; ri++)
            {
                for (int i = 0; i < num_cells; i++)
                {
                    w_avg_global[i]     += weight*w_avg_realizations[ri][i];
                    rho_w_avg_global[i] += weight*rho_w_avg_realizations[ri][i];
                    ddz_p_avg_global[i] += weight*ddz_p_avg_realizations[ri][i];
                }
            }
            
            for (int i = 0; i < num_cells; i++)
            {
                a3[i] = (rho_w_avg_global[i] - rho_avg_global[i]*w_avg_global[i])/rho_avg_global[i];
            }
        }
        
        double TKE_prod_sum = double(0);
        int count = 0;
        
        if (d_dim == tbox::Dimension(1))
        {
            for (int i = 0; i < num_cells; i++)
            {
                const double mixing_metric = double(4)*Y_0_avg_global[i]*(double(1) - Y_0_avg_global[i]);
                if (mixing_metric > double(9)/double(10))
                {
                    TKE_prod_sum += (a1[i]*ddx_p_avg_global[i]);
                    count++;
                }
            }
        }
        else if (d_dim == tbox::Dimension(2))
        {
            for (int i = 0; i < num_cells; i++)
            {
                const double mixing_metric = double(4)*Y_0_avg_global[i]*(double(1) - Y_0_avg_global[i]);
                if (mixing_metric > double(9)/double(10))
                {
                    TKE_prod_sum += (a1[i]*ddx_p_avg_global[i] + a2[i]*ddy_p_avg_global[i]);
                    count++;
                }
            }
        }
        else if (d_dim == tbox::Dimension(3))
        {
            for (int i = 0; i < num_cells; i++)
            {
                const double mixing_metric = double(4)*Y_0_avg_global[i]*(double(1) - Y_0_avg_global[i]);
                if (mixing_metric > double(9)/double(10))
                {
                    TKE_prod_sum += (a1[i]*ddx_p_avg_global[i] + a2[i]*ddy_p_avg_global[i] + a3[i]*ddz_p_avg_global[i]);
                    count++;
                }
            }
        }
        
        const double TKE_prod_mean = TKE_prod_sum/count;
        
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << TKE_prod_mean;
        
        f_out.close();
    }
}


/*
 * Output ensemble TKE production term 2 integrated to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleTKEProduction2Integrated(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
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
        
        const std::vector<std::vector<double> >& ddx_tau11_avg_realizations =
            d_ensemble_statistics->ddx_tau11_avg_realizations;
        
        const int num_realizations = d_ensemble_statistics->getNumberOfEnsembles();
        
        TBOX_ASSERT(num_realizations > 0);
        TBOX_ASSERT(num_realizations == static_cast<int>(rho_avg_realizations.size()));
        TBOX_ASSERT(num_realizations == static_cast<int>(u_avg_realizations.size()));
        TBOX_ASSERT(num_realizations == static_cast<int>(rho_u_avg_realizations.size()));
        TBOX_ASSERT(num_realizations == static_cast<int>(ddx_tau11_avg_realizations.size()));
        
        const int num_cells = static_cast<int>(rho_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> rho_avg_global(num_cells, double(0));
        std::vector<double> u_avg_global(num_cells, double(0));
        std::vector<double> rho_u_avg_global(num_cells, double(0));
        std::vector<double> ddx_tau11_avg_global(num_cells, double(0));
        
        std::vector<double> a1(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                rho_avg_global[i]       += weight*rho_avg_realizations[ri][i];
                u_avg_global[i]         += weight*u_avg_realizations[ri][i];
                rho_u_avg_global[i]     += weight*rho_u_avg_realizations[ri][i];
                ddx_tau11_avg_global[i] += weight*ddx_tau11_avg_realizations[ri][i];
            }
        }
        
        for (int i = 0; i < num_cells; i++)
        {
            a1[i] = (rho_u_avg_global[i] - rho_avg_global[i]*u_avg_global[i])/rho_avg_global[i];
        }
        
        double TKE_prod_integrated_global = double(0);
        for (int i = 0; i < num_cells; i++)
        {
            TKE_prod_integrated_global += (a1[i]*ddx_tau11_avg_global[i])*dx_finest[0];
        }
        
        const double* x_lo = d_grid_geometry->getXLower();
        const double* x_hi = d_grid_geometry->getXUpper();
        if (d_dim == tbox::Dimension(2))
        {
            const double L_y = x_hi[1] - x_lo[1];
            TKE_prod_integrated_global *= L_y;
        }
        else if (d_dim == tbox::Dimension(3))
        {
            const double L_y = x_hi[1] - x_lo[1];
            const double L_z = x_hi[2] - x_lo[2];
            TKE_prod_integrated_global *= (L_y*L_z);
        }
        
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << TKE_prod_integrated_global;
        
        f_out.close();
    }
}


/*
 * Output mean of ensemble TKE production term 2 in mixing layer with assumed homogeneity in y-direction (2D) or
 * yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleTKEProduction2MeanWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Output the spatial profile (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_out;
        
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
        
        const std::vector<std::vector<double> >& ddx_tau11_avg_realizations =
            d_ensemble_statistics->ddx_tau11_avg_realizations;
        
        const std::vector<std::vector<double> >& Y_0_avg_realizations =
            d_ensemble_statistics->Y_0_avg_realizations;
        
        const int num_realizations = d_ensemble_statistics->getNumberOfEnsembles();
        
        TBOX_ASSERT(num_realizations > 0);
        TBOX_ASSERT(num_realizations == static_cast<int>(rho_avg_realizations.size()));
        TBOX_ASSERT(num_realizations == static_cast<int>(u_avg_realizations.size()));
        TBOX_ASSERT(num_realizations == static_cast<int>(rho_u_avg_realizations.size()));
        TBOX_ASSERT(num_realizations == static_cast<int>(ddx_tau11_avg_realizations.size()));
        TBOX_ASSERT(num_realizations == static_cast<int>(Y_0_avg_realizations.size()));
        
        const int num_cells = static_cast<int>(rho_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> rho_avg_global(num_cells, double(0));
        std::vector<double> u_avg_global(num_cells, double(0));
        std::vector<double> rho_u_avg_global(num_cells, double(0));
        std::vector<double> ddx_tau11_avg_global(num_cells, double(0));
        std::vector<double> Y_0_avg_global(num_cells, double(0));
        
        std::vector<double> a1(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                rho_avg_global[i]       += weight*rho_avg_realizations[ri][i];
                u_avg_global[i]         += weight*u_avg_realizations[ri][i];
                rho_u_avg_global[i]     += weight*rho_u_avg_realizations[ri][i];
                ddx_tau11_avg_global[i] += weight*ddx_tau11_avg_realizations[ri][i];
                Y_0_avg_global[i]       += weight*Y_0_avg_realizations[ri][i];
            }
        }
        
        for (int i = 0; i < num_cells; i++)
        {
            a1[i] = (rho_u_avg_global[i] - rho_avg_global[i]*u_avg_global[i])/rho_avg_global[i];
        }
        
        double TKE_prod_sum = double(0);
        int count = 0;
        
        for (int i = 0; i < num_cells; i++)
        {
            const double mixing_metric = double(4)*Y_0_avg_global[i]*(double(1) - Y_0_avg_global[i]);
            if (mixing_metric > double(9)/double(10))
            {
                TKE_prod_sum += (a1[i]*ddx_tau11_avg_global[i]);
                count++;
            }
        }
        
        const double TKE_prod_mean = TKE_prod_sum/count;
        
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << TKE_prod_mean;
        
        f_out.close();
    }
}


/*
 * Output ensemble TKE production term 3 integrated to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleTKEProduction3Integrated(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
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
            d_ensemble_statistics->rho_u_avg_realizations;
        
        const int num_realizations = d_ensemble_statistics->getNumberOfEnsembles();
        
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
        
        std::vector<double> u_tilde(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            u_tilde[i] = rho_u_avg_global[i]/rho_avg_global[i];
        }
        
        std::vector<double> ddx_u_tilde = computeDerivativeOfVector1D(
            u_tilde,
            dx_finest[0]);
        
        std::vector<double> rho_R11(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            rho_R11[i] = rho_u_u_avg_global[i] - rho_u_avg_global[i]*u_tilde[i];
        }
        
        double TKE_prod_integrated_global = double(0);
        for (int i = 0; i < num_cells; i++)
        {
            TKE_prod_integrated_global -= rho_R11[i]*ddx_u_tilde[i];
        }
        
        const double* x_lo = d_grid_geometry->getXLower();
        const double* x_hi = d_grid_geometry->getXUpper();
        if (d_dim == tbox::Dimension(2))
        {
            const double L_y = x_hi[1] - x_lo[1];
            TKE_prod_integrated_global *= L_y;
        }
        else if (d_dim == tbox::Dimension(3))
        {
            const double L_y = x_hi[1] - x_lo[1];
            const double L_z = x_hi[2] - x_lo[2];
            TKE_prod_integrated_global *= (L_y*L_z);
        }
        
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << TKE_prod_integrated_global;
        
        f_out.close();
    }
}


/*
 * Output mean of ensemble TKE production term 3 in mixing layer with assumed homogeneity in y-direction (2D) or
 * yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleTKEProduction3MeanWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
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
            d_ensemble_statistics->rho_u_avg_realizations;
        
        const std::vector<std::vector<double> >& Y_0_avg_realizations =
            d_ensemble_statistics->Y_0_avg_realizations;
        
        const int num_realizations = d_ensemble_statistics->getNumberOfEnsembles();
        
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
        
        std::vector<double> u_tilde(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            u_tilde[i] = rho_u_avg_global[i]/rho_avg_global[i];
        }
        
        std::vector<double> ddx_u_tilde = computeDerivativeOfVector1D(
            u_tilde,
            dx_finest[0]);
        
        std::vector<double> rho_R11(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            rho_R11[i] = rho_u_u_avg_global[i] - rho_u_avg_global[i]*u_tilde[i];
        }
        
        double TKE_prod_sum = double(0);
        int count = 0;
        
        for (int i = 0; i < num_cells; i++)
        {
            const double mixing_metric = double(4)*Y_0_avg_global[i]*(double(1) - Y_0_avg_global[i]);
            if (mixing_metric > double(9)/double(10))
            {
                TKE_prod_sum -= (rho_R11[i]*ddx_u_tilde[i]);
                count++;
            }
        }
        
        const double TKE_prod_mean = TKE_prod_sum/count;
        
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << TKE_prod_mean;
        
        f_out.close();
    }
}


/*
 * Output ensemble TKE dissipation integrated to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleTKEDissipationIntegrated(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
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
        
        f_out.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        const std::vector<std::vector<double> >& ddx_u_avg_realizations =
            d_ensemble_statistics->ddx_u_avg_realizations;
        
        const std::vector<std::vector<double> >& ddy_u_avg_realizations =
            d_ensemble_statistics->ddy_u_avg_realizations;
        
        const std::vector<std::vector<double> >& ddz_u_avg_realizations =
            d_ensemble_statistics->ddz_u_avg_realizations;
        
        const std::vector<std::vector<double> >& ddx_v_avg_realizations =
            d_ensemble_statistics->ddx_v_avg_realizations;
        
        const std::vector<std::vector<double> >& ddy_v_avg_realizations =
            d_ensemble_statistics->ddy_v_avg_realizations;
        
        const std::vector<std::vector<double> >& ddz_v_avg_realizations =
            d_ensemble_statistics->ddz_v_avg_realizations;
        
        const std::vector<std::vector<double> >& ddx_w_avg_realizations =
            d_ensemble_statistics->ddx_w_avg_realizations;
        
        const std::vector<std::vector<double> >& ddy_w_avg_realizations =
            d_ensemble_statistics->ddy_w_avg_realizations;
        
        const std::vector<std::vector<double> >& ddz_w_avg_realizations =
            d_ensemble_statistics->ddz_w_avg_realizations;
        
        const std::vector<std::vector<double> >& tau11_avg_realizations =
            d_ensemble_statistics->tau11_avg_realizations;
        
        const std::vector<std::vector<double> >& tau12_avg_realizations =
            d_ensemble_statistics->tau12_avg_realizations;
        
        const std::vector<std::vector<double> >& tau13_avg_realizations =
            d_ensemble_statistics->tau13_avg_realizations;
        
        const std::vector<std::vector<double> >& tau22_avg_realizations =
            d_ensemble_statistics->tau22_avg_realizations;
        
        const std::vector<std::vector<double> >& tau23_avg_realizations =
            d_ensemble_statistics->tau23_avg_realizations;
        
        const std::vector<std::vector<double> >& tau33_avg_realizations =
            d_ensemble_statistics->tau33_avg_realizations;
        
        const std::vector<std::vector<double> >& tau11_ddx_u_avg_realizations =
            d_ensemble_statistics->tau11_ddx_u_avg_realizations;
        
        const std::vector<std::vector<double> >& tau12_ddy_u_avg_realizations =
            d_ensemble_statistics->tau12_ddy_u_avg_realizations;
        
        const std::vector<std::vector<double> >& tau13_ddz_u_avg_realizations =
            d_ensemble_statistics->tau13_ddz_u_avg_realizations;
        
        const std::vector<std::vector<double> >& tau12_ddx_v_avg_realizations =
            d_ensemble_statistics->tau12_ddx_v_avg_realizations;
        
        const std::vector<std::vector<double> >& tau22_ddy_v_avg_realizations =
            d_ensemble_statistics->tau22_ddy_v_avg_realizations;
        
        const std::vector<std::vector<double> >& tau23_ddz_v_avg_realizations =
            d_ensemble_statistics->tau23_ddz_v_avg_realizations;
        
        const std::vector<std::vector<double> >& tau13_ddx_w_avg_realizations =
            d_ensemble_statistics->tau13_ddx_w_avg_realizations;
        
        const std::vector<std::vector<double> >& tau23_ddy_w_avg_realizations =
            d_ensemble_statistics->tau23_ddy_w_avg_realizations;
        
        const std::vector<std::vector<double> >& tau33_ddz_w_avg_realizations =
            d_ensemble_statistics->tau33_ddz_w_avg_realizations;
        
        const int num_realizations = d_ensemble_statistics->getNumberOfEnsembles();
        
        TBOX_ASSERT(num_realizations > 0);
        TBOX_ASSERT(num_realizations == static_cast<int>(ddx_u_avg_realizations.size()));
        TBOX_ASSERT(num_realizations == static_cast<int>(tau11_avg_realizations.size()));
        TBOX_ASSERT(num_realizations == static_cast<int>(tau11_ddx_u_avg_realizations.size()));
        
        if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
        {
            TBOX_ASSERT(num_realizations == static_cast<int>(ddy_u_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(ddx_v_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(ddy_v_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(tau12_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(tau22_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(tau12_ddy_u_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(tau12_ddx_v_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(tau22_ddy_v_avg_realizations.size()));
        }
        if (d_dim == tbox::Dimension(3))
        {
            TBOX_ASSERT(num_realizations == static_cast<int>(ddz_u_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(ddz_v_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(ddx_w_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(ddy_w_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(ddz_w_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(tau13_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(tau23_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(tau33_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(tau13_ddz_u_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(tau23_ddz_v_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(tau13_ddx_w_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(tau23_ddy_w_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(tau33_ddz_w_avg_realizations.size()));
        }
        
        const int num_cells = static_cast<int>(ddx_u_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> ddx_u_avg_global(num_cells, double(0));
        std::vector<double> ddy_u_avg_global(num_cells, double(0));
        std::vector<double> ddz_u_avg_global(num_cells, double(0));
        std::vector<double> ddx_v_avg_global(num_cells, double(0));
        std::vector<double> ddy_v_avg_global(num_cells, double(0));
        std::vector<double> ddz_v_avg_global(num_cells, double(0));
        std::vector<double> ddx_w_avg_global(num_cells, double(0));
        std::vector<double> ddy_w_avg_global(num_cells, double(0));
        std::vector<double> ddz_w_avg_global(num_cells, double(0));
        std::vector<double> tau11_avg_global(num_cells, double(0));
        std::vector<double> tau12_avg_global(num_cells, double(0));
        std::vector<double> tau13_avg_global(num_cells, double(0));
        std::vector<double> tau22_avg_global(num_cells, double(0));
        std::vector<double> tau23_avg_global(num_cells, double(0));
        std::vector<double> tau33_avg_global(num_cells, double(0));
        std::vector<double> tau11_ddx_u_avg_global(num_cells, double(0));
        std::vector<double> tau12_ddy_u_avg_global(num_cells, double(0));
        std::vector<double> tau13_ddz_u_avg_global(num_cells, double(0));
        std::vector<double> tau12_ddx_v_avg_global(num_cells, double(0));
        std::vector<double> tau22_ddy_v_avg_global(num_cells, double(0));
        std::vector<double> tau23_ddz_v_avg_global(num_cells, double(0));
        std::vector<double> tau13_ddx_w_avg_global(num_cells, double(0));
        std::vector<double> tau23_ddy_w_avg_global(num_cells, double(0));
        std::vector<double> tau33_ddz_w_avg_global(num_cells, double(0));
        
        std::vector<double> tau11_p_ddx_u_p(num_cells, double(0));
        std::vector<double> tau12_p_ddy_u_p(num_cells, double(0));
        std::vector<double> tau13_p_ddz_u_p(num_cells, double(0));
        std::vector<double> tau12_p_ddx_v_p(num_cells, double(0));
        std::vector<double> tau22_p_ddy_v_p(num_cells, double(0));
        std::vector<double> tau23_p_ddz_v_p(num_cells, double(0));
        std::vector<double> tau13_p_ddx_w_p(num_cells, double(0));
        std::vector<double> tau23_p_ddy_w_p(num_cells, double(0));
        std::vector<double> tau33_p_ddz_w_p(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                ddx_u_avg_global[i] += weight*ddx_u_avg_realizations[ri][i];
                tau11_avg_global[i] += weight*tau11_avg_realizations[ri][i];
                tau11_ddx_u_avg_global[i] += weight*tau11_ddx_u_avg_realizations[ri][i];
            }
        }
        
        for (int i = 0; i < num_cells; i++)
        {
            tau11_p_ddx_u_p[i] = tau11_ddx_u_avg_global[i] - tau11_avg_global[i]*ddx_u_avg_global[i];
        }
        
        if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
        {
            for (int ri = 0; ri < num_realizations; ri++)
            {
                for (int i = 0; i < num_cells; i++)
                {
                    ddy_u_avg_global[i] += weight*ddy_u_avg_realizations[ri][i];
                    ddx_v_avg_global[i] += weight*ddx_v_avg_realizations[ri][i];
                    ddy_v_avg_global[i] += weight*ddy_v_avg_realizations[ri][i];
                    tau12_avg_global[i] += weight*tau12_avg_realizations[ri][i];
                    tau22_avg_global[i] += weight*tau22_avg_realizations[ri][i];
                    tau12_ddy_u_avg_global[i] += weight*tau12_ddy_u_avg_realizations[ri][i];
                    tau12_ddx_v_avg_global[i] += weight*tau12_ddx_v_avg_realizations[ri][i];
                    tau22_ddy_v_avg_global[i] += weight*tau22_ddy_v_avg_realizations[ri][i];
                }
            }
            
            for (int i = 0; i < num_cells; i++)
            {
                tau12_p_ddy_u_p[i] = tau12_ddy_u_avg_global[i] - tau12_avg_global[i]*ddy_u_avg_global[i];
                tau12_p_ddx_v_p[i] = tau12_ddx_v_avg_global[i] - tau12_avg_global[i]*ddx_v_avg_global[i];
                tau22_p_ddy_v_p[i] = tau22_ddy_v_avg_global[i] - tau22_avg_global[i]*ddy_v_avg_global[i];
            }
        }
        
        if (d_dim == tbox::Dimension(3))
        {
            for (int ri = 0; ri < num_realizations; ri++)
            {
                for (int i = 0; i < num_cells; i++)
                {
                    ddz_u_avg_global[i] += weight*ddz_u_avg_realizations[ri][i];
                    ddz_v_avg_global[i] += weight*ddz_v_avg_realizations[ri][i];
                    ddx_w_avg_global[i] += weight*ddx_w_avg_realizations[ri][i];
                    ddy_w_avg_global[i] += weight*ddy_w_avg_realizations[ri][i];
                    ddz_w_avg_global[i] += weight*ddz_w_avg_realizations[ri][i];
                    tau13_avg_global[i] += weight*tau13_avg_realizations[ri][i];
                    tau23_avg_global[i] += weight*tau23_avg_realizations[ri][i];
                    tau33_avg_global[i] += weight*tau33_avg_realizations[ri][i];
                    tau13_ddz_u_avg_global[i] += weight*tau13_ddz_u_avg_realizations[ri][i];
                    tau23_ddz_v_avg_global[i] += weight*tau23_ddz_v_avg_realizations[ri][i];
                    tau13_ddx_w_avg_global[i] += weight*tau13_ddx_w_avg_realizations[ri][i];
                    tau23_ddy_w_avg_global[i] += weight*tau23_ddy_w_avg_realizations[ri][i];
                    tau33_ddz_w_avg_global[i] += weight*tau33_ddz_w_avg_realizations[ri][i];
                }
            }
            
            for (int i = 0; i < num_cells; i++)
            {
                tau13_p_ddz_u_p[i] = tau13_ddz_u_avg_global[i] - tau13_avg_global[i]*ddz_u_avg_global[i];
                tau23_p_ddz_v_p[i] = tau23_ddz_v_avg_global[i] - tau23_avg_global[i]*ddz_v_avg_global[i];
                tau13_p_ddx_w_p[i] = tau13_ddx_w_avg_global[i] - tau13_avg_global[i]*ddx_w_avg_global[i];
                tau23_p_ddy_w_p[i] = tau23_ddy_w_avg_global[i] - tau23_avg_global[i]*ddy_w_avg_global[i];
                tau33_p_ddz_w_p[i] = tau33_ddz_w_avg_global[i] - tau33_avg_global[i]*ddz_w_avg_global[i];
            }
        }
        
        double TKE_diss_integrated_global = double(0);
        if (d_dim == tbox::Dimension(1))
        {
            for (int i = 0; i < num_cells; i++)
            {
                TKE_diss_integrated_global += (tau11_p_ddx_u_p[i])*dx_finest[0];
            }
        }
        else if (d_dim == tbox::Dimension(2))
        {
            for (int i = 0; i < num_cells; i++)
            {
                TKE_diss_integrated_global += (tau11_p_ddx_u_p[i] + tau22_p_ddy_v_p[i] +
                    tau12_p_ddx_v_p[i] + tau12_p_ddy_u_p[i]
                    )*dx_finest[0];
            }
        }
        else if (d_dim == tbox::Dimension(3))
        {
            for (int i = 0; i < num_cells; i++)
            {
                TKE_diss_integrated_global += (tau11_p_ddx_u_p[i] + tau22_p_ddy_v_p[i] + + tau33_p_ddz_w_p[i] +
                    tau12_p_ddx_v_p[i] + tau12_p_ddy_u_p[i] +
                    tau13_p_ddx_w_p[i] + tau13_p_ddz_u_p[i] +
                    tau23_p_ddy_w_p[i] + tau23_p_ddz_v_p[i]
                    )*dx_finest[0];
            }
        }
        
        const double* x_lo = d_grid_geometry->getXLower();
        const double* x_hi = d_grid_geometry->getXUpper();
        if (d_dim == tbox::Dimension(2))
        {
            const double L_y = x_hi[1] - x_lo[1];
            TKE_diss_integrated_global *= L_y;
        }
        else if (d_dim == tbox::Dimension(3))
        {
            const double L_y = x_hi[1] - x_lo[1];
            const double L_z = x_hi[2] - x_lo[2];
            TKE_diss_integrated_global *= (L_y*L_z);
        }
        
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << TKE_diss_integrated_global;
        
        f_out.close();
    }
}


/*
 * Output mean of ensemble TKE dissipation in mixing layer with assumed homogeneity in y-direction (2D) or
 * yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleTKEDissipationMeanWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const
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
        
        f_out.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        const std::vector<std::vector<double> >& ddx_u_avg_realizations =
            d_ensemble_statistics->ddx_u_avg_realizations;
        
        const std::vector<std::vector<double> >& ddy_u_avg_realizations =
            d_ensemble_statistics->ddy_u_avg_realizations;
        
        const std::vector<std::vector<double> >& ddz_u_avg_realizations =
            d_ensemble_statistics->ddz_u_avg_realizations;
        
        const std::vector<std::vector<double> >& ddx_v_avg_realizations =
            d_ensemble_statistics->ddx_v_avg_realizations;
        
        const std::vector<std::vector<double> >& ddy_v_avg_realizations =
            d_ensemble_statistics->ddy_v_avg_realizations;
        
        const std::vector<std::vector<double> >& ddz_v_avg_realizations =
            d_ensemble_statistics->ddz_v_avg_realizations;
        
        const std::vector<std::vector<double> >& ddx_w_avg_realizations =
            d_ensemble_statistics->ddx_w_avg_realizations;
        
        const std::vector<std::vector<double> >& ddy_w_avg_realizations =
            d_ensemble_statistics->ddy_w_avg_realizations;
        
        const std::vector<std::vector<double> >& ddz_w_avg_realizations =
            d_ensemble_statistics->ddz_w_avg_realizations;
        
        const std::vector<std::vector<double> >& tau11_avg_realizations =
            d_ensemble_statistics->tau11_avg_realizations;
        
        const std::vector<std::vector<double> >& tau12_avg_realizations =
            d_ensemble_statistics->tau12_avg_realizations;
        
        const std::vector<std::vector<double> >& tau13_avg_realizations =
            d_ensemble_statistics->tau13_avg_realizations;
        
        const std::vector<std::vector<double> >& tau22_avg_realizations =
            d_ensemble_statistics->tau22_avg_realizations;
        
        const std::vector<std::vector<double> >& tau23_avg_realizations =
            d_ensemble_statistics->tau23_avg_realizations;
        
        const std::vector<std::vector<double> >& tau33_avg_realizations =
            d_ensemble_statistics->tau33_avg_realizations;
        
        const std::vector<std::vector<double> >& tau11_ddx_u_avg_realizations =
            d_ensemble_statistics->tau11_ddx_u_avg_realizations;
        
        const std::vector<std::vector<double> >& tau12_ddy_u_avg_realizations =
            d_ensemble_statistics->tau12_ddy_u_avg_realizations;
        
        const std::vector<std::vector<double> >& tau13_ddz_u_avg_realizations =
            d_ensemble_statistics->tau13_ddz_u_avg_realizations;
        
        const std::vector<std::vector<double> >& tau12_ddx_v_avg_realizations =
            d_ensemble_statistics->tau12_ddx_v_avg_realizations;
        
        const std::vector<std::vector<double> >& tau22_ddy_v_avg_realizations =
            d_ensemble_statistics->tau22_ddy_v_avg_realizations;
        
        const std::vector<std::vector<double> >& tau23_ddz_v_avg_realizations =
            d_ensemble_statistics->tau23_ddz_v_avg_realizations;
        
        const std::vector<std::vector<double> >& tau13_ddx_w_avg_realizations =
            d_ensemble_statistics->tau13_ddx_w_avg_realizations;
        
        const std::vector<std::vector<double> >& tau23_ddy_w_avg_realizations =
            d_ensemble_statistics->tau23_ddy_w_avg_realizations;
        
        const std::vector<std::vector<double> >& tau33_ddz_w_avg_realizations =
            d_ensemble_statistics->tau33_ddz_w_avg_realizations;
        
        const std::vector<std::vector<double> >& Y_0_avg_realizations =
            d_ensemble_statistics->Y_0_avg_realizations;
        
        const int num_realizations = d_ensemble_statistics->getNumberOfEnsembles();
        
        TBOX_ASSERT(num_realizations > 0);
        TBOX_ASSERT(num_realizations == static_cast<int>(ddx_u_avg_realizations.size()));
        TBOX_ASSERT(num_realizations == static_cast<int>(tau11_avg_realizations.size()));
        TBOX_ASSERT(num_realizations == static_cast<int>(tau11_ddx_u_avg_realizations.size()));
        
        if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
        {
            TBOX_ASSERT(num_realizations == static_cast<int>(ddy_u_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(ddx_v_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(ddy_v_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(tau12_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(tau22_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(tau12_ddy_u_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(tau12_ddx_v_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(tau22_ddy_v_avg_realizations.size()));
        }
        if (d_dim == tbox::Dimension(3))
        {
            TBOX_ASSERT(num_realizations == static_cast<int>(ddz_u_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(ddz_v_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(ddx_w_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(ddy_w_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(ddz_w_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(tau13_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(tau23_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(tau33_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(tau13_ddz_u_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(tau23_ddz_v_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(tau13_ddx_w_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(tau23_ddy_w_avg_realizations.size()));
            TBOX_ASSERT(num_realizations == static_cast<int>(tau33_ddz_w_avg_realizations.size()));
        }
        TBOX_ASSERT(num_realizations == static_cast<int>(Y_0_avg_realizations.size()));
        
        const int num_cells = static_cast<int>(ddx_u_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> ddx_u_avg_global(num_cells, double(0));
        std::vector<double> ddy_u_avg_global(num_cells, double(0));
        std::vector<double> ddz_u_avg_global(num_cells, double(0));
        std::vector<double> ddx_v_avg_global(num_cells, double(0));
        std::vector<double> ddy_v_avg_global(num_cells, double(0));
        std::vector<double> ddz_v_avg_global(num_cells, double(0));
        std::vector<double> ddx_w_avg_global(num_cells, double(0));
        std::vector<double> ddy_w_avg_global(num_cells, double(0));
        std::vector<double> ddz_w_avg_global(num_cells, double(0));
        std::vector<double> tau11_avg_global(num_cells, double(0));
        std::vector<double> tau12_avg_global(num_cells, double(0));
        std::vector<double> tau13_avg_global(num_cells, double(0));
        std::vector<double> tau22_avg_global(num_cells, double(0));
        std::vector<double> tau23_avg_global(num_cells, double(0));
        std::vector<double> tau33_avg_global(num_cells, double(0));
        std::vector<double> tau11_ddx_u_avg_global(num_cells, double(0));
        std::vector<double> tau12_ddy_u_avg_global(num_cells, double(0));
        std::vector<double> tau13_ddz_u_avg_global(num_cells, double(0));
        std::vector<double> tau12_ddx_v_avg_global(num_cells, double(0));
        std::vector<double> tau22_ddy_v_avg_global(num_cells, double(0));
        std::vector<double> tau23_ddz_v_avg_global(num_cells, double(0));
        std::vector<double> tau13_ddx_w_avg_global(num_cells, double(0));
        std::vector<double> tau23_ddy_w_avg_global(num_cells, double(0));
        std::vector<double> tau33_ddz_w_avg_global(num_cells, double(0));
        std::vector<double> Y_0_avg_global(num_cells, double(0));
        
        std::vector<double> tau11_p_ddx_u_p(num_cells, double(0));
        std::vector<double> tau12_p_ddy_u_p(num_cells, double(0));
        std::vector<double> tau13_p_ddz_u_p(num_cells, double(0));
        std::vector<double> tau12_p_ddx_v_p(num_cells, double(0));
        std::vector<double> tau22_p_ddy_v_p(num_cells, double(0));
        std::vector<double> tau23_p_ddz_v_p(num_cells, double(0));
        std::vector<double> tau13_p_ddx_w_p(num_cells, double(0));
        std::vector<double> tau23_p_ddy_w_p(num_cells, double(0));
        std::vector<double> tau33_p_ddz_w_p(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                ddx_u_avg_global[i]       += weight*ddx_u_avg_realizations[ri][i];
                tau11_avg_global[i]       += weight*tau11_avg_realizations[ri][i];
                tau11_ddx_u_avg_global[i] += weight*tau11_ddx_u_avg_realizations[ri][i];
                Y_0_avg_global[i]         += weight*Y_0_avg_realizations[ri][i];
            }
        }
        
        for (int i = 0; i < num_cells; i++)
        {
            tau11_p_ddx_u_p[i] = tau11_ddx_u_avg_global[i] - tau11_avg_global[i]*ddx_u_avg_global[i];
        }
        
        if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
        {
            for (int ri = 0; ri < num_realizations; ri++)
            {
                for (int i = 0; i < num_cells; i++)
                {
                    ddy_u_avg_global[i]       += weight*ddy_u_avg_realizations[ri][i];
                    ddx_v_avg_global[i]       += weight*ddx_v_avg_realizations[ri][i];
                    ddy_v_avg_global[i]       += weight*ddy_v_avg_realizations[ri][i];
                    tau12_avg_global[i]       += weight*tau12_avg_realizations[ri][i];
                    tau22_avg_global[i]       += weight*tau22_avg_realizations[ri][i];
                    tau12_ddy_u_avg_global[i] += weight*tau12_ddy_u_avg_realizations[ri][i];
                    tau12_ddx_v_avg_global[i] += weight*tau12_ddx_v_avg_realizations[ri][i];
                    tau22_ddy_v_avg_global[i] += weight*tau22_ddy_v_avg_realizations[ri][i];
                }
            }
            
            for (int i = 0; i < num_cells; i++)
            {
                tau12_p_ddy_u_p[i] = tau12_ddy_u_avg_global[i] - tau12_avg_global[i]*ddy_u_avg_global[i];
                tau12_p_ddx_v_p[i] = tau12_ddx_v_avg_global[i] - tau12_avg_global[i]*ddx_v_avg_global[i];
                tau22_p_ddy_v_p[i] = tau22_ddy_v_avg_global[i] - tau22_avg_global[i]*ddy_v_avg_global[i];
            }
        }
        
        if (d_dim == tbox::Dimension(3))
        {
            for (int ri = 0; ri < num_realizations; ri++)
            {
                for (int i = 0; i < num_cells; i++)
                {
                    ddz_u_avg_global[i]       += weight*ddz_u_avg_realizations[ri][i];
                    ddz_v_avg_global[i]       += weight*ddz_v_avg_realizations[ri][i];
                    ddx_w_avg_global[i]       += weight*ddx_w_avg_realizations[ri][i];
                    ddy_w_avg_global[i]       += weight*ddy_w_avg_realizations[ri][i];
                    ddz_w_avg_global[i]       += weight*ddz_w_avg_realizations[ri][i];
                    tau13_avg_global[i]       += weight*tau13_avg_realizations[ri][i];
                    tau23_avg_global[i]       += weight*tau23_avg_realizations[ri][i];
                    tau33_avg_global[i]       += weight*tau33_avg_realizations[ri][i];
                    tau13_ddz_u_avg_global[i] += weight*tau13_ddz_u_avg_realizations[ri][i];
                    tau23_ddz_v_avg_global[i] += weight*tau23_ddz_v_avg_realizations[ri][i];
                    tau13_ddx_w_avg_global[i] += weight*tau13_ddx_w_avg_realizations[ri][i];
                    tau23_ddy_w_avg_global[i] += weight*tau23_ddy_w_avg_realizations[ri][i];
                    tau33_ddz_w_avg_global[i] += weight*tau33_ddz_w_avg_realizations[ri][i];
                }
            }
            
            for (int i = 0; i < num_cells; i++)
            {
                tau13_p_ddz_u_p[i] = tau13_ddz_u_avg_global[i] - tau13_avg_global[i]*ddz_u_avg_global[i];
                tau23_p_ddz_v_p[i] = tau23_ddz_v_avg_global[i] - tau23_avg_global[i]*ddz_v_avg_global[i];
                tau13_p_ddx_w_p[i] = tau13_ddx_w_avg_global[i] - tau13_avg_global[i]*ddx_w_avg_global[i];
                tau23_p_ddy_w_p[i] = tau23_ddy_w_avg_global[i] - tau23_avg_global[i]*ddy_w_avg_global[i];
                tau33_p_ddz_w_p[i] = tau33_ddz_w_avg_global[i] - tau33_avg_global[i]*ddz_w_avg_global[i];
            }
        }
        
        double TKE_diss_sum = double(0);
        int count = 0;
        
        if (d_dim == tbox::Dimension(1))
        {
            for (int i = 0; i < num_cells; i++)
            {
                const double mixing_metric = double(4)*Y_0_avg_global[i]*(double(1) - Y_0_avg_global[i]);
                if (mixing_metric > double(9)/double(10))
                {
                    TKE_diss_sum += tau11_p_ddx_u_p[i];
                    count++;
                }
            }
        }
        else if (d_dim == tbox::Dimension(2))
        {
            for (int i = 0; i < num_cells; i++)
            {
                const double mixing_metric = double(4)*Y_0_avg_global[i]*(double(1) - Y_0_avg_global[i]);
                if (mixing_metric > double(9)/double(10))
                {
                    TKE_diss_sum += tau11_p_ddx_u_p[i] + tau22_p_ddy_v_p[i] +
                        tau12_p_ddx_v_p[i] + tau12_p_ddy_u_p[i];
                    count++;
                }
            }
        }
        else if (d_dim == tbox::Dimension(3))
        {
            for (int i = 0; i < num_cells; i++)
            {
                const double mixing_metric = double(4)*Y_0_avg_global[i]*(double(1) - Y_0_avg_global[i]);
                if (mixing_metric > double(9)/double(10))
                {
                    TKE_diss_sum += tau11_p_ddx_u_p[i] + tau22_p_ddy_v_p[i] + + tau33_p_ddz_w_p[i] +
                        tau12_p_ddx_v_p[i] + tau12_p_ddy_u_p[i] +
                        tau13_p_ddx_w_p[i] + tau13_p_ddz_u_p[i] +
                        tau23_p_ddy_w_p[i] + tau23_p_ddz_v_p[i];
                    count++;
                }
            }
        }
        
        const double TKE_diss_mean = TKE_diss_sum/count;
        
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << TKE_diss_mean;
        
        f_out.close();
    }
}


/*
 * Output mean of ensemble dynamic shear viscosity in mixing layer with assumed homogeneity in y-direction (2D)
 * or yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleShearViscosityMeanWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const
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
        
        const std::vector<std::vector<double> >& mu_avg_realizations =
            d_ensemble_statistics->mu_avg_realizations;
        
        const std::vector<std::vector<double> >& Y_0_avg_realizations =
            d_ensemble_statistics->Y_0_avg_realizations;
        
        const int num_realizations = d_ensemble_statistics->getNumberOfEnsembles();
        
        TBOX_ASSERT(num_realizations > 0);
        TBOX_ASSERT(num_realizations == static_cast<int>(mu_avg_realizations.size()));
        TBOX_ASSERT(num_realizations == static_cast<int>(Y_0_avg_realizations.size()));
        
        const int num_cells = static_cast<int>(Y_0_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> mu_avg_global(num_cells, double(0));
        std::vector<double> Y_0_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                mu_avg_global[i]  += weight*mu_avg_realizations[ri][i];
                Y_0_avg_global[i] += weight*Y_0_avg_realizations[ri][i];
            }
        }
        
        double mu_sum = double(0);
        int count = 0;
        
        for (int i = 0; i < num_cells; i++)
        {
            const double mixing_metric = double(4)*Y_0_avg_global[i]*(double(1) - Y_0_avg_global[i]);
            if (mixing_metric > double(9)/double(10))
            {
                mu_sum += mu_avg_global[i];
                count++;
            }
        }
        
        const double mu_mean = mu_sum/count;
        
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << mu_mean;
        
        f_out.close();
    }
}


/*
 * Output mean of ensemble mass diffusivity in mixing layer with assumed homogeneity in y-direction (2D) or
 * yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleMassDiffusivityMeanWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const
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
        
        const std::vector<std::vector<double> >& D_avg_realizations =
            d_ensemble_statistics->D_avg_realizations;
        
        const std::vector<std::vector<double> >& Y_0_avg_realizations =
            d_ensemble_statistics->Y_0_avg_realizations;
        
        const int num_realizations = d_ensemble_statistics->getNumberOfEnsembles();
        
        TBOX_ASSERT(num_realizations > 0);
        TBOX_ASSERT(num_realizations == static_cast<int>(D_avg_realizations.size()));
        TBOX_ASSERT(num_realizations == static_cast<int>(Y_0_avg_realizations.size()));
        
        const int num_cells = static_cast<int>(Y_0_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> D_avg_global(num_cells, double(0));
        std::vector<double> Y_0_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                D_avg_global[i]   += weight*D_avg_realizations[ri][i];
                Y_0_avg_global[i] += weight*Y_0_avg_realizations[ri][i];
            }
        }
        
        double D_sum = double(0);
        int count = 0;
        
        for (int i = 0; i < num_cells; i++)
        {
            const double mixing_metric = double(4)*Y_0_avg_global[i]*(double(1) - Y_0_avg_global[i]);
            if (mixing_metric > double(9)/double(10))
            {
                D_sum += D_avg_global[i];
                count++;
            }
        }
        
        const double D_mean = D_sum/count;
        
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << D_mean;
        
        f_out.close();
    }
}


/*
 * Compute averaged shear stress component with only x direction as inhomogeneous direction.
 * component_idx:
 * 0: tau11
 * 1: tau12
 * 2: tau13
 * 3: tau22
 * 4: tau23
 * 5: tau33
 */
std::vector<double>
RTIRMIStatisticsUtilities::getAveragedShearStressComponentWithInhomogeneousXDirection(
    const int component_idx,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
    TBOX_ASSERT(d_num_ghosts_derivative == 3);
    
    std::vector<double> averaged_tau_ij;
    
    HAMERS_SHARED_PTR<FlowModel> d_flow_model_tmp = d_flow_model.lock();
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Get the refinement ratio from the finest level to the coarsest level.
     */
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    hier::IntVector ratio_finest_level_to_coarsest_level =
        patch_hierarchy->getRatioToCoarserLevel(num_levels - 1);
    for (int li = num_levels - 2; li > 0 ; li--)
    {
        ratio_finest_level_to_coarsest_level *= patch_hierarchy->getRatioToCoarserLevel(li);
    }
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    HAMERS_SHARED_PTR<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *patch_hierarchy,
            0,
            num_levels - 1));
    
    /*
     * Get the number of cells of physical domain refined to the finest level.
     */
    
    const hier::BoxContainer& physical_domain = d_grid_geometry->getPhysicalDomain();
    const hier::Box& physical_domain_box = physical_domain.front();
    const hier::IntVector& physical_domain_dims = physical_domain_box.numberCells();
    const hier::IntVector finest_level_dims = physical_domain_dims*ratio_finest_level_to_coarsest_level;
    
    /*
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        hier::IntVector num_ghosts_der = hier::IntVector::getOne(d_dim)*d_num_ghosts_derivative;
        
        const int finest_level_dim_0 = finest_level_dims[0];
        
        double* tau_ij_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_tau_ij.resize(finest_level_dim_0);
        double* tau_ij_avg_global = averaged_tau_ij.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            tau_ij_avg_local[i] = double(0);
            tau_ij_avg_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratio_to_coarsest_level =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarsest_level *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = ratio_finest_level_to_coarsest_level/ratio_to_coarsest_level;
            
            const int ratio_to_finest_level_0 = ratio_to_finest_level[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", num_ghosts_der));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTIONS", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("PRESSURE", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("TEMPERATURE", hier::IntVector::getZero(d_dim)));
                
                d_flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                d_flow_model_tmp->allocateMemoryForDerivedCellData();
                
                d_flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_velocity =
                    d_flow_model_tmp->getCellData("VELOCITY");
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_mass_fraction =
                    d_flow_model_tmp->getCellData("MASS_FRACTIONS");
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_pressure =
                    d_flow_model_tmp->getCellData("PRESSURE");
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_temperature =
                    d_flow_model_tmp->getCellData("TEMPERATURE");
                
                double* u = data_velocity->getPointer(0);
                std::vector<double*> Y;
                Y.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    Y.push_back(data_mass_fraction->getPointer(si));
                }
                double* p = data_pressure->getPointer(0);
                double* T = data_temperature->getPointer(0);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_velocity = data_velocity->getGhostCellWidth();
                const hier::IntVector num_ghosts_mass_fraction = data_mass_fraction->getGhostCellWidth();
                const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
                const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
                
                const int num_ghosts_0_velocity = num_ghosts_velocity[0];
                const int num_ghosts_0_mass_fraction = num_ghosts_mass_fraction[0];
                const int num_ghosts_0_pressure = num_ghosts_pressure[0];
                const int num_ghosts_0_temperature = num_ghosts_temperature[0];
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        /*
                         * Compute the index of the data point and count how many times the data is repeated.
                         */
                        
                        const hier::Index idx_pt(tbox::Dimension(1), idx_lo_0 + i);
                        
                        int n_overlapped = 1;
                        
                        for (hier::BoxContainer::BoxContainerConstIterator iob(
                                patch_overlapped_visible_boxes.begin());
                             iob != patch_overlapped_visible_boxes.end();
                             iob++)
                        {
                            const hier::Box& patch_overlapped_visible_box = *iob;
                            
                            if (patch_overlapped_visible_box.contains(idx_pt))
                            {
                                n_overlapped++;
                            }
                        }
                        
                        /*
                         * Compute the shear stress component.
                         */
                        
                        double value_to_add = double(0);
                        
                        // Compute linear indices of mass fraction, pressure and temperature.
                        const int idx_mass_fraction = relative_idx_lo_0 + i + num_ghosts_0_mass_fraction;
                        const int idx_pressure = relative_idx_lo_0 + i + num_ghosts_0_pressure;
                        const int idx_temperature = relative_idx_lo_0 + i + num_ghosts_0_temperature;
                        
                        if (component_idx == 0)
                        {
                            std::vector<const double*> Y_ptr;
                            Y_ptr.resize(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_ptr[si] = &Y[si][idx_mass_fraction];
                            }
                            const double mu = d_equation_of_shear_viscosity_mixing_rules->
                                getShearViscosity(
                                    &p[idx_pressure],
                                    &T[idx_temperature],
                                    Y_ptr);
                            const double mu_v = d_equation_of_bulk_viscosity_mixing_rules->
                                getBulkViscosity(
                                    &p[idx_pressure],
                                    &T[idx_temperature],
                                    Y_ptr);
                            
                            const int idx_vel_x_LLL = relative_idx_lo_0 + (i - 3) + num_ghosts_0_velocity;
                            const int idx_vel_x_LL  = relative_idx_lo_0 + (i - 2) + num_ghosts_0_velocity;
                            const int idx_vel_x_L   = relative_idx_lo_0 + (i - 1) + num_ghosts_0_velocity;
                            const int idx_vel_x_R   = relative_idx_lo_0 + (i + 1) + num_ghosts_0_velocity;
                            const int idx_vel_x_RR  = relative_idx_lo_0 + (i + 2) + num_ghosts_0_velocity;
                            const int idx_vel_x_RRR = relative_idx_lo_0 + (i + 3) + num_ghosts_0_velocity;
                            
                            const double dudx = (double(1)/double(60)*(u[idx_vel_x_RRR] - u[idx_vel_x_LLL])
                                - double(3)/double(20)*(u[idx_vel_x_RR] - u[idx_vel_x_LL])
                                + double(3)/double(4)*(u[idx_vel_x_R] - u[idx_vel_x_L]))/dx[0];
                            
                            value_to_add = (double(4)/double(3)*mu + mu_v)*dudx / ((double) n_overlapped);
                        }
                        else
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "Cannot compute shear stress component for one-dimensional problem!\n"
                                << "component_idx = " << component_idx << " given!\n"
                                << std::endl);
                        }
                        
                        /*
                         * Add the data.
                         */
                        
                        for (int ii = 0; ii < ratio_to_finest_level_0; ii++)
                        {
                            const int idx_fine = (idx_lo_0 + i)*ratio_to_finest_level_0 + ii;
                            
                            tau_ij_avg_local[idx_fine] += value_to_add;
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average of shear stress component.
         */
        
        mpi.Allreduce(
            tau_ij_avg_local,
            tau_ij_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(tau_ij_avg_local);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        hier::IntVector num_ghosts_der = hier::IntVector::getOne(d_dim)*d_num_ghosts_derivative;
        
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        
        double* tau_ij_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_tau_ij.resize(finest_level_dim_0);
        double* tau_ij_avg_global = averaged_tau_ij.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            tau_ij_avg_local[i] = double(0);
            tau_ij_avg_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratio_to_coarsest_level =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarsest_level *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = ratio_finest_level_to_coarsest_level/ratio_to_coarsest_level;
            
            const int ratio_to_finest_level_0 = ratio_to_finest_level[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", num_ghosts_der));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTIONS", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("PRESSURE", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("TEMPERATURE", hier::IntVector::getZero(d_dim)));
                
                d_flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                d_flow_model_tmp->allocateMemoryForDerivedCellData();
                
                d_flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_velocity =
                    d_flow_model_tmp->getCellData("VELOCITY");
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_mass_fraction =
                    d_flow_model_tmp->getCellData("MASS_FRACTIONS");
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_pressure =
                    d_flow_model_tmp->getCellData("PRESSURE");
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_temperature =
                    d_flow_model_tmp->getCellData("TEMPERATURE");
                
                double* u = data_velocity->getPointer(0);
                double* v = data_velocity->getPointer(1);
                std::vector<double*> Y;
                Y.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    Y.push_back(data_mass_fraction->getPointer(si));
                }
                double* p = data_pressure->getPointer(0);
                double* T = data_temperature->getPointer(0);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_velocity = data_velocity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_velocity = data_velocity->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_mass_fraction = data_mass_fraction->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mass_fraction = data_mass_fraction->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_pressure = data_pressure->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_temperature = data_temperature->getGhostBox().numberCells();
                
                const int num_ghosts_0_velocity = num_ghosts_velocity[0];
                const int num_ghosts_1_velocity = num_ghosts_velocity[1];
                const int ghostcell_dim_0_velocity = ghostcell_dims_velocity[0];
                
                const int num_ghosts_0_mass_fraction = num_ghosts_mass_fraction[0];
                const int num_ghosts_1_mass_fraction = num_ghosts_mass_fraction[1];
                const int ghostcell_dim_0_mass_fraction = ghostcell_dims_mass_fraction[0];
                
                const int num_ghosts_0_pressure = num_ghosts_pressure[0];
                const int num_ghosts_1_pressure = num_ghosts_pressure[1];
                const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
                
                const int num_ghosts_0_temperature = num_ghosts_temperature[0];
                const int num_ghosts_1_temperature = num_ghosts_temperature[1];
                const int ghostcell_dim_0_temperature = ghostcell_dims_temperature[0];
                
                const double weight = dx[1]/L_y;
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            /*
                             * Compute the index of the data point and count how many times the data is repeated.
                             */
                            
                            const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j);
                            
                            int n_overlapped = 1;
                            
                            for (hier::BoxContainer::BoxContainerConstIterator iob(
                                    patch_overlapped_visible_boxes.begin());
                                 iob != patch_overlapped_visible_boxes.end();
                                 iob++)
                            {
                                const hier::Box& patch_overlapped_visible_box = *iob;
                                
                                if (patch_overlapped_visible_box.contains(idx_pt))
                                {
                                    n_overlapped++;
                                }
                            }
                            
                            /*
                             * Compute the shear stress component.
                             */
                            
                            double value_to_add = double(0);
                            
                            const int idx_mass_fraction = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                            
                            const int idx_pressure = (relative_idx_lo_0 + i + num_ghosts_0_pressure) +
                                (relative_idx_lo_1 + j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure;
                            
                            const int idx_temperature = (relative_idx_lo_0 + i + num_ghosts_0_temperature) +
                                (relative_idx_lo_1 + j + num_ghosts_1_temperature)*ghostcell_dim_0_temperature;
                            
                            const int idx_vel_x_LLL = (relative_idx_lo_0 + (i - 3) + num_ghosts_0_velocity) +
                                (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                            
                            const int idx_vel_x_LL  = (relative_idx_lo_0 + (i - 2) + num_ghosts_0_velocity) +
                                (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                            
                            const int idx_vel_x_L   = (relative_idx_lo_0 + (i - 1) + num_ghosts_0_velocity) +
                                (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                            
                            const int idx_vel_x_R   = (relative_idx_lo_0 + (i + 1) + num_ghosts_0_velocity) +
                                (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                            
                            const int idx_vel_x_RR  = (relative_idx_lo_0 + (i + 2) + num_ghosts_0_velocity) +
                                (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                            
                            const int idx_vel_x_RRR = (relative_idx_lo_0 + (i + 3) + num_ghosts_0_velocity) +
                                (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                            
                            const int idx_vel_y_BBB = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                (relative_idx_lo_1 + (j - 3) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                            
                            const int idx_vel_y_BB  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                (relative_idx_lo_1 + (j - 2) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                            
                            const int idx_vel_y_B   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                (relative_idx_lo_1 + (j - 1) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                            
                            const int idx_vel_y_T   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                (relative_idx_lo_1 + (j + 1) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                            
                            const int idx_vel_y_TT  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                (relative_idx_lo_1 + (j + 2) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                            
                            const int idx_vel_y_TTT = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                (relative_idx_lo_1 + (j + 3) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                            
                            if (component_idx == 0)
                            {
                                std::vector<const double*> Y_ptr;
                                Y_ptr.resize(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_ptr[si] = &Y[si][idx_mass_fraction];
                                }
                                const double mu = d_equation_of_shear_viscosity_mixing_rules->
                                    getShearViscosity(
                                        &p[idx_pressure],
                                        &T[idx_temperature],
                                        Y_ptr);
                                const double mu_v = d_equation_of_bulk_viscosity_mixing_rules->
                                    getBulkViscosity(
                                        &p[idx_pressure],
                                        &T[idx_temperature],
                                        Y_ptr);
                                
                                const double dudx = (double(1)/double(60)*(u[idx_vel_x_RRR] - u[idx_vel_x_LLL])
                                    - double(3)/double(20)*(u[idx_vel_x_RR] - u[idx_vel_x_LL])
                                    + double(3)/double(4)*(u[idx_vel_x_R] - u[idx_vel_x_L]))/dx[0];
                                
                                const double dvdy = (double(1)/double(60)*(v[idx_vel_y_TTT] - v[idx_vel_y_BBB])
                                    - double(3)/double(20)*(v[idx_vel_y_TT] - v[idx_vel_y_BB])
                                    + double(3)/double(4)*(v[idx_vel_y_T] - v[idx_vel_y_B]))/dx[1];
                                
                                value_to_add = ((double(4)/double(3)*mu + mu_v)*dudx - (double(2)/double(3)*mu - mu_v)*dvdy)
                                    *weight/((double) n_overlapped);
                            }
                            else if (component_idx == 1)
                            {
                                std::vector<const double*> Y_ptr;
                                Y_ptr.resize(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_ptr[si] = &Y[si][idx_mass_fraction];
                                }
                                const double mu = d_equation_of_shear_viscosity_mixing_rules->
                                    getShearViscosity(
                                        &p[idx_pressure],
                                        &T[idx_temperature],
                                        Y_ptr);
                                
                                const double dvdx = (double(1)/double(60)*(v[idx_vel_x_RRR] - v[idx_vel_x_LLL])
                                    - double(3)/double(20)*(v[idx_vel_x_RR] - v[idx_vel_x_LL])
                                    + double(3)/double(4)*(v[idx_vel_x_R] - v[idx_vel_x_L]))/dx[0];
                                
                                const double dudy = (double(1)/double(60)*(u[idx_vel_y_TTT] - u[idx_vel_y_BBB])
                                    - double(3)/double(20)*(u[idx_vel_y_TT] - u[idx_vel_y_BB])
                                    + double(3)/double(4)*(u[idx_vel_y_T] - u[idx_vel_y_B]))/dx[1];
                                
                                value_to_add = mu*(dudy + dvdx)
                                    *weight/((double) n_overlapped);
                            }
                            else if (component_idx == 3)
                            {
                                std::vector<const double*> Y_ptr;
                                Y_ptr.resize(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_ptr[si] = &Y[si][idx_mass_fraction];
                                }
                                const double mu = d_equation_of_shear_viscosity_mixing_rules->
                                    getShearViscosity(
                                        &p[idx_pressure],
                                        &T[idx_temperature],
                                        Y_ptr);
                                const double mu_v = d_equation_of_bulk_viscosity_mixing_rules->
                                    getBulkViscosity(
                                        &p[idx_pressure],
                                        &T[idx_temperature],
                                        Y_ptr);
                                
                                const double dudx = (double(1)/double(60)*(u[idx_vel_x_RRR] - u[idx_vel_x_LLL])
                                    - double(3)/double(20)*(u[idx_vel_x_RR] - u[idx_vel_x_LL])
                                    + double(3)/double(4)*(u[idx_vel_x_R] - u[idx_vel_x_L]))/dx[0];
                                
                                const double dvdy = (double(1)/double(60)*(v[idx_vel_y_TTT] - v[idx_vel_y_BBB])
                                    - double(3)/double(20)*(v[idx_vel_y_TT] - v[idx_vel_y_BB])
                                    + double(3)/double(4)*(v[idx_vel_y_T] - v[idx_vel_y_B]))/dx[1];
                                
                                value_to_add = ((double(4)/double(3)*mu + mu_v)*dvdy - (double(2)/double(3)*mu - mu_v)*dudx)
                                    *weight/((double) n_overlapped);
                            }
                            else
                            {
                                TBOX_ERROR(d_object_name
                                    << ": "
                                    << "Cannot compute shear stress component for two-dimensional problem!\n"
                                    << "component_idx = " << component_idx << " given!\n"
                                    << std::endl);
                            }
                            
                            /*
                             * Add the data.
                             */
                            
                            for (int ii = 0; ii < ratio_to_finest_level_0; ii++)
                            {
                                const int idx_fine = (idx_lo_0 + i)*ratio_to_finest_level_0 + ii;
                                
                                tau_ij_avg_local[idx_fine] += value_to_add;
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average of shear stress component.
         */
        
        mpi.Allreduce(
            tau_ij_avg_local,
            tau_ij_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(tau_ij_avg_local);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        hier::IntVector num_ghosts_der = hier::IntVector::getOne(d_dim)*d_num_ghosts_derivative;
        
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        const double L_z = x_hi[2] - x_lo[2];
        
        double* tau_ij_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_tau_ij.resize(finest_level_dim_0);
        double* tau_ij_avg_global = averaged_tau_ij.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            tau_ij_avg_local[i] = double(0);
            tau_ij_avg_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratio_to_coarsest_level =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarsest_level *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = ratio_finest_level_to_coarsest_level/ratio_to_coarsest_level;
            
            const int ratio_to_finest_level_0 = ratio_to_finest_level[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", num_ghosts_der));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTIONS", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("PRESSURE", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("TEMPERATURE", hier::IntVector::getZero(d_dim)));
                
                d_flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                d_flow_model_tmp->allocateMemoryForDerivedCellData();
                
                d_flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_velocity =
                    d_flow_model_tmp->getCellData("VELOCITY");
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_mass_fraction =
                    d_flow_model_tmp->getCellData("MASS_FRACTIONS");
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_pressure =
                    d_flow_model_tmp->getCellData("PRESSURE");
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_temperature =
                    d_flow_model_tmp->getCellData("TEMPERATURE");
                
                double* u = data_velocity->getPointer(0);
                double* v = data_velocity->getPointer(1);
                double* w = data_velocity->getPointer(2);
                std::vector<double*> Y;
                Y.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    Y.push_back(data_mass_fraction->getPointer(si));
                }
                double* p = data_pressure->getPointer(0);
                double* T = data_temperature->getPointer(0);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_velocity = data_velocity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_velocity = data_velocity->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_mass_fraction = data_mass_fraction->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mass_fraction = data_mass_fraction->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_pressure = data_pressure->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_temperature = data_temperature->getGhostBox().numberCells();
                
                const int num_ghosts_0_velocity = num_ghosts_velocity[0];
                const int num_ghosts_1_velocity = num_ghosts_velocity[1];
                const int num_ghosts_2_velocity = num_ghosts_velocity[2];
                const int ghostcell_dim_0_velocity = ghostcell_dims_velocity[0];
                const int ghostcell_dim_1_velocity = ghostcell_dims_velocity[1];
                
                const int num_ghosts_0_mass_fraction = num_ghosts_mass_fraction[0];
                const int num_ghosts_1_mass_fraction = num_ghosts_mass_fraction[1];
                const int num_ghosts_2_mass_fraction = num_ghosts_mass_fraction[2];
                const int ghostcell_dim_0_mass_fraction = ghostcell_dims_mass_fraction[0];
                const int ghostcell_dim_1_mass_fraction = ghostcell_dims_mass_fraction[1];
                
                const int num_ghosts_0_pressure = num_ghosts_pressure[0];
                const int num_ghosts_1_pressure = num_ghosts_pressure[1];
                const int num_ghosts_2_pressure = num_ghosts_pressure[2];
                const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
                const int ghostcell_dim_1_pressure = ghostcell_dims_pressure[1];
                
                const int num_ghosts_0_temperature = num_ghosts_temperature[0];
                const int num_ghosts_1_temperature = num_ghosts_temperature[1];
                const int num_ghosts_2_temperature = num_ghosts_temperature[2];
                const int ghostcell_dim_0_temperature = ghostcell_dims_temperature[0];
                const int ghostcell_dim_1_temperature = ghostcell_dims_temperature[1];
                
                const double weight = (dx[1]*dx[2])/(L_y*L_z);
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    const int interior_dim_2 = interior_dims[2];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int idx_lo_2 = index_lo[2];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    const int relative_idx_lo_2 = relative_index_lo[2];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                /*
                                 * Compute the index of the data point and count how many times the data is repeated.
                                 */
                                
                                const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j, idx_lo_2 + k);
                                
                                int n_overlapped = 1;
                                
                                for (hier::BoxContainer::BoxContainerConstIterator iob(
                                        patch_overlapped_visible_boxes.begin());
                                     iob != patch_overlapped_visible_boxes.end();
                                     iob++)
                                {
                                    const hier::Box& patch_overlapped_visible_box = *iob;
                                    
                                    if (patch_overlapped_visible_box.contains(idx_pt))
                                    {
                                        n_overlapped++;
                                    }
                                }
                                
                                /*
                                 * Compute the shear stress component.
                                 */
                                
                                double value_to_add = double(0);
                                
                                const int idx_mass_fraction = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                    (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                        ghostcell_dim_1_mass_fraction;
                                
                                const int idx_pressure = (relative_idx_lo_0 + i + num_ghosts_0_pressure) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure +
                                    (relative_idx_lo_2 + k + num_ghosts_2_pressure)*ghostcell_dim_0_pressure*
                                        ghostcell_dim_1_pressure;
                                
                                const int idx_temperature = (relative_idx_lo_0 + i + num_ghosts_0_temperature) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_temperature)*ghostcell_dim_0_temperature +
                                    (relative_idx_lo_2 + k + num_ghosts_2_temperature)*ghostcell_dim_0_temperature*
                                        ghostcell_dim_1_temperature;
                                
                                if (component_idx == 0)
                                {
                                    std::vector<const double*> Y_ptr;
                                    Y_ptr.resize(d_num_species);
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        Y_ptr[si] = &Y[si][idx_mass_fraction];
                                    }
                                    const double mu = d_equation_of_shear_viscosity_mixing_rules->
                                        getShearViscosity(
                                            &p[idx_pressure],
                                            &T[idx_temperature],
                                            Y_ptr);
                                    const double mu_v = d_equation_of_bulk_viscosity_mixing_rules->
                                        getBulkViscosity(
                                            &p[idx_pressure],
                                            &T[idx_temperature],
                                            Y_ptr);
                                    
                                    const int idx_vel_x_LLL = (relative_idx_lo_0 + (i - 3) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_LL  = (relative_idx_lo_0 + (i - 2) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_L   = (relative_idx_lo_0 + (i - 1) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_R   = (relative_idx_lo_0 + (i + 1) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_RR  = (relative_idx_lo_0 + (i + 2) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_RRR = (relative_idx_lo_0 + (i + 3) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const double dudx = (double(1)/double(60)*(u[idx_vel_x_RRR] - u[idx_vel_x_LLL])
                                        - double(3)/double(20)*(u[idx_vel_x_RR] - u[idx_vel_x_LL])
                                        + double(3)/double(4)*(u[idx_vel_x_R] - u[idx_vel_x_L]))/dx[0];
                                    
                                    const int idx_vel_y_BBB = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j - 3) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_BB  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j - 2) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_B   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j - 1) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_T   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + 1) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_TT  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + 2) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_TTT = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + 3) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const double dvdy = (double(1)/double(60)*(v[idx_vel_y_TTT] - v[idx_vel_y_BBB])
                                        - double(3)/double(20)*(v[idx_vel_y_TT] - v[idx_vel_y_BB])
                                        + double(3)/double(4)*(v[idx_vel_y_T] - v[idx_vel_y_B]))/dx[1];
                                    
                                    const int idx_vel_z_BBB = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k - 3) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_BB  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k - 2) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_B   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k - 1) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_F   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k + 1) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_FF  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k + 2) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_FFF = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k + 3) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const double dwdz = (double(1)/double(60)*(w[idx_vel_z_FFF] - w[idx_vel_z_BBB])
                                        - double(3)/double(20)*(w[idx_vel_z_FF] - w[idx_vel_z_BB])
                                        + double(3)/double(4)*(w[idx_vel_z_F] - w[idx_vel_z_B]))/dx[2];
                                    
                                    value_to_add = ((double(4)/double(3)*mu + mu_v)*dudx - (double(2)/double(3)*mu - mu_v)*(dvdy + dwdz))
                                        *weight/((double) n_overlapped);
                                }
                                else if (component_idx == 1)
                                {
                                    std::vector<const double*> Y_ptr;
                                    Y_ptr.resize(d_num_species);
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        Y_ptr[si] = &Y[si][idx_mass_fraction];
                                    }
                                    const double mu = d_equation_of_shear_viscosity_mixing_rules->
                                        getShearViscosity(
                                            &p[idx_pressure],
                                            &T[idx_temperature],
                                            Y_ptr);
                                    
                                    const int idx_vel_x_LLL = (relative_idx_lo_0 + (i - 3) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_LL  = (relative_idx_lo_0 + (i - 2) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_L   = (relative_idx_lo_0 + (i - 1) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_R   = (relative_idx_lo_0 + (i + 1) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_RR  = (relative_idx_lo_0 + (i + 2) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_RRR = (relative_idx_lo_0 + (i + 3) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const double dvdx = (double(1)/double(60)*(v[idx_vel_x_RRR] - v[idx_vel_x_LLL])
                                        - double(3)/double(20)*(v[idx_vel_x_RR] - v[idx_vel_x_LL])
                                        + double(3)/double(4)*(v[idx_vel_x_R] - v[idx_vel_x_L]))/dx[0];
                                    
                                    const int idx_vel_y_BBB = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j - 3) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_BB  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j - 2) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_B   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j - 1) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_T   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + 1) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_TT  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + 2) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_TTT = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + 3) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const double dudy = (double(1)/double(60)*(u[idx_vel_y_TTT] - u[idx_vel_y_BBB])
                                        - double(3)/double(20)*(u[idx_vel_y_TT] - u[idx_vel_y_BB])
                                        + double(3)/double(4)*(u[idx_vel_y_T] - u[idx_vel_y_B]))/dx[1];
                                    
                                    value_to_add = mu*(dudy + dvdx)
                                        *weight/((double) n_overlapped);
                                }
                                else if (component_idx == 2)
                                {
                                    std::vector<const double*> Y_ptr;
                                    Y_ptr.resize(d_num_species);
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        Y_ptr[si] = &Y[si][idx_mass_fraction];
                                    }
                                    const double mu = d_equation_of_shear_viscosity_mixing_rules->
                                        getShearViscosity(
                                            &p[idx_pressure],
                                            &T[idx_temperature],
                                            Y_ptr);
                                    
                                    const int idx_vel_x_LLL = (relative_idx_lo_0 + (i - 3) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_LL  = (relative_idx_lo_0 + (i - 2) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_L   = (relative_idx_lo_0 + (i - 1) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_R   = (relative_idx_lo_0 + (i + 1) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_RR  = (relative_idx_lo_0 + (i + 2) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_RRR = (relative_idx_lo_0 + (i + 3) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const double dwdx = (double(1)/double(60)*(w[idx_vel_x_RRR] - w[idx_vel_x_LLL])
                                        - double(3)/double(20)*(w[idx_vel_x_RR] - w[idx_vel_x_LL])
                                        + double(3)/double(4)*(w[idx_vel_x_R] - w[idx_vel_x_L]))/dx[0];
                                    
                                    const int idx_vel_z_BBB = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k - 3) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_BB  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k - 2) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_B   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k - 1) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_F   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k + 1) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_FF  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k + 2) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_FFF = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k + 3) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const double dudz = (double(1)/double(60)*(u[idx_vel_z_FFF] - u[idx_vel_z_BBB])
                                        - double(3)/double(20)*(u[idx_vel_z_FF] - u[idx_vel_z_BB])
                                        + double(3)/double(4)*(u[idx_vel_z_F] - u[idx_vel_z_B]))/dx[2];
                                    
                                    value_to_add = mu*(dudz + dwdx)
                                        *weight/((double) n_overlapped);
                                }
                                else if (component_idx == 3)
                                {
                                    std::vector<const double*> Y_ptr;
                                    Y_ptr.resize(d_num_species);
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        Y_ptr[si] = &Y[si][idx_mass_fraction];
                                    }
                                    const double mu = d_equation_of_shear_viscosity_mixing_rules->
                                        getShearViscosity(
                                            &p[idx_pressure],
                                            &T[idx_temperature],
                                            Y_ptr);
                                    const double mu_v = d_equation_of_bulk_viscosity_mixing_rules->
                                        getBulkViscosity(
                                            &p[idx_pressure],
                                            &T[idx_temperature],
                                            Y_ptr);
                                    
                                    const int idx_vel_x_LLL = (relative_idx_lo_0 + (i - 3) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_LL  = (relative_idx_lo_0 + (i - 2) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_L   = (relative_idx_lo_0 + (i - 1) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_R   = (relative_idx_lo_0 + (i + 1) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_RR  = (relative_idx_lo_0 + (i + 2) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_RRR = (relative_idx_lo_0 + (i + 3) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const double dudx = (double(1)/double(60)*(u[idx_vel_x_RRR] - u[idx_vel_x_LLL])
                                        - double(3)/double(20)*(u[idx_vel_x_RR] - u[idx_vel_x_LL])
                                        + double(3)/double(4)*(u[idx_vel_x_R] - u[idx_vel_x_L]))/dx[0];
                                    
                                    const int idx_vel_y_BBB = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j - 3) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_BB  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j - 2) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_B   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j - 1) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_T   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + 1) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_TT  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + 2) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_TTT = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + 3) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const double dvdy = (double(1)/double(60)*(v[idx_vel_y_TTT] - v[idx_vel_y_BBB])
                                        - double(3)/double(20)*(v[idx_vel_y_TT] - v[idx_vel_y_BB])
                                        + double(3)/double(4)*(v[idx_vel_y_T] - v[idx_vel_y_B]))/dx[1];
                                    
                                    const int idx_vel_z_BBB = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k - 3) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_BB  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k - 2) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_B   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k - 1) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_F   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k + 1) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_FF  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k + 2) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_FFF = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k + 3) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const double dwdz = (double(1)/double(60)*(w[idx_vel_z_FFF] - w[idx_vel_z_BBB])
                                        - double(3)/double(20)*(w[idx_vel_z_FF] - w[idx_vel_z_BB])
                                        + double(3)/double(4)*(w[idx_vel_z_F] - w[idx_vel_z_B]))/dx[2];
                                    
                                    value_to_add = ((double(4)/double(3)*mu + mu_v)*dvdy - (double(2)/double(3)*mu - mu_v)*(dudx + dwdz))
                                        *weight/((double) n_overlapped);
                                }
                                else if (component_idx == 4)
                                {
                                    std::vector<const double*> Y_ptr;
                                    Y_ptr.resize(d_num_species);
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        Y_ptr[si] = &Y[si][idx_mass_fraction];
                                    }
                                    const double mu = d_equation_of_shear_viscosity_mixing_rules->
                                        getShearViscosity(
                                            &p[idx_pressure],
                                            &T[idx_temperature],
                                            Y_ptr);
                                    
                                    const int idx_vel_y_BBB = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j - 3) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_BB  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j - 2) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_B   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j - 1) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_T   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + 1) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_TT  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + 2) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_TTT = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + 3) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const double dwdy = (double(1)/double(60)*(w[idx_vel_y_TTT] - w[idx_vel_y_BBB])
                                        - double(3)/double(20)*(w[idx_vel_y_TT] - w[idx_vel_y_BB])
                                        + double(3)/double(4)*(w[idx_vel_y_T] - w[idx_vel_y_B]))/dx[1];
                                    
                                    const int idx_vel_z_BBB = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k - 3) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_BB  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k - 2) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_B   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k - 1) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_F   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k + 1) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_FF  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k + 2) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_FFF = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k + 3) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const double dvdz = (double(1)/double(60)*(v[idx_vel_z_FFF] - v[idx_vel_z_BBB])
                                        - double(3)/double(20)*(v[idx_vel_z_FF] - v[idx_vel_z_BB])
                                        + double(3)/double(4)*(v[idx_vel_z_F] - v[idx_vel_z_B]))/dx[2];
                                    
                                    value_to_add = mu*(dvdz + dwdy)
                                        *weight/((double) n_overlapped);
                                }
                                else if (component_idx == 5)
                                {
                                    std::vector<const double*> Y_ptr;
                                    Y_ptr.resize(d_num_species);
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        Y_ptr[si] = &Y[si][idx_mass_fraction];
                                    }
                                    const double mu = d_equation_of_shear_viscosity_mixing_rules->
                                        getShearViscosity(
                                            &p[idx_pressure],
                                            &T[idx_temperature],
                                            Y_ptr);
                                    const double mu_v = d_equation_of_bulk_viscosity_mixing_rules->
                                        getBulkViscosity(
                                            &p[idx_pressure],
                                            &T[idx_temperature],
                                            Y_ptr);
                                    
                                    const int idx_vel_x_LLL = (relative_idx_lo_0 + (i - 3) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_LL  = (relative_idx_lo_0 + (i - 2) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_L   = (relative_idx_lo_0 + (i - 1) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_R   = (relative_idx_lo_0 + (i + 1) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_RR  = (relative_idx_lo_0 + (i + 2) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_RRR = (relative_idx_lo_0 + (i + 3) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const double dudx = (double(1)/double(60)*(u[idx_vel_x_RRR] - u[idx_vel_x_LLL])
                                        - double(3)/double(20)*(u[idx_vel_x_RR] - u[idx_vel_x_LL])
                                        + double(3)/double(4)*(u[idx_vel_x_R] - u[idx_vel_x_L]))/dx[0];
                                    
                                    const int idx_vel_y_BBB = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j - 3) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_BB  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j - 2) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_B   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j - 1) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_T   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + 1) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_TT  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + 2) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_TTT = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + 3) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const double dvdy = (double(1)/double(60)*(v[idx_vel_y_TTT] - v[idx_vel_y_BBB])
                                        - double(3)/double(20)*(v[idx_vel_y_TT] - v[idx_vel_y_BB])
                                        + double(3)/double(4)*(v[idx_vel_y_T] - v[idx_vel_y_B]))/dx[1];
                                    
                                    const int idx_vel_z_BBB = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k - 3) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_BB  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k - 2) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_B   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k - 1) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_F   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k + 1) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_FF  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k + 2) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_FFF = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k + 3) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const double dwdz = (double(1)/double(60)*(w[idx_vel_z_FFF] - w[idx_vel_z_BBB])
                                        - double(3)/double(20)*(w[idx_vel_z_FF] - w[idx_vel_z_BB])
                                        + double(3)/double(4)*(w[idx_vel_z_F] - w[idx_vel_z_B]))/dx[2];
                                    
                                    value_to_add = ((double(4)/double(3)*mu + mu_v)*dwdz - (double(2)/double(3)*mu - mu_v)*(dudx + dvdy))
                                        *weight/((double) n_overlapped);
                                }
                                else
                                {
                                    TBOX_ERROR(d_object_name
                                        << ": "
                                        << "Cannot compute shear stress component for two-dimensional problem!\n"
                                        << "component_idx = " << component_idx << " given!\n"
                                        << std::endl);
                                }
                                
                                /*
                                 * Add the data.
                                 */
                                
                                for (int ii = 0; ii < ratio_to_finest_level_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratio_to_finest_level_0 + ii;
                                    
                                    tau_ij_avg_local[idx_fine] += value_to_add;
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average of shear stress component.
         */
        
        mpi.Allreduce(
            tau_ij_avg_local,
            tau_ij_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(tau_ij_avg_local);
    }
    
    return  averaged_tau_ij;
}



/*
 * Compute averaged derivative of shear stress component with only x direction as inhomogeneous direction.
 * component_idx:
 * 0: tau11
 * 1: tau12
 * 2: tau13
 * 3: tau22
 * 4: tau23
 * 5: tau33
 */
std::vector<double>
RTIRMIStatisticsUtilities::getAveragedDerivativeOfShearStressComponentWithInhomogeneousXDirection(
    const int component_idx,
    const int derivative_direction,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
    TBOX_ASSERT(d_num_ghosts_derivative == 3);
    
    std::vector<double> averaged_derivative_tau_ij;
    
    HAMERS_SHARED_PTR<FlowModel> d_flow_model_tmp = d_flow_model.lock();
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Get the refinement ratio from the finest level to the coarsest level.
     */
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    hier::IntVector ratio_finest_level_to_coarsest_level =
        patch_hierarchy->getRatioToCoarserLevel(num_levels - 1);
    for (int li = num_levels - 2; li > 0 ; li--)
    {
        ratio_finest_level_to_coarsest_level *= patch_hierarchy->getRatioToCoarserLevel(li);
    }
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    HAMERS_SHARED_PTR<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *patch_hierarchy,
            0,
            num_levels - 1));
    
    /*
     * Get the number of cells of physical domain refined to the finest level.
     */
    
    const hier::BoxContainer& physical_domain = d_grid_geometry->getPhysicalDomain();
    const hier::Box& physical_domain_box = physical_domain.front();
    const hier::IntVector& physical_domain_dims = physical_domain_box.numberCells();
    const hier::IntVector finest_level_dims = physical_domain_dims*ratio_finest_level_to_coarsest_level;
    
    /*
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        hier::IntVector num_ghosts_der = hier::IntVector::getOne(d_dim)*d_num_ghosts_derivative;
        
        const int finest_level_dim_0 = finest_level_dims[0];
        
        double* tau_ij_der_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_derivative_tau_ij.resize(finest_level_dim_0);
        double* tau_ij_der_avg_global = averaged_derivative_tau_ij.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            tau_ij_der_avg_local[i] = double(0);
            tau_ij_der_avg_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratio_to_coarsest_level =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarsest_level *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = ratio_finest_level_to_coarsest_level/ratio_to_coarsest_level;
            
            const int ratio_to_finest_level_0 = ratio_to_finest_level[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and the quantities in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", num_ghosts_der));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTIONS", num_ghosts_der));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("PRESSURE", num_ghosts_der));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("TEMPERATURE", num_ghosts_der));
                
                d_flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                d_flow_model_tmp->allocateMemoryForDerivedCellData();
                
                d_flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_velocity =
                    d_flow_model_tmp->getCellData("VELOCITY");
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_mass_fraction =
                    d_flow_model_tmp->getCellData("MASS_FRACTIONS");
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_pressure =
                    d_flow_model_tmp->getCellData("PRESSURE");
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_temperature =
                    d_flow_model_tmp->getCellData("TEMPERATURE");
                
                double* u = data_velocity->getPointer(0);
                std::vector<double*> Y;
                Y.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    Y.push_back(data_mass_fraction->getPointer(si));
                }
                double* p = data_pressure->getPointer(0);
                double* T = data_temperature->getPointer(0);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_velocity = data_velocity->getGhostCellWidth();
                const hier::IntVector num_ghosts_mass_fraction = data_mass_fraction->getGhostCellWidth();
                const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
                const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
                
                const int num_ghosts_0_velocity = num_ghosts_velocity[0];
                const int num_ghosts_0_mass_fraction = num_ghosts_mass_fraction[0];
                const int num_ghosts_0_pressure = num_ghosts_pressure[0];
                const int num_ghosts_0_temperature = num_ghosts_temperature[0];
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        /*
                         * Compute the index of the data point and count how many times the data is repeated.
                         */
                        
                        const hier::Index idx_pt(tbox::Dimension(1), idx_lo_0 + i);
                        
                        int n_overlapped = 1;
                        
                        for (hier::BoxContainer::BoxContainerConstIterator iob(
                                patch_overlapped_visible_boxes.begin());
                             iob != patch_overlapped_visible_boxes.end();
                             iob++)
                        {
                            const hier::Box& patch_overlapped_visible_box = *iob;
                            
                            if (patch_overlapped_visible_box.contains(idx_pt))
                            {
                                n_overlapped++;
                            }
                        }
                        
                        /*
                         * Compute the derivative of shear stress component.
                         */
                        
                        double value_to_add = double(0);
                        
                        if ((component_idx == 0) && (derivative_direction == 0))
                        {
                            double mu[7];
                            double mu_v[7];
                            
                            std::vector<const double*> Y_ptr;
                            Y_ptr.resize(d_num_species);
                            for (int count = -3; count <= 3; count++)
                            {
                                const int idx_mass_fraction = relative_idx_lo_0 + (i + count) + num_ghosts_0_mass_fraction;
                                const int idx_pressure = relative_idx_lo_0 + (i + count) + num_ghosts_0_pressure;
                                const int idx_temperature = relative_idx_lo_0 + (i + count) + num_ghosts_0_temperature;
                                
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_ptr[si] = &Y[si][idx_mass_fraction];
                                }
                                mu[count + 3] = d_equation_of_shear_viscosity_mixing_rules->
                                    getShearViscosity(
                                        &p[idx_pressure],
                                        &T[idx_temperature],
                                        Y_ptr);
                                mu_v[count + 3] = d_equation_of_bulk_viscosity_mixing_rules->
                                    getBulkViscosity(
                                        &p[idx_pressure],
                                        &T[idx_temperature],
                                        Y_ptr);
                            }
                            
                            const double dmudx = (double(1)/double(60)*(mu[6] - mu[0])
                                - double(3)/double(20)*(mu[5] - mu[1])
                                + double(3)/double(4)*(mu[4] - mu[2]))/dx[0];
                            
                            const double dmu_vdx = (double(1)/double(60)*(mu_v[6] - mu_v[0])
                                - double(3)/double(20)*(mu_v[5] - mu_v[1])
                                + double(3)/double(4)*(mu_v[4] - mu_v[2]))/dx[0];
                            
                            const int idx_vel       = relative_idx_lo_0 +  i      + num_ghosts_0_velocity;
                            
                            const int idx_vel_x_LLL = relative_idx_lo_0 + (i - 3) + num_ghosts_0_velocity;
                            const int idx_vel_x_LL  = relative_idx_lo_0 + (i - 2) + num_ghosts_0_velocity;
                            const int idx_vel_x_L   = relative_idx_lo_0 + (i - 1) + num_ghosts_0_velocity;
                            const int idx_vel_x_R   = relative_idx_lo_0 + (i + 1) + num_ghosts_0_velocity;
                            const int idx_vel_x_RR  = relative_idx_lo_0 + (i + 2) + num_ghosts_0_velocity;
                            const int idx_vel_x_RRR = relative_idx_lo_0 + (i + 3) + num_ghosts_0_velocity;
                            
                            const double dudx = (double(1)/double(60)*(u[idx_vel_x_RRR] - u[idx_vel_x_LLL])
                                - double(3)/double(20)*(u[idx_vel_x_RR] - u[idx_vel_x_LL])
                                + double(3)/double(4)*(u[idx_vel_x_R] - u[idx_vel_x_L]))/dx[0];
                            
                            const double d2udx2 = (double(1)/double(90)*(u[idx_vel_x_RRR] + u[idx_vel_x_LLL])
                                - double(3)/double(20)*(u[idx_vel_x_RR] + u[idx_vel_x_LL])
                                + double(3)/double(2)*(u[idx_vel_x_R] + u[idx_vel_x_L])
                                - double(49)/double(18)*u[idx_vel])/(dx[0]*dx[0]);
                            
                            value_to_add = ((double(4)/double(3)*dmudx + dmu_vdx)*dudx 
                                + (double(4)/double(3)*mu[3] + mu_v[3])*d2udx2) / ((double) n_overlapped);
                        }
                        else
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "Cannot compute derivative of shear stress component for one-dimensional problem!\n"
                                << "component_idx = " << component_idx << " and derivative_direction = "
                                << derivative_direction << " given!\n"
                                << std::endl);
                        }
                        
                        /*
                         * Add the data.
                         */
                        
                        for (int ii = 0; ii < ratio_to_finest_level_0; ii++)
                        {
                            const int idx_fine = (idx_lo_0 + i)*ratio_to_finest_level_0 + ii;
                            
                            tau_ij_der_avg_local[idx_fine] += value_to_add;
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average of derivative of shear stress component.
         */
        
        mpi.Allreduce(
            tau_ij_der_avg_local,
            tau_ij_der_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(tau_ij_der_avg_local);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        hier::IntVector num_ghosts_der = hier::IntVector::getOne(d_dim)*d_num_ghosts_derivative;
        
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        
        double* tau_ij_der_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_derivative_tau_ij.resize(finest_level_dim_0);
        double* tau_ij_der_avg_global = averaged_derivative_tau_ij.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            tau_ij_der_avg_local[i] = double(0);
            tau_ij_der_avg_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratio_to_coarsest_level =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarsest_level *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = ratio_finest_level_to_coarsest_level/ratio_to_coarsest_level;
            
            const int ratio_to_finest_level_0 = ratio_to_finest_level[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and the quantities in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", num_ghosts_der));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTIONS", num_ghosts_der));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("PRESSURE", num_ghosts_der));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("TEMPERATURE", num_ghosts_der));
                
                d_flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                d_flow_model_tmp->allocateMemoryForDerivedCellData();
                
                d_flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_velocity =
                    d_flow_model_tmp->getCellData("VELOCITY");
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_mass_fraction =
                    d_flow_model_tmp->getCellData("MASS_FRACTIONS");
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_pressure =
                    d_flow_model_tmp->getCellData("PRESSURE");
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_temperature =
                    d_flow_model_tmp->getCellData("TEMPERATURE");
                
                double* u = data_velocity->getPointer(0);
                double* v = data_velocity->getPointer(1);
                std::vector<double*> Y;
                Y.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    Y.push_back(data_mass_fraction->getPointer(si));
                }
                double* p = data_pressure->getPointer(0);
                double* T = data_temperature->getPointer(0);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_velocity = data_velocity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_velocity = data_velocity->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_mass_fraction = data_mass_fraction->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mass_fraction = data_mass_fraction->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_pressure = data_pressure->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_temperature = data_temperature->getGhostBox().numberCells();
                
                const int num_ghosts_0_velocity = num_ghosts_velocity[0];
                const int num_ghosts_1_velocity = num_ghosts_velocity[1];
                const int ghostcell_dim_0_velocity = ghostcell_dims_velocity[0];
                
                const int num_ghosts_0_mass_fraction = num_ghosts_mass_fraction[0];
                const int num_ghosts_1_mass_fraction = num_ghosts_mass_fraction[1];
                const int ghostcell_dim_0_mass_fraction = ghostcell_dims_mass_fraction[0];
                
                const int num_ghosts_0_pressure = num_ghosts_pressure[0];
                const int num_ghosts_1_pressure = num_ghosts_pressure[1];
                const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
                
                const int num_ghosts_0_temperature = num_ghosts_temperature[0];
                const int num_ghosts_1_temperature = num_ghosts_temperature[1];
                const int ghostcell_dim_0_temperature = ghostcell_dims_temperature[0];
                
                const double weight = dx[1]/L_y;
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            /*
                             * Compute the index of the data point and count how many times the data is repeated.
                             */
                            
                            const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j);
                            
                            int n_overlapped = 1;
                            
                            for (hier::BoxContainer::BoxContainerConstIterator iob(
                                    patch_overlapped_visible_boxes.begin());
                                 iob != patch_overlapped_visible_boxes.end();
                                 iob++)
                            {
                                const hier::Box& patch_overlapped_visible_box = *iob;
                                
                                if (patch_overlapped_visible_box.contains(idx_pt))
                                {
                                    n_overlapped++;
                                }
                            }
                            
                            /*
                             * Compute the derivative shear stress component.
                             */
                            
                            double value_to_add = double(0);
                            
                            if ((component_idx == 0) && (derivative_direction == 0))
                            {
                                double mu[7];
                                double mu_v[7];
                                
                                std::vector<const double*> Y_ptr;
                                Y_ptr.resize(d_num_species);
                                for (int count = -3; count <= 3; count++)
                                {
                                    const int idx_mass_fraction = (relative_idx_lo_0 + (i + count) + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                    
                                    const int idx_pressure = (relative_idx_lo_0 + (i + count) + num_ghosts_0_pressure) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure;
                                    
                                    const int idx_temperature = (relative_idx_lo_0 + (i + count) + num_ghosts_0_temperature) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_temperature)*ghostcell_dim_0_temperature;
                                    
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        Y_ptr[si] = &Y[si][idx_mass_fraction];
                                    }
                                    mu[count + 3] = d_equation_of_shear_viscosity_mixing_rules->
                                        getShearViscosity(
                                            &p[idx_pressure],
                                            &T[idx_temperature],
                                            Y_ptr);
                                    mu_v[count + 3] = d_equation_of_bulk_viscosity_mixing_rules->
                                        getBulkViscosity(
                                            &p[idx_pressure],
                                            &T[idx_temperature],
                                            Y_ptr);
                                }
                                
                                const double dmudx = (double(1)/double(60)*(mu[6] - mu[0])
                                    - double(3)/double(20)*(mu[5] - mu[1])
                                    + double(3)/double(4)*(mu[4] - mu[2]))/dx[0];
                                
                                const double dmu_vdx = (double(1)/double(60)*(mu_v[6] - mu_v[0])
                                    - double(3)/double(20)*(mu_v[5] - mu_v[1])
                                    + double(3)/double(4)*(mu_v[4] - mu_v[2]))/dx[0];
                                
                                const int idx_vel       = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                
                                const int idx_vel_x_LLL = (relative_idx_lo_0 + (i - 3) + num_ghosts_0_velocity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                
                                const int idx_vel_x_LL  = (relative_idx_lo_0 + (i - 2) + num_ghosts_0_velocity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                
                                const int idx_vel_x_L   = (relative_idx_lo_0 + (i - 1) + num_ghosts_0_velocity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                
                                const int idx_vel_x_R   = (relative_idx_lo_0 + (i + 1) + num_ghosts_0_velocity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                
                                const int idx_vel_x_RR  = (relative_idx_lo_0 + (i + 2) + num_ghosts_0_velocity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                
                                const int idx_vel_x_RRR = (relative_idx_lo_0 + (i + 3) + num_ghosts_0_velocity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                
                                const double dudx = (double(1)/double(60)*(u[idx_vel_x_RRR] - u[idx_vel_x_LLL])
                                    - double(3)/double(20)*(u[idx_vel_x_RR] - u[idx_vel_x_LL])
                                    + double(3)/double(4)*(u[idx_vel_x_R] - u[idx_vel_x_L]))/dx[0];
                                
                                const double d2udx2 = (double(1)/double(90)*(u[idx_vel_x_RRR] + u[idx_vel_x_LLL])
                                    - double(3)/double(20)*(u[idx_vel_x_RR] + u[idx_vel_x_LL])
                                    + double(3)/double(2)*(u[idx_vel_x_R] + u[idx_vel_x_L])
                                    - double(49)/double(18)*u[idx_vel])/(dx[0]*dx[0]);
                                
                                double visc_dvdy[7];
                                
                                for (int count = -3; count <= 3; count++)
                                {
                                    const int idx_vel_y_BBB = (relative_idx_lo_0 + (i + count) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j - 3) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                    
                                    const int idx_vel_y_BB  = (relative_idx_lo_0 + (i + count) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j - 2) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                    
                                    const int idx_vel_y_B   = (relative_idx_lo_0 + (i + count) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j - 1) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                    
                                    const int idx_vel_y_T   = (relative_idx_lo_0 + (i + count) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + 1) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                    
                                    const int idx_vel_y_TT  = (relative_idx_lo_0 + (i + count) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + 2) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                    
                                    const int idx_vel_y_TTT = (relative_idx_lo_0 + (i + count) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + 3) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                    
                                    visc_dvdy[count + 3] = (double(2)/double(3)*mu[count + 3] - mu_v[count + 3])*
                                        (double(1)/double(60)*(v[idx_vel_y_TTT] - v[idx_vel_y_BBB])
                                        - double(3)/double(20)*(v[idx_vel_y_TT] - v[idx_vel_y_BB])
                                        + double(3)/double(4)*(v[idx_vel_y_T] - v[idx_vel_y_B]))/dx[1];
                                }
                                
                                const double dvisc_dvdydx = (double(1)/double(60)*(visc_dvdy[6] - visc_dvdy[0])
                                    - double(3)/double(20)*(visc_dvdy[5] - visc_dvdy[1])
                                    + double(3)/double(4)*(visc_dvdy[4] - visc_dvdy[2]))/dx[0];
                                
                                value_to_add = (
                                    (double(4)/double(3)*dmudx + dmu_vdx)*dudx
                                    + (double(4)/double(3)*mu[3] + mu_v[3])*d2udx2
                                    - dvisc_dvdydx
                                    )*weight/((double) n_overlapped);
                            }
                            else if ((component_idx == 0) && (derivative_direction == 1))
                            {
                                double mu[7];
                                double mu_v[7];
                                
                                std::vector<const double*> Y_ptr;
                                Y_ptr.resize(d_num_species);
                                for (int count = -3; count <= 3; count++)
                                {
                                    const int idx_mass_fraction = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + (j + count) + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                    
                                    const int idx_pressure = (relative_idx_lo_0 + i + num_ghosts_0_pressure) +
                                        (relative_idx_lo_1 + (j + count) + num_ghosts_1_pressure)*ghostcell_dim_0_pressure;
                                    
                                    const int idx_temperature = (relative_idx_lo_0 + i + num_ghosts_0_temperature) +
                                        (relative_idx_lo_1 + (j + count) + num_ghosts_1_temperature)*ghostcell_dim_0_temperature;
                                    
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        Y_ptr[si] = &Y[si][idx_mass_fraction];
                                    }
                                    mu[count + 3] = d_equation_of_shear_viscosity_mixing_rules->
                                        getShearViscosity(
                                            &p[idx_pressure],
                                            &T[idx_temperature],
                                            Y_ptr);
                                    mu_v[count + 3] = d_equation_of_bulk_viscosity_mixing_rules->
                                        getBulkViscosity(
                                            &p[idx_pressure],
                                            &T[idx_temperature],
                                            Y_ptr);
                                }
                                
                                const double dmudy = (double(1)/double(60)*(mu[6] - mu[0])
                                    - double(3)/double(20)*(mu[5] - mu[1])
                                    + double(3)/double(4)*(mu[4] - mu[2]))/dx[1];
                                
                                const double dmu_vdy = (double(1)/double(60)*(mu_v[6] - mu_v[0])
                                    - double(3)/double(20)*(mu_v[5] - mu_v[1])
                                    + double(3)/double(4)*(mu_v[4] - mu_v[2]))/dx[1];
                                
                                const int idx_vel       = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                
                                const int idx_vel_y_BBB = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                    (relative_idx_lo_1 + (j - 3) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                
                                const int idx_vel_y_BB  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                    (relative_idx_lo_1 + (j - 2) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                
                                const int idx_vel_y_B   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                    (relative_idx_lo_1 + (j - 1) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                
                                const int idx_vel_y_T   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                    (relative_idx_lo_1 + (j + 1) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                
                                const int idx_vel_y_TT  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                    (relative_idx_lo_1 + (j + 2) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                
                                const int idx_vel_y_TTT = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                    (relative_idx_lo_1 + (j + 3) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                
                                const double dvdy = (double(1)/double(60)*(v[idx_vel_y_TTT] - v[idx_vel_y_BBB])
                                    - double(3)/double(20)*(v[idx_vel_y_TT] - v[idx_vel_y_BB])
                                    + double(3)/double(4)*(v[idx_vel_y_T] - v[idx_vel_y_B]))/dx[1];
                                
                                const double d2vdy2 = (double(1)/double(90)*(v[idx_vel_y_TTT] + v[idx_vel_y_BBB])
                                    - double(3)/double(20)*(v[idx_vel_y_TT] + v[idx_vel_y_BB])
                                    + double(3)/double(2)*(v[idx_vel_y_T] + v[idx_vel_y_B])
                                    - double(49)/double(18)*v[idx_vel])/(dx[1]*dx[1]);
                                
                                double visc_dudx[7];
                                
                                for (int count = -3; count <= 3; count++)
                                {
                                    const int idx_vel_x_LLL = (relative_idx_lo_0 + (i - 3) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + count) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                    
                                    const int idx_vel_x_LL  = (relative_idx_lo_0 + (i - 2) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + count) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                    
                                    const int idx_vel_x_L   = (relative_idx_lo_0 + (i - 1) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + count) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                    
                                    const int idx_vel_x_R   = (relative_idx_lo_0 + (i + 1) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + count) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                    
                                    const int idx_vel_x_RR  = (relative_idx_lo_0 + (i + 2) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + count) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                    
                                    const int idx_vel_x_RRR = (relative_idx_lo_0 + (i + 3) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + count) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                    
                                    visc_dudx[count + 3] = (double(4)/double(3)*mu[count + 3] + mu_v[count + 3])*
                                        (double(1)/double(60)*(u[idx_vel_x_RRR] - u[idx_vel_x_LLL])
                                        - double(3)/double(20)*(u[idx_vel_x_RR] - u[idx_vel_x_LL])
                                        + double(3)/double(4)*(u[idx_vel_x_R] - u[idx_vel_x_L]))/dx[0];
                                }
                                
                                const double dvisc_dudxdy = (double(1)/double(60)*(visc_dudx[6] - visc_dudx[0])
                                    - double(3)/double(20)*(visc_dudx[5] - visc_dudx[1])
                                    + double(3)/double(4)*(visc_dudx[4] - visc_dudx[2]))/dx[1];
                                
                                value_to_add = (
                                    dvisc_dudxdy
                                    - (double(2)/double(3)*dmudy - dmu_vdy)*dvdy
                                    - (double(2)/double(3)*mu[3] - mu_v[3])*d2vdy2
                                    )*weight/((double) n_overlapped);
                            }
                            else if ((component_idx == 1) && (derivative_direction == 0))
                            {
                                double mu[7];
                                
                                std::vector<const double*> Y_ptr;
                                Y_ptr.resize(d_num_species);
                                for (int count = -3; count <= 3; count++)
                                {
                                    const int idx_mass_fraction = (relative_idx_lo_0 + (i + count) + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                    
                                    const int idx_pressure = (relative_idx_lo_0 + (i + count) + num_ghosts_0_pressure) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure;
                                    
                                    const int idx_temperature = (relative_idx_lo_0 + (i + count) + num_ghosts_0_temperature) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_temperature)*ghostcell_dim_0_temperature;
                                    
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        Y_ptr[si] = &Y[si][idx_mass_fraction];
                                    }
                                    mu[count + 3] = d_equation_of_shear_viscosity_mixing_rules->
                                        getShearViscosity(
                                            &p[idx_pressure],
                                            &T[idx_temperature],
                                            Y_ptr);
                                }
                                
                                const double dmudx = (double(1)/double(60)*(mu[6] - mu[0])
                                    - double(3)/double(20)*(mu[5] - mu[1])
                                    + double(3)/double(4)*(mu[4] - mu[2]))/dx[0];
                                
                                const int idx_vel       = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                
                                const int idx_vel_x_LLL = (relative_idx_lo_0 + (i - 3) + num_ghosts_0_velocity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                
                                const int idx_vel_x_LL  = (relative_idx_lo_0 + (i - 2) + num_ghosts_0_velocity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                
                                const int idx_vel_x_L   = (relative_idx_lo_0 + (i - 1) + num_ghosts_0_velocity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                
                                const int idx_vel_x_R   = (relative_idx_lo_0 + (i + 1) + num_ghosts_0_velocity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                
                                const int idx_vel_x_RR  = (relative_idx_lo_0 + (i + 2) + num_ghosts_0_velocity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                
                                const int idx_vel_x_RRR = (relative_idx_lo_0 + (i + 3) + num_ghosts_0_velocity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                
                                const double dvdx = (double(1)/double(60)*(v[idx_vel_x_RRR] - v[idx_vel_x_LLL])
                                    - double(3)/double(20)*(v[idx_vel_x_RR] - v[idx_vel_x_LL])
                                    + double(3)/double(4)*(v[idx_vel_x_R] - v[idx_vel_x_L]))/dx[0];
                                
                                const double d2vdx2 = (double(1)/double(90)*(v[idx_vel_x_RRR] + v[idx_vel_x_LLL])
                                    - double(3)/double(20)*(v[idx_vel_x_RR] + v[idx_vel_x_LL])
                                    + double(3)/double(2)*(v[idx_vel_x_R] + v[idx_vel_x_L])
                                    - double(49)/double(18)*v[idx_vel])/(dx[0]*dx[0]);
                                
                                double visc_dudy[7];
                                
                                for (int count = -3; count <= 3; count++)
                                {
                                    const int idx_vel_y_BBB = (relative_idx_lo_0 + (i + count) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j - 3) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                    
                                    const int idx_vel_y_BB  = (relative_idx_lo_0 + (i + count) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j - 2) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                    
                                    const int idx_vel_y_B   = (relative_idx_lo_0 + (i + count) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j - 1) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                    
                                    const int idx_vel_y_T   = (relative_idx_lo_0 + (i + count) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + 1) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                    
                                    const int idx_vel_y_TT  = (relative_idx_lo_0 + (i + count) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + 2) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                    
                                    const int idx_vel_y_TTT = (relative_idx_lo_0 + (i + count) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + 3) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                    
                                    visc_dudy[count + 3] = mu[count + 3]*
                                        (double(1)/double(60)*(u[idx_vel_y_TTT] - u[idx_vel_y_BBB])
                                        - double(3)/double(20)*(u[idx_vel_y_TT] - u[idx_vel_y_BB])
                                        + double(3)/double(4)*(u[idx_vel_y_T] - u[idx_vel_y_B]))/dx[1];
                                }
                                
                                const double dvisc_dudydx = (double(1)/double(60)*(visc_dudy[6] - visc_dudy[0])
                                    - double(3)/double(20)*(visc_dudy[5] - visc_dudy[1])
                                    + double(3)/double(4)*(visc_dudy[4] - visc_dudy[2]))/dx[0];
                                
                                value_to_add = (
                                    dvisc_dudydx
                                    + dmudx*dvdx
                                    + mu[3]*d2vdx2
                                    )*weight/((double) n_overlapped);
                            }
                            else if ((component_idx == 1) && (derivative_direction == 1))
                            {
                                double mu[7];
                                
                                std::vector<const double*> Y_ptr;
                                Y_ptr.resize(d_num_species);
                                for (int count = -3; count <= 3; count++)
                                {
                                    const int idx_mass_fraction = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + (j + count) + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                    
                                    const int idx_pressure = (relative_idx_lo_0 + i + num_ghosts_0_pressure) +
                                        (relative_idx_lo_1 + (j + count) + num_ghosts_1_pressure)*ghostcell_dim_0_pressure;
                                    
                                    const int idx_temperature = (relative_idx_lo_0 + i + num_ghosts_0_temperature) +
                                        (relative_idx_lo_1 + (j + count) + num_ghosts_1_temperature)*ghostcell_dim_0_temperature;
                                    
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        Y_ptr[si] = &Y[si][idx_mass_fraction];
                                    }
                                    mu[count + 3] = d_equation_of_shear_viscosity_mixing_rules->
                                        getShearViscosity(
                                            &p[idx_pressure],
                                            &T[idx_temperature],
                                            Y_ptr);
                                }
                                
                                const double dmudy = (double(1)/double(60)*(mu[6] - mu[0])
                                    - double(3)/double(20)*(mu[5] - mu[1])
                                    + double(3)/double(4)*(mu[4] - mu[2]))/dx[1];
                                
                                const int idx_vel       = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                
                                const int idx_vel_y_BBB = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                    (relative_idx_lo_1 + (j - 3) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                
                                const int idx_vel_y_BB  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                    (relative_idx_lo_1 + (j - 2) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                
                                const int idx_vel_y_B   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                    (relative_idx_lo_1 + (j - 1) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                
                                const int idx_vel_y_T   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                    (relative_idx_lo_1 + (j + 1) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                
                                const int idx_vel_y_TT  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                    (relative_idx_lo_1 + (j + 2) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                
                                const int idx_vel_y_TTT = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                    (relative_idx_lo_1 + (j + 3) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                
                                const double dudy = (double(1)/double(60)*(u[idx_vel_y_TTT] - u[idx_vel_y_BBB])
                                    - double(3)/double(20)*(u[idx_vel_y_TT] - u[idx_vel_y_BB])
                                    + double(3)/double(4)*(u[idx_vel_y_T] - u[idx_vel_y_B]))/dx[1];
                                
                                const double d2udy2 = (double(1)/double(90)*(u[idx_vel_y_TTT] + u[idx_vel_y_BBB])
                                    - double(3)/double(20)*(u[idx_vel_y_TT] + u[idx_vel_y_BB])
                                    + double(3)/double(2)*(u[idx_vel_y_T] + u[idx_vel_y_B])
                                    - double(49)/double(18)*u[idx_vel])/(dx[1]*dx[1]);
                                
                                double visc_dvdx[7];
                                
                                for (int count = -3; count <= 3; count++)
                                {
                                    const int idx_vel_x_LLL = (relative_idx_lo_0 + (i - 3) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + count) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                    
                                    const int idx_vel_x_LL  = (relative_idx_lo_0 + (i - 2) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + count) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                    
                                    const int idx_vel_x_L   = (relative_idx_lo_0 + (i - 1) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + count) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                    
                                    const int idx_vel_x_R   = (relative_idx_lo_0 + (i + 1) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + count) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                    
                                    const int idx_vel_x_RR  = (relative_idx_lo_0 + (i + 2) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + count) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                    
                                    const int idx_vel_x_RRR = (relative_idx_lo_0 + (i + 3) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + count) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                    
                                    visc_dvdx[count + 3] = mu[count + 3]*
                                        (double(1)/double(60)*(v[idx_vel_x_RRR] - v[idx_vel_x_LLL])
                                        - double(3)/double(20)*(v[idx_vel_x_RR] - v[idx_vel_x_LL])
                                        + double(3)/double(4)*(v[idx_vel_x_R] - v[idx_vel_x_L]))/dx[0];
                                }
                                
                                const double dvisc_dvdxdy = (double(1)/double(60)*(visc_dvdx[6] - visc_dvdx[0])
                                    - double(3)/double(20)*(visc_dvdx[5] - visc_dvdx[1])
                                    + double(3)/double(4)*(visc_dvdx[4] - visc_dvdx[2]))/dx[1];
                                
                                value_to_add = (
                                    dvisc_dvdxdy
                                    + dmudy*dudy
                                    + mu[3]*d2udy2
                                    )*weight/((double) n_overlapped);
                            }
                            else
                            {
                                TBOX_ERROR(d_object_name
                                    << ": "
                                    << "Cannot compute derivative of shear stress component for two-dimensional problem!\n"
                                    << "component_idx = " << component_idx << " and derivative_direction = "
                                    << derivative_direction << " given!\n"
                                    << std::endl);
                            }
                            
                            /*
                             * Add the data.
                             */
                            
                            for (int ii = 0; ii < ratio_to_finest_level_0; ii++)
                            {
                                const int idx_fine = (idx_lo_0 + i)*ratio_to_finest_level_0 + ii;
                                
                                tau_ij_der_avg_local[idx_fine] += value_to_add;
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average of derivative of shear stress component.
         */
        
        mpi.Allreduce(
            tau_ij_der_avg_local,
            tau_ij_der_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(tau_ij_der_avg_local);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        hier::IntVector num_ghosts_der = hier::IntVector::getOne(d_dim)*d_num_ghosts_derivative;
        
       const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        const double L_z = x_hi[2] - x_lo[2];
        
        double* tau_ij_der_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_derivative_tau_ij.resize(finest_level_dim_0);
        double* tau_ij_der_avg_global = averaged_derivative_tau_ij.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            tau_ij_der_avg_local[i] = double(0);
            tau_ij_der_avg_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratio_to_coarsest_level =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarsest_level *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = ratio_finest_level_to_coarsest_level/ratio_to_coarsest_level;
            
            const int ratio_to_finest_level_0 = ratio_to_finest_level[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and the quantities in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", num_ghosts_der));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTIONS", num_ghosts_der));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("PRESSURE", num_ghosts_der));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("TEMPERATURE", num_ghosts_der));
                
                d_flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                d_flow_model_tmp->allocateMemoryForDerivedCellData();
                
                d_flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_velocity =
                    d_flow_model_tmp->getCellData("VELOCITY");
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_mass_fraction =
                    d_flow_model_tmp->getCellData("MASS_FRACTIONS");
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_pressure =
                    d_flow_model_tmp->getCellData("PRESSURE");
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_temperature =
                    d_flow_model_tmp->getCellData("TEMPERATURE");
                
                double* u = data_velocity->getPointer(0);
                double* v = data_velocity->getPointer(1);
                double* w = data_velocity->getPointer(2);
                std::vector<double*> Y;
                Y.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    Y.push_back(data_mass_fraction->getPointer(si));
                }
                double* p = data_pressure->getPointer(0);
                double* T = data_temperature->getPointer(0);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_velocity = data_velocity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_velocity = data_velocity->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_mass_fraction = data_mass_fraction->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mass_fraction = data_mass_fraction->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_pressure = data_pressure->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_temperature = data_temperature->getGhostBox().numberCells();
                
                const int num_ghosts_0_velocity = num_ghosts_velocity[0];
                const int num_ghosts_1_velocity = num_ghosts_velocity[1];
                const int num_ghosts_2_velocity = num_ghosts_velocity[2];
                const int ghostcell_dim_0_velocity = ghostcell_dims_velocity[0];
                const int ghostcell_dim_1_velocity = ghostcell_dims_velocity[1];
                
                const int num_ghosts_0_mass_fraction = num_ghosts_mass_fraction[0];
                const int num_ghosts_1_mass_fraction = num_ghosts_mass_fraction[1];
                const int num_ghosts_2_mass_fraction = num_ghosts_mass_fraction[2];
                const int ghostcell_dim_0_mass_fraction = ghostcell_dims_mass_fraction[0];
                const int ghostcell_dim_1_mass_fraction = ghostcell_dims_mass_fraction[1];
                
                const int num_ghosts_0_pressure = num_ghosts_pressure[0];
                const int num_ghosts_1_pressure = num_ghosts_pressure[1];
                const int num_ghosts_2_pressure = num_ghosts_pressure[2];
                const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
                const int ghostcell_dim_1_pressure = ghostcell_dims_pressure[1];
                
                const int num_ghosts_0_temperature = num_ghosts_temperature[0];
                const int num_ghosts_1_temperature = num_ghosts_temperature[1];
                const int num_ghosts_2_temperature = num_ghosts_temperature[2];
                const int ghostcell_dim_0_temperature = ghostcell_dims_temperature[0];
                const int ghostcell_dim_1_temperature = ghostcell_dims_temperature[1];
                
                const double weight = (dx[1]*dx[2])/(L_y*L_z);
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    const int interior_dim_2 = interior_dims[2];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int idx_lo_2 = index_lo[2];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    const int relative_idx_lo_2 = relative_index_lo[2];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                /*
                                 * Compute the index of the data point and count how many times the data is repeated.
                                 */
                                
                                const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j, idx_lo_2 + k);
                                
                                int n_overlapped = 1;
                                
                                for (hier::BoxContainer::BoxContainerConstIterator iob(
                                        patch_overlapped_visible_boxes.begin());
                                     iob != patch_overlapped_visible_boxes.end();
                                     iob++)
                                {
                                    const hier::Box& patch_overlapped_visible_box = *iob;
                                    
                                    if (patch_overlapped_visible_box.contains(idx_pt))
                                    {
                                        n_overlapped++;
                                    }
                                }
                                
                                /*
                                 * Compute the derivative of shear stress component.
                                 */
                                
                                double value_to_add = double(0);
                                
                                if ((component_idx == 0) && (derivative_direction == 0))
                                {
                                    double mu[7];
                                    double mu_v[7];
                                    
                                    std::vector<const double*> Y_ptr;
                                    Y_ptr.resize(d_num_species);
                                    for (int count = -3; count <= 3; count++)
                                    {
                                        const int idx_mass_fraction = (relative_idx_lo_0 + (i + count) + num_ghosts_0_mass_fraction) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                            (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                                ghostcell_dim_1_mass_fraction;
                                        
                                        const int idx_pressure = (relative_idx_lo_0 + (i + count) + num_ghosts_0_pressure) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure +
                                            (relative_idx_lo_2 + k + num_ghosts_2_pressure)*ghostcell_dim_0_pressure*
                                                ghostcell_dim_1_pressure;
                                        
                                        const int idx_temperature = (relative_idx_lo_0 + (i + count) + num_ghosts_0_temperature) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_temperature)*ghostcell_dim_0_temperature +
                                            (relative_idx_lo_2 + k + num_ghosts_2_temperature)*ghostcell_dim_0_temperature*
                                                ghostcell_dim_1_temperature;
                                        
                                        for (int si = 0; si < d_num_species; si++)
                                        {
                                            Y_ptr[si] = &Y[si][idx_mass_fraction];
                                        }
                                        mu[count + 3] = d_equation_of_shear_viscosity_mixing_rules->
                                            getShearViscosity(
                                                &p[idx_pressure],
                                                &T[idx_temperature],
                                                Y_ptr);
                                        mu_v[count + 3] = d_equation_of_bulk_viscosity_mixing_rules->
                                            getBulkViscosity(
                                                &p[idx_pressure],
                                                &T[idx_temperature],
                                                Y_ptr);
                                    }
                                    
                                    const double dmudx = (double(1)/double(60)*(mu[6] - mu[0])
                                        - double(3)/double(20)*(mu[5] - mu[1])
                                        + double(3)/double(4)*(mu[4] - mu[2]))/dx[0];
                                    
                                    const double dmu_vdx = (double(1)/double(60)*(mu_v[6] - mu_v[0])
                                        - double(3)/double(20)*(mu_v[5] - mu_v[1])
                                        + double(3)/double(4)*(mu_v[4] - mu_v[2]))/dx[0];
                                    
                                    const int idx_vel       = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_LLL = (relative_idx_lo_0 + (i - 3) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_LL  = (relative_idx_lo_0 + (i - 2) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_L   = (relative_idx_lo_0 + (i - 1) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_R   = (relative_idx_lo_0 + (i + 1) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_RR  = (relative_idx_lo_0 + (i + 2) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_RRR = (relative_idx_lo_0 + (i + 3) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const double dudx = (double(1)/double(60)*(u[idx_vel_x_RRR] - u[idx_vel_x_LLL])
                                        - double(3)/double(20)*(u[idx_vel_x_RR] - u[idx_vel_x_LL])
                                        + double(3)/double(4)*(u[idx_vel_x_R] - u[idx_vel_x_L]))/dx[0];
                                    
                                    const double d2udx2 = (double(1)/double(90)*(u[idx_vel_x_RRR] + u[idx_vel_x_LLL])
                                        - double(3)/double(20)*(u[idx_vel_x_RR] + u[idx_vel_x_LL])
                                        + double(3)/double(2)*(u[idx_vel_x_R] + u[idx_vel_x_L])
                                        - double(49)/double(18)*u[idx_vel])/(dx[0]*dx[0]);
                                    
                                    double visc_dvdy[7];
                                    
                                    for (int count = -3; count <= 3; count++)
                                    {
                                        const int idx_vel_y_BBB = (relative_idx_lo_0 + (i + count) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + (j - 3) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_y_BB  = (relative_idx_lo_0 + (i + count) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + (j - 2) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_y_B   = (relative_idx_lo_0 + (i + count) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + (j - 1) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_y_T   = (relative_idx_lo_0 + (i + count) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + (j + 1) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_y_TT  = (relative_idx_lo_0 + (i + count) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + (j + 2) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_y_TTT = (relative_idx_lo_0 + (i + count) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + (j + 3) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        visc_dvdy[count + 3] = (double(2)/double(3)*mu[count + 3] - mu_v[count + 3])*
                                            (double(1)/double(60)*(v[idx_vel_y_TTT] - v[idx_vel_y_BBB])
                                            - double(3)/double(20)*(v[idx_vel_y_TT] - v[idx_vel_y_BB])
                                            + double(3)/double(4)*(v[idx_vel_y_T] - v[idx_vel_y_B]))/dx[1];
                                    }
                                    
                                    const double dvisc_dvdydx = (double(1)/double(60)*(visc_dvdy[6] - visc_dvdy[0])
                                        - double(3)/double(20)*(visc_dvdy[5] - visc_dvdy[1])
                                        + double(3)/double(4)*(visc_dvdy[4] - visc_dvdy[2]))/dx[0];
                                    
                                    double visc_dwdz[7];
                                    
                                    for (int count = -3; count <= 3; count++)
                                    {
                                        const int idx_vel_z_BBB = (relative_idx_lo_0 + (i + count) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + (k - 3) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_z_BB  = (relative_idx_lo_0 + (i + count) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + (k - 2) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_z_B   = (relative_idx_lo_0 + (i + count) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + (k - 1) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_z_F   = (relative_idx_lo_0 + (i + count) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + (k + 1) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_z_FF  = (relative_idx_lo_0 + (i + count) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + (k + 2) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_z_FFF = (relative_idx_lo_0 + (i + count) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + (k + 3) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        visc_dwdz[count + 3] = (double(2)/double(3)*mu[count + 3] - mu_v[count + 3])*
                                            (double(1)/double(60)*(w[idx_vel_z_FFF] - w[idx_vel_z_BBB])
                                            - double(3)/double(20)*(w[idx_vel_z_FF] - w[idx_vel_z_BB])
                                            + double(3)/double(4)*(w[idx_vel_z_F] - w[idx_vel_z_B]))/dx[2];
                                    }
                                    
                                    const double dvisc_dwdzdx = (double(1)/double(60)*(visc_dwdz[6] - visc_dwdz[0])
                                        - double(3)/double(20)*(visc_dwdz[5] - visc_dwdz[1])
                                        + double(3)/double(4)*(visc_dwdz[4] - visc_dwdz[2]))/dx[0];
                                    
                                    value_to_add = (
                                        (double(4)/double(3)*dmudx + dmu_vdx)*dudx
                                        + (double(4)/double(3)*mu[3] + mu_v[3])*d2udx2
                                        - dvisc_dvdydx
                                        - dvisc_dwdzdx
                                        )*weight/((double) n_overlapped);
                                }
                                else if ((component_idx == 0) && (derivative_direction == 1))
                                {
                                    double mu[7];
                                    double mu_v[7];
                                    
                                    std::vector<const double*> Y_ptr;
                                    Y_ptr.resize(d_num_species);
                                    for (int count = -3; count <= 3; count++)
                                    {
                                        const int idx_mass_fraction = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                            (relative_idx_lo_1 + (j + count) + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                            (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                                ghostcell_dim_1_mass_fraction;
                                        
                                        const int idx_pressure = (relative_idx_lo_0 + i + num_ghosts_0_pressure) +
                                            (relative_idx_lo_1 + (j + count) + num_ghosts_1_pressure)*ghostcell_dim_0_pressure +
                                            (relative_idx_lo_2 + k + num_ghosts_2_pressure)*ghostcell_dim_0_pressure*
                                                ghostcell_dim_1_pressure;
                                        
                                        const int idx_temperature = (relative_idx_lo_0 + i + num_ghosts_0_temperature) +
                                            (relative_idx_lo_1 + (j + count) + num_ghosts_1_temperature)*ghostcell_dim_0_temperature +
                                            (relative_idx_lo_2 + k + num_ghosts_2_temperature)*ghostcell_dim_0_temperature*
                                                ghostcell_dim_1_temperature;
                                        
                                        for (int si = 0; si < d_num_species; si++)
                                        {
                                            Y_ptr[si] = &Y[si][idx_mass_fraction];
                                        }
                                        mu[count + 3] = d_equation_of_shear_viscosity_mixing_rules->
                                            getShearViscosity(
                                                &p[idx_pressure],
                                                &T[idx_temperature],
                                                Y_ptr);
                                        mu_v[count + 3] = d_equation_of_bulk_viscosity_mixing_rules->
                                            getBulkViscosity(
                                                &p[idx_pressure],
                                                &T[idx_temperature],
                                                Y_ptr);
                                    }
                                    
                                    const double dmudy = (double(1)/double(60)*(mu[6] - mu[0])
                                        - double(3)/double(20)*(mu[5] - mu[1])
                                        + double(3)/double(4)*(mu[4] - mu[2]))/dx[1];
                                    
                                    const double dmu_vdy = (double(1)/double(60)*(mu_v[6] - mu_v[0])
                                        - double(3)/double(20)*(mu_v[5] - mu_v[1])
                                        + double(3)/double(4)*(mu_v[4] - mu_v[2]))/dx[1];
                                    
                                    const int idx_vel       = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_BBB = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j - 3) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_BB  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j - 2) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_B   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j - 1) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_T   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + 1) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_TT  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + 2) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_TTT = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + 3) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const double dvdy = (double(1)/double(60)*(v[idx_vel_y_TTT] - v[idx_vel_y_BBB])
                                        - double(3)/double(20)*(v[idx_vel_y_TT] - v[idx_vel_y_BB])
                                        + double(3)/double(4)*(v[idx_vel_y_T] - v[idx_vel_y_B]))/dx[1];
                                    
                                    const double d2vdy2 = (double(1)/double(90)*(v[idx_vel_y_TTT] + v[idx_vel_y_BBB])
                                        - double(3)/double(20)*(v[idx_vel_y_TT] + v[idx_vel_y_BB])
                                        + double(3)/double(2)*(v[idx_vel_y_T] + v[idx_vel_y_B])
                                        - double(49)/double(18)*v[idx_vel])/(dx[1]*dx[1]);
                                    
                                    double visc_dudx[7];
                                    
                                    for (int count = -3; count <= 3; count++)
                                    {
                                        const int idx_vel_x_LLL = (relative_idx_lo_0 + (i - 3) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + (j + count) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_x_LL  = (relative_idx_lo_0 + (i - 2) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + (j + count) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_x_L   = (relative_idx_lo_0 + (i - 1) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + (j + count) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_x_R   = (relative_idx_lo_0 + (i + 1) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + (j + count) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_x_RR  = (relative_idx_lo_0 + (i + 2) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + (j + count) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_x_RRR = (relative_idx_lo_0 + (i + 3) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + (j + count) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        visc_dudx[count + 3] = (double(4)/double(3)*mu[count + 3] + mu_v[count + 3])*
                                            (double(1)/double(60)*(u[idx_vel_x_RRR] - u[idx_vel_x_LLL])
                                            - double(3)/double(20)*(u[idx_vel_x_RR] - u[idx_vel_x_LL])
                                            + double(3)/double(4)*(u[idx_vel_x_R] - u[idx_vel_x_L]))/dx[0];
                                    }
                                    
                                    const double dvisc_dudxdy = (double(1)/double(60)*(visc_dudx[6] - visc_dudx[0])
                                        - double(3)/double(20)*(visc_dudx[5] - visc_dudx[1])
                                        + double(3)/double(4)*(visc_dudx[4] - visc_dudx[2]))/dx[1];
                                    
                                    double visc_dwdz[7];
                                    
                                    for (int count = -3; count <= 3; count++)
                                    {
                                        const int idx_vel_z_BBB = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + (j + count) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + (k - 3) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_z_BB  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + (j + count) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + (k - 2) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_z_B   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + (j + count) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + (k - 1) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_z_F   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + (j + count) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + (k + 1) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_z_FF  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + (j + count) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + (k + 2) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_z_FFF = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + (j + count) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + (k + 3) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        visc_dwdz[count + 3] = (double(2)/double(3)*mu[count + 3] - mu_v[count + 3])*
                                            (double(1)/double(60)*(w[idx_vel_z_FFF] - w[idx_vel_z_BBB])
                                            - double(3)/double(20)*(w[idx_vel_z_FF] - w[idx_vel_z_BB])
                                            + double(3)/double(4)*(w[idx_vel_z_F] - w[idx_vel_z_B]))/dx[2];
                                    }
                                    
                                    const double dvisc_dwdzdy = (double(1)/double(60)*(visc_dwdz[6] - visc_dwdz[0])
                                        - double(3)/double(20)*(visc_dwdz[5] - visc_dwdz[1])
                                        + double(3)/double(4)*(visc_dwdz[4] - visc_dwdz[2]))/dx[1];
                                    
                                    value_to_add = (
                                        dvisc_dudxdy
                                        - dvisc_dwdzdy
                                        - (double(2)/double(3)*dmudy - dmu_vdy)*dvdy
                                        - (double(2)/double(3)*mu[3] - mu_v[3])*d2vdy2
                                        )*weight/((double) n_overlapped);
                                }
                                else if ((component_idx == 0) && (derivative_direction == 2))
                                {
                                    double mu[7];
                                    double mu_v[7];
                                    
                                    std::vector<const double*> Y_ptr;
                                    Y_ptr.resize(d_num_species);
                                    for (int count = -3; count <= 3; count++)
                                    {
                                        const int idx_mass_fraction = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                            (relative_idx_lo_2 + (k + count) + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                                ghostcell_dim_1_mass_fraction;
                                        
                                        const int idx_pressure = (relative_idx_lo_0 + i + num_ghosts_0_pressure) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure +
                                            (relative_idx_lo_2 + (k + count) + num_ghosts_2_pressure)*ghostcell_dim_0_pressure*
                                                ghostcell_dim_1_pressure;
                                        
                                        const int idx_temperature = (relative_idx_lo_0 + i + num_ghosts_0_temperature) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_temperature)*ghostcell_dim_0_temperature +
                                            (relative_idx_lo_2 + (k + count) + num_ghosts_2_temperature)*ghostcell_dim_0_temperature*
                                                ghostcell_dim_1_temperature;
                                        
                                        for (int si = 0; si < d_num_species; si++)
                                        {
                                            Y_ptr[si] = &Y[si][idx_mass_fraction];
                                        }
                                        mu[count + 3] = d_equation_of_shear_viscosity_mixing_rules->
                                            getShearViscosity(
                                                &p[idx_pressure],
                                                &T[idx_temperature],
                                                Y_ptr);
                                        mu_v[count + 3] = d_equation_of_bulk_viscosity_mixing_rules->
                                            getBulkViscosity(
                                                &p[idx_pressure],
                                                &T[idx_temperature],
                                                Y_ptr);
                                    }
                                    
                                    const double dmudz = (double(1)/double(60)*(mu[6] - mu[0])
                                        - double(3)/double(20)*(mu[5] - mu[1])
                                        + double(3)/double(4)*(mu[4] - mu[2]))/dx[2];
                                    
                                    const double dmu_vdz = (double(1)/double(60)*(mu_v[6] - mu_v[0])
                                        - double(3)/double(20)*(mu_v[5] - mu_v[1])
                                        + double(3)/double(4)*(mu_v[4] - mu_v[2]))/dx[2];
                                    
                                    const int idx_vel       = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_BBB = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k - 3) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_BB  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k - 2) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_B   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k - 1) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_F   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k + 1) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_FF  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k + 2) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_FFF = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k + 3) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const double dwdz = (double(1)/double(60)*(w[idx_vel_z_FFF] - w[idx_vel_z_BBB])
                                        - double(3)/double(20)*(w[idx_vel_z_FF] - w[idx_vel_z_BB])
                                        + double(3)/double(4)*(w[idx_vel_z_F] - w[idx_vel_z_B]))/dx[2];
                                    
                                    const double d2wdz2 = (double(1)/double(90)*(w[idx_vel_z_FFF] + w[idx_vel_z_BBB])
                                        - double(3)/double(20)*(w[idx_vel_z_FF] + w[idx_vel_z_BB])
                                        + double(3)/double(2)*(w[idx_vel_z_F] + w[idx_vel_z_B])
                                        - double(49)/double(18)*w[idx_vel])/(dx[2]*dx[2]);
                                    
                                    double visc_dudx[7];
                                    
                                    for (int count = -3; count <= 3; count++)
                                    {
                                        const int idx_vel_x_LLL = (relative_idx_lo_0 + (i - 3) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + (k + count) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_x_LL  = (relative_idx_lo_0 + (i - 2) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + (k + count) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_x_L   = (relative_idx_lo_0 + (i - 1) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + (k + count) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_x_R   = (relative_idx_lo_0 + (i + 1) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + (k + count) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_x_RR  = (relative_idx_lo_0 + (i + 2) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + (k + count) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_x_RRR = (relative_idx_lo_0 + (i + 3) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + (k + count) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        visc_dudx[count + 3] = (double(4)/double(3)*mu[count + 3] + mu_v[count + 3])*
                                            (double(1)/double(60)*(u[idx_vel_x_RRR] - u[idx_vel_x_LLL])
                                            - double(3)/double(20)*(u[idx_vel_x_RR] - u[idx_vel_x_LL])
                                            + double(3)/double(4)*(u[idx_vel_x_R] - u[idx_vel_x_L]))/dx[0];
                                    }
                                    
                                    const double dvisc_dudxdz = (double(1)/double(60)*(visc_dudx[6] - visc_dudx[0])
                                        - double(3)/double(20)*(visc_dudx[5] - visc_dudx[1])
                                        + double(3)/double(4)*(visc_dudx[4] - visc_dudx[2]))/dx[2];
                                    
                                    double visc_dvdy[7];
                                    
                                    for (int count = -3; count <= 3; count++)
                                    {
                                        const int idx_vel_y_BBB = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + (j - 3) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + (k + count) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_y_BB  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + (j - 2) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + (k + count) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_y_B   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + (j - 1) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + (k + count) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_y_T   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + (j + 1) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + (k + count) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_y_TT  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + (j + 2) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + (k + count) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_y_TTT = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + (j + 3) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + (k + count) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        visc_dvdy[count + 3] = (double(2)/double(3)*mu[count + 3] - mu_v[count + 3])*
                                            (double(1)/double(60)*(v[idx_vel_y_TTT] - v[idx_vel_y_BBB])
                                            - double(3)/double(20)*(v[idx_vel_y_TT] - v[idx_vel_y_BB])
                                            + double(3)/double(4)*(v[idx_vel_y_T] - v[idx_vel_y_B]))/dx[1];
                                    }
                                    
                                    const double dvisc_dvdydz = (double(1)/double(60)*(visc_dvdy[6] - visc_dvdy[0])
                                        - double(3)/double(20)*(visc_dvdy[5] - visc_dvdy[1])
                                        + double(3)/double(4)*(visc_dvdy[4] - visc_dvdy[2]))/dx[2];
                                    
                                    value_to_add = (
                                        dvisc_dudxdz
                                        - dvisc_dvdydz
                                        - (double(2)/double(3)*dmudz - dmu_vdz)*dwdz
                                        - (double(2)/double(3)*mu[3] - mu_v[3])*d2wdz2
                                        )*weight/((double) n_overlapped);
                                }
                                else if ((component_idx == 1) && (derivative_direction == 0))
                                {
                                    double mu[7];
                                    
                                    std::vector<const double*> Y_ptr;
                                    Y_ptr.resize(d_num_species);
                                    for (int count = -3; count <= 3; count++)
                                    {
                                        const int idx_mass_fraction = (relative_idx_lo_0 + (i + count) + num_ghosts_0_mass_fraction) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                            (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                                ghostcell_dim_1_mass_fraction;
                                        
                                        const int idx_pressure = (relative_idx_lo_0 + (i + count) + num_ghosts_0_pressure) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure +
                                            (relative_idx_lo_2 + k + num_ghosts_2_pressure)*ghostcell_dim_0_pressure*
                                                ghostcell_dim_1_pressure;
                                        
                                        const int idx_temperature = (relative_idx_lo_0 + (i + count) + num_ghosts_0_temperature) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_temperature)*ghostcell_dim_0_temperature +
                                            (relative_idx_lo_2 + k + num_ghosts_2_temperature)*ghostcell_dim_0_temperature*
                                                ghostcell_dim_1_temperature;
                                        
                                        for (int si = 0; si < d_num_species; si++)
                                        {
                                            Y_ptr[si] = &Y[si][idx_mass_fraction];
                                        }
                                        mu[count + 3] = d_equation_of_shear_viscosity_mixing_rules->
                                            getShearViscosity(
                                                &p[idx_pressure],
                                                &T[idx_temperature],
                                                Y_ptr);
                                    }
                                    
                                    const double dmudx = (double(1)/double(60)*(mu[6] - mu[0])
                                        - double(3)/double(20)*(mu[5] - mu[1])
                                        + double(3)/double(4)*(mu[4] - mu[2]))/dx[0];
                                    
                                    const int idx_vel       = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_LLL = (relative_idx_lo_0 + (i - 3) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_LL  = (relative_idx_lo_0 + (i - 2) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_L   = (relative_idx_lo_0 + (i - 1) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_R   = (relative_idx_lo_0 + (i + 1) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_RR  = (relative_idx_lo_0 + (i + 2) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_RRR = (relative_idx_lo_0 + (i + 3) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const double dvdx = (double(1)/double(60)*(v[idx_vel_x_RRR] - v[idx_vel_x_LLL])
                                        - double(3)/double(20)*(v[idx_vel_x_RR] - v[idx_vel_x_LL])
                                        + double(3)/double(4)*(v[idx_vel_x_R] - v[idx_vel_x_L]))/dx[0];
                                    
                                    const double d2vdx2 = (double(1)/double(90)*(v[idx_vel_x_RRR] + v[idx_vel_x_LLL])
                                        - double(3)/double(20)*(v[idx_vel_x_RR] + v[idx_vel_x_LL])
                                        + double(3)/double(2)*(v[idx_vel_x_R] + v[idx_vel_x_L])
                                        - double(49)/double(18)*v[idx_vel])/(dx[0]*dx[0]);
                                    
                                    double visc_dudy[7];
                                    
                                    for (int count = -3; count <= 3; count++)
                                    {
                                        const int idx_vel_y_BBB = (relative_idx_lo_0 + (i + count) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + (j - 3) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_y_BB  = (relative_idx_lo_0 + (i + count) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + (j - 2) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_y_B   = (relative_idx_lo_0 + (i + count) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + (j - 1) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_y_T   = (relative_idx_lo_0 + (i + count) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + (j + 1) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_y_TT  = (relative_idx_lo_0 + (i + count) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + (j + 2) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_y_TTT = (relative_idx_lo_0 + (i + count) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + (j + 3) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        visc_dudy[count + 3] = mu[count + 3]*
                                            (double(1)/double(60)*(u[idx_vel_y_TTT] - u[idx_vel_y_BBB])
                                            - double(3)/double(20)*(u[idx_vel_y_TT] - u[idx_vel_y_BB])
                                            + double(3)/double(4)*(u[idx_vel_y_T] - u[idx_vel_y_B]))/dx[1];
                                    }
                                    
                                    const double dvisc_dudydx = (double(1)/double(60)*(visc_dudy[6] - visc_dudy[0])
                                        - double(3)/double(20)*(visc_dudy[5] - visc_dudy[1])
                                        + double(3)/double(4)*(visc_dudy[4] - visc_dudy[2]))/dx[0];
                                    
                                    value_to_add = (
                                        dmudx*dvdx
                                        + mu[3]*d2vdx2
                                        + dvisc_dudydx
                                        )*weight/((double) n_overlapped);
                                }
                                else if ((component_idx == 1) && (derivative_direction == 1))
                                {
                                    double mu[7];
                                    
                                    std::vector<const double*> Y_ptr;
                                    Y_ptr.resize(d_num_species);
                                    for (int count = -3; count <= 3; count++)
                                    {
                                        const int idx_mass_fraction = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                            (relative_idx_lo_1 + (j + count) + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                            (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                                ghostcell_dim_1_mass_fraction;
                                        
                                        const int idx_pressure = (relative_idx_lo_0 + i + num_ghosts_0_pressure) +
                                            (relative_idx_lo_1 + (j + count) + num_ghosts_1_pressure)*ghostcell_dim_0_pressure +
                                            (relative_idx_lo_2 + k + num_ghosts_2_pressure)*ghostcell_dim_0_pressure*
                                                ghostcell_dim_1_pressure;
                                        
                                        const int idx_temperature = (relative_idx_lo_0 + i + num_ghosts_0_temperature) +
                                            (relative_idx_lo_1 + (j + count) + num_ghosts_1_temperature)*ghostcell_dim_0_temperature +
                                            (relative_idx_lo_2 + k + num_ghosts_2_temperature)*ghostcell_dim_0_temperature*
                                                ghostcell_dim_1_temperature;
                                        
                                        for (int si = 0; si < d_num_species; si++)
                                        {
                                            Y_ptr[si] = &Y[si][idx_mass_fraction];
                                        }
                                        mu[count + 3] = d_equation_of_shear_viscosity_mixing_rules->
                                            getShearViscosity(
                                                &p[idx_pressure],
                                                &T[idx_temperature],
                                                Y_ptr);
                                    }
                                    
                                    const double dmudy = (double(1)/double(60)*(mu[6] - mu[0])
                                        - double(3)/double(20)*(mu[5] - mu[1])
                                        + double(3)/double(4)*(mu[4] - mu[2]))/dx[1];
                                    
                                    const int idx_vel       = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_BBB = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j - 3) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_BB  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j - 2) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_B   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j - 1) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_T   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + 1) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_TT  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + 2) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_TTT = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + 3) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const double dudy = (double(1)/double(60)*(u[idx_vel_y_TTT] - u[idx_vel_y_BBB])
                                        - double(3)/double(20)*(u[idx_vel_y_TT] - u[idx_vel_y_BB])
                                        + double(3)/double(4)*(u[idx_vel_y_T] - u[idx_vel_y_B]))/dx[1];
                                    
                                    const double d2udy2 = (double(1)/double(90)*(u[idx_vel_y_TTT] + u[idx_vel_y_BBB])
                                        - double(3)/double(20)*(u[idx_vel_y_TT] + u[idx_vel_y_BB])
                                        + double(3)/double(2)*(u[idx_vel_y_T] + u[idx_vel_y_B])
                                        - double(49)/double(18)*u[idx_vel])/(dx[1]*dx[1]);
                                    
                                    double visc_dvdx[7];
                                    
                                    for (int count = -3; count <= 3; count++)
                                    {
                                        const int idx_vel_x_LLL = (relative_idx_lo_0 + (i - 3) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + (j + count) +num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_x_LL  = (relative_idx_lo_0 + (i - 2) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + (j + count) +num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_x_L   = (relative_idx_lo_0 + (i - 1) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + (j + count) +num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_x_R   = (relative_idx_lo_0 + (i + 1) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + (j + count) +num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_x_RR  = (relative_idx_lo_0 + (i + 2) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + (j + count) +num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_x_RRR = (relative_idx_lo_0 + (i + 3) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + (j + count) +num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        visc_dvdx[count + 3] = mu[count + 3]*
                                            (double(1)/double(60)*(v[idx_vel_x_RRR] - v[idx_vel_x_LLL])
                                            - double(3)/double(20)*(v[idx_vel_x_RR] - v[idx_vel_x_LL])
                                            + double(3)/double(4)*(v[idx_vel_x_R] - v[idx_vel_x_L]))/dx[0];
                                    }
                                    
                                    const double dvisc_dvdxdy = (double(1)/double(60)*(visc_dvdx[6] - visc_dvdx[0])
                                        - double(3)/double(20)*(visc_dvdx[5] - visc_dvdx[1])
                                        + double(3)/double(4)*(visc_dvdx[4] - visc_dvdx[2]))/dx[1];
                                    
                                    value_to_add = (
                                        dmudy*dudy
                                        + mu[3]*d2udy2
                                        + dvisc_dvdxdy
                                        )*weight/((double) n_overlapped);
                                }
                                else if ((component_idx == 1) && (derivative_direction == 2))
                                {
                                    double mu[7];
                                    
                                    std::vector<const double*> Y_ptr;
                                    Y_ptr.resize(d_num_species);
                                    for (int count = -3; count <= 3; count++)
                                    {
                                        const int idx_mass_fraction = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                            (relative_idx_lo_2 + (k + count) + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                                ghostcell_dim_1_mass_fraction;
                                        
                                        const int idx_pressure = (relative_idx_lo_0 + i + num_ghosts_0_pressure) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure +
                                            (relative_idx_lo_2 + (k + count) + num_ghosts_2_pressure)*ghostcell_dim_0_pressure*
                                                ghostcell_dim_1_pressure;
                                        
                                        const int idx_temperature = (relative_idx_lo_0 + i + num_ghosts_0_temperature) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_temperature)*ghostcell_dim_0_temperature +
                                            (relative_idx_lo_2 + (k + count) + num_ghosts_2_temperature)*ghostcell_dim_0_temperature*
                                                ghostcell_dim_1_temperature;
                                        
                                        for (int si = 0; si < d_num_species; si++)
                                        {
                                            Y_ptr[si] = &Y[si][idx_mass_fraction];
                                        }
                                        mu[count + 3] = d_equation_of_shear_viscosity_mixing_rules->
                                            getShearViscosity(
                                                &p[idx_pressure],
                                                &T[idx_temperature],
                                                Y_ptr);
                                    }
                                    
                                    double visc_dvdx[7];
                                    
                                    for (int count = -3; count <= 3; count++)
                                    {
                                        const int idx_vel_x_LLL = (relative_idx_lo_0 + (i - 3) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + (k + count) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_x_LL  = (relative_idx_lo_0 + (i - 2) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + (k + count) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_x_L   = (relative_idx_lo_0 + (i - 1) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + (k + count) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_x_R   = (relative_idx_lo_0 + (i + 1) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + (k + count) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_x_RR  = (relative_idx_lo_0 + (i + 2) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + (k + count) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_x_RRR = (relative_idx_lo_0 + (i + 3) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + (k + count) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        visc_dvdx[count + 3] = mu[count + 3]*
                                            (double(1)/double(60)*(v[idx_vel_x_RRR] - v[idx_vel_x_LLL])
                                            - double(3)/double(20)*(v[idx_vel_x_RR] - v[idx_vel_x_LL])
                                            + double(3)/double(4)*(v[idx_vel_x_R] - v[idx_vel_x_L]))/dx[0];
                                    }
                                    
                                    const double dvisc_dvdxdz = (double(1)/double(60)*(visc_dvdx[6] - visc_dvdx[0])
                                        - double(3)/double(20)*(visc_dvdx[5] - visc_dvdx[1])
                                        + double(3)/double(4)*(visc_dvdx[4] - visc_dvdx[2]))/dx[2];
                                    
                                    double visc_dudy[7];
                                    
                                    for (int count = -3; count <= 3; count++)
                                    {
                                        const int idx_vel_y_BBB = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + (j - 3) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + (k + count) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_y_BB  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + (j - 2) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + (k + count) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_y_B   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + (j - 1) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + (k + count) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_y_T   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + (j + 1) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + (k + count) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_y_TT  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + (j + 2) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + (k + count) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_y_TTT = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + (j + 3) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + (k + count) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        visc_dudy[count + 3] = mu[count + 3]*
                                            (double(1)/double(60)*(u[idx_vel_y_TTT] - u[idx_vel_y_BBB])
                                            - double(3)/double(20)*(u[idx_vel_y_TT] - u[idx_vel_y_BB])
                                            + double(3)/double(4)*(u[idx_vel_y_T] - u[idx_vel_y_B]))/dx[1];
                                    }
                                    
                                    const double dvisc_dudydz = (double(1)/double(60)*(visc_dudy[6] - visc_dudy[0])
                                        - double(3)/double(20)*(visc_dudy[5] - visc_dudy[1])
                                        + double(3)/double(4)*(visc_dudy[4] - visc_dudy[2]))/dx[2];
                                    
                                    value_to_add = (
                                        dvisc_dvdxdz
                                        + dvisc_dudydz
                                        )*weight/((double) n_overlapped);
                                }
                                else if ((component_idx == 2) && (derivative_direction == 0))
                                {
                                    double mu[7];
                                    
                                    std::vector<const double*> Y_ptr;
                                    Y_ptr.resize(d_num_species);
                                    for (int count = -3; count <= 3; count++)
                                    {
                                        const int idx_mass_fraction = (relative_idx_lo_0 + (i + count) + num_ghosts_0_mass_fraction) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                            (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                                ghostcell_dim_1_mass_fraction;
                                        
                                        const int idx_pressure = (relative_idx_lo_0 + (i + count) + num_ghosts_0_pressure) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure +
                                            (relative_idx_lo_2 + k + num_ghosts_2_pressure)*ghostcell_dim_0_pressure*
                                                ghostcell_dim_1_pressure;
                                        
                                        const int idx_temperature = (relative_idx_lo_0 + (i + count) + num_ghosts_0_temperature) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_temperature)*ghostcell_dim_0_temperature +
                                            (relative_idx_lo_2 + k + num_ghosts_2_temperature)*ghostcell_dim_0_temperature*
                                                ghostcell_dim_1_temperature;
                                        
                                        for (int si = 0; si < d_num_species; si++)
                                        {
                                            Y_ptr[si] = &Y[si][idx_mass_fraction];
                                        }
                                        mu[count + 3] = d_equation_of_shear_viscosity_mixing_rules->
                                            getShearViscosity(
                                                &p[idx_pressure],
                                                &T[idx_temperature],
                                                Y_ptr);
                                    }
                                    
                                    const double dmudx = (double(1)/double(60)*(mu[6] - mu[0])
                                        - double(3)/double(20)*(mu[5] - mu[1])
                                        + double(3)/double(4)*(mu[4] - mu[2]))/dx[0];
                                    
                                    const int idx_vel       = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_LLL = (relative_idx_lo_0 + (i - 3) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_LL  = (relative_idx_lo_0 + (i - 2) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_L   = (relative_idx_lo_0 + (i - 1) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_R   = (relative_idx_lo_0 + (i + 1) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_RR  = (relative_idx_lo_0 + (i + 2) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_RRR = (relative_idx_lo_0 + (i + 3) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const double dwdx = (double(1)/double(60)*(w[idx_vel_x_RRR] - w[idx_vel_x_LLL])
                                        - double(3)/double(20)*(w[idx_vel_x_RR] - w[idx_vel_x_LL])
                                        + double(3)/double(4)*(w[idx_vel_x_R] - w[idx_vel_x_L]))/dx[0];
                                    
                                    const double d2wdx2 = (double(1)/double(90)*(w[idx_vel_x_RRR] + w[idx_vel_x_LLL])
                                        - double(3)/double(20)*(w[idx_vel_x_RR] + w[idx_vel_x_LL])
                                        + double(3)/double(2)*(w[idx_vel_x_R] + w[idx_vel_x_L])
                                        - double(49)/double(18)*w[idx_vel])/(dx[0]*dx[0]);
                                    
                                    double visc_dudz[7];
                                    
                                    for (int count = -3; count <= 3; count++)
                                    {
                                        const int idx_vel_z_BBB = (relative_idx_lo_0 + (i + count) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + (k - 3) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_z_BB  = (relative_idx_lo_0 + (i + count) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + (k - 2) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_z_B   = (relative_idx_lo_0 + (i + count) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + (k - 1) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_z_F   = (relative_idx_lo_0 + (i + count) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + (k + 1) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_z_FF  = (relative_idx_lo_0 + (i + count) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + (k + 2) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_z_FFF = (relative_idx_lo_0 + (i + count) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + (k + 3) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        visc_dudz[count + 3] = mu[count + 3]*
                                            (double(1)/double(60)*(u[idx_vel_z_FFF] - u[idx_vel_z_BBB])
                                            - double(3)/double(20)*(u[idx_vel_z_FF] - u[idx_vel_z_BB])
                                            + double(3)/double(4)*(u[idx_vel_z_F] - u[idx_vel_z_B]))/dx[2];
                                    }
                                    
                                    const double dvisc_dudzdx = (double(1)/double(60)*(visc_dudz[6] - visc_dudz[0])
                                        - double(3)/double(20)*(visc_dudz[5] - visc_dudz[1])
                                        + double(3)/double(4)*(visc_dudz[4] - visc_dudz[2]))/dx[0];
                                    
                                    value_to_add = (
                                        dmudx*dwdx
                                        + mu[3]*d2wdx2
                                        + dvisc_dudzdx
                                        )*weight/((double) n_overlapped);
                                }
                                else if ((component_idx == 2) && (derivative_direction == 1))
                                {
                                    double mu[7];
                                    
                                    std::vector<const double*> Y_ptr;
                                    Y_ptr.resize(d_num_species);
                                    for (int count = -3; count <= 3; count++)
                                    {
                                        const int idx_mass_fraction = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                            (relative_idx_lo_1 + (j + count) + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                            (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                                ghostcell_dim_1_mass_fraction;
                                        
                                        const int idx_pressure = (relative_idx_lo_0 + i + num_ghosts_0_pressure) +
                                            (relative_idx_lo_1 + (j + count) + num_ghosts_1_pressure)*ghostcell_dim_0_pressure +
                                            (relative_idx_lo_2 + k + num_ghosts_2_pressure)*ghostcell_dim_0_pressure*
                                                ghostcell_dim_1_pressure;
                                        
                                        const int idx_temperature = (relative_idx_lo_0 + i + num_ghosts_0_temperature) +
                                            (relative_idx_lo_1 + (j + count) + num_ghosts_1_temperature)*ghostcell_dim_0_temperature +
                                            (relative_idx_lo_2 + k + num_ghosts_2_temperature)*ghostcell_dim_0_temperature*
                                                ghostcell_dim_1_temperature;
                                        
                                        for (int si = 0; si < d_num_species; si++)
                                        {
                                            Y_ptr[si] = &Y[si][idx_mass_fraction];
                                        }
                                        mu[count + 3] = d_equation_of_shear_viscosity_mixing_rules->
                                            getShearViscosity(
                                                &p[idx_pressure],
                                                &T[idx_temperature],
                                                Y_ptr);
                                    }
                                    
                                    double visc_dwdx[7];
                                    
                                    for (int count = -3; count <= 3; count++)
                                    {
                                        const int idx_vel_x_LLL = (relative_idx_lo_0 + (i - 3) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + (j + count) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_x_LL  = (relative_idx_lo_0 + (i - 2) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + (j + count) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_x_L   = (relative_idx_lo_0 + (i - 1) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + (j + count) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_x_R   = (relative_idx_lo_0 + (i + 1) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + (j + count) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_x_RR  = (relative_idx_lo_0 + (i + 2) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + (j + count) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_x_RRR = (relative_idx_lo_0 + (i + 3) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + (j + count) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        visc_dwdx[count + 3] = mu[count + 3]*
                                            (double(1)/double(60)*(w[idx_vel_x_RRR] - w[idx_vel_x_LLL])
                                            - double(3)/double(20)*(w[idx_vel_x_RR] - w[idx_vel_x_LL])
                                            + double(3)/double(4)*(w[idx_vel_x_R] - w[idx_vel_x_L]))/dx[0];
                                    }
                                    
                                    const double dvisc_dwdxdy = (double(1)/double(60)*(visc_dwdx[6] - visc_dwdx[0])
                                        - double(3)/double(20)*(visc_dwdx[5] - visc_dwdx[1])
                                        + double(3)/double(4)*(visc_dwdx[4] - visc_dwdx[2]))/dx[1];
                                    
                                    double visc_dudz[7];
                                    
                                    for (int count = -3; count <= 3; count++)
                                    {
                                        const int idx_vel_z_BBB = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + (j + count) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + (k - 3) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_z_BB  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + (j + count) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + (k - 2) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_z_B   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + (j + count) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + (k - 1) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_z_F   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + (j + count) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + (k + 1) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_z_FF  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + (j + count) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + (k + 2) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_z_FFF = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + (j + count) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + (k + 3) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        visc_dudz[count + 3] = mu[count + 3]*
                                            (double(1)/double(60)*(u[idx_vel_z_FFF] - u[idx_vel_z_BBB])
                                            - double(3)/double(20)*(u[idx_vel_z_FF] - u[idx_vel_z_BB])
                                            + double(3)/double(4)*(u[idx_vel_z_F] - u[idx_vel_z_B]))/dx[2];
                                    }
                                    
                                    const double dvisc_dudzdy = (double(1)/double(60)*(visc_dudz[6] - visc_dudz[0])
                                        - double(3)/double(20)*(visc_dudz[5] - visc_dudz[1])
                                        + double(3)/double(4)*(visc_dudz[4] - visc_dudz[2]))/dx[1];
                                    
                                    value_to_add = (
                                        dvisc_dwdxdy
                                        + dvisc_dudzdy
                                        )*weight/((double) n_overlapped);
                                }
                                else if ((component_idx == 2) && (derivative_direction == 2))
                                {
                                    double mu[7];
                                    
                                    std::vector<const double*> Y_ptr;
                                    Y_ptr.resize(d_num_species);
                                    for (int count = -3; count <= 3; count++)
                                    {
                                        const int idx_mass_fraction = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                            (relative_idx_lo_2 + (k + count) + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                                ghostcell_dim_1_mass_fraction;
                                        
                                        const int idx_pressure = (relative_idx_lo_0 + i + num_ghosts_0_pressure) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure +
                                            (relative_idx_lo_2 + (k + count) + num_ghosts_2_pressure)*ghostcell_dim_0_pressure*
                                                ghostcell_dim_1_pressure;
                                        
                                        const int idx_temperature = (relative_idx_lo_0 + i + num_ghosts_0_temperature) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_temperature)*ghostcell_dim_0_temperature +
                                            (relative_idx_lo_2 + (k + count) + num_ghosts_2_temperature)*ghostcell_dim_0_temperature*
                                                ghostcell_dim_1_temperature;
                                        
                                        for (int si = 0; si < d_num_species; si++)
                                        {
                                            Y_ptr[si] = &Y[si][idx_mass_fraction];
                                        }
                                        mu[count + 3] = d_equation_of_shear_viscosity_mixing_rules->
                                            getShearViscosity(
                                                &p[idx_pressure],
                                                &T[idx_temperature],
                                                Y_ptr);
                                    }
                                    
                                    const double dmudz = (double(1)/double(60)*(mu[6] - mu[0])
                                        - double(3)/double(20)*(mu[5] - mu[1])
                                        + double(3)/double(4)*(mu[4] - mu[2]))/dx[2];
                                    
                                    const int idx_vel       = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_BBB = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k - 3) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_BB  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k - 2) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_B   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k - 1) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_F   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k + 1) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_FF  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k + 2) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_FFF = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k + 3) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const double dudz = (double(1)/double(60)*(u[idx_vel_z_FFF] - u[idx_vel_z_BBB])
                                        - double(3)/double(20)*(u[idx_vel_z_FF] - u[idx_vel_z_BB])
                                        + double(3)/double(4)*(u[idx_vel_z_F] - u[idx_vel_z_B]))/dx[2];
                                    
                                    const double d2udz2 = (double(1)/double(90)*(u[idx_vel_z_FFF] + u[idx_vel_z_BBB])
                                        - double(3)/double(20)*(u[idx_vel_z_FF] + u[idx_vel_z_BB])
                                        + double(3)/double(2)*(u[idx_vel_z_F] + u[idx_vel_z_B])
                                        - double(49)/double(18)*u[idx_vel])/(dx[2]*dx[2]);
                                    
                                    double visc_dwdx[7];
                                    
                                    for (int count = -3; count <= 3; count++)
                                    {
                                        const int idx_vel_x_LLL = (relative_idx_lo_0 + (i - 3) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + (k + count) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_x_LL  = (relative_idx_lo_0 + (i - 2) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + (k + count) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_x_L   = (relative_idx_lo_0 + (i - 1) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + (k + count) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_x_R   = (relative_idx_lo_0 + (i + 1) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + (k + count) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_x_RR  = (relative_idx_lo_0 + (i + 2) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + (k + count) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        const int idx_vel_x_RRR = (relative_idx_lo_0 + (i + 3) + num_ghosts_0_velocity) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                            (relative_idx_lo_2 + (k + count) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                                ghostcell_dim_1_velocity;
                                        
                                        visc_dwdx[count + 3] = mu[count + 3]*
                                            (double(1)/double(60)*(w[idx_vel_x_RRR] - w[idx_vel_x_LLL])
                                            - double(3)/double(20)*(w[idx_vel_x_RR] - w[idx_vel_x_LL])
                                            + double(3)/double(4)*(w[idx_vel_x_R] - w[idx_vel_x_L]))/dx[0];
                                    }
                                    
                                    const double dvisc_dwdxdz = (double(1)/double(60)*(visc_dwdx[6] - visc_dwdx[0])
                                        - double(3)/double(20)*(visc_dwdx[5] - visc_dwdx[1])
                                        + double(3)/double(4)*(visc_dwdx[4] - visc_dwdx[2]))/dx[2];
                                    
                                    value_to_add = (
                                        dmudz*dudz
                                        + mu[3]*d2udz2
                                        + dvisc_dwdxdz
                                        )*weight/((double) n_overlapped);
                                }
                                else
                                {
                                    TBOX_ERROR(d_object_name
                                        << ": "
                                        << "Cannot compute derivative of shear stress component for three-dimensional problem!\n"
                                        << "component_idx = " << component_idx << " and derivative_direction = "
                                        << derivative_direction << " given!\n"
                                        << std::endl);
                                }
                                
                                /*
                                 * Add the data.
                                 */
                                
                                for (int ii = 0; ii < ratio_to_finest_level_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratio_to_finest_level_0 + ii;
                                    
                                    tau_ij_der_avg_local[idx_fine] += value_to_add;
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average of derivative of shear stress component.
         */
        
        mpi.Allreduce(
            tau_ij_der_avg_local,
            tau_ij_der_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(tau_ij_der_avg_local);
    }
    
    return averaged_derivative_tau_ij;
}


/*
 * Compute averaged value (on product of variable derivatives and shear stress component) with only x direction
 * as inhomogeneous direction.
 * component_idx:
 * 0: tau11
 * 1: tau12
 * 2: tau13
 * 3: tau22
 * 4: tau23
 * 5: tau33
 */
std::vector<double>
RTIRMIStatisticsUtilities::getAveragedQuantityWithShearStressComponentWithInhomogeneousXDirection(
    const std::vector<std::string>& quantity_names,
    const std::vector<int>& component_indices,
    const std::vector<bool>& use_derivative,
    const std::vector<int>& derivative_directions,
    const std::vector<bool>& use_reciprocal,
    const int shear_stress_component_idx,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
    TBOX_ASSERT(d_num_ghosts_derivative == 3);
    
    int num_quantities = static_cast<int>(quantity_names.size());
    
    TBOX_ASSERT(static_cast<int>(component_indices.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(use_derivative.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(derivative_directions.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(use_reciprocal.size()) == num_quantities);
    
    int num_use_derivative = 0;
    
    if (d_dim == tbox::Dimension(1))
    {
        for (int qi = 0; qi < num_quantities; qi++)
        {
            if (use_derivative[qi])
            {
                num_use_derivative++;
                if (derivative_directions[qi] < 0 || derivative_directions[qi] > 0)
                {
                    TBOX_ERROR(d_object_name
                        << ": RTIRMIStatisticsUtilities::"
                        << "getAveragedQuantityWithShearStressComponentWithInhomogeneousXDirection():\n"
                        << "Cannot take derivative for one-dimensional problem!\n"
                        << "derivative_directions[" << qi << "] = " << derivative_directions[qi] << " given!\n"
                        << std::endl);
                }
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        for (int qi = 0; qi < num_quantities; qi++)
        {
            if (use_derivative[qi])
            {
                num_use_derivative++;
                if (derivative_directions[qi] < 0 || derivative_directions[qi] > 1)
                {
                    TBOX_ERROR(d_object_name
                        << ": RTIRMIStatisticsUtilities::"
                        << "getAveragedQuantityWithShearStressComponentWithInhomogeneousXDirection():\n"
                        << "Cannot take derivative for two-dimensional problem!\n"
                        << "derivative_directions[" << qi << "] = " << derivative_directions[qi] << " given!\n"
                        << std::endl);
                }
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        for (int qi = 0; qi < num_quantities; qi++)
        {
            if (use_derivative[qi])
            {
                num_use_derivative++;
                if (derivative_directions[qi] < 0 || derivative_directions[qi] > 2)
                {
                    TBOX_ERROR(d_object_name
                        << ": RTIRMIBudgetsUtilities::"
                        << "getAveragedQuantityWithShearStressComponentWithInhomogeneousXDirection():\n"
                        << "Cannot take derivative for three-dimensional problem!\n"
                        << "derivative_directions[" << qi << "] = " << derivative_directions[qi] << " given!\n"
                        << std::endl);
                }
            }
        }
    }
    
    std::vector<double> averaged_quantity;
    
    HAMERS_SHARED_PTR<FlowModel> d_flow_model_tmp = d_flow_model.lock();
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Get the refinement ratio from the finest level to the coarsest level.
     */
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    hier::IntVector ratio_finest_level_to_coarsest_level =
        patch_hierarchy->getRatioToCoarserLevel(num_levels - 1);
    for (int li = num_levels - 2; li > 0 ; li--)
    {
        ratio_finest_level_to_coarsest_level *= patch_hierarchy->getRatioToCoarserLevel(li);
    }
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    HAMERS_SHARED_PTR<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *patch_hierarchy,
            0,
            num_levels - 1));
    
    /*
     * Get the number of cells of physical domain refined to the finest level.
     */
    
    const hier::BoxContainer& physical_domain = d_grid_geometry->getPhysicalDomain();
    const hier::Box& physical_domain_box = physical_domain.front();
    const hier::IntVector& physical_domain_dims = physical_domain_box.numberCells();
    const hier::IntVector finest_level_dims = physical_domain_dims*ratio_finest_level_to_coarsest_level;
    
    /*
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_x(
            new DerivativeFirstOrder(
                "first order derivative in x-direction",
                d_dim,
                DIRECTION::X_DIRECTION,
                d_num_ghosts_derivative));
        
        hier::IntVector num_ghosts_der = hier::IntVector::getOne(d_dim)*d_num_ghosts_derivative;
        
        const int finest_level_dim_0 = finest_level_dims[0];
        
        double* avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_quantity.resize(finest_level_dim_0);
        double* avg_global = averaged_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            avg_local[i]  = double(0);
            avg_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratio_to_coarsest_level =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarsest_level *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = ratio_finest_level_to_coarsest_level/ratio_to_coarsest_level;
            
            const int ratio_to_finest_level_0 = ratio_to_finest_level[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], num_ghosts_der));
                }
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", num_ghosts_der));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTIONS", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("PRESSURE", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("TEMPERATURE", hier::IntVector::getZero(d_dim)));
                
                d_flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                d_flow_model_tmp->allocateMemoryForDerivedCellData();
                
                d_flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_velocity =
                    d_flow_model_tmp->getCellData("VELOCITY");
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_mass_fraction =
                    d_flow_model_tmp->getCellData("MASS_FRACTIONS");
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_pressure =
                    d_flow_model_tmp->getCellData("PRESSURE");
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_temperature =
                    d_flow_model_tmp->getCellData("TEMPERATURE");
                
                double* u = data_velocity->getPointer(0);
                std::vector<double*> Y;
                Y.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    Y.push_back(data_mass_fraction->getPointer(si));
                }
                double* p = data_pressure->getPointer(0);
                double* T = data_temperature->getPointer(0);
                
                std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = d_flow_model_tmp->getCellData(quantity_names[qi]);
                }
                
                std::vector<double*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_velocity = data_velocity->getGhostCellWidth();
                const hier::IntVector num_ghosts_mass_fraction = data_mass_fraction->getGhostCellWidth();
                const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
                const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
                
                const int num_ghosts_0_velocity = num_ghosts_velocity[0];
                const int num_ghosts_0_mass_fraction = num_ghosts_mass_fraction[0];
                const int num_ghosts_0_pressure = num_ghosts_pressure[0];
                const int num_ghosts_0_temperature = num_ghosts_temperature[0];
                
                std::vector<int> num_ghosts_0_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                }
                
                /*
                 * Initialize cell data for the derivatives and get pointers to the cell data.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_derivative;
                std::vector<double*> der_qi;
                
                if (num_use_derivative > 0)
                {
                    data_derivative = HAMERS_MAKE_SHARED<pdat::CellData<double> >(
                        patch_box, num_use_derivative, hier::IntVector::getZero(d_dim));
                    
                    der_qi.resize(num_use_derivative);
                    for (int qi = 0; qi < num_use_derivative; qi++)
                    {
                        der_qi[qi] = data_derivative->getPointer(qi);
                    }
                }
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    int count_derivative = 0;
                    for (int qi = 0; qi < num_quantities; qi++)
                    {
                        if (use_derivative[qi] && derivative_directions[qi] == 0)
                        {
                            derivative_first_order_x->computeDerivative(
                                data_derivative,
                                data_quantities[qi],
                                dx[0],
                                patch_visible_box,
                                count_derivative,
                                component_indices[qi]);
                            
                            count_derivative++;
                        }
                    }
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        /*
                         * Compute the index of the data point and count how many times the data is repeated.
                         */
                        
                        const hier::Index idx_pt(tbox::Dimension(1), idx_lo_0 + i);
                        
                        int n_overlapped = 1;
                        
                        for (hier::BoxContainer::BoxContainerConstIterator iob(
                                patch_overlapped_visible_boxes.begin());
                             iob != patch_overlapped_visible_boxes.end();
                             iob++)
                        {
                            const hier::Box& patch_overlapped_visible_box = *iob;
                            
                            if (patch_overlapped_visible_box.contains(idx_pt))
                            {
                                n_overlapped++;
                            }
                        }
                        
                        /*
                         * Compute the shear stress component.
                         */
                        
                        double tau_ij = double(0);
                        
                        // Compute linear indices of mass fraction, pressure and temperature.
                        const int idx_mass_fraction = relative_idx_lo_0 + i + num_ghosts_0_mass_fraction;
                        const int idx_pressure = relative_idx_lo_0 + i + num_ghosts_0_pressure;
                        const int idx_temperature = relative_idx_lo_0 + i + num_ghosts_0_temperature;
                        
                        if (shear_stress_component_idx == 0)
                        {
                            std::vector<const double*> Y_ptr;
                            Y_ptr.resize(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_ptr[si] = &Y[si][idx_mass_fraction];
                            }
                            const double mu = d_equation_of_shear_viscosity_mixing_rules->
                                getShearViscosity(
                                    &p[idx_pressure],
                                    &T[idx_temperature],
                                    Y_ptr);
                            const double mu_v = d_equation_of_bulk_viscosity_mixing_rules->
                                getBulkViscosity(
                                    &p[idx_pressure],
                                    &T[idx_temperature],
                                    Y_ptr);
                            
                            const int idx_vel_x_LLL = relative_idx_lo_0 + (i - 3) + num_ghosts_0_velocity;
                            const int idx_vel_x_LL  = relative_idx_lo_0 + (i - 2) + num_ghosts_0_velocity;
                            const int idx_vel_x_L   = relative_idx_lo_0 + (i - 1) + num_ghosts_0_velocity;
                            const int idx_vel_x_R   = relative_idx_lo_0 + (i + 1) + num_ghosts_0_velocity;
                            const int idx_vel_x_RR  = relative_idx_lo_0 + (i + 2) + num_ghosts_0_velocity;
                            const int idx_vel_x_RRR = relative_idx_lo_0 + (i + 3) + num_ghosts_0_velocity;
                            
                            const double dudx = (double(1)/double(60)*(u[idx_vel_x_RRR] - u[idx_vel_x_LLL])
                                - double(3)/double(20)*(u[idx_vel_x_RR] - u[idx_vel_x_LL])
                                + double(3)/double(4)*(u[idx_vel_x_R] - u[idx_vel_x_L]))/dx[0];
                            
                            tau_ij = (double(4)/double(3)*mu + mu_v)*dudx;
                        }
                        else
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "Cannot compute shear stress component for one-dimensional problem!\n"
                                << "shear_stress_component_idx = " << shear_stress_component_idx << " given!\n"
                                << std::endl);
                        }
                        
                        /*
                         * Compute the linear indices and the data to add.
                         */
                        
                        double avg = double(1);
                        
                        count_derivative = 0;
                        for (int qi = 0; qi < num_quantities; qi++)
                        {
                            if (use_reciprocal[qi])
                            {
                                if (use_derivative[qi])
                                {
                                    const int idx_der = relative_idx_lo_0 + i;
                                    
                                    avg /= der_qi[count_derivative][idx_der];
                                    count_derivative++;
                                }
                                else
                                {
                                    const int idx_qi = relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi];
                                    
                                    avg /= u_qi[qi][idx_qi];
                                }
                            }
                            else
                            {
                                if (use_derivative[qi])
                                {
                                    const int idx_der = relative_idx_lo_0 + i;
                                    
                                    avg *= der_qi[count_derivative][idx_der];
                                    count_derivative++;
                                }
                                else
                                {
                                    const int idx_qi = relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi];
                                    
                                    avg *= u_qi[qi][idx_qi];
                                }
                            }
                        }
                        
                        /*
                         * Add the data.
                         */
                        
                        for (int ii = 0; ii < ratio_to_finest_level_0; ii++)
                        {
                            const int idx_fine = (idx_lo_0 + i)*ratio_to_finest_level_0 + ii;
                            
                            avg_local[idx_fine] += (avg*tau_ij/((double) n_overlapped));
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        mpi.Allreduce(
            avg_local,
            avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(avg_local);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_x(
            new DerivativeFirstOrder(
                "first order derivative in x-direction",
                d_dim,
                DIRECTION::X_DIRECTION,
                d_num_ghosts_derivative));
        
        HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_y(
            new DerivativeFirstOrder(
                "first order derivative in y-direction",
                d_dim,
                DIRECTION::Y_DIRECTION,
                d_num_ghosts_derivative));
        
        hier::IntVector num_ghosts_der = hier::IntVector::getOne(d_dim)*d_num_ghosts_derivative;
        
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        
        double* avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_quantity.resize(finest_level_dim_0);
        double* avg_global = averaged_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            avg_local[i]  = double(0);
            avg_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratio_to_coarsest_level =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarsest_level *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = ratio_finest_level_to_coarsest_level/ratio_to_coarsest_level;
            
            const int ratio_to_finest_level_0 = ratio_to_finest_level[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], num_ghosts_der));
                }
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", num_ghosts_der));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTIONS", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("PRESSURE", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("TEMPERATURE", hier::IntVector::getZero(d_dim)));
                
                d_flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                d_flow_model_tmp->allocateMemoryForDerivedCellData();
                
                d_flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_velocity =
                    d_flow_model_tmp->getCellData("VELOCITY");
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_mass_fraction =
                    d_flow_model_tmp->getCellData("MASS_FRACTIONS");
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_pressure =
                    d_flow_model_tmp->getCellData("PRESSURE");
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_temperature =
                    d_flow_model_tmp->getCellData("TEMPERATURE");
                
                double* u = data_velocity->getPointer(0);
                double* v = data_velocity->getPointer(1);
                std::vector<double*> Y;
                Y.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    Y.push_back(data_mass_fraction->getPointer(si));
                }
                double* p = data_pressure->getPointer(0);
                double* T = data_temperature->getPointer(0);
                
                std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = d_flow_model_tmp->getCellData(quantity_names[qi]);
                }
                
                std::vector<double*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_velocity = data_velocity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_velocity = data_velocity->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_mass_fraction = data_mass_fraction->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mass_fraction = data_mass_fraction->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_pressure = data_pressure->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_temperature = data_temperature->getGhostBox().numberCells();
                
                const int num_ghosts_0_velocity = num_ghosts_velocity[0];
                const int num_ghosts_1_velocity = num_ghosts_velocity[1];
                const int ghostcell_dim_0_velocity = ghostcell_dims_velocity[0];
                
                const int num_ghosts_0_mass_fraction = num_ghosts_mass_fraction[0];
                const int num_ghosts_1_mass_fraction = num_ghosts_mass_fraction[1];
                const int ghostcell_dim_0_mass_fraction = ghostcell_dims_mass_fraction[0];
                
                const int num_ghosts_0_pressure = num_ghosts_pressure[0];
                const int num_ghosts_1_pressure = num_ghosts_pressure[1];
                const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
                
                const int num_ghosts_0_temperature = num_ghosts_temperature[0];
                const int num_ghosts_1_temperature = num_ghosts_temperature[1];
                const int ghostcell_dim_0_temperature = ghostcell_dims_temperature[0];
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                }
                
                /*
                 * Initialize cell data for the derivatives and get pointers to the cell data.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_derivative;
                std::vector<double*> der_qi;
                
                if (num_use_derivative > 0)
                {
                    data_derivative = HAMERS_MAKE_SHARED<pdat::CellData<double> >(
                        patch_box, num_use_derivative, hier::IntVector::getZero(d_dim));
                    
                    der_qi.resize(num_use_derivative);
                    for (int qi = 0; qi < num_use_derivative; qi++)
                    {
                        der_qi[qi] = data_derivative->getPointer(qi);
                    }
                }
                
                const hier::IntVector patch_interior_dims = patch_box.numberCells();
                const int patch_interior_dim_0 = patch_interior_dims[0];
                
                const double weight = dx[1]/L_y;
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    int count_derivative = 0;
                    for (int qi = 0; qi < num_quantities; qi++)
                    {
                        if (use_derivative[qi] && derivative_directions[qi] == 0)
                        {
                            derivative_first_order_x->computeDerivative(
                                data_derivative,
                                data_quantities[qi],
                                dx[0],
                                patch_visible_box,
                                count_derivative,
                                component_indices[qi]);
                            
                            count_derivative++;
                        }
                        else if (use_derivative[qi] && derivative_directions[qi] == 1)
                        {
                            derivative_first_order_y->computeDerivative(
                                data_derivative,
                                data_quantities[qi],
                                dx[1],
                                patch_visible_box,
                                count_derivative,
                                component_indices[qi]);
                            
                            count_derivative++;
                        }
                    }
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            /*
                             * Compute the index of the data point and count how many times the data is repeated.
                             */
                            
                            const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j);
                            
                            int n_overlapped = 1;
                            
                            for (hier::BoxContainer::BoxContainerConstIterator iob(
                                    patch_overlapped_visible_boxes.begin());
                                 iob != patch_overlapped_visible_boxes.end();
                                 iob++)
                            {
                                const hier::Box& patch_overlapped_visible_box = *iob;
                                
                                if (patch_overlapped_visible_box.contains(idx_pt))
                                {
                                    n_overlapped++;
                                }
                            }
                            
                            /*
                             * Compute the shear stress component.
                             */
                            
                            double tau_ij = double(0);
                            
                            const int idx_mass_fraction = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                            
                            const int idx_pressure = (relative_idx_lo_0 + i + num_ghosts_0_pressure) +
                                (relative_idx_lo_1 + j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure;
                            
                            const int idx_temperature = (relative_idx_lo_0 + i + num_ghosts_0_temperature) +
                                (relative_idx_lo_1 + j + num_ghosts_1_temperature)*ghostcell_dim_0_temperature;
                            
                            const int idx_vel_x_LLL = (relative_idx_lo_0 + (i - 3) + num_ghosts_0_velocity) +
                                (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                            
                            const int idx_vel_x_LL  = (relative_idx_lo_0 + (i - 2) + num_ghosts_0_velocity) +
                                (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                            
                            const int idx_vel_x_L   = (relative_idx_lo_0 + (i - 1) + num_ghosts_0_velocity) +
                                (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                            
                            const int idx_vel_x_R   = (relative_idx_lo_0 + (i + 1) + num_ghosts_0_velocity) +
                                (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                            
                            const int idx_vel_x_RR  = (relative_idx_lo_0 + (i + 2) + num_ghosts_0_velocity) +
                                (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                            
                            const int idx_vel_x_RRR = (relative_idx_lo_0 + (i + 3) + num_ghosts_0_velocity) +
                                (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                            
                            const int idx_vel_y_BBB = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                (relative_idx_lo_1 + (j - 3) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                            
                            const int idx_vel_y_BB  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                (relative_idx_lo_1 + (j - 2) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                            
                            const int idx_vel_y_B   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                (relative_idx_lo_1 + (j - 1) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                            
                            const int idx_vel_y_T   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                (relative_idx_lo_1 + (j + 1) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                            
                            const int idx_vel_y_TT  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                (relative_idx_lo_1 + (j + 2) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                            
                            const int idx_vel_y_TTT = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                (relative_idx_lo_1 + (j + 3) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                            
                            if (shear_stress_component_idx == 0)
                            {
                                std::vector<const double*> Y_ptr;
                                Y_ptr.resize(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_ptr[si] = &Y[si][idx_mass_fraction];
                                }
                                const double mu = d_equation_of_shear_viscosity_mixing_rules->
                                    getShearViscosity(
                                        &p[idx_pressure],
                                        &T[idx_temperature],
                                        Y_ptr);
                                const double mu_v = d_equation_of_bulk_viscosity_mixing_rules->
                                    getBulkViscosity(
                                        &p[idx_pressure],
                                        &T[idx_temperature],
                                        Y_ptr);
                                
                                const double dudx = (double(1)/double(60)*(u[idx_vel_x_RRR] - u[idx_vel_x_LLL])
                                    - double(3)/double(20)*(u[idx_vel_x_RR] - u[idx_vel_x_LL])
                                    + double(3)/double(4)*(u[idx_vel_x_R] - u[idx_vel_x_L]))/dx[0];
                                
                                const double dvdy = (double(1)/double(60)*(v[idx_vel_y_TTT] - v[idx_vel_y_BBB])
                                    - double(3)/double(20)*(v[idx_vel_y_TT] - v[idx_vel_y_BB])
                                    + double(3)/double(4)*(v[idx_vel_y_T] - v[idx_vel_y_B]))/dx[1];
                                
                                tau_ij = (double(4)/double(3)*mu + mu_v)*dudx - (double(2)/double(3)*mu - mu_v)*dvdy;
                            }
                            else if (shear_stress_component_idx == 1)
                            {
                                std::vector<const double*> Y_ptr;
                                Y_ptr.resize(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_ptr[si] = &Y[si][idx_mass_fraction];
                                }
                                const double mu = d_equation_of_shear_viscosity_mixing_rules->
                                    getShearViscosity(
                                        &p[idx_pressure],
                                        &T[idx_temperature],
                                        Y_ptr);
                                
                                const double dvdx = (double(1)/double(60)*(v[idx_vel_x_RRR] - v[idx_vel_x_LLL])
                                    - double(3)/double(20)*(v[idx_vel_x_RR] - v[idx_vel_x_LL])
                                    + double(3)/double(4)*(v[idx_vel_x_R] - v[idx_vel_x_L]))/dx[0];
                                
                                const double dudy = (double(1)/double(60)*(u[idx_vel_y_TTT] - u[idx_vel_y_BBB])
                                    - double(3)/double(20)*(u[idx_vel_y_TT] - u[idx_vel_y_BB])
                                    + double(3)/double(4)*(u[idx_vel_y_T] - u[idx_vel_y_B]))/dx[1];
                                
                                tau_ij = mu*(dudy + dvdx);
                            }
                            else if (shear_stress_component_idx == 3)
                            {
                                std::vector<const double*> Y_ptr;
                                Y_ptr.resize(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_ptr[si] = &Y[si][idx_mass_fraction];
                                }
                                const double mu = d_equation_of_shear_viscosity_mixing_rules->
                                    getShearViscosity(
                                        &p[idx_pressure],
                                        &T[idx_temperature],
                                        Y_ptr);
                                const double mu_v = d_equation_of_bulk_viscosity_mixing_rules->
                                    getBulkViscosity(
                                        &p[idx_pressure],
                                        &T[idx_temperature],
                                        Y_ptr);
                                
                                const double dudx = (double(1)/double(60)*(u[idx_vel_x_RRR] - u[idx_vel_x_LLL])
                                    - double(3)/double(20)*(u[idx_vel_x_RR] - u[idx_vel_x_LL])
                                    + double(3)/double(4)*(u[idx_vel_x_R] - u[idx_vel_x_L]))/dx[0];
                                
                                const double dvdy = (double(1)/double(60)*(v[idx_vel_y_TTT] - v[idx_vel_y_BBB])
                                    - double(3)/double(20)*(v[idx_vel_y_TT] - v[idx_vel_y_BB])
                                    + double(3)/double(4)*(v[idx_vel_y_T] - v[idx_vel_y_B]))/dx[1];
                                
                                tau_ij = (double(4)/double(3)*mu + mu_v)*dvdy - (double(2)/double(3)*mu - mu_v)*dudx;
                            }
                            else
                            {
                                TBOX_ERROR(d_object_name
                                    << ": "
                                    << "Cannot compute shear stress component for two-dimensional problem!\n"
                                    << "component_idx = " << shear_stress_component_idx << " given!\n"
                                    << std::endl);
                            }
                            
                            /*
                             * Compute the linear indices and the data to add.
                             */
                            
                            double avg = double(1);
                            
                            count_derivative = 0;
                            for (int qi = 0; qi < num_quantities; qi++)
                            {
                                if (use_reciprocal[qi])
                                {
                                    if (use_derivative[qi])
                                    {
                                        const int idx_der = (relative_idx_lo_0 + i) +
                                            (relative_idx_lo_1 + j)*patch_interior_dim_0;
                                        
                                        avg /= der_qi[count_derivative][idx_der];
                                        
                                        count_derivative++;
                                    }
                                    else
                                    {
                                        const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                        
                                        avg /= u_qi[qi][idx_qi];
                                    }
                                }
                                else
                                {
                                    if (use_derivative[qi])
                                    {
                                        const int idx_der = (relative_idx_lo_0 + i) +
                                            (relative_idx_lo_1 + j)*patch_interior_dim_0;
                                        
                                        avg *= der_qi[count_derivative][idx_der];
                                        
                                        count_derivative++;
                                    }
                                    else
                                    {
                                        const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                        
                                        avg *= u_qi[qi][idx_qi];
                                    }
                                }
                            }
                            
                            /*
                             * Add the data.
                             */
                            
                            for (int ii = 0; ii < ratio_to_finest_level_0; ii++)
                            {
                                const int idx_fine = (idx_lo_0 + i)*ratio_to_finest_level_0 + ii;
                                
                                avg_local[idx_fine] += (avg*tau_ij*weight/((double) n_overlapped));
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        mpi.Allreduce(
            avg_local,
            avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(avg_local);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_x(
            new DerivativeFirstOrder(
                "first order derivative in x-direction",
                d_dim,
                DIRECTION::X_DIRECTION,
                d_num_ghosts_derivative));
        
        HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_y(
            new DerivativeFirstOrder(
                "first order derivative in y-direction",
                d_dim,
                DIRECTION::Y_DIRECTION,
                d_num_ghosts_derivative));
        
        HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_z(
            new DerivativeFirstOrder(
                "first order derivative in z-direction",
                d_dim,
                DIRECTION::Z_DIRECTION,
                d_num_ghosts_derivative));
        
        hier::IntVector num_ghosts_der = hier::IntVector::getOne(d_dim)*d_num_ghosts_derivative;
        
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        const double L_z = x_hi[2] - x_lo[2];
        
        double* avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_quantity.resize(finest_level_dim_0);
        double* avg_global = averaged_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            avg_local[i]  = double(0);
            avg_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratio_to_coarsest_level =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarsest_level *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = ratio_finest_level_to_coarsest_level/ratio_to_coarsest_level;
            
            const int ratio_to_finest_level_0 = ratio_to_finest_level[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], num_ghosts_der));
                }
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", num_ghosts_der));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTIONS", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("PRESSURE", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("TEMPERATURE", hier::IntVector::getZero(d_dim)));
                
                d_flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                d_flow_model_tmp->allocateMemoryForDerivedCellData();
                
                d_flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_velocity =
                    d_flow_model_tmp->getCellData("VELOCITY");
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_mass_fraction =
                    d_flow_model_tmp->getCellData("MASS_FRACTIONS");
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_pressure =
                    d_flow_model_tmp->getCellData("PRESSURE");
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_temperature =
                    d_flow_model_tmp->getCellData("TEMPERATURE");
                
                double* u = data_velocity->getPointer(0);
                double* v = data_velocity->getPointer(1);
                double* w = data_velocity->getPointer(2);
                std::vector<double*> Y;
                Y.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    Y.push_back(data_mass_fraction->getPointer(si));
                }
                double* p = data_pressure->getPointer(0);
                double* T = data_temperature->getPointer(0);
                
                std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = d_flow_model_tmp->getCellData(quantity_names[qi]);
                }
                
                std::vector<double*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_velocity = data_velocity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_velocity = data_velocity->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_mass_fraction = data_mass_fraction->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mass_fraction = data_mass_fraction->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_pressure = data_pressure->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_temperature = data_temperature->getGhostBox().numberCells();
                
                const int num_ghosts_0_velocity = num_ghosts_velocity[0];
                const int num_ghosts_1_velocity = num_ghosts_velocity[1];
                const int num_ghosts_2_velocity = num_ghosts_velocity[2];
                const int ghostcell_dim_0_velocity = ghostcell_dims_velocity[0];
                const int ghostcell_dim_1_velocity = ghostcell_dims_velocity[1];
                
                const int num_ghosts_0_mass_fraction = num_ghosts_mass_fraction[0];
                const int num_ghosts_1_mass_fraction = num_ghosts_mass_fraction[1];
                const int num_ghosts_2_mass_fraction = num_ghosts_mass_fraction[2];
                const int ghostcell_dim_0_mass_fraction = ghostcell_dims_mass_fraction[0];
                const int ghostcell_dim_1_mass_fraction = ghostcell_dims_mass_fraction[1];
                
                const int num_ghosts_0_pressure = num_ghosts_pressure[0];
                const int num_ghosts_1_pressure = num_ghosts_pressure[1];
                const int num_ghosts_2_pressure = num_ghosts_pressure[2];
                const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
                const int ghostcell_dim_1_pressure = ghostcell_dims_pressure[1];
                
                const int num_ghosts_0_temperature = num_ghosts_temperature[0];
                const int num_ghosts_1_temperature = num_ghosts_temperature[1];
                const int num_ghosts_2_temperature = num_ghosts_temperature[2];
                const int ghostcell_dim_0_temperature = ghostcell_dims_temperature[0];
                const int ghostcell_dim_1_temperature = ghostcell_dims_temperature[1];
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> num_ghosts_2_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                std::vector<int> ghostcell_dim_1_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                num_ghosts_2_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                ghostcell_dim_1_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    num_ghosts_2_u_qi.push_back(num_ghosts_u_qi[2]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                    ghostcell_dim_1_u_qi.push_back(ghostcell_dims_u_qi[1]);
                }
                
                /*
                 * Initialize cell data for the derivatives and get pointers to the cell data.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_derivative;
                std::vector<double*> der_qi;
                
                if (num_use_derivative > 0)
                {
                    data_derivative = HAMERS_MAKE_SHARED<pdat::CellData<double> >(
                        patch_box, num_use_derivative, hier::IntVector::getZero(d_dim));
                    
                    der_qi.resize(num_use_derivative);
                    for (int qi = 0; qi < num_use_derivative; qi++)
                    {
                        der_qi[qi] = data_derivative->getPointer(qi);
                    }
                }
                
                const hier::IntVector patch_interior_dims = patch_box.numberCells();
                const int patch_interior_dim_0 = patch_interior_dims[0];
                const int patch_interior_dim_1 = patch_interior_dims[1];
                
                const double weight = (dx[1]*dx[2])/(L_y*L_z);
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    int count_derivative = 0;
                    for (int qi = 0; qi < num_quantities; qi++)
                    {
                        if (use_derivative[qi] && derivative_directions[qi] == 0)
                        {
                            derivative_first_order_x->computeDerivative(
                                data_derivative,
                                data_quantities[qi],
                                dx[0],
                                patch_visible_box,
                                count_derivative,
                                component_indices[qi]);
                            
                            count_derivative++;
                        }
                        else if (use_derivative[qi] && derivative_directions[qi] == 1)
                        {
                            derivative_first_order_y->computeDerivative(
                                data_derivative,
                                data_quantities[qi],
                                dx[1],
                                patch_visible_box,
                                count_derivative,
                                component_indices[qi]);
                            
                            count_derivative++;
                        }
                        else if (use_derivative[qi] && derivative_directions[qi] == 2)
                        {
                            derivative_first_order_z->computeDerivative(
                                data_derivative,
                                data_quantities[qi],
                                dx[2],
                                patch_visible_box,
                                count_derivative,
                                component_indices[qi]);
                            
                            count_derivative++;
                        }
                    }
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    const int interior_dim_2 = interior_dims[2];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int idx_lo_2 = index_lo[2];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    const int relative_idx_lo_2 = relative_index_lo[2];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                /*
                                 * Compute the index of the data point and count how many times the data is repeated.
                                 */
                                
                                const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j, idx_lo_2 + k);
                                
                                int n_overlapped = 1;
                                
                                for (hier::BoxContainer::BoxContainerConstIterator iob(
                                        patch_overlapped_visible_boxes.begin());
                                     iob != patch_overlapped_visible_boxes.end();
                                     iob++)
                                {
                                    const hier::Box& patch_overlapped_visible_box = *iob;
                                    
                                    if (patch_overlapped_visible_box.contains(idx_pt))
                                    {
                                        n_overlapped++;
                                    }
                                }
                                
                                /*
                                 * Compute the shear stress component.
                                 */
                                
                                double tau_ij = double(0);
                                
                                const int idx_mass_fraction = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                    (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                        ghostcell_dim_1_mass_fraction;
                                
                                const int idx_pressure = (relative_idx_lo_0 + i + num_ghosts_0_pressure) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure +
                                    (relative_idx_lo_2 + k + num_ghosts_2_pressure)*ghostcell_dim_0_pressure*
                                        ghostcell_dim_1_pressure;
                                
                                const int idx_temperature = (relative_idx_lo_0 + i + num_ghosts_0_temperature) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_temperature)*ghostcell_dim_0_temperature +
                                    (relative_idx_lo_2 + k + num_ghosts_2_temperature)*ghostcell_dim_0_temperature*
                                        ghostcell_dim_1_temperature;
                                
                                if (shear_stress_component_idx == 0)
                                {
                                    std::vector<const double*> Y_ptr;
                                    Y_ptr.resize(d_num_species);
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        Y_ptr[si] = &Y[si][idx_mass_fraction];
                                    }
                                    const double mu = d_equation_of_shear_viscosity_mixing_rules->
                                        getShearViscosity(
                                            &p[idx_pressure],
                                            &T[idx_temperature],
                                            Y_ptr);
                                    const double mu_v = d_equation_of_bulk_viscosity_mixing_rules->
                                        getBulkViscosity(
                                            &p[idx_pressure],
                                            &T[idx_temperature],
                                            Y_ptr);
                                    
                                    const int idx_vel_x_LLL = (relative_idx_lo_0 + (i - 3) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_LL  = (relative_idx_lo_0 + (i - 2) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_L   = (relative_idx_lo_0 + (i - 1) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_R   = (relative_idx_lo_0 + (i + 1) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_RR  = (relative_idx_lo_0 + (i + 2) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_RRR = (relative_idx_lo_0 + (i + 3) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const double dudx = (double(1)/double(60)*(u[idx_vel_x_RRR] - u[idx_vel_x_LLL])
                                        - double(3)/double(20)*(u[idx_vel_x_RR] - u[idx_vel_x_LL])
                                        + double(3)/double(4)*(u[idx_vel_x_R] - u[idx_vel_x_L]))/dx[0];
                                    
                                    const int idx_vel_y_BBB = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j - 3) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_BB  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j - 2) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_B   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j - 1) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_T   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + 1) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_TT  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + 2) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_TTT = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + 3) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const double dvdy = (double(1)/double(60)*(v[idx_vel_y_TTT] - v[idx_vel_y_BBB])
                                        - double(3)/double(20)*(v[idx_vel_y_TT] - v[idx_vel_y_BB])
                                        + double(3)/double(4)*(v[idx_vel_y_T] - v[idx_vel_y_B]))/dx[1];
                                    
                                    const int idx_vel_z_BBB = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k - 3) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_BB  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k - 2) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_B   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k - 1) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_F   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k + 1) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_FF  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k + 2) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_FFF = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k + 3) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const double dwdz = (double(1)/double(60)*(w[idx_vel_z_FFF] - w[idx_vel_z_BBB])
                                        - double(3)/double(20)*(w[idx_vel_z_FF] - w[idx_vel_z_BB])
                                        + double(3)/double(4)*(w[idx_vel_z_F] - w[idx_vel_z_B]))/dx[2];
                                    
                                    tau_ij = (double(4)/double(3)*mu + mu_v)*dudx - (double(2)/double(3)*mu - mu_v)*(dvdy + dwdz);
                                }
                                else if (shear_stress_component_idx == 1)
                                {
                                    std::vector<const double*> Y_ptr;
                                    Y_ptr.resize(d_num_species);
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        Y_ptr[si] = &Y[si][idx_mass_fraction];
                                    }
                                    const double mu = d_equation_of_shear_viscosity_mixing_rules->
                                        getShearViscosity(
                                            &p[idx_pressure],
                                            &T[idx_temperature],
                                            Y_ptr);
                                    
                                    const int idx_vel_x_LLL = (relative_idx_lo_0 + (i - 3) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_LL  = (relative_idx_lo_0 + (i - 2) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_L   = (relative_idx_lo_0 + (i - 1) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_R   = (relative_idx_lo_0 + (i + 1) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_RR  = (relative_idx_lo_0 + (i + 2) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_RRR = (relative_idx_lo_0 + (i + 3) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const double dvdx = (double(1)/double(60)*(v[idx_vel_x_RRR] - v[idx_vel_x_LLL])
                                        - double(3)/double(20)*(v[idx_vel_x_RR] - v[idx_vel_x_LL])
                                        + double(3)/double(4)*(v[idx_vel_x_R] - v[idx_vel_x_L]))/dx[0];
                                    
                                    const int idx_vel_y_BBB = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j - 3) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_BB  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j - 2) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_B   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j - 1) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_T   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + 1) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_TT  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + 2) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_TTT = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + 3) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const double dudy = (double(1)/double(60)*(u[idx_vel_y_TTT] - u[idx_vel_y_BBB])
                                        - double(3)/double(20)*(u[idx_vel_y_TT] - u[idx_vel_y_BB])
                                        + double(3)/double(4)*(u[idx_vel_y_T] - u[idx_vel_y_B]))/dx[1];
                                    
                                    tau_ij = mu*(dudy + dvdx);
                                }
                                else if (shear_stress_component_idx == 2)
                                {
                                    std::vector<const double*> Y_ptr;
                                    Y_ptr.resize(d_num_species);
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        Y_ptr[si] = &Y[si][idx_mass_fraction];
                                    }
                                    const double mu = d_equation_of_shear_viscosity_mixing_rules->
                                        getShearViscosity(
                                            &p[idx_pressure],
                                            &T[idx_temperature],
                                            Y_ptr);
                                    
                                    const int idx_vel_x_LLL = (relative_idx_lo_0 + (i - 3) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_LL  = (relative_idx_lo_0 + (i - 2) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_L   = (relative_idx_lo_0 + (i - 1) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_R   = (relative_idx_lo_0 + (i + 1) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_RR  = (relative_idx_lo_0 + (i + 2) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_RRR = (relative_idx_lo_0 + (i + 3) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const double dwdx = (double(1)/double(60)*(w[idx_vel_x_RRR] - w[idx_vel_x_LLL])
                                        - double(3)/double(20)*(w[idx_vel_x_RR] - w[idx_vel_x_LL])
                                        + double(3)/double(4)*(w[idx_vel_x_R] - w[idx_vel_x_L]))/dx[0];
                                    
                                    const int idx_vel_z_BBB = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k - 3) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_BB  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k - 2) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_B   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k - 1) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_F   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k + 1) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_FF  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k + 2) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_FFF = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k + 3) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const double dudz = (double(1)/double(60)*(u[idx_vel_z_FFF] - u[idx_vel_z_BBB])
                                        - double(3)/double(20)*(u[idx_vel_z_FF] - u[idx_vel_z_BB])
                                        + double(3)/double(4)*(u[idx_vel_z_F] - u[idx_vel_z_B]))/dx[2];
                                    
                                    tau_ij = mu*(dudz + dwdx);
                                }
                                else if (shear_stress_component_idx == 3)
                                {
                                    std::vector<const double*> Y_ptr;
                                    Y_ptr.resize(d_num_species);
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        Y_ptr[si] = &Y[si][idx_mass_fraction];
                                    }
                                    const double mu = d_equation_of_shear_viscosity_mixing_rules->
                                        getShearViscosity(
                                            &p[idx_pressure],
                                            &T[idx_temperature],
                                            Y_ptr);
                                    const double mu_v = d_equation_of_bulk_viscosity_mixing_rules->
                                        getBulkViscosity(
                                            &p[idx_pressure],
                                            &T[idx_temperature],
                                            Y_ptr);
                                    
                                    const int idx_vel_x_LLL = (relative_idx_lo_0 + (i - 3) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_LL  = (relative_idx_lo_0 + (i - 2) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_L   = (relative_idx_lo_0 + (i - 1) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_R   = (relative_idx_lo_0 + (i + 1) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_RR  = (relative_idx_lo_0 + (i + 2) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_RRR = (relative_idx_lo_0 + (i + 3) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const double dudx = (double(1)/double(60)*(u[idx_vel_x_RRR] - u[idx_vel_x_LLL])
                                        - double(3)/double(20)*(u[idx_vel_x_RR] - u[idx_vel_x_LL])
                                        + double(3)/double(4)*(u[idx_vel_x_R] - u[idx_vel_x_L]))/dx[0];
                                    
                                    const int idx_vel_y_BBB = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j - 3) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_BB  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j - 2) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_B   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j - 1) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_T   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + 1) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_TT  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + 2) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_TTT = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + 3) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const double dvdy = (double(1)/double(60)*(v[idx_vel_y_TTT] - v[idx_vel_y_BBB])
                                        - double(3)/double(20)*(v[idx_vel_y_TT] - v[idx_vel_y_BB])
                                        + double(3)/double(4)*(v[idx_vel_y_T] - v[idx_vel_y_B]))/dx[1];
                                    
                                    const int idx_vel_z_BBB = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k - 3) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_BB  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k - 2) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_B   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k - 1) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_F   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k + 1) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_FF  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k + 2) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_FFF = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k + 3) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const double dwdz = (double(1)/double(60)*(w[idx_vel_z_FFF] - w[idx_vel_z_BBB])
                                        - double(3)/double(20)*(w[idx_vel_z_FF] - w[idx_vel_z_BB])
                                        + double(3)/double(4)*(w[idx_vel_z_F] - w[idx_vel_z_B]))/dx[2];
                                    
                                    tau_ij = (double(4)/double(3)*mu + mu_v)*dvdy - (double(2)/double(3)*mu - mu_v)*(dudx + dwdz);
                                }
                                else if (shear_stress_component_idx == 4)
                                {
                                    std::vector<const double*> Y_ptr;
                                    Y_ptr.resize(d_num_species);
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        Y_ptr[si] = &Y[si][idx_mass_fraction];
                                    }
                                    const double mu = d_equation_of_shear_viscosity_mixing_rules->
                                        getShearViscosity(
                                            &p[idx_pressure],
                                            &T[idx_temperature],
                                            Y_ptr);
                                    
                                    const int idx_vel_y_BBB = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j - 3) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_BB  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j - 2) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_B   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j - 1) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_T   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + 1) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_TT  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + 2) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_TTT = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + 3) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const double dwdy = (double(1)/double(60)*(w[idx_vel_y_TTT] - w[idx_vel_y_BBB])
                                        - double(3)/double(20)*(w[idx_vel_y_TT] - w[idx_vel_y_BB])
                                        + double(3)/double(4)*(w[idx_vel_y_T] - w[idx_vel_y_B]))/dx[1];
                                    
                                    const int idx_vel_z_BBB = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k - 3) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_BB  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k - 2) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_B   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k - 1) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_F   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k + 1) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_FF  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k + 2) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_FFF = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k + 3) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const double dvdz = (double(1)/double(60)*(v[idx_vel_z_FFF] - v[idx_vel_z_BBB])
                                        - double(3)/double(20)*(v[idx_vel_z_FF] - v[idx_vel_z_BB])
                                        + double(3)/double(4)*(v[idx_vel_z_F] - v[idx_vel_z_B]))/dx[2];
                                    
                                    tau_ij = mu*(dvdz + dwdy);
                                }
                                else if (shear_stress_component_idx == 5)
                                {
                                    std::vector<const double*> Y_ptr;
                                    Y_ptr.resize(d_num_species);
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        Y_ptr[si] = &Y[si][idx_mass_fraction];
                                    }
                                    const double mu = d_equation_of_shear_viscosity_mixing_rules->
                                        getShearViscosity(
                                            &p[idx_pressure],
                                            &T[idx_temperature],
                                            Y_ptr);
                                    const double mu_v = d_equation_of_bulk_viscosity_mixing_rules->
                                        getBulkViscosity(
                                            &p[idx_pressure],
                                            &T[idx_temperature],
                                            Y_ptr);
                                    
                                    const int idx_vel_x_LLL = (relative_idx_lo_0 + (i - 3) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_LL  = (relative_idx_lo_0 + (i - 2) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_L   = (relative_idx_lo_0 + (i - 1) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_R   = (relative_idx_lo_0 + (i + 1) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_RR  = (relative_idx_lo_0 + (i + 2) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_x_RRR = (relative_idx_lo_0 + (i + 3) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const double dudx = (double(1)/double(60)*(u[idx_vel_x_RRR] - u[idx_vel_x_LLL])
                                        - double(3)/double(20)*(u[idx_vel_x_RR] - u[idx_vel_x_LL])
                                        + double(3)/double(4)*(u[idx_vel_x_R] - u[idx_vel_x_L]))/dx[0];
                                    
                                    const int idx_vel_y_BBB = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j - 3) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_BB  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j - 2) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_B   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j - 1) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_T   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + 1) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_TT  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + 2) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_y_TTT = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + 3) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const double dvdy = (double(1)/double(60)*(v[idx_vel_y_TTT] - v[idx_vel_y_BBB])
                                        - double(3)/double(20)*(v[idx_vel_y_TT] - v[idx_vel_y_BB])
                                        + double(3)/double(4)*(v[idx_vel_y_T] - v[idx_vel_y_B]))/dx[1];
                                    
                                    const int idx_vel_z_BBB = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k - 3) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_BB  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k - 2) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_B   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k - 1) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_F   = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k + 1) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_FF  = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k + 2) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const int idx_vel_z_FFF = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                        (relative_idx_lo_2 + (k + 3) + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                            ghostcell_dim_1_velocity;
                                    
                                    const double dwdz = (double(1)/double(60)*(w[idx_vel_z_FFF] - w[idx_vel_z_BBB])
                                        - double(3)/double(20)*(w[idx_vel_z_FF] - w[idx_vel_z_BB])
                                        + double(3)/double(4)*(w[idx_vel_z_F] - w[idx_vel_z_B]))/dx[2];
                                    
                                    tau_ij = (double(4)/double(3)*mu + mu_v)*dwdz - (double(2)/double(3)*mu - mu_v)*(dudx + dvdy);
                                }
                                else
                                {
                                    TBOX_ERROR(d_object_name
                                        << ": "
                                        << "Cannot compute shear stress component for two-dimensional problem!\n"
                                        << "component_idx = " << shear_stress_component_idx << " given!\n"
                                        << std::endl);
                                }
                                
                                /*
                                 * Compute the linear index and the data to add.
                                 */
                                
                                double avg = double(1);
                                
                                count_derivative = 0;
                                for (int qi = 0; qi < num_quantities; qi++)
                                {
                                    if (use_reciprocal[qi])
                                    {
                                        if (use_derivative[qi])
                                        {
                                            const int idx_der = (relative_idx_lo_0 + i) +
                                                (relative_idx_lo_1 + j)*patch_interior_dim_0 +
                                                (relative_idx_lo_2 + k)*patch_interior_dim_0*
                                                    patch_interior_dim_1;
                                            
                                            avg /= der_qi[count_derivative][idx_der];
                                            
                                            count_derivative++;
                                        }
                                        else
                                        {
                                            const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                                (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                                (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                    ghostcell_dim_1_u_qi[qi];
                                            
                                            avg /= u_qi[qi][idx_qi];
                                        }
                                    }
                                    else
                                    {
                                        if (use_derivative[qi])
                                        {
                                            const int idx_der = (relative_idx_lo_0 + i) +
                                                (relative_idx_lo_1 + j)*patch_interior_dim_0 +
                                                (relative_idx_lo_2 + k)*patch_interior_dim_0*
                                                    patch_interior_dim_1;
                                            
                                            avg *= der_qi[count_derivative][idx_der];
                                            
                                            count_derivative++;
                                        }
                                        else
                                        {
                                            const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                                (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                                (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                    ghostcell_dim_1_u_qi[qi];
                                            
                                            avg *= u_qi[qi][idx_qi];
                                        }
                                    }
                                }
                                
                                /*
                                 * Add the data.
                                 */
                                
                                for (int ii = 0; ii < ratio_to_finest_level_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratio_to_finest_level_0 + ii;
                                    
                                    avg_local[idx_fine] += (avg*tau_ij*weight/((double) n_overlapped));
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        mpi.Allreduce(
            avg_local,
            avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(avg_local);
    }
    
    return averaged_quantity;
}


/*
 * Get gravity vector from the flow model database.
 */
std::vector<double>
RTIRMIStatisticsUtilities::getGravityVector() const
{
    HAMERS_SHARED_PTR<FlowModel> d_flow_model_tmp = d_flow_model.lock();
    
    const HAMERS_SHARED_PTR<tbox::Database>& flow_model_db = d_flow_model_tmp->getFlowModelDatabase();
    
    std::vector<double> gravity(d_dim.getValue(), double(0));
    
    bool has_source_terms = false;
    if (flow_model_db->keyExists("has_source_terms"))
    {
        has_source_terms = flow_model_db->getBool("has_source_terms");
    }
    else if (flow_model_db->keyExists("d_has_source_terms"))
    {
        has_source_terms = flow_model_db->getBool("d_has_source_terms");
    }
    
    bool has_gravity = false;
    if (has_source_terms)
    {
        HAMERS_SHARED_PTR<tbox::Database> source_terms_db;
        
        if (flow_model_db->keyExists("Source_terms"))
        {
            source_terms_db = flow_model_db->getDatabase("Source_terms");
        }
        else if (flow_model_db->keyExists("d_source_terms"))
        {
            source_terms_db = flow_model_db->getDatabase("d_source_terms");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "No key 'Source_terms'/'d_source_terms' found in data for flow model."
                << std::endl);
        }
        
        if (source_terms_db->keyExists("has_gravity"))
        {
            has_gravity = source_terms_db->getBool("has_gravity");
            if (has_gravity)
            {
                if (source_terms_db->keyExists("gravity"))
                {
                    source_terms_db->getVector("gravity", gravity);
                }
                else
                {
                    TBOX_ERROR(d_object_name
                        << ": "
                        << "No key 'gravity' found in data for source terms."
                        << std::endl);
                }
            }
        }
        else if (source_terms_db->keyExists("d_has_gravity"))
        {
            has_gravity = source_terms_db->getBool("d_has_gravity");
            if (has_gravity)
            {
                if (source_terms_db->keyExists("d_gravity"))
                {
                    source_terms_db->getVector("d_gravity", gravity);
                }
                else
                {
                    TBOX_ERROR(d_object_name
                        << ": "
                        << "No key 'd_gravity' found in data for source terms."
                        << std::endl);
                }
            }
        }
        
        if (has_gravity)
        {
            TBOX_ASSERT(gravity.size() == d_dim.getValue())
        }
    }
    
    return gravity;
}


/*
 * Compute the one-dimensional derivative given a vector.
 */
std::vector<double>
RTIRMIStatisticsUtilities::computeDerivativeOfVector1D(
    const std::vector<double> quantity_vector,
    const double dx) const
{
    TBOX_ASSERT(d_num_ghosts_derivative == 3);
    
    const int vector_length = quantity_vector.size();
    
    std::vector<double> derivative;
    derivative.resize(vector_length);
    
    const double* u = quantity_vector.data();
    double* dudx = derivative.data();
    
    // HARD CODE TO BE SIXTH ORDER CENTRAL SCHEME FOR DIFFERENTIATION.
    for (int i = 3; i < vector_length - 3; i++)
    {
        // Compute linear indices.
        const int idx     = i;
        
        const int idx_LLL = i - 3;
        const int idx_LL  = i - 2;
        const int idx_L   = i - 1;
        const int idx_R   = i + 1;
        const int idx_RR  = i + 2;
        const int idx_RRR = i + 3;
        
        dudx[idx] = (double(1)/double(60)*(u[idx_RRR] - u[idx_LLL])
            - double(3)/double(20)*(u[idx_RR] - u[idx_LL])
            + double(3)/double(4)*(u[idx_R] - u[idx_L]))/dx;
    }
    
    for (int i = 0; i < 3; i++)
    {
        dudx[i]                     = dudx[3];
        dudx[vector_length - i - 1] = dudx[vector_length - 4];
    }
    
    return derivative;
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
            else if (statistical_quantity_key == "MIXING_WIDTH_X_MOL_F")
            {
                f_out << "\t" << "MIXING_WIDTH_X_MOL_F ";
            }
            else if (statistical_quantity_key == "MIXING_WIDTH_X_VOL_F")
            {
                f_out << "\t" << "MIXING_WIDTH_X_VOL_F ";
            }
            else if (statistical_quantity_key == "INTERFACE_MIN_X")
            {
                f_out << "\t" << "INTERFACE_MIN_X      ";
            }
            else if (statistical_quantity_key == "INTERFACE_MIN_X_MOL_F")
            {
                f_out << "\t" << "INTERFACE_MIN_X_MOL_F";
            }
            else if (statistical_quantity_key == "INTERFACE_MIN_X_VOL_F")
            {
                f_out << "\t" << "INTERFACE_MIN_X_VOL_F";
            }
            else if (statistical_quantity_key == "INTERFACE_MAX_X")
            {
                f_out << "\t" << "INTERFACE_MAX_X      ";
            }
            else if (statistical_quantity_key == "INTERFACE_MAX_X_MOL_F")
            {
                f_out << "\t" << "INTERFACE_MAX_X_MOL_F";
            }
            else if (statistical_quantity_key == "INTERFACE_MAX_X_VOL_F")
            {
                f_out << "\t" << "INTERFACE_MAX_X_VOL_F";
            }
            else if (statistical_quantity_key == "MIXEDNESS_X")
            {
                f_out << "\t" << "MIXEDNESS_X          ";
            }
            else if (statistical_quantity_key == "MIXEDNESS_X_MOL_F")
            {
                f_out << "\t" << "MIXEDNESS_X_MOL_F    ";
            }
            else if (statistical_quantity_key == "MIXEDNESS_X_VOL_F")
            {
                f_out << "\t" << "MIXEDNESS_X_VOL_F    ";
            }
            else if (statistical_quantity_key == "DENSITY_MEAN_INHOMO_X")
            {
                f_out << "\t" << "DENSITY_MEAN_INHOMO_X";
            }
            else if (statistical_quantity_key == "SP_VOL_MEAN_INHOMO_X")
            {
                f_out << "\t" << "SP_VOL_MEAN_INHOMO_X ";
            }
            else if (statistical_quantity_key == "PRESS_MEAN_INHOMO_X")
            {
                f_out << "\t" << "PRESS_MEAN_INHOMO_X  ";
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
            else if (statistical_quantity_key == "ENSTROPHY_MEAN_INHO_X")
            {
                f_out << "\t" << "ENSTROPHY_MEAN_INHO_X";
            }
            else if (statistical_quantity_key == "S_DIS_RAT_INT")
            {
                f_out << "\t" << "S_DIS_RAT_INT        ";
            }
            else if (statistical_quantity_key == "S_DIS_RAT_MEAN_INHO_X")
            {
                f_out << "\t" << "S_DIS_RAT_MEAN_INHO_X";
            }
            else if (statistical_quantity_key == "RE_W_INHOMO_X")
            {
                f_out << "\t" << "RE_W_INHOMO_X        ";
            }
            else if (statistical_quantity_key == "MA_T_INHOMO_X")
            {
                f_out << "\t" << "MA_T_INHOMO_X        ";
            }
            else if (statistical_quantity_key == "EFF_AT_INHOMO_X")
            {
                f_out << "\t" << "EFF_AT_INHOMO_X      ";
            }
            else if (statistical_quantity_key == "BARO_TQ_Z_INT")
            {
                f_out << "\t" << "BARO_TQ_Z_INT        ";
            }
            else if (statistical_quantity_key == "BARO_TQ_Z_MEAN_INHO_X")
            {
                f_out << "\t" << "BARO_TQ_Z_MEAN_INHO_X";
            }
            else if (statistical_quantity_key == "BARO_TQ_X_INT")
            {
                f_out << "\t" << "BARO_TQ_X_INT        ";
            }
            else if (statistical_quantity_key == "BARO_TQ_X_MEAN_INHO_X")
            {
                f_out << "\t" << "BARO_TQ_X_MEAN_INHO_X";
            }
            else if (statistical_quantity_key == "BARO_TQ_Y_INT")
            {
                f_out << "\t" << "BARO_TQ_Y_INT        ";
            }
            else if (statistical_quantity_key == "BARO_TQ_Y_MEAN_INHO_X")
            {
                f_out << "\t" << "BARO_TQ_Y_MEAN_INHO_X";
            }
            else if (statistical_quantity_key == "TKE_PROD1_INT")
            {
                f_out << "\t" << "TKE_PROD1_INT        ";
            }
            else if (statistical_quantity_key == "TKE_PROD1_MEAN_INHO_X")
            {
                f_out << "\t" << "TKE_PROD1_MEAN_INHO_X";
            }
            else if (statistical_quantity_key == "TKE_PROD2_INT")
            {
                f_out << "\t" << "TKE_PROD2_INT        ";
            }
            else if (statistical_quantity_key == "TKE_PROD2_MEAN_INHO_X")
            {
                f_out << "\t" << "TKE_PROD2_MEAN_INHO_X";
            }
            else if (statistical_quantity_key == "TKE_PROD3_INT")
            {
                f_out << "\t" << "TKE_PROD3_INT        ";
            }
            else if (statistical_quantity_key == "TKE_PROD3_MEAN_INHO_X")
            {
                f_out << "\t" << "TKE_PROD3_MEAN_INHO_X";
            }
            else if (statistical_quantity_key == "TKE_DISS_INT")
            {
                f_out << "\t" << "TKE_DISS_INT         ";
            }
            else if (statistical_quantity_key == "TKE_DISS_MEAN_INHO_X")
            {
                f_out << "\t" << "TKE_DISS_MEAN_INHO_X ";
            }
            else if (statistical_quantity_key == "SHEAR_VIS_MEAN_INHO_X")
            {
                f_out << "\t" << "SHEAR_VIS_MEAN_INHO_X";
            }
            else if (statistical_quantity_key == "MASS_DIFF_MEAN_INHO_X")
            {
                f_out << "\t" << "MASS_DIFF_MEAN_INHO_X";
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
        else if (statistical_quantity_key == "MOLE_FRACTION_AVG_SP")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->X_0_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedMoleFractionWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
        }
        else if (statistical_quantity_key == "VOLUME_FRACTION_AVG_SP")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->Z_0_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedVolumeFractionWithHomogeneityInYDirectionOrInYZPlane(
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
        else if (statistical_quantity_key == "MOLE_FRACTION_VAR_SP")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->X_0_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedMoleFractionWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
            
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->X_0_sq_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeMoleFractionVarianceWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
        }
        else if (statistical_quantity_key == "VOLUME_FRACTION_VAR_SP")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->Z_0_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedVolumeFractionWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
            
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->Z_0_sq_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeVolumeFractionVarianceWithHomogeneityInYDirectionOrInYZPlane(
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
        else if (statistical_quantity_key == "PRESSURE_AVG_SP")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->p_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedPressureWithHomogeneityInYDirectionOrInYZPlane(
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
        else if (statistical_quantity_key == "BARO_TQ_Z_RMS_SP")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->B3_sq_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedSquaredBaroclinicTorqueZWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
        }
        else if (statistical_quantity_key == "BARO_TQ_X_RMS_SP")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->B1_sq_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedSquaredBaroclinicTorqueXWithHomogeneityInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
        }
        else if (statistical_quantity_key == "BARO_TQ_Y_RMS_SP")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->B2_sq_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedSquaredBaroclinicTorqueYWithHomogeneityInYZPlane(
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
        else if (statistical_quantity_key == "MIXING_WIDTH_X_MOL_F")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->X_0_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedMoleFractionWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
        }
        else if (statistical_quantity_key == "MIXING_WIDTH_X_VOL_F")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->Z_0_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedVolumeFractionWithHomogeneityInYDirectionOrInYZPlane(
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
        else if (statistical_quantity_key == "INTERFACE_MIN_X_MOL_F")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->X_0_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedMoleFractionWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
        }
        else if (statistical_quantity_key == "INTERFACE_MIN_X_VOL_F")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->Z_0_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedVolumeFractionWithHomogeneityInYDirectionOrInYZPlane(
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
        else if (statistical_quantity_key == "INTERFACE_MAX_X_MOL_F")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->X_0_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedMoleFractionWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
        }
        else if (statistical_quantity_key == "INTERFACE_MAX_X_VOL_F")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->Z_0_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedVolumeFractionWithHomogeneityInYDirectionOrInYZPlane(
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
        else if (statistical_quantity_key == "MIXEDNESS_X_MOL_F")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->X_0_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedMoleFractionWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
            
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->X_0_X_1_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeMoleFractionProductWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
        }
        else if (statistical_quantity_key == "MIXEDNESS_X_VOL_F")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->Z_0_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedVolumeFractionWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
            
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->Z_0_Z_1_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeVolumeFractionProductWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
        }
        else if (statistical_quantity_key == "DENSITY_MEAN_INHOMO_X")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->rho_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedDensityWithHomogeneityInYDirectionOrInYZPlane(
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
        else if (statistical_quantity_key == "SP_VOL_MEAN_INHOMO_X")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->rho_inv_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedSpecificVolumeWithHomogeneityInYDirectionOrInYZPlane(
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
        else if (statistical_quantity_key == "PRESS_MEAN_INHOMO_X")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->p_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedPressureWithHomogeneityInYDirectionOrInYZPlane(
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
        else if (statistical_quantity_key == "ENSTROPHY_MEAN_INHO_X")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->Omega_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedEnstrophyWithHomogeneityInYDirectionOrInYZPlane(
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
        else if (statistical_quantity_key == "S_DIS_RAT_INT")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->chi_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedScalarDissipationRateWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
        }
        else if (statistical_quantity_key == "S_DIS_RAT_MEAN_INHO_X")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->chi_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedScalarDissipationRateWithHomogeneityInYDirectionOrInYZPlane(
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
        else if (statistical_quantity_key == "EFF_AT_INHOMO_X")
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
            
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->Y_0_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedMassFractionWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
        }
        else if (statistical_quantity_key == "BARO_TQ_Z_INT")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->B3_sq_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedSquaredBaroclinicTorqueZWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
        }
        else if (statistical_quantity_key == "BARO_TQ_Z_MEAN_INHO_X")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->B3_sq_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedSquaredBaroclinicTorqueZWithHomogeneityInYDirectionOrInYZPlane(
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
        else if (statistical_quantity_key == "BARO_TQ_X_INT")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->B1_sq_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedSquaredBaroclinicTorqueXWithHomogeneityInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
        }
        else if (statistical_quantity_key == "BARO_TQ_X_MEAN_INHO_X")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->B1_sq_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedSquaredBaroclinicTorqueXWithHomogeneityInYZPlane(
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
        else if (statistical_quantity_key == "BARO_TQ_Y_INT")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->B2_sq_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedSquaredBaroclinicTorqueYWithHomogeneityInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
        }
        else if (statistical_quantity_key == "BARO_TQ_Y_MEAN_INHO_X")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->B2_sq_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedSquaredBaroclinicTorqueYWithHomogeneityInYZPlane(
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
        else if (statistical_quantity_key == "TKE_PROD1_INT")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->ddx_p_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedXDerivativeOfPressureWithHomogeneityInYDirectionOrInYZPlane(
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
            
            if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
            {
                if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->ddy_p_avg_computed))
                {
                    rti_rmi_statistics_utilities->
                        computeAveragedYDerivativeOfPressureWithHomogeneityInYDirectionOrInYZPlane(
                            patch_hierarchy,
                            data_context);
                }
                
                if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->v_avg_computed))
                {
                    rti_rmi_statistics_utilities->
                        computeAveragedVelocityYWithHomogeneityInYDirectionOrInYZPlane(
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
                if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->ddz_p_avg_computed))
                {
                    rti_rmi_statistics_utilities->
                        computeAveragedZDerivativeOfPressureWithHomogeneityInYZPlane(
                            patch_hierarchy,
                            data_context);
                }
                
                if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->w_avg_computed))
                {
                    rti_rmi_statistics_utilities->
                        computeAveragedVelocityZWithHomogeneityInYZPlane(
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
            
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->rho_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedDensityWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
        }
        else if (statistical_quantity_key == "TKE_PROD1_MEAN_INHO_X")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->ddx_p_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedXDerivativeOfPressureWithHomogeneityInYDirectionOrInYZPlane(
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
            
            if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
            {
                if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->ddy_p_avg_computed))
                {
                    rti_rmi_statistics_utilities->
                        computeAveragedYDerivativeOfPressureWithHomogeneityInYDirectionOrInYZPlane(
                            patch_hierarchy,
                            data_context);
                }
                
                if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->v_avg_computed))
                {
                    rti_rmi_statistics_utilities->
                        computeAveragedVelocityYWithHomogeneityInYDirectionOrInYZPlane(
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
                if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->ddz_p_avg_computed))
                {
                    rti_rmi_statistics_utilities->
                        computeAveragedZDerivativeOfPressureWithHomogeneityInYZPlane(
                            patch_hierarchy,
                            data_context);
                }
                
                if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->w_avg_computed))
                {
                    rti_rmi_statistics_utilities->
                        computeAveragedVelocityZWithHomogeneityInYZPlane(
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
            
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->rho_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedDensityWithHomogeneityInYDirectionOrInYZPlane(
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
        else if (statistical_quantity_key == "TKE_PROD2_INT")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->ddx_tau11_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedXDerivativeNormalShearStressXWithHomogeneityInYDirectionOrInYZPlane(
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
        else if (statistical_quantity_key == "TKE_PROD2_MEAN_INHO_X")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->ddx_tau11_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedXDerivativeNormalShearStressXWithHomogeneityInYDirectionOrInYZPlane(
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
        else if (statistical_quantity_key == "TKE_PROD3_INT")
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
        else if (statistical_quantity_key == "TKE_PROD3_MEAN_INHO_X")
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
        else if (statistical_quantity_key == "TKE_DISS_INT")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->TKE_diss_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedQuantiitesForTKEDissipationWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
        }
        else if (statistical_quantity_key == "TKE_DISS_MEAN_INHO_X")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->TKE_diss_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedQuantiitesForTKEDissipationWithHomogeneityInYDirectionOrInYZPlane(
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
        else if (statistical_quantity_key == "SHEAR_VIS_MEAN_INHO_X")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->mu_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedShearViscosityWithHomogeneityInYDirectionOrInYZPlane(
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
        else if (statistical_quantity_key == "MASS_DIFF_MEAN_INHO_X")
        {
            if (!(rti_rmi_statistics_utilities->d_ensemble_statistics->D_avg_computed))
            {
                rti_rmi_statistics_utilities->
                    computeAveragedMassDiffusivityWithHomogeneityInYDirectionOrInYZPlane(
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
                    output_time);
        }
        else if (statistical_quantity_key == "MOLE_FRACTION_AVG_SP")
        {
            rti_rmi_statistics_utilities->
                outputSpatialProfileEnsembleAveragedMoleFractionWithHomogeneityInYDirectionOrInYZPlane(
                    "X_avg.dat",
                    patch_hierarchy,
                    output_time);
        }
        else if (statistical_quantity_key == "VOLUME_FRACTION_AVG_SP")
        {
            rti_rmi_statistics_utilities->
                outputSpatialProfileEnsembleAveragedVolumeFractionWithHomogeneityInYDirectionOrInYZPlane(
                    "Z_avg.dat",
                    patch_hierarchy,
                    output_time);
        }
        else if (statistical_quantity_key == "MASS_FRACTION_VAR_SP")
        {
            rti_rmi_statistics_utilities->
                outputSpatialProfileEnsembleMassFractionVarianceWithHomogeneityInYDirectionOrInYZPlane(
                    "Y_var.dat",
                    patch_hierarchy,
                    output_time);
        }
        else if (statistical_quantity_key == "MOLE_FRACTION_VAR_SP")
        {
            rti_rmi_statistics_utilities->
                outputSpatialProfileEnsembleMoleFractionVarianceWithHomogeneityInYDirectionOrInYZPlane(
                    "X_var.dat",
                    patch_hierarchy,
                    output_time);
        }
        else if (statistical_quantity_key == "VOLUME_FRACTION_VAR_SP")
        {
            rti_rmi_statistics_utilities->
                outputSpatialProfileEnsembleVolumeFractionVarianceWithHomogeneityInYDirectionOrInYZPlane(
                    "Z_var.dat",
                    patch_hierarchy,
                    output_time);
        }
        else if (statistical_quantity_key == "DENSITY_AVG_SP")
        {
            rti_rmi_statistics_utilities->
                outputSpatialProfileEnsembleAveragedDensityWithHomogeneityInYDirectionOrInYZPlane(
                    "rho_avg.dat",
                    patch_hierarchy,
                    output_time);
        }
        else if (statistical_quantity_key == "DENSITY_VAR_SP")
        {
            rti_rmi_statistics_utilities->
                outputSpatialProfileEnsembleDensityVarianceWithHomogeneityInYDirectionOrInYZPlane(
                    "rho_var.dat",
                    patch_hierarchy,
                    output_time);
        }
        else if (statistical_quantity_key == "SPECIFIC_VOLUME_AVG_SP")
        {
            rti_rmi_statistics_utilities->
                outputSpatialProfileEnsembleAveragedSpecificVolumeWithHomogeneityInYDirectionOrInYZPlane(
                    "rho_inv_avg.dat",
                    patch_hierarchy,
                    output_time);
        }
        else if (statistical_quantity_key == "DENSITY_SPEC_VOL_COV_SP")
        {
            rti_rmi_statistics_utilities->
                outputSpatialProfileEnsembleAveragedDensitySpecificVolumeCovarianceWithHomogeneityInYDirectionOrInYZPlane(
                    "b.dat",
                    patch_hierarchy,
                    output_time);
        }
        else if (statistical_quantity_key == "PRESSURE_AVG_SP")
        {
            rti_rmi_statistics_utilities->
                outputSpatialProfileEnsembleAveragedPressureWithHomogeneityInYDirectionOrInYZPlane(
                    "p_avg.dat",
                    patch_hierarchy,
                    output_time);
        }
        else if (statistical_quantity_key == "VELOCITY_X_AVG_SP")
        {
            rti_rmi_statistics_utilities->
                outputSpatialProfileEnsembleAveragedVelocityXWithHomogeneityInYDirectionOrInYZPlane(
                    "u_avg.dat",
                    patch_hierarchy,
                    output_time);
        }
        else if (statistical_quantity_key == "VELOCITY_Y_AVG_SP")
        {
            rti_rmi_statistics_utilities->
                outputSpatialProfileEnsembleAveragedVelocityYWithHomogeneityInYDirectionOrInYZPlane(
                    "v_avg.dat",
                    patch_hierarchy,
                    output_time);
        }
        else if (statistical_quantity_key == "VELOCITY_Z_AVG_SP")
        {
            rti_rmi_statistics_utilities->
                outputSpatialProfileEnsembleAveragedVelocityZWithHomogeneityInYZPlane(
                    "w_avg.dat",
                    patch_hierarchy,
                    output_time);
        }
        else if (statistical_quantity_key == "MOMENTUM_X_AVG_SP")
        {
            rti_rmi_statistics_utilities->
                outputSpatialProfileEnsembleAveragedMomentumXWithHomogeneityInYDirectionOrInYZPlane(
                    "rho_u_avg.dat",
                    patch_hierarchy,
                    output_time);
        }
        else if (statistical_quantity_key == "MOMENTUM_Y_AVG_SP")
        {
            rti_rmi_statistics_utilities->
                outputSpatialProfileEnsembleAveragedMomentumYWithHomogeneityInYDirectionOrInYZPlane(
                    "rho_v_avg.dat",
                    patch_hierarchy,
                    output_time);
        }
        else if (statistical_quantity_key == "MOMENTUM_Z_AVG_SP")
        {
            rti_rmi_statistics_utilities->
                outputSpatialProfileEnsembleAveragedMomentumZWithHomogeneitynYZPlane(
                    "rho_w_avg.dat",
                    patch_hierarchy,
                    output_time);
        }
        else if (statistical_quantity_key == "TURB_MASS_FLUX_X_SP")
        {
            rti_rmi_statistics_utilities->
                outputSpatialProfileEnsembleTurbulentMassFluxXWithHomogeneityInYDirectionOrInYZPlane(
                    "rho_a1.dat",
                    patch_hierarchy,
                    output_time);
        }
        else if (statistical_quantity_key == "R11_SP")
        {
            rti_rmi_statistics_utilities->
                outputSpatialProfileEnsembleReynoldsNormalStressXWithHomogeneityInYDirectionOrInYZPlane(
                    "R11.dat",
                    patch_hierarchy,
                    output_time);
        }
        else if (statistical_quantity_key == "R22_SP")
        {
            rti_rmi_statistics_utilities->
                outputSpatialProfileEnsembleReynoldsNormalStressYWithHomogeneityInYDirectionOrInYZPlane(
                    "R22.dat",
                    patch_hierarchy,
                    output_time);
        }
        else if (statistical_quantity_key == "R33_SP")
        {
            rti_rmi_statistics_utilities->
                outputSpatialProfileEnsembleReynoldsNormalStressZWithHomogeneityInYZPlane(
                    "R33.dat",
                    patch_hierarchy,
                    output_time);
        }
        else if (statistical_quantity_key == "ENSTROPHY_AVG_SP")
        {
            rti_rmi_statistics_utilities->
                outputSpatialProfileEnsembleAveragedEnstrophyWithHomogeneityInYDirectionOrInYZPlane(
                    "Omega_avg.dat",
                    patch_hierarchy,
                    output_time);
        }
        else if (statistical_quantity_key == "SCAL_DISS_RAT_AVG_SP")
        {
            rti_rmi_statistics_utilities->
                outputSpatialProfileEnsembleAveragedScalarDissipationRateWithHomogeneityInYDirectionOrInYZPlane(
                    "chi_avg.dat",
                    patch_hierarchy,
                    output_time);
        }
        else if (statistical_quantity_key == "BARO_TQ_Z_RMS_SP")
        {
            rti_rmi_statistics_utilities->
                outputSpatialProfileEnsembleRMSBaroclinicTorqueZWithHomogeneityInYDirectionOrInYZPlane(
                    "B3_rms.dat",
                    patch_hierarchy,
                    output_time);
        }
        else if (statistical_quantity_key == "BARO_TQ_X_RMS_SP")
        {
            rti_rmi_statistics_utilities->
                outputSpatialProfileEnsembleRMSBaroclinicTorqueXWithHomogeneityInYZPlane(
                    "B1_rms.dat",
                    patch_hierarchy,
                    output_time);
        }
        else if (statistical_quantity_key == "BARO_TQ_Y_RMS_SP")
        {
            rti_rmi_statistics_utilities->
                outputSpatialProfileEnsembleRMSBaroclinicTorqueYWithHomogeneityInYZPlane(
                    "B2_rms.dat",
                    patch_hierarchy,
                    output_time);
        }
        // Non-spatial profiles.
        else if (statistical_quantity_key == "MIXING_WIDTH_X")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleMixingWidthInXDirection(
                    stat_dump_filename,
                    patch_hierarchy);
        }
        else if (statistical_quantity_key == "MIXING_WIDTH_X_MOL_F")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleMixingWidthInXDirectionWithMoleFractions(
                    stat_dump_filename,
                    patch_hierarchy);
        }
        else if (statistical_quantity_key == "MIXING_WIDTH_X_VOL_F")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleMixingWidthInXDirectionWithVolumeFractions(
                    stat_dump_filename,
                    patch_hierarchy);
        }
        else if (statistical_quantity_key == "INTERFACE_MIN_X")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleInterfaceMinInXDirection(
                    stat_dump_filename,
                    patch_hierarchy);
        }
        else if (statistical_quantity_key == "INTERFACE_MIN_X_MOL_F")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleInterfaceMinInXDirectionWithMoleFractions(
                    stat_dump_filename,
                    patch_hierarchy);
        }
        else if (statistical_quantity_key == "INTERFACE_MIN_X_VOL_F")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleInterfaceMinInXDirectionWithVolumeFractions(
                    stat_dump_filename,
                    patch_hierarchy);
        }
        else if (statistical_quantity_key == "INTERFACE_MAX_X")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleInterfaceMaxInXDirection(
                    stat_dump_filename,
                    patch_hierarchy);
        }
        else if (statistical_quantity_key == "INTERFACE_MAX_X_MOL_F")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleInterfaceMaxInXDirectionWithMoleFractions(
                    stat_dump_filename,
                    patch_hierarchy);
        }
        else if (statistical_quantity_key == "INTERFACE_MAX_X_VOL_F")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleInterfaceMaxInXDirectionWithVolumeFractions(
                    stat_dump_filename,
                    patch_hierarchy);
        }
        else if (statistical_quantity_key == "MIXEDNESS_X")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleMixednessInXDirection(
                    stat_dump_filename,
                    patch_hierarchy);
        }
        else if (statistical_quantity_key == "MIXEDNESS_X_MOL_F")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleMixednessInXDirectionWithMoleFractions(
                    stat_dump_filename,
                    patch_hierarchy);
        }
        else if (statistical_quantity_key == "MIXEDNESS_X_VOL_F")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleMixednessInXDirectionWithVolumeFractions(
                    stat_dump_filename,
                    patch_hierarchy);
        }
        else if (statistical_quantity_key == "DENSITY_MEAN_INHOMO_X")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleDensityMeanWithInhomogeneousXDirection(
                    stat_dump_filename,
                    patch_hierarchy);
        }
        else if (statistical_quantity_key == "SP_VOL_MEAN_INHOMO_X")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleSpecificVolumeMeanWithInhomogeneousXDirection(
                    stat_dump_filename,
                    patch_hierarchy);
        }
        else if (statistical_quantity_key == "PRESS_MEAN_INHOMO_X")
        {
            rti_rmi_statistics_utilities->
                outputEnsemblePressureMeanWithInhomogeneousXDirection(
                    stat_dump_filename,
                    patch_hierarchy);
        }
        else if (statistical_quantity_key == "TKE_INT")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleTKEIntegrateWithInhomogeneousXDirection(
                    stat_dump_filename,
                    patch_hierarchy);
        }
        else if (statistical_quantity_key == "a1_MEAN_INHOMO_X")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleTurbulentMassFluxVelocityXMeanWithInhomogeneousXDirection(
                    stat_dump_filename,
                    patch_hierarchy);
        }
        else if (statistical_quantity_key == "b_MEAN_INHOMO_X")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleDensitySpecificVolumeCovarianceMeanWithInhomogeneousXDirection(
                    stat_dump_filename,
                    patch_hierarchy);
        }
        else if (statistical_quantity_key == "R11_MEAN_INHOMO_X")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleReynoldsNormalStressXMeanWithInhomogeneousXDirection(
                    stat_dump_filename,
                    patch_hierarchy);
        }
        else if (statistical_quantity_key == "R22_MEAN_INHOMO_X")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleReynoldsNormalStressYMeanWithInhomogeneousXDirection(
                    stat_dump_filename,
                    patch_hierarchy);
        }
        else if (statistical_quantity_key == "R33_MEAN_INHOMO_X")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleReynoldsNormalStressZMeanWithInhomogeneousXDirection(
                    stat_dump_filename,
                    patch_hierarchy);
        }
        else if (statistical_quantity_key == "ENSTROPHY_INT")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleEnstrophyIntegrated(
                    stat_dump_filename,
                    patch_hierarchy);
        }
        else if (statistical_quantity_key == "ENSTROPHY_MEAN_INHO_X")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleEnstrophyMeanWithInhomogeneousXDirection(
                    stat_dump_filename,
                    patch_hierarchy);
        }
        else if (statistical_quantity_key == "S_DIS_RAT_INT")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleScalarDissipationRateIntegrated(
                    stat_dump_filename,
                    patch_hierarchy);
        }
        else if (statistical_quantity_key == "S_DIS_RAT_MEAN_INHO_X")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleScalarDissipationRateMeanWithInhomogeneousXDirection(
                    stat_dump_filename,
                    patch_hierarchy);
        }
        else if (statistical_quantity_key == "RE_W_INHOMO_X")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleMixingWidthReynoldsNumberWithInhomogeneousXDirection(
                    stat_dump_filename,
                    patch_hierarchy);
        }
        else if (statistical_quantity_key == "MA_T_INHOMO_X")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleTurbulentMachNumberWithInhomogeneousXDirection(
                    stat_dump_filename,
                    patch_hierarchy);
        }
        else if (statistical_quantity_key == "EFF_AT_INHOMO_X")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleEffectiveAtwoodNumberWithInhomogeneousXDirection(
                    stat_dump_filename,
                    patch_hierarchy);
        }
        else if (statistical_quantity_key == "BARO_TQ_Z_INT")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleRMSBaroclinicTorqueZIntegrated(
                    stat_dump_filename,
                    patch_hierarchy);
        }
        else if (statistical_quantity_key == "BARO_TQ_Z_MEAN_INHO_X")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleRMSBaroclinicTorqueZMeanWithInhomogeneousXDirection(
                    stat_dump_filename,
                    patch_hierarchy);
        }
        else if (statistical_quantity_key == "BARO_TQ_X_INT")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleRMSBaroclinicTorqueXIntegrated(
                    stat_dump_filename,
                    patch_hierarchy);
        }
        else if (statistical_quantity_key == "BARO_TQ_X_MEAN_INHO_X")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleRMSBaroclinicTorqueXMeanWithInhomogeneousXDirection(
                    stat_dump_filename,
                    patch_hierarchy);
        }
        else if (statistical_quantity_key == "BARO_TQ_Y_INT")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleRMSBaroclinicTorqueYIntegrated(
                    stat_dump_filename,
                    patch_hierarchy);
        }
        else if (statistical_quantity_key == "BARO_TQ_Y_MEAN_INHO_X")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleRMSBaroclinicTorqueYMeanWithInhomogeneousXDirection(
                    stat_dump_filename,
                    patch_hierarchy);
        }
        else if (statistical_quantity_key == "TKE_PROD1_INT")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleTKEProduction1Integrated(
                    stat_dump_filename,
                    patch_hierarchy);
        }
        else if (statistical_quantity_key == "TKE_PROD1_MEAN_INHO_X")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleTKEProduction1MeanWithInhomogeneousXDirection(
                    stat_dump_filename,
                    patch_hierarchy);
        }
        else if (statistical_quantity_key == "TKE_PROD2_INT")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleTKEProduction2Integrated(
                    stat_dump_filename,
                    patch_hierarchy);
        }
        else if (statistical_quantity_key == "TKE_PROD2_MEAN_INHO_X")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleTKEProduction2MeanWithInhomogeneousXDirection(
                    stat_dump_filename,
                    patch_hierarchy);
        }
        else if (statistical_quantity_key == "TKE_PROD3_INT")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleTKEProduction3Integrated(
                    stat_dump_filename,
                    patch_hierarchy);
        }
        else if (statistical_quantity_key == "TKE_PROD3_MEAN_INHO_X")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleTKEProduction3MeanWithInhomogeneousXDirection(
                    stat_dump_filename,
                    patch_hierarchy);
        }
        else if (statistical_quantity_key == "TKE_DISS_INT")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleTKEDissipationIntegrated(
                    stat_dump_filename,
                    patch_hierarchy);
        }
        else if (statistical_quantity_key == "TKE_DISS_MEAN_INHO_X")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleTKEDissipationMeanWithInhomogeneousXDirection(
                    stat_dump_filename,
                    patch_hierarchy);
        }
        else if (statistical_quantity_key == "SHEAR_VIS_MEAN_INHO_X")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleShearViscosityMeanWithInhomogeneousXDirection(
                    stat_dump_filename,
                    patch_hierarchy);
        }
        else if (statistical_quantity_key == "MASS_DIFF_MEAN_INHO_X")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleMassDiffusivityMeanWithInhomogeneousXDirection(
                    stat_dump_filename,
                    patch_hierarchy);
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
