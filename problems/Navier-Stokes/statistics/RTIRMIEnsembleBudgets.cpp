#include "flow/flow_models/four-eqn_conservative/FlowModelStatisticsUtilitiesFourEqnConservative.hpp"

#include "extn/patch_hierarchies/ExtendedFlattenedHierarchy.hpp"
#include "flow/flow_models/MPI_helpers/FlowModelMPIHelperAverage.hpp"
#include "flow/flow_models/MPI_helpers/FlowModelMPIHelperCorrelation.hpp"
#include "flow/flow_models/MPI_helpers/FlowModelMPIHelperMaxMin.hpp"
#include "util/MPI_helpers/MPIHelperGrid.hpp"

#include <fstream>

class EnsembleBudgetsRTIRMI: public EnsembleStatistics
{
    public:
        EnsembleBudgetsRTIRMI(const std::string& object_name):
            EnsembleStatistics(
                object_name)
        {
            setVariablesNotComputed();
        }
        
        void setVariablesNotComputed()
        {
            grid_level_num_avg_computed = false;
            
            rho_avg_coarsest_computed   = false;
            p_avg_coarsest_computed     = false;
            u_avg_coarsest_computed     = false;
            rho_u_avg_coarsest_computed = false;
            p_u_avg_coarsest_computed   = false;
            
            Y_0_avg_computed     = false;
            rho_avg_computed     = false;
            rho_inv_avg_computed = false;
            p_avg_computed       = false;
            u_avg_computed       = false;
            v_avg_computed       = false;
            w_avg_computed       = false;
            u_sq_avg_computed    = false;
            rho_u_avg_computed   = false;
            rho_v_avg_computed   = false;
            rho_w_avg_computed   = false;
            rho_u_u_avg_computed = false;
            rho_v_v_avg_computed = false;
            rho_w_w_avg_computed = false;
            
            rho_u_v_avg_computed = false;
            rho_u_w_avg_computed = false;
            
            p_u_avg_computed = false;
            
            ddx_rho_avg_computed       = false;
            ddx_p_avg_computed         = false;
            ddx_u_avg_computed         = false;
            ddx_v_avg_computed         = false;
            ddx_w_avg_computed         = false;
            ddx_u_sq_avg_computed      = false;
            ddx_rho_u_avg_computed     = false;
            ddx_rho_v_avg_computed     = false;
            ddx_rho_w_avg_computed     = false;
            ddx_rho_u_u_avg_computed   = false;
            ddx_rho_v_v_avg_computed   = false;
            ddx_rho_w_w_avg_computed   = false;
            ddx_rho_u_u_u_avg_computed = false;
            ddx_u_p_avg_computed       = false;
            
            ddx_rho_u_v_avg_computed   = false;
            ddx_rho_u_w_avg_computed   = false;
            ddx_rho_u_v_v_avg_computed = false;
            ddx_rho_u_w_w_avg_computed = false;
            
            ddy_u_avg_computed = false;
            ddy_v_avg_computed = false;
            ddy_w_avg_computed = false;
            
            ddz_u_avg_computed = false;
            ddz_v_avg_computed = false;
            ddz_w_avg_computed = false;
            
            rho_inv_ddx_p_avg_computed = false;
            
            p_ddx_u_avg_computed = false;
            p_ddy_v_avg_computed = false;
            p_ddz_w_avg_computed = false;
            
            u_ddx_u_avg_computed = false;
            u_ddy_v_avg_computed = false;
            u_ddz_w_avg_computed = false;
            
            tau11_avg_computed = false;
            tau12_avg_computed = false;
            tau13_avg_computed = false;
            tau22_avg_computed = false;
            tau23_avg_computed = false;
            tau33_avg_computed = false;
            
            u_tau11_avg_computed = false;
            v_tau12_avg_computed = false;
            w_tau13_avg_computed = false;
            
            ddx_tau11_avg_computed = false;
            ddy_tau12_avg_computed = false;
            ddz_tau13_avg_computed = false;
            
            tau11_ddx_u_avg_computed = false;
            tau12_ddy_u_avg_computed = false;
            tau13_ddz_u_avg_computed = false;
            
            tau12_ddx_v_avg_computed = false;
            tau22_ddy_v_avg_computed = false;
            tau23_ddz_v_avg_computed = false;
            tau13_ddx_w_avg_computed = false;
            tau23_ddy_w_avg_computed = false;
            tau33_ddz_w_avg_computed = false;
            
            rho_inv_ddx_tau11_avg_computed = false;
            rho_inv_ddy_tau12_avg_computed = false;
            rho_inv_ddz_tau13_avg_computed = false;
        }
        
        void clearAllData()
        {
            grid_level_num_avg_realizations.clear();
            
            rho_avg_coarsest_realizations.clear();
            p_avg_coarsest_realizations.clear();
            u_avg_coarsest_realizations.clear();
            rho_u_avg_coarsest_realizations.clear();
            p_u_avg_coarsest_realizations.clear();
            
            Y_0_avg_realizations.clear();
            rho_avg_realizations.clear();
            rho_inv_avg_realizations.clear();
            p_avg_realizations.clear();
            u_avg_realizations.clear();
            v_avg_realizations.clear();
            w_avg_realizations.clear();
            u_sq_avg_realizations.clear();
            rho_u_avg_realizations.clear();
            rho_v_avg_realizations.clear();
            rho_w_avg_realizations.clear();
            rho_u_u_avg_realizations.clear();
            rho_v_v_avg_realizations.clear();
            rho_w_w_avg_realizations.clear();
            
            rho_u_v_avg_realizations.clear();
            rho_u_w_avg_realizations.clear();
            
            p_u_avg_realizations.clear();
            
            ddx_rho_avg_realizations.clear();
            ddx_p_avg_realizations.clear();
            ddx_u_avg_realizations.clear();
            ddx_v_avg_realizations.clear();
            ddx_w_avg_realizations.clear();
            ddx_u_sq_avg_realizations.clear();
            ddx_rho_u_avg_realizations.clear();
            ddx_rho_v_avg_realizations.clear();
            ddx_rho_w_avg_realizations.clear();
            ddx_rho_u_u_avg_realizations.clear();
            ddx_rho_v_v_avg_realizations.clear();
            ddx_rho_w_w_avg_realizations.clear();
            ddx_rho_u_u_u_avg_realizations.clear();
            ddx_u_p_avg_realizations.clear();
            
            ddx_rho_u_v_avg_realizations.clear();
            ddx_rho_u_w_avg_realizations.clear();
            ddx_rho_u_v_v_avg_realizations.clear();
            ddx_rho_u_w_w_avg_realizations.clear();
            
            ddy_u_avg_realizations.clear();
            ddy_v_avg_realizations.clear();
            ddy_w_avg_realizations.clear();
            
            ddz_u_avg_realizations.clear();
            ddz_v_avg_realizations.clear();
            ddz_w_avg_realizations.clear();
            
            rho_inv_ddx_p_avg_realizations.clear();
            
            p_ddx_u_avg_realizations.clear();
            p_ddy_v_avg_realizations.clear();
            p_ddz_w_avg_realizations.clear();
            
            u_ddx_u_avg_realizations.clear();
            u_ddy_v_avg_realizations.clear();
            u_ddz_w_avg_realizations.clear();
            
            tau11_avg_realizations.clear();
            tau12_avg_realizations.clear();
            tau13_avg_realizations.clear();
            tau22_avg_realizations.clear();
            tau23_avg_realizations.clear();
            tau33_avg_realizations.clear();
            
            u_tau11_avg_realizations.clear();
            v_tau12_avg_realizations.clear();
            w_tau13_avg_realizations.clear();
            
            ddx_tau11_avg_realizations.clear();
            ddy_tau12_avg_realizations.clear();
            ddz_tau13_avg_realizations.clear();
            
            tau11_ddx_u_avg_realizations.clear();
            tau12_ddy_u_avg_realizations.clear();
            tau13_ddz_u_avg_realizations.clear();
            
            tau12_ddx_v_avg_realizations.clear();
            tau22_ddy_v_avg_realizations.clear();
            tau23_ddz_v_avg_realizations.clear();
            tau13_ddx_w_avg_realizations.clear();
            tau23_ddy_w_avg_realizations.clear();
            tau33_ddz_w_avg_realizations.clear();
            
            rho_inv_ddx_tau11_avg_realizations.clear();
            rho_inv_ddy_tau12_avg_realizations.clear();
            rho_inv_ddz_tau13_avg_realizations.clear();
            
            setVariablesNotComputed();
        }
        
        // Scratch arrays.
        // Number of realizalizations; number of cells.
        
        std::vector<std::vector<double> > grid_level_num_avg_realizations;
        
        std::vector<std::vector<double> > rho_avg_coarsest_realizations;
        std::vector<std::vector<double> > p_avg_coarsest_realizations;
        std::vector<std::vector<double> > u_avg_coarsest_realizations;
        std::vector<std::vector<double> > rho_u_avg_coarsest_realizations;
        std::vector<std::vector<double> > p_u_avg_coarsest_realizations;
        
        std::vector<std::vector<double> > Y_0_avg_realizations;
        std::vector<std::vector<double> > rho_avg_realizations;
        std::vector<std::vector<double> > rho_inv_avg_realizations;
        std::vector<std::vector<double> > p_avg_realizations;
        std::vector<std::vector<double> > u_avg_realizations;
        std::vector<std::vector<double> > v_avg_realizations;
        std::vector<std::vector<double> > w_avg_realizations;
        std::vector<std::vector<double> > u_sq_avg_realizations;
        std::vector<std::vector<double> > rho_u_avg_realizations;
        std::vector<std::vector<double> > rho_v_avg_realizations;
        std::vector<std::vector<double> > rho_w_avg_realizations;
        std::vector<std::vector<double> > rho_u_u_avg_realizations;
        std::vector<std::vector<double> > rho_v_v_avg_realizations;
        std::vector<std::vector<double> > rho_w_w_avg_realizations;
        
        std::vector<std::vector<double> > rho_u_v_avg_realizations;
        std::vector<std::vector<double> > rho_u_w_avg_realizations;
        
        std::vector<std::vector<double> > p_u_avg_realizations;
        
        std::vector<std::vector<double> > ddx_rho_avg_realizations;
        std::vector<std::vector<double> > ddx_p_avg_realizations;
        std::vector<std::vector<double> > ddx_u_avg_realizations;
        std::vector<std::vector<double> > ddx_v_avg_realizations;
        std::vector<std::vector<double> > ddx_w_avg_realizations;
        std::vector<std::vector<double> > ddx_u_sq_avg_realizations;
        std::vector<std::vector<double> > ddx_rho_u_avg_realizations;
        std::vector<std::vector<double> > ddx_rho_v_avg_realizations;
        std::vector<std::vector<double> > ddx_rho_w_avg_realizations;
        std::vector<std::vector<double> > ddx_rho_u_u_avg_realizations;
        std::vector<std::vector<double> > ddx_rho_v_v_avg_realizations;
        std::vector<std::vector<double> > ddx_rho_w_w_avg_realizations;
        std::vector<std::vector<double> > ddx_rho_u_u_u_avg_realizations;
        std::vector<std::vector<double> > ddx_u_p_avg_realizations;
        
        std::vector<std::vector<double> > ddx_rho_u_v_avg_realizations;
        std::vector<std::vector<double> > ddx_rho_u_w_avg_realizations;
        std::vector<std::vector<double> > ddx_rho_u_v_v_avg_realizations;
        std::vector<std::vector<double> > ddx_rho_u_w_w_avg_realizations;
        
        std::vector<std::vector<double> > ddy_u_avg_realizations;
        std::vector<std::vector<double> > ddy_v_avg_realizations;
        std::vector<std::vector<double> > ddy_w_avg_realizations;
        
        std::vector<std::vector<double> > ddz_u_avg_realizations;
        std::vector<std::vector<double> > ddz_v_avg_realizations;
        std::vector<std::vector<double> > ddz_w_avg_realizations;
        
        std::vector<std::vector<double> > rho_inv_ddx_p_avg_realizations;
        
        std::vector<std::vector<double> > p_ddx_u_avg_realizations;
        std::vector<std::vector<double> > p_ddy_v_avg_realizations;
        std::vector<std::vector<double> > p_ddz_w_avg_realizations;
        
        std::vector<std::vector<double> > u_ddx_u_avg_realizations;
        std::vector<std::vector<double> > u_ddy_v_avg_realizations;
        std::vector<std::vector<double> > u_ddz_w_avg_realizations;
        
        std::vector<std::vector<double> > tau11_avg_realizations;
        std::vector<std::vector<double> > tau12_avg_realizations;
        std::vector<std::vector<double> > tau13_avg_realizations;
        std::vector<std::vector<double> > tau22_avg_realizations;
        std::vector<std::vector<double> > tau23_avg_realizations;
        std::vector<std::vector<double> > tau33_avg_realizations;
        
        std::vector<std::vector<double> > u_tau11_avg_realizations;
        std::vector<std::vector<double> > v_tau12_avg_realizations;
        std::vector<std::vector<double> > w_tau13_avg_realizations;
        
        std::vector<std::vector<double> > ddx_tau11_avg_realizations;
        std::vector<std::vector<double> > ddy_tau12_avg_realizations;
        std::vector<std::vector<double> > ddz_tau13_avg_realizations;
        
        std::vector<std::vector<double> > tau11_ddx_u_avg_realizations;
        std::vector<std::vector<double> > tau12_ddy_u_avg_realizations;
        std::vector<std::vector<double> > tau13_ddz_u_avg_realizations;
        
        std::vector<std::vector<double> > tau12_ddx_v_avg_realizations;
        std::vector<std::vector<double> > tau22_ddy_v_avg_realizations;
        std::vector<std::vector<double> > tau23_ddz_v_avg_realizations;
        std::vector<std::vector<double> > tau13_ddx_w_avg_realizations;
        std::vector<std::vector<double> > tau23_ddy_w_avg_realizations;
        std::vector<std::vector<double> > tau33_ddz_w_avg_realizations;
        
        std::vector<std::vector<double> > rho_inv_ddx_tau11_avg_realizations;
        std::vector<std::vector<double> > rho_inv_ddy_tau12_avg_realizations;
        std::vector<std::vector<double> > rho_inv_ddz_tau13_avg_realizations;
        
        // Whether the scratch arrays are filled.
        
        bool grid_level_num_avg_computed;
        
        bool rho_avg_coarsest_computed;
        bool p_avg_coarsest_computed;
        bool u_avg_coarsest_computed;
        bool rho_u_avg_coarsest_computed;
        bool p_u_avg_coarsest_computed;
        
        bool Y_0_avg_computed;
        bool rho_avg_computed;
        bool rho_inv_avg_computed;
        bool p_avg_computed;
        bool u_avg_computed;
        bool v_avg_computed;
        bool w_avg_computed;
        bool u_sq_avg_computed;
        bool rho_u_avg_computed;
        bool rho_v_avg_computed;
        bool rho_w_avg_computed;
        bool rho_u_u_avg_computed;
        bool rho_v_v_avg_computed;
        bool rho_w_w_avg_computed;
        
        bool rho_u_v_avg_computed;
        bool rho_u_w_avg_computed;
        
        bool p_u_avg_computed;
        
        bool ddx_rho_avg_computed;
        bool ddx_p_avg_computed;
        bool ddx_u_avg_computed;
        bool ddx_v_avg_computed;
        bool ddx_w_avg_computed;
        bool ddx_u_sq_avg_computed;
        bool ddx_rho_u_avg_computed;
        bool ddx_rho_v_avg_computed;
        bool ddx_rho_w_avg_computed;
        bool ddx_rho_u_u_avg_computed;
        bool ddx_rho_v_v_avg_computed;
        bool ddx_rho_w_w_avg_computed;
        bool ddx_rho_u_u_u_avg_computed;
        bool ddx_u_p_avg_computed;
        
        bool ddx_rho_u_v_avg_computed;
        bool ddx_rho_u_w_avg_computed;
        bool ddx_rho_u_v_v_avg_computed;
        bool ddx_rho_u_w_w_avg_computed;
        
        bool ddy_u_avg_computed;
        bool ddy_v_avg_computed;
        bool ddy_w_avg_computed;
        
        bool ddz_u_avg_computed;
        bool ddz_v_avg_computed;
        bool ddz_w_avg_computed;
        
        bool rho_inv_ddx_p_avg_computed;
        
        bool p_ddx_u_avg_computed;
        bool p_ddy_v_avg_computed;
        bool p_ddz_w_avg_computed;
        
        bool u_ddx_u_avg_computed;
        bool u_ddy_v_avg_computed;
        bool u_ddz_w_avg_computed;
        
        bool tau11_avg_computed;
        bool tau12_avg_computed;
        bool tau13_avg_computed;
        bool tau22_avg_computed;
        bool tau23_avg_computed;
        bool tau33_avg_computed;
        
        bool u_tau11_avg_computed;
        bool v_tau12_avg_computed;
        bool w_tau13_avg_computed;
        
        bool ddx_tau11_avg_computed;
        bool ddy_tau12_avg_computed;
        bool ddz_tau13_avg_computed;
        
        bool tau11_ddx_u_avg_computed;
        bool tau12_ddy_u_avg_computed;
        bool tau13_ddz_u_avg_computed;
        
        bool tau12_ddx_v_avg_computed;
        bool tau22_ddy_v_avg_computed;
        bool tau23_ddz_v_avg_computed;
        bool tau13_ddx_w_avg_computed;
        bool tau23_ddy_w_avg_computed;
        bool tau33_ddz_w_avg_computed;
        
        bool rho_inv_ddx_tau11_avg_computed;
        bool rho_inv_ddy_tau12_avg_computed;
        bool rho_inv_ddz_tau13_avg_computed;
        
    private:
        
};


class RTIRMIBudgetsUtilities
{
    public:
        RTIRMIBudgetsUtilities(
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
            const HAMERS_SHARED_PTR<EnsembleBudgetsRTIRMI> ensemble_statistics):
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
         * Compute averaged grid level number with assumed homogeneity in y-direction (2D) or yz-plane (3D).
         */
        void
        computeAveragedGridLevelNumberWithHomogeneityInYDirectionOrInYZPlane(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy);
        
        /*
         * Compute averaged quantities with assumed homogeneity in y-direction (2D) or yz-plane (3D).
         */
        void
        computeAveragedQuantitiesWithHomogeneityInYDirectionOrInYZPlane(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Output spatial profile of ensemble averaged grid level number with assumed homogeneity in y-direction (2D) or
         * yz-plane (3D) to a file.
         */
        void
        outputSpatialProfileEnsembleAveragedGridLevelNumberWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const double output_time) const;
        
        /*
         * Output spatial profile of ensemble averaged mass fraction with assumed homogeneity in y-direction (2D) or
         * yz-plane (3D) to a file.
         */
        void
        outputSpatialProfileEnsembleAveragedMassFractionWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const double output_time) const;
        
        /*
         * Output spatial profile of ensemble averaged density with assumed homogeneity in y-direction (2D) or
         * yz-plane (3D) to a file.
         */
        void
        outputSpatialProfileEnsembleAveragedDensityWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const double output_time) const;
        
        /*
         * Output spatial profile of ensemble averaged velocity x-component with assumed homogeneity in y-direction (2D)
         * or yz-plane (3D) to a file.
         */
        void
        outputSpatialProfileEnsembleAveragedVelocityXWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const double output_time) const;
        
        /*
         * Output spatial profile of ensemble Favre averaged velocity x-component with assumed homogeneity in
         * y-direction (2D) or yz-plane (3D) to a file.
         */
        void
        outputSpatialProfileEnsembleFavreAveragedVelocityXWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const double output_time) const;
        
        /*
         * Output spatial profile of ensemble turbulent mass flux velocity in x-direction with assumed homogeneity in
         * y-direction (2D) or yz-plane (3D) to a file.
         */
        void
        outputSpatialProfileEnsembleTurbulentMassFluxVelocityXWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const double output_time) const;
        
        /*
         * Output spatial profile of ensemble averaged pressure gradient in x-direction with assumed homogeneity in
         * y-direction (2D) or yz-plane (3D) to a file.
         */
        void
        outputSpatialProfileEnsembleAveragedPressureGradientXWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const double output_time) const;
        
        /*
         * Output budget of Favre mean TKE with inhomogeneous x-direction to a file.
         */
        void
        outputBudgetFavreMeanTKEWithInhomogeneousXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const double output_time) const;
        
        /*
         * Output budget of turbulent mass flux in x-direction with inhomogeneous x-direction to a file.
         */
        void
        outputBudgetTurbMassFluxXWithInhomogeneousXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const double output_time) const;
        
        /*
         * Output budget of Reynolds normal stress in x-direction with inhomogeneous x-direction to a file.
         */
        void
        outputBudgetReynoldsNormalStressInXDirectionWithInhomogeneousXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const double output_time) const;
        
        /*
         * Output budget of Reynolds normal stress in y-direction with inhomogeneous x-direction to a file.
         */
        void
        outputBudgetReynoldsNormalStressInYDirectionWithInhomogeneousXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const double output_time) const;
        
        /*
         * Output budget of Reynolds normal stress in z-direction with inhomogeneous x-direction to a file.
         */
        void
        outputBudgetReynoldsNormalStressInZDirectionWithInhomogeneousXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const double output_time) const;
        
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
         * Compute averaged value (on product of variable derivatives and derivative of shear stress component) with only
         * x direction as inhomogeneous direction.
         * component_idx:
         * 0: tau11
         * 1: tau12
         * 2: tau13
         * 3: tau22
         * 4: tau23
         * 5: tau33
         */
        std::vector<double>
        getAveragedQuantityWithDerivativeOfShearStressComponentWithInhomogeneousXDirection(
            const std::vector<std::string>& quantity_names,
            const std::vector<int>& component_indices,
            const std::vector<bool>& use_derivative,
            const std::vector<int>& derivative_directions,
            const std::vector<bool>& use_reciprocal,
            const int shear_stress_component_idx,
            const int shear_stress_derivative_direction,
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
        HAMERS_SHARED_PTR<EnsembleBudgetsRTIRMI> d_ensemble_statistics;
        
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
 * Compute averaged grid level number with assumed homogeneity in y-direction (2D) or yz-plane (3D).
 */
void
RTIRMIBudgetsUtilities::computeAveragedGridLevelNumberWithHomogeneityInYDirectionOrInYZPlane(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy)
{
    MPIHelperGrid MPI_helper_grid = MPIHelperGrid(
        "MPI_helper_grid",
        d_dim,
        d_grid_geometry,
        patch_hierarchy);
    
    const double num_cells_global = MPI_helper_grid.getNumberOfCells();
    
    std::vector<double> grid_level_num_avg = MPI_helper_grid.getAveragedGridLevelNumberWithInhomogeneousXDirection();
    
    std::vector<std::vector<double> >& grid_level_num_avg_realizations = d_ensemble_statistics->grid_level_num_avg_realizations;
    grid_level_num_avg_realizations.push_back(grid_level_num_avg);
    
    d_ensemble_statistics->grid_level_num_avg_computed = true;
}


/*
 * Compute averaged quantities with assumed homogeneity in y-direction (2D) or yz-plane (3D).
 */
void
RTIRMIBudgetsUtilities::computeAveragedQuantitiesWithHomogeneityInYDirectionOrInYZPlane(
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
    
    // Scratch data containers to pass to MPI helper.
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    std::vector<int> derivative_directions;
    std::vector<bool> use_reciprocal;
    std::vector<bool> use_derivative;
    
    // Compute rho_avg_coarsest.
    
    std::vector<double> rho_avg_coarsest_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirectionOnCoarsestLevel(
        "DENSITY",
        0,
        data_context);
    
    std::vector<std::vector<double> >& rho_avg_coarsest_realizations = d_ensemble_statistics->rho_avg_coarsest_realizations;
    rho_avg_coarsest_realizations.push_back(rho_avg_coarsest_global);
    
    d_ensemble_statistics->rho_avg_coarsest_computed = true;
    
    // Compute p_avg_coarsest.
    
    std::vector<double> p_avg_coarsest_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirectionOnCoarsestLevel(
        "PRESSURE",
        0,
        data_context);
    
    std::vector<std::vector<double> >& p_avg_coarsest_realizations = d_ensemble_statistics->p_avg_coarsest_realizations;
    p_avg_coarsest_realizations.push_back(p_avg_coarsest_global);
    
    d_ensemble_statistics->p_avg_coarsest_computed = true;
    
    // Compute u_avg_coarsest.
    
    std::vector<double> u_avg_coarsest_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirectionOnCoarsestLevel(
        "VELOCITY",
        0,
        data_context);
    
    std::vector<std::vector<double> >& u_avg_coarsest_realizations = d_ensemble_statistics->u_avg_coarsest_realizations;
    u_avg_coarsest_realizations.push_back(u_avg_coarsest_global);
    
    d_ensemble_statistics->u_avg_coarsest_computed = true;
    
    // Compute rho_u_avg_coarsest.
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    
    std::vector<double> rho_u_avg_coarsest_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirectionOnCoarsestLevel(
        quantity_names,
        component_indices,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    
    std::vector<std::vector<double> >& rho_u_avg_coarsest_realizations = d_ensemble_statistics->rho_u_avg_coarsest_realizations;
    rho_u_avg_coarsest_realizations.push_back(rho_u_avg_coarsest_global);
    
    d_ensemble_statistics->rho_u_avg_coarsest_computed = true;
    
    // Compute p_u_avg_coarsest.
    
    quantity_names.push_back("PRESSURE");
    component_indices.push_back(0);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    
    std::vector<double> p_u_avg_coarsest_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirectionOnCoarsestLevel(
        quantity_names,
        component_indices,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    
    std::vector<std::vector<double> >& p_u_avg_coarsest_realizations = d_ensemble_statistics->p_u_avg_coarsest_realizations;
    p_u_avg_coarsest_realizations.push_back(p_u_avg_coarsest_global);
    
    d_ensemble_statistics->p_u_avg_coarsest_computed = true;
    
    // Compute Y_0_avg.
    
    std::vector<double> Y_0_avg = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        "MASS_FRACTIONS",
        0,
        data_context);
    
    std::vector<std::vector<double> >& Y_0_avg_realizations = d_ensemble_statistics->Y_0_avg_realizations;
    Y_0_avg_realizations.push_back(Y_0_avg);
    
    d_ensemble_statistics->Y_0_avg_computed = true;
    
    // Compute rho_avg.
    
    std::vector<double> rho_avg = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        data_context);
    
    std::vector<std::vector<double> >& rho_avg_realizations = d_ensemble_statistics->rho_avg_realizations;
    rho_avg_realizations.push_back(rho_avg);
    
    d_ensemble_statistics->rho_avg_computed = true;
    
    // Compute rho_inv_avg.
    
    std::vector<double> rho_inv_avg = MPI_helper_average.getAveragedReciprocalOfQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        data_context);
    
    std::vector<std::vector<double> >& rho_inv_avg_realizations = d_ensemble_statistics->rho_inv_avg_realizations;
    rho_inv_avg_realizations.push_back(rho_inv_avg);
    
    d_ensemble_statistics->rho_inv_avg_computed = true;
    
    // Compute p_avg.
    
    std::vector<double> p_avg = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        "PRESSURE",
        0,
        data_context);
    
    std::vector<std::vector<double> >& p_avg_realizations = d_ensemble_statistics->p_avg_realizations;
    p_avg_realizations.push_back(p_avg);
    
    d_ensemble_statistics->p_avg_computed = true;
    
    // Compute u_avg.
    
    std::vector<double> u_avg = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        "VELOCITY",
        0,
        data_context);
    
    std::vector<std::vector<double> >& u_avg_realizations = d_ensemble_statistics->u_avg_realizations;
    u_avg_realizations.push_back(u_avg);
    
    d_ensemble_statistics->u_avg_computed = true;
    
    // Compute v_avg.
    
    if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
    {
        std::vector<double> v_avg = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
            "VELOCITY",
            1,
            data_context);
        
        std::vector<std::vector<double> >& v_avg_realizations = d_ensemble_statistics->v_avg_realizations;
        v_avg_realizations.push_back(v_avg);
        
        d_ensemble_statistics->v_avg_computed = true;
    }
    
    // Compute w_avg.
    
    if (d_dim == tbox::Dimension(3))
    {
        std::vector<double> w_avg = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
            "VELOCITY",
            2,
            data_context);
        
        std::vector<std::vector<double> >& w_avg_realizations = d_ensemble_statistics->w_avg_realizations;
        w_avg_realizations.push_back(w_avg);
        
        d_ensemble_statistics->w_avg_computed = true;
    }
    
    // Compute u_sq_avg.
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    
    std::vector<double> u_sq_avg = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    
    std::vector<std::vector<double> >& u_sq_avg_realizations = d_ensemble_statistics->u_sq_avg_realizations;
    u_sq_avg_realizations.push_back(u_sq_avg);
    
    d_ensemble_statistics->u_sq_avg_computed = true;
    
    // Compute rho_u_avg.
    
    std::vector<double> rho_u_avg = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        "MOMENTUM",
        0,
        data_context);
    
    std::vector<std::vector<double> >& rho_u_avg_realizations = d_ensemble_statistics->rho_u_avg_realizations;
    rho_u_avg_realizations.push_back(rho_u_avg);
    
    d_ensemble_statistics->rho_u_avg_computed = true;
    
    // Compute rho_v_avg.
    
    if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
    {
        std::vector<double> rho_v_avg = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
            "MOMENTUM",
            1,
            data_context);
        
        std::vector<std::vector<double> >& rho_v_avg_realizations = d_ensemble_statistics->rho_v_avg_realizations;
        rho_v_avg_realizations.push_back(rho_v_avg);
        
        d_ensemble_statistics->rho_v_avg_computed = true;
    }
    
    // Compute rho_w_avg.
    
    if (d_dim == tbox::Dimension(3))
    {
        std::vector<double> rho_w_avg = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
            "MOMENTUM",
            2,
            data_context);
        
        std::vector<std::vector<double> >& rho_w_avg_realizations = d_ensemble_statistics->rho_w_avg_realizations;
        rho_w_avg_realizations.push_back(rho_w_avg);
        
        d_ensemble_statistics->rho_w_avg_computed = true;
    }
    
    // Compute rho_u_u_avg.
    
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
    
    // Compute rho_v_v_avg.
    
    if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
    {
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
    
    // Compute rho_w_w_avg.
    
    if (d_dim == tbox::Dimension(3))
    {
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
    
    // Compute rho_u_v_avg.
    
    if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
    {
        quantity_names.push_back("DENSITY");
        component_indices.push_back(0);
        
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(0);
        
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(1);
        
        std::vector<double> rho_u_v_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
            quantity_names,
            component_indices,
            data_context);
        
        quantity_names.clear();
        component_indices.clear();
        
        std::vector<std::vector<double> >& rho_u_v_avg_realizations = d_ensemble_statistics->rho_u_v_avg_realizations;
        rho_u_v_avg_realizations.push_back(rho_u_v_avg_global);
        
        d_ensemble_statistics->rho_u_v_avg_computed = true;
    }
    
    // Compute rho_u_w_avg.
    
    if (d_dim == tbox::Dimension(3))
    {
        quantity_names.push_back("DENSITY");
        component_indices.push_back(0);
        
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(0);
        
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(2);
        
        std::vector<double> rho_u_w_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
            quantity_names,
            component_indices,
            data_context);
        
        quantity_names.clear();
        component_indices.clear();
        
        std::vector<std::vector<double> >& rho_u_w_avg_realizations = d_ensemble_statistics->rho_u_w_avg_realizations;
        rho_u_w_avg_realizations.push_back(rho_u_w_avg_global);
        
        d_ensemble_statistics->rho_u_w_avg_computed = true;
    }
    
    // Compute p_u_avg.
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    
    quantity_names.push_back("PRESSURE");
    component_indices.push_back(0);
    
    std::vector<double> p_u_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    
    std::vector<std::vector<double> >& p_u_avg_realizations = d_ensemble_statistics->p_u_avg_realizations;
    p_u_avg_realizations.push_back(p_u_avg_global);
    
    d_ensemble_statistics->p_u_avg_computed = true;
    
    // Compute ddx_rho_avg.
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_reciprocal.push_back(false);
    
    std::vector<double> ddx_rho_avg = MPI_helper_average.getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        use_reciprocal,
        0,
        d_num_ghosts_derivative,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    use_reciprocal.clear();
    
    std::vector<std::vector<double> >& ddx_rho_avg_realizations = d_ensemble_statistics->ddx_rho_avg_realizations;
    ddx_rho_avg_realizations.push_back(ddx_rho_avg);
    
    d_ensemble_statistics->ddx_rho_avg_computed = true;
    
    // Compute ddx_p_avg.
    
    quantity_names.push_back("PRESSURE");
    component_indices.push_back(0);
    use_reciprocal.push_back(false);
    
    std::vector<double> ddx_p_avg = MPI_helper_average.getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        use_reciprocal,
        0,
        d_num_ghosts_derivative,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    use_reciprocal.clear();
    
    std::vector<std::vector<double> >& ddx_p_avg_realizations = d_ensemble_statistics->ddx_p_avg_realizations;
    ddx_p_avg_realizations.push_back(ddx_p_avg);
    
    d_ensemble_statistics->ddx_p_avg_computed = true;
    
    // Computed ddx_u_avg.
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    use_reciprocal.push_back(false);
    
    std::vector<double> ddx_u_avg = MPI_helper_average.getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        use_reciprocal,
        0,
        d_num_ghosts_derivative,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    use_reciprocal.clear();
    
    std::vector<std::vector<double> >& ddx_u_avg_realizations = d_ensemble_statistics->ddx_u_avg_realizations;
    ddx_u_avg_realizations.push_back(ddx_u_avg);
    
    d_ensemble_statistics->ddx_u_avg_computed = true;
    
    // Computed ddx_v_avg.
    
    if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
    {
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(1);
        use_reciprocal.push_back(false);
        
        std::vector<double> ddx_v_avg = MPI_helper_average.getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
            quantity_names,
            component_indices,
            use_reciprocal,
            0,
            d_num_ghosts_derivative,
            data_context);
        
        quantity_names.clear();
        component_indices.clear();
        use_reciprocal.clear();
        
        std::vector<std::vector<double> >& ddx_v_avg_realizations = d_ensemble_statistics->ddx_v_avg_realizations;
        ddx_v_avg_realizations.push_back(ddx_v_avg);
        
        d_ensemble_statistics->ddx_v_avg_computed = true;
    }
    
    // Computed ddx_w_avg.
    
    if (d_dim == tbox::Dimension(3))
    {
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(2);
        use_reciprocal.push_back(false);
        
        std::vector<double> ddx_w_avg = MPI_helper_average.getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
            quantity_names,
            component_indices,
            use_reciprocal,
            0,
            d_num_ghosts_derivative,
            data_context);
        
        quantity_names.clear();
        component_indices.clear();
        use_reciprocal.clear();
        
        std::vector<std::vector<double> >& ddx_w_avg_realizations = d_ensemble_statistics->ddx_w_avg_realizations;
        ddx_w_avg_realizations.push_back(ddx_w_avg);
        
        d_ensemble_statistics->ddx_w_avg_computed = true;
    }
    
    // Computed dx_u_sq_avg.
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    use_reciprocal.push_back(false);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    use_reciprocal.push_back(false);
    
    std::vector<double> ddx_u_sq_avg = MPI_helper_average.getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        use_reciprocal,
        0,
        d_num_ghosts_derivative,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    use_reciprocal.clear();
    
    std::vector<std::vector<double> >& ddx_u_sq_avg_realizations = d_ensemble_statistics->ddx_u_sq_avg_realizations;
    ddx_u_sq_avg_realizations.push_back(ddx_u_sq_avg);
    
    d_ensemble_statistics->ddx_u_sq_avg_computed = true;
    
    // Compute ddx_rho_u_avg.
    
    quantity_names.push_back("MOMENTUM");
    component_indices.push_back(0);
    use_reciprocal.push_back(false);
    
    std::vector<double> ddx_rho_u_avg = MPI_helper_average.getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        use_reciprocal,
        0,
        d_num_ghosts_derivative,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    use_reciprocal.clear();
    
    std::vector<std::vector<double> >& ddx_rho_u_avg_realizations = d_ensemble_statistics->ddx_rho_u_avg_realizations;
    ddx_rho_u_avg_realizations.push_back(ddx_rho_u_avg);
    
    d_ensemble_statistics->ddx_rho_u_avg_computed = true;
    
    // Compute ddx_rho_v_avg.
    
    if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
    {
        quantity_names.push_back("MOMENTUM");
        component_indices.push_back(1);
        use_reciprocal.push_back(false);
        
        std::vector<double> ddx_rho_v_avg = MPI_helper_average.getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
            quantity_names,
            component_indices,
            use_reciprocal,
            0,
            d_num_ghosts_derivative,
            data_context);
        
        quantity_names.clear();
        component_indices.clear();
        use_reciprocal.clear();
        
        std::vector<std::vector<double> >& ddx_rho_v_avg_realizations = d_ensemble_statistics->ddx_rho_v_avg_realizations;
        ddx_rho_v_avg_realizations.push_back(ddx_rho_v_avg);
        
        d_ensemble_statistics->ddx_rho_v_avg_computed = true;
    }
    
    // Compute ddx_rho_w_avg.
    
    if (d_dim == tbox::Dimension(3))
    {
        quantity_names.push_back("MOMENTUM");
        component_indices.push_back(2);
        use_reciprocal.push_back(false);
        
        std::vector<double> ddx_rho_w_avg = MPI_helper_average.getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
            quantity_names,
            component_indices,
            use_reciprocal,
            0,
            d_num_ghosts_derivative,
            data_context);
        
        quantity_names.clear();
        component_indices.clear();
        use_reciprocal.clear();
        
        std::vector<std::vector<double> >& ddx_rho_w_avg_realizations = d_ensemble_statistics->ddx_rho_w_avg_realizations;
        ddx_rho_w_avg_realizations.push_back(ddx_rho_w_avg);
        
        d_ensemble_statistics->ddx_rho_w_avg_computed = true;
    }
    
    // Compute ddx_rho_u_u_avg.
    
    quantity_names.push_back("MOMENTUM");
    component_indices.push_back(0);
    use_reciprocal.push_back(false);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    use_reciprocal.push_back(false);
    
    std::vector<double> ddx_rho_u_u_avg = MPI_helper_average.getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        use_reciprocal,
        0,
        d_num_ghosts_derivative,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    use_reciprocal.clear();
    
    std::vector<std::vector<double> >& ddx_rho_u_u_avg_realizations = d_ensemble_statistics->ddx_rho_u_u_avg_realizations;
    ddx_rho_u_u_avg_realizations.push_back(ddx_rho_u_u_avg);
    
    d_ensemble_statistics->ddx_rho_u_u_avg_computed = true;
    
    // Compute ddx_rho_v_v_avg.
    
    if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
    {
        quantity_names.push_back("MOMENTUM");
        component_indices.push_back(1);
        use_reciprocal.push_back(false);
        
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(1);
        use_reciprocal.push_back(false);
        
        std::vector<double> ddx_rho_v_v_avg = MPI_helper_average.getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
            quantity_names,
            component_indices,
            use_reciprocal,
            0,
            d_num_ghosts_derivative,
            data_context);
        
        quantity_names.clear();
        component_indices.clear();
        use_reciprocal.clear();
        
        std::vector<std::vector<double> >& ddx_rho_v_v_avg_realizations = d_ensemble_statistics->ddx_rho_v_v_avg_realizations;
        ddx_rho_v_v_avg_realizations.push_back(ddx_rho_v_v_avg);
        
        d_ensemble_statistics->ddx_rho_v_v_avg_computed = true;
    }
    
    // Compute ddx_rho_w_w_avg.
    
    if (d_dim == tbox::Dimension(3))
    {
        quantity_names.push_back("MOMENTUM");
        component_indices.push_back(2);
        use_reciprocal.push_back(false);
        
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(2);
        use_reciprocal.push_back(false);
        
        std::vector<double> ddx_rho_w_w_avg = MPI_helper_average.getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
            quantity_names,
            component_indices,
            use_reciprocal,
            0,
            d_num_ghosts_derivative,
            data_context);
        
        quantity_names.clear();
        component_indices.clear();
        use_reciprocal.clear();
        
        std::vector<std::vector<double> >& ddx_rho_w_w_avg_realizations = d_ensemble_statistics->ddx_rho_w_w_avg_realizations;
        ddx_rho_w_w_avg_realizations.push_back(ddx_rho_w_w_avg);
        
        d_ensemble_statistics->ddx_rho_w_w_avg_computed = true;
    }
    
    // Compute ddx_rho_u_u_u_avg.
    
    quantity_names.push_back("MOMENTUM");
    component_indices.push_back(0);
    use_reciprocal.push_back(false);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    use_reciprocal.push_back(false);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    use_reciprocal.push_back(false);
    
    std::vector<double> ddx_rho_u_u_u_avg = MPI_helper_average.getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        use_reciprocal,
        0,
        d_num_ghosts_derivative,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    use_reciprocal.clear();
    
    std::vector<std::vector<double> >& ddx_rho_u_u_u_avg_realizations = d_ensemble_statistics->ddx_rho_u_u_u_avg_realizations;
    ddx_rho_u_u_u_avg_realizations.push_back(ddx_rho_u_u_u_avg);
    
    // Compute ddx_u_p_avg.
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    use_reciprocal.push_back(false);
    
    quantity_names.push_back("PRESSURE");
    component_indices.push_back(0);
    use_reciprocal.push_back(false);
    
    std::vector<double> ddx_u_p_avg = MPI_helper_average.getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        use_reciprocal,
        0,
        d_num_ghosts_derivative,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    use_reciprocal.clear();
    
    std::vector<std::vector<double> >& ddx_u_p_avg_realizations = d_ensemble_statistics->ddx_u_p_avg_realizations;
    ddx_u_p_avg_realizations.push_back(ddx_u_p_avg);
    
    d_ensemble_statistics->ddx_u_p_avg_computed = true;
    
    // Compute ddx_rho_u_v_avg.
    
    if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
    {
        quantity_names.push_back("MOMENTUM");
        component_indices.push_back(0);
        use_reciprocal.push_back(false);
        
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(1);
        use_reciprocal.push_back(false);
        
        std::vector<double> ddx_rho_u_v_avg = MPI_helper_average.getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
            quantity_names,
            component_indices,
            use_reciprocal,
            0,
            d_num_ghosts_derivative,
            data_context);
        
        quantity_names.clear();
        component_indices.clear();
        use_reciprocal.clear();
        
        std::vector<std::vector<double> >& ddx_rho_u_v_avg_realizations = d_ensemble_statistics->ddx_rho_u_v_avg_realizations;
        ddx_rho_u_v_avg_realizations.push_back(ddx_rho_u_v_avg);
        
        d_ensemble_statistics->ddx_rho_u_v_avg_computed = true;
    }
    
    // Compute ddx_rho_u_w_avg.
    
    if (d_dim == tbox::Dimension(3))
    {
        quantity_names.push_back("MOMENTUM");
        component_indices.push_back(0);
        use_reciprocal.push_back(false);
        
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(2);
        use_reciprocal.push_back(false);
        
        std::vector<double> ddx_rho_u_w_avg = MPI_helper_average.getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
            quantity_names,
            component_indices,
            use_reciprocal,
            0,
            d_num_ghosts_derivative,
            data_context);
        
        quantity_names.clear();
        component_indices.clear();
        use_reciprocal.clear();
        
        std::vector<std::vector<double> >& ddx_rho_u_w_avg_realizations = d_ensemble_statistics->ddx_rho_u_w_avg_realizations;
        ddx_rho_u_w_avg_realizations.push_back(ddx_rho_u_w_avg);
        
        d_ensemble_statistics->ddx_rho_u_w_avg_computed = true;
    }
    
    // Compute ddx_rho_u_v_v_avg.
    
    if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
    {
        quantity_names.push_back("MOMENTUM");
        component_indices.push_back(0);
        use_reciprocal.push_back(false);
        
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(1);
        use_reciprocal.push_back(false);
        
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(1);
        use_reciprocal.push_back(false);
        
        std::vector<double> ddx_rho_u_v_v_avg = MPI_helper_average.getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
            quantity_names,
            component_indices,
            use_reciprocal,
            0,
            d_num_ghosts_derivative,
            data_context);
        
        quantity_names.clear();
        component_indices.clear();
        use_reciprocal.clear();
        
        std::vector<std::vector<double> >& ddx_rho_u_v_v_avg_realizations = d_ensemble_statistics->ddx_rho_u_v_v_avg_realizations;
        ddx_rho_u_v_v_avg_realizations.push_back(ddx_rho_u_v_v_avg);
        
        d_ensemble_statistics->ddx_rho_u_v_v_avg_computed = true;
    }
    
    // Compute ddx_rho_u_w_w_avg.
    
    if (d_dim == tbox::Dimension(3))
    {
        quantity_names.push_back("MOMENTUM");
        component_indices.push_back(0);
        use_reciprocal.push_back(false);
        
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(2);
        use_reciprocal.push_back(false);
        
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(2);
        use_reciprocal.push_back(false);
        
        std::vector<double> ddx_rho_u_w_w_avg = MPI_helper_average.getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
            quantity_names,
            component_indices,
            use_reciprocal,
            0,
            d_num_ghosts_derivative,
            data_context);
        
        quantity_names.clear();
        component_indices.clear();
        use_reciprocal.clear();
        
        std::vector<std::vector<double> >& ddx_rho_u_w_w_avg_realizations = d_ensemble_statistics->ddx_rho_u_w_w_avg_realizations;
        ddx_rho_u_w_w_avg_realizations.push_back(ddx_rho_u_w_w_avg);
        
        d_ensemble_statistics->ddx_rho_u_w_w_avg_computed = true;
    }
    
    // Compute ddy_u_avg.
    
    if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
    {
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(0);
        use_reciprocal.push_back(false);
        
        std::vector<double> ddy_u_avg = MPI_helper_average.getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
            quantity_names,
            component_indices,
            use_reciprocal,
            1,
            d_num_ghosts_derivative,
            data_context);
        
        quantity_names.clear();
        component_indices.clear();
        use_reciprocal.clear();
        
        std::vector<std::vector<double> >& ddy_u_avg_realizations = d_ensemble_statistics->ddy_u_avg_realizations;
        ddy_u_avg_realizations.push_back(ddy_u_avg);
        
        d_ensemble_statistics->ddy_u_avg_computed = true;
    }
    
    // Compute ddy_v_avg.
    
    if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
    {
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(1);
        use_reciprocal.push_back(false);
    
        std::vector<double> ddy_v_avg = MPI_helper_average.getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
            quantity_names,
            component_indices,
            use_reciprocal,
            1,
            d_num_ghosts_derivative,
            data_context);
        
        quantity_names.clear();
        component_indices.clear();
        use_reciprocal.clear();
        
        std::vector<std::vector<double> >& ddy_v_avg_realizations = d_ensemble_statistics->ddy_v_avg_realizations;
        ddy_v_avg_realizations.push_back(ddy_v_avg);
        
        d_ensemble_statistics->ddy_v_avg_computed = true;
    }
    
    // Compute ddy_w_avg.
    
    if (d_dim == tbox::Dimension(3))
    {
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(2);
        use_reciprocal.push_back(false);
    
        std::vector<double> ddy_w_avg = MPI_helper_average.getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
            quantity_names,
            component_indices,
            use_reciprocal,
            1,
            d_num_ghosts_derivative,
            data_context);
        
        quantity_names.clear();
        component_indices.clear();
        use_reciprocal.clear();
        
        std::vector<std::vector<double> >& ddy_w_avg_realizations = d_ensemble_statistics->ddy_w_avg_realizations;
        ddy_w_avg_realizations.push_back(ddy_w_avg);
        
        d_ensemble_statistics->ddy_w_avg_computed = true;
    }
    
    // Compute ddz_u_avg.
    
    if (d_dim == tbox::Dimension(3))
    {
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(0);
        use_reciprocal.push_back(false);
        
        std::vector<double> ddz_u_avg = MPI_helper_average.getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
            quantity_names,
            component_indices,
            use_reciprocal,
            2,
            d_num_ghosts_derivative,
            data_context);
        
        quantity_names.clear();
        component_indices.clear();
        use_reciprocal.clear();
        
        std::vector<std::vector<double> >& ddz_u_avg_realizations = d_ensemble_statistics->ddz_u_avg_realizations;
        ddz_u_avg_realizations.push_back(ddz_u_avg);
        
        d_ensemble_statistics->ddz_u_avg_computed = true;
    }
    
    // Compute ddz_v_avg.
    
    if (d_dim == tbox::Dimension(3))
    {
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(1);
        use_reciprocal.push_back(false);
        
        std::vector<double> ddz_v_avg = MPI_helper_average.getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
            quantity_names,
            component_indices,
            use_reciprocal,
            2,
            d_num_ghosts_derivative,
            data_context);
        
        quantity_names.clear();
        component_indices.clear();
        use_reciprocal.clear();
        
        std::vector<std::vector<double> >& ddz_v_avg_realizations = d_ensemble_statistics->ddz_v_avg_realizations;
        ddz_v_avg_realizations.push_back(ddz_v_avg);
        
        d_ensemble_statistics->ddz_v_avg_computed = true;
    }
    
    // Compute ddz_w_avg.
    
    if (d_dim == tbox::Dimension(3))
    {
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(2);
        use_reciprocal.push_back(false);
    
        std::vector<double> ddz_w_avg = MPI_helper_average.getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
            quantity_names,
            component_indices,
            use_reciprocal,
            2,
            d_num_ghosts_derivative,
            data_context);
        
        quantity_names.clear();
        component_indices.clear();
        use_reciprocal.clear();
        
        std::vector<std::vector<double> >& ddz_w_avg_realizations = d_ensemble_statistics->ddz_w_avg_realizations;
        ddz_w_avg_realizations.push_back(ddz_w_avg);
        
        d_ensemble_statistics->ddz_w_avg_computed = true;
    }
    
    // Compute rho_inv_ddx_p_avg.
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(false);
    derivative_directions.push_back(-1);
    use_reciprocal.push_back(true);
    
    quantity_names.push_back("PRESSURE");
    component_indices.push_back(0);
    use_derivative.push_back(true);
    derivative_directions.push_back(0);
    use_reciprocal.push_back(false);
    
    std::vector<double> rho_inv_ddx_p_avg = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
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
    
    std::vector<std::vector<double> >& rho_inv_ddx_p_avg_realizations = d_ensemble_statistics->rho_inv_ddx_p_avg_realizations;
    rho_inv_ddx_p_avg_realizations.push_back(rho_inv_ddx_p_avg);
    
    d_ensemble_statistics->rho_inv_ddx_p_avg_computed = true;
    
    // Compute p_ddx_u_avg.
    
    quantity_names.push_back("PRESSURE");
    component_indices.push_back(0);
    use_derivative.push_back(false);
    derivative_directions.push_back(-1);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    use_derivative.push_back(true);
    derivative_directions.push_back(0);
    
    std::vector<double> p_ddx_u_avg = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
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
    
    std::vector<std::vector<double> >& p_ddx_u_avg_realizations = d_ensemble_statistics->p_ddx_u_avg_realizations;
    p_ddx_u_avg_realizations.push_back(p_ddx_u_avg);
    
    d_ensemble_statistics->p_ddx_u_avg_computed = true;
    
    // Compute p_ddy_v_avg.
    
    if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
    {
        quantity_names.push_back("PRESSURE");
        component_indices.push_back(0);
        use_derivative.push_back(false);
        derivative_directions.push_back(-1);
        
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(1);
        use_derivative.push_back(true);
        derivative_directions.push_back(1);
        
        std::vector<double> p_ddy_v_avg = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
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
        
        std::vector<std::vector<double> >& p_ddy_v_avg_realizations = d_ensemble_statistics->p_ddy_v_avg_realizations;
        p_ddy_v_avg_realizations.push_back(p_ddy_v_avg);
        
        d_ensemble_statistics->p_ddy_v_avg_computed = true;
    }
    
    // Compute p_ddz_w_avg.
    
    if (d_dim == tbox::Dimension(3))
    {
        quantity_names.push_back("PRESSURE");
        component_indices.push_back(0);
        use_derivative.push_back(false);
        derivative_directions.push_back(-1);
        
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(2);
        use_derivative.push_back(true);
        derivative_directions.push_back(2);
        
        std::vector<double> p_ddz_w_avg = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
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
        
        std::vector<std::vector<double> >& p_ddz_w_avg_realizations = d_ensemble_statistics->p_ddz_w_avg_realizations;
        p_ddz_w_avg_realizations.push_back(p_ddz_w_avg);
        
        d_ensemble_statistics->p_ddz_w_avg_computed = true;
    }
    
    // Compute u_ddx_u_avg.
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    use_derivative.push_back(false);
    derivative_directions.push_back(-1);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    use_derivative.push_back(true);
    derivative_directions.push_back(0);
    
    std::vector<double> u_ddx_u_avg = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
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
    
    std::vector<std::vector<double> >& u_ddx_u_avg_realizations = d_ensemble_statistics->u_ddx_u_avg_realizations;
    u_ddx_u_avg_realizations.push_back(u_ddx_u_avg);
    
    d_ensemble_statistics->u_ddx_u_avg_computed = true;
    
    // Compute u_ddy_v_avg.
    
    if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
    {
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(0);
        use_derivative.push_back(false);
        derivative_directions.push_back(-1);
        
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(1);
        use_derivative.push_back(true);
        derivative_directions.push_back(1);
        
        std::vector<double> u_ddy_v_avg = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
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
        
        std::vector<std::vector<double> >& u_ddy_v_avg_realizations = d_ensemble_statistics->u_ddy_v_avg_realizations;
        u_ddy_v_avg_realizations.push_back(u_ddy_v_avg);
        
        d_ensemble_statistics->u_ddy_v_avg_computed = true;
    }
    
    // Compute u_ddz_w_avg.
    
    if (d_dim == tbox::Dimension(3))
    {
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(0);
        use_derivative.push_back(false);
        derivative_directions.push_back(-1);
        
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(2);
        use_derivative.push_back(true);
        derivative_directions.push_back(2);
        
        std::vector<double> u_ddz_w_avg = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
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
        
        std::vector<std::vector<double> >& u_ddz_w_avg_realizations = d_ensemble_statistics->u_ddz_w_avg_realizations;
        u_ddz_w_avg_realizations.push_back(u_ddz_w_avg);
        
        d_ensemble_statistics->u_ddz_w_avg_computed = true;
    }
    
    // Compute tau11_avg.
    
    std::vector<double> tau11_avg = getAveragedShearStressComponentWithInhomogeneousXDirection(
        0,
        patch_hierarchy,
        data_context);
    
    std::vector<std::vector<double> >& tau11_avg_realizations = d_ensemble_statistics->tau11_avg_realizations;
    tau11_avg_realizations.push_back(tau11_avg);
    
    d_ensemble_statistics->tau11_avg_computed = true;
    
    // Compute tau12_avg.
    
    if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
    {
        std::vector<double> tau12_avg = getAveragedShearStressComponentWithInhomogeneousXDirection(
            1,
            patch_hierarchy,
            data_context);
        
        std::vector<std::vector<double> >& tau12_avg_realizations = d_ensemble_statistics->tau12_avg_realizations;
        tau12_avg_realizations.push_back(tau12_avg);
        
        d_ensemble_statistics->tau12_avg_computed = true;
    }
    
    // Compute tau13_avg.
    
    if (d_dim == tbox::Dimension(3))
    {
        std::vector<double> tau13_avg = getAveragedShearStressComponentWithInhomogeneousXDirection(
            2,
            patch_hierarchy,
            data_context);
        
        std::vector<std::vector<double> >& tau13_avg_realizations = d_ensemble_statistics->tau13_avg_realizations;
        tau13_avg_realizations.push_back(tau13_avg);
        
        d_ensemble_statistics->tau13_avg_computed = true;
    }
    
    // Compute tau22_avg.
    
    if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
    {
        std::vector<double> tau22_avg = getAveragedShearStressComponentWithInhomogeneousXDirection(
            3,
            patch_hierarchy,
            data_context);
        
        std::vector<std::vector<double> >& tau22_avg_realizations = d_ensemble_statistics->tau22_avg_realizations;
        tau22_avg_realizations.push_back(tau22_avg);
        
        d_ensemble_statistics->tau22_avg_computed = true;
    }
    
    // Compute tau23_avg.
    
    if (d_dim == tbox::Dimension(3))
    {
        std::vector<double> tau23_avg = getAveragedShearStressComponentWithInhomogeneousXDirection(
            4,
            patch_hierarchy,
            data_context);
        
        std::vector<std::vector<double> >& tau23_avg_realizations = d_ensemble_statistics->tau23_avg_realizations;
        tau23_avg_realizations.push_back(tau23_avg);
        
        d_ensemble_statistics->tau23_avg_computed = true;
    }
    
    // Compute tau33_avg.
    
    if (d_dim == tbox::Dimension(3))
    {
        std::vector<double> tau33_avg = getAveragedShearStressComponentWithInhomogeneousXDirection(
            5,
            patch_hierarchy,
            data_context);
        
        std::vector<std::vector<double> >& tau33_avg_realizations = d_ensemble_statistics->tau33_avg_realizations;
        tau33_avg_realizations.push_back(tau33_avg);
        
        d_ensemble_statistics->tau33_avg_computed = true;
    }
    
    // Compute u_tau11_avg.
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    use_derivative.push_back(false);
    derivative_directions.push_back(-1);
    use_reciprocal.push_back(false);
    
    std::vector<double> u_tau11_avg = getAveragedQuantityWithShearStressComponentWithInhomogeneousXDirection(
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
    
    std::vector<std::vector<double> >& u_tau11_avg_realizations = d_ensemble_statistics->u_tau11_avg_realizations;
    u_tau11_avg_realizations.push_back(u_tau11_avg);
    
    d_ensemble_statistics->u_tau11_avg_computed = true;
    
    // Compute v_tau12_avg.
    
    if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
    {
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(1);
        use_derivative.push_back(false);
        derivative_directions.push_back(-1);
        use_reciprocal.push_back(false);
        
        std::vector<double> v_tau12_avg = getAveragedQuantityWithShearStressComponentWithInhomogeneousXDirection(
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
        
        std::vector<std::vector<double> >& v_tau12_avg_realizations = d_ensemble_statistics->v_tau12_avg_realizations;
        v_tau12_avg_realizations.push_back(v_tau12_avg);
        
        d_ensemble_statistics->v_tau12_avg_computed = true;
    }
    
    // Compute w_tau13_avg.
    
    if (d_dim == tbox::Dimension(3))
    {
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(2);
        use_derivative.push_back(false);
        derivative_directions.push_back(-1);
        use_reciprocal.push_back(false);
        
        std::vector<double> w_tau13_avg = getAveragedQuantityWithShearStressComponentWithInhomogeneousXDirection(
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
        
        std::vector<std::vector<double> >& w_tau13_avg_realizations = d_ensemble_statistics->w_tau13_avg_realizations;
        w_tau13_avg_realizations.push_back(w_tau13_avg);
        
        d_ensemble_statistics->w_tau13_avg_computed = true;
    }
    
    // Compute ddx_tau11_avg.
    
    std::vector<double> ddx_tau11_avg = getAveragedDerivativeOfShearStressComponentWithInhomogeneousXDirection(
        0,
        0,
        patch_hierarchy,
        data_context);
    
    std::vector<std::vector<double> >& ddx_tau11_avg_realizations = d_ensemble_statistics->ddx_tau11_avg_realizations;
    ddx_tau11_avg_realizations.push_back(ddx_tau11_avg);
    
    d_ensemble_statistics->ddx_tau11_avg_computed = true;
    
    // Compute ddy_tau12_avg.
    
    if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
    {
        std::vector<double> ddy_tau12_avg = getAveragedDerivativeOfShearStressComponentWithInhomogeneousXDirection(
            1,
            1,
            patch_hierarchy,
            data_context);
        
        std::vector<std::vector<double> >& ddy_tau12_avg_realizations = d_ensemble_statistics->ddy_tau12_avg_realizations;
        ddy_tau12_avg_realizations.push_back(ddy_tau12_avg);
        
        d_ensemble_statistics->ddy_tau12_avg_computed = true;
    }
    
    // Compute ddz_tau13_avg.
    
    if (d_dim == tbox::Dimension(3))
    {
        std::vector<double> ddz_tau13_avg = getAveragedDerivativeOfShearStressComponentWithInhomogeneousXDirection(
            2,
            2,
            patch_hierarchy,
            data_context);
        
        std::vector<std::vector<double> >& ddz_tau13_avg_realizations = d_ensemble_statistics->ddz_tau13_avg_realizations;
        ddz_tau13_avg_realizations.push_back(ddz_tau13_avg);
        
        d_ensemble_statistics->ddz_tau13_avg_computed = true;
    }
    
    // Compute tau11_ddx_u_avg.
    
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
    
    std::vector<std::vector<double> >& tau11_ddx_u_avg_realizations = d_ensemble_statistics->tau11_ddx_u_avg_realizations;
    tau11_ddx_u_avg_realizations.push_back(tau11_ddx_u_avg);
    
    d_ensemble_statistics->tau11_ddx_u_avg_computed = true;
    
    // Compute tau12_ddy_u_avg.
    
    if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
    {
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
        
        std::vector<std::vector<double> >& tau12_ddy_u_avg_realizations = d_ensemble_statistics->tau12_ddy_u_avg_realizations;
        tau12_ddy_u_avg_realizations.push_back(tau12_ddy_u_avg);
        
        d_ensemble_statistics->tau12_ddy_u_avg_computed = true;
    }
    
    // Compute tau13_ddz_u_avg.
    
    if (d_dim == tbox::Dimension(3))
    {
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
        
        std::vector<std::vector<double> >& tau13_ddz_u_avg_realizations = d_ensemble_statistics->tau13_ddz_u_avg_realizations;
        tau13_ddz_u_avg_realizations.push_back(tau13_ddz_u_avg);
        
        d_ensemble_statistics->tau13_ddz_u_avg_computed = true;
    }
    
    // Compute tau12_ddx_v_avg.
    
    if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
    {
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
        
        std::vector<std::vector<double> >& tau12_ddx_v_avg_realizations = d_ensemble_statistics->tau12_ddx_v_avg_realizations;
        tau12_ddx_v_avg_realizations.push_back(tau12_ddx_v_avg);
        
        d_ensemble_statistics->tau12_ddx_v_avg_computed = true;
    }
    
    // Compute tau22_ddy_v_avg.
    
    if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
    {
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
        
        std::vector<std::vector<double> >& tau22_ddy_v_avg_realizations = d_ensemble_statistics->tau22_ddy_v_avg_realizations;
        tau22_ddy_v_avg_realizations.push_back(tau22_ddy_v_avg);
        
        d_ensemble_statistics->tau22_ddy_v_avg_computed = true;
    }
    
    // Compute tau23_ddz_v_avg.
    
    if (d_dim == tbox::Dimension(3))
    {
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
        
        std::vector<std::vector<double> >& tau23_ddz_v_avg_realizations = d_ensemble_statistics->tau23_ddz_v_avg_realizations;
        tau23_ddz_v_avg_realizations.push_back(tau23_ddz_v_avg);
        
        d_ensemble_statistics->tau23_ddz_v_avg_computed = true;
    }
    
    // Compute tau13_ddx_w_avg.
    
    if (d_dim == tbox::Dimension(3))
    {
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
        
        std::vector<std::vector<double> >& tau13_ddx_w_avg_realizations = d_ensemble_statistics->tau13_ddx_w_avg_realizations;
        tau13_ddx_w_avg_realizations.push_back(tau13_ddx_w_avg);
        
        d_ensemble_statistics->tau13_ddx_w_avg_computed = true;
    }
    
    // Compute tau23_ddy_w_avg.
    
    if (d_dim == tbox::Dimension(3))
    {
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
        
        std::vector<std::vector<double> >& tau23_ddy_w_avg_realizations = d_ensemble_statistics->tau23_ddy_w_avg_realizations;
        tau23_ddy_w_avg_realizations.push_back(tau23_ddy_w_avg);
        
        d_ensemble_statistics->tau23_ddy_w_avg_computed = true;
    }
    
    // Compute tau33_ddz_w_avg.
    
    if (d_dim == tbox::Dimension(3))
    {
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
        
        std::vector<std::vector<double> >& tau33_ddz_w_avg_realizations = d_ensemble_statistics->tau33_ddz_w_avg_realizations;
        tau33_ddz_w_avg_realizations.push_back(tau33_ddz_w_avg);
        
        d_ensemble_statistics->tau33_ddz_w_avg_computed = true;
    }
    
    // Compute rho_inv_ddx_tau11_avg.
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_derivative.push_back(false);
    derivative_directions.push_back(-1);
    use_reciprocal.push_back(true);
    
    std::vector<double> rho_inv_ddx_tau11_avg =
        getAveragedQuantityWithDerivativeOfShearStressComponentWithInhomogeneousXDirection(
            quantity_names,
            component_indices,
            use_derivative,
            derivative_directions,
            use_reciprocal,
            0,
            0,
            patch_hierarchy,
            data_context);
    
    quantity_names.clear();
    component_indices.clear();
    use_derivative.clear();
    derivative_directions.clear();
    use_reciprocal.clear();
    
    std::vector<std::vector<double> >& rho_inv_ddx_tau11_avg_realizations = d_ensemble_statistics->rho_inv_ddx_tau11_avg_realizations;
    rho_inv_ddx_tau11_avg_realizations.push_back(rho_inv_ddx_tau11_avg);
    
    d_ensemble_statistics->rho_inv_ddx_tau11_avg_computed = true;
    
    // Compute rho_inv_ddy_tau12_avg.
    
    if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
    {
        quantity_names.push_back("DENSITY");
        component_indices.push_back(0);
        use_derivative.push_back(false);
        derivative_directions.push_back(-1);
        use_reciprocal.push_back(true);
        
        std::vector<double> rho_inv_ddy_tau12_avg =
            getAveragedQuantityWithDerivativeOfShearStressComponentWithInhomogeneousXDirection(
                quantity_names,
                component_indices,
                use_derivative,
                derivative_directions,
                use_reciprocal,
                1,
                1,
                patch_hierarchy,
                data_context);
        
        quantity_names.clear();
        component_indices.clear();
        use_derivative.clear();
        derivative_directions.clear();
        use_reciprocal.clear();
        
        std::vector<std::vector<double> >& rho_inv_ddy_tau12_avg_realizations = d_ensemble_statistics->rho_inv_ddy_tau12_avg_realizations;
        rho_inv_ddy_tau12_avg_realizations.push_back(rho_inv_ddy_tau12_avg);
        
        d_ensemble_statistics->rho_inv_ddy_tau12_avg_computed = true;
    }
    
    // Compute rho_inv_ddz_tau13_avg.
    
    if (d_dim == tbox::Dimension(3))
    {
        quantity_names.push_back("DENSITY");
        component_indices.push_back(0);
        use_derivative.push_back(false);
        derivative_directions.push_back(-1);
        use_reciprocal.push_back(true);
        
        std::vector<double> rho_inv_ddz_tau13_avg =
            getAveragedQuantityWithDerivativeOfShearStressComponentWithInhomogeneousXDirection(
                quantity_names,
                component_indices,
                use_derivative,
                derivative_directions,
                use_reciprocal,
                2,
                2,
                patch_hierarchy,
                data_context);
        
        quantity_names.clear();
        component_indices.clear();
        use_derivative.clear();
        derivative_directions.clear();
        use_reciprocal.clear();
        
        std::vector<std::vector<double> >& rho_inv_ddz_tau13_avg_realizations = d_ensemble_statistics->rho_inv_ddz_tau13_avg_realizations;
        rho_inv_ddz_tau13_avg_realizations.push_back(rho_inv_ddz_tau13_avg);
        
        d_ensemble_statistics->rho_inv_ddz_tau13_avg_computed = true;
    }
}


/*
 * Output spatial profile of ensemble averaged grid level number with assumed homogeneity in y-direction (2D) or
 * yz-plane (3D) to a file.
 */
void
RTIRMIBudgetsUtilities::outputSpatialProfileEnsembleAveragedGridLevelNumberWithHomogeneityInYDirectionOrInYZPlane(
    const std::string& stat_dump_filename,
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
        
        const std::vector<std::vector<double> >& grid_level_num_avg_realizations =
            d_ensemble_statistics->grid_level_num_avg_realizations;
        
        const int num_realizations = d_ensemble_statistics->getNumberOfEnsembles();
        
        TBOX_ASSERT(num_realizations > 0);
        TBOX_ASSERT(num_realizations == static_cast<int>(grid_level_num_avg_realizations.size()));
        
        const int num_cells = static_cast<int>(grid_level_num_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> grid_level_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                grid_level_avg_global[i] += weight*grid_level_num_avg_realizations[ri][i];
            }
        }
        
        f_out.write((char*)&output_time, sizeof(double));
        f_out.write((char*)&grid_level_avg_global[0], sizeof(double)*grid_level_avg_global.size());
        
        f_out.close();
    }
}


/*
 * Output spatial profile of ensemble averaged mass fraction with assumed homogeneity in y-direction (2D) or
 * yz-plane (3D) to a file.
 */
void
RTIRMIBudgetsUtilities::outputSpatialProfileEnsembleAveragedMassFractionWithHomogeneityInYDirectionOrInYZPlane(
    const std::string& stat_dump_filename,
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
        
        const int num_realizations = d_ensemble_statistics->getNumberOfEnsembles();
        
        TBOX_ASSERT(num_realizations > 0);
        TBOX_ASSERT(num_realizations == static_cast<int>(Y_0_avg_realizations.size()));
        
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
 * Output spatial profile of ensemble averaged density with assumed homogeneity in y-direction (2D) or
 * yz-plane (3D) to a file.
 */
void
RTIRMIBudgetsUtilities::outputSpatialProfileEnsembleAveragedDensityWithHomogeneityInYDirectionOrInYZPlane(
    const std::string& stat_dump_filename,
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
        
        const int num_realizations = d_ensemble_statistics->getNumberOfEnsembles();
        
        TBOX_ASSERT(num_realizations > 0);
        TBOX_ASSERT(num_realizations == static_cast<int>(rho_avg_realizations.size()));
        
        const int num_cells = static_cast<int>(rho_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> rho_avg(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                rho_avg[i] += weight*rho_avg_realizations[ri][i];
            }
        }
        
        f_out.write((char*)&output_time, sizeof(double));
        f_out.write((char*)&rho_avg[0], sizeof(double)*rho_avg.size());
        
        f_out.close();
    }
}


/*
 * Output spatial profile of ensemble averaged velocity x-component with assumed homogeneity in y-direction (2D)
 * or yz-plane (3D) to a file.
 */
void
RTIRMIBudgetsUtilities::outputSpatialProfileEnsembleAveragedVelocityXWithHomogeneityInYDirectionOrInYZPlane(
    const std::string& stat_dump_filename,
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
        
        const int num_realizations = d_ensemble_statistics->getNumberOfEnsembles();
        
        TBOX_ASSERT(num_realizations > 0);
        TBOX_ASSERT(num_realizations == static_cast<int>(u_avg_realizations.size()));
        
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
 * Output spatial profile of ensemble Favre averaged velocity x-component with assumed homogeneity in
 * y-direction (2D) or yz-plane (3D) to a file.
 */
void
RTIRMIBudgetsUtilities::outputSpatialProfileEnsembleFavreAveragedVelocityXWithHomogeneityInYDirectionOrInYZPlane(
    const std::string& stat_dump_filename,
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
        
        const int num_realizations = d_ensemble_statistics->getNumberOfEnsembles();
        
        TBOX_ASSERT(num_realizations > 0);
        TBOX_ASSERT(num_realizations == static_cast<int>(rho_avg_realizations.size()));
        TBOX_ASSERT(num_realizations == static_cast<int>(rho_u_avg_realizations.size()));
        
        const int num_cells = static_cast<int>(rho_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> rho_avg(num_cells, double(0));
        std::vector<double> rho_u_avg(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                rho_avg[i]   += weight*rho_avg_realizations[ri][i];
                rho_u_avg[i] += weight*rho_u_avg_realizations[ri][i];
            }
        }
        
        std::vector<double> u_tilde(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            u_tilde[i] = rho_u_avg[i]/rho_avg[i];
        }
        
        f_out.write((char*)&output_time, sizeof(double));
        f_out.write((char*)&u_tilde[0], sizeof(double)*u_tilde.size());
        
        f_out.close();
    }
}


/*
 * Output spatial profile of ensemble turbulent mass flux velocity in x-direction with assumed homogeneity in
 * y-direction (2D) or yz-plane (3D) to a file.
 */
void
RTIRMIBudgetsUtilities::outputSpatialProfileEnsembleTurbulentMassFluxVelocityXWithHomogeneityInYDirectionOrInYZPlane(
    const std::string& stat_dump_filename,
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
        
        const int num_realizations = d_ensemble_statistics->getNumberOfEnsembles();
        
        TBOX_ASSERT(num_realizations > 0);
        TBOX_ASSERT(num_realizations == static_cast<int>(rho_avg_realizations.size()));
        TBOX_ASSERT(num_realizations == static_cast<int>(u_avg_realizations.size()));
        TBOX_ASSERT(num_realizations == static_cast<int>(rho_u_avg_realizations.size()));
        
        const int num_cells = static_cast<int>(rho_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> rho_avg(num_cells, double(0));
        std::vector<double> u_avg(num_cells, double(0));
        std::vector<double> rho_u_avg(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                rho_avg[i]   += weight*rho_avg_realizations[ri][i];
                u_avg[i]     += weight*u_avg_realizations[ri][i];
                rho_u_avg[i] += weight*rho_u_avg_realizations[ri][i];
            }
        }
        
        std::vector<double> rho_p_u_p(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            rho_p_u_p[i] = rho_u_avg[i] - rho_avg[i]*u_avg[i];
        }
        
        std::vector<double> a1(rho_p_u_p);
        for (int i = 0; i < num_cells; i++)
        {
            a1[i] /= rho_avg[i];
        }
        
        f_out.write((char*)&output_time, sizeof(double));
        f_out.write((char*)&a1[0], sizeof(double)*a1.size());
        
        f_out.close();
    }
}


/*
 * Output spatial profile of ensemble averaged pressure gradient in x-direction with assumed homogeneity in
 * y-direction (2D) or yz-plane (3D) to a file.
 */
void
RTIRMIBudgetsUtilities::outputSpatialProfileEnsembleAveragedPressureGradientXWithHomogeneityInYDirectionOrInYZPlane(
    const std::string& stat_dump_filename,
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
        
        const std::vector<std::vector<double> >& ddx_p_avg_realizations =
            d_ensemble_statistics->ddx_p_avg_realizations;
        
        const int num_realizations = d_ensemble_statistics->getNumberOfEnsembles();
        
        TBOX_ASSERT(num_realizations > 0);
        TBOX_ASSERT(num_realizations == static_cast<int>(ddx_p_avg_realizations.size()));
        
        const int num_cells = static_cast<int>(ddx_p_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> ddx_p_avg(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                ddx_p_avg[i] += weight*ddx_p_avg_realizations[ri][i];
            }
        }
        
        f_out.write((char*)&output_time, sizeof(double));
        f_out.write((char*)&ddx_p_avg[0], sizeof(double)*ddx_p_avg.size());
        
        f_out.close();
    }
}


/*
 * Output budget of Favre mean TKE with inhomogeneous x-direction to a file.
 */
void
RTIRMIBudgetsUtilities::outputBudgetFavreMeanTKEWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    std::ofstream f_out;
    
    MPIHelper MPI_helper = MPIHelper(
        "MPI_helper",
        d_dim,
        d_grid_geometry,
        patch_hierarchy);
    
    const std::vector<double>& dx_vec = MPI_helper.getFinestRefinedDomainGridSpacing();
    const double dx = dx_vec[0];
    
    /*
     * Output the spatial profiles (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        const std::vector<std::vector<double> >& rho_avg_realizations     = d_ensemble_statistics->rho_avg_realizations;
        const std::vector<std::vector<double> >& u_avg_realizations       = d_ensemble_statistics->u_avg_realizations;
        const std::vector<std::vector<double> >& rho_u_avg_realizations   = d_ensemble_statistics->rho_u_avg_realizations;
        const std::vector<std::vector<double> >& rho_v_avg_realizations   = d_ensemble_statistics->rho_v_avg_realizations;
        const std::vector<std::vector<double> >& rho_w_avg_realizations   = d_ensemble_statistics->rho_w_avg_realizations;
        const std::vector<std::vector<double> >& rho_u_u_avg_realizations = d_ensemble_statistics->rho_u_u_avg_realizations;
        
        const std::vector<std::vector<double> >& ddx_rho_avg_realizations     = d_ensemble_statistics->ddx_rho_avg_realizations;
        const std::vector<std::vector<double> >& ddx_p_avg_realizations       = d_ensemble_statistics->ddx_p_avg_realizations;
        const std::vector<std::vector<double> >& ddx_rho_u_avg_realizations   = d_ensemble_statistics->ddx_rho_u_avg_realizations;
        const std::vector<std::vector<double> >& ddx_rho_u_u_avg_realizations = d_ensemble_statistics->ddx_rho_u_u_avg_realizations;
        
        const std::vector<std::vector<double> >& ddx_tau11_avg_realizations = d_ensemble_statistics->ddx_tau11_avg_realizations;
        
        const int num_realizations = d_ensemble_statistics->getNumberOfEnsembles();
        
        TBOX_ASSERT(num_realizations > 0);
        TBOX_ASSERT(num_realizations == static_cast<int>(rho_avg_realizations.size()));
        
        const int num_cells = static_cast<int>(rho_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> rho_avg_global(num_cells, double(0));
        std::vector<double> u_avg_global(num_cells, double(0));
        std::vector<double> rho_u_avg_global(num_cells, double(0));
        std::vector<double> rho_v_avg_global(num_cells, double(0));
        std::vector<double> rho_w_avg_global(num_cells, double(0));
        std::vector<double> rho_u_u_avg_global(num_cells, double(0));
        
        std::vector<double> ddx_rho_avg_global(num_cells, double(0));
        std::vector<double> ddx_p_avg_global(num_cells, double(0));
        std::vector<double> ddx_rho_u_avg_global(num_cells, double(0));
        std::vector<double> ddx_rho_u_u_avg_global(num_cells, double(0));
        
        std::vector<double> ddx_tau11_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                rho_avg_global[i]     += weight*rho_avg_realizations[ri][i];
                u_avg_global[i]       += weight*u_avg_realizations[ri][i];
                rho_u_avg_global[i]   += weight*rho_u_avg_realizations[ri][i];
                rho_u_u_avg_global[i] += weight*rho_u_u_avg_realizations[ri][i];
                
                ddx_rho_avg_global[i]     += weight*ddx_rho_avg_realizations[ri][i];
                ddx_p_avg_global[i]       += weight*ddx_p_avg_realizations[ri][i];
                ddx_rho_u_avg_global[i]   += weight*ddx_rho_u_avg_realizations[ri][i];
                ddx_rho_u_u_avg_global[i] += weight*ddx_rho_u_u_avg_realizations[ri][i];
                
                ddx_tau11_avg_global[i] += weight*ddx_tau11_avg_realizations[ri][i];
            }
        }
        
        if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
        {
            for (int ri = 0; ri < num_realizations; ri++)
            {
                for (int i = 0; i < num_cells; i++)
                {
                    rho_v_avg_global[i] += weight*rho_v_avg_realizations[ri][i];
                }
            }
        }
        
        if (d_dim == tbox::Dimension(3))
        {
            for (int ri = 0; ri < num_realizations; ri++)
            {
                for (int i = 0; i < num_cells; i++)
                {
                    rho_w_avg_global[i] += weight*rho_w_avg_realizations[ri][i];
                }
            }
        }
        
        /*
         * Compute u_tilde, v_tilde and w_tilde.
         */
        
        std::vector<double> u_tilde(rho_u_avg_global);
        std::vector<double> v_tilde(rho_v_avg_global);
        std::vector<double> w_tilde(rho_w_avg_global);
        
        for (int i = 0; i < num_cells; i++)
        {
            u_tilde[i] /= rho_avg_global[i];
        }
        
        if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
        {
            for (int i = 0; i < num_cells; i++)
            {
                v_tilde[i] /= rho_avg_global[i];
            }
        }
        
        if (d_dim == tbox::Dimension(3))
        {
            for (int i = 0; i < num_cells; i++)
            {
                w_tilde[i] /= rho_avg_global[i];
            }
        }
        
        /*
         * Compute a1.
         */
        
        std::vector<double> rho_p_u_p(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            rho_p_u_p[i] = rho_u_avg_global[i] - rho_avg_global[i]*u_avg_global[i];
        }
        
        std::vector<double> a1(rho_p_u_p);
        for (int i = 0; i < num_cells; i++)
        {
            a1[i] /= rho_avg_global[i];
        }
        
        /*
         * Compute R11.
         */
        
        std::vector<double> rho_R11(num_cells, double(0));
        std::vector<double> R11(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            const double u_tilde       = rho_u_avg_global[i]/rho_avg_global[i];
            const double rho_u_pp_u_pp = rho_u_u_avg_global[i] - rho_u_avg_global[i]*u_tilde;
            
            rho_R11[i] = rho_u_pp_u_pp;
            R11[i]     = rho_u_pp_u_pp/rho_avg_global[i];
        }
        
        /*
         * Compute K.
         */
        
        std::vector<double> rho_K(num_cells, double(0));
        std::vector<double> K(num_cells, double(0));
        
        for (int i = 0; i < num_cells; i++)
        {
            K[i] += u_tilde[i]*u_tilde[i];
        }
        
        if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
        {
            for (int i = 0; i < num_cells; i++)
            {
                K[i] += v_tilde[i]*v_tilde[i];
            }
        }
        
        if (d_dim == tbox::Dimension(3))
        {
            for (int i = 0; i < num_cells; i++)
            {
                K[i] += w_tilde[i]*w_tilde[i];
            }
        }
        
        for (int i = 0; i < num_cells; i++)
        {
            K[i]     *= double(1)/double(2);
            rho_K[i]  = rho_avg_global[i]*K[i]; 
        }
        
        /*
         * Compute term II.
         */
        
        std::vector<double> rho_u_tilde_K(rho_u_avg_global);
        
        for (int i = 0; i < num_cells; i++)
        {
            rho_u_tilde_K[i] *= K[i];
        }
        
        std::vector<double> ddx_rho_u_tilde_K = computeDerivativeOfVector1D(
            rho_u_tilde_K,
            dx);
        
        /*
         * Compute term II in moving frame of mixing layer.
         */
        
        std::vector<double> rho_a1_K(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            rho_a1_K[i] = rho_avg_global[i]*a1[i]*K[i];
        }
        
        std::vector<double> ddx_rho_a1_K = computeDerivativeOfVector1D(
            rho_a1_K,
            dx);
        
        /*
         * Compute term III.
         */
        
        std::vector<double> g = getGravityVector();
        
        std::vector<double> rho_ui_gi(num_cells, double(0));
        
        for (int i = 0; i < num_cells; i++)
        {
            rho_ui_gi[i] += u_tilde[i]*g[0];
        }
        
        if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
        {
            for (int i = 0; i < num_cells; i++)
            {
                rho_ui_gi[i] += v_tilde[i]*g[1];
            }
        }
        
        if (d_dim == tbox::Dimension(3))
        {
            for (int i = 0; i < num_cells; i++)
            {
                rho_ui_gi[i] += w_tilde[i]*g[2];
            }
        }
        
        for (int i = 0; i < num_cells; i++)
        {
            rho_ui_gi[i] *= rho_avg_global[i];
        }
        
        /*
         * Compute term III in moving frame of mixing layer.
         */
        
        std::vector<double> rho_a1_g1(num_cells, double(0));
        
        for (int i = 0; i < num_cells; i++)
        {
            rho_a1_g1[i] += rho_avg_global[i]*a1[i]*g[0];
        }
        
        /*
         * Compute term IV.
         */
        
        std::vector<double> m_ddx_rho_u_tilde_R11(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            const double ddx_R11_tilde = -(rho_R11[i]/(rho_avg_global[i]*rho_avg_global[i]))*ddx_rho_avg_global[i] +
                double(1)/rho_avg_global[i]*(ddx_rho_u_u_avg_global[i] - double(2)*u_tilde[i]*ddx_rho_u_avg_global[i] +
                u_tilde[i]*u_tilde[i]*ddx_rho_avg_global[i]);
            
            m_ddx_rho_u_tilde_R11[i] = -(rho_u_avg_global[i]*ddx_R11_tilde + R11[i]*ddx_rho_u_avg_global[i]);
        }
        
        /*
         * Compute term IV in moving frame of mixing layer.
         */
        
        std::vector<double> m_rho_a1_R11(rho_p_u_p);
        for (int i = 0; i < num_cells; i++)
        {
            m_rho_a1_R11[i] *= (-R11[i]);
        }
        
        std::vector<double> m_ddx_rho_a1_R11 = computeDerivativeOfVector1D(
            m_rho_a1_R11,
            dx);
        
        /*
         * Compute term V(1).
         */
        
        std::vector<double> m_u_tilde_ddx_p(ddx_p_avg_global);
        for (int i = 0; i < num_cells; i++)
        {
            m_u_tilde_ddx_p[i] *= (-u_tilde[i]);
        }
        
        /*
         * Compute term V(1) in moving frame of mixing layer
         */
        
        std::vector<double> m_a1_ddx_p(ddx_p_avg_global);
        for (int i = 0; i < num_cells; i++)
        {
            m_a1_ddx_p[i] *= (-a1[i]);
        }
        
        /*
         * Compute term V(2).
         */
    
        std::vector<double> ddx_u_tilde(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            ddx_u_tilde[i] = ddx_rho_u_avg_global[i]/rho_avg_global[i] -
                rho_u_avg_global[i]/(rho_avg_global[i]*rho_avg_global[i])*ddx_rho_avg_global[i];
        }
        
        std::vector<double> rho_R11_ddx_u_tilde(ddx_u_tilde);
        for (int i = 0; i < num_cells; i++)
        {
            rho_R11_ddx_u_tilde[i] *= rho_R11[i];
        }
        
        /*
         * Compute term VI.
         */
        
        std::vector<double> u_tilde_ddx_tau11(ddx_tau11_avg_global);
        for (int i = 0; i < num_cells; i++)
        {
            u_tilde_ddx_tau11[i] *= u_tilde[i];
        }
        
        /*
         * Output budget.
         */
        
        f_out.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        f_out.write((char*)&output_time, sizeof(double));
        f_out.write((char*)&rho_K[0], sizeof(double)*rho_K.size());
        // Term II.
        f_out.write((char*)&ddx_rho_u_tilde_K[0], sizeof(double)*ddx_rho_u_tilde_K.size());
        
        // Term III.
        f_out.write((char*)&rho_ui_gi[0], sizeof(double)*rho_ui_gi.size());
        
        // Term IV.
        f_out.write((char*)&m_ddx_rho_u_tilde_R11[0], sizeof(double)*m_ddx_rho_u_tilde_R11.size());
        
        // Term V(1).
        f_out.write((char*)&m_u_tilde_ddx_p[0], sizeof(double)*m_u_tilde_ddx_p.size());
        // Term V(2).
        f_out.write((char*)&rho_R11_ddx_u_tilde[0], sizeof(double)*rho_R11_ddx_u_tilde.size());
        
        // Term VI.
        f_out.write((char*)&u_tilde_ddx_tau11[0], sizeof(double)*u_tilde_ddx_tau11.size());
        
        // Term II in moving frame of mixing layer.
        f_out.write((char*)&ddx_rho_a1_K[0], sizeof(double)*ddx_rho_a1_K.size());
        
        // Term III in moving frame of mixing layer.
        f_out.write((char*)&rho_a1_g1[0], sizeof(double)*rho_a1_g1.size());
        
        // Term IV in moving frame of mixing layer.
        f_out.write((char*)&m_rho_a1_R11[0], sizeof(double)*m_rho_a1_R11.size());
        
        // Term V(1) in moving frame of mixing layer.
        f_out.write((char*)&m_a1_ddx_p[0], sizeof(double)*m_a1_ddx_p.size());
        
        // // a1.
        // f_out.write((char*)&a1[0], sizeof(double)*a1.size());
        
        // // u_avg.
        // f_out.write((char*)&u_avg_global[0], sizeof(double)*u_avg_global.size());
        
        // // u_tilde.
        // f_out.write((char*)&u_tilde[0], sizeof(double)*u_tilde.size());
        
        // // ddx_p_avg.
        // f_out.write((char*)&ddx_p_avg_global[0], sizeof(double)*ddx_p_avg_global.size());
        
        f_out.close();
    }
}


/*
 * Output budget of turbulent mass flux in x-direction with inhomogeneous x-direction to a file.
 */
void
RTIRMIBudgetsUtilities::outputBudgetTurbMassFluxXWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    std::ofstream f_out;
    
    MPIHelper MPI_helper = MPIHelper(
        "MPI_helper",
        d_dim,
        d_grid_geometry,
        patch_hierarchy);
    
    const std::vector<double>& dx_vec = MPI_helper.getFinestRefinedDomainGridSpacing();
    const double dx = dx_vec[0];
    
    /*
     * Output the spatial profiles (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        const std::vector<std::vector<double> >& rho_avg_realizations     = d_ensemble_statistics->rho_avg_realizations;
        const std::vector<std::vector<double> >& rho_inv_avg_realizations = d_ensemble_statistics->rho_inv_avg_realizations;
        const std::vector<std::vector<double> >& u_avg_realizations       = d_ensemble_statistics->u_avg_realizations;
        const std::vector<std::vector<double> >& u_sq_avg_realizations    = d_ensemble_statistics->u_sq_avg_realizations;
        const std::vector<std::vector<double> >& rho_u_avg_realizations   = d_ensemble_statistics->rho_u_avg_realizations;
        const std::vector<std::vector<double> >& rho_u_u_avg_realizations = d_ensemble_statistics->rho_u_u_avg_realizations;
        
        const std::vector<std::vector<double> >& ddx_rho_avg_realizations     = d_ensemble_statistics->ddx_rho_avg_realizations;
        const std::vector<std::vector<double> >& ddx_p_avg_realizations       = d_ensemble_statistics->ddx_p_avg_realizations;
        const std::vector<std::vector<double> >& ddx_u_avg_realizations       = d_ensemble_statistics->ddx_u_avg_realizations;
        const std::vector<std::vector<double> >& ddx_u_sq_avg_realizations    = d_ensemble_statistics->ddx_u_sq_avg_realizations;
        const std::vector<std::vector<double> >& ddx_rho_u_avg_realizations   = d_ensemble_statistics->ddx_rho_u_avg_realizations;
        const std::vector<std::vector<double> >& ddx_rho_u_u_avg_realizations = d_ensemble_statistics->ddx_rho_u_u_avg_realizations;
        
        const std::vector<std::vector<double> >& ddy_v_avg_realizations = d_ensemble_statistics->ddy_v_avg_realizations;
        const std::vector<std::vector<double> >& ddz_w_avg_realizations = d_ensemble_statistics->ddz_w_avg_realizations;
        
        const std::vector<std::vector<double> >& rho_inv_ddx_p_avg_realizations = d_ensemble_statistics->rho_inv_ddx_p_avg_realizations;
        
        const std::vector<std::vector<double> >& u_ddx_u_avg_realizations = d_ensemble_statistics->u_ddx_u_avg_realizations;
        const std::vector<std::vector<double> >& u_ddy_v_avg_realizations = d_ensemble_statistics->u_ddy_v_avg_realizations;
        const std::vector<std::vector<double> >& u_ddz_w_avg_realizations = d_ensemble_statistics->u_ddz_w_avg_realizations;
        
        const std::vector<std::vector<double> >& ddx_tau11_avg_realizations = d_ensemble_statistics->ddx_tau11_avg_realizations;
        const std::vector<std::vector<double> >& ddy_tau12_avg_realizations = d_ensemble_statistics->ddy_tau12_avg_realizations;
        const std::vector<std::vector<double> >& ddz_tau13_avg_realizations = d_ensemble_statistics->ddz_tau13_avg_realizations;
        
        const std::vector<std::vector<double> >& rho_inv_ddx_tau11_avg_realizations = d_ensemble_statistics->rho_inv_ddx_tau11_avg_realizations;
        const std::vector<std::vector<double> >& rho_inv_ddy_tau12_avg_realizations = d_ensemble_statistics->rho_inv_ddy_tau12_avg_realizations;
        const std::vector<std::vector<double> >& rho_inv_ddz_tau13_avg_realizations = d_ensemble_statistics->rho_inv_ddz_tau13_avg_realizations;
        
        const int num_realizations = d_ensemble_statistics->getNumberOfEnsembles();
        
        TBOX_ASSERT(num_realizations > 0);
        TBOX_ASSERT(num_realizations == static_cast<int>(rho_avg_realizations.size()));
        
        const int num_cells = static_cast<int>(rho_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> rho_avg_global(num_cells, double(0));
        std::vector<double> rho_inv_avg_global(num_cells, double(0));
        std::vector<double> u_avg_global(num_cells, double(0));
        std::vector<double> u_sq_avg_global(num_cells, double(0));
        std::vector<double> rho_u_avg_global(num_cells, double(0));
        std::vector<double> rho_u_u_avg_global(num_cells, double(0));
        
        std::vector<double> ddx_rho_avg_global(num_cells, double(0));
        std::vector<double> ddx_p_avg_global(num_cells, double(0));
        std::vector<double> ddx_u_avg_global(num_cells, double(0));
        std::vector<double> ddx_u_sq_avg_global(num_cells, double(0));
        std::vector<double> ddx_rho_u_avg_global(num_cells, double(0));
        std::vector<double> ddx_rho_u_u_avg_global(num_cells, double(0));
        
        std::vector<double> ddy_v_avg_global(num_cells, double(0));
        std::vector<double> ddz_w_avg_global(num_cells, double(0));
        
        std::vector<double> rho_inv_ddx_p_avg_global(num_cells, double(0));
        
        std::vector<double> u_ddx_u_avg_global(num_cells, double(0));
        std::vector<double> u_ddy_v_avg_global(num_cells, double(0));
        std::vector<double> u_ddz_w_avg_global(num_cells, double(0));
        
        std::vector<double> ddx_tau11_avg_global(num_cells, double(0));
        std::vector<double> ddy_tau12_avg_global(num_cells, double(0));
        std::vector<double> ddz_tau13_avg_global(num_cells, double(0));
        
        std::vector<double> rho_inv_ddx_tau11_avg_global(num_cells, double(0));
        std::vector<double> rho_inv_ddy_tau12_avg_global(num_cells, double(0));
        std::vector<double> rho_inv_ddz_tau13_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells; i++)
            {
                rho_avg_global[i]     += weight*rho_avg_realizations[ri][i];
                rho_inv_avg_global[i] += weight*rho_inv_avg_realizations[ri][i];
                u_avg_global[i]       += weight*u_avg_realizations[ri][i];
                u_sq_avg_global[i]    += weight*u_sq_avg_realizations[ri][i];
                rho_u_avg_global[i]   += weight*rho_u_avg_realizations[ri][i];
                rho_u_u_avg_global[i] += weight*rho_u_u_avg_realizations[ri][i];
                
                ddx_rho_avg_global[i]     += weight*ddx_rho_avg_realizations[ri][i];
                ddx_p_avg_global[i]       += weight*ddx_p_avg_realizations[ri][i];
                ddx_u_avg_global[i]       += weight*ddx_u_avg_realizations[ri][i];
                ddx_u_sq_avg_global[i]    += weight*ddx_u_sq_avg_realizations[ri][i];
                ddx_rho_u_avg_global[i]   += weight*ddx_rho_u_avg_realizations[ri][i];
                ddx_rho_u_u_avg_global[i] += weight*ddx_rho_u_u_avg_realizations[ri][i];
                
                
                rho_inv_ddx_p_avg_global[i] += weight*rho_inv_ddx_p_avg_realizations[ri][i];
                
                u_ddx_u_avg_global[i] += weight*u_ddx_u_avg_realizations[ri][i];
                
                ddx_tau11_avg_global[i] += weight*ddx_tau11_avg_realizations[ri][i];
                
                rho_inv_ddx_tau11_avg_global[i] += weight*rho_inv_ddx_tau11_avg_realizations[ri][i];
            }
        }
        
        if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
        {
            for (int ri = 0; ri < num_realizations; ri++)
            {
                for (int i = 0; i < num_cells; i++)
                {
                    ddy_v_avg_global[i] += weight*ddy_v_avg_realizations[ri][i];
                    
                    u_ddy_v_avg_global[i] += weight*u_ddy_v_avg_realizations[ri][i];
                    
                    ddy_tau12_avg_global[i] += weight*ddy_tau12_avg_realizations[ri][i];
                    
                    rho_inv_ddy_tau12_avg_global[i] += weight*rho_inv_ddy_tau12_avg_realizations[ri][i];
                }
            }
        }
        
        if (d_dim == tbox::Dimension(3))
        {
            for (int ri = 0; ri < num_realizations; ri++)
            {
                for (int i = 0; i < num_cells; i++)
                {
                    ddz_w_avg_global[i] += weight*ddz_w_avg_realizations[ri][i];
                    
                    u_ddz_w_avg_global[i] += weight*u_ddz_w_avg_realizations[ri][i];
                    
                    ddz_tau13_avg_global[i] += weight*ddz_tau13_avg_realizations[ri][i];
                    
                    rho_inv_ddz_tau13_avg_global[i] += weight*rho_inv_ddz_tau13_avg_realizations[ri][i];
                }
            }
        }
        
        /*
         * Compute rho_a1.
         */
        
        std::vector<double> rho_p_u_p(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            rho_p_u_p[i] = rho_u_avg_global[i] - rho_avg_global[i]*u_avg_global[i];
        }
        
        std::vector<double> a1(rho_p_u_p);
        for (int i = 0; i < num_cells; i++)
        {
            a1[i] /= rho_avg_global[i];
        }
        
        /*
         * Compute u_tilde.
         */
        
        std::vector<double> u_tilde(rho_u_avg_global);
        for (int i = 0; i < num_cells; i++)
        {
            u_tilde[i] /= rho_avg_global[i];
        }
        
        /*
         * Compute term II.
         */
        
        std::vector<double> ddx_a1(num_cells, double(0));
        std::vector<double> ddx_rho_u_tilde_a1(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            ddx_a1[i] = -rho_p_u_p[i]/(rho_avg_global[i]*rho_avg_global[i])*ddx_rho_avg_global[i] +
                double(1)/rho_avg_global[i]*(ddx_rho_u_avg_global[i] - u_avg_global[i]*ddx_rho_avg_global[i]) -
                ddx_u_avg_global[i];
            
            ddx_rho_u_tilde_a1[i] = rho_u_avg_global[i]*ddx_a1[i] + a1[i]*ddx_rho_u_avg_global[i];
        }
        
        /*
         * Compute term II in moving frame of mixing layer.
         */
        
        std::vector<double> rho_a1_a1(a1);
        for (int i = 0; i < num_cells; i++)
        {
            rho_a1_a1[i] *= rho_p_u_p[i];
        }
        
        std::vector<double> d_rho_a1_a1_dx = computeDerivativeOfVector1D(
            rho_a1_a1,
            dx);
        
        /*
         * Compute term III(1).
         */
        
        std::vector<double> b(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            b[i] = -(double(1) - rho_avg_global[i]*rho_inv_avg_global[i]);
        }
        
        std::vector<double> b_ddx_p(ddx_p_avg_global);
        for (int i = 0; i < num_cells; i++)
        {
            b_ddx_p[i] *= b[i];
        }
        
        /*
         * Compute term III(2).
         */
        
        std::vector<double> m_b_ddx_tau11(ddx_tau11_avg_global);
        for (int i = 0; i < num_cells; i++)
        {
            m_b_ddx_tau11[i] *= (-b[i]);
        }
        
        /*
         * Compute term III(3).
         */
        
        std::vector<double> R11(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            const double u_tilde = rho_u_avg_global[i]/rho_avg_global[i];
            const double rho_u_pp_u_pp = rho_u_u_avg_global[i] - rho_u_avg_global[i]*u_tilde;
            R11[i] = rho_u_pp_u_pp/rho_avg_global[i];
        }
        
        std::vector<double> m_R11_ddx_rho(ddx_rho_avg_global);
        for (int i = 0; i < num_cells; i++)
        {
            m_R11_ddx_rho[i] *= (-R11[i]);
        }
        
        /*
         * Compute term IV(1).
         */
        
        std::vector<double> rho_ddx_a1_sq(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            rho_ddx_a1_sq[i] = double(2)*rho_avg_global[i]*a1[i]*ddx_a1[i];
        }
        
        /*
         * Compute term IV(2).
         */
        
        std::vector<double> m_rho_a1_ddx_u(ddx_u_avg_global);
        for (int i = 0; i < num_cells; i++)
        {
            m_rho_a1_ddx_u[i] *= (-rho_avg_global[i]*a1[i]);
        }
        
        /*
         * Compute term V.
         */
        
        std::vector<double> rho_p_u_p_u_p(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            rho_p_u_p_u_p[i] = rho_u_u_avg_global[i] + double(2)*rho_avg_global[i]*u_avg_global[i]*u_avg_global[i] -
                rho_avg_global[i]*u_sq_avg_global[i] -
                double(2)*rho_u_avg_global[i]*u_avg_global[i];
        }
        
        std::vector<double> m_rho_ddx_rho_p_u_p_sq_over_rho(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            m_rho_ddx_rho_p_u_p_sq_over_rho[i] = -rho_avg_global[i]*(
                -rho_p_u_p_u_p[i]/(rho_avg_global[i]*rho_avg_global[i])*ddx_rho_avg_global[i] +
                double(1)/rho_avg_global[i]*(
                    ddx_rho_u_u_avg_global[i] - double(2)*rho_u_avg_global[i]*ddx_u_avg_global[i] -
                    double(2)*u_avg_global[i]*ddx_rho_u_avg_global[i] +
                    double(2)*u_avg_global[i]*u_avg_global[i]*ddx_rho_avg_global[i] +
                    double(4)*rho_avg_global[i]*u_avg_global[i]*ddx_u_avg_global[i] -
                    rho_avg_global[i]*ddx_u_sq_avg_global[i] - u_sq_avg_global[i]*ddx_rho_avg_global[i]
                ));
        }
        
        /*
         * Compute term VI(1).
         */
        
        std::vector<double> rho_rho_inv_p_ddx_p_p(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            rho_rho_inv_p_ddx_p_p[i] = rho_inv_ddx_p_avg_global[i] - rho_inv_avg_global[i]*ddx_p_avg_global[i];
        }
        
        for (int i = 0; i < num_cells; i++)
        {
            rho_rho_inv_p_ddx_p_p[i] *= rho_avg_global[i];
        }
        
        /*
         * Compute term VI(2).
         */
        
        std::vector<double> rho_inv_p_ddx_tau11_p(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            rho_inv_p_ddx_tau11_p[i] = rho_inv_ddx_tau11_avg_global[i] - rho_inv_avg_global[i]*ddx_tau11_avg_global[i];
        }
        
        std::vector<double> rho_inv_p_ddy_tau12_p(num_cells, double(0));
        if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
        {
            for (int i = 0; i < num_cells; i++)
            {
                rho_inv_p_ddy_tau12_p[i] = rho_inv_ddy_tau12_avg_global[i] - rho_inv_avg_global[i]*ddy_tau12_avg_global[i];
            }
        }
        
        std::vector<double> rho_inv_p_ddz_tau13_p(num_cells, double(0));
        if (d_dim == tbox::Dimension(3))
        {
            for (int i = 0; i < num_cells; i++)
            {
                rho_inv_p_ddz_tau13_p[i] = rho_inv_ddz_tau13_avg_global[i] - rho_inv_avg_global[i]*ddz_tau13_avg_global[i];
            }
        }
        
        std::vector<double> m_rho_rho_inv_p_ddx_tau_ij_p(num_cells, double(0));
        if (d_dim == tbox::Dimension(1))
        {
            for (int i = 0; i < num_cells; i++)
            {
                m_rho_rho_inv_p_ddx_tau_ij_p[i] = -rho_avg_global[i]*rho_inv_p_ddx_tau11_p[i];
            }
        }
        else if (d_dim == tbox::Dimension(2))
        {
            for (int i = 0; i < num_cells; i++)
            {
                m_rho_rho_inv_p_ddx_tau_ij_p[i] = -rho_avg_global[i]*(rho_inv_p_ddx_tau11_p[i] + rho_inv_p_ddy_tau12_p[i]);
            }
        }
        else if (d_dim == tbox::Dimension(3))
        {
            for (int i = 0; i < num_cells; i++)
            {
                m_rho_rho_inv_p_ddx_tau_ij_p[i] = -rho_avg_global[i]*(rho_inv_p_ddx_tau11_p[i] +
                    rho_inv_p_ddy_tau12_p[i] +
                    rho_inv_p_ddz_tau13_p[i]);
            }
        }
        
        /*
         * Compute term VI(3).
         */
        
        std::vector<double> epsilon_a1_1(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            epsilon_a1_1[i] = u_ddx_u_avg_global[i] - u_avg_global[i]*ddx_u_avg_global[i];
        }
        
        std::vector<double> epsilon_a1_2(num_cells, double(0));
        if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
        {
            for (int i = 0; i < num_cells; i++)
            {
                epsilon_a1_2[i] = u_ddy_v_avg_global[i] - u_avg_global[i]*ddy_v_avg_global[i];
            }
        }
        
        std::vector<double> epsilon_a1_3(num_cells, double(0));
        if (d_dim == tbox::Dimension(3))
        {
            for (int i = 0; i < num_cells; i++)
            {
                epsilon_a1_3[i] = u_ddz_w_avg_global[i] - u_avg_global[i]*ddz_w_avg_global[i];
            }
        }
        
        std::vector<double> rho_epsilon_a1(num_cells, double(0));
        if (d_dim == tbox::Dimension(1))
        {
            for (int i = 0; i < num_cells; i++)
            {
                rho_epsilon_a1[i] = -rho_avg_global[i]*epsilon_a1_1[i];
            }
        }
        else if (d_dim == tbox::Dimension(2))
        {
            for (int i = 0; i < num_cells; i++)
            {
                rho_epsilon_a1[i] = -rho_avg_global[i]*(epsilon_a1_1[i] + epsilon_a1_2[i]);
            }
        }
        else if (d_dim == tbox::Dimension(3))
        {
            for (int i = 0; i < num_cells; i++)
            {
                rho_epsilon_a1[i] = -rho_avg_global[i]*(epsilon_a1_1[i] + epsilon_a1_2[i] + epsilon_a1_3[i]);
            }
        }
        
        /*
         * Output budget.
         */
        
        f_out.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        f_out.write((char*)&output_time, sizeof(double));
        f_out.write((char*)&rho_p_u_p[0], sizeof(double)*rho_p_u_p.size());
        // Term II.
        f_out.write((char*)&ddx_rho_u_tilde_a1[0], sizeof(double)*ddx_rho_u_tilde_a1.size());
        
        // Term III(1).
        f_out.write((char*)&b_ddx_p[0], sizeof(double)*b_ddx_p.size());
        // Term III(2).
        f_out.write((char*)&m_b_ddx_tau11[0], sizeof(double)*m_b_ddx_tau11.size());
        // Term III(3).
        f_out.write((char*)&m_R11_ddx_rho[0], sizeof(double)*m_R11_ddx_rho.size());
        
        // Term IV(1).
        f_out.write((char*)&rho_ddx_a1_sq[0], sizeof(double)*rho_ddx_a1_sq.size());
        // Term IV(2).
        f_out.write((char*)&m_rho_a1_ddx_u[0], sizeof(double)*m_rho_a1_ddx_u.size());
        
        // Term V.
        f_out.write((char*)&m_rho_ddx_rho_p_u_p_sq_over_rho[0], sizeof(double)*m_rho_ddx_rho_p_u_p_sq_over_rho.size());
        
        // Term VI(1).
        f_out.write((char*)&rho_rho_inv_p_ddx_p_p[0], sizeof(double)*rho_rho_inv_p_ddx_p_p.size());
        // Term VI(2).
        f_out.write((char*)&m_rho_rho_inv_p_ddx_tau_ij_p[0], sizeof(double)*m_rho_rho_inv_p_ddx_tau_ij_p.size());
        // Term VI(3).
        f_out.write((char*)&rho_epsilon_a1[0], sizeof(double)*rho_epsilon_a1.size());
        
        // Term II in moving frame of mixing layer.
        f_out.write((char*)&d_rho_a1_a1_dx[0], sizeof(double)*d_rho_a1_a1_dx.size());
        
        f_out.close();
    }
}


/*
 * Output budget of Reynolds normal stress in x-direction with inhomogeneous x-direction to a file.
 */
void
RTIRMIBudgetsUtilities::outputBudgetReynoldsNormalStressInXDirectionWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    std::ofstream f_out;
    
    MPIHelper MPI_helper = MPIHelper(
        "MPI_helper",
        d_dim,
        d_grid_geometry,
        patch_hierarchy);
    
    const std::vector<double>& dx_vec_coarsest = MPI_helper.getCoarsestDomainGridSpacing();
    const double dx_coarsest = dx_vec_coarsest[0];
    
    const std::vector<double>& dx_vec = MPI_helper.getFinestRefinedDomainGridSpacing();
    const double dx = dx_vec[0];
    
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
     * Output the spatial profiles (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        const std::vector<std::vector<double> >& rho_avg_coarsest_realizations   = d_ensemble_statistics->rho_avg_coarsest_realizations;
        const std::vector<std::vector<double> >& p_avg_coarsest_realizations     = d_ensemble_statistics->p_avg_coarsest_realizations;
        const std::vector<std::vector<double> >& u_avg_coarsest_realizations     = d_ensemble_statistics->u_avg_coarsest_realizations;
        const std::vector<std::vector<double> >& rho_u_avg_coarsest_realizations = d_ensemble_statistics->rho_u_avg_coarsest_realizations;
        const std::vector<std::vector<double> >& p_u_avg_coarsest_realizations   = d_ensemble_statistics->p_u_avg_coarsest_realizations;
        
        const std::vector<std::vector<double> >& rho_avg_realizations     = d_ensemble_statistics->rho_avg_realizations;
        const std::vector<std::vector<double> >& p_avg_realizations       = d_ensemble_statistics->p_avg_realizations;
        const std::vector<std::vector<double> >& u_avg_realizations       = d_ensemble_statistics->u_avg_realizations;
        const std::vector<std::vector<double> >& rho_u_avg_realizations   = d_ensemble_statistics->rho_u_avg_realizations;
        const std::vector<std::vector<double> >& rho_u_u_avg_realizations = d_ensemble_statistics->rho_u_u_avg_realizations;
        
        const std::vector<std::vector<double> >& p_u_avg_realizations = d_ensemble_statistics->p_u_avg_realizations;
        
        const std::vector<std::vector<double> >& ddx_rho_avg_realizations       = d_ensemble_statistics->ddx_rho_avg_realizations;
        const std::vector<std::vector<double> >& ddx_p_avg_realizations         = d_ensemble_statistics->ddx_p_avg_realizations;
        const std::vector<std::vector<double> >& ddx_u_avg_realizations         = d_ensemble_statistics->ddx_u_avg_realizations;
        const std::vector<std::vector<double> >& ddx_rho_u_avg_realizations     = d_ensemble_statistics->ddx_rho_u_avg_realizations;
        const std::vector<std::vector<double> >& ddx_rho_u_u_avg_realizations   = d_ensemble_statistics->ddx_rho_u_u_avg_realizations;
        const std::vector<std::vector<double> >& ddx_rho_u_u_u_avg_realizations = d_ensemble_statistics->ddx_rho_u_u_u_avg_realizations;
        const std::vector<std::vector<double> >& ddx_u_p_avg_realizations       = d_ensemble_statistics->ddx_u_p_avg_realizations;
        
        const std::vector<std::vector<double> >& ddy_u_avg_realizations = d_ensemble_statistics->ddy_u_avg_realizations;
        
        const std::vector<std::vector<double> >& ddz_u_avg_realizations = d_ensemble_statistics->ddz_u_avg_realizations;
        
        const std::vector<std::vector<double> >& p_ddx_u_avg_realizations = d_ensemble_statistics->p_ddx_u_avg_realizations;
        
        const std::vector<std::vector<double> >& tau11_avg_realizations = d_ensemble_statistics->tau11_avg_realizations;
        const std::vector<std::vector<double> >& tau12_avg_realizations = d_ensemble_statistics->tau12_avg_realizations;
        const std::vector<std::vector<double> >& tau13_avg_realizations = d_ensemble_statistics->tau13_avg_realizations;
        
        const std::vector<std::vector<double> >& u_tau11_avg_realizations = d_ensemble_statistics->u_tau11_avg_realizations;
        
        const std::vector<std::vector<double> >& ddx_tau11_avg_realizations = d_ensemble_statistics->ddx_tau11_avg_realizations;
        const std::vector<std::vector<double> >& ddy_tau12_avg_realizations = d_ensemble_statistics->ddy_tau12_avg_realizations;
        const std::vector<std::vector<double> >& ddz_tau13_avg_realizations = d_ensemble_statistics->ddz_tau13_avg_realizations;
        
        const std::vector<std::vector<double> >& tau11_ddx_u_avg_realizations = d_ensemble_statistics->tau11_ddx_u_avg_realizations;
        const std::vector<std::vector<double> >& tau12_ddy_u_avg_realizations = d_ensemble_statistics->tau12_ddy_u_avg_realizations;
        const std::vector<std::vector<double> >& tau13_ddz_u_avg_realizations = d_ensemble_statistics->tau13_ddz_u_avg_realizations;
        
        const int num_realizations = d_ensemble_statistics->getNumberOfEnsembles();
        
        TBOX_ASSERT(num_realizations > 0);
        TBOX_ASSERT(num_realizations == static_cast<int>(rho_avg_realizations.size()));
        
        const int num_cells_coarsest = static_cast<int>(p_avg_coarsest_realizations[0].size());
        
        const int num_cells = static_cast<int>(rho_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> rho_avg_coarsest_global(num_cells_coarsest, double(0));
        std::vector<double> p_avg_coarsest_global(num_cells_coarsest, double(0));
        std::vector<double> u_avg_coarsest_global(num_cells_coarsest, double(0));
        std::vector<double> rho_u_avg_coarsest_global(num_cells_coarsest, double(0));
        std::vector<double> p_u_avg_coarsest_global(num_cells_coarsest, double(0));
        
        std::vector<double> rho_avg_global(num_cells, double(0));
        std::vector<double> p_avg_global(num_cells, double(0));
        std::vector<double> u_avg_global(num_cells, double(0));
        std::vector<double> rho_u_avg_global(num_cells, double(0));
        std::vector<double> rho_u_u_avg_global(num_cells, double(0));
        
        std::vector<double> p_u_avg_global(num_cells, double(0));
        
        std::vector<double> ddx_rho_avg_global(num_cells, double(0));
        std::vector<double> ddx_p_avg_global(num_cells, double(0));
        std::vector<double> ddx_u_avg_global(num_cells, double(0));
        std::vector<double> ddx_rho_u_avg_global(num_cells, double(0));
        std::vector<double> ddx_rho_u_u_avg_global(num_cells, double(0));
        std::vector<double> ddx_rho_u_u_u_avg_global(num_cells, double(0));
        std::vector<double> ddx_u_p_avg_global(num_cells, double(0));
        
        std::vector<double> ddy_u_avg_global(num_cells, double(0));
        
        std::vector<double> ddz_u_avg_global(num_cells, double(0));
        
        std::vector<double> p_ddx_u_avg_global(num_cells, double(0));
        
        std::vector<double> tau11_avg_global(num_cells, double(0));
        std::vector<double> tau12_avg_global(num_cells, double(0));
        std::vector<double> tau13_avg_global(num_cells, double(0));
        
        std::vector<double> u_tau11_avg_global(num_cells, double(0));
        
        std::vector<double> ddx_tau11_avg_global(num_cells, double(0));
        std::vector<double> ddy_tau12_avg_global(num_cells, double(0));
        std::vector<double> ddz_tau13_avg_global(num_cells, double(0));
        
        std::vector<double> tau11_ddx_u_avg_global(num_cells, double(0));
        std::vector<double> tau12_ddy_u_avg_global(num_cells, double(0));
        std::vector<double> tau13_ddz_u_avg_global(num_cells, double(0));
        
        for (int ri = 0; ri < num_realizations; ri++)
        {
            for (int i = 0; i < num_cells_coarsest; i++)
            {
                rho_avg_coarsest_global[i]   += weight*rho_avg_coarsest_realizations[ri][i];
                p_avg_coarsest_global[i]     += weight*p_avg_coarsest_realizations[ri][i];
                u_avg_coarsest_global[i]     += weight*u_avg_coarsest_realizations[ri][i];
                rho_u_avg_coarsest_global[i] += weight*rho_u_avg_coarsest_realizations[ri][i];
                p_u_avg_coarsest_global[i]   += weight*p_u_avg_coarsest_realizations[ri][i];
            }
            
            for (int i = 0; i < num_cells; i++)
            {
                rho_avg_global[i]     += weight*rho_avg_realizations[ri][i];
                p_avg_global[i]       += weight*p_avg_realizations[ri][i];
                u_avg_global[i]       += weight*u_avg_realizations[ri][i];
                rho_u_avg_global[i]   += weight*rho_u_avg_realizations[ri][i];
                rho_u_u_avg_global[i] += weight*rho_u_u_avg_realizations[ri][i];
                
                p_u_avg_global[i] += weight*p_u_avg_realizations[ri][i];
                
                ddx_rho_avg_global[i]       += weight*ddx_rho_avg_realizations[ri][i];
                ddx_p_avg_global[i]         += weight*ddx_p_avg_realizations[ri][i];
                ddx_u_avg_global[i]         += weight*ddx_u_avg_realizations[ri][i];
                ddx_rho_u_avg_global[i]     += weight*ddx_rho_u_avg_realizations[ri][i];
                ddx_rho_u_u_avg_global[i]   += weight*ddx_rho_u_u_avg_realizations[ri][i];
                ddx_rho_u_u_u_avg_global[i] += weight*ddx_rho_u_u_u_avg_realizations[ri][i];
                ddx_u_p_avg_global[i]       += weight*ddx_u_p_avg_realizations[ri][i];
                
                p_ddx_u_avg_global[i] += weight*p_ddx_u_avg_realizations[ri][i];
                
                tau11_avg_global[i] += weight*tau11_avg_realizations[ri][i];
                
                u_tau11_avg_global[i] += weight*u_tau11_avg_realizations[ri][i];
                
                ddx_tau11_avg_global[i] += weight*ddx_tau11_avg_realizations[ri][i];
                
                tau11_ddx_u_avg_global[i] += weight*tau11_ddx_u_avg_realizations[ri][i];
                
            }
        }
        
        if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
        {
            for (int ri = 0; ri < num_realizations; ri++)
            {
                for (int i = 0; i < num_cells; i++)
                {
                    ddy_u_avg_global[i] += weight*ddy_u_avg_realizations[ri][i];
                    
                    tau12_avg_global[i] += weight*tau12_avg_realizations[ri][i];
                    
                    ddy_tau12_avg_global[i] += weight*ddy_tau12_avg_realizations[ri][i];
                    
                    tau12_ddy_u_avg_global[i] += weight*tau12_ddy_u_avg_realizations[ri][i];
                }
            }
        }
        
        if (d_dim == tbox::Dimension(3))
        {
            for (int ri = 0; ri < num_realizations; ri++)
            {
                for (int i = 0; i < num_cells; i++)
                {
                    ddz_u_avg_global[i] += weight*ddz_u_avg_realizations[ri][i];
                    
                    tau13_avg_global[i] += weight*tau13_avg_realizations[ri][i];
                    
                    ddz_tau13_avg_global[i] += weight*ddz_tau13_avg_realizations[ri][i];
                    
                    tau13_ddz_u_avg_global[i] += weight*tau13_ddz_u_avg_realizations[ri][i];
                }
            }
        }
        
        /*
         * Compute u_tilde.
         */
        
        std::vector<double> u_tilde(rho_u_avg_global);
        for (int i = 0; i < num_cells; i++)
        {
            u_tilde[i] /= rho_avg_global[i];
        }
        
        /*
         * Compute a1.
         */
        
        std::vector<double> rho_p_u_p(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            rho_p_u_p[i] = rho_u_avg_global[i] - rho_avg_global[i]*u_avg_global[i];
        }
        
        std::vector<double> a1(rho_p_u_p);
        for (int i = 0; i < num_cells; i++)
        {
            a1[i] /= rho_avg_global[i];
        }
        
        std::vector<double> rho_p_u_p_coarsest(num_cells_coarsest, double(0));
        for (int i = 0; i < num_cells_coarsest; i++)
        {
            rho_p_u_p_coarsest[i] = rho_u_avg_coarsest_global[i] - rho_avg_coarsest_global[i]*u_avg_coarsest_global[i];
        }
        
        std::vector<double> a1_coarsest(rho_p_u_p_coarsest);
        for (int i = 0; i < num_cells_coarsest; i++)
        {
            a1_coarsest[i] /= rho_avg_coarsest_global[i];
        }
        
        /*
         * Compute R11.
         */
        
        std::vector<double> rho_R11(num_cells, double(0));
        std::vector<double> R11(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            const double u_tilde       = rho_u_avg_global[i]/rho_avg_global[i];
            const double rho_u_pp_u_pp = rho_u_u_avg_global[i] - rho_u_avg_global[i]*u_tilde;
            
            rho_R11[i] = rho_u_pp_u_pp;
            R11[i]     = rho_u_pp_u_pp/rho_avg_global[i];
        }
        
        /*
         * Compute term II.
         */
        
        std::vector<double> ddx_rho_u_tilde_R11(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            const double ddx_R11_tilde = -(rho_R11[i]/(rho_avg_global[i]*rho_avg_global[i]))*ddx_rho_avg_global[i] +
                double(1)/rho_avg_global[i]*(ddx_rho_u_u_avg_global[i] - double(2)*u_tilde[i]*ddx_rho_u_avg_global[i] +
                u_tilde[i]*u_tilde[i]*ddx_rho_avg_global[i]);
            
            ddx_rho_u_tilde_R11[i] = rho_u_avg_global[i]*ddx_R11_tilde + R11[i]*ddx_rho_u_avg_global[i];
        }
        
        /*
         * Compute term II in moving frame of mixing layer.
         */
        
        std::vector<double> rho_a1_R11(rho_R11);
        for (int i = 0; i < num_cells; i++)
        {
            rho_a1_R11[i] *= a1[i];
        }
        
        std::vector<double> ddx_rho_a1_R11 = computeDerivativeOfVector1D(
            rho_a1_R11,
            dx);
        
        /*
         * Compute term III(1).
         */
        
        std::vector<double> two_a1_ddx_p(ddx_p_avg_global);
        for (int i = 0; i < num_cells; i++)
        {
            two_a1_ddx_p[i] *= (double(2)*a1[i]);
        }
        
        std::vector<double> ddx_p_coarsest = computeDerivativeOfVector1D(
            p_avg_coarsest_global,
            dx_coarsest);
        
        std::vector<double> two_a1_ddx_p_coarsest_refined(num_cells, double(0));
        for (int i = 0; i < num_cells_coarsest; i++)
        {
            for (int ii = 0; ii < ratio_finest_level_to_coarsest_level[0]; ii++)
            {
                const int idx_fine = i*ratio_finest_level_to_coarsest_level[0] + ii;
                
                two_a1_ddx_p_coarsest_refined[idx_fine] = double(2)*a1_coarsest[i]*ddx_p_coarsest[i];
            }
        }
        
        /*
         * Compute term III(2).
         */
        
        std::vector<double> m_2a1_ddx_tau11(ddx_tau11_avg_global);
        for (int i = 0; i < num_cells; i++)
        {
            m_2a1_ddx_tau11[i] *= (-double(2)*a1[i]);
        }
        
        /*
         * Compute term III(3).
         */
    
        std::vector<double> ddx_u_tilde(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            ddx_u_tilde[i] = ddx_rho_u_avg_global[i]/rho_avg_global[i] -
                rho_u_avg_global[i]/(rho_avg_global[i]*rho_avg_global[i])*ddx_rho_avg_global[i];
        }
        
        std::vector<double> m_2rho_R11_ddx_u_tilde(ddx_u_tilde);
        for (int i = 0; i < num_cells; i++)
        {
            m_2rho_R11_ddx_u_tilde[i] *= (-double(2)*rho_R11[i]);
        }
        
        /*
         * Compute term IV(1).
         */
        
        std::vector<double> m_ddx_rho_u_pp_u_pp_u_pp(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            m_ddx_rho_u_pp_u_pp_u_pp[i] = -(ddx_rho_u_u_u_avg_global[i] - double(2)*rho_u_u_avg_global[i]*ddx_u_tilde[i] -
                double(2)*u_tilde[i]*ddx_rho_u_u_avg_global[i] + u_tilde[i]*u_tilde[i]*ddx_rho_u_avg_global[i] +
                double(2)*rho_u_avg_global[i]*u_tilde[i]*ddx_u_tilde[i] - ddx_rho_u_tilde_R11[i]);
        }
        
        /*
         * Compute term IV(2).
         */
        
        std::vector<double> m_2ddx_u_p_p_p(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            m_2ddx_u_p_p_p[i] = -double(2)*(ddx_u_p_avg_global[i] - u_avg_global[i]*ddx_p_avg_global[i] -
                p_avg_global[i]*ddx_u_avg_global[i]);
        }
        
        std::vector<double> u_p_p_p_coarsest(num_cells_coarsest, double(0));
        for (int i = 0; i < num_cells_coarsest; i++)
        {
            u_p_p_p_coarsest[i] = p_u_avg_coarsest_global[i] - u_avg_coarsest_global[i]*p_avg_coarsest_global[i];
        }
        
        std::vector<double> m_2ddx_u_p_p_p_coarsest = computeDerivativeOfVector1D(
            u_p_p_p_coarsest,
            dx_coarsest);
        
        for (int i = 0; i < num_cells_coarsest; i++)
        {
            m_2ddx_u_p_p_p_coarsest[i] *= (-double(2));
        }
        
        std::vector<double> m_2ddx_u_p_p_p_coarsest_refined(num_cells, double(0));
        for (int i = 0; i < num_cells_coarsest; i++)
        {
            for (int ii = 0; ii < ratio_finest_level_to_coarsest_level[0]; ii++)
            {
                const int idx_fine = i*ratio_finest_level_to_coarsest_level[0] + ii;
                
                m_2ddx_u_p_p_p_coarsest_refined[idx_fine] = m_2ddx_u_p_p_p_coarsest[i];
            }
        }
        
        /*
         * Compute term IV(3).
         */
        
        std::vector<double> u_p_tau11_p(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            u_p_tau11_p[i] = u_tau11_avg_global[i] - u_avg_global[i]*tau11_avg_global[i];
        }
        
        std::vector<double> two_ddx_u_p_tau11_p = computeDerivativeOfVector1D(
            u_p_tau11_p,
            dx);
        
        for (int i = 0; i < num_cells; i++)
        {
            two_ddx_u_p_tau11_p[i] *= double(2);
        }
        
        /*
         * Compute term V.
         */
        
        std::vector<double> two_p_p_ddx_u_p(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            two_p_p_ddx_u_p[i] = double(2)*(p_ddx_u_avg_global[i] - p_avg_global[i]*ddx_u_avg_global[i]);
        }
        
        /*
         * Compute term VI.
         */
        
        std::vector<double> tau11_p_ddx_u_p(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            tau11_p_ddx_u_p[i] = tau11_ddx_u_avg_global[i] - tau11_avg_global[i]*ddx_u_avg_global[i];
        }
        
        std::vector<double> tau12_p_ddy_u_p(num_cells, double(0));
        if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
        {
            for (int i = 0; i < num_cells; i++)
            {
                tau12_p_ddy_u_p[i] = tau12_ddy_u_avg_global[i] - tau12_avg_global[i]*ddy_u_avg_global[i];
            }
        }
        
        std::vector<double> tau13_p_ddz_u_p(num_cells, double(0));
        if (d_dim == tbox::Dimension(3))
        {
            for (int i = 0; i < num_cells; i++)
            {
                tau13_p_ddz_u_p[i] = tau13_ddz_u_avg_global[i] - tau13_avg_global[i]*ddz_u_avg_global[i];
            }
        }
        
        std::vector<double> m_2tau1i_p_ddxi_u_p(num_cells, double(0));
        if (d_dim == tbox::Dimension(1))
        {
            for (int i = 0; i < num_cells; i++)
            {
                m_2tau1i_p_ddxi_u_p[i] = -double(2)*tau11_p_ddx_u_p[i];
            }
        }
        else if (d_dim == tbox::Dimension(2))
        {
            for (int i = 0; i < num_cells; i++)
            {
                m_2tau1i_p_ddxi_u_p[i] = -double(2)*(tau11_p_ddx_u_p[i] + tau12_p_ddy_u_p[i]);
            }
        }
        else if (d_dim == tbox::Dimension(3))
        {
            for (int i = 0; i < num_cells; i++)
            {
                m_2tau1i_p_ddxi_u_p[i] = -double(2)*(tau11_p_ddx_u_p[i] + tau12_p_ddy_u_p[i] + tau13_p_ddz_u_p[i]);
            }
        }
        
        /*
         * Output budget.
         */
        
        f_out.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        f_out.write((char*)&output_time, sizeof(double));
        f_out.write((char*)&rho_R11[0], sizeof(double)*rho_R11.size());
        // Term II.
        f_out.write((char*)&ddx_rho_u_tilde_R11[0], sizeof(double)*ddx_rho_u_tilde_R11.size());
        
        // Term III(1).
        f_out.write((char*)&two_a1_ddx_p[0], sizeof(double)*two_a1_ddx_p.size());
        // Term III(2).
        f_out.write((char*)&m_2a1_ddx_tau11[0], sizeof(double)*m_2a1_ddx_tau11.size());
        // Term III(3).
        f_out.write((char*)&m_2rho_R11_ddx_u_tilde[0], sizeof(double)*m_2rho_R11_ddx_u_tilde.size());
        
        // Term IV(1).
        f_out.write((char*)&m_ddx_rho_u_pp_u_pp_u_pp[0], sizeof(double)*m_ddx_rho_u_pp_u_pp_u_pp.size());
        // Term IV(2).
        f_out.write((char*)&m_2ddx_u_p_p_p[0], sizeof(double)*m_2ddx_u_p_p_p.size());
        // Term IV(3).
        f_out.write((char*)&two_ddx_u_p_tau11_p[0], sizeof(double)*two_ddx_u_p_tau11_p.size());
        
        // Term V.
        f_out.write((char*)&two_p_p_ddx_u_p[0], sizeof(double)*two_p_p_ddx_u_p.size());
        
        // Term VI.
        f_out.write((char*)&m_2tau1i_p_ddxi_u_p[0], sizeof(double)*m_2tau1i_p_ddxi_u_p.size());
        
        // Term II in moving frame of mixing layer.
        f_out.write((char*)&ddx_rho_a1_R11[0], sizeof(double)*ddx_rho_a1_R11.size());
        
        // Term III(1).
        f_out.write((char*)&two_a1_ddx_p_coarsest_refined[0], sizeof(double)*two_a1_ddx_p_coarsest_refined.size());
        
        // Term IV(2) on the coarsest level.
        f_out.write((char*)&m_2ddx_u_p_p_p_coarsest_refined[0], sizeof(double)*m_2ddx_u_p_p_p_coarsest_refined.size());
        
        f_out.close();
    }
}


/*
 * Output budget of Reynolds normal stress in y-direction with inhomogeneous x-direction to a file.
 */
void
RTIRMIBudgetsUtilities::outputBudgetReynoldsNormalStressInYDirectionWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The budget of Reynolds normal stress in y-direction cannot be outputted for 1D problem!"
            << std::endl);
    }
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    std::ofstream f_out;
    
    MPIHelper MPI_helper = MPIHelper(
        "MPI_helper",
        d_dim,
        d_grid_geometry,
        patch_hierarchy);
    
    const std::vector<double>& dx_vec = MPI_helper.getFinestRefinedDomainGridSpacing();
    const double dx = dx_vec[0];
    
    /*
     * Output the spatial profiles (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        const std::vector<std::vector<double> >& rho_avg_realizations     = d_ensemble_statistics->rho_avg_realizations;
        const std::vector<std::vector<double> >& p_avg_realizations       = d_ensemble_statistics->p_avg_realizations;
        const std::vector<std::vector<double> >& u_avg_realizations       = d_ensemble_statistics->u_avg_realizations;
        const std::vector<std::vector<double> >& v_avg_realizations       = d_ensemble_statistics->v_avg_realizations;
        const std::vector<std::vector<double> >& rho_u_avg_realizations   = d_ensemble_statistics->rho_u_avg_realizations;
        const std::vector<std::vector<double> >& rho_v_avg_realizations   = d_ensemble_statistics->rho_v_avg_realizations;
        const std::vector<std::vector<double> >& rho_v_v_avg_realizations = d_ensemble_statistics->rho_v_v_avg_realizations;
        
        const std::vector<std::vector<double> >& rho_u_v_avg_realizations = d_ensemble_statistics->rho_u_v_avg_realizations;
        
        const std::vector<std::vector<double> >& ddx_rho_avg_realizations     = d_ensemble_statistics->ddx_rho_avg_realizations;
        const std::vector<std::vector<double> >& ddx_v_avg_realizations       = d_ensemble_statistics->ddx_v_avg_realizations;
        const std::vector<std::vector<double> >& ddx_rho_u_avg_realizations   = d_ensemble_statistics->ddx_rho_u_avg_realizations;
        const std::vector<std::vector<double> >& ddx_rho_v_avg_realizations   = d_ensemble_statistics->ddx_rho_v_avg_realizations;
        const std::vector<std::vector<double> >& ddx_rho_v_v_avg_realizations = d_ensemble_statistics->ddx_rho_v_v_avg_realizations;
        
        const std::vector<std::vector<double> >& ddx_rho_u_v_avg_realizations   = d_ensemble_statistics->ddx_rho_u_v_avg_realizations;
        const std::vector<std::vector<double> >& ddx_rho_u_v_v_avg_realizations = d_ensemble_statistics->ddx_rho_u_v_v_avg_realizations;
        
        const std::vector<std::vector<double> >& ddy_v_avg_realizations = d_ensemble_statistics->ddy_v_avg_realizations;
        
        const std::vector<std::vector<double> >& ddz_v_avg_realizations = d_ensemble_statistics->ddz_v_avg_realizations;
        
        const std::vector<std::vector<double> >& p_ddy_v_avg_realizations = d_ensemble_statistics->p_ddy_v_avg_realizations;
        
        const std::vector<std::vector<double> >& tau12_avg_realizations = d_ensemble_statistics->tau12_avg_realizations;
        const std::vector<std::vector<double> >& tau22_avg_realizations = d_ensemble_statistics->tau22_avg_realizations;
        const std::vector<std::vector<double> >& tau23_avg_realizations = d_ensemble_statistics->tau23_avg_realizations;
        
        const std::vector<std::vector<double> >& v_tau12_avg_realizations = d_ensemble_statistics->v_tau12_avg_realizations;
        
        const std::vector<std::vector<double> >& tau12_ddx_v_avg_realizations = d_ensemble_statistics->tau12_ddx_v_avg_realizations;
        const std::vector<std::vector<double> >& tau22_ddy_v_avg_realizations = d_ensemble_statistics->tau22_ddy_v_avg_realizations;
        const std::vector<std::vector<double> >& tau23_ddz_v_avg_realizations = d_ensemble_statistics->tau23_ddz_v_avg_realizations;
        
        const int num_realizations = d_ensemble_statistics->getNumberOfEnsembles();
        
        TBOX_ASSERT(num_realizations > 0);
        TBOX_ASSERT(num_realizations == static_cast<int>(rho_avg_realizations.size()));
        
        const int num_cells = static_cast<int>(rho_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> rho_avg_global(num_cells, double(0));
        std::vector<double> p_avg_global(num_cells, double(0));
        std::vector<double> u_avg_global(num_cells, double(0));
        std::vector<double> v_avg_global(num_cells, double(0));
        std::vector<double> rho_u_avg_global(num_cells, double(0));
        std::vector<double> rho_v_avg_global(num_cells, double(0));
        std::vector<double> rho_v_v_avg_global(num_cells, double(0));
        
        std::vector<double> rho_u_v_avg_global(num_cells, double(0));
        
        std::vector<double> ddx_rho_avg_global(num_cells, double(0));
        std::vector<double> ddx_v_avg_global(num_cells, double(0));
        std::vector<double> ddx_rho_u_avg_global(num_cells, double(0));
        std::vector<double> ddx_rho_v_avg_global(num_cells, double(0));
        std::vector<double> ddx_rho_v_v_avg_global(num_cells, double(0));
        
        std::vector<double> ddx_rho_u_v_avg_global(num_cells, double(0));
        std::vector<double> ddx_rho_u_v_v_avg_global(num_cells, double(0));
        
        std::vector<double> ddy_v_avg_global(num_cells, double(0));
        
        std::vector<double> ddz_v_avg_global(num_cells, double(0));
        
        std::vector<double> p_ddy_v_avg_global(num_cells, double(0));
        
        std::vector<double> tau12_avg_global(num_cells, double(0));
        std::vector<double> tau22_avg_global(num_cells, double(0));
        std::vector<double> tau23_avg_global(num_cells, double(0));
        
        std::vector<double> v_tau12_avg_global(num_cells, double(0));
        
        std::vector<double> tau12_ddx_v_avg_global(num_cells, double(0));
        std::vector<double> tau22_ddy_v_avg_global(num_cells, double(0));
        std::vector<double> tau23_ddz_v_avg_global(num_cells, double(0));
        
        if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
        {
            for (int ri = 0; ri < num_realizations; ri++)
            {
                for (int i = 0; i < num_cells; i++)
                {
                    rho_avg_global[i]     += weight*rho_avg_realizations[ri][i];
                    p_avg_global[i]       += weight*p_avg_realizations[ri][i];
                    u_avg_global[i]       += weight*u_avg_realizations[ri][i];
                    v_avg_global[i]       += weight*v_avg_realizations[ri][i];
                    rho_u_avg_global[i]   += weight*rho_u_avg_realizations[ri][i];
                    rho_v_avg_global[i]   += weight*rho_v_avg_realizations[ri][i];
                    rho_v_v_avg_global[i] += weight*rho_v_v_avg_realizations[ri][i];
                    
                    rho_u_v_avg_global[i] += weight*rho_u_v_avg_realizations[ri][i];
                    
                    ddx_rho_avg_global[i]     += weight*ddx_rho_avg_realizations[ri][i];
                    ddx_v_avg_global[i]       += weight*ddx_v_avg_realizations[ri][i];
                    ddx_rho_u_avg_global[i]   += weight*ddx_rho_u_avg_realizations[ri][i];
                    ddx_rho_v_avg_global[i]   += weight*ddx_rho_v_avg_realizations[ri][i];
                    ddx_rho_v_v_avg_global[i] += weight*ddx_rho_v_v_avg_realizations[ri][i];
                    
                    ddx_rho_u_v_avg_global[i]   += weight*ddx_rho_u_v_avg_realizations[ri][i];
                    ddx_rho_u_v_v_avg_global[i] += weight*ddx_rho_u_v_v_avg_realizations[ri][i];
                    
                    ddy_v_avg_global[i] += weight*ddy_v_avg_realizations[ri][i];
                    
                    p_ddy_v_avg_global[i] += weight*p_ddy_v_avg_realizations[ri][i];
                    
                    tau12_avg_global[i] += weight*tau12_avg_realizations[ri][i];
                    tau22_avg_global[i] += weight*tau22_avg_realizations[ri][i];
                    
                    v_tau12_avg_global[i] += weight*v_tau12_avg_realizations[ri][i];
                    
                    tau12_ddx_v_avg_global[i] += weight*tau12_ddx_v_avg_realizations[ri][i];
                    tau22_ddy_v_avg_global[i] += weight*tau22_ddy_v_avg_realizations[ri][i];
                }
            }
        }
        
        if (d_dim == tbox::Dimension(3))
        {
            for (int ri = 0; ri < num_realizations; ri++)
            {
                for (int i = 0; i < num_cells; i++)
                {
                    ddz_v_avg_global[i] += weight*ddz_v_avg_realizations[ri][i];
                    
                    tau23_avg_global[i] += weight*tau23_avg_realizations[ri][i];
                    
                    tau23_ddz_v_avg_global[i] += weight*tau22_ddy_v_avg_realizations[ri][i];
                }
            }
        }
        
        /*
         * Compute u_tilde.
         */
        
        std::vector<double> u_tilde(rho_u_avg_global);
        for (int i = 0; i < num_cells; i++)
        {
            u_tilde[i] /= rho_avg_global[i];
        }
        
        /*
         * Compute v_tilde.
         */
        
        std::vector<double> v_tilde(rho_v_avg_global);
        for (int i = 0; i < num_cells; i++)
        {
            v_tilde[i] /= rho_avg_global[i];
        }
        
        /*
         * Compute a1.
         */
        
        std::vector<double> rho_p_u_p(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            rho_p_u_p[i] = rho_u_avg_global[i] - rho_avg_global[i]*u_avg_global[i];
        }
        
        std::vector<double> a1(rho_p_u_p);
        for (int i = 0; i < num_cells; i++)
        {
            a1[i] /= rho_avg_global[i];
        }
        
        /*
         * Compute R22.
         */
        
        std::vector<double> rho_R22(num_cells, double(0));
        std::vector<double> R22(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            const double v_tilde       = rho_v_avg_global[i]/rho_avg_global[i];
            const double rho_v_pp_v_pp = rho_v_v_avg_global[i] - rho_v_avg_global[i]*v_tilde;
            
            rho_R22[i] = rho_v_pp_v_pp;
            R22[i]     = rho_v_pp_v_pp/rho_avg_global[i];
        }
        
        /*
         * Compute term II.
         */
        
        std::vector<double> ddx_rho_u_tilde_R22(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            const double ddx_R22_tilde = -(rho_R22[i]/(rho_avg_global[i]*rho_avg_global[i]))*ddx_rho_avg_global[i] +
                double(1)/rho_avg_global[i]*(ddx_rho_v_v_avg_global[i] - double(2)*v_tilde[i]*ddx_rho_v_avg_global[i] +
                v_tilde[i]*v_tilde[i]*ddx_rho_avg_global[i]);
            
            ddx_rho_u_tilde_R22[i] = rho_u_avg_global[i]*ddx_R22_tilde + R22[i]*ddx_rho_u_avg_global[i];
        }
        
        /*
         * Compute term II in moving frame of mixing layer.
         */
        
        std::vector<double> rho_a1_R22(rho_R22);
        for (int i = 0; i < num_cells; i++)
        {
            rho_a1_R22[i] *= a1[i];
        }
        
        std::vector<double> ddx_rho_a1_R22 = computeDerivativeOfVector1D(
            rho_a1_R22,
            dx);
        
        /*
         * Compute term IV(1).
         */
        
        std::vector<double> ddx_v_tilde(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            ddx_v_tilde[i] = ddx_rho_v_avg_global[i]/rho_avg_global[i] - rho_v_avg_global[i]/
                (rho_avg_global[i]*rho_avg_global[i])*ddx_rho_avg_global[i];
        }
        
        std::vector<double> m_ddx_rho_v_pp_v_pp_u_pp(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            m_ddx_rho_v_pp_v_pp_u_pp[i] = -(ddx_rho_u_v_v_avg_global[i] - double(2)*rho_u_v_avg_global[i]*ddx_v_tilde[i] -
                double(2)*v_tilde[i]*ddx_rho_u_v_avg_global[i] + v_tilde[i]*v_tilde[i]*ddx_rho_u_avg_global[i] +
                double(2)*rho_u_avg_global[i]*v_tilde[i]*ddx_v_tilde[i] - ddx_rho_u_tilde_R22[i]);
        }
        
        /*
         * Compute term IV(2).
         */
        
        std::vector<double> v_p_tau12_p(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            v_p_tau12_p[i] = v_tau12_avg_global[i] - v_avg_global[i]*tau12_avg_global[i];
        }
        
        std::vector<double> two_ddx_v_p_tau12_p = computeDerivativeOfVector1D(
            v_p_tau12_p,
            dx);
        
        for (int i = 0; i < num_cells; i++)
        {
            two_ddx_v_p_tau12_p[i] *= double(2);
        }
        
        /*
         * Compute term V.
         */
        
        std::vector<double> two_p_p_ddy_v_p(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            two_p_p_ddy_v_p[i] = double(2)*(p_ddy_v_avg_global[i] - p_avg_global[i]*ddy_v_avg_global[i]);
        }
        
        /*
         * Compute term VI.
         */
        
        std::vector<double> tau12_p_ddx_v_p(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            tau12_p_ddx_v_p[i] = tau12_ddx_v_avg_global[i] - tau12_avg_global[i]*ddx_v_avg_global[i];
        }
        
        std::vector<double> tau22_p_ddy_v_p(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            tau22_p_ddy_v_p[i] = tau22_ddy_v_avg_global[i] - tau22_avg_global[i]*ddy_v_avg_global[i];
        }
        
        std::vector<double> tau23_p_ddz_v_p(num_cells, double(0));
        if (d_dim == tbox::Dimension(3))
        {
            for (int i = 0; i < num_cells; i++)
            {
                tau23_p_ddz_v_p[i] = tau23_ddz_v_avg_global[i] - tau23_avg_global[i]*ddz_v_avg_global[i];
            }
        }
        
        std::vector<double> m_2tau_2i_p_ddxi_v_p(num_cells, double(0));
        if (d_dim == tbox::Dimension(2))
        {
            for (int i = 0; i < num_cells; i++)
            {
                m_2tau_2i_p_ddxi_v_p[i] = -double(2)*(tau12_p_ddx_v_p[i] + tau22_p_ddy_v_p[i]);
            }
        }
        else if (d_dim == tbox::Dimension(3))
        {
            for (int i = 0; i < num_cells; i++)
            {
                m_2tau_2i_p_ddxi_v_p[i] = -double(2)*(tau12_p_ddx_v_p[i] + tau22_p_ddy_v_p[i] + tau23_p_ddz_v_p[i]);
            }
        }
        
        /*
         * Output budget.
         */
        
        f_out.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        f_out.write((char*)&output_time, sizeof(double));
        f_out.write((char*)&rho_R22[0], sizeof(double)*rho_R22.size());
        // Term II.
        f_out.write((char*)&ddx_rho_u_tilde_R22[0], sizeof(double)*ddx_rho_u_tilde_R22.size());
        
        // Term IV(1).
        f_out.write((char*)&m_ddx_rho_v_pp_v_pp_u_pp[0], sizeof(double)*m_ddx_rho_v_pp_v_pp_u_pp.size());
        // Term IV(2).
        f_out.write((char*)&two_ddx_v_p_tau12_p[0], sizeof(double)*two_ddx_v_p_tau12_p.size());
        
        // Term V.
        f_out.write((char*)&two_p_p_ddy_v_p[0], sizeof(double)*two_p_p_ddy_v_p.size());
        
        // Term VI.
        f_out.write((char*)&m_2tau_2i_p_ddxi_v_p[0], sizeof(double)*m_2tau_2i_p_ddxi_v_p.size());
        
        // Term II in moving frame of mixing layer.
        f_out.write((char*)&ddx_rho_a1_R22[0], sizeof(double)*ddx_rho_a1_R22.size());
        
        f_out.close();
    }
}


/*
 * Output budget of Reynolds normal stress in z-direction with inhomogeneous x-direction to a file.
 */
void
RTIRMIBudgetsUtilities::outputBudgetReynoldsNormalStressInZDirectionWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_dim == tbox::Dimension(1) || d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The budget of Reynolds normal stress in z-direction cannot be outputted for 1D or 2D problem!"
            << std::endl);
    }
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    std::ofstream f_out;
    
    MPIHelper MPI_helper = MPIHelper(
        "MPI_helper",
        d_dim,
        d_grid_geometry,
        patch_hierarchy);
    
    const std::vector<double>& dx_vec = MPI_helper.getFinestRefinedDomainGridSpacing();
    const double dx = dx_vec[0];
    
    /*
     * Output the spatial profiles (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        const std::vector<std::vector<double> >& rho_avg_realizations     = d_ensemble_statistics->rho_avg_realizations;
        const std::vector<std::vector<double> >& p_avg_realizations       = d_ensemble_statistics->p_avg_realizations;
        const std::vector<std::vector<double> >& u_avg_realizations       = d_ensemble_statistics->u_avg_realizations;
        const std::vector<std::vector<double> >& w_avg_realizations       = d_ensemble_statistics->w_avg_realizations;
        const std::vector<std::vector<double> >& rho_u_avg_realizations   = d_ensemble_statistics->rho_u_avg_realizations;
        const std::vector<std::vector<double> >& rho_w_avg_realizations   = d_ensemble_statistics->rho_w_avg_realizations;
        const std::vector<std::vector<double> >& rho_w_w_avg_realizations = d_ensemble_statistics->rho_w_w_avg_realizations;
        
        const std::vector<std::vector<double> >& rho_u_w_avg_realizations = d_ensemble_statistics->rho_u_w_avg_realizations;
        
        const std::vector<std::vector<double> >& ddx_rho_avg_realizations     = d_ensemble_statistics->ddx_rho_avg_realizations;
        const std::vector<std::vector<double> >& ddx_w_avg_realizations       = d_ensemble_statistics->ddx_w_avg_realizations;
        const std::vector<std::vector<double> >& ddx_rho_u_avg_realizations   = d_ensemble_statistics->ddx_rho_u_avg_realizations;
        const std::vector<std::vector<double> >& ddx_rho_w_avg_realizations   = d_ensemble_statistics->ddx_rho_w_avg_realizations;
        const std::vector<std::vector<double> >& ddx_rho_w_w_avg_realizations = d_ensemble_statistics->ddx_rho_w_w_avg_realizations;
        
        const std::vector<std::vector<double> >& ddx_rho_u_w_avg_realizations   = d_ensemble_statistics->ddx_rho_u_w_avg_realizations;
        const std::vector<std::vector<double> >& ddx_rho_u_w_w_avg_realizations = d_ensemble_statistics->ddx_rho_u_w_w_avg_realizations;
        
        const std::vector<std::vector<double> >& ddy_w_avg_realizations = d_ensemble_statistics->ddy_w_avg_realizations;
        
        const std::vector<std::vector<double> >& ddz_w_avg_realizations = d_ensemble_statistics->ddz_w_avg_realizations;
        
        const std::vector<std::vector<double> >& p_ddz_w_avg_realizations = d_ensemble_statistics->p_ddz_w_avg_realizations;
        
        const std::vector<std::vector<double> >& tau13_avg_realizations = d_ensemble_statistics->tau13_avg_realizations;
        const std::vector<std::vector<double> >& tau23_avg_realizations = d_ensemble_statistics->tau23_avg_realizations;
        const std::vector<std::vector<double> >& tau33_avg_realizations = d_ensemble_statistics->tau33_avg_realizations;
        
        const std::vector<std::vector<double> >& w_tau13_avg_realizations = d_ensemble_statistics->w_tau13_avg_realizations;
        
        const std::vector<std::vector<double> >& tau13_ddx_w_avg_realizations = d_ensemble_statistics->tau13_ddx_w_avg_realizations;
        const std::vector<std::vector<double> >& tau23_ddy_w_avg_realizations = d_ensemble_statistics->tau23_ddy_w_avg_realizations;
        const std::vector<std::vector<double> >& tau33_ddz_w_avg_realizations = d_ensemble_statistics->tau33_ddz_w_avg_realizations;
        
        const int num_realizations = d_ensemble_statistics->getNumberOfEnsembles();
        
        TBOX_ASSERT(num_realizations > 0);
        TBOX_ASSERT(num_realizations == static_cast<int>(rho_avg_realizations.size()));
        
        const int num_cells = static_cast<int>(rho_avg_realizations[0].size());
        const double weight = double(1)/double(num_realizations);
        
        std::vector<double> rho_avg_global(num_cells, double(0));
        std::vector<double> p_avg_global(num_cells, double(0));
        std::vector<double> u_avg_global(num_cells, double(0));
        std::vector<double> w_avg_global(num_cells, double(0));
        std::vector<double> rho_u_avg_global(num_cells, double(0));
        std::vector<double> rho_w_avg_global(num_cells, double(0));
        std::vector<double> rho_w_w_avg_global(num_cells, double(0));
        
        std::vector<double> rho_u_w_avg_global(num_cells, double(0));
        
        std::vector<double> ddx_rho_avg_global(num_cells, double(0));
        std::vector<double> ddx_w_avg_global(num_cells, double(0));
        std::vector<double> ddx_rho_u_avg_global(num_cells, double(0));
        std::vector<double> ddx_rho_w_avg_global(num_cells, double(0));
        std::vector<double> ddx_rho_w_w_avg_global(num_cells, double(0));
        
        std::vector<double> ddx_rho_u_w_avg_global(num_cells, double(0));
        std::vector<double> ddx_rho_u_w_w_avg_global(num_cells, double(0));
        
        std::vector<double> ddy_w_avg_global(num_cells, double(0));
        
        std::vector<double> ddz_w_avg_global(num_cells, double(0));
        
        std::vector<double> p_ddz_w_avg_global(num_cells, double(0));
        
        std::vector<double> tau13_avg_global(num_cells, double(0));
        std::vector<double> tau23_avg_global(num_cells, double(0));
        std::vector<double> tau33_avg_global(num_cells, double(0));
        
        std::vector<double> w_tau13_avg_global(num_cells, double(0));
        
        std::vector<double> tau13_ddx_w_avg_global(num_cells, double(0));
        std::vector<double> tau23_ddy_w_avg_global(num_cells, double(0));
        std::vector<double> tau33_ddz_w_avg_global(num_cells, double(0));
        
        if (d_dim == tbox::Dimension(3))
        {
            for (int ri = 0; ri < num_realizations; ri++)
            {
                for (int i = 0; i < num_cells; i++)
                {
                    rho_avg_global[i]     += weight*rho_avg_realizations[ri][i];
                    p_avg_global[i]       += weight*p_avg_realizations[ri][i];
                    u_avg_global[i]       += weight*u_avg_realizations[ri][i];
                    w_avg_global[i]       += weight*w_avg_realizations[ri][i];
                    rho_u_avg_global[i]   += weight*rho_u_avg_realizations[ri][i];
                    rho_w_avg_global[i]   += weight*rho_w_avg_realizations[ri][i];
                    rho_w_w_avg_global[i] += weight*rho_w_w_avg_realizations[ri][i];
                    
                    rho_u_w_avg_global[i] += weight*rho_u_w_avg_realizations[ri][i];
                    
                    ddx_rho_avg_global[i]     += weight*ddx_rho_avg_realizations[ri][i];
                    ddx_w_avg_global[i]       += weight*ddx_w_avg_realizations[ri][i];
                    ddx_rho_u_avg_global[i]   += weight*ddx_rho_u_avg_realizations[ri][i];
                    ddx_rho_w_avg_global[i]   += weight*ddx_rho_w_avg_realizations[ri][i];
                    ddx_rho_w_w_avg_global[i] += weight*ddx_rho_w_w_avg_realizations[ri][i];
                    
                    ddx_rho_u_w_avg_global[i]   += weight*ddx_rho_u_w_avg_realizations[ri][i];
                    ddx_rho_u_w_w_avg_global[i] += weight*ddx_rho_u_w_w_avg_realizations[ri][i];
                    
                    ddy_w_avg_global[i] += weight*ddy_w_avg_realizations[ri][i];
                    
                    ddz_w_avg_global[i] += weight*ddz_w_avg_realizations[ri][i];
                    
                    p_ddz_w_avg_global[i] += weight*p_ddz_w_avg_realizations[ri][i];
                    
                    tau13_avg_global[i] += weight*tau13_avg_realizations[ri][i];
                    tau23_avg_global[i] += weight*tau23_avg_realizations[ri][i];
                    tau33_avg_global[i] += weight*tau33_avg_realizations[ri][i];
                    
                    w_tau13_avg_global[i] += weight*w_tau13_avg_realizations[ri][i];
                    
                    tau13_ddx_w_avg_global[i] += weight*tau13_ddx_w_avg_realizations[ri][i];
                    tau23_ddy_w_avg_global[i] += weight*tau23_ddy_w_avg_realizations[ri][i];
                    tau33_ddz_w_avg_global[i] += weight*tau33_ddz_w_avg_realizations[ri][i];
                }
            }
        }
        
        /*
         * Compute u_tilde.
         */
        
        std::vector<double> u_tilde(rho_u_avg_global);
        for (int i = 0; i < num_cells; i++)
        {
            u_tilde[i] /= rho_avg_global[i];
        }
        
        /*
         * Compute w_tilde.
         */
        
        std::vector<double> w_tilde(rho_w_avg_global);
        for (int i = 0; i < num_cells; i++)
        {
            w_tilde[i] /= rho_avg_global[i];
        }
        
        /*
         * Compute a1.
         */
        
        std::vector<double> rho_p_u_p(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            rho_p_u_p[i] = rho_u_avg_global[i] - rho_avg_global[i]*u_avg_global[i];
        }
        
        std::vector<double> a1(rho_p_u_p);
        for (int i = 0; i < num_cells; i++)
        {
            a1[i] /= rho_avg_global[i];
        }
        
        /*
         * Compute R33.
         */
        
        std::vector<double> rho_R33(num_cells, double(0));
        std::vector<double> R33(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            const double w_tilde       = rho_w_avg_global[i]/rho_avg_global[i];
            const double rho_w_pp_w_pp = rho_w_w_avg_global[i] - rho_w_avg_global[i]*w_tilde;
            
            rho_R33[i] = rho_w_pp_w_pp;
            R33[i]     = rho_w_pp_w_pp/rho_avg_global[i];
        }
        
        /*
         * Compute term II.
         */
        
        std::vector<double> ddx_rho_u_tilde_R33(num_cells, double(0));
        
        for (int i = 0; i < num_cells; i++)
        {
            const double ddx_R33_tilde = -(rho_R33[i]/(rho_avg_global[i]*rho_avg_global[i]))*ddx_rho_avg_global[i] + 
                double(1)/rho_avg_global[i]*(ddx_rho_w_w_avg_global[i] - double(2)*w_tilde[i]*ddx_rho_w_avg_global[i] + 
                w_tilde[i]*w_tilde[i]*ddx_rho_avg_global[i]);
            
            ddx_rho_u_tilde_R33[i] = rho_u_avg_global[i]*ddx_R33_tilde + R33[i]*ddx_rho_u_avg_global[i];
        }
        
        /*
         * Compute term II in moving frame of mixing layer.
         */
        
        std::vector<double> rho_a1_R33(rho_R33);
        for (int i = 0; i < num_cells; i++)
        {
            rho_a1_R33[i] *= a1[i];
        }
        
        std::vector<double> ddx_rho_a1_R33 = computeDerivativeOfVector1D(
            rho_a1_R33,
            dx);
        
        /*
         * Compute term IV(1).
         */
        
        std::vector<double> ddx_w_tilde(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            ddx_w_tilde[i] = ddx_rho_w_avg_global[i]/rho_avg_global[i] - rho_w_avg_global[i]/
                (rho_avg_global[i]*rho_avg_global[i])*ddx_rho_avg_global[i];
        }
        
        std::vector<double> m_ddx_rho_w_pp_w_pp_u_pp(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            m_ddx_rho_w_pp_w_pp_u_pp[i] = -(ddx_rho_u_w_w_avg_global[i] - double(2)*rho_u_w_avg_global[i]*ddx_w_tilde[i] -
                double(2)*w_tilde[i]*ddx_rho_u_w_avg_global[i] + w_tilde[i]*w_tilde[i]*ddx_rho_u_avg_global[i] +
                double(2)*rho_u_avg_global[i]*w_tilde[i]*ddx_w_tilde[i] - ddx_rho_u_tilde_R33[i]);
        }
        
        std::vector<double> w_p_tau13_p(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            w_p_tau13_p[i] = w_tau13_avg_global[i] - w_avg_global[i]*tau13_avg_global[i];
        }
        
        std::vector<double> two_ddx_w_p_tau13_p = computeDerivativeOfVector1D(
            w_p_tau13_p,
            dx);
        
        for (int i = 0; i < num_cells; i++)
        {
            two_ddx_w_p_tau13_p[i] *= double(2);
        }
        
        /*
         * Compute term V.
         */
        
        std::vector<double> two_p_p_ddz_w_p(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            two_p_p_ddz_w_p[i] = double(2)*(p_ddz_w_avg_global[i] - p_avg_global[i]*ddz_w_avg_global[i]);
        }
        
        /*
         * Compute term VI.
         */
        
        std::vector<double> tau13_p_ddx_w_p(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            tau13_p_ddx_w_p[i] = tau13_ddx_w_avg_global[i] - tau13_avg_global[i]*ddx_w_avg_global[i];
        }
        
        std::vector<double> tau23_p_ddy_w_p(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            tau23_p_ddy_w_p[i] = tau23_ddy_w_avg_global[i] - tau23_avg_global[i]*ddy_w_avg_global[i];
        }
        
        std::vector<double> tau33_p_ddz_w_p(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            tau33_p_ddz_w_p[i] = tau33_ddz_w_avg_global[i] - tau33_avg_global[i]*ddz_w_avg_global[i];
        }
        
        std::vector<double> m_2tau_3i_p_ddxi_w_p(num_cells, double(0));
        
        for (int i = 0; i < num_cells; i++)
        {
            m_2tau_3i_p_ddxi_w_p[i] = -double(2)*(tau13_p_ddx_w_p[i] + tau23_p_ddy_w_p[i] + tau33_p_ddz_w_p[i]);
        }
        
        /*
         * Output budget.
         */
        
        f_out.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        f_out.write((char*)&output_time, sizeof(double));
        f_out.write((char*)&rho_R33[0], sizeof(double)*rho_R33.size());
        // Term II.
        f_out.write((char*)&ddx_rho_u_tilde_R33[0], sizeof(double)*ddx_rho_u_tilde_R33.size());
        
        // Term IV(1).
        f_out.write((char*)&m_ddx_rho_w_pp_w_pp_u_pp[0], sizeof(double)*m_ddx_rho_w_pp_w_pp_u_pp.size());
        // Term IV(2).
        f_out.write((char*)&two_ddx_w_p_tau13_p[0], sizeof(double)*two_ddx_w_p_tau13_p.size());
        
        // Term V.
        f_out.write((char*)&two_p_p_ddz_w_p[0], sizeof(double)*two_p_p_ddz_w_p.size());
        
        // Term VI.
        f_out.write((char*)&m_2tau_3i_p_ddxi_w_p[0], sizeof(double)*m_2tau_3i_p_ddxi_w_p.size());
        
        // Term II in moving frame of mixing layer.
        f_out.write((char*)&ddx_rho_a1_R33[0], sizeof(double)*ddx_rho_a1_R33.size());
        
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
RTIRMIBudgetsUtilities::getAveragedShearStressComponentWithInhomogeneousXDirection(
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
RTIRMIBudgetsUtilities::getAveragedDerivativeOfShearStressComponentWithInhomogeneousXDirection(
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
RTIRMIBudgetsUtilities::getAveragedQuantityWithShearStressComponentWithInhomogeneousXDirection(
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
                        << ": RTIRMIBudgetsUtilities::"
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
                        << ": RTIRMIBudgetsUtilities::"
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
 * Compute averaged value (on product of variable derivatives and derivative of shear stress component) with only
 * x direction as inhomogeneous direction.
 * component_idx:
 * 0: tau11
 * 1: tau12
 * 2: tau13
 * 3: tau22
 * 4: tau23
 * 5: tau33
 */
std::vector<double>
RTIRMIBudgetsUtilities::getAveragedQuantityWithDerivativeOfShearStressComponentWithInhomogeneousXDirection(
    const std::vector<std::string>& quantity_names,
    const std::vector<int>& component_indices,
    const std::vector<bool>& use_derivative,
    const std::vector<int>& derivative_directions,
    const std::vector<bool>& use_reciprocal,
    const int shear_stress_component_idx,
    const int shear_stress_derivative_direction,
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
                        << ": RTIRMIBudgetsUtilities::"
                        << "getAveragedQuantityWithDerivativeOfShearStressComponentWithInhomogeneousXDirection():\n"
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
                        << ": RTIRMIBudgetsUtilities::"
                        << "getAveragedQuantityWithDerivativeOfShearStressComponentWithInhomogeneousXDirection():\n"
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
                        << "getAveragedQuantityWithDerivativeOfShearStressComponentWithInhomogeneousXDirection():\n"
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
                 * Register the patch and the quantities in the flow model and compute the
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
                         * Compute the derivative of shear stress component.
                         */
                        
                        double tau_ij_der = double(0);
                        
                        if ((shear_stress_component_idx == 0) && (shear_stress_derivative_direction == 0))
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
                            
                            tau_ij_der = (double(4)/double(3)*dmudx + dmu_vdx)*dudx 
                                + (double(4)/double(3)*mu[3] + mu_v[3])*d2udx2;
                        }
                        else
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "Cannot compute derivative of shear stress component for one-dimensional problem!\n"
                                << "shear_stress_component_idx = " << shear_stress_component_idx << " and shear_stress_derivative_direction = "
                                << shear_stress_derivative_direction << " given!\n"
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
                            
                            avg_local[idx_fine] += (avg*tau_ij_der/((double) n_overlapped));
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
                 * Register the patch and the quantities in the flow model and compute the
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
                             * Compute the derivative shear stress component.
                             */
                            
                            double tau_ij_der = double(0);
                            
                            if ((shear_stress_component_idx == 0) && (shear_stress_derivative_direction == 0))
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
                                
                                tau_ij_der = (double(4)/double(3)*dmudx + dmu_vdx)*dudx
                                    + (double(4)/double(3)*mu[3] + mu_v[3])*d2udx2
                                    - dvisc_dvdydx;
                            }
                            else if ((shear_stress_component_idx == 0) && (shear_stress_derivative_direction == 1))
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
                                
                                tau_ij_der = dvisc_dudxdy
                                    - (double(2)/double(3)*dmudy - dmu_vdy)*dvdy
                                    - (double(2)/double(3)*mu[3] - mu_v[3])*d2vdy2;
                            }
                            else if ((shear_stress_component_idx == 1) && (shear_stress_derivative_direction == 0))
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
                                
                                tau_ij_der = dvisc_dudydx
                                    + dmudx*dvdx
                                    + mu[3]*d2vdx2;
                            }
                            else if ((shear_stress_component_idx == 1) && (shear_stress_derivative_direction == 1))
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
                                
                                tau_ij_der = dvisc_dvdxdy
                                    + dmudy*dudy
                                    + mu[3]*d2udy2;
                            }
                            else
                            {
                                TBOX_ERROR(d_object_name
                                    << ": "
                                    << "Cannot compute derivative of shear stress component for two-dimensional problem!\n"
                                    << "shear_stress_component_idx = " << shear_stress_component_idx << " and shear_stress_derivative_direction = "
                                    << shear_stress_derivative_direction << " given!\n"
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
                                
                                avg_local[idx_fine] += (avg*tau_ij_der*weight/((double) n_overlapped));
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
                 * Register the patch and the quantities in the flow model and compute the
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
                                 * Compute the derivative of shear stress component.
                                 */
                                
                                double tau_ij_der = double(0);
                                
                                if ((shear_stress_component_idx == 0) && (shear_stress_derivative_direction == 0))
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
                                    
                                    tau_ij_der = (double(4)/double(3)*dmudx + dmu_vdx)*dudx
                                        + (double(4)/double(3)*mu[3] + mu_v[3])*d2udx2
                                        - dvisc_dvdydx
                                        - dvisc_dwdzdx;
                                }
                                else if ((shear_stress_component_idx == 0) && (shear_stress_derivative_direction == 1))
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
                                    
                                    tau_ij_der = dvisc_dudxdy
                                        - dvisc_dwdzdy
                                        - (double(2)/double(3)*dmudy - dmu_vdy)*dvdy
                                        - (double(2)/double(3)*mu[3] - mu_v[3])*d2vdy2;
                                }
                                else if ((shear_stress_component_idx == 0) && (shear_stress_derivative_direction == 2))
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
                                    
                                    tau_ij_der = dvisc_dudxdz
                                        - dvisc_dvdydz
                                        - (double(2)/double(3)*dmudz - dmu_vdz)*dwdz
                                        - (double(2)/double(3)*mu[3] - mu_v[3])*d2wdz2;
                                }
                                else if ((shear_stress_component_idx == 1) && (shear_stress_derivative_direction == 0))
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
                                    
                                    tau_ij_der = dmudx*dvdx
                                        + mu[3]*d2vdx2
                                        + dvisc_dudydx;
                                }
                                else if ((shear_stress_component_idx == 1) && (shear_stress_derivative_direction == 1))
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
                                    
                                    tau_ij_der = dmudy*dudy
                                        + mu[3]*d2udy2
                                        + dvisc_dvdxdy;
                                }
                                else if ((shear_stress_component_idx == 1) && (shear_stress_derivative_direction == 2))
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
                                    
                                    tau_ij_der = dvisc_dvdxdz
                                        + dvisc_dudydz;
                                }
                                else if ((shear_stress_component_idx == 2) && (shear_stress_derivative_direction == 0))
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
                                    
                                    tau_ij_der = dmudx*dwdx
                                        + mu[3]*d2wdx2
                                        + dvisc_dudzdx;
                                }
                                else if ((shear_stress_component_idx == 2) && (shear_stress_derivative_direction == 1))
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
                                    
                                    tau_ij_der = dvisc_dwdxdy
                                        + dvisc_dudzdy;
                                }
                                else if ((shear_stress_component_idx == 2) && (shear_stress_derivative_direction == 2))
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
                                    
                                    tau_ij_der = dmudz*dudz
                                        + mu[3]*d2udz2
                                        + dvisc_dwdxdz;
                                }
                                else
                                {
                                    TBOX_ERROR(d_object_name
                                        << ": "
                                        << "Cannot compute derivative of shear stress component for three-dimensional problem!\n"
                                        << "shear_stress_component_idx = " << shear_stress_component_idx << " and shear_stress_derivative_direction = "
                                        << shear_stress_derivative_direction << " given!\n"
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
                                    
                                    avg_local[idx_fine] += (avg*tau_ij_der*weight/((double) n_overlapped));
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
RTIRMIBudgetsUtilities::getGravityVector() const
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
RTIRMIBudgetsUtilities::computeDerivativeOfVector1D(
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
    // DO NOTHING.
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
        d_ensemble_statistics = HAMERS_MAKE_SHARED<EnsembleBudgetsRTIRMI>("d_ensemble_statistics");
        d_is_ensemble_statistics_initialized = true;
    }
    
    HAMERS_SHARED_PTR<RTIRMIBudgetsUtilities> rti_rmi_budgets_utilities(
        new RTIRMIBudgetsUtilities(
            "RTI RMI budgets utilities",
            d_dim,
            d_grid_geometry,
            d_num_species,
            d_flow_model,
            d_equation_of_state_mixing_rules,
            d_equation_of_mass_diffusivity_mixing_rules,
            d_equation_of_shear_viscosity_mixing_rules,
            d_equation_of_bulk_viscosity_mixing_rules,
            d_equation_of_thermal_conductivity_mixing_rules,
            HAMERS_DYNAMIC_POINTER_CAST<EnsembleBudgetsRTIRMI>(d_ensemble_statistics)));
    
    // Statistics are not computed for this realization yet.
    rti_rmi_budgets_utilities->d_ensemble_statistics->setVariablesNotComputed();
    
    // Compute the averaged grid level number no matter what.
    rti_rmi_budgets_utilities->
        computeAveragedGridLevelNumberWithHomogeneityInYDirectionOrInYZPlane(
            patch_hierarchy);
    
    // Loop over statistical quantities.
    for (int qi = 0; qi < static_cast<int>(d_statistical_quantities.size()); qi++)
    {
        // Get the key of the current variable.
        std::string statistical_quantity_key = d_statistical_quantities[qi];
        
        // Spatial profiles.
        if (statistical_quantity_key == "MASS_FRACTION_AVG_SP" ||
            statistical_quantity_key == "DENSITY_AVG_SP" ||
            statistical_quantity_key == "VELOCITY_X_AVG_SP" ||
            statistical_quantity_key == "VELOCITY_X_FAVRE_AVG_SP" ||
            statistical_quantity_key == "TURB_MASS_FLUX_VEL_X_SP" ||
            statistical_quantity_key == "PRESSURE_GRADIENT_X_AVG_SP" ||
            statistical_quantity_key == "rho_a1_budget_SP" ||
            statistical_quantity_key == "rho_K_budget_SP" ||
            statistical_quantity_key == "rho_R11_budget_SP" ||
            statistical_quantity_key == "rho_R22_budget_SP" ||
            statistical_quantity_key == "rho_R33_budget_SP")
        {
            rti_rmi_budgets_utilities->
                computeAveragedQuantitiesWithHomogeneityInYDirectionOrInYZPlane(
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
    
    HAMERS_SHARED_PTR<RTIRMIBudgetsUtilities> rti_rmi_budgets_utilities(
        new RTIRMIBudgetsUtilities(
            "RTI RMI budgets utilities",
            d_dim,
            d_grid_geometry,
            d_num_species,
            d_flow_model,
            d_equation_of_state_mixing_rules,
            d_equation_of_mass_diffusivity_mixing_rules,
            d_equation_of_shear_viscosity_mixing_rules,
            d_equation_of_bulk_viscosity_mixing_rules,
            d_equation_of_thermal_conductivity_mixing_rules,
            HAMERS_DYNAMIC_POINTER_CAST<EnsembleBudgetsRTIRMI>(d_ensemble_statistics)));
    
    // Output the averaged grid level number no matter what.
    rti_rmi_budgets_utilities->outputSpatialProfileEnsembleAveragedGridLevelNumberWithHomogeneityInYDirectionOrInYZPlane(
        "grid_level_num_avg.dat",
        output_time);
            
    // Loop over statistical quantities.
    for (int qi = 0; qi < static_cast<int>(d_statistical_quantities.size()); qi++)
    {
        // Get the key of the current variable.
        std::string statistical_quantity_key = d_statistical_quantities[qi];
        
        // Spatial profiles.
        if (statistical_quantity_key == "MASS_FRACTION_AVG_SP")
        {
            rti_rmi_budgets_utilities->
                outputSpatialProfileEnsembleAveragedMassFractionWithHomogeneityInYDirectionOrInYZPlane(
                    "Y_avg.dat",
                    output_time);
        }
        else if (statistical_quantity_key == "DENSITY_AVG_SP")
        {
            rti_rmi_budgets_utilities->
                outputSpatialProfileEnsembleAveragedDensityWithHomogeneityInYDirectionOrInYZPlane(
                    "rho_avg.dat",
                    output_time);
        }
        else if (statistical_quantity_key == "VELOCITY_X_AVG_SP")
        {
            rti_rmi_budgets_utilities->
                outputSpatialProfileEnsembleAveragedVelocityXWithHomogeneityInYDirectionOrInYZPlane(
                    "u_avg.dat",
                    output_time);
        }
        else if (statistical_quantity_key == "VELOCITY_X_FAVRE_AVG_SP")
        {
            rti_rmi_budgets_utilities->
                outputSpatialProfileEnsembleFavreAveragedVelocityXWithHomogeneityInYDirectionOrInYZPlane(
                    "u_tilde.dat",
                    output_time);
        }
        else if (statistical_quantity_key == "TURB_MASS_FLUX_VEL_X_SP")
        {
            rti_rmi_budgets_utilities->
                outputSpatialProfileEnsembleTurbulentMassFluxVelocityXWithHomogeneityInYDirectionOrInYZPlane(
                    "a1.dat",
                    output_time);
        }
        else if (statistical_quantity_key == "PRESSURE_GRADIENT_X_AVG_SP")
        {
            rti_rmi_budgets_utilities->
                outputSpatialProfileEnsembleAveragedPressureGradientXWithHomogeneityInYDirectionOrInYZPlane(
                    "ddx_p_avg.dat",
                    output_time);
        }
        else if (statistical_quantity_key == "rho_a1_budget_SP")
        {
            rti_rmi_budgets_utilities->
                outputBudgetTurbMassFluxXWithInhomogeneousXDirection(
                    "rho_a1_budget.dat",
                    patch_hierarchy,
                    output_time);
        }
        else if (statistical_quantity_key == "rho_K_budget_SP")
        {
            rti_rmi_budgets_utilities->
                outputBudgetFavreMeanTKEWithInhomogeneousXDirection(
                    "rho_K_budget.dat",
                    patch_hierarchy,
                    output_time);
        }
        else if (statistical_quantity_key == "rho_R11_budget_SP")
        {
            rti_rmi_budgets_utilities->
                outputBudgetReynoldsNormalStressInXDirectionWithInhomogeneousXDirection(
                    "rho_R11_budget.dat",
                    patch_hierarchy,
                    output_time);
        }
        else if (statistical_quantity_key == "rho_R22_budget_SP")
        {
            rti_rmi_budgets_utilities->
                outputBudgetReynoldsNormalStressInYDirectionWithInhomogeneousXDirection(
                    "rho_R22_budget.dat",
                    patch_hierarchy,
                    output_time);
        }
        else if (statistical_quantity_key == "rho_R33_budget_SP")
        {
            rti_rmi_budgets_utilities->
                outputBudgetReynoldsNormalStressInZDirectionWithInhomogeneousXDirection(
                    "rho_R33_budget.dat",
                    patch_hierarchy,
                    output_time);
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
