#ifndef FLOW_MODEL_STATISTICS_UTILITIES_FOUR_EQN_CONSERVATIVE_HPP
#define FLOW_MODEL_STATISTICS_UTILITIES_FOUR_EQN_CONSERVATIVE_HPP

#include "flow/flow_models/FlowModelStatisticsUtilities.hpp"
#include "util/mixing_rules/equations_of_bulk_viscosity/EquationOfBulkViscosityMixingRulesManager.hpp"
#include "util/mixing_rules/equations_of_mass_diffusivity/EquationOfMassDiffusivityMixingRulesManager.hpp"
#include "util/mixing_rules/equations_of_shear_viscosity/EquationOfShearViscosityMixingRulesManager.hpp"
#include "util/mixing_rules/equations_of_state/EquationOfStateMixingRulesManager.hpp"
#include "util/mixing_rules/equations_of_thermal_conductivity/EquationOfThermalConductivityMixingRulesManager.hpp"

class FlowModelStatisticsUtilitiesFourEqnConservative: public FlowModelStatisticsUtilities
{
    public:
        FlowModelStatisticsUtilitiesFourEqnConservative(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_species,
            const boost::shared_ptr<tbox::Database>& flow_model_db,
            const boost::shared_ptr<EquationOfStateMixingRules> equation_of_state_mixing_rules,
            const boost::shared_ptr<EquationOfMassDiffusivityMixingRules> equation_of_mass_diffusivity_mixing_rules,
            const boost::shared_ptr<EquationOfShearViscosityMixingRules> equation_of_shear_viscosity_mixing_rules,
            const boost::shared_ptr<EquationOfBulkViscosityMixingRules> equation_of_bulk_viscosity_mixing_rules,
            const boost::shared_ptr<EquationOfThermalConductivityMixingRules> equation_of_thermal_conductivity_mixing_rules):
                FlowModelStatisticsUtilities(
                    object_name,
                    dim,
                    grid_geometry,
                    num_species,
                    flow_model_db),
                d_equation_of_state_mixing_rules(equation_of_state_mixing_rules),
                d_equation_of_mass_diffusivity_mixing_rules(equation_of_mass_diffusivity_mixing_rules),
                d_equation_of_shear_viscosity_mixing_rules(equation_of_shear_viscosity_mixing_rules),
                d_equation_of_bulk_viscosity_mixing_rules(equation_of_bulk_viscosity_mixing_rules),
                d_equation_of_thermal_conductivity_mixing_rules(equation_of_thermal_conductivity_mixing_rules)
        {}
        
        /*
         * Output names of statistical quantities to output to a file.
         */
        void
        outputStatisticalQuantitiesNames(
            const std::string& stat_dump_filename);
        
        /*
         * Output statisitcal quantities to a file.
         */
        void
        outputStatisticalQuantities(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
    private:
        /*
         * Output mixing width in x-direction to a file.
         */
        void
        outputMixingWidthInXDirection(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mixing width in y-direction to a file.
         */
        void
        outputMixingWidthInYDirection(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mixing width in z-direction to a file.
         */
        void
        outputMixingWidthInZDirection(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mixing width with mole fraction in x-direction to a file.
         */
        void
        outputMixingWidthMoleFractionInXDirection(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mixing width with mole fraction in y-direction to a file.
         */
        void
        outputMixingWidthMoleFractionInYDirection(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mixing width with mole fraction in z-direction to a file.
         */
        void
        outputMixingWidthMoleFractionInZDirection(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mixedness in x-direction to a file.
         */
        void
        outputMixednessInXDirection(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mixedness in y-direction to a file.
         */
        void
        outputMixednessInYDirection(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mixedness in z-direction to a file.
         */
        void
        outputMixednessInZDirection(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mixedness with mole fraction in x-direction to a file.
         */
        void
        outputMixednessMoleFractionInXDirection(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mixedness with mole fraction in y-direction to a file.
         */
        void
        outputMixednessMoleFractionInYDirection(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mixedness with mole fraction in z-direction to a file.
         */
        void
        outputMixednessMoleFractionInZDirection(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output TKE integrated with assumed homogeneity in x-direction to a file.
         */
        void
        outputTKEIntegratedWithHomogeneityInXDirection(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output TKE integrated with assumed homogeneity in y-direction to a file.
         */
        void
        outputTKEIntegratedWithHomogeneityInYDirection(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output TKE integrated with assumed homogeneity in z-direction to a file.
         */
        void
        outputTKEIntegratedWithHomogeneityInZDirection(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output TKE integrated with assumed homogeneity in xy-plane to a file.
         */
        void
        outputTKEIntegratedWithHomogeneityInXYPlane(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output TKE integrated with assumed homogeneity in yz-plane to a file.
         */
        void
        outputTKEIntegratedWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output TKE integrated with assumed homogeneity in xz-plane to a file.
         */
        void
        outputTKEIntegratedWithHomogeneityInXZPlane(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output TKE in x-direction integrated with assumed homogeneity in y-direction to a file.
         */
        void
        outputTKEInXDirectionIntegratedWithHomogeneityInYDirection(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output TKE in y-direction integrated with assumed homogeneity in y-direction to a file.
         */
        void
        outputTKEInYDirectionIntegratedWithHomogeneityInYDirection(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output TKE in x-direction integrated with assumed homogeneity in yz-plane to a file.
         */
        void
        outputTKEInXDirectionIntegratedWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output TKE in y-direction integrated with assumed homogeneity in yz-plane to a file.
         */
        void
        outputTKEInYDirectionIntegratedWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output TKE in z-direction integrated with assumed homogeneity in yz-plane to a file.
         */
        void
        outputTKEInZDirectionIntegratedWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output enstrophy integrated to a file.
         */
        void
        outputEnstrophyIntegrated(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output scalar dissipation rate of first species integrated to a file.
         */
        void
        outputScalarDissipationRateIntegrated(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mean Reynolds number inside mixing layer with assumed homogeneity in y-direction to
         * a file.
         */
        void
        outputReynoldsNumberMeanInMixingLayerWithHomogeneityInYDirection(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mean Reynolds number inside mixing layer with assumed homogeneity in yz-plane to
         * a file.
         */
        void
        outputReynoldsNumberMeanInMixingLayerWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mean turbulent Mach number inside mixing layer with assumed homogeneity in y-direction to a file.
         */
        void
        outputTurbulentMachNumberMeanInMixingLayerWithHomogeneityInYDirection(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mean turbulent Mach number inside mixing layer with assumed homogeneity in yz-plane to a file.
         */
        void
        outputTurbulentMachNumberMeanInMixingLayerWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mean TKE inside mixing layer with assumed homogeneity in y-direction to a file.
         */
        void
        outputTKEMeanInMixingLayerWithHomogeneityInYDirection(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mean TKE in x-direction inside mixing layer with assumed homogeneity in y-direction
         * to a file.
         */
        void
        outputTKEInXDirectionMeanInMixingLayerWithHomogeneityInYDirection(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mean TKE in y-direction inside mixing layer with assumed homogeneity in y-direction
         * to a file.
         */
        void
        outputTKEInYDirectionMeanInMixingLayerWithHomogeneityInYDirection(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mean TKE inside mixing layer with assumed homogeneity in yz-plane to a file.
         */
        void
        outputTKEMeanInMixingLayerWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mean TKE in x-direction inside mixing layer with assumed homogeneity in yz-plane
         * to a file.
         */
        void
        outputTKEInXDirectionMeanInMixingLayerWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mean TKE in y-direction inside mixing layer with assumed homogeneity in yz-plane
         * to a file.
         */
        void
        outputTKEInYDirectionMeanInMixingLayerWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mean TKE in z-direction inside mixing layer with assumed homogeneity in yz-plane
         * to a file.
         */
        void
        outputTKEInZDirectionMeanInMixingLayerWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mean Reynolds normal stress component in x-direction inside mixing layer with
         * assumed homogeneity in y-direction to a file.
         */
        void
        outputReynoldsNormalStressInXDirectionMeanInMixingLayerWithHomogeneityInYDirection(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mean Reynolds normal stress component in y-direction inside mixing layer with
         * assumed homogeneity in y-direction to a file.
         */
        void
        outputReynoldsNormalStressInYDirectionMeanInMixingLayerWithHomogeneityInYDirection(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mean Reynolds normal stress component in x-direction inside mixing layer with
         * assumed homogeneity in yz-plane to a file.
         */
        void
        outputReynoldsNormalStressInXDirectionMeanInMixingLayerWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mean Reynolds normal stress component in y-direction inside mixing layer with
         * assumed homogeneity in yz-plane to a file.
         */
        void
        outputReynoldsNormalStressInYDirectionMeanInMixingLayerWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mean Reynolds normal stress component in z-direction inside mixing layer with
         * assumed homogeneity in yz-plane to a file.
         */
        void
        outputReynoldsNormalStressInZDirectionMeanInMixingLayerWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mean Reynolds shear stress component in x- and y-directions inside mixing layer
         * with assumed homogeneity in y-direction to a file.
         */
        void
        outputReynoldsShearStressInXYDirectionsMeanInMixingLayerWithHomogeneityInYDirection(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mean Reynolds shear stress component in x- and y-directions inside mixing layer
         * with assumed homogeneity in yz-plane to a file.
         */
        void
        outputReynoldsShearStressInXYDirectionsMeanInMixingLayerWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mean Reynolds shear stress component in x- and z-directions inside mixing layer
         * with assumed homogeneity in yz-plane to a file.
         */
        void
        outputReynoldsShearStressInXZDirectionsMeanInMixingLayerWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mean Reynolds shear stress component in y- and z-directions inside mixing layer
         * with assumed homogeneity in yz-plane to a file.
         */
        void
        outputReynoldsShearStressInYZDirectionsMeanInMixingLayerWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mean Reynolds normal stress anisotropy component in x-direction inside mixing layer with
         * assumed homogeneity in y-direction to a file.
         */
        void
        outputReynoldsNormalStressAnisotropyInXDirectionMeanInMixingLayerWithHomogeneityInYDirection(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mean Reynolds normal stress anisotropy component in y-direction inside mixing layer with
         * assumed homogeneity in y-direction to a file.
         */
        void
        outputReynoldsNormalStressAnisotropyInYDirectionMeanInMixingLayerWithHomogeneityInYDirection(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mean Reynolds normal stress anisotropy component in x-direction inside mixing layer with
         * assumed homogeneity in yz-plane to a file.
         */
        void
        outputReynoldsNormalStressAnisotropyInXDirectionMeanInMixingLayerWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mean Reynolds normal stress anisotropy component in y-direction inside mixing layer with
         * assumed homogeneity in yz-plane to a file.
         */
        void
        outputReynoldsNormalStressAnisotropyInYDirectionMeanInMixingLayerWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mean Reynolds normal stress anisotropy component in z-direction inside mixing layer with
         * assumed homogeneity in yz-plane to a file.
         */
        void
        outputReynoldsNormalStressAnisotropyInZDirectionMeanInMixingLayerWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mean Reynolds shear stress anisotropy component in x- and y-directions inside mixing layer
         * with assumed homogeneity in y-direction to a file.
         */
        void
        outputReynoldsShearStressAnisotropyInXYDirectionsMeanInMixingLayerWithHomogeneityInYDirection(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mean Reynolds shear stress anisotropy component in x- and y-directions inside mixing layer
         * with assumed homogeneity in yz-plane to a file.
         */
        void
        outputReynoldsShearStressAnisotropyInXYDirectionsMeanInMixingLayerWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mean Reynolds shear stress anisotropy component in x- and z-directions inside mixing layer
         * with assumed homogeneity in yz-plane to a file.
         */
        void
        outputReynoldsShearStressAnisotropyInXZDirectionsMeanInMixingLayerWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mean Reynolds shear stress anisotropy component in y- and z-directions inside mixing layer
         * with assumed homogeneity in yz-plane to a file.
         */
        void
        outputReynoldsShearStressAnisotropyInYZDirectionsMeanInMixingLayerWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mean turbulent mass flux in x-direction inside mixing layer with assumed homogeneity
         * in y-direction to a file.
         */
        void
        outputTurbulentMassFluxInXDirectionMeanInMixingLayerWithHomogeneityInYDirection(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mean turbulent mass flux in x-direction inside mixing layer with assumed homogeneity
         * in yz-plane to a file.
         */
        void
        outputTurbulentMassFluxInXDirectionMeanInMixingLayerWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mean density specific volume covariance inside mixing layer with assumed homogeneity
         * in y-direction to a file.
         */
        void
        outputDensitySpecificVolumeCovarianceMeanInMixingLayerWithHomogeneityInYDirection(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mean density specific volume covariance inside mixing layer with assumed homogeneity
         * in yz-plane to a file.
         */
        void
        outputDensitySpecificVolumeCovarianceMeanInMixingLayerWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mean density inside mixing layer with assumed homogeneity in y-direction to a file.
         */
        void
        outputDensityMeanInMixingLayerWithHomogeneityInYDirection(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mean density inside mixing layer with assumed homogeneity in yz-plane to a file.
         */
        void
        outputDensityMeanInMixingLayerWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mean deviation from Boussinesq approximation inside mixing layer with assumed
         * homogeneity in y-direction to a file.
         */
        void
        outputBoussinesqDeviationMeanInMixingLayerWithHomogeneityInYDirection(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mean deviation from Boussinesq approximation inside mixing layer with assumed
         * homogeneity in yz-plane to a file.
         */
        void
        outputBoussinesqDeviationMeanInMixingLayerWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mean effective Atwood number inside mixing layer with assumed homogeneity in
         * y-direction to a file.
         */
        void
        outputEffectiveAtwoodNumberMeanInMixingLayerWithHomogeneityInYDirection(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mean effective Atwood number inside mixing layer with assumed homogeneity in
         * yz-plane to a file.
         */
        void
        outputEffectiveAtwoodNumberMeanInMixingLayerWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mean mass diffusivity inside mixing layer in x-direction to a file.
         */
        void
        outputMassDiffusivityMeanInMixingLayerInXDirection(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mean dynamic shear viscosity inside mixing layer in x-direction to a file.
         */
        void
        outputDynamicShearViscosityMeanInMixingLayerInXDirection(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mean kinematic shear viscosity inside mixing layer in x-direction to a file.
         */
        void
        outputKinematicShearViscosityMeanInMixingLayerInXDirection(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mean dynamic bulk viscosity inside mixing layer in x-direction to a file.
         */
        void
        outputDynamicBulkViscosityMeanInMixingLayerInXDirection(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mean kinematic bulk viscosity inside mixing layer in x-direction to a file.
         */
        void
        outputKinematicBulkViscosityMeanInMixingLayerInXDirection(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mean thermal conductivity inside mixing layer in x-direction to a file.
         */
        void
        outputThermalConductivityMeanInMixingLayerInXDirection(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mean thermal diffusivity inside mixing layer in x-direction to a file.
         */
        void
        outputThermalDiffusivityMeanInMixingLayerInXDirection(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mixing layer width 1 in x-direction to a file.
         */
        void
        outputMixingLayerWidth1InXDirection(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mixing layer width 2 in x-direction to a file.
         */
        void
        outputMixingLayerWidth2InXDirection(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mixing layer width 3 in x-direction to a file.
         */
        void
        outputMixingLayerWidth3InXDirection(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mixing layer width 4 in x-direction to a file.
         */
        void
        outputMixingLayerWidth4InXDirection(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mixing layer width 5 in x-direction to a file.
         */
        void
        outputMixingLayerWidth5InXDirection(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mixing layer width 6 in x-direction to a file.
         */
        void
        outputMixingLayerWidth6InXDirection(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mixing layer width 7 in x-direction to a file.
         */
        void
        outputMixingLayerWidth7InXDirection(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mixing layer width 8 in x-direction to a file.
         */
        void
        outputMixingLayerWidth8InXDirection(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mixing layer width 9 in x-direction to a file.
         */
        void
        outputMixingLayerWidth9InXDirection(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output numerical interface thickness to a file.
         */
        void
        outputNumericalInterfaceThickness(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output number of cells to a file.
         */
        void
        outputNumberOfCells(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output weighted number of cells to a file.
         */
        void
        outputWeightedNumberOfCells(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
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

#endif /* FLOW_MODEL_STATISTICS_UTILITIES_FOUR_EQN_CONSERVATIVE_HPP */
