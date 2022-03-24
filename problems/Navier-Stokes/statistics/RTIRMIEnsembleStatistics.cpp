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
            Y_0_avg_computed      = false;
            Y_0_sq_avg_computed   = false;
            Y_0_Y_1_avg_computed  = false;
            rho_avg_computed      = false;
            rho_sq_avg_computed   = false;
            u_avg_computed        = false;
            v_avg_computed        = false;
            w_avg_computed        = false;
            rho_u_avg_computed    = false;
            rho_v_avg_computed    = false;
            rho_w_avg_computed    = false;
            rho_u_u_avg_computed  = false;
            rho_v_v_avg_computed  = false;
            rho_w_w_avg_computed  = false;
        }
        
        void clearAllData()
        {
            Y_0_avg_realizations.clear();
            Y_0_sq_avg_realizations.clear();
            Y_0_Y_1_avg_realizations.clear();
            rho_avg_realizations.clear();
            rho_sq_avg_realizations.clear();
            u_avg_realizations.clear();
            v_avg_realizations.clear();
            w_avg_realizations.clear();
            rho_u_avg_realizations.clear();
            rho_v_avg_realizations.clear();
            rho_w_avg_realizations.clear();
            rho_u_u_avg_realizations.clear();
            rho_v_v_avg_realizations.clear();
            rho_w_w_avg_realizations.clear();
            
            Y_0_avg_computed      = false;
            Y_0_sq_avg_computed   = false;
            Y_0_Y_1_avg_computed  = false;
            rho_avg_computed      = false;
            rho_sq_avg_computed   = false;
            u_avg_computed        = false;
            v_avg_computed        = false;
            w_avg_computed        = false;
            rho_u_avg_computed    = false;
            rho_v_avg_computed    = false;
            rho_w_avg_computed    = false;
            rho_u_u_avg_computed  = false;
            rho_v_v_avg_computed  = false;
            rho_w_w_avg_computed  = false;
        }
        
        // Scratch arrays.
        // Number of realizalizations; number of cells.
        std::vector<std::vector<double> > Y_0_avg_realizations;
        std::vector<std::vector<double> > Y_0_sq_avg_realizations;
        std::vector<std::vector<double> > Y_0_Y_1_avg_realizations;
        std::vector<std::vector<double> > rho_avg_realizations;
        std::vector<std::vector<double> > rho_sq_avg_realizations;
        std::vector<std::vector<double> > u_avg_realizations;
        std::vector<std::vector<double> > v_avg_realizations;
        std::vector<std::vector<double> > w_avg_realizations;
        std::vector<std::vector<double> > rho_u_avg_realizations;
        std::vector<std::vector<double> > rho_v_avg_realizations;
        std::vector<std::vector<double> > rho_w_avg_realizations;
        std::vector<std::vector<double> > rho_u_u_avg_realizations;
        std::vector<std::vector<double> > rho_v_v_avg_realizations;
        std::vector<std::vector<double> > rho_w_w_avg_realizations;
        
        // Whether the scratch arrays are filled.
        bool Y_0_avg_computed;
        bool Y_0_sq_avg_computed;
        bool Y_0_Y_1_avg_computed;
        bool rho_avg_computed;
        bool rho_sq_avg_computed;
        bool u_avg_computed;
        bool v_avg_computed;
        bool w_avg_computed;
        bool rho_u_avg_computed;
        bool rho_v_avg_computed;
        bool rho_w_avg_computed;
        bool rho_u_u_avg_computed;
        bool rho_v_v_avg_computed;
        bool rho_w_w_avg_computed;
        
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
         * Output ensemble averaged mass fraction with assumed homogeneity in y-direction (2D) or yz-plane (3D)
         * to a file.
         */
        void
        outputEnsembleAveragedMassFractionWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double output_time) const;
        
        /*
         * Output ensemble variance of mass fraction with assumed homogeneity in y-direction (2D) or yz-plane (3D)
         * to a file.
         */
        void
        outputEnsembleMassFractionVarianceWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double output_time) const;
        
        /*
         * Output ensemble averaged density with assumed homogeneity in y-direction (2D) or yz-plane (3D)
         * to a file.
         */
        void
        outputEnsembleAveragedDensityWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double output_time) const;
        
        /*
         * Output ensemble variance of density with assumed homogeneity in y-direction (2D) or yz-plane (3D)
         * to a file.
         */
        void
        outputEnsembleDensityVarianceWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double output_time) const;
        
        /*
         * Output ensemble averaged velocity x-component with assumed homogeneity in y-direction (2D) or yz-plane (3D)
         * to a file.
         */
        void
        outputEnsembleAveragedVelocityXWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double output_time) const;
        
        /*
         * Output ensemble averaged velocity y-component with assumed homogeneity in y-direction (2D) or yz-plane (3D)
         * to a file.
         */
        void
        outputEnsembleAveragedVelocityYWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double output_time) const;
        
        /*
         * Output ensemble averaged velocity z-component with assumed homogeneity in yz-plane (3D) to a file.
         */
        void
        outputEnsembleAveragedVelocityZWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double output_time) const;
        
        /*
         * Output ensemble averaged momentum x-component with assumed homogeneity in y-direction (2D) or yz-plane (3D)
         * to a file.
         */
        void
        outputEnsembleAveragedMomentumXWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double output_time) const;
        
        /*
         * Output ensemble averaged momentum y-component with assumed homogeneity in y-direction (2D) or yz-plane (3D)
         * to a file.
         */
        void
        outputEnsembleAveragedMomentumYWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double output_time) const;
        
        /*
         * Output ensemble averaged momentum z-component with assumed homogeneity in yz-plane (3D) to a file.
         */
        void
        outputEnsembleAveragedMomentumZWithHomogeneitynYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double output_time) const;
        
        /*
         * Output ensemble turbulent mass flux in x-direction with assumed homogeneity in y-direction (2D) or yz-plane (3D)
         * to a file.
         */
        void
        outputEnsembleTurbulentMassFluxXWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double output_time) const;
        
        /*
         * Output ensemble Reynolds normal stress x-component with assumed homogeneity in y-direction (2D) or yz-plane (3D)
         * to a file.
         */
        void
        outputEnsembleReynoldsNormalStressXWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double output_time) const;
        
        /*
         * Output ensemble Reynolds normal stress y-component with assumed homogeneity in y-direction (2D) or yz-plane (3D)
         * to a file.
         */
        void
        outputEnsembleReynoldsNormalStressYWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double output_time) const;
        
        /*
         * Output ensemble Reynolds normal stress z-component with assumed homogeneity in yz-plane (3D) to a file.
         */
        void
        outputEnsembleReynoldsNormalStressZWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double output_time) const;
        
        /*
         * Output mixing width in x-direction to a file.
         */
        void
        outputMixingWidthInXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output mixedness in x-direction to a file.
         */
        void
        outputMixednessInXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output mean velocity associated with turbulent mass flux component in x-direction with assumed homogeneity
         * in y-direction (2D) or * yz-plane (3D) to a file.
         */
        void
        outputTurbulentMassFluxVelocityXWithInhomogeneousXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output Reynolds normal stress component in x-direction with assumed homogeneity in y-direction (2D) or yz-plane (3D)
         * to a file.
         */
        void
        outputReynoldsNormalStressXWithInhomogeneousXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output Reynolds normal stress component in y-direction with assumed homogeneity in y-direction (2D) or yz-plane (3D)
         * to a file.
         */
        void
        outputReynoldsNormalStressYWithInhomogeneousXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output Reynolds normal stress component in z-direction with assumed homogeneity in yz-plane (3D) to a file.
         */
        void
        outputReynoldsNormalStressZWithInhomogeneousXDirection(
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
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
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
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
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
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
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
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
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
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
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
 * Compute averaged velocity x-component with assumed homogeneity in y-direction (2D) or yz-plane (3D).
 */
void
RTIRMIStatisticsUtilities::computeAveragedVelocityXWithHomogeneityInYDirectionOrInYZPlane(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
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
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
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
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
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
 * Compute averaged momentum x-component with assumed homogeneity in y-direction (2D) or yz-plane (3D).
 */
void
RTIRMIStatisticsUtilities::computeAveragedMomentumXWithHomogeneityInYDirectionOrInYZPlane(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
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
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
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
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
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
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
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
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
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
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
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
 * Output ensemble averaged mass fraction with assumed homogeneity in y-direction (2D) or yz-plane (3D)
 * to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleAveragedMassFractionWithHomogeneityInYDirectionOrInYZPlane(
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
 * Output ensemble variance of mass fraction with assumed homogeneity in y-direction (2D) or yz-plane (3D)
 * to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleMassFractionVarianceWithHomogeneityInYDirectionOrInYZPlane(
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
 * Output ensemble averaged density with assumed homogeneity in y-direction (2D) or yz-plane (3D)
 * to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleAveragedDensityWithHomogeneityInYDirectionOrInYZPlane(
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
 * Output ensemble variance of density with assumed homogeneity in y-direction (2D) or yz-plane (3D)
 * to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleDensityVarianceWithHomogeneityInYDirectionOrInYZPlane(
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
 * Output ensemble averaged velocity x-component with assumed homogeneity in y-direction (2D) or yz-plane (3D)
 * to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleAveragedVelocityXWithHomogeneityInYDirectionOrInYZPlane(
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
 * Output ensemble averaged velocity y-component with assumed homogeneity in y-direction (2D) or yz-plane (3D)
 * to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleAveragedVelocityYWithHomogeneityInYDirectionOrInYZPlane(
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
 * Output ensemble averaged velocity z-component with assumed homogeneity in yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleAveragedVelocityZWithHomogeneityInYZPlane(
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
 * Output ensemble averaged momentum x-component with assumed homogeneity in y-direction (2D) or yz-plane (3D)
 * to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleAveragedMomentumXWithHomogeneityInYDirectionOrInYZPlane(
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
 * Output ensemble averaged momentum y-component with assumed homogeneity in y-direction (2D) or yz-plane (3D)
 * to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleAveragedMomentumYWithHomogeneityInYDirectionOrInYZPlane(
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
 * Output ensemble averaged momentum z-component with assumed homogeneity in yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleAveragedMomentumZWithHomogeneitynYZPlane(
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
 * Output ensemble turbulent mass flux in x-direction with assumed homogeneity in y-direction (2D) or yz-plane (3D)
 * to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleTurbulentMassFluxXWithHomogeneityInYDirectionOrInYZPlane(
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
 * Output ensemble Reynolds normal stress x-component with assumed homogeneity in y-direction (2D) or yz-plane (3D)
 * to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleReynoldsNormalStressXWithHomogeneityInYDirectionOrInYZPlane(
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
 * Output ensemble Reynolds normal stress y-component with assumed homogeneity in y-direction (2D) or yz-plane (3D)
 * to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleReynoldsNormalStressYWithHomogeneityInYDirectionOrInYZPlane(
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
 * Output ensemble Reynolds normal stress z-component with assumed homogeneity in yz-plane (3D) to a file.
 * to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnsembleReynoldsNormalStressZWithHomogeneityInYZPlane(
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
 * Output mixing width in x-direction to a file.
 */
void
RTIRMIStatisticsUtilities::outputMixingWidthInXDirection(
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
 * Output mixedness in x-direction to a file.
 */
void
RTIRMIStatisticsUtilities::outputMixednessInXDirection(
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
 * Output mean velocity associated with turbulent mass flux component in x-direction with assumed homogeneity
 * in y-direction (2D) or * yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputTurbulentMassFluxVelocityXWithInhomogeneousXDirection(
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
 * Output Reynolds normal stress component in x-direction with assumed homogeneity in y-direction (2D) or yz-plane (3D)
 * to a file.
 */
void
RTIRMIStatisticsUtilities::outputReynoldsNormalStressXWithInhomogeneousXDirection(
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
 * Output Reynolds normal stress component in y-direction with assumed homogeneity in y-direction (2D) or yz-plane (3D)
 * to a file.
 */
void
RTIRMIStatisticsUtilities::outputReynoldsNormalStressYWithInhomogeneousXDirection(
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
            const double u_tilde = rho_v_avg_global[i]/rho_avg_global[i];
            const double rho_v_pp_v_pp = rho_v_v_avg_global[i] - rho_v_avg_global[i]*u_tilde;
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
 * Output Reynolds normal stress component in z-direction with assumed homogeneity in yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputReynoldsNormalStressZWithInhomogeneousXDirection(
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
            const double u_tilde = rho_w_avg_global[i]/rho_avg_global[i];
            const double rho_w_pp_w_pp = rho_w_w_avg_global[i] - rho_w_avg_global[i]*u_tilde;
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
            else if (statistical_quantity_key == "MIXEDNESS_X")
            {
                f_out << "\t" << "MIXEDNESS_X          ";
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
    rti_rmi_statistics_utilities->d_ensemble_statistics->Y_0_avg_computed      = false;
    rti_rmi_statistics_utilities->d_ensemble_statistics->Y_0_sq_avg_computed   = false;
    rti_rmi_statistics_utilities->d_ensemble_statistics->rho_avg_computed      = false;
    rti_rmi_statistics_utilities->d_ensemble_statistics->rho_sq_avg_computed   = false;
    rti_rmi_statistics_utilities->d_ensemble_statistics->u_avg_computed        = false;
    rti_rmi_statistics_utilities->d_ensemble_statistics->v_avg_computed        = false;
    rti_rmi_statistics_utilities->d_ensemble_statistics->w_avg_computed        = false;
    rti_rmi_statistics_utilities->d_ensemble_statistics->rho_u_avg_computed    = false;
    rti_rmi_statistics_utilities->d_ensemble_statistics->rho_v_avg_computed    = false;
    rti_rmi_statistics_utilities->d_ensemble_statistics->rho_w_avg_computed    = false;
    rti_rmi_statistics_utilities->d_ensemble_statistics->rho_u_u_avg_computed  = false;
    rti_rmi_statistics_utilities->d_ensemble_statistics->rho_v_v_avg_computed  = false;
    rti_rmi_statistics_utilities->d_ensemble_statistics->rho_w_w_avg_computed  = false;
    
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
            rti_rmi_statistics_utilities->
                computeAveragedDensityWithHomogeneityInYDirectionOrInYZPlane(
                    patch_hierarchy,
                    data_context);
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
        else if (statistical_quantity_key == "VELOCITY_X_AVG_SP")
        {
            rti_rmi_statistics_utilities->
                computeAveragedVelocityXWithHomogeneityInYDirectionOrInYZPlane(
                    patch_hierarchy,
                    data_context);
        }
        else if (statistical_quantity_key == "VELOCITY_Y_AVG_SP")
        {
            rti_rmi_statistics_utilities->
                computeAveragedVelocityYWithHomogeneityInYDirectionOrInYZPlane(
                    patch_hierarchy,
                    data_context);
        }
        else if (statistical_quantity_key == "VELOCITY_Z_AVG_SP")
        {
            rti_rmi_statistics_utilities->
                computeAveragedVelocityZWithHomogeneityInYZPlane(
                    patch_hierarchy,
                    data_context);
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
                outputEnsembleAveragedMassFractionWithHomogeneityInYDirectionOrInYZPlane(
                    "Y_avg.dat",
                    patch_hierarchy,
                    data_context,
                    output_time);
        }
        else if (statistical_quantity_key == "MASS_FRACTION_VAR_SP")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleMassFractionVarianceWithHomogeneityInYDirectionOrInYZPlane(
                    "Y_var.dat",
                    patch_hierarchy,
                    data_context,
                    output_time);
        }
        else if (statistical_quantity_key == "DENSITY_AVG_SP")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleAveragedDensityWithHomogeneityInYDirectionOrInYZPlane(
                    "rho_avg.dat",
                    patch_hierarchy,
                    data_context,
                    output_time);
        }
        else if (statistical_quantity_key == "DENSITY_VAR_SP")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleDensityVarianceWithHomogeneityInYDirectionOrInYZPlane(
                    "rho_var.dat",
                    patch_hierarchy,
                    data_context,
                    output_time);
        }
        else if (statistical_quantity_key == "VELOCITY_X_AVG_SP")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleAveragedVelocityXWithHomogeneityInYDirectionOrInYZPlane(
                    "u_avg.dat",
                    patch_hierarchy,
                    data_context,
                    output_time);
        }
        else if (statistical_quantity_key == "VELOCITY_Y_AVG_SP")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleAveragedVelocityYWithHomogeneityInYDirectionOrInYZPlane(
                    "v_avg.dat",
                    patch_hierarchy,
                    data_context,
                    output_time);
        }
        else if (statistical_quantity_key == "VELOCITY_Z_AVG_SP")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleAveragedVelocityZWithHomogeneityInYZPlane(
                    "w_avg.dat",
                    patch_hierarchy,
                    data_context,
                    output_time);
        }
        else if (statistical_quantity_key == "MOMENTUM_X_AVG_SP")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleAveragedMomentumXWithHomogeneityInYDirectionOrInYZPlane(
                    "rho_u_avg.dat",
                    patch_hierarchy,
                    data_context,
                    output_time);
        }
        else if (statistical_quantity_key == "MOMENTUM_Y_AVG_SP")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleAveragedMomentumYWithHomogeneityInYDirectionOrInYZPlane(
                    "rho_v_avg.dat",
                    patch_hierarchy,
                    data_context,
                    output_time);
        }
        else if (statistical_quantity_key == "MOMENTUM_Z_AVG_SP")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleAveragedMomentumZWithHomogeneitynYZPlane(
                    "rho_w_avg.dat",
                    patch_hierarchy,
                    data_context,
                    output_time);
        }
        else if (statistical_quantity_key == "TURB_MASS_FLUX_X_SP")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleTurbulentMassFluxXWithHomogeneityInYDirectionOrInYZPlane(
                    "rho_a1.dat",
                    patch_hierarchy,
                    data_context,
                    output_time);
        }
        else if (statistical_quantity_key == "R11_SP")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleReynoldsNormalStressXWithHomogeneityInYDirectionOrInYZPlane(
                    "R11.dat",
                    patch_hierarchy,
                    data_context,
                    output_time);
        }
        else if (statistical_quantity_key == "R22_SP")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleReynoldsNormalStressYWithHomogeneityInYDirectionOrInYZPlane(
                    "R22.dat",
                    patch_hierarchy,
                    data_context,
                    output_time);
        }
        else if (statistical_quantity_key == "R33_SP")
        {
            rti_rmi_statistics_utilities->
                outputEnsembleReynoldsNormalStressZWithHomogeneityInYZPlane(
                    "R33.dat",
                    patch_hierarchy,
                    data_context,
                    output_time);
        }
        // Non-spatial profiles.
        else if (statistical_quantity_key == "MIXING_WIDTH_X")
        {
            rti_rmi_statistics_utilities->
                outputMixingWidthInXDirection(
                    stat_dump_filename,
                    patch_hierarchy,
                    data_context);
        }
        else if (statistical_quantity_key == "MIXEDNESS_X")
        {
            rti_rmi_statistics_utilities->
                outputMixednessInXDirection(
                    stat_dump_filename,
                    patch_hierarchy,
                    data_context);
        }
        else if (statistical_quantity_key == "a1_MEAN_INHOMO_X")
        {
            rti_rmi_statistics_utilities->
                outputTurbulentMassFluxVelocityXWithInhomogeneousXDirection(
                    stat_dump_filename,
                    patch_hierarchy,
                    data_context);
        }
        else if (statistical_quantity_key == "R11_MEAN_INHOMO_X")
        {
            rti_rmi_statistics_utilities->
                outputReynoldsNormalStressXWithInhomogeneousXDirection(
                    stat_dump_filename,
                    patch_hierarchy,
                    data_context);
        }
        else if (statistical_quantity_key == "R22_MEAN_INHOMO_X")
        {
            rti_rmi_statistics_utilities->
                outputReynoldsNormalStressYWithInhomogeneousXDirection(
                    stat_dump_filename,
                    patch_hierarchy,
                    data_context);
        }
        else if (statistical_quantity_key == "R33_MEAN_INHOMO_X")
        {
            rti_rmi_statistics_utilities->
                outputReynoldsNormalStressZWithInhomogeneousXDirection(
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
