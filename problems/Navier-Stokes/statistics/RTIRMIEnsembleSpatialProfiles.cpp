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
            rho_avg_computed      = false;
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
            rho_avg_realizations.clear();
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
            rho_avg_computed      = false;
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
        std::vector<std::vector<double> > rho_avg_realizations;
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
        bool rho_avg_computed;
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

class RTIRMISpatialProfilesUtilities
{
    public:
        RTIRMISpatialProfilesUtilities(
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
         * Compute averaged density with assumed homogeneity in y-direction (2D) or yz-plane (3D).
         */
        void
        computeAveragedDensityWithHomogeneityInYDirectionOrInYZPlane(
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
RTIRMISpatialProfilesUtilities::computeAveragedMassFractionWithHomogeneityInYDirectionOrInYZPlane(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'MASS_FRACTION_AVG' can be computed with two species only."
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
RTIRMISpatialProfilesUtilities::computeMassFractionVarianceWithHomogeneityInYDirectionOrInYZPlane(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'MASS_FRACTION_VAR' can be computed with two species only."
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
 * Compute averaged density with assumed homogeneity in y-direction (2D) or yz-plane (3D).
 */
void
RTIRMISpatialProfilesUtilities::computeAveragedDensityWithHomogeneityInYDirectionOrInYZPlane(
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
 * Compute averaged velocity x-component with assumed homogeneity in y-direction (2D) or yz-plane (3D).
 */
void
RTIRMISpatialProfilesUtilities::computeAveragedVelocityXWithHomogeneityInYDirectionOrInYZPlane(
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
RTIRMISpatialProfilesUtilities::computeAveragedVelocityYWithHomogeneityInYDirectionOrInYZPlane(
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
RTIRMISpatialProfilesUtilities::computeAveragedVelocityZWithHomogeneityInYZPlane(
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
RTIRMISpatialProfilesUtilities::computeAveragedMomentumXWithHomogeneityInYDirectionOrInYZPlane(
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
RTIRMISpatialProfilesUtilities::computeAveragedMomentumYWithHomogeneityInYDirectionOrInYZPlane(
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
RTIRMISpatialProfilesUtilities::computeAveragedMomentumZWithHomogeneityInYZPlane(
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
RTIRMISpatialProfilesUtilities::computeReynoldsNormalStressXWithHomogeneityInYDirectionOrInYZPlane(
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
RTIRMISpatialProfilesUtilities::computeReynoldsNormalStressYWithHomogeneityInYDirectionOrInYZPlane(
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
RTIRMISpatialProfilesUtilities::computeReynoldsNormalStressZWithHomogeneityInYZPlane(
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
RTIRMISpatialProfilesUtilities::outputEnsembleAveragedMassFractionWithHomogeneityInYDirectionOrInYZPlane(
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
RTIRMISpatialProfilesUtilities::outputEnsembleMassFractionVarianceWithHomogeneityInYDirectionOrInYZPlane(
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
RTIRMISpatialProfilesUtilities::outputEnsembleAveragedDensityWithHomogeneityInYDirectionOrInYZPlane(
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
 * Output ensemble averaged velocity x-component with assumed homogeneity in y-direction (2D) or yz-plane (3D)
 * to a file.
 */
void
RTIRMISpatialProfilesUtilities::outputEnsembleAveragedVelocityXWithHomogeneityInYDirectionOrInYZPlane(
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
RTIRMISpatialProfilesUtilities::outputEnsembleAveragedVelocityYWithHomogeneityInYDirectionOrInYZPlane(
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
RTIRMISpatialProfilesUtilities::outputEnsembleAveragedVelocityZWithHomogeneityInYZPlane(
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
RTIRMISpatialProfilesUtilities::outputEnsembleAveragedMomentumXWithHomogeneityInYDirectionOrInYZPlane(
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
RTIRMISpatialProfilesUtilities::outputEnsembleAveragedMomentumYWithHomogeneityInYDirectionOrInYZPlane(
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
RTIRMISpatialProfilesUtilities::outputEnsembleAveragedMomentumZWithHomogeneitynYZPlane(
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
RTIRMISpatialProfilesUtilities::outputEnsembleTurbulentMassFluxXWithHomogeneityInYDirectionOrInYZPlane(
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
RTIRMISpatialProfilesUtilities::outputEnsembleReynoldsNormalStressXWithHomogeneityInYDirectionOrInYZPlane(
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
RTIRMISpatialProfilesUtilities::outputEnsembleReynoldsNormalStressYWithHomogeneityInYDirectionOrInYZPlane(
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
RTIRMISpatialProfilesUtilities::outputEnsembleReynoldsNormalStressZWithHomogeneityInYZPlane(
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
        // Loop over statistical quantities.
        for (int qi = 0; qi < static_cast<int>(d_statistical_quantities.size()); qi++)
        {
            // DO NOTHING.
        }
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
    
    HAMERS_SHARED_PTR<RTIRMISpatialProfilesUtilities> rti_rmi_spatial_profiles_utilities(
        new RTIRMISpatialProfilesUtilities(
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
    rti_rmi_spatial_profiles_utilities->d_ensemble_statistics->Y_0_avg_computed      = false;
    rti_rmi_spatial_profiles_utilities->d_ensemble_statistics->Y_0_sq_avg_computed   = false;
    rti_rmi_spatial_profiles_utilities->d_ensemble_statistics->rho_avg_computed      = false;
    rti_rmi_spatial_profiles_utilities->d_ensemble_statistics->u_avg_computed        = false;
    rti_rmi_spatial_profiles_utilities->d_ensemble_statistics->v_avg_computed        = false;
    rti_rmi_spatial_profiles_utilities->d_ensemble_statistics->w_avg_computed        = false;
    rti_rmi_spatial_profiles_utilities->d_ensemble_statistics->rho_u_avg_computed    = false;
    rti_rmi_spatial_profiles_utilities->d_ensemble_statistics->rho_v_avg_computed    = false;
    rti_rmi_spatial_profiles_utilities->d_ensemble_statistics->rho_w_avg_computed    = false;
    rti_rmi_spatial_profiles_utilities->d_ensemble_statistics->rho_u_u_avg_computed  = false;
    rti_rmi_spatial_profiles_utilities->d_ensemble_statistics->rho_v_v_avg_computed  = false;
    rti_rmi_spatial_profiles_utilities->d_ensemble_statistics->rho_w_w_avg_computed  = false;
    
    // Loop over statistical quantities.
    for (int qi = 0; qi < static_cast<int>(d_statistical_quantities.size()); qi++)
    {
        // Get the key of the current variable.
        std::string statistical_quantity_key = d_statistical_quantities[qi];
        
        if (statistical_quantity_key == "MASS_FRACTION_AVG")
        {
            if (!(rti_rmi_spatial_profiles_utilities->d_ensemble_statistics->Y_0_avg_computed))
            {
                rti_rmi_spatial_profiles_utilities->
                    computeAveragedMassFractionWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
        }
        else if (statistical_quantity_key == "MASS_FRACTION_VAR")
        {
            if (!(rti_rmi_spatial_profiles_utilities->d_ensemble_statistics->Y_0_avg_computed))
            {
                rti_rmi_spatial_profiles_utilities->
                    computeAveragedMassFractionWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
            
            if (!(rti_rmi_spatial_profiles_utilities->d_ensemble_statistics->Y_0_sq_avg_computed))
            {
                rti_rmi_spatial_profiles_utilities->
                    computeMassFractionVarianceWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
        }
        else if (statistical_quantity_key == "DENSITY_AVG")
        {
            rti_rmi_spatial_profiles_utilities->
                computeAveragedDensityWithHomogeneityInYDirectionOrInYZPlane(
                    patch_hierarchy,
                    data_context);
        }
        else if (statistical_quantity_key == "VELOCITY_X_AVG")
        {
            rti_rmi_spatial_profiles_utilities->
                computeAveragedVelocityXWithHomogeneityInYDirectionOrInYZPlane(
                    patch_hierarchy,
                    data_context);
        }
        else if (statistical_quantity_key == "VELOCITY_Y_AVG")
        {
            rti_rmi_spatial_profiles_utilities->
                computeAveragedVelocityYWithHomogeneityInYDirectionOrInYZPlane(
                    patch_hierarchy,
                    data_context);
        }
        else if (statistical_quantity_key == "VELOCITY_Z_AVG")
        {
            rti_rmi_spatial_profiles_utilities->
                computeAveragedVelocityZWithHomogeneityInYZPlane(
                    patch_hierarchy,
                    data_context);
        }
        else if (statistical_quantity_key == "MOMENTUM_X_AVG")
        {
            rti_rmi_spatial_profiles_utilities->
                computeAveragedMomentumXWithHomogeneityInYDirectionOrInYZPlane(
                    patch_hierarchy,
                    data_context);
        }
        else if (statistical_quantity_key == "MOMENTUM_Y_AVG")
        {
            rti_rmi_spatial_profiles_utilities->
                computeAveragedMomentumYWithHomogeneityInYDirectionOrInYZPlane(
                    patch_hierarchy,
                    data_context);
        }
        else if (statistical_quantity_key == "MOMENTUM_Z_AVG")
        {
            rti_rmi_spatial_profiles_utilities->
                computeAveragedMomentumZWithHomogeneityInYZPlane(
                    patch_hierarchy,
                    data_context);
        }
        else if (statistical_quantity_key == "TURB_MASS_FLUX_X")
        {
            if (!(rti_rmi_spatial_profiles_utilities->d_ensemble_statistics->rho_avg_computed))
            {
                rti_rmi_spatial_profiles_utilities->
                    computeAveragedDensityWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
            
            if (!(rti_rmi_spatial_profiles_utilities->d_ensemble_statistics->u_avg_computed))
            {
                rti_rmi_spatial_profiles_utilities->
                    computeAveragedVelocityXWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
            
            if (!(rti_rmi_spatial_profiles_utilities->d_ensemble_statistics->rho_u_avg_computed))
            {
                rti_rmi_spatial_profiles_utilities->
                    computeAveragedMomentumXWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
        }
        else if (statistical_quantity_key == "RE_STRESS_11")
        {
            if (!(rti_rmi_spatial_profiles_utilities->d_ensemble_statistics->rho_avg_computed))
            {
                rti_rmi_spatial_profiles_utilities->
                    computeAveragedDensityWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
            
            if (!(rti_rmi_spatial_profiles_utilities->d_ensemble_statistics->rho_u_avg_computed))
            {
                rti_rmi_spatial_profiles_utilities->
                    computeAveragedMomentumXWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
            
            if (!(rti_rmi_spatial_profiles_utilities->d_ensemble_statistics->rho_u_u_avg_computed))
            {
                rti_rmi_spatial_profiles_utilities->
                    computeReynoldsNormalStressXWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
        }
        else if (statistical_quantity_key == "RE_STRESS_22")
        {
            if (!(rti_rmi_spatial_profiles_utilities->d_ensemble_statistics->rho_avg_computed))
            {
                rti_rmi_spatial_profiles_utilities->
                    computeAveragedDensityWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
            
            if (!(rti_rmi_spatial_profiles_utilities->d_ensemble_statistics->rho_v_avg_computed))
            {
                rti_rmi_spatial_profiles_utilities->
                    computeAveragedMomentumYWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
            
            if (!(rti_rmi_spatial_profiles_utilities->d_ensemble_statistics->rho_v_v_avg_computed))
            {
                rti_rmi_spatial_profiles_utilities->
                    computeReynoldsNormalStressYWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
        }
        else if (statistical_quantity_key == "RE_STRESS_33")
        {
            if (!(rti_rmi_spatial_profiles_utilities->d_ensemble_statistics->rho_avg_computed))
            {
                rti_rmi_spatial_profiles_utilities->
                    computeAveragedDensityWithHomogeneityInYDirectionOrInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
            
            if (!(rti_rmi_spatial_profiles_utilities->d_ensemble_statistics->rho_w_avg_computed))
            {
                rti_rmi_spatial_profiles_utilities->
                    computeAveragedMomentumZWithHomogeneityInYZPlane(
                        patch_hierarchy,
                        data_context);
            }
            
            if (!(rti_rmi_spatial_profiles_utilities->d_ensemble_statistics->rho_w_w_avg_computed))
            {
                rti_rmi_spatial_profiles_utilities->
                    computeReynoldsNormalStressZWithHomogeneityInYZPlane(
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
    
    HAMERS_SHARED_PTR<RTIRMISpatialProfilesUtilities> rti_rmi_spatial_profiles_utilities(
        new RTIRMISpatialProfilesUtilities(
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
        
        if (statistical_quantity_key == "MASS_FRACTION_AVG")
        {
            rti_rmi_spatial_profiles_utilities->
                outputEnsembleAveragedMassFractionWithHomogeneityInYDirectionOrInYZPlane(
                    "Y_avg.dat",
                    patch_hierarchy,
                    data_context,
                    output_time);
        }
        else if (statistical_quantity_key == "MASS_FRACTION_VAR")
        {
            rti_rmi_spatial_profiles_utilities->
                outputEnsembleMassFractionVarianceWithHomogeneityInYDirectionOrInYZPlane(
                    "Y_var.dat",
                    patch_hierarchy,
                    data_context,
                    output_time);
        }
        else if (statistical_quantity_key == "DENSITY_AVG")
        {
            rti_rmi_spatial_profiles_utilities->
                outputEnsembleAveragedDensityWithHomogeneityInYDirectionOrInYZPlane(
                    "rho_avg.dat",
                    patch_hierarchy,
                    data_context,
                    output_time);
        }
        else if (statistical_quantity_key == "VELOCITY_X_AVG")
        {
            rti_rmi_spatial_profiles_utilities->
                outputEnsembleAveragedVelocityXWithHomogeneityInYDirectionOrInYZPlane(
                    "u_avg.dat",
                    patch_hierarchy,
                    data_context,
                    output_time);
        }
        else if (statistical_quantity_key == "VELOCITY_Y_AVG")
        {
            rti_rmi_spatial_profiles_utilities->
                outputEnsembleAveragedVelocityYWithHomogeneityInYDirectionOrInYZPlane(
                    "v_avg.dat",
                    patch_hierarchy,
                    data_context,
                    output_time);
        }
        else if (statistical_quantity_key == "VELOCITY_Z_AVG")
        {
            rti_rmi_spatial_profiles_utilities->
                outputEnsembleAveragedVelocityZWithHomogeneityInYZPlane(
                    "w_avg.dat",
                    patch_hierarchy,
                    data_context,
                    output_time);
        }
        else if (statistical_quantity_key == "MOMENTUM_X_AVG")
        {
            rti_rmi_spatial_profiles_utilities->
                outputEnsembleAveragedMomentumXWithHomogeneityInYDirectionOrInYZPlane(
                    "rho_u_avg.dat",
                    patch_hierarchy,
                    data_context,
                    output_time);
        }
        else if (statistical_quantity_key == "MOMENTUM_Y_AVG")
        {
            rti_rmi_spatial_profiles_utilities->
                outputEnsembleAveragedMomentumYWithHomogeneityInYDirectionOrInYZPlane(
                    "rho_v_avg.dat",
                    patch_hierarchy,
                    data_context,
                    output_time);
        }
        else if (statistical_quantity_key == "MOMENTUM_Z_AVG")
        {
            rti_rmi_spatial_profiles_utilities->
                outputEnsembleAveragedMomentumZWithHomogeneitynYZPlane(
                    "rho_w_avg.dat",
                    patch_hierarchy,
                    data_context,
                    output_time);
        }
        else if (statistical_quantity_key == "TURB_MASS_FLUX_X")
        {
            rti_rmi_spatial_profiles_utilities->
                outputEnsembleTurbulentMassFluxXWithHomogeneityInYDirectionOrInYZPlane(
                    "rho_a1.dat",
                    patch_hierarchy,
                    data_context,
                    output_time);
        }
        else if (statistical_quantity_key == "RE_STRESS_11")
        {
            rti_rmi_spatial_profiles_utilities->
                outputEnsembleReynoldsNormalStressXWithHomogeneityInYDirectionOrInYZPlane(
                    "R11.dat",
                    patch_hierarchy,
                    data_context,
                    output_time);
        }
        else if (statistical_quantity_key == "RE_STRESS_22")
        {
            rti_rmi_spatial_profiles_utilities->
                outputEnsembleReynoldsNormalStressYWithHomogeneityInYDirectionOrInYZPlane(
                    "R22.dat",
                    patch_hierarchy,
                    data_context,
                    output_time);
        }
        else if (statistical_quantity_key == "RE_STRESS_33")
        {
            rti_rmi_spatial_profiles_utilities->
                outputEnsembleReynoldsNormalStressZWithHomogeneityInYZPlane(
                    "R33.dat",
                    patch_hierarchy,
                    data_context,
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
