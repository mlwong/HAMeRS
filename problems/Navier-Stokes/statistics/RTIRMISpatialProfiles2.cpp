#include "flow/flow_models/four-eqn_conservative/FlowModelStatisticsUtilitiesFourEqnConservative.hpp"

#include "extn/patch_hierarchies/ExtendedFlattenedHierarchy.hpp"
#include "flow/flow_models/MPI_helpers/FlowModelMPIHelperAverage.hpp"
#include "flow/flow_models/MPI_helpers/FlowModelMPIHelperCorrelation.hpp"
#include "flow/flow_models/MPI_helpers/FlowModelMPIHelperMaxMin.hpp"
#include "util/MPI_helpers/MPIHelperGrid.hpp"

#include <fstream>

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
            const HAMERS_SHARED_PTR<EquationOfThermalConductivityMixingRules> equation_of_thermal_conductivity_mixing_rules):
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
         * Output averaged density with assumed homogeneity in y-direction (2D) or yz-plane (3D)
         * to a file.
         */
        void
        outputAveragedDensityWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double output_time) const;
        
        /*
         * Output variance of density with assumed homogeneity in y-direction (2D) or yz-plane (3D)
         * to a file.
         */
        void
        outputDensityVarianceWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double output_time) const;
        
        /*
         * Output averaged mass fraction of species 1 with assumed homogeneity in y-direction (2D) or yz-plane (3D)
         * to a file.
         */
        void
        outputAveragedMassFractionSpeciesOneWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double output_time) const;
        
        /*
         * Output averaged mass fraction of species 2 with assumed homogeneity in y-direction (2D) or yz-plane (3D)
         * to a file.
         */
        void
        outputAveragedMassFractionSpeciesTwoWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double output_time) const;
        
        /*
         * Output variance of mass fraction of species 1 with assumed homogeneity in y-direction (2D) or yz-plane (3D)
         * to a file.
         */
        void
        outputMassFractionVarianceSpeciesOneWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double output_time) const;
        
        /*
         * Output variance of mass fraction of species 2 with assumed homogeneity in y-direction (2D) or yz-plane (3D)
         * to a file.
         */
        void
        outputMassFractionVarianceSpeciesTwoWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double output_time) const;
        
        /*
         * Output covariance of mass fractions of species 1 and 2 with assumed homogeneity in y-direction (2D) or yz-plane (3D)
         * to a file.
         */
        void
        outputMassFractionCovarianceSpeciesOneTwoWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double output_time) const;
        
        /*
         * Output scalar dissipation rate with assumed homogeneity in y-direction (2D) or yz-plane (3D)
         * to a file.
         */
        void
        outputScalarDissipationRateWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double output_time) const;
        
        /*
         * Output TKE dissipation rate with assumed homogeneity in y-direction (2D) or yz-plane (3D)
         * to a file.
         */
        void
        outputTKEDissipationRateWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double output_time) const;
        
    private:
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
 * Output averaged density with assumed homogeneity in y-direction (2D) or yz-plane (3D)
 * to a file.
 */
void
RTIRMISpatialProfilesUtilities::outputAveragedDensityWithHomogeneityInYDirectionOrInYZPlane(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double output_time) const
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
    
    std::ofstream f_out;
    
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
    }
    
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
    
    /*
     * Output the spatial profile (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        f_out.write((char*)&output_time, sizeof(double));
        f_out.write((char*)&rho_avg_global[0], sizeof(double)*rho_avg_global.size());
        
        f_out.close();
    }
}


/*
 * Output variance of density with assumed homogeneity in y-direction (2D) or yz-plane (3D)
 * to a file.
 */
void
RTIRMISpatialProfilesUtilities::outputDensityVarianceWithHomogeneityInYDirectionOrInYZPlane(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
        const double output_time) const
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
    
    std::ofstream f_out;
    
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
    }
    
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
    FlowModelMPIHelperAverage MPI_helper_average = FlowModelMPIHelperAverage(
        "MPI_helper_average",
        d_dim,
        d_grid_geometry,
        patch_hierarchy,
        flow_model_tmp);
    
    FlowModelMPIHelperCorrelation MPI_helper_correlation = FlowModelMPIHelperCorrelation(
        "MPI_helper_correlation",
        d_dim,
        d_grid_geometry,
        patch_hierarchy,
        flow_model_tmp);
    
    std::vector<double> rho_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        data_context);
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    std::vector<std::vector<double> > averaged_quantities;
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(rho_avg_global);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(rho_avg_global);
    
    std::vector<double> rho_var_global = MPI_helper_correlation.getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    averaged_quantities.clear();
    
    /*
     * Output the spatial profile (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        f_out.write((char*)&output_time, sizeof(double));
        f_out.write((char*)&rho_var_global[0], sizeof(double)*rho_var_global.size());
        
        f_out.close();
    }
}


/*
 * Output averaged mass fraction of species 1 with assumed homogeneity in y-direction (2D) or yz-plane (3D)
 * to a file.
 */
void
RTIRMISpatialProfilesUtilities::
outputAveragedMassFractionSpeciesOneWithHomogeneityInYDirectionOrInYZPlane(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "MASS_FRACTION_S1_AVG' can be computed with two species only."
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
    
    std::ofstream f_out;
    
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
    
    /*
     * Output the spatial profile (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        f_out.write((char*)&output_time, sizeof(double));
        f_out.write((char*)&Y_0_avg_global[0], sizeof(double)*Y_0_avg_global.size());
        
        f_out.close();
    }
}


/*
 * Output averaged mass fraction of species 2 with assumed homogeneity in y-direction (2D) or yz-plane (3D)
 * to a file.
 */
void
RTIRMISpatialProfilesUtilities::
outputAveragedMassFractionSpeciesTwoWithHomogeneityInYDirectionOrInYZPlane(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "MASS_FRACTION_S2_AVG' can be computed with two species only."
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
    
    std::ofstream f_out;
    
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
    }
    
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
    FlowModelMPIHelperAverage MPI_helper_average = FlowModelMPIHelperAverage(
        "MPI_helper_average",
        d_dim,
        d_grid_geometry,
        patch_hierarchy,
        flow_model_tmp);
    
    std::vector<double> Y_1_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        "MASS_FRACTIONS",
        1,
        data_context);
    
    /*
     * Output the spatial profile (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        f_out.write((char*)&output_time, sizeof(double));
        f_out.write((char*)&Y_1_avg_global[0], sizeof(double)*Y_1_avg_global.size());
        
        f_out.close();
    }
}


/*
 * Output variance of mass fraction of species 1 with assumed homogeneity in y-direction (2D) or yz-plane (3D)
 * to a file.
 */
void
RTIRMISpatialProfilesUtilities::
outputMassFractionVarianceSpeciesOneWithHomogeneityInYDirectionOrInYZPlane(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'MASS_FRACTION_S1_VAR' can be computed with two species only."
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
    
    std::ofstream f_out;
    
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
    }
    
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
    FlowModelMPIHelperAverage MPI_helper_average = FlowModelMPIHelperAverage(
        "MPI_helper_average",
        d_dim,
        d_grid_geometry,
        patch_hierarchy,
        flow_model_tmp);
    
    FlowModelMPIHelperCorrelation MPI_helper_correlation = FlowModelMPIHelperCorrelation(
        "MPI_helper_correlation",
        d_dim,
        d_grid_geometry,
        patch_hierarchy,
        flow_model_tmp);
    
    std::vector<double> Y_0_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        "MASS_FRACTIONS",
        0,
        data_context);
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    std::vector<std::vector<double> > averaged_quantities;
    
    quantity_names.push_back("MASS_FRACTIONS");
    component_indices.push_back(0);
    averaged_quantities.push_back(Y_0_avg_global);
    
    quantity_names.push_back("MASS_FRACTIONS");
    component_indices.push_back(0);
    averaged_quantities.push_back(Y_0_avg_global);
    
    std::vector<double> Y_0_var_global = MPI_helper_correlation.getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    averaged_quantities.clear();
    
    /*
     * Output the spatial profile (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        f_out.write((char*)&output_time, sizeof(double));
        f_out.write((char*)&Y_0_var_global[0], sizeof(double)*Y_0_var_global.size());
        
        f_out.close();
    }
}


/*
 * Output variance of mass fraction of species 2 with assumed homogeneity in y-direction (2D) or yz-plane (3D)
 * to a file.
 */
void
RTIRMISpatialProfilesUtilities::
outputMassFractionVarianceSpeciesTwoWithHomogeneityInYDirectionOrInYZPlane(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'MASS_FRACTION_S2_VAR' can be computed with two species only."
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
    
    std::ofstream f_out;
    
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
    }
    
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
    FlowModelMPIHelperAverage MPI_helper_average = FlowModelMPIHelperAverage(
        "MPI_helper_average",
        d_dim,
        d_grid_geometry,
        patch_hierarchy,
        flow_model_tmp);
    
    FlowModelMPIHelperCorrelation MPI_helper_correlation = FlowModelMPIHelperCorrelation(
        "MPI_helper_correlation",
        d_dim,
        d_grid_geometry,
        patch_hierarchy,
        flow_model_tmp);
    
    std::vector<double> Y_1_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        "MASS_FRACTIONS",
        1,
        data_context);
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    std::vector<std::vector<double> > averaged_quantities;
    
    quantity_names.push_back("MASS_FRACTIONS");
    component_indices.push_back(1);
    averaged_quantities.push_back(Y_1_avg_global);
    
    quantity_names.push_back("MASS_FRACTIONS");
    component_indices.push_back(1);
    averaged_quantities.push_back(Y_1_avg_global);
    
    std::vector<double> Y_1_var_global = MPI_helper_correlation.getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    averaged_quantities.clear();
    
    /*
     * Output the spatial profile (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        f_out.write((char*)&output_time, sizeof(double));
        f_out.write((char*)&Y_1_var_global[0], sizeof(double)*Y_1_var_global.size());
        
        f_out.close();
    }
}


/*
 * Output covariance of mass fractions of species 1 and 2 with assumed homogeneity in y-direction (2D) or yz-plane (3D)
 * to a file.
 */
void
RTIRMISpatialProfilesUtilities::
outputMassFractionCovarianceSpeciesOneTwoWithHomogeneityInYDirectionOrInYZPlane(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'MASS_FRACTION_S12_COVAR' can be computed with two species only."
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
    
    std::ofstream f_out;
    
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
    }
    
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
    FlowModelMPIHelperAverage MPI_helper_average = FlowModelMPIHelperAverage(
        "MPI_helper_average",
        d_dim,
        d_grid_geometry,
        patch_hierarchy,
        flow_model_tmp);
    
    FlowModelMPIHelperCorrelation MPI_helper_correlation = FlowModelMPIHelperCorrelation(
        "MPI_helper_correlation",
        d_dim,
        d_grid_geometry,
        patch_hierarchy,
        flow_model_tmp);
    
    std::vector<double> Y_0_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        "MASS_FRACTIONS",
        0,
        data_context);
    
    std::vector<double> Y_1_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        "MASS_FRACTIONS",
        0,
        data_context);
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    std::vector<std::vector<double> > averaged_quantities;
    
    quantity_names.push_back("MASS_FRACTIONS");
    component_indices.push_back(0);
    averaged_quantities.push_back(Y_0_avg_global);
    
    quantity_names.push_back("MASS_FRACTIONS");
    component_indices.push_back(1);
    averaged_quantities.push_back(Y_1_avg_global);
    
    std::vector<double> Y_01_covar_global = MPI_helper_correlation.getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    averaged_quantities.clear();
    
    /*
     * Output the spatial profile (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        f_out.write((char*)&output_time, sizeof(double));
        f_out.write((char*)&Y_01_covar_global[0], sizeof(double)*Y_01_covar_global.size());
        
        f_out.close();
    }
}


/*
 * Output scalar dissipation rate with assumed homogeneity in y-direction (2D) or yz-plane (3D)
 * to a file.
 */
void
RTIRMISpatialProfilesUtilities::outputScalarDissipationRateWithHomogeneityInYDirectionOrInYZPlane(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'SCAL_DISS_RATE' can be computed with two species only."
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
    
    std::ofstream f_out;
    
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
    }
    
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
    const bool use_diffusive_flux_utilities = true;
    FlowModelMPIHelperAverage MPI_helper_average = FlowModelMPIHelperAverage(
        "MPI_helper_average",
        d_dim,
        d_grid_geometry,
        patch_hierarchy,
        flow_model_tmp,
        use_diffusive_flux_utilities);
    
    const hier::IntVector& finest_level_dims = MPI_helper_average.getFinestRefinedDomainNumberOfPoints();
    const int num_cells = finest_level_dims[0];
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    std::vector<bool> use_derivative;
    std::vector<int> derivative_directions;
    
    std::vector<double> scalar_dissipation_rate_part_1(num_cells, double(0));
    
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
    
    scalar_dissipation_rate_part_1 = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
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
    
    std::vector<double> scalar_dissipation_rate_part_2(num_cells, double(0));
    if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
    {
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
        
        scalar_dissipation_rate_part_2 = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
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
    }
    
    std::vector<double> scalar_dissipation_rate_part_3(num_cells, double(0));
    if (d_dim == tbox::Dimension(3))
    {
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
        
        scalar_dissipation_rate_part_3 = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
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
    }
    
    if (static_cast<int>(scalar_dissipation_rate_part_1.size()) != num_cells ||
        static_cast<int>(scalar_dissipation_rate_part_2.size()) != num_cells ||
        static_cast<int>(scalar_dissipation_rate_part_3.size()) != num_cells)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Size of scalar dissipation rate parts does not match the number of cells on the finest level!"
            << std::endl);
    }
    
    /*
     * Output the scalar dissipation rate integral (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        std::vector<double> scalar_dissipation_rate_global(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            scalar_dissipation_rate_global[i] =
                scalar_dissipation_rate_part_1[i] + scalar_dissipation_rate_part_2[i] + scalar_dissipation_rate_part_3[i];
        }
        
        f_out.write((char*)&output_time, sizeof(double));
        f_out.write((char*)&scalar_dissipation_rate_global[0], sizeof(double)*scalar_dissipation_rate_global.size());
        
        f_out.close();
    }
}


/*
 * Output TKE dissipation rate with assumed homogeneity in y-direction (2D) or yz-plane (3D)
 * to a file.
 */
void
RTIRMISpatialProfilesUtilities::outputTKEDissipationRateWithHomogeneityInYDirectionOrInYZPlane(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'TKE_DISS_RATE' can be computed with two species only."
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
    
    std::ofstream f_out;
    
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
    }
    
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
    const bool use_diffusive_flux_utilities = true;
    FlowModelMPIHelperAverage MPI_helper_average = FlowModelMPIHelperAverage(
        "MPI_helper_average",
        d_dim,
        d_grid_geometry,
        patch_hierarchy,
        flow_model_tmp,
        use_diffusive_flux_utilities);
    
    const hier::IntVector& finest_level_dims = MPI_helper_average.getFinestRefinedDomainNumberOfPoints();
    const int num_cells = finest_level_dims[0];
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    std::vector<bool> use_derivative;
    std::vector<int> derivative_directions;
    std::vector<bool> use_reciprocal;
    
    std::vector<double> ddx_u_avg;
    std::vector<double> ddx_v_avg;
    std::vector<double> ddx_w_avg;
    std::vector<double> ddy_u_avg;
    std::vector<double> ddy_v_avg;
    std::vector<double> ddy_w_avg;
    std::vector<double> ddz_u_avg;
    std::vector<double> ddz_v_avg;
    std::vector<double> ddz_w_avg;
    
    // Computed ddx_u_avg.
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    use_reciprocal.push_back(false);
    
    ddx_u_avg = MPI_helper_average.getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        use_reciprocal,
        0,
        d_num_ghosts_derivative,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    use_reciprocal.clear();
    
    if (static_cast<int>(ddx_u_avg.size()) != num_cells)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Size of ddx_u_avg does not match the number of cells on the finest level!"
            << std::endl);
    }
    
    if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
    {
        // Computed ddx_v_avg.
        
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(1);
        use_reciprocal.push_back(false);
        
        ddx_v_avg = MPI_helper_average.getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
            quantity_names,
            component_indices,
            use_reciprocal,
            0,
            d_num_ghosts_derivative,
            data_context);
        
        quantity_names.clear();
        component_indices.clear();
        use_reciprocal.clear();
        
        // Compute ddy_u_avg.
        
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(0);
        use_reciprocal.push_back(false);
        
        ddy_u_avg = MPI_helper_average.getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
            quantity_names,
            component_indices,
            use_reciprocal,
            1,
            d_num_ghosts_derivative,
            data_context);
        
        quantity_names.clear();
        component_indices.clear();
        use_reciprocal.clear();
        
        // Compute ddy_v_avg.
        
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(1);
        use_reciprocal.push_back(false);
    
        ddy_v_avg = MPI_helper_average.getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
            quantity_names,
            component_indices,
            use_reciprocal,
            1,
            d_num_ghosts_derivative,
            data_context);
        
        quantity_names.clear();
        component_indices.clear();
        use_reciprocal.clear();
        
        if (static_cast<int>(ddx_v_avg.size()) != num_cells ||
            static_cast<int>(ddy_u_avg.size()) != num_cells ||
            static_cast<int>(ddy_v_avg.size()) != num_cells)
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Size of ddx_v_avg, ddy_u_avg, or ddy_v_avg does not match the number of cells on the finest level!"
                << std::endl);
        }
    }
    
    if (d_dim == tbox::Dimension(3))
    {
        // Computed ddx_w_avg.
        
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(2);
        use_reciprocal.push_back(false);
        
        ddx_w_avg = MPI_helper_average.getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
            quantity_names,
            component_indices,
            use_reciprocal,
            0,
            d_num_ghosts_derivative,
            data_context);
        
        quantity_names.clear();
        component_indices.clear();
        use_reciprocal.clear();
        
        // Compute ddy_w_avg.
        
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(2);
        use_reciprocal.push_back(false);
    
        ddy_w_avg = MPI_helper_average.getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
            quantity_names,
            component_indices,
            use_reciprocal,
            1,
            d_num_ghosts_derivative,
            data_context);
        
        quantity_names.clear();
        component_indices.clear();
        use_reciprocal.clear();
        
        // Compute ddz_u_avg.
        
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(0);
        use_reciprocal.push_back(false);
        
        ddz_u_avg = MPI_helper_average.getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
            quantity_names,
            component_indices,
            use_reciprocal,
            2,
            d_num_ghosts_derivative,
            data_context);
        
        quantity_names.clear();
        component_indices.clear();
        use_reciprocal.clear();
        
        // Compute ddz_v_avg.
        
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(1);
        use_reciprocal.push_back(false);
        
        ddz_v_avg = MPI_helper_average.getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
            quantity_names,
            component_indices,
            use_reciprocal,
            2,
            d_num_ghosts_derivative,
            data_context);
        
        quantity_names.clear();
        component_indices.clear();
        use_reciprocal.clear();
        
        // Compute ddz_w_avg.
        
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(2);
        use_reciprocal.push_back(false);
    
        ddz_w_avg = MPI_helper_average.getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
            quantity_names,
            component_indices,
            use_reciprocal,
            2,
            d_num_ghosts_derivative,
            data_context);
        
        quantity_names.clear();
        component_indices.clear();
        use_reciprocal.clear();
        
        if (static_cast<int>(ddx_w_avg.size()) != num_cells ||
            static_cast<int>(ddy_w_avg.size()) != num_cells ||
            static_cast<int>(ddz_u_avg.size()) != num_cells ||
            static_cast<int>(ddz_v_avg.size()) != num_cells ||
            static_cast<int>(ddz_w_avg.size()) != num_cells)
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Size of ddx_w_avg, ddy_w_avg, ddz_u_avg, ddz_v_avg, or ddz_w_avg does not match the number of cells on the finest level!"
                << std::endl);
        }
    }
    
    std::vector<double> tau11_avg;
    std::vector<double> tau12_avg;
    std::vector<double> tau13_avg;
    std::vector<double> tau22_avg;
    std::vector<double> tau23_avg;
    std::vector<double> tau33_avg;
    
    // Compute tau11_avg.
    
    tau11_avg = getAveragedShearStressComponentWithInhomogeneousXDirection(
        0,
        patch_hierarchy,
        data_context);
    
    if (static_cast<int>(tau11_avg.size()) != num_cells)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Size of tau11_avg does not match the number of cells on the finest level!"
            << std::endl);
    }
    
    if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
    {
        // Compute tau12_avg.
        
        tau12_avg = getAveragedShearStressComponentWithInhomogeneousXDirection(
            1,
            patch_hierarchy,
            data_context);
        
        // Compute tau22_avg.
        
        tau22_avg = getAveragedShearStressComponentWithInhomogeneousXDirection(
            3,
            patch_hierarchy,
            data_context);
        
        if (static_cast<int>(tau12_avg.size()) != num_cells ||
            static_cast<int>(tau22_avg.size()) != num_cells)
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Size of tau12_avg or tau22_avg does not match the number of cells on the finest level!"
                << std::endl);
        }
    }
    
    if (d_dim == tbox::Dimension(3))
    {
        // Compute tau13_avg.
        
        tau13_avg = getAveragedShearStressComponentWithInhomogeneousXDirection(
            2,
            patch_hierarchy,
            data_context);
        
        // Compute tau23_avg.
        
        tau23_avg = getAveragedShearStressComponentWithInhomogeneousXDirection(
            4,
            patch_hierarchy,
            data_context);
        
        // Compute tau33_avg.
        
        tau33_avg = getAveragedShearStressComponentWithInhomogeneousXDirection(
            5,
            patch_hierarchy,
            data_context);
        
        if (static_cast<int>(tau13_avg.size()) != num_cells ||
            static_cast<int>(tau23_avg.size()) != num_cells ||
            static_cast<int>(tau33_avg.size()) != num_cells)
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Size of tau13_avg, tau23_avg, or tau33_avg does not match the number of cells on the finest level!"
                << std::endl);
        }
    }
    
    std::vector<double> tau11_ddx_u_avg;
    std::vector<double> tau12_ddy_u_avg;
    std::vector<double> tau13_ddz_u_avg;
    std::vector<double> tau12_ddx_v_avg;
    std::vector<double> tau22_ddy_v_avg;
    std::vector<double> tau23_ddz_v_avg;
    std::vector<double> tau13_ddx_w_avg;
    std::vector<double> tau23_ddy_w_avg;
    std::vector<double> tau33_ddz_w_avg;
    
    // Compute tau11_ddx_u_avg.
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    use_derivative.push_back(true);
    derivative_directions.push_back(0);
    use_reciprocal.push_back(false);
    
    tau11_ddx_u_avg = getAveragedQuantityWithShearStressComponentWithInhomogeneousXDirection(
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
    
    if (static_cast<int>(tau11_ddx_u_avg.size()) != num_cells)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Size of tau11_ddx_u_avg does not match the number of cells on the finest level!"
            << std::endl);
    }
    
    if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
    {
        // Compute tau12_ddy_u_avg.
        
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(0);
        use_derivative.push_back(true);
        derivative_directions.push_back(1);
        use_reciprocal.push_back(false);
        
        tau12_ddy_u_avg = getAveragedQuantityWithShearStressComponentWithInhomogeneousXDirection(
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
        
        // Compute tau12_ddx_v_avg.
        
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(1);
        use_derivative.push_back(true);
        derivative_directions.push_back(0);
        use_reciprocal.push_back(false);
        
        tau12_ddx_v_avg = getAveragedQuantityWithShearStressComponentWithInhomogeneousXDirection(
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
        
        // Compute tau22_ddy_v_avg.
        
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(1);
        use_derivative.push_back(true);
        derivative_directions.push_back(1);
        use_reciprocal.push_back(false);
        
        tau22_ddy_v_avg = getAveragedQuantityWithShearStressComponentWithInhomogeneousXDirection(
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
        
        if (static_cast<int>(tau12_ddy_u_avg.size()) != num_cells ||
            static_cast<int>(tau12_ddx_v_avg.size()) != num_cells ||
            static_cast<int>(tau22_ddy_v_avg.size()) != num_cells)
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Size of tau12_ddy_u_avg, tau12_ddx_v_avg, or tau22_ddy_v_avg does not match the number of cells on the finest level!"
                << std::endl);
        }
    }
    
    if (d_dim == tbox::Dimension(3))
    {
        // Compute tau13_ddz_u_avg.
        
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(0);
        use_derivative.push_back(true);
        derivative_directions.push_back(2);
        use_reciprocal.push_back(false);
        
        tau13_ddz_u_avg = getAveragedQuantityWithShearStressComponentWithInhomogeneousXDirection(
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
        
        // Compute tau23_ddz_v_avg.
        
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(1);
        use_derivative.push_back(true);
        derivative_directions.push_back(2);
        use_reciprocal.push_back(false);
        
        tau23_ddz_v_avg = getAveragedQuantityWithShearStressComponentWithInhomogeneousXDirection(
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
        
        // Compute tau13_ddx_w_avg.
        
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(2);
        use_derivative.push_back(true);
        derivative_directions.push_back(0);
        use_reciprocal.push_back(false);
        
        tau13_ddx_w_avg = getAveragedQuantityWithShearStressComponentWithInhomogeneousXDirection(
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
        
        // Compute tau23_ddy_w_avg.
        
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(2);
        use_derivative.push_back(true);
        derivative_directions.push_back(1);
        use_reciprocal.push_back(false);
        
        tau23_ddy_w_avg = getAveragedQuantityWithShearStressComponentWithInhomogeneousXDirection(
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
        
        // Compute tau33_ddz_w_avg.
        
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(2);
        use_derivative.push_back(true);
        derivative_directions.push_back(2);
        use_reciprocal.push_back(false);
        
        tau33_ddz_w_avg = getAveragedQuantityWithShearStressComponentWithInhomogeneousXDirection(
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
        
        if (static_cast<int>(tau13_ddz_u_avg.size()) != num_cells ||
            static_cast<int>(tau23_ddz_v_avg.size()) != num_cells ||
            static_cast<int>(tau13_ddx_w_avg.size()) != num_cells ||
            static_cast<int>(tau23_ddy_w_avg.size()) != num_cells ||
            static_cast<int>(tau33_ddz_w_avg.size()) != num_cells)
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Size of tau13_ddz_u_avg, tau23_ddz_v_avg, tau13_ddx_w_avg, tau23_ddy_w_avg, or tau33_ddz_w_avg does not match the number of cells on the finest level!"
                << std::endl);
        }
    }
    
    /*
     * Compute the dissipation rate of R11.
     */
    
    std::vector<double> m_2tau1i_p_ddxi_u_p(num_cells, double(0));
    
    std::vector<double> tau11_p_ddx_u_p(num_cells, double(0));
    for (int i = 0; i < num_cells; i++)
    {
        tau11_p_ddx_u_p[i] = tau11_ddx_u_avg[i] - tau11_avg[i]*ddx_u_avg[i];
    }
    
    std::vector<double> tau12_p_ddy_u_p(num_cells, double(0));
    if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
    {
        for (int i = 0; i < num_cells; i++)
        {
            tau12_p_ddy_u_p[i] = tau12_ddy_u_avg[i] - tau12_avg[i]*ddy_u_avg[i];
        }
    }
    
    std::vector<double> tau13_p_ddz_u_p(num_cells, double(0));
    if (d_dim == tbox::Dimension(3))
    {
        for (int i = 0; i < num_cells; i++)
        {
            tau13_p_ddz_u_p[i] = tau13_ddz_u_avg[i] - tau13_avg[i]*ddz_u_avg[i];
        }
    }
    
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
     * Compute the dissipation rate of R22.
     */
    
    std::vector<double> m_2tau_2i_p_ddxi_v_p(num_cells, double(0));
    
    if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
    {
        std::vector<double> tau12_p_ddx_v_p(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            tau12_p_ddx_v_p[i] = tau12_ddx_v_avg[i] - tau12_avg[i]*ddx_v_avg[i];
        }
        
        std::vector<double> tau22_p_ddy_v_p(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            tau22_p_ddy_v_p[i] = tau22_ddy_v_avg[i] - tau22_avg[i]*ddy_v_avg[i];
        }
        
        std::vector<double> tau23_p_ddz_v_p(num_cells, double(0));
        if (d_dim == tbox::Dimension(3))
        {
            for (int i = 0; i < num_cells; i++)
            {
                tau23_p_ddz_v_p[i] = tau23_ddz_v_avg[i] - tau23_avg[i]*ddz_v_avg[i];
            }
        }
        
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
    }
    
    /*
     * Compute the dissipation rate of R33.
     */
    
    std::vector<double> m_2tau_3i_p_ddxi_w_p(num_cells, double(0));
    
    if (d_dim == tbox::Dimension(3))
    {
        std::vector<double> tau13_p_ddx_w_p(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            tau13_p_ddx_w_p[i] = tau13_ddx_w_avg[i] - tau13_avg[i]*ddx_w_avg[i];
        }
        
        std::vector<double> tau23_p_ddy_w_p(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            tau23_p_ddy_w_p[i] = tau23_ddy_w_avg[i] - tau23_avg[i]*ddy_w_avg[i];
        }
        
        std::vector<double> tau33_p_ddz_w_p(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            tau33_p_ddz_w_p[i] = tau33_ddz_w_avg[i] - tau33_avg[i]*ddz_w_avg[i];
        }
    
        for (int i = 0; i < num_cells; i++)
        {
            m_2tau_3i_p_ddxi_w_p[i] = -double(2)*(tau13_p_ddx_w_p[i] + tau23_p_ddy_w_p[i] + tau33_p_ddz_w_p[i]);
        }
    }
    
    /*
     * Output the spatial profile (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        std::vector<double> TKE_dissipation_rate_global(num_cells, double(0));
        for (int i = 0; i < num_cells; i++)
        {
            TKE_dissipation_rate_global[i] =
                -0.5*(m_2tau1i_p_ddxi_u_p[i] + m_2tau_2i_p_ddxi_v_p[i] + m_2tau_3i_p_ddxi_w_p[i]);
        }
        
        f_out.write((char*)&output_time, sizeof(double));
        f_out.write((char*)&TKE_dissipation_rate_global[0], sizeof(double)*TKE_dissipation_rate_global.size());
        
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
RTIRMISpatialProfilesUtilities::getAveragedShearStressComponentWithInhomogeneousXDirection(
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
RTIRMISpatialProfilesUtilities::getAveragedQuantityWithShearStressComponentWithInhomogeneousXDirection(
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
                        << ": RTIRMIStatisticsUtilities::"
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
 * Compute statisitcal quantities.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::computeStatisticalQuantities(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double statistics_data_time)
{
    NULL_USE(patch_hierarchy);
    NULL_USE(data_context);
    NULL_USE(statistics_data_time);
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
            d_equation_of_thermal_conductivity_mixing_rules));
    
    // Loop over statistical quantities.
    for (int qi = 0; qi < static_cast<int>(d_statistical_quantities.size()); qi++)
    {
        // Get the key of the current variable.
        std::string statistical_quantity_key = d_statistical_quantities[qi];
        
        if (statistical_quantity_key == "DENSITY_AVG")
        {
            rti_rmi_spatial_profiles_utilities->outputAveragedDensityWithHomogeneityInYDirectionOrInYZPlane(
                "rho_avg.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "DENSITY_VAR")
        {
            rti_rmi_spatial_profiles_utilities->outputDensityVarianceWithHomogeneityInYDirectionOrInYZPlane(
                "rho_var.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "MASS_FRACTION_S1_AVG")
        {
            rti_rmi_spatial_profiles_utilities->outputAveragedMassFractionSpeciesOneWithHomogeneityInYDirectionOrInYZPlane(
                "Y1_avg.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "MASS_FRACTION_S2_AVG")
        {
            rti_rmi_spatial_profiles_utilities->outputAveragedMassFractionSpeciesTwoWithHomogeneityInYDirectionOrInYZPlane(
                "Y2_avg.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "MASS_FRACTION_S1_VAR")
        {
            rti_rmi_spatial_profiles_utilities->outputMassFractionVarianceSpeciesOneWithHomogeneityInYDirectionOrInYZPlane(
                "Y1_var.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "MASS_FRACTION_S2_VAR")
        {
            rti_rmi_spatial_profiles_utilities->outputMassFractionVarianceSpeciesTwoWithHomogeneityInYDirectionOrInYZPlane(
                "Y2_var.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "MASS_FRACTION_S12_COVAR")
        {
            rti_rmi_spatial_profiles_utilities->outputMassFractionCovarianceSpeciesOneTwoWithHomogeneityInYDirectionOrInYZPlane(
                "Y1Y2_covar.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "SCAL_DISS_RATE")
        {
            rti_rmi_spatial_profiles_utilities->outputScalarDissipationRateWithHomogeneityInYDirectionOrInYZPlane(
                "chi_avg.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "TKE_DISS_RATE")
        {
            rti_rmi_spatial_profiles_utilities->outputTKEDissipationRateWithHomogeneityInYDirectionOrInYZPlane(
                "epsilon_avg.dat",
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
