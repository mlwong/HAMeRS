#include "flow/flow_models/four-eqn_conservative/FlowModelStatisticsUtilitiesFourEqnConservative.hpp"

#include "extn/patch_hierarchies/ExtendedFlattenedHierarchy.hpp"
#include "flow/flow_models/MPI_helpers/FlowModelMPIHelperAverage.hpp"
#include "flow/flow_models/MPI_helpers/FlowModelMPIHelperCorrelation.hpp"
#include "flow/flow_models/MPI_helpers/FlowModelMPIHelperMaxMin.hpp"
#include "util/MPI_helpers/MPIHelperGrid.hpp"

#include <fstream>

class RTIThreeLayersSpatialProfilesUtilities
{
    public:
        RTIThreeLayersSpatialProfilesUtilities(
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
         * Output averaged mass fraction of species 3 with assumed homogeneity in y-direction (2D) or yz-plane (3D)
         * to a file.
         */
        void
        outputAveragedMassFractionSpeciesThreeWithHomogeneityInYDirectionOrInYZPlane(
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
         * Output variance of mass fraction of species 3 with assumed homogeneity in y-direction (2D) or yz-plane (3D)
         * to a file.
         */
        void
        outputMassFractionVarianceSpeciesThreeWithHomogeneityInYDirectionOrInYZPlane(
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
         * Output covariance of mass fractions of species 2 and 3 with assumed homogeneity in y-direction (2D) or yz-plane (3D)
         * to a file.
         */
        void
        outputMassFractionCovarianceSpeciesTwoThreeWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double output_time) const;
        
        /*
         * Output covariance of mass fractions of species 1 and 3 with assumed homogeneity in y-direction (2D) or yz-plane (3D)
         * to a file.
         */
        void
        outputMassFractionCovarianceSpeciesOneThreeWithHomogeneityInYDirectionOrInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double output_time) const;
        
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
 * Output averaged mass fraction of species 1 with assumed homogeneity in y-direction (2D) or yz-plane (3D)
 * to a file.
 */
void
RTIThreeLayersSpatialProfilesUtilities::
outputAveragedMassFractionSpeciesOneWithHomogeneityInYDirectionOrInYZPlane(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 3)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "MASS_FRACTION_S1_AVG' can be computed with three species only."
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
RTIThreeLayersSpatialProfilesUtilities::
outputAveragedMassFractionSpeciesTwoWithHomogeneityInYDirectionOrInYZPlane(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 3)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "MASS_FRACTION_S2_AVG' can be computed with three species only."
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
 * Output averaged mass fraction of species 3 with assumed homogeneity in y-direction (2D) or yz-plane (3D)
 * to a file.
 */
void
RTIThreeLayersSpatialProfilesUtilities::
outputAveragedMassFractionSpeciesThreeWithHomogeneityInYDirectionOrInYZPlane(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 3)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "MASS_FRACTION_S3_AVG' can be computed with three species only."
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
    
    std::vector<double> Y_2_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        "MASS_FRACTIONS",
        2,
        data_context);
    
    /*
     * Output the spatial profile (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        f_out.write((char*)&output_time, sizeof(double));
        f_out.write((char*)&Y_2_avg_global[0], sizeof(double)*Y_2_avg_global.size());
        
        f_out.close();
    }
}


/*
 * Output variance of mass fraction of species 1 with assumed homogeneity in y-direction (2D) or yz-plane (3D)
 * to a file.
 */
void
RTIThreeLayersSpatialProfilesUtilities::
outputMassFractionVarianceSpeciesOneWithHomogeneityInYDirectionOrInYZPlane(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 3)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'MASS_FRACTION_S1_VAR' can be computed with three species only."
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
RTIThreeLayersSpatialProfilesUtilities::
outputMassFractionVarianceSpeciesTwoWithHomogeneityInYDirectionOrInYZPlane(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 3)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'MASS_FRACTION_S2_VAR' can be computed with three species only."
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
 * Output variance of mass fraction of species 3 with assumed homogeneity in y-direction (2D) or yz-plane (3D)
 * to a file.
 */
void
RTIThreeLayersSpatialProfilesUtilities::
outputMassFractionVarianceSpeciesThreeWithHomogeneityInYDirectionOrInYZPlane(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 3)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'MASS_FRACTION_S3_VAR' can be computed with three species only."
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
    
    std::vector<double> Y_2_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        "MASS_FRACTIONS",
        2,
        data_context);
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    std::vector<std::vector<double> > averaged_quantities;
    
    quantity_names.push_back("MASS_FRACTIONS");
    component_indices.push_back(2);
    averaged_quantities.push_back(Y_2_avg_global);
    
    quantity_names.push_back("MASS_FRACTIONS");
    component_indices.push_back(2);
    averaged_quantities.push_back(Y_2_avg_global);
    
    std::vector<double> Y_2_var_global = MPI_helper_correlation.getQuantityCorrelationWithInhomogeneousXDirection(
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
        f_out.write((char*)&Y_2_var_global[0], sizeof(double)*Y_2_var_global.size());
        
        f_out.close();
    }
}


/*
 * Output covariance of mass fractions of species 1 and 2 with assumed homogeneity in y-direction (2D) or yz-plane (3D)
 * to a file.
 */
void
RTIThreeLayersSpatialProfilesUtilities::
outputMassFractionCovarianceSpeciesOneTwoWithHomogeneityInYDirectionOrInYZPlane(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 3)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'MASS_FRACTION_S12_COVAR' can be computed with three species only."
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
 * Output covariance of mass fractions of species 2 and 3 with assumed homogeneity in y-direction (2D) or yz-plane (3D)
 * to a file.
 */
void
RTIThreeLayersSpatialProfilesUtilities::
outputMassFractionCovarianceSpeciesTwoThreeWithHomogeneityInYDirectionOrInYZPlane(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 3)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'MASS_FRACTION_S23_COVAR' can be computed with three species only."
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
    
    std::vector<double> Y_2_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        "MASS_FRACTIONS",
        2,
        data_context);
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    std::vector<std::vector<double> > averaged_quantities;
    
    quantity_names.push_back("MASS_FRACTIONS");
    component_indices.push_back(1);
    averaged_quantities.push_back(Y_1_avg_global);
    
    quantity_names.push_back("MASS_FRACTIONS");
    component_indices.push_back(2);
    averaged_quantities.push_back(Y_2_avg_global);
    
    std::vector<double> Y_12_covar_global = MPI_helper_correlation.getQuantityCorrelationWithInhomogeneousXDirection(
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
        f_out.write((char*)&Y_12_covar_global[0], sizeof(double)*Y_12_covar_global.size());
        
        f_out.close();
    }
}


/*
 * Output covariance of mass fractions of species 1 and 3 with assumed homogeneity in y-direction (2D) or yz-plane (3D)
 * to a file.
 */
void
RTIThreeLayersSpatialProfilesUtilities::
outputMassFractionCovarianceSpeciesOneThreeWithHomogeneityInYDirectionOrInYZPlane(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double output_time) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 3)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'MASS_FRACTION_S13_COVAR' can be computed with three species only."
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
    
    std::vector<double> Y_2_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        "MASS_FRACTIONS",
        2,
        data_context);
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    std::vector<std::vector<double> > averaged_quantities;
    
    quantity_names.push_back("MASS_FRACTIONS");
    component_indices.push_back(0);
    averaged_quantities.push_back(Y_0_avg_global);
    
    quantity_names.push_back("MASS_FRACTIONS");
    component_indices.push_back(2);
    averaged_quantities.push_back(Y_2_avg_global);
    
    std::vector<double> Y_02_covar_global = MPI_helper_correlation.getQuantityCorrelationWithInhomogeneousXDirection(
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
        f_out.write((char*)&Y_02_covar_global[0], sizeof(double)*Y_02_covar_global.size());
        
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
    
    HAMERS_SHARED_PTR<RTIThreeLayersSpatialProfilesUtilities> rti_three_layers_spatial_profiles_utilities(
        new RTIThreeLayersSpatialProfilesUtilities(
            "RTI three layers spatial profiles utilities",
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
        
        if (statistical_quantity_key == "MASS_FRACTION_S1_AVG")
        {
            rti_three_layers_spatial_profiles_utilities->outputAveragedMassFractionSpeciesOneWithHomogeneityInYDirectionOrInYZPlane(
                "Y1_avg.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "MASS_FRACTION_S2_AVG")
        {
            rti_three_layers_spatial_profiles_utilities->outputAveragedMassFractionSpeciesTwoWithHomogeneityInYDirectionOrInYZPlane(
                "Y2_avg.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "MASS_FRACTION_S3_AVG")
        {
            rti_three_layers_spatial_profiles_utilities->outputAveragedMassFractionSpeciesThreeWithHomogeneityInYDirectionOrInYZPlane(
                "Y3_avg.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "MASS_FRACTION_S1_VAR")
        {
            rti_three_layers_spatial_profiles_utilities->outputMassFractionVarianceSpeciesOneWithHomogeneityInYDirectionOrInYZPlane(
                "Y1_var.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "MASS_FRACTION_S2_VAR")
        {
            rti_three_layers_spatial_profiles_utilities->outputMassFractionVarianceSpeciesTwoWithHomogeneityInYDirectionOrInYZPlane(
                "Y2_var.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "MASS_FRACTION_S3_VAR")
        {
            rti_three_layers_spatial_profiles_utilities->outputMassFractionVarianceSpeciesThreeWithHomogeneityInYDirectionOrInYZPlane(
                "Y3_var.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "MASS_FRACTION_S12_COVAR")
        {
            rti_three_layers_spatial_profiles_utilities->outputMassFractionCovarianceSpeciesOneTwoWithHomogeneityInYDirectionOrInYZPlane(
                "Y1Y2_covar.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "MASS_FRACTION_S23_COVAR")
        {
            rti_three_layers_spatial_profiles_utilities->outputMassFractionCovarianceSpeciesTwoThreeWithHomogeneityInYDirectionOrInYZPlane(
                "Y2Y3_covar.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "MASS_FRACTION_S13_COVAR")
        {
            rti_three_layers_spatial_profiles_utilities->outputMassFractionCovarianceSpeciesOneThreeWithHomogeneityInYDirectionOrInYZPlane(
                "Y1Y3_covar.dat",
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
