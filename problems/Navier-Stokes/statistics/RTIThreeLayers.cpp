#include "flow/flow_models/four-eqn_conservative/FlowModelStatisticsUtilitiesFourEqnConservative.hpp"

#include "extn/patch_hierarchies/ExtendedFlattenedHierarchy.hpp"
#include "flow/flow_models/MPI_helpers/FlowModelMPIHelperAverage.hpp"
#include "flow/flow_models/MPI_helpers/FlowModelMPIHelperCorrelation.hpp"
#include "flow/flow_models/MPI_helpers/FlowModelMPIHelperMaxMin.hpp"
#include "util/MPI_helpers/MPIHelperGrid.hpp"

#include <fstream>

class RTIThreeLayersStatisticsUtilities
{
    public:
        RTIThreeLayersStatisticsUtilities(
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
         * Functions for stats calculations for 3 layers RTI.
         */
        
        /*
         * Output minimum interface location of 1st species using the averaged 1D mass fraction profile
         * in x-direction to a file.
         */
        void
        outputAveragedInterfaceMinSpeciesOneInXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output maximum interface location of 1st species using the averaged 1D mass fraction profile
         * in x-direction to a file.
         */
        void
        outputAveragedInterfaceMaxSpeciesOneInXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output minimum interface location of 2nd species using the averaged 1D mass fraction profile
         * in x-direction to a file.
         */
        void
        outputAveragedInterfaceMinSpeciesTwoInXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output maximum interface location of 2nd species using the averaged 1D mass fraction profile
         * in x-direction to a file.
         */
        void
        outputAveragedInterfaceMaxSpeciesTwoInXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output minimum interface location of 3rd species using the averaged 1D mass fraction profile
         * in x-direction to a file.
         */
        void
        outputAveragedInterfaceMinSpeciesThreeInXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output maximum interface location of 3rd species using the averaged 1D mass fraction profile
         * in x-direction to a file.
         */
        void
        outputAveragedInterfaceMaxSpeciesThreeInXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output minimum interface location of 1st species in x-direction to a file.
         */
        void
        outputInterfaceMinSpeciesOneInXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output maximum interface location of 1st species in x-direction to a file.
         */
        void
        outputInterfaceMaxSpeciesOneInXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output minimum interface location of 2nd species in x-direction to a file.
         */
        void
        outputInterfaceMinSpeciesTwoInXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output maximum interface location of 2nd species in x-direction to a file.
         */
        void
        outputInterfaceMaxSpeciesTwoInXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output minimum interface location of 3rd species in x-direction to a file.
         */
        void
        outputInterfaceMinSpeciesThreeInXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output maximum interface location of 3rd species in x-direction to a file.
         */
        void
        outputInterfaceMaxSpeciesThreeInXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output mixing width in x-direction for 1st species in the domain to a file.
         */
        void
        outputMixingWidthSpeciesOneInXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output mixing width in x-direction for 2nd species in the domain to a file.
         */
        void
        outputMixingWidthSpeciesTwoInXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output mixing width in x-direction for 3rd species in the domain to a file.
         */
        void
        outputMixingWidthSpeciesThreeInXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output mixedness in x-direction for 1st species to a file.
         */
        void
        outputMixednessSpeciesOneInXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output mixedness in x-direction for 2nd species to a file.
         */
        void
        outputMixednessSpeciesTwoInXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output mixedness in x-direction for 3rd species to a file.
         */
        void
        outputMixednessSpeciesThreeInXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output integrated mixing metric between 1st and 2nd species to a file.
         */
        void
        outputMixingMetricSpeciesOneAndTwoIntegrated(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output integrated mixing metric between 2nd and 3rd species to a file.
         */
        void
        outputMixingMetricSpeciesTwoAndThreeIntegrated(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output integrated mixing metric between 1st and 3rd species to a file.
         */
        void
        outputMixingMetricSpeciesOneAndThreeIntegrated(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        // /*
        //  * Output Reynolds normal stress component in x-direction with assumed homogeneity in y-direction (2D) or yz-plane (3D)
        //  * to a file.
        //  */
        // void
        // outputReynoldsNormalStressXWithInhomogeneousXDirection(
        //     const std::string& stat_dump_filename,
        //     const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
        //     const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        // /*
        //  * Output Reynolds normal stress component in y-direction with assumed homogeneity in y-direction (2D) or yz-plane (3D)
        //  * to a file.
        //  */
        // void
        // outputReynoldsNormalStressYWithInhomogeneousXDirection(
        //     const std::string& stat_dump_filename,
        //     const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
        //     const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output enstrophy integrated to a file.
         */
        void
        outputEnstrophyIntegrated(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output TKE integrated with assumed homogeneity in y-direction to a file.
         */
        void
        outputTKEIntegratedWithHomogeneityInYDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output scalar dissipation rate (can be negative) of between 1st and 2nd species integrated to a file.
         */
        void
        outputScalarDissipationRateSpeciesOneAndTwoIntegrated(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output scalar dissipation rate (can be negative) of between 2nd and 3rd species integrated to a file.
         */
        void
        outputScalarDissipationRateSpeciesTwoAndThreeIntegrated(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output scalar dissipation rate (can be negative) of between 1st and 3rd species integrated to a file.
         */
        void
        outputScalarDissipationRateSpeciesOneAndThreeIntegrated(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output maximum value of species 1 in the domain to a file.
         */
        void
        outputMaxMassFractionSpeciesOne(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output maximum value of species 2 in the domain to a file.
         */
        void
        outputMaxMassFractionSpeciesTwo(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output maximum value of species 3 in the domain to a file.
         */
        void
        outputMaxMassFractionSpeciesThree(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
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
 * Output minimum interface location of 1st species using the averaged 1D mass fraction profile
 * in x-direction to a file.
 */
void
RTIThreeLayersStatisticsUtilities::outputAveragedInterfaceMinSpeciesOneInXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
    if (d_num_species != 3)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'INTERFACE_AVG_MIN_S1' can be computed with three species only."
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
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
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
    
    const std::vector<double>& dx_finest = MPI_helper_average.getFinestRefinedDomainGridSpacing();

    const hier::IntVector& finest_level_dims = MPI_helper_average.getFinestRefinedDomainNumberOfPoints();

    std::vector<double> Y_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        "MASS_FRACTIONS",
        0,
        data_context);
    
    /*
     * Compute and output the value (only done by process 0).
     */
    
    const double lower_bound = double(1) - double(0.99);
    const double upper_bound = double(0.99);
    
    if (mpi.getRank() == 0)
    {
        const double* x_lo = d_grid_geometry->getXLower();
        const double* x_hi = d_grid_geometry->getXUpper();
        double interface_min = x_hi[0];
        
        for (int i = finest_level_dims[0] - 1; i >= 0;  i--)
        {
            if (Y_avg_global[i] > lower_bound && Y_avg_global[i] < upper_bound)
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
 * Output maximum interface location of 1st species using the averaged 1D mass fraction profile
 * in x-direction to a file.
 */
void
RTIThreeLayersStatisticsUtilities::outputAveragedInterfaceMaxSpeciesOneInXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
    if (d_num_species != 3)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'INTERFACE_AVG_MAX_S1' can be computed with three species only."
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
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
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
    
    const std::vector<double>& dx_finest = MPI_helper_average.getFinestRefinedDomainGridSpacing();

    const hier::IntVector& finest_level_dims = MPI_helper_average.getFinestRefinedDomainNumberOfPoints();

    std::vector<double> Y_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        "MASS_FRACTIONS",
        0,
        data_context);
    
    /*
     * Compute and output the value (only done by process 0).
     */
    
    const double lower_bound = double(1) - double(0.99);
    const double upper_bound = double(0.99);
    
    if (mpi.getRank() == 0)
    {
        const double* x_lo = d_grid_geometry->getXLower();
        // const double* x_hi = d_grid_geometry->getXUpper();
        double interface_max = x_lo[0];
        
        for (int i = 0; i < finest_level_dims[0]; i++)
        {
            if (Y_avg_global[i] > lower_bound && Y_avg_global[i] < upper_bound)
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
 * Output minimum interface location of 2nd species using the averaged 1D mass fraction profile
 * in x-direction to a file.
 */
void
RTIThreeLayersStatisticsUtilities::outputAveragedInterfaceMinSpeciesTwoInXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
    if (d_num_species != 3)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'INTERFACE_AVG_MIN_S2' can be computed with three species only."
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
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
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
    
    const std::vector<double>& dx_finest = MPI_helper_average.getFinestRefinedDomainGridSpacing();

    const hier::IntVector& finest_level_dims = MPI_helper_average.getFinestRefinedDomainNumberOfPoints();

    std::vector<double> Y_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        "MASS_FRACTIONS",
        1,
        data_context);
    
    /*
     * Compute and output the value (only done by process 0).
     */
    
    const double lower_bound = double(1) - double(0.99);
    const double upper_bound = double(0.99);
    
    if (mpi.getRank() == 0)
    {
        const double* x_lo = d_grid_geometry->getXLower();
        const double* x_hi = d_grid_geometry->getXUpper();
        double interface_min = x_hi[0];
        
        for (int i = finest_level_dims[0] - 1; i >= 0;  i--)
        {
            if (Y_avg_global[i] > lower_bound && Y_avg_global[i] < upper_bound)
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
 * Output maximum interface location of 2nd species using the averaged 1D mass fraction profile
 * in x-direction to a file.
 */
void
RTIThreeLayersStatisticsUtilities::outputAveragedInterfaceMaxSpeciesTwoInXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
    if (d_num_species != 3)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'INTERFACE_AVG_MAX_S2' can be computed with three species only."
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
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
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
    
    const std::vector<double>& dx_finest = MPI_helper_average.getFinestRefinedDomainGridSpacing();

    const hier::IntVector& finest_level_dims = MPI_helper_average.getFinestRefinedDomainNumberOfPoints();

    std::vector<double> Y_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        "MASS_FRACTIONS",
        1,
        data_context);
    
    /*
     * Compute and output the value (only done by process 0).
     */
    
    const double lower_bound = double(1) - double(0.99);
    const double upper_bound = double(0.99);
    
    if (mpi.getRank() == 0)
    {
        const double* x_lo = d_grid_geometry->getXLower();
        // const double* x_hi = d_grid_geometry->getXUpper();
        double interface_max = x_lo[0];
        
        for (int i = 0; i < finest_level_dims[0]; i++)
        {
            if (Y_avg_global[i] > lower_bound && Y_avg_global[i] < upper_bound)
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
 * Output minimum interface location of 3rd species using the averaged 1D mass fraction profile
 * in x-direction to a file.
 */
void
RTIThreeLayersStatisticsUtilities::outputAveragedInterfaceMinSpeciesThreeInXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
    if (d_num_species != 3)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'INTERFACE_AVG_MIN_S3' can be computed with three species only."
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
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
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
    
    const std::vector<double>& dx_finest = MPI_helper_average.getFinestRefinedDomainGridSpacing();

    const hier::IntVector& finest_level_dims = MPI_helper_average.getFinestRefinedDomainNumberOfPoints();

    std::vector<double> Y_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        "MASS_FRACTIONS",
        2,
        data_context);
    
    /*
     * Compute and output the value (only done by process 0).
     */
    
    const double lower_bound = double(1) - double(0.99);
    const double upper_bound = double(0.99);
    
    if (mpi.getRank() == 0)
    {
        const double* x_lo = d_grid_geometry->getXLower();
        const double* x_hi = d_grid_geometry->getXUpper();
        double interface_min = x_hi[0];
        
        for (int i = finest_level_dims[0] - 1; i >= 0;  i--)
        {
            if (Y_avg_global[i] > lower_bound && Y_avg_global[i] < upper_bound)
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
 * Output maximum interface location of 3rd species using the averaged 1D mass fraction profile
 * in x-direction to a file.
 */
void
RTIThreeLayersStatisticsUtilities::outputAveragedInterfaceMaxSpeciesThreeInXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
    if (d_num_species != 3)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'INTERFACE_AVG_MAX_S3' can be computed with three species only."
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
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
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
    
    const std::vector<double>& dx_finest = MPI_helper_average.getFinestRefinedDomainGridSpacing();

    const hier::IntVector& finest_level_dims = MPI_helper_average.getFinestRefinedDomainNumberOfPoints();

    std::vector<double> Y_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        "MASS_FRACTIONS",
        2,
        data_context);
    
    /*
     * Compute and output the value (only done by process 0).
     */
    
    const double lower_bound = double(1) - double(0.99);
    const double upper_bound = double(0.99);
    
    if (mpi.getRank() == 0)
    {
        const double* x_lo = d_grid_geometry->getXLower();
        // const double* x_hi = d_grid_geometry->getXUpper();
        double interface_max = x_lo[0];
        
        for (int i = 0; i < finest_level_dims[0]; i++)
        {
            if (Y_avg_global[i] > lower_bound && Y_avg_global[i] < upper_bound)
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
 * Output minimum interface location of 1st species in x-direction to a file.
 */
void
RTIThreeLayersStatisticsUtilities::outputInterfaceMinSpeciesOneInXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 3)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'INTERFACE_MIN_S1' can be computed with three species only."
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
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
    }
    
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
    FlowModelMPIHelperMaxMin MPI_helper_max_min = FlowModelMPIHelperMaxMin(
        "MPI_helper_max_min",
        d_dim,
        d_grid_geometry,
        patch_hierarchy,
        flow_model_tmp);
    
    const double interface_x_min_global = MPI_helper_max_min.getMinLocationWithinQuantityBoundsInXDirection(
            "MASS_FRACTIONS",
            0,
            data_context,
            double(1) - double(0.99),
            double(0.99));
    
    /*
     * Compute and output the value (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << interface_x_min_global;
        
        f_out.close();
    }
}


/*
 * Output maximum interface location of 1st species in x-direction to a file.
 */
void
RTIThreeLayersStatisticsUtilities::outputInterfaceMaxSpeciesOneInXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 3)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'INTERFACE_MAX_S1' can be computed with three species only."
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
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
    }
    
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
    FlowModelMPIHelperMaxMin MPI_helper_max_min = FlowModelMPIHelperMaxMin(
        "MPI_helper_max_min",
        d_dim,
        d_grid_geometry,
        patch_hierarchy,
        flow_model_tmp);
    
    const double interface_x_max_global = MPI_helper_max_min.getMaxLocationWithinQuantityBoundsInXDirection(
            "MASS_FRACTIONS",
            0,
            data_context,
            double(1) - double(0.99),
            double(0.99));
    
    /*
     * Compute and output the value (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << interface_x_max_global;
        
        f_out.close();
    }
}


/*
 * Output minimum interface location of 2nd species in x-direction to a file.
 */
void
RTIThreeLayersStatisticsUtilities::outputInterfaceMinSpeciesTwoInXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 3)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'INTERFACE_MIN_S2' can be computed with three species only."
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
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
    }
    
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
    FlowModelMPIHelperMaxMin MPI_helper_max_min = FlowModelMPIHelperMaxMin(
        "MPI_helper_max_min",
        d_dim,
        d_grid_geometry,
        patch_hierarchy,
        flow_model_tmp);
    
    const double interface_x_min_global = MPI_helper_max_min.getMinLocationWithinQuantityBoundsInXDirection(
            "MASS_FRACTIONS",
            1,
            data_context,
            double(1) - double(0.99),
            double(0.99));
    
    /*
     * Compute and output the value (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << interface_x_min_global;
        
        f_out.close();
    }
}


/*
 * Output maximum interface location of 2nd species in x-direction to a file.
 */
void
RTIThreeLayersStatisticsUtilities::outputInterfaceMaxSpeciesTwoInXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 3)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'INTERFACE_MAX_S2' can be computed with three species only."
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
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
    }
    
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
    FlowModelMPIHelperMaxMin MPI_helper_max_min = FlowModelMPIHelperMaxMin(
        "MPI_helper_max_min",
        d_dim,
        d_grid_geometry,
        patch_hierarchy,
        flow_model_tmp);
    
    const double interface_x_max_global = MPI_helper_max_min.getMaxLocationWithinQuantityBoundsInXDirection(
            "MASS_FRACTIONS",
            1,
            data_context,
            double(1) - double(0.99),
            double(0.99));
    
    /*
     * Compute and output the value (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << interface_x_max_global;
        
        f_out.close();
    }
}


/*
 * Output minimum interface location of 3rd species in x-direction to a file.
 */
void
RTIThreeLayersStatisticsUtilities::outputInterfaceMinSpeciesThreeInXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 3)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'INTERFACE_MIN_S3' can be computed with three species only."
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
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
    }
    
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
    FlowModelMPIHelperMaxMin MPI_helper_max_min = FlowModelMPIHelperMaxMin(
        "MPI_helper_max_min",
        d_dim,
        d_grid_geometry,
        patch_hierarchy,
        flow_model_tmp);
    
    const double interface_x_min_global = MPI_helper_max_min.getMinLocationWithinQuantityBoundsInXDirection(
            "MASS_FRACTIONS",
            2,
            data_context,
            double(1) - double(0.99),
            double(0.99));
    
    /*
     * Compute and output the value (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << interface_x_min_global;
        
        f_out.close();
    }
}


/*
 * Output maximum interface location of 3rd species in x-direction to a file.
 */
void
RTIThreeLayersStatisticsUtilities::outputInterfaceMaxSpeciesThreeInXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 3)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'INTERFACE_MAX_S3' can be computed with three species only."
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
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
    }
    
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
    FlowModelMPIHelperMaxMin MPI_helper_max_min = FlowModelMPIHelperMaxMin(
        "MPI_helper_max_min",
        d_dim,
        d_grid_geometry,
        patch_hierarchy,
        flow_model_tmp);
    
    const double interface_x_max_global = MPI_helper_max_min.getMaxLocationWithinQuantityBoundsInXDirection(
            "MASS_FRACTIONS",
            2,
            data_context,
            double(1) - double(0.99),
            double(0.99));
    
    /*
     * Compute and output the value (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << interface_x_max_global;
        
        f_out.close();
    }
}


/*
 * Output mixing width in x-direction for 1st species in the domain to a file.
 */
void
RTIThreeLayersStatisticsUtilities::outputMixingWidthSpeciesOneInXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 3)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'MIXING_WIDTH_X_S1' can be computed with three species only."
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
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
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
    
    const std::vector<double>& dx_finest = MPI_helper_average.getFinestRefinedDomainGridSpacing();
    
    const hier::IntVector& finest_level_dims = MPI_helper_average.getFinestRefinedDomainNumberOfPoints();
    
    std::vector<double> Y_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        "MASS_FRACTIONS",
        0,
        data_context);
    
    /*
     * Compute and output the mixing width (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        double W = double(0);
        
        for (int i = 0; i < finest_level_dims[0]; i++)
        {
            W += Y_avg_global[i]*(double(1) - Y_avg_global[i]);
        }
        
        W = double(4)*W*dx_finest[0];
        
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << W;
        
        f_out.close();
    }
}


/*
 * Output mixing width in x-direction for 2nd species in the domain to a file.
 */
void
RTIThreeLayersStatisticsUtilities::outputMixingWidthSpeciesTwoInXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 3)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'MIXING_WIDTH_X_S2' can be computed with three species only."
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
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
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
    
    const std::vector<double>& dx_finest = MPI_helper_average.getFinestRefinedDomainGridSpacing();
    
    const hier::IntVector& finest_level_dims = MPI_helper_average.getFinestRefinedDomainNumberOfPoints();
    
    std::vector<double> Y_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        "MASS_FRACTIONS",
        1,
        data_context);
    
    /*
     * Compute and output the mixing width (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        double W = double(0);
        
        for (int i = 0; i < finest_level_dims[0]; i++)
        {
            W += Y_avg_global[i]*(double(1) - Y_avg_global[i]);
        }
        
        W = double(4)*W*dx_finest[0];
        
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << W;
        
        f_out.close();
    }
}


/*
 * Output mixing width in x-direction for 3rd species in the domain to a file.
 */
void
RTIThreeLayersStatisticsUtilities::outputMixingWidthSpeciesThreeInXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 3)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'MIXING_WIDTH_X_S3' can be computed with three species only."
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
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
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
    
    const std::vector<double>& dx_finest = MPI_helper_average.getFinestRefinedDomainGridSpacing();
    
    const hier::IntVector& finest_level_dims = MPI_helper_average.getFinestRefinedDomainNumberOfPoints();
    
    std::vector<double> Y_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        "MASS_FRACTIONS",
        2,
        data_context);
    
    /*
     * Compute and output the mixing width (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        double W = double(0);
        
        for (int i = 0; i < finest_level_dims[0]; i++)
        {
            W += Y_avg_global[i]*(double(1) - Y_avg_global[i]);
        }
        
        W = double(4)*W*dx_finest[0];
        
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << W;
        
        f_out.close();
    }
}


/*
 * Output mixedness in x-direction for 1st species to a file.
 */
void
RTIThreeLayersStatisticsUtilities::outputMixednessSpeciesOneInXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 3)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'MIXEDNESS_X_S1' can be computed with three species only."
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
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
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
    
    const hier::IntVector& finest_level_dims = MPI_helper_average.getFinestRefinedDomainNumberOfPoints();
    
    std::vector<double> Y_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        "MASS_FRACTIONS",
        0,
        data_context);
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    
    quantity_names.push_back("MASS_FRACTIONS");
    component_indices.push_back(0);
    quantity_names.push_back("MASS_FRACTIONS");
    component_indices.push_back(0);
    
    std::vector<double> Y_product_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    
    /*
     * Compute and output the mixedness (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        double num = double(0);
        double den = double(0);
        
        for (int i = 0; i < finest_level_dims[0]; i++)
        {
            num += Y_avg_global[i] - Y_product_avg_global[i];
        }
        
        for (int i = 0; i < finest_level_dims[0]; i++)
        {
            den += Y_avg_global[i]*(double(1) - Y_avg_global[i]);
        }
        
        const double Theta = num/den;
        
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << Theta;
        
        f_out.close();
    }
}


/*
 * Output mixedness in x-direction for 2nd species to a file.
 */
void
RTIThreeLayersStatisticsUtilities::outputMixednessSpeciesTwoInXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 3)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'MIXEDNESS_X_S2' can be computed with three species only."
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
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
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
    
    const hier::IntVector& finest_level_dims = MPI_helper_average.getFinestRefinedDomainNumberOfPoints();
    
    std::vector<double> Y_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        "MASS_FRACTIONS",
        1,
        data_context);
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    
    quantity_names.push_back("MASS_FRACTIONS");
    component_indices.push_back(1);
    quantity_names.push_back("MASS_FRACTIONS");
    component_indices.push_back(1);
    
    std::vector<double> Y_product_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    
    /*
     * Compute and output the mixedness (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        double num = double(0);
        double den = double(0);
        
        for (int i = 0; i < finest_level_dims[0]; i++)
        {
            num += Y_avg_global[i] - Y_product_avg_global[i];
        }
        
        for (int i = 0; i < finest_level_dims[0]; i++)
        {
            den += Y_avg_global[i]*(double(1) - Y_avg_global[i]);
        }
        
        const double Theta = num/den;
        
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << Theta;
        
        f_out.close();
    }
}


/*
 * Output mixedness in x-direction for 3rd species to a file.
 */
void
RTIThreeLayersStatisticsUtilities::outputMixednessSpeciesThreeInXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 3)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'MIXEDNESS_X_S3' can be computed with three species only."
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
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
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
    
    const hier::IntVector& finest_level_dims = MPI_helper_average.getFinestRefinedDomainNumberOfPoints();
    
    std::vector<double> Y_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        "MASS_FRACTIONS",
        2,
        data_context);
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    
    quantity_names.push_back("MASS_FRACTIONS");
    component_indices.push_back(2);
    quantity_names.push_back("MASS_FRACTIONS");
    component_indices.push_back(2);
    
    std::vector<double> Y_product_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    
    /*
     * Compute and output the mixedness (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        double num = double(0);
        double den = double(0);
        
        for (int i = 0; i < finest_level_dims[0]; i++)
        {
            num += Y_avg_global[i] - Y_product_avg_global[i];
        }
        
        for (int i = 0; i < finest_level_dims[0]; i++)
        {
            den += Y_avg_global[i]*(double(1) - Y_avg_global[i]);
        }
        
        const double Theta = num/den;
        
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << Theta;
        
        f_out.close();
    }
}


/*
 * Output integrated mixing metric between 1st and 2nd species to a file.
 */
void
RTIThreeLayersStatisticsUtilities::outputMixingMetricSpeciesOneAndTwoIntegrated(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 3)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'MIX_METRIC_INT_S12' can be computed with three species only."
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
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
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
    
    const std::vector<double>& dx_finest = MPI_helper_average.getFinestRefinedDomainGridSpacing();
    
    const hier::IntVector& finest_level_dims = MPI_helper_average.getFinestRefinedDomainNumberOfPoints();
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    
    quantity_names.push_back("MASS_FRACTIONS");
    component_indices.push_back(0);
    quantity_names.push_back("MASS_FRACTIONS");
    component_indices.push_back(1);
    
    std::vector<double> Y_product_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    
    /*
     * Compute and output the mixing metric (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        double Y1_Y2_integrated_global = double(0);
        
        for (int i = 0; i < finest_level_dims[0]; i++)
        {
            Y1_Y2_integrated_global += Y_product_avg_global[i]*dx_finest[0];
        }
        
        if (d_dim == tbox::Dimension(2))
        {
            const double* x_lo = d_grid_geometry->getXLower();
            const double* x_hi = d_grid_geometry->getXUpper();
            const double L_y = x_hi[1] - x_lo[1];
            
            Y1_Y2_integrated_global *= L_y;
        }
        else if (d_dim == tbox::Dimension(3))
        {
            const double* x_lo = d_grid_geometry->getXLower();
            const double* x_hi = d_grid_geometry->getXUpper();
            const double L_y = x_hi[1] - x_lo[1];
            const double L_z = x_hi[2] - x_lo[2];
            
            Y1_Y2_integrated_global *= (L_y*L_z);
        }
        
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << Y1_Y2_integrated_global;
        
        f_out.close();
    }
}


/*
 * Output integrated mixing metric between 2nd and 3rd species to a file.
 */
void
RTIThreeLayersStatisticsUtilities::outputMixingMetricSpeciesTwoAndThreeIntegrated(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 3)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'MIX_METRIC_INT_S23' can be computed with three species only."
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
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
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
    
    const std::vector<double>& dx_finest = MPI_helper_average.getFinestRefinedDomainGridSpacing();
    
    const hier::IntVector& finest_level_dims = MPI_helper_average.getFinestRefinedDomainNumberOfPoints();
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    
    quantity_names.push_back("MASS_FRACTIONS");
    component_indices.push_back(1);
    quantity_names.push_back("MASS_FRACTIONS");
    component_indices.push_back(2);
    
    std::vector<double> Y_product_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    
    /*
     * Compute and output the mixing metric (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        double Y2_Y3_integrated_global = double(0);
        
        for (int i = 0; i < finest_level_dims[0]; i++)
        {
            Y2_Y3_integrated_global += Y_product_avg_global[i]*dx_finest[0];
        }
        
        if (d_dim == tbox::Dimension(2))
        {
            const double* x_lo = d_grid_geometry->getXLower();
            const double* x_hi = d_grid_geometry->getXUpper();
            const double L_y = x_hi[1] - x_lo[1];
            
            Y2_Y3_integrated_global *= L_y;
        }
        else if (d_dim == tbox::Dimension(3))
        {
            const double* x_lo = d_grid_geometry->getXLower();
            const double* x_hi = d_grid_geometry->getXUpper();
            const double L_y = x_hi[1] - x_lo[1];
            const double L_z = x_hi[2] - x_lo[2];
            
            Y2_Y3_integrated_global *= (L_y*L_z);
        }
        
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << Y2_Y3_integrated_global;
        
        f_out.close();
    }
}


/*
 * Output integrated mixing metric between 1st and 3rd species to a file.
 */
void
RTIThreeLayersStatisticsUtilities::outputMixingMetricSpeciesOneAndThreeIntegrated(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 3)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'MIX_METRIC_INT_S13' can be computed with three species only."
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
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
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
    
    const std::vector<double>& dx_finest = MPI_helper_average.getFinestRefinedDomainGridSpacing();
    
    const hier::IntVector& finest_level_dims = MPI_helper_average.getFinestRefinedDomainNumberOfPoints();
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    
    quantity_names.push_back("MASS_FRACTIONS");
    component_indices.push_back(0);
    quantity_names.push_back("MASS_FRACTIONS");
    component_indices.push_back(2);
    
    std::vector<double> Y_product_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    
    /*
     * Compute and output the mixing metric (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        double Y1_Y3_integrated_global = double(0);
        
        for (int i = 0; i < finest_level_dims[0]; i++)
        {
            Y1_Y3_integrated_global += Y_product_avg_global[i]*dx_finest[0];
        }
        
        if (d_dim == tbox::Dimension(2))
        {
            const double* x_lo = d_grid_geometry->getXLower();
            const double* x_hi = d_grid_geometry->getXUpper();
            const double L_y = x_hi[1] - x_lo[1];
            
            Y1_Y3_integrated_global *= L_y;
        }
        else if (d_dim == tbox::Dimension(3))
        {
            const double* x_lo = d_grid_geometry->getXLower();
            const double* x_hi = d_grid_geometry->getXUpper();
            const double L_y = x_hi[1] - x_lo[1];
            const double L_z = x_hi[2] - x_lo[2];
            
            Y1_Y3_integrated_global *= (L_y*L_z);
        }
        
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << Y1_Y3_integrated_global;
        
        f_out.close();
    }
}


// /*
//  * Output Reynolds normal stress component in x-direction with assumed homogeneity in y-direction (2D) or yz-plane (3D)
//  * to a file.
//  */
// void
// RTIThreeLayersStatisticsUtilities::outputReynoldsNormalStressXWithInhomogeneousXDirection(
//     const std::string& stat_dump_filename,
//     const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
//     const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
// {
// #ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
//     TBOX_ASSERT(!stat_dump_filename.empty());
// #endif
//     
//     if (d_flow_model.expired())
//     {
//         TBOX_ERROR(d_object_name
//             << ": "
//             << "The object is not setup yet!"
//             << std::endl);
//     }
//     
//     const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
//     
//     std::ofstream f_out;
//     
//     if (mpi.getRank() == 0)
//     {
//         f_out.open(stat_dump_filename.c_str(), std::ios::app);
//         if (!f_out.is_open())
//         {
//             TBOX_ERROR(d_object_name
//                 << ": "
//                 << "Failed to open file to output statistics!"
//                 << std::endl);
//         }
//     }
//     
//     HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
//     
//     FlowModelMPIHelperAverage MPI_helper_average = FlowModelMPIHelperAverage(
//         "MPI_helper_average",
//         d_dim,
//         d_grid_geometry,
//         patch_hierarchy,
//         flow_model_tmp);
//     
//     FlowModelMPIHelperCorrelation MPI_helper_correlation = FlowModelMPIHelperCorrelation(
//         "MPI_helper_average",
//         d_dim,
//         d_grid_geometry,
//         patch_hierarchy,
//         flow_model_tmp);
//     
//     const hier::IntVector& finest_level_dims = MPI_helper_average.getFinestRefinedDomainNumberOfPoints();
//     
//     std::vector<double> Y_mean = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
//         "MASS_FRACTIONS",
//         0,
//         data_context);
//     
//     std::vector<double> rho_mean = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
//         "DENSITY",
//         0,
//         data_context);
//     
//     // Compute u_tilde.
//     
//     std::vector<std::string> quantity_names;
//     std::vector<int> component_indices;
//     
//     quantity_names.push_back("DENSITY");
//     component_indices.push_back(0);
//     
//     quantity_names.push_back("VELOCITY");
//     component_indices.push_back(0);
//     
//     std::vector<double> rho_u_mean = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
//         quantity_names,
//         component_indices,
//         data_context);
//     
//     quantity_names.clear();
//     component_indices.clear();
//     
//     std::vector<double> u_tilde(rho_u_mean);
//     
//     for (int i = 0; i < finest_level_dims[0]; i++)
//     {
//         u_tilde[i] /= rho_mean[i];
//     }
//     
//     // Compute R_11.
//     
//     std::vector<double> zeros(finest_level_dims[0], double(0));
//     
//     std::vector<std::vector<double> > averaged_quantities;
//     
//     quantity_names.push_back("DENSITY");
//     component_indices.push_back(0);
//     averaged_quantities.push_back(zeros);
//     
//     quantity_names.push_back("VELOCITY");
//     component_indices.push_back(0);
//     averaged_quantities.push_back(u_tilde);
//     
//     quantity_names.push_back("VELOCITY");
//     component_indices.push_back(0);
//     averaged_quantities.push_back(u_tilde);
//     
//     std::vector<double> R_11 = MPI_helper_correlation.getQuantityCorrelationWithInhomogeneousXDirection(
//         quantity_names,
//         component_indices,
//         averaged_quantities,
//         data_context);
//     
//     quantity_names.clear();
//     component_indices.clear();
//     averaged_quantities.clear();
//     
//     for (int i = 0; i < finest_level_dims[0]; i++)
//     {
//         R_11[i] /= rho_mean[i];
//     }
//     
//     /*
//      * Compute and output the mean of R_11 inside mixing layer (only done by process 0).
//      */
//     if (mpi.getRank() == 0)
//     {
//         double R_11_sum = double(0);
//         int count = 0;
//         
//         for (int i = 0; i < finest_level_dims[0]; i++)
//         {
//             const double mixing_metric = double(4)*Y_mean[i]*(double(1) - Y_mean[i]);
//             if (mixing_metric > double(9)/double(10))
//             {
//                 R_11_sum += R_11[i];
//                 count++;
//             }
//         }
//         
//         const double R_11_mean = R_11_sum/count;
//         
//         f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
//               << "\t" << R_11_mean;
//         
//         f_out.close();
//     }
// }


// /*
//  * Output Reynolds normal stress component in y-direction with assumed homogeneity in y-direction (2D) or yz-plane (3D)
//  * to a file.
//  */
// void
// RTIThreeLayersStatisticsUtilities::outputReynoldsNormalStressYWithInhomogeneousXDirection(
//     const std::string& stat_dump_filename,
//     const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
//     const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
// {
// #ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
//     TBOX_ASSERT(!stat_dump_filename.empty());
// #endif
//     
//     if (d_flow_model.expired())
//     {
//         TBOX_ERROR(d_object_name
//             << ": "
//             << "The object is not setup yet!"
//             << std::endl);
//     }
//     
//     const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
//     
//     std::ofstream f_out;
//     
//     if (mpi.getRank() == 0)
//     {
//         f_out.open(stat_dump_filename.c_str(), std::ios::app);
//         if (!f_out.is_open())
//         {
//             TBOX_ERROR(d_object_name
//                 << ": "
//                 << "Failed to open file to output statistics!"
//                 << std::endl);
//         }
//     }
//     
//     HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
//     
//     FlowModelMPIHelperAverage MPI_helper_average = FlowModelMPIHelperAverage(
//         "MPI_helper_average",
//         d_dim,
//         d_grid_geometry,
//         patch_hierarchy,
//         flow_model_tmp);
//     
//     FlowModelMPIHelperCorrelation MPI_helper_correlation = FlowModelMPIHelperCorrelation(
//         "MPI_helper_average",
//         d_dim,
//         d_grid_geometry,
//         patch_hierarchy,
//         flow_model_tmp);
//     
//     const hier::IntVector& finest_level_dims = MPI_helper_average.getFinestRefinedDomainNumberOfPoints();
//     
//     std::vector<double> Y_mean = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
//         "MASS_FRACTIONS",
//         0,
//         data_context);
//     
//     std::vector<double> rho_mean = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
//         "DENSITY",
//         0,
//         data_context);
//     
//     // Compute v_tilde.
//     
//     std::vector<std::string> quantity_names;
//     std::vector<int> component_indices;
//     
//     quantity_names.push_back("DENSITY");
//     component_indices.push_back(0);
//     
//     quantity_names.push_back("VELOCITY");
//     component_indices.push_back(1);
//     
//     std::vector<double> rho_v_mean = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
//         quantity_names,
//         component_indices,
//         data_context);
//     
//     quantity_names.clear();
//     component_indices.clear();
//     
//     std::vector<double> v_tilde(rho_v_mean);
//     
//     for (int i = 0; i < finest_level_dims[0]; i++)
//     {
//         v_tilde[i] /= rho_mean[i];
//     }
//     
//     // Compute R_22.
//     
//     std::vector<double> zeros(finest_level_dims[0], double(0));
//     
//     std::vector<std::vector<double> > averaged_quantities;
//     
//     quantity_names.push_back("DENSITY");
//     component_indices.push_back(0);
//     averaged_quantities.push_back(zeros);
//     
//     quantity_names.push_back("VELOCITY");
//     component_indices.push_back(1);
//     averaged_quantities.push_back(v_tilde);
//     
//     quantity_names.push_back("VELOCITY");
//     component_indices.push_back(1);
//     averaged_quantities.push_back(v_tilde);
//     
//     std::vector<double> R_22 = MPI_helper_correlation.getQuantityCorrelationWithInhomogeneousXDirection(
//         quantity_names,
//         component_indices,
//         averaged_quantities,
//         data_context);
//     
//     quantity_names.clear();
//     component_indices.clear();
//     averaged_quantities.clear();
//     
//     for (int i = 0; i < finest_level_dims[0]; i++)
//     {
//         R_22[i] /= rho_mean[i];
//     }
//     
//     /*
//      * Compute and output the mean of R_22 inside mixing layer (only done by process 0).
//      */
//     if (mpi.getRank() == 0)
//     {
//         double R_22_sum = double(0);
//         int count = 0;
//         
//         for (int i = 0; i < finest_level_dims[0]; i++)
//         {
//             const double mixing_metric = double(4)*Y_mean[i]*(double(1) - Y_mean[i]);
//             if (mixing_metric > double(9)/double(10))
//             {
//                 R_22_sum += R_22[i];
//                 count++;
//             }
//         }
//         
//         const double R_22_mean = R_22_sum/count;
//         
//         f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
//               << "\t" << R_22_mean;
//         
//         f_out.close();
//     }
// }


/*
 * Output enstrophy integrated to a file.
 */
void
RTIThreeLayersStatisticsUtilities::outputEnstrophyIntegrated(
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
    
    std::ofstream f_out;
    
    if (mpi.getRank() == 0)
    {
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
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
    
    const std::vector<double>& dx_finest = MPI_helper_average.getFinestRefinedDomainGridSpacing();
    
    const hier::IntVector& finest_level_dims = MPI_helper_average.getFinestRefinedDomainNumberOfPoints();
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no 'ENSTROPHY_INT' for one-dimensional problem."
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
        
        /*
         * Output the enstrophy integral (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            std::vector<double> enstrophy_full(finest_level_dims[0], double(0));
            for (int i = 0; i < finest_level_dims[0]; i++)
            {
                enstrophy_full[i] = enstrophy_part_1[i] + enstrophy_part_2[i] - double(2)*enstrophy_part_3[i];
            }
            
            double Omega_integrated_global = double(0);
            for (int i = 0; i < finest_level_dims[0]; i++)
            {
                Omega_integrated_global += enstrophy_full[i]*dx_finest[0];
            }
            
            const double* x_lo = d_grid_geometry->getXLower();
            const double* x_hi = d_grid_geometry->getXUpper();
            const double L_y = x_hi[1] - x_lo[1];
            Omega_integrated_global *= L_y;
            
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << Omega_integrated_global;
        }
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
        
        /*
         * Output the enstrophy integral (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            std::vector<double> enstrophy_full(finest_level_dims[0], double(0));
            for (int i = 0; i < finest_level_dims[0]; i++)
            {
                enstrophy_full[i] = enstrophy_part_1[i] + enstrophy_part_2[i] - double(2)*enstrophy_part_3[i] +
                    enstrophy_part_4[i] + enstrophy_part_5[i] - double(2)*enstrophy_part_6[i] +
                    enstrophy_part_7[i] + enstrophy_part_8[i] - double(2)*enstrophy_part_9[i];
            }
            
            double Omega_integrated_global = double(0);
            for (int i = 0; i < finest_level_dims[0]; i++)
            {
                Omega_integrated_global += enstrophy_full[i]*dx_finest[0];
            }
            
            const double* x_lo = d_grid_geometry->getXLower();
            const double* x_hi = d_grid_geometry->getXUpper();
            const double L_y = x_hi[1] - x_lo[1];
            const double L_z = x_hi[2] - x_lo[2];
            Omega_integrated_global *= (L_y*L_z);
            
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << Omega_integrated_global;
        }
    }
    
    if (mpi.getRank() == 0)
    {
        f_out.close();
    }
}


/*
 * Output TKE integrated with assumed homogeneity in y-direction to a file.
 */
void
RTIThreeLayersStatisticsUtilities::outputTKEIntegratedWithHomogeneityInYDirection(
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
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no 'TKE_INT_HOMO_Y' for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'TKE_INT_HOMO_Y' is not implemented for three-dimensional problem."
            << std::endl);
    }
    
    // Two-dimensional case.
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    std::ofstream f_out;
    
    if (mpi.getRank() == 0)
    {
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
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
        "MPI_helper_average",
        d_dim,
        d_grid_geometry,
        patch_hierarchy,
        flow_model_tmp);
    
    const std::vector<double>& dx_finest = MPI_helper_average.getFinestRefinedDomainGridSpacing();
    
    const hier::IntVector& finest_level_dims = MPI_helper_average.getFinestRefinedDomainNumberOfPoints();
    
    std::vector<double> zeros(finest_level_dims[0], double(0));
    
    std::vector<double> rho_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        data_context);
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    
    std::vector<double> rho_u_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(1);
    
    std::vector<double> rho_v_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    
    std::vector<double> u_tilde(rho_u_avg_global);
    std::vector<double> v_tilde(rho_v_avg_global);
    for (int i = 0; i < finest_level_dims[0]; i++)
    {
        u_tilde[i] /= rho_avg_global[i];
        v_tilde[i] /= rho_avg_global[i];
    }
    
    std::vector<std::vector<double> > averaged_quantities;
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(zeros);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(u_tilde);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(u_tilde);
    
    std::vector<double> two_TKE_x_avg_global = MPI_helper_correlation.getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    averaged_quantities.clear();
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(zeros);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(1);
    averaged_quantities.push_back(v_tilde);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(1);
    averaged_quantities.push_back(v_tilde);
    
    std::vector<double> two_TKE_y_avg_global = MPI_helper_correlation.getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    averaged_quantities.clear();
    
    double TKE_integrated_global = double(0);
    const double half = double(1)/double(2);
    for (int i = 0; i < finest_level_dims[0]; i++)
    {
        TKE_integrated_global += (half*(two_TKE_x_avg_global[i] + two_TKE_y_avg_global[i])*dx_finest[0]);
    }
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    const double L_y = x_hi[1] - x_lo[1];
    TKE_integrated_global *= L_y;
    
    if (mpi.getRank() == 0)
    {
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << TKE_integrated_global;
        
        f_out.close();
    }
}


/*
 * Output scalar dissipation rate (can be negative) of between 1st and 2nd species integrated to a file.
 */
void
RTIThreeLayersStatisticsUtilities::outputScalarDissipationRateSpeciesOneAndTwoIntegrated(
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
    
    std::ofstream f_out;
    
    if (mpi.getRank() == 0)
    {
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
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
        flow_model_tmp,
        true);
    
    const std::vector<double>& dx_finest = MPI_helper_average.getFinestRefinedDomainGridSpacing();
    
    const hier::IntVector& finest_level_dims = MPI_helper_average.getFinestRefinedDomainNumberOfPoints();
    
    if (d_dim == tbox::Dimension(1))
    {
        std::vector<std::string> quantity_names;
        std::vector<int> component_indices;
        std::vector<bool> use_derivative;
        std::vector<int> derivative_directions;
        
        quantity_names.push_back("MASS_DIFFUSIVITIES");
        component_indices.push_back(0); // Assuming effective mass diffusivities are same for all species.
        use_derivative.push_back(false);
        derivative_directions.push_back(-1);
        quantity_names.push_back("MASS_FRACTIONS");
        component_indices.push_back(0);
        use_derivative.push_back(true);
        derivative_directions.push_back(0);
        quantity_names.push_back("MASS_FRACTIONS");
        component_indices.push_back(1);
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
        
        /*
         * Output the scalar dissipation rate integral (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            double chi_integrated_global = double(0);
            for (int i = 0; i < finest_level_dims[0]; i++)
            {
                chi_integrated_global -= scalar_dissipation_rate[i]*dx_finest[0]; // Use minus sign.
            }
            
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << chi_integrated_global;
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        std::vector<std::string> quantity_names;
        std::vector<int> component_indices;
        std::vector<bool> use_derivative;
        std::vector<int> derivative_directions;
        
        quantity_names.push_back("MASS_DIFFUSIVITIES");
        component_indices.push_back(0); // Assuming effective mass diffusivities are same for all species.
        use_derivative.push_back(false);
        derivative_directions.push_back(-1);
        quantity_names.push_back("MASS_FRACTIONS");
        component_indices.push_back(0);
        use_derivative.push_back(true);
        derivative_directions.push_back(0);
        quantity_names.push_back("MASS_FRACTIONS");
        component_indices.push_back(1);
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
        component_indices.push_back(0); // Assuming effective mass diffusivities are same for all species.
        use_derivative.push_back(false);
        derivative_directions.push_back(-1);
        quantity_names.push_back("MASS_FRACTIONS");
        component_indices.push_back(0);
        use_derivative.push_back(true);
        derivative_directions.push_back(1);
        quantity_names.push_back("MASS_FRACTIONS");
        component_indices.push_back(1);
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
        
        /*
         * Output the scalar dissipation rate integral (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            std::vector<double> scalar_dissipation_rate_full(finest_level_dims[0], double(0));
            for (int i = 0; i < finest_level_dims[0]; i++)
            {
                scalar_dissipation_rate_full[i] = scalar_dissipation_rate_part_1[i] + scalar_dissipation_rate_part_2[i];
            }
            
            double chi_integrated_global = double(0);
            for (int i = 0; i < finest_level_dims[0]; i++)
            {
                // Use minus sign.
                chi_integrated_global -= scalar_dissipation_rate_full[i]*dx_finest[0];
            }
            
            const double* x_lo = d_grid_geometry->getXLower();
            const double* x_hi = d_grid_geometry->getXUpper();
            const double L_y = x_hi[1] - x_lo[1];
            chi_integrated_global *= L_y;
            
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << chi_integrated_global;
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        std::vector<std::string> quantity_names;
        std::vector<int> component_indices;
        std::vector<bool> use_derivative;
        std::vector<int> derivative_directions;
        
        quantity_names.push_back("MASS_DIFFUSIVITIES");
        component_indices.push_back(0); // Assuming effective mass diffusivities are same for all species.
        use_derivative.push_back(false);
        derivative_directions.push_back(-1);
        quantity_names.push_back("MASS_FRACTIONS");
        component_indices.push_back(0);
        use_derivative.push_back(true);
        derivative_directions.push_back(0);
        quantity_names.push_back("MASS_FRACTIONS");
        component_indices.push_back(1);
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
        component_indices.push_back(0); // Assuming effective mass diffusivities are same for all species.
        use_derivative.push_back(false);
        derivative_directions.push_back(-1);
        quantity_names.push_back("MASS_FRACTIONS");
        component_indices.push_back(0);
        use_derivative.push_back(true);
        derivative_directions.push_back(1);
        quantity_names.push_back("MASS_FRACTIONS");
        component_indices.push_back(1);
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
        component_indices.push_back(0); // Assuming effective mass diffusivities are same for all species.
        use_derivative.push_back(false);
        derivative_directions.push_back(-1);
        quantity_names.push_back("MASS_FRACTIONS");
        component_indices.push_back(0);
        use_derivative.push_back(true);
        derivative_directions.push_back(2);
        quantity_names.push_back("MASS_FRACTIONS");
        component_indices.push_back(1);
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
        
        /*
         * Output the scalar dissipation rate integral (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            std::vector<double> scalar_dissipation_rate_full(finest_level_dims[0], double(0));
            for (int i = 0; i < finest_level_dims[0]; i++)
            {
                scalar_dissipation_rate_full[i] =
                    scalar_dissipation_rate_part_1[i] + scalar_dissipation_rate_part_2[i] + scalar_dissipation_rate_part_3[i];
            }
            
            double chi_integrated_global = double(0);
            for (int i = 0; i < finest_level_dims[0]; i++)
            {
                // Use minus sign.
                chi_integrated_global -= scalar_dissipation_rate_full[i]*dx_finest[0];
            }
            
            const double* x_lo = d_grid_geometry->getXLower();
            const double* x_hi = d_grid_geometry->getXUpper();
            const double L_y = x_hi[1] - x_lo[1];
            const double L_z = x_hi[2] - x_lo[2];
            chi_integrated_global *= (L_y*L_z);
            
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << chi_integrated_global;
        }
    }
    
    if (mpi.getRank() == 0)
    {
        f_out.close();
    }
}


/*
 * Output scalar dissipation rate (can be negative) of between 2nd and 3rd species integrated to a file.
 */
void
RTIThreeLayersStatisticsUtilities::outputScalarDissipationRateSpeciesTwoAndThreeIntegrated(
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
    
    std::ofstream f_out;
    
    if (mpi.getRank() == 0)
    {
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
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
        flow_model_tmp,
        true);
    
    const std::vector<double>& dx_finest = MPI_helper_average.getFinestRefinedDomainGridSpacing();
    
    const hier::IntVector& finest_level_dims = MPI_helper_average.getFinestRefinedDomainNumberOfPoints();
    
    if (d_dim == tbox::Dimension(1))
    {
        std::vector<std::string> quantity_names;
        std::vector<int> component_indices;
        std::vector<bool> use_derivative;
        std::vector<int> derivative_directions;
        
        quantity_names.push_back("MASS_DIFFUSIVITIES");
        component_indices.push_back(0); // Assuming effective mass diffusivities are same for all species.
        use_derivative.push_back(false);
        derivative_directions.push_back(-1);
        quantity_names.push_back("MASS_FRACTIONS");
        component_indices.push_back(1);
        use_derivative.push_back(true);
        derivative_directions.push_back(0);
        quantity_names.push_back("MASS_FRACTIONS");
        component_indices.push_back(2);
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
        
        /*
         * Output the scalar dissipation rate integral (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            double chi_integrated_global = double(0);
            for (int i = 0; i < finest_level_dims[0]; i++)
            {
                chi_integrated_global -= scalar_dissipation_rate[i]*dx_finest[0]; // Use minus sign.
            }
            
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << chi_integrated_global;
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        std::vector<std::string> quantity_names;
        std::vector<int> component_indices;
        std::vector<bool> use_derivative;
        std::vector<int> derivative_directions;
        
        quantity_names.push_back("MASS_DIFFUSIVITIES");
        component_indices.push_back(0); // Assuming effective mass diffusivities are same for all species.
        use_derivative.push_back(false);
        derivative_directions.push_back(-1);
        quantity_names.push_back("MASS_FRACTIONS");
        component_indices.push_back(1);
        use_derivative.push_back(true);
        derivative_directions.push_back(0);
        quantity_names.push_back("MASS_FRACTIONS");
        component_indices.push_back(2);
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
        component_indices.push_back(0); // Assuming effective mass diffusivities are same for all species.
        use_derivative.push_back(false);
        derivative_directions.push_back(-1);
        quantity_names.push_back("MASS_FRACTIONS");
        component_indices.push_back(1);
        use_derivative.push_back(true);
        derivative_directions.push_back(1);
        quantity_names.push_back("MASS_FRACTIONS");
        component_indices.push_back(2);
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
        
        /*
         * Output the scalar dissipation rate integral (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            std::vector<double> scalar_dissipation_rate_full(finest_level_dims[0], double(0));
            for (int i = 0; i < finest_level_dims[0]; i++)
            {
                scalar_dissipation_rate_full[i] = scalar_dissipation_rate_part_1[i] + scalar_dissipation_rate_part_2[i];
            }
            
            double chi_integrated_global = double(0);
            for (int i = 0; i < finest_level_dims[0]; i++)
            {
                // Use minus sign.
                chi_integrated_global -= scalar_dissipation_rate_full[i]*dx_finest[0];
            }
            
            const double* x_lo = d_grid_geometry->getXLower();
            const double* x_hi = d_grid_geometry->getXUpper();
            const double L_y = x_hi[1] - x_lo[1];
            chi_integrated_global *= L_y;
            
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << chi_integrated_global;
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        std::vector<std::string> quantity_names;
        std::vector<int> component_indices;
        std::vector<bool> use_derivative;
        std::vector<int> derivative_directions;
        
        quantity_names.push_back("MASS_DIFFUSIVITIES");
        component_indices.push_back(0); // Assuming effective mass diffusivities are same for all species.
        use_derivative.push_back(false);
        derivative_directions.push_back(-1);
        quantity_names.push_back("MASS_FRACTIONS");
        component_indices.push_back(1);
        use_derivative.push_back(true);
        derivative_directions.push_back(0);
        quantity_names.push_back("MASS_FRACTIONS");
        component_indices.push_back(2);
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
        component_indices.push_back(0); // Assuming effective mass diffusivities are same for all species.
        use_derivative.push_back(false);
        derivative_directions.push_back(-1);
        quantity_names.push_back("MASS_FRACTIONS");
        component_indices.push_back(1);
        use_derivative.push_back(true);
        derivative_directions.push_back(1);
        quantity_names.push_back("MASS_FRACTIONS");
        component_indices.push_back(2);
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
        component_indices.push_back(0); // Assuming effective mass diffusivities are same for all species.
        use_derivative.push_back(false);
        derivative_directions.push_back(-1);
        quantity_names.push_back("MASS_FRACTIONS");
        component_indices.push_back(1);
        use_derivative.push_back(true);
        derivative_directions.push_back(2);
        quantity_names.push_back("MASS_FRACTIONS");
        component_indices.push_back(2);
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
        
        /*
         * Output the scalar dissipation rate integral (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            std::vector<double> scalar_dissipation_rate_full(finest_level_dims[0], double(0));
            for (int i = 0; i < finest_level_dims[0]; i++)
            {
                scalar_dissipation_rate_full[i] =
                    scalar_dissipation_rate_part_1[i] + scalar_dissipation_rate_part_2[i] + scalar_dissipation_rate_part_3[i];
            }
            
            double chi_integrated_global = double(0);
            for (int i = 0; i < finest_level_dims[0]; i++)
            {
                // Use minus sign.
                chi_integrated_global -= scalar_dissipation_rate_full[i]*dx_finest[0];
            }
            
            const double* x_lo = d_grid_geometry->getXLower();
            const double* x_hi = d_grid_geometry->getXUpper();
            const double L_y = x_hi[1] - x_lo[1];
            const double L_z = x_hi[2] - x_lo[2];
            chi_integrated_global *= (L_y*L_z);
            
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << chi_integrated_global;
        }
    }
    
    if (mpi.getRank() == 0)
    {
        f_out.close();
    }
}


/*
 * Output scalar dissipation rate (can be negative) of between 1st and 3rd species integrated to a file.
 */
void
RTIThreeLayersStatisticsUtilities::outputScalarDissipationRateSpeciesOneAndThreeIntegrated(
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
    
    std::ofstream f_out;
    
    if (mpi.getRank() == 0)
    {
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
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
        flow_model_tmp,
        true);
    
    const std::vector<double>& dx_finest = MPI_helper_average.getFinestRefinedDomainGridSpacing();
    
    const hier::IntVector& finest_level_dims = MPI_helper_average.getFinestRefinedDomainNumberOfPoints();
    
    if (d_dim == tbox::Dimension(1))
    {
        std::vector<std::string> quantity_names;
        std::vector<int> component_indices;
        std::vector<bool> use_derivative;
        std::vector<int> derivative_directions;
        
        quantity_names.push_back("MASS_DIFFUSIVITIES");
        component_indices.push_back(0); // Assuming effective mass diffusivities are same for all species.
        use_derivative.push_back(false);
        derivative_directions.push_back(-1);
        quantity_names.push_back("MASS_FRACTIONS");
        component_indices.push_back(0);
        use_derivative.push_back(true);
        derivative_directions.push_back(0);
        quantity_names.push_back("MASS_FRACTIONS");
        component_indices.push_back(2);
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
        
        /*
         * Output the scalar dissipation rate integral (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            double chi_integrated_global = double(0);
            for (int i = 0; i < finest_level_dims[0]; i++)
            {
                chi_integrated_global -= scalar_dissipation_rate[i]*dx_finest[0]; // Use minus sign.
            }
            
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << chi_integrated_global;
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        std::vector<std::string> quantity_names;
        std::vector<int> component_indices;
        std::vector<bool> use_derivative;
        std::vector<int> derivative_directions;
        
        quantity_names.push_back("MASS_DIFFUSIVITIES");
        component_indices.push_back(0); // Assuming effective mass diffusivities are same for all species.
        use_derivative.push_back(false);
        derivative_directions.push_back(-1);
        quantity_names.push_back("MASS_FRACTIONS");
        component_indices.push_back(0);
        use_derivative.push_back(true);
        derivative_directions.push_back(0);
        quantity_names.push_back("MASS_FRACTIONS");
        component_indices.push_back(2);
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
        component_indices.push_back(0); // Assuming effective mass diffusivities are same for all species.
        use_derivative.push_back(false);
        derivative_directions.push_back(-1);
        quantity_names.push_back("MASS_FRACTIONS");
        component_indices.push_back(0);
        use_derivative.push_back(true);
        derivative_directions.push_back(1);
        quantity_names.push_back("MASS_FRACTIONS");
        component_indices.push_back(2);
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
        
        /*
         * Output the scalar dissipation rate integral (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            std::vector<double> scalar_dissipation_rate_full(finest_level_dims[0], double(0));
            for (int i = 0; i < finest_level_dims[0]; i++)
            {
                scalar_dissipation_rate_full[i] = scalar_dissipation_rate_part_1[i] + scalar_dissipation_rate_part_2[i];
            }
            
            double chi_integrated_global = double(0);
            for (int i = 0; i < finest_level_dims[0]; i++)
            {
                // Use minus sign.
                chi_integrated_global -= scalar_dissipation_rate_full[i]*dx_finest[0];
            }
            
            const double* x_lo = d_grid_geometry->getXLower();
            const double* x_hi = d_grid_geometry->getXUpper();
            const double L_y = x_hi[1] - x_lo[1];
            chi_integrated_global *= L_y;
            
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << chi_integrated_global;
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        std::vector<std::string> quantity_names;
        std::vector<int> component_indices;
        std::vector<bool> use_derivative;
        std::vector<int> derivative_directions;
        
        quantity_names.push_back("MASS_DIFFUSIVITIES");
        component_indices.push_back(0); // Assuming effective mass diffusivities are same for all species.
        use_derivative.push_back(false);
        derivative_directions.push_back(-1);
        quantity_names.push_back("MASS_FRACTIONS");
        component_indices.push_back(0);
        use_derivative.push_back(true);
        derivative_directions.push_back(0);
        quantity_names.push_back("MASS_FRACTIONS");
        component_indices.push_back(2);
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
        component_indices.push_back(0); // Assuming effective mass diffusivities are same for all species.
        use_derivative.push_back(false);
        derivative_directions.push_back(-1);
        quantity_names.push_back("MASS_FRACTIONS");
        component_indices.push_back(0);
        use_derivative.push_back(true);
        derivative_directions.push_back(1);
        quantity_names.push_back("MASS_FRACTIONS");
        component_indices.push_back(2);
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
        component_indices.push_back(0); // Assuming effective mass diffusivities are same for all species.
        use_derivative.push_back(false);
        derivative_directions.push_back(-1);
        quantity_names.push_back("MASS_FRACTIONS");
        component_indices.push_back(0);
        use_derivative.push_back(true);
        derivative_directions.push_back(2);
        quantity_names.push_back("MASS_FRACTIONS");
        component_indices.push_back(2);
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
        
        /*
         * Output the scalar dissipation rate integral (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            std::vector<double> scalar_dissipation_rate_full(finest_level_dims[0], double(0));
            for (int i = 0; i < finest_level_dims[0]; i++)
            {
                scalar_dissipation_rate_full[i] =
                    scalar_dissipation_rate_part_1[i] + scalar_dissipation_rate_part_2[i] + scalar_dissipation_rate_part_3[i];
            }
            
            double chi_integrated_global = double(0);
            for (int i = 0; i < finest_level_dims[0]; i++)
            {
                // Use minus sign.
                chi_integrated_global -= scalar_dissipation_rate_full[i]*dx_finest[0];
            }
            
            const double* x_lo = d_grid_geometry->getXLower();
            const double* x_hi = d_grid_geometry->getXUpper();
            const double L_y = x_hi[1] - x_lo[1];
            const double L_z = x_hi[2] - x_lo[2];
            chi_integrated_global *= (L_y*L_z);
            
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << chi_integrated_global;
        }
    }
    
    if (mpi.getRank() == 0)
    {
        f_out.close();
    }
}


/*
 * Output maximum value of species 1 in the domain to a file.
 */
void
RTIThreeLayersStatisticsUtilities::outputMaxMassFractionSpeciesOne(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 3)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'Y_1_MAX' can be computed with three species only."
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
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
    }
    
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
    FlowModelMPIHelperMaxMin MPI_helper_max_min = FlowModelMPIHelperMaxMin(
        "MPI_helper_max_min",
        d_dim,
        d_grid_geometry,
        patch_hierarchy,
        flow_model_tmp);
    
    const double max_value_Y1 = MPI_helper_max_min.getMaxQuantity(
            "MASS_FRACTIONS",
            0,
            data_context);
    
    /*
     * Compute and output the value (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << max_value_Y1;
        
        f_out.close();
    }
}


/*
 * Output maximum value of species 2 in the domain to a file.
 */
void
RTIThreeLayersStatisticsUtilities::outputMaxMassFractionSpeciesTwo(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 3)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'Y_2_MAX' can be computed with three species only."
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
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
    }
    
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
    FlowModelMPIHelperMaxMin MPI_helper_max_min = FlowModelMPIHelperMaxMin(
        "MPI_helper_max_min",
        d_dim,
        d_grid_geometry,
        patch_hierarchy,
        flow_model_tmp);
    
    const double max_value_Y2 = MPI_helper_max_min.getMaxQuantity(
            "MASS_FRACTIONS",
            1,
            data_context);
    
    /*
     * Compute and output the value (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << max_value_Y2;
        
        f_out.close();
    }
}


/*
 * Output maximum value of species 3 in the domain to a file.
 */
void
RTIThreeLayersStatisticsUtilities::outputMaxMassFractionSpeciesThree(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 3)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'Y_3_MAX' can be computed with three species only."
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
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
    }
    
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
    FlowModelMPIHelperMaxMin MPI_helper_max_min = FlowModelMPIHelperMaxMin(
        "MPI_helper_max_min",
        d_dim,
        d_grid_geometry,
        patch_hierarchy,
        flow_model_tmp);
    
    const double max_value_Y3 = MPI_helper_max_min.getMaxQuantity(
            "MASS_FRACTIONS",
            2,
            data_context);
    
    /*
     * Compute and output the value (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << max_value_Y3;
        
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
            
            if (statistical_quantity_key == "INTERFACE_AVG_MIN_S1")
            {
                f_out << "\t" << "INTERFACE_AVG_MIN_S1 ";
            }
            else if (statistical_quantity_key == "INTERFACE_AVG_MAX_S1")
            {
                f_out << "\t" << "INTERFACE_AVG_MAX_S1 ";
            }
            else if (statistical_quantity_key == "INTERFACE_AVG_MIN_S2")
            {
                f_out << "\t" << "INTERFACE_AVG_MIN_S2 ";
            }
            else if (statistical_quantity_key == "INTERFACE_AVG_MAX_S2")
            {
                f_out << "\t" << "INTERFACE_AVG_MAX_S2 ";
            }
            else if (statistical_quantity_key == "INTERFACE_AVG_MIN_S3")
            {
                f_out << "\t" << "INTERFACE_AVG_MIN_S3 ";
            }
            else if (statistical_quantity_key == "INTERFACE_AVG_MAX_S3")
            {
                f_out << "\t" << "INTERFACE_AVG_MAX_S3 ";
            }
            else if (statistical_quantity_key == "INTERFACE_MIN_S1")
            {
                f_out << "\t" << "INTERFACE_MIN_S1     ";
            }
            else if (statistical_quantity_key == "INTERFACE_MAX_S1")
            {
                f_out << "\t" << "INTERFACE_MAX_S1     ";
            }
            else if (statistical_quantity_key == "INTERFACE_MIN_S2")
            {
                f_out << "\t" << "INTERFACE_MIN_S2     ";
            }
            else if (statistical_quantity_key == "INTERFACE_MAX_S2")
            {
                f_out << "\t" << "INTERFACE_MAX_S2     ";
            }
            else if (statistical_quantity_key == "INTERFACE_MIN_S3")
            {
                f_out << "\t" << "INTERFACE_MIN_S3     ";
            }
            else if (statistical_quantity_key == "INTERFACE_MAX_S3")
            {
                f_out << "\t" << "INTERFACE_MAX_S3     ";
            }
            else if (statistical_quantity_key == "MIXING_WIDTH_X_S1")
            {
                f_out << "\t" << "MIXING_WIDTH_X_S1    ";
            }
            else if (statistical_quantity_key == "MIXING_WIDTH_X_S2")
            {
                f_out << "\t" << "MIXING_WIDTH_X_S2    ";
            }
            else if (statistical_quantity_key == "MIXING_WIDTH_X_S3")
            {
                f_out << "\t" << "MIXING_WIDTH_X_S3    ";
            }
            else if (statistical_quantity_key == "MIXEDNESS_X_S1")
            {
                f_out << "\t" << "MIXEDNESS_X_S1       ";
            }
            else if (statistical_quantity_key == "MIXEDNESS_X_S2")
            {
                f_out << "\t" << "MIXEDNESS_X_S2       ";
            }
            else if (statistical_quantity_key == "MIXEDNESS_X_S3")
            {
                f_out << "\t" << "MIXEDNESS_X_S3       ";
            }
            else if (statistical_quantity_key == "MIX_METRIC_INT_S12")
            {
                f_out << "\t" << "MIX_METRIC_INT_S12    ";
            }
            else if (statistical_quantity_key == "MIX_METRIC_INT_S23")
            {
                f_out << "\t" << "MIX_METRIC_INT_S23    ";
            }
            else if (statistical_quantity_key == "MIX_METRIC_INT_S13")
            {
                f_out << "\t" << "MIX_METRIC_INT_S13    ";
            }
            // else if (statistical_quantity_key == "R11_MEAN_INHOMO_X")
            // {
            //     f_out << "\t" << "R11_MEAN_INHOMO_X    ";
            // }
            // else if (statistical_quantity_key == "R22_MEAN_INHOMO_X")
            // {
            //     f_out << "\t" << "R22_MEAN_INHOMO_X    ";
            // }
            else if (statistical_quantity_key == "ENSTROPHY_INT")
            {
                f_out << "\t" << "ENSTROPHY_INT        ";
            }
            else if (statistical_quantity_key == "TKE_INT_HOMO_Y")
            {
                f_out << "\t" << "TKE_INT_HOMO_Y       ";
            }
            else if (statistical_quantity_key == "SCAL_DISS_RAT_INT_S12")
            {
                f_out << "\t" << "SCAL_DISS_RAT_INT_S12";
            }
            else if (statistical_quantity_key == "SCAL_DISS_RAT_INT_S23")
            {
                f_out << "\t" << "SCAL_DISS_RAT_INT_S23";
            }
            else if (statistical_quantity_key == "SCAL_DISS_RAT_INT_S13")
            {
                f_out << "\t" << "SCAL_DISS_RAT_INT_S13";
            }
            else if (statistical_quantity_key == "Y_1_MAX")
            {
                f_out << "\t" << "Y_1_MAX              ";
            }
            else if (statistical_quantity_key == "Y_2_MAX")
            {
                f_out << "\t" << "Y_2_MAX              ";
            }
            else if (statistical_quantity_key == "Y_3_MAX")
            {
                f_out << "\t" << "Y_3_MAX              ";
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
    NULL_USE(patch_hierarchy);
    NULL_USE(data_context);
    NULL_USE(statistics_data_time);
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
    NULL_USE(output_time);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    HAMERS_SHARED_PTR<RTIThreeLayersStatisticsUtilities> rti_three_layers_statistics_utilities(
        new RTIThreeLayersStatisticsUtilities(
            "RTI three layers statistics utilities",
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
        
        if (statistical_quantity_key == "INTERFACE_AVG_MIN_S1")
        {
            rti_three_layers_statistics_utilities->outputAveragedInterfaceMinSpeciesOneInXDirection(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "INTERFACE_AVG_MAX_S1")
        {
            rti_three_layers_statistics_utilities->outputAveragedInterfaceMaxSpeciesOneInXDirection(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "INTERFACE_AVG_MIN_S2")
        {
            rti_three_layers_statistics_utilities->outputAveragedInterfaceMinSpeciesTwoInXDirection(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "INTERFACE_AVG_MAX_S2")
        {
            rti_three_layers_statistics_utilities->outputAveragedInterfaceMaxSpeciesTwoInXDirection(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "INTERFACE_AVG_MIN_S3")
        {
            rti_three_layers_statistics_utilities->outputAveragedInterfaceMinSpeciesThreeInXDirection(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "INTERFACE_AVG_MAX_S3")
        {
            rti_three_layers_statistics_utilities->outputAveragedInterfaceMaxSpeciesThreeInXDirection(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "INTERFACE_MIN_S1")
        {
            rti_three_layers_statistics_utilities->outputInterfaceMinSpeciesOneInXDirection(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "INTERFACE_MAX_S1")
        {
            rti_three_layers_statistics_utilities->outputInterfaceMaxSpeciesOneInXDirection(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "INTERFACE_MIN_S2")
        {
            rti_three_layers_statistics_utilities->outputInterfaceMinSpeciesTwoInXDirection(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "INTERFACE_MAX_S2")
        {
            rti_three_layers_statistics_utilities->outputInterfaceMaxSpeciesTwoInXDirection(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "INTERFACE_MIN_S3")
        {
            rti_three_layers_statistics_utilities->outputInterfaceMinSpeciesThreeInXDirection(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "INTERFACE_MAX_S3")
        {
            rti_three_layers_statistics_utilities->outputInterfaceMaxSpeciesThreeInXDirection(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "MIXING_WIDTH_X_S1")
        {
            rti_three_layers_statistics_utilities->outputMixingWidthSpeciesOneInXDirection(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "MIXING_WIDTH_X_S2")
        {
            rti_three_layers_statistics_utilities->outputMixingWidthSpeciesTwoInXDirection(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "MIXING_WIDTH_X_S3")
        {
            rti_three_layers_statistics_utilities->outputMixingWidthSpeciesThreeInXDirection(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "MIXEDNESS_X_S1")
        {
            rti_three_layers_statistics_utilities->outputMixednessSpeciesOneInXDirection(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "MIXEDNESS_X_S2")
        {
            rti_three_layers_statistics_utilities->outputMixednessSpeciesTwoInXDirection(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "MIXEDNESS_X_S3")
        {
            rti_three_layers_statistics_utilities->outputMixednessSpeciesThreeInXDirection(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "MIX_METRIC_INT_S12")
        {
            rti_three_layers_statistics_utilities->outputMixingMetricSpeciesOneAndTwoIntegrated(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "MIX_METRIC_INT_S23")
        {
            rti_three_layers_statistics_utilities->outputMixingMetricSpeciesTwoAndThreeIntegrated(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "MIX_METRIC_INT_S13")
        {
            rti_three_layers_statistics_utilities->outputMixingMetricSpeciesOneAndThreeIntegrated(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        // else if (statistical_quantity_key == "R11_MEAN_INHOMO_X")
        // {
        //     rti_three_layers_statistics_utilities->outputReynoldsNormalStressXWithInhomogeneousXDirection(
        //         stat_dump_filename,
        //         patch_hierarchy,
        //         data_context);
        // }
        // else if (statistical_quantity_key == "R22_MEAN_INHOMO_X")
        // {
        //     rti_three_layers_statistics_utilities->outputReynoldsNormalStressYWithInhomogeneousXDirection(
        //         stat_dump_filename,
        //         patch_hierarchy,
        //         data_context);
        // }
        else if (statistical_quantity_key == "ENSTROPHY_INT")
        {
            rti_three_layers_statistics_utilities->outputEnstrophyIntegrated(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "TKE_INT_HOMO_Y")
        {
            rti_three_layers_statistics_utilities->outputTKEIntegratedWithHomogeneityInYDirection(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "SCAL_DISS_RAT_INT_S12")
        {
            rti_three_layers_statistics_utilities->outputScalarDissipationRateSpeciesOneAndTwoIntegrated(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "SCAL_DISS_RAT_INT_S23")
        {
            rti_three_layers_statistics_utilities->outputScalarDissipationRateSpeciesTwoAndThreeIntegrated(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "SCAL_DISS_RAT_INT_S13")
        {
            rti_three_layers_statistics_utilities->outputScalarDissipationRateSpeciesOneAndThreeIntegrated(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "Y_1_MAX")
        {
            rti_three_layers_statistics_utilities->outputMaxMassFractionSpeciesOne(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "Y_2_MAX")
        {
            rti_three_layers_statistics_utilities->outputMaxMassFractionSpeciesTwo(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "Y_3_MAX")
        {
            rti_three_layers_statistics_utilities->outputMaxMassFractionSpeciesThree(
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
