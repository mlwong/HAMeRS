#include "flow/flow_models/four-eqn_conservative/FlowModelStatisticsUtilitiesFourEqnConservative.hpp"

#include "extn/patch_hierarchies/ExtendedFlattenedHierarchy.hpp"
#include "flow/flow_models/MPI_helpers/FlowModelMPIHelperAverage.hpp"
#include "flow/flow_models/MPI_helpers/FlowModelMPIHelperCentroid.hpp"
#include "flow/flow_models/MPI_helpers/FlowModelMPIHelperMaxMin.hpp"
#include "util/MPI_helpers/MPIHelperGrid.hpp"

#include <fstream>


class SBIStatisticsUtilities
{
    public:
        SBIStatisticsUtilities(
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
         * Output mixing width in x-direction to a file.
         */
        void
        outputMixingWidthInXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Output centroid in x-directioin to a file.
         */
        void
        outputCentroidInXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Output minimum interface location in x-direction to a file.
         */
        void
        outputInterfaceMinInXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Output maximum interface location in x-direction to a file.
         */
        void
        outputInterfaceMaxInXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Output minimum interface location in y-direction to a file.
         */
        void
        outputInterfaceMinInYDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Output maximum interface location in y-direction to a file.
         */
        void
        outputInterfaceMaxInYDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Output minimum interface location in z-direction to a file.
         */
        void
        outputInterfaceMinInZDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Output maximum interface location in z-direction to a file.
         */
        void
        outputInterfaceMaxInZDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Output circulation to a file.
         */
        void
        outputCirculation(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Output scalar dissipation rate of first species integrated to a file.
         */
        void
        outputScalarDissipationRateIntegrated(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Output enstrophy integrated to a file.
         */
        void
        outputEnstrophyIntegrated(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Output number of cells to a file.
         */
        void
        outputNumberOfCells(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Output weighted number of cells to a file.
         */
        void
        outputWeightedNumberOfCells(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
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
 * Output mixing width in x-direction to a file.
 */
void
SBIStatisticsUtilities::outputMixingWidthInXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
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
 * Output centroid in x-direction to a file.
 */
void
SBIStatisticsUtilities::outputCentroidInXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'CENTROID_X' can be computed with two species only."
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
    
    FlowModelMPIHelperCentroid MPI_helper_centroid = FlowModelMPIHelperCentroid(
        "MPI_helper_centroid",
        d_dim,
        d_grid_geometry,
        patch_hierarchy,
        flow_model_tmp);
    
    const double x_c = MPI_helper_centroid.getCentroidInXDirection(
            "MASS_FRACTIONS",
            0,
            data_context);
    
    /*
     * Compute and output the centroid (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << x_c;
        
        f_out.close();
    }
}


/*
 * Output minimum interface location in x-direction to a file.
 */
void
SBIStatisticsUtilities::outputInterfaceMinInXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
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
            double(0.01),
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
 * Output maximum interface location in x-direction to a file.
 */
void
SBIStatisticsUtilities::outputInterfaceMaxInXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
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
            double(0.01),
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
 * Output minimum interface location in y-direction to a file.
 */
void
SBIStatisticsUtilities::outputInterfaceMinInYDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'INTERFACE_MIN_Y' can be computed with two species only."
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
    
    const double interface_y_min_global = MPI_helper_max_min.getMinLocationWithinQuantityBoundsInYDirection(
            "MASS_FRACTIONS",
            0,
            data_context,
            double(0.01),
            double(0.99));
    
    /*
     * Compute and output the value (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << interface_y_min_global;
        
        f_out.close();
    }
}


/*
 * Output maximum interface location in y-direction to a file.
 */
void
SBIStatisticsUtilities::outputInterfaceMaxInYDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'INTERFACE_MAX_Y' can be computed with two species only."
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
    
    const double interface_y_max_global = MPI_helper_max_min.getMaxLocationWithinQuantityBoundsInYDirection(
            "MASS_FRACTIONS",
            0,
            data_context,
            double(0.01),
            double(0.99));
    
    /*
     * Compute and output the value (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << interface_y_max_global;
        
        f_out.close();
    }
}


/*
 * Output minimum interface location in z-direction to a file.
 */
void
SBIStatisticsUtilities::outputInterfaceMinInZDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'INTERFACE_MIN_Z' can be computed with two species only."
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
    
    const double interface_z_min_global = MPI_helper_max_min.getMinLocationWithinQuantityBoundsInZDirection(
            "MASS_FRACTIONS",
            0,
            data_context,
            double(0.01),
            double(0.99));
    
    /*
     * Compute and output the value (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << interface_z_min_global;
        
        f_out.close();
    }
}


/*
 * Output maximum interface location in z-direction to a file.
 */
void
SBIStatisticsUtilities::outputInterfaceMaxInZDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'INTERFACE_MAX_Z' can be computed with two species only."
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
    
    const double interface_z_max_global = MPI_helper_max_min.getMaxLocationWithinQuantityBoundsInZDirection(
            "MASS_FRACTIONS",
            0,
            data_context,
            double(0.01),
            double(0.99));
    
    /*
     * Compute and output the value (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << interface_z_max_global;
        
        f_out.close();
    }
}


/*
 * Output circulation to a file.
 */
void
SBIStatisticsUtilities::outputCirculation(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
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
    
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
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
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    HAMERS_SHARED_PTR<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *patch_hierarchy,
            0,
            num_levels - 1));
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no 'CIRCULATION' for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        double Circulation_local = 0.0;
        double Circulation_global = 0.0;
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch geometry.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::IntVector patch_dims = patch_box.numberCells();
                
                const int patch_dim_0 = patch_dims[0];
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and velocity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", num_ghosts));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->allocateMemoryForDerivedCellData();
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointer to velocity data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_velocity =
                    flow_model_tmp->getCellData("VELOCITY");
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_velocity = data_velocity->getGhostCellWidth();
                
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(num_ghosts_velocity >= hier::IntVector::getOne(d_dim)*d_num_ghosts_derivative);
#endif
                
                double Circulation_to_add = 0.0;
                
                /*
                 * Initialize cell data for velocity derivatives and get pointers to the derivatives.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_velocity_derivatives(
                    new pdat::CellData<double>(patch_box, d_dim.getValue()*d_dim.getValue() - d_dim.getValue(),
                        hier::IntVector::getZero(d_dim)));
                
                double* dudy = data_velocity_derivatives->getPointer(0);
                double* dvdx = data_velocity_derivatives->getPointer(1);
                
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
                    
                    // Compute dudy.
                    derivative_first_order_y->computeDerivative(
                        data_velocity_derivatives,
                        data_velocity,
                        dx[1],
                        patch_visible_box,
                        0,
                        0);
                    
                    // Compute dvdx.
                    derivative_first_order_x->computeDerivative(
                        data_velocity_derivatives,
                        data_velocity,
                        dx[0],
                        patch_visible_box,
                        1,
                        1);
                    
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
                            
                            // Compute linear indices.
                            const int idx = (relative_idx_lo_0 + i) +
                                (relative_idx_lo_1 + j)*patch_dim_0;
                            
                            const double omega = dvdx[idx] - dudy[idx];
                            Circulation_to_add += omega/((double) n_overlapped);
                        }
                    }
                }
                
                Circulation_to_add = Circulation_to_add*dx[0]*dx[1];
                Circulation_local += Circulation_to_add;
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global integral.
         */
        
        mpi.Reduce(
            &Circulation_local,
            &Circulation_global,
            1,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        /*
         * Output the circulation (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << fabs(Circulation_global);
        }
        
        if (mpi.getRank() == 0)
        {
            f_out.close();
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        double Circulation_y_local = 0.0;
        double Circulation_y_global = 0.0;
        
        double Circulation_z_local = 0.0;
        double Circulation_z_global = 0.0;
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch geometry.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::IntVector patch_dims = patch_box.numberCells();
                
                const int patch_dim_0 = patch_dims[0];
                const int patch_dim_1 = patch_dims[1];
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and velocity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", num_ghosts));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->allocateMemoryForDerivedCellData();
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointer to velocity data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_velocity = flow_model_tmp->getCellData("VELOCITY");
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_velocity = data_velocity->getGhostCellWidth();
                
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(num_ghosts_velocity >= hier::IntVector::getOne(d_dim)*d_num_ghosts_derivative);
#endif
                
                double Circulation_y_to_add = 0.0;
                double Circulation_z_to_add = 0.0;
                
                /*
                 * Initialize cell data for velocity derivatives and get pointers to the derivatives.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_velocity_derivatives(
                    new pdat::CellData<double>(patch_box, d_dim.getValue()*d_dim.getValue() - d_dim.getValue(),
                        hier::IntVector::getZero(d_dim)));
                
                double* dudy = data_velocity_derivatives->getPointer(0);
                double* dudz = data_velocity_derivatives->getPointer(1);
                double* dvdx = data_velocity_derivatives->getPointer(2);
                // double* dvdz = data_velocity_derivatives->getPointer(3);
                double* dwdx = data_velocity_derivatives->getPointer(4);
                // double* dwdy = data_velocity_derivatives->getPointer(5);
                
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
                    
                    // Compute dudy.
                    derivative_first_order_y->computeDerivative(
                        data_velocity_derivatives,
                        data_velocity,
                        dx[1],
                        patch_visible_box,
                        0,
                        0);
                    
                    // Compute dudz.
                    derivative_first_order_z->computeDerivative(
                        data_velocity_derivatives,
                        data_velocity,
                        dx[2],
                        patch_visible_box,
                        1,
                        0);
                    
                    // Compute dvdx.
                    derivative_first_order_x->computeDerivative(
                        data_velocity_derivatives,
                        data_velocity,
                        dx[0],
                        patch_visible_box,
                        2,
                        1);
                    
                    // Compute dvdz.
                    derivative_first_order_z->computeDerivative(
                        data_velocity_derivatives,
                        data_velocity,
                        dx[2],
                        patch_visible_box,
                        3,
                        1);
                    
                    // Compute dwdx.
                    derivative_first_order_x->computeDerivative(
                        data_velocity_derivatives,
                        data_velocity,
                        dx[0],
                        patch_visible_box,
                        4,
                        2);
                    
                    // Compute dwdy.
                    derivative_first_order_y->computeDerivative(
                        data_velocity_derivatives,
                        data_velocity,
                        dx[1],
                        patch_visible_box,
                        5,
                        2);
                    
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
                                
                                // Compute the linear indices.
                                const int idx = (relative_idx_lo_0 + i) +
                                    (relative_idx_lo_1 + j)*patch_dim_0 +
                                    (relative_idx_lo_2 + k)*patch_dim_0*
                                        patch_dim_1;
                                
                                const double omega_y = dudz[idx] - dwdx[idx];
                                const double omega_z = dvdx[idx] - dudy[idx];
                                Circulation_y_to_add += omega_y/((double) n_overlapped);
                                Circulation_z_to_add += omega_z/((double) n_overlapped);
                            }
                        }
                    }
                }
                
                Circulation_y_to_add = Circulation_y_to_add*dx[0]*dx[1]*dx[2];
                Circulation_z_to_add = Circulation_z_to_add*dx[0]*dx[1]*dx[2];
                
                Circulation_y_local += Circulation_y_to_add;
                Circulation_z_local += Circulation_z_to_add;
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global integral.
         */
        
        mpi.Reduce(
            &Circulation_y_local,
            &Circulation_y_global,
            1,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        mpi.Reduce(
            &Circulation_z_local,
            &Circulation_z_global,
            1,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        /*
         * Output the circulation (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << 0.5*(fabs(Circulation_y_global) + fabs(Circulation_z_global));
        }
        
        if (mpi.getRank() == 0)
        {
            f_out.close();
        }
    }
}


/*
 * Output scalar dissipation rate of first species integrated to a file.
 */
void
SBIStatisticsUtilities::outputScalarDissipationRateIntegrated(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
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
        
        /*
         * Output the scalar dissipation rate integral (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            double Chi_integrated_global = double(0);
            for (int i = 0; i < finest_level_dims[0]; i++)
            {
                Chi_integrated_global += scalar_dissipation_rate[i]*dx_finest[0];
            }
            
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << Chi_integrated_global;
        }
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
            
            double Chi_integrated_global = double(0);
            for (int i = 0; i < finest_level_dims[0]; i++)
            {
                Chi_integrated_global += scalar_dissipation_rate_full[i]*dx_finest[0];
            }
            
            const double* x_lo = d_grid_geometry->getXLower();
            const double* x_hi = d_grid_geometry->getXUpper();
            const double L_y = x_hi[1] - x_lo[1];
            Chi_integrated_global *= L_y;
            
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << Chi_integrated_global;
        }
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
            
            double Chi_integrated_global = double(0);
            for (int i = 0; i < finest_level_dims[0]; i++)
            {
                Chi_integrated_global += scalar_dissipation_rate_full[i]*dx_finest[0];
            }
            
            const double* x_lo = d_grid_geometry->getXLower();
            const double* x_hi = d_grid_geometry->getXUpper();
            const double L_y = x_hi[1] - x_lo[1];
            const double L_z = x_hi[2] - x_lo[2];
            Chi_integrated_global *= (L_y*L_z);
            
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << Chi_integrated_global;
        }
    }
    
    if (mpi.getRank() == 0)
    {
        f_out.close();
    }
}


/*
 * Output enstrophy integrated to a file.
 */
void
SBIStatisticsUtilities::outputEnstrophyIntegrated(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
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
 * Output number of cells to a file.
 */
void
SBIStatisticsUtilities::outputNumberOfCells(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
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
    
    MPIHelperGrid MPI_helper_num_cells = MPIHelperGrid(
        "MPI_helper_num_cells",
        d_dim,
        d_grid_geometry,
        patch_hierarchy);
    
    const double num_cells_global = MPI_helper_num_cells.getNumberOfCells();
    
    /*
     * Output the number of cells (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << num_cells_global;
        
        f_out.close();
    }
}


/*
 * Output weighted number of cells to a file.
 */
void
SBIStatisticsUtilities::outputWeightedNumberOfCells(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
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
    
    MPIHelperGrid MPI_helper_num_cells = MPIHelperGrid(
        "MPI_helper_num_cells",
        d_dim,
        d_grid_geometry,
        patch_hierarchy);
    
    const double weighted_num_cells_global = MPI_helper_num_cells.getWeightedNumberOfCells();
    
    /*
     * Output the number of cells (only done by process 0).
     */
    
    if (mpi.getRank() == 0)
    {
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << weighted_num_cells_global;
        
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
            else if (statistical_quantity_key == "CENTROID_X")
            {
                f_out << "\t" << "CENTROID_X           ";
            }
            else if (statistical_quantity_key == "INTERFACE_MIN_X")
            {
                f_out << "\t" << "INTERFACE_MIN_X      ";
            }
            else if (statistical_quantity_key == "INTERFACE_MAX_X")
            {
                f_out << "\t" << "INTERFACE_MAX_X      ";
            }
            else if (statistical_quantity_key == "INTERFACE_MIN_Y")
            {
                f_out << "\t" << "INTERFACE_MIN_Y      ";
            }
            else if (statistical_quantity_key == "INTERFACE_MAX_Y")
            {
                f_out << "\t" << "INTERFACE_MAX_Y      ";
            }
            else if (statistical_quantity_key == "INTERFACE_MIN_Z")
            {
                f_out << "\t" << "INTERFACE_MIN_Z      ";
            }
            else if (statistical_quantity_key == "INTERFACE_MAX_Z")
            {
                f_out << "\t" << "INTERFACE_MAX_Z      ";
            }
            else if (statistical_quantity_key == "CIRCULATION")
            {
                f_out << "\t" << "CIRCULATION          ";
            }
            else if (statistical_quantity_key == "ENSTROPHY_INT")
            {
                f_out << "\t" << "ENSTROPHY_INT        ";
            }
            else if (statistical_quantity_key == "SCAL_DISS_RAT_INT")
            {
                f_out << "\t" << "SCAL_DISS_RAT_INT    ";
            }
            else if (statistical_quantity_key == "NUM_CELLS")
            {
                f_out << "\t" << "NUM_CELLS            ";
            }
            else if (statistical_quantity_key == "WEIGHTED_NUM_CELLS")
            {
                f_out << "\t" << "WEIGHTED_NUM_CELLS   ";
            }
        }
        
        f_out.close();
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
    NULL_USE(output_time);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    HAMERS_SHARED_PTR<SBIStatisticsUtilities> sbi_statistics_utilities(
        new SBIStatisticsUtilities(
            "SBI statistics utilities",
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
        
        if (statistical_quantity_key == "MIXING_WIDTH_X")
        {
            sbi_statistics_utilities->outputMixingWidthInXDirection(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "CENTROID_X")
        {
            sbi_statistics_utilities->outputCentroidInXDirection(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "INTERFACE_MIN_X")
        {
            sbi_statistics_utilities->outputInterfaceMinInXDirection(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "INTERFACE_MAX_X")
        {
            sbi_statistics_utilities->outputInterfaceMaxInXDirection(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "INTERFACE_MIN_Y")
        {
            sbi_statistics_utilities->outputInterfaceMinInYDirection(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "INTERFACE_MAX_Y")
        {
            sbi_statistics_utilities->outputInterfaceMaxInYDirection(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "INTERFACE_MIN_Z")
        {
            sbi_statistics_utilities->outputInterfaceMinInZDirection(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "INTERFACE_MAX_Z")
        {
            sbi_statistics_utilities->outputInterfaceMaxInZDirection(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "CIRCULATION")
        {
            sbi_statistics_utilities->outputCirculation(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "ENSTROPHY_INT")
        {
            sbi_statistics_utilities->outputEnstrophyIntegrated(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "SCAL_DISS_RAT_INT")
        {
            sbi_statistics_utilities->outputScalarDissipationRateIntegrated(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "NUM_CELLS")
        {
            sbi_statistics_utilities->outputNumberOfCells(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "WEIGHTED_NUM_CELLS")
        {
            sbi_statistics_utilities->outputWeightedNumberOfCells(
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
