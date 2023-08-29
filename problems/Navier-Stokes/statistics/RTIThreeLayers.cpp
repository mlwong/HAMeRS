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
            Start new stats calculations for 3 layers RTI
        */
        
        void
        outputInterfaceMinSpeciesOneInXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;

        void
        outputInterfaceMaxSpeciesOneInXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;

        void
        outputInterfaceMinSpeciesTwoInXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;

        void
        outputInterfaceMaxSpeciesTwoInXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;

        void
        outputInterfaceMinSpeciesThreeInXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;

        void
        outputInterfaceMaxSpeciesThreeInXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;

        void
        outputMaxMassFractionSpeciesTwo(
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
 * Output minimum interface location in x-direction to a file.
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
 * Output maximum interface location in x-direction to a file.
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

//Second Species

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

//

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
            
            if (statistical_quantity_key == "INTERFACE_MIN_S1")
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
            else if (statistical_quantity_key == "Y_2_MAX")
            {
                f_out << "\t" << "Y_2_MAX              ";
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
        
        if (statistical_quantity_key == "INTERFACE_MIN_S1")
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
        else if (statistical_quantity_key == "Y_2_MAX")
        {
            rti_three_layers_statistics_utilities->outputMaxMassFractionSpeciesTwo(
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
