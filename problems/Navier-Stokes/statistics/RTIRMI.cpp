#include "flow/flow_models/four-eqn_conservative/FlowModelStatisticsUtilitiesFourEqnConservative.hpp"

#include "extn/patch_hierarchies/ExtendedFlattenedHierarchy.hpp"
#include "flow/flow_models/MPI_helpers/FlowModelMPIHelperAverage.hpp"
#include "flow/flow_models/MPI_helpers/FlowModelMPIHelperCorrelation.hpp"
#include "flow/flow_models/MPI_helpers/FlowModelMPIHelperMaxMin.hpp"
#include "util/MPI_helpers/MPIHelperGrid.hpp"

#include <fstream>

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
         * Output TKE integrated with assumed homogeneity in y-direction to a file.
         */
        void
        outputTKEIntegratedWithHomogeneityInYDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output TKE integrated with assumed homogeneity in yz-plane to a file.
         */
        void
        outputTKEIntegratedWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output TKE in x-direction integrated with assumed homogeneity in yz-plane to a file.
         */
        void
        outputTKEInXDirectionIntegratedWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output TKE in y-direction integrated with assumed homogeneity in yz-plane to a file.
         */
        void
        outputTKEInYDirectionIntegratedWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output TKE in z-direction integrated with assumed homogeneity in yz-plane to a file.
         */
        void
        outputTKEInZDirectionIntegratedWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output mean velocity associated with turbulent mass flux component in x-direction with assumed homogeneity
         * in y-direction (2D) or * yz-plane (3D) to a file.
         */
        void
        outputTurbMassFluxVelocityXWithInhomogeneousXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output density-specific volume covariance with assumed homogeneity in y-direction (2D) or yz-plane (3D)
         * to a file.
         */
        void
        outputDensitySpecificVolumeCovarianceWithInhomogeneousXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output Boussinesq approximation deviation with assumed homogeneity in y-direction (2D) or yz-plane (3D)
         * to a file.
         */
        void
        outputBoussinesqDeviationWithInhomogeneousXDirection(
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
         * Output enstrophy integrated to a file.
         */
        void
        outputEnstrophyIntegrated(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output scalar dissipation rate of first species integrated to a file.
         */
        void
        outputScalarDissipationRateIntegrated(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output numerical interface thickness to a file.
         */
        void
        outputNumericalInterfaceThickness(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output number of cells to a file.
         */
        void
        outputNumberOfCells(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output weighted number of cells to a file.
         */
        void
        outputWeightedNumberOfCells(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;

        /*
         * Output minimum interface location in x-direction to a file.
         */
        void
        outputInterfaceMinInXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Output maximum interface location in x-direction to a file.
         */
        void
        outputInterfaceMaxInXDirection(
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
            num += Y_product_avg_global[i];
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
 * Output TKE integrated with assumed homogeneity in y-direction to a file.
 */
void
RTIRMIStatisticsUtilities::outputTKEIntegratedWithHomogeneityInYDirection(
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
 * Output TKE integrated with assumed homogeneity in yz-plane to a file.
 */
void
RTIRMIStatisticsUtilities::outputTKEIntegratedWithHomogeneityInYZPlane(
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
            << "There is no 'TKE_INT_HOMO_YZ' for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no 'TKE_INT_HOMO_YZ' for two-dimensional problem."
            << std::endl);
    }
    
    // Three-dimensional case.
    
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
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(2);
    
    std::vector<double> rho_w_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    
    std::vector<double> u_tilde(rho_u_avg_global);
    std::vector<double> v_tilde(rho_v_avg_global);
    std::vector<double> w_tilde(rho_w_avg_global);
    for (int i = 0; i < finest_level_dims[0]; i++)
    {
        u_tilde[i] /= rho_avg_global[i];
        v_tilde[i] /= rho_avg_global[i];
        w_tilde[i] /= rho_avg_global[i];
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
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(zeros);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(2);
    averaged_quantities.push_back(w_tilde);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(2);
    averaged_quantities.push_back(w_tilde);
    
    std::vector<double> two_TKE_z_avg_global = MPI_helper_correlation.getQuantityCorrelationWithInhomogeneousXDirection(
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
        TKE_integrated_global += (half*(two_TKE_x_avg_global[i] + two_TKE_y_avg_global[i] + two_TKE_z_avg_global[i])*
            dx_finest[0]);
    }
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    const double L_y = x_hi[1] - x_lo[1];
    const double L_z = x_hi[2] - x_lo[2];
    TKE_integrated_global *= (L_y*L_z);
    
    if (mpi.getRank() == 0)
    {
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << TKE_integrated_global;
        
        f_out.close();
    }
}


/*
 * Output TKE in x-direction integrated with assumed homogeneity in yz-plane to a file.
 */
void
RTIRMIStatisticsUtilities::
outputTKEInXDirectionIntegratedWithHomogeneityInYZPlane(
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
            << "There is no 'TKE_X_INT_HOMO_YZ' for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no 'TKE_X_INT_HOMO_YZ' for two-dimensional problem."
            << std::endl);
    }
    
    // Three-dimensional case.
    
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
    
    std::vector<double> u_tilde(rho_u_avg_global);
    for (int i = 0; i < finest_level_dims[0]; i++)
    {
        u_tilde[i] /= rho_avg_global[i];
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
    
    double TKE_x_integrated_global = double(0);
    const double half = double(1)/double(2);
    for (int i = 0; i < finest_level_dims[0]; i++)
    {
        TKE_x_integrated_global += (half*two_TKE_x_avg_global[i]*dx_finest[0]);
    }
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    const double L_y = x_hi[1] - x_lo[1];
    const double L_z = x_hi[2] - x_lo[2];
    TKE_x_integrated_global *= (L_y*L_z);
    
    if (mpi.getRank() == 0)
    {
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << TKE_x_integrated_global;
        
        f_out.close();
    }
}


/*
 * Output TKE in y-direction integrated with assumed homogeneity in yz-plane to a file.
 */
void
RTIRMIStatisticsUtilities::
outputTKEInYDirectionIntegratedWithHomogeneityInYZPlane(
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
            << "There is no 'TKE_Y_INT_HOMO_YZ' for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no 'TKE_Y_INT_HOMO_YZ' for two-dimensional problem."
            << std::endl);
    }
    
    // Three-dimensional case.
    
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
    component_indices.push_back(1);
    
    std::vector<double> rho_v_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    
    std::vector<double> v_tilde(rho_v_avg_global);
    for (int i = 0; i < finest_level_dims[0]; i++)
    {
        v_tilde[i] /= rho_avg_global[i];
    }
    
    std::vector<std::vector<double> > averaged_quantities;
    
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
    
    double TKE_y_integrated_global = double(0);
    const double half = double(1)/double(2);
    for (int i = 0; i < finest_level_dims[0]; i++)
    {
        TKE_y_integrated_global += (half*two_TKE_y_avg_global[i]*dx_finest[0]);
    }
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    const double L_y = x_hi[1] - x_lo[1];
    const double L_z = x_hi[2] - x_lo[2];
    TKE_y_integrated_global *= (L_y*L_z);
    
    if (mpi.getRank() == 0)
    {
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << TKE_y_integrated_global;
        
        f_out.close();
    }
}


/*
 * Output TKE in z-direction integrated with assumed homogeneity in yz-plane to a file.
 */
void
RTIRMIStatisticsUtilities::
outputTKEInZDirectionIntegratedWithHomogeneityInYZPlane(
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
            << "There is no 'TKE_Z_INT_HOMO_YZ' for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no 'TKE_Z_INT_HOMO_YZ' for two-dimensional problem."
            << std::endl);
    }
    
    // Three-dimensional case.
    
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
    component_indices.push_back(2);
    
    std::vector<double> rho_w_avg_global = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    
    std::vector<double> w_tilde(rho_w_avg_global);
    for (int i = 0; i < finest_level_dims[0]; i++)
    {
        w_tilde[i] /= rho_avg_global[i];
    }
    
    std::vector<std::vector<double> > averaged_quantities;
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(zeros);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(2);
    averaged_quantities.push_back(w_tilde);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(2);
    averaged_quantities.push_back(w_tilde);
    
    std::vector<double> two_TKE_w_avg_global = MPI_helper_correlation.getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    averaged_quantities.clear();
    
    double TKE_w_integrated_global = double(0);
    const double half = double(1)/double(2);
    for (int i = 0; i < finest_level_dims[0]; i++)
    {
        TKE_w_integrated_global += (half*two_TKE_w_avg_global[i]*dx_finest[0]);
    }
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    const double L_y = x_hi[1] - x_lo[1];
    const double L_z = x_hi[2] - x_lo[2];
    TKE_w_integrated_global *= (L_y*L_z);
    
    if (mpi.getRank() == 0)
    {
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << TKE_w_integrated_global;
        
        f_out.close();
    }
}


/*
 * Output mean velocity associated with turbulent mass flux component in x-direction with assumed homogeneity
 * in y-direction (2D) or * yz-plane (3D) to a file.
 */
void
RTIRMIStatisticsUtilities::outputTurbMassFluxVelocityXWithInhomogeneousXDirection(
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
    
    FlowModelMPIHelperCorrelation MPI_helper_correlation = FlowModelMPIHelperCorrelation(
        "MPI_helper_average",
        d_dim,
        d_grid_geometry,
        patch_hierarchy,
        flow_model_tmp);
    
    const hier::IntVector& finest_level_dims = MPI_helper_average.getFinestRefinedDomainNumberOfPoints();
    
    std::vector<double> Y_mean = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        "MASS_FRACTIONS",
        0,
        data_context);
    
    std::vector<double> rho_mean = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        data_context);
    
    std::vector<double> u_mean = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        "VELOCITY",
        0,
        data_context);
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    std::vector<std::vector<double> > averaged_quantities;
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(rho_mean);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(u_mean);
    
    std::vector<double> rho_p_u_p = MPI_helper_correlation.getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    averaged_quantities.clear();
    
    /*
     * Compute and output the mean of a_1 inside mixing layer (only done by process 0).
     */
    if (mpi.getRank() == 0)
    {
        double a_sum = double(0);
        int count = 0;
        
        for (int i = 0; i < finest_level_dims[0]; i++)
        {
            const double mixing_metric = double(4)*Y_mean[i]*(double(1) - Y_mean[i]);
            if (mixing_metric > double(9)/double(10))
            {
                a_sum += rho_p_u_p[i]/rho_mean[i];
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
 * Output density-specific volume covariance with assumed homogeneity in y-direction (2D) or yz-plane (3D)
 * to a file.
 */
void
RTIRMIStatisticsUtilities::outputDensitySpecificVolumeCovarianceWithInhomogeneousXDirection(
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
    
    FlowModelMPIHelperCorrelation MPI_helper_correlation = FlowModelMPIHelperCorrelation(
        "MPI_helper_average",
        d_dim,
        d_grid_geometry,
        patch_hierarchy,
        flow_model_tmp);
    
    const hier::IntVector& finest_level_dims = MPI_helper_average.getFinestRefinedDomainNumberOfPoints();
    
    std::vector<double> Y_mean = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        "MASS_FRACTIONS",
        0,
        data_context);
    
    std::vector<double> rho_mean = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        data_context);
    
    std::vector<double> v_mean = MPI_helper_average.getAveragedReciprocalOfQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        data_context);
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    std::vector<bool> use_reciprocal;
    std::vector<std::vector<double> > averaged_quantities;
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_reciprocal.push_back(false);
    averaged_quantities.push_back(rho_mean);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_reciprocal.push_back(true);
    averaged_quantities.push_back(v_mean);
    
    std::vector<double> rho_p_v_p = MPI_helper_correlation.getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        use_reciprocal,
        averaged_quantities,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    use_reciprocal.clear();
    averaged_quantities.clear();
    
    /*
     * Compute and output the mean of b inside mixing layer (only done by process 0).
     */
    if (mpi.getRank() == 0)
    {
        double b_sum = double(0);
        int count = 0;
        
        for (int i = 0; i < finest_level_dims[0]; i++)
        {
            const double mixing_metric = double(4)*Y_mean[i]*(double(1) - Y_mean[i]);
            if (mixing_metric > double(9)/double(10))
            {
                b_sum += (-rho_p_v_p[i]);
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
 * Output Boussinesq approximation deviation with assumed homogeneity in y-direction (2D) or yz-plane (3D)
 * to a file.
 */
void
RTIRMIStatisticsUtilities::outputBoussinesqDeviationWithInhomogeneousXDirection(
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
    
    FlowModelMPIHelperCorrelation MPI_helper_correlation = FlowModelMPIHelperCorrelation(
        "MPI_helper_average",
        d_dim,
        d_grid_geometry,
        patch_hierarchy,
        flow_model_tmp);
    
    const hier::IntVector& finest_level_dims = MPI_helper_average.getFinestRefinedDomainNumberOfPoints();
    
    std::vector<double> Y_mean = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        "MASS_FRACTIONS",
        0,
        data_context);
    
    std::vector<double> rho_mean = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        data_context);
    
    std::vector<double> v_mean = MPI_helper_average.getAveragedReciprocalOfQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        data_context);
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    std::vector<bool> use_reciprocal;
    std::vector<std::vector<double> > averaged_quantities;
    
    // Compute rho_p_v_p.
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_reciprocal.push_back(false);
    averaged_quantities.push_back(rho_mean);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_reciprocal.push_back(true);
    averaged_quantities.push_back(v_mean);
    
    std::vector<double> rho_p_v_p = MPI_helper_correlation.getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        use_reciprocal,
        averaged_quantities,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    use_reciprocal.clear();
    averaged_quantities.clear();
    
    // Compute rho_p_rho_p.
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(rho_mean);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(rho_mean);
    
    std::vector<double> rho_p_rho_p = MPI_helper_correlation.getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    averaged_quantities.clear();
    
    /*
     * Compute and output the mean of Boussinesq approximation deviation inside mixing layer (only done by process 0).
     */
    if (mpi.getRank() == 0)
    {
        double Boussinesq_dev_sum = double(0);
        int count = 0;
        
        for (int i = 0; i < finest_level_dims[0]; i++)
        {
            const double mixing_metric = double(4)*Y_mean[i]*(double(1) - Y_mean[i]);
            if (mixing_metric > double(9)/double(10))
            {
                Boussinesq_dev_sum += (-rho_mean[i]*rho_mean[i]*rho_p_v_p[i]/rho_p_rho_p[i]);
                count++;
            }
        }
        
        const double Boussinesq_dev_mean = Boussinesq_dev_sum/count;
        
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << Boussinesq_dev_mean;
        
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
    
    const hier::IntVector& finest_level_dims = MPI_helper_average.getFinestRefinedDomainNumberOfPoints();
    
    std::vector<double> Y_mean = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        "MASS_FRACTIONS",
        0,
        data_context);
    
    std::vector<double> rho_mean = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        data_context);
    
    // Compute u_tilde.
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    
    std::vector<double> rho_u_mean = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    
    std::vector<double> u_tilde(rho_u_mean);
    
    for (int i = 0; i < finest_level_dims[0]; i++)
    {
        u_tilde[i] /= rho_mean[i];
    }
    
    // Compute R_11.
    
    std::vector<double> zeros(finest_level_dims[0], double(0));
    
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
    
    std::vector<double> R_11 = MPI_helper_correlation.getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    averaged_quantities.clear();
    
    for (int i = 0; i < finest_level_dims[0]; i++)
    {
        R_11[i] /= rho_mean[i];
    }
    
    /*
     * Compute and output the mean of R_11 inside mixing layer (only done by process 0).
     */
    if (mpi.getRank() == 0)
    {
        double R_11_sum = double(0);
        int count = 0;
        
        for (int i = 0; i < finest_level_dims[0]; i++)
        {
            const double mixing_metric = double(4)*Y_mean[i]*(double(1) - Y_mean[i]);
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
    
    const hier::IntVector& finest_level_dims = MPI_helper_average.getFinestRefinedDomainNumberOfPoints();
    
    std::vector<double> Y_mean = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        "MASS_FRACTIONS",
        0,
        data_context);
    
    std::vector<double> rho_mean = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        data_context);
    
    // Compute v_tilde.
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(1);
    
    std::vector<double> rho_v_mean = MPI_helper_average.getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    
    std::vector<double> v_tilde(rho_v_mean);
    
    for (int i = 0; i < finest_level_dims[0]; i++)
    {
        v_tilde[i] /= rho_mean[i];
    }
    
    // Compute R_22.
    
    std::vector<double> zeros(finest_level_dims[0], double(0));
    
    std::vector<std::vector<double> > averaged_quantities;
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(zeros);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(1);
    averaged_quantities.push_back(v_tilde);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(1);
    averaged_quantities.push_back(v_tilde);
    
    std::vector<double> R_22 = MPI_helper_correlation.getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    averaged_quantities.clear();
    
    for (int i = 0; i < finest_level_dims[0]; i++)
    {
        R_22[i] /= rho_mean[i];
    }
    
    /*
     * Compute and output the mean of R_22 inside mixing layer (only done by process 0).
     */
    if (mpi.getRank() == 0)
    {
        double R_22_sum = double(0);
        int count = 0;
        
        for (int i = 0; i < finest_level_dims[0]; i++)
        {
            const double mixing_metric = double(4)*Y_mean[i]*(double(1) - Y_mean[i]);
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
 * Output enstrophy integrated to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnstrophyIntegrated(
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
 * Output scalar dissipation rate of first species integrated to a file.
 */
void
RTIRMIStatisticsUtilities::outputScalarDissipationRateIntegrated(
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
            double chi_integrated_global = double(0);
            for (int i = 0; i < finest_level_dims[0]; i++)
            {
                chi_integrated_global += scalar_dissipation_rate[i]*dx_finest[0];
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
            
            double chi_integrated_global = double(0);
            for (int i = 0; i < finest_level_dims[0]; i++)
            {
                chi_integrated_global += scalar_dissipation_rate_full[i]*dx_finest[0];
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
            
            double chi_integrated_global = double(0);
            for (int i = 0; i < finest_level_dims[0]; i++)
            {
                chi_integrated_global += scalar_dissipation_rate_full[i]*dx_finest[0];
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
 * Output numerical interface thickness to a file.
 */
void
RTIRMIStatisticsUtilities::outputNumericalInterfaceThickness(
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
            << "'NUM_INTEF_THICK' can be computed with two species only."
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
    
    const std::vector<double>& dx_finest = MPI_helper_max_min.getFinestRefinedDomainGridSpacing();
    
    const hier::IntVector& finest_level_dims = MPI_helper_max_min.getFinestRefinedDomainNumberOfPoints();
    
    const std::vector<double> mag_grad_max_vec = MPI_helper_max_min.getMaxMagnitudeGradientWithInhomogeneousXDirection(
        "MASS_FRACTIONS",
        0,
        d_num_ghosts_derivative,
        data_context);
    
    if (mpi.getRank() == 0)
    {
        double mag_grad_max = double(0);
        for (int i = 0; i < finest_level_dims[0]; i++)
        {
            mag_grad_max = fmax(mag_grad_max, mag_grad_max_vec[i]);
        }
        
        double dx_finest_dir = dx_finest[0];
        for (int di = 1; di < d_dim.getValue(); di++)
        {
            dx_finest_dir = fmax(dx_finest_dir, dx_finest[di]);
        }
        
        f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
              << "\t" << double(1)/(dx_finest_dir*mag_grad_max);
        
        f_out.close();
    }
}


/*
 * Output number of cells to a file.
 */
void
RTIRMIStatisticsUtilities::outputNumberOfCells(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
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
    
    MPIHelperGrid MPI_helper_grid = MPIHelperGrid(
        "MPI_helper_grid",
        d_dim,
        d_grid_geometry,
        patch_hierarchy);
    
    const double num_cells_global = MPI_helper_grid.getNumberOfCells();
    
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
RTIRMIStatisticsUtilities::outputWeightedNumberOfCells(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
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
    
    MPIHelperGrid MPI_helper_grid = MPIHelperGrid(
        "MPI_helper_grid",
        d_dim,
        d_grid_geometry,
        patch_hierarchy);
    
    const double weighted_num_cells_global = MPI_helper_grid.getWeightedNumberOfCells();
    
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
 * Output minimum interface location in x-direction to a file.
 */
void
RTIRMIStatisticsUtilities::outputInterfaceMinInXDirection(
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
    
    if (mpi.getRank() == 0)
    {
        const double* x_lo = d_grid_geometry->getXLower();
        const double* x_hi = d_grid_geometry->getXUpper();
        double interface_min = x_hi[0];
        
        for (int i = finest_level_dims[0] - 1; i >= 0;  i--)
        {
            if (Y_avg_global[i] > 0.01 && Y_avg_global[i] < 0.99)
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
 * Output maximum interface location in x-direction to a file.
 */
void
RTIRMIStatisticsUtilities::outputInterfaceMaxInXDirection(
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
    
    if (mpi.getRank() == 0)
    {
        const double* x_lo = d_grid_geometry->getXLower();
        // const double* x_hi = d_grid_geometry->getXUpper();
        double interface_max = x_lo[0];
        
        for (int i = 0; i < finest_level_dims[0]; i++)
        {
            if (Y_avg_global[i] > 0.01 && Y_avg_global[i] < 0.99)
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
            else if (statistical_quantity_key == "TKE_INT_HOMO_Y")
            {
                f_out << "\t" << "TKE_INT_HOMO_Y       ";
            }
            else if (statistical_quantity_key == "TKE_INT_HOMO_YZ")
            {
                f_out << "\t" << "TKE_INT_HOMO_YZ      ";
            }
            else if (statistical_quantity_key == "TKE_X_INT_HOMO_YZ")
            {
                f_out << "\t" << "TKE_X_INT_HOME_YZ    ";
            }
            else if (statistical_quantity_key == "TKE_Y_INT_HOMO_YZ")
            {
                f_out << "\t" << "TKE_Y_INT_HOME_YZ    ";
            }
            else if (statistical_quantity_key == "TKE_Z_INT_HOMO_YZ")
            {
                f_out << "\t" << "TKE_Z_INT_HOME_YZ    ";
            }
            else if (statistical_quantity_key == "a1_MEAN_INHOMO_X")
            {
                f_out << "\t" << "a1_MEAN_INHOMO_X     ";
            }
            else if (statistical_quantity_key == "b_MEAN_INHOMO_X")
            {
                f_out << "\t" << "b_MEAN_INHOMO_X      ";
            }
            else if (statistical_quantity_key == "BOUSS_MEAN_INHOMO_X")
            {
                f_out << "\t" << "BOUSS_MEAN_INHOMO_X  ";
            }
            else if (statistical_quantity_key == "R11_MEAN_INHOMO_X")
            {
                f_out << "\t" << "R11_MEAN_INHOMO_X    ";
            }
            else if (statistical_quantity_key == "R22_MEAN_INHOMO_X")
            {
                f_out << "\t" << "R22_MEAN_INHOMO_X    ";
            }
            else if (statistical_quantity_key == "ENSTROPHY_INT")
            {
                f_out << "\t" << "ENSTROPHY_INT        ";
            }
            else if (statistical_quantity_key == "SCAL_DISS_RAT_INT")
            {
                f_out << "\t" << "SCAL_DISS_RAT_INT    ";
            }
            else if (statistical_quantity_key == "NUM_INTEF_THICK")
            {
                f_out << "\t" << "NUM_INTEF_THICK      ";
            }
            else if (statistical_quantity_key == "NUM_CELLS")
            {
                f_out << "\t" << "NUM_CELLS            ";
            }
            else if (statistical_quantity_key == "WEIGHTED_NUM_CELLS")
            {
                f_out << "\t" << "WEIGHTED_NUM_CELLS   ";
            }
            else if (statistical_quantity_key == "INTERFACE_MIN_X")
            {
                f_out << "\t" << "INTERFACE_MIN_X      ";
            }
            else if (statistical_quantity_key == "INTERFACE_MAX_X")
            {
                f_out << "\t" << "INTERFACE_MAX_X      ";
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
    
    HAMERS_SHARED_PTR<RTIRMIStatisticsUtilities> rti_rmi_statistics_utilities(
        new RTIRMIStatisticsUtilities(
            "RTI RMI statistics utilities",
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
            rti_rmi_statistics_utilities->outputMixingWidthInXDirection(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "MIXEDNESS_X")
        {
            rti_rmi_statistics_utilities->outputMixednessInXDirection(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "TKE_INT_HOMO_Y")
        {
            rti_rmi_statistics_utilities->outputTKEIntegratedWithHomogeneityInYDirection(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "TKE_INT_HOMO_YZ")
        {
            rti_rmi_statistics_utilities->outputTKEIntegratedWithHomogeneityInYZPlane(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "TKE_X_INT_HOMO_YZ")
        {
            rti_rmi_statistics_utilities->outputTKEInXDirectionIntegratedWithHomogeneityInYZPlane(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "TKE_Y_INT_HOMO_YZ")
        {
            rti_rmi_statistics_utilities->outputTKEInYDirectionIntegratedWithHomogeneityInYZPlane(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "TKE_Z_INT_HOMO_YZ")
        {
            rti_rmi_statistics_utilities->outputTKEInZDirectionIntegratedWithHomogeneityInYZPlane(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "a1_MEAN_INHOMO_X")
        {
            rti_rmi_statistics_utilities->outputTurbMassFluxVelocityXWithInhomogeneousXDirection(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "b_MEAN_INHOMO_X")
        {
            rti_rmi_statistics_utilities->outputDensitySpecificVolumeCovarianceWithInhomogeneousXDirection(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "BOUSS_MEAN_INHOMO_X")
        {
            rti_rmi_statistics_utilities->outputBoussinesqDeviationWithInhomogeneousXDirection(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "R11_MEAN_INHOMO_X")
        {
            rti_rmi_statistics_utilities->outputReynoldsNormalStressXWithInhomogeneousXDirection(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "R22_MEAN_INHOMO_X")
        {
            rti_rmi_statistics_utilities->outputReynoldsNormalStressYWithInhomogeneousXDirection(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "ENSTROPHY_INT")
        {
            rti_rmi_statistics_utilities->outputEnstrophyIntegrated(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "SCAL_DISS_RAT_INT")
        {
            rti_rmi_statistics_utilities->outputScalarDissipationRateIntegrated(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "NUM_INTEF_THICK")
        {
            rti_rmi_statistics_utilities->outputNumericalInterfaceThickness(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "NUM_CELLS")
        {
            rti_rmi_statistics_utilities->outputNumberOfCells(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "WEIGHTED_NUM_CELLS")
        {
            rti_rmi_statistics_utilities->outputWeightedNumberOfCells(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "INTERFACE_MIN_X")
        {
            rti_rmi_statistics_utilities->outputInterfaceMinInXDirection(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "INTERFACE_MAX_X")
        {
            rti_rmi_statistics_utilities->outputInterfaceMaxInXDirection(
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
