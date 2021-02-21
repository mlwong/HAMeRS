#include "flow/flow_models/four-eqn_conservative/FlowModelStatisticsUtilitiesFourEqnConservative.hpp"

#include "extn/patch_hierarchies/ExtendedFlattenedHierarchy.hpp"
#include "flow/flow_models/MPI_helpers/FlowModelMPIHelperAverage.hpp"
#include "flow/flow_models/MPI_helpers/FlowModelMPIHelperCorrelation.hpp"
#include "flow/flow_models/MPI_helpers/FlowModelMPIHelperMaxMin.hpp"
#include "util/MPI_helpers/MPIHelperNumberOfCells.hpp"

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
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Output mixedness in x-direction to a file.
         */
        void
        outputMixednessInXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Output TKE integrated with assumed homogeneity in y-direction to a file.
         */
        void
        outputTKEIntegratedWithHomogeneityInYDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Output TKE integrated with assumed homogeneity in yz-plane to a file.
         */
        void
        outputTKEIntegratedWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Output TKE in x-direction integrated with assumed homogeneity in yz-plane to a file.
         */
        void
        outputTKEInXDirectionIntegratedWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Output TKE in y-direction integrated with assumed homogeneity in yz-plane to a file.
         */
        void
        outputTKEInYDirectionIntegratedWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Output TKE in z-direction integrated with assumed homogeneity in yz-plane to a file.
         */
        void
        outputTKEInZDirectionIntegratedWithHomogeneityInYZPlane(
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
         * Output scalar dissipation rate of first species integrated to a file.
         */
        void
        outputScalarDissipationRateIntegrated(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Output numerical interface thickness to a file.
         */
        void
        outputNumericalInterfaceThickness(
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
RTIRMIStatisticsUtilities::outputMixingWidthInXDirection(
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
 * Output mixedness in x-direction to a file.
 */
void
RTIRMIStatisticsUtilities::outputMixednessInXDirection(
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
    
    const std::vector<double>& dx_finest = MPI_helper_average.getFinestRefinedDomainGridSpacing();
    
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
 * Output enstrophy integrated to a file.
 */
void
RTIRMIStatisticsUtilities::outputEnstrophyIntegrated(
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
            << "There is no 'ENSTROPHY_INT' for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        double Omega_integrated_local  = double(0);
        double Omega_integrated_global = double(0);
        
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
                 * Register the patch, density and velocity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("DENSITY", num_ghosts));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", num_ghosts));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->allocateMemoryForDerivedCellData();
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointers to density and velocity data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_density =
                    flow_model_tmp->getCellData("DENSITY");
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_velocity =
                    flow_model_tmp->getCellData("VELOCITY");
                
                double* rho = data_density->getPointer(0);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_density = data_density->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_velocity = data_velocity->getGhostCellWidth();
                
                const int num_ghosts_0_density = num_ghosts_density[0];
                const int num_ghosts_1_density = num_ghosts_density[1];
                const int ghostcell_dim_0_density = ghostcell_dims_density[0];
                
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(num_ghosts_velocity >= hier::IntVector::getOne(d_dim)*d_num_ghosts_derivative);
#endif
                
                double Omega_to_add = double(0);
                
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
                            
                            const int idx_density = (relative_idx_lo_0 + i + num_ghosts_0_density) +
                                (relative_idx_lo_1 + j + num_ghosts_1_density)*ghostcell_dim_0_density;
                            
                            const double omega = dvdx[idx] - dudy[idx];
                            Omega_to_add += rho[idx_density]*omega*omega/((double) n_overlapped);
                        }
                    }
                }
                
                Omega_to_add = Omega_to_add*dx[0]*dx[1];
                Omega_integrated_local += Omega_to_add;
                
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
            &Omega_integrated_local,
            &Omega_integrated_global,
            1,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        /*
         * Output the enstrophy integral (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << Omega_integrated_global;
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        double Omega_integrated_local  = double(0);
        double Omega_integrated_global = double(0);
        
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
                 * Register the patch, density and velocity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("DENSITY", num_ghosts));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", num_ghosts));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->allocateMemoryForDerivedCellData();
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointers to density and velocity data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_density =
                    flow_model_tmp->getCellData("DENSITY");
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_velocity =
                    flow_model_tmp->getCellData("VELOCITY");
                
                double* rho = data_density->getPointer(0);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_density = data_density->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_velocity = data_velocity->getGhostCellWidth();
                
                const int num_ghosts_0_density = num_ghosts_density[0];
                const int num_ghosts_1_density = num_ghosts_density[1];
                const int num_ghosts_2_density = num_ghosts_density[2];
                const int ghostcell_dim_0_density = ghostcell_dims_density[0];
                const int ghostcell_dim_1_density = ghostcell_dims_density[1];
                
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(num_ghosts_velocity >= hier::IntVector::getOne(d_dim)*d_num_ghosts_derivative);
#endif
                
                double Omega_to_add = double(0);
                
                /*
                 * Initialize cell data for velocity derivatives and get pointers to the derivatives.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_velocity_derivatives(
                    new pdat::CellData<double>(patch_box, d_dim.getValue()*d_dim.getValue() - d_dim.getValue(),
                        hier::IntVector::getZero(d_dim)));
                
                double* dudy = data_velocity_derivatives->getPointer(0);
                double* dudz = data_velocity_derivatives->getPointer(1);
                double* dvdx = data_velocity_derivatives->getPointer(2);
                double* dvdz = data_velocity_derivatives->getPointer(3);
                double* dwdx = data_velocity_derivatives->getPointer(4);
                double* dwdy = data_velocity_derivatives->getPointer(5);
                
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
                                
                                const int idx_density = (relative_idx_lo_0 + i + num_ghosts_0_density) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_density)*ghostcell_dim_0_density +
                                    (relative_idx_lo_2 + k + num_ghosts_2_density)*ghostcell_dim_0_density*
                                        ghostcell_dim_1_density;
                                
                                const double omega_x = dwdy[idx] - dvdz[idx];
                                const double omega_y = dudz[idx] - dwdx[idx];
                                const double omega_z = dvdx[idx] - dudy[idx];
                                
                                Omega_to_add += rho[idx_density]*(
                                    omega_x*omega_x + omega_y*omega_y + omega_z*omega_z)/((double) n_overlapped);
                            }
                        }
                    }
                }
                
                Omega_to_add = Omega_to_add*dx[0]*dx[1]*dx[2];
                Omega_integrated_local += Omega_to_add;
                
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
            &Omega_integrated_local,
            &Omega_integrated_global,
            1,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        /*
         * Output the enstrophy integral (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
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
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (!d_equation_of_mass_diffusivity_mixing_rules)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Mixing rule of mass diffusivity is not initialized yet."
            << std::endl);
    }
    
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
        double Chi_integrated_local  = double(0);
        double Chi_integrated_global = double(0);
        
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
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch, mass fractions, pressure and temperature in the flow model
                 * and compute the corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTIONS", num_ghosts));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("PRESSURE", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("TEMPERATURE", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->allocateMemoryForDerivedCellData();
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointers to mass fraction, pressure and temperature data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_mass_fractions =
                    flow_model_tmp->getCellData("MASS_FRACTIONS");
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_pressure =
                    flow_model_tmp->getCellData("PRESSURE");
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_temperature =
                    flow_model_tmp->getCellData("TEMPERATURE");
                
                std::vector<double*> Y;
                Y.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    Y.push_back(data_mass_fractions->getPointer(si));
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
                
                const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
                const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
                const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
                
                const int num_ghosts_0_mass_fractions = num_ghosts_mass_fractions[0];
                const int num_ghosts_0_pressure = num_ghosts_pressure[0];
                const int num_ghosts_0_temperature = num_ghosts_temperature[0];
                
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(num_ghosts_mass_fractions >= hier::IntVector::getOne(d_dim)*d_num_ghosts_derivative);
#endif
                
                double Chi_to_add = double(0);
                
                /*
                 * Initialize cell data for mass fraction derivatives and get pointers to the derivatives.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_mass_fraction_derivatives(
                    new pdat::CellData<double>(patch_box, d_dim.getValue(), hier::IntVector::getZero(d_dim)));
                
                double* dYdx = data_mass_fraction_derivatives->getPointer(0);
                
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
                    
                    HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_x(
                        new DerivativeFirstOrder(
                            "first order derivative in x-direction",
                            d_dim,
                            DIRECTION::X_DIRECTION,
                            d_num_ghosts_derivative));
                    
                    // Compute dYdx.
                    derivative_first_order_x->computeDerivative(
                        data_mass_fraction_derivatives,
                        data_mass_fractions,
                        dx[0],
                        patch_visible_box,
                        0,
                        0);
                    
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
                        
                        // Compute linear indices of derivatives, mass fractions, pressure and temperature.
                        const int idx = relative_idx_lo_0 + i;
                        const int idx_mass_fractions = relative_idx_lo_0 + i + num_ghosts_0_mass_fractions;
                        const int idx_pressure = relative_idx_lo_0 + i + num_ghosts_0_pressure;
                        const int idx_temperature = relative_idx_lo_0 + i + num_ghosts_0_temperature;
                        
                        std::vector<double> D;
                        D.resize(d_num_species);
                        
                        std::vector<double*> D_ptr;
                        D_ptr.reserve(d_num_species);
                        for (int si = 0; si < d_num_species; si++)
                        {
                            D_ptr.push_back(&D[si]);
                        }
                        
                        std::vector<const double*> Y_ptr;
                        Y_ptr.reserve(d_num_species);
                        for (int si = 0; si < d_num_species; si++)
                        {
                            Y_ptr.push_back(&Y[si][idx_mass_fractions]);
                        }
                        
                        d_equation_of_mass_diffusivity_mixing_rules->
                            getMassDiffusivities(
                                D_ptr,
                                &p[idx_pressure],
                                &T[idx_temperature],
                                Y_ptr);
                        
                        Chi_to_add += D[0]*dYdx[idx]*dYdx[idx]/((double) n_overlapped);
                    }
                }
                
                Chi_to_add = Chi_to_add*dx[0];
                Chi_integrated_local += Chi_to_add;
                
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
            &Chi_integrated_local,
            &Chi_integrated_global,
            1,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        /*
         * Output the scalar dissipation integral (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << Chi_integrated_global;
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        double Chi_integrated_local  = double(0);
        double Chi_integrated_global = double(0);
        
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
                 * Register the patch, mass fractions, pressure and temperature in the flow model
                 * and compute the corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTIONS", num_ghosts));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("PRESSURE", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("TEMPERATURE", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->allocateMemoryForDerivedCellData();
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointers to mass fraction, pressure and temperature data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_mass_fractions =
                    flow_model_tmp->getCellData("MASS_FRACTIONS");
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_pressure =
                    flow_model_tmp->getCellData("PRESSURE");
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_temperature =
                    flow_model_tmp->getCellData("TEMPERATURE");
                
                std::vector<double*> Y;
                Y.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    Y.push_back(data_mass_fractions->getPointer(si));
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
                
                const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mass_fractions = data_mass_fractions->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_pressure = data_pressure->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_temperature = data_temperature->getGhostBox().numberCells();
                
                const int num_ghosts_0_mass_fractions = num_ghosts_mass_fractions[0];
                const int num_ghosts_1_mass_fractions = num_ghosts_mass_fractions[1];
                const int ghostcell_dim_0_mass_fractions = ghostcell_dims_mass_fractions[0];
                
                const int num_ghosts_0_pressure = num_ghosts_pressure[0];
                const int num_ghosts_1_pressure = num_ghosts_pressure[1];
                const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
                
                const int num_ghosts_0_temperature = num_ghosts_temperature[0];
                const int num_ghosts_1_temperature = num_ghosts_temperature[1];
                const int ghostcell_dim_0_temperature = ghostcell_dims_temperature[0];
                
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(num_ghosts_mass_fractions >= hier::IntVector::getOne(d_dim)*d_num_ghosts_derivative);
#endif
                
                double Chi_to_add = double(0);
                
                /*
                 * Initialize cell data for mass fraction derivatives and get pointers to the derivatives.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_mass_fraction_derivatives(
                    new pdat::CellData<double>(patch_box, d_dim.getValue(), hier::IntVector::getZero(d_dim)));
                
                double* dYdx = data_mass_fraction_derivatives->getPointer(0);
                double* dYdy = data_mass_fraction_derivatives->getPointer(1);
                
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
                    
                    // Compute dYdx.
                    derivative_first_order_x->computeDerivative(
                        data_mass_fraction_derivatives,
                        data_mass_fractions,
                        dx[0],
                        patch_visible_box,
                        0,
                        0);
                    
                    // Compute dYdy.
                    derivative_first_order_y->computeDerivative(
                        data_mass_fraction_derivatives,
                        data_mass_fractions,
                        dx[1],
                        patch_visible_box,
                        1,
                        0);
                    
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
                            
                            // Compute linear indices of derivatives, mass fractions, pressure and temperature.
                            const int idx = (relative_idx_lo_0 + i) +
                                (relative_idx_lo_1 + j)*patch_dim_0;
                            
                            const int idx_mass_fractions = (relative_idx_lo_0 + i + num_ghosts_0_mass_fractions) +
                                (relative_idx_lo_1 + j + num_ghosts_1_mass_fractions)*ghostcell_dim_0_mass_fractions;
                            
                            const int idx_pressure = (relative_idx_lo_0 + i + num_ghosts_0_pressure) +
                                (relative_idx_lo_1 + j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure;
                            
                            const int idx_temperature = (relative_idx_lo_0 + i + num_ghosts_0_temperature) +
                                (relative_idx_lo_1 + j + num_ghosts_1_temperature)*ghostcell_dim_0_temperature;
                            
                            std::vector<double> D;
                            D.resize(d_num_species);
                            
                            std::vector<double*> D_ptr;
                            D_ptr.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                D_ptr.push_back(&D[si]);
                            }
                            
                            std::vector<const double*> Y_ptr;
                            Y_ptr.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_ptr.push_back(&Y[si][idx_mass_fractions]);
                            }
                            
                            d_equation_of_mass_diffusivity_mixing_rules->
                                getMassDiffusivities(
                                    D_ptr,
                                    &p[idx_pressure],
                                    &T[idx_temperature],
                                    Y_ptr);
                            
                            Chi_to_add += D[0]*(dYdx[idx]*dYdx[idx] + dYdy[idx]*dYdy[idx])/((double) n_overlapped);
                        }
                    }
                }
                
                Chi_to_add = Chi_to_add*dx[0]*dx[1];
                Chi_integrated_local += Chi_to_add;
                
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
            &Chi_integrated_local,
            &Chi_integrated_global,
            1,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        /*
         * Output the scalar dissipation integral (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << Chi_integrated_global;
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        double Chi_integrated_local  = double(0);
        double Chi_integrated_global = double(0);
        
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
                 * Register the patch, mass fractions, pressure and temperature in the flow model
                 * and compute the corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTIONS", num_ghosts));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("PRESSURE", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("TEMPERATURE", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->allocateMemoryForDerivedCellData();
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointers to mass fraction, pressure and temperature data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_mass_fractions =
                    flow_model_tmp->getCellData("MASS_FRACTIONS");
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_pressure =
                    flow_model_tmp->getCellData("PRESSURE");
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_temperature =
                    flow_model_tmp->getCellData("TEMPERATURE");
                
                std::vector<double*> Y;
                Y.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    Y.push_back(data_mass_fractions->getPointer(si));
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
                
                const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mass_fractions = data_mass_fractions->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_pressure = data_pressure->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_temperature = data_temperature->getGhostBox().numberCells();
                
                const int num_ghosts_0_mass_fractions = num_ghosts_mass_fractions[0];
                const int num_ghosts_1_mass_fractions = num_ghosts_mass_fractions[1];
                const int num_ghosts_2_mass_fractions = num_ghosts_mass_fractions[2];
                const int ghostcell_dim_0_mass_fractions = ghostcell_dims_mass_fractions[0];
                const int ghostcell_dim_1_mass_fractions = ghostcell_dims_mass_fractions[1];
                
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
                
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(num_ghosts_mass_fractions >= hier::IntVector::getOne(d_dim)*d_num_ghosts_derivative);
#endif
                
                double Chi_to_add = double(0);
                
                /*
                 * Initialize cell data for mass fraction derivatives and get pointers to the derivatives.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_mass_fraction_derivatives(
                    new pdat::CellData<double>(patch_box, d_dim.getValue(), hier::IntVector::getZero(d_dim)));
                
                double* dYdx = data_mass_fraction_derivatives->getPointer(0);
                double* dYdy = data_mass_fraction_derivatives->getPointer(1);
                double* dYdz = data_mass_fraction_derivatives->getPointer(2);
                
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
                    
                    // Compute dYdx.
                    derivative_first_order_x->computeDerivative(
                        data_mass_fraction_derivatives,
                        data_mass_fractions,
                        dx[0],
                        patch_visible_box,
                        0,
                        0);
                    
                    // Compute dYdy.
                    derivative_first_order_y->computeDerivative(
                        data_mass_fraction_derivatives,
                        data_mass_fractions,
                        dx[1],
                        patch_visible_box,
                        1,
                        0);
                    
                    // Compute dYdz.
                    derivative_first_order_z->computeDerivative(
                        data_mass_fraction_derivatives,
                        data_mass_fractions,
                        dx[2],
                        patch_visible_box,
                        2,
                        0);
                    
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
                                
                                // Compute linear indices of derivatives, mass fractions, pressure and temperature.
                                const int idx = (relative_idx_lo_0 + i) +
                                    (relative_idx_lo_1 + j)*patch_dim_0 +
                                    (relative_idx_lo_2 + k)*patch_dim_0*
                                        patch_dim_1;
                                
                                const int idx_mass_fractions = (relative_idx_lo_0 + i + num_ghosts_0_mass_fractions) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fractions)*ghostcell_dim_0_mass_fractions +
                                    (relative_idx_lo_2 + k + num_ghosts_2_mass_fractions)*ghostcell_dim_0_mass_fractions*
                                        ghostcell_dim_1_mass_fractions;
                                
                                const int idx_pressure = (relative_idx_lo_0 + i + num_ghosts_0_pressure) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure +
                                    (relative_idx_lo_2 + k + num_ghosts_2_pressure)*ghostcell_dim_0_pressure*
                                        ghostcell_dim_1_pressure;
                                
                                const int idx_temperature = (relative_idx_lo_0 + i + num_ghosts_0_temperature) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_temperature)*ghostcell_dim_0_temperature +
                                    (relative_idx_lo_2 + k + num_ghosts_2_temperature)*ghostcell_dim_0_temperature*
                                        ghostcell_dim_1_temperature;
                                
                                std::vector<double> D;
                                D.resize(d_num_species);
                                
                                std::vector<double*> D_ptr;
                                D_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    D_ptr.push_back(&D[si]);
                                }
                                
                                std::vector<const double*> Y_ptr;
                                Y_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_ptr.push_back(&Y[si][idx_mass_fractions]);
                                }
                                
                                d_equation_of_mass_diffusivity_mixing_rules->
                                    getMassDiffusivities(
                                        D_ptr,
                                        &p[idx_pressure],
                                        &T[idx_temperature],
                                        Y_ptr);
                                
                                Chi_to_add += D[0]*(dYdx[idx]*dYdx[idx] + dYdy[idx]*dYdy[idx] + dYdz[idx]*dYdz[idx])/
                                    ((double) n_overlapped);
                            }
                        }
                    }
                }
                
                Chi_to_add = Chi_to_add*dx[0]*dx[1]*dx[2];
                Chi_integrated_local += Chi_to_add;
                
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
            &Chi_integrated_local,
            &Chi_integrated_global,
            1,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        /*
         * Output the scalar dissipation integral (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
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
 * Output numerical interface thickness to a file.
 */
void
RTIRMIStatisticsUtilities::outputNumericalInterfaceThickness(
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
    
    const int num_ghosts_derivative = 3; // assume we have at least three ghost cells
    
    const std::vector<double> mag_grad_max_vec = MPI_helper_max_min.getMaxMagnitudeGradientWithInhomogeneousXDirection(
        "MASS_FRACTIONS",
        0,
        num_ghosts_derivative,
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
    
    MPIHelperNumberOfCells MPI_helper_num_cells = MPIHelperNumberOfCells(
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
RTIRMIStatisticsUtilities::outputWeightedNumberOfCells(
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
    
    MPIHelperNumberOfCells MPI_helper_num_cells = MPIHelperNumberOfCells(
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
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    HAMERS_SHARED_PTR<RTIRMIStatisticsUtilities> rti_rmi_statistics_utilities(
        new RTIRMIStatisticsUtilities(
            "RMI statistics utilities",
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
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Unknown statistical quantity key = '"
                << statistical_quantity_key
                << " found."
                << std::endl);
        }
    }
}
