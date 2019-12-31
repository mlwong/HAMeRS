#include "flow/flow_models/four-eqn_conservative/FlowModelStatisticsUtilitiesFourEqnConservative.hpp"

#include "extn/patch_hierarchies/ExtendedFlattenedHierarchy.hpp"

#include <fstream>



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
            // Get the key of the current variable.
            std::string statistical_quantity_key = d_statistical_quantities[qi];
            
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
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double output_time)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    // Loop over statistical quantities.
    for (int qi = 0; qi < static_cast<int>(d_statistical_quantities.size()); qi++)
    {
        // Get the key of the current variable.
        std::string statistical_quantity_key = d_statistical_quantities[qi];
        
        if (statistical_quantity_key == "DENSITY")
        {
            outputAveragedDensityWithInhomogeneousXDirection(
                "rho_mean.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "MASS_FRACTION")
        {
            outputAveragedMassFractionWithInhomogeneousXDirection(
                "Y_mean.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "MOLE_FRACTION")
        {
            outputAveragedMoleFractionWithInhomogeneousXDirection(
                "X_mean.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "SPECIFIC_VOLUME")
        {
            outputAveragedSpecificVolumeWithInhomogeneousXDirection(
                "rho_inv_mean.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "DENSITY_VARIANCE")
        {
            outputDensityVarianceWithInhomogeneousXDirection(
                "rho_variance.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "TKE")
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Not implemented!"
                << std::endl);
        }
        else if (statistical_quantity_key == "rR11")
        {
            outputReynoldsNormalStressInXDirectionWithInhomogeneousXDirection(
                "rR11.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "rR22")
        {
            outputReynoldsNormalStressInYDirectionWithInhomogeneousXDirection(
                "rR22.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "rR33")
        {
            outputReynoldsNormalStressInZDirectionWithInhomogeneousXDirection(
                "rR33.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "rR12")
        {
            outputReynoldsShearStressInXYDirectionWithInhomogeneousXDirection(
                "rR12.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "rR13")
        {
            outputReynoldsShearStressInXZDirectionWithInhomogeneousXDirection(
                "rR13.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "rR23")
        {
            outputReynoldsShearStressInYZDirectionWithInhomogeneousXDirection(
                "rR23.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "ra1")
        {
            outputAveragedTurbMassFluxXWithInhomogeneousXDirection(
                "ra1.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "ra2")
        {
            outputAveragedTurbMassFluxYWithInhomogeneousXDirection(
                "ra2.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "ra3")
        {
            outputAveragedTurbMassFluxZWithInhomogeneousXDirection(
                "ra3.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "u_p_u_p")
        {
            outputVelocityComponentInXDirectionVarianceWithInhomogeneousXDirection(
                "u_p_u_p.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "v_p_v_p")
        {
            outputVelocityComponentInYDirectionVarianceWithInhomogeneousXDirection(
                "v_p_v_p.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "w_p_w_p")
        {
            outputVelocityComponentInZDirectionVarianceWithInhomogeneousXDirection(
                "w_p_w_p.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "rho_p_u_p_u_p")
        {
            outputDensityVelocityComponentSquareInXDirectionCorrelationWithInhomogeneousXDirection(
                "rho_p_u_p_u_p.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "rho_p_v_p_v_p")
        {
            outputDensityVelocityComponentSquareInYDirectionCorrelationWithInhomogeneousXDirection(
                "rho_p_v_p_v_p.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "rho_p_w_p_w_p")
        {
            outputDensityVelocityComponentSquareInZDirectionCorrelationWithInhomogeneousXDirection(
                "rho_p_w_p_w_p.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "b")
        {
            outputDensitySpecificVolumeCovarianceWithInhomogeneousXDirection(
                "b.dat",
                patch_hierarchy,
                data_context,
                output_time);
        }
        else if (statistical_quantity_key == "ra1_budget")
        {
            outputBudgetTurbMassFluxXWithInhomogeneousXDirection(
                "ra1_budget.dat",
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
                << " found."
                << std::endl);
        }
    }
}


/*
 * Output averaged density with inhomogeneous x-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputAveragedDensityWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double output_time) const
{
    std::vector<double> rho_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        patch_hierarchy,
        data_context);
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_output;
        f_output.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        
        f_output.write((char*)&output_time, sizeof(double));
        f_output.write((char*)&rho_mean[0], sizeof(double)*rho_mean.size());
        
        f_output.close();
    }
}


/*
 * Output averaged mass fraction (first species) with inhomogeneous x-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputAveragedMassFractionWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double output_time) const
{
    std::vector<double> Y_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "MASS_FRACTION",
        0,
        patch_hierarchy,
        data_context);
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_output;
        f_output.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        
        f_output.write((char*)&output_time, sizeof(double));
        f_output.write((char*)&Y_mean[0], sizeof(double)*Y_mean.size());
        
        f_output.close();
    }
}


/*
 * Output averaged mole fraction (first species) with inhomogeneous x-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputAveragedMoleFractionWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double output_time) const
{
    std::vector<double> X_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "MOLE_FRACTION",
        0,
        patch_hierarchy,
        data_context);
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_output;
        f_output.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        
        f_output.write((char*)&output_time, sizeof(double));
        f_output.write((char*)&X_mean[0], sizeof(double)*X_mean.size());
        
        f_output.close();
    }
}


/*
 * Output averaged specific volume with inhomogeneous x-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputAveragedSpecificVolumeWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double output_time) const
{
    std::vector<double> v_mean = getAveragedReciprocalOfQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        patch_hierarchy,
        data_context);
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_output;
        f_output.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        
        f_output.write((char*)&output_time, sizeof(double));
        f_output.write((char*)&v_mean[0], sizeof(double)*v_mean.size());
        
        f_output.close();
    }
}



/*
 * Output turbulent mass flux in x-direction with inhomogeneous x-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputAveragedTurbMassFluxXWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double output_time) const
{
    std::vector<double> rho_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        patch_hierarchy,
        data_context);
    
    std::vector<double> u_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "VELOCITY",
        0,
        patch_hierarchy,
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
    
    std::vector<double> rho_p_u_p = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_output;
        f_output.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        
        f_output.write((char*)&output_time, sizeof(double));
        f_output.write((char*)&rho_p_u_p[0], sizeof(double)*rho_p_u_p.size());
        
        f_output.close();
    }
}


/*
 * Output turbulent mass flux in y-direction with inhomogeneous x-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputAveragedTurbMassFluxYWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double output_time) const
{
    std::vector<double> rho_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        patch_hierarchy,
        data_context);
    
    std::vector<double> v_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "VELOCITY",
        1,
        patch_hierarchy,
        data_context);
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    std::vector<std::vector<double> > averaged_quantities;
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(rho_mean);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(1);
    averaged_quantities.push_back(v_mean);
    
    std::vector<double> rho_p_v_p = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_output;
        f_output.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        
        f_output.write((char*)&output_time, sizeof(double));
        f_output.write((char*)&rho_p_v_p[0], sizeof(double)*rho_p_v_p.size());
        
        f_output.close();
    }
}


/*
 * Output turbulent mass flux in z-direction with inhomogeneous x-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputAveragedTurbMassFluxZWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double output_time) const
{
    std::vector<double> rho_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        patch_hierarchy,
        data_context);
    
    std::vector<double> w_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "VELOCITY",
        2,
        patch_hierarchy,
        data_context);
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    std::vector<std::vector<double> > averaged_quantities;
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(rho_mean);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(2);
    averaged_quantities.push_back(w_mean);
    
    std::vector<double> rho_p_w_p = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_output;
        f_output.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        
        f_output.write((char*)&output_time, sizeof(double));
        f_output.write((char*)&rho_p_w_p[0], sizeof(double)*rho_p_w_p.size());
        
        f_output.close();
    }
}


/*
 * Output Reynolds normal stress in x-direction with inhomogeneous x-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputReynoldsNormalStressInXDirectionWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double output_time) const
{
    const int finest_level_dim_0 = getRefinedDomainNumberOfPointsX(patch_hierarchy);
    
    std::vector<double> rho_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        patch_hierarchy,
        data_context);
    
    // Compute u_tilde.
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    
    std::vector<double> rho_u_mean = getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        patch_hierarchy,
        data_context);
    
    std::vector<double> u_tilde(rho_u_mean);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        u_tilde[i] /= rho_mean[i];
    }
    
    // Compute R_11.
    
    std::vector<double> zeros(finest_level_dim_0, double(0));
    
    quantity_names.clear();
    component_indices.clear();
    
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
    
    std::vector<double> R_11 = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    // Output R_11.
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_output;
        f_output.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        
        f_output.write((char*)&output_time, sizeof(double));
        f_output.write((char*)&R_11[0], sizeof(double)*R_11.size());
        
        f_output.close();
    }
}


/*
 * Output Reynolds normal stress in y-direction with inhomogeneous x-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputReynoldsNormalStressInYDirectionWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double output_time) const
{
    const int finest_level_dim_0 = getRefinedDomainNumberOfPointsX(patch_hierarchy);
    
    std::vector<double> rho_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        patch_hierarchy,
        data_context);
    
    // Compute v_tilde.
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(1);
    
    std::vector<double> rho_v_mean = getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        patch_hierarchy,
        data_context);
    
    std::vector<double> v_tilde(rho_v_mean);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        v_tilde[i] /= rho_mean[i];
    }
    
    // Compute R_22.
    
    std::vector<double> zeros(finest_level_dim_0, double(0));
    
    quantity_names.clear();
    component_indices.clear();
    
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
    
    std::vector<double> R_22 = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    // Output R_22.
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_output;
        f_output.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        
        f_output.write((char*)&output_time, sizeof(double));
        f_output.write((char*)&R_22[0], sizeof(double)*R_22.size());
        
        f_output.close();
    }
}


/*
 * Output Reynolds normal stress in z-direction with inhomogeneous x-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputReynoldsNormalStressInZDirectionWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double output_time) const
{
    const int finest_level_dim_0 = getRefinedDomainNumberOfPointsX(patch_hierarchy);
    
    std::vector<double> rho_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        patch_hierarchy,
        data_context);
    
    // Compute w_tilde.
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(2);
    
    std::vector<double> rho_w_mean = getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        patch_hierarchy,
        data_context);
    
    std::vector<double> w_tilde(rho_w_mean);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        w_tilde[i] /= rho_mean[i];
    }
    
    // Compute R_33.
    
    std::vector<double> zeros(finest_level_dim_0, double(0));
    
    quantity_names.clear();
    component_indices.clear();
    
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
    
    std::vector<double> R_33 = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    // Output R_33.
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_output;
        f_output.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        
        f_output.write((char*)&output_time, sizeof(double));
        f_output.write((char*)&R_33[0], sizeof(double)*R_33.size());
        
        f_output.close();
    }
}


/*
 * Output Reynolds shear stress in xy-direction with inhomogeneous x-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputReynoldsShearStressInXYDirectionWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double output_time) const
{
    const int finest_level_dim_0 = getRefinedDomainNumberOfPointsX(patch_hierarchy);
    
    std::vector<double> rho_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        patch_hierarchy,
        data_context);
    
    // Compute u_tilde.
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    
    std::vector<double> rho_u_mean = getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        patch_hierarchy,
        data_context);
    
    std::vector<double> u_tilde(rho_u_mean);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        u_tilde[i] /= rho_mean[i];
    }
    
    // Compute v_tilde.
    
    quantity_names.clear();
    component_indices.clear();
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(1);
    
    std::vector<double> rho_v_mean = getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        patch_hierarchy,
        data_context);
    
    std::vector<double> v_tilde(rho_v_mean);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        v_tilde[i] /= rho_mean[i];
    }
    
    // Compute R_12.
    
    std::vector<double> zeros(finest_level_dim_0, double(0));
    
    quantity_names.clear();
    component_indices.clear();
    
    std::vector<std::vector<double> > averaged_quantities;
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(zeros);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(u_tilde);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(1);
    averaged_quantities.push_back(v_tilde);
    
    std::vector<double> R_12 = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    // Output R_12.
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_output;
        f_output.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        
        f_output.write((char*)&output_time, sizeof(double));
        f_output.write((char*)&R_12[0], sizeof(double)*R_12.size());
        
        f_output.close();
    }
}


/*
 * Output Reynolds shear stress in xz-direction with inhomogeneous x-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputReynoldsShearStressInXZDirectionWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double output_time) const
{
    const int finest_level_dim_0 = getRefinedDomainNumberOfPointsX(patch_hierarchy);
    
    std::vector<double> rho_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        patch_hierarchy,
        data_context);
    
    // Compute u_tilde.
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    
    std::vector<double> rho_u_mean = getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        patch_hierarchy,
        data_context);
    
    std::vector<double> u_tilde(rho_u_mean);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        u_tilde[i] /= rho_mean[i];
    }
    
    // Compute w_tilde.
    
    quantity_names.clear();
    component_indices.clear();
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(2);
    
    std::vector<double> rho_w_mean = getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        patch_hierarchy,
        data_context);
    
    std::vector<double> w_tilde(rho_w_mean);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        w_tilde[i] /= rho_mean[i];
    }
    
    // Compute R_13.
    
    std::vector<double> zeros(finest_level_dim_0, double(0));
    
    quantity_names.clear();
    component_indices.clear();
    
    std::vector<std::vector<double> > averaged_quantities;
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(zeros);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(u_tilde);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(2);
    averaged_quantities.push_back(w_tilde);
    
    std::vector<double> R_13 = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    // Output R_13.
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_output;
        f_output.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        
        f_output.write((char*)&output_time, sizeof(double));
        f_output.write((char*)&R_13[0], sizeof(double)*R_13.size());
        
        f_output.close();
    }
}


/*
 * Output Reynolds shear stress in yz-direction with inhomogeneous x-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputReynoldsShearStressInYZDirectionWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double output_time) const
{
    const int finest_level_dim_0 = getRefinedDomainNumberOfPointsX(patch_hierarchy);
    
    std::vector<double> rho_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        patch_hierarchy,
        data_context);
    
    // Compute v_tilde.
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(1);
    
    std::vector<double> rho_v_mean = getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        patch_hierarchy,
        data_context);
    
    std::vector<double> v_tilde(rho_v_mean);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        v_tilde[i] /= rho_mean[i];
    }
    
    // Compute w_tilde.
    
    quantity_names.clear();
    component_indices.clear();
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(2);
    
    std::vector<double> rho_w_mean = getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        patch_hierarchy,
        data_context);
    
    std::vector<double> w_tilde(rho_w_mean);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        w_tilde[i] /= rho_mean[i];
    }
    
    // Compute R_23.
    
    std::vector<double> zeros(finest_level_dim_0, double(0));
    
    quantity_names.clear();
    component_indices.clear();
    
    std::vector<std::vector<double> > averaged_quantities;
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(zeros);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(1);
    averaged_quantities.push_back(v_tilde);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(2);
    averaged_quantities.push_back(w_tilde);
    
    std::vector<double> R_23 = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    // Output R_23.
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_output;
        f_output.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        
        f_output.write((char*)&output_time, sizeof(double));
        f_output.write((char*)&R_23[0], sizeof(double)*R_23.size());
        
        f_output.close();
    }
}


/*
 * Output variance of velocity component in x-direction with inhomogeneous x-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputVelocityComponentInXDirectionVarianceWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double output_time) const
{
    std::vector<double> u_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "VELOCITY",
        0,
        patch_hierarchy,
        data_context);
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    std::vector<std::vector<double> > averaged_quantities;
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(u_mean);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(u_mean);
    
    std::vector<double> u_p_u_p = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_output;
        f_output.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        
        f_output.write((char*)&output_time, sizeof(double));
        f_output.write((char*)&u_p_u_p[0], sizeof(double)*u_p_u_p.size());
        
        f_output.close();
    }
}


/*
 * Output variance of velocity component in y-direction with inhomogeneous x-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputVelocityComponentInYDirectionVarianceWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double output_time) const
{
    std::vector<double> v_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "VELOCITY",
        1,
        patch_hierarchy,
        data_context);
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    std::vector<std::vector<double> > averaged_quantities;
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(1);
    averaged_quantities.push_back(v_mean);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(1);
    averaged_quantities.push_back(v_mean);
    
    std::vector<double> v_p_v_p = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_output;
        f_output.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        
        f_output.write((char*)&output_time, sizeof(double));
        f_output.write((char*)&v_p_v_p[0], sizeof(double)*v_p_v_p.size());
        
        f_output.close();
    }
}


/*
 * Output variance of velocity component in z-direction with inhomogeneous x-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputVelocityComponentInZDirectionVarianceWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double output_time) const
{
    std::vector<double> w_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "VELOCITY",
        2,
        patch_hierarchy,
        data_context);
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    std::vector<std::vector<double> > averaged_quantities;
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(2);
    averaged_quantities.push_back(w_mean);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(2);
    averaged_quantities.push_back(w_mean);
    
    std::vector<double> w_p_w_p = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_output;
        f_output.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        
        f_output.write((char*)&output_time, sizeof(double));
        f_output.write((char*)&w_p_w_p[0], sizeof(double)*w_p_w_p.size());
        
        f_output.close();
    }
}


/*
 * Output correlation of density and velocity component square in x-direction with inhomogeneous x-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::
outputDensityVelocityComponentSquareInXDirectionCorrelationWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double output_time) const
{
    std::vector<double> rho_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        patch_hierarchy,
        data_context);
    
    std::vector<double> u_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "VELOCITY",
        0,
        patch_hierarchy,
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
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(u_mean);
    
    std::vector<double> rho_p_u_p_u_p = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_output;
        f_output.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        
        f_output.write((char*)&output_time, sizeof(double));
        f_output.write((char*)&rho_p_u_p_u_p[0], sizeof(double)*rho_p_u_p_u_p.size());
        
        f_output.close();
    }
}


/*
 * Output correlation of density and velocity component square in y-direction with inhomogeneous x-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::
outputDensityVelocityComponentSquareInYDirectionCorrelationWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double output_time) const
{
    std::vector<double> rho_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        patch_hierarchy,
        data_context);
    
    std::vector<double> v_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "VELOCITY",
        1,
        patch_hierarchy,
        data_context);
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    std::vector<std::vector<double> > averaged_quantities;
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(rho_mean);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(1);
    averaged_quantities.push_back(v_mean);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(1);
    averaged_quantities.push_back(v_mean);
    
    std::vector<double> rho_p_v_p_v_p = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_output;
        f_output.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        
        f_output.write((char*)&output_time, sizeof(double));
        f_output.write((char*)&rho_p_v_p_v_p[0], sizeof(double)*rho_p_v_p_v_p.size());
        
        f_output.close();
    }
}


/*
 * Output correlation of density and velocity component square in z-direction with inhomogeneous x-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::
outputDensityVelocityComponentSquareInZDirectionCorrelationWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double output_time) const
{
    std::vector<double> rho_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        patch_hierarchy,
        data_context);
    
    std::vector<double> w_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "VELOCITY",
        2,
        patch_hierarchy,
        data_context);
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    std::vector<std::vector<double> > averaged_quantities;
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(rho_mean);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(2);
    averaged_quantities.push_back(w_mean);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(2);
    averaged_quantities.push_back(w_mean);
    
    std::vector<double> rho_p_w_p_w_p = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_output;
        f_output.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        
        f_output.write((char*)&output_time, sizeof(double));
        f_output.write((char*)&rho_p_w_p_w_p[0], sizeof(double)*rho_p_w_p_w_p.size());
        
        f_output.close();
    }
}



/*
 * Output density variance with inhomogeneous x-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputDensityVarianceWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double output_time) const
{
    std::vector<double> rho_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        patch_hierarchy,
        data_context);
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    std::vector<std::vector<double> > averaged_quantities;
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(rho_mean);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(rho_mean);
    
    std::vector<double> rho_p_rho_p = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_output;
        f_output.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        
        f_output.write((char*)&output_time, sizeof(double));
        f_output.write((char*)&rho_p_rho_p[0], sizeof(double)*rho_p_rho_p.size());
        
        f_output.close();
    }
}


/*
 * Output density specific volume covariance with inhomogeneous x-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputDensitySpecificVolumeCovarianceWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double output_time) const
{
    std::vector<double> rho_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        patch_hierarchy,
        data_context);
    
    std::vector<double> v_mean = getAveragedReciprocalOfQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        patch_hierarchy,
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
    
    std::vector<double> rho_p_v_p = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        use_reciprocal,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_output;
        f_output.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        
        f_output.write((char*)&output_time, sizeof(double));
        f_output.write((char*)&rho_p_v_p[0], sizeof(double)*rho_p_v_p.size());
        
        f_output.close();
    }
}


/**
 ** Function to compute budgets.
 **/        

/*
 * Output turbulent mass flux in x-direction with inhomogeneous x-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputBudgetTurbMassFluxXWithInhomogeneousXDirection(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double output_time) const
{
    const int finest_level_dim_0 = getRefinedDomainNumberOfPointsX(patch_hierarchy);
    
    const double dx = getRefinedDomainGridSpacingX(patch_hierarchy);
    
    std::vector<std::string> quantity_names;
    std::vector<int> component_indices;
    std::vector<std::vector<double> > averaged_quantities;
    std::vector<bool> use_reciprocal;
    std::vector<int> derivative_directions;
    
    /*
     * Compute rho_a1.
     */
    
    std::vector<double> rho_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        patch_hierarchy,
        data_context);
    
    std::vector<double> u_mean = getAveragedQuantityWithInhomogeneousXDirection(
        "VELOCITY",
        0,
        patch_hierarchy,
        data_context);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(rho_mean);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(u_mean);
    
    std::vector<double> rho_p_u_p = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    averaged_quantities.clear();
    
    std::vector<double> a_1(rho_p_u_p);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        a_1[i] /= rho_mean[i];
    }
    
    /*
     * Compute term II.
     */
    
    // Compute u_tilde.
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    
    std::vector<double> rho_u_mean = getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        patch_hierarchy,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    
    std::vector<double> u_tilde(rho_u_mean);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        u_tilde[i] /= rho_mean[i];
    }
    
    std::vector<double> rho_u_tilde_a1(u_tilde);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        rho_u_tilde_a1[i] *= rho_p_u_p[i];
    }
    
    std::vector<double> d_rho_u_tilde_a1_dx = computeDerivativeOfVector1D(
        rho_u_tilde_a1,
        dx);
    
    /*
     * Compute term III(1).
     */
    
    std::vector<double> rho_inv_mean = getAveragedReciprocalOfQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        patch_hierarchy,
        data_context);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_reciprocal.push_back(false);
    averaged_quantities.push_back(rho_mean);
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_reciprocal.push_back(true);
    averaged_quantities.push_back(rho_inv_mean);
    
    std::vector<double> b = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        use_reciprocal,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    use_reciprocal.clear();
    averaged_quantities.clear();
    
    std::vector<double> dp_dx_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        "PRESSURE",
        0,
        0,
        patch_hierarchy,
        data_context);
    
    std::vector<double> b_dp_dx(dp_dx_mean);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        b_dp_dx[i] *= b[i];
    }
    
    /*
     * Compute term III(2).
     */
    
    // ...
    
    /*
     * Compute term III(3).
     */
    
    std::vector<double> drho_dx_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        "DENSITY",
        0,
        0,
        patch_hierarchy,
        data_context);
    
    // Compute R_11.
    
    std::vector<double> zeros(finest_level_dim_0, double(0));
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(zeros);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(u_tilde);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(u_tilde);
    
    std::vector<double> R_11 = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    averaged_quantities.clear();
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        R_11[i] /= rho_mean[i];
    }
    
    std::vector<double> m_R_11_drho_dx(drho_dx_mean);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        m_R_11_drho_dx[i] *= (-R_11[i]);
    }
    
    /*
     * Compute term IV(1).
     */
    
    std::vector<double> a_1_sq(a_1);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        a_1_sq[i] *= a_1[i];
    }
    
    std::vector<double> da_1_sq_dx = computeDerivativeOfVector1D(
        a_1_sq,
        dx);
    
    std::vector<double> rho_da_1_sq_dx(da_1_sq_dx);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        rho_da_1_sq_dx[i] *= rho_mean[i];
    }
    
    /*
     * Compute term IV(2).
     */
    
    std::vector<double> du_dx_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        "VELOCITY",
        0,
        0,
        patch_hierarchy,
        data_context);
    
    std::vector<double> m_rho_a_1_du_dx(du_dx_mean);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        m_rho_a_1_du_dx[i] *= (-rho_mean[i]*a_1[i]);
    }
    
    /*
     * Compute term V.
     */
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(rho_mean);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(u_mean);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    averaged_quantities.push_back(u_mean);
    
    std::vector<double> rho_p_u_p_sq_over_rho = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    averaged_quantities.clear();
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        rho_p_u_p_sq_over_rho[i] /= rho_mean[i];
    }
    
    std::vector<double> drho_p_u_p_sq_over_rho_dx = computeDerivativeOfVector1D(
        rho_p_u_p_sq_over_rho,
        dx);
    
    std::vector<double> m_rho_drho_p_u_p_sq_over_rho_dx(drho_p_u_p_sq_over_rho_dx);
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        m_rho_drho_p_u_p_sq_over_rho_dx[i] *= (-rho_mean[i]);
    }
    
    /*
     * Compute term VI(1).
     */
    
    quantity_names.push_back("DENSITY");
    component_indices.push_back(0);
    use_reciprocal.push_back(true);
    derivative_directions.push_back(-1);
    averaged_quantities.push_back(rho_inv_mean);
    
    quantity_names.push_back("PRESSURE");
    component_indices.push_back(0);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(0);
    averaged_quantities.push_back(dp_dx_mean);
    
    std::vector<double> rho_rho_inv_p_dp_dx_p = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        use_reciprocal,
        derivative_directions,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    use_reciprocal.clear();
    derivative_directions.clear();
    averaged_quantities.clear();
    
    for (int i = 0; i < finest_level_dim_0; i++)
    {
        rho_rho_inv_p_dp_dx_p[i] *= rho_mean[i];
    }
    
    /*
     * Compute term VI(3).
     */
    
    // Repeated.
    // std::vector<double> du_dx_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
    //     "VELOCITY",
    //     0,
    //     0,
    //     patch_hierarchy,
    //     data_context);
    
    std::vector<double> dv_dy_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
        "VELOCITY",
        1,
        1,
        patch_hierarchy,
        data_context);
    
    std::vector<double> dw_dz_mean;
    if (d_dim == tbox::Dimension(3))
    {
        dw_dz_mean = getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
            "VELOCITY",
            2,
            2,
            patch_hierarchy,
            data_context);
    }
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(-1);
    averaged_quantities.push_back(u_mean);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(0);
    averaged_quantities.push_back(du_dx_mean);
    
    std::vector<double> epsilon_a1_1 = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        use_reciprocal,
        derivative_directions,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    use_reciprocal.clear();
    derivative_directions.clear();
    averaged_quantities.clear();
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(0);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(-1);
    averaged_quantities.push_back(u_mean);
    
    quantity_names.push_back("VELOCITY");
    component_indices.push_back(1);
    use_reciprocal.push_back(false);
    derivative_directions.push_back(1);
    averaged_quantities.push_back(dv_dy_mean);
    
    std::vector<double> epsilon_a1_2 = getQuantityCorrelationWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        use_reciprocal,
        derivative_directions,
        averaged_quantities,
        patch_hierarchy,
        data_context);
    
    quantity_names.clear();
    component_indices.clear();
    use_reciprocal.clear();
    derivative_directions.clear();
    averaged_quantities.clear();
    
    std::vector<double> epsilon_a1_3;
    if (d_dim == tbox::Dimension(3))
    {
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(0);
        use_reciprocal.push_back(false);
        derivative_directions.push_back(-1);
        averaged_quantities.push_back(u_mean);
        
        quantity_names.push_back("VELOCITY");
        component_indices.push_back(2);
        use_reciprocal.push_back(false);
        derivative_directions.push_back(2);
        averaged_quantities.push_back(dw_dz_mean);
        
        epsilon_a1_3 = getQuantityCorrelationWithInhomogeneousXDirection(
            quantity_names,
            component_indices,
            use_reciprocal,
            derivative_directions,
            averaged_quantities,
            patch_hierarchy,
            data_context);
        
        quantity_names.clear();
        component_indices.clear();
        use_reciprocal.clear();
        derivative_directions.clear();
        averaged_quantities.clear();
    }
    
    std::vector<double> rho_epsilon_a1(finest_level_dim_0, double(0));
    
    if (d_dim == tbox::Dimension(2))
    {
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            rho_epsilon_a1[i] = -rho_mean[i]*(epsilon_a1_1[i] + epsilon_a1_2[i]);
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            rho_epsilon_a1[i] = -rho_mean[i]*(epsilon_a1_1[i] + epsilon_a1_2[i] + epsilon_a1_3[i]);
        }
    }
    
    
    /*
     * Output budget.
     */
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_output;
        f_output.open(stat_dump_filename, std::ios_base::app | std::ios::out | std::ios::binary);
        
        f_output.write((char*)&output_time, sizeof(double));
        f_output.write((char*)&rho_p_u_p[0], sizeof(double)*rho_p_u_p.size());
        // Term II.
        f_output.write((char*)&d_rho_u_tilde_a1_dx[0], sizeof(double)*d_rho_u_tilde_a1_dx.size());
        
        // Term III(1).
        f_output.write((char*)&b_dp_dx[0], sizeof(double)*b_dp_dx.size());
        // Term III(3).
        f_output.write((char*)&m_R_11_drho_dx[0], sizeof(double)*m_R_11_drho_dx.size());
        
        // Term IV(1).
        f_output.write((char*)&rho_da_1_sq_dx[0], sizeof(double)*rho_da_1_sq_dx.size());
        // Term IV(2).
        f_output.write((char*)&m_rho_a_1_du_dx[0], sizeof(double)*m_rho_a_1_du_dx.size());
        
        // Term V.
        f_output.write((char*)&m_rho_drho_p_u_p_sq_over_rho_dx[0], sizeof(double)*m_rho_drho_p_u_p_sq_over_rho_dx.size());
        
        // Term VI(1).
        f_output.write((char*)&rho_rho_inv_p_dp_dx_p[0], sizeof(double)*rho_rho_inv_p_dp_dx_p.size());
        // Term VI(3).
        f_output.write((char*)&rho_epsilon_a1[0], sizeof(double)*rho_epsilon_a1.size());
        
        f_output.close();
    }
}


/**
 ** Helper functions.
 **/        


/*
 * Get number of points in the x-direction of the refined domain.
 */
const int
FlowModelStatisticsUtilitiesFourEqnConservative::getRefinedDomainNumberOfPointsX(
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy) const
{
    /*
     * Get the refinement ratio from the finest level to the coarest level.
     */
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    hier::IntVector ratioFinestLevelToCoarestLevel =
        patch_hierarchy->getRatioToCoarserLevel(num_levels - 1);
    for (int li = num_levels - 2; li > 0 ; li--)
    {
        ratioFinestLevelToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(li);
    }
    
    /*
     * Get the number of cells of physical domain refined to the finest level.
     */
    
    const hier::BoxContainer& physical_domain = d_grid_geometry->getPhysicalDomain();
    const hier::Box& physical_domain_box = physical_domain.front();
    const hier::IntVector& physical_domain_dims = physical_domain_box.numberCells();
    const hier::IntVector finest_level_dims = physical_domain_dims*ratioFinestLevelToCoarestLevel;
    
    const int finest_level_dim_0 = finest_level_dims[0];
    return finest_level_dim_0;
}


/*
 * Get grid spacing in the x-direction of the refined domain.
 */
const double
FlowModelStatisticsUtilitiesFourEqnConservative::getRefinedDomainGridSpacingX(
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy) const
{
    /*
     * Get the refinement ratio from the finest level to the coarest level.
     */
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    hier::IntVector ratioFinestLevelToCoarestLevel =
        patch_hierarchy->getRatioToCoarserLevel(num_levels - 1);
    for (int li = num_levels - 2; li > 0 ; li--)
    {
        ratioFinestLevelToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(li);
    }
    
    const double* dx = d_grid_geometry->getDx();
    
    return dx[0]/ratioFinestLevelToCoarestLevel[0];
}


/*
 * Compute the one-dimensional derivative given a vector.
 */
std::vector<double> FlowModelStatisticsUtilitiesFourEqnConservative::computeDerivativeOfVector1D(
    const std::vector<double> quantity_vector,
    const double dx) const
{
    const int vector_length = quantity_vector.size();
    
    std::vector<double> derivative;
    derivative.resize(vector_length);
    
    const double* u = quantity_vector.data();
    double* dudx = derivative.data();
    
    // HARD CODE TO BE SIXTH ORDER CENTRAL SCHEME FOR DIFFERENTIATION.
    for (int i = 3; i < vector_length - 3; i++)
    {
        // Compute linear indices.
        const int idx     = i;
        
        const int idx_LLL = i - 3;
        const int idx_LL  = i - 2;
        const int idx_L   = i - 1;
        const int idx_R   = i + 1;
        const int idx_RR  = i + 2;
        const int idx_RRR = i + 3;
        
        dudx[idx] = (double(1)/double(60)*(u[idx_RRR] - u[idx_LLL])
            - double(3)/double(20)*(u[idx_RR] - u[idx_LL])
            + double(3)/double(4)*(u[idx_R] - u[idx_L]))/dx;
    }
    
    for (int i = 0; i < 3; i++)
    {
        dudx[i]                     = dudx[3];
        dudx[vector_length - i - 1] = dudx[vector_length - 4];
    }
    
    return derivative;
}


/*
 * Compute averaged value with only x direction as inhomogeneous direction.
 */
std::vector<double>
FlowModelStatisticsUtilitiesFourEqnConservative::getAveragedQuantityWithInhomogeneousXDirection(
    const std::string quantity_name,
    const int component_idx,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context) const
{
    std::vector<double> averaged_quantity;
    
    boost::shared_ptr<FlowModel> d_flow_model_tmp = d_flow_model.lock();
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Get the refinement ratio from the finest level to the coarest level.
     */
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    hier::IntVector ratioFinestLevelToCoarestLevel =
        patch_hierarchy->getRatioToCoarserLevel(num_levels - 1);
    for (int li = num_levels - 2; li > 0 ; li--)
    {
        ratioFinestLevelToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(li);
    }
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    boost::shared_ptr<ExtendedFlattenedHierarchy> flattened_hierarchy(
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
    const hier::IntVector finest_level_dims = physical_domain_dims*ratioFinestLevelToCoarestLevel;
    
    /*
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        double* u_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_quantity.resize(finest_level_dim_0);
        double* u_avg_global = averaged_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_avg_local[i] = 0.0;
            u_avg_global[i] = 0.0;
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_quantity =
                    d_flow_model_tmp->getGlobalCellData(quantity_name);
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                
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
                         * Compute the linear index and the data to add.
                         */
                        
                        const int idx = relative_idx_lo_0 + i + num_ghosts_0_quantity;
                        
                        const double value_to_add = u[idx]/((double) n_overlapped);
                        
                        for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                        {
                            const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                            
                            u_avg_local[idx_fine] += value_to_add;
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
            u_avg_local,
            u_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_avg_local);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        
        double* u_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_quantity.resize(finest_level_dim_0);
        double* u_avg_global = averaged_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_avg_local[i] = 0.0;
            u_avg_global[i] = 0.0;
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_quantity =
                    d_flow_model_tmp->getGlobalCellData(quantity_name);
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                
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
                             * Compute the linear index and the data to add.
                             */
                            
                            const int idx = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                            
                            const double value_to_add = u[idx]*weight/((double) n_overlapped);
                            
                            for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                            {
                                const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                
                                u_avg_local[idx_fine] += value_to_add;
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
            u_avg_local,
            u_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_avg_local);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        const double L_z = x_hi[2] - x_lo[2];
        
        double* u_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_quantity.resize(finest_level_dim_0);
        double* u_avg_global = averaged_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_avg_local[i] = 0.0;
            u_avg_global[i] = 0.0;
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_quantity =
                    d_flow_model_tmp->getGlobalCellData(quantity_name);
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int num_ghosts_2_quantity = num_ghosts_quantity[2];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                const int ghostcell_dim_1_quantity = ghostcell_dims_quantity[1];
                
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
                                 * Compute the linear index and the data to add.
                                 */
                                
                                const int idx = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                    (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                        ghostcell_dim_1_quantity;
                                
                                const double value_to_add = u[idx]*weight/((double) n_overlapped);
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    u_avg_local[idx_fine] += value_to_add;
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
            u_avg_local,
            u_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_avg_local);
    }
    
    return averaged_quantity;
}


/*
 * Compute averaged reciprocal of value with only x direction as inhomogeneous direction.
 */
std::vector<double>
FlowModelStatisticsUtilitiesFourEqnConservative::getAveragedReciprocalOfQuantityWithInhomogeneousXDirection(
    const std::string quantity_name,
    const int component_idx,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context) const
{
    std::vector<double> averaged_reciprocal_quantity;
    
    boost::shared_ptr<FlowModel> d_flow_model_tmp = d_flow_model.lock();
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Get the refinement ratio from the finest level to the coarest level.
     */
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    hier::IntVector ratioFinestLevelToCoarestLevel =
        patch_hierarchy->getRatioToCoarserLevel(num_levels - 1);
    for (int li = num_levels - 2; li > 0 ; li--)
    {
        ratioFinestLevelToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(li);
    }
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    boost::shared_ptr<ExtendedFlattenedHierarchy> flattened_hierarchy(
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
    const hier::IntVector finest_level_dims = physical_domain_dims*ratioFinestLevelToCoarestLevel;
    
    /*
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        double* u_inv_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_reciprocal_quantity.resize(finest_level_dim_0);
        double* u_inv_avg_global = averaged_reciprocal_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_inv_avg_local[i] = 0.0;
            u_inv_avg_global[i] = 0.0;
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_quantity =
                    d_flow_model_tmp->getGlobalCellData(quantity_name);
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                
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
                         * Compute the linear index and the data to add.
                         */
                        
                        const int idx = relative_idx_lo_0 + i + num_ghosts_0_quantity;
                        
                        const double value_to_add = (double(1)/u[idx])/((double) n_overlapped);
                        
                        for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                        {
                            const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                            
                            u_inv_avg_local[idx_fine] += value_to_add;
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
         * Reduction to get the global average of reciprocal.
         */
        
        mpi.Allreduce(
            u_inv_avg_local,
            u_inv_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_inv_avg_local);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        
        double* u_inv_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_reciprocal_quantity.resize(finest_level_dim_0);
        double* u_inv_avg_global = averaged_reciprocal_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_inv_avg_local[i] = 0.0;
            u_inv_avg_global[i] = 0.0;
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_quantity =
                    d_flow_model_tmp->getGlobalCellData(quantity_name);
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                
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
                             * Compute the linear index and the data to add.
                             */
                            
                            const int idx = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                            
                            const double value_to_add = (double(1)/u[idx])*weight/((double) n_overlapped);
                            
                            for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                            {
                                const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                
                                u_inv_avg_local[idx_fine] += value_to_add;
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
         * Reduction to get the global average of reciprocal.
         */
        
        mpi.Allreduce(
            u_inv_avg_local,
            u_inv_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_inv_avg_local);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        const double L_z = x_hi[2] - x_lo[2];
        
        double* u_inv_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_reciprocal_quantity.resize(finest_level_dim_0);
        double* u_inv_avg_global = averaged_reciprocal_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_inv_avg_local[i] = 0.0;
            u_inv_avg_global[i] = 0.0;
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_quantity =
                    d_flow_model_tmp->getGlobalCellData(quantity_name);
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int num_ghosts_2_quantity = num_ghosts_quantity[2];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                const int ghostcell_dim_1_quantity = ghostcell_dims_quantity[1];
                
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
                                 * Compute the linear index and the data to add.
                                 */
                                
                                const int idx = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                    (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                        ghostcell_dim_1_quantity;
                                
                                const double value_to_add = (double(1)/u[idx])*weight/((double) n_overlapped);
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    u_inv_avg_local[idx_fine] += value_to_add;
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
         * Reduction to get the global average of reciprocal.
         */
        
        mpi.Allreduce(
            u_inv_avg_local,
            u_inv_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_inv_avg_local);
    }
    
    return averaged_reciprocal_quantity;
}


/*
 * Compute averaged derivative of value with only x direction as inhomogeneous direction.
 */
std::vector<double>
FlowModelStatisticsUtilitiesFourEqnConservative::getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
    const std::string quantity_name,
    const int component_idx,
    const int derivative_direction,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context) const
{
    std::vector<double> averaged_derivative_quantity;
    
    boost::shared_ptr<FlowModel> d_flow_model_tmp = d_flow_model.lock();
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Get the refinement ratio from the finest level to the coarest level.
     */
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    hier::IntVector ratioFinestLevelToCoarestLevel =
        patch_hierarchy->getRatioToCoarserLevel(num_levels - 1);
    for (int li = num_levels - 2; li > 0 ; li--)
    {
        ratioFinestLevelToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(li);
    }
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    boost::shared_ptr<ExtendedFlattenedHierarchy> flattened_hierarchy(
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
    const hier::IntVector finest_level_dims = physical_domain_dims*ratioFinestLevelToCoarestLevel;
    
    /*
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        double* u_der_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_derivative_quantity.resize(finest_level_dim_0);
        double* u_der_avg_global = averaged_derivative_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_der_avg_local[i] = 0.0;
            u_der_avg_global[i] = 0.0;
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                // HARD CODE TO BE SIXTH ORDER CENTRAL SCHEME FOR DIFFERENTIATION.
                TBOX_ASSERT(num_ghosts >= hier::IntVector::getOne(d_dim)*3);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, num_ghosts));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_quantity =
                    d_flow_model_tmp->getGlobalCellData(quantity_name);
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                
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
                         * Compute the derivative.
                         */
                        
                        double value_to_add = double(0);
                        
                        if (derivative_direction == 0)
                        {
                            // Compute linear indices.
                            const int idx_LLL = relative_idx_lo_0 + (i - 3) + num_ghosts_0_quantity;
                            const int idx_LL  = relative_idx_lo_0 + (i - 2) + num_ghosts_0_quantity;
                            const int idx_L   = relative_idx_lo_0 + (i - 1) + num_ghosts_0_quantity;
                            const int idx_R   = relative_idx_lo_0 + (i + 1) + num_ghosts_0_quantity;
                            const int idx_RR  = relative_idx_lo_0 + (i + 2) + num_ghosts_0_quantity;
                            const int idx_RRR = relative_idx_lo_0 + (i + 3) + num_ghosts_0_quantity;
                            
                            const double dudx = (double(1)/double(60)*(u[idx_RRR] - u[idx_LLL])
                                - double(3)/double(20)*(u[idx_RR] - u[idx_LL])
                                + double(3)/double(4)*(u[idx_R] - u[idx_L]))/dx[0];
                            
                            value_to_add = dudx/((double) n_overlapped);
                        }
                        else
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "Cannot take derivative for one-dimensional problem!\n"
                                << "derivative_direction = " << derivative_direction << " given!\n"
                                << std::endl);
                        }
                        
                        /*
                         * Add the data.
                         */
                        
                        for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                        {
                            const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                            
                            u_der_avg_local[idx_fine] += value_to_add;
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
         * Reduction to get the global average of derivative.
         */
        
        mpi.Allreduce(
            u_der_avg_local,
            u_der_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_der_avg_local);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        
        double* u_der_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_derivative_quantity.resize(finest_level_dim_0);
        double* u_der_avg_global = averaged_derivative_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_der_avg_local[i] = 0.0;
            u_der_avg_global[i] = 0.0;
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                // HARD CODE TO BE SIXTH ORDER CENTRAL SCHEME FOR DIFFERENTIATION.
                TBOX_ASSERT(num_ghosts >= hier::IntVector::getOne(d_dim)*3);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, num_ghosts));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_quantity =
                    d_flow_model_tmp->getGlobalCellData(quantity_name);
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                
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
                             * Compute the derivative.
                             */
                            
                            double value_to_add = double(0);
                            
                            if (derivative_direction == 0)
                            {
                                const int idx_LLL = (relative_idx_lo_0 + (i - 3) + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const int idx_LL  = (relative_idx_lo_0 + (i - 2) + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const int idx_L   = (relative_idx_lo_0 + (i - 1) + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const int idx_R   = (relative_idx_lo_0 + (i + 1) + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const int idx_RR  = (relative_idx_lo_0 + (i + 2) + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const int idx_RRR = (relative_idx_lo_0 + (i + 3) + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const double dudx = (double(1)/double(60)*(u[idx_RRR] - u[idx_LLL])
                                    - double(3)/double(20)*(u[idx_RR] - u[idx_LL])
                                    + double(3)/double(4)*(u[idx_R] - u[idx_L]))/dx[0];
                                
                                value_to_add = dudx*weight/((double) n_overlapped);
                            }
                            else if (derivative_direction == 1)
                            {
                                const int idx_BBB = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + (j - 3) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const int idx_BB  = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + (j - 2) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const int idx_B   = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + (j - 1) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const int idx_T   = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + (j + 1) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const int idx_TT  = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + (j + 2) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const int idx_TTT = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + (j + 3) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const double dudy = (double(1)/double(60)*(u[idx_TTT] - u[idx_BBB])
                                    - double(3)/double(20)*(u[idx_TT] - u[idx_BB])
                                    + double(3)/double(4)*(u[idx_T] - u[idx_B]))/dx[1];
                                
                                value_to_add = dudy*weight/((double) n_overlapped);
                            }
                            else
                            {
                                TBOX_ERROR(d_object_name
                                    << ": "
                                    << "Cannot take derivative for two-dimensional problem!\n"
                                    << "derivative_direction = " << derivative_direction << " given!\n"
                                    << std::endl);
                            }
                            
                            /*
                             * Add the data.
                             */
                            
                            for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                            {
                                const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                
                                u_der_avg_local[idx_fine] += value_to_add;
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
         * Reduction to get the global average of derivative.
         */
        
        mpi.Allreduce(
            u_der_avg_local,
            u_der_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_der_avg_local);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        const double L_z = x_hi[2] - x_lo[2];
        
        double* u_der_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_derivative_quantity.resize(finest_level_dim_0);
        double* u_der_avg_global = averaged_derivative_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_der_avg_local[i] = 0.0;
            u_der_avg_global[i] = 0.0;
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                // HARD CODE TO BE SIXTH ORDER CENTRAL SCHEME FOR DIFFERENTIATION.
                TBOX_ASSERT(num_ghosts >= hier::IntVector::getOne(d_dim)*3);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, num_ghosts));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_quantity =
                    d_flow_model_tmp->getGlobalCellData(quantity_name);
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int num_ghosts_2_quantity = num_ghosts_quantity[2];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                const int ghostcell_dim_1_quantity = ghostcell_dims_quantity[1];
                
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
                                 * Compute the derivative.
                                 */
                                
                                double value_to_add = double(0);
                                
                                if (derivative_direction == 0)
                                {
                                    const int idx_LLL = (relative_idx_lo_0 + (i - 3) + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_LL  = (relative_idx_lo_0 + (i - 2) + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_L   = (relative_idx_lo_0 + (i - 1) + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_R   = (relative_idx_lo_0 + (i + 1) + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_RR  = (relative_idx_lo_0 + (i + 2) + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_RRR = (relative_idx_lo_0 + (i + 3) + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const double dudx = (double(1)/double(60)*(u[idx_RRR] - u[idx_LLL])
                                        - double(3)/double(20)*(u[idx_RR] - u[idx_LL])
                                        + double(3)/double(4)*(u[idx_R] - u[idx_L]))/dx[0];
                                    
                                    value_to_add = dudx*weight/((double) n_overlapped);
                                }
                                else if (derivative_direction == 1)
                                {
                                    const int idx_BBB = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + (j - 3) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_BB  = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + (j - 2) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_B   = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + (j - 1) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_T   = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + (j + 1) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_TT  = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + (j + 2) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_TTT = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + (j + 3) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const double dudy = (double(1)/double(60)*(u[idx_TTT] - u[idx_BBB])
                                        - double(3)/double(20)*(u[idx_TT] - u[idx_BB])
                                        + double(3)/double(4)*(u[idx_T] - u[idx_B]))/dx[1];
                                    
                                    value_to_add = dudy*weight/((double) n_overlapped);
                                }
                                else if (derivative_direction == 2)
                                {
                                    const int idx_BBB = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + (k - 3) + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_BB  = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + (k - 2) + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_B   = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + (k - 1) + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_F   = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + (k + 1) + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_FF  = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + (k + 2) + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_FFF = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + (k + 3) + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const double dudz = (double(1)/double(60)*(u[idx_FFF] - u[idx_BBB])
                                        - double(3)/double(20)*(u[idx_FF] - u[idx_BB])
                                        + double(3)/double(4)*(u[idx_F] - u[idx_B]))/dx[2];
                                    
                                    value_to_add = dudz*weight/((double) n_overlapped);
                                }
                                else
                                {
                                    TBOX_ERROR(d_object_name
                                        << ": "
                                        << "Cannot take derivative for three-dimensional problem!\n"
                                        << "derivative_direction = " << derivative_direction << " given!\n"
                                        << std::endl);
                                }
                                
                                /*
                                 * Add the data.
                                 */
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    u_der_avg_local[idx_fine] += value_to_add;
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
         * Reduction to get the global average of derivative.
         */
        
        mpi.Allreduce(
            u_der_avg_local,
            u_der_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_der_avg_local);
    }
    
    return averaged_derivative_quantity;
}


/*
 * Compute averaged derivative of reciprocal of value with only x direction as inhomogeneous direction.
 */
std::vector<double>
FlowModelStatisticsUtilitiesFourEqnConservative::
getAveragedDerivativeOfReciprocalOfQuantityWithInhomogeneousXDirection(
    const std::string quantity_name,
    const int component_idx,
    const int derivative_direction,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context) const
{
    std::vector<double> averaged_derivative_reciprocal_quantity;
    
    boost::shared_ptr<FlowModel> d_flow_model_tmp = d_flow_model.lock();
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Get the refinement ratio from the finest level to the coarest level.
     */
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    hier::IntVector ratioFinestLevelToCoarestLevel =
        patch_hierarchy->getRatioToCoarserLevel(num_levels - 1);
    for (int li = num_levels - 2; li > 0 ; li--)
    {
        ratioFinestLevelToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(li);
    }
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    boost::shared_ptr<ExtendedFlattenedHierarchy> flattened_hierarchy(
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
    const hier::IntVector finest_level_dims = physical_domain_dims*ratioFinestLevelToCoarestLevel;
    
    /*
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        double* u_inv_der_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_derivative_reciprocal_quantity.resize(finest_level_dim_0);
        double* u_inv_der_avg_global = averaged_derivative_reciprocal_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_inv_der_avg_local[i] = 0.0;
            u_inv_der_avg_global[i] = 0.0;
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                // HARD CODE TO BE SIXTH ORDER CENTRAL SCHEME FOR DIFFERENTIATION.
                TBOX_ASSERT(num_ghosts >= hier::IntVector::getOne(d_dim)*3);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, num_ghosts));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_quantity =
                    d_flow_model_tmp->getGlobalCellData(quantity_name);
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                
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
                         * Compute the derivative.
                         */
                        
                        double value_to_add = double(0);
                        
                        if (derivative_direction == 0)
                        {
                            // Compute linear indices.
                            const int idx_LLL = relative_idx_lo_0 + (i - 3) + num_ghosts_0_quantity;
                            const int idx_LL  = relative_idx_lo_0 + (i - 2) + num_ghosts_0_quantity;
                            const int idx_L   = relative_idx_lo_0 + (i - 1) + num_ghosts_0_quantity;
                            const int idx_R   = relative_idx_lo_0 + (i + 1) + num_ghosts_0_quantity;
                            const int idx_RR  = relative_idx_lo_0 + (i + 2) + num_ghosts_0_quantity;
                            const int idx_RRR = relative_idx_lo_0 + (i + 3) + num_ghosts_0_quantity;
                            
                            const double dudx = (double(1)/double(60)*(double(1)/u[idx_RRR] - double(1)/u[idx_LLL])
                                - double(3)/double(20)*(double(1)/u[idx_RR] - double(1)/u[idx_LL])
                                + double(3)/double(4)*(double(1)/u[idx_R] - double(1)/u[idx_L]))/dx[0];
                            
                            value_to_add = dudx/((double) n_overlapped);
                        }
                        else
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "Cannot take derivative for one-dimensional problem!\n"
                                << "derivative_direction = " << derivative_direction << " given!\n"
                                << std::endl);
                        }
                        
                        /*
                         * Add the data.
                         */
                        
                        for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                        {
                            const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                            
                            u_inv_der_avg_local[idx_fine] += value_to_add;
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
         * Reduction to get the global average of derivative of reciprocal.
         */
        
        mpi.Allreduce(
            u_inv_der_avg_local,
            u_inv_der_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_inv_der_avg_local);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        
        double* u_inv_der_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_derivative_reciprocal_quantity.resize(finest_level_dim_0);
        double* u_inv_der_avg_global = averaged_derivative_reciprocal_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_inv_der_avg_local[i] = 0.0;
            u_inv_der_avg_global[i] = 0.0;
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                // HARD CODE TO BE SIXTH ORDER CENTRAL SCHEME FOR DIFFERENTIATION.
                TBOX_ASSERT(num_ghosts >= hier::IntVector::getOne(d_dim)*3);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, num_ghosts));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_quantity =
                    d_flow_model_tmp->getGlobalCellData(quantity_name);
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                
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
                             * Compute the derivative.
                             */
                            
                            double value_to_add = double(0);
                            
                            if (derivative_direction == 0)
                            {
                                const int idx_LLL = (relative_idx_lo_0 + (i - 3) + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const int idx_LL  = (relative_idx_lo_0 + (i - 2) + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const int idx_L   = (relative_idx_lo_0 + (i - 1) + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const int idx_R   = (relative_idx_lo_0 + (i + 1) + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const int idx_RR  = (relative_idx_lo_0 + (i + 2) + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const int idx_RRR = (relative_idx_lo_0 + (i + 3) + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const double dudx = (double(1)/double(60)*(double(1)/u[idx_RRR] - double(1)/u[idx_LLL])
                                    - double(3)/double(20)*(double(1)/u[idx_RR] - double(1)/u[idx_LL])
                                    + double(3)/double(4)*(double(1)/u[idx_R] - double(1)/u[idx_L]))/dx[0];
                                
                                value_to_add = dudx*weight/((double) n_overlapped);
                            }
                            else if (derivative_direction == 1)
                            {
                                const int idx_BBB = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + (j - 3) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const int idx_BB  = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + (j - 2) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const int idx_B   = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + (j - 1) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const int idx_T   = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + (j + 1) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const int idx_TT  = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + (j + 2) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const int idx_TTT = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + (j + 3) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                                
                                const double dudy = (double(1)/double(60)*(double(1)/u[idx_TTT] - double(1)/u[idx_BBB])
                                    - double(3)/double(20)*(double(1)/u[idx_TT] - double(1)/u[idx_BB])
                                    + double(3)/double(4)*(double(1)/u[idx_T] - double(1)/u[idx_B]))/dx[1];
                                
                                value_to_add = dudy*weight/((double) n_overlapped);
                            }
                            else
                            {
                                TBOX_ERROR(d_object_name
                                    << ": "
                                    << "Cannot take derivative for two-dimensional problem!\n"
                                    << "derivative_direction = " << derivative_direction << " given!\n"
                                    << std::endl);
                            }
                            
                            /*
                             * Add the data.
                             */
                            
                            for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                            {
                                const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                
                                u_inv_der_avg_local[idx_fine] += value_to_add;
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
         * Reduction to get the global average of derivative of reciprocal.
         */
        
        mpi.Allreduce(
            u_inv_der_avg_local,
            u_inv_der_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_inv_der_avg_local);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        const double L_z = x_hi[2] - x_lo[2];
        
        double* u_inv_der_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_derivative_reciprocal_quantity.resize(finest_level_dim_0);
        double* u_inv_der_avg_global = averaged_derivative_reciprocal_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_inv_der_avg_local[i] = 0.0;
            u_inv_der_avg_global[i] = 0.0;
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                // HARD CODE TO BE SIXTH ORDER CENTRAL SCHEME FOR DIFFERENTIATION.
                TBOX_ASSERT(num_ghosts >= hier::IntVector::getOne(d_dim)*3);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, num_ghosts));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_quantity =
                    d_flow_model_tmp->getGlobalCellData(quantity_name);
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int num_ghosts_2_quantity = num_ghosts_quantity[2];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                const int ghostcell_dim_1_quantity = ghostcell_dims_quantity[1];
                
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
                                 * Compute the derivative.
                                 */
                                
                                double value_to_add = double(0);
                                
                                if (derivative_direction == 0)
                                {
                                    const int idx_LLL = (relative_idx_lo_0 + (i - 3) + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_LL  = (relative_idx_lo_0 + (i - 2) + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_L   = (relative_idx_lo_0 + (i - 1) + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_R   = (relative_idx_lo_0 + (i + 1) + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_RR  = (relative_idx_lo_0 + (i + 2) + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_RRR = (relative_idx_lo_0 + (i + 3) + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const double dudx = (double(1)/double(60)*(double(1)/u[idx_RRR] - double(1)/u[idx_LLL])
                                        - double(3)/double(20)*(double(1)/u[idx_RR] - double(1)/u[idx_LL])
                                        + double(3)/double(4)*(double(1)/u[idx_R] - double(1)/u[idx_L]))/dx[0];
                                    
                                    value_to_add = dudx*weight/((double) n_overlapped);
                                }
                                else if (derivative_direction == 1)
                                {
                                    const int idx_BBB = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + (j - 3) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_BB  = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + (j - 2) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_B   = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + (j - 1) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_T   = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + (j + 1) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_TT  = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + (j + 2) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_TTT = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + (j + 3) + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const double dudy = (double(1)/double(60)*(double(1)/u[idx_TTT] - double(1)/u[idx_BBB])
                                        - double(3)/double(20)*(double(1)/u[idx_TT] - double(1)/u[idx_BB])
                                        + double(3)/double(4)*(double(1)/u[idx_T] - double(1)/u[idx_B]))/dx[1];
                                    
                                    value_to_add = dudy*weight/((double) n_overlapped);
                                }
                                else if (derivative_direction == 2)
                                {
                                    const int idx_BBB = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + (k - 3) + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_BB  = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + (k - 2) + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_B   = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + (k - 1) + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_F   = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + (k + 1) + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_FF  = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + (k + 2) + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const int idx_FFF = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                        (relative_idx_lo_2 + (k + 3) + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                            ghostcell_dim_1_quantity;
                                    
                                    const double dudz = (double(1)/double(60)*(double(1)/u[idx_FFF] - double(1)/u[idx_BBB])
                                        - double(3)/double(20)*(double(1)/u[idx_FF] - double(1)/u[idx_BB])
                                        + double(3)/double(4)*(double(1)/u[idx_F] - double(1)/u[idx_B]))/dx[2];
                                    
                                    value_to_add = dudz*weight/((double) n_overlapped);
                                }
                                else
                                {
                                    TBOX_ERROR(d_object_name
                                        << ": "
                                        << "Cannot take derivative for three-dimensional problem!\n"
                                        << "derivative_direction = " << derivative_direction << " given!\n"
                                        << std::endl);
                                }
                                
                                /*
                                 * Add the data.
                                 */
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    u_inv_der_avg_local[idx_fine] += value_to_add;
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
         * Reduction to get the global average of derivative of reciprocal.
         */
        
        mpi.Allreduce(
            u_inv_der_avg_local,
            u_inv_der_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_inv_der_avg_local);
    }
    
    return averaged_derivative_reciprocal_quantity;
}


/*
 * Compute averaged derivative of shear stress component with only x direction as inhomogeneous direction.
 * Component index:
 * 0: tau_11
 * 1: tau_12
 * 2: tau_13
 * 3: tau_22
 * 4: tau_23
 * 5: tau_33
 */
std::vector<double> 
FlowModelStatisticsUtilitiesFourEqnConservative::getAveragedShearStressComponentWithInhomogeneousXDirection(
    const int component_idx,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context) const
{
    std::vector<double> averaged_tau_ij;
    
    boost::shared_ptr<FlowModel> d_flow_model_tmp = d_flow_model.lock();
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Get the refinement ratio from the finest level to the coarest level.
     */
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    hier::IntVector ratioFinestLevelToCoarestLevel =
        patch_hierarchy->getRatioToCoarserLevel(num_levels - 1);
    for (int li = num_levels - 2; li > 0 ; li--)
    {
        ratioFinestLevelToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(li);
    }
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    boost::shared_ptr<ExtendedFlattenedHierarchy> flattened_hierarchy(
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
    const hier::IntVector finest_level_dims = physical_domain_dims*ratioFinestLevelToCoarestLevel;
    
    /*
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        double* tau_ij_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_tau_ij.resize(finest_level_dim_0);
        double* tau_ij_avg_global = averaged_tau_ij.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            tau_ij_avg_local[i] = 0.0;
            tau_ij_avg_global[i] = 0.0;
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                // HARD CODE TO BE SIXTH ORDER CENTRAL SCHEME FOR DIFFERENTIATION.
                TBOX_ASSERT(num_ghosts >= hier::IntVector::getOne(d_dim)*3);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", num_ghosts));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTION", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("PRESSURE", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("TEMPERATURE", hier::IntVector::getZero(d_dim)));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_velocity =
                    d_flow_model_tmp->getGlobalCellData("VELOCITY");
                
                boost::shared_ptr<pdat::CellData<double> > data_mass_fraction =
                    d_flow_model_tmp->getGlobalCellData("MASS_FRACTION");
                
                boost::shared_ptr<pdat::CellData<double> > data_pressure =
                    d_flow_model_tmp->getGlobalCellData("PRESSURE");
                
                boost::shared_ptr<pdat::CellData<double> > data_temperature =
                    d_flow_model_tmp->getGlobalCellData("TEMPERATURE");
                
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
                        
                        for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                        {
                            const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                            
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
            tau_ij_avg_local[i] = 0.0;
            tau_ij_avg_global[i] = 0.0;
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                // HARD CODE TO BE SIXTH ORDER CENTRAL SCHEME FOR DIFFERENTIATION.
                TBOX_ASSERT(num_ghosts >= hier::IntVector::getOne(d_dim)*3);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", num_ghosts));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTION", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("PRESSURE", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("TEMPERATURE", hier::IntVector::getZero(d_dim)));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_velocity =
                    d_flow_model_tmp->getGlobalCellData("VELOCITY");
                
                boost::shared_ptr<pdat::CellData<double> > data_mass_fraction =
                    d_flow_model_tmp->getGlobalCellData("MASS_FRACTION");
                
                boost::shared_ptr<pdat::CellData<double> > data_pressure =
                    d_flow_model_tmp->getGlobalCellData("PRESSURE");
                
                boost::shared_ptr<pdat::CellData<double> > data_temperature =
                    d_flow_model_tmp->getGlobalCellData("TEMPERATURE");
                
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
                             * Compute the the shear stress component.
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
                            
                            for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                            {
                                const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                
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
            tau_ij_avg_local[i] = 0.0;
            tau_ij_avg_global[i] = 0.0;
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                // HARD CODE TO BE SIXTH ORDER CENTRAL SCHEME FOR DIFFERENTIATION.
                TBOX_ASSERT(num_ghosts >= hier::IntVector::getOne(d_dim)*3);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", num_ghosts));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTION", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("PRESSURE", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("TEMPERATURE", hier::IntVector::getZero(d_dim)));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_velocity =
                    d_flow_model_tmp->getGlobalCellData("VELOCITY");
                
                boost::shared_ptr<pdat::CellData<double> > data_mass_fraction =
                    d_flow_model_tmp->getGlobalCellData("MASS_FRACTION");
                
                boost::shared_ptr<pdat::CellData<double> > data_pressure =
                    d_flow_model_tmp->getGlobalCellData("PRESSURE");
                
                boost::shared_ptr<pdat::CellData<double> > data_temperature =
                    d_flow_model_tmp->getGlobalCellData("TEMPERATURE");
                
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
                                 * Compute the the shear stress component.
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
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
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
    
    return averaged_tau_ij;
}


/*
 * Compute averaged derivative of shear stress component with only x direction as inhomogeneous direction.
 * Component index:
 * 0: tau_11
 * 1: tau_12
 * 2: tau_13
 * 3: tau_22
 * 4: tau_23
 * 5: tau_33
 */
std::vector<double>
FlowModelStatisticsUtilitiesFourEqnConservative::
getAveragedDerivativeOfShearStressComponentWithInhomogeneousXDirection(
    const int component_idx,
    const int derivative_direction,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context) const
{
    std::vector<double> averaged_derivative_tau_ij;
    
    boost::shared_ptr<FlowModel> d_flow_model_tmp = d_flow_model.lock();
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Get the refinement ratio from the finest level to the coarest level.
     */
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    hier::IntVector ratioFinestLevelToCoarestLevel =
        patch_hierarchy->getRatioToCoarserLevel(num_levels - 1);
    for (int li = num_levels - 2; li > 0 ; li--)
    {
        ratioFinestLevelToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(li);
    }
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    boost::shared_ptr<ExtendedFlattenedHierarchy> flattened_hierarchy(
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
    const hier::IntVector finest_level_dims = physical_domain_dims*ratioFinestLevelToCoarestLevel;
    
    /*
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        double* tau_ij_der_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_derivative_tau_ij.resize(finest_level_dim_0);
        double* tau_ij_der_avg_global = averaged_derivative_tau_ij.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            tau_ij_der_avg_local[i] = 0.0;
            tau_ij_der_avg_global[i] = 0.0;
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and the quantities in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                // HARD CODE TO BE SIXTH ORDER CENTRAL SCHEME FOR DIFFERENTIATION.
                TBOX_ASSERT(num_ghosts >= hier::IntVector::getOne(d_dim)*3);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", num_ghosts));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTION", num_ghosts));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("PRESSURE", num_ghosts));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("TEMPERATURE", num_ghosts));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_velocity =
                    d_flow_model_tmp->getGlobalCellData("VELOCITY");
                
                boost::shared_ptr<pdat::CellData<double> > data_mass_fraction =
                    d_flow_model_tmp->getGlobalCellData("MASS_FRACTION");
                
                boost::shared_ptr<pdat::CellData<double> > data_pressure =
                    d_flow_model_tmp->getGlobalCellData("PRESSURE");
                
                boost::shared_ptr<pdat::CellData<double> > data_temperature =
                    d_flow_model_tmp->getGlobalCellData("TEMPERATURE");
                
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
                        
                        if ((component_idx == 0) && (derivative_direction == 0))
                        {
                            double mu[7];
                            double mu_v[7];
                            
                            std::vector<const double*> Y_ptr;
                            Y_ptr.resize(d_num_species);
                            for (int count = -3; count <=3; count++)
                            {
                                const int idx_mass_fraction = relative_idx_lo_0 + (i + count) + num_ghosts_0_mass_fraction;
                                const int idx_pressure = relative_idx_lo_0 + (i + count) + num_ghosts_0_pressure;
                                const int idx_temperature = relative_idx_lo_0 + (i + count) + num_ghosts_0_temperature;
                                
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_ptr[si] = &Y[si][idx_mass_fraction];
                                }
                                mu[count + 3] = d_equation_of_shear_viscosity_mixing_rules->
                                    getShearViscosity(
                                        &p[idx_pressure],
                                        &T[idx_temperature],
                                        Y_ptr);
                                mu_v[count + 3] = d_equation_of_bulk_viscosity_mixing_rules->
                                    getBulkViscosity(
                                        &p[idx_pressure],
                                        &T[idx_temperature],
                                        Y_ptr);
                            }
                            
                            const double dmudx = (double(1)/double(60)*(mu[6] - mu[0])
                                - double(3)/double(20)*(mu[5] - mu[1])
                                + double(3)/double(4)*(mu[4] - mu[2]))/dx[0];
                            
                            const double dmu_vdx = (double(1)/double(60)*(mu_v[6] - mu_v[0])
                                - double(3)/double(20)*(mu_v[5] - mu_v[1])
                                + double(3)/double(4)*(mu_v[4] - mu_v[2]))/dx[0];
                            
                            const int idx_vel       = relative_idx_lo_0 + (i - 1) + num_ghosts_0_velocity;
                            
                            const int idx_vel_x_LLL = relative_idx_lo_0 + (i - 3) + num_ghosts_0_velocity;
                            const int idx_vel_x_LL  = relative_idx_lo_0 + (i - 2) + num_ghosts_0_velocity;
                            const int idx_vel_x_L   = relative_idx_lo_0 + (i - 1) + num_ghosts_0_velocity;
                            const int idx_vel_x_R   = relative_idx_lo_0 + (i + 1) + num_ghosts_0_velocity;
                            const int idx_vel_x_RR  = relative_idx_lo_0 + (i + 2) + num_ghosts_0_velocity;
                            const int idx_vel_x_RRR = relative_idx_lo_0 + (i + 3) + num_ghosts_0_velocity;
                            
                            const double dudx = (double(1)/double(60)*(u[idx_vel_x_RRR] - u[idx_vel_x_LLL])
                                - double(3)/double(20)*(u[idx_vel_x_RR] - u[idx_vel_x_LL])
                                + double(3)/double(4)*(u[idx_vel_x_R] - u[idx_vel_x_L]))/dx[0];
                            
                            const double d2udx2 = (double(1)/double(90)*(u[idx_vel_x_RRR] + u[idx_vel_x_LLL])
                                - double(3)/double(20)*(u[idx_vel_x_RR] + u[idx_vel_x_LL])
                                + double(3)/double(2)*(u[idx_vel_x_R] + u[idx_vel_x_L])
                                - double(49)/double(18)*u[idx_vel])/(dx[0]*dx[0]);
                            
                            value_to_add = ((double(4)/double(3)*dmudx + dmu_vdx)*dudx 
                                + (double(4)/double(3)*mu[3] + mu_v[3])*d2udx2) / ((double) n_overlapped);
                        }
                        else
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "Cannot compute derivative of shear stress component for one-dimensional problem!\n"
                                << "component_idx = " << component_idx << " and derivative_direction = "
                                << derivative_direction << " given!\n"
                                << std::endl);
                        }
                        
                        /*
                         * Add the data.
                         */
                        
                        for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                        {
                            const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                            
                            tau_ij_der_avg_local[idx_fine] += value_to_add;
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
         * Reduction to get the global average of derivative of shear stress component.
         */
        
        mpi.Allreduce(
            tau_ij_der_avg_local,
            tau_ij_der_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(tau_ij_der_avg_local);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        
        double* tau_ij_der_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_derivative_tau_ij.resize(finest_level_dim_0);
        double* tau_ij_der_avg_global = averaged_derivative_tau_ij.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            tau_ij_der_avg_local[i] = 0.0;
            tau_ij_der_avg_global[i] = 0.0;
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and the quantities in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                // HARD CODE TO BE SIXTH ORDER CENTRAL SCHEME FOR DIFFERENTIATION.
                TBOX_ASSERT(num_ghosts >= hier::IntVector::getOne(d_dim)*3);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", num_ghosts));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTION", num_ghosts));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("PRESSURE", num_ghosts));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("TEMPERATURE", num_ghosts));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_velocity =
                    d_flow_model_tmp->getGlobalCellData("VELOCITY");
                
                boost::shared_ptr<pdat::CellData<double> > data_mass_fraction =
                    d_flow_model_tmp->getGlobalCellData("MASS_FRACTION");
                
                boost::shared_ptr<pdat::CellData<double> > data_pressure =
                    d_flow_model_tmp->getGlobalCellData("PRESSURE");
                
                boost::shared_ptr<pdat::CellData<double> > data_temperature =
                    d_flow_model_tmp->getGlobalCellData("TEMPERATURE");
                
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
                             * Compute the the shear stress component.
                             */
                            
                            double value_to_add = double(0);
                            
                            if ((component_idx == 0) && (derivative_direction == 0))
                            {
                                double mu[7];
                                double mu_v[7];
                                
                                std::vector<const double*> Y_ptr;
                                Y_ptr.resize(d_num_species);
                                for (int count = -3; count <=3; count++)
                                {
                                    const int idx_mass_fraction = (relative_idx_lo_0 + (i + count) + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                    
                                    const int idx_pressure = (relative_idx_lo_0 + (i + count) + num_ghosts_0_pressure) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure;
                                    
                                    const int idx_temperature = (relative_idx_lo_0 + (i + count) + num_ghosts_0_temperature) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_temperature)*ghostcell_dim_0_temperature;
                                    
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        Y_ptr[si] = &Y[si][idx_mass_fraction];
                                    }
                                    mu[count + 3] = d_equation_of_shear_viscosity_mixing_rules->
                                        getShearViscosity(
                                            &p[idx_pressure],
                                            &T[idx_temperature],
                                            Y_ptr);
                                    mu_v[count + 3] = d_equation_of_bulk_viscosity_mixing_rules->
                                        getBulkViscosity(
                                            &p[idx_pressure],
                                            &T[idx_temperature],
                                            Y_ptr);
                                }
                                
                                const double dmudx = (double(1)/double(60)*(mu[6] - mu[0])
                                    - double(3)/double(20)*(mu[5] - mu[1])
                                    + double(3)/double(4)*(mu[4] - mu[2]))/dx[0];
                                
                                const double dmu_vdx = (double(1)/double(60)*(mu_v[6] - mu_v[0])
                                    - double(3)/double(20)*(mu_v[5] - mu_v[1])
                                    + double(3)/double(4)*(mu_v[4] - mu_v[2]))/dx[0];
                                
                                const int idx_vel       = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                
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
                                
                                const double dudx = (double(1)/double(60)*(u[idx_vel_x_RRR] - u[idx_vel_x_LLL])
                                    - double(3)/double(20)*(u[idx_vel_x_RR] - u[idx_vel_x_LL])
                                    + double(3)/double(4)*(u[idx_vel_x_R] - u[idx_vel_x_L]))/dx[0];
                                
                                const double d2udx2 = (double(1)/double(90)*(u[idx_vel_x_RRR] + u[idx_vel_x_LLL])
                                    - double(3)/double(20)*(u[idx_vel_x_RR] + u[idx_vel_x_LL])
                                    + double(3)/double(2)*(u[idx_vel_x_R] + u[idx_vel_x_L])
                                    - double(49)/double(18)*u[idx_vel])/(dx[0]*dx[0]);
                                
                                double visc_dvdy[7];
                                
                                for (int count = -3; count <=3; count++)
                                {
                                    const int idx_vel_y_BBB = (relative_idx_lo_0 + (i + count) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j - 3) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                    
                                    const int idx_vel_y_BB  = (relative_idx_lo_0 + (i + count) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j - 2) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                    
                                    const int idx_vel_y_B   = (relative_idx_lo_0 + (i + count) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j - 1) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                    
                                    const int idx_vel_y_T   = (relative_idx_lo_0 + (i + count) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + 1) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                    
                                    const int idx_vel_y_TT  = (relative_idx_lo_0 + (i + count) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + 2) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                    
                                    const int idx_vel_y_TTT = (relative_idx_lo_0 + (i + count) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + 3) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                    
                                    visc_dvdy[count + 3] = (double(2)/double(3)*mu[count + 3] - mu_v[count + 3])*
                                        (double(1)/double(60)*(v[idx_vel_y_TTT] - v[idx_vel_y_BBB])
                                        - double(3)/double(20)*(v[idx_vel_y_TT] - v[idx_vel_y_BB])
                                        + double(3)/double(4)*(v[idx_vel_y_T] - v[idx_vel_y_B]))/dx[1];
                                }
                                
                                const double dvisc_dvdydx = (double(1)/double(60)*(visc_dvdy[6] - visc_dvdy[0])
                                    - double(3)/double(20)*(visc_dvdy[5] - visc_dvdy[1])
                                    + double(3)/double(4)*(visc_dvdy[4] - visc_dvdy[2]))/dx[0];
                                
                                value_to_add = (
                                    (double(4)/double(3)*dmudx + dmu_vdx)*dudx
                                    + (double(4)/double(3)*mu[3] + mu_v[3])*d2udx2
                                    - dvisc_dvdydx
                                    )*weight/((double) n_overlapped);
                            }
                            else if ((component_idx == 0) && (derivative_direction == 1))
                            {
                                double mu[7];
                                double mu_v[7];
                                
                                std::vector<const double*> Y_ptr;
                                Y_ptr.resize(d_num_species);
                                for (int count = -3; count <=3; count++)
                                {
                                    const int idx_mass_fraction = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + (j + count) + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                    
                                    const int idx_pressure = (relative_idx_lo_0 + i + num_ghosts_0_pressure) +
                                        (relative_idx_lo_1 + (j + count) + num_ghosts_1_pressure)*ghostcell_dim_0_pressure;
                                    
                                    const int idx_temperature = (relative_idx_lo_0 + i + num_ghosts_0_temperature) +
                                        (relative_idx_lo_1 + (j + count) + num_ghosts_1_temperature)*ghostcell_dim_0_temperature;
                                    
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        Y_ptr[si] = &Y[si][idx_mass_fraction];
                                    }
                                    mu[count + 3] = d_equation_of_shear_viscosity_mixing_rules->
                                        getShearViscosity(
                                            &p[idx_pressure],
                                            &T[idx_temperature],
                                            Y_ptr);
                                    mu_v[count + 3] = d_equation_of_bulk_viscosity_mixing_rules->
                                        getBulkViscosity(
                                            &p[idx_pressure],
                                            &T[idx_temperature],
                                            Y_ptr);
                                }
                                
                                const double dmudy = (double(1)/double(60)*(mu[6] - mu[0])
                                    - double(3)/double(20)*(mu[5] - mu[1])
                                    + double(3)/double(4)*(mu[4] - mu[2]))/dx[1];
                                
                                const double dmu_vdy = (double(1)/double(60)*(mu_v[6] - mu_v[0])
                                    - double(3)/double(20)*(mu_v[5] - mu_v[1])
                                    + double(3)/double(4)*(mu_v[4] - mu_v[2]))/dx[1];
                                
                                const int idx_vel       = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
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
                                
                                const double dvdy = (double(1)/double(60)*(v[idx_vel_y_TTT] - v[idx_vel_y_BBB])
                                    - double(3)/double(20)*(v[idx_vel_y_TT] - v[idx_vel_y_BB])
                                    + double(3)/double(4)*(v[idx_vel_y_T] - v[idx_vel_y_B]))/dx[1];
                                
                                const double d2vdy2 = (double(1)/double(90)*(v[idx_vel_y_TTT] + v[idx_vel_y_BBB])
                                    - double(3)/double(20)*(v[idx_vel_y_TT] + v[idx_vel_y_BB])
                                    + double(3)/double(2)*(v[idx_vel_y_T] + v[idx_vel_y_B])
                                    - double(49)/double(18)*v[idx_vel])/(dx[1]*dx[1]);
                                
                                double visc_dudx[7];
                                
                                for (int count = -3; count <=3; count++)
                                {
                                    const int idx_vel_x_LLL = (relative_idx_lo_0 + (i - 3) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + count) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                    
                                    const int idx_vel_x_LL  = (relative_idx_lo_0 + (i - 2) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + count) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                    
                                    const int idx_vel_x_L   = (relative_idx_lo_0 + (i - 1) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + count) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                    
                                    const int idx_vel_x_R   = (relative_idx_lo_0 + (i + 1) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + count) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                    
                                    const int idx_vel_x_RR  = (relative_idx_lo_0 + (i + 2) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + count) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                    
                                    const int idx_vel_x_RRR = (relative_idx_lo_0 + (i + 3) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + count) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                    
                                    visc_dudx[count + 3] = (double(4)/double(3)*mu[count + 3] + mu_v[count + 3])*
                                        (double(1)/double(60)*(u[idx_vel_x_RRR] - u[idx_vel_x_LLL])
                                        - double(3)/double(20)*(u[idx_vel_x_RR] - u[idx_vel_x_LL])
                                        + double(3)/double(4)*(u[idx_vel_x_R] - u[idx_vel_x_L]))/dx[0];
                                }
                                
                                const double dvisc_dudxdy = (double(1)/double(60)*(visc_dudx[6] - visc_dudx[0])
                                    - double(3)/double(20)*(visc_dudx[5] - visc_dudx[1])
                                    + double(3)/double(4)*(visc_dudx[4] - visc_dudx[2]))/dx[1];
                                
                                value_to_add = (
                                    dvisc_dudxdy
                                    - (double(2)/double(3)*dmudy - dmu_vdy)*dvdy
                                    - (double(2)/double(3)*mu[3] - mu_v[3])*d2vdy2
                                    )*weight/((double) n_overlapped);
                            }
                            else if ((component_idx == 1) && (derivative_direction == 0))
                            {
                                double mu[7];
                                
                                std::vector<const double*> Y_ptr;
                                Y_ptr.resize(d_num_species);
                                for (int count = -3; count <=3; count++)
                                {
                                    const int idx_mass_fraction = (relative_idx_lo_0 + (i + count) + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                    
                                    const int idx_pressure = (relative_idx_lo_0 + (i + count) + num_ghosts_0_pressure) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure;
                                    
                                    const int idx_temperature = (relative_idx_lo_0 + (i + count) + num_ghosts_0_temperature) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_temperature)*ghostcell_dim_0_temperature;
                                    
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        Y_ptr[si] = &Y[si][idx_mass_fraction];
                                    }
                                    mu[count + 3] = d_equation_of_shear_viscosity_mixing_rules->
                                        getShearViscosity(
                                            &p[idx_pressure],
                                            &T[idx_temperature],
                                            Y_ptr);
                                }
                                
                                const double dmudx = (double(1)/double(60)*(mu[6] - mu[0])
                                    - double(3)/double(20)*(mu[5] - mu[1])
                                    + double(3)/double(4)*(mu[4] - mu[2]))/dx[0];
                                
                                const int idx_vel       = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                
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
                                
                                const double dvdx = (double(1)/double(60)*(v[idx_vel_x_RRR] - v[idx_vel_x_LLL])
                                    - double(3)/double(20)*(v[idx_vel_x_RR] - v[idx_vel_x_LL])
                                    + double(3)/double(4)*(v[idx_vel_x_R] - v[idx_vel_x_L]))/dx[0];
                                
                                const double d2vdx2 = (double(1)/double(90)*(v[idx_vel_x_RRR] + v[idx_vel_x_LLL])
                                    - double(3)/double(20)*(v[idx_vel_x_RR] + v[idx_vel_x_LL])
                                    + double(3)/double(2)*(v[idx_vel_x_R] + v[idx_vel_x_L])
                                    - double(49)/double(18)*v[idx_vel])/(dx[0]*dx[0]);
                                
                                double visc_dudy[7];
                                
                                for (int count = -3; count <=3; count++)
                                {
                                    const int idx_vel_y_BBB = (relative_idx_lo_0 + (i + count) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j - 3) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                    
                                    const int idx_vel_y_BB  = (relative_idx_lo_0 + (i + count) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j - 2) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                    
                                    const int idx_vel_y_B   = (relative_idx_lo_0 + (i + count) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j - 1) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                    
                                    const int idx_vel_y_T   = (relative_idx_lo_0 + (i + count) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + 1) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                    
                                    const int idx_vel_y_TT  = (relative_idx_lo_0 + (i + count) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + 2) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                    
                                    const int idx_vel_y_TTT = (relative_idx_lo_0 + (i + count) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + 3) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                    
                                    visc_dudy[count + 3] = mu[count + 3]*
                                        (double(1)/double(60)*(u[idx_vel_y_TTT] - u[idx_vel_y_BBB])
                                        - double(3)/double(20)*(u[idx_vel_y_TT] - u[idx_vel_y_BB])
                                        + double(3)/double(4)*(u[idx_vel_y_T] - u[idx_vel_y_B]))/dx[1];
                                }
                                
                                const double dvisc_dudydx = (double(1)/double(60)*(visc_dudy[6] - visc_dudy[0])
                                    - double(3)/double(20)*(visc_dudy[5] - visc_dudy[1])
                                    + double(3)/double(4)*(visc_dudy[4] - visc_dudy[2]))/dx[0];
                                
                                value_to_add = (
                                    dvisc_dudydx
                                    + dmudx*dvdx
                                    + mu[3]*d2vdx2
                                    )*weight/((double) n_overlapped);
                            }
                            else if ((component_idx == 1) && (derivative_direction == 1))
                            {
                                double mu[7];
                                
                                std::vector<const double*> Y_ptr;
                                Y_ptr.resize(d_num_species);
                                for (int count = -3; count <=3; count++)
                                {
                                    const int idx_mass_fraction = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + (j + count) + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                    
                                    const int idx_pressure = (relative_idx_lo_0 + i + num_ghosts_0_pressure) +
                                        (relative_idx_lo_1 + (j + count) + num_ghosts_1_pressure)*ghostcell_dim_0_pressure;
                                    
                                    const int idx_temperature = (relative_idx_lo_0 + i + num_ghosts_0_temperature) +
                                        (relative_idx_lo_1 + (j + count) + num_ghosts_1_temperature)*ghostcell_dim_0_temperature;
                                    
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        Y_ptr[si] = &Y[si][idx_mass_fraction];
                                    }
                                    mu[count + 3] = d_equation_of_shear_viscosity_mixing_rules->
                                        getShearViscosity(
                                            &p[idx_pressure],
                                            &T[idx_temperature],
                                            Y_ptr);
                                }
                                
                                const double dmudy = (double(1)/double(60)*(mu[6] - mu[0])
                                    - double(3)/double(20)*(mu[5] - mu[1])
                                    + double(3)/double(4)*(mu[4] - mu[2]))/dx[1];
                                
                                const int idx_vel       = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
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
                                
                                const double dudy = (double(1)/double(60)*(u[idx_vel_y_TTT] - u[idx_vel_y_BBB])
                                    - double(3)/double(20)*(u[idx_vel_y_TT] - u[idx_vel_y_BB])
                                    + double(3)/double(4)*(u[idx_vel_y_T] - u[idx_vel_y_B]))/dx[1];
                                
                                const double d2udy2 = (double(1)/double(90)*(u[idx_vel_y_TTT] + u[idx_vel_y_BBB])
                                    - double(3)/double(20)*(u[idx_vel_y_TT] + u[idx_vel_y_BB])
                                    + double(3)/double(2)*(u[idx_vel_y_T] + u[idx_vel_y_B])
                                    - double(49)/double(18)*u[idx_vel])/(dx[1]*dx[1]);
                                
                                double visc_dvdx[7];
                                
                                for (int count = -3; count <=3; count++)
                                {
                                    const int idx_vel_x_LLL = (relative_idx_lo_0 + (i - 3) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + count) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                    
                                    const int idx_vel_x_LL  = (relative_idx_lo_0 + (i - 2) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + count) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                    
                                    const int idx_vel_x_L   = (relative_idx_lo_0 + (i - 1) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + count) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                    
                                    const int idx_vel_x_R   = (relative_idx_lo_0 + (i + 1) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + count) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                    
                                    const int idx_vel_x_RR  = (relative_idx_lo_0 + (i + 2) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + count) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                    
                                    const int idx_vel_x_RRR = (relative_idx_lo_0 + (i + 3) + num_ghosts_0_velocity) +
                                        (relative_idx_lo_1 + (j + count) + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                                    
                                    visc_dvdx[count + 3] = mu[count + 3]*
                                        (double(1)/double(60)*(v[idx_vel_x_RRR] - v[idx_vel_x_LLL])
                                        - double(3)/double(20)*(v[idx_vel_x_RR] - v[idx_vel_x_LL])
                                        + double(3)/double(4)*(v[idx_vel_x_R] - v[idx_vel_x_L]))/dx[0];
                                }
                                
                                const double dvisc_dvdxdy = (double(1)/double(60)*(visc_dvdx[6] - visc_dvdx[0])
                                    - double(3)/double(20)*(visc_dvdx[5] - visc_dvdx[1])
                                    + double(3)/double(4)*(visc_dvdx[4] - visc_dvdx[2]))/dx[1];
                                
                                value_to_add = (
                                    dvisc_dvdxdy
                                    + dmudy*dudy
                                    + mu[3]*d2udy2
                                    )*weight/((double) n_overlapped);
                            }
                            else
                            {
                                TBOX_ERROR(d_object_name
                                    << ": "
                                    << "Cannot compute derivative of shear stress component for two-dimensional problem!\n"
                                    << "component_idx = " << component_idx << " and derivative_direction = "
                                    << derivative_direction << " given!\n"
                                    << std::endl);
                            }
                            
                            /*
                             * Add the data.
                             */
                            
                            for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                            {
                                const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                
                                tau_ij_der_avg_local[idx_fine] += value_to_add;
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
         * Reduction to get the global average of derivative of shear stress component.
         */
        
        mpi.Allreduce(
            tau_ij_der_avg_local,
            tau_ij_der_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(tau_ij_der_avg_local);
    }
    else if (d_dim == tbox::Dimension(3))
    {
       const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        const double L_z = x_hi[2] - x_lo[2];
        
        double* tau_ij_der_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_derivative_tau_ij.resize(finest_level_dim_0);
        double* tau_ij_der_avg_global = averaged_derivative_tau_ij.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            tau_ij_der_avg_local[i] = 0.0;
            tau_ij_der_avg_global[i] = 0.0;
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and the quantities in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                // HARD CODE TO BE SIXTH ORDER CENTRAL SCHEME FOR DIFFERENTIATION.
                TBOX_ASSERT(num_ghosts >= hier::IntVector::getOne(d_dim)*3);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", num_ghosts));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTION", num_ghosts));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("PRESSURE", num_ghosts));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("TEMPERATURE", num_ghosts));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_velocity =
                    d_flow_model_tmp->getGlobalCellData("VELOCITY");
                
                boost::shared_ptr<pdat::CellData<double> > data_mass_fraction =
                    d_flow_model_tmp->getGlobalCellData("MASS_FRACTION");
                
                boost::shared_ptr<pdat::CellData<double> > data_pressure =
                    d_flow_model_tmp->getGlobalCellData("PRESSURE");
                
                boost::shared_ptr<pdat::CellData<double> > data_temperature =
                    d_flow_model_tmp->getGlobalCellData("TEMPERATURE");
                
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
                                 * Compute the the shear stress component.
                                 */
                                
                                double value_to_add = double(0);
                                
                                if ((component_idx == 0) && (derivative_direction == 0))
                                {
                                }
                                else if ((component_idx == 0) && (derivative_direction == 1))
                                {
                                }
                                else if ((component_idx == 0) && (derivative_direction == 2))
                                {
                                }
                                else if ((component_idx == 1) && (derivative_direction == 0))
                                {
                                }
                                else if ((component_idx == 1) && (derivative_direction == 1))
                                {
                                }
                                else if ((component_idx == 1) && (derivative_direction == 2))
                                {
                                }
                                else if ((component_idx == 2) && (derivative_direction == 0))
                                {
                                }
                                else if ((component_idx == 2) && (derivative_direction == 1))
                                {
                                }
                                else if ((component_idx == 2) && (derivative_direction == 2))
                                {
                                }
                                else
                                {
                                    TBOX_ERROR(d_object_name
                                        << ": "
                                        << "Cannot compute derivative of shear stress component for three-dimensional problem!\n"
                                        << "component_idx = " << component_idx << " and derivative_direction = "
                                        << derivative_direction << " given!\n"
                                        << std::endl);
                                }
                                
                                /*
                                 * Add the data.
                                 */
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    tau_ij_der_avg_local[idx_fine] += value_to_add;
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
         * Reduction to get the global average of derivative of shear stress component.
         */
        
        mpi.Allreduce(
            tau_ij_der_avg_local,
            tau_ij_der_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(tau_ij_der_avg_local);
    }
    
    return averaged_derivative_tau_ij;
}


/*
 * Compute averaged value (on product of variables) with only x direction as inhomogeneous direction.
 */
std::vector<double>
FlowModelStatisticsUtilitiesFourEqnConservative::getAveragedQuantityWithInhomogeneousXDirection(
    const std::vector<std::string>& quantity_names,
    const std::vector<int>& component_indices,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context) const
{
    int num_quantities = static_cast<int>(quantity_names.size());
    
    TBOX_ASSERT(static_cast<int>(component_indices.size()) == num_quantities);
    
    std::vector<double> averaged_quantity;
    
    boost::shared_ptr<FlowModel> d_flow_model_tmp = d_flow_model.lock();
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Get the refinement ratio from the finest level to the coarest level.
     */
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    hier::IntVector ratioFinestLevelToCoarestLevel =
        patch_hierarchy->getRatioToCoarserLevel(num_levels - 1);
    for (int li = num_levels - 2; li > 0 ; li--)
    {
        ratioFinestLevelToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(li);
    }
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    boost::shared_ptr<ExtendedFlattenedHierarchy> flattened_hierarchy(
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
    const hier::IntVector finest_level_dims = physical_domain_dims*ratioFinestLevelToCoarestLevel;
    
    /*
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        double* avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        averaged_quantity.resize(finest_level_dim_0);
        double* avg_global = averaged_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            avg_local[i] = 0.0;
            avg_global[i] = 0.0;
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                /*
                 * Register the patch and data in the flow model and compute the corresponding
                 * correlation.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
                }
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<boost::shared_ptr<pdat::CellData<double> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = d_flow_model_tmp->getGlobalCellData(quantity_names[qi]);
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
                
                std::vector<int> num_ghosts_0_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                }
                
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
                         * Compute the linear indices and the data to add.
                         */
                        
                        for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                        {
                            const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                            
                            const int idx_q0 = relative_idx_lo_0 + i + num_ghosts_0_u_qi[0];
                            
                            double avg = u_qi[0][idx_q0];
                            
                            for (int qi = 1; qi < num_quantities; qi++)
                            {
                                const int idx_qi = relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi];
                                
                                avg *= u_qi[qi][idx_qi];
                            }
                            
                            avg_local[idx_fine] += (avg/((double) n_overlapped));
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
            avg_local[i] = 0.0;
            avg_global[i] = 0.0;
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and data in the flow model and compute the corresponding
                 * correlation.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
                }
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<boost::shared_ptr<pdat::CellData<double> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = d_flow_model_tmp->getGlobalCellData(quantity_names[qi]);
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
                             * Compute the linear indices and the data to add.
                             */
                            
                            for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                            {
                                const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                
                                const int idx_q0 = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[0]) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0];
                                
                                double avg = u_qi[0][idx_q0];
                                
                                for (int qi = 1; qi < num_quantities; qi++)
                                {
                                    const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    avg *= u_qi[qi][idx_qi];
                                }
                                
                                avg_local[idx_fine] += (avg*weight/((double) n_overlapped));
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
            avg_local[i] = 0.0;
            avg_global[i] = 0.0;
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and data in the flow model and compute the corresponding
                 * correlation.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
                }
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<boost::shared_ptr<pdat::CellData<double> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = d_flow_model_tmp->getGlobalCellData(quantity_names[qi]);
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
                
                const double weight = dx[1]*dx[2]/(L_y*L_z);
                
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
                                 * Compute the linear index and the data to add.
                                 */
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    const int idx_q0 = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[0]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0] +
                                        (relative_idx_lo_2 + k + num_ghosts_2_u_qi[0])*ghostcell_dim_0_u_qi[0]*
                                            ghostcell_dim_1_u_qi[0];
                                    
                                    double avg = u_qi[0][idx_q0];
                                    
                                    for (int qi = 1; qi < num_quantities; qi++)
                                    {
                                        const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        avg *= u_qi[qi][idx_qi];
                                    }
                                    
                                    avg_local[idx_fine] += (avg*weight/((double) n_overlapped));
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
 * Compute correlation with only x direction as inhomogeneous direction.
 */
std::vector<double>
FlowModelStatisticsUtilitiesFourEqnConservative::getQuantityCorrelationWithInhomogeneousXDirection(
    const std::vector<std::string>& quantity_names,
    const std::vector<int>& component_indices,
    const std::vector<std::vector<double> >& averaged_quantities,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context) const
{
    int num_quantities = static_cast<int>(quantity_names.size());
    
    TBOX_ASSERT(static_cast<int>(component_indices.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(averaged_quantities.size()) == num_quantities);
    
    std::vector<double> correlation;
    
    boost::shared_ptr<FlowModel> d_flow_model_tmp = d_flow_model.lock();
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Get the refinement ratio from the finest level to the coarest level.
     */
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    hier::IntVector ratioFinestLevelToCoarestLevel =
        patch_hierarchy->getRatioToCoarserLevel(num_levels - 1);
    for (int li = num_levels - 2; li > 0 ; li--)
    {
        ratioFinestLevelToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(li);
    }
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    boost::shared_ptr<ExtendedFlattenedHierarchy> flattened_hierarchy(
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
    const hier::IntVector finest_level_dims = physical_domain_dims*ratioFinestLevelToCoarestLevel;
    
    /*
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << "FlowModelStatisticsUtilitiesFourEqnConservative::\n"
            << "getQuantityCorrelationWithInhomogeneousXDirection:\n"
            << "Not implemented for 1D problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        
        std::vector<const double*> u_qi_avg_global;
        u_qi_avg_global.reserve(num_quantities);
        for (int qi = 0; qi < num_quantities; qi++)
        {
            u_qi_avg_global.push_back(averaged_quantities[qi].data());
        }
        
        double* corr_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        correlation.resize(finest_level_dim_0);
        double* corr_global = correlation.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            corr_local[i] = 0.0;
            corr_global[i] = 0.0;
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                // Get the dimensions of box that covers the interior of patch.
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and data in the flow model and compute the corresponding
                 * correlation.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
                }
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<boost::shared_ptr<pdat::CellData<double> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = d_flow_model_tmp->getGlobalCellData(quantity_names[qi]);
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
                             * Compute the linear indices and the data to add.
                             */
                            
                            for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                            {
                                const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                
                                const int idx_q0 = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[0]) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0];
                                
                                double corr = u_qi[0][idx_q0] - u_qi_avg_global[0][idx_fine];
                                
                                for (int qi = 1; qi < num_quantities; qi++)
                                {
                                    const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    corr *= (u_qi[qi][idx_qi] - u_qi_avg_global[qi][idx_fine]);
                                }
                                
                                corr_local[idx_fine] += (corr*weight/((double) n_overlapped));
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
        
        mpi.Allreduce(
            corr_local,
            corr_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(corr_local);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        const double L_z = x_hi[2] - x_lo[2];
        
        std::vector<const double*> u_qi_avg_global;
        u_qi_avg_global.reserve(num_quantities);
        for (int qi = 0; qi < num_quantities; qi++)
        {
            u_qi_avg_global.push_back(averaged_quantities[qi].data());
        }
        
        double* corr_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        correlation.resize(finest_level_dim_0);
        double* corr_global = correlation.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            corr_local[i] = 0.0;
            corr_global[i] = 0.0;
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                // Get the dimensions of box that covers the interior of patch.
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and data in the flow model and compute the corresponding
                 * correlation.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
                }
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<boost::shared_ptr<pdat::CellData<double> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = d_flow_model_tmp->getGlobalCellData(quantity_names[qi]);
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
                
                const double weight = dx[1]*dx[2]/(L_y*L_z);
                
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
                                 * Compute the linear index and the data to add.
                                 */
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    const int idx_q0 = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[0]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0] +
                                        (relative_idx_lo_2 + k + num_ghosts_2_u_qi[0])*ghostcell_dim_0_u_qi[0]*
                                            ghostcell_dim_1_u_qi[0];
                                    
                                    double corr = u_qi[0][idx_q0] - u_qi_avg_global[0][idx_fine];
                                    
                                    for (int qi = 1; qi < num_quantities; qi++)
                                    {
                                        const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        corr *= (u_qi[qi][idx_qi] - u_qi_avg_global[qi][idx_fine]);
                                    }
                                    
                                    corr_local[idx_fine] += (corr*weight/((double) n_overlapped));
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
        
        mpi.Allreduce(
            corr_local,
            corr_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(corr_local);
    }
    
    return correlation;
}


/*
 * Compute correlation with only x direction as inhomogeneous direction.
 */
std::vector<double>
FlowModelStatisticsUtilitiesFourEqnConservative::getQuantityCorrelationWithInhomogeneousXDirection(
    const std::vector<std::string>& quantity_names,
    const std::vector<int>& component_indices,
    const std::vector<bool>& use_reciprocal,
    const std::vector<std::vector<double> >& averaged_quantities,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context) const
{
    int num_quantities = static_cast<int>(quantity_names.size());
    
    TBOX_ASSERT(static_cast<int>(component_indices.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(use_reciprocal.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(averaged_quantities.size()) == num_quantities);
    
    std::vector<double> correlation;
    
    boost::shared_ptr<FlowModel> d_flow_model_tmp = d_flow_model.lock();
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Get the refinement ratio from the finest level to the coarest level.
     */
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    hier::IntVector ratioFinestLevelToCoarestLevel =
        patch_hierarchy->getRatioToCoarserLevel(num_levels - 1);
    for (int li = num_levels - 2; li > 0 ; li--)
    {
        ratioFinestLevelToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(li);
    }
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    boost::shared_ptr<ExtendedFlattenedHierarchy> flattened_hierarchy(
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
    const hier::IntVector finest_level_dims = physical_domain_dims*ratioFinestLevelToCoarestLevel;
    
    /*
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << "FlowModelStatisticsUtilitiesFourEqnConservative::\n"
            << "getQuantityCorrelationWithInhomogeneousXDirection:\n"
            << "Not implemented for 1D problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        
        std::vector<const double*> u_qi_avg_global;
        u_qi_avg_global.reserve(num_quantities);
        for (int qi = 0; qi < num_quantities; qi++)
        {
            u_qi_avg_global.push_back(averaged_quantities[qi].data());
        }
        
        double* corr_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        correlation.resize(finest_level_dim_0);
        double* corr_global = correlation.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            corr_local[i] = 0.0;
            corr_global[i] = 0.0;
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                // Get the dimensions of box that covers the interior of patch.
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and data in the flow model and compute the corresponding
                 * correlation.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
                }
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<boost::shared_ptr<pdat::CellData<double> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = d_flow_model_tmp->getGlobalCellData(quantity_names[qi]);
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
                             * Compute the linear indices and the data to add.
                             */
                            
                            for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                            {
                                const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                
                                const int idx_q0 = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[0]) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0];
                                
                                double corr = double(0);
                                if (use_reciprocal[0])
                                {
                                    corr = double(1)/(u_qi[0][idx_q0]) - u_qi_avg_global[0][idx_fine];
                                }
                                else
                                {
                                    corr = u_qi[0][idx_q0] - u_qi_avg_global[0][idx_fine];
                                }
                                
                                for (int qi = 1; qi < num_quantities; qi++)
                                {
                                    const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    if (use_reciprocal[qi])
                                    {
                                        corr *= (double(1)/(u_qi[qi][idx_qi]) - u_qi_avg_global[qi][idx_fine]);
                                    }
                                    else
                                    {
                                        corr *= (u_qi[qi][idx_qi] - u_qi_avg_global[qi][idx_fine]);
                                    }
                                }
                                
                                corr_local[idx_fine] += (corr*weight/((double) n_overlapped));
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
        
        mpi.Allreduce(
            corr_local,
            corr_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(corr_local);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        const double L_z = x_hi[2] - x_lo[2];
        
        std::vector<const double*> u_qi_avg_global;
        u_qi_avg_global.reserve(num_quantities);
        for (int qi = 0; qi < num_quantities; qi++)
        {
            u_qi_avg_global.push_back(averaged_quantities[qi].data());
        }
        
        double* corr_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        correlation.resize(finest_level_dim_0);
        double* corr_global = correlation.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            corr_local[i] = 0.0;
            corr_global[i] = 0.0;
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                // Get the dimensions of box that covers the interior of patch.
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and data in the flow model and compute the corresponding
                 * correlation.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
                }
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<boost::shared_ptr<pdat::CellData<double> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = d_flow_model_tmp->getGlobalCellData(quantity_names[qi]);
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
                
                const double weight = dx[1]*dx[2]/(L_y*L_z);
                
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
                                 * Compute the linear index and the data to add.
                                 */
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    const int idx_q0 = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[0]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0] +
                                        (relative_idx_lo_2 + k + num_ghosts_2_u_qi[0])*ghostcell_dim_0_u_qi[0]*
                                            ghostcell_dim_1_u_qi[0];
                                    
                                    double corr = double(0);
                                    if (use_reciprocal[0])
                                    {
                                        corr = double(1)/(u_qi[0][idx_q0]) - u_qi_avg_global[0][idx_fine];
                                    }
                                    else
                                    {
                                        corr = u_qi[0][idx_q0] - u_qi_avg_global[0][idx_fine];
                                    }
                                    
                                    for (int qi = 1; qi < num_quantities; qi++)
                                    {
                                        const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        if (use_reciprocal[qi])
                                        {
                                            corr *= (double(1)/(u_qi[qi][idx_qi]) - u_qi_avg_global[qi][idx_fine]);
                                        }
                                        else
                                        {
                                            corr *= (u_qi[qi][idx_qi] - u_qi_avg_global[qi][idx_fine]);
                                        }
                                    }
                                    
                                    corr_local[idx_fine] += (corr*weight/((double) n_overlapped));
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
        
        mpi.Allreduce(
            corr_local,
            corr_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(corr_local);
    }
    
    return correlation;
}


/*
 * Compute correlation with only x direction as inhomogeneous direction.
 */
std::vector<double>
FlowModelStatisticsUtilitiesFourEqnConservative::getQuantityCorrelationWithInhomogeneousXDirection(
    const std::vector<std::string>& quantity_names,
    const std::vector<int>& component_indices,
    const std::vector<bool>& use_reciprocal,
    const std::vector<int>& derivative_directions,
    const std::vector<std::vector<double> >& averaged_quantities,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context) const
{
    int num_quantities = static_cast<int>(quantity_names.size());
    
    TBOX_ASSERT(static_cast<int>(component_indices.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(use_reciprocal.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(derivative_directions.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(averaged_quantities.size()) == num_quantities);
    
    std::vector<double> correlation;
    
    boost::shared_ptr<FlowModel> d_flow_model_tmp = d_flow_model.lock();
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Get the refinement ratio from the finest level to the coarest level.
     */
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    hier::IntVector ratioFinestLevelToCoarestLevel =
        patch_hierarchy->getRatioToCoarserLevel(num_levels - 1);
    for (int li = num_levels - 2; li > 0 ; li--)
    {
        ratioFinestLevelToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(li);
    }
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    boost::shared_ptr<ExtendedFlattenedHierarchy> flattened_hierarchy(
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
    const hier::IntVector finest_level_dims = physical_domain_dims*ratioFinestLevelToCoarestLevel;
    
    /*
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << "FlowModelStatisticsUtilitiesFourEqnConservative::\n"
            << "getQuantityCorrelationWithInhomogeneousXDirection:\n"
            << "Not implemented for 1D problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        
        std::vector<const double*> u_qi_avg_global;
        u_qi_avg_global.reserve(num_quantities);
        for (int qi = 0; qi < num_quantities; qi++)
        {
            u_qi_avg_global.push_back(averaged_quantities[qi].data());
        }
        
        double* corr_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        correlation.resize(finest_level_dim_0);
        double* corr_global = correlation.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            corr_local[i] = 0.0;
            corr_global[i] = 0.0;
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                // Get the dimensions of box that covers the interior of patch.
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and data in the flow model and compute the corresponding
                 * correlation.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                // HARD CODE TO BE SIXTH ORDER CENTRAL SCHEME FOR DIFFERENTIATION.
                TBOX_ASSERT(num_ghosts >= hier::IntVector::getOne(d_dim)*3);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    std::unordered_map<std::string, hier::IntVector>::const_iterator quantity_it = num_subghosts_of_data.find(quantity_names[qi]);
                    
                    if (quantity_it == num_subghosts_of_data.end())
                    {
                        if (derivative_directions[qi] >= 0)
                        {
                            num_subghosts_of_data.insert(
                                std::pair<std::string, hier::IntVector>(quantity_names[qi], num_ghosts));
                        }
                        else
                        {
                            num_subghosts_of_data.insert(
                                std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
                        }
                    }
                    else
                    {
                        if (derivative_directions[qi] >= 0 && quantity_it->second < num_ghosts)
                        {
                            num_subghosts_of_data.erase(quantity_it);
                            
                            num_subghosts_of_data.insert(
                                std::pair<std::string, hier::IntVector>(quantity_names[qi], num_ghosts));
                        }
                    }
                }
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<boost::shared_ptr<pdat::CellData<double> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = d_flow_model_tmp->getGlobalCellData(quantity_names[qi]);
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
                             * Compute the linear indices and the data to add.
                             */
                            
                            for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                            {
                                const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                
                                double corr = double(0);
                                if (derivative_directions[0] == -1)
                                {
                                    const int idx_q0 = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[0]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0];
                                    
                                    if (use_reciprocal[0])
                                    {
                                        corr = double(1)/(u_qi[0][idx_q0]) - u_qi_avg_global[0][idx_fine];
                                    }
                                    else
                                    {
                                        corr = u_qi[0][idx_q0] - u_qi_avg_global[0][idx_fine];
                                    }
                                }
                                else if (derivative_directions[0] == 0)
                                {
                                    const int idx_LLL = (relative_idx_lo_0 + (i - 3) + num_ghosts_0_u_qi[0]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0];
                                    
                                    const int idx_LL  = (relative_idx_lo_0 + (i - 2) + num_ghosts_0_u_qi[0]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0];
                                    
                                    const int idx_L   = (relative_idx_lo_0 + (i - 1) + num_ghosts_0_u_qi[0]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0];
                                    
                                    const int idx_R   = (relative_idx_lo_0 + (i + 1) + num_ghosts_0_u_qi[0]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0];
                                    
                                    const int idx_RR  = (relative_idx_lo_0 + (i + 2) + num_ghosts_0_u_qi[0]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0];
                                    
                                    const int idx_RRR = (relative_idx_lo_0 + (i + 3) + num_ghosts_0_u_qi[0]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0];
                                    
                                    double dudx = double(0);
                                    if (use_reciprocal[0])
                                    {
                                        dudx = (double(1)/double(60)*(double(1)/u_qi[0][idx_RRR] - double(1)/u_qi[0][idx_LLL])
                                            - double(3)/double(20)*(double(1)/u_qi[0][idx_RR] - double(1)/u_qi[0][idx_LL])
                                            + double(3)/double(4)*(double(1)/u_qi[0][idx_R] - double(1)/u_qi[0][idx_L]))/dx[0];
                                    }
                                    else
                                    {
                                        dudx = (double(1)/double(60)*(u_qi[0][idx_RRR] - u_qi[0][idx_LLL])
                                            - double(3)/double(20)*(u_qi[0][idx_RR] - u_qi[0][idx_LL])
                                            + double(3)/double(4)*(u_qi[0][idx_R] - u_qi[0][idx_L]))/dx[0];
                                    }
                                    corr = dudx - u_qi_avg_global[0][idx_fine];
                                }
                                else if (derivative_directions[0] == 1)
                                {
                                    const int idx_BBB = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[0]) +
                                        (relative_idx_lo_1 + (j - 3) + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0];
                                    
                                    const int idx_BB  = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[0]) +
                                        (relative_idx_lo_1 + (j - 2) + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0];
                                    
                                    const int idx_B   = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[0]) +
                                        (relative_idx_lo_1 + (j - 1) + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0];
                                    
                                    const int idx_T   = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[0]) +
                                        (relative_idx_lo_1 + (j + 1) + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0];
                                    
                                    const int idx_TT  = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[0]) +
                                        (relative_idx_lo_1 + (j + 2) + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0];
                                    
                                    const int idx_TTT = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[0]) +
                                        (relative_idx_lo_1 + (j + 3) + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0];
                                    
                                    double dudy = double(0);
                                    if (use_reciprocal[0])
                                    {
                                        dudy = (double(1)/double(60)*(double(1)/u_qi[0][idx_TTT] - double(1)/u_qi[0][idx_BBB])
                                            - double(3)/double(20)*(double(1)/u_qi[0][idx_TT] - double(1)/u_qi[0][idx_BB])
                                            + double(3)/double(4)*(double(1)/u_qi[0][idx_T] - double(1)/u_qi[0][idx_B]))/dx[1];
                                    }
                                    else
                                    {
                                        dudy = (double(1)/double(60)*(u_qi[0][idx_TTT] - u_qi[0][idx_BBB])
                                            - double(3)/double(20)*(u_qi[0][idx_TT] - u_qi[0][idx_BB])
                                            + double(3)/double(4)*(u_qi[0][idx_T] - u_qi[0][idx_B]))/dx[1];
                                    }
                                    corr = dudy - u_qi_avg_global[0][idx_fine];
                                }
                                else
                                {
                                    TBOX_ERROR(d_object_name
                                        << ": "
                                        << "Cannot take derivative for two-dimensional problem!\n"
                                        << "derivative_direction = " << derivative_directions[0] << " given!\n"
                                        << std::endl);
                                }
                                
                                for (int qi = 1; qi < num_quantities; qi++)
                                {
                                    if (derivative_directions[qi] == -1)
                                    {
                                        const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                        
                                        if (use_reciprocal[qi])
                                        {
                                            corr *= (double(1)/(u_qi[qi][idx_qi]) - u_qi_avg_global[qi][idx_fine]);
                                        }
                                        else
                                        {
                                            corr *= (u_qi[qi][idx_qi] - u_qi_avg_global[qi][idx_fine]);
                                        }
                                    }
                                    else if (derivative_directions[qi] == 0)
                                    {
                                        const int idx_LLL = (relative_idx_lo_0 + (i - 3) + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                        
                                        const int idx_LL  = (relative_idx_lo_0 + (i - 2) + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                        
                                        const int idx_L   = (relative_idx_lo_0 + (i - 1) + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                        
                                        const int idx_R   = (relative_idx_lo_0 + (i + 1) + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                        
                                        const int idx_RR  = (relative_idx_lo_0 + (i + 2) + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                        
                                        const int idx_RRR = (relative_idx_lo_0 + (i + 3) + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                        
                                        double dudx = double(0);
                                        if (use_reciprocal[qi])
                                        {
                                            dudx = (double(1)/double(60)*(double(1)/u_qi[qi][idx_RRR] - double(1)/u_qi[qi][idx_LLL])
                                                - double(3)/double(20)*(double(1)/u_qi[qi][idx_RR] - double(1)/u_qi[qi][idx_LL])
                                                + double(3)/double(4)*(double(1)/u_qi[qi][idx_R] - double(1)/u_qi[qi][idx_L]))/dx[0];
                                        }
                                        else
                                        {
                                            dudx = (double(1)/double(60)*(u_qi[qi][idx_RRR] - u_qi[qi][idx_LLL])
                                                - double(3)/double(20)*(u_qi[qi][idx_RR] - u_qi[qi][idx_LL])
                                                + double(3)/double(4)*(u_qi[qi][idx_R] - u_qi[qi][idx_L]))/dx[0];
                                        }
                                        corr *= (dudx - u_qi_avg_global[qi][idx_fine]);
                                    }
                                    else if (derivative_directions[qi] == 1)
                                    {
                                        const int idx_BBB = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + (j - 3) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                        
                                        const int idx_BB  = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + (j - 2) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                        
                                        const int idx_B   = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + (j - 1) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                        
                                        const int idx_T   = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + (j + 1) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                        
                                        const int idx_TT  = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + (j + 2) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                        
                                        const int idx_TTT = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + (j + 3) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                        
                                        double dudy = double(0);
                                        if (use_reciprocal[qi])
                                        {
                                            dudy = (double(1)/double(60)*(double(1)/u_qi[qi][idx_TTT] - double(1)/u_qi[qi][idx_BBB])
                                                - double(3)/double(20)*(double(1)/u_qi[qi][idx_TT] - double(1)/u_qi[qi][idx_BB])
                                                + double(3)/double(4)*(double(1)/u_qi[qi][idx_T] - double(1)/u_qi[qi][idx_B]))/dx[1];
                                        }
                                        else
                                        {
                                            dudy = (double(1)/double(60)*(u_qi[qi][idx_TTT] - u_qi[qi][idx_BBB])
                                                - double(3)/double(20)*(u_qi[qi][idx_TT] - u_qi[qi][idx_BB])
                                                + double(3)/double(4)*(u_qi[qi][idx_T] - u_qi[qi][idx_B]))/dx[1];
                                        }
                                        corr *= (dudy - u_qi_avg_global[qi][idx_fine]);
                                    }
                                    else
                                    {
                                        TBOX_ERROR(d_object_name
                                            << ": "
                                            << "Cannot take derivative for two-dimensional problem!\n"
                                            << "derivative_direction = " << derivative_directions[qi] << " given!\n"
                                            << std::endl);
                                    }
                                }
                                
                                corr_local[idx_fine] += (corr*weight/((double) n_overlapped));
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
        
        mpi.Allreduce(
            corr_local,
            corr_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(corr_local);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        const double L_z = x_hi[2] - x_lo[2];
        
        std::vector<const double*> u_qi_avg_global;
        u_qi_avg_global.reserve(num_quantities);
        for (int qi = 0; qi < num_quantities; qi++)
        {
            u_qi_avg_global.push_back(averaged_quantities[qi].data());
        }
        
        double* corr_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        correlation.resize(finest_level_dim_0);
        double* corr_global = correlation.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            corr_local[i] = 0.0;
            corr_global[i] = 0.0;
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                // Get the dimensions of box that covers the interior of patch.
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and data in the flow model and compute the corresponding
                 * correlation.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                // HARD CODE TO BE SIXTH ORDER CENTRAL SCHEME FOR DIFFERENTIATION.
                TBOX_ASSERT(num_ghosts >= hier::IntVector::getOne(d_dim)*3);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    std::unordered_map<std::string, hier::IntVector>::const_iterator quantity_it = num_subghosts_of_data.find(quantity_names[qi]);
                    
                    if (quantity_it == num_subghosts_of_data.end())
                    {
                        if (derivative_directions[qi] >= 0)
                        {
                            num_subghosts_of_data.insert(
                                std::pair<std::string, hier::IntVector>(quantity_names[qi], num_ghosts));
                        }
                        else
                        {
                            num_subghosts_of_data.insert(
                                std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
                        }
                    }
                    else
                    {
                        if (derivative_directions[qi] >= 0 && quantity_it->second < num_ghosts)
                        {
                            num_subghosts_of_data.erase(quantity_it);
                            
                            num_subghosts_of_data.insert(
                                std::pair<std::string, hier::IntVector>(quantity_names[qi], num_ghosts));
                        }
                    }
                }
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<boost::shared_ptr<pdat::CellData<double> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = d_flow_model_tmp->getGlobalCellData(quantity_names[qi]);
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
                
                const double weight = dx[1]*dx[2]/(L_y*L_z);
                
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
                                 * Compute the linear index and the data to add.
                                 */
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    double corr = double(0);
                                    if (derivative_directions[0] == -1)
                                    {
                                        const int idx_q0 = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[0]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[0])*ghostcell_dim_0_u_qi[0]*
                                                ghostcell_dim_1_u_qi[0];
                                        
                                        if (use_reciprocal[0])
                                        {
                                            corr = double(1)/(u_qi[0][idx_q0]) - u_qi_avg_global[0][idx_fine];
                                        }
                                        else
                                        {
                                            corr = u_qi[0][idx_q0] - u_qi_avg_global[0][idx_fine];
                                        }
                                    }
                                    else if (derivative_directions[0] == 0)
                                    {
                                        const int idx_LLL = (relative_idx_lo_0 + (i - 3) + num_ghosts_0_u_qi[0]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[0])*ghostcell_dim_0_u_qi[0]*
                                                ghostcell_dim_1_u_qi[0];
                                        
                                        const int idx_LL  = (relative_idx_lo_0 + (i - 2) + num_ghosts_0_u_qi[0]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[0])*ghostcell_dim_0_u_qi[0]*
                                                ghostcell_dim_1_u_qi[0];
                                        
                                        const int idx_L   = (relative_idx_lo_0 + (i - 1) + num_ghosts_0_u_qi[0]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[0])*ghostcell_dim_0_u_qi[0]*
                                                ghostcell_dim_1_u_qi[0];
                                        
                                        const int idx_R   = (relative_idx_lo_0 + (i + 1) + num_ghosts_0_u_qi[0]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[0])*ghostcell_dim_0_u_qi[0]*
                                                ghostcell_dim_1_u_qi[0];
                                        
                                        const int idx_RR  = (relative_idx_lo_0 + (i + 2) + num_ghosts_0_u_qi[0]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[0])*ghostcell_dim_0_u_qi[0]*
                                                ghostcell_dim_1_u_qi[0];
                                        
                                        const int idx_RRR = (relative_idx_lo_0 + (i + 3) + num_ghosts_0_u_qi[0]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[0])*ghostcell_dim_0_u_qi[0]*
                                                ghostcell_dim_1_u_qi[0];
                                        
                                        double dudx = double(0);
                                        if (use_reciprocal[0])
                                        {
                                            dudx = (double(1)/double(60)*(double(1)/u_qi[0][idx_RRR] - double(1)/u_qi[0][idx_LLL])
                                                - double(3)/double(20)*(double(1)/u_qi[0][idx_RR] - double(1)/u_qi[0][idx_LL])
                                                + double(3)/double(4)*(double(1)/u_qi[0][idx_R] - double(1)/u_qi[0][idx_L]))/dx[0];
                                        }
                                        else
                                        {
                                            dudx = (double(1)/double(60)*(u_qi[0][idx_RRR] - u_qi[0][idx_LLL])
                                                - double(3)/double(20)*(u_qi[0][idx_RR] - u_qi[0][idx_LL])
                                                + double(3)/double(4)*(u_qi[0][idx_R] - u_qi[0][idx_L]))/dx[0];
                                        }
                                        corr = dudx - u_qi_avg_global[0][idx_fine];
                                    }
                                    else if (derivative_directions[0] == 1)
                                    {
                                        const int idx_BBB = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[0]) +
                                            (relative_idx_lo_1 + (j - 3) + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[0])*ghostcell_dim_0_u_qi[0]*
                                                ghostcell_dim_1_u_qi[0];
                                        
                                        const int idx_BB  = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[0]) +
                                            (relative_idx_lo_1 + (j - 2) + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[0])*ghostcell_dim_0_u_qi[0]*
                                                ghostcell_dim_1_u_qi[0];
                                        
                                        const int idx_B   = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[0]) +
                                            (relative_idx_lo_1 + (j - 1) + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[0])*ghostcell_dim_0_u_qi[0]*
                                                ghostcell_dim_1_u_qi[0];
                                        
                                        const int idx_T   = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[0]) +
                                            (relative_idx_lo_1 + (j + 1) + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[0])*ghostcell_dim_0_u_qi[0]*
                                                ghostcell_dim_1_u_qi[0];
                                        
                                        const int idx_TT  = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[0]) +
                                            (relative_idx_lo_1 + (j + 2) + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[0])*ghostcell_dim_0_u_qi[0]*
                                                ghostcell_dim_1_u_qi[0];
                                        
                                        const int idx_TTT = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[0]) +
                                            (relative_idx_lo_1 + (j + 3) + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[0])*ghostcell_dim_0_u_qi[0]*
                                                ghostcell_dim_1_u_qi[0];
                                        
                                        double dudy = double(0);
                                        if (use_reciprocal[0])
                                        {
                                            dudy = (double(1)/double(60)*(double(1)/u_qi[0][idx_TTT] - double(1)/u_qi[0][idx_BBB])
                                                - double(3)/double(20)*(double(1)/u_qi[0][idx_TT] - double(1)/u_qi[0][idx_BB])
                                                + double(3)/double(4)*(double(1)/u_qi[0][idx_T] - double(1)/u_qi[0][idx_B]))/dx[1];
                                        }
                                        else
                                        {
                                            dudy = (double(1)/double(60)*(u_qi[0][idx_TTT] - u_qi[0][idx_BBB])
                                                - double(3)/double(20)*(u_qi[0][idx_TT] - u_qi[0][idx_BB])
                                                + double(3)/double(4)*(u_qi[0][idx_T] - u_qi[0][idx_B]))/dx[1];
                                        }
                                        corr = dudy - u_qi_avg_global[0][idx_fine];
                                    }
                                    else if (derivative_directions[0] == 2)
                                    {
                                        const int idx_BBB = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[0]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0] +
                                            (relative_idx_lo_2 + (k - 3) + num_ghosts_2_u_qi[0])*ghostcell_dim_0_u_qi[0]*
                                                ghostcell_dim_1_u_qi[0];
                                        
                                        const int idx_BB  = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[0]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0] +
                                            (relative_idx_lo_2 + (k - 2) + num_ghosts_2_u_qi[0])*ghostcell_dim_0_u_qi[0]*
                                                ghostcell_dim_1_u_qi[0];
                                        
                                        const int idx_B   = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[0]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0] +
                                            (relative_idx_lo_2 + (k - 1) + num_ghosts_2_u_qi[0])*ghostcell_dim_0_u_qi[0]*
                                                ghostcell_dim_1_u_qi[0];
                                        
                                        const int idx_F   = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[0]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0] +
                                            (relative_idx_lo_2 + (k + 1) + num_ghosts_2_u_qi[0])*ghostcell_dim_0_u_qi[0]*
                                                ghostcell_dim_1_u_qi[0];
                                        
                                        const int idx_FF  = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[0]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0] +
                                            (relative_idx_lo_2 + (k + 2) + num_ghosts_2_u_qi[0])*ghostcell_dim_0_u_qi[0]*
                                                ghostcell_dim_1_u_qi[0];
                                        
                                        const int idx_FFF = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[0]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0] +
                                            (relative_idx_lo_2 + (k + 3) + num_ghosts_2_u_qi[0])*ghostcell_dim_0_u_qi[0]*
                                                ghostcell_dim_1_u_qi[0];
                                        
                                        double dudz = double(0);
                                        if (use_reciprocal[0])
                                        {
                                            dudz = (double(1)/double(60)*(double(1)/u_qi[0][idx_FFF] - double(1)/u_qi[0][idx_BBB])
                                                - double(3)/double(20)*(double(1)/u_qi[0][idx_FF] - double(1)/u_qi[0][idx_BB])
                                                + double(3)/double(4)*(double(1)/u_qi[0][idx_F] - double(1)/u_qi[0][idx_B]))/dx[2];
                                        }
                                        else
                                        {
                                            dudz = (double(1)/double(60)*(u_qi[0][idx_FFF] - u_qi[0][idx_BBB])
                                                - double(3)/double(20)*(u_qi[0][idx_FF] - u_qi[0][idx_BB])
                                                + double(3)/double(4)*(u_qi[0][idx_F] - u_qi[0][idx_B]))/dx[2];
                                        }
                                        corr = dudz - u_qi_avg_global[0][idx_fine];
                                    }
                                    else
                                    {
                                        TBOX_ERROR(d_object_name
                                            << ": "
                                            << "Cannot take derivative for three-dimensional problem!\n"
                                            << "derivative_direction = " << derivative_directions[0] << " given!\n"
                                            << std::endl);
                                    }
                                    
                                    for (int qi = 1; qi < num_quantities; qi++)
                                    {
                                        if (derivative_directions[qi] == -1)
                                        {
                                            const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                                (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                                (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                    ghostcell_dim_1_u_qi[qi];
                                            
                                            if (use_reciprocal[qi])
                                            {
                                                corr *= (double(1)/(u_qi[qi][idx_qi]) - u_qi_avg_global[qi][idx_fine]);
                                            }
                                            else
                                            {
                                                corr *= (u_qi[qi][idx_qi] - u_qi_avg_global[qi][idx_fine]);
                                            }
                                        }
                                        else if (derivative_directions[qi] == 0)
                                        {
                                            const int idx_LLL = (relative_idx_lo_0 + (i - 3) + num_ghosts_0_u_qi[qi]) +
                                                (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                                (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                    ghostcell_dim_1_u_qi[qi];
                                            
                                            const int idx_LL  = (relative_idx_lo_0 + (i - 2) + num_ghosts_0_u_qi[qi]) +
                                                (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                                (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                    ghostcell_dim_1_u_qi[qi];
                                            
                                            const int idx_L   = (relative_idx_lo_0 + (i - 1) + num_ghosts_0_u_qi[qi]) +
                                                (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                                (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                    ghostcell_dim_1_u_qi[qi];
                                            
                                            const int idx_R   = (relative_idx_lo_0 + (i + 1) + num_ghosts_0_u_qi[qi]) +
                                                (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                                (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                    ghostcell_dim_1_u_qi[qi];
                                            
                                            const int idx_RR  = (relative_idx_lo_0 + (i + 2) + num_ghosts_0_u_qi[qi]) +
                                                (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                                (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                    ghostcell_dim_1_u_qi[qi];
                                            
                                            const int idx_RRR = (relative_idx_lo_0 + (i + 3) + num_ghosts_0_u_qi[qi]) +
                                                (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                                (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                    ghostcell_dim_1_u_qi[qi];
                                            
                                            double dudx = double(0);
                                            if (use_reciprocal[qi])
                                            {
                                                dudx = (double(1)/double(60)*(double(1)/u_qi[qi][idx_RRR] - double(1)/u_qi[qi][idx_LLL])
                                                    - double(3)/double(20)*(double(1)/u_qi[qi][idx_RR] - double(1)/u_qi[qi][idx_LL])
                                                    + double(3)/double(4)*(double(1)/u_qi[qi][idx_R] - double(1)/u_qi[qi][idx_L]))/dx[0];
                                            }
                                            else
                                            {
                                                dudx = (double(1)/double(60)*(u_qi[qi][idx_RRR] - u_qi[qi][idx_LLL])
                                                    - double(3)/double(20)*(u_qi[qi][idx_RR] - u_qi[qi][idx_LL])
                                                    + double(3)/double(4)*(u_qi[qi][idx_R] - u_qi[qi][idx_L]))/dx[0];
                                            }
                                            corr *= (dudx - u_qi_avg_global[qi][idx_fine]);
                                        }
                                        else if (derivative_directions[qi] == 1)
                                        {
                                            const int idx_BBB = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                                (relative_idx_lo_1 + (j - 3) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                                (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                    ghostcell_dim_1_u_qi[qi];
                                            
                                            const int idx_BB  = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                                (relative_idx_lo_1 + (j - 2) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                                (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                    ghostcell_dim_1_u_qi[qi];
                                            
                                            const int idx_B   = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                                (relative_idx_lo_1 + (j - 1) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                                (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                    ghostcell_dim_1_u_qi[qi];
                                            
                                            const int idx_T   = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                                (relative_idx_lo_1 + (j + 1) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                                (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                    ghostcell_dim_1_u_qi[qi];
                                            
                                            const int idx_TT  = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                                (relative_idx_lo_1 + (j + 2) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                                (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                    ghostcell_dim_1_u_qi[qi];
                                            
                                            const int idx_TTT = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                                (relative_idx_lo_1 + (j + 3) + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                                (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                    ghostcell_dim_1_u_qi[qi];
                                            
                                            double dudy = double(0);
                                            if (use_reciprocal[qi])
                                            {
                                                dudy = (double(1)/double(60)*(double(1)/u_qi[qi][idx_TTT] - double(1)/u_qi[qi][idx_BBB])
                                                    - double(3)/double(20)*(double(1)/u_qi[qi][idx_TT] - double(1)/u_qi[qi][idx_BB])
                                                    + double(3)/double(4)*(double(1)/u_qi[qi][idx_T] - double(1)/u_qi[qi][idx_B]))/dx[1];
                                            }
                                            else
                                            {
                                                dudy = (double(1)/double(60)*(u_qi[qi][idx_TTT] - u_qi[qi][idx_BBB])
                                                    - double(3)/double(20)*(u_qi[qi][idx_TT] - u_qi[qi][idx_BB])
                                                    + double(3)/double(4)*(u_qi[qi][idx_T] - u_qi[qi][idx_B]))/dx[1];
                                            }
                                            corr *= (dudy - u_qi_avg_global[qi][idx_fine]);
                                        }
                                        else if (derivative_directions[qi] == 2)
                                        {
                                            const int idx_BBB = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                                (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                                (relative_idx_lo_2 + (k - 3) + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                    ghostcell_dim_1_u_qi[qi];
                                            
                                            const int idx_BB  = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                                (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                                (relative_idx_lo_2 + (k - 2) + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                    ghostcell_dim_1_u_qi[qi];
                                            
                                            const int idx_B   = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                                (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                                (relative_idx_lo_2 + (k - 1) + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                    ghostcell_dim_1_u_qi[qi];
                                            
                                            const int idx_F   = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                                (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                                (relative_idx_lo_2 + (k + 1) + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                    ghostcell_dim_1_u_qi[qi];
                                            
                                            const int idx_FF  = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                                (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                                (relative_idx_lo_2 + (k + 2) + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                    ghostcell_dim_1_u_qi[qi];
                                            
                                            const int idx_FFF = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                                (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                                (relative_idx_lo_2 + (k + 3) + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                    ghostcell_dim_1_u_qi[qi];
                                            
                                            double dudz = double(0);
                                            if (use_reciprocal[qi])
                                            {
                                                dudz = (double(1)/double(60)*(double(1)/u_qi[qi][idx_FFF] - double(1)/u_qi[qi][idx_BBB])
                                                    - double(3)/double(20)*(double(1)/u_qi[qi][idx_FF] - double(1)/u_qi[qi][idx_BB])
                                                    + double(3)/double(4)*(double(1)/u_qi[qi][idx_F] - double(1)/u_qi[qi][idx_B]))/dx[2];
                                            }
                                            else
                                            {
                                                dudz = (double(1)/double(60)*(u_qi[qi][idx_FFF] - u_qi[qi][idx_BBB])
                                                    - double(3)/double(20)*(u_qi[qi][idx_FF] - u_qi[qi][idx_BB])
                                                    + double(3)/double(4)*(u_qi[qi][idx_F] - u_qi[qi][idx_B]))/dx[2];
                                            }
                                            corr *= (dudz - u_qi_avg_global[qi][idx_fine]);
                                        }
                                        else
                                        {
                                            TBOX_ERROR(d_object_name
                                                << ": "
                                                << "Cannot take derivative for three-dimensional problem!\n"
                                                << "derivative_direction = " << derivative_directions[qi] << " given!\n"
                                                << std::endl);
                                        }
                                    }
                                    
                                    corr_local[idx_fine] += (corr*weight/((double) n_overlapped));
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
        
        mpi.Allreduce(
            corr_local,
            corr_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(corr_local);
    }
    
    return correlation;
}

