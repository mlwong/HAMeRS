#include "equation_of_state/EquationOfStateIdealGas.hpp"

EquationOfStateIdealGas::EquationOfStateIdealGas(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const int& num_species,
    const boost::shared_ptr<tbox::Database>& equation_of_state_db,
    const THERMAL_PROCESS_ASSUMPTION& thermal_pro_assum):
        EquationOfState(
            object_name,
            dim,
            num_species,
            equation_of_state_db,
            thermal_pro_assum),
        has_Cp(false),
        has_Cv(false)
{
    if (d_equation_of_state_db->keyExists("species_gamma"))
    {
        size_t species_gamma_array_size = d_equation_of_state_db->getArraySize("species_gamma");
        if (static_cast<int>(species_gamma_array_size) == d_num_species)
        {
            d_species_gamma = d_equation_of_state_db->getDoubleVector("species_gamma");
        }
        else
        {
            TBOX_ERROR(d_object_name
                       << ": "
                       << "number of 'species_gamma' entries must be equal to 'num_species'."
                       << std::endl);
        }
    }
    else if (d_equation_of_state_db->keyExists("d_species_gamma"))
    {
        size_t species_gamma_array_size = d_equation_of_state_db->getArraySize("d_species_gamma");
        if (static_cast<int>(species_gamma_array_size) == d_num_species)
        {
            d_species_gamma = d_equation_of_state_db->getDoubleVector("d_species_gamma");
        }
        else
        {
            TBOX_ERROR(d_object_name
                       << ": "
                       << "number of 'd_species_gamma' entries must be equal to 'd_num_species'."
                       << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "Key data 'd_species_gamma' or 'species_gamma' not found in data for Equation_of_state."
                   << std::endl);
    }
    
    if (d_equation_of_state_db->keyExists("species_Cp"))
    {
        size_t species_Cp_array_size = d_equation_of_state_db->getArraySize("species_Cp");
        if (static_cast<int>(species_Cp_array_size) == d_num_species)
        {
            d_species_Cp = d_equation_of_state_db->getDoubleVector("species_Cp");
        }
        else
        {
            TBOX_ERROR(d_object_name
                       << ": "
                       << "number of 'species_Cp' entries must be equal to 'num_species'."
                       << std::endl);
        }
    }
    else if (d_equation_of_state_db->keyExists("d_species_Cp"))
    {
        size_t species_Cp_array_size = d_equation_of_state_db->getArraySize("d_species_Cp");
        if (static_cast<int>(species_Cp_array_size) == d_num_species)
        {
            d_species_Cp = d_equation_of_state_db->getDoubleVector("d_species_Cp");
        }
        else
        {
            TBOX_ERROR(d_object_name
                       << ": "
                       << "number of 'd_species_Cp' entries must be equal to 'd_num_species'."
                       << std::endl);
        }
    }
    
    if (d_equation_of_state_db->keyExists("species_Cv"))
    {
        size_t species_Cv_array_size = d_equation_of_state_db->getArraySize("species_Cv");
        if (static_cast<int>(species_Cv_array_size) == d_num_species)
        {
            d_species_Cv = d_equation_of_state_db->getDoubleVector("species_Cv");
        }
        else
        {
            TBOX_ERROR(d_object_name
                       << ": "
                       << "number of 'species_Cv' entries must be equal to 'num_species'."
                       << std::endl);
        }
    }
    else if (d_equation_of_state_db->keyExists("d_species_Cv"))
    {
        size_t species_Cv_array_size = d_equation_of_state_db->getArraySize("d_species_Cv");
        if (static_cast<int>(species_Cv_array_size) == d_num_species)
        {
            d_species_Cv = d_equation_of_state_db->getDoubleVector("d_species_Cv");
        }
        else
        {
            TBOX_ERROR(d_object_name
                       << ": "
                       << "number of 'd_species_Cv' entries must be equal to 'd_num_species'."
                       << std::endl);
        }
    }
    
    if ((thermal_pro_assum == ISOTHERMAL) &&
        (static_cast<int>(d_species_Cp.size()) != num_species) &&
        (static_cast<int>(d_species_Cv.size()) != num_species))
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "Please provide the correct number of Cp/Cv"
                   << " while isothermal process is assumed.\n"
                   << "Number of Cp/Cv provided is not equal to the total number of species."
                   << std::endl);
    }
    
    if (static_cast<int>(d_species_Cp.size()) == num_species)
    {
        has_Cp = true;
    }
    
    if (static_cast<int>(d_species_Cv.size()) == num_species)
    {
        has_Cv = true;
    }
}


/*
 * Print all characteristics of the equation of state class.
 */
void
EquationOfStateIdealGas::printClassData(
    std::ostream& os) const
{
    os << "\nPrint EquationOfStateIdealGas object..."
       << std::endl;
       
    os << std::endl;
    
    os << "EquationOfStateIdealGas: this = "
       << (EquationOfStateIdealGas *)this
       << std::endl;
    os << "d_object_name = "
       << d_object_name
       << std::endl;
    
    os << "d_thermal_pro_assum = "
       << d_thermal_pro_assum
       << std::endl;
    
    os << "d_species_gamma = "
       << d_species_gamma[0];
    for (int si = 1; si < d_num_species; si++)
    {
        os << " , "
           << d_species_gamma[si];
    }
    os << std::endl;
    
    if (has_Cp)
    {
        os << "d_species_Cp = "
           << d_species_Cp[0];
        for (int si = 1; si < d_num_species; si++)
        {
            os << " , "
               << d_species_Cp[si];
        }
        os << std::endl;
    }
    
    if (has_Cv)
    {
        os << "d_species_Cv = "
           << d_species_Cv[0];
        for (int si = 1; si < d_num_species; si++)
        {
            os << " , "
               << d_species_Cv[si];
        }
        os << std::endl;
    }
}


/*
 * Put the characteristics of the convective flux reconstruction class
 * into the restart database.
 */
void
EquationOfStateIdealGas::putToRestart(
   const boost::shared_ptr<tbox::Database>& restart_db) const
{
    restart_db->putString("d_equation_of_state", "IDEAL_GAS");

    restart_db->putDoubleVector("d_species_gamma", d_species_gamma);
    
    if (has_Cp)
    {
        restart_db->putDoubleVector("d_species_Cp", d_species_Cp);
    }
    
    if (has_Cv)
    {
        restart_db->putDoubleVector("d_species_Cv", d_species_Cv);
    }
}


/*
 * Compute the pressure for single-species flow.
 */
double
EquationOfStateIdealGas::getPressure(
   const double* const density,
   const std::vector<const double*> momentum,
   const double* const total_energy)
{
    double p = 0.0;
    
    if (d_num_species == 1)
    {
        const double& rho = *density;
        const double& E = *total_energy;
        
        if (d_dim == tbox::Dimension(1))
        {
            // Get the reference to vector time-dependent variables.
            const double& rho_u = *(momentum[0]);
            
            p = (d_species_gamma[0] - 1)*
                (E - 0.5*(rho_u*rho_u)/rho);
        }
        else if (d_dim == tbox::Dimension(2))
        {
            // Get the reference to vector time-dependent variables.
            const double& rho_u = *(momentum[0]);
            const double& rho_v = *(momentum[1]);
            
            p = (d_species_gamma[0] - 1)*
                (E - 0.5*(rho_u*rho_u + rho_v*rho_v)/rho);
        }
        else if (d_dim == tbox::Dimension(3))
        {
            // Get the pointer of vector time-dependent variables.
            const double& rho_u = *(momentum[0]);
            const double& rho_v = *(momentum[1]);
            const double& rho_w = *(momentum[2]);
            
            p = (d_species_gamma[0] - 1)*
                (E - 0.5*(rho_u*rho_u + rho_v*rho_v + rho_w*rho_w)/rho);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "Mass fractions/volume fractions are required to"
                   << " compute pressure for multi-species flow."
                   << std::endl);
    }
    
    return p;
}


/*
 * Compute the pressure for multi-species flow with total density and
 * mass fractions.
 */
double
EquationOfStateIdealGas::getPressureWithMassFraction(
   const double* const density,
   const std::vector<const double*> momentum,
   const double* const total_energy,
   const std::vector<const double*> mass_fraction)
{
    double p = 0.0;
    
    if (d_num_species > 1)
    {
        if (d_thermal_pro_assum == ISOTHERMAL)
        {
            const double& rho = *density;
            const double& E = *total_energy;
            
            double mixture_gamma = getMixtureGammaWithMassFraction(
                mass_fraction);
            
            if (d_dim == tbox::Dimension(1))
            {
                // Get the reference to vector time-dependent variables.
                const double& rho_u = *(momentum[0]);
                
                p = (mixture_gamma - 1)*
                    (E - 0.5*(rho_u*rho_u)/rho);
            }
            else if (d_dim == tbox::Dimension(2))
            {
                // Get the reference to vector time-dependent variables.
                const double& rho_u = *(momentum[0]);
                const double& rho_v = *(momentum[1]);
                
                p = (mixture_gamma - 1)*
                    (E - 0.5*(rho_u*rho_u + rho_v*rho_v)/rho);
            }
            else if (d_dim == tbox::Dimension(3))
            {
                // Get the pointer of vector time-dependent variables.
                const double& rho_u = *(momentum[0]);
                const double& rho_v = *(momentum[1]);
                const double& rho_w = *(momentum[2]);
                
                p = (mixture_gamma - 1)*
                    (E - 0.5*(rho_u*rho_u + rho_v*rho_v + rho_w*rho_w)/rho);
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                       << ": "
                       << "getPressureWithMassFraction() can only be used"
                       << " with isothermal assumption."
                       << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "Mass fractions are not required to compute"
                   << " pressure for single-species flow"
                   << std::endl);
    }
    
    return p;
}


/*
 * Compute the pressure for multi-species flow with partial densities and
 * mass fractions.
 */
double
EquationOfStateIdealGas::getPressureWithMassFraction(
   const std::vector<const double*> partial_density,
   const std::vector<const double*> momentum,
   const double* const total_energy,
   const std::vector<const double*> mass_fraction)
{
    double p = 0.0;
    
    if (d_num_species > 1)
    {
        if (d_thermal_pro_assum == ISOTHERMAL)
        {
            const double& E = *total_energy;
            const double rho = getTotalDensity(partial_density);
            
            double mixture_gamma = getMixtureGammaWithMassFraction(
                mass_fraction);
            
            if (d_dim == tbox::Dimension(1))
            {
                // Get the reference to vector time-dependent variables.
                const double& rho_u = *(momentum[0]);
                
                p = (mixture_gamma - 1)*
                    (E - 0.5*(rho_u*rho_u)/rho);
            }
            else if (d_dim == tbox::Dimension(2))
            {
                // Get the reference to vector time-dependent variables.
                const double& rho_u = *(momentum[0]);
                const double& rho_v = *(momentum[1]);
                
                p = (mixture_gamma - 1)*
                    (E - 0.5*(rho_u*rho_u + rho_v*rho_v)/rho);
            }
            else if (d_dim == tbox::Dimension(3))
            {
                // Get the pointer of vector time-dependent variables.
                const double& rho_u = *(momentum[0]);
                const double& rho_v = *(momentum[1]);
                const double& rho_w = *(momentum[2]);
                
                p = (mixture_gamma - 1)*
                    (E - 0.5*(rho_u*rho_u + rho_v*rho_v + rho_w*rho_w)/rho);
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                       << ": "
                       << "getPressureWithMassFraction() can only be used"
                       << " with isothermal assumption."
                       << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "Mass fractions are not required to compute"
                   << " pressure for single-species flow"
                   << std::endl);
    }
    
    return p;
}


/*
 * Compute the pressure for multi-species flow with total density and
 * volume fractions.
 */
double
EquationOfStateIdealGas::getPressureWithVolumeFraction(
   const double* const density,
   const std::vector<const double*> momentum,
   const double* const total_energy,
   const std::vector<const double*> volume_fraction)
{
    double p = 0.0;
    
    if (d_num_species > 1)
    {
        if (d_thermal_pro_assum == ISOBARIC)
        {
            const double& rho = *density;
            const double& E = *total_energy;
            
            double mixture_gamma = getMixtureGammaWithVolumeFraction(
                volume_fraction);
            
            if (d_dim == tbox::Dimension(1))
            {
                // Get the reference to vector time-dependent variables.
                const double& rho_u = *(momentum[0]);
                
                p = (mixture_gamma - 1)*
                    (E - 0.5*(rho_u*rho_u)/rho);
            }
            else if (d_dim == tbox::Dimension(2))
            {
                // Get the reference to vector time-dependent variables.
                const double& rho_u = *(momentum[0]);
                const double& rho_v = *(momentum[1]);
                
                p = (mixture_gamma - 1)*
                    (E - 0.5*(rho_u*rho_u + rho_v*rho_v)/rho);
            }
            else if (d_dim == tbox::Dimension(3))
            {
                // Get the pointer of vector time-dependent variables.
                const double& rho_u = *(momentum[0]);
                const double& rho_v = *(momentum[1]);
                const double& rho_w = *(momentum[2]);
                
                p = (mixture_gamma- 1)*
                    (E - 0.5*(rho_u*rho_u + rho_v*rho_v + rho_w*rho_w)/rho);
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                       << ": "
                       << "getPressureWithVolumeFraction() without partial densities"
                       << " can only be used with isobaric assumption."
                       << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "Volume fractions are not required to compute"
                   << " pressure for single-species flow"
                   << std::endl);
    }
    
    return p;
}


/*
 * Compute the pressure for multi-species flow with partial densities and
 * volume fractions.
 */
double
EquationOfStateIdealGas::getPressureWithVolumeFraction(
   const std::vector<const double*> partial_density,
   const std::vector<const double*> momentum,
   const double* const total_energy,
   const std::vector<const double*> volume_fraction)
{
    double p = 0.0;
    
    if (d_num_species > 1)
    {
        if (d_thermal_pro_assum == ISOTHERMAL)
        {
            const double& E = *total_energy;
            const double rho = getTotalDensity(partial_density);
            
            double mixture_gamma = getMixtureGammaWithPartialDensity(
                partial_density);
            
            if (d_dim == tbox::Dimension(1))
            {
                // Get the reference to vector time-dependent variables.
                const double& rho_u = *(momentum[0]);
                
                p = (mixture_gamma - 1)*
                    (E - 0.5*(rho_u*rho_u)/rho);
            }
            else if (d_dim == tbox::Dimension(2))
            {
                // Get the reference to vector time-dependent variables.
                const double& rho_u = *(momentum[0]);
                const double& rho_v = *(momentum[1]);
                
                p = (mixture_gamma - 1)*
                    (E - 0.5*(rho_u*rho_u + rho_v*rho_v)/rho);
            }
            else if (d_dim == tbox::Dimension(3))
            {
                // Get the pointer of vector time-dependent variables.
                const double& rho_u = *(momentum[0]);
                const double& rho_v = *(momentum[1]);
                const double& rho_w = *(momentum[2]);
                
                p = (mixture_gamma- 1)*
                    (E - 0.5*(rho_u*rho_u + rho_v*rho_v + rho_w*rho_w)/rho);
            }
        }
        else if (d_thermal_pro_assum == ISOBARIC)
        {
            const double& E = *total_energy;
            const double rho = getTotalDensity(partial_density);
            
            double mixture_gamma = getMixtureGammaWithVolumeFraction(
                volume_fraction);
            
            if (d_dim == tbox::Dimension(1))
            {
                // Get the reference to vector time-dependent variables.
                const double& rho_u = *(momentum[0]);
                
                p = (mixture_gamma - 1)*
                    (E - 0.5*(rho_u*rho_u)/rho);
            }
            else if (d_dim == tbox::Dimension(2))
            {
                // Get the reference to vector time-dependent variables.
                const double& rho_u = *(momentum[0]);
                const double& rho_v = *(momentum[1]);
                
                p = (mixture_gamma - 1)*
                    (E - 0.5*(rho_u*rho_u + rho_v*rho_v)/rho);
            }
            else if (d_dim == tbox::Dimension(3))
            {
                // Get the pointer of vector time-dependent variables.
                const double& rho_u = *(momentum[0]);
                const double& rho_v = *(momentum[1]);
                const double& rho_w = *(momentum[2]);
                
                p = (mixture_gamma- 1)*
                    (E - 0.5*(rho_u*rho_u + rho_v*rho_v + rho_w*rho_w)/rho);
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                       << ": "
                       << "getPressureWithVolumeFraction() with partial densities"
                       << " can only be used with isothermal or isobaric assumptions."
                       << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "Volume fractions are not required to compute"
                   << " pressure for single-species flow"
                   << std::endl);
    }
    
    return p;
}


/*
 * Compute the sound speed for single-species flow.
 */
double
EquationOfStateIdealGas::getSoundSpeed(
   const double* const density,
   const std::vector<const double*> momentum,
   const double* const total_energy)
{
    double c = 0.0;
    
    if (d_num_species == 1)
    {
        const double& rho = *density;
        
        // Compute the pressure
        double p = getPressure(
            density,
            momentum,
            total_energy);
        
        c = sqrt(d_species_gamma[0]*p/rho);
    }
    else
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "Mass fractions/volume fractions are required to"
                   << " compute sound speed for multi-species flow."
                   << std::endl);
    }
    
    return c;
}


/*
 * Compute the sound speed for single-species flow with pressure.
 * 
 * Recommended to use with getPressure() if both pressure and sound speed
 * data is required.
 */
double
EquationOfStateIdealGas::getSoundSpeedWithPressure(
    const double* const density,
    const double* const pressure)
{
    double c = 0.0;
    
    if (d_num_species == 1)
    {
        const double& rho = *density;
        const double& p = *pressure;
        
        c = sqrt(d_species_gamma[0]*p/rho);
    }
    else
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "Mass fractions/volume fractions are required to"
                   << " compute sound speed for multi-species flow."
                   << std::endl);
    }
    
    return c;
}


/*
 * Compute the sound speed for multi-species flow with total density and
 * mass fractions.
 */
double
EquationOfStateIdealGas::getSoundSpeedWithMassFraction(
    const double* const density,
    const std::vector<const double*> momentum,
    const double* const total_energy,
    const std::vector<const double*> mass_fraction)
{
    double c = 0.0;
    
    if (d_num_species > 1)
    {
        if (d_thermal_pro_assum == ISOTHERMAL)
        {
            const double& rho = *density;
            
            double mixture_gamma = getMixtureGammaWithMassFraction(
                mass_fraction);
            
            // Compute the pressure
            double p = getPressureWithMassFraction(
                density,
                momentum,
                total_energy,
                mass_fraction);
            
            c = sqrt(mixture_gamma*p/rho);
        }
        else
        {
            TBOX_ERROR(d_object_name
                       << ": "
                       << "getSoundSpeedWithMassFraction() can only be used"
                       << " with isothermal assumption."
                       << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "Mass fractions are not required to"
                   << " compute sound speed for multi-species flow."
                   << std::endl);
    }
    
    return c;
}


/*
 * Compute the sound speed for multi-species flow with partial densities and
 * mass fractions.
 */
double
EquationOfStateIdealGas::getSoundSpeedWithMassFraction(
    const std::vector<const double*> partial_density,
    const std::vector<const double*> momentum,
    const double* const total_energy,
    const std::vector<const double*> mass_fraction)
{
    double c = 0.0;
    
    if (d_num_species > 1)
    {
        if (d_thermal_pro_assum == ISOTHERMAL)
        {
            const double rho = getTotalDensity(partial_density);
            const double* const density = &rho;
            
            double mixture_gamma = getMixtureGammaWithMassFraction(
                mass_fraction);
            
            // Compute the pressure
            double p = getPressureWithMassFraction(
                density,
                momentum,
                total_energy,
                mass_fraction);
            
            c = sqrt(mixture_gamma*p/rho);
        }
        else
        {
            TBOX_ERROR(d_object_name
                       << ": "
                       << "getSoundSpeedWithMassFraction() can only be used"
                       << " with isothermal assumption."
                       << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "Mass fractions are not required to"
                   << " compute sound speed for multi-species flow."
                   << std::endl);
    }
    
    return c;
}


/*
 * Compute the sound speed for multi-species flow with total density,
 * mass fractions and pressure.
 * 
 * Recommended to use with getPressureWithMassFraction() if both pressure
 * and sound speed data is required.
 */
double
EquationOfStateIdealGas::getSoundSpeedWithMassFractionAndPressure(
    const double* const density,
    const std::vector<const double*> mass_fraction,
    const double* const pressure)
{
    double c = 0.0;
    
    if (d_num_species > 1)
    {
        if (d_thermal_pro_assum == ISOTHERMAL)
        {
            const double& rho = *density;
            const double& p = *pressure;
            
            double mixture_gamma = getMixtureGammaWithMassFraction(
                mass_fraction);
            
            c = sqrt(mixture_gamma*p/rho);
        }
        else
        {
            TBOX_ERROR(d_object_name
                       << ": "
                       << "getSoundSpeedWithMassFractionAndPressure() can only be used"
                       << " with isothermal assumption."
                       << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "Mass fractions are not required to"
                   << " compute sound speed for multi-species flow."
                   << std::endl);
    }
    
    return c;
}


/*
 * Compute the sound speed for multi-species flow with partial densities,
 * mass fractions and pressure.
 * 
 * Recommended to use with getPressureWithMassFraction() if both pressure
 * and sound speed data is required.
 */
double
EquationOfStateIdealGas::getSoundSpeedWithMassFractionAndPressure(
    const std::vector<const double*> partial_density,
    const std::vector<const double*> mass_fraction,
    const double* const pressure)
{
    double c = 0.0;
    
    if (d_num_species > 1)
    {
        if (d_thermal_pro_assum == ISOTHERMAL)
        {
            const double rho = getTotalDensity(partial_density);
            const double& p = *pressure;
            
            double mixture_gamma = getMixtureGammaWithMassFraction(
                mass_fraction);
            
            c = sqrt(mixture_gamma*p/rho);
        }
        else
        {
            TBOX_ERROR(d_object_name
                       << ": "
                       << "getSoundSpeedWithMassFractionAndPressure() can only be used"
                       << " with isothermal assumption."
                       << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "Mass fractions are not required to"
                   << " compute sound speed for multi-species flow."
                   << std::endl);
    }
    
    return c;
}


/*
 * Compute the sound speed for multi-species flow with total density and
 * volume fractions.
 */
double
EquationOfStateIdealGas::getSoundSpeedWithVolumeFraction(
    const double* const density,
    const std::vector<const double*> momentum,
    const double* const total_energy,
    const std::vector<const double*> volume_fraction)
{
    double c = 0.0;
    
    if (d_num_species > 1)
    {
        if (d_thermal_pro_assum == ISOBARIC)
        {
            const double& rho = *density;
            
            double mixture_gamma = getMixtureGammaWithVolumeFraction(
                volume_fraction);
            
            // Compute the pressure
            double p = getPressureWithVolumeFraction(
                density,
                momentum,
                total_energy,
                volume_fraction);
            
            c = sqrt(mixture_gamma*p/rho);
        }
        else
        {
            TBOX_ERROR(d_object_name
                       << ": "
                       << "getSoundSpeedWithVolumeFraction() without partial densities"
                       << " can only be used with isobaric assumption."
                       << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "Volume fractions are not required to"
                   << " compute sound speed for multi-species flow."
                   << std::endl);
    }
    
    return c;
}


/*
 * Compute the sound speed for multi-species flow with partial densities and
 * volume fractions.
 */
double
EquationOfStateIdealGas::getSoundSpeedWithVolumeFraction(
    const std::vector<const double*> partial_density,
    const std::vector<const double*> momentum,
    const double* const total_energy,
    const std::vector<const double*> volume_fraction)
{
    double c = 0.0;
    
    if (d_num_species > 1)
    {
        if (d_thermal_pro_assum == ISOTHERMAL)
        {
            const double rho = getTotalDensity(partial_density);
            const double* const density = &rho;
            
            double mixture_gamma = getMixtureGammaWithPartialDensity(
                partial_density);
            
            // Compute the pressure
            double p = getPressureWithVolumeFraction(
                density,
                momentum,
                total_energy,
                volume_fraction);
            
            c = sqrt(mixture_gamma*p/rho);
        }
        else if (d_thermal_pro_assum == ISOBARIC)
        {
            const double rho = getTotalDensity(partial_density);
            const double* const density = &rho;
            
            double mixture_gamma = getMixtureGammaWithVolumeFraction(
                volume_fraction);
            
            // Compute the pressure
            double p = getPressureWithVolumeFraction(
                density,
                momentum,
                total_energy,
                volume_fraction);
            
            c = sqrt(mixture_gamma*p/rho);
        }
        else
        {
            TBOX_ERROR(d_object_name
                       << ": "
                       << "getSoundSpeedWithVolumeFraction() with partial densities"
                       << " can only be used with isothermal or isobaric assumptions."
                       << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "Volume fractions are not required to"
                   << " compute sound speed for multi-species flow."
                   << std::endl);
    }
    
    return c;
}


/*
 * Compute the sound speed for multi-species flow with total density,
 * volume fractions and pressure.
 * 
 * Recommended to use with getPressureWithVolumeFraction() if both pressure
 * and sound speed data is required.
 */
double
EquationOfStateIdealGas::getSoundSpeedWithVolumeFractionAndPressure(
    const double* const density,
    const std::vector<const double*> volume_fraction,
    const double* const pressure)
{
    double c = 0.0;
    
    if (d_num_species > 1)
    {
        if (d_thermal_pro_assum == ISOBARIC)
        {
            const double& rho = *density;
            const double& p = *pressure;
            
            double mixture_gamma = getMixtureGammaWithVolumeFraction(
                volume_fraction);
            
            c = sqrt(mixture_gamma*p/rho);
        }
        else
        {
            TBOX_ERROR(d_object_name
                       << ": "
                       << "getSoundSpeedWithVolumeFractionAndPressure() without partial densities"
                       << " can only be used with isobaric assumption."
                       << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "Volume fractions are not required to"
                   << " compute sound speed for multi-species flow."
                   << std::endl);
    }
    
    return c;
}


/*
 * Compute the sound speed for multi-species flow with partial densities,
 * volume fractions and pressure.
 * Recommended to use with getPressureWithVolumeFraction() if both pressure
 * and sound speed data is required.
 */
double
EquationOfStateIdealGas::getSoundSpeedWithVolumeFractionAndPressure(
    const std::vector<const double*> partial_density,
    const std::vector<const double*> volume_fraction,
    const double* const pressure)
{
    double c = 0.0;
    
    if (d_num_species > 1)
    {
        if (d_thermal_pro_assum == ISOTHERMAL)
        {
            const double rho = getTotalDensity(partial_density);
            const double& p = *pressure;
            
            double mixture_gamma = getMixtureGammaWithPartialDensity(
                partial_density);
            
            c = sqrt(mixture_gamma*p/rho);
        }
        else if (d_thermal_pro_assum == ISOBARIC)
        {
            const double rho = getTotalDensity(partial_density);
            const double& p = *pressure;
            
            double mixture_gamma = getMixtureGammaWithVolumeFraction(
                volume_fraction);
            
            c = sqrt(mixture_gamma*p/rho);
        }
        else
        {
            TBOX_ERROR(d_object_name
                       << ": "
                       << "getSoundSpeedWithVolumeFractionAndPressure() with partial densities"
                       << " can only be used with isothermal or isobaric assumptions."
                       << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "Volume fractions are not required to"
                   << " compute sound speed for multi-species flow."
                   << std::endl);
    }
    
    return c;
}


/*
 * Compute the total energy for single-species flow.
 */
double
EquationOfStateIdealGas::getTotalEnergy(
    const double* const density,
    const std::vector<const double*> velocity,
    const double* const pressure)
{
    double E = 0.0;
    
    if (d_num_species == 1)
    {
        const double& rho = *density;
        const double& p = *pressure;
        
        if (d_dim == tbox::Dimension(1))
        {
            // Get the reference to vector time-dependent variables.
            const double& u = *(velocity[0]);
            
            E = p/(d_species_gamma[0] - 1.0) + 0.5*rho*u*u;
        }
        else if (d_dim == tbox::Dimension(2))
        {
            // Get the reference to vector time-dependent variables.
            const double& u = *(velocity[0]);
            const double& v = *(velocity[1]);
            
            E = p/(d_species_gamma[0] - 1.0) + 0.5*rho*(u*u + v*v);
        }
        else if (d_dim == tbox::Dimension(3))
        {
            // Get the pointer of vector time-dependent variables.
            const double& u = *(velocity[0]);
            const double& v = *(velocity[1]);
            const double& w = *(velocity[2]);
            
            E = p/(d_species_gamma[0] - 1.0) + 0.5*rho*(u*u + v*v + w*w);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "Mass fractions/volume fractions are required to"
                   << " compute total energy for multi-species flow."
                   << std::endl);
    }
    
    return E;
}


/*
 * Compute the total energy for multi-species flow with total density and
 * mass fractions.
 */
double
EquationOfStateIdealGas::getTotalEnergyWithMassFraction(
   const double* const density,
   const std::vector<const double*> velocity,
   const double* const pressure,
   const std::vector<const double*> mass_fraction)
{
    double E = 0.0;
    
    if (d_num_species > 1)
    {
        if (d_thermal_pro_assum == ISOTHERMAL)
        {
            const double& rho = *density;
            const double& p = *pressure;
            
            double mixture_gamma = getMixtureGammaWithMassFraction(
                mass_fraction);
            
            if (d_dim == tbox::Dimension(1))
            {
                // Get the reference to vector time-dependent variables.
                const double& u = *(velocity[0]);
                
                E = p/(mixture_gamma - 1.0) + 0.5*rho*u*u;
            }
            else if (d_dim == tbox::Dimension(2))
            {
                // Get the reference to vector time-dependent variables.
                const double& u = *(velocity[0]);
                const double& v = *(velocity[1]);
                
                E = p/(mixture_gamma - 1.0) + 0.5*rho*(u*u + v*v);
            }
            else if (d_dim == tbox::Dimension(3))
            {
                // Get the pointer of vector time-dependent variables.
                const double& u = *(velocity[0]);
                const double& v = *(velocity[1]);
                const double& w = *(velocity[2]);
                
                E = p/(mixture_gamma - 1.0) + 0.5*rho*(u*u + v*v + w*w);
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                       << ": "
                       << "getTotalEnergyWithMassFraction() can only be used"
                       << " with isothermal assumption."
                       << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "Mass fractions are not required to compute"
                   << " total energy for single-species flow"
                   << std::endl);
    }
    
    return E;
}


/*
 * Compute the total energy for multi-species flow with partial densities and
 * mass fractions.
 */
double
EquationOfStateIdealGas::getTotalEnergyWithMassFraction(
   const std::vector<const double*> partial_density,
   const std::vector<const double*> velocity,
   const double* const pressure,
   const std::vector<const double*> mass_fraction)
{
    double E = 0.0;
    
    if (d_num_species > 1)
    {
        if (d_thermal_pro_assum == ISOTHERMAL)
        {
            const double& p = *pressure;
            const double rho = getTotalDensity(partial_density);
            
            double mixture_gamma = getMixtureGammaWithMassFraction(
                mass_fraction);
            
            if (d_dim == tbox::Dimension(1))
            {
                // Get the reference to vector time-dependent variables.
                const double& u = *(velocity[0]);
                
                E = p/(mixture_gamma - 1.0) + 0.5*rho*u*u;
            }
            else if (d_dim == tbox::Dimension(2))
            {
                // Get the reference to vector time-dependent variables.
                const double& u = *(velocity[0]);
                const double& v = *(velocity[1]);
                
                E = p/(mixture_gamma - 1.0) + 0.5*rho*(u*u + v*v);
            }
            else if (d_dim == tbox::Dimension(3))
            {
                // Get the pointer of vector time-dependent variables.
                const double& u = *(velocity[0]);
                const double& v = *(velocity[1]);
                const double& w = *(velocity[2]);
                
                E = p/(mixture_gamma - 1.0) + 0.5*rho*(u*u + v*v + w*w);
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                       << ": "
                       << "getTotalEnergyWithMassFraction() can only be used"
                       << " with isothermal assumption."
                       << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "Mass fractions are not required to compute"
                   << " total energy for single-species flow"
                   << std::endl);
    }
    
    return E;
}


/*
 * Compute the total energy for multi-species flow with total density and
 * volume fractions.
 */
double
EquationOfStateIdealGas::getTotalEnergyWithVolumeFraction(
   const double* const density,
   const std::vector<const double*> velocity,
   const double* const pressure,
   const std::vector<const double*> volume_fraction)
{
    double E = 0.0;
    
    if (d_num_species > 1)
    {
        if (d_thermal_pro_assum == ISOBARIC)
        {
            const double& rho = *density;
            const double& p = *pressure;
            
            double mixture_gamma = getMixtureGammaWithVolumeFraction(
                volume_fraction);
            
            if (d_dim == tbox::Dimension(1))
            {
                // Get the reference to vector time-dependent variables.
                const double& u = *(velocity[0]);
                
                E = p/(mixture_gamma - 1.0) + 0.5*rho*u*u;
            }
            else if (d_dim == tbox::Dimension(2))
            {
                // Get the reference to vector time-dependent variables.
                const double& u = *(velocity[0]);
                const double& v = *(velocity[1]);
                
                E = p/(mixture_gamma - 1.0) + 0.5*rho*(u*u + v*v);
            }
            else if (d_dim == tbox::Dimension(3))
            {
                // Get the pointer of vector time-dependent variables.
                const double& u = *(velocity[0]);
                const double& v = *(velocity[1]);
                const double& w = *(velocity[2]);
                
                E = p/(mixture_gamma - 1.0) + 0.5*rho*(u*u + v*v + w*w);
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                       << ": "
                       << "getTotalEnergyWithVolumeFraction() without partial densities"
                       << " can only be used with isobaric assumption."
                       << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "Volume fractions are not required to compute"
                   << " total energy for single-species flow"
                   << std::endl);
    }
    
    return E;
}


/*
 * Compute the total energy for multi-species flow with partial densities and
 * volume fractions.
 */
double
EquationOfStateIdealGas::getTotalEnergyWithVolumeFraction(
   const std::vector<const double*> partial_density,
   const std::vector<const double*> velocity,
   const double* const pressure,
   const std::vector<const double*> volume_fraction)
{
    double E = 0.0;
    
    if (d_num_species > 1)
    {
        if (d_thermal_pro_assum == ISOTHERMAL)
        {
            const double& p = *pressure;
            const double rho = getTotalDensity(partial_density);
            
            double mixture_gamma = getMixtureGammaWithPartialDensity(
                partial_density);
            
            if (d_dim == tbox::Dimension(1))
            {
                // Get the reference to vector time-dependent variables.
                const double& u = *(velocity[0]);
                
                E = p/(mixture_gamma - 1.0) + 0.5*rho*u*u;
            }
            else if (d_dim == tbox::Dimension(2))
            {
                // Get the reference to vector time-dependent variables.
                const double& u = *(velocity[0]);
                const double& v = *(velocity[1]);
                
                E = p/(mixture_gamma - 1.0) + 0.5*rho*(u*u + v*v);
            }
            else if (d_dim == tbox::Dimension(3))
            {
                // Get the pointer of vector time-dependent variables.
                const double& u = *(velocity[0]);
                const double& v = *(velocity[1]);
                const double& w = *(velocity[2]);
                
                E = p/(mixture_gamma - 1.0) + 0.5*rho*(u*u + v*v + w*w);
            }
        }
        else if (d_thermal_pro_assum == ISOBARIC)
        {
            const double& p = *pressure;
            const double rho = getTotalDensity(partial_density);
            
            double mixture_gamma = getMixtureGammaWithVolumeFraction(
                volume_fraction);
            
            if (d_dim == tbox::Dimension(1))
            {
                // Get the reference to vector time-dependent variables.
                const double& u = *(velocity[0]);
                
                E = p/(mixture_gamma - 1.0) + 0.5*rho*u*u;
            }
            else if (d_dim == tbox::Dimension(2))
            {
                // Get the reference to vector time-dependent variables.
                const double& u = *(velocity[0]);
                const double& v = *(velocity[1]);
                
                E = p/(mixture_gamma - 1.0) + 0.5*rho*(u*u + v*v);
            }
            else if (d_dim == tbox::Dimension(3))
            {
                // Get the pointer of vector time-dependent variables.
                const double& u = *(velocity[0]);
                const double& v = *(velocity[1]);
                const double& w = *(velocity[2]);
                
                E = p/(mixture_gamma - 1.0) + 0.5*rho*(u*u + v*v + w*w);
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                       << ": "
                       << "getTotalEnergyWithVolumeFraction() with partial densities"
                       << " can only be used with isothermal or isobaric assumptions."
                       << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "Volume fractions are not required to compute"
                   << " total energy for single-species flow"
                   << std::endl);
    }
    
    return E;
}


/*
 * Compute a thermodynamic property of a particular species
 * with mass fractions.
 */
double
EquationOfStateIdealGas::getMixtureThermodynamicPropertyWithMassFraction(
    const std::string& property_name,
    const std::vector<const double*> mass_fraction)
{
    double therm_propty_value = 0.0;
    
    if (property_name == "gamma")
    {
        therm_propty_value =  getMixtureGammaWithMassFraction(mass_fraction);
    }
    else if (property_name == "Cp")
    {
        if (has_Cp)
        {
            // NOT YET IMPLEMENTED
            therm_propty_value =  0;
        }
        else
        {
            TBOX_ERROR(d_object_name
                       << ": "
                       << "Cp is not registered."
                       << std::endl);
        }
    }
    else if (property_name == "Cv")
    {
        if (has_Cv)
        {
            // NOT YET IMPLEMENTED
            therm_propty_value =  0;
        }
        else
        {
            TBOX_ERROR(d_object_name
                       << ": "
                       << "Cv is not registered."
                       << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "Unknown thermodynamic property name = '"
                   << property_name
                   << "' provided."
                   << std::endl);
    }
    
    return therm_propty_value;
}


/*
 * Compute a thermodynamic property of a particular species
 * with volume fractions.
 */
double
EquationOfStateIdealGas::getMixtureThermodynamicPropertyWithVolumeFraction(
    const std::string& property_name,
    const std::vector<const double*> volume_fraction)
{
    double therm_propty_value = 0.0;
    
    if (property_name == "gamma")
    {
        therm_propty_value =  getMixtureGammaWithVolumeFraction(volume_fraction);
    }
    else if (property_name == "Cp")
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "Mixture's Cp cannot be computed from volume fractions."
                   << std::endl);
    }
    else if (property_name == "Cv")
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "Mixture's Cv cannot be computed from volume fractions."
                   << std::endl);
    }
    else
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "Unknown thermodynamic property name = '"
                   << property_name
                   << "' provided."
                   << std::endl);
    }
    
    return therm_propty_value;
}


/*
 * Get a thermodynamic property of a particular species.
 */
double&
EquationOfStateIdealGas::getSpeciesThermodynamicProperty(
    const std::string& property_name,
    const int species_index)
{
    double* therm_propty_value = NULL;
    if (property_name == "gamma")
    {
        therm_propty_value =  &(d_species_gamma[species_index]);
    }
    else if (property_name == "Cp")
    {
        if (has_Cp)
        {
            therm_propty_value =  &(d_species_Cp[species_index]);
        }
        else
        {
            TBOX_ERROR(d_object_name
                       << ": "
                       << "Cp is not registered."
                       << std::endl);
        }
    }
    else if (property_name == "Cv")
    {
        if (has_Cv)
        {
            therm_propty_value = &(d_species_Cv[species_index]);
        }
        else
        {
            TBOX_ERROR(d_object_name
                       << ": "
                       << "Cv is not registered."
                       << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "Unknown thermodynamic property name = '"
                   << property_name
                   << "' provided."
                   << std::endl);
    }
    
    return *therm_propty_value;
}


/*
 * Determine whether a thermodynamic property is registered.
 */
bool
EquationOfStateIdealGas::hasThermodynamicProperty(
    const std::string& property_name)
{
    bool has_therm_propty_value = false;
    if (property_name == "gamma")
    {
        has_therm_propty_value = true;
    }
    else if (property_name == "Cp")
    {
        if (has_Cp)
        {
            has_therm_propty_value = true;
        }
    }
    else if (property_name == "Cv")
    {
        if (has_Cv)
        {
            has_therm_propty_value = true;
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "Unknown thermodynamic property name = '"
                   << property_name
                   << "' provided."
                   << std::endl);
    }
    
    return has_therm_propty_value;
}


/*
 * Compute the mixture's gamma for multi-species flow with partial densities.
 */
double
EquationOfStateIdealGas::getMixtureGammaWithPartialDensity(
    const std::vector<const double*> partial_density)
{
    double mixture_gamma = 0.0;
    
    if (d_thermal_pro_assum == ISOTHERMAL)
    {
        if (static_cast<int>(partial_density.size()) == d_num_species)
        {
            const double rho = getTotalDensity(partial_density);
            double mixture_Cp = 0.0;
            double mixture_Cv = 0.0;
            
            for (int si = 0; si < d_num_species; si++)
            {
                const double& Y = (*(partial_density[si]))/rho;
                mixture_Cp += d_species_Cp[si]*Y;
                mixture_Cv += d_species_Cv[si]*Y;
            }
            
            mixture_gamma = mixture_Cp/mixture_Cv;
        }
        else
        {
            TBOX_ERROR(d_object_name
                       << ": "
                       << "Number of partial densities provided is not"
                       << " equal to the total number of species."
                       << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "getMixtureGammaWithPartialDensity() can only be used"
                   << " with isothermal assumption."
                   << std::endl);
    }
    
    return mixture_gamma;
}


/*
 * Compute the mixture's gamma for multi-species flow with mass fractions.
 */
double
EquationOfStateIdealGas::getMixtureGammaWithMassFraction(
    const std::vector<const double*> mass_fraction)
{
    double mixture_gamma = 0.0;
    
    if (d_thermal_pro_assum == ISOTHERMAL)
    {
        if (static_cast<int>(mass_fraction.size()) == d_num_species)
        {
            double mixture_Cp = 0.0;
            double mixture_Cv = 0.0;
            
            for (int si = 0; si < d_num_species; si++)
            {
                const double& Y = *(mass_fraction[si]);
                mixture_Cp += d_species_Cp[si]*Y;
                mixture_Cv += d_species_Cv[si]*Y;
            }
            
            mixture_gamma = mixture_Cp/mixture_Cv;
        }
        else if (static_cast<int>(mass_fraction.size()) == d_num_species - 1)
        {
            double mixture_Cp = 0.0;
            double mixture_Cv = 0.0;
            double Y_last = 1.0;
            
            for (int si = 0; si < d_num_species - 1; si++)
            {
                const double& Y = *(mass_fraction[si]);
                mixture_Cp += d_species_Cp[si]*Y;
                mixture_Cv += d_species_Cv[si]*Y;
                
                // Compute the mass fraction of the last species.
                Y_last -= Y;
            }
            
            // Add the contribution from the last species.
            mixture_Cp += d_species_Cp.back()*Y_last;
            mixture_Cv += d_species_Cv.back()*Y_last;
            
            mixture_gamma = mixture_Cp/mixture_Cv;
        }
        else
        {
            TBOX_ERROR(d_object_name
                       << ": "
                       << "Number of mass fractions provided is not"
                       << " equal to the total number of species."
                       << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "getMixtureGammaWithMassFraction() can only be used"
                   << " with isothermal assumption."
                   << std::endl);
    }
    
    return mixture_gamma;
}


/*
 * Compute the mixture's gamma for multi-species flow with volume fractions.
 */
double
EquationOfStateIdealGas::getMixtureGammaWithVolumeFraction(
    const std::vector<const double*> volume_fraction)
{
    double mixture_gamma = 0.0;
    
    if (d_thermal_pro_assum == ISOBARIC)
    {
        if (static_cast<int>(volume_fraction.size()) == d_num_species)
        {
            double sum = 0.0;
            
            for (int si = 0; si < d_num_species; si++)
            {
                const double& Z = *(volume_fraction[si]);
                sum += Z/(d_species_gamma[si] - 1.0);
            }
            
            mixture_gamma = 1.0/sum + 1.0;
        }
        else if (static_cast<int>(volume_fraction.size()) == d_num_species - 1)
        {
            double Z_last = 1.0;
            double sum = 0.0;
            
            for (int si = 0; si < d_num_species - 1; si++)
            {
                const double& Z = *(volume_fraction[si]);
                sum += Z/(d_species_gamma[si] - 1.0);
                
                // Compute the mass fraction of the last species.
                Z_last -= Z;
            }
            
            // Add the contribution from the last species.
            sum += Z_last/(d_species_gamma.back() - 1.0);
            
            mixture_gamma = 1.0/sum + 1.0;
        }
        else
        {
            TBOX_ERROR(d_object_name
                       << ": "
                       << "Number of volume fractions provided is not"
                       << " equal to the total number of species."
                       << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
                       << ": "
                       << "getMixtureGammaWithVolumeFraction() can only be used"
                       << " with isobaric assumption."
                       << std::endl);
    }
    
    return mixture_gamma;
}
