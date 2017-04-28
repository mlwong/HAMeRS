#include "util/mixing_rules/equations_of_state/ideal_gas/EquationOfStateMixingRulesIdealGas.hpp"

EquationOfStateMixingRulesIdealGas::EquationOfStateMixingRulesIdealGas(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const int& num_species,
    const MIXING_CLOSURE_MODEL::TYPE& mixing_closure_model,
    const boost::shared_ptr<tbox::Database>& equation_of_state_mixing_rules_db):
        EquationOfStateMixingRules(
            object_name,
            dim,
            num_species,
            mixing_closure_model,
            equation_of_state_mixing_rules_db)
{
    d_equation_of_state.reset(new EquationOfStateIdealGas(
        "d_equation_of_state",
        dim));
    
    /*
     * Get the ratio of specific heats of each species from the database.
     */
    
    if (equation_of_state_mixing_rules_db->keyExists("species_gamma"))
    {
        size_t species_gamma_array_size = equation_of_state_mixing_rules_db->getArraySize("species_gamma");
        if (static_cast<int>(species_gamma_array_size) == d_num_species)
        {
            d_species_gamma = equation_of_state_mixing_rules_db->getDoubleVector("species_gamma");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'species_gamma' entries must be equal to 'num_species'."
                << std::endl);
        }
    }
    else if (equation_of_state_mixing_rules_db->keyExists("d_species_gamma"))
    {
        size_t species_gamma_array_size = equation_of_state_mixing_rules_db->getArraySize("d_species_gamma");
        if (static_cast<int>(species_gamma_array_size) == d_num_species)
        {
            d_species_gamma = equation_of_state_mixing_rules_db->getDoubleVector("d_species_gamma");
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
            << "Key data 'species_gamma'/'d_species_gamma'"
            << "not found in data for Equation_of_state_mixing_rules."
            << std::endl);
    }
    
    /*
     * Get the gas constant of each species from the database.
     */
    
    if (equation_of_state_mixing_rules_db->keyExists("species_R"))
    {
        size_t species_R_array_size = equation_of_state_mixing_rules_db->getArraySize("species_R");
        if (static_cast<int>(species_R_array_size) == d_num_species)
        {
            d_species_R = equation_of_state_mixing_rules_db->getDoubleVector("species_R");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'species_R' entries must be equal to 'num_species'."
                << std::endl);
        }
    }
    else if (equation_of_state_mixing_rules_db->keyExists("d_species_R"))
    {
        size_t species_R_array_size = equation_of_state_mixing_rules_db->getArraySize("d_species_R");
        if (static_cast<int>(species_R_array_size) == d_num_species)
        {
            d_species_R = equation_of_state_mixing_rules_db->getDoubleVector("d_species_R");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'd_species_R' entries must be equal to 'd_num_species'."
                << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Key data 'species_R'/'d_species_R'"
            << "not found in data for Equation_of_state_mixing_rules."
            << std::endl);
    }
    
    /*
     * Compute the specific heats of each species.
     */
    
    d_species_c_p.reserve(d_num_species);
    for (int si = 0; si < d_num_species; si++)
    {
        d_species_c_p.push_back(d_species_gamma[si]/(d_species_gamma[si] - 1.0)*d_species_R[si]);
    }
    
    d_species_c_v.reserve(d_num_species);
    for (int si = 0; si < d_num_species; si++)
    {
        d_species_c_v.push_back(1.0/(d_species_gamma[si] - 1.0)*d_species_R[si]);
    }
}


/*
 * Print all characteristics of the equation of state class.
 */
void
EquationOfStateMixingRulesIdealGas::printClassData(
    std::ostream& os) const
{
    os << "\nPrint EquationOfStateMixingRulesIdealGas object..."
       << std::endl;
    
    os << std::endl;
    os << "EquationOfStateMixingRulesIdealGas: this = "
       << (EquationOfStateMixingRulesIdealGas *)this
       << std::endl;
    
    os << "d_object_name = "
       << d_object_name
       << std::endl;
    
    os << "d_mixing_closure_model = "
       << d_mixing_closure_model
       << std::endl;
    
    /*
     * Print the ratio of specific heats of each species.
     */
    
    os << "d_species_gamma = ";
    for (int si = 0; si < d_num_species - 1; si++)
    {
        os << d_species_gamma[si] << ", ";
    }
    os << d_species_gamma[d_num_species - 1];
    os << std::endl;
    
    /*
     * Print the gas constant of each species.
     */
    
    os << "d_species_R = ";
    for (int si = 0; si < d_num_species - 1; si++)
    {
        os << d_species_R[si] << ", ";
    }
    os << d_species_R[d_num_species - 1];
    os << std::endl;
}


/*
 * Put the characteristics of the equation of state mixing rules class into the restart
 * database.
 */
void
EquationOfStateMixingRulesIdealGas::putToRestart(
    const boost::shared_ptr<tbox::Database>& restart_db) const
{
    restart_db->putDoubleVector("d_species_gamma", d_species_gamma);
    restart_db->putDoubleVector("d_species_R", d_species_R);
}


/*
 * Compute the pressure of the mixture with isothermal and isobaric assumptions.
 */
double
EquationOfStateMixingRulesIdealGas::getPressure(
    const double* const density,
    const double* const internal_energy,
    const std::vector<const double*>& mass_fraction) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT((d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC) ||
                (d_mixing_closure_model == MIXING_CLOSURE_MODEL::NO_MODEL && d_num_species == 1));
    TBOX_ASSERT((static_cast<int>(mass_fraction.size()) == d_num_species) ||
                (static_cast<int>(mass_fraction.size()) == d_num_species - 1));
#endif
    
    // Get the mixture thermodynamic properties.
    std::vector<double> mixture_thermo_properties;
    std::vector<double*> mixture_thermo_properties_ptr;
    std::vector<const double*> mixture_thermo_properties_const_ptr;
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    mixture_thermo_properties.resize(num_thermo_properties);
    mixture_thermo_properties_ptr.reserve(num_thermo_properties);
    mixture_thermo_properties_const_ptr.reserve(num_thermo_properties);
    
    for (int ti = 0; ti < num_thermo_properties; ti++)
    {
        mixture_thermo_properties_ptr.push_back(&mixture_thermo_properties[ti]);
        mixture_thermo_properties_const_ptr.push_back(&mixture_thermo_properties[ti]);
    }
    
    getMixtureThermodynamicProperties(
        mixture_thermo_properties_ptr,
        mass_fraction);
    
    return d_equation_of_state->getPressure(
        density,
        internal_energy,
        mixture_thermo_properties_const_ptr);
}


/*
 * Compute the pressure of the mixture with isobaric assumption.
 */
double
EquationOfStateMixingRulesIdealGas::getPressure(
    const double* const density,
    const double* const internal_energy,
    const std::vector<const double*>& mass_fraction,
    const std::vector<const double*>& volume_fraction) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOBARIC);
    TBOX_ASSERT((static_cast<int>(mass_fraction.size()) == d_num_species) ||
                (static_cast<int>(mass_fraction.size()) == d_num_species - 1));
    TBOX_ASSERT((static_cast<int>(volume_fraction.size()) == d_num_species) ||
                (static_cast<int>(volume_fraction.size()) == d_num_species - 1));
#endif
    
    NULL_USE(mass_fraction);
    
    // Get the mixture thermodynamic properties.
    std::vector<double> mixture_thermo_properties;
    std::vector<double*> mixture_thermo_properties_ptr;
    std::vector<const double*> mixture_thermo_properties_const_ptr;
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    mixture_thermo_properties.resize(num_thermo_properties);
    mixture_thermo_properties_ptr.reserve(num_thermo_properties);
    mixture_thermo_properties_const_ptr.reserve(num_thermo_properties);
    
    for (int ti = 0; ti < num_thermo_properties; ti++)
    {
        mixture_thermo_properties_ptr.push_back(&mixture_thermo_properties[ti]);
        mixture_thermo_properties_const_ptr.push_back(&mixture_thermo_properties[ti]);
    }
    
    getMixtureThermodynamicProperties(
        mixture_thermo_properties_ptr,
        volume_fraction);
    
    return d_equation_of_state->getPressure(
        density,
        internal_energy,
        mixture_thermo_properties_const_ptr);
}


/*
 * Compute the sound speed of the mixture with isothermal and isobaric assumptions.
 */
double
EquationOfStateMixingRulesIdealGas::getSoundSpeed(
    const double* const density,
    const double* const pressure,
    const std::vector<const double*>& mass_fraction) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT((d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC) ||
                (d_mixing_closure_model == MIXING_CLOSURE_MODEL::NO_MODEL && d_num_species == 1));
    TBOX_ASSERT((static_cast<int>(mass_fraction.size()) == d_num_species) ||
                (static_cast<int>(mass_fraction.size()) == d_num_species - 1));
#endif
    
    // Get the mixture thermodynamic properties.
    std::vector<double> mixture_thermo_properties;
    std::vector<double*> mixture_thermo_properties_ptr;
    std::vector<const double*> mixture_thermo_properties_const_ptr;
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    mixture_thermo_properties.resize(num_thermo_properties);
    mixture_thermo_properties_ptr.reserve(num_thermo_properties);
    mixture_thermo_properties_const_ptr.reserve(num_thermo_properties);
    
    for (int ti = 0; ti < num_thermo_properties; ti++)
    {
        mixture_thermo_properties_ptr.push_back(&mixture_thermo_properties[ti]);
        mixture_thermo_properties_const_ptr.push_back(&mixture_thermo_properties[ti]);
    }
    
    getMixtureThermodynamicProperties(
        mixture_thermo_properties_ptr,
        mass_fraction);
    
    return d_equation_of_state->getSoundSpeed(
        density,
        pressure,
        mixture_thermo_properties_const_ptr);
}


/*
 * Compute the sound speed of the mixture with isobaric assumption.
 */
double
EquationOfStateMixingRulesIdealGas::getSoundSpeed(
    const double* const density,
    const double* const pressure,
    const std::vector<const double*>& mass_fraction,
    const std::vector<const double*>& volume_fraction) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOBARIC);
    TBOX_ASSERT((static_cast<int>(mass_fraction.size()) == d_num_species) ||
                (static_cast<int>(mass_fraction.size()) == d_num_species - 1));
    TBOX_ASSERT((static_cast<int>(volume_fraction.size()) == d_num_species) ||
                (static_cast<int>(volume_fraction.size()) == d_num_species - 1));
#endif
    
    NULL_USE(mass_fraction);
    
    // Get the mixture thermodynamic properties.
    std::vector<double> mixture_thermo_properties;
    std::vector<double*> mixture_thermo_properties_ptr;
    std::vector<const double*> mixture_thermo_properties_const_ptr;
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    mixture_thermo_properties.resize(num_thermo_properties);
    mixture_thermo_properties_ptr.reserve(num_thermo_properties);
    mixture_thermo_properties_const_ptr.reserve(num_thermo_properties);
    
    for (int ti = 0; ti < num_thermo_properties; ti++)
    {
        mixture_thermo_properties_ptr.push_back(&mixture_thermo_properties[ti]);
        mixture_thermo_properties_const_ptr.push_back(&mixture_thermo_properties[ti]);
    }
    
    getMixtureThermodynamicProperties(
        mixture_thermo_properties_ptr,
        volume_fraction);
    
    return d_equation_of_state->getSoundSpeed(
        density,
        pressure,
        mixture_thermo_properties_const_ptr);
}


/*
 * Compute the specific internal energy of the mixture with isothermal and isobaric assumptions.
 */
double
EquationOfStateMixingRulesIdealGas::getInternalEnergy(
    const double* const density,
    const double* const pressure,
    const std::vector<const double*>& mass_fraction) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT((d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC) ||
                (d_mixing_closure_model == MIXING_CLOSURE_MODEL::NO_MODEL && d_num_species == 1));
    TBOX_ASSERT((static_cast<int>(mass_fraction.size()) == d_num_species) ||
                (static_cast<int>(mass_fraction.size()) == d_num_species - 1));
#endif
    
    // Get the mixture thermodynamic properties.
    std::vector<double> mixture_thermo_properties;
    std::vector<double*> mixture_thermo_properties_ptr;
    std::vector<const double*> mixture_thermo_properties_const_ptr;
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    mixture_thermo_properties.resize(num_thermo_properties);
    mixture_thermo_properties_ptr.reserve(num_thermo_properties);
    mixture_thermo_properties_const_ptr.reserve(num_thermo_properties);
    
    for (int ti = 0; ti < num_thermo_properties; ti++)
    {
        mixture_thermo_properties_ptr.push_back(&mixture_thermo_properties[ti]);
        mixture_thermo_properties_const_ptr.push_back(&mixture_thermo_properties[ti]);
    }
    
    getMixtureThermodynamicProperties(
        mixture_thermo_properties_ptr,
        mass_fraction);
    
    return d_equation_of_state->getInternalEnergy(
        density,
        pressure,
        mixture_thermo_properties_const_ptr);
}


/*
 * Compute the specific internal energy of the mixture with isobaric assumption.
 */
double
EquationOfStateMixingRulesIdealGas::getInternalEnergy(
    const double* const density,
    const double* const pressure,
    const std::vector<const double*>& mass_fraction,
    const std::vector<const double*>& volume_fraction) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOBARIC);
    TBOX_ASSERT((static_cast<int>(mass_fraction.size()) == d_num_species) ||
                (static_cast<int>(mass_fraction.size()) == d_num_species - 1));
    TBOX_ASSERT((static_cast<int>(volume_fraction.size()) == d_num_species) ||
                (static_cast<int>(volume_fraction.size()) == d_num_species - 1));
#endif
    
    NULL_USE(mass_fraction);
    
    // Get the mixture thermodynamic properties.
    std::vector<double> mixture_thermo_properties;
    std::vector<double*> mixture_thermo_properties_ptr;
    std::vector<const double*> mixture_thermo_properties_const_ptr;
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    mixture_thermo_properties.resize(num_thermo_properties);
    mixture_thermo_properties_ptr.reserve(num_thermo_properties);
    mixture_thermo_properties_const_ptr.reserve(num_thermo_properties);
    
    for (int ti = 0; ti < num_thermo_properties; ti++)
    {
        mixture_thermo_properties_ptr.push_back(&mixture_thermo_properties[ti]);
        mixture_thermo_properties_const_ptr.push_back(&mixture_thermo_properties[ti]);
    }
    
    getMixtureThermodynamicProperties(
        mixture_thermo_properties_ptr,
        volume_fraction);
    
    return d_equation_of_state->getInternalEnergy(
        density,
        pressure,
        mixture_thermo_properties_const_ptr);
}


/*
 * Compute the temperature of the mixture with isothermal and isobaric assumptions.
 */
double
EquationOfStateMixingRulesIdealGas::getTemperature(
    const double* const density,
    const double* const pressure,
    const std::vector<const double*>& mass_fraction) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT((d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC) ||
                (d_mixing_closure_model == MIXING_CLOSURE_MODEL::NO_MODEL && d_num_species == 1));
    TBOX_ASSERT((static_cast<int>(mass_fraction.size()) == d_num_species) ||
                (static_cast<int>(mass_fraction.size()) == d_num_species - 1));
#endif
    
    // Get the mixture thermodynamic properties.
    std::vector<double> mixture_thermo_properties;
    std::vector<double*> mixture_thermo_properties_ptr;
    std::vector<const double*> mixture_thermo_properties_const_ptr;
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    mixture_thermo_properties.resize(num_thermo_properties);
    mixture_thermo_properties_ptr.reserve(num_thermo_properties);
    mixture_thermo_properties_const_ptr.reserve(num_thermo_properties);
    
    for (int ti = 0; ti < num_thermo_properties; ti++)
    {
        mixture_thermo_properties_ptr.push_back(&mixture_thermo_properties[ti]);
        mixture_thermo_properties_const_ptr.push_back(&mixture_thermo_properties[ti]);
    }
    
    getMixtureThermodynamicProperties(
        mixture_thermo_properties_ptr,
        mass_fraction);
    
    return d_equation_of_state->getTemperature(
        density,
        pressure,
        mixture_thermo_properties_const_ptr);
}


/*
 * Compute the specific internal energy of the mixture from temperature with isothermal
 * and isobaric assumptions.
 */
double
EquationOfStateMixingRulesIdealGas::getInternalEnergyFromTemperature(
    const double* const density,
    const double* const temperature,
    const std::vector<const double*>& mass_fraction) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT((d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC) ||
                (d_mixing_closure_model == MIXING_CLOSURE_MODEL::NO_MODEL && d_num_species == 1));
    TBOX_ASSERT((static_cast<int>(mass_fraction.size()) == d_num_species) ||
                (static_cast<int>(mass_fraction.size()) == d_num_species - 1));
#endif
    
    // Get the mixture thermodynamic properties.
    std::vector<double> mixture_thermo_properties;
    std::vector<double*> mixture_thermo_properties_ptr;
    std::vector<const double*> mixture_thermo_properties_const_ptr;
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    mixture_thermo_properties.resize(num_thermo_properties);
    mixture_thermo_properties_ptr.reserve(num_thermo_properties);
    mixture_thermo_properties_const_ptr.reserve(num_thermo_properties);
    
    for (int ti = 0; ti < num_thermo_properties; ti++)
    {
        mixture_thermo_properties_ptr.push_back(&mixture_thermo_properties[ti]);
        mixture_thermo_properties_const_ptr.push_back(&mixture_thermo_properties[ti]);
    }
    
    getMixtureThermodynamicProperties(
        mixture_thermo_properties_ptr,
        mass_fraction);
    
    return d_equation_of_state->getInternalEnergyFromTemperature(
        density,
        temperature,
        mixture_thermo_properties_const_ptr);
}


/*
 * Compute the isochoric specific heat capacity of mixture with isothermal assumption.
 */
double
EquationOfStateMixingRulesIdealGas::getIsochoricSpecificHeatCapacity(
    const double* const density,
    const double* const pressure,
    const std::vector<const double*>& mass_fraction) const
{
    double c_v = 0.0;
    
    if (static_cast<int>(mass_fraction.size()) == d_num_species)
    {
        for (int si = 0; si < d_num_species; si++)
        {
            c_v += *(mass_fraction[si])*d_species_c_v[si];
        }
    }
    else if (static_cast<int>(mass_fraction.size()) == d_num_species - 1)
    {
        double Y_last = 1.0;
        
        for (int si = 0; si < d_num_species - 1; si++)
        {
            c_v += *(mass_fraction[si])*d_species_c_v[si];
            
            // Compute the mass fraction of the last species.
            Y_last -= *(mass_fraction[si]);
        }
        
        // Add the contribution from the last species.
        c_v += Y_last*d_species_c_v.back();
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Number of mass fractions provided is not"
            << " equal to the total number of species or (total number of species - 1)."
            << std::endl);
    }
    
    return c_v;
}


/*
 * Compute the isobaric specific heat capacity of mixture with isothermal assumption.
 */
double
EquationOfStateMixingRulesIdealGas::getIsobaricSpecificHeatCapacity(
    const double* const density,
    const double* const pressure,
    const std::vector<const double*>& mass_fraction) const
{
    double c_p = 0.0;
    
    if (static_cast<int>(mass_fraction.size()) == d_num_species)
    {
        for (int si = 0; si < d_num_species; si++)
        {
            c_p += *(mass_fraction[si])*d_species_c_p[si];
        }
    }
    else if (static_cast<int>(mass_fraction.size()) == d_num_species - 1)
    {
        double Y_last = 1.0;
        
        for (int si = 0; si < d_num_species - 1; si++)
        {
            c_p += *(mass_fraction[si])*d_species_c_p[si];
            
            // Compute the mass fraction of the last species.
            Y_last -= *(mass_fraction[si]);
        }
        
        // Add the contribution from the last species.
        c_p += Y_last*d_species_c_p.back();
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Number of mass fractions provided is not"
            << " equal to the total number of species or (total number of species - 1)."
            << std::endl);
    }
    
    return c_p;
}


/*
 * Get the thermodynamic properties of a species.
 */
void
EquationOfStateMixingRulesIdealGas::getSpeciesThermodynamicProperties(
    std::vector<double*>& species_thermo_properties,
    const int& species_index) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(species_thermo_properties.size()) == 4);
    TBOX_ASSERT(species_index >= 0);
    TBOX_ASSERT(species_index < d_num_species);
#endif
    
    // Get references to gamma and R of species.
    double& gamma = *(species_thermo_properties[0]);
    double& R = *(species_thermo_properties[1]);
    
    // Get references to c_p and c_v of species.
    double& c_p = *(species_thermo_properties[2]);
    double& c_v = *(species_thermo_properties[3]);
    
    gamma = d_species_gamma[species_index];
    R = d_species_R[species_index];
    
    c_p = d_species_c_p[species_index];
    c_v = d_species_c_v[species_index];
}


/*
 * Get the number of thermodynamic properties of the mixture.
 */
int
EquationOfStateMixingRulesIdealGas::getNumberOfMixtureThermodynamicProperties() const
{
    int num_of_thermo_properties = 0;
    
    switch (d_mixing_closure_model)
    {
        case MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC:
        {
            num_of_thermo_properties = 4;
            
            break;
        }
        case MIXING_CLOSURE_MODEL::ISOBARIC:
        {
            num_of_thermo_properties = 1;
            
            break;
        }
        case MIXING_CLOSURE_MODEL::NO_MODEL:
        {
            num_of_thermo_properties = 0;
            
            break;
        }
    }
    
    return num_of_thermo_properties;
}


/*
 * Get the thermodynamic properties of the mixture.
 */
void
EquationOfStateMixingRulesIdealGas::getMixtureThermodynamicProperties(
    std::vector<double*>& mixture_thermo_properties,
    const std::vector<const double*>& species_fraction) const
{
    switch (d_mixing_closure_model)
    {
        case MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC:
        {
            getMixtureThermodynamicPropertiesWithMassFraction(
                mixture_thermo_properties,
                species_fraction);
            
            break;
        }
        case MIXING_CLOSURE_MODEL::ISOBARIC:
        {
            getMixtureThermodynamicPropertiesWithVolumeFraction(
                mixture_thermo_properties,
                species_fraction);
            
            break;
        }
        case MIXING_CLOSURE_MODEL::NO_MODEL:
        {
            TBOX_WARNING(d_object_name
                << ": "
                << "No thermodynamic properties of mixture are returned!\n"
                << "d_mixing_closure_model = "
                << d_mixing_closure_model
                << std::endl);
            
            break;
        }
    }
}


/*
 * Compute the thermodynamic properties of the mixture with mass fraction.
 */
void
EquationOfStateMixingRulesIdealGas::getMixtureThermodynamicPropertiesWithMassFraction(
    std::vector<double*>& mixture_thermo_properties,
    const std::vector<const double*>& mass_fraction) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(mixture_thermo_properties.size()) == 4);
#endif
    
    // Get references to gamma and R of mixture.
    double& gamma = *(mixture_thermo_properties[0]);
    double& R = *(mixture_thermo_properties[1]);
    
    // Get references to c_p and c_v of mixture.
    double& c_p = *(mixture_thermo_properties[2]);
    double& c_v = *(mixture_thermo_properties[3]);
    
    c_p = 0.0;
    c_v = 0.0;
    
    if (static_cast<int>(mass_fraction.size()) == d_num_species)
    {
        for (int si = 0; si < d_num_species; si++)
        {
            c_p += *(mass_fraction[si])*d_species_c_p[si];
            c_v += *(mass_fraction[si])*d_species_c_v[si];
        }
    }
    else if (static_cast<int>(mass_fraction.size()) == d_num_species - 1)
    {
        double Y_last = 1.0;
        
        for (int si = 0; si < d_num_species - 1; si++)
        {
            c_p += *(mass_fraction[si])*d_species_c_p[si];
            c_v += *(mass_fraction[si])*d_species_c_v[si];
            
            // Compute the mass fraction of the last species.
            Y_last -= *(mass_fraction[si]);
        }
        
        // Add the contribution from the last species.
        c_p += Y_last*d_species_c_p.back();
        c_v += Y_last*d_species_c_v.back();
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Number of mass fractions provided is not"
            << " equal to the total number of species or (total number of species - 1)."
            << std::endl);
    }
    
    gamma = c_p/c_v;
    R = c_p - c_v;
}


/*
 * Compute the thermodynamic properties of the mixture with volume fraction.
 */
void
EquationOfStateMixingRulesIdealGas::getMixtureThermodynamicPropertiesWithVolumeFraction(
    std::vector<double*>& mixture_thermo_properties,
    const std::vector<const double*>& volume_fraction) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(mixture_thermo_properties.size()) == 1);
#endif
    
    // Get reference to gamma of mixture.
    double& gamma = *(mixture_thermo_properties[0]);
    double xi = 0.0;
    
    if (static_cast<int>(volume_fraction.size()) == d_num_species)
    {
        for (int si = 0; si < d_num_species; si++)
        {
            xi += *(volume_fraction[si])/(d_species_gamma[si] - 1.0);
        }
    }
    else if (static_cast<int>(volume_fraction.size()) == d_num_species - 1)
    {
        double Z_last = 1.0;
        
        for (int si = 0; si < d_num_species - 1; si++)
        {
            xi += *(volume_fraction[si])/(d_species_gamma[si] - 1.0);
            
            // Compute the volume fraction of the last species.
            Z_last -= *(volume_fraction[si]);
        }
        
        // Add the contribution from the last species.
        xi += Z_last/(d_species_gamma.back() - 1.0);
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Number of volume fractions provided is not"
            << " equal to the total number of species or (total number of species - 1)."
            << std::endl);
    }
    
    gamma = 1.0/xi + 1.0;
}
