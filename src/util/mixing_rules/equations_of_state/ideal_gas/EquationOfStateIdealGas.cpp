#include "util/mixing_rules/equations_of_state/ideal_gas/EquationOfStateIdealGas.hpp"

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
}


/*
 * Compute the pressure.
 */
double
EquationOfStateIdealGas::getPressure(
    const double* const density,
    const double* const internal_energy,
    const std::vector<const double*>& thermo_properties) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 1);
#endif
    
    const double& gamma = *(thermo_properties[0]);
    
    const double& rho = *density;
    const double& epsilon = *internal_energy;
    
    return (gamma - 1.0)*rho*epsilon; // Return p.
}


/*
 * Compute the sound speed.
 */
double
EquationOfStateIdealGas::getSoundSpeed(
    const double* const density,
    const double* const pressure,
    const std::vector<const double*>& thermo_properties) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 1);
#endif
    
    const double& gamma = *(thermo_properties[0]);
    
    const double& rho = *density;
    const double& p = *pressure;
    
    return sqrt(gamma*p/rho); // Return c.
}


/*
 * Compute the specific internal energy.
 */
double
EquationOfStateIdealGas::getInternalEnergy(
    const double* const density,
    const double* const pressure,
    const std::vector<const double*>& thermo_properties) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 1);
#endif
    
    const double& gamma = *(thermo_properties[0]);
    
    const double& rho = *density;
    const double& p = *pressure;
    
    return p/((gamma - 1.0)*rho); // Return epsilon.
}


/*
 * Compute the specific enthalpy.
 */
double
EquationOfStateIdealGas::getEnthalpy(
    const double* const density,
    const double* const pressure,
    const std::vector<const double*>& thermo_properties) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 1);
#endif
    
    const double& gamma = *(thermo_properties[0]);
    
    const double& rho = *density;
    const double& p = *pressure;
    
    return gamma*p/((gamma - 1.0)*rho); // Return h.
}


/*
 * Compute the temperature.
 */
double
EquationOfStateIdealGas::getTemperature(
    const double* const density,
    const double* const pressure,
    const std::vector<const double*>& thermo_properties) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) == 4);
#endif
    
    const double& gamma = *(thermo_properties[0]);
    const double& c_v = *(thermo_properties[3]);
    
    const double& rho = *density;
    const double& p = *pressure;
    
    return p/((gamma - 1.0)*c_v*rho); // Return T.
}


/*
 * Compute the specific internal energy from temperature.
 */
double
EquationOfStateIdealGas::getInternalEnergyFromTemperature(
    const double* const density,
    const double* const temperature,
    const std::vector<const double*>& thermo_properties) const
{
    NULL_USE(density);
    
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) == 4);
#endif
    
    const double& c_v = *(thermo_properties[3]);
    
    const double& T = *temperature;
    
    return c_v*T; // Return epsilon.
}


/*
 * Compute the partial derivative of internal energy w.r.t. pressure under constant density.
 */
double
EquationOfStateIdealGas::getIsochoricPartialInternalEnergyPartialPressure(
    const double* const density,
    const double* const pressure,
    const std::vector<const double*>& thermo_properties) const
{
    NULL_USE(density);
    NULL_USE(pressure);
    
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 1);
#endif
    
    const double& gamma = *(thermo_properties[0]);
    
    return 1.0/(gamma - 1.0); // Return xi.
}


/*
 * Compute the partial derivative of internal energy w.r.t. density under constant pressure.
 */
double
EquationOfStateIdealGas::getIsobaricPartialInternalEnergyPartialDensity(
    const double* const density,
    const double* const pressure,
    const std::vector<const double*>& thermo_properties) const
{
    NULL_USE(density);
    NULL_USE(pressure);
    NULL_USE(thermo_properties);
    
    return 0.0; // Return delta.
}


/*
 * Compute the partial derivative of internal energy w.r.t. density under constant pressure.
 */
double
EquationOfStateIdealGas::getDensity(
    const double* const pressure,
    const double* const temperature,
    const std::vector<const double*>& thermo_properties) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 2);
#endif
    
    const double& R = *(thermo_properties[1]);
    
    const double& p = *pressure;
    const double& T = *temperature;
    
    return p/(R*T); // Return rho.
}
