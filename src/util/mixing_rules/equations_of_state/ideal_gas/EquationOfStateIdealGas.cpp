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
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 1);
#endif
    
    const double& gamma = *(thermo_properties[0]);
    
    const double& rho = *density;
    const double& epsilon = *internal_energy;
    
    double p = (gamma - 1.0)*rho*epsilon;
    
    return p;
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
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 1);
#endif
    
    const double& gamma = *(thermo_properties[0]);
    
    const double& rho = *density;
    const double& p = *pressure;
    
    double c = sqrt(gamma*p/rho);
    
    return c;
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
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 1);
#endif
    
    const double& gamma = *(thermo_properties[0]);
    
    const double& rho = *density;
    const double& p = *pressure;
    
    double epsilon = p/((gamma - 1.0)*rho);
    
    return epsilon;
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
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 1);
#endif
    
    const double& gamma = *(thermo_properties[0]);
    
    const double& rho = *density;
    const double& p = *pressure;
    
    double h = gamma*p/((gamma - 1.0)*rho);
    
    return h;
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
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) == 4);
#endif
    
    const double& gamma = *(thermo_properties[0]);
    const double& c_v = *(thermo_properties[3]);
    
    const double& rho = *density;
    const double& p = *pressure;
    
    double T = p/((gamma - 1.0)*c_v*rho);
    
    return T;
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
    
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 1);
#endif
    
    const double& gamma = *(thermo_properties[0]);
    
    double xi = 1.0/(gamma - 1.0);
    
    return xi;
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
    
    return 0.0;
}
