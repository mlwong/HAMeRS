#include "util/equations_of_state/ideal_gas/EquationOfStateIdealGas.hpp"

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
    const std::vector<const double*>& momentum,
    const double* const total_energy,
    const std::vector<const double*>& thermo_properties) const
{
    double p = 0.0;
    
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 1);
#endif
    
    const double& gamma = *(thermo_properties[0]);
    
    const double& rho = *density;
    const double& E = *total_energy;
    
    if (d_dim == tbox::Dimension(1))
    {
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
        TBOX_ASSERT(static_cast<int>(momentum.size()) == 1);
#endif
        const double& rho_u = *(momentum[0]);
        
        p = (gamma - 1.0)*(E - 0.5*(rho_u*rho_u)/rho);
    }
    else if (d_dim == tbox::Dimension(2))
    {
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
       TBOX_ASSERT(static_cast<int>(momentum.size()) == 2);
#endif
        const double& rho_u = *(momentum[0]);
        const double& rho_v = *(momentum[1]);
        
        p = (gamma - 1.0)*(E - 0.5*(rho_u*rho_u + rho_v*rho_v)/rho);
    }
    else if (d_dim == tbox::Dimension(3))
    {
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
        TBOX_ASSERT(static_cast<int>(momentum.size()) == 3);
#endif
        const double& rho_u = *(momentum[0]);
        const double& rho_v = *(momentum[1]);
        const double& rho_w = *(momentum[2]);
        
        p = (gamma - 1.0)*(E - 0.5*(rho_u*rho_u + rho_v*rho_v + rho_w*rho_w)/rho);
    }
    
    return p;
}


/*
 * Compute the sound speed.
 */
double
EquationOfStateIdealGas::getSoundSpeed(
    const double* const density,
    const std::vector<const double*>& velocity,
    const double* const pressure,
    const std::vector<const double*>& thermo_properties) const
{
    NULL_USE(velocity);
    
    double c = 0.0;
    
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 1);
#endif
    
    const double& gamma = *(thermo_properties[0]);
    
    const double& rho = *density;
    const double& p = *pressure;
    
    c = sqrt(gamma*p/rho);
    
    return c;
}


/*
 * Compute the total energy per unit volume.
 */
double
EquationOfStateIdealGas::getTotalEnergy(
    const double* const density,
    const std::vector<const double*>& velocity,
    const double* const pressure,
    const std::vector<const double*>& thermo_properties) const
{
    double E = 0.0;
    
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 1);
#endif
    
    const double& gamma = *(thermo_properties[0]);
    
    const double& rho = *density;
    const double& p = *pressure;
    
    if (d_dim == tbox::Dimension(1))
    {
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
        TBOX_ASSERT(static_cast<int>(velocity.size()) == 1);
#endif
        const double& u = *(velocity[0]);
        
        E = p/(gamma - 1.0) + 0.5*rho*u*u;
    }
    else if (d_dim == tbox::Dimension(2))
    {
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
        TBOX_ASSERT(static_cast<int>(velocity.size()) == 2);
#endif
        const double& u = *(velocity[0]);
        const double& v = *(velocity[1]);
        
        E = p/(gamma - 1.0) + 0.5*rho*(u*u + v*v);
    }
    else if (d_dim == tbox::Dimension(3))
    {
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
        TBOX_ASSERT(static_cast<int>(velocity.size()) == 3);
#endif
        const double& u = *(velocity[0]);
        const double& v = *(velocity[1]);
        const double& w = *(velocity[2]);
        
        E = p/(gamma - 1.0) + 0.5*rho*(u*u + v*v + w*w);
    }
    
    return E;
}


/*
 * Compute the specific total enthalpy.
 */
double
EquationOfStateIdealGas::getTotalEnthalpy(
    const double* const density,
    const double* const total_energy,
    const double* const pressure) const
{
    const double& rho = *density;
    const double& E = *total_energy;
    const double& p = *pressure;
    
    double H = E + p/rho;
    
    return H;
}


/*
 * Compute the temperature.
 */
double
EquationOfStateIdealGas::getTemperature(
    const double* const density,
    const std::vector<const double*>& momentum,
    const double* const total_energy,
    const std::vector<const double*>& thermo_properties) const
{
    double T = 0.0;
    
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) == 4);
#endif
    
    const double& c_v = *(thermo_properties[3]);
    
    const double& rho = *density;
    const double& E = *total_energy;
    
    if (d_dim == tbox::Dimension(1))
    {
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
        TBOX_ASSERT(static_cast<int>(momentum.size()) == 1);
#endif
        const double& rho_u = *(momentum[0]);
        
        T = (E - 0.5*rho_u*rho_u/rho)/(rho*c_v);
    }
    else if (d_dim == tbox::Dimension(2))
    {
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
        TBOX_ASSERT(static_cast<int>(momentum.size()) == 2);
#endif
        const double& rho_u = *(momentum[0]);
        const double& rho_v = *(momentum[1]);
        
        T = (E - 0.5*(rho_u*rho_u + rho_v*rho_v)/rho)/(rho*c_v);
    }
    else if (d_dim == tbox::Dimension(3))
    {
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
        TBOX_ASSERT(static_cast<int>(momentum.size()) == 3);
#endif
        const double& rho_u = *(momentum[0]);
        const double& rho_v = *(momentum[1]);
        const double& rho_w = *(momentum[2]);
        
        T = (E - 0.5*(rho_u*rho_u + rho_v*rho_v + rho_w*rho_w)/rho)/(rho*c_v);
    }
    
    return T;
}


/*
 * Compute the partial derivative of internal energy w.r.t. pressure under constant density.
 */
double
EquationOfStateIdealGas::getIsochoricPartialInternalEnergyPartialPressure(
    const double* const density,
    const std::vector<const double*>& momentum,
    const double* const total_energy,
    const std::vector<const double*>& thermo_properties) const
{
    NULL_USE(density);
    NULL_USE(momentum);
    NULL_USE(total_energy);
    
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
    const std::vector<const double*>& momentum,
    const double* const total_energy,
    const std::vector<const double*>& thermo_properties) const
{
    NULL_USE(density);
    NULL_USE(momentum);
    NULL_USE(total_energy);
    NULL_USE(thermo_properties);
    
    return 0.0;
}
