#include "util/mixing_rules/equations_of_thermal_conductivity/Prandtl/EquationOfThermalConductivityPrandtl.hpp"

/*
 * Print all characteristics of the equation of thermal conductivity class.
 */
void
EquationOfThermalConductivityPrandtl::printClassData(
    std::ostream& os) const
{
    os << "\nPrint EquationOfThermalConductivityPrandtl object..."
       << std::endl;
       
    os << std::endl;
    
    os << "EquationOfThermalConductivityPrandtl: this = "
       << (EquationOfThermalConductivityPrandtl *)this
       << std::endl;
    
    os << "d_object_name = "
       << d_object_name
       << std::endl;
}


/*
 * Compute the thermal conductivity.
 */
double
EquationOfThermalConductivityPrandtl::getThermalConductivity(
    const double* const pressure,
    const double* const temperature,
    const std::vector<const double*>& molecular_properties) const
{
    NULL_USE(pressure);
    NULL_USE(temperature);
    
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(molecular_properties.size()) >= 3);
#endif
    
    const double& c_p = *(molecular_properties[0]);
    const double& Pr = *(molecular_properties[1]);
    const double& mu = *(molecular_properties[3]);
    
    return c_p*mu/Pr;
}
