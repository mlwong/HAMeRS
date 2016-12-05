#include "util/mixing_rules/equations_of_thermal_conductivity/constant/EquationOfThermalConductivityConstant.hpp"

/*
 * Print all characteristics of the equation of thermal conductivity class.
 */
void
EquationOfThermalConductivityConstant::printClassData(
    std::ostream& os) const
{
    os << "\nPrint EquationOfThermalConductivityConstant object..."
       << std::endl;
       
    os << std::endl;
    
    os << "EquationOfThermalConductivityConstant: this = "
       << (EquationOfThermalConductivityConstant *)this
       << std::endl;
    
    os << "d_object_name = "
       << d_object_name
       << std::endl;
}


/*
 * Compute the thermal conductivity.
 */
double
EquationOfThermalConductivityConstant::getThermalConductivity(
    const double* const pressure,
    const double* const temperature,
    const std::vector<const double*>& molecular_properties) const
{
    NULL_USE(pressure);
    NULL_USE(temperature);
    
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(molecular_properties.size()) >= 1);
#endif
    
    const double& kappa = *(molecular_properties[0]);
    
    return kappa;
}
