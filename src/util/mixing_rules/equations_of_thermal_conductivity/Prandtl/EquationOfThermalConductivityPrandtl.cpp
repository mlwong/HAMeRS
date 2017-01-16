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
    
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(molecular_properties.size()) >= 3);
#endif
    
    const double& c_p = *(molecular_properties[0]);
    const double& Pr = *(molecular_properties[1]);
    
    std::vector<double> mu_molecular_properties;
    std::vector<const double*> mu_molecular_properties_const_ptr;
    
    mu_molecular_properties.reserve(static_cast<int>(molecular_properties.size()) - 3);
    mu_molecular_properties_const_ptr.reserve(static_cast<int>(molecular_properties.size()) - 3);
    
    for (int mi = 3; mi < static_cast<int>(molecular_properties.size()); mi++)
    {
        mu_molecular_properties.push_back(*molecular_properties[mi]);
    }
    
    for (int mi = 0; mi < (static_cast<int>(molecular_properties.size()) - 3); mi++)
    {
        mu_molecular_properties_const_ptr.push_back(&mu_molecular_properties[mi]);
    }
    
    const double mu = d_equation_of_shear_viscosity->
        getShearViscosity(
            pressure,
            temperature,
            mu_molecular_properties_const_ptr);
    
    return c_p*mu/Pr;
}
