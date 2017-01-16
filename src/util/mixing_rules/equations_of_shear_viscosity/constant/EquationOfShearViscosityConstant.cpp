#include "util/mixing_rules/equations_of_shear_viscosity/constant/EquationOfShearViscosityConstant.hpp"

/*
 * Print all characteristics of the equation of shear viscosity class.
 */
void
EquationOfShearViscosityConstant::printClassData(
    std::ostream& os) const
{
    os << "\nPrint EquationOfShearViscosityConstant object..."
       << std::endl;
       
    os << std::endl;
    
    os << "EquationOfShearViscosityConstant: this = "
       << (EquationOfShearViscosityConstant *)this
       << std::endl;
    
    os << "d_object_name = "
       << d_object_name
       << std::endl;
}


/*
 * Compute the shear viscosity.
 */
double
EquationOfShearViscosityConstant::getShearViscosity(
    const double* const pressure,
    const double* const temperature,
    const std::vector<const double*>& molecular_properties) const
{
    NULL_USE(pressure);
    NULL_USE(temperature);
    
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(molecular_properties.size()) >= 1);
#endif
    
    return *(molecular_properties[0]);
}
