#include "util/mixing_rules/equations_of_shear_viscosity/Chapman-Enskog/EquationOfShearViscosityChapmanEnskog.hpp"

#include <cmath>

/*
 * Print all characteristics of the equation of shear viscosity class.
 */
void
EquationOfShearViscosityChapmanEnskog::printClassData(
    std::ostream& os) const
{
    os << "\nPrint EquationOfShearViscosityChapmanEnskog object..."
       << std::endl;
       
    os << std::endl;
    
    os << "EquationOfShearViscosityChapmanEnskog: this = "
       << (EquationOfShearViscosityChapmanEnskog *)this
       << std::endl;
    
    os << "d_object_name = "
       << d_object_name
       << std::endl;
}


/*
 * Compute the shear viscosity.
 */
double
EquationOfShearViscosityChapmanEnskog::getShearViscosity(
    const double* const pressure,
    const double* const temperature,
    const std::vector<const double*>& molecular_properties) const
{
    NULL_USE(pressure);
    
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(molecular_properties.size()) >= 3);
#endif
    
    double mu = 0.0;
    
    const double& epsilon_by_k = *(molecular_properties[0]);
    const double& sigma = *(molecular_properties[1]);
    const double& M = *(molecular_properties[2]);
    
    const double& T = *temperature;
    
    const double A = 1.16145;
    const double B = -0.14874;
    const double C = 0.52487;
    const double D = -0.7732;
    const double E = 2.16178;
    const double F = -2.43787;
    
    const double T_star = T/epsilon_by_k;
    const double Omega = A*pow(T_star, B) + C*exp(D*T_star) + E*exp(F*T_star);
    
    mu = 2.6693e-6*sqrt(M*T)/(Omega*sigma*sigma);
    
    return mu;
}
