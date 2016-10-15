#ifndef EQUATION_OF_THERMAL_CONDUCTIVITY_PRANDTL_HPP
#define EQUATION_OF_THERMAL_CONDUCTIVITY_PRANDTL_HPP

#include "util/mixing_rules/equations_of_thermal_conductivity/EquationOfThermalConductivity.hpp"

class EquationOfThermalConductivityPrandtl: public EquationOfThermalConductivity
{
    public:
        EquationOfThermalConductivityPrandtl(
            const std::string& object_name,
            const tbox::Dimension& dim):
                EquationOfThermalConductivity(
                    object_name,
                    dim)
        {}
        
        /*
         * Print all characteristics of the equation of thermal conductivity class.
         */
        void
        printClassData(std::ostream& os) const;
        
        /*
         * Compute the thermal conductivity.
         */
        double
        getThermalConductivity(
            const double* const pressure,
            const double* const temperature,
            const std::vector<const double*>& molecular_properties) const;
        
    private:
        
};

#endif /* EQUATION_OF_THERMAL_CONDUCTIVITY_PRANDTL_HPP */