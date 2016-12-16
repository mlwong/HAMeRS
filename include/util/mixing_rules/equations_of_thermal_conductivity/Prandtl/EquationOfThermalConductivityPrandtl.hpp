#ifndef EQUATION_OF_THERMAL_CONDUCTIVITY_PRANDTL_HPP
#define EQUATION_OF_THERMAL_CONDUCTIVITY_PRANDTL_HPP

#include "util/mixing_rules/equations_of_thermal_conductivity/EquationOfThermalConductivity.hpp"

#include "util/mixing_rules/equations_of_shear_viscosity/EquationOfShearViscosityMixingRulesManager.hpp"

class EquationOfThermalConductivityPrandtl: public EquationOfThermalConductivity
{
    public:
        EquationOfThermalConductivityPrandtl(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const boost::shared_ptr<EquationOfShearViscosity>& equation_of_shear_viscosity):
                EquationOfThermalConductivity(
                    object_name,
                    dim),
                d_equation_of_shear_viscosity(equation_of_shear_viscosity)
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
        /*
         * Boost shared pointer to equation of shear viscosity.
         */
        const boost::shared_ptr<EquationOfShearViscosity> d_equation_of_shear_viscosity;
        
};

#endif /* EQUATION_OF_THERMAL_CONDUCTIVITY_PRANDTL_HPP */