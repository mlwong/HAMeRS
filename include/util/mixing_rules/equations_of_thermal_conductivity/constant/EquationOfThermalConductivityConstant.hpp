#ifndef EQUATION_OF_THERMAL_CONDUCTIVITY_CONSTANT_HPP
#define EQUATION_OF_THERMAL_CONDUCTIVITY_CONSTANT_HPP

#include "util/mixing_rules/equations_of_thermal_conductivity/EquationOfThermalConductivity.hpp"

class EquationOfThermalConductivityConstant: public EquationOfThermalConductivity
{
    public:
        EquationOfThermalConductivityConstant(
            const std::string& object_name,
            const tbox::Dimension& dim):
                EquationOfThermalConductivity(
                    object_name,
                    dim)
        {}
        
        ~EquationOfThermalConductivityConstant() {}
        
        /*
         * Print all characteristics of the equation of thermal conductivity class.
         */
        void
        printClassData(std::ostream& os) const;
        
        /*
         * Compute the thermal conductivity.
         */
        Real
        getThermalConductivity(
            const Real* const pressure,
            const Real* const temperature,
            const std::vector<const Real*>& molecular_properties) const;
        
        /*
         * Compute the thermal conductivity.
         */
        void
        computeThermalConductivity(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_thermal_conductivity,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_temperature,
            const std::vector<const Real*>& molecular_properties,
            const hier::Box& domain) const;
        
        /*
         * Compute the thermal conductivity.
         */
        void
        computeThermalConductivity(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_thermal_conductivity,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_temperature,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_molecular_properties,
            const hier::Box& domain) const;
        
    private:
        
};

#endif /* EQUATION_OF_THERMAL_CONDUCTIVITY_CONSTANT_HPP */
