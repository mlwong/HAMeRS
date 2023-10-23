#ifndef EQUATION_OF_SHEAR_VISCOSITY_POWER_LAW_HPP
#define EQUATION_OF_SHEAR_VISCOSITY_POWER_LAW_HPP

#include "util/mixing_rules/equations_of_shear_viscosity/EquationOfShearViscosity.hpp"

class EquationOfShearViscosityPowerLaw: public EquationOfShearViscosity
{
    public:
        EquationOfShearViscosityPowerLaw(
            const std::string& object_name,
            const tbox::Dimension& dim):
                EquationOfShearViscosity(
                    object_name,
                    dim)
        {}
        
        ~EquationOfShearViscosityPowerLaw() {}
        
        /*
         * Print all characteristics of the equation of shear viscosity class.
         */
        void
        printClassData(std::ostream& os) const;
        
        /*
         * Compute the shear viscosity.
         */
        Real
        getShearViscosity(
            const Real* const pressure,
            const Real* const temperature,
            const std::vector<const Real*>& molecular_properties) const;
        
        /*
         * Compute the shear viscosity.
         */
        void
        computeShearViscosity(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_shear_viscosity,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_temperature,
            const std::vector<const Real*>& molecular_properties,
            const hier::Box& domain) const;
        
        /*
         * Compute the shear viscosity.
         */
        void
        computeShearViscosity(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_shear_viscosity,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_temperature,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_molecular_properties,
            const hier::Box& domain) const;
        
    private:
        
};

#endif /* EQUATION_OF_SHEAR_VISCOSITY_POWER_LAW_HPP */
