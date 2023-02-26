#ifndef EQUATION_OF_BULK_VISCOSITY_CONSTANT_HPP
#define EQUATION_OF_BULK_VISCOSITY_CONSTANT_HPP

#include "util/mixing_rules/equations_of_bulk_viscosity/EquationOfBulkViscosity.hpp"

class EquationOfBulkViscosityConstant: public EquationOfBulkViscosity
{
    public:
        EquationOfBulkViscosityConstant(
            const std::string& object_name,
            const tbox::Dimension& dim):
                EquationOfBulkViscosity(
                    object_name,
                    dim)
        {}
        
        ~EquationOfBulkViscosityConstant() {}
        
        /*
         * Print all characteristics of the equation of thermal conductivity class.
         */
        void
        printClassData(std::ostream& os) const;
        
        /*
         * Compute the thermal conductivity.
         */
        Real
        getBulkViscosity(
            const Real* const pressure,
            const Real* const temperature,
            const std::vector<const Real*>& molecular_properties) const;
        
        /*
         * Compute the bulk viscosity.
         */
        void
        computeBulkViscosity(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_bulk_viscosity,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_temperature,
            const std::vector<const Real*>& molecular_properties,
            const hier::Box& domain) const;
        
        /*
         * Compute the bulk viscosity.
         */
        void
        computeBulkViscosity(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_bulk_viscosity,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_temperature,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_molecular_properties,
            const hier::Box& domain) const;
        
    private:
        
};

#endif /* EQUATION_OF_BULK_VISCOSITY_CONSTANT_HPP */
