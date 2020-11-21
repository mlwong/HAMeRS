#ifndef EQUATION_OF_SHEAR_VISCOSITY_CHAPMAN_ENSKOG_HPP
#define EQUATION_OF_SHEAR_VISCOSITY_CHAPMAN_ENSKOG_HPP

#include "util/mixing_rules/equations_of_shear_viscosity/EquationOfShearViscosity.hpp"

#include <cmath>

class EquationOfShearViscosityChapmanEnskog: public EquationOfShearViscosity
{
    public:
        EquationOfShearViscosityChapmanEnskog(
            const std::string& object_name,
            const tbox::Dimension& dim):
                EquationOfShearViscosity(
                    object_name,
                    dim)
        {}
        
        ~EquationOfShearViscosityChapmanEnskog() {}
        
        /*
         * Print all characteristics of the equation of shear viscosity class.
         */
        void
        printClassData(std::ostream& os) const;
        
        /*
         * Compute the shear viscosity.
         */
        double
        getShearViscosity(
            const double* const pressure,
            const double* const temperature,
            const std::vector<const double*>& molecular_properties) const;
        
        /*
         * Compute the shear viscosity.
         */
        void
        computeShearViscosity(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& data_shear_viscosity,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_temperature,
            const std::vector<const double*>& molecular_properties,
            const hier::Box& domain) const;
        
        /*
         * Compute the shear viscosity.
         */
        void
        computeShearViscosity(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& data_shear_viscosity,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_temperature,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_molecular_properties,
            const hier::Box& domain) const;
        
    private:
        
};

#endif /* EQUATION_OF_SHEAR_VISCOSITY_CHAPMAN_ENSKOG_HPP */
