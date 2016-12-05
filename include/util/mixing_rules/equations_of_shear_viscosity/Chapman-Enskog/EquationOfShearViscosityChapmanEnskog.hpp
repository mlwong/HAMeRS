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
        
    private:
        
};

#endif /* EQUATION_OF_SHEAR_VISCOSITY_CHAPMAN_ENSKOG_HPP */
