#ifndef EQUATION_OF_SHEAR_VISCOSITY_CONSTANT_HPP
#define EQUATION_OF_SHEAR_VISCOSITY_CONSTANT_HPP

#include "util/mixing_rules/equations_of_shear_viscosity/EquationOfShearViscosity.hpp"

class EquationOfShearViscosityConstant: public EquationOfShearViscosity
{
    public:
        EquationOfShearViscosityConstant(
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

#endif /* #ifndef EQUATION_OF_SHEAR_VISCOSITY_CONSTANT_HPP */