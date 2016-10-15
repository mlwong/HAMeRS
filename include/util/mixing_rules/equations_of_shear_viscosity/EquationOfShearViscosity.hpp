#ifndef EQUATION_OF_SHEAR_VISCOSITY_HPP
#define EQUATION_OF_SHEAR_VISCOSITY_HPP

#include "SAMRAI/tbox/Dimension.h"

#include <string>
#include <vector>

using namespace SAMRAI;

class EquationOfShearViscosity
{
    public:
        EquationOfShearViscosity(
            const std::string& object_name,
            const tbox::Dimension& dim):
                d_object_name(object_name),
                d_dim(dim)
        {}
        
        /*
         * Print all characteristics of the equation of shear viscosity class.
         */
        virtual void
        printClassData(std::ostream& os) const = 0;
        
        /*
         * Compute the shear viscosity.
         */
        virtual double
        getShearViscosity(
            const double* const pressure,
            const double* const temperature,
            const std::vector<const double*>& molecular_properties) const = 0;
        
    protected:
        /*
         * The object name is used for error/warning reporting.
         */
        const std::string d_object_name;
        
        /*
         * Problem dimension.
         */
        const tbox::Dimension d_dim;
        
};

#endif /* EQUATION_OF_SHEAR_VISCOSITY_HPP */