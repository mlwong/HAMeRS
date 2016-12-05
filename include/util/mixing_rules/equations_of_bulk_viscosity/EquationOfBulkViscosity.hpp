#ifndef EQUATION_OF_BULK_VISCOSITY_HPP
#define EQUATION_OF_BULK_VISCOSITY_HPP

#include "SAMRAI/tbox/Dimension.h"

#include <string>
#include <vector>

using namespace SAMRAI;

class EquationOfBulkViscosity
{
    public:
        EquationOfBulkViscosity(
            const std::string& object_name,
            const tbox::Dimension& dim):
                d_object_name(object_name),
                d_dim(dim)
        {}
        
        /*
         * Print all characteristics of the equation of bulk viscosity class.
         */
        virtual void
        printClassData(std::ostream& os) const = 0;
        
        /*
         * Compute the bulk viscosity.
         */
        virtual double
        getBulkViscosity(
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

#endif /* EQUATION_OF_BULK_VISCOSITY_HPP */
