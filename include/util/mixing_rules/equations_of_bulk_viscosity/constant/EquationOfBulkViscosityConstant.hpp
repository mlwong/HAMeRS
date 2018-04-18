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
        double
        getBulkViscosity(
            const double* const pressure,
            const double* const temperature,
            const std::vector<const double*>& molecular_properties) const;
        
        /*
         * Compute the bulk viscosity.
         */
        void
        computeBulkViscosity(
            boost::shared_ptr<pdat::CellData<double> >& data_bulk_viscosity,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
            const std::vector<const double*>& molecular_properties,
            const hier::Box& domain) const;
        
        /*
         * Compute the bulk viscosity.
         */
        void
        computeBulkViscosity(
            boost::shared_ptr<pdat::CellData<double> >& data_bulk_viscosity,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
            const boost::shared_ptr<pdat::CellData<double> >& data_molecular_properties,
            const hier::Box& domain) const;
        
    private:
        
};

#endif /* EQUATION_OF_BULK_VISCOSITY_CONSTANT_HPP */
