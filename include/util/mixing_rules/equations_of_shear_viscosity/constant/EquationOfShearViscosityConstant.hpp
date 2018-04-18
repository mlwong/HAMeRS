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
        
        ~EquationOfShearViscosityConstant() {}
        
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
            boost::shared_ptr<pdat::CellData<double> >& data_shear_viscosity,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
            const std::vector<const double*>& molecular_properties,
            const hier::Box& domain) const;
        
        /*
         * Compute the shear viscosity.
         */
        void
        computeShearViscosity(
            boost::shared_ptr<pdat::CellData<double> >& data_shear_viscosity,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
            const boost::shared_ptr<pdat::CellData<double> >& data_molecular_properties,
            const hier::Box& domain) const;
        
    private:
        
};

#endif /* #ifndef EQUATION_OF_SHEAR_VISCOSITY_CONSTANT_HPP */
