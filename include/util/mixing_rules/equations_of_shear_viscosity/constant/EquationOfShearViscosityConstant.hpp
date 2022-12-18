#ifndef EQUATION_OF_SHEAR_VISCOSITY_CONSTANT_HPP
#define EQUATION_OF_SHEAR_VISCOSITY_CONSTANT_HPP

#include "util/mixing_rules/equations_of_shear_viscosity/EquationOfShearViscosity.hpp"

class EquationOfShearViscosityConstant: public EquationOfShearViscosity
{
    public:
        EquationOfShearViscosityConstant(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const bool use_constant_kinematic_viscosity_and_ideal_gas_assumptions,
            const double R_u = double(8.314462618153240)): // Universal gas constant in SI units.
                EquationOfShearViscosity(
                    object_name,
                    dim),
                d_use_constant_kinematic_viscosity_and_ideal_gas_assumptions(use_constant_kinematic_viscosity_and_ideal_gas_assumptions),
                d_R_u(R_u)
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
        /*
         * Whether to assume the kinematic viscosity (instead of dynamic viscosity) is constant and the species is ideal gas.
         */
        bool d_use_constant_kinematic_viscosity_and_ideal_gas_assumptions;
        
        /*
         * Universal gas constant.
         */
        double d_R_u;
        
};

#endif /* #ifndef EQUATION_OF_SHEAR_VISCOSITY_CONSTANT_HPP */
