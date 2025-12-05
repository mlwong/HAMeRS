#ifndef EQUATION_OF_BULK_VISCOSITY_CONSTANT_RATIO_TO_SHEAR_VISCOSITY_HPP
#define EQUATION_OF_BULK_VISCOSITY_CONSTANT_RATIO_TO_SHEAR_VISCOSITY_HPP

#include "util/mixing_rules/equations_of_bulk_viscosity/EquationOfBulkViscosity.hpp"

#include "util/mixing_rules/equations_of_shear_viscosity/EquationOfShearViscosityMixingRulesManager.hpp"

class EquationOfBulkViscosityConstantRatioToShearViscosity: public EquationOfBulkViscosity
{
    public:
        EquationOfBulkViscosityConstantRatioToShearViscosity(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<EquationOfShearViscosity>& equation_of_shear_viscosity):
                EquationOfBulkViscosity(
                    object_name,
                    dim),
                d_equation_of_shear_viscosity(equation_of_shear_viscosity)
        {}
        
        ~EquationOfBulkViscosityConstantRatioToShearViscosity() {}
        
        /*
         * Print all characteristics of the equation of bulk viscosity class.
         */
        void
        printClassData(std::ostream& os) const;
        
        /*
         * Compute the bulk viscosity.
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
        /*
         * Boost shared pointer to equation of shear viscosity.
         */
        const HAMERS_SHARED_PTR<EquationOfShearViscosity> d_equation_of_shear_viscosity;
        
};

#endif /* EQUATION_OF_BULK_VISCOSITY_CONSTANT_RATIO_TO_SHEAR_VISCOSITY_HPP */
