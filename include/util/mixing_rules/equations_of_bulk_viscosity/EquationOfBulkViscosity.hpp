#ifndef EQUATION_OF_BULK_VISCOSITY_HPP
#define EQUATION_OF_BULK_VISCOSITY_HPP

#include "HAMeRS_config.hpp"

#include "HAMeRS_memory.hpp"

#include "SAMRAI/pdat/CellData.h"
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
        
        virtual ~EquationOfBulkViscosity() {}
        
        /*
         * Print all characteristics of the equation of bulk viscosity class.
         */
        virtual void
        printClassData(std::ostream& os) const = 0;
        
        /*
         * Compute the bulk viscosity.
         */
        virtual Real
        getBulkViscosity(
            const Real* const pressure,
            const Real* const temperature,
            const std::vector<const Real*>& molecular_properties) const = 0;
        
        /*
         * Compute the bulk viscosity.
         */
        void
        computeBulkViscosity(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_bulk_viscosity,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_temperature,
            const std::vector<const Real*>& molecular_properties) const
        {
            const hier::Box empty_box(d_dim);
            computeBulkViscosity(
                data_bulk_viscosity,
                data_pressure,
                data_temperature,
                molecular_properties,
                empty_box);
        }
        
        /*
         * Compute the bulk viscosity.
         */
        virtual void
        computeBulkViscosity(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_bulk_viscosity,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_temperature,
            const std::vector<const Real*>& molecular_properties,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the bulk viscosity.
         */
        void
        computeBulkViscosity(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_bulk_viscosity,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_temperature,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_molecular_properties) const
        {
            const hier::Box empty_box(d_dim);
            computeBulkViscosity(
                data_bulk_viscosity,
                data_pressure,
                data_temperature,
                data_molecular_properties,
                empty_box);
        }
        
        /*
         * Compute the bulk viscosity.
         */
        virtual void
        computeBulkViscosity(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_bulk_viscosity,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_temperature,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_molecular_properties,
            const hier::Box& domain) const = 0;
        
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
