#ifndef EQUATION_OF_THERMAL_CONDUCTIVITY_HPP
#define EQUATION_OF_THERMAL_CONDUCTIVITY_HPP

#include "HAMeRS_config.hpp"

#include "HAMeRS_memory.hpp"

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/tbox/Dimension.h"

#include <string>
#include <vector>

using namespace SAMRAI;

class EquationOfThermalConductivity
{
    public:
        EquationOfThermalConductivity(
            const std::string& object_name,
            const tbox::Dimension& dim):
                d_object_name(object_name),
                d_dim(dim)
        {}
        
        virtual ~EquationOfThermalConductivity() {}
        
        /*
         * Print all characteristics of the equation of thermal conductivity class.
         */
        virtual void
        printClassData(std::ostream& os) const = 0;
        
        /*
         * Compute the thermal conductivity.
         */
        virtual double
        getThermalConductivity(
            const double* const pressure,
            const double* const temperature,
            const std::vector<const double*>& molecular_properties) const = 0;
        
        /*
         * Compute the thermal conductivity.
         */
        void
        computeThermalConductivity(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& data_thermal_conductivity,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_temperature,
            const std::vector<const double*>& molecular_properties) const
        {
            const hier::Box empty_box(d_dim);
            computeThermalConductivity(
                data_thermal_conductivity,
                data_pressure,
                data_temperature,
                molecular_properties,
                empty_box);
        }
        
        /*
         * Compute the thermal conductivity.
         */
        virtual void
        computeThermalConductivity(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& data_thermal_conductivity,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_temperature,
            const std::vector<const double*>& molecular_properties,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the thermal conductivity.
         */
        void
        computeThermalConductivity(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& data_thermal_conductivity,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_temperature,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_molecular_properties) const
        {
            const hier::Box empty_box(d_dim);
            computeThermalConductivity(
                data_thermal_conductivity,
                data_pressure,
                data_temperature,
                data_molecular_properties,
                empty_box);
        }
        
        /*
         * Compute the thermal conductivity.
         */
        virtual void
        computeThermalConductivity(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& data_thermal_conductivity,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_temperature,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_molecular_properties,
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

#endif /* EQUATION_OF_THERMAL_CONDUCTIVITY_HPP */
