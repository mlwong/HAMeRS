#ifndef EQUATION_OF_STATE_HPP
#define EQUATION_OF_STATE_HPP

#include "SAMRAI/tbox/Dimension.h"

#include <string>
#include <vector>

using namespace SAMRAI;

class EquationOfState
{
    public:
        EquationOfState(
            const std::string& object_name,
            const tbox::Dimension& dim):
                d_object_name(object_name),
                d_dim(dim)
        {}
        
        /*
         * Print all characteristics of the equation of state class.
         */
        virtual void
        printClassData(std::ostream& os) const = 0;
        
        /*
         * Compute the pressure.
         */
        virtual double
        getPressure(
            const double* const density,
            const double* const internal_energy,
            const std::vector<const double*>& thermo_properties) const = 0;
        
        /*
         * Compute the sound speed.
         */
        virtual double
        getSoundSpeed(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& thermo_properties) const = 0;
        
        /*
         * Compute the specific internal energy.
         */
        virtual double
        getInternalEnergy(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& thermo_properties) const = 0;
        
        /*
         * Compute the specific enthalpy.
         */
        virtual double
        getEnthalpy(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& thermo_properties) const = 0;
        
        /*
         * Compute the temperature.
         */
        virtual double
        getTemperature(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& thermo_properties) const = 0;
        
        /*
         * Compute the specific internal energy from temperature.
         */
        virtual double
        getInternalEnergyFromTemperature(
            const double* const density,
            const double* const temperature,
            const std::vector<const double*>& thermo_properties) const = 0;
        
        /*
         * Compute the isochoric specific heat capacity.
         */
        virtual double
        getIsochoricSpecificHeatCapacity(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& thermo_properties) const = 0;
        
        /*
         * Compute the isobaric specific heat capacity.
         */
        virtual double
        getIsobaricSpecificHeatCapacity(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& thermo_properties) const = 0;
        
        /*
         * Compute the partial derivative of internal energy w.r.t. pressure under constant density.
         */
        virtual double
        getIsochoricPartialInternalEnergyPartialPressure(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& thermo_properties) const = 0;
        
        /*
         * Compute the partial derivative of internal energy w.r.t. density under constant pressure.
         */
        virtual double
        getIsobaricPartialInternalEnergyPartialDensity(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& thermo_properties) const = 0;
        
        /*
         * Compute the density.
         */
        virtual double
        getDensity(
            const double* const pressure,
            const double* const temperature,
            const std::vector<const double*>& thermo_properties) const = 0;
        
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

#endif /* EQUATION_OF_STATE_HPP */
