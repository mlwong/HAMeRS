#ifndef EQUATION_OF_STATE_IDEAL_GAS_HPP
#define EQUATION_OF_STATE_IDEAL_GAS_HPP

#include "util/mixing_rules/equations_of_state/EquationOfState.hpp"

class EquationOfStateIdealGas: public EquationOfState
{
    public:        
        EquationOfStateIdealGas(
            const std::string& object_name,
            const tbox::Dimension& dim):
                EquationOfState(
                    object_name,
                    dim)
        {}
        
        /*
         * Print all characteristics of the equation of state class.
         */
        void
        printClassData(std::ostream& os) const;
        
        /*
         * Compute the pressure.
         */
        double
        getPressure(
            const double* const density,
            const double* const internal_energy,
            const std::vector<const double*>& thermo_properties) const;
        
        /*
         * Compute the sound speed.
         */
        double
        getSoundSpeed(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& thermo_properties) const;
        
        /*
         * Compute the specific internal energy.
         */
        double
        getInternalEnergy(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& thermo_properties) const;
        
        /*
         * Compute the specific enthalpy.
         */
        double
        getEnthalpy(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& thermo_properties) const;
        
        /*
         * Compute the temperature.
         */
        double
        getTemperature(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& thermo_properties) const;
        
        /*
         * Compute the partial derivative of internal energy w.r.t. pressure under constant density.
         */
        double
        getIsochoricPartialInternalEnergyPartialPressure(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& thermo_properties) const;
        
        /*
         * Compute the partial derivative of internal energy w.r.t. density under constant pressure.
         */
        double
        getIsobaricPartialInternalEnergyPartialDensity(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& thermo_properties) const;
        
    private:
        
};

#endif /* EQUATION_OF_STATE_IDEAL_GAS_HPP */
