#ifndef EQUATION_OF_STATE_IDEAL_GAS_HPP
#define EQUATION_OF_STATE_IDEAL_GAS_HPP

#include "util/equations_of_state/EquationOfState.hpp"

class EquationOfStateIdealGas: public EquationOfState
{
    public:        
        EquationOfStateIdealGas(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const int& num_species,
            const boost::shared_ptr<tbox::Database>& equation_of_state_db);
        
        /*
         * Print all characteristics of the equation of state class.
         */
        void
        printClassData(std::ostream& os) const;
        
        /*
         * Put the characteristics of the class into the restart database.
         */
        void
        putToRestart(
            const boost::shared_ptr<tbox::Database>& restart_db) const;
        
        /*
         * Compute the pressure for single-species flow.
         */
        double
        getPressure(
            const double* const density,
            const std::vector<const double*>& momentum,
            const double* const total_energy);
        
        /*
         * Compute the pressure for multi-species flow with total density and
         * mass fractions.
         */
        double
        getPressureWithMassFraction(
            const double* const density,
            const std::vector<const double*>& momentum,
            const double* const total_energy,
            const std::vector<const double*>& mass_fraction);
        
        /*
         * Compute the pressure for multi-species flow with partial densities and
         * mass fractions.
         */
        double
        getPressureWithMassFraction(
            const std::vector<const double*>& partial_density,
            const std::vector<const double*>& momentum,
            const double* const total_energy,
            const std::vector<const double*>& mass_fraction);
        
        /*
         * Compute the pressure for multi-species flow with total density and
         * volume fractions.
         */
        double
        getPressureWithVolumeFraction(
            const double* const density,
            const std::vector<const double*>& momentum,
            const double* const total_energy,
            const std::vector<const double*>& volume_fraction);
        
        /*
         * Compute the pressure for multi-species flow with partial densities and
         * volume fractions.
         */
        double
        getPressureWithVolumeFraction(
            const std::vector<const double*>& partial_density,
            const std::vector<const double*>& momentum,
            const double* const total_energy,
            const std::vector<const double*>& volume_fraction);
        
        /*
         * Compute the sound speed for single-species flow.
         */
        double
        getSoundSpeed(
            const double* const density,
            const std::vector<const double*>& momentum,
            const double* const total_energy);
        
        /*
         * Compute the sound speed for single-species flow with pressure.
         * 
         * Recommended to use with getPressure() if both pressure and sound speed
         * data is required.
         */
        double
        getSoundSpeedWithPressure(
            const double* const density,
            const double* const pressure);
        
        /*
         * Compute the sound speed for multi-species flow with total density and
         * mass fractions.
         */
        double
        getSoundSpeedWithMassFraction(
            const double* const density,
            const std::vector<const double*>& momentum,
            const double* const total_energy,
            const std::vector<const double*>& mass_fraction);
        
        /*
         * Compute the sound speed for multi-species flow with partial densities and
         * mass fractions.
         */
        double
        getSoundSpeedWithMassFraction(
            const std::vector<const double*>& partial_density,
            const std::vector<const double*>& momentum,
            const double* const total_energy,
            const std::vector<const double*>& mass_fraction);
        
        /*
         * Compute the sound speed for multi-species flow with total density,
         * mass fractions and pressure.
         * 
         * Recommended to use with getPressureWithMassFraction() if both pressure
         * and sound speed data is required.
         */
        double
        getSoundSpeedWithMassFractionAndPressure(
            const double* const density,
            const std::vector<const double*>& mass_fraction,
            const double* const pressure);
        
        /*
         * Compute the sound speed for multi-species flow with partial densities,
         * mass fractions and pressure.
         * 
         * Recommended to use with getPressureWithMassFraction() if both pressure
         * and sound speed data is required.
         */
        double
        getSoundSpeedWithMassFractionAndPressure(
            const std::vector<const double*>& partial_density,
            const std::vector<const double*>& mass_fraction,
            const double* const pressure);
        
        /*
         * Compute the sound speed for multi-species flow with total density and
         * volume fractions.
         */
        double
        getSoundSpeedWithVolumeFraction(
            const double* const density,
            const std::vector<const double*>& momentum,
            const double* const total_energy,
            const std::vector<const double*>& volume_fraction);
        
        /*
         * Compute the sound speed for multi-species flow with partial densities and
         * volume fractions.
         */
        double
        getSoundSpeedWithVolumeFraction(
            const std::vector<const double*>& partial_density,
            const std::vector<const double*>& momentum,
            const double* const total_energy,
            const std::vector<const double*>& volume_fraction);
        
        /*
         * Compute the sound speed for multi-species flow with total density,
         * volume fractions and pressure.
         * 
         * Recommended to use with getPressureWithVolumeFraction() if both pressure
         * and sound speed data is required.
         */
        double
        getSoundSpeedWithVolumeFractionAndPressure(
            const double* const density,
            const std::vector<const double*>& volume_fraction,
            const double* const pressure);
        
        /*
         * Compute the sound speed for multi-species flow with partial densities,
         * volume fractions and pressure.
         * Recommended to use with getPressureWithVolumeFraction() if both pressure
         * and sound speed data is required.
         */
        double
        getSoundSpeedWithVolumeFractionAndPressure(
            const std::vector<const double*>& partial_density,
            const std::vector<const double*>& volume_fraction,
            const double* const pressure);
        
        /*
         * Compute the total energy for single-species flow.
         */
        double
        getTotalEnergy(
            const double* const density,
            const std::vector<const double*>& velocity,
            const double* const pressure);
        
        /*
         * Compute the total energy for multi-species flow with total density and
         * mass fractions.
         */
        double
        getTotalEnergyWithMassFraction(
            const double* const density,
            const std::vector<const double*>& velocity,
            const double* const pressure,
            const std::vector<const double*>& mass_fraction);
        
        /*
         * Compute the total energy for multi-species flow with partial densities and
         * mass fractions.
         */
        double
        getTotalEnergyWithMassFraction(
           const std::vector<const double*>& partial_density,
           const std::vector<const double*>& velocity,
           const double* const pressure,
           const std::vector<const double*>& mass_fraction);
        
        /*
         * Compute the total energy for multi-species flow with total density and
         * volume fractions.
         */
        double
        getTotalEnergyWithVolumeFraction(
           const double* const density,
           const std::vector<const double*>& velocity,
           const double* const pressure,
           const std::vector<const double*>& volume_fraction);
        
        /*
         * Compute the total energy for multi-species flow with partial densities and
         * volume fractions.
         */
        double
        getTotalEnergyWithVolumeFraction(
           const std::vector<const double*>& partial_density,
           const std::vector<const double*>& velocity,
           const double* const pressure,
           const std::vector<const double*>& volume_fraction);
        
        /*
         * Compute a thermodynamic property of a particular species
         * with mass fractions.
         */
        double
        getMixtureThermodynamicPropertyWithMassFraction(
            const std::string& property_name,
            const std::vector<const double*>& mass_fraction);
        
        /*
         * Compute a thermodynamic property of a particular species
         * with volume fractions.
         */
        double
        getMixtureThermodynamicPropertyWithVolumeFraction(
            const std::string& property_name,
            const std::vector<const double*>& volume_fraction);
        
        /*
         * Get a thermodynamic property of a particular species.
         */
        double&
        getSpeciesThermodynamicProperty(
            const std::string& property_name,
            const int species_index);
        
        /*
         * Determine whether a thermodynamic property is registered.
         */
        bool
        hasThermodynamicProperty(
            const std::string& property_name);
        
    private:
        /*
         * Compute the mixture's gamma for multi-species flow with partial densities.
         */
        double
        getMixtureGammaWithPartialDensity(
            const std::vector<const double*>& partial_density);
        
        /*
         * Compute the mixture's gamma for multi-species flow with mass fractions.
         */
        double
        getMixtureGammaWithMassFraction(
            const std::vector<const double*>& mass_fraction);
        
        /*
         * Compute the mixture's gamma for multi-species flow with volume fractions.
         */
        double
        getMixtureGammaWithVolumeFraction(
            const std::vector<const double*>& volume_fraction);
        
        /*
         * Ratio of specific heats of different species.
         */
        std::vector<double> d_species_gamma;
        
        /*
         * Specific heats at constant pressure of different species.
         */
        std::vector<double> d_species_Cp;
        
        /*
         * Specific heats at constant volume of different species.
         */
        std::vector<double> d_species_Cv;
        
};

#endif /* EQUATION_OF_STATE_IDEAL_GAS_HPP */
