#ifndef EQUATION_OF_STATE_MIXING_RULES_IDEAL_GAS_HPP
#define EQUATION_OF_STATE_MIXING_RULES_IDEAL_GAS_HPP

#include "util/mixing_rules/equations_of_state/EquationOfStateMixingRules.hpp"

#include "util/mixing_rules/equations_of_state/ideal_gas/EquationOfStateIdealGas.hpp"

class EquationOfStateMixingRulesIdealGas: public EquationOfStateMixingRules
{
    public:        
        EquationOfStateMixingRulesIdealGas(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const int& num_species,
            const MIXING_CLOSURE_MODEL::TYPE& mixing_closure_model,
            const boost::shared_ptr<tbox::Database>& equation_of_state_mixing_rules_db);
        
        /*
         * Print all characteristics of the equation of state mixing rules class.
         */
        void
        printClassData(std::ostream& os) const;
        
        /*
         * Put the characteristics of the equation of state mixing rules class into the restart
         * database.
         */
        void
        putToRestart(
            const boost::shared_ptr<tbox::Database>& restart_db) const;
        
        /*
         * Return the boost::shared_ptr to the equation of state.
         */
        const boost::shared_ptr<EquationOfState>&
        getEquationOfState(const int species_index = 0) const
        {
            NULL_USE(species_index);
            return d_equation_of_state;
        }
        
        /*
         * Compute the pressure of the mixture with isothermal and isobaric assumptions.
         */
        double
        getPressure(
            const double* const density,
            const double* const internal_energy,
            const std::vector<const double*>& mass_fraction) const;
        
        /*
         * Compute the pressure of the mixture with isothermal and isobaric assumptions.
         */
        void
        computePressure(
            boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fraction,
            const hier::Box& domain) const;
        
        /*
         * Compute the pressure of the mixture with isobaric assumption.
         */
        double
        getPressure(
            const double* const density,
            const double* const internal_energy,
            const std::vector<const double*>& mass_fraction,
            const std::vector<const double*>& volume_fraction) const;
        
        /*
         * Compute the pressure of the mixture with isobaric assumption.
         */
        void
        computePressure(
            boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fraction,
            const boost::shared_ptr<pdat::CellData<double> >& data_volume_fraction,
            const hier::Box& domain) const;
        
        /*
         * Compute the sound speed of the mixture with isothermal and isobaric assumptions.
         */
        double
        getSoundSpeed(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& mass_fraction) const;
        
        /*
         * Compute the sound speed of the mixture with isothermal and isobaric assumptions.
         */
        void
        computeSoundSpeed(
            boost::shared_ptr<pdat::CellData<double> >& data_sound_speed,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fraction,
            const hier::Box& domain) const;
        
        /*
         * Compute the sound speed of the mixture with isobaric assumption.
         */
        double
        getSoundSpeed(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& mass_fraction,
            const std::vector<const double*>& volume_fraction) const;
        
        /*
         * Compute the sound speed of the mixture with isobaric assumption.
         */
        void
        computeSoundSpeed(
            boost::shared_ptr<pdat::CellData<double> >& data_sound_speed,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fraction,
            const boost::shared_ptr<pdat::CellData<double> >& data_volume_fraction,
            const hier::Box& domain) const;
        
        /*
         * Compute the specific internal energy of the mixture with isothermal and isobaric assumptions.
         */
        double
        getInternalEnergy(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& mass_fraction) const;
        
        /*
         * Compute the specific internal energy of the mixture with isothermal and isobaric assumptions.
         */
        void
        computeInternalEnergy(
            boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fraction,
            const hier::Box& domain) const;
        
        /*
         * Compute the specific internal energy of the mixture with isobaric assumption.
         */
        double
        getInternalEnergy(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& mass_fraction,
            const std::vector<const double*>& volume_fraction) const;
        
        /*
         * Compute the specific internal energy of the mixture with isobaric assumption.
         */
        void
        computeInternalEnergy(
            boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fraction,
            const boost::shared_ptr<pdat::CellData<double> >& data_volume_fraction,
            const hier::Box& domain) const;
        
        /*
         * Compute the temperature of the mixture with isothermal and isobaric assumptions.
         */
        double
        getTemperature(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& mass_fraction) const;
        
        /*
         * Compute the temperature of the mixture with isothermal and isobaric assumptions.
         */
        void
        computeTemperature(
            boost::shared_ptr<pdat::CellData<double> >& data_temperature,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fraction,
            const hier::Box& domain) const;
        
        /*
         * Compute the specific internal energy of the mixture from temperature with isothermal
         * and isobaric assumptions.
         */
        double
        getInternalEnergyFromTemperature(
            const double* const density,
            const double* const temperature,
            const std::vector<const double*>& mass_fraction) const;
        
        /*
         * Compute the specific internal energy of the mixture from temperature with isothermal
         * and isobaric assumptions.
         */
        void
        computeInternalEnergyFromTemperature(
            boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fraction,
            const hier::Box& domain) const;
        
        /*
         * Compute the isochoric specific heat capacity of mixture with isothermal assumption.
         */
        double
        getIsochoricSpecificHeatCapacity(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& mass_fraction) const;
        
        /*
         * Compute the isochoric specific heat capacity of mixture with isothermal assumption.
         */
        void
        computeIsochoricSpecificHeatCapacity(
            boost::shared_ptr<pdat::CellData<double> >& data_isochoric_specific_heat_capacity,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fraction,
            const hier::Box& domain) const;
        
        /*
         * Compute the isobaric specific heat capacity of mixture with isothermal assumption.
         */
        double
        getIsobaricSpecificHeatCapacity(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& mass_fraction) const;
        
        /*
         * Compute the isobaric specific heat capacity of mixture with isothermal assumption.
         */
        void
        computeIsobaricSpecificHeatCapacity(
            boost::shared_ptr<pdat::CellData<double> >& data_isobaric_specific_heat_capacity,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fraction,
            const hier::Box& domain) const;
        
        /*
         * Compute the density of mixture with isothermal and isobaric assumptions.
         */
        double
        getMixtureDensity(
            const double* const pressure,
            const double* const temperature,
            const std::vector<const double*>& mass_fraction) const;
        
        /*
         * Compute the density of mixture with isothermal and isobaric assumptions.
         */
        void
        computeMixtureDensity(
            boost::shared_ptr<pdat::CellData<double> >& data_mixture_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fraction,
            const hier::Box& domain) const;
        
        /*
         * Get the number of thermodynamic properties of a species.
         */
        int
        getNumberOfSpeciesThermodynamicProperties(const int species_index = 0) const
        {
            NULL_USE(species_index);
            return 4;
        }
        
        /*
         * Get the thermodynamic properties of a species.
         */
        void
        getSpeciesThermodynamicProperties(
            std::vector<double*>& species_thermo_properties,
            const int species_index = 0) const;
        
    private:
        /*
         * Get the number of thermodynamic properties of the mixture.
         */
        int
        getNumberOfMixtureThermodynamicProperties() const;
        
        /*
         * Get the thermodynamic properties of the mixture.
         */
        void
        getMixtureThermodynamicProperties(
            std::vector<double*>& mixture_thermo_properties,
            const std::vector<const double*>& species_fraction) const;
        
        /*
         * Get the thermodynamic properties of the mixture.
         */
        void
        computeMixtureThermodynamicProperties(
            boost::shared_ptr<pdat::CellData<double> >& data_mixture_thermo_properties,
            const boost::shared_ptr<pdat::CellData<double> >& data_species_fraction,
            const hier::Box& domain) const;
        
        /*
         * Compute the thermodynamic properties of the mixture with mass fraction.
         */
        void
        getMixtureThermodynamicPropertiesWithMassFraction(
            std::vector<double*>& mixture_thermo_properties,
            const std::vector<const double*>& mass_fraction) const;
        
        /*
         * Compute the thermodynamic properties of the mixture with mass fraction.
         */
        void
        computeMixtureThermodynamicPropertiesWithMassFraction(
            boost::shared_ptr<pdat::CellData<double> >& data_mixture_thermo_properties,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fraction,
            const hier::Box& domain) const;
        
        /*
         * Compute the thermodynamic properties of the mixture with volume fraction.
         */
        void
        getMixtureThermodynamicPropertiesWithVolumeFraction(
            std::vector<double*>& mixture_thermo_properties,
            const std::vector<const double*>& volume_fraction) const;
        
        /*
         * Compute the thermodynamic properties of the mixture with volume fraction.
         */
        void
        computeMixtureThermodynamicPropertiesWithVolumeFraction(
            boost::shared_ptr<pdat::CellData<double> >& data_mixture_thermo_properties,
            const boost::shared_ptr<pdat::CellData<double> >& data_volume_fraction,
            const hier::Box& domain) const;
        
        /*
         * Ratio of specific heats of different species.
         */
        std::vector<double> d_species_gamma;
        
        /*
         * Gas constants of different species.
         */
        std::vector<double> d_species_R;
        
        /*
         * Specific heats of different species.
         */
        std::vector<double> d_species_c_p;
        std::vector<double> d_species_c_v;
        
        /*
         * boost::shared_ptr to EquationOfState.
         */
        boost::shared_ptr<EquationOfState> d_equation_of_state;
        
};

#endif /* EQUATION_OF_STATE_MIXING_RULES_IDEAL_GAS_HPP */
