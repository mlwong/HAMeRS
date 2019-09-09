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
        
        ~EquationOfStateMixingRulesIdealGas() {}
        
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
         * Compute the pressure of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        double
        getPressure(
            const double* const density,
            const double* const internal_energy,
            const std::vector<const double*>& mass_fractions) const;
        
        /*
         * Compute the pressure of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        void
        computePressure(
            boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
            const hier::Box& domain) const;
        
        /*
         * Compute the pressure of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        void
        computePressure(
            boost::shared_ptr<pdat::SideData<double> >& data_pressure,
            const boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the pressure of the mixture with isobaric equilibrium assumption.
         */
        double
        getPressure(
            const double* const density,
            const double* const internal_energy,
            const std::vector<const double*>& mass_fractions,
            const std::vector<const double*>& volume_fractions) const;
        
        /*
         * Compute the pressure of the mixture with isobaric equilibrium assumption.
         */
        void
        computePressure(
            boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
            const boost::shared_ptr<pdat::CellData<double> >& data_volume_fractions,
            const hier::Box& domain) const;
        
        /*
         * Compute the pressure of the mixture with isobaric equilibrium assumption.
         */
        void
        computePressure(
            boost::shared_ptr<pdat::SideData<double> >& data_pressure,
            const boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
            const boost::shared_ptr<pdat::SideData<double> >& data_volume_fractions,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the sound speed of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        double
        getSoundSpeed(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& mass_fractions) const;
        
        /*
         * Compute the sound speed of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        void
        computeSoundSpeed(
            boost::shared_ptr<pdat::CellData<double> >& data_sound_speed,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
            const hier::Box& domain) const;
        
        /*
         * Compute the sound speed of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        void
        computeSoundSpeed(
            boost::shared_ptr<pdat::SideData<double> >& data_sound_speed,
            const boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
            const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the sound speed of the mixture with isobaric equilibrium assumption.
         */
        double
        getSoundSpeed(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& mass_fractions,
            const std::vector<const double*>& volume_fractions) const;
        
        /*
         * Compute the sound speed of the mixture with isobaric equilibrium assumption.
         */
        void
        computeSoundSpeed(
            boost::shared_ptr<pdat::CellData<double> >& data_sound_speed,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
            const boost::shared_ptr<pdat::CellData<double> >& data_volume_fractions,
            const hier::Box& domain) const;
        
        /*
         * Compute the sound speed of the mixture with isobaric equilibrium assumption.
         */
        void
        computeSoundSpeed(
            boost::shared_ptr<pdat::SideData<double> >& data_sound_speed,
            const boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
            const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
            const boost::shared_ptr<pdat::SideData<double> >& data_volume_fractions,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the specific internal energy of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        double
        getInternalEnergy(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& mass_fractions) const;
        
        /*
         * Compute the specific internal energy of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        void
        computeInternalEnergy(
            boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
            const hier::Box& domain) const;
        
        /*
         * Compute the specific internal energy of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        void
        computeInternalEnergy(
            boost::shared_ptr<pdat::SideData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
            const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the specific internal energy of the mixture with isobaric equilibrium assumption.
         */
        double
        getInternalEnergy(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& mass_fractions,
            const std::vector<const double*>& volume_fractions) const;
        
        /*
         * Compute the specific internal energy of the mixture with isobaric equilibrium assumption.
         */
        void
        computeInternalEnergy(
            boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
            const boost::shared_ptr<pdat::CellData<double> >& data_volume_fractions,
            const hier::Box& domain) const;
        
        /*
         * Compute the specific internal energy of the mixture with isobaric equilibrium assumption.
         */
        void
        computeInternalEnergy(
            boost::shared_ptr<pdat::SideData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
            const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
            const boost::shared_ptr<pdat::SideData<double> >& data_volume_fractions,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the temperature of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        double
        getTemperature(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& mass_fractions) const;
        
        /*
         * Compute the temperature of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        void
        computeTemperature(
            boost::shared_ptr<pdat::CellData<double> >& data_temperature,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
            const hier::Box& domain) const;
        
        /*
         * Compute the temperature of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        void
        computeTemperature(
            boost::shared_ptr<pdat::SideData<double> >& data_temperature,
            const boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
            const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the specific internal energy of the mixture from temperature with isothermal
         * and isobaric equilibrium assumptions.
         */
        double
        getInternalEnergyFromTemperature(
            const double* const density,
            const double* const temperature,
            const std::vector<const double*>& mass_fractions) const;
        
        /*
         * Compute the specific internal energy of the mixture from temperature with isothermal
         * and isobaric equilibrium assumptions.
         */
        void
        computeInternalEnergyFromTemperature(
            boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
            const hier::Box& domain) const;
        
        /*
         * Compute the specific internal energy of the mixture from temperature with isothermal
         * and isobaric equilibrium assumptions.
         */
        void
        computeInternalEnergyFromTemperature(
            boost::shared_ptr<pdat::SideData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_temperature,
            const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the isochoric specific heat capacity of mixture with isothermal and isobaric
         * equilibrium assumptions.
         */
        double
        getIsochoricSpecificHeatCapacity(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& mass_fractions) const;
        
        /*
         * Compute the isochoric specific heat capacity of mixture with isothermal and isobaric
         * equilibrium assumptions.
         */
        void
        computeIsochoricSpecificHeatCapacity(
            boost::shared_ptr<pdat::CellData<double> >& data_isochoric_specific_heat_capacity,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
            const hier::Box& domain) const;
        
        /*
         * Compute the isochoric specific heat capacity of mixture with isothermal and isobaric
         * equilibrium assumptions.
         */
        void
        computeIsochoricSpecificHeatCapacity(
            boost::shared_ptr<pdat::SideData<double> >& data_isochoric_specific_heat_capacity,
            const boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
            const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the isobaric specific heat capacity of mixture with isothermal and isobaric
         * equilibrium assumptions.
         */
        double
        getIsobaricSpecificHeatCapacity(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& mass_fractions) const;
        
        /*
         * Compute the isobaric specific heat capacity of mixture with isothermal and isobaric
         * equilibrium assumptions.
         */
        void
        computeIsobaricSpecificHeatCapacity(
            boost::shared_ptr<pdat::CellData<double> >& data_isobaric_specific_heat_capacity,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
            const hier::Box& domain) const;
        
        /*
         * Compute the isobaric specific heat capacity of mixture with isothermal and isobaric
         * equilibrium assumptions.
         */
        void
        computeIsobaricSpecificHeatCapacity(
            boost::shared_ptr<pdat::SideData<double> >& data_isobaric_specific_heat_capacity,
            const boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
            const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the Gruneisen parameter of the mixture with isothermal and isobaric equilibrium assumptions
         * (partial derivative of pressure w.r.t. specific internal energy under constant partial densities
         * divided by mixture density).
         */
        double
        getGruneisenParameter(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& mass_fractions) const;
        
        /*
         * Compute the Gruneisen parameter of the mixture with isothermal and isobaric equilibrium assumptions
         * (partial derivative of pressure w.r.t. specific internal energy under constant partial densities
         * divided by mixture density).
         */
        void
        computeGruneisenParameter(
            boost::shared_ptr<pdat::CellData<double> >& data_gruneisen_parameter,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
            const hier::Box& domain) const;
        
        /*
         * Compute the Gruneisen parameter of the mixture with isothermal and isobaric equilibrium assumptions
         * (partial derivative of pressure w.r.t. specific internal energy under constant partial densities
         * divided by mixture density).
         */
        void
        computeGruneisenParameter(
            boost::shared_ptr<pdat::SideData<double> >& data_gruneisen_parameter,
            const boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
            const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the Gruneisen parameter of the mixture with isobaric equilibrium assumption
         * (partial derivative of pressure w.r.t. specific internal energy under constant partial densities
         * and volume fractions divided by mixture density).
         */
        double
        getGruneisenParameter(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& mass_fractions,
            const std::vector<const double*>& volume_fractions) const;
        
        /*
         * Compute the Gruneisen parameter of the mixture with isobaric equilibrium assumption
         * (partial derivative of pressure w.r.t. specific internal energy under constant partial densities
         * and volume fractions divided by mixture density).
         */
        void
        computeGruneisenParameter(
            boost::shared_ptr<pdat::CellData<double> >& data_gruneisen_parameter,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
            const boost::shared_ptr<pdat::CellData<double> >& data_volume_fractions,
            const hier::Box& domain) const;
        
        /*
         * Compute the Gruneisen parameter of the mixture with isobaric equilibrium assumption
         * (partial derivative of pressure w.r.t. specific internal energy under constant partial densities
         * and volume fractions divided by mixture density).
         */
        void
        computeGruneisenParameter(
            boost::shared_ptr<pdat::SideData<double> >& data_gruneisen_parameter,
            const boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
            const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
            const boost::shared_ptr<pdat::SideData<double> >& data_volume_fractions,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the mixture partial derivative of pressure w.r.t. partial densities under constant specific
         * internal energy with isothermal and isobaric equilibrium assumptions.
         */
        std::vector<double>
        getPressureDerivativeWithPartialDensities(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& mass_fractions) const;
        
        /*
         * Compute the mixture partial derivative of pressure w.r.t. partial densities under constant specific
         * internal energy with isothermal and isobaric equilibrium assumptions.
         */
        void
        computePressureDerivativeWithPartialDensities(
            boost::shared_ptr<pdat::CellData<double> >& data_partial_pressure_partial_partial_densities,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
            const hier::Box& domain) const;
        
        /*
         * Compute the mixture partial derivative of pressure w.r.t. partial densities under constant specific
         * internal energy with isothermal and isobaric equilibrium assumptions.
         */
        void
        computePressureDerivativeWithPartialDensities(
            boost::shared_ptr<pdat::SideData<double> >& data_partial_pressure_partial_partial_densities,
            const boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
            const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the mixture partial derivative of pressure w.r.t. partial densities under constant specific
         * internal energy and volume fractions with isobaric equilibrium assumption.
         */
        std::vector<double>
        getPressureDerivativeWithPartialDensities(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& mass_fractions,
            const std::vector<const double*>& volume_fractions) const;
        
        /*
         * Compute the mixture partial derivative of pressure w.r.t. partial densities under constant specific
         * internal energy and volume fractions with isobaric equilibrium assumption.
         */
        void
        computePressureDerivativeWithPartialDensities(
            boost::shared_ptr<pdat::CellData<double> >& data_partial_pressure_partial_partial_densities,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
            const boost::shared_ptr<pdat::CellData<double> >& data_volume_fractions,
            const hier::Box& domain) const;
        
        /*
         * Compute the mixture partial derivative of pressure w.r.t. partial densities under constant specific
         * internal energy and volume fractions with isobaric equilibrium assumption.
         */
        void
        computePressureDerivativeWithPartialDensities(
            boost::shared_ptr<pdat::SideData<double> >& data_partial_pressure_partial_partial_densities,
            const boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
            const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
            const boost::shared_ptr<pdat::SideData<double> >& data_volume_fractions,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the mixture partial derivative of pressure w.r.t. volume fractions under constant specific
         * internal energy and partial densities with isobaric equilibrium assumption.
         */
        std::vector<double>
        getPressureDerivativeWithVolumeFractions(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& mass_fractions,
            const std::vector<const double*>& volume_fractions) const;
        
        /*
         * Compute the mixture partial derivative of pressure w.r.t. volume fractions under constant specific
         * internal energy and partial densities with isobaric equilibrium assumption.
         */
        void
        computePressureDerivativeWithVolumeFractions(
            boost::shared_ptr<pdat::CellData<double> >& data_partial_pressure_partial_volume_fractions,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
            const boost::shared_ptr<pdat::CellData<double> >& data_volume_fractions,
            const hier::Box& domain) const;
        
        /*
         * Compute the mixture partial derivative of pressure w.r.t. volume fractions under constant specific
         * internal energy and partial densities with isobaric equilibrium assumption.
         */
        void
        computePressureDerivativeWithVolumeFractions(
            boost::shared_ptr<pdat::SideData<double> >& data_partial_pressure_partial_volume_fractions,
            const boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
            const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
            const boost::shared_ptr<pdat::SideData<double> >& data_volume_fractions,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the density of mixture with isothermal and isobaric equilibrium assumptions.
         */
        double
        getMixtureDensity(
            const double* const pressure,
            const double* const temperature,
            const std::vector<const double*>& mass_fractions) const;
        
        /*
         * Compute the density of mixture with isothermal and isobaric equilibrium assumptions.
         */
        void
        computeMixtureDensity(
            boost::shared_ptr<pdat::CellData<double> >& data_mixture_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
            const hier::Box& domain) const;
        
        /*
         * Compute the density of mixture with isothermal and isobaric equilibrium assumptions.
         */
        void
        computeMixtureDensity(
            boost::shared_ptr<pdat::SideData<double> >& data_mixture_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
            const boost::shared_ptr<pdat::SideData<double> >& data_temperature,
            const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
            int side_normal,
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
         * Get the thermodynamic properties of the mixture.
         */
        void
        computeMixtureThermodynamicProperties(
            boost::shared_ptr<pdat::SideData<double> >& data_mixture_thermo_properties,
            const boost::shared_ptr<pdat::SideData<double> >& data_species_fraction,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the thermodynamic properties of the mixture with mass fractions.
         */
        void
        getMixtureThermodynamicPropertiesWithMassFraction(
            std::vector<double*>& mixture_thermo_properties,
            const std::vector<const double*>& mass_fractions) const;
        
        /*
         * Compute the thermodynamic properties of the mixture with mass fractions.
         */
        void
        computeMixtureThermodynamicPropertiesWithMassFraction(
            boost::shared_ptr<pdat::CellData<double> >& data_mixture_thermo_properties,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
            const hier::Box& domain) const;
        
        /*
         * Compute the thermodynamic properties of the mixture with mass fractions.
         */
        void
        computeMixtureThermodynamicPropertiesWithMassFraction(
            boost::shared_ptr<pdat::SideData<double> >& data_mixture_thermo_properties,
            const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the thermodynamic properties of the mixture with volume fractions.
         */
        void
        getMixtureThermodynamicPropertiesWithVolumeFraction(
            std::vector<double*>& mixture_thermo_properties,
            const std::vector<const double*>& volume_fractions) const;
        
        /*
         * Compute the thermodynamic properties of the mixture with volume fractions.
         */
        void
        computeMixtureThermodynamicPropertiesWithVolumeFraction(
            boost::shared_ptr<pdat::CellData<double> >& data_mixture_thermo_properties,
            const boost::shared_ptr<pdat::CellData<double> >& data_volume_fractions,
            const hier::Box& domain) const;
        
        /*
         * Compute the thermodynamic properties of the mixture with volume fractions.
         */
        void
        computeMixtureThermodynamicPropertiesWithVolumeFraction(
            boost::shared_ptr<pdat::SideData<double> >& data_mixture_thermo_properties,
            const boost::shared_ptr<pdat::SideData<double> >& data_volume_fractions,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the isochoric specific heat capacity of mixture with isothermal and isobaric
         * equilibrium assumptions.
         */
        void
        computeIsochoricSpecificHeatCapacity(
            double* const c_v,
            const std::vector<const double*> Y,
            const hier::IntVector& num_ghosts_isochoric_specific_heat_capacity,
            const hier::IntVector& num_ghosts_mass_fractions,
            const hier::IntVector& ghostcell_dims_isochoric_specific_heat_capacity,
            const hier::IntVector& ghostcell_dims_mass_fractions,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims) const;
        
        /*
         * Compute the isochoric specific heat capacity of mixture with isothermal and isobaric
         * equilibrium assumptions.
         */
        void
        computeIsochoricSpecificHeatCapacity(
            double* const c_v,
            double* const Y_last,
            const std::vector<const double*> Y,
            const hier::IntVector& num_ghosts_isochoric_specific_heat_capacity,
            const hier::IntVector& num_ghosts_mass_fractions,
            const hier::IntVector& ghostcell_dims_isochoric_specific_heat_capacity,
            const hier::IntVector& ghostcell_dims_mass_fractions,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims) const;
        
        /*
         * Compute the isobaric specific heat capacity of mixture with isothermal and isobaric equilibrium
         * assumptions.
         */
        void
        computeIsobaricSpecificHeatCapacity(
            double* const c_p,
            const std::vector<const double*> Y,
            const hier::IntVector& num_ghosts_isobaric_specific_heat_capacity,
            const hier::IntVector& num_ghosts_mass_fractions,
            const hier::IntVector& ghostcell_dims_isobaric_specific_heat_capacity,
            const hier::IntVector& ghostcell_dims_mass_fractions,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims) const;
        
        /*
         * Compute the isobaric specific heat capacity of mixture with isothermal and isobaric equilibrium
         * assumptions.
         */
        void
        computeIsobaricSpecificHeatCapacity(
            double* const c_p,
            double* const Y_last,
            const std::vector<const double*> Y,
            const hier::IntVector& num_ghosts_isobaric_specific_heat_capacity,
            const hier::IntVector& num_ghosts_mass_fractions,
            const hier::IntVector& ghostcell_dims_isobaric_specific_heat_capacity,
            const hier::IntVector& ghostcell_dims_mass_fractions,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims) const;
        
        /*
         * Compute the mixture partial derivative of pressure w.r.t. partial densities under constant specific
         * internal energy with isothermal and isobaric equilibrium assumptions.
         */
        void
        computePressureDerivativeWithPartialDensities(
            std::vector<double*> Psi,
            const double* const epsilon,
            const double* const gamma,
            const double* const c_v,
            const hier::IntVector& num_ghosts_partial_pressure_partial_partial_densities,
            const hier::IntVector& num_ghosts_internal_energy,
            const hier::IntVector& num_ghosts_mixture_thermo_properties,
            const hier::IntVector& ghostcell_dims_partial_pressure_partial_partial_densities,
            const hier::IntVector& ghostcell_dims_internal_energy,
            const hier::IntVector& ghostcell_dims_mixture_thermo_properties,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims) const;
        
        /*
         * Compute the mixture partial derivative of pressure w.r.t. partial densities under constant specific
         * internal energy and volume fractions with isobaric equilibrium assumption.
         */
        void
        computePressureDerivativeWithPartialDensities(
            std::vector<double*> Psi,
            const double* const rho,
            const double* const p,
            const hier::IntVector& num_ghosts_partial_pressure_partial_partial_densities,
            const hier::IntVector& num_ghosts_density,
            const hier::IntVector& num_ghosts_pressure,
            const hier::IntVector& ghostcell_dims_partial_pressure_partial_partial_densities,
            const hier::IntVector& ghostcell_dims_density,
            const hier::IntVector& ghostcell_dims_pressure,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims) const;
        
        /*
         * Compute the mixture partial derivative of pressure w.r.t. volume fractions under constant specific
         * internal energy and partial densities with isobaric equilibrium assumption.
         */
        void
        computePressureDerivativeWithVolumeFractions(
            std::vector<double*> M,
            const double* const p,
            const double* const gamma,
            const hier::IntVector& num_ghosts_partial_pressure_partial_volume_fractions,
            const hier::IntVector& num_ghosts_pressure,
            const hier::IntVector& num_ghosts_mixture_thermo_properties,
            const hier::IntVector& ghostcell_dims_partial_pressure_partial_volume_fractions,
            const hier::IntVector& ghostcell_dims_pressure,
            const hier::IntVector& ghostcell_dims_mixture_thermo_properties,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims) const;
        
        /*
         * Compute the thermodynamic properties of the mixture with mass fractions.
         */
        void
        computeMixtureThermodynamicPropertiesWithMassFraction(
            double* const gamma,
            double* const R,
            double* const c_p,
            double* const c_v,
            const std::vector<const double*> Y,
            const hier::IntVector& num_ghosts_mixture_thermo_properties,
            const hier::IntVector& num_ghosts_mass_fractions,
            const hier::IntVector& ghostcell_dims_mixture_thermo_properties,
            const hier::IntVector& ghostcell_dims_mass_fractions,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims) const;
        
        /*
         * Compute the thermodynamic properties of the mixture with mass fractions.
         */
        void
        computeMixtureThermodynamicPropertiesWithMassFraction(
            double* const gamma,
            double* const R,
            double* const c_p,
            double* const c_v,
            double* const Y_last,
            const std::vector<const double*> Y,
            const hier::IntVector& num_ghosts_mixture_thermo_properties,
            const hier::IntVector& num_ghosts_mass_fractions,
            const hier::IntVector& ghostcell_dims_mixture_thermo_properties,
            const hier::IntVector& ghostcell_dims_mass_fractions,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims) const;
        
        /*
         * Compute the thermodynamic properties of the mixture with volume fractions.
         */
        void
        getMixtureThermodynamicPropertiesWithVolumeFraction(
            double* const gamma,
            const std::vector<const double*> Z,
            const hier::IntVector& num_ghosts_mixture_thermo_properties,
            const hier::IntVector& num_ghosts_volume_fractions,
            const hier::IntVector& ghostcell_dims_mixture_thermo_properties,
            const hier::IntVector& ghostcell_dims_volume_fractions,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims) const;
        
        /*
         * Compute the thermodynamic properties of the mixture with volume fractions.
         */
        void
        getMixtureThermodynamicPropertiesWithVolumeFraction(
            double* const gamma,
            double* const Z_last,
            const std::vector<const double*> Z,
            const hier::IntVector& num_ghosts_mixture_thermo_properties,
            const hier::IntVector& num_ghosts_volume_fractions,
            const hier::IntVector& ghostcell_dims_mixture_thermo_properties,
            const hier::IntVector& ghostcell_dims_volume_fractions,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims) const;
        
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
