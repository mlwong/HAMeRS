#ifndef EQUATION_OF_STATE_MIXING_RULES_STIFFENED_GAS_HPP
#define EQUATION_OF_STATE_MIXING_RULES_STIFFENED_GAS_HPP

#include "util/mixing_rules/equations_of_state/EquationOfStateMixingRules.hpp"

#include "util/mixing_rules/equations_of_state/stiffened_gas/EquationOfStateStiffenedGas.hpp"

class EquationOfStateMixingRulesStiffenedGas: public EquationOfStateMixingRules
{
    public:        
        EquationOfStateMixingRulesStiffenedGas(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const int& num_species,
            const MIXING_CLOSURE_MODEL::TYPE& mixing_closure_model,
            const HAMERS_SHARED_PTR<tbox::Database>& equation_of_state_mixing_rules_db);
        
        ~EquationOfStateMixingRulesStiffenedGas() {}
        
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
            const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const;
        
        /*
         * Return the HAMERS_SHARED_PTR to the equation of state.
         */
        const HAMERS_SHARED_PTR<EquationOfState>&
        getEquationOfState(const int species_index = 0) const
        {
            NULL_USE(species_index);
            return d_equation_of_state;
        }
        
        /*
         * Compute the pressure of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        Real
        getPressure(
            const Real* const density,
            const Real* const internal_energy,
            const std::vector<const Real*>& mass_fractions) const;
        
        /*
         * Compute the pressure of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        void
        computePressure(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_internal_energy,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_mass_fractions,
            const hier::Box& domain) const;
        
        /*
         * Compute the pressure of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        void
        computePressure(
            HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_internal_energy,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_mass_fractions,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the pressure of the mixture with isobaric equilibrium assumption.
         */
        Real
        getPressure(
            const Real* const density,
            const Real* const internal_energy,
            const std::vector<const Real*>& mass_fractions,
            const std::vector<const Real*>& volume_fractions) const;
        
        /*
         * Compute the pressure of the mixture with isobaric equilibrium assumption.
         */
        void
        computePressure(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_internal_energy,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_mass_fractions,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_volume_fractions,
            const hier::Box& domain) const;
        
        /*
         * Compute the pressure of the mixture with isobaric equilibrium assumption.
         */
        void
        computePressure(
            HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_internal_energy,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_mass_fractions,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_volume_fractions,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the specific internal energy of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        Real
        getInternalEnergy(
            const Real* const density,
            const Real* const pressure,
            const std::vector<const Real*>& mass_fractions) const;
        
        /*
         * Compute the specific internal energy of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        void
        computeInternalEnergy(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_internal_energy,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_mass_fractions,
            const hier::Box& domain) const;
        
        /*
         * Compute the specific internal energy of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        void
        computeInternalEnergy(
            HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_internal_energy,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_mass_fractions,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the specific internal energy of the mixture with isobaric equilibrium assumption.
         */
        Real
        getInternalEnergy(
            const Real* const density,
            const Real* const pressure,
            const std::vector<const Real*>& mass_fractions,
            const std::vector<const Real*>& volume_fractions) const;
        
        /*
         * Compute the specific internal energy of the mixture with isobaric equilibrium assumption.
         */
        void
        computeInternalEnergy(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_internal_energy,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_mass_fractions,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_volume_fractions,
            const hier::Box& domain) const;
        
        /*
         * Compute the specific internal energy of the mixture with isobaric equilibrium assumption.
         */
        void
        computeInternalEnergy(
            HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_internal_energy,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_mass_fractions,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_volume_fractions,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the temperature of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        Real
        getTemperature(
            const Real* const density,
            const Real* const pressure,
            const std::vector<const Real*>& mass_fractions) const;
        
        /*
         * Compute the temperature of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        void
        computeTemperature(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_temperature,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_mass_fractions,
            const hier::Box& domain) const;
        
        /*
         * Compute the temperature of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        void
        computeTemperature(
            HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_temperature,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_mass_fractions,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the specific internal energy of the mixture from temperature with isothermal
         * and isobaric equilibrium assumptions.
         */
        Real
        getInternalEnergyFromTemperature(
            const Real* const density,
            const Real* const temperature,
            const std::vector<const Real*>& mass_fractions) const;
        
        /*
         * Compute the specific internal energy of the mixture from temperature with isothermal
         * and isobaric equilibrium assumptions.
         */
        void
        computeInternalEnergyFromTemperature(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_internal_energy,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_temperature,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_mass_fractions,
            const hier::Box& domain) const;
        
        /*
         * Compute the specific internal energy of the mixture from temperature with isothermal
         * and isobaric equilibrium assumptions.
         */
        void
        computeInternalEnergyFromTemperature(
            HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_internal_energy,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_temperature,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_mass_fractions,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the isochoric specific heat capacity of mixture with isothermal and isobaric
         * equilibrium assumptions.
         */
        Real
        getIsochoricSpecificHeatCapacity(
            const Real* const density,
            const Real* const pressure,
            const std::vector<const Real*>& mass_fractions) const;
        
        /*
         * Compute the isochoric specific heat capacity of mixture with isothermal and isobaric
         * equilibrium assumptions.
         */
        void
        computeIsochoricSpecificHeatCapacity(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_isochoric_specific_heat_capacity,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_mass_fractions,
            const hier::Box& domain) const;
        
        /*
         * Compute the isochoric specific heat capacity of mixture with isothermal and isobaric
         * equilibrium assumptions.
         */
        void
        computeIsochoricSpecificHeatCapacity(
            HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_isochoric_specific_heat_capacity,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_mass_fractions,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the isobaric specific heat capacity of mixture with isothermal and isobaric
         * equilibrium assumptions.
         */
        Real
        getIsobaricSpecificHeatCapacity(
            const Real* const density,
            const Real* const pressure,
            const std::vector<const Real*>& mass_fractions) const;
        
        /*
         * Compute the isobaric specific heat capacity of mixture with isothermal and isobaric
         * equilibrium assumptions.
         */
        void
        computeIsobaricSpecificHeatCapacity(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_isobaric_specific_heat_capacity,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_mass_fractions,
            const hier::Box& domain) const;
        
        /*
         * Compute the isobaric specific heat capacity of mixture with isothermal and isobaric
         * equilibrium assumptions.
         */
        void
        computeIsobaricSpecificHeatCapacity(
            HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_isobaric_specific_heat_capacity,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_mass_fractions,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the Gruneisen parameter of the mixture with isothermal and isobaric equilibrium assumptions
         * (partial derivative of pressure w.r.t. specific internal energy under constant partial densities
         * divided by mixture density).
         */
        Real
        getGruneisenParameter(
            const Real* const density,
            const Real* const pressure,
            const std::vector<const Real*>& mass_fractions) const;
        
        /*
         * Compute the Gruneisen parameter of the mixture with isothermal and isobaric equilibrium assumptions
         * (partial derivative of pressure w.r.t. specific internal energy under constant partial densities
         * divided by mixture density).
         */
        void
        computeGruneisenParameter(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_gruneisen_parameter,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_mass_fractions,
            const hier::Box& domain) const;
        
        /*
         * Compute the Gruneisen parameter of the mixture with isothermal and isobaric equilibrium assumptions
         * (partial derivative of pressure w.r.t. specific internal energy under constant partial densities
         * divided by mixture density).
         */
        void
        computeGruneisenParameter(
            HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_gruneisen_parameter,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_mass_fractions,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the Gruneisen parameter of the mixture with isobaric equilibrium assumption
         * (partial derivative of pressure w.r.t. specific internal energy under constant partial densities
         * and volume fractions divided by mixture density).
         */
        Real
        getGruneisenParameter(
            const Real* const density,
            const Real* const pressure,
            const std::vector<const Real*>& mass_fractions,
            const std::vector<const Real*>& volume_fractions) const;
        
        /*
         * Compute the Gruneisen parameter of the mixture with isobaric equilibrium assumption
         * (partial derivative of pressure w.r.t. specific internal energy under constant partial densities
         * and volume fractions divided by mixture density).
         */
        void
        computeGruneisenParameter(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_gruneisen_parameter,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_mass_fractions,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_volume_fractions,
            const hier::Box& domain) const;
        
        /*
         * Compute the Gruneisen parameter of the mixture with isobaric equilibrium assumption
         * (partial derivative of pressure w.r.t. specific internal energy under constant partial densities
         * and volume fractions divided by mixture density).
         */
        void
        computeGruneisenParameter(
            HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_gruneisen_parameter,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_mass_fractions,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_volume_fractions,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the mixture partial derivative of pressure w.r.t. partial densities under constant specific
         * internal energy with isothermal and isobaric equilibrium assumptions.
         */
        std::vector<Real>
        getPressureDerivativeWithPartialDensities(
            const Real* const density,
            const Real* const pressure,
            const std::vector<const Real*>& mass_fractions) const;
        
        /*
         * Compute the mixture partial derivative of pressure w.r.t. partial densities under constant specific
         * internal energy with isothermal and isobaric equilibrium assumptions.
         */
        void
        computePressureDerivativeWithPartialDensities(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_partial_pressure_partial_partial_densities,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_mass_fractions,
            const hier::Box& domain) const;
        
        /*
         * Compute the mixture partial derivative of pressure w.r.t. partial densities under constant specific
         * internal energy with isothermal and isobaric equilibrium assumptions.
         */
        void
        computePressureDerivativeWithPartialDensities(
            HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_partial_pressure_partial_partial_densities,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_mass_fractions,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the mixture partial derivative of pressure w.r.t. partial densities under constant specific
         * internal energy and volume fractions with isobaric equilibrium assumption.
         */
        std::vector<Real>
        getPressureDerivativeWithPartialDensities(
            const Real* const density,
            const Real* const pressure,
            const std::vector<const Real*>& mass_fractions,
            const std::vector<const Real*>& volume_fractions) const;
        
        /*
         * Compute the mixture partial derivative of pressure w.r.t. partial densities under constant specific
         * internal energy and volume fractions with isobaric equilibrium assumption.
         */
        void
        computePressureDerivativeWithPartialDensities(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_partial_pressure_partial_partial_densities,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_mass_fractions,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_volume_fractions,
            const hier::Box& domain) const;
        
        /*
         * Compute the mixture partial derivative of pressure w.r.t. partial densities under constant specific
         * internal energy and volume fractions with isobaric equilibrium assumption.
         */
        void
        computePressureDerivativeWithPartialDensities(
            HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_partial_pressure_partial_partial_densities,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_mass_fractions,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_volume_fractions,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the mixture partial derivative of pressure w.r.t. volume fractions under constant specific
         * internal energy and partial densities with isobaric equilibrium assumption.
         */
        std::vector<Real>
        getPressureDerivativeWithVolumeFractions(
            const Real* const density,
            const Real* const pressure,
            const std::vector<const Real*>& mass_fractions,
            const std::vector<const Real*>& volume_fractions) const;
        
        /*
         * Compute the mixture partial derivative of pressure w.r.t. volume fractions under constant specific
         * internal energy and partial densities with isobaric equilibrium assumption.
         */
        void
        computePressureDerivativeWithVolumeFractions(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_partial_pressure_partial_volume_fractions,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_mass_fractions,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_volume_fractions,
            const hier::Box& domain) const;
        
        /*
         * Compute the mixture partial derivative of pressure w.r.t. volume fractions under constant specific
         * internal energy and partial densities with isobaric equilibrium assumption.
         */
        void
        computePressureDerivativeWithVolumeFractions(
            HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_partial_pressure_partial_volume_fractions,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_mass_fractions,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_volume_fractions,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the density of mixture with isothermal and isobaric equilibrium assumptions.
         */
        Real
        getMixtureDensity(
            const Real* const pressure,
            const Real* const temperature,
            const std::vector<const Real*>& mass_fractions) const;
        
        /*
         * Compute the density of mixture with isothermal and isobaric equilibrium assumptions.
         */
        void
        computeMixtureDensity(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_mixture_density,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_temperature,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_mass_fractions,
            const hier::Box& domain) const;
        
        /*
         * Compute the density of mixture with isothermal and isobaric equilibrium assumptions.
         */
        void
        computeMixtureDensity(
            HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_mixture_density,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_temperature,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_mass_fractions,
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
            std::vector<Real*>& species_thermo_properties,
            const int species_index = 0) const;
        
    private:
        /*
         * Get the number of thermodynamic properties of the mixture.
         */
        int
        getNumberOfMixtureThermodynamicProperties() const;
        
        /*
         * Compute the thermodynamic properties of the mixture with volume fractions.
         */
        void
        getMixtureThermodynamicPropertiesWithVolumeFractions(
            std::vector<Real*>& mixture_thermo_properties,
            const std::vector<const Real*>& volume_fractions) const;
        
        /*
         * Compute the thermodynamic properties of the mixture with volume fractions.
         */
        void
        computeMixtureThermodynamicPropertiesWithVolumeFractions(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_mixture_thermo_properties,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_volume_fractions,
            const hier::Box& domain) const;
        
        /*
         * Compute the thermodynamic properties of the mixture with volume fractions.
         */
        void
        computeMixtureThermodynamicPropertiesWithVolumeFractions(
            HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_mixture_thermo_properties,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_volume_fractions,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the mixture partial derivative of pressure w.r.t. partial densities under constant specific
         * internal energy and volume fractions with isobaric equilibrium assumption.
         */
        void
        computePressureDerivativeWithPartialDensities(
            std::vector<Real*>& Psi,
            const Real* const rho,
            const Real* const p,
            const Real* const gamma,
            const Real* const p_inf,
            const hier::IntVector& offset_partial_pressure_partial_partial_densities,
            const hier::IntVector& offset_density,
            const hier::IntVector& offset_pressure,
            const hier::IntVector& offset_mixture_thermo_properties,
            const hier::IntVector& ghostcell_dims_partial_pressure_partial_partial_densities,
            const hier::IntVector& ghostcell_dims_density,
            const hier::IntVector& ghostcell_dims_pressure,
            const hier::IntVector& ghostcell_dims_mixture_thermo_properties,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims) const;
        
        /*
         * Compute the mixture partial derivative of pressure w.r.t. volume fractions under constant specific
         * internal energy and partial densities with isobaric equilibrium assumption.
         */
        void
        computePressureDerivativeWithVolumeFractions(
            std::vector<Real*>& M,
            const Real* const p,
            const Real* const gamma,
            const hier::IntVector& offset_partial_pressure_partial_volume_fractions,
            const hier::IntVector& offset_pressure,
            const hier::IntVector& offset_mixture_thermo_properties,
            const hier::IntVector& ghostcell_dims_partial_pressure_partial_volume_fractions,
            const hier::IntVector& ghostcell_dims_pressure,
            const hier::IntVector& ghostcell_dims_mixture_thermo_properties,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims) const;
        
        /*
         * Compute the thermodynamic properties of the mixture with volume fractions.
         */
        void
        computeMixtureThermodynamicPropertiesWithVolumeFractions(
            Real* const gamma,
            Real* const p_inf,
            const std::vector<const Real*>& Z,
            const hier::IntVector& offset_mixture_thermo_properties,
            const hier::IntVector& offset_volume_fractions,
            const hier::IntVector& ghostcell_dims_mixture_thermo_properties,
            const hier::IntVector& ghostcell_dims_volume_fractions,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims) const;
        
        /*
         * Compute the thermodynamic properties of the mixture with volume fractions.
         */
        void
        computeMixtureThermodynamicPropertiesWithVolumeFractions(
            Real* const gamma,
            Real* const p_inf,
            Real* const Z_last,
            const std::vector<const Real*>& Z,
            const hier::IntVector& offset_mixture_thermo_properties,
            const hier::IntVector& offset_volume_fractions_last,
            const hier::IntVector& offset_volume_fractions,
            const hier::IntVector& ghostcell_dims_mixture_thermo_properties,
            const hier::IntVector& ghostcell_dims_volume_fractions_last,
            const hier::IntVector& ghostcell_dims_volume_fractions,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims) const;
        
        /*
         * Ratio of specific heats of different species.
         */
        std::vector<Real> d_species_gamma;
        
        /*
         * Reference pressure.
         */
        std::vector<Real> d_species_p_inf;
        
        /*
         * Specific heats of different species.
         */
        std::vector<Real> d_species_c_p;
        std::vector<Real> d_species_c_v;
        
        /*
         * HAMERS_SHARED_PTR to EquationOfState.
         */
        HAMERS_SHARED_PTR<EquationOfState> d_equation_of_state;
        
};

#endif /* EQUATION_OF_STATE_MIXING_RULES_STIFFENED_GAS_HPP */
