#ifndef EQUATION_OF_STATE_IDEAL_GAS_HPP
#define EQUATION_OF_STATE_IDEAL_GAS_HPP

#include "util/mixing_rules/equations_of_state/EquationOfState.hpp"

#include <cmath>

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
        
        ~EquationOfStateIdealGas() {}
        
        /*
         * Print all characteristics of the equation of state class.
         */
        void
        printClassData(std::ostream& os) const;
        
        /*
         * Compute the pressure.
         */
        Real
        getPressure(
            const Real* const density,
            const Real* const internal_energy,
            const std::vector<const Real*>& thermo_properties) const;
        
        /*
         * Compute the pressure.
         */
        void
        computePressure(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_internal_energy,
            const std::vector<const Real*>& thermo_properties,
            const hier::Box& domain) const;
        
        /*
         * Compute the pressure.
         */
        void
        computePressure(
            HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_internal_energy,
            const std::vector<const Real*>& thermo_properties,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the pressure.
         */
        void
        computePressure(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_internal_energy,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_thermo_properties,
            const hier::Box& domain) const;
        
        /*
         * Compute the pressure.
         */
        void
        computePressure(
            HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_internal_energy,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_thermo_properties,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the sound speed.
         */
        Real
        getSoundSpeed(
            const Real* const density,
            const Real* const pressure,
            const std::vector<const Real*>& thermo_properties) const;
        
        /*
         * Compute the sound speed.
         */
        void
        computeSoundSpeed(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_sound_speed,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
            const std::vector<const Real*>& thermo_properties,
            const hier::Box& domain) const;
        
        /*
         * Compute the sound speed.
         */
        void
        computeSoundSpeed(
            HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_sound_speed,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
            const std::vector<const Real*>& thermo_properties,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the sound speed.
         */
        void
        computeSoundSpeed(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_sound_speed,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_thermo_properties,
            const hier::Box& domain) const;
        
        /*
         * Compute the sound speed.
         */
        void
        computeSoundSpeed(
            HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_sound_speed,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_thermo_properties,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the specific internal energy.
         */
        Real
        getInternalEnergy(
            const Real* const density,
            const Real* const pressure,
            const std::vector<const Real*>& thermo_properties) const;
        
        /*
         * Compute the specific internal energy.
         */
        void
        computeInternalEnergy(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_internal_energy,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
            const std::vector<const Real*>& thermo_properties,
            const hier::Box& domain) const;
        
        /*
         * Compute the specific internal energy.
         */
        void
        computeInternalEnergy(
            HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_internal_energy,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
            const std::vector<const Real*>& thermo_properties,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the specific internal energy.
         */
        void
        computeInternalEnergy(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_internal_energy,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_thermo_properties,
            const hier::Box& domain) const;
        
        /*
         * Compute the specific internal energy.
         */
        void
        computeInternalEnergy(
            HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_internal_energy,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_thermo_properties,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the specific enthalpy.
         */
        Real
        getEnthalpy(
            const Real* const density,
            const Real* const pressure,
            const std::vector<const Real*>& thermo_properties) const;
        
        /*
         * Compute the specific enthalpy.
         */
        void
        computeEnthalpy(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_enthalpy,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
            const std::vector<const Real*>& thermo_properties,
            const hier::Box& domain) const;
        
        /*
         * Compute the specific enthalpy.
         */
        void
        computeEnthalpy(
            HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_enthalpy,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
            const std::vector<const Real*>& thermo_properties,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the specific enthalpy.
         */
        void
        computeEnthalpy(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_enthalpy,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_thermo_properties,
            const hier::Box& domain) const;
        
        /*
         * Compute the specific enthalpy.
         */
        void
        computeEnthalpy(
            HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_enthalpy,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_thermo_properties,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the temperature.
         */
        Real
        getTemperature(
            const Real* const density,
            const Real* const pressure,
            const std::vector<const Real*>& thermo_properties) const;
        
        /*
         * Compute the temperature.
         */
        void
        computeTemperature(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_temperature,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
            const std::vector<const Real*>& thermo_properties,
            const hier::Box& domain) const;
        
        /*
         * Compute the temperature.
         */
        void
        computeTemperature(
            HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_temperature,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
            const std::vector<const Real*>& thermo_properties,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the temperature.
         */
        void
        computeTemperature(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_temperature,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_thermo_properties,
            const hier::Box& domain) const;
        
        /*
         * Compute the temperature.
         */
        void
        computeTemperature(
            HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_temperature,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_thermo_properties,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the specific internal energy from temperature.
         */
        Real
        getInternalEnergyFromTemperature(
            const Real* const density,
            const Real* const temperature,
            const std::vector<const Real*>& thermo_properties) const;
        
        /*
         * Compute the specific internal energy from temperature.
         */
        void
        computeInternalEnergyFromTemperature(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_internal_energy,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_temperature,
            const std::vector<const Real*>& thermo_properties,
            const hier::Box& domain) const;
        
        /*
         * Compute the specific internal energy from temperature.
         */
        void
        computeInternalEnergyFromTemperature(
            HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_internal_energy,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_temperature,
            const std::vector<const Real*>& thermo_properties,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the specific internal energy from temperature.
         */
        void
        computeInternalEnergyFromTemperature(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_internal_energy,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_temperature,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_thermo_properties,
            const hier::Box& domain) const;
        
        /*
         * Compute the specific internal energy from temperature.
         */
        void
        computeInternalEnergyFromTemperature(
            HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_internal_energy,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_temperature,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_thermo_properties,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the isochoric specific heat capacity.
         */
        Real
        getIsochoricSpecificHeatCapacity(
            const Real* const density,
            const Real* const pressure,
            const std::vector<const Real*>& thermo_properties) const;
        
        /*
         * Compute the isochoric specific heat capacity.
         */
        void
        computeIsochoricSpecificHeatCapacity(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_isochoric_specific_heat_capacity,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
            const std::vector<const Real*>& thermo_properties,
            const hier::Box& domain) const;
        
        /*
         * Compute the isochoric specific heat capacity.
         */
        void
        computeIsochoricSpecificHeatCapacity(
            HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_isochoric_specific_heat_capacity,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
            const std::vector<const Real*>& thermo_properties,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the isochoric specific heat capacity.
         */
        void
        computeIsochoricSpecificHeatCapacity(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_isochoric_specific_heat_capacity,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_thermo_properties,
            const hier::Box& domain) const;
        
        /*
         * Compute the isochoric specific heat capacity.
         */
        void
        computeIsochoricSpecificHeatCapacity(
            HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_isochoric_specific_heat_capacity,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_thermo_properties,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the isobaric specific heat capacity.
         */
        Real
        getIsobaricSpecificHeatCapacity(
            const Real* const density,
            const Real* const pressure,
            const std::vector<const Real*>& thermo_properties) const;
        
        /*
         * Compute the isobaric specific heat capacity.
         */
        void
        computeIsobaricSpecificHeatCapacity(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_isobaric_specific_heat_capacity,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
            const std::vector<const Real*>& thermo_properties,
            const hier::Box& domain) const;
        
        /*
         * Compute the isobaric specific heat capacity.
         */
        void
        computeIsobaricSpecificHeatCapacity(
            HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_isobaric_specific_heat_capacity,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
            const std::vector<const Real*>& thermo_properties,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the isobaric specific heat capacity.
         */
        void
        computeIsobaricSpecificHeatCapacity(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_isobaric_specific_heat_capacity,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_thermo_properties,
            const hier::Box& domain) const;
        
        /*
         * Compute the isobaric specific heat capacity.
         */
        void
        computeIsobaricSpecificHeatCapacity(
            HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_isobaric_specific_heat_capacity,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_thermo_properties,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the Gruneisen parameter (partial derivative of pressure w.r.t. specific internal energy under
         * constant density divided by density).
         */
        Real
        getGruneisenParameter(
            const Real* const density,
            const Real* const pressure,
            const std::vector<const Real*>& thermo_properties) const;
        
        /*
         * Compute the Gruneisen parameter (partial derivative of pressure w.r.t. specific internal energy under
         * constant density divided by density).
         */
        void
        computeGruneisenParameter(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_gruneisen_parameter,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
            const std::vector<const Real*>& thermo_properties,
            const hier::Box& domain) const;
        
        /*
         * Compute the Gruneisen parameter (partial derivative of pressure w.r.t. specific internal energy under
         * constant density divided by density).
         */
        void
        computeGruneisenParameter(
            HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_gruneisen_parameter,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
            const std::vector<const Real*>& thermo_properties,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the Gruneisen parameter (partial derivative of pressure w.r.t. specific internal energy under
         * constant density divided by density).
         */
        void
        computeGruneisenParameter(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_gruneisen_parameter,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_thermo_properties,
            const hier::Box& domain) const;
        
        /*
         * Compute the Gruneisen parameter (partial derivative of pressure w.r.t. specific internal energy under
         * constant density divided by density).
         */
        void
        computeGruneisenParameter(
            HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_gruneisen_parameter,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_thermo_properties,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the partial derivative of pressure w.r.t. density under constant specific internal energy.
         */
        Real
        getPressureDerivativeWithDensity(
            const Real* const density,
            const Real* const pressure,
            const std::vector<const Real*>& thermo_properties) const;
        
        /*
         * Compute the partial derivative of pressure w.r.t. density under constant specific internal energy.
         */
        void
        computePressureDerivativeWithDensity(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_partial_pressure_partial_density,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
            const std::vector<const Real*>& thermo_properties,
            const hier::Box& domain) const;
        
        /*
         * Compute the partial derivative of pressure w.r.t. density under constant specific internal energy.
         */
        void
        computePressureDerivativeWithDensity(
            HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_partial_pressure_partial_density,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
            const std::vector<const Real*>& thermo_properties,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the partial derivative of pressure w.r.t. density under constant specific internal energy.
         */
        void
        computePressureDerivativeWithDensity(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_partial_pressure_partial_density,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_thermo_properties,
            const hier::Box& domain) const;
        
        /*
         * Compute the partial derivative of pressure w.r.t. density under constant specific internal energy.
         */
        void
        computePressureDerivativeWithDensity(
            HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_partial_pressure_partial_density,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_thermo_properties,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the density.
         */
        Real
        getDensity(
            const Real* const pressure,
            const Real* const temperature,
            const std::vector<const Real*>& thermo_properties) const;
        
        /*
         * Compute the density.
         */
        void
        computeDensity(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_temperature,
            const std::vector<const Real*>& thermo_properties,
            const hier::Box& domain) const;
        
        /*
         * Compute the density.
         */
        void
        computeDensity(
            HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_temperature,
            const std::vector<const Real*>& thermo_properties,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the density.
         */
        void
        computeDensity(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_temperature,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_thermo_properties,
            const hier::Box& domain) const;
        
        /*
         * Compute the density.
         */
        void
        computeDensity(
            HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_temperature,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_thermo_properties,
            int side_normal,
            const hier::Box& domain) const;
        
    private:
        /*
         * Compute the pressure.
         */
        void
        computePressure(
            Real* const p,
            const Real* const rho,
            const Real* const epsilon,
            const Real& gamma,
            const hier::IntVector& offset_pressure,
            const hier::IntVector& offset_density,
            const hier::IntVector& offset_internal_energy,
            const hier::IntVector& ghostcell_dims_pressure,
            const hier::IntVector& ghostcell_dims_density,
            const hier::IntVector& ghostcell_dims_internal_energy,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims) const;
        
        /*
         * Compute the pressure.
         */
        void
        computePressure(
            Real* const p,
            const Real* const rho,
            const Real* const epsilon,
            const Real* const gamma,
            const hier::IntVector& offset_pressure,
            const hier::IntVector& offset_density,
            const hier::IntVector& offset_internal_energy,
            const hier::IntVector& offset_thermo_properties,
            const hier::IntVector& ghostcell_dims_pressure,
            const hier::IntVector& ghostcell_dims_density,
            const hier::IntVector& ghostcell_dims_internal_energy,
            const hier::IntVector& ghostcell_dims_thermo_properties,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims) const;
        
        /*
         * Compute the sound speed.
         */
        void
        computeSoundSpeed(
            Real* const c,
            const Real* const rho,
            const Real* const p,
            const Real& gamma,
            const hier::IntVector& offset_sound_speed,
            const hier::IntVector& offset_density,
            const hier::IntVector& offset_pressure,
            const hier::IntVector& ghostcell_dims_sound_speed,
            const hier::IntVector& ghostcell_dims_density,
            const hier::IntVector& ghostcell_dims_pressure,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims) const;
        
        /*
         * Compute the sound speed.
         */
        void
        computeSoundSpeed(
            Real* const c,
            const Real* const rho,
            const Real* const p,
            const Real* const gamma,
            const hier::IntVector& offset_sound_speed,
            const hier::IntVector& offset_density,
            const hier::IntVector& offset_pressure,
            const hier::IntVector& offset_thermo_properties,
            const hier::IntVector& ghostcell_dims_sound_speed,
            const hier::IntVector& ghostcell_dims_density,
            const hier::IntVector& ghostcell_dims_pressure,
            const hier::IntVector& ghostcell_dims_thermo_properties,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims) const;
        
        /*
         * Compute the specific internal energy.
         */
        void
        computeInternalEnergy(
            Real* const epsilon,
            const Real* const rho,
            const Real* const p,
            const Real& gamma,
            const hier::IntVector& offset_internal_energy,
            const hier::IntVector& offset_density,
            const hier::IntVector& offset_pressure,
            const hier::IntVector& ghostcell_dims_internal_energy,
            const hier::IntVector& ghostcell_dims_density,
            const hier::IntVector& ghostcell_dims_pressure,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims) const;
        
        /*
         * Compute the specific internal energy.
         */
        void
        computeInternalEnergy(
            Real* const epsilon,
            const Real* const rho,
            const Real* const p,
            const Real* const gamma,
            const hier::IntVector& offset_internal_energy,
            const hier::IntVector& offset_density,
            const hier::IntVector& offset_pressure,
            const hier::IntVector& offset_thermo_properties,
            const hier::IntVector& ghostcell_dims_internal_energy,
            const hier::IntVector& ghostcell_dims_density,
            const hier::IntVector& ghostcell_dims_pressure,
            const hier::IntVector& ghostcell_dims_thermo_properties,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims) const;
        
        /*
         * Compute the specific enthalpy.
         */
        void
        computeEnthalpy(
            Real* const h,
            const Real* const rho,
            const Real* const p,
            const Real& gamma,
            const hier::IntVector& offset_enthalpy,
            const hier::IntVector& offset_density,
            const hier::IntVector& offset_pressure,
            const hier::IntVector& ghostcell_dims_enthalpy,
            const hier::IntVector& ghostcell_dims_density,
            const hier::IntVector& ghostcell_dims_pressure,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims) const;
        
        /*
         * Compute the specific enthalpy.
         */
        void
        computeEnthalpy(
            Real* const h,
            const Real* const rho,
            const Real* const p,
            const Real* const gamma,
            const hier::IntVector& offset_enthalpy,
            const hier::IntVector& offset_density,
            const hier::IntVector& offset_pressure,
            const hier::IntVector& offset_thermo_properties,
            const hier::IntVector& ghostcell_dims_enthalpy,
            const hier::IntVector& ghostcell_dims_density,
            const hier::IntVector& ghostcell_dims_pressure,
            const hier::IntVector& ghostcell_dims_thermo_properties,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims) const;
        
        /*
         * Compute the temperature.
         */
        void
        computeTemperature(
            Real* const T,
            const Real* const rho,
            const Real* const p,
            const Real& gamma,
            const Real& c_v,
            const hier::IntVector& offset_temperature,
            const hier::IntVector& offset_density,
            const hier::IntVector& offset_pressure,
            const hier::IntVector& ghostcell_dims_temperature,
            const hier::IntVector& ghostcell_dims_density,
            const hier::IntVector& ghostcell_dims_pressure,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims) const;
        
        /*
         * Compute the temperature.
         */
        void
        computeTemperature(
            Real* const T,
            const Real* const rho,
            const Real* const p,
            const Real* const gamma,
            const Real* const c_v,
            const hier::IntVector& offset_temperature,
            const hier::IntVector& offset_density,
            const hier::IntVector& offset_pressure,
            const hier::IntVector& offset_thermo_properties,
            const hier::IntVector& ghostcell_dims_temperature,
            const hier::IntVector& ghostcell_dims_density,
            const hier::IntVector& ghostcell_dims_pressure,
            const hier::IntVector& ghostcell_dims_thermo_properties,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims) const;
        
        /*
         * Compute the specific internal energy from temperature.
         */
        void
        computeInternalEnergyFromTemperature(
            Real* const epsilon,
            const Real* const T,
            const Real& c_v,
            const hier::IntVector& offset_internal_energy,
            const hier::IntVector& offset_temperature,
            const hier::IntVector& ghostcell_dims_internal_energy,
            const hier::IntVector& ghostcell_dims_temperature,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims) const;
        
        /*
         * Compute the specific internal energy from temperature.
         */
        void
        computeInternalEnergyFromTemperature(
            Real* const epsilon,
            const Real* const T,
            const Real* const c_v,
            const hier::IntVector& offset_internal_energy,
            const hier::IntVector& offset_temperature,
            const hier::IntVector& offset_thermo_properties,
            const hier::IntVector& ghostcell_dims_internal_energy,
            const hier::IntVector& ghostcell_dims_temperature,
            const hier::IntVector& ghostcell_dims_thermo_properties,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims) const;
        
        /*
         * Compute the isochoric specific heat capacity.
         */
        void
        computeIsochoricSpecificHeatCapacity(
            Real* const c_v,
            const Real& c_v_src,
            const hier::IntVector& offset_isochoric_specific_heat_capacity,
            const hier::IntVector& ghostcell_dims_isochoric_specific_heat_capacity,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims) const;
        
        /*
         * Compute the isochoric specific heat capacity.
         */
        void
        computeIsochoricSpecificHeatCapacity(
            Real* const c_v,
            const Real* const c_v_src,
            const hier::IntVector& offset_isochoric_specific_heat_capacity,
            const hier::IntVector& offset_thermo_properties,
            const hier::IntVector& ghostcell_dims_isochoric_specific_heat_capacity,
            const hier::IntVector& ghostcell_dims_thermo_properties,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims) const;
        
        /*
         * Compute the isobaric specific heat capacity.
         */
        void
        computeIsobaricSpecificHeatCapacity(
            Real* const c_p,
            const Real& c_p_src,
            const hier::IntVector& offset_isobaric_specific_heat_capacity,
            const hier::IntVector& ghostcell_dims_isobaric_specific_heat_capacity,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims) const;
        
        /*
         * Compute the isobaric specific heat capacity.
         */
        void
        computeIsobaricSpecificHeatCapacity(
            Real* const c_p,
            const Real* const c_p_src,
            const hier::IntVector& offset_isobaric_specific_heat_capacity,
            const hier::IntVector& offset_thermo_properties,
            const hier::IntVector& ghostcell_dims_isobaric_specific_heat_capacity,
            const hier::IntVector& ghostcell_dims_thermo_properties,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims) const;
        
        /*
         * Compute the Gruneisen parameter (partial derivative of pressure w.r.t. specific internal energy under
         * constant density divided by density).
         */
        void
        computeGruneisenParameter(
            Real* const Gamma,
            const Real& gamma,
            const hier::IntVector& offset_gruneisen_parameter,
            const hier::IntVector& ghostcell_dims_gruneisen_parameter,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims) const;
        
        /*
         * Compute the Gruneisen parameter (partial derivative of pressure w.r.t. specific internal energy under
         * constant density divided by density).
         */
        void
        computeGruneisenParameter(
            Real* const Gamma,
            const Real* const gamma,
            const hier::IntVector& offset_gruneisen_parameter,
            const hier::IntVector& offset_thermo_properties,
            const hier::IntVector& ghostcell_dims_gruneisen_parameter,
            const hier::IntVector& ghostcell_dims_thermo_properties,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims) const;
        
        /*
         * Compute the partial derivative of pressure w.r.t. density under constant specific internal energy.
         */
        void
        computePressureDerivativeWithDensity(
            Real* const Psi,
            const Real* const rho,
            const Real* const p,
            const hier::IntVector& offset_partial_pressure_partial_density,
            const hier::IntVector& offset_density,
            const hier::IntVector& offset_pressure,
            const hier::IntVector& ghostcell_dims_partial_pressure_partial_density,
            const hier::IntVector& ghostcell_dims_density,
            const hier::IntVector& ghostcell_dims_pressure,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims) const;
        
        /*
         * Compute the density.
         */
        void
        computeDensity(
            Real* const rho,
            const Real* const p,
            const Real* const T,
            const Real& R,
            const hier::IntVector& offset_density,
            const hier::IntVector& offset_pressure,
            const hier::IntVector& offset_temperature,
            const hier::IntVector& ghostcell_dims_density,
            const hier::IntVector& ghostcell_dims_pressure,
            const hier::IntVector& ghostcell_dims_temperature,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims) const;
        
        /*
         * Compute the density.
         */
        void
        computeDensity(
            Real* const rho,
            const Real* const p,
            const Real* const T,
            const Real* const R,
            const hier::IntVector& offset_density,
            const hier::IntVector& offset_pressure,
            const hier::IntVector& offset_temperature,
            const hier::IntVector& offset_thermo_properties,
            const hier::IntVector& ghostcell_dims_density,
            const hier::IntVector& ghostcell_dims_pressure,
            const hier::IntVector& ghostcell_dims_temperature,
            const hier::IntVector& ghostcell_dims_thermo_properties,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims) const;
        
};

#endif /* EQUATION_OF_STATE_IDEAL_GAS_HPP */
