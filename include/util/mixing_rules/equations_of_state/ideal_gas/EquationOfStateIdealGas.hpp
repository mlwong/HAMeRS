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
         * Compute the pressure.
         */
        void
        computePressure(
            boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
            const std::vector<const double*>& thermo_properties,
            const hier::Box& domain) const;
        
        /*
         * Compute the pressure.
         */
        void
        computePressure(
            boost::shared_ptr<pdat::SideData<double> >& data_pressure,
            const boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_internal_energy,
            const std::vector<const double*>& thermo_properties,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the pressure.
         */
        void
        computePressure(
            boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties,
            const hier::Box& domain) const;
        
        /*
         * Compute the pressure.
         */
        void
        computePressure(
            boost::shared_ptr<pdat::SideData<double> >& data_pressure,
            const boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::SideData<double> >& data_thermo_properties,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the sound speed.
         */
        double
        getSoundSpeed(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& thermo_properties) const;
        
        /*
         * Compute the sound speed.
         */
        void
        computeSoundSpeed(
            boost::shared_ptr<pdat::CellData<double> >& data_sound_speed,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const std::vector<const double*>& thermo_properties,
            const hier::Box& domain) const;
        
        /*
         * Compute the sound speed.
         */
        void
        computeSoundSpeed(
            boost::shared_ptr<pdat::SideData<double> >& data_sound_speed,
            const boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
            const std::vector<const double*>& thermo_properties,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the sound speed.
         */
        void
        computeSoundSpeed(
            boost::shared_ptr<pdat::CellData<double> >& data_sound_speed,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties,
            const hier::Box& domain) const;
        
        /*
         * Compute the sound speed.
         */
        void
        computeSoundSpeed(
            boost::shared_ptr<pdat::SideData<double> >& data_sound_speed,
            const boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
            const boost::shared_ptr<pdat::SideData<double> >& data_thermo_properties,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the specific internal energy.
         */
        double
        getInternalEnergy(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& thermo_properties) const;
        
        /*
         * Compute the specific internal energy.
         */
        void
        computeInternalEnergy(
            boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const std::vector<const double*>& thermo_properties,
            const hier::Box& domain) const;
        
        /*
         * Compute the specific internal energy.
         */
        void
        computeInternalEnergy(
            boost::shared_ptr<pdat::SideData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
            const std::vector<const double*>& thermo_properties,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the specific internal energy.
         */
        void
        computeInternalEnergy(
            boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties,
            const hier::Box& domain) const;
        
        /*
         * Compute the specific internal energy.
         */
        void
        computeInternalEnergy(
            boost::shared_ptr<pdat::SideData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
            const boost::shared_ptr<pdat::SideData<double> >& data_thermo_properties,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the specific enthalpy.
         */
        double
        getEnthalpy(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& thermo_properties) const;
        
        /*
         * Compute the specific enthalpy.
         */
        void
        computeEnthalpy(
            boost::shared_ptr<pdat::CellData<double> >& data_enthalpy,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const std::vector<const double*>& thermo_properties,
            const hier::Box& domain) const;
        
        /*
         * Compute the specific enthalpy.
         */
        void
        computeEnthalpy(
            boost::shared_ptr<pdat::SideData<double> >& data_enthalpy,
            const boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
            const std::vector<const double*>& thermo_properties,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the specific enthalpy.
         */
        void
        computeEnthalpy(
            boost::shared_ptr<pdat::CellData<double> >& data_enthalpy,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties,
            const hier::Box& domain) const;
        
        /*
         * Compute the specific enthalpy.
         */
        void
        computeEnthalpy(
            boost::shared_ptr<pdat::SideData<double> >& data_enthalpy,
            const boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
            const boost::shared_ptr<pdat::SideData<double> >& data_thermo_properties,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the temperature.
         */
        double
        getTemperature(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& thermo_properties) const;
        
        /*
         * Compute the temperature.
         */
        void
        computeTemperature(
            boost::shared_ptr<pdat::CellData<double> >& data_temperature,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const std::vector<const double*>& thermo_properties,
            const hier::Box& domain) const;
        
        /*
         * Compute the temperature.
         */
        void
        computeTemperature(
            boost::shared_ptr<pdat::SideData<double> >& data_temperature,
            const boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
            const std::vector<const double*>& thermo_properties,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the temperature.
         */
        void
        computeTemperature(
            boost::shared_ptr<pdat::CellData<double> >& data_temperature,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties,
            const hier::Box& domain) const;
        
        /*
         * Compute the temperature.
         */
        void
        computeTemperature(
            boost::shared_ptr<pdat::SideData<double> >& data_temperature,
            const boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
            const boost::shared_ptr<pdat::SideData<double> >& data_thermo_properties,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the specific internal energy from temperature.
         */
        double
        getInternalEnergyFromTemperature(
            const double* const density,
            const double* const temperature,
            const std::vector<const double*>& thermo_properties) const;
        
        /*
         * Compute the specific internal energy from temperature.
         */
        void
        computeInternalEnergyFromTemperature(
            boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
            const std::vector<const double*>& thermo_properties,
            const hier::Box& domain) const;
        
        /*
         * Compute the specific internal energy from temperature.
         */
        void
        computeInternalEnergyFromTemperature(
            boost::shared_ptr<pdat::SideData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_temperature,
            const std::vector<const double*>& thermo_properties,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the specific internal energy from temperature.
         */
        void
        computeInternalEnergyFromTemperature(
            boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
            const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties,
            const hier::Box& domain) const;
        
        /*
         * Compute the specific internal energy from temperature.
         */
        void
        computeInternalEnergyFromTemperature(
            boost::shared_ptr<pdat::SideData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_temperature,
            const boost::shared_ptr<pdat::SideData<double> >& data_thermo_properties,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the isochoric specific heat capacity.
         */
        double
        getIsochoricSpecificHeatCapacity(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& thermo_properties) const;
        
        /*
         * Compute the isochoric specific heat capacity.
         */
        void
        computeIsochoricSpecificHeatCapacity(
            boost::shared_ptr<pdat::CellData<double> >& data_isochoric_specific_heat_capacity,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const std::vector<const double*>& thermo_properties,
            const hier::Box& domain) const;
        
        /*
         * Compute the isochoric specific heat capacity.
         */
        void
        computeIsochoricSpecificHeatCapacity(
            boost::shared_ptr<pdat::SideData<double> >& data_isochoric_specific_heat_capacity,
            const boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
            const std::vector<const double*>& thermo_properties,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the isochoric specific heat capacity.
         */
        void
        computeIsochoricSpecificHeatCapacity(
            boost::shared_ptr<pdat::CellData<double> >& data_isochoric_specific_heat_capacity,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties,
            const hier::Box& domain) const;
        
        /*
         * Compute the isochoric specific heat capacity.
         */
        void
        computeIsochoricSpecificHeatCapacity(
            boost::shared_ptr<pdat::SideData<double> >& data_isochoric_specific_heat_capacity,
            const boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
            const boost::shared_ptr<pdat::SideData<double> >& data_thermo_properties,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the isobaric specific heat capacity.
         */
        double
        getIsobaricSpecificHeatCapacity(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& thermo_properties) const;
        
        /*
         * Compute the isobaric specific heat capacity.
         */
        void
        computeIsobaricSpecificHeatCapacity(
            boost::shared_ptr<pdat::CellData<double> >& data_isobaric_specific_heat_capacity,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const std::vector<const double*>& thermo_properties,
            const hier::Box& domain) const;
        
        /*
         * Compute the isobaric specific heat capacity.
         */
        void
        computeIsobaricSpecificHeatCapacity(
            boost::shared_ptr<pdat::SideData<double> >& data_isobaric_specific_heat_capacity,
            const boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
            const std::vector<const double*>& thermo_properties,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the isobaric specific heat capacity.
         */
        void
        computeIsobaricSpecificHeatCapacity(
            boost::shared_ptr<pdat::CellData<double> >& data_isobaric_specific_heat_capacity,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties,
            const hier::Box& domain) const;
        
        /*
         * Compute the isobaric specific heat capacity.
         */
        void
        computeIsobaricSpecificHeatCapacity(
            boost::shared_ptr<pdat::SideData<double> >& data_isobaric_specific_heat_capacity,
            const boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
            const boost::shared_ptr<pdat::SideData<double> >& data_thermo_properties,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the partial derivative of internal energy w.r.t. pressure under constant density.
         */
        double
        getIsochoricPartialInternalEnergyPartialPressure(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& thermo_properties) const;
        
        /*
         * Compute the partial derivative of internal energy w.r.t. pressure under constant density.
         */
        void
        computeIsochoricPartialInternalEnergyPartialPressure(
            boost::shared_ptr<pdat::CellData<double> >& data_partial_internal_energy_partial_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const std::vector<const double*>& thermo_properties,
            const hier::Box& domain) const;
        
        /*
         * Compute the partial derivative of internal energy w.r.t. pressure under constant density.
         */
        void
        computeIsochoricPartialInternalEnergyPartialPressure(
            boost::shared_ptr<pdat::SideData<double> >& data_partial_internal_energy_partial_pressure,
            const boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
            const std::vector<const double*>& thermo_properties,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the partial derivative of internal energy w.r.t. pressure under constant density.
         */
        void
        computeIsochoricPartialInternalEnergyPartialPressure(
            boost::shared_ptr<pdat::CellData<double> >& data_partial_internal_energy_partial_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties,
            const hier::Box& domain) const;
        
        /*
         * Compute the partial derivative of internal energy w.r.t. pressure under constant density.
         */
        void
        computeIsochoricPartialInternalEnergyPartialPressure(
            boost::shared_ptr<pdat::SideData<double> >& data_partial_internal_energy_partial_pressure,
            const boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
            const boost::shared_ptr<pdat::SideData<double> >& data_thermo_properties,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the partial derivative of internal energy w.r.t. density under constant pressure.
         */
        double
        getIsobaricPartialInternalEnergyPartialDensity(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& thermo_properties) const;
        
        /*
         * Compute the partial derivative of internal energy w.r.t. density under constant pressure.
         */
        void
        computeIsobaricPartialInternalEnergyPartialDensity(
            boost::shared_ptr<pdat::CellData<double> >& data_partial_internal_energy_partial_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const std::vector<const double*>& thermo_properties,
            const hier::Box& domain) const;
        
        /*
         * Compute the partial derivative of internal energy w.r.t. density under constant pressure.
         */
        void
        computeIsobaricPartialInternalEnergyPartialDensity(
            boost::shared_ptr<pdat::SideData<double> >& data_partial_internal_energy_partial_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
            const std::vector<const double*>& thermo_properties,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the partial derivative of internal energy w.r.t. density under constant pressure.
         */
        void
        computeIsobaricPartialInternalEnergyPartialDensity(
            boost::shared_ptr<pdat::CellData<double> >& data_partial_internal_energy_partial_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties,
            const hier::Box& domain) const;
        
        /*
         * Compute the partial derivative of internal energy w.r.t. density under constant pressure.
         */
        void
        computeIsobaricPartialInternalEnergyPartialDensity(
            boost::shared_ptr<pdat::SideData<double> >& data_partial_internal_energy_partial_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
            const boost::shared_ptr<pdat::SideData<double> >& data_thermo_properties,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the density.
         */
        double
        getDensity(
            const double* const pressure,
            const double* const temperature,
            const std::vector<const double*>& thermo_properties) const;
        
        /*
         * Compute the density.
         */
        void
        computeDensity(
            boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
            const std::vector<const double*>& thermo_properties,
            const hier::Box& domain) const;
        
        /*
         * Compute the density.
         */
        void
        computeDensity(
            boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
            const boost::shared_ptr<pdat::SideData<double> >& data_temperature,
            const std::vector<const double*>& thermo_properties,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the density.
         */
        void
        computeDensity(
            boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
            const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties,
            const hier::Box& domain) const;
        
        /*
         * Compute the density.
         */
        void
        computeDensity(
            boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
            const boost::shared_ptr<pdat::SideData<double> >& data_temperature,
            const boost::shared_ptr<pdat::SideData<double> >& data_thermo_properties,
            int side_normal,
            const hier::Box& domain) const;
        
    private:
        
};

#endif /* EQUATION_OF_STATE_IDEAL_GAS_HPP */
