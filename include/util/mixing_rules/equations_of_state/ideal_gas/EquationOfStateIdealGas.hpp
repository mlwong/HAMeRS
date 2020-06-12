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
         * Compute the Gruneisen parameter (partial derivative of pressure w.r.t. specific internal energy under
         * constant density divided by density).
         */
        double
        getGruneisenParameter(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& thermo_properties) const;
        
        /*
         * Compute the Gruneisen parameter (partial derivative of pressure w.r.t. specific internal energy under
         * constant density divided by density).
         */
        void
        computeGruneisenParameter(
            boost::shared_ptr<pdat::CellData<double> >& data_gruneisen_parameter,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const std::vector<const double*>& thermo_properties,
            const hier::Box& domain) const;
        
        /*
         * Compute the Gruneisen parameter (partial derivative of pressure w.r.t. specific internal energy under
         * constant density divided by density).
         */
        void
        computeGruneisenParameter(
            boost::shared_ptr<pdat::SideData<double> >& data_gruneisen_parameter,
            const boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
            const std::vector<const double*>& thermo_properties,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the Gruneisen parameter (partial derivative of pressure w.r.t. specific internal energy under
         * constant density divided by density).
         */
        void
        computeGruneisenParameter(
            boost::shared_ptr<pdat::CellData<double> >& data_gruneisen_parameter,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties,
            const hier::Box& domain) const;
        
        /*
         * Compute the Gruneisen parameter (partial derivative of pressure w.r.t. specific internal energy under
         * constant density divided by density).
         */
        void
        computeGruneisenParameter(
            boost::shared_ptr<pdat::SideData<double> >& data_gruneisen_parameter,
            const boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
            const boost::shared_ptr<pdat::SideData<double> >& data_thermo_properties,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the partial derivative of pressure w.r.t. density under constant specific internal energy.
         */
        double
        getPressureDerivativeWithDensity(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& thermo_properties) const;
        
        /*
         * Compute the partial derivative of pressure w.r.t. density under constant specific internal energy.
         */
        void
        computePressureDerivativeWithDensity(
            boost::shared_ptr<pdat::CellData<double> >& data_partial_pressure_partial_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const std::vector<const double*>& thermo_properties,
            const hier::Box& domain) const;
        
        /*
         * Compute the partial derivative of pressure w.r.t. density under constant specific internal energy.
         */
        void
        computePressureDerivativeWithDensity(
            boost::shared_ptr<pdat::SideData<double> >& data_partial_pressure_partial_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
            const std::vector<const double*>& thermo_properties,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Compute the partial derivative of pressure w.r.t. density under constant specific internal energy.
         */
        void
        computePressureDerivativeWithDensity(
            boost::shared_ptr<pdat::CellData<double> >& data_partial_pressure_partial_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties,
            const hier::Box& domain) const;
        
        /*
         * Compute the partial derivative of pressure w.r.t. density under constant specific internal energy.
         */
        void
        computePressureDerivativeWithDensity(
            boost::shared_ptr<pdat::SideData<double> >& data_partial_pressure_partial_density,
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
        /*
         * Compute the pressure.
         */
        void
        computePressure(
            double* const p,
            const double* const rho,
            const double* const epsilon,
            const double& gamma,
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
            double* const p,
            const double* const rho,
            const double* const epsilon,
            const double* const gamma,
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
            double* const c,
            const double* const rho,
            const double* const p,
            const double& gamma,
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
            double* const c,
            const double* const rho,
            const double* const p,
            const double* const gamma,
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
            double* const epsilon,
            const double* const rho,
            const double* const p,
            const double& gamma,
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
            double* const epsilon,
            const double* const rho,
            const double* const p,
            const double* const gamma,
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
            double* const h,
            const double* const rho,
            const double* const p,
            const double& gamma,
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
            double* const h,
            const double* const rho,
            const double* const p,
            const double* const gamma,
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
            double* const T,
            const double* const rho,
            const double* const p,
            const double& gamma,
            const double& c_v,
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
            double* const T,
            const double* const rho,
            const double* const p,
            const double* const gamma,
            const double* const c_v,
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
            double* const epsilon,
            const double* const T,
            const double& c_v,
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
            double* const epsilon,
            const double* const T,
            const double* const c_v,
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
            double* const c_v,
            const double& c_v_src,
            const hier::IntVector& offset_isochoric_specific_heat_capacity,
            const hier::IntVector& ghostcell_dims_isochoric_specific_heat_capacity,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims) const;
        
        /*
         * Compute the isochoric specific heat capacity.
         */
        void
        computeIsochoricSpecificHeatCapacity(
            double* const c_v,
            const double* const c_v_src,
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
            double* const c_p,
            const double& c_p_src,
            const hier::IntVector& offset_isobaric_specific_heat_capacity,
            const hier::IntVector& ghostcell_dims_isobaric_specific_heat_capacity,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims) const;
        
        /*
         * Compute the isobaric specific heat capacity.
         */
        void
        computeIsobaricSpecificHeatCapacity(
            double* const c_p,
            const double* const c_p_src,
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
            double* const Gamma,
            const double& gamma,
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
            double* const Gamma,
            const double* const gamma,
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
            double* const Psi,
            const double* const rho,
            const double* const p,
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
            double* const rho,
            const double* const p,
            const double* const T,
            const double& R,
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
            double* const rho,
            const double* const p,
            const double* const T,
            const double* const R,
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
