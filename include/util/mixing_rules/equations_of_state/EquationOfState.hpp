#ifndef EQUATION_OF_STATE_HPP
#define EQUATION_OF_STATE_HPP

#include "HAMeRS_config.hpp"

#include "HAMeRS_memory.hpp"

#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/SideData.h"

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
        
        virtual ~EquationOfState() {}
        
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
         * Compute the pressure.
         */
        void
        computePressure(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_internal_energy,
            const std::vector<const double*>& thermo_properties) const
        {
            const hier::Box empty_box(d_dim);
            computePressure(
                data_pressure,
                data_density,
                data_internal_energy,
                thermo_properties,
                empty_box);
        }
        
        /*
         * Compute the pressure.
         */
        virtual void
        computePressure(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_internal_energy,
            const std::vector<const double*>& thermo_properties,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the pressure.
         */
        void
        computePressure(
            HAMERS_SHARED_PTR<pdat::SideData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_internal_energy,
            const std::vector<const double*>& thermo_properties,
            int side_normal) const
        {
            const hier::Box empty_box(d_dim);
            computePressure(
                data_pressure,
                data_density,
                data_internal_energy,
                thermo_properties,
                side_normal,
                empty_box);
        }
        
        /*
         * Compute the pressure.
         */
        virtual void
        computePressure(
            HAMERS_SHARED_PTR<pdat::SideData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_internal_energy,
            const std::vector<const double*>& thermo_properties,
            int side_normal,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the pressure.
         */
        void
        computePressure(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_internal_energy,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_thermo_properties) const
        {
            const hier::Box empty_box(d_dim);
            computePressure(
                data_pressure,
                data_density,
                data_internal_energy,
                data_thermo_properties,
                empty_box);
        }
        
        /*
         * Compute the pressure.
         */
        virtual void
        computePressure(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_internal_energy,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_thermo_properties,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the pressure.
         */
        void
        computePressure(
            HAMERS_SHARED_PTR<pdat::SideData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_internal_energy,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_thermo_properties,
            int side_normal) const
        {
            const hier::Box empty_box(d_dim);
            computePressure(
                data_pressure,
                data_density,
                data_internal_energy,
                data_thermo_properties,
                side_normal,
                empty_box);
        }
        
        /*
         * Compute the pressure.
         */
        virtual void
        computePressure(
            HAMERS_SHARED_PTR<pdat::SideData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_internal_energy,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_thermo_properties,
            int side_normal,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the sound speed.
         */
        virtual double
        getSoundSpeed(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& thermo_properties) const = 0;
        
        /*
         * Compute the sound speed.
         */
        void
        computeSoundSpeed(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& data_sound_speed,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
            const std::vector<const double*>& thermo_properties) const
        {
            const hier::Box empty_box(d_dim);
            computeSoundSpeed(
                data_sound_speed,
                data_density,
                data_pressure,
                thermo_properties,
                empty_box);
        }
        
        /*
         * Compute the sound speed.
         */
        virtual void
        computeSoundSpeed(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& data_sound_speed,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
            const std::vector<const double*>& thermo_properties,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the sound speed.
         */
        void
        computeSoundSpeed(
            HAMERS_SHARED_PTR<pdat::SideData<double> >& data_sound_speed,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_pressure,
            const std::vector<const double*>& thermo_properties,
            int side_normal) const
        {
            const hier::Box empty_box(d_dim);
            computeSoundSpeed(
                data_sound_speed,
                data_density,
                data_pressure,
                thermo_properties,
                side_normal,
                empty_box);
        }
        
        /*
         * Compute the sound speed.
         */
        virtual void
        computeSoundSpeed(
            HAMERS_SHARED_PTR<pdat::SideData<double> >& data_sound_speed,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_pressure,
            const std::vector<const double*>& thermo_properties,
            int side_normal,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the sound speed.
         */
        void
        computeSoundSpeed(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& data_sound_speed,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_thermo_properties) const
        {
            const hier::Box empty_box(d_dim);
            computeSoundSpeed(
                data_sound_speed,
                data_density,
                data_pressure,
                data_thermo_properties,
                empty_box);
        }
        
        /*
         * Compute the sound speed.
         */
        virtual void
        computeSoundSpeed(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& data_sound_speed,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_thermo_properties,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the sound speed.
         */
        void
        computeSoundSpeed(
            HAMERS_SHARED_PTR<pdat::SideData<double> >& data_sound_speed,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_thermo_properties,
            int side_normal) const
        {
            const hier::Box empty_box(d_dim);
            computeSoundSpeed(
                data_sound_speed,
                data_density,
                data_pressure,
                data_thermo_properties,
                side_normal,
                empty_box);
        }
        
        /*
         * Compute the sound speed.
         */
        virtual void
        computeSoundSpeed(
            HAMERS_SHARED_PTR<pdat::SideData<double> >& data_sound_speed,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_thermo_properties,
            int side_normal,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the specific internal energy.
         */
        virtual double
        getInternalEnergy(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& thermo_properties) const = 0;
        
        /*
         * Compute the specific internal energy.
         */
        void
        computeInternalEnergy(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& data_internal_energy,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
            const std::vector<const double*>& thermo_properties) const
        {
            const hier::Box empty_box(d_dim);
            computeInternalEnergy(
                data_internal_energy,
                data_density,
                data_pressure,
                thermo_properties,
                empty_box);
        }
        
        /*
         * Compute the specific internal energy.
         */
        virtual void
        computeInternalEnergy(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& data_internal_energy,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
            const std::vector<const double*>& thermo_properties,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the specific internal energy.
         */
        void
        computeInternalEnergy(
            HAMERS_SHARED_PTR<pdat::SideData<double> >& data_internal_energy,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_pressure,
            const std::vector<const double*>& thermo_properties,
            int side_normal) const
        {
            const hier::Box empty_box(d_dim);
            computeInternalEnergy(
                data_internal_energy,
                data_density,
                data_pressure,
                thermo_properties,
                side_normal,
                empty_box);
        }
        
        /*
         * Compute the specific internal energy.
         */
        virtual void
        computeInternalEnergy(
            HAMERS_SHARED_PTR<pdat::SideData<double> >& data_internal_energy,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_pressure,
            const std::vector<const double*>& thermo_properties,
            int side_normal,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the specific internal energy.
         */
        void
        computeInternalEnergy(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& data_internal_energy,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_thermo_properties) const
        {
            const hier::Box empty_box(d_dim);
            computeInternalEnergy(
                data_internal_energy,
                data_density,
                data_pressure,
                data_thermo_properties,
                empty_box);
        }
        
        /*
         * Compute the specific internal energy.
         */
        virtual void
        computeInternalEnergy(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& data_internal_energy,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_thermo_properties,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the specific internal energy.
         */
        void
        computeInternalEnergy(
            HAMERS_SHARED_PTR<pdat::SideData<double> >& data_internal_energy,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_thermo_properties,
            int side_normal) const
        {
            const hier::Box empty_box(d_dim);
            computeInternalEnergy(
                data_internal_energy,
                data_density,
                data_pressure,
                data_thermo_properties,
                side_normal,
                empty_box);
        }
        
        /*
         * Compute the specific internal energy.
         */
        virtual void
        computeInternalEnergy(
            HAMERS_SHARED_PTR<pdat::SideData<double> >& data_internal_energy,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_thermo_properties,
            int side_normal,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the specific enthalpy.
         */
        virtual double
        getEnthalpy(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& thermo_properties) const = 0;
        
        /*
         * Compute the specific enthalpy.
         */
        void
        computeEnthalpy(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& data_enthalpy,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
            const std::vector<const double*>& thermo_properties) const
        {
            const hier::Box empty_box(d_dim);
            computeEnthalpy(
                data_enthalpy,
                data_density,
                data_pressure,
                thermo_properties,
                empty_box);
        }
        
        /*
         * Compute the specific enthalpy.
         */
        virtual void
        computeEnthalpy(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& data_enthalpy,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
            const std::vector<const double*>& thermo_properties,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the specific enthalpy.
         */
        void
        computeEnthalpy(
            HAMERS_SHARED_PTR<pdat::SideData<double> >& data_enthalpy,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_pressure,
            const std::vector<const double*>& thermo_properties,
            int side_normal) const
        {
            const hier::Box empty_box(d_dim);
            computeEnthalpy(
                data_enthalpy,
                data_density,
                data_pressure,
                thermo_properties,
                side_normal,
                empty_box);
        }
        
        /*
         * Compute the specific enthalpy.
         */
        virtual void
        computeEnthalpy(
            HAMERS_SHARED_PTR<pdat::SideData<double> >& data_enthalpy,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_pressure,
            const std::vector<const double*>& thermo_properties,
            int side_normal,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the specific enthalpy.
         */
        void
        computeEnthalpy(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& data_enthalpy,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_thermo_properties) const
        {
            const hier::Box empty_box(d_dim);
            computeEnthalpy(
                data_enthalpy,
                data_density,
                data_pressure,
                data_thermo_properties,
                empty_box);
        }
        
        /*
         * Compute the specific enthalpy.
         */
        virtual void
        computeEnthalpy(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& data_enthalpy,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_thermo_properties,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the specific enthalpy.
         */
        void
        computeEnthalpy(
            HAMERS_SHARED_PTR<pdat::SideData<double> >& data_enthalpy,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_thermo_properties,
            int side_normal) const
        {
            const hier::Box empty_box(d_dim);
            computeEnthalpy(
                data_enthalpy,
                data_density,
                data_pressure,
                data_thermo_properties,
                side_normal,
                empty_box);
        }
        
        /*
         * Compute the specific enthalpy.
         */
        virtual void
        computeEnthalpy(
            HAMERS_SHARED_PTR<pdat::SideData<double> >& data_enthalpy,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_thermo_properties,
            int side_normal,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the temperature.
         */
        virtual double
        getTemperature(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& thermo_properties) const = 0;
        
        /*
         * Compute the temperature.
         */
        void
        computeTemperature(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& data_temperature,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
            const std::vector<const double*>& thermo_properties) const
        {
            const hier::Box empty_box(d_dim);
            computeTemperature(
                data_temperature,
                data_density,
                data_pressure,
                thermo_properties,
                empty_box);
        }
        
        /*
         * Compute the temperature.
         */
        virtual void
        computeTemperature(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& data_temperature,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
            const std::vector<const double*>& thermo_properties,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the temperature.
         */
        void
        computeTemperature(
            HAMERS_SHARED_PTR<pdat::SideData<double> >& data_temperature,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_pressure,
            const std::vector<const double*>& thermo_properties,
            int side_normal) const
        {
            const hier::Box empty_box(d_dim);
            computeTemperature(
                data_temperature,
                data_density,
                data_pressure,
                thermo_properties,
                side_normal,
                empty_box);
        }
        
        /*
         * Compute the temperature.
         */
        virtual void
        computeTemperature(
            HAMERS_SHARED_PTR<pdat::SideData<double> >& data_temperature,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_pressure,
            const std::vector<const double*>& thermo_properties,
            int side_normal,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the temperature.
         */
        void
        computeTemperature(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& data_temperature,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_thermo_properties) const
        {
            const hier::Box empty_box(d_dim);
            computeTemperature(
                data_temperature,
                data_density,
                data_pressure,
                data_thermo_properties,
                empty_box);
        }
        
        /*
         * Compute the temperature.
         */
        virtual void
        computeTemperature(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& data_temperature,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_thermo_properties,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the temperature.
         */
        void
        computeTemperature(
            HAMERS_SHARED_PTR<pdat::SideData<double> >& data_temperature,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_thermo_properties,
            int side_normal) const
        {
            const hier::Box empty_box(d_dim);
            computeTemperature(
                data_temperature,
                data_density,
                data_pressure,
                data_thermo_properties,
                side_normal,
                empty_box);
        }
        
        /*
         * Compute the temperature.
         */
        virtual void
        computeTemperature(
            HAMERS_SHARED_PTR<pdat::SideData<double> >& data_temperature,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_thermo_properties,
            int side_normal,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the specific internal energy from temperature.
         */
        virtual double
        getInternalEnergyFromTemperature(
            const double* const density,
            const double* const temperature,
            const std::vector<const double*>& thermo_properties) const = 0;
        
        /*
         * Compute the specific internal energy from temperature.
         */
        void
        computeInternalEnergyFromTemperature(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& data_internal_energy,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_temperature,
            const std::vector<const double*>& thermo_properties) const
        {
            const hier::Box empty_box(d_dim);
            computeInternalEnergyFromTemperature(
                data_internal_energy,
                data_density,
                data_temperature,
                thermo_properties,
                empty_box);
        }
        
        /*
         * Compute the specific internal energy from temperature.
         */
        virtual void
        computeInternalEnergyFromTemperature(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& data_internal_energy,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_temperature,
            const std::vector<const double*>& thermo_properties,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the specific internal energy from temperature.
         */
        void
        computeInternalEnergyFromTemperature(
            HAMERS_SHARED_PTR<pdat::SideData<double> >& data_internal_energy,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_temperature,
            const std::vector<const double*>& thermo_properties,
            int side_normal) const
        {
            const hier::Box empty_box(d_dim);
            computeInternalEnergyFromTemperature(
                data_internal_energy,
                data_density,
                data_temperature,
                thermo_properties,
                side_normal,
                empty_box);
        }
        
        /*
         * Compute the specific internal energy from temperature.
         */
        virtual void
        computeInternalEnergyFromTemperature(
            HAMERS_SHARED_PTR<pdat::SideData<double> >& data_internal_energy,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_temperature,
            const std::vector<const double*>& thermo_properties,
            int side_normal,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the specific internal energy from temperature.
         */
        void
        computeInternalEnergyFromTemperature(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& data_internal_energy,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_temperature,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_thermo_properties) const
        {
            const hier::Box empty_box(d_dim);
            computeInternalEnergyFromTemperature(
                data_internal_energy,
                data_density,
                data_temperature,
                data_thermo_properties,
                empty_box);
        }
        
        /*
         * Compute the specific internal energy from temperature.
         */
        virtual void
        computeInternalEnergyFromTemperature(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& data_internal_energy,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_temperature,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_thermo_properties,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the specific internal energy from temperature.
         */
        void
        computeInternalEnergyFromTemperature(
            HAMERS_SHARED_PTR<pdat::SideData<double> >& data_internal_energy,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_temperature,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_thermo_properties,
            int side_normal) const
        {
            const hier::Box empty_box(d_dim);
            computeInternalEnergyFromTemperature(
                data_internal_energy,
                data_density,
                data_temperature,
                data_thermo_properties,
                side_normal,
                empty_box);
        }
        
        /*
         * Compute the specific internal energy from temperature.
         */
        virtual void
        computeInternalEnergyFromTemperature(
            HAMERS_SHARED_PTR<pdat::SideData<double> >& data_internal_energy,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_temperature,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_thermo_properties,
            int side_normal,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the isochoric specific heat capacity.
         */
        virtual double
        getIsochoricSpecificHeatCapacity(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& thermo_properties) const = 0;
        
        /*
         * Compute the isochoric specific heat capacity.
         */
        void
        computeIsochoricSpecificHeatCapacity(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& data_isochoric_specific_heat_capacity,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
            const std::vector<const double*>& thermo_properties) const
        {
            const hier::Box empty_box(d_dim);
            computeIsochoricSpecificHeatCapacity(
                data_isochoric_specific_heat_capacity,
                data_density,
                data_pressure,
                thermo_properties,
                empty_box);
        }
        
        /*
         * Compute the isochoric specific heat capacity.
         */
        virtual void
        computeIsochoricSpecificHeatCapacity(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& data_isochoric_specific_heat_capacity,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
            const std::vector<const double*>& thermo_properties,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the isochoric specific heat capacity.
         */
        void
        computeIsochoricSpecificHeatCapacity(
            HAMERS_SHARED_PTR<pdat::SideData<double> >& data_isochoric_specific_heat_capacity,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_pressure,
            const std::vector<const double*>& thermo_properties,
            int side_normal) const
        {
            const hier::Box empty_box(d_dim);
            computeIsochoricSpecificHeatCapacity(
                data_isochoric_specific_heat_capacity,
                data_density,
                data_pressure,
                thermo_properties,
                side_normal,
                empty_box);
        }
        
        /*
         * Compute the isochoric specific heat capacity.
         */
        virtual void
        computeIsochoricSpecificHeatCapacity(
            HAMERS_SHARED_PTR<pdat::SideData<double> >& data_isochoric_specific_heat_capacity,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_pressure,
            const std::vector<const double*>& thermo_properties,
            int side_normal,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the isochoric specific heat capacity.
         */
        void
        computeIsochoricSpecificHeatCapacity(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& data_isochoric_specific_heat_capacity,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_thermo_properties) const
        {
            const hier::Box empty_box(d_dim);
            computeIsochoricSpecificHeatCapacity(
                data_isochoric_specific_heat_capacity,
                data_density,
                data_pressure,
                data_thermo_properties,
                empty_box);
        }
        
        /*
         * Compute the isochoric specific heat capacity.
         */
        virtual void
        computeIsochoricSpecificHeatCapacity(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& data_isochoric_specific_heat_capacity,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_thermo_properties,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the isochoric specific heat capacity.
         */
        void
        computeIsochoricSpecificHeatCapacity(
            HAMERS_SHARED_PTR<pdat::SideData<double> >& data_isochoric_specific_heat_capacity,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_thermo_properties,
            int side_normal) const
        {
            const hier::Box empty_box(d_dim);
            computeIsochoricSpecificHeatCapacity(
                data_isochoric_specific_heat_capacity,
                data_density,
                data_pressure,
                data_thermo_properties,
                side_normal,
                empty_box);
        }
        
        /*
         * Compute the isochoric specific heat capacity.
         */
        virtual void
        computeIsochoricSpecificHeatCapacity(
            HAMERS_SHARED_PTR<pdat::SideData<double> >& data_isochoric_specific_heat_capacity,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_thermo_properties,
            int side_normal,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the isobaric specific heat capacity.
         */
        virtual double
        getIsobaricSpecificHeatCapacity(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& thermo_properties) const = 0;
        
        /*
         * Compute the isobaric specific heat capacity.
         */
        void
        computeIsobaricSpecificHeatCapacity(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& data_isobaric_specific_heat_capacity,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
            const std::vector<const double*>& thermo_properties) const
        {
            const hier::Box empty_box(d_dim);
            computeIsobaricSpecificHeatCapacity(
                data_isobaric_specific_heat_capacity,
                data_density,
                data_pressure,
                thermo_properties,
                empty_box);
        }
        
        /*
         * Compute the isobaric specific heat capacity.
         */
        virtual void
        computeIsobaricSpecificHeatCapacity(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& data_isobaric_specific_heat_capacity,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
            const std::vector<const double*>& thermo_properties,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the isobaric specific heat capacity.
         */
        void
        computeIsobaricSpecificHeatCapacity(
            HAMERS_SHARED_PTR<pdat::SideData<double> >& data_isobaric_specific_heat_capacity,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_pressure,
            const std::vector<const double*>& thermo_properties,
            int side_normal) const
        {
            const hier::Box empty_box(d_dim);
            computeIsobaricSpecificHeatCapacity(
                data_isobaric_specific_heat_capacity,
                data_density,
                data_pressure,
                thermo_properties,
                side_normal,
                empty_box);
        }
        
        /*
         * Compute the isobaric specific heat capacity.
         */
        virtual void
        computeIsobaricSpecificHeatCapacity(
            HAMERS_SHARED_PTR<pdat::SideData<double> >& data_isobaric_specific_heat_capacity,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_pressure,
            const std::vector<const double*>& thermo_properties,
            int side_normal,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the isobaric specific heat capacity.
         */
        void
        computeIsobaricSpecificHeatCapacity(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& data_isobaric_specific_heat_capacity,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_thermo_properties) const
        {
            const hier::Box empty_box(d_dim);
            computeIsobaricSpecificHeatCapacity(
                data_isobaric_specific_heat_capacity,
                data_density,
                data_pressure,
                data_thermo_properties,
                empty_box);
        }
        
        /*
         * Compute the isobaric specific heat capacity.
         */
        virtual void
        computeIsobaricSpecificHeatCapacity(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& data_isobaric_specific_heat_capacity,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_thermo_properties,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the isobaric specific heat capacity.
         */
        void
        computeIsobaricSpecificHeatCapacity(
            HAMERS_SHARED_PTR<pdat::SideData<double> >& data_isobaric_specific_heat_capacity,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_thermo_properties,
            int side_normal) const
        {
            const hier::Box empty_box(d_dim);
            computeIsobaricSpecificHeatCapacity(
                data_isobaric_specific_heat_capacity,
                data_density,
                data_pressure,
                data_thermo_properties,
                side_normal,
                empty_box);
        }
        
        /*
         * Compute the isobaric specific heat capacity.
         */
        virtual void
        computeIsobaricSpecificHeatCapacity(
            HAMERS_SHARED_PTR<pdat::SideData<double> >& data_isobaric_specific_heat_capacity,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_thermo_properties,
            int side_normal,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the Gruneisen parameter (partial derivative of pressure w.r.t. specific internal energy under
         * constant density divided by density).
         */
        virtual double
        getGruneisenParameter(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& thermo_properties) const = 0;
        
        /*
         * Compute the Gruneisen parameter (partial derivative of pressure w.r.t. specific internal energy under
         * constant density divided by density).
         */
        void
        computeGruneisenParameter(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& data_gruneisen_parameter,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
            const std::vector<const double*>& thermo_properties) const
        {
            const hier::Box empty_box(d_dim);
            computeGruneisenParameter(
                data_gruneisen_parameter,
                data_density,
                data_pressure,
                thermo_properties,
                empty_box);
        }
        
        /*
         * Compute the Gruneisen parameter (partial derivative of pressure w.r.t. specific internal energy under
         * constant density divided by density).
         */
        virtual void
        computeGruneisenParameter(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& data_gruneisen_parameter,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
            const std::vector<const double*>& thermo_properties,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the Gruneisen parameter (partial derivative of pressure w.r.t. specific internal energy under
         * constant density divided by density).
         */
        void
        computeGruneisenParameter(
            HAMERS_SHARED_PTR<pdat::SideData<double> >& data_gruneisen_parameter,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_pressure,
            const std::vector<const double*>& thermo_properties,
            int side_normal) const
        {
            const hier::Box empty_box(d_dim);
            computeGruneisenParameter(
                data_gruneisen_parameter,
                data_density,
                data_pressure,
                thermo_properties,
                side_normal,
                empty_box);
        }
        
        /*
         * Compute the Gruneisen parameter (partial derivative of pressure w.r.t. specific internal energy under
         * constant density divided by density).
         */
        virtual void
        computeGruneisenParameter(
            HAMERS_SHARED_PTR<pdat::SideData<double> >& data_gruneisen_parameter,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_pressure,
            const std::vector<const double*>& thermo_properties,
            int side_normal,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the Gruneisen parameter (partial derivative of pressure w.r.t. specific internal energy under
         * constant density divided by density).
         */
        void
        computeGruneisenParameter(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& data_gruneisen_parameter,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_thermo_properties) const
        {
            const hier::Box empty_box(d_dim);
            computeGruneisenParameter(
                data_gruneisen_parameter,
                data_density,
                data_pressure,
                data_thermo_properties,
                empty_box);
        }
        
        /*
         * Compute the Gruneisen parameter (partial derivative of pressure w.r.t. specific internal energy under
         * constant density divided by density).
         */
        virtual void
        computeGruneisenParameter(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& data_gruneisen_parameter,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_thermo_properties,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the Gruneisen parameter (partial derivative of pressure w.r.t. specific internal energy under
         * constant density divided by density).
         */
        void
        computeGruneisenParameter(
            HAMERS_SHARED_PTR<pdat::SideData<double> >& data_gruneisen_parameter,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_thermo_properties,
            int side_normal) const
        {
            const hier::Box empty_box(d_dim);
            computeGruneisenParameter(
                data_gruneisen_parameter,
                data_density,
                data_pressure,
                data_thermo_properties,
                side_normal,
                empty_box);
        }
        
        /*
         * Compute the Gruneisen parameter (partial derivative of pressure w.r.t. specific internal energy under
         * constant density divided by density).
         */
        virtual void
        computeGruneisenParameter(
            HAMERS_SHARED_PTR<pdat::SideData<double> >& data_gruneisen_parameter,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_thermo_properties,
            int side_normal,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the partial derivative of pressure w.r.t. density under constant specific internal energy.
         */
        virtual double
        getPressureDerivativeWithDensity(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& thermo_properties) const = 0;
        
        /*
         * Compute the partial derivative of pressure w.r.t. density under constant specific internal energy.
         */
        void
        computePressureDerivativeWithDensity(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& data_partial_pressure_partial_density,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
            const std::vector<const double*>& thermo_properties) const
        {
            const hier::Box empty_box(d_dim);
            computePressureDerivativeWithDensity(
                data_partial_pressure_partial_density,
                data_density,
                data_pressure,
                thermo_properties,
                empty_box);
        }
        
        /*
         * Compute the partial derivative of pressure w.r.t. density under constant specific internal energy.
         */
        virtual void
        computePressureDerivativeWithDensity(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& data_partial_pressure_partial_density,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
            const std::vector<const double*>& thermo_properties,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the partial derivative of pressure w.r.t. density under constant specific internal energy.
         */
        void
        computePressureDerivativeWithDensity(
            HAMERS_SHARED_PTR<pdat::SideData<double> >& data_partial_pressure_partial_density,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_pressure,
            const std::vector<const double*>& thermo_properties,
            int side_normal) const
        {
            const hier::Box empty_box(d_dim);
            computePressureDerivativeWithDensity(
                data_partial_pressure_partial_density,
                data_density,
                data_pressure,
                thermo_properties,
                side_normal,
                empty_box);
        }
        
        /*
         * Compute the partial derivative of pressure w.r.t. density under constant specific internal energy.
         */
        virtual void
        computePressureDerivativeWithDensity(
            HAMERS_SHARED_PTR<pdat::SideData<double> >& data_partial_pressure_partial_density,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_pressure,
            const std::vector<const double*>& thermo_properties,
            int side_normal,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the partial derivative of pressure w.r.t. density under constant specific internal energy.
         */
        void
        computePressureDerivativeWithDensity(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& data_partial_pressure_partial_density,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_thermo_properties) const
        {
            const hier::Box empty_box(d_dim);
            computePressureDerivativeWithDensity(
                data_partial_pressure_partial_density,
                data_density,
                data_pressure,
                data_thermo_properties,
                empty_box);
        }
        
        /*
         * Compute the partial derivative of pressure w.r.t. density under constant specific internal energy.
         */
        virtual void
        computePressureDerivativeWithDensity(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& data_partial_pressure_partial_density,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_thermo_properties,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the partial derivative of pressure w.r.t. density under constant specific internal energy.
         */
        void
        computePressureDerivativeWithDensity(
            HAMERS_SHARED_PTR<pdat::SideData<double> >& data_partial_pressure_partial_density,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_thermo_properties,
            int side_normal) const
        {
            const hier::Box empty_box(d_dim);
            computePressureDerivativeWithDensity(
                data_partial_pressure_partial_density,
                data_density,
                data_pressure,
                data_thermo_properties,
                side_normal,
                empty_box);
        }
        
        /*
         * Compute the partial derivative of pressure w.r.t. density under constant specific internal energy.
         */
        virtual void
        computePressureDerivativeWithDensity(
            HAMERS_SHARED_PTR<pdat::SideData<double> >& data_partial_pressure_partial_density,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_thermo_properties,
            int side_normal,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the density.
         */
        virtual double
        getDensity(
            const double* const pressure,
            const double* const temperature,
            const std::vector<const double*>& thermo_properties) const = 0;
        
        /*
         * Compute the density.
         */
        void
        computeDensity(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_temperature,
            const std::vector<const double*>& thermo_properties) const
        {
            const hier::Box empty_box(d_dim);
            computeDensity(
                data_density,
                data_pressure,
                data_temperature,
                thermo_properties,
                empty_box);
        }
        
        /*
         * Compute the density.
         */
        virtual void
        computeDensity(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_temperature,
            const std::vector<const double*>& thermo_properties,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the density.
         */
        void
        computeDensity(
            HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_temperature,
            const std::vector<const double*>& thermo_properties,
            int side_normal) const
        {
            const hier::Box empty_box(d_dim);
            computeDensity(
                data_density,
                data_pressure,
                data_temperature,
                thermo_properties,
                side_normal,
                empty_box);
        }
        
        /*
         * Compute the density.
         */
        virtual void
        computeDensity(
            HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_temperature,
            const std::vector<const double*>& thermo_properties,
            int side_normal,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the density.
         */
        void
        computeDensity(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_temperature,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_thermo_properties) const
        {
            const hier::Box empty_box(d_dim);
            computeDensity(
                data_density,
                data_pressure,
                data_temperature,
                data_thermo_properties,
                empty_box);
        }
        
        /*
         * Compute the density.
         */
        virtual void
        computeDensity(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_temperature,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_thermo_properties,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the density.
         */
        void
        computeDensity(
            HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_temperature,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_thermo_properties,
            int side_normal) const
        {
            const hier::Box empty_box(d_dim);
            computeDensity(
                data_density,
                data_pressure,
                data_temperature,
                data_thermo_properties,
                side_normal,
                empty_box);
        }
        
        /*
         * Compute the density.
         */
        virtual void
        computeDensity(
            HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_temperature,
            const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_thermo_properties,
            int side_normal,
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

#endif /* EQUATION_OF_STATE_HPP */
