#ifndef EQUATION_OF_STATE_HPP
#define EQUATION_OF_STATE_HPP

#include "HAMeRS_config.hpp"

#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/pdat/CellData.h"

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
         * Compute the pressure.
         */
        void
        getPressure(
            boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties) const
        {
            const hier::Box empty_box(d_dim);
            getPressure(
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
        getPressure(
            boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties,
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
        getSoundSpeed(
            boost::shared_ptr<pdat::CellData<double> >& data_sound_speed,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties) const
        {
            const hier::Box empty_box(d_dim);
            getSoundSpeed(
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
        getSoundSpeed(
            boost::shared_ptr<pdat::CellData<double> >& data_sound_speed,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties,
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
        getInternalEnergy(
            boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties) const
        {
            const hier::Box empty_box(d_dim);
            getInternalEnergy(
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
        getInternalEnergy(
            boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties,
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
        getEnthalpy(
            boost::shared_ptr<pdat::CellData<double> >& data_enthalpy,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties) const
        {
            const hier::Box empty_box(d_dim);
            getEnthalpy(
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
        getEnthalpy(
            boost::shared_ptr<pdat::CellData<double> >& data_enthalpy,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties,
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
        getTemperature(
            boost::shared_ptr<pdat::CellData<double> >& data_temperature,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties) const
        {
            const hier::Box empty_box(d_dim);
            getTemperature(
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
        getTemperature(
            boost::shared_ptr<pdat::CellData<double> >& data_temperature,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties,
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
        getInternalEnergyFromTemperature(
            boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
            const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties) const
        {
            const hier::Box empty_box(d_dim);
            getInternalEnergyFromTemperature(
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
        getInternalEnergyFromTemperature(
            boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
            const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties,
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
        getIsochoricSpecificHeatCapacity(
            boost::shared_ptr<pdat::CellData<double> >& data_isochoric_specific_heat_capacity,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties) const
        {
            const hier::Box empty_box(d_dim);
            getIsochoricSpecificHeatCapacity(
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
        getIsochoricSpecificHeatCapacity(
            boost::shared_ptr<pdat::CellData<double> >& data_isochoric_specific_heat_capacity,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties,
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
        getIsobaricSpecificHeatCapacity(
            boost::shared_ptr<pdat::CellData<double> >& data_isobaric_specific_heat_capacity,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties) const
        {
            const hier::Box empty_box(d_dim);
            getIsobaricSpecificHeatCapacity(
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
        getIsobaricSpecificHeatCapacity(
            boost::shared_ptr<pdat::CellData<double> >& data_isobaric_specific_heat_capacity,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the partial derivative of internal energy w.r.t. pressure under constant density.
         */
        virtual double
        getIsochoricPartialInternalEnergyPartialPressure(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& thermo_properties) const = 0;
        
        /*
         * Compute the partial derivative of internal energy w.r.t. pressure under constant density.
         */
        void
        getIsochoricPartialInternalEnergyPartialPressure(
            boost::shared_ptr<pdat::CellData<double> >& data_partial_internal_energy_partial_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties) const
        {
            const hier::Box empty_box(d_dim);
            getIsochoricPartialInternalEnergyPartialPressure(
                data_partial_internal_energy_partial_pressure,
                data_density,
                data_pressure,
                data_thermo_properties,
                empty_box);
        }
        
        /*
         * Compute the partial derivative of internal energy w.r.t. pressure under constant density.
         */
        virtual void
        getIsochoricPartialInternalEnergyPartialPressure(
            boost::shared_ptr<pdat::CellData<double> >& data_partial_internal_energy_partial_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the partial derivative of internal energy w.r.t. density under constant pressure.
         */
        virtual double
        getIsobaricPartialInternalEnergyPartialDensity(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& thermo_properties) const = 0;
        
        /*
         * Compute the partial derivative of internal energy w.r.t. density under constant pressure.
         */
        void
        getIsobaricPartialInternalEnergyPartialDensity(
            boost::shared_ptr<pdat::CellData<double> >& data_partial_internal_energy_partial_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties) const
        {
            const hier::Box empty_box(d_dim);
            getIsobaricPartialInternalEnergyPartialDensity(
                data_partial_internal_energy_partial_density,
                data_density,
                data_pressure,
                data_thermo_properties,
                empty_box);
        }
        
        /*
         * Compute the partial derivative of internal energy w.r.t. density under constant pressure.
         */
        virtual void
        getIsobaricPartialInternalEnergyPartialDensity(
            boost::shared_ptr<pdat::CellData<double> >& data_partial_internal_energy_partial_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties,
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
        getDensity(
            boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
            const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties) const
        {
            const hier::Box empty_box(d_dim);
            getDensity(
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
        getDensity(
            boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
            const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties,
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
