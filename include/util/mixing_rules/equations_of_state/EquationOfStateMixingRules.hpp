#ifndef EQUATION_OF_STATE_MIXING_RULES_HPP
#define EQUATION_OF_STATE_MIXING_RULES_HPP

#include "util/mixing_rules/equations_of_state/EquationOfState.hpp"
#include "util/mixing_rules/MixingClosureModels.hpp"

#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Dimension.h"

#include <string>
#include <vector>

using namespace SAMRAI;

class EquationOfStateMixingRules
{
    public:
        EquationOfStateMixingRules(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const int& num_species,
            const MIXING_CLOSURE_MODEL::TYPE& mixing_closure_model,
            const boost::shared_ptr<tbox::Database>& equation_of_state_mixing_rules_db):
                d_object_name(object_name),
                d_dim(dim),
                d_num_species(num_species),
                d_mixing_closure_model(mixing_closure_model),
                d_equation_of_state_mixing_rules_db(equation_of_state_mixing_rules_db)
        {}
        
        virtual ~EquationOfStateMixingRules() {}
        
        /*
         * Return the boost::shared_ptr to the equation of state.
         */
        virtual const boost::shared_ptr<EquationOfState>&
        getEquationOfState(const int species_index = 0) const = 0;
        
        /*
         * Print all characteristics of the equation of state mixing rules class.
         */
        virtual void
        printClassData(std::ostream& os) const = 0;
        
        /*
         * Put the characteristics of the equation of state mixing rules class into the restart
         * database.
         */
        virtual void
        putToRestart(
            const boost::shared_ptr<tbox::Database>& restart_db) const = 0;
        
        /*
         * Compute the pressure of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        virtual double
        getPressure(
            const double* const density,
            const double* const internal_energy,
            const std::vector<const double*>& mass_fractions) const = 0;
        
        /*
         * Compute the pressure of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        void
        computePressure(
            boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions) const
        {
            const hier::Box empty_box(d_dim);
            computePressure(
                data_pressure,
                data_density,
                data_internal_energy,
                data_mass_fractions,
                empty_box);
        }
        
        /*
         * Compute the pressure of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        virtual void
        computePressure(
            boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the pressure of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        void
        computePressure(
            boost::shared_ptr<pdat::SideData<double> >& data_pressure,
            const boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
            int side_normal) const
        {
            const hier::Box empty_box(d_dim);
            computePressure(
                data_pressure,
                data_density,
                data_internal_energy,
                data_mass_fractions,
                side_normal,
                empty_box);
        }
        
        /*
         * Compute the pressure of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        virtual void
        computePressure(
            boost::shared_ptr<pdat::SideData<double> >& data_pressure,
            const boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
            int side_normal,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the pressure of the mixture with isobaric equilibrium assumption.
         */
        virtual double
        getPressure(
            const double* const density,
            const double* const internal_energy,
            const std::vector<const double*>& mass_fractions,
            const std::vector<const double*>& volume_fractions) const = 0;
        
        /*
         * Compute the pressure of the mixture with isobaric equilibrium assumption.
         */
        void
        computePressure(
            boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
            const boost::shared_ptr<pdat::CellData<double> >& data_volume_fractions) const
        {
            const hier::Box empty_box(d_dim);
            computePressure(
                data_pressure,
                data_density,
                data_internal_energy,
                data_mass_fractions,
                data_volume_fractions,
                empty_box);
        }
        
        /*
         * Compute the pressure of the mixture with isobaric equilibrium assumption.
         */
        virtual void
        computePressure(
            boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
            const boost::shared_ptr<pdat::CellData<double> >& data_volume_fractions,
            const hier::Box& domain) const = 0;
        
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
            int side_normal) const
        {
            const hier::Box empty_box(d_dim);
            computePressure(
                data_pressure,
                data_density,
                data_internal_energy,
                data_mass_fractions,
                data_volume_fractions,
                side_normal,
                empty_box);
        }
        
        /*
         * Compute the pressure of the mixture with isobaric equilibrium assumption.
         */
        virtual void
        computePressure(
            boost::shared_ptr<pdat::SideData<double> >& data_pressure,
            const boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
            const boost::shared_ptr<pdat::SideData<double> >& data_volume_fractions,
            int side_normal,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the specific internal energy of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        virtual double
        getInternalEnergy(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& mass_fractions) const = 0;
        
        /*
         * Compute the specific internal energy of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        void
        computeInternalEnergy(
            boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions) const
        {
            const hier::Box empty_box(d_dim);
            computeInternalEnergy(
                data_internal_energy,
                data_density,
                data_pressure,
                data_mass_fractions,
                empty_box);
        }
        
        /*
         * Compute the specific internal energy of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        virtual void
        computeInternalEnergy(
            boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the specific internal energy of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        void
        computeInternalEnergy(
            boost::shared_ptr<pdat::SideData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
            const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
            int side_normal) const
        {
            const hier::Box empty_box(d_dim);
            computeInternalEnergy(
                data_internal_energy,
                data_density,
                data_pressure,
                data_mass_fractions,
                side_normal,
                empty_box);
        }
        
        /*
         * Compute the specific internal energy of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        virtual void
        computeInternalEnergy(
            boost::shared_ptr<pdat::SideData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
            const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
            int side_normal,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the specific internal energy of the mixture with isobaric equilibrium assumption.
         */
        virtual double
        getInternalEnergy(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& mass_fractions,
            const std::vector<const double*>& volume_fractions) const = 0;
        
        /*
         * Compute the specific internal energy of the mixture with isobaric equilibrium assumption.
         */
        void
        computeInternalEnergy(
            boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
            const boost::shared_ptr<pdat::CellData<double> >& data_volume_fractions) const
        {
            const hier::Box empty_box(d_dim);
            computeInternalEnergy(
                data_internal_energy,
                data_density,
                data_pressure,
                data_mass_fractions,
                data_volume_fractions,
                empty_box);
        }
        
        /*
         * Compute the specific internal energy of the mixture with isobaric equilibrium assumption.
         */
        virtual void
        computeInternalEnergy(
            boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
            const boost::shared_ptr<pdat::CellData<double> >& data_volume_fractions,
            const hier::Box& domain) const = 0;
        
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
            int side_normal) const
        {
            const hier::Box empty_box(d_dim);
            computeInternalEnergy(
                data_internal_energy,
                data_density,
                data_pressure,
                data_mass_fractions,
                data_volume_fractions,
                side_normal,
                empty_box);
        }
        
        /*
         * Compute the specific internal energy of the mixture with isobaric equilibrium assumption.
         */
        virtual void
        computeInternalEnergy(
            boost::shared_ptr<pdat::SideData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
            const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
            const boost::shared_ptr<pdat::SideData<double> >& data_volume_fractions,
            int side_normal,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the temperature of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        virtual double
        getTemperature(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& mass_fractions) const = 0;
        
        /*
         * Compute the temperature of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        void
        computeTemperature(
            boost::shared_ptr<pdat::CellData<double> >& data_temperature,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions) const
        {
            const hier::Box empty_box(d_dim);
            computeTemperature(
                data_temperature,
                data_density,
                data_pressure,
                data_mass_fractions,
                empty_box);
        } 
        
        /*
         * Compute the temperature of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        virtual void
        computeTemperature(
            boost::shared_ptr<pdat::CellData<double> >& data_temperature,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the temperature of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        void
        computeTemperature(
            boost::shared_ptr<pdat::SideData<double> >& data_temperature,
            const boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
            const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
            int side_normal) const
        {
            const hier::Box empty_box(d_dim);
            computeTemperature(
                data_temperature,
                data_density,
                data_pressure,
                data_mass_fractions,
                side_normal,
                empty_box);
        } 
        
        /*
         * Compute the temperature of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        virtual void
        computeTemperature(
            boost::shared_ptr<pdat::SideData<double> >& data_temperature,
            const boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
            const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
            int side_normal,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the specific internal energy of the mixture from temperature with isothermal
         * and isobaric equilibrium assumptions.
         */
        virtual double
        getInternalEnergyFromTemperature(
            const double* const density,
            const double* const temperature,
            const std::vector<const double*>& mass_fractions) const = 0;
        
        /*
         * Compute the specific internal energy of the mixture from temperature with isothermal
         * and isobaric equilibrium assumptions.
         */
        void
        computeInternalEnergyFromTemperature(
            boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions) const
        {
            const hier::Box empty_box(d_dim);
            computeInternalEnergyFromTemperature(
                data_internal_energy,
                data_density,
                data_temperature,
                data_mass_fractions,
                empty_box);
        }
        
        /*
         * Compute the specific internal energy of the mixture from temperature with isothermal
         * and isobaric equilibrium assumptions.
         */
        virtual void
        computeInternalEnergyFromTemperature(
            boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
            const hier::Box& domain) const = 0;
        
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
            int side_normal) const
        {
            const hier::Box empty_box(d_dim);
            computeInternalEnergyFromTemperature(
                data_internal_energy,
                data_density,
                data_temperature,
                data_mass_fractions,
                side_normal,
                empty_box);
        }
        
        /*
         * Compute the specific internal energy of the mixture from temperature with isothermal
         * and isobaric equilibrium assumptions.
         */
        virtual void
        computeInternalEnergyFromTemperature(
            boost::shared_ptr<pdat::SideData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_temperature,
            const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
            int side_normal,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the isochoric specific heat capacity of mixture with isothermal and isobaric
         * equilibrium assumptions.
         */
        virtual double
        getIsochoricSpecificHeatCapacity(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& mass_fractions) const = 0;
        
        /*
         * Compute the isochoric specific heat capacity of mixture with isothermal and isobaric
         * equilibrium assumptions.
         */
        void
        computeIsochoricSpecificHeatCapacity(
            boost::shared_ptr<pdat::CellData<double> >& data_isochoric_specific_heat_capacity,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions) const
        {
            const hier::Box empty_box(d_dim);
            computeIsochoricSpecificHeatCapacity(
                data_isochoric_specific_heat_capacity,
                data_density,
                data_pressure,
                data_mass_fractions,
                empty_box);
        }
        
        /*
         * Compute the isochoric specific heat capacity of mixture with isothermal and isobaric
         * equilibrium assumptions.
         */
        virtual void
        computeIsochoricSpecificHeatCapacity(
            boost::shared_ptr<pdat::CellData<double> >& data_isochoric_specific_heat_capacity,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
            const hier::Box& domain) const = 0;
        
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
            int side_normal) const
        {
            const hier::Box empty_box(d_dim);
            computeIsochoricSpecificHeatCapacity(
                data_isochoric_specific_heat_capacity,
                data_density,
                data_pressure,
                data_mass_fractions,
                side_normal,
                empty_box);
        }
        
        /*
         * Compute the isochoric specific heat capacity of mixture with isothermal and isobaric
         * equilibrium assumptions.
         */
        virtual void
        computeIsochoricSpecificHeatCapacity(
            boost::shared_ptr<pdat::SideData<double> >& data_isochoric_specific_heat_capacity,
            const boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
            const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
            int side_normal,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the isobaric specific heat capacity of mixture with isothermal and isobaric
         * equilibrium assumptions.
         */
        virtual double
        getIsobaricSpecificHeatCapacity(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& mass_fractions) const = 0;
        
        /*
         * Compute the isobaric specific heat capacity of mixture with isothermal and isobaric
         * equilibrium assumptions.
         */
        void
        computeIsobaricSpecificHeatCapacity(
            boost::shared_ptr<pdat::CellData<double> >& data_isobaric_specific_heat_capacity,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions) const
        {
            const hier::Box empty_box(d_dim);
            computeIsobaricSpecificHeatCapacity(
                data_isobaric_specific_heat_capacity,
                data_density,
                data_pressure,
                data_mass_fractions,
                empty_box);
        }
        
        /*
         * Compute the isobaric specific heat capacity of mixture with isothermal and isobaric
         * equilibrium assumptions.
         */
        virtual void
        computeIsobaricSpecificHeatCapacity(
            boost::shared_ptr<pdat::CellData<double> >& data_isobaric_specific_heat_capacity,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
            const hier::Box& domain) const = 0;
        
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
            int side_normal) const
        {
            const hier::Box empty_box(d_dim);
            computeIsobaricSpecificHeatCapacity(
                data_isobaric_specific_heat_capacity,
                data_density,
                data_pressure,
                data_mass_fractions,
                side_normal,
                empty_box);
        }
        
        /*
         * Compute the isobaric specific heat capacity of mixture with isothermal and isobaric
         * equilibrium assumptions.
         */
        virtual void
        computeIsobaricSpecificHeatCapacity(
            boost::shared_ptr<pdat::SideData<double> >& data_isobaric_specific_heat_capacity,
            const boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
            const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
            int side_normal,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the Gruneisen parameter of the mixture with isothermal and isobaric equilibrium assumptions
         * (partial derivative of pressure w.r.t. specific internal energy under constant partial densities
         * divided by mixture density).
         */
        virtual double
        getGruneisenParameter(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& mass_fractions) const = 0;
        
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
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions) const
        {
            const hier::Box empty_box(d_dim);
            computeGruneisenParameter(
                data_gruneisen_parameter,
                data_density,
                data_pressure,
                data_mass_fractions,
                empty_box);
        }
        
        /*
         * Compute the Gruneisen parameter of the mixture with isothermal and isobaric equilibrium assumptions
         * (partial derivative of pressure w.r.t. specific internal energy under constant partial densities
         * divided by mixture density).
         */
        virtual void
        computeGruneisenParameter(
            boost::shared_ptr<pdat::CellData<double> >& data_gruneisen_parameter,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
            const hier::Box& domain) const = 0;
        
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
            int side_normal) const
        {
            const hier::Box empty_box(d_dim);
            computeGruneisenParameter(
                data_gruneisen_parameter,
                data_density,
                data_pressure,
                data_mass_fractions,
                side_normal,
                empty_box);
        }
        
        /*
         * Compute the Gruneisen parameter of the mixture with isothermal and isobaric equilibrium assumptions
         * (partial derivative of pressure w.r.t. specific internal energy under constant partial densities
         * divided by mixture density).
         */
        virtual void
        computeGruneisenParameter(
            boost::shared_ptr<pdat::SideData<double> >& data_gruneisen_parameter,
            const boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
            const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
            int side_normal,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the Gruneisen parameter of the mixture with isobaric equilibrium assumption
         * (partial derivative of pressure w.r.t. specific internal energy under constant partial densities
         * and volume fractions divided by mixture density).
         */
        virtual double
        getGruneisenParameter(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& mass_fractions,
            const std::vector<const double*>& volume_fractions) const = 0;
        
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
            const boost::shared_ptr<pdat::CellData<double> >& data_volume_fractions) const
        {
            const hier::Box empty_box(d_dim);
            computeGruneisenParameter(
                data_gruneisen_parameter,
                data_density,
                data_pressure,
                data_mass_fractions,
                data_volume_fractions,
                empty_box);
        }
        
        /*
         * Compute the Gruneisen parameter of the mixture with isobaric equilibrium assumption
         * (partial derivative of pressure w.r.t. specific internal energy under constant partial densities
         * and volume fractions divided by mixture density).
         */
        virtual void
        computeGruneisenParameter(
            boost::shared_ptr<pdat::CellData<double> >& data_gruneisen_parameter,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
            const boost::shared_ptr<pdat::CellData<double> >& data_volume_fractions,
            const hier::Box& domain) const = 0;
        
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
            int side_normal) const
        {
            const hier::Box empty_box(d_dim);
            computeGruneisenParameter(
                data_gruneisen_parameter,
                data_density,
                data_pressure,
                data_mass_fractions,
                data_volume_fractions,
                side_normal,
                empty_box);
        }
        
        /*
         * Compute the Gruneisen parameter of the mixture with isobaric equilibrium assumption
         * (partial derivative of pressure w.r.t. specific internal energy under constant partial densities
         * and volume fractions divided by mixture density).
         */
        virtual void
        computeGruneisenParameter(
            boost::shared_ptr<pdat::SideData<double> >& data_gruneisen_parameter,
            const boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
            const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
            const boost::shared_ptr<pdat::SideData<double> >& data_volume_fractions,
            int side_normal,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the mixture partial derivative of pressure w.r.t. partial densities under constant specific
         * internal energy with isothermal and isobaric equilibrium assumptions.
         */
        virtual std::vector<double>
        getPressureDerivativeWithPartialDensities(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& mass_fractions) const = 0;
        
        /*
         * Compute the mixture partial derivative of pressure w.r.t. partial densities under constant specific
         * internal energy with isothermal and isobaric equilibrium assumptions.
         */
        void
        computePressureDerivativeWithPartialDensities(
            boost::shared_ptr<pdat::CellData<double> >& data_partial_pressure_partial_partial_densities,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions) const
        {
            const hier::Box empty_box(d_dim);
            computePressureDerivativeWithPartialDensities(
                data_partial_pressure_partial_partial_densities,
                data_density,
                data_pressure,
                data_mass_fractions,
                empty_box);
        }
        
        /*
         * Compute the mixture partial derivative of pressure w.r.t. partial densities under constant specific
         * internal energy with isothermal and isobaric equilibrium assumptions.
         */
        virtual void
        computePressureDerivativeWithPartialDensities(
            boost::shared_ptr<pdat::CellData<double> >& data_partial_pressure_partial_partial_densities,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
            const hier::Box& domain) const = 0;
        
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
            int side_normal) const
        {
            const hier::Box empty_box(d_dim);
            computePressureDerivativeWithPartialDensities(
                data_partial_pressure_partial_partial_densities,
                data_density,
                data_pressure,
                data_mass_fractions,
                side_normal,
                empty_box);
        }
        
        /*
         * Compute the mixture partial derivative of pressure w.r.t. partial densities under constant specific
         * internal energy with isothermal and isobaric equilibrium assumptions.
         */
        virtual void
        computePressureDerivativeWithPartialDensities(
            boost::shared_ptr<pdat::SideData<double> >& data_partial_pressure_partial_partial_densities,
            const boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
            const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
            int side_normal,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the mixture partial derivative of pressure w.r.t. partial densities under constant specific
         * internal energy and volume fractions with isobaric equilibrium assumption.
         */
        virtual std::vector<double>
        getPressureDerivativeWithPartialDensities(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& mass_fractions,
            const std::vector<const double*>& volume_fractions) const = 0;
        
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
            const boost::shared_ptr<pdat::CellData<double> >& data_volume_fractions) const
        {
            const hier::Box empty_box(d_dim);
            computePressureDerivativeWithPartialDensities(
                data_partial_pressure_partial_partial_densities,
                data_density,
                data_pressure,
                data_mass_fractions,
                data_volume_fractions,
                empty_box);
        }
        
        /*
         * Compute the mixture partial derivative of pressure w.r.t. partial densities under constant specific
         * internal energy and volume fractions with isobaric equilibrium assumption.
         */
        virtual void
        computePressureDerivativeWithPartialDensities(
            boost::shared_ptr<pdat::CellData<double> >& data_partial_pressure_partial_partial_densities,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
            const boost::shared_ptr<pdat::CellData<double> >& data_volume_fractions,
            const hier::Box& domain) const = 0;
        
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
            int side_normal) const
        {
            const hier::Box empty_box(d_dim);
            computePressureDerivativeWithPartialDensities(
                data_partial_pressure_partial_partial_densities,
                data_density,
                data_pressure,
                data_mass_fractions,
                data_volume_fractions,
                side_normal,
                empty_box);
        }
        
        /*
         * Compute the mixture partial derivative of pressure w.r.t. partial densities under constant specific
         * internal energy and volume fractions with isobaric equilibrium assumption.
         */
        virtual void
        computePressureDerivativeWithPartialDensities(
            boost::shared_ptr<pdat::SideData<double> >& data_partial_pressure_partial_partial_densities,
            const boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
            const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
            const boost::shared_ptr<pdat::SideData<double> >& data_volume_fractions,
            int side_normal,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the mixture partial derivative of pressure w.r.t. volume fractions under constant specific
         * internal energy and partial densities with isobaric equilibrium assumption.
         */
        virtual std::vector<double>
        getPressureDerivativeWithVolumeFractions(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& mass_fractions,
            const std::vector<const double*>& volume_fractions) const = 0;
        
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
            const boost::shared_ptr<pdat::CellData<double> >& data_volume_fractions) const
        {
            const hier::Box empty_box(d_dim);
            computePressureDerivativeWithVolumeFractions(
                data_partial_pressure_partial_volume_fractions,
                data_density,
                data_pressure,
                data_mass_fractions,
                data_volume_fractions,
                empty_box);
        }
        
        /*
         * Compute the mixture partial derivative of pressure w.r.t. volume fractions under constant specific
         * internal energy and partial densities with isobaric equilibrium assumption.
         */
        virtual void
        computePressureDerivativeWithVolumeFractions(
            boost::shared_ptr<pdat::CellData<double> >& data_partial_pressure_partial_volume_fractions,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
            const boost::shared_ptr<pdat::CellData<double> >& data_volume_fractions,
            const hier::Box& domain) const = 0;
        
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
            int side_normal) const
        {
            const hier::Box empty_box(d_dim);
            computePressureDerivativeWithVolumeFractions(
                data_partial_pressure_partial_volume_fractions,
                data_density,
                data_pressure,
                data_mass_fractions,
                data_volume_fractions,
                side_normal,
                empty_box);
        }
        
        /*
         * Compute the mixture partial derivative of pressure w.r.t. volume fractions under constant specific
         * internal energy and partial densities with isobaric equilibrium assumption.
         */
        virtual void
        computePressureDerivativeWithVolumeFractions(
            boost::shared_ptr<pdat::SideData<double> >& data_partial_pressure_partial_volume_fractions,
            const boost::shared_ptr<pdat::SideData<double> >& data_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
            const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
            const boost::shared_ptr<pdat::SideData<double> >& data_volume_fractions,
            int side_normal,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the density of mixture with isothermal and isobaric equilibrium assumptions.
         */
        virtual double
        getMixtureDensity(
            const double* const pressure,
            const double* const temperature,
            const std::vector<const double*>& mass_fractions) const = 0;
        
        /*
         * Compute the density of mixture with isothermal and isobaric equilibrium assumptions.
         */
        void
        computeMixtureDensity(
            boost::shared_ptr<pdat::CellData<double> >& data_mixture_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions) const
        {
            const hier::Box empty_box(d_dim);
            computeMixtureDensity(
                data_mixture_density,
                data_pressure,
                data_temperature,
                data_mass_fractions,
                empty_box);
        }
        
        /*
         * Compute the density of mixture with isothermal and isobaric equilibrium assumptions.
         */
        virtual void
        computeMixtureDensity(
            boost::shared_ptr<pdat::CellData<double> >& data_mixture_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the density of mixture with isothermal and isobaric equilibrium assumptions.
         */
        void
        computeMixtureDensity(
            boost::shared_ptr<pdat::SideData<double> >& data_mixture_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
            const boost::shared_ptr<pdat::SideData<double> >& data_temperature,
            const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
            int side_normal) const
        {
            const hier::Box empty_box(d_dim);
            computeMixtureDensity(
                data_mixture_density,
                data_pressure,
                data_temperature,
                data_mass_fractions,
                side_normal,
                empty_box);
        }
        
        /*
         * Compute the density of mixture with isothermal and isobaric equilibrium assumptions.
         */
        virtual void
        computeMixtureDensity(
            boost::shared_ptr<pdat::SideData<double> >& data_mixture_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
            const boost::shared_ptr<pdat::SideData<double> >& data_temperature,
            const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
            int side_normal,
            const hier::Box& domain) const = 0;
        
        /*
         * Get the number of thermodynamic properties of a species.
         */
        virtual int
        getNumberOfSpeciesThermodynamicProperties(const int species_index = 0) const = 0;
        
        /*
         * Get the thermodynamic properties of a species.
         */
        virtual void
        getSpeciesThermodynamicProperties(
            std::vector<double*>& species_thermo_properties,
            const int species_index = 0) const = 0;
        
        /*
         * Get the molecular weight of a species.
         */
        double
        getSpeciesMolecularWeight(
            const int species_index = 0) const;
        
        /*
         * Compute the molecular weight of mixture.
         */
        double
        getMixtureMolecularWeight(
            const std::vector<const double*>& mass_fractions) const;
        
        /*
         * Compute the molecular weight of mixture.
         */
        void
        computeMixtureMolecularWeight(
            boost::shared_ptr<pdat::CellData<double> >& data_mixture_molecular_weight,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions) const
        {
            const hier::Box empty_box(d_dim);
            computeMixtureMolecularWeight(
                data_mixture_molecular_weight,
                data_mass_fractions,
                empty_box);
        }
        
        /*
         * Compute the molecular weight of mixture.
         */
        void
        computeMixtureMolecularWeight(
            boost::shared_ptr<pdat::CellData<double> >& data_mixture_molecular_weight,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
            const hier::Box& domain) const;
        
        /*
         * Compute the molecular weight of mixture.
         */
        void
        computeMixtureMolecularWeight(
            boost::shared_ptr<pdat::SideData<double> >& data_mixture_molecular_weight,
            const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
            int side_normal) const
        {
            const hier::Box empty_box(d_dim);
            computeMixtureMolecularWeight(
                data_mixture_molecular_weight,
                data_mass_fractions,
                side_normal,
                empty_box);
        }
        
        /*
         * Compute the molecular weight of mixture.
         */
        void
        computeMixtureMolecularWeight(
            boost::shared_ptr<pdat::SideData<double> >& data_mixture_molecular_weight,
            const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
            int side_normal,
            const hier::Box& domain) const;
        
        /*
         * Helper function to compute the density of mixture given the partial densities.
         */
        double
        getMixtureDensity(
            const std::vector<const double*>& partial_densities) const;
        
        /*
         * Helper function to compute the density of mixture given the partial densities.
         */
        void
        computeMixtureDensity(
            boost::shared_ptr<pdat::CellData<double> >& data_mixture_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_partial_densities) const
        {
            const hier::Box empty_box(d_dim);
            computeMixtureDensity(
                data_mixture_density,
                data_partial_densities,
                empty_box);
        }
        
        /*
         * Helper function to compute the density of mixture given the partial densities.
         */
        void
        computeMixtureDensity(
            boost::shared_ptr<pdat::CellData<double> >& data_mixture_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_partial_densities,
            const hier::Box& domain) const;
        
        /*
         * Helper function to compute the density of mixture given the partial densities.
         */
        void
        computeMixtureDensity(
            boost::shared_ptr<pdat::SideData<double> >& data_mixture_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_partial_densities,
            int side_normal) const
        {
            const hier::Box empty_box(d_dim);
            computeMixtureDensity(
                data_mixture_density,
                data_partial_densities,
                side_normal,
                empty_box);
        }
        
        /*
         * Helper function to compute the density of mixture given the partial densities.
         */
        void
        computeMixtureDensity(
            boost::shared_ptr<pdat::SideData<double> >& data_mixture_density,
            const boost::shared_ptr<pdat::SideData<double> >& data_partial_densities,
            int side_normal,
            const hier::Box& domain) const;
        
    protected:
        /*
         * Compute the molecular weight of mixture given the mass fractions.
         */
        void
        computeMixtureMolecularWeight(
            double* const M,
            const std::vector<const double*>& Y,
            const hier::IntVector& num_ghosts_mixture_molecular_weight,
            const hier::IntVector& num_ghosts_mass_fractions,
            const hier::IntVector& ghostcell_dims_mixture_molecular_weight,
            const hier::IntVector& ghostcell_dims_mass_fractions,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims) const;
        
        /*
         * Compute the density of mixture given the partial densities.
         */
        void
        computeMixtureDensity(
            double* const rho,
            const std::vector<const double*>& Z_rho,
            const hier::IntVector& num_ghosts_mixture_density,
            const hier::IntVector& num_ghosts_partial_densities,
            const hier::IntVector& ghostcell_dims_mixture_density,
            const hier::IntVector& ghostcell_dims_partial_densities,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims) const;
        
        /*
         * The object name is used for error/warning reporting.
         */
        const std::string d_object_name;
        
        /*
         * Problem dimension.
         */
        const tbox::Dimension d_dim;
        
        /*
         * Number of species.
         */
        const int d_num_species;
        
        /*
         * Type of mixing closure model.
         */
        const MIXING_CLOSURE_MODEL::TYPE d_mixing_closure_model;
        
        /*
         * boost::shared_ptr to the database of equation of state mixing rules.
         */
        const boost::shared_ptr<tbox::Database> d_equation_of_state_mixing_rules_db;
        
        /*
         * Molecular weights of different species.
         */
        std::vector<double> d_species_M;
        
};

#endif /* EQUATION_OF_STATE_MIXING_RULES_HPP */
