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
        
        /*
         * Return the boost::shared_ptr to the equation of state.
         */
        const boost::shared_ptr<EquationOfState>&
        getEquationOfState() const
        {
            return d_equation_of_state;
        }
        
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
         * Compute the pressure of the mixture with isothermal and isobaric assumptions.
         */
        virtual double
        getPressure(
            const double* const density,
            const double* const internal_energy,
            const std::vector<const double*>& mass_fraction) const = 0;
        
        /*
         * Compute the pressure of the mixture with isothermal and isobaric assumptions.
         */
        void
        getPressure(
            boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fraction) const
        {
            const hier::Box empty_box(d_dim);
            getPressure(
                data_pressure,
                data_density,
                data_internal_energy,
                data_mass_fraction,
                empty_box);
        }
        
        /*
         * Compute the pressure of the mixture with isothermal and isobaric assumptions.
         */
        virtual void
        getPressure(
            boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fraction,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the pressure of the mixture with isobaric assumption.
         */
        virtual double
        getPressure(
            const double* const density,
            const double* const internal_energy,
            const std::vector<const double*>& mass_fraction,
            const std::vector<const double*>& volume_fraction) const = 0;
        
        /*
         * Compute the pressure of the mixture with isobaric assumption.
         */
        void
        getPressure(
            boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fraction,
            const boost::shared_ptr<pdat::CellData<double> >& data_volume_fraction) const
        {
            const hier::Box empty_box(d_dim);
            getPressure(
                data_pressure,
                data_density,
                data_internal_energy,
                data_mass_fraction,
                data_volume_fraction,
                empty_box);
        }
        
        /*
         * Compute the pressure of the mixture with isobaric assumption.
         */
        virtual void
        getPressure(
            boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fraction,
            const boost::shared_ptr<pdat::CellData<double> >& data_volume_fraction,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the sound speed of the mixture with isothermal and isobaric assumptions.
         */
        virtual double
        getSoundSpeed(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& mass_fraction) const = 0;
        
        /*
         * Compute the sound speed of the mixture with isothermal and isobaric assumptions.
         */
        void
        getSoundSpeed(
            boost::shared_ptr<pdat::CellData<double> >& data_sound_speed,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fraction) const
        {
            const hier::Box empty_box(d_dim);
            getSoundSpeed(
                data_sound_speed,
                data_density,
                data_pressure,
                data_mass_fraction,
                empty_box);
        }
        
        /*
         * Compute the sound speed of the mixture with isothermal and isobaric assumptions.
         */
        virtual void
        getSoundSpeed(
            boost::shared_ptr<pdat::CellData<double> >& data_sound_speed,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fraction,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the sound speed of the mixture with isobaric assumption.
         */
        virtual double
        getSoundSpeed(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& mass_fraction,
            const std::vector<const double*>& volume_fraction) const = 0;
        
        /*
         * Compute the sound speed of the mixture with isobaric assumption.
         */
        void
        getSoundSpeed(
            boost::shared_ptr<pdat::CellData<double> >& data_sound_speed,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fraction,
            const boost::shared_ptr<pdat::CellData<double> >& data_volume_fraction) const
        {
            const hier::Box empty_box(d_dim);
            getSoundSpeed(
                data_sound_speed,
                data_density,
                data_pressure,
                data_mass_fraction,
                data_volume_fraction,
                empty_box);
        }
        
        /*
         * Compute the sound speed of the mixture with isobaric assumption.
         */
        virtual void
        getSoundSpeed(
            boost::shared_ptr<pdat::CellData<double> >& data_sound_speed,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fraction,
            const boost::shared_ptr<pdat::CellData<double> >& data_volume_fraction,
            const hier::Box& domain) const= 0;
        
        /*
         * Compute the specific internal energy of the mixture with isothermal and isobaric assumptions.
         */
        virtual double
        getInternalEnergy(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& mass_fraction) const = 0;
        
        /*
         * Compute the specific internal energy of the mixture with isothermal and isobaric assumptions.
         */
        void
        getInternalEnergy(
            boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fraction) const
        {
            const hier::Box empty_box(d_dim);
            getInternalEnergy(
                data_internal_energy,
                data_density,
                data_pressure,
                data_mass_fraction,
                empty_box);
        }
        
        /*
         * Compute the specific internal energy of the mixture with isothermal and isobaric assumptions.
         */
        virtual void
        getInternalEnergy(
            boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fraction,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the specific internal energy of the mixture with isobaric assumption.
         */
        virtual double
        getInternalEnergy(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& mass_fraction,
            const std::vector<const double*>& volume_fraction) const = 0;
        
        /*
         * Compute the specific internal energy of the mixture with isobaric assumption.
         */
        void
        getInternalEnergy(
            boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fraction,
            const boost::shared_ptr<pdat::CellData<double> >& data_volume_fraction) const
        {
            const hier::Box empty_box(d_dim);
            getInternalEnergy(
                data_internal_energy,
                data_density,
                data_pressure,
                data_mass_fraction,
                data_volume_fraction,
                empty_box);
        }
        
        /*
         * Compute the specific internal energy of the mixture with isobaric assumption.
         */
        virtual void
        getInternalEnergy(
            boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fraction,
            const boost::shared_ptr<pdat::CellData<double> >& data_volume_fraction,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the temperature of the mixture with isothermal and isobaric assumptions.
         */
        virtual double
        getTemperature(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& mass_fraction) const = 0;
        
        /*
         * Compute the temperature of the mixture with isothermal and isobaric assumptions.
         */
        void
        getTemperature(
            boost::shared_ptr<pdat::CellData<double> >& data_temperature,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fraction) const
        {
            const hier::Box empty_box(d_dim);
            getTemperature(
                data_temperature,
                data_density,
                data_pressure,
                data_mass_fraction,
                empty_box);
        } 
        
        /*
         * Compute the temperature of the mixture with isothermal and isobaric assumptions.
         */
        virtual void
        getTemperature(
            boost::shared_ptr<pdat::CellData<double> >& data_temperature,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fraction,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the specific internal energy of the mixture from temperature with isothermal
         * and isobaric assumptions.
         */
        virtual double
        getInternalEnergyFromTemperature(
            const double* const density,
            const double* const temperature,
            const std::vector<const double*>& mass_fraction) const = 0;
        
        /*
         * Compute the specific internal energy of the mixture from temperature with isothermal
         * and isobaric assumptions.
         */
        void
        getInternalEnergyFromTemperature(
            boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fraction) const
        {
            const hier::Box empty_box(d_dim);
            getInternalEnergyFromTemperature(
                data_internal_energy,
                data_density,
                data_temperature,
                data_mass_fraction,
                empty_box);
        }
        
        /*
         * Compute the specific internal energy of the mixture from temperature with isothermal
         * and isobaric assumptions.
         */
        virtual void
        getInternalEnergyFromTemperature(
            boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fraction,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the isochoric specific heat capacity of mixture with isothermal assumption.
         */
        virtual double
        getIsochoricSpecificHeatCapacity(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& mass_fraction) const = 0;
        
        /*
         * Compute the isochoric specific heat capacity of mixture with isothermal assumption.
         */
        void
        getIsochoricSpecificHeatCapacity(
            boost::shared_ptr<pdat::CellData<double> >& data_isochoric_specific_heat_capacity,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fraction) const
        {
            const hier::Box empty_box(d_dim);
            getIsochoricSpecificHeatCapacity(
                data_isochoric_specific_heat_capacity,
                data_density,
                data_pressure,
                data_mass_fraction,
                empty_box);
        }
        
        /*
         * Compute the isochoric specific heat capacity of mixture with isothermal assumption.
         */
        virtual void
        getIsochoricSpecificHeatCapacity(
            boost::shared_ptr<pdat::CellData<double> >& data_isochoric_specific_heat_capacity,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fraction,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the isobaric specific heat capacity of mixture with isothermal assumption.
         */
        virtual double
        getIsobaricSpecificHeatCapacity(
            const double* const density,
            const double* const pressure,
            const std::vector<const double*>& mass_fraction) const = 0;
        
        /*
         * Compute the isobaric specific heat capacity of mixture with isothermal assumption.
         */
        void
        getIsobaricSpecificHeatCapacity(
            boost::shared_ptr<pdat::CellData<double> >& data_isobaric_specific_heat_capacity,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fraction) const
        {
            const hier::Box empty_box(d_dim);
            getIsobaricSpecificHeatCapacity(
                data_isobaric_specific_heat_capacity,
                data_density,
                data_pressure,
                data_mass_fraction,
                empty_box);
        }
        
        /*
         * Compute the isobaric specific heat capacity of mixture with isothermal assumption.
         */
        virtual void
        getIsobaricSpecificHeatCapacity(
            boost::shared_ptr<pdat::CellData<double> >& data_isobaric_specific_heat_capacity,
            const boost::shared_ptr<pdat::CellData<double> >& data_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fraction,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the density of mixture with isothermal and isobaric assumptions.
         */
        virtual double
        getMixtureDensity(
            const double* const pressure,
            const double* const temperature,
            const std::vector<const double*>& mass_fraction) const = 0;
        
        /*
         * Compute the density of mixture with isothermal and isobaric assumptions.
         */
        void
        getMixtureDensity(
            boost::shared_ptr<pdat::CellData<double> >& data_mixture_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fraction) const
        {
            const hier::Box empty_box(d_dim);
            getMixtureDensity(
                data_mixture_density,
                data_pressure,
                data_temperature,
                data_mass_fraction,
                empty_box);
        }
        
        /*
         * Compute the density of mixture with isothermal and isobaric assumptions.
         */
        virtual void
        getMixtureDensity(
            boost::shared_ptr<pdat::CellData<double> >& data_mixture_density,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fraction,
            const hier::Box& domain) const = 0;
        
        /*
         * Get the number of thermodynamic properties of a species.
         */
        virtual int
        getNumberOfSpeciesThermodynamicProperties() const = 0;
        
        /*
         * Get the thermodynamic properties of a species.
         */
        virtual void
        getSpeciesThermodynamicProperties(
            std::vector<double*>& species_thermo_properties,
            const int& species_index) const = 0;
        
        /*
         * Helper function to compute the density of mixture given the partial densities.
         */
        double
        getMixtureDensity(
            const std::vector<const double*>& partial_density) const
        {
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
            if (static_cast<int>(partial_density.size()) != d_num_species)
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "Number of partial densities provided is not"
                    << "equal to the total number of species."
                    << std::endl);
            }
#endif
            
            double mixture_density = 0.0;
            
            for (int si = 0; si < d_num_species; si++)
            {
                const double& Z_rho = *(partial_density[si]);
                
                mixture_density += Z_rho;
            }
            
            return mixture_density;            
        }
        
    protected:
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
         * boost::shared_ptr to EquationOfState.
         */
        boost::shared_ptr<EquationOfState> d_equation_of_state;
        
};

#endif /* EQUATION_OF_STATE_MIXING_RULES_HPP */
