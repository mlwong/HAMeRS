#ifndef EQUATION_OF_STATE_HPP
#define EQUATION_OF_STATE_HPP

#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/pdat/CellData.h"

#include "boost/shared_ptr.hpp"
#include <string>
#include <vector>

#include "equation_of_state/ThermalProcessAssumptions.hpp"

using namespace SAMRAI;

enum EQUATIONS_OF_STATE { IDEAL_GAS };

class EquationOfState
{
    public:
        EquationOfState(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const int& num_species,
            const boost::shared_ptr<tbox::Database>& equation_of_state_db,
            const THERMAL_PROCESS_ASSUMPTION& thermal_pro_assum):
                d_object_name(object_name),
                d_dim(dim),
                d_num_species(num_species),
                d_equation_of_state_db(equation_of_state_db),
                d_thermal_pro_assum(thermal_pro_assum)
        {}
        
        /*
         * Get the label of the equation of state.
         */
        virtual EQUATIONS_OF_STATE getLabel() = 0;
        
        /*
         * Helper function to compute the total density of the mixture given
         * the partial densities.
         */
        double
        getTotalDensity(
            const std::vector<const double*>& partial_density)
        {
            double mixture_density = 0.0;
            
            if (static_cast<int>(partial_density.size()) == d_num_species)
            {
                for (int si = 0; si < d_num_species; si++)
                {
                    const double& Z_rho = *(partial_density[si]);
                    
                    mixture_density += Z_rho;
                }
            }
            else
            {
                TBOX_ERROR(d_object_name
                           << ": "
                           << "Number of partial densities provided is not"
                           << "equal to the total number of species."
                           << std::endl);
            }
            
            return mixture_density;            
        }
        
        /*
         * Print all characteristics of the equation of state class.
         */
        virtual void
        printClassData(std::ostream& os) const = 0;
        
        /*
         * Put the characteristics of the equation of state into the restart
         * database.
         */
        virtual void
        putToRestart(
            const boost::shared_ptr<tbox::Database>& restart_db) const = 0;
        
        /*
         * Compute the pressure for single-species flow.
         */
        virtual double
        getPressure(
            const double* const density,
            const std::vector<const double*>& momentum,
            const double* const total_energy) = 0;
        
        /*
         * Compute the pressure for multi-species flow with total density and
         * mass fractions.
         */
        virtual double
        getPressureWithMassFraction(
            const double* const density,
            const std::vector<const double*>& momentum,
            const double* const total_energy,
            const std::vector<const double*>& mass_fraction) = 0;
        
        /*
         * Compute the pressure for multi-species flow with partial densities and
         * mass fractions.
         */
        virtual double
        getPressureWithMassFraction(
            const std::vector<const double*>& partial_density,
            const std::vector<const double*>& momentum,
            const double* const total_energy,
            const std::vector<const double*>& mass_fraction) = 0;
        
        /*
         * Compute the pressure for multi-species flow with total density and
         * volume fractions.
         */
        virtual double
        getPressureWithVolumeFraction(
            const double* const density,
            const std::vector<const double*>& momentum,
            const double* const total_energy,
            const std::vector<const double*>& volume_fraction) = 0;
        
        /*
         * Compute the pressure for multi-species flow with partial densities and
         * volume fractions.
         */
        virtual double
        getPressureWithVolumeFraction(
            const std::vector<const double*>& partial_density,
            const std::vector<const double*>& momentum,
            const double* const total_energy,
            const std::vector<const double*>& volume_fraction) = 0;
        
        /*
         * Compute the sound speed for single-species flow.
         */
        virtual double
        getSoundSpeed(
            const double* const density,
            const std::vector<const double*>& momentum,
            const double* const total_energy) = 0;
        
        /*
         * Compute the sound speed for single-species flow with pressure.
         * 
         * Recommended to use with getPressure() if both pressure and sound speed
         * data is required.
         */
        virtual double
        getSoundSpeedWithPressure(
            const double* const density,
            const double* const pressure) = 0;
        
        /*
         * Compute the sound speed for multi-species flow with total density and
         * mass fractions.
         */
        virtual double
        getSoundSpeedWithMassFraction(
            const double* const density,
            const std::vector<const double*>& momentum,
            const double* const total_energy,
            const std::vector<const double*>& mass_fraction) = 0;
        
        /*
         * Compute the sound speed for multi-species flow with partial densities and
         * mass fractions.
         */
        virtual double
        getSoundSpeedWithMassFraction(
            const std::vector<const double*>& partial_density,
            const std::vector<const double*>& momentum,
            const double* const total_energy,
            const std::vector<const double*>& mass_fraction) = 0;
        
        /*
         * Compute the sound speed for multi-species flow with total density,
         * mass fractions and pressure.
         * 
         * Recommended to use with getPressureWithMassFraction() if both pressure
         * and sound speed data is required.
         */
        virtual double
        getSoundSpeedWithMassFractionAndPressure(
            const double* const density,
            const std::vector<const double*>& mass_fraction,
            const double* const pressure) = 0;
        
        /*
         * Compute the sound speed for multi-species flow with partial densities,
         * mass fractions and pressure.
         * 
         * Recommended to use with getPressureWithMassFraction() if both pressure
         * and sound speed data is required.
         */
        virtual double
        getSoundSpeedWithMassFractionAndPressure(
            const std::vector<const double*>& partial_density,
            const std::vector<const double*>& mass_fraction,
            const double* const pressure) = 0;
        
        /*
         * Compute the sound speed for multi-species flow with total density and
         * volume fractions.
         */
        virtual double
        getSoundSpeedWithVolumeFraction(
            const double* const density,
            const std::vector<const double*>& momentum,
            const double* const total_energy,
            const std::vector<const double*>& volume_fraction) = 0;
        
        /*
         * Compute the sound speed for multi-species flow with partial densities and
         * volume fractions.
         */
        virtual double
        getSoundSpeedWithVolumeFraction(
            const std::vector<const double*>& partial_density,
            const std::vector<const double*>& momentum,
            const double* const total_energy,
            const std::vector<const double*>& volume_fraction) = 0;
        
        /*
         * Compute the sound speed for multi-species flow with total density,
         * volume fractions and pressure.
         * 
         * Recommended to use with getPressureWithVolumeFraction() if both pressure
         * and sound speed data is required.
         */
        virtual double
        getSoundSpeedWithVolumeFractionAndPressure(
            const double* const density,
            const std::vector<const double*>& volume_fraction,
            const double* const pressure) = 0;
        
        /*
         * Compute the sound speed for multi-species flow with partial densities,
         * volume fractions and pressure.
         * 
         * Recommended to use with getPressureWithVolumeFraction() if both pressure
         * and sound speed data is required.
         */
        virtual double
        getSoundSpeedWithVolumeFractionAndPressure(
            const std::vector<const double*>& partial_density,
            const std::vector<const double*>& volume_fraction,
            const double* const pressure) = 0;
        
        /*
         * Compute the total energy for single-species flow.
         */
        virtual double
        getTotalEnergy(
            const double* const density,
            const std::vector<const double*>& velocity,
            const double* const pressure) = 0;
        
        /*
         * Compute the total energy for multi-species flow with total density and
         * mass fractions.
         */
        virtual double
        getTotalEnergyWithMassFraction(
            const double* const density,
            const std::vector<const double*>& velocity,
            const double* const pressure,
            const std::vector<const double*>& mass_fraction) = 0;
        
        /*
         * Compute the total energy for multi-species flow with partial densities and
         * mass fractions.
         */
        virtual double
        getTotalEnergyWithMassFraction(
           const std::vector<const double*>& partial_density,
           const std::vector<const double*>& velocity,
           const double* const pressure,
           const std::vector<const double*>& mass_fraction) = 0;
        
        /*
         * Compute the total energy for multi-species flow with total density and
         * volume fractions.
         */
        virtual double
        getTotalEnergyWithVolumeFraction(
           const double* const density,
           const std::vector<const double*>& velocity,
           const double* const pressure,
           const std::vector<const double*>& volume_fraction) = 0;
        
        /*
         * Compute the total energy for multi-species flow with partial densities and
         * volume fractions.
         */
        virtual double
        getTotalEnergyWithVolumeFraction(
           const std::vector<const double*>& partial_density,
           const std::vector<const double*>& velocity,
           const double* const pressure,
           const std::vector<const double*>& volume_fraction) = 0;
        
        /*
         * Compute a thermodynamic property of a particular species
         * with mass fractions.
         */
        virtual double
        getMixtureThermodynamicPropertyWithMassFraction(
            const std::string& property_name,
            const std::vector<const double*>& mass_fraction) = 0;
        
        /*
         * Compute a thermodynamic property of a particular species
         * with volume fractions.
         */
        virtual double
        getMixtureThermodynamicPropertyWithVolumeFraction(
            const std::string& property_name,
            const std::vector<const double*>& volume_fraction) = 0;
        
        /*
         * Get a thermodynamic property of a particular species.
         */
        virtual double&
        getSpeciesThermodynamicProperty(
            const std::string& property_name,
            const int species_index = 0) = 0;
        
        /*
         * Determine whether a thermodynamic property is registered.
         */
        virtual bool
        hasThermodynamicProperty(
            const std::string& property_name) = 0;
        
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
         * boost::shared_ptr to the grid geometry.
         */
        const int d_num_species;
        
        /*
         * boost::shared_ptr to the database of equation of state.
         */
        const boost::shared_ptr<tbox::Database> d_equation_of_state_db;
        
        /*
         * Assumption of thermal process:
         * ISOTHERMAL, ISOBARIC
         */
        const THERMAL_PROCESS_ASSUMPTION d_thermal_pro_assum;
};

#endif /* EQUATION_OF_STATE_HPP */
