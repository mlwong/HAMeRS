#ifndef EQUATION_OF_STATE_MIXING_RULES_IDEAL_GAS_HPP
#define EQUATION_OF_STATE_MIXING_RULES_IDEAL_GAS_HPP

#include "util/equations_of_state/EquationOfStateMixingRules.hpp"

class EquationOfStateMixingRulesIdealGas: public EquationOfStateMixingRules
{
    public:        
        EquationOfStateMixingRulesIdealGas(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const int& num_species,
            const MIXING_CLOSURE_MODEL& mixing_closure_model,
            const boost::shared_ptr<tbox::Database>& species_db);
        
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
         * Get the number of thermodynamic properties of the mixture.
         */
        int
        getNumberOfMixtureThermodynamicProperties() const;
        
        /*
         * Get the number of thermodynamic properties of a species.
         */
        int
        getNumberOfSpeciesThermodynamicProperties() const { return 4; }
        
        /*
         * Get the thermodynamic properties of the mixture.
         */
        void
        getMixtureThermodynamicProperties(
            std::vector<double*>& mixture_thermo_properties,
            const std::vector<const double*>& species_fraction) const;
        
        /*
         * Get the thermodynamic properties of a species.
         */
        void
        getSpeciesThermodynamicProperties(
            std::vector<double*>& species_thermo_properties,
            const int& species_index) const;
        
    private:
        /*
         * Compute the thermodynamic properties of the mixture with mass fraction.
         */
        void
        getMixtureThermodynamicPropertiesWithMassFraction(
            std::vector<double*>& mixture_thermo_properties,
            const std::vector<const double*>& mass_fraction) const;
        
        /*
         * Compute the thermodynamic properties of the mixture with volume fraction.
         */
        void
        getMixtureThermodynamicPropertiesWithVolumeFraction(
            std::vector<double*>& mixture_thermo_properties,
            const std::vector<const double*>& volume_fraction) const;
        
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
        
};

#endif /* EQUATION_OF_STATE_MIXING_RULES_IDEAL_GAS_HPP */
