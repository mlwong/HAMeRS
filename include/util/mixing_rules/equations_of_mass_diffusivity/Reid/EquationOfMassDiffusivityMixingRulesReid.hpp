#ifndef EQUATION_OF_MASS_DIFFUSIVITY_MIXING_RULES_REID_HPP
#define EQUATION_OF_MASS_DIFFUSIVITY_MIXING_RULES_REID_HPP

#include "util/mixing_rules/equations_of_mass_diffusivity/EquationOfMassDiffusivityMixingRules.hpp"

class EquationOfMassDiffusivityMixingRulesReid: public EquationOfMassDiffusivityMixingRules
{
    public:
        EquationOfMassDiffusivityMixingRulesReid(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const int& num_species,
            const MIXING_CLOSURE_MODEL::TYPE& mixing_closure_model,
            const HAMERS_SHARED_PTR<tbox::Database>& equation_of_mass_diffusivity_mixing_rules_db);
        
        ~EquationOfMassDiffusivityMixingRulesReid() {}
        
        /*
         * Print all characteristics of the equation of mass diffusivity mixing rules class.
         */
        void
        printClassData(std::ostream& os) const;
        
        /*
         * Put the characteristics of the equation of mass diffusivity mixing rules class into the restart
         * database.
         */
        void
        putToRestart(
            const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const;
        
        /*
         * Compute the mass diffusivities of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        void
        getMassDiffusivities(
            std::vector<Real*>& mass_diffusivities,
            const Real* const pressure,
            const Real* const temperature,
            const std::vector<const Real*>& mass_fractions) const;
        
        /*
         * Compute the mass diffusivities of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        void
        computeMassDiffusivities(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_mass_diffusivities,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_temperature,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_mass_fractions,
            const hier::Box& domain) const;
        
        /*
         * Get the number of molecular properties of a species.
         */
        int
        getNumberOfSpeciesMolecularProperties(const int species_index = 0) const
        {
            NULL_USE(species_index);
            return 3;
        }
        
        /*
         * Get the molecular properties of a species.
         */
        void
        getSpeciesMolecularProperties(
            std::vector<Real*>& species_molecular_properties,
            const int species_index = 0) const;
        
    private:
        /*
         * Compute the mass diffusivity of a binary mixture.
         */
        Real
        getMassDiffusivity(
            const Real* const pressure,
            const Real* const temperature,
            const std::vector<const Real*>& molecular_properties_1,
            const std::vector<const Real*>& molecular_properties_2) const;
        
        /*
         * Lennard-Jones energy parameter of different species.
         */
        std::vector<Real> d_species_epsilon_by_k;
        
        /*
         * Collision diameter of different species.
         */
        std::vector<Real> d_species_sigma;
        
        /*
         * Molecular weight of different species.
         */
        std::vector<Real> d_species_M;
        
};

#endif /* EQUATION_OF_MASS_DIFFUSIVITY_MIXING_RULES_REID_HPP */
