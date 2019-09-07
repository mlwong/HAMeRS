#ifndef EQUATION_OF_MASS_DIFFUSIVITY_MIXING_RULES_CONSTANT_HPP
#define EQUATION_OF_MASS_DIFFUSIVITY_MIXING_RULES_CONSTANT_HPP

#include "util/mixing_rules/equations_of_mass_diffusivity/EquationOfMassDiffusivityMixingRules.hpp"

class EquationOfMassDiffusivityMixingRulesConstant: public EquationOfMassDiffusivityMixingRules
{
    public:
        EquationOfMassDiffusivityMixingRulesConstant(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const int& num_species,
            const MIXING_CLOSURE_MODEL::TYPE& mixing_closure_model,
            const boost::shared_ptr<tbox::Database>& equation_of_mass_diffusivity_mixing_rules_db);
        
        ~EquationOfMassDiffusivityMixingRulesConstant() {}
        
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
            const boost::shared_ptr<tbox::Database>& restart_db) const;
        
        /*
         * Compute the mass diffusivities of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        void
        getMassDiffusivities(
            std::vector<double*>& mass_diffusivities,
            const double* const pressure,
            const double* const temperature,
            const std::vector<const double*>& mass_fractions) const;
        
        /*
         * Compute the mass diffusivities of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        void
        computeMassDiffusivities(
            boost::shared_ptr<pdat::CellData<double> >& data_mass_diffusivities,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
            const hier::Box& domain) const;
        
        /*
         * Get the number of molecular properties of a species.
         */
        int
        getNumberOfSpeciesMolecularProperties(const int species_index = 0) const
        {
            NULL_USE(species_index);
            return 1;
        }
        
        /*
         * Get the molecular properties of a species.
         */
        void
        getSpeciesMolecularProperties(
            std::vector<double*>& species_molecular_properties,
            const int species_index = 0) const;
        
    private:
        /*
         * Mass diffusivity of different species.
         */
        std::vector<double> d_species_D;
        
};

#endif /* EQUATION_OF_MASS_DIFFUSIVITY_MIXING_RULES_CONSTANT_HPP */
