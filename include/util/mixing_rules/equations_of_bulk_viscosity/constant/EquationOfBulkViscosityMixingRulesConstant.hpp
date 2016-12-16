#ifndef EQUATION_OF_BULK_VISCOSITY_MIXING_RULES_CONSTANT_HPP
#define EQUATION_OF_BULK_VISCOSITY_MIXING_RULES_CONSTANT_HPP

#include "util/mixing_rules/equations_of_bulk_viscosity/EquationOfBulkViscosityMixingRules.hpp"
#include "util/mixing_rules/equations_of_bulk_viscosity/constant/EquationOfBulkViscosityConstant.hpp"

class EquationOfBulkViscosityMixingRulesConstant: public EquationOfBulkViscosityMixingRules
{
    public:
        EquationOfBulkViscosityMixingRulesConstant(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const int& num_species,
            const MIXING_CLOSURE_MODEL::TYPE& mixing_closure_model,
            const boost::shared_ptr<tbox::Database>& equation_of_bulk_viscosity_mixing_rules_db);
        
        /*
         * Print all characteristics of the equation of bulk viscosity mixing rules class.
         */
        void
        printClassData(std::ostream& os) const;
        
        /*
         * Put the characteristics of the equation of bulk viscosity mixing rules class into the restart
         * database.
         */
        void
        putToRestart(
            const boost::shared_ptr<tbox::Database>& restart_db) const;
        
        /*
         * Compute the bulk viscosity of the mixture with isothermal and isobaric assumptions.
         */
        double
        getBulkViscosity(
            const double* const pressure,
            const double* const temperature,
            const std::vector<const double*>& mass_fraction) const;
        
        /*
         * Compute the bulk viscosity of the mixture with isobaric assumption.
         */
        double
        getBulkViscosity(
            const double* const pressure,
            const std::vector<const double*>& temperature,
            const std::vector<const double*>& mass_fraction,
            const std::vector<const double*>& volume_fraction) const;
        
        /*
         * Get the number of molecular properties of a species.
         */
        int
        getNumberOfSpeciesMolecularProperties() const { return 2; }
        
        /*
         * Get the molecular properties of a species.
         */
        void
        getSpeciesMolecularProperties(
            std::vector<double*>& species_molecular_properties,
            const int& species_index) const;
        
    private:
        /*
         * Bulk viscosities of different species.
         */
        std::vector<double> d_species_mu_v;
        
        /*
         * Molecular weights of different species.
         */
        std::vector<double> d_species_M;
        
};

#endif /* EQUATION_OF_BULK_VISCOSITY_MIXING_RULES_CONSTANT_HPP */