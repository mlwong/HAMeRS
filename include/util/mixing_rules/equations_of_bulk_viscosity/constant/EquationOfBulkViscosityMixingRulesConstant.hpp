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
        
        ~EquationOfBulkViscosityMixingRulesConstant() {}
        
        /*
         * Return the boost::shared_ptr to the equation of bulk viscosity.
         */
        const boost::shared_ptr<EquationOfBulkViscosity>&
        getEquationOfBulkViscosity(const int species_index = 0) const
        {
            NULL_USE(species_index);
            return d_equation_of_bulk_viscosity;
        }
        
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
         * Compute the bulk viscosity of the mixture with isothermal and isobaric equilibria assumptions.
         */
        double
        getBulkViscosity(
            const double* const pressure,
            const double* const temperature,
            const std::vector<const double*>& mass_fractions) const;
        
        /*
         * Compute the bulk viscosity of the mixture with isothermal and isobaric equilibria assumptions.
         */
        void
        computeBulkViscosity(
            boost::shared_ptr<pdat::CellData<double> >& data_bulk_viscosity,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
            const hier::Box& domain) const;
        
        /*
         * Compute the bulk viscosity of the mixture with isobaric equilibrium assumption.
         */
        double
        getBulkViscosity(
            const double* const pressure,
            const std::vector<const double*>& species_temperatures,
            const std::vector<const double*>& mass_fractions,
            const std::vector<const double*>& volume_fractions) const;
        
        /*
         * Compute the bulk viscosity of the mixture with isobaric equilibrium assumption.
         */
        void
        computeBulkViscosity(
            boost::shared_ptr<pdat::CellData<double> >& data_bulk_viscosity,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_species_temperatures,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
            const boost::shared_ptr<pdat::CellData<double> >& data_volume_fractions,
            const hier::Box& domain) const;
        
        /*
         * Get the number of molecular properties of a species.
         */
        int
        getNumberOfSpeciesMolecularProperties(const int species_index = 0) const
        {
            NULL_USE(species_index);
            return 2;
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
         * Bulk viscosities of different species.
         */
        std::vector<double> d_species_mu_v;
        
        /*
         * Molecular weights of different species.
         */
        std::vector<double> d_species_M;
        
        /*
         * boost::shared_ptr to EquationOfBulkViscosity.
         */
        boost::shared_ptr<EquationOfBulkViscosity> d_equation_of_bulk_viscosity;
        
};

#endif /* EQUATION_OF_BULK_VISCOSITY_MIXING_RULES_CONSTANT_HPP */
