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
            const HAMERS_SHARED_PTR<tbox::Database>& equation_of_bulk_viscosity_mixing_rules_db);
        
        ~EquationOfBulkViscosityMixingRulesConstant() {}
        
        /*
         * Return the HAMERS_SHARED_PTR to the equation of bulk viscosity.
         */
        const HAMERS_SHARED_PTR<EquationOfBulkViscosity>&
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
            const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const;
        
        /*
         * Compute the bulk viscosity of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        Real
        getBulkViscosity(
            const Real* const pressure,
            const Real* const temperature,
            const std::vector<const Real*>& mass_fractions) const;
        
        /*
         * Compute the bulk viscosity of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        void
        computeBulkViscosity(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_bulk_viscosity,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_temperature,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_mass_fractions,
            const hier::Box& domain) const;
        
        /*
         * Compute the bulk viscosity of the mixture with isobaric equilibrium assumption.
         */
        Real
        getBulkViscosity(
            const Real* const pressure,
            const std::vector<const Real*>& species_temperatures,
            const std::vector<const Real*>& mass_fractions,
            const std::vector<const Real*>& volume_fractions) const;
        
        /*
         * Compute the bulk viscosity of the mixture with isobaric equilibrium assumption.
         */
        void
        computeBulkViscosity(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_bulk_viscosity,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& data_species_temperatures,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_mass_fractions,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_volume_fractions,
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
            std::vector<Real*>& species_molecular_properties,
            const int species_index = 0) const;
        
    private:
        /*
         * Bulk viscosities of different species.
         */
        std::vector<Real> d_species_mu_v;
        
        /*
         * Molecular weights of different species.
         */
        std::vector<Real> d_species_M;
        
        /*
         * HAMERS_SHARED_PTR to EquationOfBulkViscosity.
         */
        HAMERS_SHARED_PTR<EquationOfBulkViscosity> d_equation_of_bulk_viscosity;
        
};

#endif /* EQUATION_OF_BULK_VISCOSITY_MIXING_RULES_CONSTANT_HPP */
