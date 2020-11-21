#ifndef EQUATION_OF_BULK_VISCOSITY_MIXING_RULES_CRAMER_HPP
#define EQUATION_OF_BULK_VISCOSITY_MIXING_RULES_CRAMER_HPP

#include "util/mixing_rules/equations_of_bulk_viscosity/EquationOfBulkViscosityMixingRules.hpp"
#include "util/mixing_rules/equations_of_bulk_viscosity/Cramer/EquationOfBulkViscosityCramer.hpp"

class EquationOfBulkViscosityMixingRulesCramer: public EquationOfBulkViscosityMixingRules
{
    public:
        EquationOfBulkViscosityMixingRulesCramer(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const int& num_species,
            const MIXING_CLOSURE_MODEL::TYPE& mixing_closure_model,
            const HAMERS_SHARED_PTR<tbox::Database>& equation_of_bulk_viscosity_mixing_rules_db);
        
        ~EquationOfBulkViscosityMixingRulesCramer() {}
        
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
        double
        getBulkViscosity(
            const double* const pressure,
            const double* const temperature,
            const std::vector<const double*>& mass_fractions) const;
        
        /*
         * Compute the bulk viscosity of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        void
        computeBulkViscosity(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& data_bulk_viscosity,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_temperature,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_mass_fractions,
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
            HAMERS_SHARED_PTR<pdat::CellData<double> >& data_bulk_viscosity,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& data_species_temperatures,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_mass_fractions,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_volume_fractions,
            const hier::Box& domain) const;
        
        /*
         * Get the number of molecular properties of a species.
         */
        int
        getNumberOfSpeciesMolecularProperties(const int species_index = 0) const
        {
            NULL_USE(species_index);
            return 8;
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
         * Ratio of specific heats of different species.
         */
        std::vector<double> d_species_gamma;
        
        /*
         * Parameters for rotational mode.
         */
        std::vector<double> d_species_A_r;
        std::vector<double> d_species_B_r;
        
        /*
         * Parameters for vibrational mode.
         */
        std::vector<double> d_species_c_v_v;
        std::vector<double> d_species_A_v;
        std::vector<double> d_species_B_v;
        std::vector<double> d_species_C_v;
        
        /*
         * Molecular weights of different species.
         */
        std::vector<double> d_species_M;
        
        /*
         * HAMERS_SHARED_PTR to EquationOfBulkViscosity.
         */
        HAMERS_SHARED_PTR<EquationOfBulkViscosity> d_equation_of_bulk_viscosity;
        
};

#endif /* EQUATION_OF_BULK_VISCOSITY_MIXING_RULES_CRAMER_HPP */
