#ifndef EQUATION_OF_SHEAR_VISCOSITY_MIXING_RULES_CONSTANT_HPP
#define EQUATION_OF_SHEAR_VISCOSITY_MIXING_RULES_CONSTANT_HPP

#include "util/mixing_rules/equations_of_shear_viscosity/EquationOfShearViscosityMixingRules.hpp"

#include "util/mixing_rules/equations_of_shear_viscosity/constant/EquationOfShearViscosityConstant.hpp"

class EquationOfShearViscosityMixingRulesConstant: public EquationOfShearViscosityMixingRules
{
    public:
        EquationOfShearViscosityMixingRulesConstant(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const int& num_species,
            const MIXING_CLOSURE_MODEL::TYPE& mixing_closure_model,
            const HAMERS_SHARED_PTR<tbox::Database>& equation_of_shear_viscosity_mixing_rules_db);
        
        ~EquationOfShearViscosityMixingRulesConstant() {}
        
        /*
         * Return the HAMERS_SHARED_PTR to the equation of shear viscosity.
         */
        const HAMERS_SHARED_PTR<EquationOfShearViscosity>&
        getEquationOfShearViscosity(const int species_index = 0) const
        {
            NULL_USE(species_index);
            return d_equation_of_shear_viscosity;
        }
        
        /*
         * Print all characteristics of the equation of shear viscosity mixing rules class.
         */
        void
        printClassData(std::ostream& os) const;
        
        /*
         * Put the characteristics of the equation of shear viscosity mixing rules class into the restart
         * database.
         */
        void
        putToRestart(
            const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const;
        
        /*
         * Compute the shear viscosity of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        double
        getShearViscosity(
            const double* const pressure,
            const double* const temperature,
            const std::vector<const double*>& mass_fractions) const;
        
        /*
         * Compute the shear viscosity of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        void
        computeShearViscosity(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& data_shear_viscosity,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_temperature,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_mass_fractions,
            const hier::Box& domain) const;
        
        /*
         * Compute the shear viscosity of the mixture with isobaric equilibrium assumption.
         */
        double
        getShearViscosity(
            const double* const pressure,
            const std::vector<const double*>& species_temperatures,
            const std::vector<const double*>& mass_fractions,
            const std::vector<const double*>& volume_fractions) const;
        
        /*
         * Compute the shear viscosity of the mixture with isobaric equilibrium assumption.
         */
        void
        computeShearViscosity(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& data_shear_viscosity,
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
         * Shear viscosities of different species.
         */
        std::vector<double> d_species_mu;
        
        /*
         * Molecular weights of different species.
         */
        std::vector<double> d_species_M;
        
        /*
         * HAMERS_SHARED_PTR to EquationOfShearViscosity.
         */
        HAMERS_SHARED_PTR<EquationOfShearViscosity> d_equation_of_shear_viscosity;
        
};

#endif /* EQUATION_OF_SHEAR_VISCOSITY_MIXING_RULES_CONSTANT_HPP */
