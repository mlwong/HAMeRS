#ifndef EQUATION_OF_SHEAR_VISCOSITY_MIXING_RULES_CHAPMAN_ENSKOG_HPP
#define EQUATION_OF_SHEAR_VISCOSITY_MIXING_RULES_CHAPMAN_ENSKOG_HPP

#include "util/mixing_rules/equations_of_shear_viscosity/EquationOfShearViscosityMixingRules.hpp"

#include "util/mixing_rules/equations_of_shear_viscosity/Chapman-Enskog/EquationOfShearViscosityChapmanEnskog.hpp"

class EquationOfShearViscosityMixingRulesChapmanEnskog: public EquationOfShearViscosityMixingRules
{
    public:
        EquationOfShearViscosityMixingRulesChapmanEnskog(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const int& num_species,
            const MIXING_CLOSURE_MODEL::TYPE& mixing_closure_model,
            const boost::shared_ptr<tbox::Database>& equation_of_shear_viscosity_mixing_rules_db);
        
        ~EquationOfShearViscosityMixingRulesChapmanEnskog() {}
        
        /*
         * Return the boost::shared_ptr to the equation of shear viscosity.
         */
        const boost::shared_ptr<EquationOfShearViscosity>&
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
            const boost::shared_ptr<tbox::Database>& restart_db) const;
        
        /*
         * Compute the shear viscosity of the mixture with isothermal and isobaric equilibria assumptions.
         */
        double
        getShearViscosity(
            const double* const pressure,
            const double* const temperature,
            const std::vector<const double*>& mass_fractions) const;
        
        /*
         * Compute the shear viscosity of the mixture with isothermal and isobaric equilibria assumptions.
         */
        void
        computeShearViscosity(
            boost::shared_ptr<pdat::CellData<double> >& data_shear_viscosity,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
            const hier::Box& domain) const;
        
        /*
         * Compute the shear viscosity of the mixture with isobaric equilibrium assumption.
         */
        double
        getShearViscosity(
            const double* const pressure,
            const std::vector<const double*>& species_temperature,
            const std::vector<const double*>& mass_fractions,
            const std::vector<const double*>& volume_fractions) const;
        
        /*
         * Compute the shear viscosity of the mixture with isobaric equilibrium assumption.
         */
        void
        computeShearViscosity(
            boost::shared_ptr<pdat::CellData<double> >& data_shear_viscosity,
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
            return 3;
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
         * Lennard-Jones energy parameter of different species.
         */
        std::vector<double> d_species_epsilon_by_k;
        
        /*
         * Collision diameter of different species.
         */
        std::vector<double> d_species_sigma;
        
        /*
         * Molecular weight of different species.
         */
        std::vector<double> d_species_M;
        
        /*
         * boost::shared_ptr to EquationOfShearViscosity.
         */
        boost::shared_ptr<EquationOfShearViscosity> d_equation_of_shear_viscosity;
        
};
    

#endif /* EQUATION_OF_SHEAR_VISCOSITY_MIXING_RULES_CHAPMAN_ENSKOG_HPP */
