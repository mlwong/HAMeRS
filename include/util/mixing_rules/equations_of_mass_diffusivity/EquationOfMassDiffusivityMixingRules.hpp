#ifndef EQUATION_OF_MASS_DIFFUSIVITY_MIXING_RULES_HPP
#define EQUATION_OF_MASS_DIFFUSIVITY_MIXING_RULES_HPP

#include "HAMeRS_config.hpp"

#include "util/mixing_rules/MixingClosureModels.hpp"

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Dimension.h"

#include <string>
#include <vector>

using namespace SAMRAI;

class EquationOfMassDiffusivityMixingRules
{
    public:
        EquationOfMassDiffusivityMixingRules(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const int& num_species,
            const MIXING_CLOSURE_MODEL::TYPE& mixing_closure_model,
            const boost::shared_ptr<tbox::Database>& equation_of_mass_diffusivity_mixing_rules_db):
                d_object_name(object_name),
                d_dim(dim),
                d_num_species(num_species),
                d_mixing_closure_model(mixing_closure_model),
                d_equation_of_mass_diffusivity_mixing_rules_db(equation_of_mass_diffusivity_mixing_rules_db)
        {}
        
        /*
         * Print all characteristics of the equation of mass diffusivity mixing rules class.
         */
        virtual void
        printClassData(std::ostream& os) const = 0;
        
        /*
         * Put the characteristics of the equation of mass diffusivity mixing rules class into the restart
         * database.
         */
        virtual void
        putToRestart(
            const boost::shared_ptr<tbox::Database>& restart_db) const = 0;
        
        /*
         * Compute the mass diffusivities of the mixture with isothermal and isobaric equilibria assumptions.
         */
        virtual void
        getMassDiffusivities(
            std::vector<double*>& mass_diffusivities,
            const double* const pressure,
            const double* const temperature,
            const std::vector<const double*>& mass_fractions) const = 0;
        
        /*
         * Compute the mass diffusivities of the mixture with isothermal and isobaric equilibria assumptions.
         */
        virtual void
        computeMassDiffusivities(
            boost::shared_ptr<pdat::CellData<double> >& data_mass_diffusivities,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
            const hier::Box& domain) const = 0;
        
        /*
         * Get the number of molecular properties of a species.
         */
        virtual int
        getNumberOfSpeciesMolecularProperties(const int species_index = 0) const = 0;
        
        /*
         * Get the molecular properties of a species.
         */
        virtual void
        getSpeciesMolecularProperties(
            std::vector<double*>& species_molecular_properties,
            const int species_index = 0) const = 0;
        
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
         * Number species.
         */
        const int d_num_species;
        
        /*
         * Type of mixing closure model.
         */
        const MIXING_CLOSURE_MODEL::TYPE d_mixing_closure_model;
        
        /*
         * boost::shared_ptr to the database of equation of mass diffusivity mixing rules.
         */
        const boost::shared_ptr<tbox::Database> d_equation_of_mass_diffusivity_mixing_rules_db;
        
};
    

#endif /* EQUATION_OF_MASS_DIFFUSIVITY_MIXING_RULES_HPP */
