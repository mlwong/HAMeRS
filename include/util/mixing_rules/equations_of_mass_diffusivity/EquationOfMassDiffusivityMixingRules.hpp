#ifndef EQUATION_OF_MASS_DIFFUSIVITY_MIXING_RULES_HPP
#define EQUATION_OF_MASS_DIFFUSIVITY_MIXING_RULES_HPP

#include "HAMeRS_config.hpp"

#include "HAMeRS_memory.hpp"

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
            const HAMERS_SHARED_PTR<tbox::Database>& equation_of_mass_diffusivity_mixing_rules_db):
                d_object_name(object_name),
                d_dim(dim),
                d_num_species(num_species),
                d_mixing_closure_model(mixing_closure_model),
                d_equation_of_mass_diffusivity_mixing_rules_db(equation_of_mass_diffusivity_mixing_rules_db)
        {}
        
        virtual ~EquationOfMassDiffusivityMixingRules() {}
        
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
            const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const = 0;
        
        /*
         * Compute the mass diffusivities of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        virtual void
        getMassDiffusivities(
            std::vector<Real*>& mass_diffusivities,
            const Real* const pressure,
            const Real* const temperature,
            const std::vector<const Real*>& mass_fractions) const = 0;
        
        /*
         * Compute the mass diffusivities of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        void
        computeMassDiffusivities(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_mass_diffusivities,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_temperature,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_mass_fractions) const
        {
            const hier::Box empty_box(d_dim);
            computeMassDiffusivities(
                data_mass_diffusivities,
                data_pressure,
                data_temperature,
                data_mass_fractions,
                empty_box);
        }
        
        /*
         * Compute the mass diffusivities of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        virtual void
        computeMassDiffusivities(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_mass_diffusivities,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_temperature,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_mass_fractions,
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
            std::vector<Real*>& species_molecular_properties,
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
         * HAMERS_SHARED_PTR to the database of equation of mass diffusivity mixing rules.
         */
        const HAMERS_SHARED_PTR<tbox::Database> d_equation_of_mass_diffusivity_mixing_rules_db;
        
};
    

#endif /* EQUATION_OF_MASS_DIFFUSIVITY_MIXING_RULES_HPP */
