#ifndef EQUATION_OF_SHEAR_VISCOSITY_MIXING_RULES_HPP
#define EQUATION_OF_SHEAR_VISCOSITY_MIXING_RULES_HPP

#include "util/mixing_rules/equations_of_shear_viscosity/EquationOfShearViscosity.hpp"
#include "util/mixing_rules/MixingClosureModels.hpp"

#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Dimension.h"

#include <string>
#include <vector>

using namespace SAMRAI;

class EquationOfShearViscosityMixingRules
{
    public:
        EquationOfShearViscosityMixingRules(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const int& num_species,
            const MIXING_CLOSURE_MODEL::TYPE& mixing_closure_model,
            const boost::shared_ptr<tbox::Database>& equation_of_shear_viscosity_mixing_rules_db):
                d_object_name(object_name),
                d_dim(dim),
                d_num_species(num_species),
                d_mixing_closure_model(mixing_closure_model),
                d_equation_of_shear_viscosity_mixing_rules_db(equation_of_shear_viscosity_mixing_rules_db)
        {}
        
        /*
         * Return the boost::shared_ptr to the equation of shear viscosity.
         */
        const boost::shared_ptr<EquationOfShearViscosity>&
        getEquationOfShearViscosity() const
        {
            return d_equation_of_shear_viscosity;
        }
        
        /*
         * Print all characteristics of the equation of shear viscosity mixing rules class.
         */
        virtual void
        printClassData(std::ostream& os) const = 0;
        
        /*
         * Put the characteristics of the equation of shear viscosity mixing rules class into the restart
         * database.
         */
        virtual void
        putToRestart(
            const boost::shared_ptr<tbox::Database>& restart_db) const = 0;
        
        /*
         * Compute the shear viscosity of the mixture with isothermal and isobaric assumptions.
         */
        virtual double
        getShearViscosity(
            const double* const pressure,
            const double* const temperature,
            const std::vector<const double*>& mass_fraction) const = 0;
        
        /*
         * Compute the shear viscosity of the mixture with isobaric assumption.
         */
        virtual double
        getShearViscosity(
            const double* const pressure,
            const std::vector<const double*>& temperature,
            const std::vector<const double*>& mass_fraction,
            const std::vector<const double*>& volume_fraction) const = 0;
        
        /*
         * Get the number of molecular properties of a species.
         */
        virtual int
        getNumberOfSpeciesMolecularProperties() const = 0;
        
        /*
         * Get the molecular properties of a species.
         */
        virtual void
        getSpeciesMolecularProperties(
            std::vector<double*>& species_molecular_properties,
            const int& species_index) const = 0;
        
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
         * boost::shared_ptr to the database of equation of shear viscosity mixing rules.
         */
        const boost::shared_ptr<tbox::Database> d_equation_of_shear_viscosity_mixing_rules_db;
        
        /*
         * boost::shared_ptr to EquationOfShearViscosity.
         */
        boost::shared_ptr<EquationOfShearViscosity> d_equation_of_shear_viscosity;
        
};
    

#endif /* EQUATION_OF_SHEAR_VISCOSITY_MIXING_RULES_HPP */
