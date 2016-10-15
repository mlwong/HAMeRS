#ifndef EQUATION_OF_THERMAL_CONDUCTIVITY_MIXING_RULES_HPP
#define EQUATION_OF_THERMAL_CONDUCTIVITY_MIXING_RULES_HPP

#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Dimension.h"

#include <string>
#include <vector>

#include "util/mixing_rules/equations_of_thermal_conductivity/EquationOfThermalConductivity.hpp"
#include "util/mixing_rules/MixingClosureModels.hpp"

using namespace SAMRAI;

class EquationOfThermalConductivityMixingRules
{
    public:
        EquationOfThermalConductivityMixingRules(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const int& num_species,
            const MIXING_CLOSURE_MODEL& mixing_closure_model,
            const boost::shared_ptr<tbox::Database>& equation_of_thermal_conductivity_mixing_rules_db):
                d_object_name(object_name),
                d_dim(dim),
                d_num_species(num_species),
                d_mixing_closure_model(mixing_closure_model),
                d_equation_of_thermal_conductivity_mixing_rules_db(equation_of_thermal_conductivity_mixing_rules_db)
        {}
        
        /*
         * Print all characteristics of the equation of thermal conductivity mixing rules class.
         */
        virtual void
        printClassData(std::ostream& os) const = 0;
        
        /*
         * Put the characteristics of the equation of thermal conductivity mixing rules class into the restart
         * database.
         */
        virtual void
        putToRestart(
            const boost::shared_ptr<tbox::Database>& restart_db) const = 0;
        
        /*
         * Compute the thermal conductivity of the mixture with isothermal and isobaric assumptions.
         */
        virtual double
        getThermalConductivity(
            const double* const pressure,
            const double* const temperature,
            const std::vector<const double*>& mass_fraction) const = 0;
        
        /*
         * Compute the thermal conductivity of the mixture with isobaric assumption.
         */
        virtual double
        getThermalConductivity(
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
        
        /*
         * boost::shared_ptr to EquationOfThermalConductivity.
         */
        boost::shared_ptr<EquationOfThermalConductivity> d_equation_of_thermal_conductivity;
        
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
         * Mixing closure model.
         */
        const MIXING_CLOSURE_MODEL d_mixing_closure_model;
        
        /*
         * boost::shared_ptr to the database of equation of thermal conductivity mixing rules.
         */
        const boost::shared_ptr<tbox::Database> d_equation_of_thermal_conductivity_mixing_rules_db;
        
};
    

#endif /* EQUATION_OF_THERMAL_CONDUCTIVITY_MIXING_RULES_HPP */