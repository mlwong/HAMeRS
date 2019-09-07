#ifndef EQUATION_OF_THERMAL_CONDUCTIVITY_MIXING_RULES_HPP
#define EQUATION_OF_THERMAL_CONDUCTIVITY_MIXING_RULES_HPP

#include "util/mixing_rules/equations_of_thermal_conductivity/EquationOfThermalConductivity.hpp"
#include "util/mixing_rules/MixingClosureModels.hpp"

#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Dimension.h"

#include <string>
#include <vector>

using namespace SAMRAI;

class EquationOfThermalConductivityMixingRules
{
    public:
        EquationOfThermalConductivityMixingRules(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const int& num_species,
            const MIXING_CLOSURE_MODEL::TYPE& mixing_closure_model,
            const boost::shared_ptr<tbox::Database>& equation_of_thermal_conductivity_mixing_rules_db):
                d_object_name(object_name),
                d_dim(dim),
                d_num_species(num_species),
                d_mixing_closure_model(mixing_closure_model),
                d_equation_of_thermal_conductivity_mixing_rules_db(equation_of_thermal_conductivity_mixing_rules_db)
        {}
        
        virtual ~EquationOfThermalConductivityMixingRules() {}
        
        /*
         * Return the boost::shared_ptr to the equation of thermal conductivity.
         */
        virtual const boost::shared_ptr<EquationOfThermalConductivity>&
        getEquationOfThermalConductivity(const int species_index = 0) const = 0;
        
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
         * Compute the thermal conductivity of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        virtual double
        getThermalConductivity(
            const double* const pressure,
            const double* const temperature,
            const std::vector<const double*>& mass_fractions) const = 0;
        
        /*
         * Compute the thermal conductivity of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        void
        computeThermalConductivity(
            boost::shared_ptr<pdat::CellData<double> >& data_thermal_conductivity,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
            const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions) const
        {
            const hier::Box empty_box(d_dim);
            computeThermalConductivity(
                data_thermal_conductivity,
                data_pressure,
                data_temperature,
                data_mass_fractions,
                empty_box);
        }
        
        /*
         * Compute the thermal conductivity of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        virtual void
        computeThermalConductivity(
            boost::shared_ptr<pdat::CellData<double> >& data_thermal_conductivity,
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
         * boost::shared_ptr to the database of equation of thermal conductivity mixing rules.
         */
        const boost::shared_ptr<tbox::Database> d_equation_of_thermal_conductivity_mixing_rules_db;
        
};
    

#endif /* EQUATION_OF_THERMAL_CONDUCTIVITY_MIXING_RULES_HPP */
