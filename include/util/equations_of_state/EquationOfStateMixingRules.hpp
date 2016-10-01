#ifndef EQUATION_OF_STATE_MIXING_RULES_HPP
#define EQUATION_OF_STATE_MIXING_RULES_HPP

#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Dimension.h"

#include <string>
#include <vector>

#include "util/equations_of_state/MixingClosureModels.hpp"

using namespace SAMRAI;

class EquationOfStateMixingRules
{
    public:
        EquationOfStateMixingRules(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const int& num_species,
            const MIXING_CLOSURE_MODEL& mixing_closure_model,
            const boost::shared_ptr<tbox::Database>& equation_of_state_mixing_rules_db):
                d_object_name(object_name),
                d_dim(dim),
                d_num_species(num_species),
                d_mixing_closure_model(mixing_closure_model),
                d_equation_of_state_mixing_rules_db(equation_of_state_mixing_rules_db)
        {}
        
        /*
         * Print all characteristics of the equation of state mixing rules class.
         */
        virtual void
        printClassData(std::ostream& os) const = 0;
        
        /*
         * Put the characteristics of the equation of state mixing rules class into the restart
         * database.
         */
        virtual void
        putToRestart(
            const boost::shared_ptr<tbox::Database>& restart_db) const = 0;
        
        /*
         * Get the number of thermodynamic properties of the mixture.
         */
        virtual int
        getNumberOfMixtureThermodynamicProperties() const = 0;
        
        /*
         * Get the number of thermodynamic properties of a species.
         */
        virtual int
        getNumberOfSpeciesThermodynamicProperties() const = 0;
        
        /*
         * Get the thermodynamic properties of the mixture.
         */
        virtual void
        getMixtureThermodynamicProperties(
            std::vector<double*>& mixture_thermo_properties,
            const std::vector<const double*>& species_fraction) const = 0;
        
        /*
         * Get the thermodynamic properties of a species.
         */
        virtual void
        getSpeciesThermodynamicProperties(
            std::vector<double*>& species_thermo_properties,
            const int& species_index) const = 0;
        
        /*
         * Helper function to compute the mixture density given the partial densities.
         */
        double
        getMixtureDensity(
            const std::vector<const double*>& partial_density) const
        {
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
            if (static_cast<int>(partial_density.size()) != d_num_species)
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "Number of partial densities provided is not"
                    << "equal to the total number of species."
                    << std::endl);
            }
#endif
            
            double mixture_density = 0.0;
            
            for (int si = 0; si < d_num_species; si++)
            {
                const double& Z_rho = *(partial_density[si]);
                
                mixture_density += Z_rho;
            }
            
            return mixture_density;            
        }
        
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
         * boost::shared_ptr to the database of equation of state mixing rules.
         */
        const boost::shared_ptr<tbox::Database> d_equation_of_state_mixing_rules_db;
        
};

#endif /* EQUATION_OF_STATE_MIXING_RULES_HPP */
