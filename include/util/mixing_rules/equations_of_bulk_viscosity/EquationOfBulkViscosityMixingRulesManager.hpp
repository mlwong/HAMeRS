#ifndef EQUATION_OF_BULK_VISCOSITY_MIXING_RULES_MANAGER_HPP
#define EQUATION_OF_BULK_VISCOSITY_MIXING_RULES_MANAGER_HPP

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/Dimension.h"

#include "util/mixing_rules/equations_of_bulk_viscosity/EquationsOfBulkViscosity.hpp"

#include "boost/shared_ptr.hpp"
#include <string>

using namespace SAMRAI;

class EquationOfBulkViscosityMixingRulesManager
{
    public:
        EquationOfBulkViscosityMixingRulesManager(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const int& num_species,
            const MIXING_CLOSURE_MODEL::TYPE& mixing_closure_model,
            const boost::shared_ptr<tbox::Database>& equation_of_bulk_viscosity_mixing_rules_db,
            const std::string& equation_of_bulk_viscosity_str);
        
        /*
         * Get the type of equation of bulk viscosity.
         */
        const EQN_BULK_VISCOSITY::TYPE&
        getEquationOfBulkViscosityType() const
        {
            return d_equation_of_bulk_viscosity_type;
        }
        
        /*
         * Get the equation of bulk viscosity mixing rules.
         */
        boost::shared_ptr<EquationOfBulkViscosityMixingRules>
        getEquationOfBulkViscosityMixingRules() const
        {
            return d_equation_of_bulk_viscosity_mixing_rules;
        }
        
        /*
         * Print all characteristics of equation of bulk viscosity manager.
         */
        void
        printClassData(std::ostream& os) const;
        
    private:
        /*
         * The object name is used for error/warning reporting.
         */
        const std::string d_object_name;
        
        /*
         * Type of equation of bulk viscosity.
         */
        EQN_BULK_VISCOSITY::TYPE d_equation_of_bulk_viscosity_type;
        
        /*
         * boost::shared_ptr to the equation of bulk viscosity.
         */
        boost::shared_ptr<EquationOfBulkViscosityMixingRules> d_equation_of_bulk_viscosity_mixing_rules;
        
};

#endif /* EQUATION_OF_BULK_VISCOSITY_MIXING_RULES_MANAGER_HPP */
