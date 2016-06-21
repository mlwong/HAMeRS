#ifndef EQUATION_OF_STATE_MANAGER_HPP
#define EQUATION_OF_STATE_MANAGER_HPP

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/Dimension.h"

#include "util/equations_of_state/EquationsOfState.hpp"

#include "boost/shared_ptr.hpp"
#include <string>

using namespace SAMRAI;

class EquationOfStateManager
{
    public:
        EquationOfStateManager(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const int& num_species,
            const boost::shared_ptr<tbox::Database>& equation_of_state_db,
            const std::string& equation_of_state_str);
        
        /*
         * Get the label of equation of state.
         */
        const EQUATION_OF_STATE_LABEL&
        getEquationOfStateLabel() const
        {
            return d_equation_of_state_label;
        }
        
        /*
         * Get the equation of state.
         */
        boost::shared_ptr<EquationOfState>
        getEquationOfState() const
        {
            return d_equation_of_state;
        }
        
        /*
         * Print all characteristics of equation of state manager.
         */
        void
        printClassData(std::ostream& os) const;
        
    private:
        /*
         * The object name is used for error/warning reporting.
         */
        const std::string d_object_name;
        
        /*
         * Label of equation of state.
         */
        EQUATION_OF_STATE_LABEL d_equation_of_state_label;
        
        /*
         * boost::shared_ptr to the equation of state.
         */
        boost::shared_ptr<EquationOfState> d_equation_of_state;
        
};

#endif /* EQUATION_OF_STATE_MANAGER_HPP */
