#include "util/equations_of_state/EquationOfStateManager.hpp"

EquationOfStateManager::EquationOfStateManager(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const int& num_species,
    const boost::shared_ptr<tbox::Database>& equation_of_state_db,
    const std::string& equation_of_state_str):
        d_object_name(object_name)
{
    TBOX_ASSERT(!object_name.empty());
    TBOX_ASSERT(equation_of_state_db);
    
    if (equation_of_state_str == "IDEAL_GAS")
    {
        d_equation_of_state_label = IDEAL_GAS;
        
        d_equation_of_state.reset(new EquationOfStateIdealGas(
            "d_equation_of_state",
            dim,
            num_species,
            equation_of_state_db));
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Unknown equation_of_state/d_equation_of_state string = '"
            << equation_of_state_str
            << "' found in input/restart file."
            << std::endl);        
    }
}


/*
 * Print all characteristics of equation of state manager.
 */
void
EquationOfStateManager::printClassData(std::ostream& os) const
{
    os << "\nPrint EquationOfStateManager object..."
       << std::endl;
    
    os << std::endl;
    
    os << "EquationOfStateManager: this = "
       << (EquationOfStateManager *)this
       << std::endl;
    
    os << "d_object_name = "
       << d_object_name
       << std::endl;
    
    os << "d_equation_of_state_label = "
       << d_equation_of_state_label
       << std::endl;
}
