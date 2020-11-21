#include "util/mixing_rules/equations_of_bulk_viscosity/EquationOfBulkViscosityMixingRulesManager.hpp"

EquationOfBulkViscosityMixingRulesManager::EquationOfBulkViscosityMixingRulesManager(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const int& num_species,
    const MIXING_CLOSURE_MODEL::TYPE& mixing_closure_model,
    const HAMERS_SHARED_PTR<tbox::Database>& equation_of_bulk_viscosity_mixing_rules_db,
    const std::string& equation_of_bulk_viscosity_str):
        d_object_name(object_name)
{
    TBOX_ASSERT(!object_name.empty());
    
    if (equation_of_bulk_viscosity_str == "CONSTANT")
    {
        d_equation_of_bulk_viscosity_type = EQN_BULK_VISCOSITY::CONSTANT;
        
        d_equation_of_bulk_viscosity_mixing_rules.reset(new EquationOfBulkViscosityMixingRulesConstant(
                "d_equation_of_bulk_viscosity_mixing_rules",
                dim,
                num_species,
                mixing_closure_model,
                equation_of_bulk_viscosity_mixing_rules_db));
    }
    else if (equation_of_bulk_viscosity_str == "CRAMER")
    {
        d_equation_of_bulk_viscosity_type = EQN_BULK_VISCOSITY::CRAMER;
        
        d_equation_of_bulk_viscosity_mixing_rules.reset(new EquationOfBulkViscosityMixingRulesCramer(
                "d_equation_of_bulk_viscosity_mixing_rules",
                dim,
                num_species,
                mixing_closure_model,
                equation_of_bulk_viscosity_mixing_rules_db));
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Unknown equation_of_bulk_viscosity/d_equation_of_bulk_viscosity string = '"
            << equation_of_bulk_viscosity_str
            << "' found in input/restart file."
            << std::endl);        
    }
}


/*
 * Print all characteristics of equation of bulk viscosity manager.
 */
void
EquationOfBulkViscosityMixingRulesManager::printClassData(std::ostream& os) const
{
    os << "\nPrint EquationOfBulkViscosityMixingRulesManager object..."
       << std::endl;
    
    os << std::endl;
    
    os << "EquationOfBulkViscosityMixingRulesManager: this = "
       << (EquationOfBulkViscosityMixingRulesManager *)this
       << std::endl;
    
    os << "d_object_name = "
       << d_object_name
       << std::endl;
    
    os << "d_equation_of_bulk_viscosity_type = "
       << d_equation_of_bulk_viscosity_type
       << std::endl;
    
    os << "................................................................................";
    d_equation_of_bulk_viscosity_mixing_rules->printClassData(os);
}
