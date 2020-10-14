#include "util/mixing_rules/equations_of_thermal_conductivity/EquationOfThermalConductivityMixingRulesManager.hpp"

EquationOfThermalConductivityMixingRulesManager::EquationOfThermalConductivityMixingRulesManager(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const int& num_species,
    const MIXING_CLOSURE_MODEL::TYPE& mixing_closure_model,
    const HAMERS_SHARED_PTR<tbox::Database>& equation_of_thermal_conductivity_mixing_rules_db,
    const std::string& equation_of_thermal_conductivity_str):
        d_object_name(object_name)
{
    TBOX_ASSERT(!object_name.empty());
    
    if (equation_of_thermal_conductivity_str == "CONSTANT")
    {
        d_equation_of_thermal_conductivity_type = EQN_THERMAL_CONDUCTIVITY::CONSTANT;
        
        d_equation_of_thermal_conductivity_mixing_rules.reset(new EquationOfThermalConductivityMixingRulesConstant(
                "d_equation_of_thermal_conductivity_mixing_rules",
                dim,
                num_species,
                mixing_closure_model,
                equation_of_thermal_conductivity_mixing_rules_db));
    }
    else if (equation_of_thermal_conductivity_str == "PRANDTL")
    {
        d_equation_of_thermal_conductivity_type = EQN_THERMAL_CONDUCTIVITY::PRANDTL;
        
        d_equation_of_thermal_conductivity_mixing_rules.reset(new EquationOfThermalConductivityMixingRulesPrandtl(
                "d_equation_of_thermal_conductivity_mixing_rules",
                dim,
                num_species,
                mixing_closure_model,
                equation_of_thermal_conductivity_mixing_rules_db));
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Unknown equation_of_thermal_conductivity/d_equation_of_thermal_conductivity string = '"
            << equation_of_thermal_conductivity_str
            << "' found in input/restart file."
            << std::endl);        
    }
}


/*
 * Print all characteristics of equation of thermal conductivity manager.
 */
void
EquationOfThermalConductivityMixingRulesManager::printClassData(std::ostream& os) const
{
    os << "\nPrint EquationOfThermalConductivityMixingRulesManager object..."
       << std::endl;
    
    os << std::endl;
    
    os << "EquationOfThermalConductivityMixingRulesManager: this = "
       << (EquationOfThermalConductivityMixingRulesManager *)this
       << std::endl;
    
    os << "d_object_name = "
       << d_object_name
       << std::endl;
    
    os << "d_equation_of_thermal_conductivity_type = "
       << d_equation_of_thermal_conductivity_type
       << std::endl;
    
    os << "................................................................................";
    d_equation_of_thermal_conductivity_mixing_rules->printClassData(os);
}
