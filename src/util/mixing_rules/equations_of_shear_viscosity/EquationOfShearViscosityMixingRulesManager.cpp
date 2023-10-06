#include "util/mixing_rules/equations_of_shear_viscosity/EquationOfShearViscosityMixingRulesManager.hpp"

EquationOfShearViscosityMixingRulesManager::EquationOfShearViscosityMixingRulesManager(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const int& num_species,
    const MIXING_CLOSURE_MODEL::TYPE& mixing_closure_model,
    const HAMERS_SHARED_PTR<tbox::Database>& equation_of_shear_viscosity_mixing_rules_db,
    const std::string& equation_of_shear_viscosity_str):
        d_object_name(object_name)
{
    TBOX_ASSERT(!object_name.empty());
    
    if (equation_of_shear_viscosity_str == "CONSTANT")
    {
        d_equation_of_shear_viscosity_type = EQN_SHEAR_VISCOSITY::CONSTANT;
        
        d_equation_of_shear_viscosity_mixing_rules.reset(new EquationOfShearViscosityMixingRulesConstant(
                "d_equation_of_shear_viscosity_mixing_rules",
                dim,
                num_species,
                mixing_closure_model,
                equation_of_shear_viscosity_mixing_rules_db));
    }
    else if (equation_of_shear_viscosity_str == "POWER_LAW")
    {
        d_equation_of_shear_viscosity_type = EQN_SHEAR_VISCOSITY::POWER_LAW;
        
        d_equation_of_shear_viscosity_mixing_rules.reset(new EquationOfShearViscosityMixingRulesPowerLaw(
                "d_equation_of_shear_viscosity_mixing_rules",
                dim,
                num_species,
                mixing_closure_model,
                equation_of_shear_viscosity_mixing_rules_db));
    }
    else if (equation_of_shear_viscosity_str == "CHAPMAN_ENSKOG")
    {
        d_equation_of_shear_viscosity_type = EQN_SHEAR_VISCOSITY::CHAPMAN_ENSKOG;
        
        d_equation_of_shear_viscosity_mixing_rules.reset(new EquationOfShearViscosityMixingRulesChapmanEnskog(
                "d_equation_of_shear_viscosity_mixing_rules",
                dim,
                num_species,
                mixing_closure_model,
                equation_of_shear_viscosity_mixing_rules_db));
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Unknown equation_of_shear_viscosity/d_equation_of_shear_viscosity string = '"
            << equation_of_shear_viscosity_str
            << "' found in input/restart file."
            << std::endl);
    }
}


/*
 * Print all characteristics of equation of shear viscosity manager.
 */
void
EquationOfShearViscosityMixingRulesManager::printClassData(std::ostream& os) const
{
    os << "\nPrint EquationOfShearViscosityMixingRulesManager object..."
       << std::endl;
    
    os << std::endl;
    
    os << "EquationOfShearViscosityMixingRulesManager: this = "
       << (EquationOfShearViscosityMixingRulesManager *)this
       << std::endl;
    
    os << "d_object_name = "
       << d_object_name
       << std::endl;
    
    os << "d_equation_of_shear_viscosity_type = "
       << d_equation_of_shear_viscosity_type
       << std::endl;
    
    os << "................................................................................";
    d_equation_of_shear_viscosity_mixing_rules->printClassData(os);
}
