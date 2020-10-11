#include "util/mixing_rules/equations_of_mass_diffusivity/EquationOfMassDiffusivityMixingRulesManager.hpp"

EquationOfMassDiffusivityMixingRulesManager::EquationOfMassDiffusivityMixingRulesManager(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const int& num_species,
    const MIXING_CLOSURE_MODEL::TYPE& mixing_closure_model,
    const HAMERS_SHARED_PTR<tbox::Database>& equation_of_mass_diffusivity_mixing_rules_db,
    const std::string& equation_of_mass_diffusivity_str):
        d_object_name(object_name)
{
    TBOX_ASSERT(!object_name.empty());
    
    if (equation_of_mass_diffusivity_str == "CONSTANT")
    {
        d_equation_of_mass_diffusivity_type = EQN_MASS_DIFFUSIVITY::CONSTANT;
        
        d_equation_of_mass_diffusivity_mixing_rules.reset(new EquationOfMassDiffusivityMixingRulesConstant(
                "d_equation_of_mass_diffusivity_mixing_rules",
                dim,
                num_species,
                mixing_closure_model,
                equation_of_mass_diffusivity_mixing_rules_db));
    }
    else if (equation_of_mass_diffusivity_str == "REID")
    {
        d_equation_of_mass_diffusivity_type = EQN_MASS_DIFFUSIVITY::REID;
        
        d_equation_of_mass_diffusivity_mixing_rules.reset(new EquationOfMassDiffusivityMixingRulesReid(
                "d_equation_of_mass_diffusivity_mixing_rules",
                dim,
                num_species,
                mixing_closure_model,
                equation_of_mass_diffusivity_mixing_rules_db));
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Unknown equation_of_mass_diffusivity/d_equation_of_mass_diffusivity string = '"
            << equation_of_mass_diffusivity_str
            << "' found in input/restart file."
            << std::endl);        
    }
}


/*
 * Print all characteristics of equation of shear viscosity manager.
 */
void
EquationOfMassDiffusivityMixingRulesManager::printClassData(std::ostream& os) const
{
    os << "\nPrint EquationOfMassDiffusivityMixingRulesManager object..."
       << std::endl;
    
    os << std::endl;
    
    os << "EquationOfMassDiffusivityMixingRulesManager: this = "
       << (EquationOfMassDiffusivityMixingRulesManager *)this
       << std::endl;
    
    os << "d_object_name = "
       << d_object_name
       << std::endl;
    
    os << "d_equation_of_mass_diffusivity_type = "
       << d_equation_of_mass_diffusivity_type
       << std::endl;
    
    os << "................................................................................";
    d_equation_of_mass_diffusivity_mixing_rules->printClassData(os);
}
