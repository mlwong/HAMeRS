#include "util/mixing_rules/equations_of_bulk_viscosity/constant/EquationOfBulkViscosityConstant.hpp"
#include "util/mixing_rules/equations_of_bulk_viscosity/Cramer/EquationOfBulkViscosityCramer.hpp"
#include "util/mixing_rules/equations_of_mass_diffusivity/constant/EquationOfMassDiffusivityMixingRulesConstant.hpp"
#include "util/mixing_rules/equations_of_mass_diffusivity/Reid/EquationOfMassDiffusivityMixingRulesReid.hpp"
#include "util/mixing_rules/equations_of_shear_viscosity/constant/EquationOfShearViscosityConstant.hpp"
#include "util/mixing_rules/equations_of_shear_viscosity/Chapman-Enskog/EquationOfShearViscosityChapmanEnskog.hpp"
#include "util/mixing_rules/equations_of_thermal_conductivity/constant/EquationOfThermalConductivityConstant.hpp"
#include "util/mixing_rules/equations_of_thermal_conductivity/Prandtl/EquationOfThermalConductivityPrandtl.hpp"

#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/InputDatabase.h"

using namespace SAMRAI;

int main(int argc, char *argv[])
{
    /*
     * Set the dimension.
     */
    tbox::Dimension dim(3);
    
    Real p = Real(0);
    Real T = Real(0);
    
    std::vector<Real> molecular_properties;
    std::vector<const Real*> molecular_properties_const_ptr;
    
    /*
     * Verify that the equations of shear viscosity are implemented correctly.
     */
    
    Real mu = Real(0);
    
    HAMERS_SHARED_PTR<EquationOfShearViscosity> equation_of_shear_viscosity;
    
    const bool use_constant_kinematic_viscosity_and_ideal_gas_assumptions = false;
    
    equation_of_shear_viscosity.reset(new EquationOfShearViscosityConstant(
        "equation_of_shear_viscosity",
        dim,
        use_constant_kinematic_viscosity_and_ideal_gas_assumptions));
    
    molecular_properties.resize(2);
    molecular_properties_const_ptr.reserve(2);
    for (int i = 0; i < 2; i++)
    {
        molecular_properties_const_ptr.push_back(&molecular_properties[i]);
    }
    
    molecular_properties[0] = Real(1.81e-5);
    molecular_properties[1] = Real(29.0);
    
    p = Real(1.0e5);
    T = Real(300);
    
    mu = equation_of_shear_viscosity->getShearViscosity(
        &p,
        &T,
        molecular_properties_const_ptr);
    
    if (std::abs(mu - Real(1.81e-5)) < Real(1.0e-8))
    {
        std::cout << "EquationOfShearViscosityConstant is implemented correctly!" << std::endl;
    }
    else
    {
        std::cout << "EquationOfShearViscosityConstant is not implemented correctly!" << std::endl;
    }
    
    molecular_properties.clear();
    molecular_properties_const_ptr.clear();
    
    equation_of_shear_viscosity.reset();
    
    equation_of_shear_viscosity.reset(new EquationOfShearViscosityChapmanEnskog(
        "equation_of_shear_viscosity",
        dim));
    
    molecular_properties.resize(3);
    molecular_properties_const_ptr.reserve(3);
    for (int i = 0; i < 3; i++)
    {
        molecular_properties_const_ptr.push_back(&molecular_properties[i]);
    }
    
    molecular_properties[0] = Real(78.6);
    molecular_properties[1] = Real(3.711);
    molecular_properties[2] = Real(29);
    
    p = Real(1.0e5);
    T = Real(300);
    
    mu = equation_of_shear_viscosity->getShearViscosity(
        &p,
        &T,
        molecular_properties_const_ptr);
    
    if (std::abs(mu - Real(1.8461366680483000e-5)) < Real(1.0e-8))
    {
        std::cout << "EquationOfShearViscosityChapmanEnskog is implemented correctly!" << std::endl;
    }
    else
    {
        std::cout << "EquationOfShearViscosityChapmanEnskog is not implemented correctly!" << std::endl;
    }
    
    molecular_properties.clear();
    molecular_properties_const_ptr.clear();
    
    equation_of_shear_viscosity.reset();
    
    /*
     * Verify that the equations of bulk viscosity are implemented correctly.
     */
    
    Real mu_v = Real(0);
    
    HAMERS_SHARED_PTR<EquationOfBulkViscosity> equation_of_bulk_viscosity;
    
    equation_of_bulk_viscosity.reset(new EquationOfBulkViscosityConstant(
        "equation_of_bulk_viscosity",
        dim));
    
    molecular_properties.resize(2);
    molecular_properties_const_ptr.reserve(2);
    for (int i = 0; i < 2; i++)
    {
        molecular_properties_const_ptr.push_back(&molecular_properties[i]);
    }
    
    molecular_properties[0] = Real(1.59e-5);
    molecular_properties[1] = Real(29);
    
    p = Real(1.0e5);
    T = Real(300);
    
    mu_v = equation_of_bulk_viscosity->getBulkViscosity(
        &p,
        &T,
        molecular_properties_const_ptr);
    
    if (fabs(mu_v - 1.59e-5) < 1.0e-8)
    {
        std::cout << "EquationOfBulkViscosityConstant is implemented correctly!" << std::endl;
    }
    else
    {
        std::cout << "EquationOfBulkViscosityConstant is not implemented correctly!" << std::endl;
    }
    
    molecular_properties.clear();
    molecular_properties_const_ptr.clear();
    
    equation_of_bulk_viscosity.reset();
    
    equation_of_bulk_viscosity.reset(new EquationOfBulkViscosityCramer(
        "equation_of_bulk_viscosity",
        dim));
    
    molecular_properties.resize(8);
    molecular_properties_const_ptr.reserve(8);
    for (int i = 0; i < 8; i++)
    {
        molecular_properties_const_ptr.push_back(&molecular_properties[i]);
    }
    
    molecular_properties[0] = Real(1.4);
    
    molecular_properties[1] = Real(-3.15e-5);
    molecular_properties[2] = Real(1.58e-7);
    
    molecular_properties[3] = Real(0);
    molecular_properties[4] = Real(0);
    molecular_properties[5] = Real(0);
    molecular_properties[6] = Real(0);
    
    molecular_properties[7] = Real(29);
    
    p = Real(1.0e5);
    T = Real(300);
    
    mu_v = equation_of_bulk_viscosity->getBulkViscosity(
        &p,
        &T,
        molecular_properties_const_ptr);
    
    if (std::abs(mu_v - Real(1.59e-5)) < Real(1.0e-8))
    {
        std::cout << "EquationOfBulkViscosityCramer is implemented correctly!" << std::endl;
    }
    else
    {
        std::cout << "EquationOfBulkViscosityCramer is not implemented correctly!" << std::endl;
    }
    
    molecular_properties[0] = Real(1.1);
    
    molecular_properties[1] = Real(0);
    molecular_properties[2] = Real(0);
    
    molecular_properties[3] = Real(1)/(molecular_properties[0] - Real(1)) - (Real(3) + Real(3))/Real(2);
    molecular_properties[4] = Real(0.2064e-5);
    molecular_properties[5] = Real(121);
    molecular_properties[6] = -Real(339);
    
    molecular_properties[7] = Real(146.0);
    
    p = Real(1.0e5);
    T = Real(300);
    
    mu_v = equation_of_bulk_viscosity->getBulkViscosity(
        &p,
        &T,
        molecular_properties_const_ptr);
    
    if (std::abs(mu_v - Real(5.30175204375121e-3)) < Real(1.0e-8))
    {
        std::cout << "EquationOfBulkViscosityCramer is implemented correctly!" << std::endl;
    }
    else
    {
        std::cout << "EquationOfBulkViscosityCramer is not implemented correctly!" << std::endl;
    }
    
    molecular_properties.clear();
    molecular_properties_const_ptr.clear();
    
    equation_of_bulk_viscosity.reset();
    
    /*
     * Verify that the equations of thermal conductivity are implemented correctly.
     */
    
    Real kappa = Real(0);
    
    HAMERS_SHARED_PTR<EquationOfThermalConductivity> equation_of_thermal_conductivity;
    
    equation_of_thermal_conductivity.reset(new EquationOfThermalConductivityConstant(
        "equation_of_thermal_conductivity",
        dim));
    
    molecular_properties.resize(2);
    molecular_properties_const_ptr.reserve(2);
    for (int i = 0; i < 2; i++)
    {
        molecular_properties_const_ptr.push_back(&molecular_properties[i]);
    }
    
    molecular_properties[0] = Real(2.71e-2);
    molecular_properties[1] = Real(29);
    
    p = Real(1.0e5);
    T = Real(300);
    
    kappa = equation_of_thermal_conductivity->getThermalConductivity(
        &p,
        &T,
        molecular_properties_const_ptr);
    
    if (std::abs(kappa - Real(2.71e-2)) < Real(1.0e-8))
    {
        std::cout << "EquationOfThermalConductivityConstant is implemented correctly!" << std::endl;
    }
    else
    {
        std::cout << "EquationOfThermalConductivityConstant is not implemented correctly!" << std::endl;
    }
    
    molecular_properties.clear();
    molecular_properties_const_ptr.clear();
    
    equation_of_shear_viscosity.reset();
    equation_of_thermal_conductivity.reset();
    
    equation_of_shear_viscosity.reset(new EquationOfShearViscosityChapmanEnskog(
        "equation_of_shear_viscosity",
        dim));
    
    equation_of_thermal_conductivity.reset(new EquationOfThermalConductivityPrandtl(
        "equation_of_thermal_conductivity",
        dim,
        equation_of_shear_viscosity));
    
    molecular_properties.resize(6);
    molecular_properties_const_ptr.reserve(6);
    for (int i = 0; i < 6; i++)
    {
        molecular_properties_const_ptr.push_back(&molecular_properties[i]);
    }
    
    molecular_properties[0] = Real(1005);
    molecular_properties[1] = Real(0.72);
    molecular_properties[2] = Real(29);
    molecular_properties[3] = Real(78.6);
    molecular_properties[4] = Real(3.711);
    molecular_properties[5] = Real(29);
    
    p = Real(1.0e5);
    T = Real(300);
    
    kappa = equation_of_thermal_conductivity->getThermalConductivity(
        &p,
        &T,
        molecular_properties_const_ptr);
    
    if (std::abs(kappa - Real(2.5768990991507500e-2)) < Real(1.0e-8))
    {
        std::cout << "EquationOfThermalConductivityPrandtl is implemented correctly!" << std::endl;
    }
    else
    {
        std::cout << "EquationOfThermalConductivityPrandtl is not implemented correctly!" << std::endl;
    }
    
    molecular_properties.clear();
    molecular_properties_const_ptr.clear();
    
    equation_of_shear_viscosity.reset();
    equation_of_thermal_conductivity.reset();
    
    /*
     * Verify that the equations of mass diffusivity are implemented correctly.
     */
    
    std::vector<Real> D;
    std::vector<Real*> D_ptr;
    
    std::vector<Real> Y;
    std::vector<const Real*> Y_const_ptr;
    
    HAMERS_SHARED_PTR<EquationOfMassDiffusivityMixingRules> equation_of_mass_diffusivity_mixing_rules;
    
    HAMERS_SHARED_PTR<tbox::Database> equation_of_mass_diffusivity_mixing_rules_db;
    
    equation_of_mass_diffusivity_mixing_rules_db.reset(new
        tbox::MemoryDatabase("equation_of_mass_diffusivity_mixing_rules_db"));
    
    std::vector<Real> species_D;
    species_D.reserve(2);
    species_D.push_back(Real(1.85e-5));
    species_D.push_back(Real(1.85e-5));
    
    equation_of_mass_diffusivity_mixing_rules_db->putRealVector(
        "species_D",
        species_D);
    
    equation_of_mass_diffusivity_mixing_rules.reset(new EquationOfMassDiffusivityMixingRulesConstant(
        "equation_of_mass_diffusivity_mixing_rules",
        dim,
        2,
        MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC,
        equation_of_mass_diffusivity_mixing_rules_db));
    
    D.resize(2);
    D_ptr.reserve(2);
    for (int i = 0; i < 2; i++)
    {
        D_ptr.push_back(&D[i]);
    }
    
    Y.resize(2);
    Y_const_ptr.reserve(2);
    for (int i = 0; i < 2; i++)
    {
        Y_const_ptr.push_back(&Y[i]);
    }
    
    p = Real(23000);
    T = Real(298);
    
    Y[0] = Real(0.1);
    Y[1] = Real(0.9);
    
    equation_of_mass_diffusivity_mixing_rules->getMassDiffusivities(
        D_ptr,
        &p,
        &T,
        Y_const_ptr);
    
    if (fabs(D[1] - Real(1.85e-5)) < Real(1.0e-8))
    {
        std::cout << "EquationOfMassDiffusivityMixingRulesConstant is implemented correctly!" << std::endl;
    }
    else
    {
        std::cout << "EquationOfMassDiffusivityMixingRulesConstant is not implemented correctly!" << std::endl;
    }
    
    D.clear();
    D_ptr.clear();
    Y.clear();
    Y_const_ptr.clear();
    
    equation_of_mass_diffusivity_mixing_rules_db.reset();
    equation_of_mass_diffusivity_mixing_rules.reset();
    
    equation_of_mass_diffusivity_mixing_rules_db.reset(new
        tbox::MemoryDatabase("equation_of_mass_diffusivity_mixing_rules_db"));
    
    std::vector<Real> species_epsilon_by_k;
    species_epsilon_by_k.reserve(2);
    species_epsilon_by_k.push_back(Real(212));
    species_epsilon_by_k.push_back(Real(458));
    
    std::vector<Real> species_sigma;
    species_sigma.reserve(2);
    species_sigma.push_back(Real(5.199));
    species_sigma.push_back(Real(4.599));
    
    std::vector<Real> species_M;
    species_M.reserve(2);
    species_M.push_back(Real(146.057));
    species_M.push_back(Real(58.0805));
    
    equation_of_mass_diffusivity_mixing_rules_db->putRealVector(
        "species_epsilon_by_k",
        species_epsilon_by_k);
    
    equation_of_mass_diffusivity_mixing_rules_db->putRealVector(
        "species_sigma",
        species_sigma);
    
    equation_of_mass_diffusivity_mixing_rules_db->putRealVector(
        "species_M",
        species_M);
    
    equation_of_mass_diffusivity_mixing_rules.reset(new EquationOfMassDiffusivityMixingRulesReid(
        "equation_of_mass_diffusivity_mixing_rules",
        dim,
        2,
        MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC,
        equation_of_mass_diffusivity_mixing_rules_db));
    
    D.resize(2);
    D_ptr.reserve(2);
    for (int i = 0; i < 2; i++)
    {
        D_ptr.push_back(&D[i]);
    }
    
    Y.resize(2);
    Y_const_ptr.reserve(2);
    for (int i = 0; i < 2; i++)
    {
        Y_const_ptr.push_back(&Y[i]);
    }
    
    p = Real(23000);
    T = Real(298);
    
    Y[0] = Real(0.1);
    Y[1] = Real(0.9);
    
    equation_of_mass_diffusivity_mixing_rules->getMassDiffusivities(
        D_ptr,
        &p,
        &T,
        Y_const_ptr);
    
    if (fabs(D[1] - Real(1.84656024212537e-5)) < Real(1.0e-8))
    {
        std::cout << "EquationOfMassDiffusivityMixingRulesReid is implemented correctly!" << std::endl;
    }
    else
    {
        std::cout << "EquationOfMassDiffusivityMixingRulesReid is not implemented correctly!" << std::endl;
    }
    
    D.clear();
    Y.clear();
    
    equation_of_mass_diffusivity_mixing_rules_db.reset();
    equation_of_mass_diffusivity_mixing_rules.reset();
    
    return 0;
}
