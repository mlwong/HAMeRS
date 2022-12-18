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
    
    double p = 0.0;
    double T = 0.0;
    
    std::vector<double> molecular_properties;
    std::vector<const double*> molecular_properties_const_ptr;
    
    /*
     * Verify that the equations of shear viscosity are implemented correctly.
     */
    
    double mu = 0.0;
    
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
    
    molecular_properties[0] = 1.81e-5;
    molecular_properties[1] = 29.0;
    
    p = 1.0e5;
    T = 300.0;
    
    mu = equation_of_shear_viscosity->getShearViscosity(
        &p,
        &T,
        molecular_properties_const_ptr);
    
    if (fabs(mu - 1.81e-5) < 1.0e-8)
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
    
    molecular_properties[0] = 78.6;
    molecular_properties[1] = 3.711;
    molecular_properties[2] = 29.0;
    
    p = 1.0e5;
    T = 300.0;
    
    mu = equation_of_shear_viscosity->getShearViscosity(
        &p,
        &T,
        molecular_properties_const_ptr);
    
    if (fabs(mu - 1.8461366680483000e-5) < 1.0e-8)
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
    
    double mu_v = 0.0;
    
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
    
    molecular_properties[0] = 1.59e-5;
    molecular_properties[1] = 29.0;
    
    p = 1.0e5;
    T = 300.0;
    
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
    
    molecular_properties[0] = 1.4;
    
    molecular_properties[1] = -3.15e-5;
    molecular_properties[2] = 1.58e-7;
    
    molecular_properties[3] = 0.0;
    molecular_properties[4] = 0.0;
    molecular_properties[5] = 0.0;
    molecular_properties[6] = 0.0;
    
    molecular_properties[7] = 29.0;
    
    p = 1.0e5;
    T = 300.0;
    
    mu_v = equation_of_bulk_viscosity->getBulkViscosity(
        &p,
        &T,
        molecular_properties_const_ptr);
    
    if (fabs(mu_v - 1.59e-5) < 1.0e-8)
    {
        std::cout << "EquationOfBulkViscosityCramer is implemented correctly!" << std::endl;
    }
    else
    {
        std::cout << "EquationOfBulkViscosityCramer is not implemented correctly!" << std::endl;
    }
    
    molecular_properties[0] = 1.1;
    
    molecular_properties[1] = 0.0;
    molecular_properties[2] = 0.0;
    
    molecular_properties[3] = 1.0/(molecular_properties[0] - 1.0) - (3.0 + 3.0)/2.0;
    molecular_properties[4] = 0.2064e-5;
    molecular_properties[5] = 121.0;
    molecular_properties[6] = -339.0;
    
    molecular_properties[7] = 146.0;
    
    p = 1.0e5;
    T = 300.0;
    
    mu_v = equation_of_bulk_viscosity->getBulkViscosity(
        &p,
        &T,
        molecular_properties_const_ptr);
    
    if (fabs(mu_v - 5.30175204375121e-3) < 1.0e-8)
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
    
    double kappa = 0.0;
    
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
    
    molecular_properties[0] = 2.71e-2;
    molecular_properties[1] = 29.0;
    
    p = 1.0e5;
    T = 300.0;
    
    kappa = equation_of_thermal_conductivity->getThermalConductivity(
        &p,
        &T,
        molecular_properties_const_ptr);
    
    if (fabs(kappa - 2.71e-2) < 1.0e-8)
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
    
    molecular_properties[0] = 1005.0;
    molecular_properties[1] = 0.72;
    molecular_properties[2] = 29.0;
    molecular_properties[3] = 78.6;
    molecular_properties[4] = 3.711;
    molecular_properties[5] = 29.0;
    
    p = 1.0e5;
    T = 300.0;
    
    kappa = equation_of_thermal_conductivity->getThermalConductivity(
        &p,
        &T,
        molecular_properties_const_ptr);
    
    if (fabs(kappa - 2.5768990991507500e-2) < 1.0e-8)
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
    
    std::vector<double> D;
    std::vector<double*> D_ptr;
    
    std::vector<double> Y;
    std::vector<const double*> Y_const_ptr;
    
    HAMERS_SHARED_PTR<EquationOfMassDiffusivityMixingRules> equation_of_mass_diffusivity_mixing_rules;
    
    HAMERS_SHARED_PTR<tbox::Database> equation_of_mass_diffusivity_mixing_rules_db;
    
    equation_of_mass_diffusivity_mixing_rules_db.reset(new
        tbox::MemoryDatabase("equation_of_mass_diffusivity_mixing_rules_db"));
    
    std::vector<double> species_D;
    species_D.reserve(2);
    species_D.push_back(1.85e-5);
    species_D.push_back(1.85e-5);
    
    equation_of_mass_diffusivity_mixing_rules_db->putDoubleVector(
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
    
    p = 23000;
    T = 298;
    
    Y[0] = 0.1;
    Y[1] = 0.9;
    
    equation_of_mass_diffusivity_mixing_rules->getMassDiffusivities(
        D_ptr,
        &p,
        &T,
        Y_const_ptr);
    
    if (fabs(D[1] - 1.85e-5) < 1.0e-8)
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
    
    std::vector<double> species_epsilon_by_k;
    species_epsilon_by_k.reserve(2);
    species_epsilon_by_k.push_back(212.0);
    species_epsilon_by_k.push_back(458.0);
    
    std::vector<double> species_sigma;
    species_sigma.reserve(2);
    species_sigma.push_back(5.199);
    species_sigma.push_back(4.599);
    
    std::vector<double> species_M;
    species_M.reserve(2);
    species_M.push_back(146.057);
    species_M.push_back(58.0805);
    
    equation_of_mass_diffusivity_mixing_rules_db->putDoubleVector(
        "species_epsilon_by_k",
        species_epsilon_by_k);
    
    equation_of_mass_diffusivity_mixing_rules_db->putDoubleVector(
        "species_sigma",
        species_sigma);
    
    equation_of_mass_diffusivity_mixing_rules_db->putDoubleVector(
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
    
    p = 23000;
    T = 298;
    
    Y[0] = 0.1;
    Y[1] = 0.9;
    
    equation_of_mass_diffusivity_mixing_rules->getMassDiffusivities(
        D_ptr,
        &p,
        &T,
        Y_const_ptr);
    
    if (fabs(D[1] - 1.84656024212537e-5) < 1.0e-8)
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
