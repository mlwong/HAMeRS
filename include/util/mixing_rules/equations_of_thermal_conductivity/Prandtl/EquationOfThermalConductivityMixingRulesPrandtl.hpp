#ifndef EQUATION_OF_THERMAL_CONDUCTIVITY_MIXING_RULES_PRANDTL_HPP
#define EQUATION_OF_THERMAL_CONDUCTIVITY_MIXING_RULES_PRANDTL_HPP

#include "util/mixing_rules/equations_of_thermal_conductivity/EquationOfThermalConductivityMixingRules.hpp"

#include "util/mixing_rules/equations_of_thermal_conductivity/Prandtl/EquationOfThermalConductivityPrandtl.hpp"
#include "util/mixing_rules/equations_of_shear_viscosity/EquationOfShearViscosityMixingRulesManager.hpp"

class EquationOfThermalConductivityMixingRulesPrandtl: public EquationOfThermalConductivityMixingRules
{
    public:
        EquationOfThermalConductivityMixingRulesPrandtl(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const int& num_species,
            const MIXING_CLOSURE_MODEL::TYPE& mixing_closure_model,
            const HAMERS_SHARED_PTR<tbox::Database>& equation_of_thermal_conductivity_mixing_rules_db);
        
        ~EquationOfThermalConductivityMixingRulesPrandtl() {}
        
        /*
         * Return the HAMERS_SHARED_PTR to the equation of thermal conductivity.
         */
        const HAMERS_SHARED_PTR<EquationOfThermalConductivity>&
        getEquationOfThermalConductivity(const int species_index = 0) const
        {
            NULL_USE(species_index);
            return d_equation_of_thermal_conductivity;
        }
        
        /*
         * Print all characteristics of the equation of thermal conductivity mixing rules class.
         */
        void
        printClassData(std::ostream& os) const;
        
        /*
         * Put the characteristics of the equation of thermal conductivity mixing rules class into the restart
         * database.
         */
        void
        putToRestart(
            const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const;
        
        /*
         * Compute the thermal conductivity of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        Real
        getThermalConductivity(
            const Real* const pressure,
            const Real* const temperature,
            const std::vector<const Real*>& mass_fractions) const;
        
        /*
         * Compute the thermal conductivity of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        void
        computeThermalConductivity(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_thermal_conductivity,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_temperature,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_mass_fractions,
            const hier::Box& domain) const;
        
        /*
         * Get the number of molecular properties of a species.
         */
        int
        getNumberOfSpeciesMolecularProperties(const int species_index = 0) const
        {
            NULL_USE(species_index);
            return (3 +
                d_equation_of_shear_viscosity_mixing_rules->
                    getNumberOfSpeciesMolecularProperties());
        }
        
        /*
         * Get the molecular properties of a species.
         */
        void
        getSpeciesMolecularProperties(
            std::vector<Real*>& species_molecular_properties,
            const int species_index = 0) const;
        
    private:
        /*
         * Specific heats at constant pressure of different species.
         */
        std::vector<Real> d_species_c_p;
        
        /*
         * Prandtl number of different species.
         */
        std::vector<Real> d_species_Pr;
        
        /*
         * Molecular weight of different species.
         */
        std::vector<Real> d_species_M;
        
        /*
         * A string variable to describe the equation of shear viscosity used.
         */
        std::string d_equation_of_shear_viscosity_str;
        
        /*
         * HAMERS_SHARED_PTR to EquationOfShearViscosityMixingRules.
         */
        HAMERS_SHARED_PTR<EquationOfShearViscosityMixingRules> d_equation_of_shear_viscosity_mixing_rules;
        
        /*
         * HAMERS_SHARED_PTR to EquationOfShearViscosityMixingRulesManager.
         */
        HAMERS_SHARED_PTR<EquationOfShearViscosityMixingRulesManager> d_equation_of_shear_viscosity_mixing_rules_manager;
        
        /*
         * HAMERS_SHARED_PTR to EquationOfThermalConductivity.
         */
        HAMERS_SHARED_PTR<EquationOfThermalConductivity> d_equation_of_thermal_conductivity;
        
};
    

#endif /* EQUATION_OF_THERMAL_CONDUCTIVITY_MIXING_RULES_PRANDTL_HPP */
