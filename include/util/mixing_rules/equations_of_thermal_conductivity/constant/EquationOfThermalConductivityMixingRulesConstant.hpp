#ifndef EQUATION_OF_THERMAL_CONDUCTIVITY_MIXING_RULES_CONSTANT_HPP
#define EQUATION_OF_THERMAL_CONDUCTIVITY_MIXING_RULES_CONSTANT_HPP

#include "util/mixing_rules/equations_of_thermal_conductivity/EquationOfThermalConductivityMixingRules.hpp"

#include "util/mixing_rules/equations_of_thermal_conductivity/constant/EquationOfThermalConductivityConstant.hpp"

class EquationOfThermalConductivityMixingRulesConstant: public EquationOfThermalConductivityMixingRules
{
    public:
        EquationOfThermalConductivityMixingRulesConstant(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const int& num_species,
            const MIXING_CLOSURE_MODEL::TYPE& mixing_closure_model,
            const HAMERS_SHARED_PTR<tbox::Database>& equation_of_thermal_conductivity_mixing_rules_db);
        
        ~EquationOfThermalConductivityMixingRulesConstant() {}
        
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
            return 2;
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
         * Thermal conductivities of different species.
         */
        std::vector<Real> d_species_kappa;
        
        /*
         * Molecular weight of different species.
         */
        std::vector<Real> d_species_M;
        
        /*
         * HAMERS_SHARED_PTR to EquationOfThermalConductivity.
         */
        HAMERS_SHARED_PTR<EquationOfThermalConductivity> d_equation_of_thermal_conductivity;
        
};

#endif /* EQUATION_OF_THERMAL_CONDUCTIVITY_MIXING_RULES_CONSTANT_HPP */
