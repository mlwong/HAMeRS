#ifndef EQUATION_OF_SHEAR_VISCOSITY_MIXING_RULES_HPP
#define EQUATION_OF_SHEAR_VISCOSITY_MIXING_RULES_HPP

#include "util/mixing_rules/equations_of_shear_viscosity/EquationOfShearViscosity.hpp"
#include "util/mixing_rules/MixingClosureModels.hpp"

#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Dimension.h"

#include <string>
#include <vector>

using namespace SAMRAI;

class EquationOfShearViscosityMixingRules
{
    public:
        EquationOfShearViscosityMixingRules(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const int& num_species,
            const MIXING_CLOSURE_MODEL::TYPE& mixing_closure_model,
            const HAMERS_SHARED_PTR<tbox::Database>& equation_of_shear_viscosity_mixing_rules_db):
                d_object_name(object_name),
                d_dim(dim),
                d_num_species(num_species),
                d_mixing_closure_model(mixing_closure_model),
                d_equation_of_shear_viscosity_mixing_rules_db(equation_of_shear_viscosity_mixing_rules_db)
        {}
        
        virtual ~EquationOfShearViscosityMixingRules() {}
        
        /*
         * Return the HAMERS_SHARED_PTR to the equation of shear viscosity.
         */
        virtual const HAMERS_SHARED_PTR<EquationOfShearViscosity>&
        getEquationOfShearViscosity(const int species_index = 0) const = 0;
        
        /*
         * Print all characteristics of the equation of shear viscosity mixing rules class.
         */
        virtual void
        printClassData(std::ostream& os) const = 0;
        
        /*
         * Put the characteristics of the equation of shear viscosity mixing rules class into the restart
         * database.
         */
        virtual void
        putToRestart(
            const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const = 0;
        
        /*
         * Compute the shear viscosity of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        virtual Real
        getShearViscosity(
            const Real* const pressure,
            const Real* const temperature,
            const std::vector<const Real*>& mass_fractions) const = 0;
        
        /*
         * Compute the shear viscosity of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        void
        computeShearViscosity(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_shear_viscosity,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_temperature,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_mass_fractions) const
        {
            const hier::Box empty_box(d_dim);
            computeShearViscosity(
                data_shear_viscosity,
                data_pressure,
                data_temperature,
                data_mass_fractions,
                empty_box);
        }
        
        /*
         * Compute the shear viscosity of the mixture with isothermal and isobaric equilibrium assumptions.
         */
        virtual void
        computeShearViscosity(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_shear_viscosity,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_temperature,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_mass_fractions,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the shear viscosity of the mixture with isobaric equilibrium assumption.
         */
        virtual Real
        getShearViscosity(
            const Real* const pressure,
            const std::vector<const Real*>& species_temperatures,
            const std::vector<const Real*>& mass_fractions,
            const std::vector<const Real*>& volume_fractions) const = 0;
        
        /*
         * Compute the shear viscosity of the mixture with isobaric equilibrium assumption.
         */
        void
        computeShearViscosity(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_shear_viscosity,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& data_species_temperatures,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_mass_fractions,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_volume_fractions) const
        {
            const hier::Box empty_box(d_dim);
            computeShearViscosity(
                data_shear_viscosity,
                data_pressure,
                data_species_temperatures,
                data_mass_fractions,
                data_volume_fractions,
                empty_box);
            
        }
        
        /*
         * Compute the shear viscosity of the mixture with isobaric equilibrium assumption.
         */
        virtual void
        computeShearViscosity(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_shear_viscosity,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& data_species_temperatures,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_mass_fractions,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_volume_fractions,
            const hier::Box& domain) const = 0;
        
        /*
         * Get the number of molecular properties of a species.
         */
        virtual int
        getNumberOfSpeciesMolecularProperties(const int species_index = 0) const = 0;
        
        /*
         * Get the molecular properties of a species.
         */
        virtual void
        getSpeciesMolecularProperties(
            std::vector<Real*>& species_molecular_properties,
            const int species_index = 0) const = 0;
        
    protected:
        /*
         * The object name is used for error/warning reporting.
         */
        const std::string d_object_name;
        
        /*
         * Problem dimension.
         */
        const tbox::Dimension d_dim;
        
        /*
         * Number species.
         */
        const int d_num_species;
        
        /*
         * Type of mixing closure model.
         */
        const MIXING_CLOSURE_MODEL::TYPE d_mixing_closure_model;
        
        /*
         * HAMERS_SHARED_PTR to the database of equation of shear viscosity mixing rules.
         */
        const HAMERS_SHARED_PTR<tbox::Database> d_equation_of_shear_viscosity_mixing_rules_db;
        
};
    

#endif /* EQUATION_OF_SHEAR_VISCOSITY_MIXING_RULES_HPP */
