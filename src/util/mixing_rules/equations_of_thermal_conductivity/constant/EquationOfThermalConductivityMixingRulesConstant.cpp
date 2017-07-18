#include "util/mixing_rules/equations_of_thermal_conductivity/constant/EquationOfThermalConductivityMixingRulesConstant.hpp"

EquationOfThermalConductivityMixingRulesConstant::EquationOfThermalConductivityMixingRulesConstant(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const int& num_species,
    const MIXING_CLOSURE_MODEL::TYPE& mixing_closure_model,
    const boost::shared_ptr<tbox::Database>& equation_of_thermal_conductivity_mixing_rules_db):
        EquationOfThermalConductivityMixingRules(
            object_name,
            dim,
            num_species,
            mixing_closure_model,
            equation_of_thermal_conductivity_mixing_rules_db)
{
    /*
     * Get the thermal conductivity of each species from the database.
     */
    
    if (equation_of_thermal_conductivity_mixing_rules_db->keyExists("species_kappa"))
    {
        size_t species_kappa_array_size =
            equation_of_thermal_conductivity_mixing_rules_db->getArraySize("species_kappa");
        if (static_cast<int>(species_kappa_array_size) == d_num_species)
        {
            d_species_kappa =
                equation_of_thermal_conductivity_mixing_rules_db->getDoubleVector("species_kappa");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'species_kappa' entries must be equal to 'num_species'."
                << std::endl);
        }
    }
    else if (equation_of_thermal_conductivity_mixing_rules_db->keyExists("d_species_kappa"))
    {
        size_t species_kappa_array_size =
            equation_of_thermal_conductivity_mixing_rules_db->getArraySize("d_species_kappa");
        if (static_cast<int>(species_kappa_array_size) == d_num_species)
        {
            d_species_kappa =
                equation_of_thermal_conductivity_mixing_rules_db->getDoubleVector("d_species_kappa");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'd_species_kappa' entries must be equal to 'd_num_species'."
                << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Key data 'species_kappa'/'d_species_kappa'"
            << "not found in data for equation of thermal conductivity mixing rules."
            << std::endl);
    }
    
    /*
     * Get the molecular weight of each species from the database.
     */
    
    if (equation_of_thermal_conductivity_mixing_rules_db->keyExists("species_M"))
    {
        size_t species_M_array_size =
            equation_of_thermal_conductivity_mixing_rules_db->getArraySize("species_M");
        if (static_cast<int>(species_M_array_size) == d_num_species)
        {
            d_species_M =
                equation_of_thermal_conductivity_mixing_rules_db->getDoubleVector("species_M");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'species_M' entries must be equal to 'num_species'."
                << std::endl);
        }
    }
    else if (equation_of_thermal_conductivity_mixing_rules_db->keyExists("d_species_M"))
    {
        size_t species_M_array_size =
            equation_of_thermal_conductivity_mixing_rules_db->getArraySize("d_species_M");
        if (static_cast<int>(species_M_array_size) == d_num_species)
        {
            d_species_M =
                equation_of_thermal_conductivity_mixing_rules_db->getDoubleVector("d_species_M");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'd_species_M' entries must be equal to 'num_species'."
                << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Key data 'species_M'/'d_species_M'"
            << "not found in data for equation of thermal conductivity mixing rules."
            << std::endl);
    }
}


/*
 * Print all characteristics of the equation of thermal conductivity class.
 */
void
EquationOfThermalConductivityMixingRulesConstant::printClassData(
    std::ostream& os) const
{
    os << "\nPrint EquationOfThermalConductivityMixingRulesConstant object..."
       << std::endl;
    
    os << std::endl;
    os << "EquationOfThermalConductivityMixingRulesConstant: this = "
       << (EquationOfThermalConductivityMixingRulesConstant *)this
       << std::endl;
    
    os << "d_object_name = "
       << d_object_name
       << std::endl;
    
    os << "d_mixing_closure_model = "
       << d_mixing_closure_model
       << std::endl;
    
    /*
     * Print the thermal conductivity of each species.
     */
    
    os << "d_species_kappa = ";
    for (int si = 0; si < d_num_species - 1; si++)
    {
        os << d_species_kappa[si] << ", ";
    }
    os << d_species_kappa[d_num_species - 1];
    os << std::endl;
    
    /*
     * Print the molecular weight of each species.
     */
    
    os << "d_species_M = ";
    for (int si = 0; si < d_num_species - 1; si++)
    {
        os << d_species_M[si] << ", ";
    }
    os << d_species_M[d_num_species - 1];
    os << std::endl;
}


/*
 * Put the characteristics of the equation of thermal conductivity mixing rules class into the restart
 * database.
 */
void
EquationOfThermalConductivityMixingRulesConstant::putToRestart(
    const boost::shared_ptr<tbox::Database>& restart_db) const
{
    restart_db->putDoubleVector("d_species_kappa", d_species_kappa);
    restart_db->putDoubleVector("d_species_M", d_species_M);
}


/*
 * Compute the thermal conductivity of the mixture with isothermal and isobaric equilibria assumptions.
 */
double
EquationOfThermalConductivityMixingRulesConstant::getThermalConductivity(
    const double* const pressure,
    const double* const temperature,
    const std::vector<const double*>& mass_fractions) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT((d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC) ||
                (d_mixing_closure_model == MIXING_CLOSURE_MODEL::NO_MODEL && d_num_species == 1));
    TBOX_ASSERT((static_cast<int>(mass_fractions.size()) == d_num_species) ||
                (static_cast<int>(mass_fractions.size()) == d_num_species - 1));
#endif
    
    double kappa = 0.0;
    
    double num = 0.0;
    double den = 0.0;
    
    /*
     * Initialize the container and pointers to the container for the molecular properties
     * of a species.
     */
    
    std::vector<double> species_molecular_properties;
    std::vector<double*> species_molecular_properties_ptr;
    std::vector<const double*> species_molecular_properties_const_ptr;
    
    const int num_molecular_properties = getNumberOfSpeciesMolecularProperties();
    
    species_molecular_properties.resize(num_molecular_properties);
    species_molecular_properties_ptr.reserve(num_molecular_properties);
    species_molecular_properties_const_ptr.reserve(num_molecular_properties);
    
    for (int mi = 0; mi < num_molecular_properties; mi++)
    {
        species_molecular_properties_ptr.push_back(&species_molecular_properties[mi]);
        species_molecular_properties_const_ptr.push_back(&species_molecular_properties[mi]);
    }
    
    if (static_cast<int>(mass_fractions.size()) == d_num_species)
    {
        for (int si = 0; si < d_num_species; si++)
        {
            getSpeciesMolecularProperties(species_molecular_properties_ptr, si);
            
            const double kappa_i = d_equation_of_thermal_conductivity->
                getThermalConductivity(
                    pressure,
                    temperature,
                    species_molecular_properties_const_ptr);
            
            num += kappa_i*(*(mass_fractions[si]))/(sqrt(species_molecular_properties[1]));
            den += *(mass_fractions[si])/(sqrt(species_molecular_properties[1]));
        }
    }
    else if (static_cast<int>(mass_fractions.size()) == d_num_species - 1)
    {
        double Y_last = 1.0;
        
        for (int si = 0; si < d_num_species - 1; si++)
        {
            getSpeciesMolecularProperties(species_molecular_properties_ptr, si);
            
            const double kappa_i = d_equation_of_thermal_conductivity->
                getThermalConductivity(
                    pressure,
                    temperature,
                    species_molecular_properties_const_ptr);
            
            num += kappa_i*(*(mass_fractions[si]))/(sqrt(species_molecular_properties[1]));
            den += *(mass_fractions[si])/(sqrt(species_molecular_properties[1]));
            
            // Compute the mass fraction of the last species.
            Y_last -= *(mass_fractions[si]);
        }
        
        /*
         * Add the contribution from the last species.
         */
        getSpeciesMolecularProperties(species_molecular_properties_ptr, d_num_species - 1);
        
        const double kappa_last = d_equation_of_thermal_conductivity->
            getThermalConductivity(
                pressure,
                temperature,
                species_molecular_properties_const_ptr);
        
        num += kappa_last*Y_last/(sqrt(species_molecular_properties[1]));
        den += Y_last/(sqrt(species_molecular_properties[1]));
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Number of mass fractions provided is not"
            << " equal to the total number of species or (total number of species - 1)."
            << std::endl);
    }
    
    kappa = num/den;
    
    return kappa;
}


/*
 * Compute the thermal conductivity of the mixture with isothermal and isobaric equilibria assumptions.
 */
void
EquationOfThermalConductivityMixingRulesConstant::computeThermalConductivity(
    boost::shared_ptr<pdat::CellData<double> >& data_thermal_conductivity,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
    const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
    const hier::Box& domain) const
{

}


/*
 * Get the molecular properties of a species.
 */
void
EquationOfThermalConductivityMixingRulesConstant::getSpeciesMolecularProperties(
    std::vector<double*>& species_molecular_properties,
    const int species_index) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(species_molecular_properties.size()) >= 2);
    TBOX_ASSERT(species_index >= 0);
    TBOX_ASSERT(species_index < d_num_species);
#endif
    
    *(species_molecular_properties[0]) = d_species_kappa[species_index];
    *(species_molecular_properties[1]) = d_species_M[species_index];
}

