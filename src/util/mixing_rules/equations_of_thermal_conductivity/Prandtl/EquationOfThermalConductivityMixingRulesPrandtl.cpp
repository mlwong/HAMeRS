#include "util/mixing_rules/equations_of_thermal_conductivity/Prandtl/EquationOfThermalConductivityMixingRulesPrandtl.hpp"

EquationOfThermalConductivityMixingRulesPrandtl::EquationOfThermalConductivityMixingRulesPrandtl(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const int& num_species,
    const MIXING_CLOSURE_MODEL& mixing_closure_model,
    const boost::shared_ptr<tbox::Database>& equation_of_thermal_conductivity_mixing_rules_db):
        EquationOfThermalConductivityMixingRules(
            object_name,
            dim,
            num_species,
            mixing_closure_model,
            equation_of_thermal_conductivity_mixing_rules_db)
{
    d_equation_of_thermal_conductivity.reset(new EquationOfThermalConductivityPrandtl(
        "d_equation_of_thermal_conductivity",
        dim));
    
    /*
     * Get the specific heat at constant pressure of each species from the database.
     */
    
    if (equation_of_thermal_conductivity_mixing_rules_db->keyExists("species_c_p"))
    {
        size_t species_c_p_array_size =
            equation_of_thermal_conductivity_mixing_rules_db->getArraySize("species_c_p");
        if (static_cast<int>(species_c_p_array_size) == d_num_species)
        {
            d_species_c_p =
                equation_of_thermal_conductivity_mixing_rules_db->getDoubleVector("species_c_p");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'species_c_p' entries must be equal to 'num_species'."
                << std::endl);
        }
    }
    else if (equation_of_thermal_conductivity_mixing_rules_db->keyExists("d_species_c_p"))
    {
        size_t species_c_p_array_size =
            equation_of_thermal_conductivity_mixing_rules_db->getArraySize("d_species_c_p");
        if (static_cast<int>(species_c_p_array_size) == d_num_species)
        {
            d_species_c_p =
                equation_of_thermal_conductivity_mixing_rules_db->getDoubleVector("d_species_c_p");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'd_species_c_p' entries must be equal to 'd_num_species'."
                << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Key data 'species_c_p'/'d_species_c_p'"
            << "not found in data for equation of thermal conductivity mixing rules."
            << std::endl);
    }
    
    /*
     * Get the Prandtl number of each species from the database.
     */
    
    if (equation_of_thermal_conductivity_mixing_rules_db->keyExists("species_Pr"))
    {
        size_t species_Pr_array_size =
            equation_of_thermal_conductivity_mixing_rules_db->getArraySize("species_Pr");
        if (static_cast<int>(species_Pr_array_size) == d_num_species)
        {
            d_species_Pr =
                equation_of_thermal_conductivity_mixing_rules_db->getDoubleVector("species_Pr");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'species_Pr' entries must be equal to 'num_species'."
                << std::endl);
        }
    }
    else if (equation_of_thermal_conductivity_mixing_rules_db->keyExists("d_species_Pr"))
    {
        size_t species_Pr_array_size =
            equation_of_thermal_conductivity_mixing_rules_db->getArraySize("d_species_Pr");
        if (static_cast<int>(species_Pr_array_size) == d_num_species)
        {
            d_species_Pr =
                equation_of_thermal_conductivity_mixing_rules_db->getDoubleVector("d_species_Pr");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'd_species_Pr' entries must be equal to 'num_species'."
                << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Key data 'species_Pr'/'d_species_Pr'"
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
    
    /*
     * Initialize d_equation_of_shear_viscosity_mixing_rules_manager and get the equation of shear viscosity
     * mixing rules object.
     */
    
    if (equation_of_thermal_conductivity_mixing_rules_db->keyExists("equation_of_shear_viscosity"))
    {
        d_equation_of_shear_viscosity_str =
            equation_of_thermal_conductivity_mixing_rules_db->getString("equation_of_shear_viscosity");
    }
    else if (equation_of_thermal_conductivity_mixing_rules_db->keyExists("d_equation_of_shear_viscosity_str"))
    {
        d_equation_of_shear_viscosity_str =
            equation_of_thermal_conductivity_mixing_rules_db->getString("d_equation_of_shear_viscosity_str");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "No key 'equation_of_shear_viscosity'/'d_equation_of_shear_viscosity_str' found in data"
            << " for equation of thermal conductivity mixing rules."
            << std::endl);
    }
    
    d_equation_of_shear_viscosity_mixing_rules_manager.reset(new EquationOfShearViscosityMixingRulesManager(
        "d_equation_of_shear_viscosity_mixing_rules_manager",
        d_dim,
        d_num_species,
        mixing_closure_model,
        equation_of_thermal_conductivity_mixing_rules_db,
        d_equation_of_shear_viscosity_str));
    
    d_equation_of_shear_viscosity_mixing_rules =
        d_equation_of_shear_viscosity_mixing_rules_manager->getEquationOfShearViscosityMixingRules();
}


/*
 * Print all characteristics of the equation of thermal conductivity class.
 */
void
EquationOfThermalConductivityMixingRulesPrandtl::printClassData(
    std::ostream& os) const
{
    os << "\nPrint EquationOfThermalConductivityMixingRulesPrandtl object..."
       << std::endl;
    
    os << std::endl;
    os << "EquationOfThermalConductivityMixingRulesPrandtl: this = "
       << (EquationOfThermalConductivityMixingRulesPrandtl *)this
       << std::endl;
    
    os << "d_object_name = "
       << d_object_name
       << std::endl;
    
    os << "d_mixing_closure_model = "
       << d_mixing_closure_model
       << std::endl;
    
    /*
     * Print the specific heats at constant pressure of each species.
     */
    
    os << "d_species_c_p = ";
    for (int si = 0; si < d_num_species - 1; si++)
    {
        os << d_species_c_p[si] << ", ";
    }
    os << d_species_c_p[d_num_species - 1];
    os << std::endl;
    
    /*
     * Print the Prandtl number of each species.
     */
    
    os << "d_species_Pr = ";
    for (int si = 0; si < d_num_species - 1; si++)
    {
        os << d_species_Pr[si] << ", ";
    }
    os << d_species_Pr[d_num_species - 1];
    os << std::endl;
    
    os << "................................................................................";
    
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
    
    os << "................................................................................";
    
    d_equation_of_shear_viscosity_mixing_rules->printClassData(os);
}


/*
 * Put the characteristics of the equation of thermal conductivity mixing rules class into the restart
 * database.
 */
void
EquationOfThermalConductivityMixingRulesPrandtl::putToRestart(
    const boost::shared_ptr<tbox::Database>& restart_db) const
{
    restart_db->putDoubleVector("d_species_c_p", d_species_c_p);
    restart_db->putDoubleVector("d_species_Pr", d_species_Pr);
    restart_db->putDoubleVector("d_species_M", d_species_M);
    
    restart_db->putString("d_equation_of_shear_viscosity_str", d_equation_of_shear_viscosity_str);
    d_equation_of_shear_viscosity_mixing_rules->putToRestart(restart_db);
}


/*
 * Compute the thermal conductivity of the mixture with isothermal and isobaric assumptions.
 */
double
EquationOfThermalConductivityMixingRulesPrandtl::getThermalConductivity(
    const double* const pressure,
    const double* const temperature,
    const std::vector<const double*>& mass_fraction) const
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT((d_mixing_closure_model == ISOTHERMAL_AND_ISOBARIC) ||
                (d_mixing_closure_model == NO_MODEL && d_num_species == 1));
    TBOX_ASSERT((static_cast<int>(mass_fraction.size()) == d_num_species) ||
                (static_cast<int>(mass_fraction.size()) == d_num_species - 1));
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
    
    const int num_molecular_properties = getNumberOfSpeciesMolecularProperties() + 1;
    
    species_molecular_properties.resize(num_molecular_properties);
    species_molecular_properties_ptr.reserve(num_molecular_properties);
    species_molecular_properties_const_ptr.reserve(num_molecular_properties);
    
    for (int mi = 0; mi < num_molecular_properties; mi++)
    {
        species_molecular_properties_ptr.push_back(&species_molecular_properties[mi]);
        species_molecular_properties_const_ptr.push_back(&species_molecular_properties[mi]);
    }
    
    std::vector<double> mu_molecular_properties;
    std::vector<double*> mu_molecular_properties_ptr;
    std::vector<const double*> mu_molecular_properties_const_ptr;
    
    const int num_mu_molecular_properties = d_equation_of_shear_viscosity_mixing_rules->
        getNumberOfSpeciesMolecularProperties();
    
    mu_molecular_properties.resize(num_mu_molecular_properties);
    mu_molecular_properties_ptr.reserve(num_mu_molecular_properties);
    mu_molecular_properties_const_ptr.reserve(num_mu_molecular_properties);
    
    for (int mi = 0; mi < num_mu_molecular_properties; mi++)
    {
        mu_molecular_properties_ptr.push_back(&mu_molecular_properties[mi]);
        mu_molecular_properties_const_ptr.push_back(&mu_molecular_properties[mi]);
    }
    
    if (static_cast<int>(mass_fraction.size()) == d_num_species - 1)
    {
        double Y_last = 1.0;
        
        for (int si = 0; si < d_num_species - 1; si++)
        {
            getSpeciesMolecularProperties(species_molecular_properties_ptr, si);
            
            d_equation_of_shear_viscosity_mixing_rules->
                getSpeciesMolecularProperties(mu_molecular_properties_ptr, si);
            
            species_molecular_properties[num_molecular_properties - 1] =
                d_equation_of_shear_viscosity_mixing_rules->
                    d_equation_of_shear_viscosity->
                        getShearViscosity(
                            pressure,
                            temperature,
                            mu_molecular_properties_const_ptr);
            
            const double kappa_i = d_equation_of_thermal_conductivity->
                getThermalConductivity(
                    pressure,
                    temperature,
                    species_molecular_properties_const_ptr);
            
            num += kappa_i*(*(mass_fraction[si]))/(sqrt(species_molecular_properties[2]));
            den += *(mass_fraction[si])/(sqrt(species_molecular_properties[2]));
            
            // Compute the mass fraction of the last species.
            Y_last -= *(mass_fraction[si]);
        }
        
        /*
         * Add the contribution from the last species.
         */
        
        getSpeciesMolecularProperties(species_molecular_properties_ptr, d_num_species - 1);
        
        d_equation_of_shear_viscosity_mixing_rules->
            getSpeciesMolecularProperties(mu_molecular_properties_ptr, d_num_species - 1);
        
        species_molecular_properties[num_molecular_properties - 1] =
            d_equation_of_shear_viscosity_mixing_rules->
                d_equation_of_shear_viscosity->
                    getShearViscosity(
                        pressure,
                        temperature,
                        mu_molecular_properties_const_ptr);
        
        const double kappa_last = d_equation_of_thermal_conductivity->
            getThermalConductivity(
                pressure,
                temperature,
                species_molecular_properties_const_ptr);
        
        num += kappa_last*Y_last/(sqrt(species_molecular_properties[2]));
        den += Y_last/(sqrt(species_molecular_properties[2]));
    }
    else
    {
        for (int si = 0; si < d_num_species; si++)
        {
            getSpeciesMolecularProperties(species_molecular_properties_ptr, si);
            
            d_equation_of_shear_viscosity_mixing_rules->
                getSpeciesMolecularProperties(mu_molecular_properties_ptr, si);
            
            species_molecular_properties[num_molecular_properties - 1] =
                d_equation_of_shear_viscosity_mixing_rules->
                    d_equation_of_shear_viscosity->
                        getShearViscosity(
                            pressure,
                            temperature,
                            mu_molecular_properties_const_ptr);
            
            const double kappa_i = d_equation_of_thermal_conductivity->
                getThermalConductivity(
                    pressure,
                    temperature,
                    species_molecular_properties_const_ptr);
            
            num += kappa_i*(*(mass_fraction[si]))/(sqrt(species_molecular_properties[2]));
            den += *(mass_fraction[si])/(sqrt(species_molecular_properties[2]));
        }
    }
    
    kappa = num/den;
    
    return kappa;
}


/*
 * Compute the thermal conductivity of the mixture with isobaric assumption.
 */
double
EquationOfThermalConductivityMixingRulesPrandtl::getThermalConductivity(
    const double* const pressure,
    const std::vector<const double*>& temperature,
    const std::vector<const double*>& mass_fraction,
    const std::vector<const double*>& volume_fraction) const
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(d_mixing_closure_model == ISOBARIC);
    TBOX_ASSERT((static_cast<int>(temperature.size()) == d_num_species);
    TBOX_ASSERT((static_cast<int>(volume_fraction.size()) == d_num_species) ||
                (static_cast<int>(volume_fraction.size()) == d_num_species - 1));
#endif
    
    NULL_USE(mass_fraction);
    
    double kappa = 0.0;
    
    /*
     * Initialize the container and pointers to the container for the molecular properties
     * of a species.
     */
    
    std::vector<double> species_molecular_properties;
    std::vector<double*> species_molecular_properties_ptr;
    std::vector<const double*> species_molecular_properties_const_ptr;
    
    const int num_molecular_properties = getNumberOfSpeciesMolecularProperties() + 1;
    
    species_molecular_properties.resize(num_molecular_properties);
    species_molecular_properties_ptr.reserve(num_molecular_properties);
    species_molecular_properties_const_ptr.reserve(num_molecular_properties);
    
    for (int mi = 0; mi < num_molecular_properties; mi++)
    {
        species_molecular_properties_ptr.push_back(&species_molecular_properties[mi]);
        species_molecular_properties_const_ptr.push_back(&species_molecular_properties[mi]);
    }
    
    std::vector<double> mu_molecular_properties;
    std::vector<double*> mu_molecular_properties_ptr;
    std::vector<const double*> mu_molecular_properties_const_ptr;
    
    const int num_mu_molecular_properties = d_equation_of_shear_viscosity_mixing_rules->
        getNumberOfSpeciesMolecularProperties();
    
    mu_molecular_properties.resize(num_mu_molecular_properties);
    mu_molecular_properties_ptr.reserve(num_mu_molecular_properties);
    mu_molecular_properties_const_ptr.reserve(num_mu_molecular_properties);
    
    for (int mi = 0; mi < num_mu_molecular_properties; mi++)
    {
        mu_molecular_properties_ptr.push_back(&mu_molecular_properties[mi]);
        mu_molecular_properties_const_ptr.push_back(&mu_molecular_properties[mi]);
    }
    
    if (static_cast<int>(mass_fraction.size()) == d_num_species - 1)
    {
        double Z_last = 1.0;
        
        for (int si = 0; si < d_num_species - 1; si++)
        {
            getSpeciesMolecularProperties(species_molecular_properties_ptr, si);
            
            d_equation_of_shear_viscosity_mixing_rules->
                getSpeciesMolecularProperties(mu_molecular_properties_ptr, si);
            
            species_molecular_properties[num_molecular_properties - 1] =
                d_equation_of_shear_viscosity_mixing_rules->
                    d_equation_of_shear_viscosity->
                        getShearViscosity(
                            pressure,
                            temperature[si],
                            mu_molecular_properties_const_ptr);
            
            const double kappa_i = d_equation_of_thermal_conductivity->
                getThermalConductivity(
                    pressure,
                    temperature[si],
                    species_molecular_properties_const_ptr);
            
            kappa += *(volume_fraction[si])*kappa_i;
            
            // Compute the volume fraction of the last species.
            Z_last -= *(volume_fraction[si]);
        }
        
        /*
         * Add the contribution from the last species.
         */
        
        getSpeciesMolecularProperties(species_molecular_properties_ptr, d_num_species - 1);
        
        d_equation_of_shear_viscosity_mixing_rules->
            getSpeciesMolecularProperties(mu_molecular_properties_ptr, d_num_species - 1);
        
        species_molecular_properties[num_molecular_properties - 1] =
            d_equation_of_shear_viscosity_mixing_rules->
                d_equation_of_shear_viscosity->
                    getShearViscosity(
                        pressure,
                        temperature[d_num_species - 1],
                        mu_molecular_properties_const_ptr);
        
        const double kappa_last = d_equation_of_thermal_conductivity->
            getThermalConductivity(
                pressure,
                temperature[d_num_species - 1],
                species_molecular_properties_const_ptr);
        
        kappa += Z_last*kappa_last;
    }
    else
    {
        for (int si = 0; si < d_num_species; si++)
        {
            getSpeciesMolecularProperties(species_molecular_properties_ptr, si);
            
            d_equation_of_shear_viscosity_mixing_rules->
                getSpeciesMolecularProperties(mu_molecular_properties_ptr, si);
            
            species_molecular_properties[num_molecular_properties - 1] =
                d_equation_of_shear_viscosity_mixing_rules->
                    d_equation_of_shear_viscosity->
                        getShearViscosity(
                            pressure,
                            temperature[si],
                            mu_molecular_properties_const_ptr);
            
            const double kappa_i = d_equation_of_thermal_conductivity->
                getThermalConductivity(
                    pressure,
                    temperature[si],
                    species_molecular_properties_const_ptr);
            
            kappa += *(volume_fraction[si])*kappa_i;
        }
    }
    
    return kappa;
}


/*
 * Get the molecular properties of a species.
 */
void
EquationOfThermalConductivityMixingRulesPrandtl::getSpeciesMolecularProperties(
    std::vector<double*>& species_molecular_properties,
    const int& species_index) const
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(species_molecular_properties.size()) >= 3);
    TBOX_ASSERT(species_index >= 0);
    TBOX_ASSERT(species_index < d_num_species);
#endif
    
    *(species_molecular_properties[0]) = d_species_c_p[species_index];
    *(species_molecular_properties[1]) = d_species_Pr[species_index];
    *(species_molecular_properties[2]) = d_species_M[species_index];
}
