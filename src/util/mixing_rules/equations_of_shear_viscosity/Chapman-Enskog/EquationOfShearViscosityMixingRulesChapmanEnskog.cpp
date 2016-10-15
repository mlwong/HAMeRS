#include "util/mixing_rules/equations_of_shear_viscosity/Chapman-Enskog/EquationOfShearViscosityMixingRulesChapmanEnskog.hpp"

EquationOfShearViscosityMixingRulesChapmanEnskog::EquationOfShearViscosityMixingRulesChapmanEnskog(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const int& num_species,
    const MIXING_CLOSURE_MODEL& mixing_closure_model,
    const boost::shared_ptr<tbox::Database>& equation_of_shear_viscosity_mixing_rules_db):
        EquationOfShearViscosityMixingRules(
            object_name,
            dim,
            num_species,
            mixing_closure_model,
            equation_of_shear_viscosity_mixing_rules_db)
{
    d_equation_of_shear_viscosity.reset(new EquationOfShearViscosityChapmanEnskog(
        "d_equation_of_shear_viscosity",
        dim));
    
    /*
     * Get the Lennard-Jones energy parameter of each species from the database.
     */
    
    if (equation_of_shear_viscosity_mixing_rules_db->keyExists("species_epsilon_by_k"))
    {
        size_t species_epsilon_by_k_array_size =
            equation_of_shear_viscosity_mixing_rules_db->getArraySize("species_epsilon_by_k");
        if (static_cast<int>(species_epsilon_by_k_array_size) == d_num_species)
        {
            d_species_epsilon_by_k =
                equation_of_shear_viscosity_mixing_rules_db->getDoubleVector("species_epsilon_by_k");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'species_epsilon_by_k' entries must be equal to 'num_species'."
                << std::endl);
        }
    }
    else if (equation_of_shear_viscosity_mixing_rules_db->keyExists("d_species_epsilon_by_k"))
    {
        size_t species_epsilon_by_k_array_size =
            equation_of_shear_viscosity_mixing_rules_db->getArraySize("d_species_epsilon_by_k");
        if (static_cast<int>(species_epsilon_by_k_array_size) == d_num_species)
        {
            d_species_epsilon_by_k =
                equation_of_shear_viscosity_mixing_rules_db->getDoubleVector("d_species_epsilon_by_k");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'd_species_epsilon_by_k' entries must be equal to 'd_num_species'."
                << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Key data 'species_epsilon_by_k'/'d_species_epsilon_by_k'"
            << "not found in data for equation of shear viscosity mixing rules."
            << std::endl);
    }
    
    /*
     * Get the collision diameter of each species from the database.
     */
    
    if (equation_of_shear_viscosity_mixing_rules_db->keyExists("species_sigma"))
    {
        size_t species_sigma_array_size =
            equation_of_shear_viscosity_mixing_rules_db->getArraySize("species_sigma");
        if (static_cast<int>(species_sigma_array_size) == d_num_species)
        {
            d_species_sigma =
                equation_of_shear_viscosity_mixing_rules_db->getDoubleVector("species_sigma");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'species_sigma' entries must be equal to 'num_species'."
                << std::endl);
        }
    }
    else if (equation_of_shear_viscosity_mixing_rules_db->keyExists("d_species_sigma"))
    {
        size_t species_sigma_array_size =
            equation_of_shear_viscosity_mixing_rules_db->getArraySize("d_species_sigma");
        if (static_cast<int>(species_sigma_array_size) == d_num_species)
        {
            d_species_sigma =
                equation_of_shear_viscosity_mixing_rules_db->getDoubleVector("d_species_sigma");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'd_species_sigma' entries must be equal to 'num_species'."
                << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Key data 'species_sigma'/'d_species_sigma'"
            << "not found in data for equation of shear viscosity mixing rules."
            << std::endl);
    }
    
    /*
     * Get the molecular weight of each species from the database.
     */
    
    if (equation_of_shear_viscosity_mixing_rules_db->keyExists("species_M"))
    {
        size_t species_M_array_size =
            equation_of_shear_viscosity_mixing_rules_db->getArraySize("species_M");
        if (static_cast<int>(species_M_array_size) == d_num_species)
        {
            d_species_M =
                equation_of_shear_viscosity_mixing_rules_db->getDoubleVector("species_M");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'species_M' entries must be equal to 'num_species'."
                << std::endl);
        }
    }
    else if (equation_of_shear_viscosity_mixing_rules_db->keyExists("d_species_M"))
    {
        size_t species_M_array_size =
            equation_of_shear_viscosity_mixing_rules_db->getArraySize("d_species_M");
        if (static_cast<int>(species_M_array_size) == d_num_species)
        {
            d_species_M =
                equation_of_shear_viscosity_mixing_rules_db->getDoubleVector("d_species_M");
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
            << "not found in data for equation of shear viscosity mixing rules."
            << std::endl);
    }
}


/*
 * Print all characteristics of the equation of shear viscosity class.
 */
void
EquationOfShearViscosityMixingRulesChapmanEnskog::printClassData(
    std::ostream& os) const
{
    os << "\nPrint EquationOfShearViscosityMixingRulesChapmanEnskog object..."
       << std::endl;
    
    os << std::endl;
    os << "EquationOfShearViscosityMixingRulesChapmanEnskog: this = "
       << (EquationOfShearViscosityMixingRulesChapmanEnskog *)this
       << std::endl;
    
    os << "d_object_name = "
       << d_object_name
       << std::endl;
    
    os << "d_mixing_closure_model = "
       << d_mixing_closure_model
       << std::endl;
    
    /*
     * Print the Lennard-Jones energy parameter of each species.
     */
    
    os << "d_species_epsilon_by_k = ";
    for (int si = 0; si < d_num_species - 1; si++)
    {
        os << d_species_epsilon_by_k[si] << ", ";
    }
    os << d_species_epsilon_by_k[d_num_species - 1];
    os << std::endl;
    
    /*
     * Print the collision diameter of each species.
     */
    
    os << "d_species_sigma = ";
    for (int si = 0; si < d_num_species - 1; si++)
    {
        os << d_species_sigma[si] << ", ";
    }
    os << d_species_sigma[d_num_species - 1];
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
 * Put the characteristics of the equation of shear viscosity mixing rules class into the restart
 * database.
 */
void
EquationOfShearViscosityMixingRulesChapmanEnskog::putToRestart(
    const boost::shared_ptr<tbox::Database>& restart_db) const
{
    restart_db->putDoubleVector("d_species_epsilon_by_k", d_species_epsilon_by_k);
    restart_db->putDoubleVector("d_species_sigma", d_species_sigma);
    restart_db->putDoubleVector("d_species_M", d_species_M);
}


/*
 * Compute the shear viscosity of the mixture with isothermal and isobaric assumptions.
 */
double
EquationOfShearViscosityMixingRulesChapmanEnskog::getShearViscosity(
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
    
    double mu = 0.0;
    
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
    
    if (static_cast<int>(mass_fraction.size()) == d_num_species - 1)
    {
        double Y_last = 1.0;
        
        for (int si = 0; si < d_num_species - 1; si++)
        {
            getSpeciesMolecularProperties(species_molecular_properties_ptr, si);
            
            const double mu_i = d_equation_of_shear_viscosity->
                getShearViscosity(
                    pressure,
                    temperature,
                    species_molecular_properties_const_ptr);
            
            num += mu_i*(*(mass_fraction[si]))/(sqrt(species_molecular_properties[2]));
            den += *(mass_fraction[si])/(sqrt(species_molecular_properties[2]));
            
            // Compute the mass fraction of the last species.
            Y_last -= *(mass_fraction[si]);
        }
        
        /*
         * Add the contribution from the last species.
         */
        
        getSpeciesMolecularProperties(species_molecular_properties_ptr, d_num_species - 1);
            
        const double mu_last = d_equation_of_shear_viscosity->
            getShearViscosity(
                pressure,
                temperature,
                species_molecular_properties_const_ptr);
        
        num += mu_last*Y_last/(sqrt(species_molecular_properties[2]));
        den += Y_last/(sqrt(species_molecular_properties[2]));
    }
    else
    {
        for (int si = 0; si < d_num_species; si++)
        {
            getSpeciesMolecularProperties(species_molecular_properties_ptr, si);
            
            const double mu_i = d_equation_of_shear_viscosity->
                getShearViscosity(
                    pressure,
                    temperature,
                    species_molecular_properties_const_ptr);
            
            num += mu_i*(*(mass_fraction[si]))/(sqrt(species_molecular_properties[2]));
            den += *(mass_fraction[si])/(sqrt(species_molecular_properties[2]));
        }
    }
    
    mu = num/den;
    
    return mu;
}


/*
 * Compute the shear viscosity of the mixture with isobaric assumption.
 */
double
EquationOfShearViscosityMixingRulesChapmanEnskog::getShearViscosity(
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
    
    double mu = 0.0;
    
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
    
    if (static_cast<int>(volume_fraction.size()) == d_num_species - 1)
    {
        double Z_last = 1.0;
        
        for (int si = 0; si < d_num_species - 1; si++)
        {
            getSpeciesMolecularProperties(species_molecular_properties_ptr, si);
            
            const double mu_i = d_equation_of_shear_viscosity->
                getShearViscosity(
                    pressure,
                    temperature[si],
                    species_molecular_properties_const_ptr);
            
            mu += *(volume_fraction[si])*mu_i;
            
            // Compute the volume fraction of the last species.
            Z_last -= *(volume_fraction[si]);
        }
        
        /*
         * Add the contribution from the last species.
         */
        
        getSpeciesMolecularProperties(species_molecular_properties_ptr, d_num_species - 1);
            
        const double mu_last = d_equation_of_shear_viscosity->
            getShearViscosity(
                pressure,
                temperature[d_num_species - 1],
                species_molecular_properties_const_ptr);
        
        mu += Z_last*mu_last;
    }
    else
    {
        for (int si = 0; si < d_num_species; si++)
        {
            getSpeciesMolecularProperties(species_molecular_properties_ptr, si);
            
            const double mu_i = d_equation_of_shear_viscosity->
                getShearViscosity(
                    pressure,
                    temperature[si],
                    species_molecular_properties_const_ptr);
            
            mu += *(volume_fraction[si])*mu_i;
        }
    }
    
    return mu;
}


/*
 * Get the molecular properties of a species.
 */
void
EquationOfShearViscosityMixingRulesChapmanEnskog::getSpeciesMolecularProperties(
    std::vector<double*>& species_molecular_properties,
    const int& species_index) const
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(species_molecular_properties.size()) >= 3);
    TBOX_ASSERT(species_index >= 0);
    TBOX_ASSERT(species_index < d_num_species);
#endif
    
    *(species_molecular_properties[0]) = d_species_epsilon_by_k[species_index];
    *(species_molecular_properties[1]) = d_species_sigma[species_index];
    *(species_molecular_properties[2]) = d_species_M[species_index];
}
