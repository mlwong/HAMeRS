#include "util/mixing_rules/equations_of_bulk_viscosity/Cramer/EquationOfBulkViscosityMixingRulesCramer.hpp"

EquationOfBulkViscosityMixingRulesCramer::EquationOfBulkViscosityMixingRulesCramer(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const int& num_species,
    const MIXING_CLOSURE_MODEL::TYPE& mixing_closure_model,
    const boost::shared_ptr<tbox::Database>& equation_of_bulk_viscosity_mixing_rules_db):
        EquationOfBulkViscosityMixingRules(
            object_name,
            dim,
            num_species,
            mixing_closure_model,
            equation_of_bulk_viscosity_mixing_rules_db)
{
    d_equation_of_bulk_viscosity.reset(new EquationOfBulkViscosityCramer(
        "d_equation_of_bulk_viscosity",
        dim));
    
    /*
     * Get the ratio of specific heats of each species from the database.
     */
    
    if (equation_of_bulk_viscosity_mixing_rules_db->keyExists("species_gamma"))
    {
        size_t species_gamma_array_size =
            equation_of_bulk_viscosity_mixing_rules_db->getArraySize("species_gamma");
        if (static_cast<int>(species_gamma_array_size) == d_num_species)
        {
            d_species_gamma =
                equation_of_bulk_viscosity_mixing_rules_db->getDoubleVector("species_gamma");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'species_gamma' entries must be equal to 'num_species'."
                << std::endl);
        }
    }
    else if (equation_of_bulk_viscosity_mixing_rules_db->keyExists("d_species_gamma"))
    {
        size_t species_gamma_array_size =
            equation_of_bulk_viscosity_mixing_rules_db->getArraySize("d_species_gamma");
        if (static_cast<int>(species_gamma_array_size) == d_num_species)
        {
            d_species_gamma =
                equation_of_bulk_viscosity_mixing_rules_db->getDoubleVector("d_species_gamma");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'd_species_gamma' entries must be equal to 'd_num_species'."
                << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Key data 'species_gamma'/'d_species_gamma'"
            << "not found in data for equation of bulk viscosity mixing rules."
            << std::endl);
    }
    
    /*
     * Get parameters for rotational mode of each species from the database.
     */
    
    if (equation_of_bulk_viscosity_mixing_rules_db->keyExists("species_A_r"))
    {
        size_t species_A_r_array_size =
            equation_of_bulk_viscosity_mixing_rules_db->getArraySize("species_A_r");
        if (static_cast<int>(species_A_r_array_size) == d_num_species)
        {
            d_species_A_r =
                equation_of_bulk_viscosity_mixing_rules_db->getDoubleVector("species_A_r");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'species_A_r' entries must be equal to 'num_species'."
                << std::endl);
        }
    }
    else if (equation_of_bulk_viscosity_mixing_rules_db->keyExists("d_species_A_r"))
    {
        size_t species_A_r_array_size =
            equation_of_bulk_viscosity_mixing_rules_db->getArraySize("d_species_A_r");
        if (static_cast<int>(species_A_r_array_size) == d_num_species)
        {
            d_species_A_r =
                equation_of_bulk_viscosity_mixing_rules_db->getDoubleVector("d_species_A_r");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'd_species_A_r' entries must be equal to 'num_species'."
                << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Key data 'species_A_r'/'d_species_A_r'"
            << "not found in data for equation of bulk viscosity mixing rules."
            << std::endl);
    }
    
    if (equation_of_bulk_viscosity_mixing_rules_db->keyExists("species_B_r"))
    {
        size_t species_B_r_array_size =
            equation_of_bulk_viscosity_mixing_rules_db->getArraySize("species_B_r");
        if (static_cast<int>(species_B_r_array_size) == d_num_species)
        {
            d_species_B_r =
                equation_of_bulk_viscosity_mixing_rules_db->getDoubleVector("species_B_r");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'species_B_r' entries must be equal to 'num_species'."
                << std::endl);
        }
    }
    else if (equation_of_bulk_viscosity_mixing_rules_db->keyExists("d_species_B_r"))
    {
        size_t species_B_r_array_size =
            equation_of_bulk_viscosity_mixing_rules_db->getArraySize("d_species_B_r");
        if (static_cast<int>(species_B_r_array_size) == d_num_species)
        {
            d_species_B_r =
                equation_of_bulk_viscosity_mixing_rules_db->getDoubleVector("d_species_B_r");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'd_species_B_r' entries must be equal to 'num_species'."
                << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Key data 'species_B_r'/'d_species_B_r'"
            << "not found in data for equation of bulk viscosity mixing rules."
            << std::endl);
    }
    
    /*
     * Get parameters for vibrational mode of each species from the database.
     */
    
    if (equation_of_bulk_viscosity_mixing_rules_db->keyExists("species_c_v_v"))
    {
        size_t species_c_v_v_array_size =
            equation_of_bulk_viscosity_mixing_rules_db->getArraySize("species_c_v_v");
        if (static_cast<int>(species_c_v_v_array_size) == d_num_species)
        {
            d_species_c_v_v =
                equation_of_bulk_viscosity_mixing_rules_db->getDoubleVector("species_c_v_v");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'species_c_v_v' entries must be equal to 'num_species'."
                << std::endl);
        }
    }
    else if (equation_of_bulk_viscosity_mixing_rules_db->keyExists("d_species_c_v_v"))
    {
        size_t species_c_v_v_array_size =
            equation_of_bulk_viscosity_mixing_rules_db->getArraySize("d_species_c_v_v");
        if (static_cast<int>(species_c_v_v_array_size) == d_num_species)
        {
            d_species_c_v_v =
                equation_of_bulk_viscosity_mixing_rules_db->getDoubleVector("d_species_c_v_v");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'd_species_c_v_v' entries must be equal to 'num_species'."
                << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Key data 'species_c_v_v'/'d_species_c_v_v'"
            << "not found in data for equation of bulk viscosity mixing rules."
            << std::endl);
    }
    
    if (equation_of_bulk_viscosity_mixing_rules_db->keyExists("species_A_v"))
    {
        size_t species_A_v_array_size =
            equation_of_bulk_viscosity_mixing_rules_db->getArraySize("species_A_v");
        if (static_cast<int>(species_A_v_array_size) == d_num_species)
        {
            d_species_A_v =
                equation_of_bulk_viscosity_mixing_rules_db->getDoubleVector("species_A_v");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'species_A_v' entries must be equal to 'num_species'."
                << std::endl);
        }
    }
    else if (equation_of_bulk_viscosity_mixing_rules_db->keyExists("d_species_A_v"))
    {
        size_t species_A_v_array_size =
            equation_of_bulk_viscosity_mixing_rules_db->getArraySize("d_species_A_v");
        if (static_cast<int>(species_A_v_array_size) == d_num_species)
        {
            d_species_A_v =
                equation_of_bulk_viscosity_mixing_rules_db->getDoubleVector("d_species_A_v");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'd_species_A_v' entries must be equal to 'num_species'."
                << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Key data 'species_A_v'/'d_species_A_v'"
            << "not found in data for equation of bulk viscosity mixing rules."
            << std::endl);
    }
    
    if (equation_of_bulk_viscosity_mixing_rules_db->keyExists("species_B_v"))
    {
        size_t species_B_v_array_size =
            equation_of_bulk_viscosity_mixing_rules_db->getArraySize("species_B_v");
        if (static_cast<int>(species_B_v_array_size) == d_num_species)
        {
            d_species_B_v =
                equation_of_bulk_viscosity_mixing_rules_db->getDoubleVector("species_B_v");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'species_B_v' entries must be equal to 'num_species'."
                << std::endl);
        }
    }
    else if (equation_of_bulk_viscosity_mixing_rules_db->keyExists("d_species_B_v"))
    {
        size_t species_B_v_array_size =
            equation_of_bulk_viscosity_mixing_rules_db->getArraySize("d_species_B_v");
        if (static_cast<int>(species_B_v_array_size) == d_num_species)
        {
            d_species_B_v =
                equation_of_bulk_viscosity_mixing_rules_db->getDoubleVector("d_species_B_v");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'd_species_B_v' entries must be equal to 'num_species'."
                << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Key data 'species_B_v'/'d_species_B_v'"
            << "not found in data for equation of bulk viscosity mixing rules."
            << std::endl);
    }
    
    if (equation_of_bulk_viscosity_mixing_rules_db->keyExists("species_C_v"))
    {
        size_t species_C_v_array_size =
            equation_of_bulk_viscosity_mixing_rules_db->getArraySize("species_C_v");
        if (static_cast<int>(species_C_v_array_size) == d_num_species)
        {
            d_species_C_v =
                equation_of_bulk_viscosity_mixing_rules_db->getDoubleVector("species_C_v");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'species_C_v' entries must be equal to 'num_species'."
                << std::endl);
        }
    }
    else if (equation_of_bulk_viscosity_mixing_rules_db->keyExists("d_species_C_v"))
    {
        size_t species_C_v_array_size =
            equation_of_bulk_viscosity_mixing_rules_db->getArraySize("d_species_C_v");
        if (static_cast<int>(species_C_v_array_size) == d_num_species)
        {
            d_species_C_v =
                equation_of_bulk_viscosity_mixing_rules_db->getDoubleVector("d_species_C_v");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'd_species_C_v' entries must be equal to 'num_species'."
                << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Key data 'species_C_v'/'d_species_C_v'"
            << "not found in data for equation of bulk viscosity mixing rules."
            << std::endl);
    }
    
    /*
     * Get the molecular weight of each species from the database.
     */
    
    if (equation_of_bulk_viscosity_mixing_rules_db->keyExists("species_M"))
    {
        size_t species_M_array_size =
            equation_of_bulk_viscosity_mixing_rules_db->getArraySize("species_M");
        if (static_cast<int>(species_M_array_size) == d_num_species)
        {
            d_species_M =
                equation_of_bulk_viscosity_mixing_rules_db->getDoubleVector("species_M");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'species_M' entries must be equal to 'num_species'."
                << std::endl);
        }
    }
    else if (equation_of_bulk_viscosity_mixing_rules_db->keyExists("d_species_M"))
    {
        size_t species_M_array_size =
            equation_of_bulk_viscosity_mixing_rules_db->getArraySize("d_species_M");
        if (static_cast<int>(species_M_array_size) == d_num_species)
        {
            d_species_M =
                equation_of_bulk_viscosity_mixing_rules_db->getDoubleVector("d_species_M");
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
            << "not found in data for equation of bulk viscosity mixing rules."
            << std::endl);
    }
}


/*
 * Print all characteristics of the equation of bulk viscosity class.
 */
void
EquationOfBulkViscosityMixingRulesCramer::printClassData(
    std::ostream& os) const
{
    os << "\nPrint EquationOfBulkViscosityMixingRulesCramer object..."
       << std::endl;
    
    os << std::endl;
    os << "EquationOfBulkViscosityMixingRulesCramer: this = "
       << (EquationOfBulkViscosityMixingRulesCramer *)this
       << std::endl;
    
    os << "d_object_name = "
       << d_object_name
       << std::endl;
    
    os << "d_mixing_closure_model = "
       << d_mixing_closure_model
       << std::endl;
    
    /*
     * Print the ratio of specific heats of each species.
     */
    
    os << "d_species_gamma = ";
    for (int si = 0; si < d_num_species - 1; si++)
    {
        os << d_species_gamma[si] << ", ";
    }
    os << d_species_gamma[d_num_species - 1];
    os << std::endl;
    
    /*
     * Print the parameters for rotational mode of each species.
     */
    
    os << "d_species_A_r = ";
    for (int si = 0; si < d_num_species - 1; si++)
    {
        os << d_species_A_r[si] << ", ";
    }
    os << d_species_A_r[d_num_species - 1];
    os << std::endl;
    
    os << "d_species_B_r = ";
    for (int si = 0; si < d_num_species - 1; si++)
    {
        os << d_species_B_r[si] << ", ";
    }
    os << d_species_B_r[d_num_species - 1];
    os << std::endl;
    
    /*
     * Print the parameters for vibrational mode of each species.
     */
    
    os << "d_species_c_v_v = ";
    for (int si = 0; si < d_num_species - 1; si++)
    {
        os << d_species_c_v_v[si] << ", ";
    }
    os << d_species_c_v_v[d_num_species - 1];
    os << std::endl;
    
    os << "d_species_A_v = ";
    for (int si = 0; si < d_num_species - 1; si++)
    {
        os << d_species_A_v[si] << ", ";
    }
    os << d_species_A_v[d_num_species - 1];
    os << std::endl;
    
    os << "d_species_B_v = ";
    for (int si = 0; si < d_num_species - 1; si++)
    {
        os << d_species_B_v[si] << ", ";
    }
    os << d_species_B_v[d_num_species - 1];
    os << std::endl;
    
    os << "d_species_C_v = ";
    for (int si = 0; si < d_num_species - 1; si++)
    {
        os << d_species_C_v[si] << ", ";
    }
    os << d_species_C_v[d_num_species - 1];
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
 * Put the characteristics of the equation of bulk viscosity mixing rules class into the restart
 * database.
 */
void
EquationOfBulkViscosityMixingRulesCramer::putToRestart(
    const boost::shared_ptr<tbox::Database>& restart_db) const
{
    restart_db->putDoubleVector("d_species_gamma", d_species_gamma);
    
    restart_db->putDoubleVector("d_species_A_r", d_species_A_r);
    restart_db->putDoubleVector("d_species_B_r", d_species_B_r);
    
    restart_db->putDoubleVector("d_species_c_v_v", d_species_c_v_v);
    restart_db->putDoubleVector("d_species_A_v", d_species_A_v);
    restart_db->putDoubleVector("d_species_B_v", d_species_B_v);
    restart_db->putDoubleVector("d_species_C_v", d_species_C_v);
    
    restart_db->putDoubleVector("d_species_M", d_species_M);
}


/*
 * Compute the bulk viscosity of the mixture with isothermal and isobaric assumptions.
 */
double
EquationOfBulkViscosityMixingRulesCramer::getBulkViscosity(
    const double* const pressure,
    const double* const temperature,
    const std::vector<const double*>& mass_fraction) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT((d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC) ||
                (d_mixing_closure_model == MIXING_CLOSURE_MODEL::NO_MODEL && d_num_species == 1));
    TBOX_ASSERT((static_cast<int>(mass_fraction.size()) == d_num_species) ||
                (static_cast<int>(mass_fraction.size()) == d_num_species - 1));
#endif
    
    double mu_v = 0.0;
    
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
            
            const double mu_v_i = d_equation_of_bulk_viscosity->
                getBulkViscosity(
                    pressure,
                    temperature,
                    species_molecular_properties_const_ptr);
            
            num += mu_v_i*(*(mass_fraction[si]))/(sqrt(species_molecular_properties[7]));
            den += *(mass_fraction[si])/(sqrt(species_molecular_properties[7]));
            
            // Compute the mass fraction of the last species.
            Y_last -= *(mass_fraction[si]);
        }
        
        /*
         * Add the contribution from the last species.
         */
        
        getSpeciesMolecularProperties(species_molecular_properties_ptr, d_num_species - 1);
            
        const double mu_v_last = d_equation_of_bulk_viscosity->
            getBulkViscosity(
                pressure,
                temperature,
                species_molecular_properties_const_ptr);
        
        num += mu_v_last*Y_last/(sqrt(species_molecular_properties[7]));
        den += Y_last/(sqrt(species_molecular_properties[7]));
    }
    else
    {
        for (int si = 0; si < d_num_species; si++)
        {
            getSpeciesMolecularProperties(species_molecular_properties_ptr, si);
            
            const double mu_v_i = d_equation_of_bulk_viscosity->
                getBulkViscosity(
                    pressure,
                    temperature,
                    species_molecular_properties_const_ptr);
            
            num += mu_v_i*(*(mass_fraction[si]))/(sqrt(species_molecular_properties[7]));
            den += *(mass_fraction[si])/(sqrt(species_molecular_properties[7]));
        }
    }
    
    mu_v = num/den;
    
    return mu_v;
}


/*
 * Compute the bulk viscosity of the mixture with isobaric assumption.
 */
double
EquationOfBulkViscosityMixingRulesCramer::getBulkViscosity(
    const double* const pressure,
    const std::vector<const double*>& temperature,
    const std::vector<const double*>& mass_fraction,
    const std::vector<const double*>& volume_fraction) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOBARIC);
    TBOX_ASSERT((static_cast<int>(temperature.size()) == d_num_species);
    TBOX_ASSERT((static_cast<int>(volume_fraction.size()) == d_num_species) ||
                (static_cast<int>(volume_fraction.size()) == d_num_species - 1));
#endif
    
    NULL_USE(mass_fraction);
    
    double mu_v = 0.0;
    
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
            
            const double mu_v_i = d_equation_of_bulk_viscosity->
                getBulkViscosity(
                    pressure,
                    temperature[si],
                    species_molecular_properties_const_ptr);
            
            mu_v += *(volume_fraction[si])*mu_v_i;
            
            // Compute the volume fraction of the last species.
            Z_last -= *(volume_fraction[si]);
        }
        
        /*
         * Add the contribution from the last species.
         */
        
        getSpeciesMolecularProperties(species_molecular_properties_ptr, d_num_species - 1);
            
        const double mu_v_last = d_equation_of_bulk_viscosity->
            getBulkViscosity(
                pressure,
                temperature[d_num_species - 1],
                species_molecular_properties_const_ptr);
        
        mu_v += Z_last*mu_v_last;
    }
    else
    {
        for (int si = 0; si < d_num_species; si++)
        {
            getSpeciesMolecularProperties(species_molecular_properties_ptr, si);
            
            const double mu_v_i = d_equation_of_bulk_viscosity->
                getBulkViscosity(
                    pressure,
                    temperature[si],
                    species_molecular_properties_const_ptr);
            
            mu_v += *(volume_fraction[si])*mu_v_i;
        }
    }
    
    return mu_v;
}


/*
 * Get the molecular properties of a species.
 */
void
EquationOfBulkViscosityMixingRulesCramer::getSpeciesMolecularProperties(
    std::vector<double*>& species_molecular_properties,
    const int& species_index) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(species_molecular_properties.size()) >= 8);
    TBOX_ASSERT(species_index >= 0);
    TBOX_ASSERT(species_index < d_num_species);
#endif
    
    *(species_molecular_properties[0]) = d_species_gamma[species_index];
    
    *(species_molecular_properties[1]) = d_species_A_r[species_index];
    *(species_molecular_properties[2]) = d_species_B_r[species_index];
    
    *(species_molecular_properties[3]) = d_species_c_v_v[species_index];
    *(species_molecular_properties[4]) = d_species_A_v[species_index];
    *(species_molecular_properties[5]) = d_species_B_v[species_index];
    *(species_molecular_properties[6]) = d_species_C_v[species_index];
    
    *(species_molecular_properties[7]) = d_species_M[species_index];
}
