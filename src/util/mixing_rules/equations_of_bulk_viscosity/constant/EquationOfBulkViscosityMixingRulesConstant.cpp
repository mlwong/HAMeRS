#include "util/mixing_rules/equations_of_bulk_viscosity/constant/EquationOfBulkViscosityMixingRulesConstant.hpp"

EquationOfBulkViscosityMixingRulesConstant::EquationOfBulkViscosityMixingRulesConstant(
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
    d_equation_of_bulk_viscosity.reset(new EquationOfBulkViscosityConstant(
        "d_equation_of_bulk_viscosity",
        dim));
    
    /*
     * Get the bulk viscosity of each species from the database.
     */
    
    if (equation_of_bulk_viscosity_mixing_rules_db->keyExists("species_mu_v"))
    {
        size_t species_mu_v_array_size =
            equation_of_bulk_viscosity_mixing_rules_db->getArraySize("species_mu_v");
        if (static_cast<int>(species_mu_v_array_size) == d_num_species)
        {
            d_species_mu_v =
                equation_of_bulk_viscosity_mixing_rules_db->getDoubleVector("species_mu_v");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'species_mu_v' entries must be equal to 'num_species'."
                << std::endl);
        }
    }
    else if (equation_of_bulk_viscosity_mixing_rules_db->keyExists("d_species_mu_v"))
    {
        size_t species_mu_v_array_size =
            equation_of_bulk_viscosity_mixing_rules_db->getArraySize("d_species_mu_v");
        if (static_cast<int>(species_mu_v_array_size) == d_num_species)
        {
            d_species_mu_v =
                equation_of_bulk_viscosity_mixing_rules_db->getDoubleVector("d_species_mu_v");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'd_species_mu_v' entries must be equal to 'd_num_species'."
                << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Key data 'species_mu_v'/'d_species_mu_v'"
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
EquationOfBulkViscosityMixingRulesConstant::printClassData(
    std::ostream& os) const
{
    os << "\nPrint EquationOfBulkViscosityMixingRulesConstant object..."
       << std::endl;
    
    os << std::endl;
    os << "EquationOfBulkViscosityMixingRulesConstant: this = "
       << (EquationOfBulkViscosityMixingRulesConstant *)this
       << std::endl;
    
    os << "d_object_name = "
       << d_object_name
       << std::endl;
    
    os << "d_mixing_closure_model = "
       << d_mixing_closure_model
       << std::endl;
    
    /*
     * Print the bulk viscosity of each species.
     */
    
    os << "d_species_mu_v = ";
    for (int si = 0; si < d_num_species - 1; si++)
    {
        os << d_species_mu_v[si] << ", ";
    }
    os << d_species_mu_v[d_num_species - 1];
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
EquationOfBulkViscosityMixingRulesConstant::putToRestart(
    const boost::shared_ptr<tbox::Database>& restart_db) const
{
    restart_db->putDoubleVector("d_species_mu_v", d_species_mu_v);
    restart_db->putDoubleVector("d_species_M", d_species_M);
}


/*
 * Compute the bulk viscosity of the mixture with isothermal and isobaric assumptions.
 */
double
EquationOfBulkViscosityMixingRulesConstant::getBulkViscosity(
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
            
            num += mu_v_i*(*(mass_fraction[si]))/(sqrt(species_molecular_properties[1]));
            den += *(mass_fraction[si])/(sqrt(species_molecular_properties[1]));
            
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
        
        num += mu_v_last*Y_last/(sqrt(species_molecular_properties[1]));
        den += Y_last/(sqrt(species_molecular_properties[1]));
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
            
            num += mu_v_i*(*(mass_fraction[si]))/(sqrt(species_molecular_properties[1]));
            den += *(mass_fraction[si])/(sqrt(species_molecular_properties[1]));
        }
    }
    
    mu_v = num/den;
    
    return mu_v;
}


/*
 * Compute the bulk viscosity of the mixture with isobaric assumption.
 */
double
EquationOfBulkViscosityMixingRulesConstant::getBulkViscosity(
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
EquationOfBulkViscosityMixingRulesConstant::getSpeciesMolecularProperties(
    std::vector<double*>& species_molecular_properties,
    const int& species_index) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(species_molecular_properties.size()) >= 2);
    TBOX_ASSERT(species_index >= 0);
    TBOX_ASSERT(species_index < d_num_species);
#endif
    
    *(species_molecular_properties[0]) = d_species_mu_v[species_index];
    *(species_molecular_properties[1]) = d_species_M[species_index];
}
