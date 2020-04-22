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
 * Compute the bulk viscosity of the mixture with isothermal and isobaric equilibrium assumptions.
 */
double
EquationOfBulkViscosityMixingRulesCramer::getBulkViscosity(
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
    
    double mu_v = double(0);
    
    double num = double(0);
    double den = double(0);
    
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
            
            const double mu_v_i = d_equation_of_bulk_viscosity->
                getBulkViscosity(
                    pressure,
                    temperature,
                    species_molecular_properties_const_ptr);
            
            const double weight = *(mass_fractions[si])/(sqrt(species_molecular_properties[7]));
            
            num += mu_v_i*weight;
            den += weight;
        }
    }
    else if (static_cast<int>(mass_fractions.size()) == d_num_species - 1)
    {
        double Y_last = double(1);
        
        for (int si = 0; si < d_num_species - 1; si++)
        {
            getSpeciesMolecularProperties(species_molecular_properties_ptr, si);
            
            const double mu_v_i = d_equation_of_bulk_viscosity->
                getBulkViscosity(
                    pressure,
                    temperature,
                    species_molecular_properties_const_ptr);
            
            const double weight = *(mass_fractions[si])/(sqrt(species_molecular_properties[7]));
            
            num += mu_v_i*weight;
            den += weight;
            
            // Compute the mass fraction of the last species.
            Y_last -= *(mass_fractions[si]);
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
        
        const double weight = Y_last/(sqrt(species_molecular_properties[7]));
        
        num += mu_v_last*weight;
        den += weight;
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Number of mass fractions provided is not"
            << " equal to the total number of species or (total number of species - 1)."
            << std::endl);
    }
    
    mu_v = num/den;
    
    return mu_v;
}


/*
 * Compute the bulk viscosity of the mixture with isothermal and isobaric equilibrium assumptions.
 */
void
EquationOfBulkViscosityMixingRulesCramer::computeBulkViscosity(
    boost::shared_ptr<pdat::CellData<double> >& data_bulk_viscosity,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
    const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT((d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC) ||
                (d_mixing_closure_model == MIXING_CLOSURE_MODEL::NO_MODEL && d_num_species == 1));
    
    TBOX_ASSERT(data_bulk_viscosity);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_temperature);
    TBOX_ASSERT(data_mass_fractions);
    
    TBOX_ASSERT((data_mass_fractions->getDepth() == d_num_species) ||
                (data_mass_fractions->getDepth() == d_num_species - 1));
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_bulk_viscosity = data_bulk_viscosity->getGhostBox();
    const hier::IntVector ghostcell_dims_bulk_viscosity = ghost_box_bulk_viscosity.numberCells();
    
    const hier::Box ghost_box_mass_fractions = data_mass_fractions->getGhostBox();
    const hier::IntVector ghostcell_dims_mass_fractions = ghost_box_mass_fractions.numberCells();
    
    // Delcare data containers for bulk viscosity of a species, denominator and numerator.
    boost::shared_ptr<pdat::CellData<double> > data_bulk_viscosity_species;
    boost::shared_ptr<pdat::CellData<double> > data_den;
    boost::shared_ptr<pdat::CellData<double> > data_num;
    
    // Declare data container for last mass fraction.
    boost::shared_ptr<pdat::CellData<double> > data_mass_fractions_last;
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets of all data and dimensions of the ghost cell box for denominator,
     * numerator and last mass fraction and allocate memory.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_bulk_viscosity(d_dim);
    hier::IntVector offset_mass_fractions(d_dim);
    hier::IntVector offset_min(d_dim);
    
    hier::IntVector ghostcell_dims_min(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_bulk_viscosity = data_bulk_viscosity->getGhostCellWidth();
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
        const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
        
        // Get the interior box and the dimensions of box that covers the interior of patch.
        const hier::Box interior_box = data_bulk_viscosity->getBox();
        const hier::IntVector interior_dims = interior_box.numberCells();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_pressure->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_temperature->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_mass_fractions->getBox().isSpatiallyEqual(interior_box));
#endif
        
        /*
         * Get the minimum number of ghost cells and the dimensions of the ghost cell box for denominator,
         * numerator and last mass fraction.
         */
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_bulk_viscosity;
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_temperature, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_mass_fractions, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_min = num_ghosts_min;
        offset_bulk_viscosity = num_ghosts_bulk_viscosity;
        offset_mass_fractions = num_ghosts_mass_fractions;
        
        ghostcell_dims_min = interior_dims + num_ghosts_min*2;
        
        data_bulk_viscosity_species = boost::make_shared<pdat::CellData<double> >(interior_box, 1, num_ghosts_min);
        data_den = boost::make_shared<pdat::CellData<double> >(interior_box, 1, num_ghosts_min);
        data_num = boost::make_shared<pdat::CellData<double> >(interior_box, 1, num_ghosts_min);
        
        if (data_mass_fractions->getDepth() == d_num_species - 1)
        {
            data_mass_fractions_last = boost::make_shared<pdat::CellData<double> >(interior_box, 1, num_ghosts_min);
        }
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
        TBOX_ASSERT(data_bulk_viscosity->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_temperature->getGhostBox().contains(domain));
        TBOX_ASSERT(data_mass_fractions->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_min = hier::IntVector::getZero(d_dim);
        offset_bulk_viscosity = domain.lower() - ghost_box_bulk_viscosity.lower();
        offset_mass_fractions = domain.lower() - ghost_box_mass_fractions.lower();
        
        ghostcell_dims_min = domain_dims;
        
        data_bulk_viscosity_species =
            boost::make_shared<pdat::CellData<double> >(domain, 1, hier::IntVector::getZero(d_dim));
        data_den = boost::make_shared<pdat::CellData<double> >(domain, 1, hier::IntVector::getZero(d_dim));
        data_num = boost::make_shared<pdat::CellData<double> >(domain, 1, hier::IntVector::getZero(d_dim));
        
        if (data_mass_fractions->getDepth() == d_num_species - 1)
        {
            data_mass_fractions_last =
                boost::make_shared<pdat::CellData<double> >(domain, 1, hier::IntVector::getZero(d_dim));
        }
    }
    
    // Declare data containers for species molecular properties.
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
    
    /*
     * Get the pointers to the cell data of mixture bulk viscosity, species bulk viscosity, denominator
     * and numerator.
     */
    
    double* mu_v = data_bulk_viscosity->getPointer(0);
    double* mu_v_i = data_bulk_viscosity_species->getPointer(0);
    double* den = data_den->getPointer(0);
    double* num = data_num->getPointer(0);
    
    /*
     * Fill zeros for denominator and numerator.
     */
    
    data_den->fillAll(double(0));
    data_num->fillAll(double(0));
    
    if (data_mass_fractions->getDepth() == d_num_species)
    {
        /*
         * Get the pointers to the cell data of mass fractions.
         */
        
        std::vector<double*> Y;
        Y.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            Y.push_back(data_mass_fractions->getPointer(si));
        }
        
        if (d_dim == tbox::Dimension(1))
        {
            /*
             * Get the local lower index, numbers of cells in each dimension and offsets.
             */
            
            const int domain_lo_0 = domain_lo[0];
            const int domain_dim_0 = domain_dims[0];
            
            const int offset0_bulk_viscosity = offset_bulk_viscosity[0];
            const int offset0_min = offset_min[0];
            const int offset0_mass_fractions = offset_mass_fractions[0];
            
            // Compute the mixture bulk viscosity field.
            for (int si = 0; si < d_num_species; si++)
            {
                getSpeciesMolecularProperties(species_molecular_properties_ptr, si);
                
                d_equation_of_bulk_viscosity->
                    computeBulkViscosity(
                        data_bulk_viscosity_species,
                        data_pressure,
                        data_temperature,
                        species_molecular_properties_const_ptr,
                        domain);
                
                const double factor = double(1)/(sqrt(species_molecular_properties[7]));
                
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_min = i + offset0_min;
                    const int idx_mass_fractions = i + offset0_mass_fractions;
                    
                    const double weight = Y[si][idx_mass_fractions]*factor;
                    
                    num[idx_min] += mu_v_i[idx_min]*weight;
                    den[idx_min] += weight;
                }
            }
            
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_bulk_viscosity = i + offset0_bulk_viscosity;
                const int idx_min = i + offset0_min;
                
                mu_v[idx_bulk_viscosity] = num[idx_min]/den[idx_min];
            }
        }
        else if (d_dim == tbox::Dimension(2))
        {
            /*
             * Get the local lower indices, numbers of cells in each dimension and offsets.
             */
            
            const int domain_lo_0 = domain_lo[0];
            const int domain_lo_1 = domain_lo[1];
            const int domain_dim_0 = domain_dims[0];
            const int domain_dim_1 = domain_dims[1];
            
            const int offset0_bulk_viscosity = offset_bulk_viscosity[0];
            const int offset1_bulk_viscosity = offset_bulk_viscosity[1];
            const int ghostcell_dim_0_bulk_viscosity = ghostcell_dims_bulk_viscosity[0];
            
            const int offset0_min = offset_min[0];
            const int offset1_min = offset_min[1];
            const int ghostcell_dim_0_min = ghostcell_dims_min[0];
            
            const int offset0_mass_fractions = offset_mass_fractions[0];
            const int offset1_mass_fractions = offset_mass_fractions[1];
            const int ghostcell_dim_0_mass_fractions = ghostcell_dims_mass_fractions[0];
            
            // Compute the mixture bulk viscosity field.
            for (int si = 0; si < d_num_species; si++)
            {
                getSpeciesMolecularProperties(species_molecular_properties_ptr, si);
                
                d_equation_of_bulk_viscosity->
                    computeBulkViscosity(
                        data_bulk_viscosity_species,
                        data_pressure,
                        data_temperature,
                        species_molecular_properties_const_ptr,
                        domain);
                
                const double factor = double(1)/(sqrt(species_molecular_properties[7]));
                
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_min = (i + offset0_min) +
                            (j + offset1_min)*ghostcell_dim_0_min;
                        
                        const int idx_mass_fractions = (i + offset0_mass_fractions) +
                            (j + offset1_mass_fractions)*ghostcell_dim_0_mass_fractions;
                        
                        const double weight = Y[si][idx_mass_fractions]*factor;
                        
                        num[idx_min] += mu_v_i[idx_min]*weight;
                        den[idx_min] += weight;
                    }
                }
            }
            
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_bulk_viscosity = (i + offset0_bulk_viscosity) +
                        (j + offset1_bulk_viscosity)*ghostcell_dim_0_bulk_viscosity;
                    
                    const int idx_min = (i + offset0_min) +
                        (j + offset1_min)*ghostcell_dim_0_min;
                    
                    mu_v[idx_bulk_viscosity] = num[idx_min]/den[idx_min];
                }
            }
        }
        else if (d_dim == tbox::Dimension(3))
        {
            /*
             * Get the local lower indices, numbers of cells in each dimension and offsets.
             */
            
            const int domain_lo_0 = domain_lo[0];
            const int domain_lo_1 = domain_lo[1];
            const int domain_lo_2 = domain_lo[2];
            const int domain_dim_0 = domain_dims[0];
            const int domain_dim_1 = domain_dims[1];
            const int domain_dim_2 = domain_dims[2];
            
            const int offset0_bulk_viscosity = offset_bulk_viscosity[0];
            const int offset1_bulk_viscosity = offset_bulk_viscosity[1];
            const int offset2_bulk_viscosity = offset_bulk_viscosity[2];
            const int ghostcell_dim_0_bulk_viscosity = ghostcell_dims_bulk_viscosity[0];
            const int ghostcell_dim_1_bulk_viscosity = ghostcell_dims_bulk_viscosity[1];
            
            const int offset0_min = offset_min[0];
            const int offset1_min = offset_min[1];
            const int offset2_min = offset_min[2];
            const int ghostcell_dim_0_min = ghostcell_dims_min[0];
            const int ghostcell_dim_1_min = ghostcell_dims_min[1];
            
            const int offset0_mass_fractions = offset_mass_fractions[0];
            const int offset1_mass_fractions = offset_mass_fractions[1];
            const int offset2_mass_fractions = offset_mass_fractions[2];
            const int ghostcell_dim_0_mass_fractions = ghostcell_dims_mass_fractions[0];
            const int ghostcell_dim_1_mass_fractions = ghostcell_dims_mass_fractions[1];
            
            // Compute the mixture bulk viscosity field.
            for (int si = 0; si < d_num_species; si++)
            {
                getSpeciesMolecularProperties(species_molecular_properties_ptr, si);
                
                d_equation_of_bulk_viscosity->
                    computeBulkViscosity(
                        data_bulk_viscosity_species,
                        data_pressure,
                        data_temperature,
                        species_molecular_properties_const_ptr,
                        domain);
                
                const double factor = double(1)/(sqrt(species_molecular_properties[7]));
                
                for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                {
                    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_min = (i + offset0_min) +
                                (j + offset1_min)*ghostcell_dim_0_min +
                                (k + offset2_min)*ghostcell_dim_0_min*
                                    ghostcell_dim_1_min;
                            
                            const int idx_mass_fractions = (i + offset0_mass_fractions) +
                                (j + offset1_mass_fractions)*ghostcell_dim_0_mass_fractions +
                                (k + offset2_mass_fractions)*ghostcell_dim_0_mass_fractions*
                                    ghostcell_dim_1_mass_fractions;
                            
                            const double weight = Y[si][idx_mass_fractions]*factor;
                            
                            num[idx_min] += mu_v_i[idx_min]*weight;
                            den[idx_min] += weight;
                        }
                    }
                }
            }
            
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_bulk_viscosity = (i + offset0_bulk_viscosity) +
                            (j + offset1_bulk_viscosity)*ghostcell_dim_0_bulk_viscosity +
                            (k + offset2_bulk_viscosity)*ghostcell_dim_0_bulk_viscosity*
                                ghostcell_dim_1_bulk_viscosity;
                        
                        const int idx_min = (i + offset0_min) +
                            (j + offset1_min)*ghostcell_dim_0_min +
                            (k + offset2_min)*ghostcell_dim_0_min*
                                ghostcell_dim_1_min;
                        
                        mu_v[idx_bulk_viscosity] = num[idx_min]/den[idx_min];
                    }
                }
            }
        }
    }
    else if (data_mass_fractions->getDepth() == d_num_species - 1)
    {
        data_mass_fractions_last->fillAll(double(1));
        
        /*
         * Get the pointers to the cell data of mass fractions.
         */
        
        std::vector<double*> Y;
        Y.reserve(d_num_species - 1);
        for (int si = 0; si < d_num_species - 1; si++)
        {
            Y.push_back(data_mass_fractions->getPointer(si));
        }
        
        double* Y_last = data_mass_fractions_last->getPointer(0);
        
        if (d_dim == tbox::Dimension(1))
        {
            /*
             * Get the local lower index, numbers of cells in each dimension and offsets.
             */
            
            const int domain_lo_0 = domain_lo[0];
            const int domain_dim_0 = domain_dims[0];
            
            const int offset0_bulk_viscosity = offset_bulk_viscosity[0];
            const int offset0_min = offset_min[0];
            const int offset0_mass_fractions = offset_mass_fractions[0];
            
            // Compute the mixture bulk viscosity field.
            for (int si = 0; si < d_num_species - 1; si++)
            {
                getSpeciesMolecularProperties(species_molecular_properties_ptr, si);
                
                d_equation_of_bulk_viscosity->
                    computeBulkViscosity(
                        data_bulk_viscosity_species,
                        data_pressure,
                        data_temperature,
                        species_molecular_properties_const_ptr,
                        domain);
                
                const double factor = double(1)/(sqrt(species_molecular_properties[7]));
                
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_min = i + offset0_min;
                    const int idx_mass_fractions = i + offset0_mass_fractions;
                    
                    const double weight = Y[si][idx_mass_fractions]*factor;
                    
                    num[idx_min] += mu_v_i[idx_min]*weight;
                    den[idx_min] += weight;
                    
                    // Compute the mass fraction of the last species.
                    Y_last[idx_min] -= Y[si][idx_mass_fractions];
                }
            }
            
            getSpeciesMolecularProperties(species_molecular_properties_ptr, d_num_species - 1);
            
            d_equation_of_bulk_viscosity->
                computeBulkViscosity(
                    data_bulk_viscosity_species,
                    data_pressure,
                    data_temperature,
                    species_molecular_properties_const_ptr,
                    domain);
            
            const double factor = double(1)/(sqrt(species_molecular_properties[7]));
            
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_bulk_viscosity = i + offset0_bulk_viscosity;
                const int idx_min = i + offset0_min;
                
                const double weight = Y_last[idx_min]*factor;
                
                num[idx_min] += mu_v_i[idx_min]*weight;
                den[idx_min] += weight;
                
                mu_v[idx_bulk_viscosity] = num[idx_min]/den[idx_min];
            }
        }
        else if (d_dim == tbox::Dimension(2))
        {
            /*
             * Get the local lower indices, numbers of cells in each dimension and offsets.
             */
            
            const int domain_lo_0 = domain_lo[0];
            const int domain_lo_1 = domain_lo[1];
            const int domain_dim_0 = domain_dims[0];
            const int domain_dim_1 = domain_dims[1];
            
            const int offset0_bulk_viscosity = offset_bulk_viscosity[0];
            const int offset1_bulk_viscosity = offset_bulk_viscosity[1];
            const int ghostcell_dim_0_bulk_viscosity = ghostcell_dims_bulk_viscosity[0];
            
            const int offset0_min = offset_min[0];
            const int offset1_min = offset_min[1];
            const int ghostcell_dim_0_min = ghostcell_dims_min[0];
            
            const int offset0_mass_fractions = offset_mass_fractions[0];
            const int offset1_mass_fractions = offset_mass_fractions[1];
            const int ghostcell_dim_0_mass_fractions = ghostcell_dims_mass_fractions[0];
            
            // Compute the mixture bulk viscosity field.
            for (int si = 0; si < d_num_species - 1; si++)
            {
                getSpeciesMolecularProperties(species_molecular_properties_ptr, si);
                
                d_equation_of_bulk_viscosity->
                    computeBulkViscosity(
                        data_bulk_viscosity_species,
                        data_pressure,
                        data_temperature,
                        species_molecular_properties_const_ptr,
                        domain);
                
                const double factor = double(1)/(sqrt(species_molecular_properties[7]));
                
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_min = (i + offset0_min) +
                            (j + offset1_min)*ghostcell_dim_0_min;
                        
                        const int idx_mass_fractions = (i + offset0_mass_fractions) +
                            (j + offset1_mass_fractions)*ghostcell_dim_0_mass_fractions;
                        
                        const double weight = Y[si][idx_mass_fractions]*factor;
                        
                        num[idx_min] += mu_v_i[idx_min]*weight;
                        den[idx_min] += weight;
                        
                        // Compute the mass fraction of the last species.
                        Y_last[idx_min] -= Y[si][idx_mass_fractions];
                    }
                }
            }
            
            getSpeciesMolecularProperties(species_molecular_properties_ptr, d_num_species - 1);
            
            d_equation_of_bulk_viscosity->
                computeBulkViscosity(
                    data_bulk_viscosity_species,
                    data_pressure,
                    data_temperature,
                    species_molecular_properties_const_ptr,
                    domain);
            
            const double factor = double(1)/(sqrt(species_molecular_properties[7]));
            
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_bulk_viscosity = (i + offset0_bulk_viscosity) +
                        (j + offset1_bulk_viscosity)*ghostcell_dim_0_bulk_viscosity;
                    
                    const int idx_min = (i + offset0_min) +
                        (j + offset1_min)*ghostcell_dim_0_min;
                    
                    const double weight = Y_last[idx_min]*factor;
                    
                    num[idx_min] += mu_v_i[idx_min]*weight;
                    den[idx_min] += weight;
                    
                    mu_v[idx_bulk_viscosity] = num[idx_min]/den[idx_min];
                }
            }
        }
        else if (d_dim == tbox::Dimension(3))
        {
            /*
             * Get the local lower indices, numbers of cells in each dimension and offsets.
             */
            
            const int domain_lo_0 = domain_lo[0];
            const int domain_lo_1 = domain_lo[1];
            const int domain_lo_2 = domain_lo[2];
            const int domain_dim_0 = domain_dims[0];
            const int domain_dim_1 = domain_dims[1];
            const int domain_dim_2 = domain_dims[2];
            
            const int offset0_bulk_viscosity = offset_bulk_viscosity[0];
            const int offset1_bulk_viscosity = offset_bulk_viscosity[1];
            const int offset2_bulk_viscosity = offset_bulk_viscosity[2];
            const int ghostcell_dim_0_bulk_viscosity = ghostcell_dims_bulk_viscosity[0];
            const int ghostcell_dim_1_bulk_viscosity = ghostcell_dims_bulk_viscosity[1];
            
            const int offset0_min = offset_min[0];
            const int offset1_min = offset_min[1];
            const int offset2_min = offset_min[2];
            const int ghostcell_dim_0_min = ghostcell_dims_min[0];
            const int ghostcell_dim_1_min = ghostcell_dims_min[1];
            
            const int offset0_mass_fractions = offset_mass_fractions[0];
            const int offset1_mass_fractions = offset_mass_fractions[1];
            const int offset2_mass_fractions = offset_mass_fractions[2];
            const int ghostcell_dim_0_mass_fractions = ghostcell_dims_mass_fractions[0];
            const int ghostcell_dim_1_mass_fractions = ghostcell_dims_mass_fractions[1];
            
            // Compute the mixture bulk viscosity field.
            for (int si = 0; si < d_num_species - 1; si++)
            {
                getSpeciesMolecularProperties(species_molecular_properties_ptr, si);
                
                d_equation_of_bulk_viscosity->
                    computeBulkViscosity(
                        data_bulk_viscosity_species,
                        data_pressure,
                        data_temperature,
                        species_molecular_properties_const_ptr,
                        domain);
                
                const double factor = double(1)/(sqrt(species_molecular_properties[7]));
                
                for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                {
                    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_min = (i + offset0_min) +
                                (j + offset1_min)*ghostcell_dim_0_min +
                                (k + offset2_min)*ghostcell_dim_0_min*
                                    ghostcell_dim_1_min;
                            
                            const int idx_mass_fractions = (i + offset0_mass_fractions) +
                                (j + offset1_mass_fractions)*ghostcell_dim_0_mass_fractions +
                                (k + offset2_mass_fractions)*ghostcell_dim_0_mass_fractions*
                                    ghostcell_dim_1_mass_fractions;
                            
                            const double weight = Y[si][idx_mass_fractions]*factor;
                            
                            num[idx_min] += mu_v_i[idx_min]*weight;
                            den[idx_min] += weight;
                            
                            // Compute the mass fraction of the last species.
                            Y_last[idx_min] -= Y[si][idx_mass_fractions];
                        }
                    }
                }
            }
            
            getSpeciesMolecularProperties(species_molecular_properties_ptr, d_num_species - 1);
            
            d_equation_of_bulk_viscosity->
                computeBulkViscosity(
                    data_bulk_viscosity_species,
                    data_pressure,
                    data_temperature,
                    species_molecular_properties_const_ptr,
                    domain);
            
            const double factor = double(1)/(sqrt(species_molecular_properties[7]));
            
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_bulk_viscosity = (i + offset0_bulk_viscosity) +
                            (j + offset1_bulk_viscosity)*ghostcell_dim_0_bulk_viscosity +
                            (k + offset2_bulk_viscosity)*ghostcell_dim_0_bulk_viscosity*
                                ghostcell_dim_1_bulk_viscosity;
                        
                        const int idx_min = (i + offset0_min) +
                            (j + offset1_min)*ghostcell_dim_0_min +
                            (k + offset2_min)*ghostcell_dim_0_min*
                                ghostcell_dim_1_min;
                        
                        const double weight = Y_last[idx_min]*factor;
                        
                        num[idx_min] += mu_v_i[idx_min]*weight;
                        den[idx_min] += weight;
                        
                        mu_v[idx_bulk_viscosity] = num[idx_min]/den[idx_min];
                    }
                }
            }
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Number of components in the data of mass fractions provided is not"
            << " equal to the total number of species or (total number of species - 1)."
            << std::endl);
    }
}


/*
 * Compute the bulk viscosity of the mixture with isobaric equilibrium assumption.
 */
double
EquationOfBulkViscosityMixingRulesCramer::getBulkViscosity(
    const double* const pressure,
    const std::vector<const double*>& species_temperatures,
    const std::vector<const double*>& mass_fractions,
    const std::vector<const double*>& volume_fractions) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOBARIC);
    TBOX_ASSERT((static_cast<int>(species_temperatures.size()) == d_num_species));
    TBOX_ASSERT((static_cast<int>(volume_fractions.size()) == d_num_species) ||
                (static_cast<int>(volume_fractions.size()) == d_num_species - 1));
#endif
    
    NULL_USE(mass_fractions);
    
    double mu_v = double(0);
    
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
    
    if (static_cast<int>(volume_fractions.size()) == d_num_species - 1)
    {
        for (int si = 0; si < d_num_species; si++)
        {
            getSpeciesMolecularProperties(species_molecular_properties_ptr, si);
            
            const double mu_v_i = d_equation_of_bulk_viscosity->
                getBulkViscosity(
                    pressure,
                    species_temperatures[si],
                    species_molecular_properties_const_ptr);
            
            mu_v += *(volume_fractions[si])*mu_v_i;
        }
    }
    else if (static_cast<int>(volume_fractions.size()) == d_num_species - 1)
    {
        double Z_last = double(1);
        
        for (int si = 0; si < d_num_species - 1; si++)
        {
            getSpeciesMolecularProperties(species_molecular_properties_ptr, si);
            
            const double mu_v_i = d_equation_of_bulk_viscosity->
                getBulkViscosity(
                    pressure,
                    species_temperatures[si],
                    species_molecular_properties_const_ptr);
            
            mu_v += *(volume_fractions[si])*mu_v_i;
            
            // Compute the volume fraction of the last species.
            Z_last -= *(volume_fractions[si]);
        }
        
        /*
         * Add the contribution from the last species.
         */
        
        getSpeciesMolecularProperties(species_molecular_properties_ptr, d_num_species - 1);
            
        const double mu_v_last = d_equation_of_bulk_viscosity->
            getBulkViscosity(
                pressure,
                species_temperatures[d_num_species - 1],
                species_molecular_properties_const_ptr);
        
        mu_v += Z_last*mu_v_last;
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Number of volume fractions provided is not"
            << " equal to the total number of species or (total number of species - 1)."
            << std::endl);
    }
    
    return mu_v;
}


/*
 * Compute the bulk viscosity of the mixture with isobaric equilibrium assumption.
 */
void
EquationOfBulkViscosityMixingRulesCramer::computeBulkViscosity(
    boost::shared_ptr<pdat::CellData<double> >& data_bulk_viscosity,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const std::vector<boost::shared_ptr<pdat::CellData<double> > >& data_species_temperatures,
    const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
    const boost::shared_ptr<pdat::CellData<double> >& data_volume_fractions,
    const hier::Box& domain) const
{
    // NEED TO IMPLEMENT DUE TO VECTOR OF SPECIES TEMPERATURES
}


/*
 * Compute the bulk viscosity of the mixture with isobaric equilibrium assumption.
 */
// NEED TO REMOVE DUE TO VECTOR OF SPECIES TEMPERATURES
void
EquationOfBulkViscosityMixingRulesCramer::computeBulkViscosity(
    boost::shared_ptr<pdat::CellData<double> >& data_bulk_viscosity,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_species_temperatures,
    const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
    const boost::shared_ptr<pdat::CellData<double> >& data_volume_fractions,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOBARIC);
    
    TBOX_ASSERT(data_bulk_viscosity);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_species_temperatures);
    TBOX_ASSERT(data_volume_fractions);
    
    TBOX_ASSERT(data_species_temperatures->getDepth() == d_num_species);
    TBOX_ASSERT((data_volume_fractions->getDepth() == d_num_species) ||
                (data_volume_fractions->getDepth() == d_num_species - 1));
#endif
    
    NULL_USE(data_mass_fractions);
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_bulk_viscosity = data_bulk_viscosity->getGhostBox();
    const hier::IntVector ghostcell_dims_bulk_viscosity = ghost_box_bulk_viscosity.numberCells();
    
    const hier::Box ghost_box_volume_fractions = data_volume_fractions->getGhostBox();
    const hier::IntVector ghostcell_dims_volume_fractions = ghost_box_volume_fractions.numberCells();
    
    // Delcare data containers for bulk viscosity and temperature of a species.
    boost::shared_ptr<pdat::CellData<double> > data_bulk_viscosity_species;
    boost::shared_ptr<pdat::CellData<double> > data_temperature_species;
    
    // Declare data container for last volume fraction.
    boost::shared_ptr<pdat::CellData<double> > data_volume_fractions_last;
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets of all data and dimensions of the ghost cell box for bulk viscosity,
     * temperature of a species and last volume fraction and allocate memory.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_bulk_viscosity(d_dim);
    hier::IntVector offset_volume_fractions(d_dim);
    hier::IntVector offset_min(d_dim);
    
    hier::IntVector ghostcell_dims_min(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_bulk_viscosity = data_bulk_viscosity->getGhostCellWidth();
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        const hier::IntVector num_ghosts_species_temperatures = data_species_temperatures->getGhostCellWidth();
        const hier::IntVector num_ghosts_volume_fractions = data_volume_fractions->getGhostCellWidth();
        
        // Get the interior box and the dimensions of box that covers the interior of patch.
        const hier::Box interior_box = data_bulk_viscosity->getBox();
        const hier::IntVector interior_dims = interior_box.numberCells();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_pressure->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_species_temperatures->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_volume_fractions->getBox().isSpatiallyEqual(interior_box));
#endif
        
        /*
         * Get the minimum number of ghost cells and the dimensions of the ghost cell box for bulk viscosity,
         * temperature of a species and last volume fraction.
         */
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_bulk_viscosity;
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_species_temperatures, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_volume_fractions, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_min = num_ghosts_min;
        offset_bulk_viscosity = num_ghosts_bulk_viscosity;
        offset_volume_fractions = num_ghosts_volume_fractions;
        
        ghostcell_dims_min = interior_dims + num_ghosts_min*2;
        
        data_bulk_viscosity_species = boost::make_shared<pdat::CellData<double> >(interior_box, 1, num_ghosts_min);
        data_temperature_species = boost::make_shared<pdat::CellData<double> >(interior_box, 1, num_ghosts_min);
        
        if (data_volume_fractions->getDepth() == d_num_species - 1)
        {
            data_volume_fractions_last = boost::make_shared<pdat::CellData<double> >(interior_box, 1, num_ghosts_min);
        }
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
        TBOX_ASSERT(data_bulk_viscosity->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_species_temperatures->getGhostBox().contains(domain));
        TBOX_ASSERT(data_volume_fractions->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_min = hier::IntVector::getZero(d_dim);
        offset_bulk_viscosity = domain.lower() - ghost_box_bulk_viscosity.lower();
        offset_volume_fractions = domain.lower() - ghost_box_volume_fractions.lower();
        
        ghostcell_dims_min = domain_dims;
        
        data_bulk_viscosity_species =
            boost::make_shared<pdat::CellData<double> >(domain, 1, hier::IntVector::getZero(d_dim));
        data_temperature_species =
            boost::make_shared<pdat::CellData<double> >(domain, 1, hier::IntVector::getZero(d_dim));
        
        if (data_volume_fractions->getDepth() == d_num_species - 1)
        {
            data_volume_fractions_last =
                boost::make_shared<pdat::CellData<double> >(domain, 1, hier::IntVector::getZero(d_dim));
        }
    }
    
    // Delcare data containers for species molecular properties.
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
    
    /*
     * Get the pointers to the cell data of mixture bulk viscosity, species bulk viscosity, denominator
     * and numerator.
     */
    
    double* mu_v = data_bulk_viscosity->getPointer(0);
    double* mu_v_i = data_bulk_viscosity_species->getPointer(0);
    
    /*
     * Fill zeros for mixture bulk viscosity.
     */
    
    if (domain.empty())
    {
        data_bulk_viscosity->fillAll(double(0));
    }
    else
    {
        data_bulk_viscosity->fillAll(double(0), domain);
    }
    
    if (data_volume_fractions->getDepth() == d_num_species)
    {
        /*
         * Get the pointers to the cell data of volume fractions.
         */
        
        std::vector<double*> Z;
        Z.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            Z.push_back(data_volume_fractions->getPointer(si));
        }
        
        if (d_dim == tbox::Dimension(1))
        {
            /*
             * Get the local lower index, numbers of cells in each dimension and offsets.
             */
            
            const int domain_lo_0 = domain_lo[0];
            const int domain_dim_0 = domain_dims[0];
            
            const int offset_0_bulk_viscosity = offset_bulk_viscosity[0];
            const int offset_0_min = offset_min[0];
            const int offset_0_volume_fractions = offset_volume_fractions[0];
            
            for (int si = 0; si < d_num_species; si++)
            {
                data_temperature_species->copyDepth(0, *data_species_temperatures, si);
                
                getSpeciesMolecularProperties(species_molecular_properties_ptr, si);
                
                d_equation_of_bulk_viscosity->
                    computeBulkViscosity(
                        data_bulk_viscosity_species,
                        data_pressure,
                        data_temperature_species,
                        species_molecular_properties_const_ptr,
                        domain);
                
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_bulk_viscosity = i + offset_0_bulk_viscosity;
                    const int idx_min = i + offset_0_min;
                    const int idx_volume_fractions = i + offset_0_volume_fractions;
                    
                    mu_v[idx_bulk_viscosity] += mu_v_i[idx_min]*Z[si][idx_volume_fractions];
                }
            }
        }
        else if (d_dim == tbox::Dimension(2))
        {
            /*
             * Get the local lower indices, numbers of cells in each dimension and offsets.
             */
            
            const int domain_lo_0 = domain_lo[0];
            const int domain_lo_1 = domain_lo[1];
            const int domain_dim_0 = domain_dims[0];
            const int domain_dim_1 = domain_dims[1];
            
            const int offset_0_bulk_viscosity = offset_bulk_viscosity[0];
            const int offset_1_bulk_viscosity = offset_bulk_viscosity[1];
            const int ghostcell_dim_0_bulk_viscosity = ghostcell_dims_bulk_viscosity[0];
            
            const int offset_0_min = offset_min[0];
            const int offset_1_min = offset_min[1];
            const int ghostcell_dim_0_min = ghostcell_dims_min[0];
            
            const int offset_0_volume_fractions = offset_volume_fractions[0];
            const int offset_1_volume_fractions = offset_volume_fractions[1];
            const int ghostcell_dim_0_volume_fractions = ghostcell_dims_volume_fractions[0];
            
            for (int si = 0; si < d_num_species; si++)
            {
                data_temperature_species->copyDepth(0, *data_species_temperatures, si);
                
                getSpeciesMolecularProperties(species_molecular_properties_ptr, si);
                
                d_equation_of_bulk_viscosity->
                    computeBulkViscosity(
                        data_bulk_viscosity_species,
                        data_pressure,
                        data_temperature_species,
                        species_molecular_properties_const_ptr,
                        domain);
                
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_bulk_viscosity = (i + offset_0_bulk_viscosity) +
                            (j + offset_1_bulk_viscosity)*ghostcell_dim_0_bulk_viscosity;
                        
                        const int idx_min = (i + offset_0_min) +
                            (j + offset_1_min)*ghostcell_dim_0_min;
                        
                        const int idx_volume_fractions = (i + offset_0_volume_fractions) +
                            (j + offset_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        mu_v[idx_bulk_viscosity] += mu_v_i[idx_min]*Z[si][idx_volume_fractions];
                    }
                }
            }
        }
        else if (d_dim == tbox::Dimension(3))
        {
            /*
             * Get the local lower indices, numbers of cells in each dimension and offsets.
             */
            
            const int domain_lo_0 = domain_lo[0];
            const int domain_lo_1 = domain_lo[1];
            const int domain_lo_2 = domain_lo[2];
            const int domain_dim_0 = domain_dims[0];
            const int domain_dim_1 = domain_dims[1];
            const int domain_dim_2 = domain_dims[2];
            
            const int offset_0_bulk_viscosity = offset_bulk_viscosity[0];
            const int offset_1_bulk_viscosity = offset_bulk_viscosity[1];
            const int offset_2_bulk_viscosity = offset_bulk_viscosity[2];
            const int ghostcell_dim_0_bulk_viscosity = ghostcell_dims_bulk_viscosity[0];
            const int ghostcell_dim_1_bulk_viscosity = ghostcell_dims_bulk_viscosity[1];
            
            const int offset_0_min = offset_min[0];
            const int offset_1_min = offset_min[1];
            const int offset_2_min = offset_min[2];
            const int ghostcell_dim_0_min = ghostcell_dims_min[0];
            const int ghostcell_dim_1_min = ghostcell_dims_min[1];
            
            const int offset_0_volume_fractions = offset_volume_fractions[0];
            const int offset_1_volume_fractions = offset_volume_fractions[1];
            const int offset_2_volume_fractions = offset_volume_fractions[2];
            const int ghostcell_dim_0_volume_fractions = ghostcell_dims_volume_fractions[0];
            const int ghostcell_dim_1_volume_fractions = ghostcell_dims_volume_fractions[1];
            
            for (int si = 0; si < d_num_species; si++)
            {
                data_temperature_species->copyDepth(0, *data_species_temperatures, si);
                
                getSpeciesMolecularProperties(species_molecular_properties_ptr, si);
                
                d_equation_of_bulk_viscosity->
                    computeBulkViscosity(
                        data_bulk_viscosity_species,
                        data_pressure,
                        data_temperature_species,
                        species_molecular_properties_const_ptr,
                        domain);
                
                for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                {
                    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_bulk_viscosity = (i + offset_0_bulk_viscosity) +
                                (j + offset_1_bulk_viscosity)*ghostcell_dim_0_bulk_viscosity +
                                (k + offset_2_bulk_viscosity)*ghostcell_dim_0_bulk_viscosity*
                                    ghostcell_dim_1_bulk_viscosity;
                            
                            const int idx_min = (i + offset_0_min) +
                                (j + offset_1_min)*ghostcell_dim_0_min +
                                (k + offset_2_min)*ghostcell_dim_0_min*
                                    ghostcell_dim_1_min;
                            
                            const int idx_volume_fractions = (i + offset_0_volume_fractions) +
                                (j + offset_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + offset_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            mu_v[idx_bulk_viscosity] += mu_v_i[idx_min]*Z[si][idx_volume_fractions];
                        }
                    }
                }
            }
        }
    }
    else if (data_volume_fractions->getDepth() == d_num_species - 1)
    {
        data_volume_fractions_last->fillAll(double(1));
        
        /*
         * Get the pointers to the cell data of volume fractions.
         */
        
        std::vector<double*> Z;
        Z.reserve(d_num_species - 1);
        for (int si = 0; si < d_num_species - 1; si++)
        {
            Z.push_back(data_volume_fractions->getPointer(si));
        }
        
        double* Z_last = data_volume_fractions_last->getPointer(0);
        
        if (d_dim == tbox::Dimension(1))
        {
            /*
             * Get the local lower index, numbers of cells in each dimension and offsets.
             */
            
            const int domain_lo_0 = domain_lo[0];
            const int domain_dim_0 = domain_dims[0];
            
            const int offset_0_bulk_viscosity = offset_bulk_viscosity[0];
            const int offset_0_min = offset_min[0];
            const int offset_0_volume_fractions = offset_volume_fractions[0];
            
            for (int si = 0; si < d_num_species - 1; si++)
            {
                data_temperature_species->copyDepth(0, *data_species_temperatures, si);
                
                getSpeciesMolecularProperties(species_molecular_properties_ptr, si);
                
                d_equation_of_bulk_viscosity->
                    computeBulkViscosity(
                        data_bulk_viscosity_species,
                        data_pressure,
                        data_temperature_species,
                        species_molecular_properties_const_ptr,
                        domain);
                
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_bulk_viscosity = i + offset_0_bulk_viscosity;
                    const int idx_min = i + offset_0_min;
                    const int idx_volume_fractions = i + offset_0_volume_fractions;
                    
                    mu_v[idx_bulk_viscosity] += mu_v_i[idx_min]*Z[si][idx_volume_fractions];
                    
                    // Compute the volume fraction of the last species.
                    Z_last[idx_min] -= Z[si][idx_volume_fractions];
                }
            }
            
            data_temperature_species->copyDepth(0, *data_species_temperatures, d_num_species - 1);
            
            getSpeciesMolecularProperties(species_molecular_properties_ptr, d_num_species - 1);
            
            d_equation_of_bulk_viscosity->
                computeBulkViscosity(
                    data_bulk_viscosity_species,
                    data_pressure,
                    data_temperature_species,
                    species_molecular_properties_const_ptr,
                    domain);
            
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_bulk_viscosity = i + offset_0_bulk_viscosity;
                const int idx_min = i + offset_0_min;
                
                mu_v[idx_bulk_viscosity] += mu_v_i[idx_min]*Z_last[idx_min];
            }
        }
        else if (d_dim == tbox::Dimension(2))
        {
            /*
             * Get the local lower indices, numbers of cells in each dimension and offsets.
             */
            
            const int domain_lo_0 = domain_lo[0];
            const int domain_lo_1 = domain_lo[1];
            const int domain_dim_0 = domain_dims[0];
            const int domain_dim_1 = domain_dims[1];
            
            const int offset_0_bulk_viscosity = offset_bulk_viscosity[0];
            const int offset_1_bulk_viscosity = offset_bulk_viscosity[1];
            const int ghostcell_dim_0_bulk_viscosity = ghostcell_dims_bulk_viscosity[0];
            
            const int offset_0_min = offset_min[0];
            const int offset_1_min = offset_min[1];
            const int ghostcell_dim_0_min = ghostcell_dims_min[0];
            
            const int offset_0_volume_fractions = offset_volume_fractions[0];
            const int offset_1_volume_fractions = offset_volume_fractions[1];
            const int ghostcell_dim_0_volume_fractions = ghostcell_dims_volume_fractions[0];
            
            for (int si = 0; si < d_num_species - 1; si++)
            {
                data_temperature_species->copyDepth(0, *data_species_temperatures, si);
                
                getSpeciesMolecularProperties(species_molecular_properties_ptr, si);
                
                d_equation_of_bulk_viscosity->
                    computeBulkViscosity(
                        data_bulk_viscosity_species,
                        data_pressure,
                        data_temperature_species,
                        species_molecular_properties_const_ptr,
                        domain);
                
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_bulk_viscosity = (i + offset_0_bulk_viscosity) +
                            (j + offset_1_bulk_viscosity)*ghostcell_dim_0_bulk_viscosity;
                        
                        const int idx_min = (i + offset_0_min) +
                            (j + offset_1_min)*ghostcell_dim_0_min;
                        
                        const int idx_volume_fractions = (i + offset_0_volume_fractions) +
                            (j + offset_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        mu_v[idx_bulk_viscosity] += mu_v_i[idx_min]*Z[si][idx_volume_fractions];
                        
                        // Compute the volume fraction of the last species.
                        Z_last[idx_min] -= Z[si][idx_volume_fractions];
                    }
                }
            }
            
            data_temperature_species->copyDepth(0, *data_species_temperatures, d_num_species - 1);
            
            getSpeciesMolecularProperties(species_molecular_properties_ptr, d_num_species - 1);
            
            d_equation_of_bulk_viscosity->
                computeBulkViscosity(
                    data_bulk_viscosity_species,
                    data_pressure,
                    data_temperature_species,
                    species_molecular_properties_const_ptr,
                    domain);
            
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_bulk_viscosity = (i + offset_0_bulk_viscosity) +
                        (j + offset_1_bulk_viscosity)*ghostcell_dim_0_bulk_viscosity;
                    
                    const int idx_min = (i + offset_0_min) +
                        (j + offset_1_min)*ghostcell_dim_0_min;
                    
                    mu_v[idx_bulk_viscosity] += mu_v_i[idx_min]*Z_last[idx_min];
                }
            }
        }
        else if (d_dim == tbox::Dimension(3))
        {
            /*
             * Get the local lower indices, numbers of cells in each dimension and offsets.
             */
            
            const int domain_lo_0 = domain_lo[0];
            const int domain_lo_1 = domain_lo[1];
            const int domain_lo_2 = domain_lo[2];
            const int domain_dim_0 = domain_dims[0];
            const int domain_dim_1 = domain_dims[1];
            const int domain_dim_2 = domain_dims[2];
            
            const int offset_0_bulk_viscosity = offset_bulk_viscosity[0];
            const int offset_1_bulk_viscosity = offset_bulk_viscosity[1];
            const int offset_2_bulk_viscosity = offset_bulk_viscosity[2];
            const int ghostcell_dim_0_bulk_viscosity = ghostcell_dims_bulk_viscosity[0];
            const int ghostcell_dim_1_bulk_viscosity = ghostcell_dims_bulk_viscosity[1];
            
            const int offset_0_min = offset_min[0];
            const int offset_1_min = offset_min[1];
            const int offset_2_min = offset_min[2];
            const int ghostcell_dim_0_min = ghostcell_dims_min[0];
            const int ghostcell_dim_1_min = ghostcell_dims_min[1];
            
            const int offset_0_volume_fractions = offset_volume_fractions[0];
            const int offset_1_volume_fractions = offset_volume_fractions[1];
            const int offset_2_volume_fractions = offset_volume_fractions[2];
            const int ghostcell_dim_0_volume_fractions = ghostcell_dims_volume_fractions[0];
            const int ghostcell_dim_1_volume_fractions = ghostcell_dims_volume_fractions[1];
            
            for (int si = 0; si < d_num_species - 1; si++)
            {
                data_temperature_species->copyDepth(0, *data_species_temperatures, si);
                
                getSpeciesMolecularProperties(species_molecular_properties_ptr, si);
                
                d_equation_of_bulk_viscosity->
                    computeBulkViscosity(
                        data_bulk_viscosity_species,
                        data_pressure,
                        data_temperature_species,
                        species_molecular_properties_const_ptr,
                        domain);
                
                for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                {
                    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_bulk_viscosity = (i + offset_0_bulk_viscosity) +
                                (j + offset_1_bulk_viscosity)*ghostcell_dim_0_bulk_viscosity +
                                (k + offset_2_bulk_viscosity)*ghostcell_dim_0_bulk_viscosity*
                                    ghostcell_dim_1_bulk_viscosity;
                            
                            const int idx_min = (i + offset_0_min) +
                                (j + offset_1_min)*ghostcell_dim_0_min +
                                (k + offset_2_min)*ghostcell_dim_0_min*
                                    ghostcell_dim_1_min;
                            
                            const int idx_volume_fractions = (i + offset_0_volume_fractions) +
                                (j + offset_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + offset_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            mu_v[idx_bulk_viscosity] += mu_v_i[idx_min]*Z[si][idx_volume_fractions];
                            
                            // Compute the volume fraction of the last species.
                            Z_last[idx_min] -= Z[si][idx_volume_fractions];
                        }
                    }
                }
            }
            
            data_temperature_species->copyDepth(0, *data_species_temperatures, d_num_species - 1);
            
            getSpeciesMolecularProperties(species_molecular_properties_ptr, d_num_species - 1);
            
            d_equation_of_bulk_viscosity->
                computeBulkViscosity(
                    data_bulk_viscosity_species,
                    data_pressure,
                    data_temperature_species,
                    species_molecular_properties_const_ptr,
                    domain);
            
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_bulk_viscosity = (i + offset_0_bulk_viscosity) +
                            (j + offset_1_bulk_viscosity)*ghostcell_dim_0_bulk_viscosity +
                            (k + offset_2_bulk_viscosity)*ghostcell_dim_0_bulk_viscosity*
                                ghostcell_dim_1_bulk_viscosity;
                        
                        const int idx_min = (i + offset_0_min) +
                            (j + offset_1_min)*ghostcell_dim_0_min +
                            (k + offset_2_min)*ghostcell_dim_0_min*
                                ghostcell_dim_1_min;
                        
                        mu_v[idx_bulk_viscosity] += mu_v_i[idx_min]*Z_last[idx_min];
                    }
                }
            }
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Number of components in the data of volume fractions provided is not"
            << " equal to the total number of species or (total number of species - 1)."
            << std::endl);
    }
}


/*
 * Get the molecular properties of a species.
 */
void
EquationOfBulkViscosityMixingRulesCramer::getSpeciesMolecularProperties(
    std::vector<double*>& species_molecular_properties,
    const int species_index) const
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
