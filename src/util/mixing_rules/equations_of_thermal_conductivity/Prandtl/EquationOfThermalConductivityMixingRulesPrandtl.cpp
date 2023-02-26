#include "util/mixing_rules/equations_of_thermal_conductivity/Prandtl/EquationOfThermalConductivityMixingRulesPrandtl.hpp"

EquationOfThermalConductivityMixingRulesPrandtl::EquationOfThermalConductivityMixingRulesPrandtl(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const int& num_species,
    const MIXING_CLOSURE_MODEL::TYPE& mixing_closure_model,
    const HAMERS_SHARED_PTR<tbox::Database>& equation_of_thermal_conductivity_mixing_rules_db):
        EquationOfThermalConductivityMixingRules(
            object_name,
            dim,
            num_species,
            mixing_closure_model,
            equation_of_thermal_conductivity_mixing_rules_db)
{
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
                equation_of_thermal_conductivity_mixing_rules_db->getRealVector("species_c_p");
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
                equation_of_thermal_conductivity_mixing_rules_db->getRealVector("d_species_c_p");
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
                equation_of_thermal_conductivity_mixing_rules_db->getRealVector("species_Pr");
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
                equation_of_thermal_conductivity_mixing_rules_db->getRealVector("d_species_Pr");
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
                equation_of_thermal_conductivity_mixing_rules_db->getRealVector("species_M");
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
                equation_of_thermal_conductivity_mixing_rules_db->getRealVector("d_species_M");
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
    
    d_equation_of_thermal_conductivity.reset(new EquationOfThermalConductivityPrandtl(
        "d_equation_of_thermal_conductivity",
        dim,
        d_equation_of_shear_viscosity_mixing_rules->
            getEquationOfShearViscosity()));
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
    const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const
{
    restart_db->putString("d_equation_of_shear_viscosity_str", d_equation_of_shear_viscosity_str);
    d_equation_of_shear_viscosity_mixing_rules->putToRestart(restart_db);
    
    if (!restart_db->keyExists("d_species_c_p"))
    {
        restart_db->putRealVector("d_species_c_p", d_species_c_p);
    }
    
    if (!restart_db->keyExists("d_species_Pr"))
    {
        restart_db->putRealVector("d_species_Pr", d_species_Pr);
    }
    
    if (!restart_db->keyExists("d_species_M"))
    {
        restart_db->putRealVector("d_species_M", d_species_M);
    }
}


/*
 * Compute the thermal conductivity of the mixture with isothermal and isobaric equilibrium assumptions.
 */
Real
EquationOfThermalConductivityMixingRulesPrandtl::getThermalConductivity(
    const Real* const pressure,
    const Real* const temperature,
    const std::vector<const Real*>& mass_fractions) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT((d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC) ||
                (d_mixing_closure_model == MIXING_CLOSURE_MODEL::NO_MODEL && d_num_species == 1));
    TBOX_ASSERT((static_cast<int>(mass_fractions.size()) == d_num_species) ||
                (static_cast<int>(mass_fractions.size()) == d_num_species - 1));
#endif
    
    Real kappa = Real(0);
    
    Real num = Real(0);
    Real den = Real(0);
    
    /*
     * Initialize the container and pointers to the container for the molecular properties
     * of a species.
     */
    
    std::vector<Real> species_molecular_properties;
    std::vector<Real*> species_molecular_properties_ptr;
    std::vector<const Real*> species_molecular_properties_const_ptr;
    
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
            
            const Real kappa_i = d_equation_of_thermal_conductivity->
                getThermalConductivity(
                    pressure,
                    temperature,
                    species_molecular_properties_const_ptr);
            
            const Real weight = *(mass_fractions[si])/(std::sqrt(species_molecular_properties[2]));
            
            num += kappa_i*weight;
            den += weight;
        }
    }
    else if (static_cast<int>(mass_fractions.size()) == d_num_species - 1)
    {
        Real Y_last = Real(1);
        
        for (int si = 0; si < d_num_species - 1; si++)
        {
            getSpeciesMolecularProperties(species_molecular_properties_ptr, si);
            
            const Real kappa_i = d_equation_of_thermal_conductivity->
                getThermalConductivity(
                    pressure,
                    temperature,
                    species_molecular_properties_const_ptr);
            
            const Real weight = *(mass_fractions[si])/(std::sqrt(species_molecular_properties[2]));
            
            num += kappa_i*weight;
            den += weight;
            
            // Compute the mass fraction of the last species.
            Y_last -= *(mass_fractions[si]);
        }
        
        /*
         * Add the contribution from the last species.
         */
        getSpeciesMolecularProperties(species_molecular_properties_ptr, d_num_species - 1);
        
        const Real kappa_last = d_equation_of_thermal_conductivity->
            getThermalConductivity(
                pressure,
                temperature,
                species_molecular_properties_const_ptr);
        
        const Real weight = Y_last/(std::sqrt(species_molecular_properties[2]));
        
        num += kappa_last*weight;
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
    
    kappa = num/den;
    
    return kappa;
}


/*
 * Compute the thermal conductivity of the mixture with isothermal and isobaric equilibrium assumptions.
 */
void
EquationOfThermalConductivityMixingRulesPrandtl::computeThermalConductivity(
    HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_thermal_conductivity,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_temperature,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_mass_fractions,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT((d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC) ||
                (d_mixing_closure_model == MIXING_CLOSURE_MODEL::NO_MODEL && d_num_species == 1));
    
    TBOX_ASSERT(data_thermal_conductivity);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_temperature);
    TBOX_ASSERT(data_mass_fractions);
    
    TBOX_ASSERT((data_mass_fractions->getDepth() == d_num_species) ||
                (data_mass_fractions->getDepth() == d_num_species - 1));
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_thermal_conductivity = data_thermal_conductivity->getGhostBox();
    const hier::IntVector ghostcell_dims_thermal_conductivity = ghost_box_thermal_conductivity.numberCells();
    
    const hier::Box ghost_box_mass_fractions = data_mass_fractions->getGhostBox();
    const hier::IntVector ghostcell_dims_mass_fractions = ghost_box_mass_fractions.numberCells();
    
    // Delcare data containers for thermal conductivity of a species, denominator and numerator.
    HAMERS_SHARED_PTR<pdat::CellData<Real> > data_thermal_conductivity_species;;
    HAMERS_SHARED_PTR<pdat::CellData<Real> > data_den;
    HAMERS_SHARED_PTR<pdat::CellData<Real> > data_num;
    
    // Declare data container for last mass fraction.
    HAMERS_SHARED_PTR<pdat::CellData<Real> > data_mass_fractions_last;
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets of all data and dimensions of the ghost cell box for denominator,
     * numerator and last mass fraction and allocate memory.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_thermal_conductivity(d_dim);
    hier::IntVector offset_mass_fractions(d_dim);
    hier::IntVector offset_min(d_dim);
    
    hier::IntVector ghostcell_dims_min(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_thermal_conductivity = data_thermal_conductivity->getGhostCellWidth();
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
        const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
        
        // Get the dimensions of box that covers the interior of patch.
        const hier::Box interior_box = data_thermal_conductivity->getBox();
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
        
        num_ghosts_min = num_ghosts_thermal_conductivity;
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_temperature, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_mass_fractions, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_min = num_ghosts_min;
        offset_thermal_conductivity = num_ghosts_thermal_conductivity;
        offset_mass_fractions = num_ghosts_mass_fractions;
        
        ghostcell_dims_min = interior_dims + num_ghosts_min*2;
        
        data_thermal_conductivity_species = HAMERS_MAKE_SHARED<pdat::CellData<Real> >(interior_box, 1, num_ghosts_min);
        data_den = HAMERS_MAKE_SHARED<pdat::CellData<Real> >(interior_box, 1, num_ghosts_min);
        data_num = HAMERS_MAKE_SHARED<pdat::CellData<Real> >(interior_box, 1, num_ghosts_min);
        
        if (data_mass_fractions->getDepth() == d_num_species - 1)
        {
            data_mass_fractions_last = HAMERS_MAKE_SHARED<pdat::CellData<Real> >(interior_box, 1, num_ghosts_min);
        }
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
        TBOX_ASSERT(data_thermal_conductivity->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_temperature->getGhostBox().contains(domain));
        TBOX_ASSERT(data_mass_fractions->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_min = hier::IntVector::getZero(d_dim);
        offset_thermal_conductivity = domain.lower() - ghost_box_thermal_conductivity.lower();
        offset_mass_fractions = domain.lower() - ghost_box_mass_fractions.lower();
        
        ghostcell_dims_min = domain_dims;
        
        data_thermal_conductivity_species =
            HAMERS_MAKE_SHARED<pdat::CellData<Real> >(domain, 1, hier::IntVector::getZero(d_dim));
        data_den = HAMERS_MAKE_SHARED<pdat::CellData<Real> >(domain, 1, hier::IntVector::getZero(d_dim));
        data_num = HAMERS_MAKE_SHARED<pdat::CellData<Real> >(domain, 1, hier::IntVector::getZero(d_dim));
        
        if (data_mass_fractions->getDepth() == d_num_species - 1)
        {
            data_mass_fractions_last =
                HAMERS_MAKE_SHARED<pdat::CellData<Real> >(domain, 1, hier::IntVector::getZero(d_dim));
        }
    }
    
    // Declare data containers for species molecular properties.
    std::vector<Real> species_molecular_properties;
    std::vector<Real*> species_molecular_properties_ptr;
    std::vector<const Real*> species_molecular_properties_const_ptr;
    
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
     * Get the pointers to the cell data of mixture thermal conductivity, species thermal conductivity,
     * denominator and numerator.
     */
    
    Real* kappa = data_thermal_conductivity->getPointer(0);
    Real* kappa_i = data_thermal_conductivity_species->getPointer(0);
    Real* den = data_den->getPointer(0);
    Real* num = data_num->getPointer(0);
    
    /*
     * Fill zeros for denominator and numerator.
     */
    
    data_den->fillAll(Real(0));
    data_num->fillAll(Real(0));
    
    if (data_mass_fractions->getDepth() == d_num_species)
    {
        /*
         * Get the pointers to the cell data of mass fractions.
         */
        
        std::vector<Real*> Y;
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
            
            const int offset_0_thermal_conductivity = offset_thermal_conductivity[0];
            const int offset_0_min = offset_min[0];
            const int offset_0_mass_fractions = offset_mass_fractions[0];
            
            // Compute the mixture thermal conductivity field.
            for (int si = 0; si < d_num_species; si++)
            {
                getSpeciesMolecularProperties(species_molecular_properties_ptr, si);
                
                d_equation_of_thermal_conductivity->
                    computeThermalConductivity(
                        data_thermal_conductivity_species,
                        data_pressure,
                        data_temperature,
                        species_molecular_properties_const_ptr,
                        domain);
                
                const Real factor = Real(1)/(std::sqrt(species_molecular_properties[2]));
                
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_min = i + offset_0_min;
                    const int idx_mass_fractions = i + offset_0_mass_fractions;
                    
                    const Real weight = Y[si][idx_mass_fractions]*factor;
                    
                    num[idx_min] += kappa_i[idx_min]*weight;
                    den[idx_min] += weight;
                }
            }
            
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_thermal_conductivity = i + offset_0_thermal_conductivity;
                const int idx_min = i + offset_0_min;
                
                kappa[idx_thermal_conductivity] = num[idx_min]/den[idx_min];
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
            
            const int offset_0_thermal_conductivity = offset_thermal_conductivity[0];
            const int offset_1_thermal_conductivity = offset_thermal_conductivity[1];
            const int ghostcell_dim_0_thermal_conductivity = ghostcell_dims_thermal_conductivity[0];
            
            const int offset_0_min = offset_min[0];
            const int offset_1_min = offset_min[1];
            const int ghostcell_dim_0_min = ghostcell_dims_min[0];
            
            const int offset_0_mass_fractions = offset_mass_fractions[0];
            const int offset_1_mass_fractions = offset_mass_fractions[1];
            const int ghostcell_dim_0_mass_fractions = ghostcell_dims_mass_fractions[0];
            
            // Compute the mixture thermal conductivity field.
            for (int si = 0; si < d_num_species; si++)
            {
                getSpeciesMolecularProperties(species_molecular_properties_ptr, si);
                
                d_equation_of_thermal_conductivity->
                    computeThermalConductivity(
                        data_thermal_conductivity_species,
                        data_pressure,
                        data_temperature,
                        species_molecular_properties_const_ptr,
                        domain);
                
                const Real factor = Real(1)/(std::sqrt(species_molecular_properties[2]));
                
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_min = (i + offset_0_min) +
                            (j + offset_1_min)*ghostcell_dim_0_min;
                        
                        const int idx_mass_fractions = (i + offset_0_mass_fractions) +
                            (j + offset_1_mass_fractions)*ghostcell_dim_0_mass_fractions;
                        
                        const Real weight = Y[si][idx_mass_fractions]*factor;
                        
                        num[idx_min] += kappa_i[idx_min]*weight;
                        den[idx_min] += weight;
                    }
                }
            }
            
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_thermal_conductivity = (i + offset_0_thermal_conductivity) +
                        (j + offset_1_thermal_conductivity)*ghostcell_dim_0_thermal_conductivity;
                    
                    const int idx_min = (i + offset_0_min) +
                        (j + offset_1_min)*ghostcell_dim_0_min;
                    
                    kappa[idx_thermal_conductivity] = num[idx_min]/den[idx_min];
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
            
            const int offset_0_thermal_conductivity = offset_thermal_conductivity[0];
            const int offset_1_thermal_conductivity = offset_thermal_conductivity[1];
            const int offset_2_thermal_conductivity = offset_thermal_conductivity[2];
            const int ghostcell_dim_0_thermal_conductivity = ghostcell_dims_thermal_conductivity[0];
            const int ghostcell_dim_1_thermal_conductivity = ghostcell_dims_thermal_conductivity[1];
            
            const int offset_0_min = offset_min[0];
            const int offset_1_min = offset_min[1];
            const int offset_2_min = offset_min[2];
            const int ghostcell_dim_0_min = ghostcell_dims_min[0];
            const int ghostcell_dim_1_min = ghostcell_dims_min[1];
            
            const int offset_0_mass_fractions = offset_mass_fractions[0];
            const int offset_1_mass_fractions = offset_mass_fractions[1];
            const int offset_2_mass_fractions = offset_mass_fractions[2];
            const int ghostcell_dim_0_mass_fractions = ghostcell_dims_mass_fractions[0];
            const int ghostcell_dim_1_mass_fractions = ghostcell_dims_mass_fractions[1];
            
            // Compute the mixture thermal conductivity field.
            for (int si = 0; si < d_num_species; si++)
            {
                getSpeciesMolecularProperties(species_molecular_properties_ptr, si);
                
                d_equation_of_thermal_conductivity->
                    computeThermalConductivity(
                        data_thermal_conductivity_species,
                        data_pressure,
                        data_temperature,
                        species_molecular_properties_const_ptr,
                        domain);
                
                const Real factor = Real(1)/(std::sqrt(species_molecular_properties[2]));
                
                for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                {
                    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_min = (i + offset_0_min) +
                                (j + offset_1_min)*ghostcell_dim_0_min +
                                (k + offset_2_min)*ghostcell_dim_0_min*
                                    ghostcell_dim_1_min;
                            
                            const int idx_mass_fractions = (i + offset_0_mass_fractions) +
                                (j + offset_1_mass_fractions)*ghostcell_dim_0_mass_fractions +
                                (k + offset_2_mass_fractions)*ghostcell_dim_0_mass_fractions*
                                    ghostcell_dim_1_mass_fractions;
                            
                            const Real weight = Y[si][idx_mass_fractions]*factor;
                            
                            num[idx_min] += kappa_i[idx_min]*weight;
                            den[idx_min] += weight;
                        }
                    }
                }
            }
            
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_thermal_conductivity = (i + offset_0_thermal_conductivity) +
                            (j + offset_1_thermal_conductivity)*ghostcell_dim_0_thermal_conductivity +
                            (k + offset_2_thermal_conductivity)*ghostcell_dim_0_thermal_conductivity*
                                ghostcell_dim_1_thermal_conductivity;
                        
                        const int idx_min = (i + offset_0_min) +
                            (j + offset_1_min)*ghostcell_dim_0_min +
                            (k + offset_2_min)*ghostcell_dim_0_min*
                                ghostcell_dim_1_min;
                        
                        kappa[idx_thermal_conductivity] = num[idx_min]/den[idx_min];
                    }
                }
            }
        }
    }
    else if (data_mass_fractions->getDepth() == d_num_species - 1)
    {
        data_mass_fractions_last->fillAll(Real(1));
        
        /*
         * Get the pointers to the cell data of mass fractions.
         */
        
        std::vector<Real*> Y;
        Y.reserve(d_num_species - 1);
        for (int si = 0; si < d_num_species - 1; si++)
        {
            Y.push_back(data_mass_fractions->getPointer(si));
        }
        
        Real* Y_last = data_mass_fractions_last->getPointer(0);
        
        if (d_dim == tbox::Dimension(1))
        {
            /*
             * Get the local lower index, numbers of cells in each dimension and offsets.
             */
            
            const int domain_lo_0 = domain_lo[0];
            const int domain_dim_0 = domain_dims[0];
            
            const int offset_0_thermal_conductivity = offset_thermal_conductivity[0];
            const int offset_0_min = offset_min[0];
            const int offset_0_mass_fractions = offset_mass_fractions[0];
            
            // Compute the mixture thermal conductivity field.
            for (int si = 0; si < d_num_species - 1; si++)
            {
                getSpeciesMolecularProperties(species_molecular_properties_ptr, si);
                
                d_equation_of_thermal_conductivity->
                    computeThermalConductivity(
                        data_thermal_conductivity_species,
                        data_pressure,
                        data_temperature,
                        species_molecular_properties_const_ptr,
                        domain);
                
                const Real factor = Real(1)/(std::sqrt(species_molecular_properties[2]));
                
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_min = i + offset_0_min;
                    const int idx_mass_fractions = i + offset_0_mass_fractions;
                    
                    const Real weight = Y[si][idx_mass_fractions]*factor;
                    
                    num[idx_min] += kappa_i[idx_min]*weight;
                    den[idx_min] += weight;
                    
                    // Compute the mass fraction of the last species.
                    Y_last[idx_min] -= Y[si][idx_mass_fractions];
                }
            }
            
            getSpeciesMolecularProperties(species_molecular_properties_ptr, d_num_species - 1);
            
            d_equation_of_thermal_conductivity->
                computeThermalConductivity(
                    data_thermal_conductivity_species,
                    data_pressure,
                    data_temperature,
                    species_molecular_properties_const_ptr,
                    domain);
            
            const Real factor = Real(1)/(std::sqrt(species_molecular_properties[2]));
            
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_thermal_conductivity = i + offset_0_thermal_conductivity;
                const int idx_min = i + offset_0_min;
                
                const Real weight = Y_last[idx_min]*factor;
                
                num[idx_min] += kappa_i[idx_min]*weight;
                den[idx_min] += weight;
                
                kappa[idx_thermal_conductivity] = num[idx_min]/den[idx_min];
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
            
            const int offset_0_thermal_conductivity = offset_thermal_conductivity[0];
            const int offset_1_thermal_conductivity = offset_thermal_conductivity[1];
            const int ghostcell_dim_0_thermal_conductivity = ghostcell_dims_thermal_conductivity[0];
            
            const int offset_0_min = offset_min[0];
            const int offset_1_min = offset_min[1];
            const int ghostcell_dim_0_min = ghostcell_dims_min[0];
            
            const int offset_0_mass_fractions = offset_mass_fractions[0];
            const int offset_1_mass_fractions = offset_mass_fractions[1];
            const int ghostcell_dim_0_mass_fractions = ghostcell_dims_mass_fractions[0];
            
            // Compute the mixture thermal conductivity field.
            for (int si = 0; si < d_num_species - 1; si++)
            {
                getSpeciesMolecularProperties(species_molecular_properties_ptr, si);
                
                d_equation_of_thermal_conductivity->
                    computeThermalConductivity(
                        data_thermal_conductivity_species,
                        data_pressure,
                        data_temperature,
                        species_molecular_properties_const_ptr,
                        domain);
                
                const Real factor = Real(1)/(std::sqrt(species_molecular_properties[2]));
                
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_min = (i + offset_0_min) +
                            (j + offset_1_min)*ghostcell_dim_0_min;
                        
                        const int idx_mass_fractions = (i + offset_0_mass_fractions) +
                            (j + offset_1_mass_fractions)*ghostcell_dim_0_mass_fractions;
                        
                        const Real weight = Y[si][idx_mass_fractions]*factor;
                        
                        num[idx_min] += kappa_i[idx_min]*weight;
                        den[idx_min] += weight;
                        
                        // Compute the mass fraction of the last species.
                        Y_last[idx_min] -= Y[si][idx_mass_fractions];
                    }
                }
            }
            
            getSpeciesMolecularProperties(species_molecular_properties_ptr, d_num_species - 1);
            
            d_equation_of_thermal_conductivity->
                computeThermalConductivity(
                    data_thermal_conductivity_species,
                    data_pressure,
                    data_temperature,
                    species_molecular_properties_const_ptr,
                    domain);
            
            const Real factor = Real(1)/(std::sqrt(species_molecular_properties[2]));
            
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_thermal_conductivity = (i + offset_0_thermal_conductivity) +
                        (j + offset_1_thermal_conductivity)*ghostcell_dim_0_thermal_conductivity;
                    
                    const int idx_min = (i + offset_0_min) +
                        (j + offset_1_min)*ghostcell_dim_0_min;
                    
                    const Real weight = Y_last[idx_min]*factor;
                    
                    num[idx_min] += kappa_i[idx_min]*weight;
                    den[idx_min] += weight;
                    
                    kappa[idx_thermal_conductivity] = num[idx_min]/den[idx_min];
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
            
            const int offset_0_thermal_conductivity = offset_thermal_conductivity[0];
            const int offset_1_thermal_conductivity = offset_thermal_conductivity[1];
            const int offset_2_thermal_conductivity = offset_thermal_conductivity[2];
            const int ghostcell_dim_0_thermal_conductivity = ghostcell_dims_thermal_conductivity[0];
            const int ghostcell_dim_1_thermal_conductivity = ghostcell_dims_thermal_conductivity[1];
            
            const int offset_0_min = offset_min[0];
            const int offset_1_min = offset_min[1];
            const int offset_2_min = offset_min[2];
            const int ghostcell_dim_0_min = ghostcell_dims_min[0];
            const int ghostcell_dim_1_min = ghostcell_dims_min[1];
            
            const int offset_0_mass_fractions = offset_mass_fractions[0];
            const int offset_1_mass_fractions = offset_mass_fractions[1];
            const int offset_2_mass_fractions = offset_mass_fractions[2];
            const int ghostcell_dim_0_mass_fractions = ghostcell_dims_mass_fractions[0];
            const int ghostcell_dim_1_mass_fractions = ghostcell_dims_mass_fractions[1];
            
            // Compute the mixture thermal conductivity field.
            for (int si = 0; si < d_num_species - 1; si++)
            {
                getSpeciesMolecularProperties(species_molecular_properties_ptr, si);
                
                d_equation_of_thermal_conductivity->
                    computeThermalConductivity(
                        data_thermal_conductivity_species,
                        data_pressure,
                        data_temperature,
                        species_molecular_properties_const_ptr,
                        domain);
                
                const Real factor = Real(1)/(std::sqrt(species_molecular_properties[2]));
                
                for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                {
                    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_min = (i + offset_0_min) +
                                (j + offset_1_min)*ghostcell_dim_0_min +
                                (k + offset_2_min)*ghostcell_dim_0_min*
                                    ghostcell_dim_1_min;
                            
                            const int idx_mass_fractions = (i + offset_0_mass_fractions) +
                                (j + offset_1_mass_fractions)*ghostcell_dim_0_mass_fractions +
                                (k + offset_2_mass_fractions)*ghostcell_dim_0_mass_fractions*
                                    ghostcell_dim_1_mass_fractions;
                            
                            const Real weight = Y[si][idx_mass_fractions]*factor;
                            
                            num[idx_min] += kappa_i[idx_min]*weight;
                            den[idx_min] += weight;
                            
                            // Compute the mass fraction of the last species.
                            Y_last[idx_min] -= Y[si][idx_mass_fractions];
                        }
                    }
                }
            }
            
            getSpeciesMolecularProperties(species_molecular_properties_ptr, d_num_species - 1);
            
            d_equation_of_thermal_conductivity->
                computeThermalConductivity(
                    data_thermal_conductivity_species,
                    data_pressure,
                    data_temperature,
                    species_molecular_properties_const_ptr,
                    domain);
            
            const Real factor = Real(1)/(std::sqrt(species_molecular_properties[2]));
            
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_thermal_conductivity = (i + offset_0_thermal_conductivity) +
                            (j + offset_1_thermal_conductivity)*ghostcell_dim_0_thermal_conductivity +
                            (k + offset_2_thermal_conductivity)*ghostcell_dim_0_thermal_conductivity*
                                ghostcell_dim_1_thermal_conductivity;
                        
                        const int idx_min = (i + offset_0_min) +
                            (j + offset_1_min)*ghostcell_dim_0_min +
                            (k + offset_2_min)*ghostcell_dim_0_min*
                                ghostcell_dim_1_min;
                        
                        const Real weight = Y_last[idx_min]*factor;
                        
                        num[idx_min] += kappa_i[idx_min]*weight;
                        den[idx_min] += weight;
                        
                        kappa[idx_thermal_conductivity] = num[idx_min]/den[idx_min];
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
 * Get the molecular properties of a species.
 */
void
EquationOfThermalConductivityMixingRulesPrandtl::getSpeciesMolecularProperties(
    std::vector<Real*>& species_molecular_properties,
    const int species_index) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(species_molecular_properties.size()) >=
                3 + d_equation_of_shear_viscosity_mixing_rules->
                    getNumberOfSpeciesMolecularProperties());
    TBOX_ASSERT(species_index >= 0);
    TBOX_ASSERT(species_index < d_num_species);
#endif
    
    *(species_molecular_properties[0]) = d_species_c_p[species_index];
    *(species_molecular_properties[1]) = d_species_Pr[species_index];
    *(species_molecular_properties[2]) = d_species_M[species_index];
    
    /*
     * Add the molecular properties of the species in the shear viscosity mixing rule
     * object.
     */
    
    std::vector<Real> mu_molecular_properties;
    std::vector<Real*> mu_molecular_properties_ptr;
    
    int num_mu_molecular_properties = d_equation_of_shear_viscosity_mixing_rules->
        getNumberOfSpeciesMolecularProperties();
    
    mu_molecular_properties.resize(num_mu_molecular_properties);
    mu_molecular_properties_ptr.reserve(num_mu_molecular_properties);
    
    for (int mi = 0; mi < num_mu_molecular_properties; mi++)
    {
        mu_molecular_properties_ptr.push_back(&mu_molecular_properties[mi]);
    }
    
    d_equation_of_shear_viscosity_mixing_rules->getSpeciesMolecularProperties(
        mu_molecular_properties_ptr,
        species_index);
    
    for (int mi = 0; mi < num_mu_molecular_properties; mi++)
    {
        *(species_molecular_properties[3 + mi]) = mu_molecular_properties[mi];
    }
}
