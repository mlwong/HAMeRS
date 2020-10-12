#include "util/mixing_rules/equations_of_state/ideal_gas/EquationOfStateMixingRulesIdealGas.hpp"

EquationOfStateMixingRulesIdealGas::EquationOfStateMixingRulesIdealGas(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const int& num_species,
    const MIXING_CLOSURE_MODEL::TYPE& mixing_closure_model,
    const HAMERS_SHARED_PTR<tbox::Database>& equation_of_state_mixing_rules_db):
        EquationOfStateMixingRules(
            object_name,
            dim,
            num_species,
            mixing_closure_model,
            equation_of_state_mixing_rules_db)
{
    d_equation_of_state.reset(new EquationOfStateIdealGas(
        "d_equation_of_state",
        dim));
    
    /*
     * Get the ratio of specific heats of each species from the database.
     */
    
    if (equation_of_state_mixing_rules_db->keyExists("species_gamma"))
    {
        size_t species_gamma_array_size = equation_of_state_mixing_rules_db->getArraySize("species_gamma");
        if (static_cast<int>(species_gamma_array_size) == d_num_species)
        {
            d_species_gamma = equation_of_state_mixing_rules_db->getDoubleVector("species_gamma");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'species_gamma' entries must be equal to 'num_species'."
                << std::endl);
        }
    }
    else if (equation_of_state_mixing_rules_db->keyExists("d_species_gamma"))
    {
        size_t species_gamma_array_size = equation_of_state_mixing_rules_db->getArraySize("d_species_gamma");
        if (static_cast<int>(species_gamma_array_size) == d_num_species)
        {
            d_species_gamma = equation_of_state_mixing_rules_db->getDoubleVector("d_species_gamma");
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
            << "not found in data for Equation_of_state_mixing_rules."
            << std::endl);
    }
    
    /*
     * Get the gas constant of each species from the database.
     */
    
    if (equation_of_state_mixing_rules_db->keyExists("species_R"))
    {
        size_t species_R_array_size = equation_of_state_mixing_rules_db->getArraySize("species_R");
        if (static_cast<int>(species_R_array_size) == d_num_species)
        {
            d_species_R = equation_of_state_mixing_rules_db->getDoubleVector("species_R");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'species_R' entries must be equal to 'num_species'."
                << std::endl);
        }
    }
    else if (equation_of_state_mixing_rules_db->keyExists("d_species_R"))
    {
        size_t species_R_array_size = equation_of_state_mixing_rules_db->getArraySize("d_species_R");
        if (static_cast<int>(species_R_array_size) == d_num_species)
        {
            d_species_R = equation_of_state_mixing_rules_db->getDoubleVector("d_species_R");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'd_species_R' entries must be equal to 'd_num_species'."
                << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Key data 'species_R'/'d_species_R'"
            << "not found in data for Equation_of_state_mixing_rules."
            << std::endl);
    }
    
    /*
     * Compute the specific heats of each species.
     */
    
    d_species_c_p.reserve(d_num_species);
    for (int si = 0; si < d_num_species; si++)
    {
        d_species_c_p.push_back(d_species_gamma[si]/(d_species_gamma[si] - double(1))*d_species_R[si]);
    }
    
    d_species_c_v.reserve(d_num_species);
    for (int si = 0; si < d_num_species; si++)
    {
        d_species_c_v.push_back(double(1)/(d_species_gamma[si] - double(1))*d_species_R[si]);
    }
    
    /*
     * Compute the molecular weight of each species.
     */
    
    d_species_M.reserve(d_num_species);
    for (int si = 0; si < d_num_species; si++)
    {
        d_species_M.push_back(double(8.3144598)/d_species_R[si]);
    }
}


/*
 * Print all characteristics of the equation of state class.
 */
void
EquationOfStateMixingRulesIdealGas::printClassData(
    std::ostream& os) const
{
    os << "\nPrint EquationOfStateMixingRulesIdealGas object..."
       << std::endl;
    
    os << std::endl;
    os << "EquationOfStateMixingRulesIdealGas: this = "
       << (EquationOfStateMixingRulesIdealGas *)this
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
     * Print the gas constant of each species.
     */
    
    os << "d_species_R = ";
    for (int si = 0; si < d_num_species - 1; si++)
    {
        os << d_species_R[si] << ", ";
    }
    os << d_species_R[d_num_species - 1];
    os << std::endl;
}


/*
 * Put the characteristics of the equation of state mixing rules class into the restart
 * database.
 */
void
EquationOfStateMixingRulesIdealGas::putToRestart(
    const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const
{
    restart_db->putDoubleVector("d_species_gamma", d_species_gamma);
    restart_db->putDoubleVector("d_species_R", d_species_R);
}


/*
 * Compute the pressure of the mixture with isothermal and isobaric equilibrium assumptions.
 */
double
EquationOfStateMixingRulesIdealGas::getPressure(
    const double* const density,
    const double* const internal_energy,
    const std::vector<const double*>& mass_fractions) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT((d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC) ||
                (d_mixing_closure_model == MIXING_CLOSURE_MODEL::NO_MODEL && d_num_species == 1));
    TBOX_ASSERT((static_cast<int>(mass_fractions.size()) == d_num_species) ||
                (static_cast<int>(mass_fractions.size()) == d_num_species - 1));
#endif
    
    // Get the mixture thermodynamic properties.
    std::vector<double> mixture_thermo_properties;
    std::vector<double*> mixture_thermo_properties_ptr;
    std::vector<const double*> mixture_thermo_properties_const_ptr;
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    mixture_thermo_properties.resize(num_thermo_properties);
    mixture_thermo_properties_ptr.reserve(num_thermo_properties);
    mixture_thermo_properties_const_ptr.reserve(num_thermo_properties);
    
    for (int ti = 0; ti < num_thermo_properties; ti++)
    {
        mixture_thermo_properties_ptr.push_back(&mixture_thermo_properties[ti]);
        mixture_thermo_properties_const_ptr.push_back(&mixture_thermo_properties[ti]);
    }
    
    getMixtureThermodynamicProperties(
        mixture_thermo_properties_ptr,
        mass_fractions);
    
    return d_equation_of_state->getPressure(
        density,
        internal_energy,
        mixture_thermo_properties_const_ptr);
}


/*
 * Compute the pressure of the mixture with isothermal and isobaric equilibrium assumptions.
 */
void
EquationOfStateMixingRulesIdealGas::computePressure(
    HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_internal_energy,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_mass_fractions,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT((d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC) ||
                (d_mixing_closure_model == MIXING_CLOSURE_MODEL::NO_MODEL && d_num_species == 1));
    
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_internal_energy);
    TBOX_ASSERT(data_mass_fractions);
    
    TBOX_ASSERT((data_mass_fractions->getDepth() == d_num_species) ||
                (data_mass_fractions->getDepth() == d_num_species - 1));
#endif
    
    /*
     * Get the mixture thermodyanmic properties.
     */
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    // Declare data container for mixture thermodyanmic properties.
    HAMERS_SHARED_PTR<pdat::CellData<double> > data_mixture_thermo_properties;
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_internal_energy = data_internal_energy->getGhostCellWidth();
        const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_pressure->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_internal_energy->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_mass_fractions->getBox().isSpatiallyEqual(interior_box));
#endif
        
        /*
         * Get the minimum number of ghost cells for mixture thermodynamic properties.
         */
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_pressure;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_internal_energy, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_mass_fractions, num_ghosts_min);
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::CellData<double> >(
            interior_box, num_thermo_properties, num_ghosts_min);
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_internal_energy->getGhostBox().contains(domain));
        TBOX_ASSERT(data_mass_fractions->getGhostBox().contains(domain));
#endif
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::CellData<double> >(
            domain, num_thermo_properties, hier::IntVector::getZero(d_dim));
    }
    
    computeMixtureThermodynamicProperties(
        data_mixture_thermo_properties,
        data_mass_fractions,
        domain);
    
    d_equation_of_state->computePressure(
        data_pressure,
        data_density,
        data_internal_energy,
        data_mixture_thermo_properties,
        domain);
}


/*
 * Compute the pressure of the mixture with isothermal and isobaric equilibrium assumptions.
 */
void
EquationOfStateMixingRulesIdealGas::computePressure(
    HAMERS_SHARED_PTR<pdat::SideData<double> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
    const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_internal_energy,
    const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_mass_fractions,
    int side_normal,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT((d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC) ||
                (d_mixing_closure_model == MIXING_CLOSURE_MODEL::NO_MODEL && d_num_species == 1));
    
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_internal_energy);
    TBOX_ASSERT(data_mass_fractions);
    
    TBOX_ASSERT((data_mass_fractions->getDepth() == d_num_species) ||
                (data_mass_fractions->getDepth() == d_num_species - 1));
#endif
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(side_normal < d_dim.getValue());
    
    TBOX_ASSERT(data_pressure->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_density->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_internal_energy->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_mass_fractions->getDirectionVector()[side_normal] > 0);
#endif
    
    hier::IntVector direction = hier::IntVector::getZero(d_dim);
    direction[side_normal] = 1;
    
    /*
     * Get the mixture thermodyanmic properties.
     */
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    // Declare data container for mixture thermodyanmic properties.
    HAMERS_SHARED_PTR<pdat::SideData<double> > data_mixture_thermo_properties;
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_internal_energy = data_internal_energy->getGhostCellWidth();
        const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_pressure->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_internal_energy->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_mass_fractions->getBox().isSpatiallyEqual(interior_box));
#endif
        
        /*
         * Get the minimum number of ghost cells for mixture thermodynamic properties.
         */
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_pressure;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_internal_energy, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_mass_fractions, num_ghosts_min);
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::SideData<double> >(
            interior_box, num_thermo_properties, num_ghosts_min, direction);
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_internal_energy->getGhostBox().contains(domain));
        TBOX_ASSERT(data_mass_fractions->getGhostBox().contains(domain));
#endif
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::SideData<double> >(
            domain, num_thermo_properties, hier::IntVector::getZero(d_dim), direction);
    }
    
    computeMixtureThermodynamicProperties(
        data_mixture_thermo_properties,
        data_mass_fractions,
        side_normal,
        domain);
    
    d_equation_of_state->computePressure(
        data_pressure,
        data_density,
        data_internal_energy,
        data_mixture_thermo_properties,
        side_normal,
        domain);
}


/*
 * Compute the pressure of the mixture with isobaric equilibrium assumption.
 */
double
EquationOfStateMixingRulesIdealGas::getPressure(
    const double* const density,
    const double* const internal_energy,
    const std::vector<const double*>& mass_fractions,
    const std::vector<const double*>& volume_fractions) const
{
    NULL_USE(mass_fractions);
    
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOBARIC);
    TBOX_ASSERT((static_cast<int>(volume_fractions.size()) == d_num_species) ||
                (static_cast<int>(volume_fractions.size()) == d_num_species - 1));
#endif
    
    // Get the mixture thermodynamic properties.
    std::vector<double> mixture_thermo_properties;
    std::vector<double*> mixture_thermo_properties_ptr;
    std::vector<const double*> mixture_thermo_properties_const_ptr;
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    mixture_thermo_properties.resize(num_thermo_properties);
    mixture_thermo_properties_ptr.reserve(num_thermo_properties);
    mixture_thermo_properties_const_ptr.reserve(num_thermo_properties);
    
    for (int ti = 0; ti < num_thermo_properties; ti++)
    {
        mixture_thermo_properties_ptr.push_back(&mixture_thermo_properties[ti]);
        mixture_thermo_properties_const_ptr.push_back(&mixture_thermo_properties[ti]);
    }
    
    getMixtureThermodynamicProperties(
        mixture_thermo_properties_ptr,
        volume_fractions);
    
    return d_equation_of_state->getPressure(
        density,
        internal_energy,
        mixture_thermo_properties_const_ptr);
}


/*
 * Compute the pressure of the mixture with isobaric equilibrium assumption.
 */
void
EquationOfStateMixingRulesIdealGas::computePressure(
    HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_internal_energy,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_mass_fractions,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_volume_fractions,
    const hier::Box& domain) const
{
    NULL_USE(data_mass_fractions);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOBARIC);
    
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_internal_energy);
    TBOX_ASSERT(data_volume_fractions);
    
    TBOX_ASSERT((data_volume_fractions->getDepth() == d_num_species) ||
                (data_volume_fractions->getDepth() == d_num_species - 1));
#endif
    
    /*
     * Get the mixture thermodyanmic properties.
     */
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    // Declare data container for mixture thermodyanmic properties.
    HAMERS_SHARED_PTR<pdat::CellData<double> > data_mixture_thermo_properties;
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_internal_energy = data_internal_energy->getGhostCellWidth();
        const hier::IntVector num_ghosts_volume_fractions = data_volume_fractions->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_pressure->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_internal_energy->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_volume_fractions->getBox().isSpatiallyEqual(interior_box));
#endif
        
        /*
         * Get the minimum number of ghost cells for mixture thermodynamic properties.
         */
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_pressure;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_internal_energy, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_volume_fractions, num_ghosts_min);
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::CellData<double> >(
            interior_box, num_thermo_properties, num_ghosts_min);
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_internal_energy->getGhostBox().contains(domain));
        TBOX_ASSERT(data_volume_fractions->getGhostBox().contains(domain));
#endif
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::CellData<double> >(
            domain, num_thermo_properties, hier::IntVector::getZero(d_dim));
    }
    
    computeMixtureThermodynamicProperties(
        data_mixture_thermo_properties,
        data_volume_fractions,
        domain);
    
    d_equation_of_state->computePressure(
        data_pressure,
        data_density,
        data_internal_energy,
        data_mixture_thermo_properties,
        domain);
}


/*
 * Compute the pressure of the mixture with isobaric equilibrium assumption.
 */
void
EquationOfStateMixingRulesIdealGas::computePressure(
    HAMERS_SHARED_PTR<pdat::SideData<double> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
    const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_internal_energy,
    const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_mass_fractions,
    const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_volume_fractions,
    int side_normal,
    const hier::Box& domain) const
{
    NULL_USE(data_mass_fractions);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOBARIC);
    
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_internal_energy);
    TBOX_ASSERT(data_volume_fractions);
    
    TBOX_ASSERT((data_volume_fractions->getDepth() == d_num_species) ||
                (data_volume_fractions->getDepth() == d_num_species - 1));
#endif
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(side_normal < d_dim.getValue());
    
    TBOX_ASSERT(data_pressure->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_density->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_internal_energy->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_volume_fractions->getDirectionVector()[side_normal] > 0);
#endif
    
    hier::IntVector direction = hier::IntVector::getZero(d_dim);
    direction[side_normal] = 1;
    
    /*
     * Get the mixture thermodyanmic properties.
     */
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    // Declare data container for mixture thermodyanmic properties.
    HAMERS_SHARED_PTR<pdat::SideData<double> > data_mixture_thermo_properties;
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_internal_energy = data_internal_energy->getGhostCellWidth();
        const hier::IntVector num_ghosts_volume_fractions = data_volume_fractions->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_pressure->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_internal_energy->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_volume_fractions->getBox().isSpatiallyEqual(interior_box));
#endif
        
        /*
         * Get the minimum number of ghost cells for mixture thermodynamic properties.
         */
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_pressure;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_internal_energy, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_volume_fractions, num_ghosts_min);
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::SideData<double> >(
            interior_box, num_thermo_properties, num_ghosts_min, direction);
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_internal_energy->getGhostBox().contains(domain));
        TBOX_ASSERT(data_volume_fractions->getGhostBox().contains(domain));
#endif
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::SideData<double> >(
            domain, num_thermo_properties, hier::IntVector::getZero(d_dim), direction);
    }
    
    computeMixtureThermodynamicProperties(
        data_mixture_thermo_properties,
        data_volume_fractions,
        side_normal,
        domain);
    
    d_equation_of_state->computePressure(
        data_pressure,
        data_density,
        data_internal_energy,
        data_mixture_thermo_properties,
        side_normal,
        domain);
}


/*
 * Compute the specific internal energy of the mixture with isothermal and isobaric equilibrium assumptions.
 */
double
EquationOfStateMixingRulesIdealGas::getInternalEnergy(
    const double* const density,
    const double* const pressure,
    const std::vector<const double*>& mass_fractions) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT((d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC) ||
                (d_mixing_closure_model == MIXING_CLOSURE_MODEL::NO_MODEL && d_num_species == 1));
    TBOX_ASSERT((static_cast<int>(mass_fractions.size()) == d_num_species) ||
                (static_cast<int>(mass_fractions.size()) == d_num_species - 1));
#endif
    
    // Get the mixture thermodynamic properties.
    std::vector<double> mixture_thermo_properties;
    std::vector<double*> mixture_thermo_properties_ptr;
    std::vector<const double*> mixture_thermo_properties_const_ptr;
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    mixture_thermo_properties.resize(num_thermo_properties);
    mixture_thermo_properties_ptr.reserve(num_thermo_properties);
    mixture_thermo_properties_const_ptr.reserve(num_thermo_properties);
    
    for (int ti = 0; ti < num_thermo_properties; ti++)
    {
        mixture_thermo_properties_ptr.push_back(&mixture_thermo_properties[ti]);
        mixture_thermo_properties_const_ptr.push_back(&mixture_thermo_properties[ti]);
    }
    
    getMixtureThermodynamicProperties(
        mixture_thermo_properties_ptr,
        mass_fractions);
    
    return d_equation_of_state->getInternalEnergy(
        density,
        pressure,
        mixture_thermo_properties_const_ptr);
}


/*
 * Compute the specific internal energy of the mixture with isothermal and isobaric equilibrium assumptions.
 */
void
EquationOfStateMixingRulesIdealGas::computeInternalEnergy(
    HAMERS_SHARED_PTR<pdat::CellData<double> >& data_internal_energy,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_mass_fractions,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT((d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC) ||
                (d_mixing_closure_model == MIXING_CLOSURE_MODEL::NO_MODEL && d_num_species == 1));
    
    TBOX_ASSERT(data_internal_energy);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_mass_fractions);
    
    TBOX_ASSERT((data_mass_fractions->getDepth() == d_num_species) ||
                (data_mass_fractions->getDepth() == d_num_species - 1));
#endif
    
    /*
     * Get the mixture thermodyanmic properties.
     */
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    // Declare data container for mixture thermodyanmic properties.
    HAMERS_SHARED_PTR<pdat::CellData<double> > data_mixture_thermo_properties;
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_internal_energy = data_internal_energy->getGhostCellWidth();
        const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_internal_energy->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_pressure->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_mass_fractions->getBox().isSpatiallyEqual(interior_box));
#endif
        
        /*
         * Get the minimum number of ghost cells for mixture thermodynamic properties.
         */
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_internal_energy;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_mass_fractions, num_ghosts_min);
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::CellData<double> >(
            interior_box, num_thermo_properties, num_ghosts_min);
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_internal_energy->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_mass_fractions->getGhostBox().contains(domain));
#endif
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::CellData<double> >(
            domain, num_thermo_properties, hier::IntVector::getZero(d_dim));
    }
    
    computeMixtureThermodynamicProperties(
        data_mixture_thermo_properties,
        data_mass_fractions,
        domain);
    
    d_equation_of_state->computeInternalEnergy(
        data_internal_energy,
        data_density,
        data_pressure,
        data_mixture_thermo_properties,
        domain);
}


/*
 * Compute the specific internal energy of the mixture with isothermal and isobaric equilibrium assumptions.
 */
void
EquationOfStateMixingRulesIdealGas::computeInternalEnergy(
    HAMERS_SHARED_PTR<pdat::SideData<double> >& data_internal_energy,
    const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
    const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_mass_fractions,
    int side_normal,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT((d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC) ||
                (d_mixing_closure_model == MIXING_CLOSURE_MODEL::NO_MODEL && d_num_species == 1));
    
    TBOX_ASSERT(data_internal_energy);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_mass_fractions);
    
    TBOX_ASSERT((data_mass_fractions->getDepth() == d_num_species) ||
                (data_mass_fractions->getDepth() == d_num_species - 1));
#endif
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(side_normal < d_dim.getValue());
    
    TBOX_ASSERT(data_internal_energy->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_density->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_pressure->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_mass_fractions->getDirectionVector()[side_normal] > 0);
#endif
    
    hier::IntVector direction = hier::IntVector::getZero(d_dim);
    direction[side_normal] = 1;
    
    /*
     * Get the mixture thermodyanmic properties.
     */
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    // Declare data container for mixture thermodyanmic properties.
    HAMERS_SHARED_PTR<pdat::SideData<double> > data_mixture_thermo_properties;
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_internal_energy = data_internal_energy->getGhostCellWidth();
        const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_internal_energy->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_pressure->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_mass_fractions->getBox().isSpatiallyEqual(interior_box));
#endif
        
        /*
         * Get the minimum number of ghost cells for mixture thermodynamic properties.
         */
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_internal_energy;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_mass_fractions, num_ghosts_min);
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::SideData<double> >(
            interior_box, num_thermo_properties, num_ghosts_min, direction);
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_internal_energy->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_mass_fractions->getGhostBox().contains(domain));
#endif
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::SideData<double> >(
            domain, num_thermo_properties, hier::IntVector::getZero(d_dim), direction);
    }
    
    computeMixtureThermodynamicProperties(
        data_mixture_thermo_properties,
        data_mass_fractions,
        side_normal,
        domain);
    
    d_equation_of_state->computeInternalEnergy(
        data_internal_energy,
        data_density,
        data_pressure,
        data_mixture_thermo_properties,
        side_normal,
        domain);
}


/*
 * Compute the specific internal energy of the mixture with isobaric equilibrium assumption.
 */
double
EquationOfStateMixingRulesIdealGas::getInternalEnergy(
    const double* const density,
    const double* const pressure,
    const std::vector<const double*>& mass_fractions,
    const std::vector<const double*>& volume_fractions) const
{
    NULL_USE(mass_fractions);
    
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOBARIC);
    TBOX_ASSERT((static_cast<int>(volume_fractions.size()) == d_num_species) ||
                (static_cast<int>(volume_fractions.size()) == d_num_species - 1));
#endif
    
    // Get the mixture thermodynamic properties.
    std::vector<double> mixture_thermo_properties;
    std::vector<double*> mixture_thermo_properties_ptr;
    std::vector<const double*> mixture_thermo_properties_const_ptr;
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    mixture_thermo_properties.resize(num_thermo_properties);
    mixture_thermo_properties_ptr.reserve(num_thermo_properties);
    mixture_thermo_properties_const_ptr.reserve(num_thermo_properties);
    
    for (int ti = 0; ti < num_thermo_properties; ti++)
    {
        mixture_thermo_properties_ptr.push_back(&mixture_thermo_properties[ti]);
        mixture_thermo_properties_const_ptr.push_back(&mixture_thermo_properties[ti]);
    }
    
    getMixtureThermodynamicProperties(
        mixture_thermo_properties_ptr,
        volume_fractions);
    
    return d_equation_of_state->getInternalEnergy(
        density,
        pressure,
        mixture_thermo_properties_const_ptr);
}


/*
 * Compute the specific internal energy of the mixture with isobaric equilibrium assumption.
 */
void
EquationOfStateMixingRulesIdealGas::computeInternalEnergy(
    HAMERS_SHARED_PTR<pdat::CellData<double> >& data_internal_energy,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_mass_fractions,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_volume_fractions,
    const hier::Box& domain) const
{
    NULL_USE(data_mass_fractions);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOBARIC);
    
    TBOX_ASSERT(data_internal_energy);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_volume_fractions);
    
    TBOX_ASSERT((data_volume_fractions->getDepth() == d_num_species) ||
                (data_volume_fractions->getDepth() == d_num_species - 1));
#endif
    
    /*
     * Get the mixture thermodyanmic properties.
     */
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    // Declare data container for mixture thermodyanmic properties.
    HAMERS_SHARED_PTR<pdat::CellData<double> > data_mixture_thermo_properties;
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_internal_energy = data_internal_energy->getGhostCellWidth();
        const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        const hier::IntVector num_ghosts_volume_fractions = data_volume_fractions->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_internal_energy->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_pressure->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_volume_fractions->getBox().isSpatiallyEqual(interior_box));
#endif
        
        /*
         * Get the minimum number of ghost cells for mixture thermodynamic properties.
         */
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_internal_energy;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_volume_fractions, num_ghosts_min);
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::CellData<double> >(
            interior_box, num_thermo_properties, num_ghosts_min);
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_internal_energy->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_volume_fractions->getGhostBox().contains(domain));
#endif
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::CellData<double> >(
            domain, num_thermo_properties, hier::IntVector::getZero(d_dim));
    }
    
    computeMixtureThermodynamicProperties(
        data_mixture_thermo_properties,
        data_volume_fractions,
        domain);
    
    d_equation_of_state->computeInternalEnergy(
        data_internal_energy,
        data_density,
        data_pressure,
        data_mixture_thermo_properties,
        domain);
}


/*
 * Compute the specific internal energy of the mixture with isobaric equilibrium assumption.
 */
void
EquationOfStateMixingRulesIdealGas::computeInternalEnergy(
    HAMERS_SHARED_PTR<pdat::SideData<double> >& data_internal_energy,
    const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
    const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_mass_fractions,
    const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_volume_fractions,
    int side_normal,
    const hier::Box& domain) const
{
    NULL_USE(data_mass_fractions);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOBARIC);
    
    TBOX_ASSERT(data_internal_energy);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_volume_fractions);
    
    TBOX_ASSERT((data_volume_fractions->getDepth() == d_num_species) ||
                (data_volume_fractions->getDepth() == d_num_species - 1));
#endif
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(side_normal < d_dim.getValue());
    
    TBOX_ASSERT(data_internal_energy->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_density->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_pressure->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_volume_fractions->getDirectionVector()[side_normal] > 0);
#endif
    
    hier::IntVector direction = hier::IntVector::getZero(d_dim);
    direction[side_normal] = 1;
    
    /*
     * Get the mixture thermodyanmic properties.
     */
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    // Declare data container for mixture thermodyanmic properties.
    HAMERS_SHARED_PTR<pdat::SideData<double> > data_mixture_thermo_properties;
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_internal_energy = data_internal_energy->getGhostCellWidth();
        const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        const hier::IntVector num_ghosts_volume_fractions = data_volume_fractions->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_internal_energy->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_pressure->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_volume_fractions->getBox().isSpatiallyEqual(interior_box));
#endif
        
        /*
         * Get the minimum number of ghost cells for mixture thermodynamic properties.
         */
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_internal_energy;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_volume_fractions, num_ghosts_min);
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::SideData<double> >(
            interior_box, num_thermo_properties, num_ghosts_min, direction);
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_internal_energy->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_volume_fractions->getGhostBox().contains(domain));
#endif
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::SideData<double> >(
            domain, num_thermo_properties, hier::IntVector::getZero(d_dim), direction);
    }
    
    computeMixtureThermodynamicProperties(
        data_mixture_thermo_properties,
        data_volume_fractions,
        side_normal,
        domain);
    
    d_equation_of_state->computeInternalEnergy(
        data_internal_energy,
        data_density,
        data_pressure,
        data_mixture_thermo_properties,
        side_normal,
        domain);
}


/*
 * Compute the temperature of the mixture with isothermal and isobaric equilibrium assumptions.
 */
double
EquationOfStateMixingRulesIdealGas::getTemperature(
    const double* const density,
    const double* const pressure,
    const std::vector<const double*>& mass_fractions) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT((d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC) ||
                (d_mixing_closure_model == MIXING_CLOSURE_MODEL::NO_MODEL && d_num_species == 1));
    TBOX_ASSERT((static_cast<int>(mass_fractions.size()) == d_num_species) ||
                (static_cast<int>(mass_fractions.size()) == d_num_species - 1));
#endif
    
    // Get the mixture thermodynamic properties.
    std::vector<double> mixture_thermo_properties;
    std::vector<double*> mixture_thermo_properties_ptr;
    std::vector<const double*> mixture_thermo_properties_const_ptr;
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    mixture_thermo_properties.resize(num_thermo_properties);
    mixture_thermo_properties_ptr.reserve(num_thermo_properties);
    mixture_thermo_properties_const_ptr.reserve(num_thermo_properties);
    
    for (int ti = 0; ti < num_thermo_properties; ti++)
    {
        mixture_thermo_properties_ptr.push_back(&mixture_thermo_properties[ti]);
        mixture_thermo_properties_const_ptr.push_back(&mixture_thermo_properties[ti]);
    }
    
    getMixtureThermodynamicProperties(
        mixture_thermo_properties_ptr,
        mass_fractions);
    
    return d_equation_of_state->getTemperature(
        density,
        pressure,
        mixture_thermo_properties_const_ptr);
}


/*
 * Compute the temperature of the mixture with isothermal and isobaric equilibrium assumptions.
 */
void
EquationOfStateMixingRulesIdealGas::computeTemperature(
    HAMERS_SHARED_PTR<pdat::CellData<double> >& data_temperature,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_mass_fractions,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT((d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC) ||
                (d_mixing_closure_model == MIXING_CLOSURE_MODEL::NO_MODEL && d_num_species == 1));
    
    TBOX_ASSERT(data_temperature);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_mass_fractions);
    
    TBOX_ASSERT((data_mass_fractions->getDepth() == d_num_species) ||
                (data_mass_fractions->getDepth() == d_num_species - 1));
#endif
    
    /*
     * Get the mixture thermodyanmic properties.
     */
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    // Declare data container for mixture thermodyanmic properties.
    HAMERS_SHARED_PTR<pdat::CellData<double> > data_mixture_thermo_properties;
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
        const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_temperature->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_pressure->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_mass_fractions->getBox().isSpatiallyEqual(interior_box));
#endif
        
        /*
         * Get the minimum number of ghost cells for mixture thermodynamic properties.
         */
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_temperature;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_mass_fractions, num_ghosts_min);
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::CellData<double> >(
            interior_box, num_thermo_properties, num_ghosts_min);
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_temperature->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_mass_fractions->getGhostBox().contains(domain));
#endif
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::CellData<double> >(
            domain, num_thermo_properties, hier::IntVector::getZero(d_dim));
    }
    
    computeMixtureThermodynamicProperties(
        data_mixture_thermo_properties,
        data_mass_fractions,
        domain);
    
    d_equation_of_state->computeTemperature(
        data_temperature,
        data_density,
        data_pressure,
        data_mixture_thermo_properties,
        domain);
}


/*
 * Compute the temperature of the mixture with isothermal and isobaric equilibrium assumptions.
 */
void
EquationOfStateMixingRulesIdealGas::computeTemperature(
    HAMERS_SHARED_PTR<pdat::SideData<double> >& data_temperature,
    const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
    const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_mass_fractions,
    int side_normal,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT((d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC) ||
                (d_mixing_closure_model == MIXING_CLOSURE_MODEL::NO_MODEL && d_num_species == 1));
    
    TBOX_ASSERT(data_temperature);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_mass_fractions);
    
    TBOX_ASSERT((data_mass_fractions->getDepth() == d_num_species) ||
                (data_mass_fractions->getDepth() == d_num_species - 1));
#endif
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(side_normal < d_dim.getValue());
    
    TBOX_ASSERT(data_temperature->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_density->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_pressure->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_mass_fractions->getDirectionVector()[side_normal] > 0);
#endif
    
    hier::IntVector direction = hier::IntVector::getZero(d_dim);
    direction[side_normal] = 1;
    
    /*
     * Get the mixture thermodyanmic properties.
     */
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    // Declare data container for mixture thermodyanmic properties.
    HAMERS_SHARED_PTR<pdat::SideData<double> > data_mixture_thermo_properties;
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
        const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_temperature->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_pressure->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_mass_fractions->getBox().isSpatiallyEqual(interior_box));
#endif
        
        /*
         * Get the minimum number of ghost cells for mixture thermodynamic properties.
         */
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_temperature;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_mass_fractions, num_ghosts_min);
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::SideData<double> >(
            interior_box, num_thermo_properties, num_ghosts_min, direction);
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_temperature->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_mass_fractions->getGhostBox().contains(domain));
#endif
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::SideData<double> >(
            domain, num_thermo_properties, hier::IntVector::getZero(d_dim), direction);
    }
    
    computeMixtureThermodynamicProperties(
        data_mixture_thermo_properties,
        data_mass_fractions,
        side_normal,
        domain);
    
    d_equation_of_state->computeTemperature(
        data_temperature,
        data_density,
        data_pressure,
        data_mixture_thermo_properties,
        side_normal,
        domain);
}


/*
 * Compute the specific internal energy of the mixture from temperature with isothermal
 * and isobaric equilibrium assumptions.
 */
double
EquationOfStateMixingRulesIdealGas::getInternalEnergyFromTemperature(
    const double* const density,
    const double* const temperature,
    const std::vector<const double*>& mass_fractions) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT((d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC) ||
                (d_mixing_closure_model == MIXING_CLOSURE_MODEL::NO_MODEL && d_num_species == 1));
    TBOX_ASSERT((static_cast<int>(mass_fractions.size()) == d_num_species) ||
                (static_cast<int>(mass_fractions.size()) == d_num_species - 1));
#endif
    
    // Get the mixture thermodynamic properties.
    std::vector<double> mixture_thermo_properties;
    std::vector<double*> mixture_thermo_properties_ptr;
    std::vector<const double*> mixture_thermo_properties_const_ptr;
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    mixture_thermo_properties.resize(num_thermo_properties);
    mixture_thermo_properties_ptr.reserve(num_thermo_properties);
    mixture_thermo_properties_const_ptr.reserve(num_thermo_properties);
    
    for (int ti = 0; ti < num_thermo_properties; ti++)
    {
        mixture_thermo_properties_ptr.push_back(&mixture_thermo_properties[ti]);
        mixture_thermo_properties_const_ptr.push_back(&mixture_thermo_properties[ti]);
    }
    
    getMixtureThermodynamicProperties(
        mixture_thermo_properties_ptr,
        mass_fractions);
    
    return d_equation_of_state->getInternalEnergyFromTemperature(
        density,
        temperature,
        mixture_thermo_properties_const_ptr);
}


/*
 * Compute the specific internal energy of the mixture from temperature with isothermal
 * and isobaric equilibrium assumptions.
 */
void
EquationOfStateMixingRulesIdealGas::computeInternalEnergyFromTemperature(
    HAMERS_SHARED_PTR<pdat::CellData<double> >& data_internal_energy,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_temperature,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_mass_fractions,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT((d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC) ||
                (d_mixing_closure_model == MIXING_CLOSURE_MODEL::NO_MODEL && d_num_species == 1));
    
    TBOX_ASSERT(data_internal_energy);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_temperature);
    TBOX_ASSERT(data_mass_fractions);
    
    TBOX_ASSERT((data_mass_fractions->getDepth() == d_num_species) ||
                (data_mass_fractions->getDepth() == d_num_species - 1));
#endif
    
    /*
     * Get the mixture thermodyanmic properties.
     */
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    // Declare data container for mixture thermodyanmic properties.
    HAMERS_SHARED_PTR<pdat::CellData<double> > data_mixture_thermo_properties;
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_internal_energy = data_internal_energy->getGhostCellWidth();
        const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
        const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_internal_energy->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_temperature->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_mass_fractions->getBox().isSpatiallyEqual(interior_box));
#endif
        
        /*
         * Get the minimum number of ghost cells for mixture thermodynamic properties.
         */
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_internal_energy;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_temperature, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_mass_fractions, num_ghosts_min);
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::CellData<double> >(
            interior_box, num_thermo_properties, num_ghosts_min);
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_internal_energy->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_temperature->getGhostBox().contains(domain));
        TBOX_ASSERT(data_mass_fractions->getGhostBox().contains(domain));
#endif
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::CellData<double> >(
            domain, num_thermo_properties, hier::IntVector::getZero(d_dim));
    }
    
    computeMixtureThermodynamicProperties(
        data_mixture_thermo_properties,
        data_mass_fractions,
        domain);
    
    d_equation_of_state->computeInternalEnergyFromTemperature(
        data_internal_energy,
        data_density,
        data_temperature,
        data_mixture_thermo_properties,
        domain);
}


/*
 * Compute the specific internal energy of the mixture from temperature with isothermal
 * and isobaric equilibrium assumptions.
 */
void
EquationOfStateMixingRulesIdealGas::computeInternalEnergyFromTemperature(
    HAMERS_SHARED_PTR<pdat::SideData<double> >& data_internal_energy,
    const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
    const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_temperature,
    const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_mass_fractions,
    int side_normal,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT((d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC) ||
                (d_mixing_closure_model == MIXING_CLOSURE_MODEL::NO_MODEL && d_num_species == 1));
    
    TBOX_ASSERT(data_internal_energy);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_temperature);
    TBOX_ASSERT(data_mass_fractions);
    
    TBOX_ASSERT((data_mass_fractions->getDepth() == d_num_species) ||
                (data_mass_fractions->getDepth() == d_num_species - 1));
#endif
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(side_normal < d_dim.getValue());
    
    TBOX_ASSERT(data_internal_energy->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_density->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_temperature->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_mass_fractions->getDirectionVector()[side_normal] > 0);
#endif
    
    hier::IntVector direction = hier::IntVector::getZero(d_dim);
    direction[side_normal] = 1;
    
    /*
     * Get the mixture thermodyanmic properties.
     */
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    // Declare data container for mixture thermodyanmic properties.
    HAMERS_SHARED_PTR<pdat::SideData<double> > data_mixture_thermo_properties;
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_internal_energy = data_internal_energy->getGhostCellWidth();
        const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
        const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_internal_energy->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_temperature->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_mass_fractions->getBox().isSpatiallyEqual(interior_box));
#endif
        
        /*
         * Get the minimum number of ghost cells for mixture thermodynamic properties.
         */
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_internal_energy;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_temperature, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_mass_fractions, num_ghosts_min);
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::SideData<double> >(
            interior_box, num_thermo_properties, num_ghosts_min, direction);
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_internal_energy->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_temperature->getGhostBox().contains(domain));
        TBOX_ASSERT(data_mass_fractions->getGhostBox().contains(domain));
#endif
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::SideData<double> >(
            domain, num_thermo_properties, hier::IntVector::getZero(d_dim), direction);
    }
    
    computeMixtureThermodynamicProperties(
        data_mixture_thermo_properties,
        data_mass_fractions,
        side_normal,
        domain);
    
    d_equation_of_state->computeInternalEnergyFromTemperature(
        data_internal_energy,
        data_density,
        data_temperature,
        data_mixture_thermo_properties,
        side_normal,
        domain);
}


/*
 * Compute the isochoric specific heat capacity of mixture with isothermal and isobaric equilibrium
 * assumptions.
 */
double
EquationOfStateMixingRulesIdealGas::getIsochoricSpecificHeatCapacity(
    const double* const density,
    const double* const pressure,
    const std::vector<const double*>& mass_fractions) const
{
    NULL_USE(density);
    NULL_USE(pressure);
    
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT((d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC) ||
                (d_mixing_closure_model == MIXING_CLOSURE_MODEL::NO_MODEL && d_num_species == 1));
    TBOX_ASSERT((static_cast<int>(mass_fractions.size()) == d_num_species) ||
                (static_cast<int>(mass_fractions.size()) == d_num_species - 1));
#endif
    
    double c_v = double(0);
    
    if (static_cast<int>(mass_fractions.size()) == d_num_species)
    {
        for (int si = 0; si < d_num_species; si++)
        {
            c_v += *(mass_fractions[si])*d_species_c_v[si];
        }
    }
    else if (static_cast<int>(mass_fractions.size()) == d_num_species - 1)
    {
        double Y_last = double(1);
        
        for (int si = 0; si < d_num_species - 1; si++)
        {
            c_v += *(mass_fractions[si])*d_species_c_v[si];
            
            // Compute the mass fraction of the last species.
            Y_last -= *(mass_fractions[si]);
        }
        
        // Add the contribution from the last species.
        c_v += Y_last*d_species_c_v.back();
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Number of mass fractions provided is not"
            << " equal to the total number of species or (total number of species - 1)."
            << std::endl);
    }
    
    return c_v;
}


/*
 * Compute the isochoric specific heat capacity of mixture with isothermal and isobaric equilibrium
 * assumptions.
 */
void
EquationOfStateMixingRulesIdealGas::computeIsochoricSpecificHeatCapacity(
    HAMERS_SHARED_PTR<pdat::CellData<double> >& data_isochoric_specific_heat_capacity,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_mass_fractions,
    const hier::Box& domain) const
{
    NULL_USE(data_density);
    NULL_USE(data_pressure);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT((d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC) ||
                (d_mixing_closure_model == MIXING_CLOSURE_MODEL::NO_MODEL && d_num_species == 1));
    
    TBOX_ASSERT(data_isochoric_specific_heat_capacity);
    TBOX_ASSERT(data_mass_fractions);
    
    TBOX_ASSERT((data_mass_fractions->getDepth() == d_num_species) ||
                (data_mass_fractions->getDepth() == d_num_species - 1));
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_isochoric_specific_heat_capacity =
        data_isochoric_specific_heat_capacity->getGhostBox();
    const hier::IntVector ghostcell_dims_isochoric_specific_heat_capacity =
        ghost_box_isochoric_specific_heat_capacity.numberCells();
    
    const hier::Box ghost_box_mass_fractions = data_mass_fractions->getGhostBox();
    const hier::IntVector ghostcell_dims_mass_fractions = ghost_box_mass_fractions.numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_isochoric_specific_heat_capacity(d_dim);
    hier::IntVector offset_mass_fractions(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_isochoric_specific_heat_capacity =
            data_isochoric_specific_heat_capacity->getGhostCellWidth();
        
        const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_isochoric_specific_heat_capacity->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_mass_fractions->getBox().isSpatiallyEqual(interior_box));
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_isochoric_specific_heat_capacity;
        num_ghosts_min = hier::IntVector::min(num_ghosts_mass_fractions, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_isochoric_specific_heat_capacity = num_ghosts_isochoric_specific_heat_capacity;
        offset_mass_fractions = num_ghosts_mass_fractions;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
        TBOX_ASSERT(data_isochoric_specific_heat_capacity->getGhostBox().contains(domain));
        TBOX_ASSERT(data_mass_fractions->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_isochoric_specific_heat_capacity = domain.lower() - ghost_box_isochoric_specific_heat_capacity.lower();
        offset_mass_fractions = domain.lower() - ghost_box_mass_fractions.lower();
    }
    
    /*
     * Get the pointer to the cell data of isochoric specific heat capacity.
     */
    
    double* const c_v = data_isochoric_specific_heat_capacity->getPointer(0);
    
    /*
     * Fill zeros for c_v.
     */
    
    if (domain.empty())
    {
        data_isochoric_specific_heat_capacity->fillAll(double(0));
    }
    else
    {
        data_isochoric_specific_heat_capacity->fillAll(double(0), domain);
    }
    
    if (data_mass_fractions->getDepth() == d_num_species)
    {
        /*
         * Get the pointers to the cell data of mass fractions.
         */
        
        std::vector<const double*> Y;
        Y.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            Y.push_back(data_mass_fractions->getPointer(si));
        }
        
        computeIsochoricSpecificHeatCapacity(
            c_v,
            Y,
            offset_isochoric_specific_heat_capacity,
            offset_mass_fractions,
            ghostcell_dims_isochoric_specific_heat_capacity,
            ghostcell_dims_mass_fractions,
            domain_lo,
            domain_dims);
    }
    else if (data_mass_fractions->getDepth() == d_num_species - 1)
    {
        HAMERS_SHARED_PTR<pdat::CellData<double> > data_mass_fractions_last;
        
        hier::IntVector offset_mass_fractions_last(d_dim);
        
        if (domain.empty())
        {
            const hier::Box interior_box = data_mass_fractions->getBox();
            
            offset_mass_fractions_last = data_mass_fractions->getGhostCellWidth();
            
            data_mass_fractions_last =
                HAMERS_MAKE_SHARED<pdat::CellData<double> >(interior_box, 1, data_mass_fractions->getGhostCellWidth());
            
            data_mass_fractions_last->fillAll(double(1));
        }
        else
        {
            offset_mass_fractions_last = hier::IntVector::getZero(d_dim);
            
            data_mass_fractions_last =
                HAMERS_MAKE_SHARED<pdat::CellData<double> >(domain, 1, hier::IntVector::getZero(d_dim));
            
            data_mass_fractions_last->fillAll(double(1), domain);
        }
        
        /*
         * Get the pointers to the cell data of mass fractions and the dimensions of the ghost cell box of
         * last mass fraction.
         */
        
        std::vector<const double*> Y;
        Y.reserve(d_num_species - 1);
        for (int si = 0; si < d_num_species - 1; si++)
        {
            Y.push_back(data_mass_fractions->getPointer(si));
        }
        
        double* const Y_last = data_mass_fractions_last->getPointer(0);
        
        const hier::Box ghost_box_mass_fractions_last = data_mass_fractions_last->getGhostBox();
        const hier::IntVector ghostcell_dims_mass_fractions_last = ghost_box_mass_fractions_last.numberCells();
        
        computeIsochoricSpecificHeatCapacity(
            c_v,
            Y_last,
            Y,
            offset_isochoric_specific_heat_capacity,
            offset_mass_fractions_last,
            offset_mass_fractions,
            ghostcell_dims_isochoric_specific_heat_capacity,
            ghostcell_dims_mass_fractions_last,
            ghostcell_dims_mass_fractions,
            domain_lo,
            domain_dims);
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
 * Compute the isochoric specific heat capacity of mixture with isothermal and isobaric equilibrium
 * assumptions.
 */
void
EquationOfStateMixingRulesIdealGas::computeIsochoricSpecificHeatCapacity(
    HAMERS_SHARED_PTR<pdat::SideData<double> >& data_isochoric_specific_heat_capacity,
    const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
    const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_mass_fractions,
    int side_normal,
    const hier::Box& domain) const
{
    NULL_USE(data_density);
    NULL_USE(data_pressure);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT((d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC) ||
                (d_mixing_closure_model == MIXING_CLOSURE_MODEL::NO_MODEL && d_num_species == 1));
    
    TBOX_ASSERT(data_isochoric_specific_heat_capacity);
    TBOX_ASSERT(data_mass_fractions);
    
    TBOX_ASSERT((data_mass_fractions->getDepth() == d_num_species) ||
                (data_mass_fractions->getDepth() == d_num_species - 1));
#endif
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(side_normal < d_dim.getValue());
    
    TBOX_ASSERT(data_isochoric_specific_heat_capacity->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_mass_fractions->getDirectionVector()[side_normal] > 0);
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_isochoric_specific_heat_capacity =
        data_isochoric_specific_heat_capacity->getGhostBox();
    hier::IntVector ghostcell_dims_isochoric_specific_heat_capacity =
        ghost_box_isochoric_specific_heat_capacity.numberCells();
    
    const hier::Box ghost_box_mass_fractions = data_mass_fractions->getGhostBox();
    hier::IntVector ghostcell_dims_mass_fractions = ghost_box_mass_fractions.numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_isochoric_specific_heat_capacity(d_dim);
    hier::IntVector offset_mass_fractions(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_isochoric_specific_heat_capacity =
            data_isochoric_specific_heat_capacity->getGhostCellWidth();
        
        const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_isochoric_specific_heat_capacity->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_mass_fractions->getBox().isSpatiallyEqual(interior_box));
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_isochoric_specific_heat_capacity;
        num_ghosts_min = hier::IntVector::min(num_ghosts_mass_fractions, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_isochoric_specific_heat_capacity = num_ghosts_isochoric_specific_heat_capacity;
        offset_mass_fractions = num_ghosts_mass_fractions;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
        TBOX_ASSERT(data_isochoric_specific_heat_capacity->getGhostBox().contains(domain));
        TBOX_ASSERT(data_mass_fractions->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_isochoric_specific_heat_capacity = domain.lower() - ghost_box_isochoric_specific_heat_capacity.lower();
        offset_mass_fractions = domain.lower() - ghost_box_mass_fractions.lower();
    }
    
    ghostcell_dims_isochoric_specific_heat_capacity[side_normal]++;
    ghostcell_dims_mass_fractions[side_normal]++;
    domain_dims[side_normal]++;
    
    /*
     * Get the pointer to the cell data of isochoric specific heat capacity.
     */
    
    double* const c_v = data_isochoric_specific_heat_capacity->getPointer(side_normal, 0);
    
    /*
     * Fill zeros for c_v.
     */
    
    if (domain.empty())
    {
        data_isochoric_specific_heat_capacity->fillAll(double(0));
    }
    else
    {
        data_isochoric_specific_heat_capacity->fillAll(double(0), domain);
    }
    
    if (data_mass_fractions->getDepth() == d_num_species)
    {
        /*
         * Get the pointers to the cell data of mass fractions.
         */
        
        std::vector<const double*> Y;
        Y.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            Y.push_back(data_mass_fractions->getPointer(side_normal, si));
        }
        
        computeIsochoricSpecificHeatCapacity(
            c_v,
            Y,
            offset_isochoric_specific_heat_capacity,
            offset_mass_fractions,
            ghostcell_dims_isochoric_specific_heat_capacity,
            ghostcell_dims_mass_fractions,
            domain_lo,
            domain_dims);
    }
    else if (data_mass_fractions->getDepth() == d_num_species - 1)
    {
        hier::IntVector direction = hier::IntVector::getZero(d_dim);
        direction[side_normal] = 1;
        
        HAMERS_SHARED_PTR<pdat::SideData<double> > data_mass_fractions_last;
        
        hier::IntVector offset_mass_fractions_last(d_dim);
        
        if (domain.empty())
        {
            const hier::Box interior_box = data_mass_fractions->getBox();
            
            offset_mass_fractions_last = data_mass_fractions->getGhostCellWidth();
            
            HAMERS_SHARED_PTR<pdat::SideData<double> > data_mass_fractions_last =
                HAMERS_MAKE_SHARED<pdat::SideData<double> >(
                    interior_box, 1, data_mass_fractions->getGhostCellWidth(), direction);
            
            data_mass_fractions_last->fillAll(double(1));
        }
        else
        {
            offset_mass_fractions_last = hier::IntVector::getZero(d_dim);
            
            data_mass_fractions_last =
                HAMERS_MAKE_SHARED<pdat::SideData<double> >(domain, 1, hier::IntVector::getZero(d_dim), direction);
            
            data_mass_fractions_last->fillAll(double(1), domain);
        }
        
        /*
         * Get the pointers to the cell data of mass fractions and the dimensions of the ghost cell box of
         * last mass fraction.
         */
        
        std::vector<const double*> Y;
        Y.reserve(d_num_species - 1);
        for (int si = 0; si < d_num_species - 1; si++)
        {
            Y.push_back(data_mass_fractions->getPointer(side_normal, si));
        }
        
        double* const Y_last = data_mass_fractions_last->getPointer(side_normal, 0);
        
        const hier::Box ghost_box_mass_fractions_last = data_mass_fractions_last->getGhostBox();
        hier::IntVector ghostcell_dims_mass_fractions_last = ghost_box_mass_fractions_last.numberCells();
        ghostcell_dims_mass_fractions_last[side_normal]++;
        
        computeIsochoricSpecificHeatCapacity(
            c_v,
            Y_last,
            Y,
            offset_isochoric_specific_heat_capacity,
            offset_mass_fractions_last,
            offset_mass_fractions,
            ghostcell_dims_isochoric_specific_heat_capacity,
            ghostcell_dims_mass_fractions_last,
            ghostcell_dims_mass_fractions,
            domain_lo,
            domain_dims);
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
 * Compute the isobaric specific heat capacity of mixture with isothermal and isobaric equilibrium
 * assumptions.
 */
double
EquationOfStateMixingRulesIdealGas::getIsobaricSpecificHeatCapacity(
    const double* const density,
    const double* const pressure,
    const std::vector<const double*>& mass_fractions) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT((d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC) ||
                (d_mixing_closure_model == MIXING_CLOSURE_MODEL::NO_MODEL && d_num_species == 1));
    TBOX_ASSERT((static_cast<int>(mass_fractions.size()) == d_num_species) ||
                (static_cast<int>(mass_fractions.size()) == d_num_species - 1));
#endif
    
    double c_p = double(0);
    
    if (static_cast<int>(mass_fractions.size()) == d_num_species)
    {
        for (int si = 0; si < d_num_species; si++)
        {
            c_p += *(mass_fractions[si])*d_species_c_p[si];
        }
    }
    else if (static_cast<int>(mass_fractions.size()) == d_num_species - 1)
    {
        double Y_last = double(1);
        
        for (int si = 0; si < d_num_species - 1; si++)
        {
            c_p += *(mass_fractions[si])*d_species_c_p[si];
            
            // Compute the mass fraction of the last species.
            Y_last -= *(mass_fractions[si]);
        }
        
        // Add the contribution from the last species.
        c_p += Y_last*d_species_c_p.back();
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Number of mass fractions provided is not"
            << " equal to the total number of species or (total number of species - 1)."
            << std::endl);
    }
    
    return c_p;
}


/*
 * Compute the isobaric specific heat capacity of mixture with isothermal and isobaric equilibrium
 * assumptions.
 */
void
EquationOfStateMixingRulesIdealGas::computeIsobaricSpecificHeatCapacity(
    HAMERS_SHARED_PTR<pdat::CellData<double> >& data_isobaric_specific_heat_capacity,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_mass_fractions,
    const hier::Box& domain) const
{
    NULL_USE(data_density);
    NULL_USE(data_pressure);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT((d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC) ||
                (d_mixing_closure_model == MIXING_CLOSURE_MODEL::NO_MODEL && d_num_species == 1));
    
    TBOX_ASSERT(data_isobaric_specific_heat_capacity);
    TBOX_ASSERT(data_mass_fractions);
    
    TBOX_ASSERT((data_mass_fractions->getDepth() == d_num_species) ||
                (data_mass_fractions->getDepth() == d_num_species - 1));
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_isobaric_specific_heat_capacity =
        data_isobaric_specific_heat_capacity->getGhostBox();
    const hier::IntVector ghostcell_dims_isobaric_specific_heat_capacity =
        ghost_box_isobaric_specific_heat_capacity.numberCells();
    
    const hier::Box ghost_box_mass_fractions = data_mass_fractions->getGhostBox();
    const hier::IntVector ghostcell_dims_mass_fractions = ghost_box_mass_fractions.numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_isobaric_specific_heat_capacity(d_dim);
    hier::IntVector offset_mass_fractions(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_isobaric_specific_heat_capacity =
            data_isobaric_specific_heat_capacity->getGhostCellWidth();
        
        const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_isobaric_specific_heat_capacity->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_mass_fractions->getBox().isSpatiallyEqual(interior_box));
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_isobaric_specific_heat_capacity;
        num_ghosts_min = hier::IntVector::min(num_ghosts_mass_fractions, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_isobaric_specific_heat_capacity = num_ghosts_isobaric_specific_heat_capacity;
        offset_mass_fractions = num_ghosts_mass_fractions;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
        TBOX_ASSERT(data_isobaric_specific_heat_capacity->getGhostBox().contains(domain));
        TBOX_ASSERT(data_mass_fractions->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_isobaric_specific_heat_capacity = domain.lower() - ghost_box_isobaric_specific_heat_capacity.lower();
        offset_mass_fractions = domain.lower() - ghost_box_mass_fractions.lower();
    }
    
    /*
     * Get the pointer to the cell data of isobaric specific heat capacity.
     */
    
    double* const c_p = data_isobaric_specific_heat_capacity->getPointer(0);
    
    /*
     * Fill zeros for c_p.
     */
    
    if (domain.empty())
    {
        data_isobaric_specific_heat_capacity->fillAll(double(0));
    }
    else
    {
        data_isobaric_specific_heat_capacity->fillAll(double(0), domain);
    }
    
    if (data_mass_fractions->getDepth() == d_num_species)
    {
        /*
         * Get the pointers to the cell data of mass fractions.
         */
        
        std::vector<const double*> Y;
        Y.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            Y.push_back(data_mass_fractions->getPointer(si));
        }
        
        computeIsobaricSpecificHeatCapacity(
            c_p,
            Y,
            offset_isobaric_specific_heat_capacity,
            offset_mass_fractions,
            ghostcell_dims_isobaric_specific_heat_capacity,
            ghostcell_dims_mass_fractions,
            domain_lo,
            domain_dims);
    }
    else if (data_mass_fractions->getDepth() == d_num_species - 1)
    {
        HAMERS_SHARED_PTR<pdat::CellData<double> > data_mass_fractions_last;
        
        hier::IntVector offset_mass_fractions_last(d_dim);
        
        if (domain.empty())
        {
            const hier::Box interior_box = data_mass_fractions->getBox();
            
            offset_mass_fractions_last = data_mass_fractions->getGhostCellWidth();
            
            data_mass_fractions_last =
                HAMERS_MAKE_SHARED<pdat::CellData<double> >(interior_box, 1, data_mass_fractions->getGhostCellWidth());
            
            data_mass_fractions_last->fillAll(double(1));
        }
        else
        {
            offset_mass_fractions_last = hier::IntVector::getZero(d_dim);
            
            data_mass_fractions_last =
                HAMERS_MAKE_SHARED<pdat::CellData<double> >(domain, 1, hier::IntVector::getZero(d_dim));
            
            data_mass_fractions_last->fillAll(double(1), domain);
        }
        
        /*
         * Get the pointers to the cell data of mass fractions and the dimensions of the ghost cell box of
         * last mass fraction.
         */
        
        std::vector<const double*> Y;
        Y.reserve(d_num_species - 1);
        for (int si = 0; si < d_num_species - 1; si++)
        {
            Y.push_back(data_mass_fractions->getPointer(si));
        }
        
        double* const Y_last = data_mass_fractions_last->getPointer(0);
        
        const hier::Box ghost_box_mass_fractions_last = data_mass_fractions_last->getGhostBox();
        const hier::IntVector ghostcell_dims_mass_fractions_last = ghost_box_mass_fractions_last.numberCells();
        
        computeIsobaricSpecificHeatCapacity(
            c_p,
            Y_last,
            Y,
            offset_isobaric_specific_heat_capacity,
            offset_mass_fractions_last,
            offset_mass_fractions,
            ghostcell_dims_isobaric_specific_heat_capacity,
            ghostcell_dims_mass_fractions_last,
            ghostcell_dims_mass_fractions,
            domain_lo,
            domain_dims);
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
 * Compute the isobaric specific heat capacity of mixture with isothermal and isobaric equilibrium
 * assumptions.
 */
void
EquationOfStateMixingRulesIdealGas::computeIsobaricSpecificHeatCapacity(
    HAMERS_SHARED_PTR<pdat::SideData<double> >& data_isobaric_specific_heat_capacity,
    const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
    const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_mass_fractions,
    int side_normal,
    const hier::Box& domain) const
{
    NULL_USE(data_density);
    NULL_USE(data_pressure);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT((d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC) ||
                (d_mixing_closure_model == MIXING_CLOSURE_MODEL::NO_MODEL && d_num_species == 1));
    
    TBOX_ASSERT(data_isobaric_specific_heat_capacity);
    TBOX_ASSERT(data_mass_fractions);
    
    TBOX_ASSERT((data_mass_fractions->getDepth() == d_num_species) ||
                (data_mass_fractions->getDepth() == d_num_species - 1));
#endif
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(side_normal < d_dim.getValue());
    
    TBOX_ASSERT(data_isobaric_specific_heat_capacity->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_mass_fractions->getDirectionVector()[side_normal] > 0);
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_isobaric_specific_heat_capacity =
        data_isobaric_specific_heat_capacity->getGhostBox();
    hier::IntVector ghostcell_dims_isobaric_specific_heat_capacity =
        ghost_box_isobaric_specific_heat_capacity.numberCells();
    
    const hier::Box ghost_box_mass_fractions = data_mass_fractions->getGhostBox();
    hier::IntVector ghostcell_dims_mass_fractions = ghost_box_mass_fractions.numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_isobaric_specific_heat_capacity(d_dim);
    hier::IntVector offset_mass_fractions(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_isobaric_specific_heat_capacity =
            data_isobaric_specific_heat_capacity->getGhostCellWidth();
        
        const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_isobaric_specific_heat_capacity->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_mass_fractions->getBox().isSpatiallyEqual(interior_box));
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_isobaric_specific_heat_capacity;
        num_ghosts_min = hier::IntVector::min(num_ghosts_mass_fractions, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_isobaric_specific_heat_capacity = num_ghosts_isobaric_specific_heat_capacity;
        offset_mass_fractions = num_ghosts_mass_fractions;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
        TBOX_ASSERT(data_isobaric_specific_heat_capacity->getGhostBox().contains(domain));
        TBOX_ASSERT(data_mass_fractions->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_isobaric_specific_heat_capacity = domain.lower() - ghost_box_isobaric_specific_heat_capacity.lower();
        offset_mass_fractions = domain.lower() - ghost_box_mass_fractions.lower();
    }
    
    ghostcell_dims_isobaric_specific_heat_capacity[side_normal]++;
    ghostcell_dims_mass_fractions[side_normal]++;
    domain_dims[side_normal]++;
    
    /*
     * Get the pointer to the cell data of isobaric specific heat capacity.
     */
    
    double* const c_p = data_isobaric_specific_heat_capacity->getPointer(side_normal, 0);
    
    /*
     * Fill zeros for c_p.
     */
    
    if (domain.empty())
    {
        data_isobaric_specific_heat_capacity->fillAll(double(0));
    }
    else
    {
        data_isobaric_specific_heat_capacity->fillAll(double(0), domain);
    }
    
    if (data_mass_fractions->getDepth() == d_num_species)
    {
        /*
         * Get the pointers to the cell data of mass fractions.
         */
        
        std::vector<const double*> Y;
        Y.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            Y.push_back(data_mass_fractions->getPointer(side_normal, si));
        }
        
        computeIsobaricSpecificHeatCapacity(
            c_p,
            Y,
            offset_isobaric_specific_heat_capacity,
            offset_mass_fractions,
            ghostcell_dims_isobaric_specific_heat_capacity,
            ghostcell_dims_mass_fractions,
            domain_lo,
            domain_dims);
    }
    else if (data_mass_fractions->getDepth() == d_num_species - 1)
    {
        hier::IntVector direction = hier::IntVector::getZero(d_dim);
        direction[side_normal] = 1;
        
        HAMERS_SHARED_PTR<pdat::SideData<double> > data_mass_fractions_last;
        
        hier::IntVector offset_mass_fractions_last(d_dim);
        
        if (domain.empty())
        {
            const hier::Box interior_box = data_mass_fractions->getBox();
            
            offset_mass_fractions_last = data_mass_fractions->getGhostCellWidth();
            
            HAMERS_SHARED_PTR<pdat::SideData<double> > data_mass_fractions_last =
                HAMERS_MAKE_SHARED<pdat::SideData<double> >(
                    interior_box, 1, data_mass_fractions->getGhostCellWidth(), direction);
            
            data_mass_fractions_last->fillAll(double(1));
        }
        else
        {
            offset_mass_fractions_last = hier::IntVector::getZero(d_dim);
            
            data_mass_fractions_last =
                HAMERS_MAKE_SHARED<pdat::SideData<double> >(domain, 1, hier::IntVector::getZero(d_dim), direction);
            
            data_mass_fractions_last->fillAll(double(1), domain);
        }
        
        /*
         * Get the pointers to the cell data of mass fractions and the dimensions of the ghost cell box of
         * last mass fraction.
         */
        
        std::vector<const double*> Y;
        Y.reserve(d_num_species - 1);
        for (int si = 0; si < d_num_species - 1; si++)
        {
            Y.push_back(data_mass_fractions->getPointer(side_normal, si));
        }
        
        double* const Y_last = data_mass_fractions_last->getPointer(side_normal, 0);
        
        const hier::Box ghost_box_mass_fractions_last = data_mass_fractions_last->getGhostBox();
        hier::IntVector ghostcell_dims_mass_fractions_last = ghost_box_mass_fractions_last.numberCells();
        ghostcell_dims_mass_fractions_last[side_normal]++;
        
        computeIsobaricSpecificHeatCapacity(
            c_p,
            Y_last,
            Y,
            offset_isobaric_specific_heat_capacity,
            offset_mass_fractions_last,
            offset_mass_fractions,
            ghostcell_dims_isobaric_specific_heat_capacity,
            ghostcell_dims_mass_fractions_last,
            ghostcell_dims_mass_fractions,
            domain_lo,
            domain_dims);
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
 * Compute the Gruneisen parameter of the mixture with isothermal and isobaric equilibrium assumptions
 * (partial derivative of pressure w.r.t. specific internal energy under constant partial densities
 * divided by mixture density).
 */
double
EquationOfStateMixingRulesIdealGas::getGruneisenParameter(
    const double* const density,
    const double* const pressure,
    const std::vector<const double*>& mass_fractions) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT((d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC) ||
                (d_mixing_closure_model == MIXING_CLOSURE_MODEL::NO_MODEL && d_num_species == 1));
    TBOX_ASSERT((static_cast<int>(mass_fractions.size()) == d_num_species) ||
                (static_cast<int>(mass_fractions.size()) == d_num_species - 1));
#endif
    
    // Get the mixture thermodynamic properties.
    std::vector<double> mixture_thermo_properties;
    std::vector<double*> mixture_thermo_properties_ptr;
    std::vector<const double*> mixture_thermo_properties_const_ptr;
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    mixture_thermo_properties.resize(num_thermo_properties);
    mixture_thermo_properties_ptr.reserve(num_thermo_properties);
    mixture_thermo_properties_const_ptr.reserve(num_thermo_properties);
    
    for (int ti = 0; ti < num_thermo_properties; ti++)
    {
        mixture_thermo_properties_ptr.push_back(&mixture_thermo_properties[ti]);
        mixture_thermo_properties_const_ptr.push_back(&mixture_thermo_properties[ti]);
    }
    
    getMixtureThermodynamicProperties(
        mixture_thermo_properties_ptr,
        mass_fractions);
    
    return d_equation_of_state->getGruneisenParameter(
        density,
        pressure,
        mixture_thermo_properties_const_ptr);
}


/*
 * Compute the Gruneisen parameter of the mixture with isothermal and isobaric equilibrium assumptions
 * (partial derivative of pressure w.r.t. specific internal energy under constant partial densities
 * divided by mixture density).
 */
void
EquationOfStateMixingRulesIdealGas::computeGruneisenParameter(
    HAMERS_SHARED_PTR<pdat::CellData<double> >& data_gruneisen_parameter,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_mass_fractions,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT((d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC) ||
                (d_mixing_closure_model == MIXING_CLOSURE_MODEL::NO_MODEL && d_num_species == 1));
    
    TBOX_ASSERT(data_gruneisen_parameter);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_mass_fractions);
    
    TBOX_ASSERT((data_mass_fractions->getDepth() == d_num_species) ||
                (data_mass_fractions->getDepth() == d_num_species - 1));
#endif
    
    /*
     * Get the mixture thermodyanmic properties.
     */
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    // Declare data container for mixture thermodyanmic properties.
    HAMERS_SHARED_PTR<pdat::CellData<double> > data_mixture_thermo_properties;
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_gruneisen_parameter = data_gruneisen_parameter->getGhostCellWidth();
        const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_gruneisen_parameter->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_pressure->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_mass_fractions->getBox().isSpatiallyEqual(interior_box));
#endif
        
        /*
         * Get the minimum number of ghost cells for mixture thermodynamic properties.
         */
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_gruneisen_parameter;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_mass_fractions, num_ghosts_min);
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::CellData<double> >(
            interior_box, num_thermo_properties, num_ghosts_min);
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_gruneisen_parameter->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_mass_fractions->getGhostBox().contains(domain));
#endif
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::CellData<double> >(
            domain, num_thermo_properties, hier::IntVector::getZero(d_dim));
    }
    
    computeMixtureThermodynamicProperties(
        data_mixture_thermo_properties,
        data_mass_fractions,
        domain);
    
    d_equation_of_state->computeGruneisenParameter(
        data_gruneisen_parameter,
        data_density,
        data_pressure,
        data_mixture_thermo_properties,
        domain);
}


/*
 * Compute the Gruneisen parameter of the mixture with isothermal and isobaric equilibrium assumptions
 * (partial derivative of pressure w.r.t. specific internal energy under constant partial densities
 * divided by mixture density).
 */
void
EquationOfStateMixingRulesIdealGas::computeGruneisenParameter(
    HAMERS_SHARED_PTR<pdat::SideData<double> >& data_gruneisen_parameter,
    const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
    const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_mass_fractions,
    int side_normal,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT((d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC) ||
                (d_mixing_closure_model == MIXING_CLOSURE_MODEL::NO_MODEL && d_num_species == 1));
    
    TBOX_ASSERT(data_gruneisen_parameter);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_mass_fractions);
    
    TBOX_ASSERT((data_mass_fractions->getDepth() == d_num_species) ||
                (data_mass_fractions->getDepth() == d_num_species - 1));
#endif
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(side_normal < d_dim.getValue());
    
    TBOX_ASSERT(data_gruneisen_parameter->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_density->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_pressure->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_mass_fractions->getDirectionVector()[side_normal] > 0);
#endif
    
    hier::IntVector direction = hier::IntVector::getZero(d_dim);
    direction[side_normal] = 1;
    
    /*
     * Get the mixture thermodyanmic properties.
     */
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    // Declare data container for mixture thermodyanmic properties.
    HAMERS_SHARED_PTR<pdat::SideData<double> > data_mixture_thermo_properties;
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_gruneisen_parameter = data_gruneisen_parameter->getGhostCellWidth();
        const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_gruneisen_parameter->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_pressure->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_mass_fractions->getBox().isSpatiallyEqual(interior_box));
#endif
        
        /*
         * Get the minimum number of ghost cells for mixture thermodynamic properties.
         */
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_gruneisen_parameter;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_mass_fractions, num_ghosts_min);
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::SideData<double> >(
            interior_box, num_thermo_properties, num_ghosts_min, direction);
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_gruneisen_parameter->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_mass_fractions->getGhostBox().contains(domain));
#endif
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::SideData<double> >(
            domain, num_thermo_properties, hier::IntVector::getZero(d_dim), direction);
    }
    
    computeMixtureThermodynamicProperties(
        data_mixture_thermo_properties,
        data_mass_fractions,
        side_normal,
        domain);
    
    d_equation_of_state->computeGruneisenParameter(
        data_gruneisen_parameter,
        data_density,
        data_pressure,
        data_mixture_thermo_properties,
        side_normal,
        domain);
}


/*
 * Compute the Gruneisen parameter of the mixture with isobaric equilibrium assumption
 * (partial derivative of pressure w.r.t. specific internal energy under constant partial densities
 * and volume fractions divided by mixture density).
 */
double
EquationOfStateMixingRulesIdealGas::getGruneisenParameter(
    const double* const density,
    const double* const pressure,
    const std::vector<const double*>& mass_fractions,
    const std::vector<const double*>& volume_fractions) const
{
    NULL_USE(mass_fractions);
    
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOBARIC);
    TBOX_ASSERT((static_cast<int>(volume_fractions.size()) == d_num_species) ||
                (static_cast<int>(volume_fractions.size()) == d_num_species - 1));
#endif
    
    // Get the mixture thermodynamic properties.
    std::vector<double> mixture_thermo_properties;
    std::vector<double*> mixture_thermo_properties_ptr;
    std::vector<const double*> mixture_thermo_properties_const_ptr;
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    mixture_thermo_properties.resize(num_thermo_properties);
    mixture_thermo_properties_ptr.reserve(num_thermo_properties);
    mixture_thermo_properties_const_ptr.reserve(num_thermo_properties);
    
    for (int ti = 0; ti < num_thermo_properties; ti++)
    {
        mixture_thermo_properties_ptr.push_back(&mixture_thermo_properties[ti]);
        mixture_thermo_properties_const_ptr.push_back(&mixture_thermo_properties[ti]);
    }
    
    getMixtureThermodynamicProperties(
        mixture_thermo_properties_ptr,
        volume_fractions);
    
    return d_equation_of_state->getGruneisenParameter(
        density,
        pressure,
        mixture_thermo_properties_const_ptr);
}


/*
 * Compute the Gruneisen parameter of the mixture with isobaric equilibrium assumption
 * (partial derivative of pressure w.r.t. specific internal energy under constant partial densities
 * and volume fractions divided by mixture density).
 */
void
EquationOfStateMixingRulesIdealGas::computeGruneisenParameter(
    HAMERS_SHARED_PTR<pdat::CellData<double> >& data_gruneisen_parameter,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_mass_fractions,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_volume_fractions,
    const hier::Box& domain) const
{
    NULL_USE(data_mass_fractions);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOBARIC);
    
    TBOX_ASSERT(data_gruneisen_parameter);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_volume_fractions);
    
    TBOX_ASSERT((data_volume_fractions->getDepth() == d_num_species) ||
                (data_volume_fractions->getDepth() == d_num_species - 1));
#endif
    
    /*
     * Get the mixture thermodyanmic properties.
     */
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    // Declare data container for mixture thermodyanmic properties.
    HAMERS_SHARED_PTR<pdat::CellData<double> > data_mixture_thermo_properties;
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_gruneisen_parameter = data_gruneisen_parameter->getGhostCellWidth();
        const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        const hier::IntVector num_ghosts_volume_fractions = data_volume_fractions->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_gruneisen_parameter->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_pressure->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_volume_fractions->getBox().isSpatiallyEqual(interior_box));
#endif
        
        /*
         * Get the minimum number of ghost cells for mixture thermodynamic properties.
         */
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_gruneisen_parameter;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_volume_fractions, num_ghosts_min);
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::CellData<double> >(
            interior_box, num_thermo_properties, num_ghosts_min);
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_gruneisen_parameter->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_volume_fractions->getGhostBox().contains(domain));
#endif
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::CellData<double> >(
            domain, num_thermo_properties, hier::IntVector::getZero(d_dim));
    }
    
    computeMixtureThermodynamicProperties(
        data_mixture_thermo_properties,
        data_volume_fractions,
        domain);
    
    d_equation_of_state->computeGruneisenParameter(
        data_gruneisen_parameter,
        data_density,
        data_pressure,
        data_mixture_thermo_properties,
        domain);
}


/*
 * Compute the Gruneisen parameter of the mixture with isobaric equilibrium assumption
 * (partial derivative of pressure w.r.t. specific internal energy under constant partial densities
 * and volume fractions divided by mixture density).
 */
void
EquationOfStateMixingRulesIdealGas::computeGruneisenParameter(
    HAMERS_SHARED_PTR<pdat::SideData<double> >& data_gruneisen_parameter,
    const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
    const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_mass_fractions,
    const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_volume_fractions,
    int side_normal,
    const hier::Box& domain) const
{
    NULL_USE(data_mass_fractions);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOBARIC);
    
    TBOX_ASSERT(data_gruneisen_parameter);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_volume_fractions);
    
    TBOX_ASSERT((data_volume_fractions->getDepth() == d_num_species) ||
                (data_volume_fractions->getDepth() == d_num_species - 1));
#endif
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(side_normal < d_dim.getValue());
    
    TBOX_ASSERT(data_gruneisen_parameter->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_density->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_pressure->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_volume_fractions->getDirectionVector()[side_normal] > 0);
#endif
    
    hier::IntVector direction = hier::IntVector::getZero(d_dim);
    direction[side_normal] = 1;
    
    /*
     * Get the mixture thermodyanmic properties.
     */
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    // Declare data container for mixture thermodyanmic properties.
    HAMERS_SHARED_PTR<pdat::SideData<double> > data_mixture_thermo_properties;
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_gruneisen_parameter = data_gruneisen_parameter->getGhostCellWidth();
        const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        const hier::IntVector num_ghosts_volume_fractions = data_volume_fractions->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_gruneisen_parameter->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_pressure->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_volume_fractions->getBox().isSpatiallyEqual(interior_box));
#endif
        
        /*
         * Get the minimum number of ghost cells for mixture thermodynamic properties.
         */
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_gruneisen_parameter;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_volume_fractions, num_ghosts_min);
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::SideData<double> >(
            interior_box, num_thermo_properties, num_ghosts_min, direction);
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_gruneisen_parameter->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_volume_fractions->getGhostBox().contains(domain));
#endif
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::SideData<double> >(
            domain, num_thermo_properties, hier::IntVector::getZero(d_dim), direction);
    }
    
    computeMixtureThermodynamicProperties(
        data_mixture_thermo_properties,
        data_volume_fractions,
        side_normal,
        domain);
    
    d_equation_of_state->computeGruneisenParameter(
        data_gruneisen_parameter,
        data_density,
        data_pressure,
        data_mixture_thermo_properties,
        side_normal,
        domain);
}


/*
 * Compute the mixture partial derivative of pressure w.r.t. partial densities under constant specific
 * internal energy with isothermal and isobaric equilibrium assumptions.
 */
std::vector<double>
EquationOfStateMixingRulesIdealGas::getPressureDerivativeWithPartialDensities(
    const double* const density,
    const double* const pressure,
    const std::vector<const double*>& mass_fractions) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT((d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC) ||
                (d_mixing_closure_model == MIXING_CLOSURE_MODEL::NO_MODEL && d_num_species == 1));
    TBOX_ASSERT((static_cast<int>(mass_fractions.size()) == d_num_species) ||
                (static_cast<int>(mass_fractions.size()) == d_num_species - 1));
#endif
    
    // Get the mixture thermodynamic properties.
    std::vector<double> mixture_thermo_properties;
    std::vector<double*> mixture_thermo_properties_ptr;
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    mixture_thermo_properties.resize(num_thermo_properties);
    mixture_thermo_properties_ptr.reserve(num_thermo_properties);
    
    for (int ti = 0; ti < num_thermo_properties; ti++)
    {
        mixture_thermo_properties_ptr.push_back(&mixture_thermo_properties[ti]);
    }
    
    getMixtureThermodynamicProperties(
        mixture_thermo_properties_ptr,
        mass_fractions);
    
    const double& gamma = mixture_thermo_properties[0];
    const double& c_v   = mixture_thermo_properties[3];
    
    const double& rho = *density;
    const double& p   = *pressure;
    
    const double epsilon = p/((gamma - double(1))*rho);
    
    std::vector<double> Psi;
    Psi.reserve(d_num_species);
    for (int si = 0; si < d_num_species; si++)
    {
        double Psi_i = ((d_species_c_p[si] - gamma*d_species_c_v[si])/c_v + gamma - double(1))*epsilon;
        
        Psi.push_back(Psi_i);
    }
    
    return Psi;
}


/*
 * Compute the mixture partial derivative of pressure w.r.t. partial densities under constant specific
 * internal energy with isothermal and isobaric equilibrium assumptions.
 */
void
EquationOfStateMixingRulesIdealGas::computePressureDerivativeWithPartialDensities(
    HAMERS_SHARED_PTR<pdat::CellData<double> >& data_partial_pressure_partial_partial_densities,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_mass_fractions,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT((d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC) ||
                (d_mixing_closure_model == MIXING_CLOSURE_MODEL::NO_MODEL && d_num_species == 1));
    
    TBOX_ASSERT(data_partial_pressure_partial_partial_densities);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_mass_fractions);
    
    TBOX_ASSERT(data_partial_pressure_partial_partial_densities->getDepth() == d_num_species);
    
    TBOX_ASSERT((data_mass_fractions->getDepth() == d_num_species) ||
                (data_mass_fractions->getDepth() == d_num_species - 1));
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_partial_pressure_partial_partial_densities =
        data_partial_pressure_partial_partial_densities->getGhostBox();
    const hier::IntVector ghostcell_dims_partial_pressure_partial_partial_densities =
        ghost_box_partial_pressure_partial_partial_densities.numberCells();
    
    hier::IntVector ghostcell_dims_min(d_dim);
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets and allocate memory for the mixture thermodyanmic properties.
     * and the specific internal energy.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_partial_pressure_partial_partial_densities(d_dim);
    hier::IntVector offset_min(d_dim);
    
    HAMERS_SHARED_PTR<pdat::CellData<double> > data_mixture_thermo_properties;
    HAMERS_SHARED_PTR<pdat::CellData<double> > data_internal_energy;
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_partial_pressure_partial_partial_densities =
            data_partial_pressure_partial_partial_densities->getGhostCellWidth();
        
        const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_partial_pressure_partial_partial_densities->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_pressure->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_mass_fractions->getBox().isSpatiallyEqual(interior_box));
#endif
        
        /*
         * Get the minimum number of ghost cells and the dimension of the ghost cell box for
         * mixture thermodynamic properties and specific internal energy.
         */
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_partial_pressure_partial_partial_densities;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_mass_fractions, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        ghostcell_dims_min = ghost_box.numberCells();
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_partial_pressure_partial_partial_densities = num_ghosts_partial_pressure_partial_partial_densities;
        offset_min = num_ghosts_min;
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::CellData<double> >(
            interior_box, num_thermo_properties, num_ghosts_min);
        
        data_internal_energy = HAMERS_MAKE_SHARED<pdat::CellData<double> >(
            interior_box, 1, num_ghosts_min);
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
        TBOX_ASSERT(data_partial_pressure_partial_partial_densities->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_mass_fractions->getGhostBox().contains(domain));
#endif
        
        ghostcell_dims_min = domain.numberCells();
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_partial_pressure_partial_partial_densities =
            domain.lower() - ghost_box_partial_pressure_partial_partial_densities.lower();
        offset_min = hier::IntVector::getZero(d_dim);
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::CellData<double> >(
            domain, num_thermo_properties, hier::IntVector::getZero(d_dim));
        
        data_internal_energy = HAMERS_MAKE_SHARED<pdat::CellData<double> >(
            domain, 1, hier::IntVector::getZero(d_dim));
    }
    
    /*
     * Compute the mixture thermodyanmic properties and the specific internal energy.
     */
    
    computeMixtureThermodynamicProperties(
        data_mixture_thermo_properties,
        data_mass_fractions,
        domain);
    
    d_equation_of_state->computeInternalEnergy(
        data_internal_energy,
        data_density,
        data_pressure,
        data_mixture_thermo_properties,
        domain);
    
    /*
     * Get the pointers to the cell data.
     */
    
    double* const epsilon = data_internal_energy->getPointer(0);
    double* const gamma   = data_mixture_thermo_properties->getPointer(0);
    double* const c_v     = data_mixture_thermo_properties->getPointer(3);
    
    /*
     * Get the partial derivative.
     */
    
    std::vector<double*> Psi;
    Psi.reserve(d_num_species);
    for (int si = 0; si < d_num_species; si++)
    {
        Psi.push_back(data_partial_pressure_partial_partial_densities->getPointer(si));
    }
    
    computePressureDerivativeWithPartialDensities(
        Psi,
        epsilon,
        gamma,
        c_v,
        offset_partial_pressure_partial_partial_densities,
        offset_min,
        offset_min,
        ghostcell_dims_partial_pressure_partial_partial_densities,
        ghostcell_dims_min,
        ghostcell_dims_min,
        domain_lo,
        domain_dims);
}


/*
 * Compute the mixture partial derivative of pressure w.r.t. partial densities under constant specific
 * internal energy with isothermal and isobaric equilibrium assumptions.
 */
void
EquationOfStateMixingRulesIdealGas::computePressureDerivativeWithPartialDensities(
    HAMERS_SHARED_PTR<pdat::SideData<double> >& data_partial_pressure_partial_partial_densities,
    const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
    const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_mass_fractions,
    int side_normal,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT((d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC) ||
                (d_mixing_closure_model == MIXING_CLOSURE_MODEL::NO_MODEL && d_num_species == 1));
    
    TBOX_ASSERT(data_partial_pressure_partial_partial_densities);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_mass_fractions);
    
    TBOX_ASSERT(data_partial_pressure_partial_partial_densities->getDepth() == d_num_species);
    
    TBOX_ASSERT((data_mass_fractions->getDepth() == d_num_species) ||
                (data_mass_fractions->getDepth() == d_num_species - 1));
#endif
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(side_normal < d_dim.getValue());
    
    TBOX_ASSERT(data_partial_pressure_partial_partial_densities->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_density->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_pressure->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_mass_fractions->getDirectionVector()[side_normal] > 0);
#endif
    
    hier::IntVector direction = hier::IntVector::getZero(d_dim);
    direction[side_normal] = 1;
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_partial_pressure_partial_partial_densities =
        data_partial_pressure_partial_partial_densities->getGhostBox();
    hier::IntVector ghostcell_dims_partial_pressure_partial_partial_densities =
        ghost_box_partial_pressure_partial_partial_densities.numberCells();
    
    hier::IntVector ghostcell_dims_min(d_dim);
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets and allocate memory for the mixture thermodyanmic properties.
     * and the specific internal energy.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_partial_pressure_partial_partial_densities(d_dim);
    hier::IntVector offset_min(d_dim);
    
    HAMERS_SHARED_PTR<pdat::SideData<double> > data_mixture_thermo_properties;
    HAMERS_SHARED_PTR<pdat::SideData<double> > data_internal_energy;
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_partial_pressure_partial_partial_densities =
            data_partial_pressure_partial_partial_densities->getGhostCellWidth();
        
        const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_partial_pressure_partial_partial_densities->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_pressure->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_mass_fractions->getBox().isSpatiallyEqual(interior_box));
#endif
        
        /*
         * Get the minimum number of ghost cells and the dimension of the ghost cell box for
         * mixture thermodynamic properties and specific internal energy.
         */
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_partial_pressure_partial_partial_densities;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_mass_fractions, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        ghostcell_dims_min = ghost_box.numberCells();
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_partial_pressure_partial_partial_densities = num_ghosts_partial_pressure_partial_partial_densities;
        offset_min = num_ghosts_min;
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::SideData<double> >(
            interior_box, num_thermo_properties, num_ghosts_min, direction);
        
        data_internal_energy = HAMERS_MAKE_SHARED<pdat::SideData<double> >(
            interior_box, 1, num_ghosts_min, direction);
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
        TBOX_ASSERT(data_partial_pressure_partial_partial_densities->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_mass_fractions->getGhostBox().contains(domain));
#endif
        
        ghostcell_dims_min = domain.numberCells();
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_partial_pressure_partial_partial_densities =
            domain.lower() - ghost_box_partial_pressure_partial_partial_densities.lower();
        offset_min = hier::IntVector::getZero(d_dim);
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::SideData<double> >(
                domain, num_thermo_properties, hier::IntVector::getZero(d_dim), direction);
        
        data_internal_energy = HAMERS_MAKE_SHARED<pdat::SideData<double> >(
            domain, 1, hier::IntVector::getZero(d_dim), direction);
    }
    
    ghostcell_dims_partial_pressure_partial_partial_densities[side_normal]++;
    ghostcell_dims_min[side_normal]++;
    domain_dims[side_normal]++;
    
    /*
     * Compute the mixture thermodyanmic properties and the specific internal energy.
     */
    
    computeMixtureThermodynamicProperties(
        data_mixture_thermo_properties,
        data_mass_fractions,
        side_normal,
        domain);
    
    d_equation_of_state->computeInternalEnergy(
        data_internal_energy,
        data_density,
        data_pressure,
        data_mixture_thermo_properties,
        side_normal,
        domain);
    
    /*
     * Get the pointers to the cell data.
     */
    
    double* const epsilon = data_internal_energy->getPointer(side_normal, 0);
    double* const gamma   = data_mixture_thermo_properties->getPointer(side_normal, 0);
    double* const c_v     = data_mixture_thermo_properties->getPointer(side_normal, 3);
    
    /*
     * Get the partial derivative.
     */
    
    std::vector<double*> Psi;
    Psi.reserve(d_num_species);
    for (int si = 0; si < d_num_species; si++)
    {
        Psi.push_back(data_partial_pressure_partial_partial_densities->getPointer(side_normal, si));
    }
    
    computePressureDerivativeWithPartialDensities(
        Psi,
        epsilon,
        gamma,
        c_v,
        offset_partial_pressure_partial_partial_densities,
        offset_min,
        offset_min,
        ghostcell_dims_partial_pressure_partial_partial_densities,
        ghostcell_dims_min,
        ghostcell_dims_min,
        domain_lo,
        domain_dims);
}


/*
 * Compute the mixture partial derivative of pressure w.r.t. partial densities under constant specific
 * internal energy and volume fractions with isobaric equilibrium assumption.
 */
std::vector<double>
EquationOfStateMixingRulesIdealGas::getPressureDerivativeWithPartialDensities(
    const double* const density,
    const double* const pressure,
    const std::vector<const double*>& mass_fractions,
    const std::vector<const double*>& volume_fractions) const
{
    NULL_USE(mass_fractions);
    NULL_USE(volume_fractions);
    
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOBARIC);
    TBOX_ASSERT((static_cast<int>(volume_fractions.size()) == d_num_species) ||
                (static_cast<int>(volume_fractions.size()) == d_num_species - 1));
#endif
    
    const double& rho = *density;
    const double& p   = *pressure;
    
    double Psi_i = p/rho;
    
    std::vector<double> Psi;
    Psi.reserve(d_num_species);
    for (int si = 0; si < d_num_species; si++)
    {
        Psi.push_back(Psi_i);
    }
    
    return Psi;
}


/*
 * Compute the mixture partial derivative of pressure w.r.t. partial densities under constant specific
 * internal energy and volume fractions with isobaric equilibrium assumption.
 */
void
EquationOfStateMixingRulesIdealGas::computePressureDerivativeWithPartialDensities(
    HAMERS_SHARED_PTR<pdat::CellData<double> >& data_partial_pressure_partial_partial_densities,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_mass_fractions,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_volume_fractions,
    const hier::Box& domain) const
{
    NULL_USE(data_mass_fractions);
    NULL_USE(data_volume_fractions);
    
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOBARIC);
    
    TBOX_ASSERT(data_partial_pressure_partial_partial_densities);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    
    TBOX_ASSERT(data_partial_pressure_partial_partial_densities->getDepth() == d_num_species);
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_partial_pressure_partial_partial_densities =
        data_partial_pressure_partial_partial_densities->getGhostBox();
    const hier::IntVector ghostcell_dims_partial_pressure_partial_partial_densities =
        ghost_box_partial_pressure_partial_partial_densities.numberCells();
    
    const hier::Box ghost_box_density = data_density->getGhostBox();
    const hier::IntVector ghostcell_dims_density = ghost_box_density.numberCells();
    
    const hier::Box ghost_box_pressure = data_pressure->getGhostBox();
    const hier::IntVector ghostcell_dims_pressure = data_pressure->getGhostBox().numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_partial_pressure_partial_partial_densities(d_dim);
    hier::IntVector offset_density(d_dim);
    hier::IntVector offset_pressure(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_partial_pressure_partial_partial_densities =
            data_partial_pressure_partial_partial_densities->getGhostCellWidth();
        
        const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_partial_pressure_partial_partial_densities->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_pressure->getBox().isSpatiallyEqual(interior_box));
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_partial_pressure_partial_partial_densities;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_partial_pressure_partial_partial_densities = num_ghosts_partial_pressure_partial_partial_densities;
        offset_density = num_ghosts_density;
        offset_pressure = num_ghosts_pressure;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
        TBOX_ASSERT(data_partial_pressure_partial_partial_densities->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_partial_pressure_partial_partial_densities = domain.lower() -
            ghost_box_partial_pressure_partial_partial_densities.lower();
        
        offset_density = domain.lower() - ghost_box_density.lower();
        offset_pressure = domain.lower() - ghost_box_pressure.lower();
    }
    
    /*
     * Get the pointers to the cell data.
     */
    
    double* const rho = data_density->getPointer(0);
    double* const p   = data_pressure->getPointer(0);
    
    /*
     * Get the partial derivative.
     */
    
    std::vector<double*> Psi;
    Psi.reserve(d_num_species);
    for (int si = 0; si < d_num_species; si++)
    {
        Psi.push_back(data_partial_pressure_partial_partial_densities->getPointer(si));
    }
    
    computePressureDerivativeWithPartialDensities(
        Psi,
        rho,
        p,
        offset_partial_pressure_partial_partial_densities,
        offset_density,
        offset_pressure,
        ghostcell_dims_partial_pressure_partial_partial_densities,
        ghostcell_dims_density,
        ghostcell_dims_pressure,
        domain_lo,
        domain_dims);
}


/*
 * Compute the mixture partial derivative of pressure w.r.t. partial densities under constant specific
 * internal energy and volume fractions with isobaric equilibrium assumption.
 */
void
EquationOfStateMixingRulesIdealGas::computePressureDerivativeWithPartialDensities(
    HAMERS_SHARED_PTR<pdat::SideData<double> >& data_partial_pressure_partial_partial_densities,
    const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
    const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_mass_fractions,
    const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_volume_fractions,
    int side_normal,
    const hier::Box& domain) const
{
    NULL_USE(data_mass_fractions);
    NULL_USE(data_volume_fractions);
    
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOBARIC);
    
    TBOX_ASSERT(data_partial_pressure_partial_partial_densities);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    
    TBOX_ASSERT(data_partial_pressure_partial_partial_densities->getDepth() == d_num_species);
#endif
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(side_normal < d_dim.getValue());
    
    TBOX_ASSERT(data_partial_pressure_partial_partial_densities->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_density->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_pressure->getDirectionVector()[side_normal] > 0);
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_partial_pressure_partial_partial_densities =
        data_partial_pressure_partial_partial_densities->getGhostBox();
    hier::IntVector ghostcell_dims_partial_pressure_partial_partial_densities =
        ghost_box_partial_pressure_partial_partial_densities.numberCells();
    
    const hier::Box ghost_box_density = data_density->getGhostBox();
    hier::IntVector ghostcell_dims_density = ghost_box_density.numberCells();
    
    const hier::Box ghost_box_pressure = data_pressure->getGhostBox();
    hier::IntVector ghostcell_dims_pressure = ghost_box_pressure.numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_partial_pressure_partial_partial_densities(d_dim);
    hier::IntVector offset_density(d_dim);
    hier::IntVector offset_pressure(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_partial_pressure_partial_partial_densities =
            data_partial_pressure_partial_partial_densities->getGhostCellWidth();
        
        const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_partial_pressure_partial_partial_densities->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_pressure->getBox().isSpatiallyEqual(interior_box));
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_partial_pressure_partial_partial_densities;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_partial_pressure_partial_partial_densities = num_ghosts_partial_pressure_partial_partial_densities;
        offset_density = num_ghosts_density;
        offset_pressure = num_ghosts_pressure;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
        TBOX_ASSERT(data_partial_pressure_partial_partial_densities->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_partial_pressure_partial_partial_densities =
            domain.lower() - ghost_box_partial_pressure_partial_partial_densities.lower();
        
        offset_density = domain.lower() - ghost_box_density.lower();
        offset_pressure = domain.lower() - ghost_box_pressure.lower();
    }
    
    ghostcell_dims_partial_pressure_partial_partial_densities[side_normal]++;
    ghostcell_dims_density[side_normal]++;
    ghostcell_dims_pressure[side_normal]++;
    domain_dims[side_normal]++;
    
    /*
     * Get the pointers to the cell data.
     */
    
    double* const rho = data_density->getPointer(side_normal, 0);
    double* const p   = data_pressure->getPointer(side_normal, 0);
    
    /*
     * Get the partial derivative.
     */
    
    std::vector<double*> Psi;
    Psi.reserve(d_num_species);
    for (int si = 0; si < d_num_species; si++)
    {
        Psi.push_back(data_partial_pressure_partial_partial_densities->getPointer(side_normal, si));
    }
    
    computePressureDerivativeWithPartialDensities(
        Psi,
        rho,
        p,
        offset_partial_pressure_partial_partial_densities,
        offset_density,
        offset_pressure,
        ghostcell_dims_partial_pressure_partial_partial_densities,
        ghostcell_dims_density,
        ghostcell_dims_pressure,
        domain_lo,
        domain_dims);
}


/*
 * Compute the mixture partial derivative of pressure w.r.t. volume fractions under constant specific
 * internal energy and partial densities with isobaric equilibrium assumption.
 */
std::vector<double>
EquationOfStateMixingRulesIdealGas::getPressureDerivativeWithVolumeFractions(
    const double* const density,
    const double* const pressure,
    const std::vector<const double*>& mass_fractions,
    const std::vector<const double*>& volume_fractions) const
{
    NULL_USE(density);
    NULL_USE(mass_fractions);
    
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOBARIC);
    TBOX_ASSERT((static_cast<int>(volume_fractions.size()) == d_num_species) ||
                (static_cast<int>(volume_fractions.size()) == d_num_species - 1));
#endif
    
    // Get the mixture thermodynamic properties.
    std::vector<double> mixture_thermo_properties;
    std::vector<double*> mixture_thermo_properties_ptr;
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    mixture_thermo_properties.resize(num_thermo_properties);
    mixture_thermo_properties_ptr.reserve(num_thermo_properties);
    
    for (int ti = 0; ti < num_thermo_properties; ti++)
    {
        mixture_thermo_properties_ptr.push_back(&mixture_thermo_properties[ti]);
    }
    
    getMixtureThermodynamicProperties(
        mixture_thermo_properties_ptr,
        volume_fractions);
    
    const double& gamma = mixture_thermo_properties[0];
    
    const double& p   = *pressure;
    
    const double tmp = double(1)/(d_species_gamma[d_num_species - 1] - double(1));
    
    std::vector<double> M;
    M.reserve(d_num_species - 1);
    for (int si = 0; si < d_num_species - 1; si++)
    {
        double M_i = (tmp - double(1)/(d_species_gamma[si] - double(1)))*(gamma - double(1))*p;
        M.push_back(M_i);
    }
    
    return M;
}


/*
 * Compute the mixture partial derivative of pressure w.r.t. volume fractions under constant specific
 * internal energy and partial densities with isobaric equilibrium assumption.
 */
void
EquationOfStateMixingRulesIdealGas::computePressureDerivativeWithVolumeFractions(
    HAMERS_SHARED_PTR<pdat::CellData<double> >& data_partial_pressure_partial_volume_fractions,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_density,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_mass_fractions,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_volume_fractions,
    const hier::Box& domain) const
{
    NULL_USE(data_density);
    NULL_USE(data_mass_fractions);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOBARIC);
    
    TBOX_ASSERT(data_partial_pressure_partial_volume_fractions);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_volume_fractions);
    
    TBOX_ASSERT(data_partial_pressure_partial_volume_fractions->getDepth() == d_num_species - 1);
    
    TBOX_ASSERT((data_volume_fractions->getDepth() == d_num_species) ||
                (data_volume_fractions->getDepth() == d_num_species - 1));
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_partial_pressure_partial_volume_fractions =
        data_partial_pressure_partial_volume_fractions->getGhostBox();
    const hier::IntVector ghostcell_dims_partial_pressure_partial_volume_fractions =
        ghost_box_partial_pressure_partial_volume_fractions.numberCells();
    
    const hier::Box ghost_box_pressure = data_pressure->getGhostBox();
    const hier::IntVector ghostcell_dims_pressure = ghost_box_pressure.numberCells();
    
    hier::IntVector ghostcell_dims_min(d_dim);
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets and allocate memory for the mixture thermodyanmic properties.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_partial_pressure_partial_volume_fractions(d_dim);
    hier::IntVector offset_pressure(d_dim);
    hier::IntVector offset_min(d_dim);
    
    HAMERS_SHARED_PTR<pdat::CellData<double> > data_mixture_thermo_properties;
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_partial_pressure_partial_volume_fractions =
            data_partial_pressure_partial_volume_fractions->getGhostCellWidth();
        
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        const hier::IntVector num_ghosts_volume_fractions = data_volume_fractions->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_partial_pressure_partial_volume_fractions->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_pressure->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_volume_fractions->getBox().isSpatiallyEqual(interior_box));
#endif
        
        /*
         * Get the minimum number of ghost cells and the dimension of the ghost cell box for
         * mixture thermodynamic properties.
         */
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_partial_pressure_partial_volume_fractions;
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_volume_fractions, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        ghostcell_dims_min = ghost_box.numberCells();
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_partial_pressure_partial_volume_fractions = num_ghosts_partial_pressure_partial_volume_fractions;
        offset_pressure = num_ghosts_pressure;
        offset_min = num_ghosts_min;
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::CellData<double> >(
            interior_box, num_thermo_properties, num_ghosts_min);
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
        TBOX_ASSERT(data_partial_pressure_partial_volume_fractions->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_volume_fractions->getGhostBox().contains(domain));
#endif
        
        ghostcell_dims_min = domain.numberCells();
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_partial_pressure_partial_volume_fractions =
            domain.lower() - ghost_box_partial_pressure_partial_volume_fractions.lower();
        
        offset_pressure = domain.lower() - ghost_box_pressure.lower();
        offset_min = hier::IntVector::getZero(d_dim);
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::CellData<double> >(
            domain, num_thermo_properties, hier::IntVector::getZero(d_dim));
    }
    
    // Compute the mixture thermodyanmic properties.
    computeMixtureThermodynamicProperties(
        data_mixture_thermo_properties,
        data_volume_fractions,
        domain);
    
    /*
     * Get the pointers to the cell data.
     */
    
    double* const p     = data_pressure->getPointer(0);
    double* const gamma = data_mixture_thermo_properties->getPointer(0);
    
    /*
     * Get the partial derivative.
     */
    
    std::vector<double*> M;
    M.reserve(d_num_species - 1);
    for (int si = 0; si < d_num_species - 1; si++)
    {
        M.push_back(data_partial_pressure_partial_volume_fractions->getPointer(si));
    }
    
    computePressureDerivativeWithVolumeFractions(
        M,
        p,
        gamma,
        offset_partial_pressure_partial_volume_fractions,
        offset_pressure,
        offset_min,
        ghostcell_dims_partial_pressure_partial_volume_fractions,
        ghostcell_dims_pressure,
        ghostcell_dims_min,
        domain_lo,
        domain_dims);
}


/*
 * Compute the mixture partial derivative of pressure w.r.t. volume fractions under constant specific
 * internal energy and partial densities with isobaric equilibrium assumption.
 */
void
EquationOfStateMixingRulesIdealGas::computePressureDerivativeWithVolumeFractions(
    HAMERS_SHARED_PTR<pdat::SideData<double> >& data_partial_pressure_partial_volume_fractions,
    const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_density,
    const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_mass_fractions,
    const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_volume_fractions,
    int side_normal,
    const hier::Box& domain) const
{
    NULL_USE(data_density);
    NULL_USE(data_mass_fractions);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOBARIC);
    
    TBOX_ASSERT(data_partial_pressure_partial_volume_fractions);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_volume_fractions);
    
    TBOX_ASSERT(data_partial_pressure_partial_volume_fractions->getDepth() == d_num_species - 1);
    
    TBOX_ASSERT((data_volume_fractions->getDepth() == d_num_species) ||
                (data_volume_fractions->getDepth() == d_num_species - 1));
#endif
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(side_normal < d_dim.getValue());
    
    TBOX_ASSERT(data_partial_pressure_partial_volume_fractions->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_pressure->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_volume_fractions->getDirectionVector()[side_normal] > 0);
#endif
    
    hier::IntVector direction = hier::IntVector::getZero(d_dim);
    direction[side_normal] = 1;
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_partial_pressure_partial_volume_fractions =
        data_partial_pressure_partial_volume_fractions->getGhostBox();
    hier::IntVector ghostcell_dims_partial_pressure_partial_volume_fractions =
        ghost_box_partial_pressure_partial_volume_fractions.numberCells();
    
    const hier::Box ghost_box_pressure = data_pressure->getGhostBox();
    hier::IntVector ghostcell_dims_pressure = ghost_box_pressure.numberCells();
    
    hier::IntVector ghostcell_dims_min(d_dim);
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets and allocate memory for the mixture thermodyanmic properties.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_partial_pressure_partial_volume_fractions(d_dim);
    hier::IntVector offset_pressure(d_dim);
    hier::IntVector offset_min(d_dim);
    
    HAMERS_SHARED_PTR<pdat::SideData<double> > data_mixture_thermo_properties;
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_partial_pressure_partial_volume_fractions =
            data_partial_pressure_partial_volume_fractions->getGhostCellWidth();
        
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        const hier::IntVector num_ghosts_volume_fractions = data_volume_fractions->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_partial_pressure_partial_volume_fractions->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_pressure->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_volume_fractions->getBox().isSpatiallyEqual(interior_box));
#endif
        
        /*
         * Get the minimum number of ghost cells and the dimension of the ghost cell box for
         * mixture thermodynamic properties.
         */
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_partial_pressure_partial_volume_fractions;
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_volume_fractions, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        ghostcell_dims_min = ghost_box.numberCells();
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_partial_pressure_partial_volume_fractions = num_ghosts_partial_pressure_partial_volume_fractions;
        offset_pressure = num_ghosts_pressure;
        offset_min = num_ghosts_min;
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::SideData<double> >(
            interior_box, num_thermo_properties, num_ghosts_min, direction);
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
        TBOX_ASSERT(data_partial_pressure_partial_volume_fractions->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_volume_fractions->getGhostBox().contains(domain));
#endif
        
        ghostcell_dims_min = domain.numberCells();
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_partial_pressure_partial_volume_fractions =
            domain.lower() - ghost_box_partial_pressure_partial_volume_fractions.lower();
        
        offset_pressure = domain.lower() - ghost_box_pressure.lower();
        offset_min = hier::IntVector::getZero(d_dim);
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::SideData<double> >(
            domain, num_thermo_properties, hier::IntVector::getZero(d_dim), direction);
    }
    
    ghostcell_dims_partial_pressure_partial_volume_fractions[side_normal]++;
    ghostcell_dims_pressure[side_normal]++;
    ghostcell_dims_min[side_normal]++;
    domain_dims[side_normal]++;
    
    // Compute the mixture thermodyanmic properties.
    computeMixtureThermodynamicProperties(
        data_mixture_thermo_properties,
        data_volume_fractions,
        side_normal,
        domain);
    
    /*
     * Get the pointers to the cell data.
     */
    
    double* const p     = data_pressure->getPointer(side_normal, 0);
    double* const gamma = data_mixture_thermo_properties->getPointer(side_normal, 0);
    
    /*
     * Get the partial derivative.
     */
    
    std::vector<double*> M;
    M.reserve(d_num_species - 1);
    for (int si = 0; si < d_num_species - 1; si++)
    {
        M.push_back(data_partial_pressure_partial_volume_fractions->getPointer(side_normal, si));
    }
    
    computePressureDerivativeWithVolumeFractions(
        M,
        p,
        gamma,
        offset_partial_pressure_partial_volume_fractions,
        offset_pressure,
        offset_min,
        ghostcell_dims_partial_pressure_partial_volume_fractions,
        ghostcell_dims_pressure,
        ghostcell_dims_min,
        domain_lo,
        domain_dims);
}


/*
 * Compute the density of mixture with isothermal and isobaric equilibrium assumptions.
 */
double
EquationOfStateMixingRulesIdealGas::getMixtureDensity(
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
    
    // Get the mixture thermodynamic properties.
    std::vector<double> mixture_thermo_properties;
    std::vector<double*> mixture_thermo_properties_ptr;
    std::vector<const double*> mixture_thermo_properties_const_ptr;
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    mixture_thermo_properties.resize(num_thermo_properties);
    mixture_thermo_properties_ptr.reserve(num_thermo_properties);
    mixture_thermo_properties_const_ptr.reserve(num_thermo_properties);
    
    for (int ti = 0; ti < num_thermo_properties; ti++)
    {
        mixture_thermo_properties_ptr.push_back(&mixture_thermo_properties[ti]);
        mixture_thermo_properties_const_ptr.push_back(&mixture_thermo_properties[ti]);
    }
    
    getMixtureThermodynamicProperties(
        mixture_thermo_properties_ptr,
        mass_fractions);
    
    return d_equation_of_state->getDensity(
        pressure,
        temperature,
        mixture_thermo_properties_const_ptr);
}


/*
 * Compute the density of mixture with isothermal and isobaric equilibrium assumptions.
 */
void
EquationOfStateMixingRulesIdealGas::computeMixtureDensity(
    HAMERS_SHARED_PTR<pdat::CellData<double> >& data_mixture_density,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_temperature,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_mass_fractions,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT((d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC) ||
                (d_mixing_closure_model == MIXING_CLOSURE_MODEL::NO_MODEL && d_num_species == 1));
    
    TBOX_ASSERT(data_mixture_density);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_temperature);
    TBOX_ASSERT(data_mass_fractions);
    
    TBOX_ASSERT((data_mass_fractions->getDepth() == d_num_species) ||
                (data_mass_fractions->getDepth() == d_num_species - 1));
#endif
    
    /*
     * Get the mixture thermodyanmic properties.
     */
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    // Declare data container for mixture thermodyanmic properties.
    HAMERS_SHARED_PTR<pdat::CellData<double> > data_mixture_thermo_properties;
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_mixture_density = data_mixture_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
        const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_mixture_density->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_pressure->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_temperature->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_mass_fractions->getBox().isSpatiallyEqual(interior_box));
#endif
        
        /*
         * Get the minimum number of ghost cells for mixture thermodynamic properties.
         */
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_mixture_density;
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_temperature, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_mass_fractions, num_ghosts_min);
        
        HAMERS_SHARED_PTR<pdat::CellData<double> > data_mixture_thermo_properties(
            new pdat::CellData<double>(interior_box, num_thermo_properties, num_ghosts_min));
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_mixture_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_temperature->getGhostBox().contains(domain));
        TBOX_ASSERT(data_mass_fractions->getGhostBox().contains(domain));
#endif
        HAMERS_SHARED_PTR<pdat::CellData<double> > data_mixture_thermo_properties(
            new pdat::CellData<double>(domain, num_thermo_properties, hier::IntVector::getZero(d_dim)));
    }
    
    computeMixtureThermodynamicProperties(
        data_mixture_thermo_properties,
        data_mass_fractions,
        domain);
    
    d_equation_of_state->computeDensity(
        data_mixture_density,
        data_pressure,
        data_temperature,
        data_mixture_thermo_properties,
        domain);
}


/*
 * Compute the density of mixture with isothermal and isobaric equilibrium assumptions.
 */
void
EquationOfStateMixingRulesIdealGas::computeMixtureDensity(
    HAMERS_SHARED_PTR<pdat::SideData<double> >& data_mixture_density,
    const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_temperature,
    const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_mass_fractions,
    int side_normal,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT((d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC) ||
                (d_mixing_closure_model == MIXING_CLOSURE_MODEL::NO_MODEL && d_num_species == 1));

    TBOX_ASSERT(data_mixture_density);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_temperature);
    TBOX_ASSERT(data_mass_fractions);
    
    TBOX_ASSERT((data_mass_fractions->getDepth() == d_num_species) ||
                (data_mass_fractions->getDepth() == d_num_species - 1));
#endif
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(side_normal < d_dim.getValue());
    
    TBOX_ASSERT(data_mixture_density->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_pressure->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_temperature->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_mass_fractions->getDirectionVector()[side_normal] > 0);
#endif
    
    hier::IntVector direction = hier::IntVector::getZero(d_dim);
    direction[side_normal] = 1;
    
    /*
     * Get the mixture thermodyanmic properties.
     */
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    HAMERS_SHARED_PTR<pdat::SideData<double> > data_mixture_thermo_properties;
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_mixture_density = data_mixture_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
        const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_mixture_density->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_pressure->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_temperature->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_mass_fractions->getBox().isSpatiallyEqual(interior_box));
#endif
        
        /*
         * Get the minimum number of ghost cells for mixture thermodynamic properties.
         */
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_mixture_density;
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_temperature, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_mass_fractions, num_ghosts_min);
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::SideData<double> >(
            interior_box, num_thermo_properties, num_ghosts_min, direction);
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_mixture_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_temperature->getGhostBox().contains(domain));
        TBOX_ASSERT(data_mass_fractions->getGhostBox().contains(domain));
#endif
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::SideData<double> >(
            domain, num_thermo_properties, hier::IntVector::getZero(d_dim), direction);
    }
    
    computeMixtureThermodynamicProperties(
        data_mixture_thermo_properties,
        data_mass_fractions,
        side_normal,
        domain);
    
    d_equation_of_state->computeDensity(
        data_mixture_density,
        data_pressure,
        data_temperature,
        data_mixture_thermo_properties,
        side_normal,
        domain);
}


/*
 * Get the thermodynamic properties of a species.
 */
void
EquationOfStateMixingRulesIdealGas::getSpeciesThermodynamicProperties(
    std::vector<double*>& species_thermo_properties,
    const int species_index) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(species_thermo_properties.size()) == 4);
    TBOX_ASSERT(species_index >= 0);
    TBOX_ASSERT(species_index < d_num_species);
#endif
    
    // Get references to gamma and R of species.
    double& gamma = *(species_thermo_properties[0]);
    double& R = *(species_thermo_properties[1]);
    
    // Get references to c_p and c_v of species.
    double& c_p = *(species_thermo_properties[2]);
    double& c_v = *(species_thermo_properties[3]);
    
    gamma = d_species_gamma[species_index];
    R = d_species_R[species_index];
    
    c_p = d_species_c_p[species_index];
    c_v = d_species_c_v[species_index];
}


/*
 * Get the number of thermodynamic properties of the mixture.
 */
int
EquationOfStateMixingRulesIdealGas::getNumberOfMixtureThermodynamicProperties() const
{
    int num_of_thermo_properties = 0;
    
    switch (d_mixing_closure_model)
    {
        case MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC:
        {
            num_of_thermo_properties = 4;
            
            break;
        }
        case MIXING_CLOSURE_MODEL::ISOBARIC:
        {
            num_of_thermo_properties = 1;
            
            break;
        }
        case MIXING_CLOSURE_MODEL::NO_MODEL:
        {
            num_of_thermo_properties = 0;
            
            break;
        }
    }
    
    return num_of_thermo_properties;
}


/*
 * Get the thermodynamic properties of the mixture.
 */
void
EquationOfStateMixingRulesIdealGas::getMixtureThermodynamicProperties(
    std::vector<double*>& mixture_thermo_properties,
    const std::vector<const double*>& species_fraction) const
{
    switch (d_mixing_closure_model)
    {
        case MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC:
        {
            getMixtureThermodynamicPropertiesWithMassFractions(
                mixture_thermo_properties,
                species_fraction);
            
            break;
        }
        case MIXING_CLOSURE_MODEL::ISOBARIC:
        {
            getMixtureThermodynamicPropertiesWithVolumeFractions(
                mixture_thermo_properties,
                species_fraction);
            
            break;
        }
        case MIXING_CLOSURE_MODEL::NO_MODEL:
        {
            TBOX_WARNING(d_object_name
                << ": "
                << "No thermodynamic properties of mixture are returned!\n"
                << "d_mixing_closure_model = "
                << d_mixing_closure_model
                << std::endl);
            
            break;
        }
    }
}


/*
 * Compute the thermodynamic properties of the mixture.
 */
void
EquationOfStateMixingRulesIdealGas::computeMixtureThermodynamicProperties(
    HAMERS_SHARED_PTR<pdat::CellData<double> >& data_mixture_thermo_properties,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_species_fraction,
    const hier::Box& domain) const
{
    switch (d_mixing_closure_model)
    {
        case MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC:
        {
            computeMixtureThermodynamicPropertiesWithMassFractions(
                data_mixture_thermo_properties,
                data_species_fraction,
                domain);
            
            break;
        }
        case MIXING_CLOSURE_MODEL::ISOBARIC:
        {
            computeMixtureThermodynamicPropertiesWithVolumeFractions(
                data_mixture_thermo_properties,
                data_species_fraction,
                domain);
            
            break;
        }
        case MIXING_CLOSURE_MODEL::NO_MODEL:
        {
            TBOX_WARNING(d_object_name
                << ": "
                << "No thermodynamic properties of mixture are returned!\n"
                << "d_mixing_closure_model = "
                << d_mixing_closure_model
                << std::endl);
            
            break;
        }
    }
}


/*
 * Get the thermodynamic properties of the mixture.
 */
void
EquationOfStateMixingRulesIdealGas::computeMixtureThermodynamicProperties(
    HAMERS_SHARED_PTR<pdat::SideData<double> >& data_mixture_thermo_properties,
    const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_species_fraction,
    int side_normal,
    const hier::Box& domain) const
{
    switch (d_mixing_closure_model)
    {
        case MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC:
        {
            computeMixtureThermodynamicPropertiesWithMassFractions(
                data_mixture_thermo_properties,
                data_species_fraction,
                side_normal,
                domain);
            
            break;
        }
        case MIXING_CLOSURE_MODEL::ISOBARIC:
        {
            computeMixtureThermodynamicPropertiesWithVolumeFractions(
                data_mixture_thermo_properties,
                data_species_fraction,
                side_normal,
                domain);
            
            break;
        }
        case MIXING_CLOSURE_MODEL::NO_MODEL:
        {
            TBOX_WARNING(d_object_name
                << ": "
                << "No thermodynamic properties of mixture are returned!\n"
                << "d_mixing_closure_model = "
                << d_mixing_closure_model
                << std::endl);
            
            break;
        }
    }
}


/*
 * Compute the thermodynamic properties of the mixture with mass fractions.
 */
void
EquationOfStateMixingRulesIdealGas::getMixtureThermodynamicPropertiesWithMassFractions(
    std::vector<double*>& mixture_thermo_properties,
    const std::vector<const double*>& mass_fractions) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(mixture_thermo_properties.size()) == 4);
#endif
    
    // Get references to gamma and R of mixture.
    double& gamma = *(mixture_thermo_properties[0]);
    double& R = *(mixture_thermo_properties[1]);
    
    // Get references to c_p and c_v of mixture.
    double& c_p = *(mixture_thermo_properties[2]);
    double& c_v = *(mixture_thermo_properties[3]);
    
    c_p = double(0);
    c_v = double(0);
    
    if (static_cast<int>(mass_fractions.size()) == d_num_species)
    {
        for (int si = 0; si < d_num_species; si++)
        {
            c_p += *(mass_fractions[si])*d_species_c_p[si];
            c_v += *(mass_fractions[si])*d_species_c_v[si];
        }
    }
    else if (static_cast<int>(mass_fractions.size()) == d_num_species - 1)
    {
        double Y_last = double(1);
        
        for (int si = 0; si < d_num_species - 1; si++)
        {
            c_p += *(mass_fractions[si])*d_species_c_p[si];
            c_v += *(mass_fractions[si])*d_species_c_v[si];
            
            // Compute the mass fraction of the last species.
            Y_last -= *(mass_fractions[si]);
        }
        
        // Add the contribution from the last species.
        c_p += Y_last*d_species_c_p.back();
        c_v += Y_last*d_species_c_v.back();
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Number of mass fractions provided is not"
            << " equal to the total number of species or (total number of species - 1)."
            << std::endl);
    }
    
    gamma = c_p/c_v;
    R = c_p - c_v;
}


/*
 * Compute the thermodynamic properties of the mixture with mass fractions.
 */
void
EquationOfStateMixingRulesIdealGas::computeMixtureThermodynamicPropertiesWithMassFractions(
    HAMERS_SHARED_PTR<pdat::CellData<double> >& data_mixture_thermo_properties,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_mass_fractions,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(data_mixture_thermo_properties);
    TBOX_ASSERT(data_mixture_thermo_properties->getDepth() == 4);
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_mixture_thermo_properties = data_mixture_thermo_properties->getGhostBox();
    const hier::IntVector ghostcell_dims_mixture_thermo_properties = ghost_box_mixture_thermo_properties.numberCells();
    
    const hier::Box ghost_box_mass_fractions = data_mass_fractions->getGhostBox();
    const hier::IntVector ghostcell_dims_mass_fractions = ghost_box_mass_fractions.numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_mixture_thermo_properties(d_dim);
    hier::IntVector offset_mass_fractions(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_mixture_thermo_properties = data_mixture_thermo_properties->getGhostCellWidth();
        const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_mixture_thermo_properties->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
        TBOX_ASSERT(data_mass_fractions->getBox().isSpatiallyEqual(interior_box));
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_mixture_thermo_properties;
        num_ghosts_min = hier::IntVector::min(num_ghosts_mass_fractions, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_mixture_thermo_properties = num_ghosts_mixture_thermo_properties;
        offset_mass_fractions = num_ghosts_mass_fractions;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
        TBOX_ASSERT(data_mixture_thermo_properties->getGhostBox().contains(domain));
        TBOX_ASSERT(data_mass_fractions->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_mixture_thermo_properties = domain.lower() - ghost_box_mixture_thermo_properties.lower();
        offset_mass_fractions = domain.lower() - ghost_box_mass_fractions.lower();
    }
    
    /*
     * Get the pointers to the cell data of mixture thermodynamic properties.
     */
    
    double* const gamma = data_mixture_thermo_properties->getPointer(0);
    double* const R = data_mixture_thermo_properties->getPointer(1);
    double* const c_p = data_mixture_thermo_properties->getPointer(2);
    double* const c_v = data_mixture_thermo_properties->getPointer(3);
    
    /*
     * Fill zeros for c_p and c_v.
     */
    
    if (domain.empty())
    {
        data_mixture_thermo_properties->fill(double(0), 2);
        data_mixture_thermo_properties->fill(double(0), 3);
    }
    else
    {
        data_mixture_thermo_properties->fill(double(0), domain, 2);
        data_mixture_thermo_properties->fill(double(0), domain, 3);
    }
    
    if (data_mass_fractions->getDepth() == d_num_species)
    {
        /*
         * Get the pointers to the cell data of mass fractions.
         */
        
        std::vector<const double*> Y;
        Y.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            Y.push_back(data_mass_fractions->getPointer(si));
        }
        
        computeMixtureThermodynamicPropertiesWithMassFractions(
            gamma,
            R,
            c_p,
            c_v,
            Y,
            offset_mixture_thermo_properties,
            offset_mass_fractions,
            ghostcell_dims_mixture_thermo_properties,
            ghostcell_dims_mass_fractions,
            domain_lo,
            domain_dims);
    }
    else if (data_mass_fractions->getDepth() == d_num_species - 1)
    {
        HAMERS_SHARED_PTR<pdat::CellData<double> > data_mass_fractions_last;
        
        hier::IntVector offset_mass_fractions_last(d_dim);
        
        if (domain.empty())
        {
            const hier::Box interior_box = data_mass_fractions->getBox();
            
            offset_mass_fractions_last = data_mass_fractions->getGhostCellWidth();
            
            data_mass_fractions_last = HAMERS_MAKE_SHARED<pdat::CellData<double> >(
                interior_box, 1, data_mass_fractions->getGhostCellWidth());
            
            data_mass_fractions_last->fillAll(double(1));
        }
        else
        {
            offset_mass_fractions_last = hier::IntVector::getZero(d_dim);
            
            data_mass_fractions_last = HAMERS_MAKE_SHARED<pdat::CellData<double> >(
                domain, 1, hier::IntVector::getZero(d_dim));
            
            data_mass_fractions_last->fillAll(double(1), domain);
        }
        
        /*
         * Get the pointers to the cell data of mass fractions and the dimensions of the ghost cell box of
         * last mass fraction.
         */
        
        std::vector<const double*> Y;
        Y.reserve(d_num_species - 1);
        for (int si = 0; si < d_num_species - 1; si++)
        {
            Y.push_back(data_mass_fractions->getPointer(si));
        }
        
        double* const Y_last = data_mass_fractions_last->getPointer(0);
        
        const hier::Box ghost_box_mass_fractions_last = data_mass_fractions_last->getGhostBox();
        const hier::IntVector ghostcell_dims_mass_fractions_last = ghost_box_mass_fractions_last.numberCells();
        
        computeMixtureThermodynamicPropertiesWithMassFractions(
            gamma,
            R,
            c_p,
            c_v,
            Y_last,
            Y,
            offset_mixture_thermo_properties,
            offset_mass_fractions_last,
            offset_mass_fractions,
            ghostcell_dims_mixture_thermo_properties,
            ghostcell_dims_mass_fractions_last,
            ghostcell_dims_mass_fractions,
            domain_lo,
            domain_dims);
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
 * Compute the thermodynamic properties of the mixture with mass fractions.
 */
void
EquationOfStateMixingRulesIdealGas::computeMixtureThermodynamicPropertiesWithMassFractions(
    HAMERS_SHARED_PTR<pdat::SideData<double> >& data_mixture_thermo_properties,
    const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_mass_fractions,
    int side_normal,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(data_mixture_thermo_properties);
    TBOX_ASSERT(data_mixture_thermo_properties->getDepth() == 4);
#endif
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(side_normal < d_dim.getValue());
    
    TBOX_ASSERT(data_mixture_thermo_properties->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_mass_fractions->getDirectionVector()[side_normal] > 0);
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_mixture_thermo_properties = data_mixture_thermo_properties->getGhostBox();
    hier::IntVector ghostcell_dims_mixture_thermo_properties = ghost_box_mixture_thermo_properties.numberCells();
    
    const hier::Box ghost_box_mass_fractions = data_mass_fractions->getGhostBox();
    hier::IntVector ghostcell_dims_mass_fractions = ghost_box_mass_fractions.numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_mixture_thermo_properties(d_dim);
    hier::IntVector offset_mass_fractions(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_mixture_thermo_properties = data_mixture_thermo_properties->getGhostCellWidth();
        const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_mixture_thermo_properties->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
        TBOX_ASSERT(data_mass_fractions->getBox().isSpatiallyEqual(interior_box));
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_mixture_thermo_properties;
        num_ghosts_min = hier::IntVector::min(num_ghosts_mass_fractions, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_mixture_thermo_properties = num_ghosts_mixture_thermo_properties;
        offset_mass_fractions = num_ghosts_mass_fractions;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
        TBOX_ASSERT(data_mixture_thermo_properties->getGhostBox().contains(domain));
        TBOX_ASSERT(data_mass_fractions->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_mixture_thermo_properties = domain.lower() - ghost_box_mixture_thermo_properties.lower();
        offset_mass_fractions = domain.lower() - ghost_box_mass_fractions.lower();
    }
    
    ghostcell_dims_mixture_thermo_properties[side_normal]++;
    ghostcell_dims_mass_fractions[side_normal]++;
    domain_dims[side_normal]++;
    
    /*
     * Get the pointers to the cell data of mixture thermodynamic properties.
     */
    
    double* const gamma = data_mixture_thermo_properties->getPointer(side_normal, 0);
    double* const R = data_mixture_thermo_properties->getPointer(side_normal, 1);
    double* const c_p = data_mixture_thermo_properties->getPointer(side_normal, 2);
    double* const c_v = data_mixture_thermo_properties->getPointer(side_normal, 3);
    
    /*
     * Fill zeros for c_p and c_v.
     */
    
    if (domain.empty())
    {
        data_mixture_thermo_properties->fill(double(0), 2);
        data_mixture_thermo_properties->fill(double(0), 3);
    }
    else
    {
        data_mixture_thermo_properties->fill(double(0), domain, 2);
        data_mixture_thermo_properties->fill(double(0), domain, 3);
    }
    
    if (data_mass_fractions->getDepth() == d_num_species)
    {
        /*
         * Get the pointers to the cell data of mass fractions.
         */
        
        std::vector<const double*> Y;
        Y.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            Y.push_back(data_mass_fractions->getPointer(side_normal, si));
        }
        
        computeMixtureThermodynamicPropertiesWithMassFractions(
            gamma,
            R,
            c_p,
            c_v,
            Y,
            offset_mixture_thermo_properties,
            offset_mass_fractions,
            ghostcell_dims_mixture_thermo_properties,
            ghostcell_dims_mass_fractions,
            domain_lo,
            domain_dims);
    }
    else if (data_mass_fractions->getDepth() == d_num_species - 1)
    {
        hier::IntVector direction = hier::IntVector::getZero(d_dim);
        direction[side_normal] = 1;
        
        HAMERS_SHARED_PTR<pdat::SideData<double> > data_mass_fractions_last;
        
        hier::IntVector offset_mass_fractions_last(d_dim);
        
        if (domain.empty())
        {
            const hier::Box interior_box = data_mass_fractions->getBox();
            
            offset_mass_fractions_last = data_mass_fractions->getGhostCellWidth();
            
            data_mass_fractions_last = HAMERS_MAKE_SHARED<pdat::SideData<double> >(
                interior_box, 1, data_mass_fractions->getGhostCellWidth(), direction);
            
            data_mass_fractions_last->fillAll(double(1));
        }
        else
        {
            offset_mass_fractions_last = hier::IntVector::getZero(d_dim);
            
            data_mass_fractions_last = HAMERS_MAKE_SHARED<pdat::SideData<double> >(
                domain, 1, hier::IntVector::getZero(d_dim), direction);
            
            data_mass_fractions_last->fillAll(double(1), domain);
        }
        
        /*
         * Get the pointers to the cell data of mass fractions and the dimensions of the ghost cell box of
         * last mass fraction.
         */
        
        std::vector<const double*> Y;
        Y.reserve(d_num_species - 1);
        for (int si = 0; si < d_num_species - 1; si++)
        {
            Y.push_back(data_mass_fractions->getPointer(side_normal, si));
        }
        
        double* const Y_last = data_mass_fractions_last->getPointer(side_normal, 0);
        
        const hier::Box ghost_box_mass_fractions_last = data_mass_fractions_last->getGhostBox();
        hier::IntVector ghostcell_dims_mass_fractions_last = ghost_box_mass_fractions_last.numberCells();
        ghostcell_dims_mass_fractions_last[side_normal]++;
               
        computeMixtureThermodynamicPropertiesWithMassFractions(
            gamma,
            R,
            c_p,
            c_v,
            Y_last,
            Y,
            offset_mixture_thermo_properties,
            offset_mass_fractions_last,
            offset_mass_fractions,
            ghostcell_dims_mixture_thermo_properties,
            ghostcell_dims_mass_fractions_last,
            ghostcell_dims_mass_fractions,
            domain_lo,
            domain_dims);
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
 * Compute the thermodynamic properties of the mixture with volume fractions.
 */
void
EquationOfStateMixingRulesIdealGas::getMixtureThermodynamicPropertiesWithVolumeFractions(
    std::vector<double*>& mixture_thermo_properties,
    const std::vector<const double*>& volume_fractions) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(mixture_thermo_properties.size()) == 1);
#endif
    
    // Get reference to gamma of mixture.
    double& gamma = *(mixture_thermo_properties[0]);
    double xi = double(0);
    
    if (static_cast<int>(volume_fractions.size()) == d_num_species)
    {
        for (int si = 0; si < d_num_species; si++)
        {
            xi += *(volume_fractions[si])/(d_species_gamma[si] - double(1));
        }
    }
    else if (static_cast<int>(volume_fractions.size()) == d_num_species - 1)
    {
        double Z_last = double(1);
        
        for (int si = 0; si < d_num_species - 1; si++)
        {
            xi += *(volume_fractions[si])/(d_species_gamma[si] - double(1));
            
            // Compute the volume fraction of the last species.
            Z_last -= *(volume_fractions[si]);
        }
        
        // Add the contribution from the last species.
        xi += Z_last/(d_species_gamma.back() - double(1));
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Number of volume fractions provided is not"
            << " equal to the total number of species or (total number of species - 1)."
            << std::endl);
    }
    
    gamma = double(1)/xi + double(1);
}


/*
 * Compute the thermodynamic properties of the mixture with volume fractions.
 */
void
EquationOfStateMixingRulesIdealGas::computeMixtureThermodynamicPropertiesWithVolumeFractions(
    HAMERS_SHARED_PTR<pdat::CellData<double> >& data_mixture_thermo_properties,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_volume_fractions,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(data_mixture_thermo_properties);
    TBOX_ASSERT(data_mixture_thermo_properties->getDepth() == 1);
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_mixture_thermo_properties = data_mixture_thermo_properties->getGhostBox();
    const hier::IntVector ghostcell_dims_mixture_thermo_properties = ghost_box_mixture_thermo_properties.numberCells();
    
    const hier::Box ghost_box_volume_fractions = data_volume_fractions->getGhostBox();
    const hier::IntVector ghostcell_dims_volume_fractions = ghost_box_volume_fractions.numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_mixture_thermo_properties(d_dim);
    hier::IntVector offset_volume_fractions(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_mixture_thermo_properties = data_mixture_thermo_properties->getGhostCellWidth();
        const hier::IntVector num_ghosts_volume_fractions = data_volume_fractions->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_mixture_thermo_properties->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
        TBOX_ASSERT(data_volume_fractions->getBox().isSpatiallyEqual(interior_box));
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_mixture_thermo_properties;
        num_ghosts_min = hier::IntVector::min(num_ghosts_volume_fractions, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_mixture_thermo_properties = num_ghosts_mixture_thermo_properties;
        offset_volume_fractions = num_ghosts_volume_fractions;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
        TBOX_ASSERT(data_mixture_thermo_properties->getGhostBox().contains(domain));
        TBOX_ASSERT(data_volume_fractions->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_mixture_thermo_properties = domain.lower() - ghost_box_mixture_thermo_properties.lower();
        offset_volume_fractions = domain.lower() - ghost_box_volume_fractions.lower();
    }
    
    /*
     * Get the pointer to the cell data of mixture thermodynamic property.
     */
    
    double* const gamma = data_mixture_thermo_properties->getPointer(0);
    
    /*
     * Fill zeros for gamma.
     */
    
    if (domain.empty())
    {
        data_mixture_thermo_properties->fill(double(0), 0);
    }
    else
    {
        data_mixture_thermo_properties->fill(double(0), domain, 0);
    }
    
    if (data_volume_fractions->getDepth() == d_num_species)
    {
        /*
         * Get the pointers to the cell data of volume fractions.
         */
        
        std::vector<const double*> Z;
        Z.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            Z.push_back(data_volume_fractions->getPointer(si));
        }
        
        computeMixtureThermodynamicPropertiesWithVolumeFractions(
            gamma,
            Z,
            offset_mixture_thermo_properties,
            offset_volume_fractions,
            ghostcell_dims_mixture_thermo_properties,
            ghostcell_dims_volume_fractions,
            domain_lo,
            domain_dims);
    }
    else if (data_volume_fractions->getDepth() == d_num_species - 1)
    {
        HAMERS_SHARED_PTR<pdat::CellData<double> > data_volume_fractions_last;
        
        hier::IntVector offset_volume_fractions_last(d_dim);
        
        if (domain.empty())
        {
            const hier::Box interior_box = data_volume_fractions->getBox();
            
            offset_volume_fractions_last = data_volume_fractions->getGhostCellWidth();
            
            data_volume_fractions_last = HAMERS_MAKE_SHARED<pdat::CellData<double> >(
                interior_box, 1, data_volume_fractions->getGhostCellWidth());
            
            data_volume_fractions_last->fillAll(double(1));
        }
        else
        {
            offset_volume_fractions_last = hier::IntVector::getZero(d_dim);
            
            data_volume_fractions_last = HAMERS_MAKE_SHARED<pdat::CellData<double> >(
                domain, 1, hier::IntVector::getZero(d_dim));
            
            data_volume_fractions_last->fillAll(double(1), domain);
        }
        
        /*
         * Get the pointers to the cell data of volume fractions and the dimensions of the ghost cell box of
         * last volume fraction.
         */
        
        std::vector<const double*> Z;
        Z.reserve(d_num_species - 1);
        for (int si = 0; si < d_num_species - 1; si++)
        {
            Z.push_back(data_volume_fractions->getPointer(si));
        }
        
        double* const Z_last = data_volume_fractions_last->getPointer(0);
        
        const hier::Box ghost_box_volume_fractions_last = data_volume_fractions_last->getGhostBox();
        const hier::IntVector ghostcell_dims_volume_fractions_last = ghost_box_volume_fractions_last.numberCells();
        
        computeMixtureThermodynamicPropertiesWithVolumeFractions(
            gamma,
            Z_last,
            Z,
            offset_mixture_thermo_properties,
            offset_volume_fractions_last,
            offset_volume_fractions,
            ghostcell_dims_mixture_thermo_properties,
            ghostcell_dims_volume_fractions_last,
            ghostcell_dims_volume_fractions,
            domain_lo,
            domain_dims);
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
 * Compute the thermodynamic properties of the mixture with volume fractions.
 */
void
EquationOfStateMixingRulesIdealGas::computeMixtureThermodynamicPropertiesWithVolumeFractions(
    HAMERS_SHARED_PTR<pdat::SideData<double> >& data_mixture_thermo_properties,
    const HAMERS_SHARED_PTR<pdat::SideData<double> >& data_volume_fractions,
    int side_normal,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(data_mixture_thermo_properties);
    TBOX_ASSERT(data_mixture_thermo_properties->getDepth() == 1);
#endif
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(side_normal < d_dim.getValue());
    
    TBOX_ASSERT(data_mixture_thermo_properties->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_volume_fractions->getDirectionVector()[side_normal] > 0);
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_mixture_thermo_properties = data_mixture_thermo_properties->getGhostBox();
    hier::IntVector ghostcell_dims_mixture_thermo_properties = ghost_box_mixture_thermo_properties.numberCells();
    
    const hier::Box ghost_box_volume_fractions = data_volume_fractions->getGhostBox();
    hier::IntVector ghostcell_dims_volume_fractions = ghost_box_volume_fractions.numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_mixture_thermo_properties(d_dim);
    hier::IntVector offset_volume_fractions(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_mixture_thermo_properties = data_mixture_thermo_properties->getGhostCellWidth();
        const hier::IntVector num_ghosts_volume_fractions = data_volume_fractions->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_mixture_thermo_properties->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
        TBOX_ASSERT(data_volume_fractions->getBox().isSpatiallyEqual(interior_box));
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_mixture_thermo_properties;
        num_ghosts_min = hier::IntVector::min(num_ghosts_volume_fractions, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_mixture_thermo_properties = num_ghosts_mixture_thermo_properties;
        offset_volume_fractions = num_ghosts_volume_fractions;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
        TBOX_ASSERT(data_mixture_thermo_properties->getGhostBox().contains(domain));
        TBOX_ASSERT(data_volume_fractions->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_mixture_thermo_properties = domain.lower() - ghost_box_mixture_thermo_properties.lower();
        offset_volume_fractions = domain.lower() - ghost_box_volume_fractions.lower();
    }
    
    ghostcell_dims_mixture_thermo_properties[side_normal]++;
    ghostcell_dims_volume_fractions[side_normal]++;
    domain_dims[side_normal]++;
    
    /*
     * Get the pointer to the cell data of mixture thermodynamic property.
     */
    
    double* const gamma = data_mixture_thermo_properties->getPointer(side_normal, 0);
    
    /*
     * Fill zeros for gamma.
     */
    
    if (domain.empty())
    {
        data_mixture_thermo_properties->fill(double(0), 0);
    }
    else
    {
        data_mixture_thermo_properties->fill(double(0), domain, 0);
    }
    
    if (data_volume_fractions->getDepth() == d_num_species)
    {
        /*
         * Get the pointers to the cell data of volume fractions.
         */
        
        std::vector<const double*> Z;
        Z.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            Z.push_back(data_volume_fractions->getPointer(side_normal, si));
        }
        
        computeMixtureThermodynamicPropertiesWithVolumeFractions(
            gamma,
            Z,
            offset_mixture_thermo_properties,
            offset_volume_fractions,
            ghostcell_dims_mixture_thermo_properties,
            ghostcell_dims_volume_fractions,
            domain_lo,
            domain_dims);
    }
    else if (data_volume_fractions->getDepth() == d_num_species - 1)
    {
        hier::IntVector direction = hier::IntVector::getZero(d_dim);
        direction[side_normal] = 1;
        
        HAMERS_SHARED_PTR<pdat::SideData<double> > data_volume_fractions_last;
        
        hier::IntVector offset_volume_fractions_last(d_dim);
        
        if (domain.empty())
        {
            const hier::Box interior_box = data_volume_fractions->getBox();
            
            offset_volume_fractions_last = data_volume_fractions->getGhostCellWidth();
            
            data_volume_fractions_last = HAMERS_MAKE_SHARED<pdat::SideData<double> >(
                interior_box, 1, data_volume_fractions->getGhostCellWidth(), direction);
            
            data_volume_fractions_last->fillAll(double(1));
        }
        else
        {
            offset_volume_fractions_last = hier::IntVector::getZero(d_dim);
            
            data_volume_fractions_last = HAMERS_MAKE_SHARED<pdat::SideData<double> >(
                domain, 1, hier::IntVector::getZero(d_dim), direction);
            
            data_volume_fractions_last->fillAll(double(1), domain);
        }
        
        /*
         * Get the pointers to the cell data of volume fractions and the dimensions of the ghost cell box of
         * last volume fraction.
         */
        
        std::vector<const double*> Z;
        Z.reserve(d_num_species - 1);
        for (int si = 0; si < d_num_species - 1; si++)
        {
            Z.push_back(data_volume_fractions->getPointer(side_normal, si));
        }
        
        double* const Z_last = data_volume_fractions_last->getPointer(side_normal, 0);
        
        const hier::Box ghost_box_volume_fractions_last = data_volume_fractions_last->getGhostBox();
        hier::IntVector ghostcell_dims_volume_fractions_last = ghost_box_volume_fractions_last.numberCells();
        ghostcell_dims_volume_fractions_last[side_normal]++;
        
        computeMixtureThermodynamicPropertiesWithVolumeFractions(
            gamma,
            Z_last,
            Z,
            offset_mixture_thermo_properties,
            offset_volume_fractions_last,
            offset_volume_fractions,
            ghostcell_dims_mixture_thermo_properties,
            ghostcell_dims_volume_fractions_last,
            ghostcell_dims_volume_fractions,
            domain_lo,
            domain_dims);
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
 * Compute the isochoric specific heat capacity of mixture with isothermal and isobaric
 * equilibrium assumptions.
 */
void
EquationOfStateMixingRulesIdealGas::computeIsochoricSpecificHeatCapacity(
    double* const c_v,
    const std::vector<const double*>& Y,
    const hier::IntVector& offset_isochoric_specific_heat_capacity,
    const hier::IntVector& offset_mass_fractions,
    const hier::IntVector& ghostcell_dims_isochoric_specific_heat_capacity,
    const hier::IntVector& ghostcell_dims_mass_fractions,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims) const
{
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and offsets.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int offset_0_isochoric_specific_heat_capacity =
            offset_isochoric_specific_heat_capacity[0];
        const int offset_0_mass_fractions = offset_mass_fractions[0];
        
        // Compute c_v.
        for (int si = 0; si < d_num_species; si++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_isochoric_specific_heat_capacity = i + offset_0_isochoric_specific_heat_capacity;
                const int idx_mass_fractions = i + offset_0_mass_fractions;
                
                c_v[idx_isochoric_specific_heat_capacity] += Y[si][idx_mass_fractions]*d_species_c_v[si];
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
        
        const int offset_0_isochoric_specific_heat_capacity = offset_isochoric_specific_heat_capacity[0];
        const int offset_1_isochoric_specific_heat_capacity = offset_isochoric_specific_heat_capacity[1];
        const int ghostcell_dim_0_isochoric_specific_heat_capacity =
            ghostcell_dims_isochoric_specific_heat_capacity[0];
        
        const int offset_0_mass_fractions = offset_mass_fractions[0];
        const int offset_1_mass_fractions = offset_mass_fractions[1];
        const int ghostcell_dim_0_mass_fractions = ghostcell_dims_mass_fractions[0];
        
        // Compute c_v.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_isochoric_specific_heat_capacity =
                        (i + offset_0_isochoric_specific_heat_capacity) +
                        (j + offset_1_isochoric_specific_heat_capacity)*
                            ghostcell_dim_0_isochoric_specific_heat_capacity;
                    
                    const int idx_mass_fractions = (i + offset_0_mass_fractions) +
                        (j + offset_1_mass_fractions)*ghostcell_dim_0_mass_fractions;
                    
                    c_v[idx_isochoric_specific_heat_capacity] += Y[si][idx_mass_fractions]*d_species_c_v[si];
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
        
        const int offset_0_isochoric_specific_heat_capacity = offset_isochoric_specific_heat_capacity[0];
        const int offset_1_isochoric_specific_heat_capacity = offset_isochoric_specific_heat_capacity[1];
        const int offset_2_isochoric_specific_heat_capacity = offset_isochoric_specific_heat_capacity[2];
        const int ghostcell_dim_0_isochoric_specific_heat_capacity =
            ghostcell_dims_isochoric_specific_heat_capacity[0];
        const int ghostcell_dim_1_isochoric_specific_heat_capacity =
            ghostcell_dims_isochoric_specific_heat_capacity[1];
        
        const int offset_0_mass_fractions = offset_mass_fractions[0];
        const int offset_1_mass_fractions = offset_mass_fractions[1];
        const int offset_2_mass_fractions = offset_mass_fractions[2];
        const int ghostcell_dim_0_mass_fractions = ghostcell_dims_mass_fractions[0];
        const int ghostcell_dim_1_mass_fractions = ghostcell_dims_mass_fractions[1];
        
        // Compute c_v.
        for (int si = 0; si < d_num_species; si++)
        {
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
                        const int idx_isochoric_specific_heat_capacity =
                            (i + offset_0_isochoric_specific_heat_capacity) +
                            (j + offset_1_isochoric_specific_heat_capacity)*
                                ghostcell_dim_0_isochoric_specific_heat_capacity +
                            (k + offset_2_isochoric_specific_heat_capacity)*
                                ghostcell_dim_0_isochoric_specific_heat_capacity*
                                    ghostcell_dim_1_isochoric_specific_heat_capacity;
                        
                        const int idx_mass_fractions = (i + offset_0_mass_fractions) +
                            (j + offset_1_mass_fractions)*ghostcell_dim_0_mass_fractions +
                            (k + offset_2_mass_fractions)*ghostcell_dim_0_mass_fractions*
                                ghostcell_dim_1_mass_fractions;
                        
                        c_v[idx_isochoric_specific_heat_capacity] += Y[si][idx_mass_fractions]*d_species_c_v[si];
                    }
                }
            }
        }
    }
}


/*
 * Compute the isochoric specific heat capacity of mixture with isothermal and isobaric
 * equilibrium assumptions.
 */
void
EquationOfStateMixingRulesIdealGas::computeIsochoricSpecificHeatCapacity(
    double* const c_v,
    double* const Y_last,
    const std::vector<const double*>& Y,
    const hier::IntVector& offset_isochoric_specific_heat_capacity,
    const hier::IntVector& offset_mass_fractions_last,
    const hier::IntVector& offset_mass_fractions,
    const hier::IntVector& ghostcell_dims_isochoric_specific_heat_capacity,
    const hier::IntVector& ghostcell_dims_mass_fractions_last,
    const hier::IntVector& ghostcell_dims_mass_fractions,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims) const
{
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and offsets.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int offset_0_isochoric_specific_heat_capacity =
            offset_isochoric_specific_heat_capacity[0];
        const int offset_0_mass_fractions_last = offset_mass_fractions_last[0];
        const int offset_0_mass_fractions = offset_mass_fractions[0];
        
        // Compute c_v.
        for (int si = 0; si < d_num_species - 1; si++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_isochoric_specific_heat_capacity = i + offset_0_isochoric_specific_heat_capacity;
                const int idx_mass_fractions_last = i + offset_0_mass_fractions_last;
                const int idx_mass_fractions = i + offset_0_mass_fractions;
                
                c_v[idx_isochoric_specific_heat_capacity] += Y[si][idx_mass_fractions]*d_species_c_v[si];
                
                // Compute the mass fraction of the last species.
                Y_last[idx_mass_fractions_last] -= Y[si][idx_mass_fractions];
            }
        }
        
        // Add the contribution from the last species.
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_isochoric_specific_heat_capacity = i + offset_0_isochoric_specific_heat_capacity;
            const int idx_mass_fractions_last = i + offset_0_mass_fractions_last;
            
            c_v[idx_isochoric_specific_heat_capacity] += Y_last[idx_mass_fractions_last]*d_species_c_v.back();
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
        
        const int offset_0_isochoric_specific_heat_capacity = offset_isochoric_specific_heat_capacity[0];
        const int offset_1_isochoric_specific_heat_capacity = offset_isochoric_specific_heat_capacity[1];
        const int ghostcell_dim_0_isochoric_specific_heat_capacity =
            ghostcell_dims_isochoric_specific_heat_capacity[0];
        
        const int offset_0_mass_fractions_last = offset_mass_fractions_last[0];
        const int offset_1_mass_fractions_last = offset_mass_fractions_last[1];
        const int ghostcell_dim_0_mass_fractions_last = ghostcell_dims_mass_fractions_last[0];
        
        const int offset_0_mass_fractions = offset_mass_fractions[0];
        const int offset_1_mass_fractions = offset_mass_fractions[1];
        const int ghostcell_dim_0_mass_fractions = ghostcell_dims_mass_fractions[0];
        
        // Compute c_v.
        for (int si = 0; si < d_num_species - 1; si++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_isochoric_specific_heat_capacity =
                        (i + offset_0_isochoric_specific_heat_capacity) +
                        (j + offset_1_isochoric_specific_heat_capacity)*
                            ghostcell_dim_0_isochoric_specific_heat_capacity;
                    
                    const int idx_mass_fractions_last = (i + offset_0_mass_fractions_last) +
                        (j + offset_1_mass_fractions_last)*ghostcell_dim_0_mass_fractions_last;
                    
                    const int idx_mass_fractions = (i + offset_0_mass_fractions) +
                        (j + offset_1_mass_fractions)*ghostcell_dim_0_mass_fractions;
                    
                    c_v[idx_isochoric_specific_heat_capacity] += Y[si][idx_mass_fractions]*d_species_c_v[si];
                    
                    // Compute the mass fraction of the last species.
                    Y_last[idx_mass_fractions_last] -= Y[si][idx_mass_fractions];
                }
            }
        }
        
        // Add the contribution from the last species.
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_isochoric_specific_heat_capacity =
                    (i + offset_0_isochoric_specific_heat_capacity) +
                    (j + offset_1_isochoric_specific_heat_capacity)*
                        ghostcell_dim_0_isochoric_specific_heat_capacity;
                
                const int idx_mass_fractions_last = (i + offset_0_mass_fractions_last) +
                    (j + offset_1_mass_fractions_last)*ghostcell_dim_0_mass_fractions_last;
                
                c_v[idx_isochoric_specific_heat_capacity] += Y_last[idx_mass_fractions_last]*d_species_c_v.back();
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
        
        const int offset_0_isochoric_specific_heat_capacity = offset_isochoric_specific_heat_capacity[0];
        const int offset_1_isochoric_specific_heat_capacity = offset_isochoric_specific_heat_capacity[1];
        const int offset_2_isochoric_specific_heat_capacity = offset_isochoric_specific_heat_capacity[2];
        const int ghostcell_dim_0_isochoric_specific_heat_capacity =
            ghostcell_dims_isochoric_specific_heat_capacity[0];
        const int ghostcell_dim_1_isochoric_specific_heat_capacity =
            ghostcell_dims_isochoric_specific_heat_capacity[1];
        
        const int offset_0_mass_fractions_last = offset_mass_fractions_last[0];
        const int offset_1_mass_fractions_last = offset_mass_fractions_last[1];
        const int offset_2_mass_fractions_last = offset_mass_fractions_last[2];
        const int ghostcell_dim_0_mass_fractions_last = ghostcell_dims_mass_fractions_last[0];
        const int ghostcell_dim_1_mass_fractions_last = ghostcell_dims_mass_fractions_last[1];
        
        const int offset_0_mass_fractions = offset_mass_fractions[0];
        const int offset_1_mass_fractions = offset_mass_fractions[1];
        const int offset_2_mass_fractions = offset_mass_fractions[2];
        const int ghostcell_dim_0_mass_fractions = ghostcell_dims_mass_fractions[0];
        const int ghostcell_dim_1_mass_fractions = ghostcell_dims_mass_fractions[1];
        
        // Compute c_v.
        for (int si = 0; si < d_num_species - 1; si++)
        {
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
                        const int idx_isochoric_specific_heat_capacity =
                            (i + offset_0_isochoric_specific_heat_capacity) +
                            (j + offset_1_isochoric_specific_heat_capacity)*
                                ghostcell_dim_0_isochoric_specific_heat_capacity +
                            (k + offset_2_isochoric_specific_heat_capacity)*
                                ghostcell_dim_0_isochoric_specific_heat_capacity*
                                    ghostcell_dim_1_isochoric_specific_heat_capacity;
                        
                        const int idx_mass_fractions_last = (i + offset_0_mass_fractions_last) +
                            (j + offset_1_mass_fractions_last)*ghostcell_dim_0_mass_fractions_last +
                            (k + offset_2_mass_fractions_last)*ghostcell_dim_0_mass_fractions_last*
                                ghostcell_dim_1_mass_fractions_last;
                        
                        const int idx_mass_fractions = (i + offset_0_mass_fractions) +
                            (j + offset_1_mass_fractions)*ghostcell_dim_0_mass_fractions +
                            (k + offset_2_mass_fractions)*ghostcell_dim_0_mass_fractions*
                                ghostcell_dim_1_mass_fractions;
                        
                        c_v[idx_isochoric_specific_heat_capacity] += Y[si][idx_mass_fractions]*d_species_c_v[si];
                        
                        // Compute the mass fraction of the last species.
                        Y_last[idx_mass_fractions_last] -= Y[si][idx_mass_fractions];
                    }
                }
            }
        }
        
        // Add the contribution from the last species.
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
                    const int idx_isochoric_specific_heat_capacity =
                        (i + offset_0_isochoric_specific_heat_capacity) +
                        (j + offset_1_isochoric_specific_heat_capacity)*
                            ghostcell_dim_0_isochoric_specific_heat_capacity +
                        (k + offset_2_isochoric_specific_heat_capacity)*
                            ghostcell_dim_0_isochoric_specific_heat_capacity*
                                ghostcell_dim_1_isochoric_specific_heat_capacity;
                    
                    const int idx_mass_fractions_last = (i + offset_0_mass_fractions_last) +
                        (j + offset_1_mass_fractions_last)*ghostcell_dim_0_mass_fractions_last +
                        (k + offset_2_mass_fractions_last)*ghostcell_dim_0_mass_fractions_last*
                            ghostcell_dim_1_mass_fractions_last;
                    
                    c_v[idx_isochoric_specific_heat_capacity] += Y_last[idx_mass_fractions_last]*d_species_c_v.back();
                }
            }
        }
    }
}


/*
 * Compute the isobaric specific heat capacity of mixture with isothermal and isobaric equilibrium
 * assumptions.
 */
void
EquationOfStateMixingRulesIdealGas::computeIsobaricSpecificHeatCapacity(
    double* const c_p,
    const std::vector<const double*>& Y,
    const hier::IntVector& offset_isobaric_specific_heat_capacity,
    const hier::IntVector& offset_mass_fractions,
    const hier::IntVector& ghostcell_dims_isobaric_specific_heat_capacity,
    const hier::IntVector& ghostcell_dims_mass_fractions,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims) const
{
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and offsets.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int offset_0_isobaric_specific_heat_capacity =
            offset_isobaric_specific_heat_capacity[0];
        const int offset_0_mass_fractions = offset_mass_fractions[0];
        
        // Compute c_p.
        for (int si = 0; si < d_num_species; si++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_isobaric_specific_heat_capacity = i + offset_0_isobaric_specific_heat_capacity;
                const int idx_mass_fractions = i + offset_0_mass_fractions;
                
                c_p[idx_isobaric_specific_heat_capacity] += Y[si][idx_mass_fractions]*d_species_c_p[si];
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
        
        const int offset_0_isobaric_specific_heat_capacity = offset_isobaric_specific_heat_capacity[0];
        const int offset_1_isobaric_specific_heat_capacity = offset_isobaric_specific_heat_capacity[1];
        const int ghostcell_dim_0_isobaric_specific_heat_capacity =
            ghostcell_dims_isobaric_specific_heat_capacity[0];
        
        const int offset_0_mass_fractions = offset_mass_fractions[0];
        const int offset_1_mass_fractions = offset_mass_fractions[1];
        const int ghostcell_dim_0_mass_fractions = ghostcell_dims_mass_fractions[0];
        
        // Compute c_p.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_isobaric_specific_heat_capacity =
                        (i + offset_0_isobaric_specific_heat_capacity) +
                        (j + offset_1_isobaric_specific_heat_capacity)*
                            ghostcell_dim_0_isobaric_specific_heat_capacity;
                    
                    const int idx_mass_fractions = (i + offset_0_mass_fractions) +
                        (j + offset_1_mass_fractions)*ghostcell_dim_0_mass_fractions;
                    
                    c_p[idx_isobaric_specific_heat_capacity] += Y[si][idx_mass_fractions]*d_species_c_p[si];
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
        
        const int offset_0_isobaric_specific_heat_capacity = offset_isobaric_specific_heat_capacity[0];
        const int offset_1_isobaric_specific_heat_capacity = offset_isobaric_specific_heat_capacity[1];
        const int offset_2_isobaric_specific_heat_capacity = offset_isobaric_specific_heat_capacity[2];
        const int ghostcell_dim_0_isobaric_specific_heat_capacity =
            ghostcell_dims_isobaric_specific_heat_capacity[0];
        const int ghostcell_dim_1_isobaric_specific_heat_capacity =
            ghostcell_dims_isobaric_specific_heat_capacity[1];
        
        const int offset_0_mass_fractions = offset_mass_fractions[0];
        const int offset_1_mass_fractions = offset_mass_fractions[1];
        const int offset_2_mass_fractions = offset_mass_fractions[2];
        const int ghostcell_dim_0_mass_fractions = ghostcell_dims_mass_fractions[0];
        const int ghostcell_dim_1_mass_fractions = ghostcell_dims_mass_fractions[1];
        
        // Compute c_p.
        for (int si = 0; si < d_num_species; si++)
        {
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
                        const int idx_isobaric_specific_heat_capacity =
                            (i + offset_0_isobaric_specific_heat_capacity) +
                            (j + offset_1_isobaric_specific_heat_capacity)*
                                ghostcell_dim_0_isobaric_specific_heat_capacity +
                            (k + offset_2_isobaric_specific_heat_capacity)*
                                ghostcell_dim_0_isobaric_specific_heat_capacity*
                                    ghostcell_dim_1_isobaric_specific_heat_capacity;
                        
                        const int idx_mass_fractions = (i + offset_0_mass_fractions) +
                            (j + offset_1_mass_fractions)*ghostcell_dim_0_mass_fractions +
                            (k + offset_2_mass_fractions)*ghostcell_dim_0_mass_fractions*
                                ghostcell_dim_1_mass_fractions;
                        
                        c_p[idx_isobaric_specific_heat_capacity] += Y[si][idx_mass_fractions]*d_species_c_p[si];
                    }
                }
            }
        }
    }
}


/*
 * Compute the isobaric specific heat capacity of mixture with isothermal and isobaric equilibrium
 * assumptions.
 */
void
EquationOfStateMixingRulesIdealGas::computeIsobaricSpecificHeatCapacity(
    double* const c_p,
    double* const Y_last,
    const std::vector<const double*>& Y,
    const hier::IntVector& offset_isobaric_specific_heat_capacity,
    const hier::IntVector& offset_mass_fractions_last,
    const hier::IntVector& offset_mass_fractions,
    const hier::IntVector& ghostcell_dims_isobaric_specific_heat_capacity,
    const hier::IntVector& ghostcell_dims_mass_fractions_last,
    const hier::IntVector& ghostcell_dims_mass_fractions,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims) const
{
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and offsets.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int offset_0_isobaric_specific_heat_capacity =
            offset_isobaric_specific_heat_capacity[0];
        const int offset_0_mass_fractions_last = offset_mass_fractions_last[0];
        const int offset_0_mass_fractions = offset_mass_fractions[0];
        
        // Compute c_p.
        for (int si = 0; si < d_num_species - 1; si++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_isobaric_specific_heat_capacity = i + offset_0_isobaric_specific_heat_capacity;
                const int idx_mass_fractions_last = i + offset_0_mass_fractions_last;
                const int idx_mass_fractions = i + offset_0_mass_fractions;
                
                c_p[idx_isobaric_specific_heat_capacity] += Y[si][idx_mass_fractions]*d_species_c_p[si];
                
                // Compute the mass fraction of the last species.
                Y_last[idx_mass_fractions_last] -= Y[si][idx_mass_fractions];
            }
        }
        
        // Add the contribution from the last species.
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_isobaric_specific_heat_capacity = i + offset_0_isobaric_specific_heat_capacity;
            const int idx_mass_fractions_last = i + offset_0_mass_fractions_last;
            
            c_p[idx_isobaric_specific_heat_capacity] += Y_last[idx_mass_fractions_last]*d_species_c_p.back();
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
        
        const int offset_0_isobaric_specific_heat_capacity = offset_isobaric_specific_heat_capacity[0];
        const int offset_1_isobaric_specific_heat_capacity = offset_isobaric_specific_heat_capacity[1];
        const int ghostcell_dim_0_isobaric_specific_heat_capacity =
            ghostcell_dims_isobaric_specific_heat_capacity[0];
        
        const int offset_0_mass_fractions_last = offset_mass_fractions_last[0];
        const int offset_1_mass_fractions_last = offset_mass_fractions_last[1];
        const int ghostcell_dim_0_mass_fractions_last = ghostcell_dims_mass_fractions_last[0];
        
        const int offset_0_mass_fractions = offset_mass_fractions[0];
        const int offset_1_mass_fractions = offset_mass_fractions[1];
        const int ghostcell_dim_0_mass_fractions = ghostcell_dims_mass_fractions[0];
        
        // Compute c_p.
        for (int si = 0; si < d_num_species - 1; si++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_isobaric_specific_heat_capacity =
                        (i + offset_0_isobaric_specific_heat_capacity) +
                        (j + offset_1_isobaric_specific_heat_capacity)*
                            ghostcell_dim_0_isobaric_specific_heat_capacity;
                    
                    const int idx_mass_fractions_last = (i + offset_0_mass_fractions_last) +
                        (j + offset_1_mass_fractions_last)*ghostcell_dim_0_mass_fractions_last;
                    
                    const int idx_mass_fractions = (i + offset_0_mass_fractions) +
                        (j + offset_1_mass_fractions)*ghostcell_dim_0_mass_fractions;
                    
                    c_p[idx_isobaric_specific_heat_capacity] += Y[si][idx_mass_fractions]*d_species_c_p[si];
                    
                    // Compute the mass fraction of the last species.
                    Y_last[idx_mass_fractions_last] -= Y[si][idx_mass_fractions];
                }
            }
        }
        
        // Add the contribution from the last species.
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_isobaric_specific_heat_capacity =
                    (i + offset_0_isobaric_specific_heat_capacity) +
                    (j + offset_1_isobaric_specific_heat_capacity)*
                        ghostcell_dim_0_isobaric_specific_heat_capacity;
                
                const int idx_mass_fractions_last = (i + offset_0_mass_fractions_last) +
                    (j + offset_1_mass_fractions_last)*ghostcell_dim_0_mass_fractions_last;
                
                c_p[idx_isobaric_specific_heat_capacity] += Y_last[idx_mass_fractions_last]*d_species_c_p.back();
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
        
        const int offset_0_isobaric_specific_heat_capacity = offset_isobaric_specific_heat_capacity[0];
        const int offset_1_isobaric_specific_heat_capacity = offset_isobaric_specific_heat_capacity[1];
        const int offset_2_isobaric_specific_heat_capacity = offset_isobaric_specific_heat_capacity[2];
        const int ghostcell_dim_0_isobaric_specific_heat_capacity =
            ghostcell_dims_isobaric_specific_heat_capacity[0];
        const int ghostcell_dim_1_isobaric_specific_heat_capacity =
            ghostcell_dims_isobaric_specific_heat_capacity[1];
        
        const int offset_0_mass_fractions_last = offset_mass_fractions_last[0];
        const int offset_1_mass_fractions_last = offset_mass_fractions_last[1];
        const int offset_2_mass_fractions_last = offset_mass_fractions_last[2];
        const int ghostcell_dim_0_mass_fractions_last = ghostcell_dims_mass_fractions_last[0];
        const int ghostcell_dim_1_mass_fractions_last = ghostcell_dims_mass_fractions_last[1];
        
        const int offset_0_mass_fractions = offset_mass_fractions[0];
        const int offset_1_mass_fractions = offset_mass_fractions[1];
        const int offset_2_mass_fractions = offset_mass_fractions[2];
        const int ghostcell_dim_0_mass_fractions = ghostcell_dims_mass_fractions[0];
        const int ghostcell_dim_1_mass_fractions = ghostcell_dims_mass_fractions[1];
        
        // Compute c_p.
        for (int si = 0; si < d_num_species - 1; si++)
        {
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
                        const int idx_isobaric_specific_heat_capacity =
                            (i + offset_0_isobaric_specific_heat_capacity) +
                            (j + offset_1_isobaric_specific_heat_capacity)*
                                ghostcell_dim_0_isobaric_specific_heat_capacity +
                            (k + offset_2_isobaric_specific_heat_capacity)*
                                ghostcell_dim_0_isobaric_specific_heat_capacity*
                                    ghostcell_dim_1_isobaric_specific_heat_capacity;
                        
                        const int idx_mass_fractions_last = (i + offset_0_mass_fractions_last) +
                            (j + offset_1_mass_fractions_last)*ghostcell_dim_0_mass_fractions_last +
                            (k + offset_2_mass_fractions_last)*ghostcell_dim_0_mass_fractions_last*
                                ghostcell_dim_1_mass_fractions_last;
                        
                        const int idx_mass_fractions = (i + offset_0_mass_fractions) +
                            (j + offset_1_mass_fractions)*ghostcell_dim_0_mass_fractions +
                            (k + offset_2_mass_fractions)*ghostcell_dim_0_mass_fractions*
                                ghostcell_dim_1_mass_fractions;
                        
                        c_p[idx_isobaric_specific_heat_capacity] += Y[si][idx_mass_fractions]*d_species_c_p[si];
                        
                        // Compute the mass fraction of the last species.
                        Y_last[idx_mass_fractions_last] -= Y[si][idx_mass_fractions];
                    }
                }
            }
        }
        
        // Add the contribution from the last species.
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
                    const int idx_isobaric_specific_heat_capacity =
                        (i + offset_0_isobaric_specific_heat_capacity) +
                        (j + offset_1_isobaric_specific_heat_capacity)*
                            ghostcell_dim_0_isobaric_specific_heat_capacity +
                        (k + offset_2_isobaric_specific_heat_capacity)*
                            ghostcell_dim_0_isobaric_specific_heat_capacity*
                                ghostcell_dim_1_isobaric_specific_heat_capacity;
                    
                    const int idx_mass_fractions_last = (i + offset_0_mass_fractions_last) +
                        (j + offset_1_mass_fractions_last)*ghostcell_dim_0_mass_fractions_last +
                        (k + offset_2_mass_fractions_last)*ghostcell_dim_0_mass_fractions_last*
                            ghostcell_dim_1_mass_fractions_last;
                    
                    c_p[idx_isobaric_specific_heat_capacity] += Y_last[idx_mass_fractions_last]*d_species_c_p.back();
                }
            }
        }
    }
}


/*
 * Compute the mixture partial derivative of pressure w.r.t. partial densities under constant specific
 * internal energy with isothermal and isobaric equilibrium assumptions.
 */
void
EquationOfStateMixingRulesIdealGas::computePressureDerivativeWithPartialDensities(
    std::vector<double*>& Psi,
    const double* const epsilon,
    const double* const gamma,
    const double* const c_v,
    const hier::IntVector& offset_partial_pressure_partial_partial_densities,
    const hier::IntVector& offset_internal_energy,
    const hier::IntVector& offset_mixture_thermo_properties,
    const hier::IntVector& ghostcell_dims_partial_pressure_partial_partial_densities,
    const hier::IntVector& ghostcell_dims_internal_energy,
    const hier::IntVector& ghostcell_dims_mixture_thermo_properties,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims) const
{
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and offsets.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int offset_0_partial_pressure_partial_partial_densities =
            offset_partial_pressure_partial_partial_densities[0];
        const int offset_0_internal_energy = offset_internal_energy[0];
        const int offset_0_mixture_thermo_properties = offset_mixture_thermo_properties[0];
        
        // Compute Psi.
        for (int si = 0; si < d_num_species; si++)
        {
            double* Psi_i = Psi[si];
            
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_partial_pressure_partial_partial_densities =
                    i + offset_0_partial_pressure_partial_partial_densities;
                
                const int idx_internal_energy = i + offset_0_internal_energy;
                
                const int idx_mixture_thermo_properties = i + offset_0_mixture_thermo_properties;
                
                Psi_i[idx_partial_pressure_partial_partial_densities] =
                    ((d_species_c_p[si] - gamma[idx_mixture_thermo_properties]*d_species_c_v[si])/
                    c_v[idx_mixture_thermo_properties] + gamma[idx_mixture_thermo_properties] - double(1))*
                    epsilon[idx_internal_energy];
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
        
        const int offset_0_partial_pressure_partial_partial_densities =
            offset_partial_pressure_partial_partial_densities[0];
        const int offset_1_partial_pressure_partial_partial_densities =
            offset_partial_pressure_partial_partial_densities[1];
        const int ghostcell_dim_0_partial_pressure_partial_partial_densities =
            ghostcell_dims_partial_pressure_partial_partial_densities[0];
        
        const int offset_0_internal_energy = offset_internal_energy[0];
        const int offset_1_internal_energy = offset_internal_energy[1];
        const int ghostcell_dim_0_internal_energy = ghostcell_dims_internal_energy[0];
        
        const int offset_0_mixture_thermo_properties = offset_mixture_thermo_properties[0];
        const int offset_1_mixture_thermo_properties = offset_mixture_thermo_properties[1];
        const int ghostcell_dim_0_mixture_thermo_properties = ghostcell_dims_mixture_thermo_properties[0];
        
        // Compute Psi.
        for (int si = 0; si < d_num_species; si++)
        {
            double* Psi_i = Psi[si];
            
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_partial_pressure_partial_partial_densities =
                        (i + offset_0_partial_pressure_partial_partial_densities) +
                        (j + offset_1_partial_pressure_partial_partial_densities)*
                            ghostcell_dim_0_partial_pressure_partial_partial_densities;
                    
                    const int idx_internal_energy = (i + offset_0_internal_energy) +
                        (j + offset_1_internal_energy)*ghostcell_dim_0_internal_energy;
                    
                    const int idx_mixture_thermo_properties = (i + offset_0_mixture_thermo_properties) +
                        (j + offset_1_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties;
                    
                    Psi_i[idx_partial_pressure_partial_partial_densities] =
                        ((d_species_c_p[si] - gamma[idx_mixture_thermo_properties]*d_species_c_v[si])/
                        c_v[idx_mixture_thermo_properties] + gamma[idx_mixture_thermo_properties] - double(1))*
                        epsilon[idx_internal_energy];
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
        
        const int offset_0_partial_pressure_partial_partial_densities =
            offset_partial_pressure_partial_partial_densities[0];
        const int offset_1_partial_pressure_partial_partial_densities =
            offset_partial_pressure_partial_partial_densities[1];
        const int offset_2_partial_pressure_partial_partial_densities =
            offset_partial_pressure_partial_partial_densities[2];
        const int ghostcell_dim_0_partial_pressure_partial_partial_densities =
            ghostcell_dims_partial_pressure_partial_partial_densities[0];
        const int ghostcell_dim_1_partial_pressure_partial_partial_densities =
            ghostcell_dims_partial_pressure_partial_partial_densities[1];
        
        const int offset_0_internal_energy = offset_internal_energy[0];
        const int offset_1_internal_energy = offset_internal_energy[1];
        const int offset_2_internal_energy = offset_internal_energy[2];
        const int ghostcell_dim_0_internal_energy = ghostcell_dims_internal_energy[0];
        const int ghostcell_dim_1_internal_energy = ghostcell_dims_internal_energy[1];
        
        const int offset_0_mixture_thermo_properties = offset_mixture_thermo_properties[0];
        const int offset_1_mixture_thermo_properties = offset_mixture_thermo_properties[1];
        const int offset_2_mixture_thermo_properties = offset_mixture_thermo_properties[2];
        const int ghostcell_dim_0_mixture_thermo_properties = ghostcell_dims_mixture_thermo_properties[0];
        const int ghostcell_dim_1_mixture_thermo_properties = ghostcell_dims_mixture_thermo_properties[1];
        
        // Compute Psi.
        for (int si = 0; si < d_num_species; si++)
        {
            double* Psi_i = Psi[si];
            
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
                        const int idx_partial_pressure_partial_partial_densities =
                            (i + offset_0_partial_pressure_partial_partial_densities) +
                            (j + offset_1_partial_pressure_partial_partial_densities)*
                                ghostcell_dim_0_partial_pressure_partial_partial_densities +
                            (k + offset_2_partial_pressure_partial_partial_densities)*
                                ghostcell_dim_0_partial_pressure_partial_partial_densities*
                                ghostcell_dim_1_partial_pressure_partial_partial_densities;
                        
                        const int idx_internal_energy = (i + offset_0_internal_energy) +
                            (j + offset_1_internal_energy)*ghostcell_dim_0_internal_energy +
                            (k + offset_2_internal_energy)*ghostcell_dim_0_internal_energy*
                                ghostcell_dim_1_internal_energy;
                        
                        const int idx_mixture_thermo_properties = (i + offset_0_mixture_thermo_properties) +
                            (j + offset_1_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties +
                            (k + offset_2_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties*
                                ghostcell_dim_1_mixture_thermo_properties;
                        
                        Psi_i[idx_partial_pressure_partial_partial_densities] =
                            ((d_species_c_p[si] - gamma[idx_mixture_thermo_properties]*d_species_c_v[si])/
                            c_v[idx_mixture_thermo_properties] + gamma[idx_mixture_thermo_properties] - double(1))*
                            epsilon[idx_internal_energy];
                    }
                }
            }
        }
    }
}


/*
 * Compute the mixture partial derivative of pressure w.r.t. partial densities under constant specific
 * internal energy and volume fractions with isobaric equilibrium assumption.
 */
void
EquationOfStateMixingRulesIdealGas::computePressureDerivativeWithPartialDensities(
    std::vector<double*>& Psi,
    const double* const rho,
    const double* const p,
    const hier::IntVector& offset_partial_pressure_partial_partial_densities,
    const hier::IntVector& offset_density,
    const hier::IntVector& offset_pressure,
    const hier::IntVector& ghostcell_dims_partial_pressure_partial_partial_densities,
    const hier::IntVector& ghostcell_dims_density,
    const hier::IntVector& ghostcell_dims_pressure,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims) const
{
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and offsets.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int offset_0_partial_pressure_partial_partial_densities =
            offset_partial_pressure_partial_partial_densities[0];
        const int offset_0_density = offset_density[0];
        const int offset_0_pressure = offset_pressure[0];
        
        // Compute Psi.
        for (int si = 0; si < d_num_species; si++)
        {
            double* Psi_i = Psi[si];
            
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_partial_pressure_partial_partial_densities =
                    i + offset_0_partial_pressure_partial_partial_densities;
                
                const int idx_density = i + offset_0_density;
                
                const int idx_pressure = i + offset_0_pressure;
                
                Psi_i[idx_partial_pressure_partial_partial_densities] = p[idx_pressure]/rho[idx_density];
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
        
        const int offset_0_partial_pressure_partial_partial_densities =
            offset_partial_pressure_partial_partial_densities[0];
        const int offset_1_partial_pressure_partial_partial_densities =
            offset_partial_pressure_partial_partial_densities[1];
        const int ghostcell_dim_0_partial_pressure_partial_partial_densities =
            ghostcell_dims_partial_pressure_partial_partial_densities[0];
        
        const int offset_0_density = offset_density[0];
        const int offset_1_density = offset_density[1];
        const int ghostcell_dim_0_density = ghostcell_dims_density[0];
        
        const int offset_0_pressure = offset_pressure[0];
        const int offset_1_pressure = offset_pressure[1];
        const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
        
        // Compute Psi.
        for (int si = 0; si < d_num_species; si++)
        {
            double* Psi_i = Psi[si];
            
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_partial_pressure_partial_partial_densities =
                        (i + offset_0_partial_pressure_partial_partial_densities) +
                        (j + offset_1_partial_pressure_partial_partial_densities)*
                            ghostcell_dim_0_partial_pressure_partial_partial_densities;
                    
                    const int idx_density = (i + offset_0_density) +
                        (j + offset_1_density)*ghostcell_dim_0_density;
                    
                    const int idx_pressure = (i + offset_0_pressure) +
                        (j + offset_1_pressure)*ghostcell_dim_0_pressure;
                    
                    Psi_i[idx_partial_pressure_partial_partial_densities] = p[idx_pressure]/rho[idx_density];
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
        
        const int offset_0_partial_pressure_partial_partial_densities =
            offset_partial_pressure_partial_partial_densities[0];
        const int offset_1_partial_pressure_partial_partial_densities =
            offset_partial_pressure_partial_partial_densities[1];
        const int offset_2_partial_pressure_partial_partial_densities =
            offset_partial_pressure_partial_partial_densities[2];
        const int ghostcell_dim_0_partial_pressure_partial_partial_densities =
            ghostcell_dims_partial_pressure_partial_partial_densities[0];
        const int ghostcell_dim_1_partial_pressure_partial_partial_densities =
            ghostcell_dims_partial_pressure_partial_partial_densities[1];
        
        const int offset_0_density = offset_density[0];
        const int offset_1_density = offset_density[1];
        const int offset_2_density = offset_density[2];
        const int ghostcell_dim_0_density = ghostcell_dims_density[0];
        const int ghostcell_dim_1_density = ghostcell_dims_density[1];
        
        const int offset_0_pressure = offset_pressure[0];
        const int offset_1_pressure = offset_pressure[1];
        const int offset_2_pressure = offset_pressure[2];
        const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
        const int ghostcell_dim_1_pressure = ghostcell_dims_pressure[1];
        
        // Compute Psi.
        for (int si = 0; si < d_num_species; si++)
        {
            double* Psi_i = Psi[si];
            
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
                        const int idx_partial_pressure_partial_partial_densities =
                            (i + offset_0_partial_pressure_partial_partial_densities) +
                            (j + offset_1_partial_pressure_partial_partial_densities)*
                                ghostcell_dim_0_partial_pressure_partial_partial_densities +
                            (k + offset_2_partial_pressure_partial_partial_densities)*
                                ghostcell_dim_0_partial_pressure_partial_partial_densities*
                                ghostcell_dim_1_partial_pressure_partial_partial_densities;
                        
                        const int idx_density = (i + offset_0_density) +
                            (j + offset_1_density)*ghostcell_dim_0_density +
                            (k + offset_2_density)*ghostcell_dim_0_density*
                                ghostcell_dim_1_density;
                        
                        const int idx_pressure = (i + offset_0_pressure) +
                            (j + offset_1_pressure)*ghostcell_dim_0_pressure +
                            (k + offset_2_pressure)*ghostcell_dim_0_pressure*
                                ghostcell_dim_1_pressure;
                        
                        Psi_i[idx_partial_pressure_partial_partial_densities] = p[idx_pressure]/rho[idx_density];
                    }
                }
            }
        }
    }
}


/*
 * Compute the mixture partial derivative of pressure w.r.t. volume fractions under constant specific
 * internal energy and partial densities with isobaric equilibrium assumption.
 */
void
EquationOfStateMixingRulesIdealGas::computePressureDerivativeWithVolumeFractions(
    std::vector<double*>& M,
    const double* const p,
    const double* const gamma,
    const hier::IntVector& offset_partial_pressure_partial_volume_fractions,
    const hier::IntVector& offset_pressure,
    const hier::IntVector& offset_mixture_thermo_properties,
    const hier::IntVector& ghostcell_dims_partial_pressure_partial_volume_fractions,
    const hier::IntVector& ghostcell_dims_pressure,
    const hier::IntVector& ghostcell_dims_mixture_thermo_properties,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims) const
{
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and offsets.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int offset_0_partial_pressure_partial_volume_fractions =
            offset_partial_pressure_partial_volume_fractions[0];
        const int offset_0_pressure = offset_pressure[0];
        const int offset_0_mixture_thermo_properties = offset_mixture_thermo_properties[0];
        
        // Compute M.
        for (int si = 0; si < d_num_species - 1; si++)
        {
            double* M_i = M[si];
            
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_partial_pressure_partial_volume_fractions =
                    i + offset_0_partial_pressure_partial_volume_fractions;
                
                const int idx_pressure = i + offset_0_pressure;
                
                const int idx_mixture_thermo_properties = i + offset_0_mixture_thermo_properties;
                
                M_i[idx_partial_pressure_partial_volume_fractions] =
                    (double(1)/(d_species_gamma[d_num_species - 1] - double(1)) -
                     double(1)/(d_species_gamma[si] - double(1)))*
                    (gamma[idx_mixture_thermo_properties] - double(1))*p[idx_pressure];
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
        
        const int offset_0_partial_pressure_partial_volume_fractions =
            offset_partial_pressure_partial_volume_fractions[0];
        const int offset_1_partial_pressure_partial_volume_fractions =
            offset_partial_pressure_partial_volume_fractions[1];
        const int ghostcell_dim_0_partial_pressure_partial_volume_fractions =
            ghostcell_dims_partial_pressure_partial_volume_fractions[0];
        
        const int offset_0_pressure = offset_pressure[0];
        const int offset_1_pressure = offset_pressure[1];
        const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
        
        const int offset_0_mixture_thermo_properties = offset_mixture_thermo_properties[0];
        const int offset_1_mixture_thermo_properties = offset_mixture_thermo_properties[1];
        const int ghostcell_dim_0_mixture_thermo_properties = ghostcell_dims_mixture_thermo_properties[0];
        
        // Compute M.
        for (int si = 0; si < d_num_species - 1; si++)
        {
            double* M_i = M[si];
            
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_partial_pressure_partial_volume_fractions =
                        (i + offset_0_partial_pressure_partial_volume_fractions) +
                        (j + offset_1_partial_pressure_partial_volume_fractions)*
                            ghostcell_dim_0_partial_pressure_partial_volume_fractions;
                    
                    const int idx_pressure = (i + offset_0_pressure) +
                        (j + offset_1_pressure)*ghostcell_dim_0_pressure;
                    
                    const int idx_mixture_thermo_properties = (i + offset_0_mixture_thermo_properties) +
                        (j + offset_1_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties;
                    
                    M_i[idx_partial_pressure_partial_volume_fractions] =
                        (double(1)/(d_species_gamma[d_num_species - 1] - double(1)) -
                         double(1)/(d_species_gamma[si] - double(1)))*
                        (gamma[idx_mixture_thermo_properties] - double(1))*p[idx_pressure];
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
        
        const int offset_0_partial_pressure_partial_volume_fractions =
            offset_partial_pressure_partial_volume_fractions[0];
        const int offset_1_partial_pressure_partial_volume_fractions =
            offset_partial_pressure_partial_volume_fractions[1];
        const int offset_2_partial_pressure_partial_volume_fractions =
            offset_partial_pressure_partial_volume_fractions[2];
        const int ghostcell_dim_0_partial_pressure_partial_volume_fractions =
            ghostcell_dims_partial_pressure_partial_volume_fractions[0];
        const int ghostcell_dim_1_partial_pressure_partial_volume_fractions =
            ghostcell_dims_partial_pressure_partial_volume_fractions[1];
        
        const int offset_0_pressure = offset_pressure[0];
        const int offset_1_pressure = offset_pressure[1];
        const int offset_2_pressure = offset_pressure[2];
        const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
        const int ghostcell_dim_1_pressure = ghostcell_dims_pressure[1];
        
        const int offset_0_mixture_thermo_properties = offset_mixture_thermo_properties[0];
        const int offset_1_mixture_thermo_properties = offset_mixture_thermo_properties[1];
        const int offset_2_mixture_thermo_properties = offset_mixture_thermo_properties[2];
        const int ghostcell_dim_0_mixture_thermo_properties = ghostcell_dims_mixture_thermo_properties[0];
        const int ghostcell_dim_1_mixture_thermo_properties = ghostcell_dims_mixture_thermo_properties[1];
        
        // Compute M.
        for (int si = 0; si < d_num_species - 1; si++)
        {
            double* M_i = M[si];
            
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
                        const int idx_partial_pressure_partial_volume_fractions =
                            (i + offset_0_partial_pressure_partial_volume_fractions) +
                            (j + offset_1_partial_pressure_partial_volume_fractions)*
                                ghostcell_dim_0_partial_pressure_partial_volume_fractions +
                            (k + offset_2_partial_pressure_partial_volume_fractions)*
                                ghostcell_dim_0_partial_pressure_partial_volume_fractions*
                                ghostcell_dim_1_partial_pressure_partial_volume_fractions;
                        
                        const int idx_pressure = (i + offset_0_pressure) +
                            (j + offset_1_pressure)*ghostcell_dim_0_pressure +
                            (k + offset_2_pressure)*ghostcell_dim_0_pressure*
                                ghostcell_dim_1_pressure;
                        
                        const int idx_mixture_thermo_properties = (i + offset_0_mixture_thermo_properties) +
                            (j + offset_1_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties +
                            (k + offset_2_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties*
                                ghostcell_dim_1_mixture_thermo_properties;
                        
                        M_i[idx_partial_pressure_partial_volume_fractions] =
                            (double(1)/(d_species_gamma[d_num_species - 1] - double(1)) -
                             double(1)/(d_species_gamma[si] - double(1)))*
                            (gamma[idx_mixture_thermo_properties] - double(1))*p[idx_pressure];
                    }
                }
            }
        }
    }
}


/*
 * Compute the thermodynamic properties of the mixture with mass fractions.
 */
void
EquationOfStateMixingRulesIdealGas::computeMixtureThermodynamicPropertiesWithMassFractions(
    double* const gamma,
    double* const R,
    double* const c_p,
    double* const c_v,
    const std::vector<const double*>& Y,
    const hier::IntVector& offset_mixture_thermo_properties,
    const hier::IntVector& offset_mass_fractions,
    const hier::IntVector& ghostcell_dims_mixture_thermo_properties,
    const hier::IntVector& ghostcell_dims_mass_fractions,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims) const
{
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and offsets.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int offset_0_mixture_thermo_properties = offset_mixture_thermo_properties[0];
        const int offset_0_mass_fractions = offset_mass_fractions[0];
        
        // Compute c_p and c_v.
        for (int si = 0; si < d_num_species; si++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_mixture_thermo_properties = i + offset_0_mixture_thermo_properties;
                const int idx_mass_fractions = i + offset_0_mass_fractions;
                
                c_p[idx_mixture_thermo_properties] += Y[si][idx_mass_fractions]*d_species_c_p[si];
                c_v[idx_mixture_thermo_properties] += Y[si][idx_mass_fractions]*d_species_c_v[si];
            }
        }
        
        // Compute gamma and R.
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_mixture_thermo_properties = i + offset_0_mixture_thermo_properties;
            
            gamma[idx_mixture_thermo_properties] = c_p[idx_mixture_thermo_properties]/
                c_v[idx_mixture_thermo_properties];
            R[idx_mixture_thermo_properties] = c_p[idx_mixture_thermo_properties] -
                c_v[idx_mixture_thermo_properties];
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
        
        const int offset_0_mixture_thermo_properties = offset_mixture_thermo_properties[0];
        const int offset_1_mixture_thermo_properties = offset_mixture_thermo_properties[1];
        const int ghostcell_dim_0_mixture_thermo_properties = ghostcell_dims_mixture_thermo_properties[0];
        
        const int offset_0_mass_fractions = offset_mass_fractions[0];
        const int offset_1_mass_fractions = offset_mass_fractions[1];
        const int ghostcell_dim_0_mass_fractions = ghostcell_dims_mass_fractions[0];
        
        // Compute c_p and c_v.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_mixture_thermo_properties = (i + offset_0_mixture_thermo_properties) +
                        (j + offset_1_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties;
                    
                    const int idx_mass_fractions = (i + offset_0_mass_fractions) +
                        (j + offset_1_mass_fractions)*ghostcell_dim_0_mass_fractions;
                    
                    c_p[idx_mixture_thermo_properties] += Y[si][idx_mass_fractions]*d_species_c_p[si];
                    c_v[idx_mixture_thermo_properties] += Y[si][idx_mass_fractions]*d_species_c_v[si];
                }
            }
        }
        
        // Compute gamma and R.
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_mixture_thermo_properties = (i + offset_0_mixture_thermo_properties) +
                    (j + offset_1_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties;
                
                gamma[idx_mixture_thermo_properties] = c_p[idx_mixture_thermo_properties]/
                    c_v[idx_mixture_thermo_properties];
                R[idx_mixture_thermo_properties] = c_p[idx_mixture_thermo_properties] -
                    c_v[idx_mixture_thermo_properties];
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
        
        const int offset_0_mixture_thermo_properties = offset_mixture_thermo_properties[0];
        const int offset_1_mixture_thermo_properties = offset_mixture_thermo_properties[1];
        const int offset_2_mixture_thermo_properties = offset_mixture_thermo_properties[2];
        const int ghostcell_dim_0_mixture_thermo_properties = ghostcell_dims_mixture_thermo_properties[0];
        const int ghostcell_dim_1_mixture_thermo_properties = ghostcell_dims_mixture_thermo_properties[1];
        
        const int offset_0_mass_fractions = offset_mass_fractions[0];
        const int offset_1_mass_fractions = offset_mass_fractions[1];
        const int offset_2_mass_fractions = offset_mass_fractions[2];
        const int ghostcell_dim_0_mass_fractions = ghostcell_dims_mass_fractions[0];
        const int ghostcell_dim_1_mass_fractions = ghostcell_dims_mass_fractions[1];
        
        // Compute c_p and c_v.
        for (int si = 0; si < d_num_species; si++)
        {
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
                        const int idx_mixture_thermo_properties = (i + offset_0_mixture_thermo_properties) +
                            (j + offset_1_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties +
                            (k + offset_2_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties*
                                ghostcell_dim_1_mixture_thermo_properties;
                        
                        const int idx_mass_fractions = (i + offset_0_mass_fractions) +
                            (j + offset_1_mass_fractions)*ghostcell_dim_0_mass_fractions +
                            (k + offset_2_mass_fractions)*ghostcell_dim_0_mass_fractions*
                                ghostcell_dim_1_mass_fractions;
                        
                        c_p[idx_mixture_thermo_properties] += Y[si][idx_mass_fractions]*d_species_c_p[si];
                        c_v[idx_mixture_thermo_properties] += Y[si][idx_mass_fractions]*d_species_c_v[si];
                    }
                }
            }
        }
        
        // Compute gamma and R.
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
                    const int idx_mixture_thermo_properties = (i + offset_0_mixture_thermo_properties) +
                        (j + offset_1_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties +
                        (k + offset_2_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties*
                            ghostcell_dim_1_mixture_thermo_properties;
                    
                    gamma[idx_mixture_thermo_properties] = c_p[idx_mixture_thermo_properties]/
                        c_v[idx_mixture_thermo_properties];
                    R[idx_mixture_thermo_properties] = c_p[idx_mixture_thermo_properties] -
                        c_v[idx_mixture_thermo_properties];
                }
            }
        }
    }
}


/*
 * Compute the thermodynamic properties of the mixture with mass fractions.
 */
void
EquationOfStateMixingRulesIdealGas::computeMixtureThermodynamicPropertiesWithMassFractions(
    double* const gamma,
    double* const R,
    double* const c_p,
    double* const c_v,
    double* const Y_last,
    const std::vector<const double*>& Y,
    const hier::IntVector& offset_mixture_thermo_properties,
    const hier::IntVector& offset_mass_fractions_last,
    const hier::IntVector& offset_mass_fractions,
    const hier::IntVector& ghostcell_dims_mixture_thermo_properties,
    const hier::IntVector& ghostcell_dims_mass_fractions_last,
    const hier::IntVector& ghostcell_dims_mass_fractions,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims) const
{
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and offsets.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int offset_0_mixture_thermo_properties = offset_mixture_thermo_properties[0];
        const int offset_0_mass_fractions_last = offset_mass_fractions_last[0];
        const int offset_0_mass_fractions = offset_mass_fractions[0];
        
        // Compute c_p and c_v.
        for (int si = 0; si < d_num_species - 1; si++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_mixture_thermo_properties = i + offset_0_mixture_thermo_properties;
                const int idx_mass_fractions_last = i + offset_0_mass_fractions_last;
                const int idx_mass_fractions = i + offset_0_mass_fractions;
                
                c_p[idx_mixture_thermo_properties] += Y[si][idx_mass_fractions]*d_species_c_p[si];
                c_v[idx_mixture_thermo_properties] += Y[si][idx_mass_fractions]*d_species_c_v[si];
                
                // Compute the mass fraction of the last species.
                Y_last[idx_mass_fractions_last] -= Y[si][idx_mass_fractions];
            }
        }
        
        // Add the contribution from the last species and compute gamma and R.
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_mixture_thermo_properties = i + offset_0_mixture_thermo_properties;
            const int idx_mass_fractions_last = i + offset_0_mass_fractions_last;
            
            c_p[idx_mixture_thermo_properties] += Y_last[idx_mass_fractions_last]*d_species_c_p.back();
            c_v[idx_mixture_thermo_properties] += Y_last[idx_mass_fractions_last]*d_species_c_v.back();
            
            gamma[idx_mixture_thermo_properties] = c_p[idx_mixture_thermo_properties]/
                c_v[idx_mixture_thermo_properties];
            R[idx_mixture_thermo_properties] = c_p[idx_mixture_thermo_properties] -
                c_v[idx_mixture_thermo_properties];
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
        
        const int offset_0_mixture_thermo_properties = offset_mixture_thermo_properties[0];
        const int offset_1_mixture_thermo_properties = offset_mixture_thermo_properties[1];
        const int ghostcell_dim_0_mixture_thermo_properties = ghostcell_dims_mixture_thermo_properties[0];
        
        const int offset_0_mass_fractions_last = offset_mass_fractions_last[0];
        const int offset_1_mass_fractions_last = offset_mass_fractions_last[1];
        const int ghostcell_dim_0_mass_fractions_last = ghostcell_dims_mass_fractions_last[0];
        
        const int offset_0_mass_fractions = offset_mass_fractions[0];
        const int offset_1_mass_fractions = offset_mass_fractions[1];
        const int ghostcell_dim_0_mass_fractions = ghostcell_dims_mass_fractions[0];
        
        // Compute c_p and c_v.
        for (int si = 0; si < d_num_species - 1; si++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_mixture_thermo_properties = (i + offset_0_mixture_thermo_properties) +
                        (j + offset_1_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties;
                    
                    const int idx_mass_fractions_last = (i + offset_0_mass_fractions_last) +
                        (j + offset_1_mass_fractions_last)*ghostcell_dim_0_mass_fractions_last;
                    
                    const int idx_mass_fractions = (i + offset_0_mass_fractions) +
                        (j + offset_1_mass_fractions)*ghostcell_dim_0_mass_fractions;
                    
                    c_p[idx_mixture_thermo_properties] += Y[si][idx_mass_fractions]*d_species_c_p[si];
                    c_v[idx_mixture_thermo_properties] += Y[si][idx_mass_fractions]*d_species_c_v[si];
                    
                    // Compute the mass fraction of the last species.
                    Y_last[idx_mass_fractions_last] -= Y[si][idx_mass_fractions];
                }
            }
        }
        
        // Add the contribution from the last species and compute gamma and R.
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_mixture_thermo_properties = (i + offset_0_mixture_thermo_properties) +
                    (j + offset_1_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties;
                
                const int idx_mass_fractions_last = (i + offset_0_mass_fractions_last) +
                    (j + offset_1_mass_fractions_last)*ghostcell_dim_0_mass_fractions_last;
                
                c_p[idx_mixture_thermo_properties] += Y_last[idx_mass_fractions_last]*d_species_c_p.back();
                c_v[idx_mixture_thermo_properties] += Y_last[idx_mass_fractions_last]*d_species_c_v.back();
                
                gamma[idx_mixture_thermo_properties] = c_p[idx_mixture_thermo_properties]/
                    c_v[idx_mixture_thermo_properties];
                R[idx_mixture_thermo_properties] = c_p[idx_mixture_thermo_properties] -
                    c_v[idx_mixture_thermo_properties];
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
        
        const int offset_0_mixture_thermo_properties = offset_mixture_thermo_properties[0];
        const int offset_1_mixture_thermo_properties = offset_mixture_thermo_properties[1];
        const int offset_2_mixture_thermo_properties = offset_mixture_thermo_properties[2];
        const int ghostcell_dim_0_mixture_thermo_properties = ghostcell_dims_mixture_thermo_properties[0];
        const int ghostcell_dim_1_mixture_thermo_properties = ghostcell_dims_mixture_thermo_properties[1];
        
        const int offset_0_mass_fractions_last = offset_mass_fractions_last[0];
        const int offset_1_mass_fractions_last = offset_mass_fractions_last[1];
        const int offset_2_mass_fractions_last = offset_mass_fractions_last[2];
        const int ghostcell_dim_0_mass_fractions_last = ghostcell_dims_mass_fractions_last[0];
        const int ghostcell_dim_1_mass_fractions_last = ghostcell_dims_mass_fractions_last[1];
        
        const int offset_0_mass_fractions = offset_mass_fractions[0];
        const int offset_1_mass_fractions = offset_mass_fractions[1];
        const int offset_2_mass_fractions = offset_mass_fractions[2];
        const int ghostcell_dim_0_mass_fractions = ghostcell_dims_mass_fractions[0];
        const int ghostcell_dim_1_mass_fractions = ghostcell_dims_mass_fractions[1];
        
        // Compute c_p and c_v.
        for (int si = 0; si < d_num_species - 1; si++)
        {
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
                        const int idx_mixture_thermo_properties = (i + offset_0_mixture_thermo_properties) +
                            (j + offset_1_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties +
                            (k + offset_2_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties*
                                ghostcell_dim_1_mixture_thermo_properties;
                        
                        const int idx_mass_fractions_last = (i + offset_0_mass_fractions_last) +
                            (j + offset_1_mass_fractions_last)*ghostcell_dim_0_mass_fractions_last +
                            (k + offset_2_mass_fractions_last)*ghostcell_dim_0_mass_fractions_last*
                                ghostcell_dim_1_mass_fractions_last;
                        
                        const int idx_mass_fractions = (i + offset_0_mass_fractions) +
                            (j + offset_1_mass_fractions)*ghostcell_dim_0_mass_fractions +
                            (k + offset_2_mass_fractions)*ghostcell_dim_0_mass_fractions*
                                ghostcell_dim_1_mass_fractions;
                        
                        c_p[idx_mixture_thermo_properties] += Y[si][idx_mass_fractions]*d_species_c_p[si];
                        c_v[idx_mixture_thermo_properties] += Y[si][idx_mass_fractions]*d_species_c_v[si];
                        
                        // Compute the mass fraction of the last species.
                        Y_last[idx_mass_fractions_last] -= Y[si][idx_mass_fractions];
                    }
                }
            }
        }
        
        // Add the contribution from the last species and compute gamma and R.
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
                    const int idx_mixture_thermo_properties = (i + offset_0_mixture_thermo_properties) +
                        (j + offset_1_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties +
                        (k + offset_2_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties*
                            ghostcell_dim_1_mixture_thermo_properties;
                    
                    const int idx_mass_fractions_last = (i + offset_0_mass_fractions_last) +
                        (j + offset_1_mass_fractions_last)*ghostcell_dim_0_mass_fractions_last +
                        (k + offset_2_mass_fractions_last)*ghostcell_dim_0_mass_fractions_last*
                            ghostcell_dim_1_mass_fractions_last;
                    
                    c_p[idx_mixture_thermo_properties] += Y_last[idx_mass_fractions_last]*d_species_c_p.back();
                    c_v[idx_mixture_thermo_properties] += Y_last[idx_mass_fractions_last]*d_species_c_v.back();
                    
                    gamma[idx_mixture_thermo_properties] = c_p[idx_mixture_thermo_properties]/
                        c_v[idx_mixture_thermo_properties];
                    R[idx_mixture_thermo_properties] = c_p[idx_mixture_thermo_properties] -
                        c_v[idx_mixture_thermo_properties];
                }
            }
        }
    }
}


/*
 * Compute the thermodynamic properties of the mixture with volume fractions.
 */
void
EquationOfStateMixingRulesIdealGas::computeMixtureThermodynamicPropertiesWithVolumeFractions(
    double* const gamma,
    const std::vector<const double*>& Z,
    const hier::IntVector& offset_mixture_thermo_properties,
    const hier::IntVector& offset_volume_fractions,
    const hier::IntVector& ghostcell_dims_mixture_thermo_properties,
    const hier::IntVector& ghostcell_dims_volume_fractions,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims) const
{
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and offsets.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int offset_0_mixture_thermo_properties = offset_mixture_thermo_properties[0];
        const int offset_0_volume_fractions = offset_volume_fractions[0];
        
        // Compute xi and store it in the data of gamma temporarily.
        for (int si = 0; si < d_num_species; si++)
        {
            const double one_over_denominator = double(1)/(d_species_gamma[si] - double(1));
            
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_mixture_thermo_properties = i + offset_0_mixture_thermo_properties;
                const int idx_volume_fractions = i + offset_0_volume_fractions;
                
                gamma[idx_mixture_thermo_properties] += Z[si][idx_volume_fractions]*one_over_denominator;
            }
        }
        
        // Compute gamma.
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear index.
            const int idx_mixture_thermo_properties = i + offset_0_mixture_thermo_properties;
            
            gamma[idx_mixture_thermo_properties] = double(1)/gamma[idx_mixture_thermo_properties] + double(1);
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
        
        const int offset_0_mixture_thermo_properties = offset_mixture_thermo_properties[0];
        const int offset_1_mixture_thermo_properties = offset_mixture_thermo_properties[1];
        const int ghostcell_dim_0_mixture_thermo_properties = ghostcell_dims_mixture_thermo_properties[0];
        
        const int offset_0_volume_fractions = offset_volume_fractions[0];
        const int offset_1_volume_fractions = offset_volume_fractions[1];
        const int ghostcell_dim_0_volume_fractions = ghostcell_dims_volume_fractions[0];
        
        // Compute xi and store it in the data of gamma temporarily.
        for (int si = 0; si < d_num_species; si++)
        {
            const double one_over_denominator = double(1)/(d_species_gamma[si] - double(1));
            
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_mixture_thermo_properties = (i + offset_0_mixture_thermo_properties) +
                        (j + offset_1_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties;
                    
                    const int idx_volume_fractions = (i + offset_0_volume_fractions) +
                        (j + offset_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                    
                    gamma[idx_mixture_thermo_properties] += Z[si][idx_volume_fractions]*one_over_denominator;
                }
            }
        }
        
        // Compute gamma.
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear index.
                const int idx_mixture_thermo_properties = (i + offset_0_mixture_thermo_properties) +
                    (j + offset_1_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties;
                
                gamma[idx_mixture_thermo_properties] = double(1)/gamma[idx_mixture_thermo_properties] + double(1);
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
        
        const int offset_0_mixture_thermo_properties = offset_mixture_thermo_properties[0];
        const int offset_1_mixture_thermo_properties = offset_mixture_thermo_properties[1];
        const int offset_2_mixture_thermo_properties = offset_mixture_thermo_properties[2];
        const int ghostcell_dim_0_mixture_thermo_properties = ghostcell_dims_mixture_thermo_properties[0];
        const int ghostcell_dim_1_mixture_thermo_properties = ghostcell_dims_mixture_thermo_properties[1];
        
        const int offset_0_volume_fractions = offset_volume_fractions[0];
        const int offset_1_volume_fractions = offset_volume_fractions[1];
        const int offset_2_volume_fractions = offset_volume_fractions[2];
        const int ghostcell_dim_0_volume_fractions = ghostcell_dims_volume_fractions[0];
        const int ghostcell_dim_1_volume_fractions = ghostcell_dims_volume_fractions[1];
        
        // Compute xi and store it in the data of gamma temporarily.
        for (int si = 0; si < d_num_species; si++)
        {
            const double one_over_denominator = double(1)/(d_species_gamma[si] - double(1));
            
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
                        const int idx_mixture_thermo_properties = (i + offset_0_mixture_thermo_properties) +
                            (j + offset_1_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties +
                            (k + offset_2_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties*
                                ghostcell_dim_1_mixture_thermo_properties;
                        
                        const int idx_volume_fractions = (i + offset_0_volume_fractions) +
                            (j + offset_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                            (k + offset_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                ghostcell_dim_1_volume_fractions;
                        
                        gamma[idx_mixture_thermo_properties] += Z[si][idx_volume_fractions]*one_over_denominator;
                    }
                }
            }
        }
        
        // Compute gamma.
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_mixture_thermo_properties = (i + offset_0_mixture_thermo_properties) +
                        (j + offset_1_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties +
                        (k + offset_2_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties*
                            ghostcell_dim_1_mixture_thermo_properties;
                    
                    gamma[idx_mixture_thermo_properties] = double(1)/gamma[idx_mixture_thermo_properties] + double(1);
                }
            }
        }
    }
}


/*
 * Compute the thermodynamic properties of the mixture with volume fractions.
 */
void
EquationOfStateMixingRulesIdealGas::computeMixtureThermodynamicPropertiesWithVolumeFractions(
    double* const gamma,
    double* const Z_last,
    const std::vector<const double*>& Z,
    const hier::IntVector& offset_mixture_thermo_properties,
    const hier::IntVector& offset_volume_fractions_last,
    const hier::IntVector& offset_volume_fractions,
    const hier::IntVector& ghostcell_dims_mixture_thermo_properties,
    const hier::IntVector& ghostcell_dims_volume_fractions_last,
    const hier::IntVector& ghostcell_dims_volume_fractions,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims) const
{
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and offsets.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int offset_0_mixture_thermo_properties = offset_mixture_thermo_properties[0];
        const int offset_0_volume_fractions_last = offset_volume_fractions_last[0];
        const int offset_0_volume_fractions = offset_volume_fractions[0];
        
        // Compute xi and store it in the data of gamma temporarily.
        for (int si = 0; si < d_num_species - 1; si++)
        {
            const double one_over_denominator = double(1)/(d_species_gamma[si] - double(1));
            
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_mixture_thermo_properties = i + offset_0_mixture_thermo_properties;
                const int idx_volume_fractions_last = i + offset_0_volume_fractions_last;
                const int idx_volume_fractions = i + offset_0_volume_fractions;
                
                gamma[idx_mixture_thermo_properties] += Z[si][idx_volume_fractions]*one_over_denominator;
                
                // Compute the volume fraction of the last species.
                Z_last[idx_volume_fractions_last] -= Z[si][idx_volume_fractions];
            }
        }
        
        // Add the contribution from the last species and compute gamma.
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_mixture_thermo_properties = i + offset_0_mixture_thermo_properties;
            const int idx_volume_fractions_last = i + offset_0_volume_fractions_last;
            
            gamma[idx_mixture_thermo_properties] += Z_last[idx_volume_fractions_last]/
                (d_species_gamma.back() - double(1));
            gamma[idx_mixture_thermo_properties] = double(1)/gamma[idx_mixture_thermo_properties] + double(1);
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
        
        const int offset_0_mixture_thermo_properties = offset_mixture_thermo_properties[0];
        const int offset_1_mixture_thermo_properties = offset_mixture_thermo_properties[1];
        const int ghostcell_dim_0_mixture_thermo_properties = ghostcell_dims_mixture_thermo_properties[0];
        
        const int offset_0_volume_fractions_last = offset_volume_fractions_last[0];
        const int offset_1_volume_fractions_last = offset_volume_fractions_last[1];
        const int ghostcell_dim_0_volume_fractions_last = ghostcell_dims_volume_fractions_last[0];
        
        const int offset_0_volume_fractions = offset_volume_fractions[0];
        const int offset_1_volume_fractions = offset_volume_fractions[1];
        const int ghostcell_dim_0_volume_fractions = ghostcell_dims_volume_fractions[0];
        
        // Compute xi and store it in the data of gamma temporarily.
        for (int si = 0; si < d_num_species - 1; si++)
        {
            const double one_over_denominator = double(1)/(d_species_gamma[si] - double(1));
            
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_mixture_thermo_properties = (i + offset_0_mixture_thermo_properties) +
                        (j + offset_1_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties;
                    
                    const int idx_volume_fractions_last = (i + offset_0_volume_fractions_last) +
                        (j + offset_1_volume_fractions_last)*ghostcell_dim_0_volume_fractions_last;
                    
                    const int idx_volume_fractions = (i + offset_0_volume_fractions) +
                        (j + offset_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                    
                    gamma[idx_mixture_thermo_properties] += Z[si][idx_volume_fractions]*one_over_denominator;
                    
                    // Compute the volume fraction of the last species.
                    Z_last[idx_volume_fractions_last] -= Z[si][idx_volume_fractions];
                }
            }
        }
        
        // Add the contribution from the last species and compute gamma.
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_mixture_thermo_properties = (i + offset_0_mixture_thermo_properties) +
                    (j + offset_1_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties;
                
                const int idx_volume_fractions_last = (i + offset_0_volume_fractions_last) +
                    (j + offset_1_volume_fractions_last)*ghostcell_dim_0_volume_fractions_last;
                
                gamma[idx_mixture_thermo_properties] += Z_last[idx_volume_fractions_last]/
                    (d_species_gamma.back() - double(1));
                gamma[idx_mixture_thermo_properties] = double(1)/gamma[idx_mixture_thermo_properties] + double(1);
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
        
        const int offset_0_mixture_thermo_properties = offset_mixture_thermo_properties[0];
        const int offset_1_mixture_thermo_properties = offset_mixture_thermo_properties[1];
        const int offset_2_mixture_thermo_properties = offset_mixture_thermo_properties[2];
        const int ghostcell_dim_0_mixture_thermo_properties = ghostcell_dims_mixture_thermo_properties[0];
        const int ghostcell_dim_1_mixture_thermo_properties = ghostcell_dims_mixture_thermo_properties[1];
        
        const int offset_0_volume_fractions_last = offset_volume_fractions_last[0];
        const int offset_1_volume_fractions_last = offset_volume_fractions_last[1];
        const int offset_2_volume_fractions_last = offset_volume_fractions_last[2];
        const int ghostcell_dim_0_volume_fractions_last = ghostcell_dims_volume_fractions_last[0];
        const int ghostcell_dim_1_volume_fractions_last = ghostcell_dims_volume_fractions_last[1];
        
        const int offset_0_volume_fractions = offset_volume_fractions[0];
        const int offset_1_volume_fractions = offset_volume_fractions[1];
        const int offset_2_volume_fractions = offset_volume_fractions[2];
        const int ghostcell_dim_0_volume_fractions = ghostcell_dims_volume_fractions[0];
        const int ghostcell_dim_1_volume_fractions = ghostcell_dims_volume_fractions[1];
        
        // Compute xi and store it in the data of gamma temporarily.
        for (int si = 0; si < d_num_species - 1; si++)
        {
            const double one_over_denominator = double(1)/(d_species_gamma[si] - double(1));
            
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
                        const int idx_mixture_thermo_properties = (i + offset_0_mixture_thermo_properties) +
                            (j + offset_1_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties +
                            (k + offset_2_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties*
                                ghostcell_dim_1_mixture_thermo_properties;
                        
                        const int idx_volume_fractions_last = (i + offset_0_volume_fractions_last) +
                            (j + offset_1_volume_fractions_last)*ghostcell_dim_0_volume_fractions_last +
                            (k + offset_2_volume_fractions_last)*ghostcell_dim_0_volume_fractions_last*
                                ghostcell_dim_1_volume_fractions_last;
                        
                        const int idx_volume_fractions = (i + offset_0_volume_fractions) +
                            (j + offset_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                            (k + offset_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                ghostcell_dim_1_volume_fractions;
                        
                        gamma[idx_mixture_thermo_properties] += Z[si][idx_volume_fractions]*one_over_denominator;
                        
                        // Compute the volume fraction of the last species.
                        Z_last[idx_volume_fractions_last] -= Z[si][idx_volume_fractions];
                    }
                }
            }
        }
        
        // Add the contribution from the last species and compute gamma.
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
                    const int idx_mixture_thermo_properties = (i + offset_0_mixture_thermo_properties) +
                        (j + offset_1_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties +
                        (k + offset_2_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties*
                            ghostcell_dim_1_mixture_thermo_properties;
                    
                    const int idx_volume_fractions_last = (i + offset_0_volume_fractions_last) +
                        (j + offset_1_volume_fractions_last)*ghostcell_dim_0_volume_fractions_last +
                        (k + offset_2_volume_fractions_last)*ghostcell_dim_0_volume_fractions_last*
                            ghostcell_dim_1_volume_fractions_last;
                    
                    gamma[idx_mixture_thermo_properties] += Z_last[idx_volume_fractions_last]/
                        (d_species_gamma.back() - double(1));
                    gamma[idx_mixture_thermo_properties] = double(1)/gamma[idx_mixture_thermo_properties] + double(1);
                }
            }
        }
    }
}
