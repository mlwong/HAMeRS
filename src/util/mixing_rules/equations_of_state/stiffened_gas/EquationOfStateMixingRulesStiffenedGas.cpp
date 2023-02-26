#include "util/mixing_rules/equations_of_state/stiffened_gas/EquationOfStateMixingRulesStiffenedGas.hpp"

EquationOfStateMixingRulesStiffenedGas::EquationOfStateMixingRulesStiffenedGas(
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
    d_equation_of_state.reset(new EquationOfStateStiffenedGas(
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
            d_species_gamma = equation_of_state_mixing_rules_db->getRealVector("species_gamma");
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
            d_species_gamma = equation_of_state_mixing_rules_db->getRealVector("d_species_gamma");
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
     * Get the reference pressure of each species from the database.
     */
    
    if (equation_of_state_mixing_rules_db->keyExists("species_p_inf"))
    {
        size_t species_p_inf_array_size = equation_of_state_mixing_rules_db->getArraySize("species_p_inf");
        if (static_cast<int>(species_p_inf_array_size) == d_num_species)
        {
            d_species_p_inf = equation_of_state_mixing_rules_db->getRealVector("species_p_inf");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'species_p_inf' entries must be equal to 'num_species'."
                << std::endl);
        }
    }
    else if (equation_of_state_mixing_rules_db->keyExists("d_species_p_inf"))
    {
        size_t species_p_inf_array_size = equation_of_state_mixing_rules_db->getArraySize("d_species_p_inf");
        if (static_cast<int>(species_p_inf_array_size) == d_num_species)
        {
            d_species_p_inf = equation_of_state_mixing_rules_db->getRealVector("d_species_p_inf");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'd_species_p_inf' entries must be equal to 'd_num_species'."
                << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Key data 'species_p_inf'/'d_species_p_inf'"
            << "not found in data for Equation_of_state_mixing_rules."
            << std::endl);
    }
    
    /*
     * Get the specific heat at constant volume of each species from the database.
     */
    
    if (equation_of_state_mixing_rules_db->keyExists("species_c_v"))
    {
        size_t species_c_v_array_size = equation_of_state_mixing_rules_db->getArraySize("species_c_v");
        if (static_cast<int>(species_c_v_array_size) == d_num_species)
        {
            d_species_c_v = equation_of_state_mixing_rules_db->getRealVector("species_c_v");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'species_c_v' entries must be equal to 'num_species'."
                << std::endl);
        }
    }
    else if (equation_of_state_mixing_rules_db->keyExists("d_species_c_v"))
    {
        size_t species_c_v_array_size = equation_of_state_mixing_rules_db->getArraySize("d_species_c_v");
        if (static_cast<int>(species_c_v_array_size) == d_num_species)
        {
            d_species_c_v = equation_of_state_mixing_rules_db->getRealVector("d_species_c_v");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'd_species_c_v' entries must be equal to 'd_num_species'."
                << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Key data 'species_c_v'/'d_species_c_v'"
            << "not found in data for Equation_of_state_mixing_rules."
            << std::endl);
    }
    
    /*
     * Get the molcular weight of each species from the database.
     */
    
    if (equation_of_state_mixing_rules_db->keyExists("species_M"))
    {
        size_t species_M_array_size = equation_of_state_mixing_rules_db->getArraySize("species_M");
        if (static_cast<int>(species_M_array_size) == d_num_species)
        {
            d_species_M = equation_of_state_mixing_rules_db->getRealVector("species_M");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'species_M' entries must be equal to 'num_species'."
                << std::endl);
        }
    }
    else if (equation_of_state_mixing_rules_db->keyExists("d_species_M"))
    {
        size_t species_M_array_size = equation_of_state_mixing_rules_db->getArraySize("d_species_M");
        if (static_cast<int>(species_M_array_size) == d_num_species)
        {
            d_species_M = equation_of_state_mixing_rules_db->getRealVector("d_species_M");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'd_species_M' entries must be equal to 'd_num_species'."
                << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Key data 'species_M'/'d_species_M'"
            << "not found in data for Equation_of_state_mixing_rules."
            << std::endl);
    }
    
    /*
     * Compute the specific heat at constant pressure of each species.
     */
    
    d_species_c_p.reserve(d_num_species);
    for (int si = 0; si < d_num_species; si++)
    {
        d_species_c_p.push_back(d_species_gamma[si]*d_species_c_v[si]);
    }
}


/*
 * Print all characteristics of the equation of state class.
 */
void
EquationOfStateMixingRulesStiffenedGas::printClassData(
    std::ostream& os) const
{
    os << "\nPrint EquationOfStateMixingRulesStiffenedGas object..."
       << std::endl;
    
    os << std::endl;
    os << "EquationOfStateMixingRulesStiffenedGas: this = "
       << (EquationOfStateMixingRulesStiffenedGas *)this
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
     * Print the reference pressure of each species.
     */
    
    os << "d_species_p_inf = ";
    for (int si = 0; si < d_num_species - 1; si++)
    {
        os << d_species_p_inf[si] << ", ";
    }
    os << d_species_p_inf[d_num_species - 1];
    os << std::endl;
    
    /*
     * Print the specific heat at constant volume of each species.
     */
    
    os << "d_species_c_v = ";
    for (int si = 0; si < d_num_species - 1; si++)
    {
        os << d_species_c_v[si] << ", ";
    }
    os << d_species_c_v[d_num_species - 1];
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
 * Put the characteristics of the equation of state mixing rules class into the restart
 * database.
 */
void
EquationOfStateMixingRulesStiffenedGas::putToRestart(
    const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const
{
    restart_db->putRealVector("d_species_gamma", d_species_gamma);
    restart_db->putRealVector("d_species_p_inf", d_species_p_inf);
    restart_db->putRealVector("d_species_c_v", d_species_c_v);
    restart_db->putRealVector("d_species_M", d_species_M);
}


/*
 * Compute the pressure of the mixture with isothermal and isobaric equilibrium assumptions.
 */
Real
EquationOfStateMixingRulesStiffenedGas::getPressure(
    const Real* const density,
    const Real* const internal_energy,
    const std::vector<const Real*>& mass_fractions) const
{
    TBOX_ERROR(d_object_name
        << ": EquationOfStateMixingRulesStiffenedGas::getPressure()\n"
        << "Method getPressure() for mixture with isothermal and isobaric equilibrium assumptions"
        << " is not yet implemented."
        << std::endl);
    
    return Real(0);
}


/*
 * Compute the pressure of the mixture with isothermal and isobaric equilibrium assumptions.
 */
void
EquationOfStateMixingRulesStiffenedGas::computePressure(
    HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_internal_energy,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_mass_fractions,
    const hier::Box& domain) const
{
    TBOX_ERROR(d_object_name
        << ": EquationOfStateMixingRulesStiffenedGas::computePressure()\n"
        << "Method computePressure() for mixture with isothermal and isobaric equilibrium assumptions"
        << " is not yet implemented."
        << std::endl);
}


/*
 * Compute the pressure of the mixture with isothermal and isobaric equilibrium assumptions.
 */
void
EquationOfStateMixingRulesStiffenedGas::computePressure(
    HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_internal_energy,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_mass_fractions,
    int side_normal,
    const hier::Box& domain) const
{
    TBOX_ERROR(d_object_name
        << ": EquationOfStateMixingRulesStiffenedGas::computePressure()\n"
        << "Method computePressure() for mixture with isothermal and isobaric equilibrium assumptions"
        << " is not yet implemented."
        << std::endl);
}


/*
 * Compute the pressure of the mixture with isobaric equilibrium assumption.
 */
Real
EquationOfStateMixingRulesStiffenedGas::getPressure(
    const Real* const density,
    const Real* const internal_energy,
    const std::vector<const Real*>& mass_fractions,
    const std::vector<const Real*>& volume_fractions) const
{
    NULL_USE(mass_fractions);
    
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOBARIC);
    TBOX_ASSERT((static_cast<int>(volume_fractions.size()) == d_num_species) ||
                (static_cast<int>(volume_fractions.size()) == d_num_species - 1));
#endif
    
    // Get the mixture thermodynamic properties.
    std::vector<Real> mixture_thermo_properties;
    std::vector<Real*> mixture_thermo_properties_ptr;
    std::vector<const Real*> mixture_thermo_properties_const_ptr;
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    mixture_thermo_properties.resize(num_thermo_properties);
    mixture_thermo_properties_ptr.reserve(num_thermo_properties);
    mixture_thermo_properties_const_ptr.reserve(num_thermo_properties);
    
    for (int ti = 0; ti < num_thermo_properties; ti++)
    {
        mixture_thermo_properties_ptr.push_back(&mixture_thermo_properties[ti]);
        mixture_thermo_properties_const_ptr.push_back(&mixture_thermo_properties[ti]);
    }
    
    getMixtureThermodynamicPropertiesWithVolumeFractions(
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
EquationOfStateMixingRulesStiffenedGas::computePressure(
    HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_internal_energy,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_mass_fractions,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_volume_fractions,
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
    HAMERS_SHARED_PTR<pdat::CellData<Real> > data_mixture_thermo_properties;
    
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
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::CellData<Real> >(
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
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::CellData<Real> >(
            domain, num_thermo_properties, hier::IntVector::getZero(d_dim));
    }
    
    computeMixtureThermodynamicPropertiesWithVolumeFractions(
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
EquationOfStateMixingRulesStiffenedGas::computePressure(
    HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_internal_energy,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_mass_fractions,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_volume_fractions,
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
    HAMERS_SHARED_PTR<pdat::SideData<Real> > data_mixture_thermo_properties;
    
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
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::SideData<Real> >(
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
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::SideData<Real> >(
            domain, num_thermo_properties, hier::IntVector::getZero(d_dim), direction);
    }
    
    computeMixtureThermodynamicPropertiesWithVolumeFractions(
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
Real
EquationOfStateMixingRulesStiffenedGas::getInternalEnergy(
    const Real* const density,
    const Real* const pressure,
    const std::vector<const Real*>& mass_fractions) const
{
    TBOX_ERROR(d_object_name
        << ": EquationOfStateMixingRulesStiffenedGas::getInternalEnergy()\n"
        << "Method getInternalEnergy() for mixture with isothermal and isobaric equilibrium assumptions"
        << " is not yet implemented."
        << std::endl);
    
    return Real(0);
}


/*
 * Compute the specific internal energy of the mixture with isothermal and isobaric equilibrium assumptions.
 */
void
EquationOfStateMixingRulesStiffenedGas::computeInternalEnergy(
    HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_internal_energy,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_mass_fractions,
    const hier::Box& domain) const
{
    TBOX_ERROR(d_object_name
        << ": EquationOfStateMixingRulesStiffenedGas::computeInternalEnergy()\n"
        << "Method computeInternalEnergy() for mixture with isothermal and isobaric equilibrium assumptions"
        << " is not yet implemented."
        << std::endl);
}


/*
 * Compute the specific internal energy of the mixture with isothermal and isobaric equilibrium assumptions.
 */
void
EquationOfStateMixingRulesStiffenedGas::computeInternalEnergy(
    HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_internal_energy,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_mass_fractions,
    int side_normal,
    const hier::Box& domain) const
{
    TBOX_ERROR(d_object_name
        << ": EquationOfStateMixingRulesStiffenedGas::computeInternalEnergy()\n"
        << "Method computeInternalEnergy() for mixture with isothermal and isobaric equilibrium assumptions"
        << " is not yet implemented."
        << std::endl);
}


/*
 * Compute the specific internal energy of the mixture with isobaric equilibrium assumption.
 */
Real
EquationOfStateMixingRulesStiffenedGas::getInternalEnergy(
    const Real* const density,
    const Real* const pressure,
    const std::vector<const Real*>& mass_fractions,
    const std::vector<const Real*>& volume_fractions) const
{
    NULL_USE(mass_fractions);
    
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOBARIC);
    TBOX_ASSERT((static_cast<int>(volume_fractions.size()) == d_num_species) ||
                (static_cast<int>(volume_fractions.size()) == d_num_species - 1));
#endif
    
    // Get the mixture thermodynamic properties.
    std::vector<Real> mixture_thermo_properties;
    std::vector<Real*> mixture_thermo_properties_ptr;
    std::vector<const Real*> mixture_thermo_properties_const_ptr;
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    mixture_thermo_properties.resize(num_thermo_properties);
    mixture_thermo_properties_ptr.reserve(num_thermo_properties);
    mixture_thermo_properties_const_ptr.reserve(num_thermo_properties);
    
    for (int ti = 0; ti < num_thermo_properties; ti++)
    {
        mixture_thermo_properties_ptr.push_back(&mixture_thermo_properties[ti]);
        mixture_thermo_properties_const_ptr.push_back(&mixture_thermo_properties[ti]);
    }
    
    getMixtureThermodynamicPropertiesWithVolumeFractions(
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
EquationOfStateMixingRulesStiffenedGas::computeInternalEnergy(
    HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_internal_energy,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_mass_fractions,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_volume_fractions,
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
    HAMERS_SHARED_PTR<pdat::CellData<Real> > data_mixture_thermo_properties;
    
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
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::CellData<Real> >(
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
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::CellData<Real> >(
            domain, num_thermo_properties, hier::IntVector::getZero(d_dim));
    }
    
    computeMixtureThermodynamicPropertiesWithVolumeFractions(
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
EquationOfStateMixingRulesStiffenedGas::computeInternalEnergy(
    HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_internal_energy,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_mass_fractions,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_volume_fractions,
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
    HAMERS_SHARED_PTR<pdat::SideData<Real> > data_mixture_thermo_properties;
    
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
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::SideData<Real> >(
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
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::SideData<Real> >(
            domain, num_thermo_properties, hier::IntVector::getZero(d_dim), direction);
    }
    
    computeMixtureThermodynamicPropertiesWithVolumeFractions(
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
Real
EquationOfStateMixingRulesStiffenedGas::getTemperature(
    const Real* const density,
    const Real* const pressure,
    const std::vector<const Real*>& mass_fractions) const
{
    TBOX_ERROR(d_object_name
        << ": EquationOfStateMixingRulesStiffenedGas::getTemperature()\n"
        << "Method getTemperature() for mixture with isothermal and isobaric equilibrium assumptions"
        << " is not yet implemented."
        << std::endl);
    
    return Real(0);
}


/*
 * Compute the temperature of the mixture with isothermal and isobaric equilibrium assumptions.
 */
void
EquationOfStateMixingRulesStiffenedGas::computeTemperature(
    HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_temperature,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_mass_fractions,
    const hier::Box& domain) const
{
    TBOX_ERROR(d_object_name
        << ": EquationOfStateMixingRulesStiffenedGas::computeTemperature()\n"
        << "Method computeTemperature() for mixture with isothermal and isobaric equilibrium assumptions"
        << " is not yet implemented."
        << std::endl);
}


/*
 * Compute the temperature of the mixture with isothermal and isobaric equilibrium assumptions.
 */
void
EquationOfStateMixingRulesStiffenedGas::computeTemperature(
    HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_temperature,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_mass_fractions,
    int side_normal,
    const hier::Box& domain) const
{
    TBOX_ERROR(d_object_name
        << ": EquationOfStateMixingRulesStiffenedGas::computeTemperature()\n"
        << "Method computeTemperature() for mixture with isothermal and isobaric equilibrium assumptions"
        << " is not yet implemented."
        << std::endl);
}


/*
 * Compute the specific internal energy of the mixture from temperature with isothermal
 * and isobaric equilibrium assumptions.
 */
Real
EquationOfStateMixingRulesStiffenedGas::getInternalEnergyFromTemperature(
    const Real* const density,
    const Real* const temperature,
    const std::vector<const Real*>& mass_fractions) const
{
    TBOX_ERROR(d_object_name
        << ": EquationOfStateMixingRulesStiffenedGas::getInternalEnergyFromTemperature()\n"
        << "Method getInternalEnergyFromTemperature() for mixture with isothermal and isobaric equilibrium assumptions"
        << " is not yet implemented."
        << std::endl);
    
    return Real(0);
}


/*
 * Compute the specific internal energy of the mixture from temperature with isothermal
 * and isobaric equilibrium assumptions.
 */
void
EquationOfStateMixingRulesStiffenedGas::computeInternalEnergyFromTemperature(
    HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_internal_energy,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_temperature,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_mass_fractions,
    const hier::Box& domain) const
{
    TBOX_ERROR(d_object_name
        << ": EquationOfStateMixingRulesStiffenedGas::computeInternalEnergyFromTemperature()\n"
        << "Method computeInternalEnergyFromTemperature() for mixture with isothermal and isobaric equilibrium"
        << " assumptions is not yet implemented."
        << std::endl);
}


/*
 * Compute the specific internal energy of the mixture from temperature with isothermal
 * and isobaric equilibrium assumptions.
 */
void
EquationOfStateMixingRulesStiffenedGas::computeInternalEnergyFromTemperature(
    HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_internal_energy,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_temperature,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_mass_fractions,
    int side_normal,
    const hier::Box& domain) const
{
    TBOX_ERROR(d_object_name
        << ": EquationOfStateMixingRulesStiffenedGas::computeInternalEnergyFromTemperature()\n"
        << "Method computeInternalEnergyFromTemperature() for mixture with isothermal and isobaric equilibrium"
        << " assumptions is not yet implemented."
        << std::endl);
}


/*
 * Compute the isochoric specific heat capacity of mixture with isothermal and isobaric equilibrium
 * assumptions.
 */
Real
EquationOfStateMixingRulesStiffenedGas::getIsochoricSpecificHeatCapacity(
    const Real* const density,
    const Real* const pressure,
    const std::vector<const Real*>& mass_fractions) const
{
    TBOX_ERROR(d_object_name
        << ": EquationOfStateMixingRulesStiffenedGas::getIsochoricSpecificHeatCapacity()\n"
        << "Method getIsochoricSpecificHeatCapacity() for mixture with isothermal and isobaric equilibrium assumptions"
        << " is not yet implemented."
        << std::endl);
    
    return Real(0);
}


/*
 * Compute the isochoric specific heat capacity of mixture with isothermal and isobaric equilibrium
 * assumptions.
 */
void
EquationOfStateMixingRulesStiffenedGas::computeIsochoricSpecificHeatCapacity(
    HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_isochoric_specific_heat_capacity,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_mass_fractions,
    const hier::Box& domain) const
{
    TBOX_ERROR(d_object_name
        << ": EquationOfStateMixingRulesStiffenedGas::computeIsochoricSpecificHeatCapacity()\n"
        << "Method computeIsochoricSpecificHeatCapacity() for mixture with isothermal and isobaric equilibrium"
        << " assumptions is not yet implemented."
        << std::endl);
}


/*
 * Compute the isochoric specific heat capacity of mixture with isothermal and isobaric equilibrium
 * assumptions.
 */
void
EquationOfStateMixingRulesStiffenedGas::computeIsochoricSpecificHeatCapacity(
    HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_isochoric_specific_heat_capacity,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_mass_fractions,
    int side_normal,
    const hier::Box& domain) const
{
    TBOX_ERROR(d_object_name
        << ": EquationOfStateMixingRulesStiffenedGas::computeIsochoricSpecificHeatCapacity()\n"
        << "Method computeIsochoricSpecificHeatCapacity() for mixture with isothermal and isobaric equilibrium"
        << " assumptions is not yet implemented."
        << std::endl);
}


/*
 * Compute the isobaric specific heat capacity of mixture with isothermal and isobaric equilibrium
 * assumptions.
 */
Real
EquationOfStateMixingRulesStiffenedGas::getIsobaricSpecificHeatCapacity(
    const Real* const density,
    const Real* const pressure,
    const std::vector<const Real*>& mass_fractions) const
{
    TBOX_ERROR(d_object_name
        << ": EquationOfStateMixingRulesStiffenedGas::getIsobaricSpecificHeatCapacity()\n"
        << "Method getIsobaricSpecificHeatCapacity() for mixture with isothermal and isobaric equilibrium assumptions"
        << " is not yet implemented."
        << std::endl);
    
    return Real(0);
}


/*
 * Compute the isobaric specific heat capacity of mixture with isothermal and isobaric equilibrium
 * assumptions.
 */
void
EquationOfStateMixingRulesStiffenedGas::computeIsobaricSpecificHeatCapacity(
    HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_isobaric_specific_heat_capacity,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_mass_fractions,
    const hier::Box& domain) const
{
    TBOX_ERROR(d_object_name
        << ": EquationOfStateMixingRulesStiffenedGas::computeIsobaricSpecificHeatCapacity()\n"
        << "Method computeIsobaricSpecificHeatCapacity() for mixture with isothermal and isobaric equilibrium"
        << " assumptions is not yet implemented."
        << std::endl);
}


/*
 * Compute the isobaric specific heat capacity of mixture with isothermal and isobaric equilibrium
 * assumptions.
 */
void
EquationOfStateMixingRulesStiffenedGas::computeIsobaricSpecificHeatCapacity(
    HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_isobaric_specific_heat_capacity,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_mass_fractions,
    int side_normal,
    const hier::Box& domain) const
{
    TBOX_ERROR(d_object_name
        << ": EquationOfStateMixingRulesStiffenedGas::computeIsobaricSpecificHeatCapacity()\n"
        << "Method computeIsobaricSpecificHeatCapacity() for mixture with isothermal and isobaric equilibrium"
        << " assumptions is not yet implemented."
        << std::endl);
}


/*
 * Compute the Gruneisen parameter of the mixture with isothermal and isobaric equilibrium assumptions
 * (partial derivative of pressure w.r.t. specific internal energy under constant partial densities
 * divided by mixture density).
 */
Real
EquationOfStateMixingRulesStiffenedGas::getGruneisenParameter(
    const Real* const density,
    const Real* const pressure,
    const std::vector<const Real*>& mass_fractions) const
{
    TBOX_ERROR(d_object_name
        << ": EquationOfStateMixingRulesStiffenedGas::getGruneisenParameter()\n"
        << "Method getGruneisenParameter() for mixture with isothermal and isobaric equilibrium assumptions"
        << " is not yet implemented."
        << std::endl);
    
    return Real(0);
}


/*
 * Compute the Gruneisen parameter of the mixture with isothermal and isobaric equilibrium assumptions
 * (partial derivative of pressure w.r.t. specific internal energy under constant partial densities
 * divided by mixture density).
 */
void
EquationOfStateMixingRulesStiffenedGas::computeGruneisenParameter(
    HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_gruneisen_parameter,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_mass_fractions,
    const hier::Box& domain) const
{
    TBOX_ERROR(d_object_name
        << ": EquationOfStateMixingRulesStiffenedGas::computeGruneisenParameter()\n"
        << "Method computeGruneisenParameter() for mixture with isothermal and isobaric equilibrium assumptions"
        << " is not yet implemented."
        << std::endl);
}


/*
 * Compute the Gruneisen parameter of the mixture with isothermal and isobaric equilibrium assumptions
 * (partial derivative of pressure w.r.t. specific internal energy under constant partial densities
 * divided by mixture density).
 */
void
EquationOfStateMixingRulesStiffenedGas::computeGruneisenParameter(
    HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_gruneisen_parameter,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_mass_fractions,
    int side_normal,
    const hier::Box& domain) const
{
    TBOX_ERROR(d_object_name
        << ": EquationOfStateMixingRulesStiffenedGas::computeGruneisenParameter()\n"
        << "Method computeGruneisenParameter() for mixture with isothermal and isobaric equilibrium assumptions"
        << " is not yet implemented."
        << std::endl);
}


/*
 * Compute the Gruneisen parameter of the mixture with isobaric equilibrium assumption
 * (partial derivative of pressure w.r.t. specific internal energy under constant partial densities
 * and volume fractions divided by mixture density).
 */
Real
EquationOfStateMixingRulesStiffenedGas::getGruneisenParameter(
    const Real* const density,
    const Real* const pressure,
    const std::vector<const Real*>& mass_fractions,
    const std::vector<const Real*>& volume_fractions) const
{
    NULL_USE(mass_fractions);
    
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOBARIC);
    TBOX_ASSERT((static_cast<int>(volume_fractions.size()) == d_num_species) ||
                (static_cast<int>(volume_fractions.size()) == d_num_species - 1));
#endif
    
    // Get the mixture thermodynamic properties.
    std::vector<Real> mixture_thermo_properties;
    std::vector<Real*> mixture_thermo_properties_ptr;
    std::vector<const Real*> mixture_thermo_properties_const_ptr;
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    mixture_thermo_properties.resize(num_thermo_properties);
    mixture_thermo_properties_ptr.reserve(num_thermo_properties);
    mixture_thermo_properties_const_ptr.reserve(num_thermo_properties);
    
    for (int ti = 0; ti < num_thermo_properties; ti++)
    {
        mixture_thermo_properties_ptr.push_back(&mixture_thermo_properties[ti]);
        mixture_thermo_properties_const_ptr.push_back(&mixture_thermo_properties[ti]);
    }
    
    getMixtureThermodynamicPropertiesWithVolumeFractions(
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
EquationOfStateMixingRulesStiffenedGas::computeGruneisenParameter(
    HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_gruneisen_parameter,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_mass_fractions,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_volume_fractions,
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
    HAMERS_SHARED_PTR<pdat::CellData<Real> > data_mixture_thermo_properties;
    
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
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::CellData<Real> >(
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
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::CellData<Real> >(
            domain, num_thermo_properties, hier::IntVector::getZero(d_dim));
    }
    
    computeMixtureThermodynamicPropertiesWithVolumeFractions(
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
EquationOfStateMixingRulesStiffenedGas::computeGruneisenParameter(
    HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_gruneisen_parameter,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_mass_fractions,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_volume_fractions,
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
    
    // Get the box that covers the interior of patch.
    HAMERS_SHARED_PTR<pdat::SideData<Real> > data_mixture_thermo_properties;
    
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
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::SideData<Real> >(
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
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::SideData<Real> >(
            domain, num_thermo_properties, hier::IntVector::getZero(d_dim), direction);
    }
    
    computeMixtureThermodynamicPropertiesWithVolumeFractions(
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
std::vector<Real>
EquationOfStateMixingRulesStiffenedGas::getPressureDerivativeWithPartialDensities(
    const Real* const density,
    const Real* const pressure,
    const std::vector<const Real*>& mass_fractions) const
{
    TBOX_ERROR(d_object_name
        << ": EquationOfStateMixingRulesStiffenedGas::getPressureDerivativeWithPartialDensities()\n"
        << "Method getPressureDerivativeWithPartialDensities() for mixture with isothermal and isobaric equilibrium"
        << " assumptions is not yet implemented."
        << std::endl);
    
    std::vector<Real> Psi;
    Psi.reserve(d_num_species);
    for (int si = 0; si < d_num_species; si++)
    {
        Psi.push_back(Real(0));
    }
    
    return Psi;
}


/*
 * Compute the mixture partial derivative of pressure w.r.t. partial densities under constant specific
 * internal energy with isothermal and isobaric equilibrium assumptions.
 */
void
EquationOfStateMixingRulesStiffenedGas::computePressureDerivativeWithPartialDensities(
    HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_partial_pressure_partial_partial_densities,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_mass_fractions,
    const hier::Box& domain) const
{
    TBOX_ERROR(d_object_name
        << ": EquationOfStateMixingRulesStiffenedGas::computePressureDerivativeWithPartialDensities()\n"
        << "Method computePressureDerivativeWithPartialDensities() for mixture with isothermal and isobaric"
        << " equilibrium assumptions is not yet implemented."
        << std::endl);
}


/*
 * Compute the mixture partial derivative of pressure w.r.t. partial densities under constant specific
 * internal energy with isothermal and isobaric equilibrium assumptions.
 */
void
EquationOfStateMixingRulesStiffenedGas::computePressureDerivativeWithPartialDensities(
    HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_partial_pressure_partial_partial_densities,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_mass_fractions,
    int side_normal,
    const hier::Box& domain) const
{
    TBOX_ERROR(d_object_name
        << ": EquationOfStateMixingRulesStiffenedGas::computePressureDerivativeWithPartialDensities()\n"
        << "Method computePressureDerivativeWithPartialDensities() for mixture with isothermal and isobaric"
        << " equilibrium assumptions is not yet implemented."
        << std::endl);
}


/*
 * Compute the mixture partial derivative of pressure w.r.t. partial densities under constant specific
 * internal energy and volume fractions with isobaric equilibrium assumption.
 */
std::vector<Real>
EquationOfStateMixingRulesStiffenedGas::getPressureDerivativeWithPartialDensities(
    const Real* const density,
    const Real* const pressure,
    const std::vector<const Real*>& mass_fractions,
    const std::vector<const Real*>& volume_fractions) const
{
    NULL_USE(mass_fractions);
    
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOBARIC);
    TBOX_ASSERT((static_cast<int>(volume_fractions.size()) == d_num_species) ||
                (static_cast<int>(volume_fractions.size()) == d_num_species - 1));
#endif
    
    // Get the mixture thermodynamic properties.
    std::vector<Real> mixture_thermo_properties;
    std::vector<Real*> mixture_thermo_properties_ptr;
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    mixture_thermo_properties.resize(num_thermo_properties);
    mixture_thermo_properties_ptr.reserve(num_thermo_properties);
    
    for (int ti = 0; ti < num_thermo_properties; ti++)
    {
        mixture_thermo_properties_ptr.push_back(&mixture_thermo_properties[ti]);
    }
    
    getMixtureThermodynamicPropertiesWithVolumeFractions(
        mixture_thermo_properties_ptr,
        volume_fractions);
    
    const Real& gamma = mixture_thermo_properties[0];
    const Real& p_inf = mixture_thermo_properties[1];
    
    const Real& rho = *density;
    const Real& p   = *pressure;
    
    std::vector<Real> Psi;
    Psi.reserve(d_num_species);
    for (int si = 0; si < d_num_species; si++)
    {
        Real Psi_i = (p + gamma*p_inf)/rho;
        Psi.push_back(Psi_i);
    }
    
    return Psi;
}


/*
 * Compute the mixture partial derivative of pressure w.r.t. partial densities under constant specific
 * internal energy and volume fractions with isobaric equilibrium assumption.
 */
void
EquationOfStateMixingRulesStiffenedGas::computePressureDerivativeWithPartialDensities(
    HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_partial_pressure_partial_partial_densities,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_mass_fractions,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_volume_fractions,
    const hier::Box& domain) const
{
    NULL_USE(data_mass_fractions);
    
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOBARIC);
    
    TBOX_ASSERT(data_partial_pressure_partial_partial_densities);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_volume_fractions);
    
    TBOX_ASSERT(data_partial_pressure_partial_partial_densities->getDepth() == d_num_species);
    
    TBOX_ASSERT((data_volume_fractions->getDepth() == d_num_species) ||
                (data_volume_fractions->getDepth() == d_num_species - 1));
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
    
    hier::IntVector ghostcell_dims_min(d_dim);
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets and allocate memory for the mixture thermodyanmic properties.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_partial_pressure_partial_partial_densities(d_dim);
    hier::IntVector offset_density(d_dim);
    hier::IntVector offset_pressure(d_dim);
    hier::IntVector offset_min(d_dim);
    
    HAMERS_SHARED_PTR<pdat::CellData<Real> > data_mixture_thermo_properties;
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_partial_pressure_partial_partial_densities =
            data_partial_pressure_partial_partial_densities->getGhostCellWidth();
        
        const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        const hier::IntVector num_ghosts_volume_fractions = data_volume_fractions->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_partial_pressure_partial_partial_densities->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_pressure->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_volume_fractions->getBox().isSpatiallyEqual(interior_box));
#endif
        
        /*
         * Get the minimum number of ghost cells and the dimension of the ghost cell box for
         * mixture thermodynamic properties.
         */
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_partial_pressure_partial_partial_densities;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_volume_fractions, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        ghostcell_dims_min = ghost_box.numberCells();
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_partial_pressure_partial_partial_densities = num_ghosts_partial_pressure_partial_partial_densities;
        offset_density = num_ghosts_density;
        offset_pressure = num_ghosts_pressure;
        offset_min = num_ghosts_min;
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::CellData<Real> >(
            interior_box, num_thermo_properties, num_ghosts_min);
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
        TBOX_ASSERT(data_partial_pressure_partial_partial_densities->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_volume_fractions->getGhostBox().contains(domain));
#endif
        
        ghostcell_dims_min = domain.numberCells();
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_partial_pressure_partial_partial_densities =
            domain.lower() - ghost_box_partial_pressure_partial_partial_densities.lower();
        
        offset_density = domain.lower() - ghost_box_density.lower();
        offset_pressure = domain.lower() - ghost_box_pressure.lower();
        offset_min = hier::IntVector::getZero(d_dim);
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::CellData<Real> >(
            domain, num_thermo_properties, hier::IntVector::getZero(d_dim));
    }
    
    computeMixtureThermodynamicPropertiesWithVolumeFractions(
        data_mixture_thermo_properties,
        data_volume_fractions,
        domain);
    
    /*
     * Get the pointers to the cell data.
     */
    
    Real* const rho = data_density->getPointer(0);
    Real* const p   = data_pressure->getPointer(0);
    Real* const gamma = data_mixture_thermo_properties->getPointer(0);
    Real* const p_inf = data_mixture_thermo_properties->getPointer(1);
    
    /*
     * Get the partial derivative.
     */
    
    std::vector<Real*> Psi;
    Psi.reserve(d_num_species);
    for (int si = 0; si < d_num_species; si++)
    {
        Psi.push_back(data_partial_pressure_partial_partial_densities->getPointer(si));
    }
    
    computePressureDerivativeWithPartialDensities(
        Psi,
        rho,
        p,
        gamma,
        p_inf,
        offset_partial_pressure_partial_partial_densities,
        offset_density,
        offset_pressure,
        offset_min,
        ghostcell_dims_partial_pressure_partial_partial_densities,
        ghostcell_dims_density,
        ghostcell_dims_pressure,
        ghostcell_dims_min,
        domain_lo,
        domain_dims);
}


/*
 * Compute the mixture partial derivative of pressure w.r.t. partial densities under constant specific
 * internal energy and volume fractions with isobaric equilibrium assumption.
 */
void
EquationOfStateMixingRulesStiffenedGas::computePressureDerivativeWithPartialDensities(
    HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_partial_pressure_partial_partial_densities,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_mass_fractions,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_volume_fractions,
    int side_normal,
    const hier::Box& domain) const
{
    NULL_USE(data_mass_fractions);
    
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOBARIC);
    
    TBOX_ASSERT(data_partial_pressure_partial_partial_densities);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_volume_fractions);
    
    TBOX_ASSERT(data_partial_pressure_partial_partial_densities->getDepth() == d_num_species);
    
    TBOX_ASSERT((data_volume_fractions->getDepth() == d_num_species) ||
                (data_volume_fractions->getDepth() == d_num_species - 1));
#endif
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(side_normal < d_dim.getValue());
    
    TBOX_ASSERT(data_partial_pressure_partial_partial_densities->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_density->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_pressure->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_volume_fractions->getDirectionVector()[side_normal] > 0);
#endif
    
    hier::IntVector direction = hier::IntVector::getZero(d_dim);
    direction[side_normal] = 1;
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_partial_pressure_partial_partial_densities =
        data_partial_pressure_partial_partial_densities->getGhostBox();
    hier::IntVector ghostcell_dims_partial_pressure_partial_partial_densities =
        ghost_box_partial_pressure_partial_partial_densities.numberCells();
    
    const hier::Box ghost_box_density = data_density->getGhostBox();
    hier::IntVector ghostcell_dims_density = ghost_box_density.numberCells();
    
    const hier::Box ghost_box_pressure = data_pressure->getGhostBox();
    hier::IntVector ghostcell_dims_pressure = ghost_box_pressure.numberCells();
    
    hier::IntVector ghostcell_dims_min(d_dim);
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets and allocate memory for the mixture thermodyanmic properties.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_partial_pressure_partial_partial_densities(d_dim);
    hier::IntVector offset_density(d_dim);
    hier::IntVector offset_pressure(d_dim);
    hier::IntVector offset_min(d_dim);
    
    HAMERS_SHARED_PTR<pdat::SideData<Real> > data_mixture_thermo_properties;
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_partial_pressure_partial_partial_densities =
            data_partial_pressure_partial_partial_densities->getGhostCellWidth();
        
        const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        const hier::IntVector num_ghosts_volume_fractions = data_volume_fractions->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_partial_pressure_partial_partial_densities->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_pressure->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_volume_fractions->getBox().isSpatiallyEqual(interior_box));
#endif
        
        /*
         * Get the minimum number of ghost cells and the dimension of the ghost cell box for
         * mixture thermodynamic properties.
         */
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_partial_pressure_partial_partial_densities;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_volume_fractions, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        ghostcell_dims_min = ghost_box.numberCells();
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_partial_pressure_partial_partial_densities = num_ghosts_partial_pressure_partial_partial_densities;
        offset_density = num_ghosts_density;
        offset_pressure = num_ghosts_pressure;
        offset_min = num_ghosts_min;
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::SideData<Real> >(
            interior_box, num_thermo_properties, num_ghosts_min, direction);
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
        TBOX_ASSERT(data_partial_pressure_partial_partial_densities->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_volume_fractions->getGhostBox().contains(domain));
#endif
        
        ghostcell_dims_min = domain.numberCells();
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_partial_pressure_partial_partial_densities =
            domain.lower() - ghost_box_partial_pressure_partial_partial_densities.lower();
        
        offset_density = domain.lower() - ghost_box_density.lower();
        offset_pressure = domain.lower() - ghost_box_pressure.lower();
        offset_min = hier::IntVector::getZero(d_dim);
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::SideData<Real> >(
            domain, num_thermo_properties, hier::IntVector::getZero(d_dim), direction);
    }
    
    ghostcell_dims_partial_pressure_partial_partial_densities[side_normal]++;
    ghostcell_dims_density[side_normal]++;
    ghostcell_dims_pressure[side_normal]++;
    ghostcell_dims_min[side_normal]++;
    domain_dims[side_normal]++;
    
    computeMixtureThermodynamicPropertiesWithVolumeFractions(
        data_mixture_thermo_properties,
        data_volume_fractions,
        side_normal,
        domain);
    
    /*
     * Get the pointers to the cell data.
     */
    
    Real* const rho = data_density->getPointer(side_normal, 0);
    Real* const p   = data_pressure->getPointer(side_normal, 0);
    Real* const gamma = data_mixture_thermo_properties->getPointer(side_normal, 0);
    Real* const p_inf = data_mixture_thermo_properties->getPointer(side_normal, 1);
    
    /*
     * Get the partial derivative.
     */
    
    std::vector<Real*> Psi;
    Psi.reserve(d_num_species);
    for (int si = 0; si < d_num_species; si++)
    {
        Psi.push_back(data_partial_pressure_partial_partial_densities->getPointer(side_normal, si));
    }
    
    computePressureDerivativeWithPartialDensities(
        Psi,
        rho,
        p,
        gamma,
        p_inf,
        offset_partial_pressure_partial_partial_densities,
        offset_density,
        offset_pressure,
        offset_min,
        ghostcell_dims_partial_pressure_partial_partial_densities,
        ghostcell_dims_density,
        ghostcell_dims_pressure,
        ghostcell_dims_min,
        domain_lo,
        domain_dims);
}


/*
 * Compute the mixture partial derivative of pressure w.r.t. volume fractions under constant specific
 * internal energy and partial densities with isobaric equilibrium assumption.
 */
std::vector<Real>
EquationOfStateMixingRulesStiffenedGas::getPressureDerivativeWithVolumeFractions(
    const Real* const density,
    const Real* const pressure,
    const std::vector<const Real*>& mass_fractions,
    const std::vector<const Real*>& volume_fractions) const
{
    NULL_USE(mass_fractions);
    
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOBARIC);
    TBOX_ASSERT((static_cast<int>(volume_fractions.size()) == d_num_species) ||
                (static_cast<int>(volume_fractions.size()) == d_num_species - 1));
#endif
    
    // Get the mixture thermodynamic properties.
    std::vector<Real> mixture_thermo_properties;
    std::vector<Real*> mixture_thermo_properties_ptr;
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    mixture_thermo_properties.resize(num_thermo_properties);
    mixture_thermo_properties_ptr.reserve(num_thermo_properties);
    
    for (int ti = 0; ti < num_thermo_properties; ti++)
    {
        mixture_thermo_properties_ptr.push_back(&mixture_thermo_properties[ti]);
    }
    
    getMixtureThermodynamicPropertiesWithVolumeFractions(
        mixture_thermo_properties_ptr,
        volume_fractions);
    
    const Real& gamma = mixture_thermo_properties[0];
    
    const Real& p   = *pressure;
    
    const Real tmp = (p + d_species_gamma[d_num_species - 1]*d_species_p_inf[d_num_species - 1])/
        (d_species_gamma[d_num_species - 1] - Real(1));
    
    std::vector<Real> M;
    M.reserve(d_num_species - 1);
    for (int si = 0; si < d_num_species - 1; si++)
    {
        Real M_i = (tmp - (p + d_species_gamma[si]*d_species_p_inf[si])/(d_species_gamma[si] - Real(1)))*
            (gamma - Real(1));
        M.push_back(M_i);
    }
    
    return M;
}


/*
 * Compute the mixture partial derivative of pressure w.r.t. volume fractions under constant specific
 * internal energy and partial densities with isobaric equilibrium assumption.
 */
void
EquationOfStateMixingRulesStiffenedGas::computePressureDerivativeWithVolumeFractions(
    HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_partial_pressure_partial_volume_fractions,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_mass_fractions,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_volume_fractions,
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
    
    HAMERS_SHARED_PTR<pdat::CellData<Real> > data_mixture_thermo_properties;
    
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
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::CellData<Real> >(
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
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::CellData<Real> >(
            domain, num_thermo_properties, hier::IntVector::getZero(d_dim));
    }
    
    // Compute the mixture thermodyanmic properties.
    computeMixtureThermodynamicPropertiesWithVolumeFractions(
        data_mixture_thermo_properties,
        data_volume_fractions,
        domain);
    
    /*
     * Get the pointers to the cell data.
     */
    
    Real* const p     = data_pressure->getPointer(0);
    Real* const gamma = data_mixture_thermo_properties->getPointer(0);
    
    /*
     * Get the partial derivative.
     */
    
    std::vector<Real*> M;
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
EquationOfStateMixingRulesStiffenedGas::computePressureDerivativeWithVolumeFractions(
    HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_partial_pressure_partial_volume_fractions,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_mass_fractions,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_volume_fractions,
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
    
    HAMERS_SHARED_PTR<pdat::SideData<Real> > data_mixture_thermo_properties;
    
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
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::SideData<Real> >(
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
        
        data_mixture_thermo_properties = HAMERS_MAKE_SHARED<pdat::SideData<Real> >(
            domain, num_thermo_properties, hier::IntVector::getZero(d_dim), direction);
    }
    
    ghostcell_dims_partial_pressure_partial_volume_fractions[side_normal]++;
    ghostcell_dims_pressure[side_normal]++;
    ghostcell_dims_min[side_normal]++;
    domain_dims[side_normal]++;
    
    // Compute the mixture thermodyanmic properties.
    computeMixtureThermodynamicPropertiesWithVolumeFractions(
        data_mixture_thermo_properties,
        data_volume_fractions,
        side_normal,
        domain);
    
    /*
     * Get the pointers to the cell data.
     */
    
    Real* const p     = data_pressure->getPointer(side_normal, 0);
    Real* const gamma = data_mixture_thermo_properties->getPointer(side_normal, 0);
    
    /*
     * Get the partial derivative.
     */
    
    std::vector<Real*> M;
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
Real
EquationOfStateMixingRulesStiffenedGas::getMixtureDensity(
    const Real* const pressure,
    const Real* const temperature,
    const std::vector<const Real*>& mass_fractions) const
{
    TBOX_ERROR(d_object_name
        << ": EquationOfStateMixingRulesStiffenedGas::getMixtureDensity()\n"
        << "Method getMixtureDensity() for mixture with isothermal and isobaric equilibrium assumptions"
        << " is not yet implemented."
        << std::endl);
    
    return Real(0);
}


/*
 * Compute the density of mixture with isothermal and isobaric equilibrium assumptions.
 */
void
EquationOfStateMixingRulesStiffenedGas::computeMixtureDensity(
    HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_mixture_density,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_temperature,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_mass_fractions,
    const hier::Box& domain) const
{
    TBOX_ERROR(d_object_name
        << ": EquationOfStateMixingRulesStiffenedGas::computeMixtureDensity()\n"
        << "Method computeMixtureDensity() for mixture with isothermal and isobaric equilibrium assumptions"
        << " is not yet implemented."
        << std::endl);
}


/*
 * Compute the density of mixture with isothermal and isobaric equilibrium assumptions.
 */
void
EquationOfStateMixingRulesStiffenedGas::computeMixtureDensity(
    HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_mixture_density,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_temperature,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_mass_fractions,
    int side_normal,
    const hier::Box& domain) const
{
    TBOX_ERROR(d_object_name
        << ": EquationOfStateMixingRulesStiffenedGas::computeMixtureDensity()\n"
        << "Method computeMixtureDensity() for mixture with isothermal and isobaric equilibrium assumptions"
        << " is not yet implemented."
        << std::endl);
}


/*
 * Get the thermodynamic properties of a species.
 */
void
EquationOfStateMixingRulesStiffenedGas::getSpeciesThermodynamicProperties(
    std::vector<Real*>& species_thermo_properties,
    const int species_index) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(species_thermo_properties.size()) == 4);
    TBOX_ASSERT(species_index >= 0);
    TBOX_ASSERT(species_index < d_num_species);
#endif
    
    // Get references to gamma and p_inf of species.
    Real& gamma = *(species_thermo_properties[0]);
    Real& p_inf = *(species_thermo_properties[1]);
    
    // Get references to c_p and c_v of species.
    Real& c_p = *(species_thermo_properties[2]);
    Real& c_v = *(species_thermo_properties[3]);
    
    gamma = d_species_gamma[species_index];
    p_inf = d_species_p_inf[species_index];
    
    c_p = d_species_c_p[species_index];
    c_v = d_species_c_v[species_index];
}


/*
 * Get the number of thermodynamic properties of the mixture.
 */
int
EquationOfStateMixingRulesStiffenedGas::getNumberOfMixtureThermodynamicProperties() const
{
    int num_of_thermo_properties = 0;
    
    switch (d_mixing_closure_model)
    {
        case MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC:
        {
            num_of_thermo_properties = 0;
            
            break;
        }
        case MIXING_CLOSURE_MODEL::ISOBARIC:
        {
            num_of_thermo_properties = 2;
            
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
 * Compute the thermodynamic properties of the mixture with volume fractions.
 */
void
EquationOfStateMixingRulesStiffenedGas::getMixtureThermodynamicPropertiesWithVolumeFractions(
    std::vector<Real*>& mixture_thermo_properties,
    const std::vector<const Real*>& volume_fractions) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(mixture_thermo_properties.size()) == 2);
#endif
    
    // Get reference to gamma and p_inf of mixture.
    Real& gamma = *(mixture_thermo_properties[0]);
    Real& p_inf = *(mixture_thermo_properties[1]);
    
    // Compute gamma and p_inf of mixture.
    
    Real tmp_1 = Real(0);
    Real tmp_2 = Real(0);
    
    if (static_cast<int>(volume_fractions.size()) == d_num_species)
    {
        for (int si = 0; si < d_num_species; si++)
        {
            tmp_1 += *(volume_fractions[si])/(d_species_gamma[si] - Real(1));
            
            tmp_2 += (*(volume_fractions[si])*d_species_gamma[si]*d_species_p_inf[si])/
                (d_species_gamma[si] - Real(1));
        }
    }
    else if (static_cast<int>(volume_fractions.size()) == d_num_species - 1)
    {
        Real Z_last = Real(1);
        
        for (int si = 0; si < d_num_species - 1; si++)
        {
            tmp_1 += *(volume_fractions[si])/(d_species_gamma[si] - Real(1));
            
            tmp_2 += (*(volume_fractions[si])*d_species_gamma[si]*d_species_p_inf[si])/
                (d_species_gamma[si] - Real(1));
            
            // Compute the volume fraction of the last species.
            Z_last -= *(volume_fractions[si]);
        }
        
        // Add the contributions from the last species.
        tmp_1 += Z_last/(d_species_gamma.back() - Real(1));
        tmp_2 += (Z_last*d_species_gamma.back()*d_species_p_inf.back())/(d_species_gamma.back() - Real(1));
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Number of volume fractions provided is not"
            << " equal to the total number of species or (total number of species - 1)."
            << std::endl);
    }
    
    gamma = Real(1)/tmp_1 + Real(1);
    p_inf = (gamma - Real(1))/gamma*tmp_2;
}


/*
 * Compute the thermodynamic properties of the mixture with volume fractions.
 */
void
EquationOfStateMixingRulesStiffenedGas::computeMixtureThermodynamicPropertiesWithVolumeFractions(
    HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_mixture_thermo_properties,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_volume_fractions,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(data_mixture_thermo_properties);
    TBOX_ASSERT(data_mixture_thermo_properties->getDepth() == 2);
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
     * Get the pointers to the cell data of mixture thermodynamic properties.
     */
    
    Real* const gamma = data_mixture_thermo_properties->getPointer(0);
    Real* const p_inf = data_mixture_thermo_properties->getPointer(1);
    
    /*
     * Fill zeros for gamma and p_inf.
     */
    
    if (domain.empty())
    {
        data_mixture_thermo_properties->fill(Real(0), 0);
        data_mixture_thermo_properties->fill(Real(0), 1);
    }
    else
    {
        data_mixture_thermo_properties->fill(Real(0), domain, 0);
        data_mixture_thermo_properties->fill(Real(0), domain, 1);
    }
    
    if (data_volume_fractions->getDepth() == d_num_species)
    {
        /*
         * Get the pointers to the cell data of volume fractions.
         */
        
        std::vector<const Real*> Z;
        Z.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            Z.push_back(data_volume_fractions->getPointer(si));
        }
        
        computeMixtureThermodynamicPropertiesWithVolumeFractions(
            gamma,
            p_inf,
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
        HAMERS_SHARED_PTR<pdat::CellData<Real> > data_volume_fractions_last;
        
        hier::IntVector offset_volume_fractions_last(d_dim);
        
        if (domain.empty())
        {
            const hier::Box interior_box = data_volume_fractions->getBox();
            
            offset_volume_fractions_last = data_volume_fractions->getGhostCellWidth();
            
            data_volume_fractions_last = HAMERS_MAKE_SHARED<pdat::CellData<Real> >(
                interior_box, 1, data_volume_fractions->getGhostCellWidth());
            
            data_volume_fractions_last->fillAll(Real(1));
        }
        else
        {
            offset_volume_fractions_last = hier::IntVector::getZero(d_dim);
            
            data_volume_fractions_last = HAMERS_MAKE_SHARED<pdat::CellData<Real> >(
                domain, 1, hier::IntVector::getZero(d_dim));
            
            data_volume_fractions_last->fillAll(Real(1), domain);
        }
        
        /*
         * Get the pointers to the cell data of volume fractions and the dimensions of the ghost cell box of
         * last volume fraction.
         */
        
        std::vector<const Real*> Z;
        Z.reserve(d_num_species - 1);
        for (int si = 0; si < d_num_species - 1; si++)
        {
            Z.push_back(data_volume_fractions->getPointer(si));
        }
        
        Real* const Z_last = data_volume_fractions_last->getPointer(0);
        
        const hier::Box ghost_box_volume_fractions_last = data_volume_fractions_last->getGhostBox();
        const hier::IntVector ghostcell_dims_volume_fractions_last = ghost_box_volume_fractions_last.numberCells();
        
        computeMixtureThermodynamicPropertiesWithVolumeFractions(
            gamma,
            p_inf,
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
EquationOfStateMixingRulesStiffenedGas::computeMixtureThermodynamicPropertiesWithVolumeFractions(
    HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_mixture_thermo_properties,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_volume_fractions,
    int side_normal,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(data_mixture_thermo_properties);
    TBOX_ASSERT(data_mixture_thermo_properties->getDepth() == 2);
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
     * Get the pointers to the cell data of mixture thermodynamic properties.
     */
    
    Real* const gamma = data_mixture_thermo_properties->getPointer(side_normal, 0);
    Real* const p_inf = data_mixture_thermo_properties->getPointer(side_normal, 1);
    
    /*
     * Fill zeros for gamma and p_inf.
     */
    
    if (domain.empty())
    {
        data_mixture_thermo_properties->fill(Real(0), 0);
        data_mixture_thermo_properties->fill(Real(0), 1);
    }
    else
    {
        data_mixture_thermo_properties->fill(Real(0), domain, 0);
        data_mixture_thermo_properties->fill(Real(0), domain, 1);
    }
    
    if (data_volume_fractions->getDepth() == d_num_species)
    {
        /*
         * Get the pointers to the cell data of volume fractions.
         */
        
        std::vector<const Real*> Z;
        Z.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            Z.push_back(data_volume_fractions->getPointer(side_normal, si));
        }
        
        computeMixtureThermodynamicPropertiesWithVolumeFractions(
            gamma,
            p_inf,
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
        
        HAMERS_SHARED_PTR<pdat::SideData<Real> > data_volume_fractions_last;
        
        hier::IntVector offset_volume_fractions_last(d_dim);
        
        if (domain.empty())
        {
            const hier::Box interior_box = data_volume_fractions->getBox();
            
            offset_volume_fractions_last = data_volume_fractions->getGhostCellWidth();
            
            data_volume_fractions_last = HAMERS_MAKE_SHARED<pdat::SideData<Real> >(
                interior_box, 1, data_volume_fractions->getGhostCellWidth(), direction);
        
            data_volume_fractions_last->fillAll(Real(1));
        }
        else
        {
            offset_volume_fractions_last = hier::IntVector::getZero(d_dim);
            
            data_volume_fractions_last = HAMERS_MAKE_SHARED<pdat::SideData<Real> >(
                domain, 1, hier::IntVector::getZero(d_dim), direction);
            
            data_volume_fractions_last->fillAll(Real(1), domain);
        }
        
        /*
         * Get the pointers to the cell data of volume fractions and the dimensions of the ghost cell box of
         * last volume fraction.
         */
        
        std::vector<const Real*> Z;
        Z.reserve(d_num_species - 1);
        for (int si = 0; si < d_num_species - 1; si++)
        {
            Z.push_back(data_volume_fractions->getPointer(side_normal, si));
        }
        
        Real* const Z_last = data_volume_fractions_last->getPointer(side_normal, 0);
        
        const hier::Box ghost_box_volume_fractions_last = data_volume_fractions_last->getGhostBox();
        hier::IntVector ghostcell_dims_volume_fractions_last = ghost_box_volume_fractions_last.numberCells();
        ghostcell_dims_volume_fractions_last[side_normal]++;
        
        computeMixtureThermodynamicPropertiesWithVolumeFractions(
            gamma,
            p_inf,
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
 * Compute the mixture partial derivative of pressure w.r.t. partial densities under constant specific
 * internal energy and volume fractions with isobaric equilibrium assumption.
 */
void
EquationOfStateMixingRulesStiffenedGas::computePressureDerivativeWithPartialDensities(
    std::vector<Real*>& Psi,
    const Real* const rho,
    const Real* const p,
    const Real* const gamma,
    const Real* const p_inf,
    const hier::IntVector& offset_partial_pressure_partial_partial_densities,
    const hier::IntVector& offset_density,
    const hier::IntVector& offset_pressure,
    const hier::IntVector& offset_mixture_thermo_properties,
    const hier::IntVector& ghostcell_dims_partial_pressure_partial_partial_densities,
    const hier::IntVector& ghostcell_dims_density,
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
        
        const int offset_0_partial_pressure_partial_partial_densities =
            offset_partial_pressure_partial_partial_densities[0];
        const int offset_0_density = offset_density[0];
        const int offset_0_pressure = offset_pressure[0];
        const int offset_0_mixture_thermo_properties = offset_mixture_thermo_properties[0];
        
        // Compute Psi.
        for (int si = 0; si < d_num_species; si++)
        {
            Real* Psi_i = Psi[si];
            
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_partial_pressure_partial_partial_densities =
                    i + offset_0_partial_pressure_partial_partial_densities;
                
                const int idx_density = i + offset_0_density;
                
                const int idx_pressure = i + offset_0_pressure;
                
                const int idx_mixture_thermo_properties = i + offset_0_mixture_thermo_properties;
                
                Psi_i[idx_partial_pressure_partial_partial_densities] = (p[idx_pressure] +
                    gamma[idx_mixture_thermo_properties]*p_inf[idx_mixture_thermo_properties])/
                    rho[idx_density];
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
        
        const int offset_0_mixture_thermo_properties = offset_mixture_thermo_properties[0];
        const int offset_1_mixture_thermo_properties = offset_mixture_thermo_properties[1];
        const int ghostcell_dim_0_mixture_thermo_properties = ghostcell_dims_mixture_thermo_properties[0];
        
        // Compute Psi.
        for (int si = 0; si < d_num_species; si++)
        {
            Real* Psi_i = Psi[si];
            
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
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
                    
                    const int idx_mixture_thermo_properties = (i + offset_0_mixture_thermo_properties) +
                        (j + offset_1_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties;
                    
                    Psi_i[idx_partial_pressure_partial_partial_densities] = (p[idx_pressure] +
                        gamma[idx_mixture_thermo_properties]*p_inf[idx_mixture_thermo_properties])/
                        rho[idx_density];
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
        
        const int offset_0_mixture_thermo_properties = offset_mixture_thermo_properties[0];
        const int offset_1_mixture_thermo_properties = offset_mixture_thermo_properties[1];
        const int offset_2_mixture_thermo_properties = offset_mixture_thermo_properties[2];
        const int ghostcell_dim_0_mixture_thermo_properties = ghostcell_dims_mixture_thermo_properties[0];
        const int ghostcell_dim_1_mixture_thermo_properties = ghostcell_dims_mixture_thermo_properties[1];
        
        // Compute Psi.
        for (int si = 0; si < d_num_species; si++)
        {
            Real* Psi_i = Psi[si];
            
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
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
                        
                        const int idx_mixture_thermo_properties = (i + offset_0_mixture_thermo_properties) +
                            (j + offset_1_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties +
                            (k + offset_2_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties*
                                ghostcell_dim_1_mixture_thermo_properties;
                        
                        Psi_i[idx_partial_pressure_partial_partial_densities] = (p[idx_pressure] +
                            gamma[idx_mixture_thermo_properties]*p_inf[idx_mixture_thermo_properties])/
                            rho[idx_density];
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
EquationOfStateMixingRulesStiffenedGas::computePressureDerivativeWithVolumeFractions(
    std::vector<Real*>& M,
    const Real* const p,
    const Real* const gamma,
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
            Real* M_i = M[si];
            
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_partial_pressure_partial_volume_fractions =
                    i + offset_0_partial_pressure_partial_volume_fractions;
                
                const int idx_pressure = i + offset_0_pressure;
                
                const int idx_mixture_thermo_properties = i + offset_0_mixture_thermo_properties;
                
                M_i[idx_partial_pressure_partial_volume_fractions] =
                    (gamma[idx_mixture_thermo_properties] - Real(1))*
                    ((p[idx_pressure] + d_species_gamma[d_num_species - 1]*d_species_p_inf[d_num_species - 1])/
                        (d_species_gamma[d_num_species - 1] - Real(1)) -
                    (p[idx_pressure] + d_species_gamma[si]*d_species_p_inf[si])/
                        (d_species_gamma[si] - Real(1)));
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
            Real* M_i = M[si];
            
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
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
                        (gamma[idx_mixture_thermo_properties] - Real(1))*
                        ((p[idx_pressure] + d_species_gamma[d_num_species - 1]*d_species_p_inf[d_num_species - 1])/
                            (d_species_gamma[d_num_species - 1] - Real(1)) -
                        (p[idx_pressure] + d_species_gamma[si]*d_species_p_inf[si])/
                            (d_species_gamma[si] - Real(1)));
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
            Real* M_i = M[si];
            
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
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
                            (gamma[idx_mixture_thermo_properties] - Real(1))*
                            ((p[idx_pressure] + d_species_gamma[d_num_species - 1]*d_species_p_inf[d_num_species - 1])/
                                (d_species_gamma[d_num_species - 1] - Real(1)) -
                            (p[idx_pressure] + d_species_gamma[si]*d_species_p_inf[si])/
                                (d_species_gamma[si] - Real(1)));
                    }
                }
            }
        }
    }
}


/*
 * Compute the thermodynamic properties of the mixture with volume fractions.
 */
void
EquationOfStateMixingRulesStiffenedGas::computeMixtureThermodynamicPropertiesWithVolumeFractions(
    Real* const gamma,
    Real* const p_inf,
    const std::vector<const Real*>& Z,
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
        
        // Compute temporary variables and store them in the data of gamma of p_inf temporarily.
        for (int si = 0; si < d_num_species; si++)
        {
            const Real tmp_1 = Real(1)/(d_species_gamma[si] - Real(1));
            const Real tmp_2 = d_species_gamma[si]*d_species_p_inf[si]/(d_species_gamma[si] - Real(1));
            
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_mixture_thermo_properties = i + offset_0_mixture_thermo_properties;
                const int idx_volume_fractions = i + offset_0_volume_fractions;
                
                gamma[idx_mixture_thermo_properties] += Z[si][idx_volume_fractions]*tmp_1;
                p_inf[idx_mixture_thermo_properties] += Z[si][idx_volume_fractions]*tmp_2;
            }
        }
        
        // Compute gamma and p_inf.
        HAMERS_PRAGMA_SIMD
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear index.
            const int idx_mixture_thermo_properties = i + offset_0_mixture_thermo_properties;
            
            gamma[idx_mixture_thermo_properties] = Real(1)/gamma[idx_mixture_thermo_properties] + Real(1);
            
            p_inf[idx_mixture_thermo_properties] = (gamma[idx_mixture_thermo_properties] - Real(1))/
                gamma[idx_mixture_thermo_properties]*p_inf[idx_mixture_thermo_properties];
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
        
        // Compute temporary variables and store them in the data of gamma of p_inf temporarily.
        for (int si = 0; si < d_num_species; si++)
        {
            const Real tmp_1 = Real(1)/(d_species_gamma[si] - Real(1));
            const Real tmp_2 = d_species_gamma[si]*d_species_p_inf[si]/(d_species_gamma[si] - Real(1));
            
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_mixture_thermo_properties = (i + offset_0_mixture_thermo_properties) +
                        (j + offset_1_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties;
                    
                    const int idx_volume_fractions = (i + offset_0_volume_fractions) +
                        (j + offset_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                    
                    gamma[idx_mixture_thermo_properties] += Z[si][idx_volume_fractions]*tmp_1;
                    p_inf[idx_mixture_thermo_properties] += Z[si][idx_volume_fractions]*tmp_2;
                }
            }
        }
        
        // Compute gamma and p_inf.
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear index.
                const int idx_mixture_thermo_properties = (i + offset_0_mixture_thermo_properties) +
                    (j + offset_1_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties;
                
                gamma[idx_mixture_thermo_properties] = Real(1)/gamma[idx_mixture_thermo_properties] + Real(1);
                
                p_inf[idx_mixture_thermo_properties] = (gamma[idx_mixture_thermo_properties] - Real(1))/
                    gamma[idx_mixture_thermo_properties]*p_inf[idx_mixture_thermo_properties];
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
        
        // Compute temporary variables and store them in the data of gamma of p_inf temporarily.
        for (int si = 0; si < d_num_species; si++)
        {
            const Real tmp_1 = Real(1)/(d_species_gamma[si] - Real(1));
            const Real tmp_2 = d_species_gamma[si]*d_species_p_inf[si]/(d_species_gamma[si] - Real(1));
            
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
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
                        
                        gamma[idx_mixture_thermo_properties] += Z[si][idx_volume_fractions]*tmp_1;
                        p_inf[idx_mixture_thermo_properties] += Z[si][idx_volume_fractions]*tmp_2;
                    }
                }
            }
        }
        
        // Compute gamma and p_inf.
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_mixture_thermo_properties = (i + offset_0_mixture_thermo_properties) +
                        (j + offset_1_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties +
                        (k + offset_2_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties*
                            ghostcell_dim_1_mixture_thermo_properties;
                    
                    gamma[idx_mixture_thermo_properties] = Real(1)/gamma[idx_mixture_thermo_properties] + Real(1);
                    
                    p_inf[idx_mixture_thermo_properties] = (gamma[idx_mixture_thermo_properties] - Real(1))/
                        gamma[idx_mixture_thermo_properties]*p_inf[idx_mixture_thermo_properties];
                }
            }
        }
    }
}


/*
 * Compute the thermodynamic properties of the mixture with volume fractions.
 */
void
EquationOfStateMixingRulesStiffenedGas::computeMixtureThermodynamicPropertiesWithVolumeFractions(
    Real* const gamma,
    Real* const p_inf,
    Real* const Z_last,
    const std::vector<const Real*>& Z,
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
        
        // Compute temporary variables and store them in the data of gamma of p_inf temporarily.
        for (int si = 0; si < d_num_species - 1; si++)
        {
            const Real tmp_1 = Real(1)/(d_species_gamma[si] - Real(1));
            const Real tmp_2 = d_species_gamma[si]*d_species_p_inf[si]/(d_species_gamma[si] - Real(1));
            
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_mixture_thermo_properties = i + offset_0_mixture_thermo_properties;
                const int idx_volume_fractions_last = i + offset_0_volume_fractions_last;
                const int idx_volume_fractions = i + offset_0_volume_fractions;
                
                gamma[idx_mixture_thermo_properties] += Z[si][idx_volume_fractions]*tmp_1;
                p_inf[idx_mixture_thermo_properties] += Z[si][idx_volume_fractions]*tmp_2;
                
                // Compute the volume fraction of the last species.
                Z_last[idx_volume_fractions_last] -= Z[si][idx_volume_fractions];
            }
        }
        
        // Add the contribution from the last species and compute gamma and p_inf.
        
        const Real tmp_1 = Real(1)/(d_species_gamma.back() - Real(1));
        const Real tmp_2 = d_species_gamma.back()*d_species_p_inf.back()/(d_species_gamma.back() - Real(1));
        
        HAMERS_PRAGMA_SIMD
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_mixture_thermo_properties = i + offset_0_mixture_thermo_properties;
            const int idx_volume_fractions_last = i + offset_0_volume_fractions_last;
            
            gamma[idx_mixture_thermo_properties] += Z_last[idx_volume_fractions_last]*tmp_1;
            p_inf[idx_mixture_thermo_properties] += Z_last[idx_volume_fractions_last]*tmp_2;
            
            gamma[idx_mixture_thermo_properties] = Real(1)/gamma[idx_mixture_thermo_properties] + Real(1);
            
            p_inf[idx_mixture_thermo_properties] = (gamma[idx_mixture_thermo_properties] - Real(1))/
                gamma[idx_mixture_thermo_properties]*p_inf[idx_mixture_thermo_properties];
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
        
        // Compute temporary variables and store them in the data of gamma of p_inf temporarily.
        for (int si = 0; si < d_num_species - 1; si++)
        {
            const Real tmp_1 = Real(1)/(d_species_gamma[si] - Real(1));
            const Real tmp_2 = d_species_gamma[si]*d_species_p_inf[si]/(d_species_gamma[si] - Real(1));
            
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_mixture_thermo_properties = (i + offset_0_mixture_thermo_properties) +
                        (j + offset_1_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties;
                    
                    const int idx_volume_fractions_last = (i + offset_0_volume_fractions_last) +
                        (j + offset_1_volume_fractions_last)*ghostcell_dim_0_volume_fractions_last;
                    
                    const int idx_volume_fractions = (i + offset_0_volume_fractions) +
                        (j + offset_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                    
                    gamma[idx_mixture_thermo_properties] += Z[si][idx_volume_fractions]*tmp_1;
                    p_inf[idx_mixture_thermo_properties] += Z[si][idx_volume_fractions]*tmp_2;
                    
                    // Compute the volume fraction of the last species.
                    Z_last[idx_volume_fractions_last] -= Z[si][idx_volume_fractions];
                }
            }
        }
        
        // Add the contribution from the last species and compute gamma and p_inf.
        
        const Real tmp_1 = Real(1)/(d_species_gamma.back() - Real(1));
        const Real tmp_2 = d_species_gamma.back()*d_species_p_inf.back()/(d_species_gamma.back() - Real(1));
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_mixture_thermo_properties = (i + offset_0_mixture_thermo_properties) +
                    (j + offset_1_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties;
                
                const int idx_volume_fractions_last = (i + offset_0_volume_fractions_last) +
                    (j + offset_1_volume_fractions_last)*ghostcell_dim_0_volume_fractions_last;
                
                gamma[idx_mixture_thermo_properties] += Z_last[idx_volume_fractions_last]*tmp_1;
                p_inf[idx_mixture_thermo_properties] += Z_last[idx_volume_fractions_last]*tmp_2;
                
                gamma[idx_mixture_thermo_properties] = Real(1)/gamma[idx_mixture_thermo_properties] + Real(1);
                
                p_inf[idx_mixture_thermo_properties] = (gamma[idx_mixture_thermo_properties] - Real(1))/
                    gamma[idx_mixture_thermo_properties]*p_inf[idx_mixture_thermo_properties];
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
        
        // Compute temporary variables and store them in the data of gamma of p_inf temporarily.
        for (int si = 0; si < d_num_species - 1; si++)
        {
            const Real tmp_1 = Real(1)/(d_species_gamma[si] - Real(1));
            const Real tmp_2 = d_species_gamma[si]*d_species_p_inf[si]/(d_species_gamma[si] - Real(1));
            
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
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
                        
                        gamma[idx_mixture_thermo_properties] += Z[si][idx_volume_fractions]*tmp_1;
                        p_inf[idx_mixture_thermo_properties] += Z[si][idx_volume_fractions]*tmp_2;
                        
                        // Compute the volume fraction of the last species.
                        Z_last[idx_volume_fractions_last] -= Z[si][idx_volume_fractions];
                    }
                }
            }
        }
        
        // Add the contribution from the last species and compute gamma and p_inf.
        
        const Real tmp_1 = Real(1)/(d_species_gamma.back() - Real(1));
        const Real tmp_2 = d_species_gamma.back()*d_species_p_inf.back()/(d_species_gamma.back() - Real(1));
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
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
                        
                    gamma[idx_mixture_thermo_properties] += Z_last[idx_volume_fractions_last]*tmp_1;
                    p_inf[idx_mixture_thermo_properties] += Z_last[idx_volume_fractions_last]*tmp_2;
                    
                    gamma[idx_mixture_thermo_properties] = Real(1)/gamma[idx_mixture_thermo_properties] + Real(1);
                    
                    p_inf[idx_mixture_thermo_properties] = (gamma[idx_mixture_thermo_properties] - Real(1))/
                        gamma[idx_mixture_thermo_properties]*p_inf[idx_mixture_thermo_properties];
                }
            }
        }
    }
}
