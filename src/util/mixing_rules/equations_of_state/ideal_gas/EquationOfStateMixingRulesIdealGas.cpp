#include "util/mixing_rules/equations_of_state/ideal_gas/EquationOfStateMixingRulesIdealGas.hpp"

EquationOfStateMixingRulesIdealGas::EquationOfStateMixingRulesIdealGas(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const int& num_species,
    const MIXING_CLOSURE_MODEL::TYPE& mixing_closure_model,
    const boost::shared_ptr<tbox::Database>& equation_of_state_mixing_rules_db):
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
    const boost::shared_ptr<tbox::Database>& restart_db) const
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
    boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
    const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
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
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_pressure->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_density->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_internal_energy->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_mass_fractions->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells.
     */
    
    const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
    const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
    const hier::IntVector num_ghosts_internal_energy = data_internal_energy->getGhostCellWidth();
    const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
    
    /*
     * Get the minimum number of ghost cells for mixture thermodynamic properties.
     */
    
    hier::IntVector num_ghosts_min(d_dim);
    
    num_ghosts_min = num_ghosts_pressure;
    num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
    num_ghosts_min = hier::IntVector::min(num_ghosts_internal_energy, num_ghosts_min);
    num_ghosts_min = hier::IntVector::min(num_ghosts_mass_fractions, num_ghosts_min);
    
    /*
     * Get the mixture thermodyanmic properties.
     */
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    boost::shared_ptr<pdat::CellData<double> > data_mixture_thermo_properties(
        new pdat::CellData<double>(interior_box, num_thermo_properties, num_ghosts_min));
    
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
    boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_internal_energy,
    const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
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
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_pressure->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_density->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_internal_energy->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_mass_fractions->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells.
     */
    
    const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
    const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
    const hier::IntVector num_ghosts_internal_energy = data_internal_energy->getGhostCellWidth();
    const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
    
    /*
     * Get the minimum number of ghost cells for mixture thermodynamic properties.
     */
    
    hier::IntVector num_ghosts_min(d_dim);
    
    num_ghosts_min = num_ghosts_pressure;
    num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
    num_ghosts_min = hier::IntVector::min(num_ghosts_internal_energy, num_ghosts_min);
    num_ghosts_min = hier::IntVector::min(num_ghosts_mass_fractions, num_ghosts_min);
    
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
    
    boost::shared_ptr<pdat::SideData<double> > data_mixture_thermo_properties(
        new pdat::SideData<double>(interior_box, num_thermo_properties, num_ghosts_min,
            direction));
    
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
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOBARIC);
    TBOX_ASSERT((static_cast<int>(volume_fractions.size()) == d_num_species) ||
                (static_cast<int>(volume_fractions.size()) == d_num_species - 1));
#endif
    
    NULL_USE(mass_fractions);
    
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
    boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
    const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
    const boost::shared_ptr<pdat::CellData<double> >& data_volume_fractions,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOBARIC);
    
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_internal_energy);
    TBOX_ASSERT(data_volume_fractions);
    
    TBOX_ASSERT((data_volume_fractions->getDepth() == d_num_species) ||
                (data_volume_fractions->getDepth() == d_num_species - 1));
#endif
    
    NULL_USE(data_mass_fractions);
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_pressure->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_density->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_internal_energy->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_volume_fractions->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells.
     */
    
    const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
    const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
    const hier::IntVector num_ghosts_internal_energy = data_internal_energy->getGhostCellWidth();
    const hier::IntVector num_ghosts_volume_fractions = data_volume_fractions->getGhostCellWidth();
    
    /*
     * Get the minimum number of ghost cells for mixture thermodynamic properties.
     */
    
    hier::IntVector num_ghosts_min(d_dim);
    
    num_ghosts_min = num_ghosts_pressure;
    num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
    num_ghosts_min = hier::IntVector::min(num_ghosts_internal_energy, num_ghosts_min);
    num_ghosts_min = hier::IntVector::min(num_ghosts_volume_fractions, num_ghosts_min);
    
    /*
     * Get the mixture thermodyanmic properties.
     */
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    boost::shared_ptr<pdat::CellData<double> > data_mixture_thermo_properties(
        new pdat::CellData<double>(interior_box, num_thermo_properties, num_ghosts_min));
    
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
    boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_internal_energy,
    const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
    const boost::shared_ptr<pdat::SideData<double> >& data_volume_fractions,
    int side_normal,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOBARIC);
    
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_internal_energy);
    TBOX_ASSERT(data_volume_fractions);
    
    TBOX_ASSERT((data_volume_fractions->getDepth() == d_num_species) ||
                (data_volume_fractions->getDepth() == d_num_species - 1));
#endif
    
    NULL_USE(data_mass_fractions);
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_pressure->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_density->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_internal_energy->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_volume_fractions->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells.
     */
    
    const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
    const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
    const hier::IntVector num_ghosts_internal_energy = data_internal_energy->getGhostCellWidth();
    const hier::IntVector num_ghosts_volume_fractions = data_volume_fractions->getGhostCellWidth();
    
    /*
     * Get the minimum number of ghost cells for mixture thermodynamic properties.
     */
    
    hier::IntVector num_ghosts_min(d_dim);
    
    num_ghosts_min = num_ghosts_pressure;
    num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
    num_ghosts_min = hier::IntVector::min(num_ghosts_internal_energy, num_ghosts_min);
    num_ghosts_min = hier::IntVector::min(num_ghosts_volume_fractions, num_ghosts_min);
    
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
    
    boost::shared_ptr<pdat::SideData<double> > data_mixture_thermo_properties(
        new pdat::SideData<double>(interior_box, num_thermo_properties, num_ghosts_min,
            direction));
    
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
 * Compute the sound speed of the mixture with isothermal and isobaric equilibrium assumptions.
 */
double
EquationOfStateMixingRulesIdealGas::getSoundSpeed(
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
    
    return d_equation_of_state->getSoundSpeed(
        density,
        pressure,
        mixture_thermo_properties_const_ptr);
}


/*
 * Compute the sound speed of the mixture with isothermal and isobaric equilibrium assumptions.
 */
void
EquationOfStateMixingRulesIdealGas::computeSoundSpeed(
    boost::shared_ptr<pdat::CellData<double> >& data_sound_speed,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT((d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC) ||
                (d_mixing_closure_model == MIXING_CLOSURE_MODEL::NO_MODEL && d_num_species == 1));
    
    TBOX_ASSERT(data_sound_speed);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_mass_fractions);
    
    TBOX_ASSERT((data_mass_fractions->getDepth() == d_num_species) ||
                (data_mass_fractions->getDepth() == d_num_species - 1));
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_sound_speed->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_density->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_pressure->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_mass_fractions->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells.
     */
    
    const hier::IntVector num_ghosts_sound_speed = data_sound_speed->getGhostCellWidth();
    const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
    const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
    const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
    
    /*
     * Get the minimum number of ghost cells for mixture thermodynamic properties.
     */
    
    hier::IntVector num_ghosts_min(d_dim);
    
    num_ghosts_min = num_ghosts_sound_speed;
    num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
    num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
    num_ghosts_min = hier::IntVector::min(num_ghosts_mass_fractions, num_ghosts_min);
    
    /*
     * Get the mixture thermodyanmic properties.
     */
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    boost::shared_ptr<pdat::CellData<double> > data_mixture_thermo_properties(
        new pdat::CellData<double>(interior_box, num_thermo_properties, num_ghosts_min));
    
    computeMixtureThermodynamicProperties(
        data_mixture_thermo_properties,
        data_mass_fractions,
        domain);
    
    d_equation_of_state->computeSoundSpeed(
        data_sound_speed,
        data_density,
        data_pressure,
        data_mixture_thermo_properties,
        domain);
}


/*
 * Compute the sound speed of the mixture with isothermal and isobaric equilibrium assumptions.
 */
void
EquationOfStateMixingRulesIdealGas::computeSoundSpeed(
    boost::shared_ptr<pdat::SideData<double> >& data_sound_speed,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
    int side_normal,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT((d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC) ||
                (d_mixing_closure_model == MIXING_CLOSURE_MODEL::NO_MODEL && d_num_species == 1));
    
    TBOX_ASSERT(data_sound_speed);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_mass_fractions);
    
    TBOX_ASSERT((data_mass_fractions->getDepth() == d_num_species) ||
                (data_mass_fractions->getDepth() == d_num_species - 1));
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_sound_speed->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_density->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_pressure->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_mass_fractions->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells.
     */
    
    const hier::IntVector num_ghosts_sound_speed = data_sound_speed->getGhostCellWidth();
    const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
    const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
    const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
    
    /*
     * Get the minimum number of ghost cells for mixture thermodynamic properties.
     */
    
    hier::IntVector num_ghosts_min(d_dim);
    
    num_ghosts_min = num_ghosts_sound_speed;
    num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
    num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
    num_ghosts_min = hier::IntVector::min(num_ghosts_mass_fractions, num_ghosts_min);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(side_normal < d_dim.getValue());
    
    TBOX_ASSERT(data_sound_speed->getDirectionVector()[side_normal] > 0);
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
    
    boost::shared_ptr<pdat::SideData<double> > data_mixture_thermo_properties(
        new pdat::SideData<double>(interior_box, num_thermo_properties, num_ghosts_min,
            direction));
    
    computeMixtureThermodynamicProperties(
        data_mixture_thermo_properties,
        data_mass_fractions,
        side_normal,
        domain);
    
    d_equation_of_state->computeSoundSpeed(
        data_sound_speed,
        data_density,
        data_pressure,
        data_mixture_thermo_properties,
        side_normal,
        domain);
}


/*
 * Compute the sound speed of the mixture with isobaric equilibrium assumption.
 */
double
EquationOfStateMixingRulesIdealGas::getSoundSpeed(
    const double* const density,
    const double* const pressure,
    const std::vector<const double*>& mass_fractions,
    const std::vector<const double*>& volume_fractions) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOBARIC);
    TBOX_ASSERT((static_cast<int>(volume_fractions.size()) == d_num_species) ||
                (static_cast<int>(volume_fractions.size()) == d_num_species - 1));
#endif
    
    NULL_USE(mass_fractions);
    
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
    
    return d_equation_of_state->getSoundSpeed(
        density,
        pressure,
        mixture_thermo_properties_const_ptr);
}


/*
 * Compute the sound speed of the mixture with isobaric equilibrium assumption.
 */
void
EquationOfStateMixingRulesIdealGas::computeSoundSpeed(
    boost::shared_ptr<pdat::CellData<double> >& data_sound_speed,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
    const boost::shared_ptr<pdat::CellData<double> >& data_volume_fractions,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOBARIC);
    
    TBOX_ASSERT(data_sound_speed);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_volume_fractions);
    
    TBOX_ASSERT((data_volume_fractions->getDepth() == d_num_species) ||
                (data_volume_fractions->getDepth() == d_num_species - 1));
#endif
    
    NULL_USE(data_mass_fractions);
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_sound_speed->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_density->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_pressure->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_volume_fractions->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells.
     */
    
    const hier::IntVector num_ghosts_sound_speed = data_sound_speed->getGhostCellWidth();
    const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
    const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
    const hier::IntVector num_ghosts_volume_fractions = data_volume_fractions->getGhostCellWidth();
    
    /*
     * Get the minimum number of ghost cells for mixture thermodynamic properties.
     */
    
    hier::IntVector num_ghosts_min(d_dim);
    
    num_ghosts_min = num_ghosts_sound_speed;
    num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
    num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
    num_ghosts_min = hier::IntVector::min(num_ghosts_volume_fractions, num_ghosts_min);
    
    /*
     * Get the mixture thermodyanmic properties.
     */
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    boost::shared_ptr<pdat::CellData<double> > data_mixture_thermo_properties(
        new pdat::CellData<double>(interior_box, num_thermo_properties, num_ghosts_min));
    
    computeMixtureThermodynamicProperties(
        data_mixture_thermo_properties,
        data_volume_fractions,
        domain);
    
    d_equation_of_state->computeSoundSpeed(
        data_sound_speed,
        data_density,
        data_pressure,
        data_mixture_thermo_properties,
        domain);
}


/*
 * Compute the sound speed of the mixture with isobaric equilibrium assumption.
 */
void
EquationOfStateMixingRulesIdealGas::computeSoundSpeed(
    boost::shared_ptr<pdat::SideData<double> >& data_sound_speed,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
    const boost::shared_ptr<pdat::SideData<double> >& data_volume_fractions,
    int side_normal,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOBARIC);
    
    TBOX_ASSERT(data_sound_speed);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_volume_fractions);
    
    TBOX_ASSERT((data_volume_fractions->getDepth() == d_num_species) ||
                (data_volume_fractions->getDepth() == d_num_species - 1));
#endif
    
    NULL_USE(data_mass_fractions);
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_sound_speed->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_density->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_pressure->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_volume_fractions->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells.
     */
    
    const hier::IntVector num_ghosts_sound_speed = data_sound_speed->getGhostCellWidth();
    const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
    const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
    const hier::IntVector num_ghosts_volume_fractions = data_volume_fractions->getGhostCellWidth();
    
    /*
     * Get the minimum number of ghost cells for mixture thermodynamic properties.
     */
    
    hier::IntVector num_ghosts_min(d_dim);
    
    num_ghosts_min = num_ghosts_sound_speed;
    num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
    num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
    num_ghosts_min = hier::IntVector::min(num_ghosts_volume_fractions, num_ghosts_min);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(side_normal < d_dim.getValue());
    
    TBOX_ASSERT(data_sound_speed->getDirectionVector()[side_normal] > 0);
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
    
    boost::shared_ptr<pdat::SideData<double> > data_mixture_thermo_properties(
        new pdat::SideData<double>(interior_box, num_thermo_properties, num_ghosts_min,
            direction));
    
    computeMixtureThermodynamicProperties(
        data_mixture_thermo_properties,
        data_volume_fractions,
        side_normal,
        domain);
    
    d_equation_of_state->computeSoundSpeed(
        data_sound_speed,
        data_density,
        data_pressure,
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
    boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
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
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_internal_energy->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_density->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_pressure->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_mass_fractions->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells.
     */
    
    const hier::IntVector num_ghosts_internal_energy = data_internal_energy->getGhostCellWidth();
    const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
    const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
    const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
    
    /*
     * Get the minimum number of ghost cells for mixture thermodynamic properties.
     */
    
    hier::IntVector num_ghosts_min(d_dim);
    
    num_ghosts_min = num_ghosts_internal_energy;
    num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
    num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
    num_ghosts_min = hier::IntVector::min(num_ghosts_mass_fractions, num_ghosts_min);
    
    /*
     * Get the mixture thermodyanmic properties.
     */
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    boost::shared_ptr<pdat::CellData<double> > data_mixture_thermo_properties(
        new pdat::CellData<double>(interior_box, num_thermo_properties, num_ghosts_min));
    
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
    boost::shared_ptr<pdat::SideData<double> >& data_internal_energy,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
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
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_internal_energy->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_density->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_pressure->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_mass_fractions->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells.
     */
    
    const hier::IntVector num_ghosts_internal_energy = data_internal_energy->getGhostCellWidth();
    const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
    const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
    const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
    
    /*
     * Get the minimum number of ghost cells for mixture thermodynamic properties.
     */
    
    hier::IntVector num_ghosts_min(d_dim);
    
    num_ghosts_min = num_ghosts_internal_energy;
    num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
    num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
    num_ghosts_min = hier::IntVector::min(num_ghosts_mass_fractions, num_ghosts_min);
    
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
    
    boost::shared_ptr<pdat::SideData<double> > data_mixture_thermo_properties(
        new pdat::SideData<double>(interior_box, num_thermo_properties, num_ghosts_min,
            direction));
    
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
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOBARIC);
    TBOX_ASSERT((static_cast<int>(volume_fractions.size()) == d_num_species) ||
                (static_cast<int>(volume_fractions.size()) == d_num_species - 1));
#endif
    
    NULL_USE(mass_fractions);
    
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
    boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
    const boost::shared_ptr<pdat::CellData<double> >& data_volume_fractions,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOBARIC);
    
    TBOX_ASSERT(data_internal_energy);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_volume_fractions);
    
    TBOX_ASSERT((data_volume_fractions->getDepth() == d_num_species) ||
                (data_volume_fractions->getDepth() == d_num_species - 1));
#endif
    
    NULL_USE(data_mass_fractions);
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_internal_energy->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_density->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_pressure->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_volume_fractions->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells.
     */
    
    const hier::IntVector num_ghosts_internal_energy = data_internal_energy->getGhostCellWidth();
    const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
    const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
    const hier::IntVector num_ghosts_volume_fractions = data_volume_fractions->getGhostCellWidth();
    
    /*
     * Get the minimum number of ghost cells for mixture thermodynamic properties.
     */
    
    hier::IntVector num_ghosts_min(d_dim);
    
    num_ghosts_min = num_ghosts_internal_energy;
    num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
    num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
    num_ghosts_min = hier::IntVector::min(num_ghosts_volume_fractions, num_ghosts_min);
    
    /*
     * Get the mixture thermodyanmic properties.
     */
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    boost::shared_ptr<pdat::CellData<double> > data_mixture_thermo_properties(
        new pdat::CellData<double>(interior_box, num_thermo_properties, num_ghosts_min));
    
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
    boost::shared_ptr<pdat::SideData<double> >& data_internal_energy,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
    const boost::shared_ptr<pdat::SideData<double> >& data_volume_fractions,
    int side_normal,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOBARIC);
    
    TBOX_ASSERT(data_internal_energy);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_volume_fractions);
    
    TBOX_ASSERT((data_volume_fractions->getDepth() == d_num_species) ||
                (data_volume_fractions->getDepth() == d_num_species - 1));
#endif
    
    NULL_USE(data_mass_fractions);
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_internal_energy->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_density->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_pressure->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_volume_fractions->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells.
     */
    
    const hier::IntVector num_ghosts_internal_energy = data_internal_energy->getGhostCellWidth();
    const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
    const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
    const hier::IntVector num_ghosts_volume_fractions = data_volume_fractions->getGhostCellWidth();
    
    /*
     * Get the minimum number of ghost cells for mixture thermodynamic properties.
     */
    
    hier::IntVector num_ghosts_min(d_dim);
    
    num_ghosts_min = num_ghosts_internal_energy;
    num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
    num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
    num_ghosts_min = hier::IntVector::min(num_ghosts_volume_fractions, num_ghosts_min);
    
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
    
    boost::shared_ptr<pdat::SideData<double> > data_mixture_thermo_properties(
        new pdat::SideData<double>(interior_box, num_thermo_properties, num_ghosts_min,
            direction));
    
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
    boost::shared_ptr<pdat::CellData<double> >& data_temperature,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
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
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_temperature->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_density->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_pressure->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_mass_fractions->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells.
     */
    
    const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
    const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
    const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
    const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
    
    /*
     * Get the minimum number of ghost cells for mixture thermodynamic properties.
     */
    
    hier::IntVector num_ghosts_min(d_dim);
    
    num_ghosts_min = num_ghosts_temperature;
    num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
    num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
    num_ghosts_min = hier::IntVector::min(num_ghosts_mass_fractions, num_ghosts_min);
    
    /*
     * Get the mixture thermodyanmic properties.
     */
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    boost::shared_ptr<pdat::CellData<double> > data_mixture_thermo_properties(
        new pdat::CellData<double>(interior_box, num_thermo_properties, num_ghosts_min));
    
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
    boost::shared_ptr<pdat::SideData<double> >& data_temperature,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
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
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_temperature->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_density->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_pressure->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_mass_fractions->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells.
     */
    
    const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
    const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
    const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
    const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
    
    /*
     * Get the minimum number of ghost cells for mixture thermodynamic properties.
     */
    
    hier::IntVector num_ghosts_min(d_dim);
    
    num_ghosts_min = num_ghosts_temperature;
    num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
    num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
    num_ghosts_min = hier::IntVector::min(num_ghosts_mass_fractions, num_ghosts_min);
    
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
    
    boost::shared_ptr<pdat::SideData<double> > data_mixture_thermo_properties(
        new pdat::SideData<double>(interior_box, num_thermo_properties, num_ghosts_min,
            direction));
    
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
    boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
    const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
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
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_internal_energy->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_density->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_temperature->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_mass_fractions->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells.
     */
    
    const hier::IntVector num_ghosts_internal_energy = data_internal_energy->getGhostCellWidth();
    const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
    const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
    const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
    
    /*
     * Get the minimum number of ghost cells for mixture thermodynamic properties.
     */
    
    hier::IntVector num_ghosts_min(d_dim);
    
    num_ghosts_min = num_ghosts_internal_energy;
    num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
    num_ghosts_min = hier::IntVector::min(num_ghosts_temperature, num_ghosts_min);
    num_ghosts_min = hier::IntVector::min(num_ghosts_mass_fractions, num_ghosts_min);
    
    /*
     * Get the mixture thermodyanmic properties.
     */
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    boost::shared_ptr<pdat::CellData<double> > data_mixture_thermo_properties(
        new pdat::CellData<double>(interior_box, num_thermo_properties, num_ghosts_min));
    
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
    boost::shared_ptr<pdat::SideData<double> >& data_internal_energy,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_temperature,
    const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
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
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_internal_energy->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_density->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_temperature->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_mass_fractions->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells.
     */
    
    const hier::IntVector num_ghosts_internal_energy = data_internal_energy->getGhostCellWidth();
    const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
    const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
    const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
    
    /*
     * Get the minimum number of ghost cells for mixture thermodynamic properties.
     */
    
    hier::IntVector num_ghosts_min(d_dim);
    
    num_ghosts_min = num_ghosts_internal_energy;
    num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
    num_ghosts_min = hier::IntVector::min(num_ghosts_temperature, num_ghosts_min);
    num_ghosts_min = hier::IntVector::min(num_ghosts_mass_fractions, num_ghosts_min);
    
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
    
    boost::shared_ptr<pdat::SideData<double> > data_mixture_thermo_properties(
        new pdat::SideData<double>(interior_box, num_thermo_properties, num_ghosts_min));
    
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
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT((d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC) ||
                (d_mixing_closure_model == MIXING_CLOSURE_MODEL::NO_MODEL && d_num_species == 1));
    TBOX_ASSERT((static_cast<int>(mass_fractions.size()) == d_num_species) ||
                (static_cast<int>(mass_fractions.size()) == d_num_species - 1));
#endif
    
    NULL_USE(density);
    NULL_USE(pressure);
    
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
    boost::shared_ptr<pdat::CellData<double> >& data_isochoric_specific_heat_capacity,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT((d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC) ||
                (d_mixing_closure_model == MIXING_CLOSURE_MODEL::NO_MODEL && d_num_species == 1));
    
    TBOX_ASSERT(data_isochoric_specific_heat_capacity);
    TBOX_ASSERT(data_mass_fractions);
    
    TBOX_ASSERT((data_mass_fractions->getDepth() == d_num_species) ||
                (data_mass_fractions->getDepth() == d_num_species - 1));
#endif
    
    NULL_USE(data_density);
    NULL_USE(data_pressure);
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_isochoric_specific_heat_capacity->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_mass_fractions->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_isochoric_specific_heat_capacity =
        data_isochoric_specific_heat_capacity->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_isochoric_specific_heat_capacity =
        data_isochoric_specific_heat_capacity->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_mass_fractions =
        data_mass_fractions->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_isochoric_specific_heat_capacity;
        num_ghosts_min = hier::IntVector::min(num_ghosts_mass_fractions, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
        TBOX_ASSERT(data_isochoric_specific_heat_capacity->getGhostBox().contains(domain));
        TBOX_ASSERT(data_mass_fractions->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
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
            num_ghosts_isochoric_specific_heat_capacity,
            num_ghosts_mass_fractions,
            ghostcell_dims_isochoric_specific_heat_capacity,
            ghostcell_dims_mass_fractions,
            domain_lo,
            domain_dims);
    }
    else if (data_mass_fractions->getDepth() == d_num_species - 1)
    {
        boost::shared_ptr<pdat::CellData<double> > data_mass_fractions_last(
            new pdat::CellData<double>(interior_box, 1, num_ghosts_mass_fractions));
        
        if (domain.empty())
        {
            data_mass_fractions_last->fillAll(double(1));
        }
        else
        {
            data_mass_fractions_last->fillAll(double(1), domain);
        }
        
        /*
         * Get the pointers to the cell data of mass fractions.
         */
        
        std::vector<const double*> Y;
        Y.reserve(d_num_species - 1);
        for (int si = 0; si < d_num_species - 1; si++)
        {
            Y.push_back(data_mass_fractions->getPointer(si));
        }
        
        double* const Y_last = data_mass_fractions_last->getPointer(0);
        
        computeIsochoricSpecificHeatCapacity(
            c_v,
            Y_last,
            Y,
            num_ghosts_isochoric_specific_heat_capacity,
            num_ghosts_mass_fractions,
            ghostcell_dims_isochoric_specific_heat_capacity,
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
    boost::shared_ptr<pdat::SideData<double> >& data_isochoric_specific_heat_capacity,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
    int side_normal,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT((d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC) ||
                (d_mixing_closure_model == MIXING_CLOSURE_MODEL::NO_MODEL && d_num_species == 1));
    
    TBOX_ASSERT(data_isochoric_specific_heat_capacity);
    TBOX_ASSERT(data_mass_fractions);
    
    TBOX_ASSERT((data_mass_fractions->getDepth() == d_num_species) ||
                (data_mass_fractions->getDepth() == d_num_species - 1));
#endif
    
    NULL_USE(data_density);
    NULL_USE(data_pressure);
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_isochoric_specific_heat_capacity->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_mass_fractions->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_isochoric_specific_heat_capacity =
        data_isochoric_specific_heat_capacity->getGhostCellWidth();
    hier::IntVector ghostcell_dims_isochoric_specific_heat_capacity =
        data_isochoric_specific_heat_capacity->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
    hier::IntVector ghostcell_dims_mass_fractions =
        data_mass_fractions->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_isochoric_specific_heat_capacity;
        num_ghosts_min = hier::IntVector::min(num_ghosts_mass_fractions, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
        TBOX_ASSERT(data_isochoric_specific_heat_capacity->getGhostBox().contains(domain));
        TBOX_ASSERT(data_mass_fractions->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(side_normal < d_dim.getValue());
    
    TBOX_ASSERT(data_isochoric_specific_heat_capacity->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_mass_fractions->getDirectionVector()[side_normal] > 0);
#endif
    
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
            num_ghosts_isochoric_specific_heat_capacity,
            num_ghosts_mass_fractions,
            ghostcell_dims_isochoric_specific_heat_capacity,
            ghostcell_dims_mass_fractions,
            domain_lo,
            domain_dims);
    }
    else if (data_mass_fractions->getDepth() == d_num_species - 1)
    {
        hier::IntVector direction = hier::IntVector::getZero(d_dim);
        direction[side_normal] = 1;
        
        boost::shared_ptr<pdat::SideData<double> > data_mass_fractions_last(
            new pdat::SideData<double>(interior_box, 1, num_ghosts_mass_fractions,
                direction));
        
        if (domain.empty())
        {
            data_mass_fractions_last->fillAll(double(1));
        }
        else
        {
            data_mass_fractions_last->fillAll(double(1), domain);
        }
        
        /*
         * Get the pointers to the cell data of mass fractions.
         */
        
        std::vector<const double*> Y;
        Y.reserve(d_num_species - 1);
        for (int si = 0; si < d_num_species - 1; si++)
        {
            Y.push_back(data_mass_fractions->getPointer(side_normal, si));
        }
        
        double* const Y_last = data_mass_fractions_last->getPointer(side_normal, 0);
        
        computeIsochoricSpecificHeatCapacity(
            c_v,
            Y_last,
            Y,
            num_ghosts_isochoric_specific_heat_capacity,
            num_ghosts_mass_fractions,
            ghostcell_dims_isochoric_specific_heat_capacity,
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
    boost::shared_ptr<pdat::CellData<double> >& data_isobaric_specific_heat_capacity,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT((d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC) ||
                (d_mixing_closure_model == MIXING_CLOSURE_MODEL::NO_MODEL && d_num_species == 1));
    
    TBOX_ASSERT(data_isobaric_specific_heat_capacity);
    TBOX_ASSERT(data_mass_fractions);
    
    TBOX_ASSERT((data_mass_fractions->getDepth() == d_num_species) ||
                (data_mass_fractions->getDepth() == d_num_species - 1));
#endif
    
    NULL_USE(data_density);
    NULL_USE(data_pressure);
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_isobaric_specific_heat_capacity->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_mass_fractions->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_isobaric_specific_heat_capacity =
        data_isobaric_specific_heat_capacity->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_isobaric_specific_heat_capacity =
        data_isobaric_specific_heat_capacity->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_mass_fractions =
        data_mass_fractions->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_isobaric_specific_heat_capacity;
        num_ghosts_min = hier::IntVector::min(num_ghosts_mass_fractions, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
        TBOX_ASSERT(data_isobaric_specific_heat_capacity->getGhostBox().contains(domain));
        TBOX_ASSERT(data_mass_fractions->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
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
            num_ghosts_isobaric_specific_heat_capacity,
            num_ghosts_mass_fractions,
            ghostcell_dims_isobaric_specific_heat_capacity,
            ghostcell_dims_mass_fractions,
            domain_lo,
            domain_dims);
    }
    else if (data_mass_fractions->getDepth() == d_num_species - 1)
    {
        boost::shared_ptr<pdat::CellData<double> > data_mass_fractions_last(
            new pdat::CellData<double>(interior_box, 1, num_ghosts_mass_fractions));
        
        if (domain.empty())
        {
            data_mass_fractions_last->fillAll(double(1));
        }
        else
        {
            data_mass_fractions_last->fillAll(double(1), domain);
        }
        
        /*
         * Get the pointers to the cell data of mass fractions.
         */
        
        std::vector<const double*> Y;
        Y.reserve(d_num_species - 1);
        for (int si = 0; si < d_num_species - 1; si++)
        {
            Y.push_back(data_mass_fractions->getPointer(si));
        }
        
        double* const Y_last = data_mass_fractions_last->getPointer(0);
        
        computeIsobaricSpecificHeatCapacity(
            c_p,
            Y_last,
            Y,
            num_ghosts_isobaric_specific_heat_capacity,
            num_ghosts_mass_fractions,
            ghostcell_dims_isobaric_specific_heat_capacity,
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
    boost::shared_ptr<pdat::SideData<double> >& data_isobaric_specific_heat_capacity,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
    int side_normal,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT((d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC) ||
                (d_mixing_closure_model == MIXING_CLOSURE_MODEL::NO_MODEL && d_num_species == 1));
    
    TBOX_ASSERT(data_isobaric_specific_heat_capacity);
    TBOX_ASSERT(data_mass_fractions);
    
    TBOX_ASSERT((data_mass_fractions->getDepth() == d_num_species) ||
                (data_mass_fractions->getDepth() == d_num_species - 1));
#endif
    
    NULL_USE(data_density);
    NULL_USE(data_pressure);
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_isobaric_specific_heat_capacity->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_mass_fractions->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_isobaric_specific_heat_capacity =
        data_isobaric_specific_heat_capacity->getGhostCellWidth();
    hier::IntVector ghostcell_dims_isobaric_specific_heat_capacity =
        data_isobaric_specific_heat_capacity->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
    hier::IntVector ghostcell_dims_mass_fractions =
        data_mass_fractions->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_isobaric_specific_heat_capacity;
        num_ghosts_min = hier::IntVector::min(num_ghosts_mass_fractions, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
        TBOX_ASSERT(data_isobaric_specific_heat_capacity->getGhostBox().contains(domain));
        TBOX_ASSERT(data_mass_fractions->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(side_normal < d_dim.getValue());
    
    TBOX_ASSERT(data_isobaric_specific_heat_capacity->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_mass_fractions->getDirectionVector()[side_normal] > 0);
#endif
    
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
            num_ghosts_isobaric_specific_heat_capacity,
            num_ghosts_mass_fractions,
            ghostcell_dims_isobaric_specific_heat_capacity,
            ghostcell_dims_mass_fractions,
            domain_lo,
            domain_dims);
    }
    else if (data_mass_fractions->getDepth() == d_num_species - 1)
    {
        hier::IntVector direction = hier::IntVector::getZero(d_dim);
        direction[side_normal] = 1;
        
        boost::shared_ptr<pdat::SideData<double> > data_mass_fractions_last(
            new pdat::SideData<double>(interior_box, 1, num_ghosts_mass_fractions,
                direction));
        
        if (domain.empty())
        {
            data_mass_fractions_last->fillAll(double(1));
        }
        else
        {
            data_mass_fractions_last->fillAll(double(1), domain);
        }
        
        /*
         * Get the pointers to the cell data of mass fractions.
         */
        
        std::vector<const double*> Y;
        Y.reserve(d_num_species - 1);
        for (int si = 0; si < d_num_species - 1; si++)
        {
            Y.push_back(data_mass_fractions->getPointer(side_normal, si));
        }
        
        double* const Y_last = data_mass_fractions_last->getPointer(side_normal, 0);
        
        computeIsobaricSpecificHeatCapacity(
            c_p,
            Y_last,
            Y,
            num_ghosts_isobaric_specific_heat_capacity,
            num_ghosts_mass_fractions,
            ghostcell_dims_isobaric_specific_heat_capacity,
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


///////////////////////////////////////////////////////////////


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
    
    NULL_USE(density);
    NULL_USE(pressure);
    
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
    
    return (gamma - double(1));
}


/*
 * Compute the Gruneisen parameter of the mixture with isothermal and isobaric equilibrium assumptions
 * (partial derivative of pressure w.r.t. specific internal energy under constant partial densities
 * divided by mixture density).
 */
void
EquationOfStateMixingRulesIdealGas::computeGruneisenParameter(
    boost::shared_ptr<pdat::CellData<double> >& data_gruneisen_parameter,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
    const hier::Box& domain) const
{
}


/*
 * Compute the Gruneisen parameter of the mixture with isothermal and isobaric equilibrium assumptions
 * (partial derivative of pressure w.r.t. specific internal energy under constant partial densities
 * divided by mixture density).
 */
void
EquationOfStateMixingRulesIdealGas::computeGruneisenParameter(
    boost::shared_ptr<pdat::SideData<double> >& data_gruneisen_parameter,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
    int side_normal,
    const hier::Box& domain) const
{
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
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOBARIC);
    TBOX_ASSERT((static_cast<int>(volume_fractions.size()) == d_num_species) ||
                (static_cast<int>(volume_fractions.size()) == d_num_species - 1));
#endif
    
    NULL_USE(density);
    NULL_USE(pressure);
    NULL_USE(mass_fractions);
    
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
    
    return (gamma - double(1));
}


/*
 * Compute the Gruneisen parameter of the mixture with isobaric equilibrium assumption
 * (partial derivative of pressure w.r.t. specific internal energy under constant partial densities
 * and volume fractions divided by mixture density).
 */
void
EquationOfStateMixingRulesIdealGas::computeGruneisenParameter(
    boost::shared_ptr<pdat::CellData<double> >& data_gruneisen_parameter,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
    const boost::shared_ptr<pdat::CellData<double> >& data_volume_fractions,
    const hier::Box& domain) const
{
}


/*
 * Compute the Gruneisen parameter of the mixture with isobaric equilibrium assumption
 * (partial derivative of pressure w.r.t. specific internal energy under constant partial densities
 * and volume fractions divided by mixture density).
 */
void
EquationOfStateMixingRulesIdealGas::computeGruneisenParameter(
    boost::shared_ptr<pdat::SideData<double> >& data_gruneisen_parameter,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
    const boost::shared_ptr<pdat::SideData<double> >& data_volume_fractions,
    int side_normal,
    const hier::Box& domain) const
{
}


///////////////////////////////////////////////////////////////


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
    boost::shared_ptr<pdat::CellData<double> >& data_partial_pressure_partial_densities,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
    const hier::Box& domain) const
{
}


/*
 * Compute the mixture partial derivative of pressure w.r.t. partial densities under constant specific
 * internal energy with isothermal and isobaric equilibrium assumptions.
 */
void
EquationOfStateMixingRulesIdealGas::computePressureDerivativeWithPartialDensities(
    boost::shared_ptr<pdat::SideData<double> >& data_partial_pressure_partial_densities,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
    int side_normal,
    const hier::Box& domain) const
{
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
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOBARIC);
    TBOX_ASSERT((static_cast<int>(volume_fractions.size()) == d_num_species) ||
                (static_cast<int>(volume_fractions.size()) == d_num_species - 1));
#endif
    
    NULL_USE(mass_fractions);
    
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
    
    const double& rho = *density;
    const double& p   = *pressure;
    
    const double epsilon = p/((gamma - double(1))*rho);
    
    double Psi_i = (gamma - double(1))*epsilon;
    
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
    boost::shared_ptr<pdat::CellData<double> >& data_partial_pressure_partial_densities,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
    const boost::shared_ptr<pdat::CellData<double> >& data_volume_fractions,
    const hier::Box& domain) const
{
}


/*
 * Compute the mixture partial derivative of pressure w.r.t. partial densities under constant specific
 * internal energy and volume fractions with isobaric equilibrium assumption.
 */
void
EquationOfStateMixingRulesIdealGas::computePressureDerivativeWithPartialDensities(
    boost::shared_ptr<pdat::SideData<double> >& data_partial_pressure_partial_densities,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
    const boost::shared_ptr<pdat::SideData<double> >& data_volume_fractions,
    int side_normal,
    const hier::Box& domain) const
{
}


///////////////////////////////////////////////////////////////


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
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOBARIC);
    TBOX_ASSERT((static_cast<int>(volume_fractions.size()) == d_num_species) ||
                (static_cast<int>(volume_fractions.size()) == d_num_species - 1));
#endif
    
    NULL_USE(mass_fractions);
    
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
    
    const double& rho = *density;
    const double& p   = *pressure;
    
    const double epsilon = p/((gamma - double(1))*rho);
    const double temp_1 = (gamma - double(1))*(gamma - double(1))*rho*epsilon;
    const double temp_2 = double(1)/(d_species_gamma[d_num_species - 1] - double(1));
    
    std::vector<double> M;
    M.reserve(d_num_species - 1);
    for (int si = 0; si < d_num_species - 1; si++)
    {
        double M_i = temp_1*(temp_2 - double(1)/(d_species_gamma[si] - double(1)));
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
    boost::shared_ptr<pdat::CellData<double> >& data_partial_pressure_partial_volume_fractions,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
    const boost::shared_ptr<pdat::CellData<double> >& data_volume_fractions,
    const hier::Box& domain) const
{
}


/*
 * Compute the mixture partial derivative of pressure w.r.t. volume fractions under constant specific
 * internal energy and partial densities with isobaric equilibrium assumption.
 */
void
EquationOfStateMixingRulesIdealGas::computePressureDerivativeWithVolumeFractions(
    boost::shared_ptr<pdat::SideData<double> >& data_partial_pressure_partial_volume_fractions,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
    const boost::shared_ptr<pdat::SideData<double> >& data_volume_fractions,
    int side_normal,
    const hier::Box& domain) const
{
}


///////////////////////////////////////////////////////////////


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
    boost::shared_ptr<pdat::CellData<double> >& data_mixture_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
    const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
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
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_mixture_density->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_pressure->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_temperature->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_mass_fractions->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells.
     */
    
    const hier::IntVector num_ghosts_mixture_density = data_mixture_density->getGhostCellWidth();
    const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
    const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
    const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
    
    /*
     * Get the minimum number of ghost cells for mixture thermodynamic properties.
     */
    
    hier::IntVector num_ghosts_min(d_dim);
    
    num_ghosts_min = num_ghosts_mixture_density;
    num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
    num_ghosts_min = hier::IntVector::min(num_ghosts_temperature, num_ghosts_min);
    num_ghosts_min = hier::IntVector::min(num_ghosts_mass_fractions, num_ghosts_min);
    
    /*
     * Get the mixture thermodyanmic properties.
     */
    
    const int num_thermo_properties = getNumberOfMixtureThermodynamicProperties();
    
    boost::shared_ptr<pdat::CellData<double> > data_mixture_thermo_properties(
        new pdat::CellData<double>(interior_box, num_thermo_properties, num_ghosts_min));
    
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
    boost::shared_ptr<pdat::SideData<double> >& data_mixture_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const boost::shared_ptr<pdat::SideData<double> >& data_temperature,
    const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
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
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_mixture_density->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_pressure->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_temperature->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_mass_fractions->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells.
     */
    
    const hier::IntVector num_ghosts_mixture_density = data_mixture_density->getGhostCellWidth();
    const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
    const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
    const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
    
    /*
     * Get the minimum number of ghost cells for mixture thermodynamic properties.
     */
    
    hier::IntVector num_ghosts_min(d_dim);
    
    num_ghosts_min = num_ghosts_mixture_density;
    num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
    num_ghosts_min = hier::IntVector::min(num_ghosts_temperature, num_ghosts_min);
    num_ghosts_min = hier::IntVector::min(num_ghosts_mass_fractions, num_ghosts_min);
    
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
    
    boost::shared_ptr<pdat::SideData<double> > data_mixture_thermo_properties(
        new pdat::SideData<double>(interior_box, num_thermo_properties, num_ghosts_min,
            direction));
    
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
            getMixtureThermodynamicPropertiesWithMassFraction(
                mixture_thermo_properties,
                species_fraction);
            
            break;
        }
        case MIXING_CLOSURE_MODEL::ISOBARIC:
        {
            getMixtureThermodynamicPropertiesWithVolumeFraction(
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
    boost::shared_ptr<pdat::CellData<double> >& data_mixture_thermo_properties,
    const boost::shared_ptr<pdat::CellData<double> >& data_species_fraction,
    const hier::Box& domain) const
{
    switch (d_mixing_closure_model)
    {
        case MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC:
        {
            computeMixtureThermodynamicPropertiesWithMassFraction(
                data_mixture_thermo_properties,
                data_species_fraction,
                domain);
            
            break;
        }
        case MIXING_CLOSURE_MODEL::ISOBARIC:
        {
            computeMixtureThermodynamicPropertiesWithVolumeFraction(
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
    boost::shared_ptr<pdat::SideData<double> >& data_mixture_thermo_properties,
    const boost::shared_ptr<pdat::SideData<double> >& data_species_fraction,
    int side_normal,
    const hier::Box& domain) const
{
    switch (d_mixing_closure_model)
    {
        case MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC:
        {
            computeMixtureThermodynamicPropertiesWithMassFraction(
                data_mixture_thermo_properties,
                data_species_fraction,
                side_normal,
                domain);
            
            break;
        }
        case MIXING_CLOSURE_MODEL::ISOBARIC:
        {
            computeMixtureThermodynamicPropertiesWithVolumeFraction(
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
EquationOfStateMixingRulesIdealGas::getMixtureThermodynamicPropertiesWithMassFraction(
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
EquationOfStateMixingRulesIdealGas::computeMixtureThermodynamicPropertiesWithMassFraction(
    boost::shared_ptr<pdat::CellData<double> >& data_mixture_thermo_properties,
    const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(data_mixture_thermo_properties);
    TBOX_ASSERT(data_mixture_thermo_properties->getDepth() == 4);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_mixture_thermo_properties->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(data_mass_fractions->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_mixture_thermo_properties = data_mixture_thermo_properties->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_mixture_thermo_properties =
        data_mixture_thermo_properties->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_mass_fractions =
        data_mass_fractions->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_mixture_thermo_properties;
        num_ghosts_min = hier::IntVector::min(num_ghosts_mass_fractions, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
        TBOX_ASSERT(data_mixture_thermo_properties->getGhostBox().contains(domain));
        TBOX_ASSERT(data_mass_fractions->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
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
        
        computeMixtureThermodynamicPropertiesWithMassFraction(
            gamma,
            R,
            c_p,
            c_v,
            Y,
            num_ghosts_mixture_thermo_properties,
            num_ghosts_mass_fractions,
            ghostcell_dims_mixture_thermo_properties,
            ghostcell_dims_mass_fractions,
            domain_lo,
            domain_dims);
    }
    else if (data_mass_fractions->getDepth() == d_num_species - 1)
    {
        boost::shared_ptr<pdat::CellData<double> > data_mass_fractions_last(
            new pdat::CellData<double>(interior_box, 1, num_ghosts_mass_fractions));
        
        if (domain.empty())
        {
            data_mass_fractions_last->fillAll(double(1));
        }
        else
        {
            data_mass_fractions_last->fillAll(double(1), domain);
        }
        
        /*
         * Get the pointers to the cell data of mass fractions.
         */
        
        std::vector<const double*> Y;
        Y.reserve(d_num_species - 1);
        for (int si = 0; si < d_num_species - 1; si++)
        {
            Y.push_back(data_mass_fractions->getPointer(si));
        }
        
        double* const Y_last = data_mass_fractions_last->getPointer(0);
        
        computeMixtureThermodynamicPropertiesWithMassFraction(
            gamma,
            R,
            c_p,
            c_v,
            Y_last,
            Y,
            num_ghosts_mixture_thermo_properties,
            num_ghosts_mass_fractions,
            ghostcell_dims_mixture_thermo_properties,
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
EquationOfStateMixingRulesIdealGas::computeMixtureThermodynamicPropertiesWithMassFraction(
    boost::shared_ptr<pdat::SideData<double> >& data_mixture_thermo_properties,
    const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
    int side_normal,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(data_mixture_thermo_properties);
    TBOX_ASSERT(data_mixture_thermo_properties->getDepth() == 4);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_mixture_thermo_properties->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(data_mass_fractions->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_mixture_thermo_properties = data_mixture_thermo_properties->getGhostCellWidth();
    hier::IntVector ghostcell_dims_mixture_thermo_properties =
        data_mixture_thermo_properties->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
    hier::IntVector ghostcell_dims_mass_fractions =
        data_mass_fractions->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_mixture_thermo_properties;
        num_ghosts_min = hier::IntVector::min(num_ghosts_mass_fractions, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
        TBOX_ASSERT(data_mixture_thermo_properties->getGhostBox().contains(domain));
        TBOX_ASSERT(data_mass_fractions->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(side_normal < d_dim.getValue());
    
    TBOX_ASSERT(data_mixture_thermo_properties->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_mass_fractions->getDirectionVector()[side_normal] > 0);
#endif
    
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
        
        computeMixtureThermodynamicPropertiesWithMassFraction(
            gamma,
            R,
            c_p,
            c_v,
            Y,
            num_ghosts_mixture_thermo_properties,
            num_ghosts_mass_fractions,
            ghostcell_dims_mixture_thermo_properties,
            ghostcell_dims_mass_fractions,
            domain_lo,
            domain_dims);
    }
    else if (data_mass_fractions->getDepth() == d_num_species - 1)
    {
        hier::IntVector direction = hier::IntVector::getZero(d_dim);
        direction[side_normal] = 1;
        
        boost::shared_ptr<pdat::SideData<double> > data_mass_fractions_last(
            new pdat::SideData<double>(interior_box, 1, num_ghosts_mass_fractions,
                direction));
        
        if (domain.empty())
        {
            data_mass_fractions_last->fillAll(double(1));
        }
        else
        {
            data_mass_fractions_last->fillAll(double(1), domain);
        }
        
        /*
         * Get the pointers to the cell data of mass fractions.
         */
        
        std::vector<const double*> Y;
        Y.reserve(d_num_species - 1);
        for (int si = 0; si < d_num_species - 1; si++)
        {
            Y.push_back(data_mass_fractions->getPointer(side_normal, si));
        }
        
        double* const Y_last = data_mass_fractions_last->getPointer(side_normal, 0);
        
        computeMixtureThermodynamicPropertiesWithMassFraction(
            gamma,
            R,
            c_p,
            c_v,
            Y_last,
            Y,
            num_ghosts_mixture_thermo_properties,
            num_ghosts_mass_fractions,
            ghostcell_dims_mixture_thermo_properties,
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
EquationOfStateMixingRulesIdealGas::getMixtureThermodynamicPropertiesWithVolumeFraction(
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
EquationOfStateMixingRulesIdealGas::computeMixtureThermodynamicPropertiesWithVolumeFraction(
    boost::shared_ptr<pdat::CellData<double> >& data_mixture_thermo_properties,
    const boost::shared_ptr<pdat::CellData<double> >& data_volume_fractions,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(data_mixture_thermo_properties);
    TBOX_ASSERT(data_mixture_thermo_properties->getDepth() == 1);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_mixture_thermo_properties->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(data_volume_fractions->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_mixture_thermo_properties = data_mixture_thermo_properties->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_mixture_thermo_properties =
        data_mixture_thermo_properties->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_volume_fractions = data_volume_fractions->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_volume_fractions =
        data_volume_fractions->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_mixture_thermo_properties;
        num_ghosts_min = hier::IntVector::min(num_ghosts_volume_fractions, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
        TBOX_ASSERT(data_mixture_thermo_properties->getGhostBox().contains(domain));
        TBOX_ASSERT(data_volume_fractions->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
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
        
        getMixtureThermodynamicPropertiesWithVolumeFraction(
            gamma,
            Z,
            num_ghosts_mixture_thermo_properties,
            num_ghosts_volume_fractions,
            ghostcell_dims_mixture_thermo_properties,
            ghostcell_dims_volume_fractions,
            domain_lo,
            domain_dims);
    }
    else if (data_volume_fractions->getDepth() == d_num_species - 1)
    {
        boost::shared_ptr<pdat::CellData<double> > data_volume_fractions_last(
            new pdat::CellData<double>(interior_box, 1, num_ghosts_volume_fractions));
        
        if (domain.empty())
        {
            data_volume_fractions_last->fillAll(double(1));
        }
        else
        {
            data_volume_fractions_last->fillAll(double(1), domain);
        }
        
        /*
         * Get the pointers to the cell data of volume fractions.
         */
        
        std::vector<const double*> Z;
        Z.reserve(d_num_species - 1);
        for (int si = 0; si < d_num_species - 1; si++)
        {
            Z.push_back(data_volume_fractions->getPointer(si));
        }
        
        double* const Z_last = data_volume_fractions_last->getPointer(0);
        
        getMixtureThermodynamicPropertiesWithVolumeFraction(
            gamma,
            Z_last,
            Z,
            num_ghosts_mixture_thermo_properties,
            num_ghosts_volume_fractions,
            ghostcell_dims_mixture_thermo_properties,
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
EquationOfStateMixingRulesIdealGas::computeMixtureThermodynamicPropertiesWithVolumeFraction(
    boost::shared_ptr<pdat::SideData<double> >& data_mixture_thermo_properties,
    const boost::shared_ptr<pdat::SideData<double> >& data_volume_fractions,
    int side_normal,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(data_mixture_thermo_properties);
    TBOX_ASSERT(data_mixture_thermo_properties->getDepth() == 1);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_mixture_thermo_properties->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(data_volume_fractions->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_mixture_thermo_properties = data_mixture_thermo_properties->getGhostCellWidth();
    hier::IntVector ghostcell_dims_mixture_thermo_properties =
        data_mixture_thermo_properties->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_volume_fractions = data_volume_fractions->getGhostCellWidth();
    hier::IntVector ghostcell_dims_volume_fractions =
        data_volume_fractions->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_mixture_thermo_properties;
        num_ghosts_min = hier::IntVector::min(num_ghosts_volume_fractions, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
        TBOX_ASSERT(data_mixture_thermo_properties->getGhostBox().contains(domain));
        TBOX_ASSERT(data_volume_fractions->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(side_normal < d_dim.getValue());
    
    TBOX_ASSERT(data_mixture_thermo_properties->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_volume_fractions->getDirectionVector()[side_normal] > 0);
#endif
    
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
        
        getMixtureThermodynamicPropertiesWithVolumeFraction(
            gamma,
            Z,
            num_ghosts_mixture_thermo_properties,
            num_ghosts_volume_fractions,
            ghostcell_dims_mixture_thermo_properties,
            ghostcell_dims_volume_fractions,
            domain_lo,
            domain_dims);
    }
    else if (data_volume_fractions->getDepth() == d_num_species - 1)
    {
        hier::IntVector direction = hier::IntVector::getZero(d_dim);
        direction[side_normal] = 1;
        
        boost::shared_ptr<pdat::SideData<double> > data_volume_fractions_last(
            new pdat::SideData<double>(interior_box, 1, num_ghosts_volume_fractions,
                direction));
        
        if (domain.empty())
        {
            data_volume_fractions_last->fillAll(double(1));
        }
        else
        {
            data_volume_fractions_last->fillAll(double(1), domain);
        }
        
        /*
         * Get the pointers to the cell data of volume fractions.
         */
        
        std::vector<const double*> Z;
        Z.reserve(d_num_species - 1);
        for (int si = 0; si < d_num_species - 1; si++)
        {
            Z.push_back(data_volume_fractions->getPointer(side_normal, si));
        }
        
        double* const Z_last = data_volume_fractions_last->getPointer(side_normal, 0);
        
        getMixtureThermodynamicPropertiesWithVolumeFraction(
            gamma,
            Z_last,
            Z,
            num_ghosts_mixture_thermo_properties,
            num_ghosts_volume_fractions,
            ghostcell_dims_mixture_thermo_properties,
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
    const std::vector<const double*> Y,
    const hier::IntVector& num_ghosts_isochoric_specific_heat_capacity,
    const hier::IntVector& num_ghosts_mass_fractions,
    const hier::IntVector& ghostcell_dims_isochoric_specific_heat_capacity,
    const hier::IntVector& ghostcell_dims_mass_fractions,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims) const
{
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and numbers of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int num_ghosts_0_isochoric_specific_heat_capacity =
            num_ghosts_isochoric_specific_heat_capacity[0];
        const int num_ghosts_0_mass_fractions = num_ghosts_mass_fractions[0];
        
        // Compute c_v.
        for (int si = 0; si < d_num_species; si++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_isochoric_specific_heat_capacity = i + num_ghosts_0_isochoric_specific_heat_capacity;
                const int idx_mass_fractions = i + num_ghosts_0_mass_fractions;
                
                c_v[idx_isochoric_specific_heat_capacity] += Y[si][idx_mass_fractions]*d_species_c_v[si];
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the local lower indices, numbers of cells in each dimension and numbers of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        
        const int num_ghosts_0_isochoric_specific_heat_capacity = num_ghosts_isochoric_specific_heat_capacity[0];
        const int num_ghosts_1_isochoric_specific_heat_capacity = num_ghosts_isochoric_specific_heat_capacity[1];
        const int ghostcell_dim_0_isochoric_specific_heat_capacity =
            ghostcell_dims_isochoric_specific_heat_capacity[0];
        
        const int num_ghosts_0_mass_fractions = num_ghosts_mass_fractions[0];
        const int num_ghosts_1_mass_fractions = num_ghosts_mass_fractions[1];
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
                        (i + num_ghosts_0_isochoric_specific_heat_capacity) +
                        (j + num_ghosts_1_isochoric_specific_heat_capacity)*
                            ghostcell_dim_0_isochoric_specific_heat_capacity;
                    
                    const int idx_mass_fractions = (i + num_ghosts_0_mass_fractions) +
                        (j + num_ghosts_1_mass_fractions)*ghostcell_dim_0_mass_fractions;
                    
                    c_v[idx_isochoric_specific_heat_capacity] += Y[si][idx_mass_fractions]*d_species_c_v[si];
                }
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the local lower indices, numbers of cells in each dimension and numbers of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_lo_2 = domain_lo[2];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        const int domain_dim_2 = domain_dims[2];
        
        const int num_ghosts_0_isochoric_specific_heat_capacity = num_ghosts_isochoric_specific_heat_capacity[0];
        const int num_ghosts_1_isochoric_specific_heat_capacity = num_ghosts_isochoric_specific_heat_capacity[1];
        const int num_ghosts_2_isochoric_specific_heat_capacity = num_ghosts_isochoric_specific_heat_capacity[2];
        const int ghostcell_dim_0_isochoric_specific_heat_capacity =
            ghostcell_dims_isochoric_specific_heat_capacity[0];
        const int ghostcell_dim_1_isochoric_specific_heat_capacity =
            ghostcell_dims_isochoric_specific_heat_capacity[1];
        
        const int num_ghosts_0_mass_fractions = num_ghosts_mass_fractions[0];
        const int num_ghosts_1_mass_fractions = num_ghosts_mass_fractions[1];
        const int num_ghosts_2_mass_fractions = num_ghosts_mass_fractions[2];
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
                            (i + num_ghosts_0_isochoric_specific_heat_capacity) +
                            (j + num_ghosts_1_isochoric_specific_heat_capacity)*
                                ghostcell_dim_0_isochoric_specific_heat_capacity +
                            (k + num_ghosts_2_isochoric_specific_heat_capacity)*
                                ghostcell_dim_0_isochoric_specific_heat_capacity*
                                    ghostcell_dim_1_isochoric_specific_heat_capacity;
                        
                        const int idx_mass_fractions = (i + num_ghosts_0_mass_fractions) +
                            (j + num_ghosts_1_mass_fractions)*ghostcell_dim_0_mass_fractions +
                            (k + num_ghosts_2_mass_fractions)*ghostcell_dim_0_mass_fractions*
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
    const std::vector<const double*> Y,
    const hier::IntVector& num_ghosts_isochoric_specific_heat_capacity,
    const hier::IntVector& num_ghosts_mass_fractions,
    const hier::IntVector& ghostcell_dims_isochoric_specific_heat_capacity,
    const hier::IntVector& ghostcell_dims_mass_fractions,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims) const
{
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and numbers of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int num_ghosts_0_isochoric_specific_heat_capacity =
            num_ghosts_isochoric_specific_heat_capacity[0];
        const int num_ghosts_0_mass_fractions = num_ghosts_mass_fractions[0];
        
        // Compute c_v.
        for (int si = 0; si < d_num_species - 1; si++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_isochoric_specific_heat_capacity = i + num_ghosts_0_isochoric_specific_heat_capacity;
                const int idx_mass_fractions = i + num_ghosts_0_mass_fractions;
                
                c_v[idx_isochoric_specific_heat_capacity] += Y[si][idx_mass_fractions]*d_species_c_v[si];
                
                // Compute the mass fraction of the last species.
                Y_last[idx_mass_fractions] -= Y[si][idx_mass_fractions];
            }
        }
        
        // Add the contribution from the last species.
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_isochoric_specific_heat_capacity = i + num_ghosts_0_isochoric_specific_heat_capacity;
            const int idx_mass_fractions = i + num_ghosts_0_mass_fractions;
            
            c_v[idx_isochoric_specific_heat_capacity] += Y_last[idx_mass_fractions]*d_species_c_v.back();
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the local lower indices, numbers of cells in each dimension and numbers of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        
        const int num_ghosts_0_isochoric_specific_heat_capacity = num_ghosts_isochoric_specific_heat_capacity[0];
        const int num_ghosts_1_isochoric_specific_heat_capacity = num_ghosts_isochoric_specific_heat_capacity[1];
        const int ghostcell_dim_0_isochoric_specific_heat_capacity =
            ghostcell_dims_isochoric_specific_heat_capacity[0];
        
        const int num_ghosts_0_mass_fractions = num_ghosts_mass_fractions[0];
        const int num_ghosts_1_mass_fractions = num_ghosts_mass_fractions[1];
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
                        (i + num_ghosts_0_isochoric_specific_heat_capacity) +
                        (j + num_ghosts_1_isochoric_specific_heat_capacity)*
                            ghostcell_dim_0_isochoric_specific_heat_capacity;
                    
                    const int idx_mass_fractions = (i + num_ghosts_0_mass_fractions) +
                        (j + num_ghosts_1_mass_fractions)*ghostcell_dim_0_mass_fractions;
                    
                    c_v[idx_isochoric_specific_heat_capacity] += Y[si][idx_mass_fractions]*d_species_c_v[si];
                    
                    // Compute the mass fraction of the last species.
                    Y_last[idx_mass_fractions] -= Y[si][idx_mass_fractions];
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
                    (i + num_ghosts_0_isochoric_specific_heat_capacity) +
                    (j + num_ghosts_1_isochoric_specific_heat_capacity)*
                        ghostcell_dim_0_isochoric_specific_heat_capacity;
                
                const int idx_mass_fractions = (i + num_ghosts_0_mass_fractions) +
                    (j + num_ghosts_1_mass_fractions)*ghostcell_dim_0_mass_fractions;
                
                c_v[idx_isochoric_specific_heat_capacity] += Y_last[idx_mass_fractions]*d_species_c_v.back();
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the local lower indices, numbers of cells in each dimension and numbers of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_lo_2 = domain_lo[2];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        const int domain_dim_2 = domain_dims[2];
        
        const int num_ghosts_0_isochoric_specific_heat_capacity = num_ghosts_isochoric_specific_heat_capacity[0];
        const int num_ghosts_1_isochoric_specific_heat_capacity = num_ghosts_isochoric_specific_heat_capacity[1];
        const int num_ghosts_2_isochoric_specific_heat_capacity = num_ghosts_isochoric_specific_heat_capacity[2];
        const int ghostcell_dim_0_isochoric_specific_heat_capacity =
            ghostcell_dims_isochoric_specific_heat_capacity[0];
        const int ghostcell_dim_1_isochoric_specific_heat_capacity =
            ghostcell_dims_isochoric_specific_heat_capacity[1];
        
        const int num_ghosts_0_mass_fractions = num_ghosts_mass_fractions[0];
        const int num_ghosts_1_mass_fractions = num_ghosts_mass_fractions[1];
        const int num_ghosts_2_mass_fractions = num_ghosts_mass_fractions[2];
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
                            (i + num_ghosts_0_isochoric_specific_heat_capacity) +
                            (j + num_ghosts_1_isochoric_specific_heat_capacity)*
                                ghostcell_dim_0_isochoric_specific_heat_capacity +
                            (k + num_ghosts_2_isochoric_specific_heat_capacity)*
                                ghostcell_dim_0_isochoric_specific_heat_capacity*
                                    ghostcell_dim_1_isochoric_specific_heat_capacity;
                        
                        const int idx_mass_fractions = (i + num_ghosts_0_mass_fractions) +
                            (j + num_ghosts_1_mass_fractions)*ghostcell_dim_0_mass_fractions +
                            (k + num_ghosts_2_mass_fractions)*ghostcell_dim_0_mass_fractions*
                                ghostcell_dim_1_mass_fractions;
                        
                        c_v[idx_isochoric_specific_heat_capacity] += Y[si][idx_mass_fractions]*d_species_c_v[si];
                        
                        // Compute the mass fraction of the last species.
                        Y_last[idx_mass_fractions] -= Y[si][idx_mass_fractions];
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
                        (i + num_ghosts_0_isochoric_specific_heat_capacity) +
                        (j + num_ghosts_1_isochoric_specific_heat_capacity)*
                            ghostcell_dim_0_isochoric_specific_heat_capacity +
                        (k + num_ghosts_2_isochoric_specific_heat_capacity)*
                            ghostcell_dim_0_isochoric_specific_heat_capacity*
                                ghostcell_dim_1_isochoric_specific_heat_capacity;
                    
                    const int idx_mass_fractions = (i + num_ghosts_0_mass_fractions) +
                        (j + num_ghosts_1_mass_fractions)*ghostcell_dim_0_mass_fractions +
                        (k + num_ghosts_2_mass_fractions)*ghostcell_dim_0_mass_fractions*
                            ghostcell_dim_1_mass_fractions;
                        
                    c_v[idx_isochoric_specific_heat_capacity] += Y_last[idx_mass_fractions]*d_species_c_v.back();
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
    const std::vector<const double*> Y,
    const hier::IntVector& num_ghosts_isobaric_specific_heat_capacity,
    const hier::IntVector& num_ghosts_mass_fractions,
    const hier::IntVector& ghostcell_dims_isobaric_specific_heat_capacity,
    const hier::IntVector& ghostcell_dims_mass_fractions,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims) const
{
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and numbers of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int num_ghosts_0_isobaric_specific_heat_capacity =
            num_ghosts_isobaric_specific_heat_capacity[0];
        const int num_ghosts_0_mass_fractions = num_ghosts_mass_fractions[0];
        
        // Compute c_p.
        for (int si = 0; si < d_num_species; si++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_isobaric_specific_heat_capacity = i + num_ghosts_0_isobaric_specific_heat_capacity;
                const int idx_mass_fractions = i + num_ghosts_0_mass_fractions;
                
                c_p[idx_isobaric_specific_heat_capacity] += Y[si][idx_mass_fractions]*d_species_c_p[si];
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the local lower indices, numbers of cells in each dimension and numbers of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        
        const int num_ghosts_0_isobaric_specific_heat_capacity = num_ghosts_isobaric_specific_heat_capacity[0];
        const int num_ghosts_1_isobaric_specific_heat_capacity = num_ghosts_isobaric_specific_heat_capacity[1];
        const int ghostcell_dim_0_isobaric_specific_heat_capacity =
            ghostcell_dims_isobaric_specific_heat_capacity[0];
        
        const int num_ghosts_0_mass_fractions = num_ghosts_mass_fractions[0];
        const int num_ghosts_1_mass_fractions = num_ghosts_mass_fractions[1];
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
                        (i + num_ghosts_0_isobaric_specific_heat_capacity) +
                        (j + num_ghosts_1_isobaric_specific_heat_capacity)*
                            ghostcell_dim_0_isobaric_specific_heat_capacity;
                    
                    const int idx_mass_fractions = (i + num_ghosts_0_mass_fractions) +
                        (j + num_ghosts_1_mass_fractions)*ghostcell_dim_0_mass_fractions;
                    
                    c_p[idx_isobaric_specific_heat_capacity] += Y[si][idx_mass_fractions]*d_species_c_p[si];
                }
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the local lower indices, numbers of cells in each dimension and numbers of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_lo_2 = domain_lo[2];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        const int domain_dim_2 = domain_dims[2];
        
        const int num_ghosts_0_isobaric_specific_heat_capacity = num_ghosts_isobaric_specific_heat_capacity[0];
        const int num_ghosts_1_isobaric_specific_heat_capacity = num_ghosts_isobaric_specific_heat_capacity[1];
        const int num_ghosts_2_isobaric_specific_heat_capacity = num_ghosts_isobaric_specific_heat_capacity[2];
        const int ghostcell_dim_0_isobaric_specific_heat_capacity =
            ghostcell_dims_isobaric_specific_heat_capacity[0];
        const int ghostcell_dim_1_isobaric_specific_heat_capacity =
            ghostcell_dims_isobaric_specific_heat_capacity[1];
        
        const int num_ghosts_0_mass_fractions = num_ghosts_mass_fractions[0];
        const int num_ghosts_1_mass_fractions = num_ghosts_mass_fractions[1];
        const int num_ghosts_2_mass_fractions = num_ghosts_mass_fractions[2];
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
                            (i + num_ghosts_0_isobaric_specific_heat_capacity) +
                            (j + num_ghosts_1_isobaric_specific_heat_capacity)*
                                ghostcell_dim_0_isobaric_specific_heat_capacity +
                            (k + num_ghosts_2_isobaric_specific_heat_capacity)*
                                ghostcell_dim_0_isobaric_specific_heat_capacity*
                                    ghostcell_dim_1_isobaric_specific_heat_capacity;
                        
                        const int idx_mass_fractions = (i + num_ghosts_0_mass_fractions) +
                            (j + num_ghosts_1_mass_fractions)*ghostcell_dim_0_mass_fractions +
                            (k + num_ghosts_2_mass_fractions)*ghostcell_dim_0_mass_fractions*
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
    const std::vector<const double*> Y,
    const hier::IntVector& num_ghosts_isobaric_specific_heat_capacity,
    const hier::IntVector& num_ghosts_mass_fractions,
    const hier::IntVector& ghostcell_dims_isobaric_specific_heat_capacity,
    const hier::IntVector& ghostcell_dims_mass_fractions,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims) const
{
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and numbers of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int num_ghosts_0_isobaric_specific_heat_capacity =
            num_ghosts_isobaric_specific_heat_capacity[0];
        const int num_ghosts_0_mass_fractions = num_ghosts_mass_fractions[0];
        
        // Compute c_p.
        for (int si = 0; si < d_num_species - 1; si++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_isobaric_specific_heat_capacity = i + num_ghosts_0_isobaric_specific_heat_capacity;
                const int idx_mass_fractions = i + num_ghosts_0_mass_fractions;
                
                c_p[idx_isobaric_specific_heat_capacity] += Y[si][idx_mass_fractions]*d_species_c_p[si];
                
                // Compute the mass fraction of the last species.
                Y_last[idx_mass_fractions] -= Y[si][idx_mass_fractions];
            }
        }
        
        // Add the contribution from the last species.
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_isobaric_specific_heat_capacity = i + num_ghosts_0_isobaric_specific_heat_capacity;
            const int idx_mass_fractions = i + num_ghosts_0_mass_fractions;
            
            c_p[idx_isobaric_specific_heat_capacity] += Y_last[idx_mass_fractions]*d_species_c_p.back();
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the local lower indices, numbers of cells in each dimension and numbers of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        
        const int num_ghosts_0_isobaric_specific_heat_capacity = num_ghosts_isobaric_specific_heat_capacity[0];
        const int num_ghosts_1_isobaric_specific_heat_capacity = num_ghosts_isobaric_specific_heat_capacity[1];
        const int ghostcell_dim_0_isobaric_specific_heat_capacity =
            ghostcell_dims_isobaric_specific_heat_capacity[0];
        
        const int num_ghosts_0_mass_fractions = num_ghosts_mass_fractions[0];
        const int num_ghosts_1_mass_fractions = num_ghosts_mass_fractions[1];
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
                        (i + num_ghosts_0_isobaric_specific_heat_capacity) +
                        (j + num_ghosts_1_isobaric_specific_heat_capacity)*
                            ghostcell_dim_0_isobaric_specific_heat_capacity;
                    
                    const int idx_mass_fractions = (i + num_ghosts_0_mass_fractions) +
                        (j + num_ghosts_1_mass_fractions)*ghostcell_dim_0_mass_fractions;
                    
                    c_p[idx_isobaric_specific_heat_capacity] += Y[si][idx_mass_fractions]*d_species_c_p[si];
                    
                    // Compute the mass fraction of the last species.
                    Y_last[idx_mass_fractions] -= Y[si][idx_mass_fractions];
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
                    (i + num_ghosts_0_isobaric_specific_heat_capacity) +
                    (j + num_ghosts_1_isobaric_specific_heat_capacity)*
                        ghostcell_dim_0_isobaric_specific_heat_capacity;
                
                const int idx_mass_fractions = (i + num_ghosts_0_mass_fractions) +
                    (j + num_ghosts_1_mass_fractions)*ghostcell_dim_0_mass_fractions;
                
                c_p[idx_isobaric_specific_heat_capacity] += Y_last[idx_mass_fractions]*d_species_c_p.back();
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the local lower indices, numbers of cells in each dimension and numbers of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_lo_2 = domain_lo[2];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        const int domain_dim_2 = domain_dims[2];
        
        const int num_ghosts_0_isobaric_specific_heat_capacity = num_ghosts_isobaric_specific_heat_capacity[0];
        const int num_ghosts_1_isobaric_specific_heat_capacity = num_ghosts_isobaric_specific_heat_capacity[1];
        const int num_ghosts_2_isobaric_specific_heat_capacity = num_ghosts_isobaric_specific_heat_capacity[2];
        const int ghostcell_dim_0_isobaric_specific_heat_capacity =
            ghostcell_dims_isobaric_specific_heat_capacity[0];
        const int ghostcell_dim_1_isobaric_specific_heat_capacity =
            ghostcell_dims_isobaric_specific_heat_capacity[1];
        
        const int num_ghosts_0_mass_fractions = num_ghosts_mass_fractions[0];
        const int num_ghosts_1_mass_fractions = num_ghosts_mass_fractions[1];
        const int num_ghosts_2_mass_fractions = num_ghosts_mass_fractions[2];
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
                            (i + num_ghosts_0_isobaric_specific_heat_capacity) +
                            (j + num_ghosts_1_isobaric_specific_heat_capacity)*
                                ghostcell_dim_0_isobaric_specific_heat_capacity +
                            (k + num_ghosts_2_isobaric_specific_heat_capacity)*
                                ghostcell_dim_0_isobaric_specific_heat_capacity*
                                    ghostcell_dim_1_isobaric_specific_heat_capacity;
                        
                        const int idx_mass_fractions = (i + num_ghosts_0_mass_fractions) +
                            (j + num_ghosts_1_mass_fractions)*ghostcell_dim_0_mass_fractions +
                            (k + num_ghosts_2_mass_fractions)*ghostcell_dim_0_mass_fractions*
                                ghostcell_dim_1_mass_fractions;
                        
                        c_p[idx_isobaric_specific_heat_capacity] += Y[si][idx_mass_fractions]*d_species_c_p[si];
                        
                        // Compute the mass fraction of the last species.
                        Y_last[idx_mass_fractions] -= Y[si][idx_mass_fractions];
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
                        (i + num_ghosts_0_isobaric_specific_heat_capacity) +
                        (j + num_ghosts_1_isobaric_specific_heat_capacity)*
                            ghostcell_dim_0_isobaric_specific_heat_capacity +
                        (k + num_ghosts_2_isobaric_specific_heat_capacity)*
                            ghostcell_dim_0_isobaric_specific_heat_capacity*
                                ghostcell_dim_1_isobaric_specific_heat_capacity;
                    
                    const int idx_mass_fractions = (i + num_ghosts_0_mass_fractions) +
                        (j + num_ghosts_1_mass_fractions)*ghostcell_dim_0_mass_fractions +
                        (k + num_ghosts_2_mass_fractions)*ghostcell_dim_0_mass_fractions*
                            ghostcell_dim_1_mass_fractions;
                        
                    c_p[idx_isobaric_specific_heat_capacity] += Y_last[idx_mass_fractions]*d_species_c_p.back();
                }
            }
        }
    }
}


/*
 * Compute the thermodynamic properties of the mixture with mass fractions.
 */
void
EquationOfStateMixingRulesIdealGas::computeMixtureThermodynamicPropertiesWithMassFraction(
    double* const gamma,
    double* const R,
    double* const c_p,
    double* const c_v,
    const std::vector<const double*> Y,
    const hier::IntVector& num_ghosts_mixture_thermo_properties,
    const hier::IntVector& num_ghosts_mass_fractions,
    const hier::IntVector& ghostcell_dims_mixture_thermo_properties,
    const hier::IntVector& ghostcell_dims_mass_fractions,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims) const
{
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and numbers of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int num_ghosts_0_mixture_thermo_properties = num_ghosts_mixture_thermo_properties[0];
        const int num_ghosts_0_mass_fractions = num_ghosts_mass_fractions[0];
        
        // Compute c_p and c_v.
        for (int si = 0; si < d_num_species; si++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_mixture_thermo_properties = i + num_ghosts_0_mixture_thermo_properties;
                const int idx_mass_fractions = i + num_ghosts_0_mass_fractions;
                
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
            const int idx_mixture_thermo_properties = i + num_ghosts_0_mixture_thermo_properties;
            
            gamma[idx_mixture_thermo_properties] = c_p[idx_mixture_thermo_properties]/
                c_v[idx_mixture_thermo_properties];
            R[idx_mixture_thermo_properties] = c_p[idx_mixture_thermo_properties] -
                c_v[idx_mixture_thermo_properties];
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the local lower indices, numbers of cells in each dimension and numbers of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        
        const int num_ghosts_0_mixture_thermo_properties = num_ghosts_mixture_thermo_properties[0];
        const int num_ghosts_1_mixture_thermo_properties = num_ghosts_mixture_thermo_properties[1];
        const int ghostcell_dim_0_mixture_thermo_properties = ghostcell_dims_mixture_thermo_properties[0];
        
        const int num_ghosts_0_mass_fractions = num_ghosts_mass_fractions[0];
        const int num_ghosts_1_mass_fractions = num_ghosts_mass_fractions[1];
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
                    const int idx_mixture_thermo_properties = (i + num_ghosts_0_mixture_thermo_properties) +
                        (j + num_ghosts_1_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties;
                    
                    const int idx_mass_fractions = (i + num_ghosts_0_mass_fractions) +
                        (j + num_ghosts_1_mass_fractions)*ghostcell_dim_0_mass_fractions;
                    
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
                const int idx_mixture_thermo_properties = (i + num_ghosts_0_mixture_thermo_properties) +
                    (j + num_ghosts_1_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties;
                
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
         * Get the local lower indices, numbers of cells in each dimension and numbers of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_lo_2 = domain_lo[2];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        const int domain_dim_2 = domain_dims[2];
        
        const int num_ghosts_0_mixture_thermo_properties = num_ghosts_mixture_thermo_properties[0];
        const int num_ghosts_1_mixture_thermo_properties = num_ghosts_mixture_thermo_properties[1];
        const int num_ghosts_2_mixture_thermo_properties = num_ghosts_mixture_thermo_properties[2];
        const int ghostcell_dim_0_mixture_thermo_properties = ghostcell_dims_mixture_thermo_properties[0];
        const int ghostcell_dim_1_mixture_thermo_properties = ghostcell_dims_mixture_thermo_properties[1];
        
        const int num_ghosts_0_mass_fractions = num_ghosts_mass_fractions[0];
        const int num_ghosts_1_mass_fractions = num_ghosts_mass_fractions[1];
        const int num_ghosts_2_mass_fractions = num_ghosts_mass_fractions[2];
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
                        const int idx_mixture_thermo_properties = (i + num_ghosts_0_mixture_thermo_properties) +
                            (j + num_ghosts_1_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties +
                            (k + num_ghosts_2_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties*
                                ghostcell_dim_1_mixture_thermo_properties;
                        
                        const int idx_mass_fractions = (i + num_ghosts_0_mass_fractions) +
                            (j + num_ghosts_1_mass_fractions)*ghostcell_dim_0_mass_fractions +
                            (k + num_ghosts_2_mass_fractions)*ghostcell_dim_0_mass_fractions*
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
                    const int idx_mixture_thermo_properties = (i + num_ghosts_0_mixture_thermo_properties) +
                        (j + num_ghosts_1_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties +
                        (k + num_ghosts_2_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties*
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
EquationOfStateMixingRulesIdealGas::computeMixtureThermodynamicPropertiesWithMassFraction(
    double* const gamma,
    double* const R,
    double* const c_p,
    double* const c_v,
    double* const Y_last,
    const std::vector<const double*> Y,
    const hier::IntVector& num_ghosts_mixture_thermo_properties,
    const hier::IntVector& num_ghosts_mass_fractions,
    const hier::IntVector& ghostcell_dims_mixture_thermo_properties,
    const hier::IntVector& ghostcell_dims_mass_fractions,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims) const
{
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and numbers of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int num_ghosts_0_mixture_thermo_properties = num_ghosts_mixture_thermo_properties[0];
        const int num_ghosts_0_mass_fractions = num_ghosts_mass_fractions[0];
        
        // Compute c_p and c_v.
        for (int si = 0; si < d_num_species - 1; si++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_mixture_thermo_properties = i + num_ghosts_0_mixture_thermo_properties;
                const int idx_mass_fractions = i + num_ghosts_0_mass_fractions;
                
                c_p[idx_mixture_thermo_properties] += Y[si][idx_mass_fractions]*d_species_c_p[si];
                c_v[idx_mixture_thermo_properties] += Y[si][idx_mass_fractions]*d_species_c_v[si];
                
                // Compute the mass fraction of the last species.
                Y_last[idx_mass_fractions] -= Y[si][idx_mass_fractions];
            }
        }
        
        // Add the contribution from the last species and compute gamma and R.
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_mixture_thermo_properties = i + num_ghosts_0_mixture_thermo_properties;
            const int idx_mass_fractions = i + num_ghosts_0_mass_fractions;
            
            c_p[idx_mixture_thermo_properties] += Y_last[idx_mass_fractions]*d_species_c_p.back();
            c_v[idx_mixture_thermo_properties] += Y_last[idx_mass_fractions]*d_species_c_v.back();
            
            gamma[idx_mixture_thermo_properties] = c_p[idx_mixture_thermo_properties]/
                c_v[idx_mixture_thermo_properties];
            R[idx_mixture_thermo_properties] = c_p[idx_mixture_thermo_properties] -
                c_v[idx_mixture_thermo_properties];
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the local lower indices, numbers of cells in each dimension and numbers of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        
        const int num_ghosts_0_mixture_thermo_properties = num_ghosts_mixture_thermo_properties[0];
        const int num_ghosts_1_mixture_thermo_properties = num_ghosts_mixture_thermo_properties[1];
        const int ghostcell_dim_0_mixture_thermo_properties = ghostcell_dims_mixture_thermo_properties[0];
        
        const int num_ghosts_0_mass_fractions = num_ghosts_mass_fractions[0];
        const int num_ghosts_1_mass_fractions = num_ghosts_mass_fractions[1];
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
                    const int idx_mixture_thermo_properties = (i + num_ghosts_0_mixture_thermo_properties) +
                        (j + num_ghosts_1_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties;
                    
                    const int idx_mass_fractions = (i + num_ghosts_0_mass_fractions) +
                        (j + num_ghosts_1_mass_fractions)*ghostcell_dim_0_mass_fractions;
                    
                    c_p[idx_mixture_thermo_properties] += Y[si][idx_mass_fractions]*d_species_c_p[si];
                    c_v[idx_mixture_thermo_properties] += Y[si][idx_mass_fractions]*d_species_c_v[si];
                    
                    // Compute the mass fraction of the last species.
                    Y_last[idx_mass_fractions] -= Y[si][idx_mass_fractions];
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
                const int idx_mixture_thermo_properties = (i + num_ghosts_0_mixture_thermo_properties) +
                    (j + num_ghosts_1_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties;
                
                const int idx_mass_fractions = (i + num_ghosts_0_mass_fractions) +
                    (j + num_ghosts_1_mass_fractions)*ghostcell_dim_0_mass_fractions;
                
                c_p[idx_mixture_thermo_properties] += Y_last[idx_mass_fractions]*d_species_c_p.back();
                c_v[idx_mixture_thermo_properties] += Y_last[idx_mass_fractions]*d_species_c_v.back();
                
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
         * Get the local lower indices, numbers of cells in each dimension and numbers of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_lo_2 = domain_lo[2];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        const int domain_dim_2 = domain_dims[2];
        
        const int num_ghosts_0_mixture_thermo_properties = num_ghosts_mixture_thermo_properties[0];
        const int num_ghosts_1_mixture_thermo_properties = num_ghosts_mixture_thermo_properties[1];
        const int num_ghosts_2_mixture_thermo_properties = num_ghosts_mixture_thermo_properties[2];
        const int ghostcell_dim_0_mixture_thermo_properties = ghostcell_dims_mixture_thermo_properties[0];
        const int ghostcell_dim_1_mixture_thermo_properties = ghostcell_dims_mixture_thermo_properties[1];
        
        const int num_ghosts_0_mass_fractions = num_ghosts_mass_fractions[0];
        const int num_ghosts_1_mass_fractions = num_ghosts_mass_fractions[1];
        const int num_ghosts_2_mass_fractions = num_ghosts_mass_fractions[2];
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
                        const int idx_mixture_thermo_properties = (i + num_ghosts_0_mixture_thermo_properties) +
                            (j + num_ghosts_1_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties +
                            (k + num_ghosts_2_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties*
                                ghostcell_dim_1_mixture_thermo_properties;
                        
                        const int idx_mass_fractions = (i + num_ghosts_0_mass_fractions) +
                            (j + num_ghosts_1_mass_fractions)*ghostcell_dim_0_mass_fractions +
                            (k + num_ghosts_2_mass_fractions)*ghostcell_dim_0_mass_fractions*
                                ghostcell_dim_1_mass_fractions;
                        
                        c_p[idx_mixture_thermo_properties] += Y[si][idx_mass_fractions]*d_species_c_p[si];
                        c_v[idx_mixture_thermo_properties] += Y[si][idx_mass_fractions]*d_species_c_v[si];
                        
                        // Compute the mass fraction of the last species.
                        Y_last[idx_mass_fractions] -= Y[si][idx_mass_fractions];
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
                    const int idx_mixture_thermo_properties = (i + num_ghosts_0_mixture_thermo_properties) +
                        (j + num_ghosts_1_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties +
                        (k + num_ghosts_2_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties*
                            ghostcell_dim_1_mixture_thermo_properties;
                    
                    const int idx_mass_fractions = (i + num_ghosts_0_mass_fractions) +
                        (j + num_ghosts_1_mass_fractions)*ghostcell_dim_0_mass_fractions +
                        (k + num_ghosts_2_mass_fractions)*ghostcell_dim_0_mass_fractions*
                            ghostcell_dim_1_mass_fractions;
                        
                    c_p[idx_mixture_thermo_properties] += Y_last[idx_mass_fractions]*d_species_c_p.back();
                    c_v[idx_mixture_thermo_properties] += Y_last[idx_mass_fractions]*d_species_c_v.back();
                    
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
EquationOfStateMixingRulesIdealGas::getMixtureThermodynamicPropertiesWithVolumeFraction(
    double* const gamma,
    const std::vector<const double*> Z,
    const hier::IntVector& num_ghosts_mixture_thermo_properties,
    const hier::IntVector& num_ghosts_volume_fractions,
    const hier::IntVector& ghostcell_dims_mixture_thermo_properties,
    const hier::IntVector& ghostcell_dims_volume_fractions,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims) const
{
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and numbers of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int num_ghosts_0_mixture_thermo_properties = num_ghosts_mixture_thermo_properties[0];
        const int num_ghosts_0_volume_fractions = num_ghosts_volume_fractions[0];
        
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
                const int idx_mixture_thermo_properties = i + num_ghosts_0_mixture_thermo_properties;
                const int idx_volume_fractions = i + num_ghosts_0_volume_fractions;
                
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
            const int idx_mixture_thermo_properties = i + num_ghosts_0_mixture_thermo_properties;
            
            gamma[idx_mixture_thermo_properties] = double(1)/gamma[idx_mixture_thermo_properties] + double(1);
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the local lower indices, numbers of cells in each dimension and numbers of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        
        const int num_ghosts_0_mixture_thermo_properties = num_ghosts_mixture_thermo_properties[0];
        const int num_ghosts_1_mixture_thermo_properties = num_ghosts_mixture_thermo_properties[1];
        const int ghostcell_dim_0_mixture_thermo_properties = ghostcell_dims_mixture_thermo_properties[0];
        
        const int num_ghosts_0_volume_fractions = num_ghosts_volume_fractions[0];
        const int num_ghosts_1_volume_fractions = num_ghosts_volume_fractions[1];
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
                    const int idx_mixture_thermo_properties = (i + num_ghosts_0_mixture_thermo_properties) +
                        (j + num_ghosts_1_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties;
                    
                    const int idx_volume_fractions = (i + num_ghosts_0_volume_fractions) +
                        (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                    
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
                const int idx_mixture_thermo_properties = (i + num_ghosts_0_mixture_thermo_properties) +
                    (j + num_ghosts_1_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties;
                
                gamma[idx_mixture_thermo_properties] = double(1)/gamma[idx_mixture_thermo_properties] + double(1);
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the local lower indices, numbers of cells in each dimension and numbers of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_lo_2 = domain_lo[2];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        const int domain_dim_2 = domain_dims[2];
        
        const int num_ghosts_0_mixture_thermo_properties = num_ghosts_mixture_thermo_properties[0];
        const int num_ghosts_1_mixture_thermo_properties = num_ghosts_mixture_thermo_properties[1];
        const int num_ghosts_2_mixture_thermo_properties = num_ghosts_mixture_thermo_properties[2];
        const int ghostcell_dim_0_mixture_thermo_properties = ghostcell_dims_mixture_thermo_properties[0];
        const int ghostcell_dim_1_mixture_thermo_properties = ghostcell_dims_mixture_thermo_properties[1];
        
        const int num_ghosts_0_volume_fractions = num_ghosts_volume_fractions[0];
        const int num_ghosts_1_volume_fractions = num_ghosts_volume_fractions[1];
        const int num_ghosts_2_volume_fractions = num_ghosts_volume_fractions[2];
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
                        const int idx_mixture_thermo_properties = (i + num_ghosts_0_mixture_thermo_properties) +
                            (j + num_ghosts_1_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties +
                            (k + num_ghosts_2_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties*
                                ghostcell_dim_1_mixture_thermo_properties;
                        
                        const int idx_volume_fractions = (i + num_ghosts_0_volume_fractions) +
                            (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                            (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
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
                    const int idx_mixture_thermo_properties = (i + num_ghosts_0_mixture_thermo_properties) +
                        (j + num_ghosts_1_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties +
                        (k + num_ghosts_2_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties*
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
EquationOfStateMixingRulesIdealGas::getMixtureThermodynamicPropertiesWithVolumeFraction(
    double* const gamma,
    double* const Z_last,
    const std::vector<const double*> Z,
    const hier::IntVector& num_ghosts_mixture_thermo_properties,
    const hier::IntVector& num_ghosts_volume_fractions,
    const hier::IntVector& ghostcell_dims_mixture_thermo_properties,
    const hier::IntVector& ghostcell_dims_volume_fractions,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims) const
{
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and numbers of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int num_ghosts_0_mixture_thermo_properties = num_ghosts_mixture_thermo_properties[0];
        const int num_ghosts_0_volume_fractions = num_ghosts_volume_fractions[0];
        
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
                const int idx_mixture_thermo_properties = i + num_ghosts_0_mixture_thermo_properties;
                const int idx_volume_fractions = i + num_ghosts_0_volume_fractions;
                
                gamma[idx_mixture_thermo_properties] += Z[si][idx_volume_fractions]*one_over_denominator;
                
                // Compute the volume fraction of the last species.
                Z_last[idx_volume_fractions] -= Z[si][idx_volume_fractions];
            }
        }
        
        // Add the contribution from the last species and compute gamma.
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_mixture_thermo_properties = i + num_ghosts_0_mixture_thermo_properties;
            const int idx_volume_fractions = i + num_ghosts_0_volume_fractions;
            
            gamma[idx_mixture_thermo_properties] += Z_last[idx_volume_fractions]/
                (d_species_gamma.back() - double(1));
            gamma[idx_mixture_thermo_properties] = double(1)/gamma[idx_mixture_thermo_properties] + double(1);
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the local lower indices, numbers of cells in each dimension and numbers of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        
        const int num_ghosts_0_mixture_thermo_properties = num_ghosts_mixture_thermo_properties[0];
        const int num_ghosts_1_mixture_thermo_properties = num_ghosts_mixture_thermo_properties[1];
        const int ghostcell_dim_0_mixture_thermo_properties = ghostcell_dims_mixture_thermo_properties[0];
        
        const int num_ghosts_0_volume_fractions = num_ghosts_volume_fractions[0];
        const int num_ghosts_1_volume_fractions = num_ghosts_volume_fractions[1];
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
                    const int idx_mixture_thermo_properties = (i + num_ghosts_0_mixture_thermo_properties) +
                        (j + num_ghosts_1_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties;
                    
                    const int idx_volume_fractions = (i + num_ghosts_0_volume_fractions) +
                        (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                    
                    gamma[idx_mixture_thermo_properties] += Z[si][idx_volume_fractions]*one_over_denominator;
                    
                    // Compute the volume fraction of the last species.
                    Z_last[idx_volume_fractions] -= Z[si][idx_volume_fractions];
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
                const int idx_mixture_thermo_properties = (i + num_ghosts_0_mixture_thermo_properties) +
                    (j + num_ghosts_1_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties;
                
                const int idx_volume_fractions = (i + num_ghosts_0_volume_fractions) +
                    (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                
                gamma[idx_mixture_thermo_properties] += Z_last[idx_volume_fractions]/
                    (d_species_gamma.back() - double(1));
                gamma[idx_mixture_thermo_properties] = double(1)/gamma[idx_mixture_thermo_properties] + double(1);
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the local lower indices, numbers of cells in each dimension and numbers of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_lo_2 = domain_lo[2];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        const int domain_dim_2 = domain_dims[2];
        
        const int num_ghosts_0_mixture_thermo_properties = num_ghosts_mixture_thermo_properties[0];
        const int num_ghosts_1_mixture_thermo_properties = num_ghosts_mixture_thermo_properties[1];
        const int num_ghosts_2_mixture_thermo_properties = num_ghosts_mixture_thermo_properties[2];
        const int ghostcell_dim_0_mixture_thermo_properties = ghostcell_dims_mixture_thermo_properties[0];
        const int ghostcell_dim_1_mixture_thermo_properties = ghostcell_dims_mixture_thermo_properties[1];
        
        const int num_ghosts_0_volume_fractions = num_ghosts_volume_fractions[0];
        const int num_ghosts_1_volume_fractions = num_ghosts_volume_fractions[1];
        const int num_ghosts_2_volume_fractions = num_ghosts_volume_fractions[2];
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
                        const int idx_mixture_thermo_properties = (i + num_ghosts_0_mixture_thermo_properties) +
                            (j + num_ghosts_1_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties +
                            (k + num_ghosts_2_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties*
                                ghostcell_dim_1_mixture_thermo_properties;
                        
                        const int idx_volume_fractions = (i + num_ghosts_0_volume_fractions) +
                            (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                            (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                ghostcell_dim_1_volume_fractions;
                        
                        gamma[idx_mixture_thermo_properties] += Z[si][idx_volume_fractions]*one_over_denominator;
                        
                        // Compute the volume fraction of the last species.
                        Z_last[idx_volume_fractions] -= Z[si][idx_volume_fractions];
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
                    const int idx_mixture_thermo_properties = (i + num_ghosts_0_mixture_thermo_properties) +
                        (j + num_ghosts_1_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties +
                        (k + num_ghosts_2_mixture_thermo_properties)*ghostcell_dim_0_mixture_thermo_properties*
                            ghostcell_dim_1_mixture_thermo_properties;
                    
                    const int idx_volume_fractions = (i + num_ghosts_0_volume_fractions) +
                        (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                        (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                            ghostcell_dim_1_volume_fractions;
                        
                    gamma[idx_mixture_thermo_properties] += Z_last[idx_volume_fractions]/
                        (d_species_gamma.back() - double(1));
                    gamma[idx_mixture_thermo_properties] = double(1)/gamma[idx_mixture_thermo_properties] + double(1);
                }
            }
        }
    }
}
