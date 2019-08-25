#include "util/mixing_rules/equations_of_state/ideal_gas/EquationOfStateIdealGas.hpp"

/*
 * Print all characteristics of the equation of state class.
 */
void
EquationOfStateIdealGas::printClassData(
    std::ostream& os) const
{
    os << "\nPrint EquationOfStateIdealGas object..."
       << std::endl;
       
    os << std::endl;
    
    os << "EquationOfStateIdealGas: this = "
       << (EquationOfStateIdealGas *)this
       << std::endl;
    
    os << "d_object_name = "
       << d_object_name
       << std::endl;
}


/*
 * Compute the pressure.
 */
double
EquationOfStateIdealGas::getPressure(
    const double* const density,
    const double* const internal_energy,
    const std::vector<const double*>& thermo_properties) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 1);
#endif
    
    const double& gamma = *(thermo_properties[0]);
    
    const double& rho = *density;
    const double& epsilon = *internal_energy;
    
    return (gamma - double(1))*rho*epsilon; // Return p.
}


/*
 * Compute the pressure.
 */
void
EquationOfStateIdealGas::computePressure(
    boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
    const std::vector<const double*>& thermo_properties,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_internal_energy);
    
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 1);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_pressure->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_density->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_internal_energy->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_pressure =
        data_pressure->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_density =
        data_density->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_internal_energy = data_internal_energy->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_internal_energy =
        data_internal_energy->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_pressure;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_internal_energy, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_internal_energy->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
    /*
     * Get the pointers to the cell data.
     */
    
    double* const p = data_pressure->getPointer(0);
    const double* const rho = data_density->getPointer(0);
    const double* const epsilon = data_internal_energy->getPointer(0);
    
    const double& gamma = *(thermo_properties[0]);
    
    computePressure(
        p,
        rho,
        epsilon,
        gamma,
        num_ghosts_pressure,
        num_ghosts_density,
        num_ghosts_internal_energy,
        ghostcell_dims_pressure,
        ghostcell_dims_density,
        ghostcell_dims_internal_energy,
        domain_lo,
        domain_dims);
}


/*
 * Compute the pressure.
 */
void
EquationOfStateIdealGas::computePressure(
    boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_internal_energy,
    const std::vector<const double*>& thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_internal_energy);
    
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 1);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_pressure->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_density->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_internal_energy->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
    hier::IntVector ghostcell_dims_pressure = data_pressure->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
    hier::IntVector ghostcell_dims_density = data_density->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_internal_energy = data_internal_energy->getGhostCellWidth();
    hier::IntVector ghostcell_dims_internal_energy = data_internal_energy->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_pressure;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_internal_energy, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_internal_energy->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(side_normal < d_dim.getValue());
    
    TBOX_ASSERT(data_pressure->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_density->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_internal_energy->getDirectionVector()[side_normal] > 0);
#endif
    
    ghostcell_dims_pressure[side_normal]++;
    ghostcell_dims_density[side_normal]++;
    ghostcell_dims_internal_energy[side_normal]++;
    domain_dims[side_normal]++;
    
    /*
     * Get the pointers to the cell data.
     */
    
    double* const p = data_pressure->getPointer(side_normal, 0);
    const double* const rho = data_density->getPointer(side_normal, 0);
    const double* const epsilon = data_internal_energy->getPointer(side_normal, 0);
    
    const double& gamma = *(thermo_properties[0]);
    
    computePressure(
        p,
        rho,
        epsilon,
        gamma,
        num_ghosts_pressure,
        num_ghosts_density,
        num_ghosts_internal_energy,
        ghostcell_dims_pressure,
        ghostcell_dims_density,
        ghostcell_dims_internal_energy,
        domain_lo,
        domain_dims);
}


/*
 * Compute the pressure.
 */
void
EquationOfStateIdealGas::computePressure(
    boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
    const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_internal_energy);
    TBOX_ASSERT(data_thermo_properties);
    
    TBOX_ASSERT(data_thermo_properties->getDepth() >= 1);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_pressure->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_density->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_internal_energy->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_thermo_properties->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_pressure =
        data_pressure->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_density =
        data_density->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_internal_energy = data_internal_energy->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_internal_energy =
        data_internal_energy->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_thermo_properties = data_thermo_properties->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_thermo_properties =
        data_thermo_properties->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_pressure;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_internal_energy, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_thermo_properties, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_internal_energy->getGhostBox().contains(domain));
        TBOX_ASSERT(data_thermo_properties->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
    /*
     * Get the pointers to the cell data.
     */
    
    double* p = data_pressure->getPointer(0);
    const double* const rho = data_density->getPointer(0);
    const double* const epsilon = data_internal_energy->getPointer(0);
    const double* const gamma = data_thermo_properties->getPointer(0);
    
    computePressure(
        p,
        rho,
        epsilon,
        gamma,
        num_ghosts_pressure,
        num_ghosts_density,
        num_ghosts_internal_energy,
        num_ghosts_thermo_properties,
        ghostcell_dims_pressure,
        ghostcell_dims_density,
        ghostcell_dims_internal_energy,
        ghostcell_dims_thermo_properties,
        domain_lo,
        domain_dims);
}


/*
 * Compute the pressure.
 */
void
EquationOfStateIdealGas::computePressure(
    boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_internal_energy,
    const boost::shared_ptr<pdat::SideData<double> >& data_thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_internal_energy);
    TBOX_ASSERT(data_thermo_properties);
    
    TBOX_ASSERT(data_thermo_properties->getDepth() >= 1);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_pressure->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_density->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_internal_energy->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_thermo_properties->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
    hier::IntVector ghostcell_dims_pressure = data_pressure->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
    hier::IntVector ghostcell_dims_density = data_density->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_internal_energy = data_internal_energy->getGhostCellWidth();
    hier::IntVector ghostcell_dims_internal_energy = data_internal_energy->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_thermo_properties = data_thermo_properties->getGhostCellWidth();
    hier::IntVector ghostcell_dims_thermo_properties = data_thermo_properties->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_pressure;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_internal_energy, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_thermo_properties, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_internal_energy->getGhostBox().contains(domain));
        TBOX_ASSERT(data_thermo_properties->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(side_normal < d_dim.getValue());
    
    TBOX_ASSERT(data_pressure->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_density->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_internal_energy->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_thermo_properties->getDirectionVector()[side_normal] > 0);
#endif
    
    ghostcell_dims_pressure[side_normal]++;
    ghostcell_dims_density[side_normal]++;
    ghostcell_dims_internal_energy[side_normal]++;
    ghostcell_dims_thermo_properties[side_normal]++;
    domain_dims[side_normal]++;
    
    /*
     * Get the pointers to the cell data.
     */
    
    double* p = data_pressure->getPointer(side_normal, 0);
    const double* const rho = data_density->getPointer(side_normal, 0);
    const double* const epsilon = data_internal_energy->getPointer(side_normal, 0);
    const double* const gamma = data_thermo_properties->getPointer(side_normal, 0);
    
    computePressure(
        p,
        rho,
        epsilon,
        gamma,
        num_ghosts_pressure,
        num_ghosts_density,
        num_ghosts_internal_energy,
        num_ghosts_thermo_properties,
        ghostcell_dims_pressure,
        ghostcell_dims_density,
        ghostcell_dims_internal_energy,
        ghostcell_dims_thermo_properties,
        domain_lo,
        domain_dims);
}


/*
 * Compute the sound speed.
 */
double
EquationOfStateIdealGas::getSoundSpeed(
    const double* const density,
    const double* const pressure,
    const std::vector<const double*>& thermo_properties) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 1);
#endif
    
    const double& gamma = *(thermo_properties[0]);
    
    const double& rho = *density;
    const double& p = *pressure;
    
    return sqrt(gamma*p/rho); // Return c.
}


/*
 * Compute the sound speed.
 */
void
EquationOfStateIdealGas::computeSoundSpeed(
    boost::shared_ptr<pdat::CellData<double> >& data_sound_speed,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const std::vector<const double*>& thermo_properties,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_sound_speed);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 1);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_sound_speed->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_density->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_pressure->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_sound_speed = data_sound_speed->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_sound_speed =
        data_sound_speed->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_density =
        data_density->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_pressure =
        data_pressure->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_sound_speed;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_sound_speed->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
    /*
     * Get the pointers to the cell data.
     */
    
    double* const c = data_sound_speed->getPointer(0);
    const double* const rho = data_density->getPointer(0);
    const double* const p = data_pressure->getPointer(0);
    
    const double& gamma = *(thermo_properties[0]);
    
    computeSoundSpeed(
        c,
        rho,
        p,
        gamma,
        num_ghosts_sound_speed,
        num_ghosts_density,
        num_ghosts_pressure,
        ghostcell_dims_sound_speed,
        ghostcell_dims_density,
        ghostcell_dims_pressure,
        domain_lo,
        domain_dims);
}


/*
 * Compute the sound speed.
 */
void
EquationOfStateIdealGas::computeSoundSpeed(
    boost::shared_ptr<pdat::SideData<double> >& data_sound_speed,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const std::vector<const double*>& thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_sound_speed);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 1);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_sound_speed->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_density->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_pressure->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_sound_speed = data_sound_speed->getGhostCellWidth();
    hier::IntVector ghostcell_dims_sound_speed = data_sound_speed->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
    hier::IntVector ghostcell_dims_density = data_density->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
    hier::IntVector ghostcell_dims_pressure = data_pressure->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_sound_speed;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_sound_speed->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(side_normal < d_dim.getValue());
    
    TBOX_ASSERT(data_sound_speed->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_density->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_pressure->getDirectionVector()[side_normal] > 0);
#endif
    
    ghostcell_dims_sound_speed[side_normal]++;
    ghostcell_dims_density[side_normal]++;
    ghostcell_dims_pressure[side_normal]++;
    domain_dims[side_normal]++;
    
    /*
     * Get the pointers to the cell data.
     */
    
    double* const c = data_sound_speed->getPointer(side_normal, 0);
    const double* const rho = data_density->getPointer(side_normal, 0);
    const double* const p = data_pressure->getPointer(side_normal, 0);
    
    const double& gamma = *(thermo_properties[0]);
    
    computeSoundSpeed(
        c,
        rho,
        p,
        gamma,
        num_ghosts_sound_speed,
        num_ghosts_density,
        num_ghosts_pressure,
        ghostcell_dims_sound_speed,
        ghostcell_dims_density,
        ghostcell_dims_pressure,
        domain_lo,
        domain_dims);
}


/*
 * Compute the sound speed.
 */
void
EquationOfStateIdealGas::computeSoundSpeed(
    boost::shared_ptr<pdat::CellData<double> >& data_sound_speed,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_sound_speed);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_thermo_properties);
    
    TBOX_ASSERT(data_thermo_properties->getDepth() >= 1);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_sound_speed->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_density->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_pressure->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_thermo_properties->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_sound_speed = data_sound_speed->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_sound_speed =
        data_sound_speed->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_density =
        data_density->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_pressure =
        data_pressure->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_thermo_properties = data_thermo_properties->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_thermo_properties =
        data_thermo_properties->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_sound_speed;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_thermo_properties, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_sound_speed->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_thermo_properties->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
    /*
     * Get the pointers to the cell data.
     */
    
    double* const c = data_sound_speed->getPointer(0);
    const double* const rho = data_density->getPointer(0);
    const double* const p = data_pressure->getPointer(0);
    const double* const gamma = data_thermo_properties->getPointer(0);
    
    computeSoundSpeed(
        c,
        rho,
        p,
        gamma,
        num_ghosts_sound_speed,
        num_ghosts_density,
        num_ghosts_pressure,
        num_ghosts_thermo_properties,
        ghostcell_dims_sound_speed,
        ghostcell_dims_density,
        ghostcell_dims_pressure,
        ghostcell_dims_thermo_properties,
        domain_lo,
        domain_dims);
}


/*
 * Compute the sound speed.
 */
void
EquationOfStateIdealGas::computeSoundSpeed(
    boost::shared_ptr<pdat::SideData<double> >& data_sound_speed,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const boost::shared_ptr<pdat::SideData<double> >& data_thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_sound_speed);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_thermo_properties);
    
    TBOX_ASSERT(data_thermo_properties->getDepth() >= 1);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_sound_speed->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_density->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_pressure->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_thermo_properties->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_sound_speed = data_sound_speed->getGhostCellWidth();
    hier::IntVector ghostcell_dims_sound_speed = data_sound_speed->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
    hier::IntVector ghostcell_dims_density = data_density->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
    hier::IntVector ghostcell_dims_pressure = data_pressure->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_thermo_properties = data_thermo_properties->getGhostCellWidth();
    hier::IntVector ghostcell_dims_thermo_properties = data_thermo_properties->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_sound_speed;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_thermo_properties, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_sound_speed->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_thermo_properties->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(side_normal < d_dim.getValue());
    
    TBOX_ASSERT(data_sound_speed->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_density->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_pressure->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_thermo_properties->getDirectionVector()[side_normal] > 0);
#endif
    
    ghostcell_dims_sound_speed[side_normal]++;
    ghostcell_dims_density[side_normal]++;
    ghostcell_dims_pressure[side_normal]++;
    ghostcell_dims_thermo_properties[side_normal]++;
    domain_dims[side_normal]++;
    
    /*
     * Get the pointers to the cell data.
     */
    
    double* const c = data_sound_speed->getPointer(side_normal, 0);
    const double* const rho = data_density->getPointer(side_normal, 0);
    const double* const p = data_pressure->getPointer(side_normal, 0);
    const double* const gamma = data_thermo_properties->getPointer(side_normal, 0);
    
    computeSoundSpeed(
        c,
        rho,
        p,
        gamma,
        num_ghosts_sound_speed,
        num_ghosts_density,
        num_ghosts_pressure,
        num_ghosts_thermo_properties,
        ghostcell_dims_sound_speed,
        ghostcell_dims_density,
        ghostcell_dims_pressure,
        ghostcell_dims_thermo_properties,
        domain_lo,
        domain_dims);
}


/*
 * Compute the specific internal energy.
 */
double
EquationOfStateIdealGas::getInternalEnergy(
    const double* const density,
    const double* const pressure,
    const std::vector<const double*>& thermo_properties) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 1);
#endif
    
    const double& gamma = *(thermo_properties[0]);
    
    const double& rho = *density;
    const double& p = *pressure;
    
    return p/((gamma - double(1))*rho); // Return epsilon.
}


/*
 * Compute the specific internal energy.
 */
void
EquationOfStateIdealGas::computeInternalEnergy(
    boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const std::vector<const double*>& thermo_properties,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_internal_energy);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 1);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_internal_energy->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_density->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_pressure->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_internal_energy = data_internal_energy->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_internal_energy =
        data_internal_energy->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_density =
        data_density->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_pressure =
        data_pressure->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_internal_energy;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_internal_energy->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
    /*
     * Get the pointers to the cell data.
     */
    
    double* const epsilon = data_internal_energy->getPointer(0);
    const double* const rho = data_density->getPointer(0);
    const double* const p = data_pressure->getPointer(0);
    
    const double& gamma = *(thermo_properties[0]);
    
    computeInternalEnergy(
        epsilon,
        rho,
        p,
        gamma,
        num_ghosts_internal_energy,
        num_ghosts_density,
        num_ghosts_pressure,
        ghostcell_dims_internal_energy,
        ghostcell_dims_density,
        ghostcell_dims_pressure,
        domain_lo,
        domain_dims);
}


/*
 * Compute the specific internal energy.
 */
void
EquationOfStateIdealGas::computeInternalEnergy(
    boost::shared_ptr<pdat::SideData<double> >& data_internal_energy,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const std::vector<const double*>& thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_internal_energy);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 1);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_internal_energy->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_density->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_pressure->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_internal_energy = data_internal_energy->getGhostCellWidth();
    hier::IntVector ghostcell_dims_internal_energy = data_internal_energy->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
    hier::IntVector ghostcell_dims_density = data_density->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
    hier::IntVector ghostcell_dims_pressure = data_pressure->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_internal_energy;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_internal_energy->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(side_normal < d_dim.getValue());
    
    TBOX_ASSERT(data_internal_energy->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_density->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_pressure->getDirectionVector()[side_normal] > 0);
#endif
    
    ghostcell_dims_internal_energy[side_normal]++;
    ghostcell_dims_density[side_normal]++;
    ghostcell_dims_pressure[side_normal]++;
    domain_dims[side_normal]++;
    
    /*
     * Get the pointers to the cell data.
     */
    
    double* const epsilon = data_internal_energy->getPointer(side_normal, 0);
    const double* const rho = data_density->getPointer(side_normal, 0);
    const double* const p = data_pressure->getPointer(side_normal, 0);
    
    const double& gamma = *(thermo_properties[0]);
    
    computeInternalEnergy(
        epsilon,
        rho,
        p,
        gamma,
        num_ghosts_internal_energy,
        num_ghosts_density,
        num_ghosts_pressure,
        ghostcell_dims_internal_energy,
        ghostcell_dims_density,
        ghostcell_dims_pressure,
        domain_lo,
        domain_dims);
}


/*
 * Compute the specific internal energy.
 */
void
EquationOfStateIdealGas::computeInternalEnergy(
    boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_internal_energy);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_thermo_properties);
    
    TBOX_ASSERT(data_thermo_properties->getDepth() >= 1);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_internal_energy->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_density->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_pressure->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_thermo_properties->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_internal_energy = data_internal_energy->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_internal_energy =
        data_internal_energy->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_density =
        data_density->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_pressure =
        data_pressure->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_thermo_properties = data_thermo_properties->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_thermo_properties =
        data_thermo_properties->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_internal_energy;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_thermo_properties, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_internal_energy->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_thermo_properties->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
    /*
     * Get the pointers to the cell data.
     */
    
    double* const epsilon = data_internal_energy->getPointer(0);
    const double* const rho = data_density->getPointer(0);
    const double* const p = data_pressure->getPointer(0);
    const double* const gamma = data_thermo_properties->getPointer(0);
    
    computeInternalEnergy(
        epsilon,
        rho,
        p,
        gamma,
        num_ghosts_internal_energy,
        num_ghosts_density,
        num_ghosts_pressure,
        num_ghosts_thermo_properties,
        ghostcell_dims_internal_energy,
        ghostcell_dims_density,
        ghostcell_dims_pressure,
        ghostcell_dims_thermo_properties,
        domain_lo,
        domain_dims);
}


/*
 * Compute the specific internal energy.
 */
void
EquationOfStateIdealGas::computeInternalEnergy(
    boost::shared_ptr<pdat::SideData<double> >& data_internal_energy,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const boost::shared_ptr<pdat::SideData<double> >& data_thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_internal_energy);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_thermo_properties);
    
    TBOX_ASSERT(data_thermo_properties->getDepth() >= 1);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_internal_energy->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_density->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_pressure->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_thermo_properties->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_internal_energy = data_internal_energy->getGhostCellWidth();
    hier::IntVector ghostcell_dims_internal_energy = data_internal_energy->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
    hier::IntVector ghostcell_dims_density = data_density->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
    hier::IntVector ghostcell_dims_pressure = data_pressure->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_thermo_properties = data_thermo_properties->getGhostCellWidth();
    hier::IntVector ghostcell_dims_thermo_properties = data_thermo_properties->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_internal_energy;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_thermo_properties, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_internal_energy->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_thermo_properties->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(side_normal < d_dim.getValue());
    
    TBOX_ASSERT(data_internal_energy->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_density->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_pressure->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_thermo_properties->getDirectionVector()[side_normal] > 0);
#endif
    
    ghostcell_dims_internal_energy[side_normal]++;
    ghostcell_dims_density[side_normal]++;
    ghostcell_dims_pressure[side_normal]++;
    ghostcell_dims_thermo_properties[side_normal]++;
    domain_dims[side_normal]++;
    
    /*
     * Get the pointers to the cell data.
     */
    
    double* const epsilon = data_internal_energy->getPointer(side_normal, 0);
    const double* const rho = data_density->getPointer(side_normal, 0);
    const double* const p = data_pressure->getPointer(side_normal, 0);
    const double* const gamma = data_thermo_properties->getPointer(side_normal, 0);
    
    computeInternalEnergy(
        epsilon,
        rho,
        p,
        gamma,
        num_ghosts_internal_energy,
        num_ghosts_density,
        num_ghosts_pressure,
        num_ghosts_thermo_properties,
        ghostcell_dims_internal_energy,
        ghostcell_dims_density,
        ghostcell_dims_pressure,
        ghostcell_dims_thermo_properties,
        domain_lo,
        domain_dims);
}


/*
 * Compute the specific enthalpy.
 */
double
EquationOfStateIdealGas::getEnthalpy(
    const double* const density,
    const double* const pressure,
    const std::vector<const double*>& thermo_properties) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 1);
#endif
    
    const double& gamma = *(thermo_properties[0]);
    
    const double& rho = *density;
    const double& p = *pressure;
    
    return gamma*p/((gamma - double(1))*rho); // Return h.
}


/*
 * Compute the specific enthalpy.
 */
void
EquationOfStateIdealGas::computeEnthalpy(
    boost::shared_ptr<pdat::CellData<double> >& data_enthalpy,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const std::vector<const double*>& thermo_properties,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_enthalpy);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 1);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_enthalpy->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_density->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_pressure->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_enthalpy = data_enthalpy->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_enthalpy =
        data_enthalpy->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_density =
        data_density->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_pressure =
        data_pressure->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_enthalpy;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_enthalpy->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
    /*
     * Get the pointers to the cell data.
     */
    
    double* const h = data_enthalpy->getPointer(0);
    const double* const rho = data_density->getPointer(0);
    const double* const p = data_pressure->getPointer(0);
    
    const double& gamma = *(thermo_properties[0]);
    
    computeEnthalpy(
        h,
        rho,
        p,
        gamma,
        num_ghosts_enthalpy,
        num_ghosts_density,
        num_ghosts_pressure,
        ghostcell_dims_enthalpy,
        ghostcell_dims_density,
        ghostcell_dims_pressure,
        domain_lo,
        domain_dims);
}


/*
 * Compute the specific enthalpy.
 */
void
EquationOfStateIdealGas::computeEnthalpy(
    boost::shared_ptr<pdat::SideData<double> >& data_enthalpy,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const std::vector<const double*>& thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_enthalpy);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 1);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_enthalpy->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_density->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_pressure->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_enthalpy = data_enthalpy->getGhostCellWidth();
    hier::IntVector ghostcell_dims_enthalpy = data_enthalpy->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
    hier::IntVector ghostcell_dims_density = data_density->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
    hier::IntVector ghostcell_dims_pressure = data_pressure->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_enthalpy;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_enthalpy->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(side_normal < d_dim.getValue());
    
    TBOX_ASSERT(data_enthalpy->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_density->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_pressure->getDirectionVector()[side_normal] > 0);
#endif
    
    ghostcell_dims_enthalpy[side_normal]++;
    ghostcell_dims_density[side_normal]++;
    ghostcell_dims_pressure[side_normal]++;
    domain_dims[side_normal]++;
    
    /*
     * Get the pointers to the cell data.
     */
    
    double* const h = data_enthalpy->getPointer(side_normal, 0);
    const double* const rho = data_density->getPointer(side_normal, 0);
    const double* const p = data_pressure->getPointer(side_normal, 0);
    
    const double& gamma = *(thermo_properties[0]);
    
    computeEnthalpy(
        h,
        rho,
        p,
        gamma,
        num_ghosts_enthalpy,
        num_ghosts_density,
        num_ghosts_pressure,
        ghostcell_dims_enthalpy,
        ghostcell_dims_density,
        ghostcell_dims_pressure,
        domain_lo,
        domain_dims);
}


/*
 * Compute the specific enthalpy.
 */
void
EquationOfStateIdealGas::computeEnthalpy(
    boost::shared_ptr<pdat::CellData<double> >& data_enthalpy,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_enthalpy);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_thermo_properties);
    
    TBOX_ASSERT(data_thermo_properties->getDepth() >= 1);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_enthalpy->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_density->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_pressure->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_thermo_properties->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_enthalpy = data_enthalpy->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_enthalpy =
        data_enthalpy->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_density =
        data_density->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_pressure =
        data_pressure->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_thermo_properties = data_thermo_properties->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_thermo_properties =
        data_thermo_properties->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_enthalpy;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_thermo_properties, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_enthalpy->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_thermo_properties->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
    /*
     * Get the pointers to the cell data.
     */
    
    double* const h = data_enthalpy->getPointer(0);
    const double* const rho = data_density->getPointer(0);
    const double* const p = data_pressure->getPointer(0);
    const double* const gamma = data_thermo_properties->getPointer(0);
    
    computeEnthalpy(
        h,
        rho,
        p,
        gamma,
        num_ghosts_enthalpy,
        num_ghosts_density,
        num_ghosts_pressure,
        num_ghosts_thermo_properties,
        ghostcell_dims_enthalpy,
        ghostcell_dims_density,
        ghostcell_dims_pressure,
        ghostcell_dims_thermo_properties,
        domain_lo,
        domain_dims);
}


/*
 * Compute the specific enthalpy.
 */
void
EquationOfStateIdealGas::computeEnthalpy(
    boost::shared_ptr<pdat::SideData<double> >& data_enthalpy,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const boost::shared_ptr<pdat::SideData<double> >& data_thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_enthalpy);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_thermo_properties);
    
    TBOX_ASSERT(data_thermo_properties->getDepth() >= 1);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_enthalpy->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_density->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_pressure->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_thermo_properties->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_enthalpy = data_enthalpy->getGhostCellWidth();
    hier::IntVector ghostcell_dims_enthalpy = data_enthalpy->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
    hier::IntVector ghostcell_dims_density = data_density->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
    hier::IntVector ghostcell_dims_pressure = data_pressure->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_thermo_properties = data_thermo_properties->getGhostCellWidth();
    hier::IntVector ghostcell_dims_thermo_properties = data_thermo_properties->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_enthalpy;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_thermo_properties, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_enthalpy->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_thermo_properties->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(side_normal < d_dim.getValue());
    
    TBOX_ASSERT(data_enthalpy->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_density->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_pressure->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_thermo_properties->getDirectionVector()[side_normal] > 0);
#endif
    
    ghostcell_dims_enthalpy[side_normal]++;
    ghostcell_dims_density[side_normal]++;
    ghostcell_dims_pressure[side_normal]++;
    ghostcell_dims_thermo_properties[side_normal]++;
    domain_dims[side_normal]++;
    
    /*
     * Get the pointers to the cell data.
     */
    
    double* const h = data_enthalpy->getPointer(side_normal, 0);
    const double* const rho = data_density->getPointer(side_normal, 0);
    const double* const p = data_pressure->getPointer(side_normal, 0);
    const double* const gamma = data_thermo_properties->getPointer(side_normal, 0);
    
    computeEnthalpy(
        h,
        rho,
        p,
        gamma,
        num_ghosts_enthalpy,
        num_ghosts_density,
        num_ghosts_pressure,
        num_ghosts_thermo_properties,
        ghostcell_dims_enthalpy,
        ghostcell_dims_density,
        ghostcell_dims_pressure,
        ghostcell_dims_thermo_properties,
        domain_lo,
        domain_dims);
}


/*
 * Compute the temperature.
 */
double
EquationOfStateIdealGas::getTemperature(
    const double* const density,
    const double* const pressure,
    const std::vector<const double*>& thermo_properties) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) == 4);
#endif
    
    const double& gamma = *(thermo_properties[0]);
    const double& c_v = *(thermo_properties[3]);
    
    const double& rho = *density;
    const double& p = *pressure;
    
    return p/((gamma - double(1))*c_v*rho); // Return T.
}


/*
 * Compute the temperature.
 */
void
EquationOfStateIdealGas::computeTemperature(
    boost::shared_ptr<pdat::CellData<double> >& data_temperature,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const std::vector<const double*>& thermo_properties,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_temperature);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) == 4);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_temperature->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_density->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_pressure->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_temperature =
        data_temperature->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_density =
        data_density->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_pressure =
        data_pressure->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_temperature;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_temperature->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
    /*
     * Get the pointers to the cell data.
     */
    
    double* const T = data_temperature->getPointer(0);
    const double* const rho = data_density->getPointer(0);
    const double* const p = data_pressure->getPointer(0);
    
    const double& gamma = *(thermo_properties[0]);
    const double& c_v = *(thermo_properties[3]);
    
    computeTemperature(
        T,
        rho,
        p,
        gamma,
        c_v,
        num_ghosts_temperature,
        num_ghosts_density,
        num_ghosts_pressure,
        ghostcell_dims_temperature,
        ghostcell_dims_density,
        ghostcell_dims_pressure,
        domain_lo,
        domain_dims);
}


/*
 * Compute the temperature.
 */
void
EquationOfStateIdealGas::computeTemperature(
    boost::shared_ptr<pdat::SideData<double> >& data_temperature,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const std::vector<const double*>& thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_temperature);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) == 4);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_temperature->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_density->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_pressure->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
    hier::IntVector ghostcell_dims_temperature = data_temperature->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
    hier::IntVector ghostcell_dims_density = data_density->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
    hier::IntVector ghostcell_dims_pressure = data_pressure->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_temperature;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_temperature->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(side_normal < d_dim.getValue());
    
    TBOX_ASSERT(data_temperature->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_density->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_pressure->getDirectionVector()[side_normal] > 0);
#endif
    
    ghostcell_dims_temperature[side_normal]++;
    ghostcell_dims_density[side_normal]++;
    ghostcell_dims_pressure[side_normal]++;
    domain_dims[side_normal]++;
    
    /*
     * Get the pointers to the cell data.
     */
    
    double* const T = data_temperature->getPointer(side_normal, 0);
    const double* const rho = data_density->getPointer(side_normal, 0);
    const double* const p = data_pressure->getPointer(side_normal, 0);
    
    const double& gamma = *(thermo_properties[0]);
    const double& c_v = *(thermo_properties[3]);
    
    computeTemperature(
        T,
        rho,
        p,
        gamma,
        c_v,
        num_ghosts_temperature,
        num_ghosts_density,
        num_ghosts_pressure,
        ghostcell_dims_temperature,
        ghostcell_dims_density,
        ghostcell_dims_pressure,
        domain_lo,
        domain_dims);
}


/*
 * Compute the temperature.
 */
void
EquationOfStateIdealGas::computeTemperature(
    boost::shared_ptr<pdat::CellData<double> >& data_temperature,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_temperature);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_thermo_properties);
    
    TBOX_ASSERT(data_thermo_properties->getDepth() == 4);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_temperature->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_density->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_pressure->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_thermo_properties->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_temperature =
        data_temperature->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_density =
        data_density->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_pressure =
        data_pressure->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_thermo_properties = data_thermo_properties->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_thermo_properties =
        data_thermo_properties->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_temperature;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_thermo_properties, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_temperature->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_thermo_properties->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
    /*
     * Get the pointers to the cell data.
     */
    
    double* const T = data_temperature->getPointer(0);
    const double* const rho = data_density->getPointer(0);
    const double* const p = data_pressure->getPointer(0);
    const double* const gamma = data_thermo_properties->getPointer(0);
    const double* const c_v = data_thermo_properties->getPointer(3);
    
    computeTemperature(
        T,
        rho,
        p,
        gamma,
        c_v,
        num_ghosts_temperature,
        num_ghosts_density,
        num_ghosts_pressure,
        num_ghosts_thermo_properties,
        ghostcell_dims_temperature,
        ghostcell_dims_density,
        ghostcell_dims_pressure,
        ghostcell_dims_thermo_properties,
        domain_lo,
        domain_dims);
}


/*
 * Compute the temperature.
 */
void
EquationOfStateIdealGas::computeTemperature(
    boost::shared_ptr<pdat::SideData<double> >& data_temperature,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const boost::shared_ptr<pdat::SideData<double> >& data_thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_temperature);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_thermo_properties);
    
    TBOX_ASSERT(data_thermo_properties->getDepth() == 4);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_temperature->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_density->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_pressure->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_thermo_properties->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
    hier::IntVector ghostcell_dims_temperature = data_temperature->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
    hier::IntVector ghostcell_dims_density = data_density->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
    hier::IntVector ghostcell_dims_pressure = data_pressure->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_thermo_properties = data_thermo_properties->getGhostCellWidth();
    hier::IntVector ghostcell_dims_thermo_properties = data_thermo_properties->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_temperature;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_thermo_properties, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_temperature->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_thermo_properties->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(side_normal < d_dim.getValue());
    
    TBOX_ASSERT(data_temperature->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_density->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_pressure->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_thermo_properties->getDirectionVector()[side_normal] > 0);
#endif
    
    ghostcell_dims_temperature[side_normal]++;
    ghostcell_dims_density[side_normal]++;
    ghostcell_dims_pressure[side_normal]++;
    ghostcell_dims_thermo_properties[side_normal]++;
    domain_dims[side_normal]++;
    
    /*
     * Get the pointers to the cell data.
     */
    
    double* const T = data_temperature->getPointer(side_normal, 0);
    const double* const rho = data_density->getPointer(side_normal, 0);
    const double* const p = data_pressure->getPointer(side_normal, 0);
    const double* const gamma = data_thermo_properties->getPointer(side_normal, 0);
    const double* const c_v = data_thermo_properties->getPointer(side_normal, 3);
    
    computeTemperature(
        T,
        rho,
        p,
        gamma,
        c_v,
        num_ghosts_temperature,
        num_ghosts_density,
        num_ghosts_pressure,
        num_ghosts_thermo_properties,
        ghostcell_dims_temperature,
        ghostcell_dims_density,
        ghostcell_dims_pressure,
        ghostcell_dims_thermo_properties,
        domain_lo,
        domain_dims);
}


/*
 * Compute the specific internal energy from temperature.
 */
double
EquationOfStateIdealGas::getInternalEnergyFromTemperature(
    const double* const density,
    const double* const temperature,
    const std::vector<const double*>& thermo_properties) const
{
    NULL_USE(density);
    
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) == 4);
#endif
    
    const double& c_v = *(thermo_properties[3]);
    
    const double& T = *temperature;
    
    return c_v*T; // Return epsilon.
}


/*
 * Compute the specific internal energy from temperature.
 */
void
EquationOfStateIdealGas::computeInternalEnergyFromTemperature(
    boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
    const std::vector<const double*>& thermo_properties,
    const hier::Box& domain) const
{
    NULL_USE(data_density);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_internal_energy);
    TBOX_ASSERT(data_temperature);
    
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) == 4);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_internal_energy->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_temperature->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_internal_energy = data_internal_energy->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_internal_energy =
        data_internal_energy->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_temperature =
        data_temperature->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_internal_energy;
        num_ghosts_min = hier::IntVector::min(num_ghosts_temperature, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_internal_energy->getGhostBox().contains(domain));
        TBOX_ASSERT(data_temperature->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
    /*
     * Get the pointers to the cell data.
     */
    
    double* const epsilon = data_internal_energy->getPointer(0);
    const double* const T = data_temperature->getPointer(0);
    
    const double& c_v = *(thermo_properties[3]);
    
    computeInternalEnergyFromTemperature(
        epsilon,
        T,
        c_v,
        num_ghosts_internal_energy,
        num_ghosts_temperature,
        ghostcell_dims_internal_energy,
        ghostcell_dims_temperature,
        domain_lo,
        domain_dims);
}


/*
 * Compute the specific internal energy from temperature.
 */
void
EquationOfStateIdealGas::computeInternalEnergyFromTemperature(
    boost::shared_ptr<pdat::SideData<double> >& data_internal_energy,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_temperature,
    const std::vector<const double*>& thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
    NULL_USE(data_density);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_internal_energy);
    TBOX_ASSERT(data_temperature);
    
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) == 4);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_internal_energy->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_temperature->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_internal_energy = data_internal_energy->getGhostCellWidth();
    hier::IntVector ghostcell_dims_internal_energy = data_internal_energy->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
    hier::IntVector ghostcell_dims_temperature = data_temperature->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_internal_energy;
        num_ghosts_min = hier::IntVector::min(num_ghosts_temperature, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_internal_energy->getGhostBox().contains(domain));
        TBOX_ASSERT(data_temperature->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(side_normal < d_dim.getValue());
    
    TBOX_ASSERT(data_internal_energy->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_temperature->getDirectionVector()[side_normal] > 0);
#endif
    
    ghostcell_dims_internal_energy[side_normal]++;
    ghostcell_dims_temperature[side_normal]++;
    domain_dims[side_normal]++;
    
    /*
     * Get the pointers to the cell data.
     */
    
    double* const epsilon = data_internal_energy->getPointer(side_normal, 0);
    const double* const T = data_temperature->getPointer(side_normal, 0);
    
    const double& c_v = *(thermo_properties[3]);
    
    computeInternalEnergyFromTemperature(
        epsilon,
        T,
        c_v,
        num_ghosts_internal_energy,
        num_ghosts_temperature,
        ghostcell_dims_internal_energy,
        ghostcell_dims_temperature,
        domain_lo,
        domain_dims);
}


/*
 * Compute the specific internal energy from temperature.
 */
void
EquationOfStateIdealGas::computeInternalEnergyFromTemperature(
    boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
    const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties,
    const hier::Box& domain) const
{
    NULL_USE(data_density);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_internal_energy);
    TBOX_ASSERT(data_temperature);
    TBOX_ASSERT(data_thermo_properties);
    
    TBOX_ASSERT(data_thermo_properties->getDepth() == 4);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_internal_energy->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_temperature->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_thermo_properties->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_internal_energy = data_internal_energy->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_internal_energy =
        data_internal_energy->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_temperature =
        data_temperature->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_thermo_properties = data_thermo_properties->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_thermo_properties =
        data_thermo_properties->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_internal_energy;
        num_ghosts_min = hier::IntVector::min(num_ghosts_temperature, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_thermo_properties, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_internal_energy->getGhostBox().contains(domain));
        TBOX_ASSERT(data_temperature->getGhostBox().contains(domain));
        TBOX_ASSERT(data_thermo_properties->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
    /*
     * Get the pointers to the cell data.
     */
    
    double* const epsilon = data_internal_energy->getPointer(0);
    const double* const T = data_temperature->getPointer(0);
    const double* const c_v = data_thermo_properties->getPointer(3);
    
    computeInternalEnergyFromTemperature(
        epsilon,
        T,
        c_v,
        num_ghosts_internal_energy,
        num_ghosts_temperature,
        num_ghosts_thermo_properties,
        ghostcell_dims_internal_energy,
        ghostcell_dims_temperature,
        ghostcell_dims_thermo_properties,
        domain_lo,
        domain_dims);
}


/*
 * Compute the specific internal energy from temperature.
 */
void
EquationOfStateIdealGas::computeInternalEnergyFromTemperature(
    boost::shared_ptr<pdat::SideData<double> >& data_internal_energy,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_temperature,
    const boost::shared_ptr<pdat::SideData<double> >& data_thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
    NULL_USE(data_density);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_internal_energy);
    TBOX_ASSERT(data_temperature);
    TBOX_ASSERT(data_thermo_properties);
    
    TBOX_ASSERT(data_thermo_properties->getDepth() == 4);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_internal_energy->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_temperature->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_thermo_properties->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_internal_energy = data_internal_energy->getGhostCellWidth();
    hier::IntVector ghostcell_dims_internal_energy = data_internal_energy->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
    hier::IntVector ghostcell_dims_temperature = data_temperature->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_thermo_properties = data_thermo_properties->getGhostCellWidth();
    hier::IntVector ghostcell_dims_thermo_properties = data_thermo_properties->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_internal_energy;
        num_ghosts_min = hier::IntVector::min(num_ghosts_temperature, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_thermo_properties, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_internal_energy->getGhostBox().contains(domain));
        TBOX_ASSERT(data_temperature->getGhostBox().contains(domain));
        TBOX_ASSERT(data_thermo_properties->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(side_normal < d_dim.getValue());
    
    TBOX_ASSERT(data_internal_energy->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_temperature->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_thermo_properties->getDirectionVector()[side_normal] > 0);
#endif
    
    ghostcell_dims_internal_energy[side_normal]++;
    ghostcell_dims_temperature[side_normal]++;
    ghostcell_dims_thermo_properties[side_normal]++;
    domain_dims[side_normal]++;
    
    /*
     * Get the pointers to the cell data.
     */
    
    double* const epsilon = data_internal_energy->getPointer(side_normal, 0);
    const double* const T = data_temperature->getPointer(side_normal, 0);
    const double* const c_v = data_thermo_properties->getPointer(side_normal, 3);
    
    computeInternalEnergyFromTemperature(
        epsilon,
        T,
        c_v,
        num_ghosts_internal_energy,
        num_ghosts_temperature,
        num_ghosts_thermo_properties,
        ghostcell_dims_internal_energy,
        ghostcell_dims_temperature,
        ghostcell_dims_thermo_properties,
        domain_lo,
        domain_dims);
}


/*
 * Compute the isochoric specific heat capacity.
 */
double
EquationOfStateIdealGas::getIsochoricSpecificHeatCapacity(
    const double* const density,
    const double* const pressure,
    const std::vector<const double*>& thermo_properties) const
{
    NULL_USE(density);
    NULL_USE(pressure);
    
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) == 4);
#endif
    
    return *(thermo_properties[3]);
}


/*
 * Compute the isochoric specific heat capacity.
 */
void
EquationOfStateIdealGas::computeIsochoricSpecificHeatCapacity(
    boost::shared_ptr<pdat::CellData<double> >& data_isochoric_specific_heat_capacity,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const std::vector<const double*>& thermo_properties,
    const hier::Box& domain) const
{
    NULL_USE(data_density);
    NULL_USE(data_pressure);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_isochoric_specific_heat_capacity);
    
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) == 4);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_isochoric_specific_heat_capacity->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell box.
     */
    
    const hier::IntVector num_ghosts_isochoric_specific_heat_capacity =
        data_isochoric_specific_heat_capacity->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_isochoric_specific_heat_capacity =
        data_isochoric_specific_heat_capacity->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_isochoric_specific_heat_capacity);
        
        domain_lo = -num_ghosts_isochoric_specific_heat_capacity;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_isochoric_specific_heat_capacity->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
    /*
     * Get the pointer to the cell data.
     */
    
    double* const c_v = data_isochoric_specific_heat_capacity->getPointer(0);
    
    const double& c_v_src = *(thermo_properties[3]);
    
    computeIsochoricSpecificHeatCapacity(
        c_v,
        c_v_src,
        num_ghosts_isochoric_specific_heat_capacity,
        ghostcell_dims_isochoric_specific_heat_capacity,
        domain_lo,
        domain_dims);
}


/*
 * Compute the isochoric specific heat capacity.
 */
void
EquationOfStateIdealGas::computeIsochoricSpecificHeatCapacity(
    boost::shared_ptr<pdat::SideData<double> >& data_isochoric_specific_heat_capacity,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const std::vector<const double*>& thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
    NULL_USE(data_density);
    NULL_USE(data_pressure);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_isochoric_specific_heat_capacity);
    
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) == 4);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_isochoric_specific_heat_capacity->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell box.
     */
    
    const hier::IntVector num_ghosts_isochoric_specific_heat_capacity =
        data_isochoric_specific_heat_capacity->getGhostCellWidth();
    hier::IntVector ghostcell_dims_isochoric_specific_heat_capacity =
        data_isochoric_specific_heat_capacity->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_isochoric_specific_heat_capacity);
        
        domain_lo = -num_ghosts_isochoric_specific_heat_capacity;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_isochoric_specific_heat_capacity->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(side_normal < d_dim.getValue());
    
    TBOX_ASSERT(data_isochoric_specific_heat_capacity->getDirectionVector()[side_normal] > 0);
#endif
    
    ghostcell_dims_isochoric_specific_heat_capacity[side_normal]++;
    domain_dims[side_normal]++;
    
    /*
     * Get the pointer to the cell data.
     */
    
    double* const c_v = data_isochoric_specific_heat_capacity->getPointer(side_normal, 0);
    
    const double& c_v_src = *(thermo_properties[3]);
    
    computeIsochoricSpecificHeatCapacity(
        c_v,
        c_v_src,
        num_ghosts_isochoric_specific_heat_capacity,
        ghostcell_dims_isochoric_specific_heat_capacity,
        domain_lo,
        domain_dims);
}


/*
 * Compute the isochoric specific heat capacity.
 */
void
EquationOfStateIdealGas::computeIsochoricSpecificHeatCapacity(
    boost::shared_ptr<pdat::CellData<double> >& data_isochoric_specific_heat_capacity,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties,
    const hier::Box& domain) const
{
    NULL_USE(data_density);
    NULL_USE(data_pressure);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_isochoric_specific_heat_capacity);
    TBOX_ASSERT(data_thermo_properties);
    
    TBOX_ASSERT(data_thermo_properties->getDepth() == 4);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_isochoric_specific_heat_capacity->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_thermo_properties->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_isochoric_specific_heat_capacity =
        data_isochoric_specific_heat_capacity->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_isochoric_specific_heat_capacity =
        data_isochoric_specific_heat_capacity->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_thermo_properties = data_thermo_properties->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_thermo_properties =
        data_thermo_properties->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_isochoric_specific_heat_capacity;
        num_ghosts_min = hier::IntVector::min(num_ghosts_thermo_properties, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_isochoric_specific_heat_capacity->getGhostBox().contains(domain));
        TBOX_ASSERT(data_thermo_properties->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
    /*
     * Get the pointers to the cell data.
     */
    
    double* const c_v = data_isochoric_specific_heat_capacity->getPointer(0);
    const double* const c_v_src = data_thermo_properties->getPointer(3);
    
    computeIsochoricSpecificHeatCapacity(
        c_v,
        c_v_src,
        num_ghosts_isochoric_specific_heat_capacity,
        num_ghosts_thermo_properties,
        ghostcell_dims_isochoric_specific_heat_capacity,
        ghostcell_dims_thermo_properties,
        domain_lo,
        domain_dims);
}


/*
 * Compute the isochoric specific heat capacity.
 */
void
EquationOfStateIdealGas::computeIsochoricSpecificHeatCapacity(
    boost::shared_ptr<pdat::SideData<double> >& data_isochoric_specific_heat_capacity,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const boost::shared_ptr<pdat::SideData<double> >& data_thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
    NULL_USE(data_density);
    NULL_USE(data_pressure);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_isochoric_specific_heat_capacity);
    TBOX_ASSERT(data_thermo_properties);
    
    TBOX_ASSERT(data_thermo_properties->getDepth() == 4);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_isochoric_specific_heat_capacity->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_thermo_properties->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_isochoric_specific_heat_capacity =
        data_isochoric_specific_heat_capacity->getGhostCellWidth();
    hier::IntVector ghostcell_dims_isochoric_specific_heat_capacity =
        data_isochoric_specific_heat_capacity->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_thermo_properties = data_thermo_properties->getGhostCellWidth();
    hier::IntVector ghostcell_dims_thermo_properties = data_thermo_properties->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_isochoric_specific_heat_capacity;
        num_ghosts_min = hier::IntVector::min(num_ghosts_thermo_properties, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_isochoric_specific_heat_capacity->getGhostBox().contains(domain));
        TBOX_ASSERT(data_thermo_properties->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(side_normal < d_dim.getValue());
    
    TBOX_ASSERT(data_isochoric_specific_heat_capacity->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_thermo_properties->getDirectionVector()[side_normal] > 0);
#endif
    
    ghostcell_dims_isochoric_specific_heat_capacity[side_normal]++;
    ghostcell_dims_thermo_properties[side_normal]++;
    domain_dims[side_normal]++;
    
    /*
     * Get the pointers to the cell data.
     */
    
    double* const c_v = data_isochoric_specific_heat_capacity->getPointer(side_normal, 0);
    const double* const c_v_src = data_thermo_properties->getPointer(side_normal, 3);
    
    computeIsochoricSpecificHeatCapacity(
        c_v,
        c_v_src,
        num_ghosts_isochoric_specific_heat_capacity,
        num_ghosts_thermo_properties,
        ghostcell_dims_isochoric_specific_heat_capacity,
        ghostcell_dims_thermo_properties,
        domain_lo,
        domain_dims);
}


/*
 * Compute the isobaric specific heat capacity.
 */
double
EquationOfStateIdealGas::getIsobaricSpecificHeatCapacity(
    const double* const density,
    const double* const pressure,
    const std::vector<const double*>& thermo_properties) const
{
    NULL_USE(density);
    NULL_USE(pressure);
    
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 3);
#endif
    
    return *(thermo_properties[2]);
}


/*
 * Compute the isobaric specific heat capacity.
 */
void
EquationOfStateIdealGas::computeIsobaricSpecificHeatCapacity(
    boost::shared_ptr<pdat::CellData<double> >& data_isobaric_specific_heat_capacity,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const std::vector<const double*>& thermo_properties,
    const hier::Box& domain) const
{
    NULL_USE(data_density);
    NULL_USE(data_pressure);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_isobaric_specific_heat_capacity);
    
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 3);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_isobaric_specific_heat_capacity->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell box.
     */
    
    const hier::IntVector num_ghosts_isobaric_specific_heat_capacity =
        data_isobaric_specific_heat_capacity->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_isobaric_specific_heat_capacity =
        data_isobaric_specific_heat_capacity->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_isobaric_specific_heat_capacity);
        
        domain_lo = -num_ghosts_isobaric_specific_heat_capacity;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_isobaric_specific_heat_capacity->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
    /*
     * Get the pointer to the cell data.
     */
    
    double* const c_p = data_isobaric_specific_heat_capacity->getPointer(0);
    
    const double& c_p_src = *(thermo_properties[2]);
    
    computeIsobaricSpecificHeatCapacity(
        c_p,
        c_p_src,
        num_ghosts_isobaric_specific_heat_capacity,
        ghostcell_dims_isobaric_specific_heat_capacity,
        domain_lo,
        domain_dims);
}


/*
 * Compute the isobaric specific heat capacity.
 */
void
EquationOfStateIdealGas::computeIsobaricSpecificHeatCapacity(
    boost::shared_ptr<pdat::SideData<double> >& data_isobaric_specific_heat_capacity,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const std::vector<const double*>& thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
    NULL_USE(data_density);
    NULL_USE(data_pressure);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_isobaric_specific_heat_capacity);
    
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 3);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_isobaric_specific_heat_capacity->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell box.
     */
    
    const hier::IntVector num_ghosts_isobaric_specific_heat_capacity =
        data_isobaric_specific_heat_capacity->getGhostCellWidth();
    hier::IntVector ghostcell_dims_isobaric_specific_heat_capacity =
        data_isobaric_specific_heat_capacity->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_isobaric_specific_heat_capacity);
        
        domain_lo = -num_ghosts_isobaric_specific_heat_capacity;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_isobaric_specific_heat_capacity->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(side_normal < d_dim.getValue());
    
    TBOX_ASSERT(data_isobaric_specific_heat_capacity->getDirectionVector()[side_normal] > 0);
#endif
    
    ghostcell_dims_isobaric_specific_heat_capacity[side_normal]++;
    domain_dims[side_normal]++;
    
    /*
     * Get the pointer to the cell data.
     */
    
    double* const c_p = data_isobaric_specific_heat_capacity->getPointer(side_normal, 0);
    
    const double& c_p_src = *(thermo_properties[2]);
    
    computeIsobaricSpecificHeatCapacity(
        c_p,
        c_p_src,
        num_ghosts_isobaric_specific_heat_capacity,
        ghostcell_dims_isobaric_specific_heat_capacity,
        domain_lo,
        domain_dims);
}


/*
 * Compute the isobaric specific heat capacity.
 */
void
EquationOfStateIdealGas::computeIsobaricSpecificHeatCapacity(
    boost::shared_ptr<pdat::CellData<double> >& data_isobaric_specific_heat_capacity,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties,
    const hier::Box& domain) const
{
    NULL_USE(data_density);
    NULL_USE(data_pressure);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_isobaric_specific_heat_capacity);
    TBOX_ASSERT(data_thermo_properties);
    
    TBOX_ASSERT(data_thermo_properties->getDepth() >= 3);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_isobaric_specific_heat_capacity->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_thermo_properties->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_isobaric_specific_heat_capacity =
        data_isobaric_specific_heat_capacity->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_isobaric_specific_heat_capacity =
        data_isobaric_specific_heat_capacity->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_thermo_properties = data_thermo_properties->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_thermo_properties =
        data_thermo_properties->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_isobaric_specific_heat_capacity;
        num_ghosts_min = hier::IntVector::min(num_ghosts_thermo_properties, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_isobaric_specific_heat_capacity->getGhostBox().contains(domain));
        TBOX_ASSERT(data_thermo_properties->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
    /*
     * Get the pointers to the cell data.
     */
    
    double* const c_p = data_isobaric_specific_heat_capacity->getPointer(0);
    const double* const c_p_src = data_thermo_properties->getPointer(2);
    
    computeIsobaricSpecificHeatCapacity(
        c_p,
        c_p_src,
        num_ghosts_isobaric_specific_heat_capacity,
        num_ghosts_thermo_properties,
        ghostcell_dims_isobaric_specific_heat_capacity,
        ghostcell_dims_thermo_properties,
        domain_lo,
        domain_dims);
}


/*
 * Compute the isobaric specific heat capacity.
 */
void
EquationOfStateIdealGas::computeIsobaricSpecificHeatCapacity(
    boost::shared_ptr<pdat::SideData<double> >& data_isobaric_specific_heat_capacity,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const boost::shared_ptr<pdat::SideData<double> >& data_thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
    NULL_USE(data_density);
    NULL_USE(data_pressure);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_isobaric_specific_heat_capacity);
    TBOX_ASSERT(data_thermo_properties);
    
    TBOX_ASSERT(data_thermo_properties->getDepth() >= 3);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_isobaric_specific_heat_capacity->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_thermo_properties->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_isobaric_specific_heat_capacity =
        data_isobaric_specific_heat_capacity->getGhostCellWidth();
    hier::IntVector ghostcell_dims_isobaric_specific_heat_capacity =
        data_isobaric_specific_heat_capacity->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_thermo_properties = data_thermo_properties->getGhostCellWidth();
    hier::IntVector ghostcell_dims_thermo_properties = data_thermo_properties->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_isobaric_specific_heat_capacity;
        num_ghosts_min = hier::IntVector::min(num_ghosts_thermo_properties, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_isobaric_specific_heat_capacity->getGhostBox().contains(domain));
        TBOX_ASSERT(data_thermo_properties->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(side_normal < d_dim.getValue());
    
    TBOX_ASSERT(data_isobaric_specific_heat_capacity->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_thermo_properties->getDirectionVector()[side_normal] > 0);
#endif
    
    ghostcell_dims_isobaric_specific_heat_capacity[side_normal]++;
    ghostcell_dims_thermo_properties[side_normal]++;
    domain_dims[side_normal]++;
    
    /*
     * Get the pointers to the cell data.
     */
    
    double* const c_p = data_isobaric_specific_heat_capacity->getPointer(side_normal, 0);
    const double* const c_p_src = data_thermo_properties->getPointer(side_normal, 2);
    
    computeIsobaricSpecificHeatCapacity(
        c_p,
        c_p_src,
        num_ghosts_isobaric_specific_heat_capacity,
        num_ghosts_thermo_properties,
        ghostcell_dims_isobaric_specific_heat_capacity,
        ghostcell_dims_thermo_properties,
        domain_lo,
        domain_dims);
}


////////////////////////////////////////////////////////////////////////////////////////////////////


/*
 * Compute the partial derivative of pressure w.r.t. density under constant specific internal energy.
 */
double
EquationOfStateIdealGas::getPartialPressurePartialDensity(
    const double* const density,
    const double* const pressure,
    const std::vector<const double*>& thermo_properties) const
{
    NULL_USE(thermo_properties);
    
    const double& rho = *density;
    const double& p = *pressure;
    
    return p/rho;
}


/*
 * Compute the partial derivative of pressure w.r.t. density under constant specific internal energy.
 */
void
EquationOfStateIdealGas::computePartialPressurePartialDensity(
    boost::shared_ptr<pdat::CellData<double> >& data_partial_pressure_partial_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const std::vector<const double*>& thermo_properties,
    const hier::Box& domain) const
{
    NULL_USE(thermo_properties);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_partial_pressure_partial_density);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_partial_pressure_partial_density->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_density->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_pressure->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_partial_pressure_partial_density = data_partial_pressure_partial_density->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_partial_pressure_partial_density =
        data_partial_pressure_partial_density->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_density =
        data_density->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_pressure =
        data_pressure->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_partial_pressure_partial_density;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_partial_pressure_partial_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
    /*
     * Get the pointers to the cell data.
     */
    
    double* const Psi = data_partial_pressure_partial_density->getPointer(0);
    const double* const rho = data_density->getPointer(0);
    const double* const p = data_pressure->getPointer(0);
    
    computePartialPressurePartialDensity(
        Psi,
        rho,
        p,
        num_ghosts_partial_pressure_partial_density,
        num_ghosts_density,
        num_ghosts_pressure,
        ghostcell_dims_partial_pressure_partial_density,
        ghostcell_dims_density,
        ghostcell_dims_pressure,
        domain_lo,
        domain_dims);
}


/*
 * Compute the partial derivative of pressure w.r.t. density under constant specific internal energy.
 */
void
EquationOfStateIdealGas::computePartialPressurePartialDensity(
    boost::shared_ptr<pdat::SideData<double> >& data_partial_pressure_partial_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const std::vector<const double*>& thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
    NULL_USE(thermo_properties);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_partial_pressure_partial_density);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_partial_pressure_partial_density->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_density->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_pressure->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_partial_pressure_partial_density = data_partial_pressure_partial_density->getGhostCellWidth();
    hier::IntVector ghostcell_dims_partial_pressure_partial_density = data_partial_pressure_partial_density->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
    hier::IntVector ghostcell_dims_density = data_density->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
    hier::IntVector ghostcell_dims_pressure = data_pressure->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_partial_pressure_partial_density;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_partial_pressure_partial_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(side_normal < d_dim.getValue());
    
    TBOX_ASSERT(data_partial_pressure_partial_density->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_density->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_pressure->getDirectionVector()[side_normal] > 0);
#endif
    
    ghostcell_dims_partial_pressure_partial_density[side_normal]++;
    ghostcell_dims_density[side_normal]++;
    ghostcell_dims_pressure[side_normal]++;
    domain_dims[side_normal]++;
    
    /*
     * Get the pointers to the cell data.
     */
    
    double* const Psi = data_partial_pressure_partial_density->getPointer(side_normal, 0);
    const double* const rho = data_density->getPointer(side_normal, 0);
    const double* const p = data_pressure->getPointer(side_normal, 0);
    
    computePartialPressurePartialDensity(
        Psi,
        rho,
        p,
        num_ghosts_partial_pressure_partial_density,
        num_ghosts_density,
        num_ghosts_pressure,
        ghostcell_dims_partial_pressure_partial_density,
        ghostcell_dims_density,
        ghostcell_dims_pressure,
        domain_lo,
        domain_dims);
}


/*
 * Compute the partial derivative of pressure w.r.t. density under constant specific internal energy.
 */
void
EquationOfStateIdealGas::computePartialPressurePartialDensity(
    boost::shared_ptr<pdat::CellData<double> >& data_partial_pressure_partial_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties,
    const hier::Box& domain) const
{
    NULL_USE(data_thermo_properties);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_partial_pressure_partial_density);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_partial_pressure_partial_density->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_density->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_pressure->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_partial_pressure_partial_density = data_partial_pressure_partial_density->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_partial_pressure_partial_density =
        data_partial_pressure_partial_density->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_density =
        data_density->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_pressure =
        data_pressure->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_partial_pressure_partial_density;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_partial_pressure_partial_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
    /*
     * Get the pointers to the cell data.
     */
    
    double* const Psi = data_partial_pressure_partial_density->getPointer(0);
    const double* const rho = data_density->getPointer(0);
    const double* const p = data_pressure->getPointer(0);
    
    computePartialPressurePartialDensity(
        Psi,
        rho,
        p,
        num_ghosts_partial_pressure_partial_density,
        num_ghosts_density,
        num_ghosts_pressure,
        ghostcell_dims_partial_pressure_partial_density,
        ghostcell_dims_density,
        ghostcell_dims_pressure,
        domain_lo,
        domain_dims);
}


/*
 * Compute the partial derivative of pressure w.r.t. density under constant specific internal energy.
 */
void
EquationOfStateIdealGas::computePartialPressurePartialDensity(
    boost::shared_ptr<pdat::SideData<double> >& data_partial_pressure_partial_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const boost::shared_ptr<pdat::SideData<double> >& data_thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
    NULL_USE(data_thermo_properties);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_partial_pressure_partial_density);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_partial_pressure_partial_density->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_density->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_pressure->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_partial_pressure_partial_density = data_partial_pressure_partial_density->getGhostCellWidth();
    hier::IntVector ghostcell_dims_partial_pressure_partial_density = data_partial_pressure_partial_density->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
    hier::IntVector ghostcell_dims_density = data_density->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
    hier::IntVector ghostcell_dims_pressure = data_pressure->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_partial_pressure_partial_density;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_partial_pressure_partial_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(side_normal < d_dim.getValue());
    
    TBOX_ASSERT(data_partial_pressure_partial_density->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_density->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_pressure->getDirectionVector()[side_normal] > 0);
#endif
    
    ghostcell_dims_partial_pressure_partial_density[side_normal]++;
    ghostcell_dims_density[side_normal]++;
    ghostcell_dims_pressure[side_normal]++;
    domain_dims[side_normal]++;
    
    /*
     * Get the pointers to the cell data.
     */
    
    double* const Psi = data_partial_pressure_partial_density->getPointer(side_normal, 0);
    const double* const rho = data_density->getPointer(side_normal, 0);
    const double* const p = data_pressure->getPointer(side_normal, 0);
    
    computePartialPressurePartialDensity(
        Psi,
        rho,
        p,
        num_ghosts_partial_pressure_partial_density,
        num_ghosts_density,
        num_ghosts_pressure,
        ghostcell_dims_partial_pressure_partial_density,
        ghostcell_dims_density,
        ghostcell_dims_pressure,
        domain_lo,
        domain_dims);
}


////////////////////////////////////////////////////////////////////////////////////////////////////


/*
 * Compute the Gruneisen parameter (partial derivative of pressure w.r.t. specific internal energy under
 * constant density divided by density).
 */
double
EquationOfStateIdealGas::getGruneisenParameter(
    const double* const density,
    const double* const pressure,
    const std::vector<const double*>& thermo_properties) const
{
    NULL_USE(density);
    NULL_USE(pressure);
    
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 1);
#endif
    
    const double& gamma = *(thermo_properties[0]);
    
    return (gamma - double(1));
}


/*
 * Compute the Gruneisen parameter (partial derivative of pressure w.r.t. specific internal energy under
 * constant density divided by density).
 */
void
EquationOfStateIdealGas::computeGruneisenParameter(
    boost::shared_ptr<pdat::CellData<double> >& data_gruneisen_parameter,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const std::vector<const double*>& thermo_properties,
    const hier::Box& domain) const
{
    NULL_USE(data_density);
    NULL_USE(data_pressure);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_gruneisen_parameter);
    
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 1);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_gruneisen_parameter->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell box.
     */
    
    const hier::IntVector num_ghosts_gruneisen_parameter =
        data_gruneisen_parameter->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_gruneisen_parameter =
        data_gruneisen_parameter->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_gruneisen_parameter);
        
        domain_lo = -num_ghosts_gruneisen_parameter;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_gruneisen_parameter->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
    /*
     * Get the pointer to the cell data.
     */
    
    double* const Gamma = data_gruneisen_parameter->getPointer(0);
    
    const double& gamma = *(thermo_properties[0]);
    
    computeGruneisenParameter(
        Gamma,
        gamma,
        num_ghosts_gruneisen_parameter,
        ghostcell_dims_gruneisen_parameter,
        domain_lo,
        domain_dims);
}


/*
 * Compute the Gruneisen parameter (partial derivative of pressure w.r.t. specific internal energy under
 * constant density divided by density).
 */
void
EquationOfStateIdealGas::computeGruneisenParameter(
    boost::shared_ptr<pdat::SideData<double> >& data_gruneisen_parameter,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const std::vector<const double*>& thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
    NULL_USE(data_density);
    NULL_USE(data_pressure);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_gruneisen_parameter);
    
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 1);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_gruneisen_parameter->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell box.
     */
    
    const hier::IntVector num_ghosts_gruneisen_parameter =
        data_gruneisen_parameter->getGhostCellWidth();
    hier::IntVector ghostcell_dims_gruneisen_parameter =
        data_gruneisen_parameter->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_gruneisen_parameter);
        
        domain_lo = -num_ghosts_gruneisen_parameter;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_gruneisen_parameter->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(side_normal < d_dim.getValue());
    
    TBOX_ASSERT(data_gruneisen_parameter->getDirectionVector()[side_normal] > 0);
#endif
    
    ghostcell_dims_gruneisen_parameter[side_normal]++;
    domain_dims[side_normal]++;
    
    /*
     * Get the pointer to the cell data.
     */
    
    double* const Gamma = data_gruneisen_parameter->getPointer(side_normal, 0);
    
    const double& gamma = *(thermo_properties[0]);
    
    computeGruneisenParameter(
        Gamma,
        gamma,
        num_ghosts_gruneisen_parameter,
        ghostcell_dims_gruneisen_parameter,
        domain_lo,
        domain_dims);
}


/*
 * Compute the Gruneisen parameter (partial derivative of pressure w.r.t. specific internal energy under
 * constant density divided by density).
 */
void
EquationOfStateIdealGas::computeGruneisenParameter(
    boost::shared_ptr<pdat::CellData<double> >& data_gruneisen_parameter,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties,
    const hier::Box& domain) const
{
    NULL_USE(data_density);
    NULL_USE(data_pressure);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_gruneisen_parameter);
    TBOX_ASSERT(data_thermo_properties);
    
    TBOX_ASSERT(data_thermo_properties->getDepth() >= 1);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_gruneisen_parameter->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_thermo_properties->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_gruneisen_parameter =
        data_gruneisen_parameter->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_gruneisen_parameter =
        data_gruneisen_parameter->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_thermo_properties = data_thermo_properties->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_thermo_properties =
        data_thermo_properties->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_gruneisen_parameter;
        num_ghosts_min = hier::IntVector::min(num_ghosts_thermo_properties, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_gruneisen_parameter->getGhostBox().contains(domain));
        TBOX_ASSERT(data_thermo_properties->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
    /*
     * Get the pointers to the cell data.
     */
    
    double* const Gamma = data_gruneisen_parameter->getPointer(0);
    const double* const gamma = data_thermo_properties->getPointer(0);
    
    computeGruneisenParameter(
        Gamma,
        gamma,
        num_ghosts_gruneisen_parameter,
        num_ghosts_thermo_properties,
        ghostcell_dims_gruneisen_parameter,
        ghostcell_dims_thermo_properties,
        domain_lo,
        domain_dims);
}


/*
 * Compute the Gruneisen parameter (partial derivative of pressure w.r.t. specific internal energy under
 * constant density divided by density).
 */
void
EquationOfStateIdealGas::computeGruneisenParameter(
    boost::shared_ptr<pdat::SideData<double> >& data_gruneisen_parameter,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const boost::shared_ptr<pdat::SideData<double> >& data_thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
    NULL_USE(data_density);
    NULL_USE(data_pressure);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_gruneisen_parameter);
    TBOX_ASSERT(data_thermo_properties);
    
    TBOX_ASSERT(data_thermo_properties->getDepth() >= 1);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_gruneisen_parameter->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_thermo_properties->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_gruneisen_parameter =
        data_gruneisen_parameter->getGhostCellWidth();
    hier::IntVector ghostcell_dims_gruneisen_parameter =
        data_gruneisen_parameter->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_thermo_properties = data_thermo_properties->getGhostCellWidth();
    hier::IntVector ghostcell_dims_thermo_properties = data_thermo_properties->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_gruneisen_parameter;
        num_ghosts_min = hier::IntVector::min(num_ghosts_thermo_properties, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_gruneisen_parameter->getGhostBox().contains(domain));
        TBOX_ASSERT(data_thermo_properties->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(side_normal < d_dim.getValue());
    
    TBOX_ASSERT(data_gruneisen_parameter->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_thermo_properties->getDirectionVector()[side_normal] > 0);
#endif
    
    ghostcell_dims_gruneisen_parameter[side_normal]++;
    ghostcell_dims_thermo_properties[side_normal]++;
    domain_dims[side_normal]++;
    
    /*
     * Get the pointers to the cell data.
     */
    
    double* const Gamma = data_gruneisen_parameter->getPointer(side_normal, 0);
    const double* const gamma = data_thermo_properties->getPointer(side_normal, 0);
    
    computeGruneisenParameter(
        Gamma,
        gamma,
        num_ghosts_gruneisen_parameter,
        num_ghosts_thermo_properties,
        ghostcell_dims_gruneisen_parameter,
        ghostcell_dims_thermo_properties,
        domain_lo,
        domain_dims);
}


////////////////////////////////////////////////////////////////////////////////////////////////////


/*
 * Compute the partial derivative of internal energy w.r.t. pressure under constant density.
 */
double
EquationOfStateIdealGas::getIsochoricPartialInternalEnergyPartialPressure(
    const double* const density,
    const double* const pressure,
    const std::vector<const double*>& thermo_properties) const
{
    NULL_USE(density);
    NULL_USE(pressure);
    
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 1);
#endif
    
    const double& gamma = *(thermo_properties[0]);
    
    return double(1)/(gamma - double(1)); // Return xi.
}


/*
 * Compute the partial derivative of internal energy w.r.t. pressure under constant density.
 */
void
EquationOfStateIdealGas::computeIsochoricPartialInternalEnergyPartialPressure(
    boost::shared_ptr<pdat::CellData<double> >& data_partial_internal_energy_partial_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const std::vector<const double*>& thermo_properties,
    const hier::Box& domain) const
{
    NULL_USE(data_density);
    NULL_USE(data_pressure);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_partial_internal_energy_partial_pressure);
    
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 1);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_partial_internal_energy_partial_pressure->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell box.
     */
    
    const hier::IntVector num_ghosts_partial_internal_energy_partial_pressure =
        data_partial_internal_energy_partial_pressure->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_partial_internal_energy_partial_pressure =
        data_partial_internal_energy_partial_pressure->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_partial_internal_energy_partial_pressure);
        
        domain_lo = -num_ghosts_partial_internal_energy_partial_pressure;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_partial_internal_energy_partial_pressure->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
    /*
     * Get the pointer to the cell data.
     */
    
    double* const xi = data_partial_internal_energy_partial_pressure->getPointer(0);
    
    const double& gamma = *(thermo_properties[0]);
    
    computeIsochoricPartialInternalEnergyPartialPressure(
        xi,
        gamma,
        num_ghosts_partial_internal_energy_partial_pressure,
        ghostcell_dims_partial_internal_energy_partial_pressure,
        domain_lo,
        domain_dims);
}


/*
 * Compute the partial derivative of internal energy w.r.t. pressure under constant density.
 */
void
EquationOfStateIdealGas::computeIsochoricPartialInternalEnergyPartialPressure(
    boost::shared_ptr<pdat::SideData<double> >& data_partial_internal_energy_partial_pressure,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const std::vector<const double*>& thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
    NULL_USE(data_density);
    NULL_USE(data_pressure);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_partial_internal_energy_partial_pressure);
    
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 1);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_partial_internal_energy_partial_pressure->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell box.
     */
    
    const hier::IntVector num_ghosts_partial_internal_energy_partial_pressure =
        data_partial_internal_energy_partial_pressure->getGhostCellWidth();
    hier::IntVector ghostcell_dims_partial_internal_energy_partial_pressure =
        data_partial_internal_energy_partial_pressure->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_partial_internal_energy_partial_pressure);
        
        domain_lo = -num_ghosts_partial_internal_energy_partial_pressure;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_partial_internal_energy_partial_pressure->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(side_normal < d_dim.getValue());
    
    TBOX_ASSERT(data_partial_internal_energy_partial_pressure->getDirectionVector()[side_normal] > 0);
#endif
    
    ghostcell_dims_partial_internal_energy_partial_pressure[side_normal]++;
    domain_dims[side_normal]++;
    
    /*
     * Get the pointer to the cell data.
     */
    
    double* const xi = data_partial_internal_energy_partial_pressure->getPointer(side_normal, 0);
    
    const double& gamma = *(thermo_properties[0]);
    
    computeIsochoricPartialInternalEnergyPartialPressure(
        xi,
        gamma,
        num_ghosts_partial_internal_energy_partial_pressure,
        ghostcell_dims_partial_internal_energy_partial_pressure,
        domain_lo,
        domain_dims);
}


/*
 * Compute the partial derivative of internal energy w.r.t. pressure under constant density.
 */
void
EquationOfStateIdealGas::computeIsochoricPartialInternalEnergyPartialPressure(
    boost::shared_ptr<pdat::CellData<double> >& data_partial_internal_energy_partial_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties,
    const hier::Box& domain) const
{
    NULL_USE(data_density);
    NULL_USE(data_pressure);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_partial_internal_energy_partial_pressure);
    TBOX_ASSERT(data_thermo_properties);
    
    TBOX_ASSERT(data_thermo_properties->getDepth() >= 1);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_partial_internal_energy_partial_pressure->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_thermo_properties->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_partial_internal_energy_partial_pressure =
        data_partial_internal_energy_partial_pressure->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_partial_internal_energy_partial_pressure =
        data_partial_internal_energy_partial_pressure->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_thermo_properties = data_thermo_properties->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_thermo_properties =
        data_thermo_properties->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_partial_internal_energy_partial_pressure;
        num_ghosts_min = hier::IntVector::min(num_ghosts_thermo_properties, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_partial_internal_energy_partial_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_thermo_properties->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
    /*
     * Get the pointers to the cell data.
     */
    
    double* const xi = data_partial_internal_energy_partial_pressure->getPointer(0);
    const double* const gamma = data_thermo_properties->getPointer(0);
    
    computeIsochoricPartialInternalEnergyPartialPressure(
        xi,
        gamma,
        num_ghosts_partial_internal_energy_partial_pressure,
        num_ghosts_thermo_properties,
        ghostcell_dims_partial_internal_energy_partial_pressure,
        ghostcell_dims_thermo_properties,
        domain_lo,
        domain_dims);
}


/*
 * Compute the partial derivative of internal energy w.r.t. pressure under constant density.
 */
void
EquationOfStateIdealGas::computeIsochoricPartialInternalEnergyPartialPressure(
    boost::shared_ptr<pdat::SideData<double> >& data_partial_internal_energy_partial_pressure,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const boost::shared_ptr<pdat::SideData<double> >& data_thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
    NULL_USE(data_density);
    NULL_USE(data_pressure);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_partial_internal_energy_partial_pressure);
    TBOX_ASSERT(data_thermo_properties);
    
    TBOX_ASSERT(data_thermo_properties->getDepth() >= 1);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_partial_internal_energy_partial_pressure->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_thermo_properties->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_partial_internal_energy_partial_pressure =
        data_partial_internal_energy_partial_pressure->getGhostCellWidth();
    hier::IntVector ghostcell_dims_partial_internal_energy_partial_pressure =
        data_partial_internal_energy_partial_pressure->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_thermo_properties = data_thermo_properties->getGhostCellWidth();
    hier::IntVector ghostcell_dims_thermo_properties = data_thermo_properties->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_partial_internal_energy_partial_pressure;
        num_ghosts_min = hier::IntVector::min(num_ghosts_thermo_properties, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_partial_internal_energy_partial_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_thermo_properties->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(side_normal < d_dim.getValue());
    
    TBOX_ASSERT(data_partial_internal_energy_partial_pressure->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_thermo_properties->getDirectionVector()[side_normal] > 0);
#endif
    
    ghostcell_dims_partial_internal_energy_partial_pressure[side_normal]++;
    ghostcell_dims_thermo_properties[side_normal]++;
    domain_dims[side_normal]++;
    
    /*
     * Get the pointers to the cell data.
     */
    
    double* const xi = data_partial_internal_energy_partial_pressure->getPointer(side_normal, 0);
    const double* const gamma = data_thermo_properties->getPointer(side_normal, 0);
    
    computeIsochoricPartialInternalEnergyPartialPressure(
        xi,
        gamma,
        num_ghosts_partial_internal_energy_partial_pressure,
        num_ghosts_thermo_properties,
        ghostcell_dims_partial_internal_energy_partial_pressure,
        ghostcell_dims_thermo_properties,
        domain_lo,
        domain_dims);
}


/*
 * Compute the partial derivative of internal energy w.r.t. density under constant pressure.
 */
double
EquationOfStateIdealGas::getIsobaricPartialInternalEnergyPartialDensity(
    const double* const density,
    const double* const pressure,
    const std::vector<const double*>& thermo_properties) const
{
    NULL_USE(density);
    NULL_USE(pressure);
    NULL_USE(thermo_properties);
    
    return double(0); // Return delta.
}


/*
 * Compute the partial derivative of internal energy w.r.t. density under constant pressure.
 */
void
EquationOfStateIdealGas::computeIsobaricPartialInternalEnergyPartialDensity(
    boost::shared_ptr<pdat::CellData<double> >& data_partial_internal_energy_partial_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const std::vector<const double*>& thermo_properties,
    const hier::Box& domain) const
{
    NULL_USE(data_density);
    NULL_USE(data_pressure);
    NULL_USE(thermo_properties);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_partial_internal_energy_partial_density);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_partial_internal_energy_partial_density->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell box.
     */
    
    const hier::IntVector num_ghosts_partial_internal_energy_partial_density =
        data_partial_internal_energy_partial_density->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_partial_internal_energy_partial_density =
        data_partial_internal_energy_partial_density->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_partial_internal_energy_partial_density);
        
        domain_lo = -num_ghosts_partial_internal_energy_partial_density;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_partial_internal_energy_partial_density->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
    /*
     * Get the pointer to the cell data.
     */
    
    double* const delta = data_partial_internal_energy_partial_density->getPointer(0);
    
    computeIsobaricPartialInternalEnergyPartialDensity(
        delta,
        num_ghosts_partial_internal_energy_partial_density,
        ghostcell_dims_partial_internal_energy_partial_density,
        domain_lo,
        domain_dims);
}


/*
 * Compute the partial derivative of internal energy w.r.t. density under constant pressure.
 */
void
EquationOfStateIdealGas::computeIsobaricPartialInternalEnergyPartialDensity(
    boost::shared_ptr<pdat::SideData<double> >& data_partial_internal_energy_partial_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const std::vector<const double*>& thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
    NULL_USE(data_density);
    NULL_USE(data_pressure);
    NULL_USE(thermo_properties);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_partial_internal_energy_partial_density);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_partial_internal_energy_partial_density->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell box.
     */
    
    const hier::IntVector num_ghosts_partial_internal_energy_partial_density =
        data_partial_internal_energy_partial_density->getGhostCellWidth();
    hier::IntVector ghostcell_dims_partial_internal_energy_partial_density =
        data_partial_internal_energy_partial_density->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_partial_internal_energy_partial_density);
        
        domain_lo = -num_ghosts_partial_internal_energy_partial_density;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_partial_internal_energy_partial_density->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(side_normal < d_dim.getValue());
    
    TBOX_ASSERT(data_partial_internal_energy_partial_density->getDirectionVector()[side_normal] > 0);
#endif
    
    ghostcell_dims_partial_internal_energy_partial_density[side_normal]++;
    domain_dims[side_normal]++;
    
    /*
     * Get the pointer to the cell data.
     */
    
    double* const delta = data_partial_internal_energy_partial_density->getPointer(side_normal, 0);
    
    computeIsobaricPartialInternalEnergyPartialDensity(
        delta,
        num_ghosts_partial_internal_energy_partial_density,
        ghostcell_dims_partial_internal_energy_partial_density,
        domain_lo,
        domain_dims);
}


/*
 * Compute the partial derivative of internal energy w.r.t. density under constant pressure.
 */
void
EquationOfStateIdealGas::computeIsobaricPartialInternalEnergyPartialDensity(
    boost::shared_ptr<pdat::CellData<double> >& data_partial_internal_energy_partial_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties,
    const hier::Box& domain) const
{
    NULL_USE(data_density);
    NULL_USE(data_pressure);
    NULL_USE(data_thermo_properties);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_partial_internal_energy_partial_density);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_partial_internal_energy_partial_density->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell box.
     */
    
    const hier::IntVector num_ghosts_partial_internal_energy_partial_density =
        data_partial_internal_energy_partial_density->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_partial_internal_energy_partial_density =
        data_partial_internal_energy_partial_density->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_partial_internal_energy_partial_density);
        
        domain_lo = -num_ghosts_partial_internal_energy_partial_density;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_partial_internal_energy_partial_density->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
    /*
     * Get the pointer to the cell data.
     */
    
    double* const delta = data_partial_internal_energy_partial_density->getPointer(0);
    
    computeIsobaricPartialInternalEnergyPartialDensity(
        delta,
        num_ghosts_partial_internal_energy_partial_density,
        ghostcell_dims_partial_internal_energy_partial_density,
        domain_lo,
        domain_dims);
}


/*
 * Compute the partial derivative of internal energy w.r.t. density under constant pressure.
 */
void
EquationOfStateIdealGas::computeIsobaricPartialInternalEnergyPartialDensity(
    boost::shared_ptr<pdat::SideData<double> >& data_partial_internal_energy_partial_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const boost::shared_ptr<pdat::SideData<double> >& data_thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
    NULL_USE(data_density);
    NULL_USE(data_pressure);
    NULL_USE(data_thermo_properties);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_partial_internal_energy_partial_density);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_partial_internal_energy_partial_density->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell box.
     */
    
    const hier::IntVector num_ghosts_partial_internal_energy_partial_density =
        data_partial_internal_energy_partial_density->getGhostCellWidth();
    hier::IntVector ghostcell_dims_partial_internal_energy_partial_density =
        data_partial_internal_energy_partial_density->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_partial_internal_energy_partial_density);
        
        domain_lo = -num_ghosts_partial_internal_energy_partial_density;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_partial_internal_energy_partial_density->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(side_normal < d_dim.getValue());
    
    TBOX_ASSERT(data_partial_internal_energy_partial_density->getDirectionVector()[side_normal] > 0);
#endif
    
    ghostcell_dims_partial_internal_energy_partial_density[side_normal]++;
    domain_dims[side_normal]++;
    
    /*
     * Get the pointer to the cell data.
     */
    
    double* const delta = data_partial_internal_energy_partial_density->getPointer(side_normal, 0);
    
    computeIsobaricPartialInternalEnergyPartialDensity(
        delta,
        num_ghosts_partial_internal_energy_partial_density,
        ghostcell_dims_partial_internal_energy_partial_density,
        domain_lo,
        domain_dims);
}


/*
 * Compute the density.
 */
double
EquationOfStateIdealGas::getDensity(
    const double* const pressure,
    const double* const temperature,
    const std::vector<const double*>& thermo_properties) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 2);
#endif
    
    const double& R = *(thermo_properties[1]);
    
    const double& p = *pressure;
    const double& T = *temperature;
    
    return p/(R*T); // Return rho.
}


/*
 * Compute the density.
 */
void
EquationOfStateIdealGas::computeDensity(
    boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
    const std::vector<const double*>& thermo_properties,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_temperature);
    
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 2);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_density->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_pressure->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_temperature->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_density =
        data_density->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_pressure =
        data_pressure->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_temperature =
        data_temperature->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_density;
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_temperature, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_temperature->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
    /*
     * Get the pointers to the cell data.
     */
    
    double* const rho = data_density->getPointer(0);
    const double* const p = data_pressure->getPointer(0);
    const double* const T = data_temperature->getPointer(0);
    
    const double& R = *(thermo_properties[1]);
    
    computeDensity(
        rho,
        p,
        T,
        R,
        num_ghosts_density,
        num_ghosts_pressure,
        num_ghosts_temperature,
        ghostcell_dims_density,
        ghostcell_dims_pressure,
        ghostcell_dims_temperature,
        domain_lo,
        domain_dims);
}


/*
 * Compute the density.
 */
void
EquationOfStateIdealGas::computeDensity(
    boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const boost::shared_ptr<pdat::SideData<double> >& data_temperature,
    const std::vector<const double*>& thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_temperature);
    
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 2);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_density->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_pressure->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_temperature->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
    hier::IntVector ghostcell_dims_density = data_density->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
    hier::IntVector ghostcell_dims_pressure = data_pressure->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
    hier::IntVector ghostcell_dims_temperature = data_temperature->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_density;
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_temperature, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_temperature->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(side_normal < d_dim.getValue());
    
    TBOX_ASSERT(data_density->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_pressure->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_temperature->getDirectionVector()[side_normal] > 0);
#endif
    
    ghostcell_dims_density[side_normal]++;
    ghostcell_dims_pressure[side_normal]++;
    ghostcell_dims_temperature[side_normal]++;
    domain_dims[side_normal]++;
    
    /*
     * Get the pointers to the cell data.
     */
    
    double* const rho = data_density->getPointer(side_normal, 0);
    const double* const p = data_pressure->getPointer(side_normal, 0);
    const double* const T = data_temperature->getPointer(side_normal, 0);
    
    const double& R = *(thermo_properties[1]);
    
    computeDensity(
        rho,
        p,
        T,
        R,
        num_ghosts_density,
        num_ghosts_pressure,
        num_ghosts_temperature,
        ghostcell_dims_density,
        ghostcell_dims_pressure,
        ghostcell_dims_temperature,
        domain_lo,
        domain_dims);
}


/*
 * Compute the density.
 */
void
EquationOfStateIdealGas::computeDensity(
    boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
    const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_temperature);
    TBOX_ASSERT(data_thermo_properties);
    
    TBOX_ASSERT(data_thermo_properties->getDepth() >= 2);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_density->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_pressure->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_temperature->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_thermo_properties->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_density =
        data_density->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_pressure =
        data_pressure->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_temperature =
        data_temperature->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_thermo_properties = data_thermo_properties->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_thermo_properties =
        data_thermo_properties->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_density;
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_temperature, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_thermo_properties, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_temperature->getGhostBox().contains(domain));
        TBOX_ASSERT(data_thermo_properties->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
    /*
     * Get the pointers to the cell data.
     */
    
    double* const rho = data_density->getPointer(0);
    const double* const p = data_pressure->getPointer(0);
    const double* const T = data_temperature->getPointer(0);
    const double* const R = data_thermo_properties->getPointer(1);
    
    computeDensity(
        rho,
        p,
        T,
        R,
        num_ghosts_density,
        num_ghosts_pressure,
        num_ghosts_temperature,
        num_ghosts_thermo_properties,
        ghostcell_dims_density,
        ghostcell_dims_pressure,
        ghostcell_dims_temperature,
        ghostcell_dims_thermo_properties,
        domain_lo,
        domain_dims);
}


/*
 * Compute the density.
 */
void
EquationOfStateIdealGas::computeDensity(
    boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const boost::shared_ptr<pdat::SideData<double> >& data_temperature,
    const boost::shared_ptr<pdat::SideData<double> >& data_thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_temperature);
    TBOX_ASSERT(data_thermo_properties);
    
    TBOX_ASSERT(data_thermo_properties->getDepth() >= 2);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_density->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_pressure->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_temperature->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_thermo_properties->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
    hier::IntVector ghostcell_dims_density = data_density->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
    hier::IntVector ghostcell_dims_pressure = data_pressure->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
    hier::IntVector ghostcell_dims_temperature = data_temperature->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_thermo_properties = data_thermo_properties->getGhostCellWidth();
    hier::IntVector ghostcell_dims_thermo_properties = data_thermo_properties->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_density;
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_temperature, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_thermo_properties, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_temperature->getGhostBox().contains(domain));
        TBOX_ASSERT(data_thermo_properties->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(side_normal < d_dim.getValue());
    
    TBOX_ASSERT(data_density->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_pressure->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_temperature->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_thermo_properties->getDirectionVector()[side_normal] > 0);
#endif
    
    ghostcell_dims_density[side_normal]++;
    ghostcell_dims_pressure[side_normal]++;
    ghostcell_dims_temperature[side_normal]++;
    ghostcell_dims_thermo_properties[side_normal]++;
    domain_dims[side_normal]++;
    
    /*
     * Get the pointers to the cell data.
     */
    
    double* const rho = data_density->getPointer(side_normal, 0);
    const double* const p = data_pressure->getPointer(side_normal, 0);
    const double* const T = data_temperature->getPointer(side_normal, 0);
    const double* const R = data_thermo_properties->getPointer(side_normal, 1);
    
    computeDensity(
        rho,
        p,
        T,
        R,
        num_ghosts_density,
        num_ghosts_pressure,
        num_ghosts_temperature,
        num_ghosts_thermo_properties,
        ghostcell_dims_density,
        ghostcell_dims_pressure,
        ghostcell_dims_temperature,
        ghostcell_dims_thermo_properties,
        domain_lo,
        domain_dims);
}


/*
 * Compute the pressure.
 */
void
EquationOfStateIdealGas::computePressure(
    double* const p,
    const double* const rho,
    const double* const epsilon,
    const double& gamma,
    const hier::IntVector& num_ghosts_pressure,
    const hier::IntVector& num_ghosts_density,
    const hier::IntVector& num_ghosts_internal_energy,
    const hier::IntVector& ghostcell_dims_pressure,
    const hier::IntVector& ghostcell_dims_density,
    const hier::IntVector& ghostcell_dims_internal_energy,
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
        
        const int num_ghosts_0_pressure = num_ghosts_pressure[0];
        const int num_ghosts_0_density = num_ghosts_density[0];
        const int num_ghosts_0_internal_energy = num_ghosts_internal_energy[0];
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_pressure = i + num_ghosts_0_pressure;
            const int idx_density = i + num_ghosts_0_density;
            const int idx_internal_energy = i + num_ghosts_0_internal_energy;
            
            p[idx_pressure] = (gamma - double(1))*rho[idx_density]*epsilon[idx_internal_energy];
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
        
        const int num_ghosts_0_pressure = num_ghosts_pressure[0];
        const int num_ghosts_1_pressure = num_ghosts_pressure[1];
        const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
        
        const int num_ghosts_0_density = num_ghosts_density[0];
        const int num_ghosts_1_density = num_ghosts_density[1];
        const int ghostcell_dim_0_density = ghostcell_dims_density[0];
        
        const int num_ghosts_0_internal_energy = num_ghosts_internal_energy[0];
        const int num_ghosts_1_internal_energy = num_ghosts_internal_energy[1];
        const int ghostcell_dim_0_internal_energy = ghostcell_dims_internal_energy[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_pressure = (i + num_ghosts_0_pressure) +
                    (j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure;
                
                const int idx_density = (i + num_ghosts_0_density) +
                    (j + num_ghosts_1_density)*ghostcell_dim_0_density;
                
                const int idx_internal_energy = (i + num_ghosts_0_internal_energy) +
                    (j + num_ghosts_1_internal_energy)*ghostcell_dim_0_internal_energy;
                
                p[idx_pressure] = (gamma - double(1))*rho[idx_density]*epsilon[idx_internal_energy];
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
        
        const int num_ghosts_0_pressure = num_ghosts_pressure[0];
        const int num_ghosts_1_pressure = num_ghosts_pressure[1];
        const int num_ghosts_2_pressure = num_ghosts_pressure[2];
        const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
        const int ghostcell_dim_1_pressure = ghostcell_dims_pressure[1];
        
        const int num_ghosts_0_density = num_ghosts_density[0];
        const int num_ghosts_1_density = num_ghosts_density[1];
        const int num_ghosts_2_density = num_ghosts_density[2];
        const int ghostcell_dim_0_density = ghostcell_dims_density[0];
        const int ghostcell_dim_1_density = ghostcell_dims_density[1];
        
        const int num_ghosts_0_internal_energy = num_ghosts_internal_energy[0];
        const int num_ghosts_1_internal_energy = num_ghosts_internal_energy[1];
        const int num_ghosts_2_internal_energy = num_ghosts_internal_energy[2];
        const int ghostcell_dim_0_internal_energy = ghostcell_dims_internal_energy[0];
        const int ghostcell_dim_1_internal_energy = ghostcell_dims_internal_energy[1];
        
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
                    const int idx_pressure = (i + num_ghosts_0_pressure) +
                        (j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure +
                        (k + num_ghosts_2_pressure)*ghostcell_dim_0_pressure*
                            ghostcell_dim_1_pressure;
                    
                    const int idx_density = (i + num_ghosts_0_density) +
                        (j + num_ghosts_1_density)*ghostcell_dim_0_density +
                        (k + num_ghosts_2_density)*ghostcell_dim_0_density*
                            ghostcell_dim_1_density;
                    
                    const int idx_internal_energy = (i + num_ghosts_0_internal_energy) +
                        (j + num_ghosts_1_internal_energy)*ghostcell_dim_0_internal_energy +
                        (k + num_ghosts_2_internal_energy)*ghostcell_dim_0_internal_energy*
                            ghostcell_dim_1_internal_energy;
                    
                    p[idx_pressure] = (gamma - double(1))*rho[idx_density]*epsilon[idx_internal_energy];
                }
            }
        }
    }
}


/*
 * Compute the pressure.
 */
void
EquationOfStateIdealGas::computePressure(
    double* const p,
    const double* const rho,
    const double* const epsilon,
    const double* const gamma,
    const hier::IntVector& num_ghosts_pressure,
    const hier::IntVector& num_ghosts_density,
    const hier::IntVector& num_ghosts_internal_energy,
    const hier::IntVector& num_ghosts_thermo_properties,
    const hier::IntVector& ghostcell_dims_pressure,
    const hier::IntVector& ghostcell_dims_density,
    const hier::IntVector& ghostcell_dims_internal_energy,
    const hier::IntVector& ghostcell_dims_thermo_properties,
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
        
        const int num_ghosts_0_pressure = num_ghosts_pressure[0];
        const int num_ghosts_0_density = num_ghosts_density[0];
        const int num_ghosts_0_internal_energy = num_ghosts_internal_energy[0];
        const int num_ghosts_0_thermo_properties = num_ghosts_thermo_properties[0];
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_pressure = i + num_ghosts_0_pressure;
            const int idx_density = i + num_ghosts_0_density;
            const int idx_internal_energy = i + num_ghosts_0_internal_energy;
            const int idx_thermo_properties = i + num_ghosts_0_thermo_properties;
            
            p[idx_pressure] = (gamma[idx_thermo_properties] - double(1))*rho[idx_density]*
                epsilon[idx_internal_energy];
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
        
        const int num_ghosts_0_pressure = num_ghosts_pressure[0];
        const int num_ghosts_1_pressure = num_ghosts_pressure[1];
        const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
        
        const int num_ghosts_0_density = num_ghosts_density[0];
        const int num_ghosts_1_density = num_ghosts_density[1];
        const int ghostcell_dim_0_density = ghostcell_dims_density[0];
        
        const int num_ghosts_0_internal_energy = num_ghosts_internal_energy[0];
        const int num_ghosts_1_internal_energy = num_ghosts_internal_energy[1];
        const int ghostcell_dim_0_internal_energy = ghostcell_dims_internal_energy[0];
        
        const int num_ghosts_0_thermo_properties = num_ghosts_thermo_properties[0];
        const int num_ghosts_1_thermo_properties = num_ghosts_thermo_properties[1];
        const int ghostcell_dim_0_thermo_properties = ghostcell_dims_thermo_properties[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_pressure = (i + num_ghosts_0_pressure) +
                    (j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure;
                
                const int idx_density = (i + num_ghosts_0_density) +
                    (j + num_ghosts_1_density)*ghostcell_dim_0_density;
                
                const int idx_internal_energy = (i + num_ghosts_0_internal_energy) +
                    (j + num_ghosts_1_internal_energy)*ghostcell_dim_0_internal_energy;
                
                const int idx_thermo_properties = (i + num_ghosts_0_thermo_properties) +
                    (j + num_ghosts_1_thermo_properties)*ghostcell_dim_0_thermo_properties;
                
                p[idx_pressure] = (gamma[idx_thermo_properties] - double(1))*rho[idx_density]*
                    epsilon[idx_internal_energy];
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
        
        const int num_ghosts_0_pressure = num_ghosts_pressure[0];
        const int num_ghosts_1_pressure = num_ghosts_pressure[1];
        const int num_ghosts_2_pressure = num_ghosts_pressure[2];
        const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
        const int ghostcell_dim_1_pressure = ghostcell_dims_pressure[1];
        
        const int num_ghosts_0_density = num_ghosts_density[0];
        const int num_ghosts_1_density = num_ghosts_density[1];
        const int num_ghosts_2_density = num_ghosts_density[2];
        const int ghostcell_dim_0_density = ghostcell_dims_density[0];
        const int ghostcell_dim_1_density = ghostcell_dims_density[1];
        
        const int num_ghosts_0_internal_energy = num_ghosts_internal_energy[0];
        const int num_ghosts_1_internal_energy = num_ghosts_internal_energy[1];
        const int num_ghosts_2_internal_energy = num_ghosts_internal_energy[2];
        const int ghostcell_dim_0_internal_energy = ghostcell_dims_internal_energy[0];
        const int ghostcell_dim_1_internal_energy = ghostcell_dims_internal_energy[1];
        
        const int num_ghosts_0_thermo_properties = num_ghosts_thermo_properties[0];
        const int num_ghosts_1_thermo_properties = num_ghosts_thermo_properties[1];
        const int num_ghosts_2_thermo_properties = num_ghosts_thermo_properties[2];
        const int ghostcell_dim_0_thermo_properties = ghostcell_dims_thermo_properties[0];
        const int ghostcell_dim_1_thermo_properties = ghostcell_dims_thermo_properties[1];
        
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
                    const int idx_pressure = (i + num_ghosts_0_pressure) +
                        (j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure +
                        (k + num_ghosts_2_pressure)*ghostcell_dim_0_pressure*
                            ghostcell_dim_1_pressure;
                    
                    const int idx_density = (i + num_ghosts_0_density) +
                        (j + num_ghosts_1_density)*ghostcell_dim_0_density +
                        (k + num_ghosts_2_density)*ghostcell_dim_0_density*
                            ghostcell_dim_1_density;
                    
                    const int idx_internal_energy = (i + num_ghosts_0_internal_energy) +
                        (j + num_ghosts_1_internal_energy)*ghostcell_dim_0_internal_energy +
                        (k + num_ghosts_2_internal_energy)*ghostcell_dim_0_internal_energy*
                            ghostcell_dim_1_internal_energy;
                    
                    const int idx_thermo_properties = (i + num_ghosts_0_thermo_properties) +
                        (j + num_ghosts_1_thermo_properties)*ghostcell_dim_0_thermo_properties +
                        (k + num_ghosts_2_thermo_properties)*ghostcell_dim_0_thermo_properties*
                            ghostcell_dim_1_thermo_properties;
                    
                    p[idx_pressure] = (gamma[idx_thermo_properties] - double(1))*rho[idx_density]*
                        epsilon[idx_internal_energy];
                }
            }
        }
    }
}


/*
 * Compute the sound speed.
 */
void
EquationOfStateIdealGas::computeSoundSpeed(
    double* const c,
    const double* const rho,
    const double* const p,
    const double& gamma,
    const hier::IntVector& num_ghosts_sound_speed,
    const hier::IntVector& num_ghosts_density,
    const hier::IntVector& num_ghosts_pressure,
    const hier::IntVector& ghostcell_dims_sound_speed,
    const hier::IntVector& ghostcell_dims_density,
    const hier::IntVector& ghostcell_dims_pressure,
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
        
        const int num_ghosts_0_sound_speed = num_ghosts_sound_speed[0];
        const int num_ghosts_0_density = num_ghosts_density[0];
        const int num_ghosts_0_pressure = num_ghosts_pressure[0];
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_sound_speed = i + num_ghosts_0_sound_speed;
            const int idx_density = i + num_ghosts_0_density;
            const int idx_pressure = i + num_ghosts_0_pressure;
            
            c[idx_sound_speed] = sqrt(gamma*p[idx_pressure]/rho[idx_density]);
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
        
        const int num_ghosts_0_sound_speed = num_ghosts_sound_speed[0];
        const int num_ghosts_1_sound_speed = num_ghosts_sound_speed[1];
        const int ghostcell_dim_0_sound_speed = ghostcell_dims_sound_speed[0];
        
        const int num_ghosts_0_density = num_ghosts_density[0];
        const int num_ghosts_1_density = num_ghosts_density[1];
        const int ghostcell_dim_0_density = ghostcell_dims_density[0];
        
        const int num_ghosts_0_pressure = num_ghosts_pressure[0];
        const int num_ghosts_1_pressure = num_ghosts_pressure[1];
        const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_sound_speed = (i + num_ghosts_0_sound_speed) +
                    (j + num_ghosts_1_sound_speed)*ghostcell_dim_0_sound_speed;
                
                const int idx_density = (i + num_ghosts_0_density) +
                    (j + num_ghosts_1_density)*ghostcell_dim_0_density;
                
                const int idx_pressure = (i + num_ghosts_0_pressure) +
                    (j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure;
                
                c[idx_sound_speed] = sqrt(gamma*p[idx_pressure]/rho[idx_density]);
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
        
        const int num_ghosts_0_sound_speed = num_ghosts_sound_speed[0];
        const int num_ghosts_1_sound_speed = num_ghosts_sound_speed[1];
        const int num_ghosts_2_sound_speed = num_ghosts_sound_speed[2];
        const int ghostcell_dim_0_sound_speed = ghostcell_dims_sound_speed[0];
        const int ghostcell_dim_1_sound_speed = ghostcell_dims_sound_speed[1];
        
        const int num_ghosts_0_density = num_ghosts_density[0];
        const int num_ghosts_1_density = num_ghosts_density[1];
        const int num_ghosts_2_density = num_ghosts_density[2];
        const int ghostcell_dim_0_density = ghostcell_dims_density[0];
        const int ghostcell_dim_1_density = ghostcell_dims_density[1];
        
        const int num_ghosts_0_pressure = num_ghosts_pressure[0];
        const int num_ghosts_1_pressure = num_ghosts_pressure[1];
        const int num_ghosts_2_pressure = num_ghosts_pressure[2];
        const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
        const int ghostcell_dim_1_pressure = ghostcell_dims_pressure[1];
        
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
                    const int idx_sound_speed = (i + num_ghosts_0_sound_speed) +
                        (j + num_ghosts_1_sound_speed)*ghostcell_dim_0_sound_speed +
                        (k + num_ghosts_2_sound_speed)*ghostcell_dim_0_sound_speed*
                            ghostcell_dim_1_sound_speed;
                    
                    const int idx_density = (i + num_ghosts_0_density) +
                        (j + num_ghosts_1_density)*ghostcell_dim_0_density +
                        (k + num_ghosts_2_density)*ghostcell_dim_0_density*
                            ghostcell_dim_1_density;
                    
                    const int idx_pressure = (i + num_ghosts_0_pressure) +
                        (j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure +
                        (k + num_ghosts_2_pressure)*ghostcell_dim_0_pressure*
                            ghostcell_dim_1_pressure;
                    
                    c[idx_sound_speed] = sqrt(gamma*p[idx_pressure]/rho[idx_density]);
                }
            }
        }
    }
}


/*
 * Compute the sound speed.
 */
void
EquationOfStateIdealGas::computeSoundSpeed(
    double* const c,
    const double* const rho,
    const double* const p,
    const double* const gamma,
    const hier::IntVector& num_ghosts_sound_speed,
    const hier::IntVector& num_ghosts_density,
    const hier::IntVector& num_ghosts_pressure,
    const hier::IntVector& num_ghosts_thermo_properties,
    const hier::IntVector& ghostcell_dims_sound_speed,
    const hier::IntVector& ghostcell_dims_density,
    const hier::IntVector& ghostcell_dims_pressure,
    const hier::IntVector& ghostcell_dims_thermo_properties,
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
        
        const int num_ghosts_0_sound_speed = num_ghosts_sound_speed[0];
        const int num_ghosts_0_density = num_ghosts_density[0];
        const int num_ghosts_0_pressure = num_ghosts_pressure[0];
        const int num_ghosts_0_thermo_properties = num_ghosts_thermo_properties[0];
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_sound_speed = i + num_ghosts_0_sound_speed;
            const int idx_density = i + num_ghosts_0_density;
            const int idx_pressure = i + num_ghosts_0_pressure;
            const int idx_thermo_properties = i + num_ghosts_0_thermo_properties;
            
            c[idx_sound_speed] = sqrt(gamma[idx_thermo_properties]*p[idx_pressure]/
                rho[idx_density]);
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
        
        const int num_ghosts_0_sound_speed = num_ghosts_sound_speed[0];
        const int num_ghosts_1_sound_speed = num_ghosts_sound_speed[1];
        const int ghostcell_dim_0_sound_speed = ghostcell_dims_sound_speed[0];
        
        const int num_ghosts_0_density = num_ghosts_density[0];
        const int num_ghosts_1_density = num_ghosts_density[1];
        const int ghostcell_dim_0_density = ghostcell_dims_density[0];
        
        const int num_ghosts_0_pressure = num_ghosts_pressure[0];
        const int num_ghosts_1_pressure = num_ghosts_pressure[1];
        const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
        
        const int num_ghosts_0_thermo_properties = num_ghosts_thermo_properties[0];
        const int num_ghosts_1_thermo_properties = num_ghosts_thermo_properties[1];
        const int ghostcell_dim_0_thermo_properties = ghostcell_dims_thermo_properties[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_sound_speed = (i + num_ghosts_0_sound_speed) +
                    (j + num_ghosts_1_sound_speed)*ghostcell_dim_0_sound_speed;
                
                const int idx_density = (i + num_ghosts_0_density) +
                    (j + num_ghosts_1_density)*ghostcell_dim_0_density;
                
                const int idx_pressure = (i + num_ghosts_0_pressure) +
                    (j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure;
                
                const int idx_thermo_properties = (i + num_ghosts_0_thermo_properties) +
                    (j + num_ghosts_1_thermo_properties)*ghostcell_dim_0_thermo_properties;
                
                c[idx_sound_speed] = sqrt(gamma[idx_thermo_properties]*p[idx_pressure]/
                    rho[idx_density]);
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
        
        const int num_ghosts_0_sound_speed = num_ghosts_sound_speed[0];
        const int num_ghosts_1_sound_speed = num_ghosts_sound_speed[1];
        const int num_ghosts_2_sound_speed = num_ghosts_sound_speed[2];
        const int ghostcell_dim_0_sound_speed = ghostcell_dims_sound_speed[0];
        const int ghostcell_dim_1_sound_speed = ghostcell_dims_sound_speed[1];
        
        const int num_ghosts_0_density = num_ghosts_density[0];
        const int num_ghosts_1_density = num_ghosts_density[1];
        const int num_ghosts_2_density = num_ghosts_density[2];
        const int ghostcell_dim_0_density = ghostcell_dims_density[0];
        const int ghostcell_dim_1_density = ghostcell_dims_density[1];
        
        const int num_ghosts_0_pressure = num_ghosts_pressure[0];
        const int num_ghosts_1_pressure = num_ghosts_pressure[1];
        const int num_ghosts_2_pressure = num_ghosts_pressure[2];
        const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
        const int ghostcell_dim_1_pressure = ghostcell_dims_pressure[1];
        
        const int num_ghosts_0_thermo_properties = num_ghosts_thermo_properties[0];
        const int num_ghosts_1_thermo_properties = num_ghosts_thermo_properties[1];
        const int num_ghosts_2_thermo_properties = num_ghosts_thermo_properties[2];
        const int ghostcell_dim_0_thermo_properties = ghostcell_dims_thermo_properties[0];
        const int ghostcell_dim_1_thermo_properties = ghostcell_dims_thermo_properties[1];
        
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
                    const int idx_sound_speed = (i + num_ghosts_0_sound_speed) +
                        (j + num_ghosts_1_sound_speed)*ghostcell_dim_0_sound_speed +
                        (k + num_ghosts_2_sound_speed)*ghostcell_dim_0_sound_speed*
                            ghostcell_dim_1_sound_speed;
                    
                    const int idx_density = (i + num_ghosts_0_density) +
                        (j + num_ghosts_1_density)*ghostcell_dim_0_density +
                        (k + num_ghosts_2_density)*ghostcell_dim_0_density*
                            ghostcell_dim_1_density;
                    
                    const int idx_pressure = (i + num_ghosts_0_pressure) +
                        (j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure +
                        (k + num_ghosts_2_pressure)*ghostcell_dim_0_pressure*
                            ghostcell_dim_1_pressure;
                    
                    const int idx_thermo_properties = (i + num_ghosts_0_thermo_properties) +
                        (j + num_ghosts_1_thermo_properties)*ghostcell_dim_0_thermo_properties +
                        (k + num_ghosts_2_thermo_properties)*ghostcell_dim_0_thermo_properties*
                            ghostcell_dim_1_thermo_properties;
                    
                    c[idx_sound_speed] = sqrt(gamma[idx_thermo_properties]*p[idx_pressure]/
                        rho[idx_density]);
                }
            }
        }
    }
}


/*
 * Compute the specific internal energy.
 */
void
EquationOfStateIdealGas::computeInternalEnergy(
    double* const epsilon,
    const double* const rho,
    const double* const p,
    const double& gamma,
    const hier::IntVector& num_ghosts_internal_energy,
    const hier::IntVector& num_ghosts_density,
    const hier::IntVector& num_ghosts_pressure,
    const hier::IntVector& ghostcell_dims_internal_energy,
    const hier::IntVector& ghostcell_dims_density,
    const hier::IntVector& ghostcell_dims_pressure,
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
        
        const int num_ghosts_0_internal_energy = num_ghosts_internal_energy[0];
        const int num_ghosts_0_density = num_ghosts_density[0];
        const int num_ghosts_0_pressure = num_ghosts_pressure[0];
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_internal_energy = i + num_ghosts_0_internal_energy;
            const int idx_density = i + num_ghosts_0_density;
            const int idx_pressure = i + num_ghosts_0_pressure;
            
            epsilon[idx_internal_energy] = p[idx_pressure]/((gamma - double(1))*rho[idx_density]);
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
        
        const int num_ghosts_0_internal_energy = num_ghosts_internal_energy[0];
        const int num_ghosts_1_internal_energy = num_ghosts_internal_energy[1];
        const int ghostcell_dim_0_internal_energy = ghostcell_dims_internal_energy[0];
        
        const int num_ghosts_0_density = num_ghosts_density[0];
        const int num_ghosts_1_density = num_ghosts_density[1];
        const int ghostcell_dim_0_density = ghostcell_dims_density[0];
        
        const int num_ghosts_0_pressure = num_ghosts_pressure[0];
        const int num_ghosts_1_pressure = num_ghosts_pressure[1];
        const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_internal_energy = (i + num_ghosts_0_internal_energy) +
                    (j + num_ghosts_1_internal_energy)*ghostcell_dim_0_internal_energy;
                
                const int idx_density = (i + num_ghosts_0_density) +
                    (j + num_ghosts_1_density)*ghostcell_dim_0_density;
                
                const int idx_pressure = (i + num_ghosts_0_pressure) +
                    (j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure;
                
                epsilon[idx_internal_energy] = p[idx_pressure]/((gamma - double(1))*rho[idx_density]);
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
        
        const int num_ghosts_0_internal_energy = num_ghosts_internal_energy[0];
        const int num_ghosts_1_internal_energy = num_ghosts_internal_energy[1];
        const int num_ghosts_2_internal_energy = num_ghosts_internal_energy[2];
        const int ghostcell_dim_0_internal_energy = ghostcell_dims_internal_energy[0];
        const int ghostcell_dim_1_internal_energy = ghostcell_dims_internal_energy[1];
        
        const int num_ghosts_0_density = num_ghosts_density[0];
        const int num_ghosts_1_density = num_ghosts_density[1];
        const int num_ghosts_2_density = num_ghosts_density[2];
        const int ghostcell_dim_0_density = ghostcell_dims_density[0];
        const int ghostcell_dim_1_density = ghostcell_dims_density[1];
        
        const int num_ghosts_0_pressure = num_ghosts_pressure[0];
        const int num_ghosts_1_pressure = num_ghosts_pressure[1];
        const int num_ghosts_2_pressure = num_ghosts_pressure[2];
        const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
        const int ghostcell_dim_1_pressure = ghostcell_dims_pressure[1];
        
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
                    const int idx_internal_energy = (i + num_ghosts_0_internal_energy) +
                        (j + num_ghosts_1_internal_energy)*ghostcell_dim_0_internal_energy +
                        (k + num_ghosts_2_internal_energy)*ghostcell_dim_0_internal_energy*
                            ghostcell_dim_1_internal_energy;
                    
                    const int idx_density = (i + num_ghosts_0_density) +
                        (j + num_ghosts_1_density)*ghostcell_dim_0_density +
                        (k + num_ghosts_2_density)*ghostcell_dim_0_density*
                            ghostcell_dim_1_density;
                    
                    const int idx_pressure = (i + num_ghosts_0_pressure) +
                        (j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure +
                        (k + num_ghosts_2_pressure)*ghostcell_dim_0_pressure*
                            ghostcell_dim_1_pressure;
                    
                    epsilon[idx_internal_energy] = p[idx_pressure]/((gamma - double(1))*rho[idx_density]);
                }
            }
        }
    }
}


/*
 * Compute the specific internal energy.
 */
void
EquationOfStateIdealGas::computeInternalEnergy(
    double* const epsilon,
    const double* const rho,
    const double* const p,
    const double* const gamma,
    const hier::IntVector& num_ghosts_internal_energy,
    const hier::IntVector& num_ghosts_density,
    const hier::IntVector& num_ghosts_pressure,
    const hier::IntVector& num_ghosts_thermo_properties,
    const hier::IntVector& ghostcell_dims_internal_energy,
    const hier::IntVector& ghostcell_dims_density,
    const hier::IntVector& ghostcell_dims_pressure,
    const hier::IntVector& ghostcell_dims_thermo_properties,
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
        
        const int num_ghosts_0_internal_energy = num_ghosts_internal_energy[0];
        const int num_ghosts_0_density = num_ghosts_density[0];
        const int num_ghosts_0_pressure = num_ghosts_pressure[0];
        const int num_ghosts_0_thermo_properties = num_ghosts_thermo_properties[0];
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_internal_energy = i + num_ghosts_0_internal_energy;
            const int idx_density = i + num_ghosts_0_density;
            const int idx_pressure = i + num_ghosts_0_pressure;
            const int idx_thermo_properties = i + num_ghosts_0_thermo_properties;
            
            epsilon[idx_internal_energy] = p[idx_pressure]/((gamma[idx_thermo_properties] - double(1))*
                rho[idx_density]);
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
        
        const int num_ghosts_0_internal_energy = num_ghosts_internal_energy[0];
        const int num_ghosts_1_internal_energy = num_ghosts_internal_energy[1];
        const int ghostcell_dim_0_internal_energy = ghostcell_dims_internal_energy[0];
        
        const int num_ghosts_0_density = num_ghosts_density[0];
        const int num_ghosts_1_density = num_ghosts_density[1];
        const int ghostcell_dim_0_density = ghostcell_dims_density[0];
        
        const int num_ghosts_0_pressure = num_ghosts_pressure[0];
        const int num_ghosts_1_pressure = num_ghosts_pressure[1];
        const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
        
        const int num_ghosts_0_thermo_properties = num_ghosts_thermo_properties[0];
        const int num_ghosts_1_thermo_properties = num_ghosts_thermo_properties[1];
        const int ghostcell_dim_0_thermo_properties = ghostcell_dims_thermo_properties[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_internal_energy = (i + num_ghosts_0_internal_energy) +
                    (j + num_ghosts_1_internal_energy)*ghostcell_dim_0_internal_energy;
                
                const int idx_density = (i + num_ghosts_0_density) +
                    (j + num_ghosts_1_density)*ghostcell_dim_0_density;
                
                const int idx_pressure = (i + num_ghosts_0_pressure) +
                    (j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure;
                
                const int idx_thermo_properties = (i + num_ghosts_0_thermo_properties) +
                    (j + num_ghosts_1_thermo_properties)*ghostcell_dim_0_thermo_properties;
                
                epsilon[idx_internal_energy] = p[idx_pressure]/((gamma[idx_thermo_properties] - double(1))*
                    rho[idx_density]);
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
        
        const int num_ghosts_0_internal_energy = num_ghosts_internal_energy[0];
        const int num_ghosts_1_internal_energy = num_ghosts_internal_energy[1];
        const int num_ghosts_2_internal_energy = num_ghosts_internal_energy[2];
        const int ghostcell_dim_0_internal_energy = ghostcell_dims_internal_energy[0];
        const int ghostcell_dim_1_internal_energy = ghostcell_dims_internal_energy[1];
        
        const int num_ghosts_0_density = num_ghosts_density[0];
        const int num_ghosts_1_density = num_ghosts_density[1];
        const int num_ghosts_2_density = num_ghosts_density[2];
        const int ghostcell_dim_0_density = ghostcell_dims_density[0];
        const int ghostcell_dim_1_density = ghostcell_dims_density[1];
        
        const int num_ghosts_0_pressure = num_ghosts_pressure[0];
        const int num_ghosts_1_pressure = num_ghosts_pressure[1];
        const int num_ghosts_2_pressure = num_ghosts_pressure[2];
        const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
        const int ghostcell_dim_1_pressure = ghostcell_dims_pressure[1];
        
        const int num_ghosts_0_thermo_properties = num_ghosts_thermo_properties[0];
        const int num_ghosts_1_thermo_properties = num_ghosts_thermo_properties[1];
        const int num_ghosts_2_thermo_properties = num_ghosts_thermo_properties[2];
        const int ghostcell_dim_0_thermo_properties = ghostcell_dims_thermo_properties[0];
        const int ghostcell_dim_1_thermo_properties = ghostcell_dims_thermo_properties[1];
        
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
                    const int idx_internal_energy = (i + num_ghosts_0_internal_energy) +
                        (j + num_ghosts_1_internal_energy)*ghostcell_dim_0_internal_energy +
                        (k + num_ghosts_2_internal_energy)*ghostcell_dim_0_internal_energy*
                            ghostcell_dim_1_internal_energy;
                    
                    const int idx_density = (i + num_ghosts_0_density) +
                        (j + num_ghosts_1_density)*ghostcell_dim_0_density +
                        (k + num_ghosts_2_density)*ghostcell_dim_0_density*
                            ghostcell_dim_1_density;
                    
                    const int idx_pressure = (i + num_ghosts_0_pressure) +
                        (j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure +
                        (k + num_ghosts_2_pressure)*ghostcell_dim_0_pressure*
                            ghostcell_dim_1_pressure;
                    
                    const int idx_thermo_properties = (i + num_ghosts_0_thermo_properties) +
                        (j + num_ghosts_1_thermo_properties)*ghostcell_dim_0_thermo_properties +
                        (k + num_ghosts_2_thermo_properties)*ghostcell_dim_0_thermo_properties*
                            ghostcell_dim_1_thermo_properties;
                    
                    epsilon[idx_internal_energy] = p[idx_pressure]/((gamma[idx_thermo_properties] - double(1))*
                        rho[idx_density]);
                }
            }
        }
    }
}


/*
 * Compute the specific enthalpy.
 */
void
EquationOfStateIdealGas::computeEnthalpy(
    double* const h,
    const double* const rho,
    const double* const p,
    const double& gamma,
    const hier::IntVector& num_ghosts_enthalpy,
    const hier::IntVector& num_ghosts_density,
    const hier::IntVector& num_ghosts_pressure,
    const hier::IntVector& ghostcell_dims_enthalpy,
    const hier::IntVector& ghostcell_dims_density,
    const hier::IntVector& ghostcell_dims_pressure,
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
        
        const int num_ghosts_0_enthalpy = num_ghosts_enthalpy[0];
        const int num_ghosts_0_density = num_ghosts_density[0];
        const int num_ghosts_0_pressure = num_ghosts_pressure[0];
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_enthalpy = i + num_ghosts_0_enthalpy;
            const int idx_density = i + num_ghosts_0_density;
            const int idx_pressure = i + num_ghosts_0_pressure;
            
            h[idx_enthalpy] = gamma*p[idx_pressure]/((gamma - double(1))*rho[idx_density]);
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
        
        const int num_ghosts_0_enthalpy = num_ghosts_enthalpy[0];
        const int num_ghosts_1_enthalpy = num_ghosts_enthalpy[1];
        const int ghostcell_dim_0_enthalpy = ghostcell_dims_enthalpy[0];
        
        const int num_ghosts_0_density = num_ghosts_density[0];
        const int num_ghosts_1_density = num_ghosts_density[1];
        const int ghostcell_dim_0_density = ghostcell_dims_density[0];
        
        const int num_ghosts_0_pressure = num_ghosts_pressure[0];
        const int num_ghosts_1_pressure = num_ghosts_pressure[1];
        const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_enthalpy = (i + num_ghosts_0_enthalpy) +
                    (j + num_ghosts_1_enthalpy)*ghostcell_dim_0_enthalpy;
                
                const int idx_density = (i + num_ghosts_0_density) +
                    (j + num_ghosts_1_density)*ghostcell_dim_0_density;
                
                const int idx_pressure = (i + num_ghosts_0_pressure) +
                    (j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure;
                
                h[idx_enthalpy] = gamma*p[idx_pressure]/((gamma - double(1))*rho[idx_density]);
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
        
        const int num_ghosts_0_enthalpy = num_ghosts_enthalpy[0];
        const int num_ghosts_1_enthalpy = num_ghosts_enthalpy[1];
        const int num_ghosts_2_enthalpy = num_ghosts_enthalpy[2];
        const int ghostcell_dim_0_enthalpy = ghostcell_dims_enthalpy[0];
        const int ghostcell_dim_1_enthalpy = ghostcell_dims_enthalpy[1];
        
        const int num_ghosts_0_density = num_ghosts_density[0];
        const int num_ghosts_1_density = num_ghosts_density[1];
        const int num_ghosts_2_density = num_ghosts_density[2];
        const int ghostcell_dim_0_density = ghostcell_dims_density[0];
        const int ghostcell_dim_1_density = ghostcell_dims_density[1];
        
        const int num_ghosts_0_pressure = num_ghosts_pressure[0];
        const int num_ghosts_1_pressure = num_ghosts_pressure[1];
        const int num_ghosts_2_pressure = num_ghosts_pressure[2];
        const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
        const int ghostcell_dim_1_pressure = ghostcell_dims_pressure[1];
        
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
                    const int idx_enthalpy = (i + num_ghosts_0_enthalpy) +
                        (j + num_ghosts_1_enthalpy)*ghostcell_dim_0_enthalpy +
                        (k + num_ghosts_2_enthalpy)*ghostcell_dim_0_enthalpy*
                            ghostcell_dim_1_enthalpy;
                    
                    const int idx_density = (i + num_ghosts_0_density) +
                        (j + num_ghosts_1_density)*ghostcell_dim_0_density +
                        (k + num_ghosts_2_density)*ghostcell_dim_0_density*
                            ghostcell_dim_1_density;
                    
                    const int idx_pressure = (i + num_ghosts_0_pressure) +
                        (j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure +
                        (k + num_ghosts_2_pressure)*ghostcell_dim_0_pressure*
                            ghostcell_dim_1_pressure;
                    
                    h[idx_enthalpy] = gamma*p[idx_pressure]/((gamma - double(1))*rho[idx_density]);
                }
            }
        }
    }
}


/*
 * Compute the specific enthalpy.
 */
void
EquationOfStateIdealGas::computeEnthalpy(
    double* const h,
    const double* const rho,
    const double* const p,
    const double* const gamma,
    const hier::IntVector& num_ghosts_enthalpy,
    const hier::IntVector& num_ghosts_density,
    const hier::IntVector& num_ghosts_pressure,
    const hier::IntVector& num_ghosts_thermo_properties,
    const hier::IntVector& ghostcell_dims_enthalpy,
    const hier::IntVector& ghostcell_dims_density,
    const hier::IntVector& ghostcell_dims_pressure,
    const hier::IntVector& ghostcell_dims_thermo_properties,
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
        
        const int num_ghosts_0_enthalpy = num_ghosts_enthalpy[0];
        const int num_ghosts_0_density = num_ghosts_density[0];
        const int num_ghosts_0_pressure = num_ghosts_pressure[0];
        const int num_ghosts_0_thermo_properties = num_ghosts_thermo_properties[0];
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_enthalpy = i + num_ghosts_0_enthalpy;
            const int idx_density = i + num_ghosts_0_density;
            const int idx_pressure = i + num_ghosts_0_pressure;
            const int idx_thermo_properties = i + num_ghosts_0_thermo_properties;
            
            h[idx_enthalpy] = gamma[idx_thermo_properties]*p[idx_pressure]/
                ((gamma[idx_thermo_properties] - double(1))*rho[idx_density]);
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
        
        const int num_ghosts_0_enthalpy = num_ghosts_enthalpy[0];
        const int num_ghosts_1_enthalpy = num_ghosts_enthalpy[1];
        const int ghostcell_dim_0_enthalpy = ghostcell_dims_enthalpy[0];
        
        const int num_ghosts_0_density = num_ghosts_density[0];
        const int num_ghosts_1_density = num_ghosts_density[1];
        const int ghostcell_dim_0_density = ghostcell_dims_density[0];
        
        const int num_ghosts_0_pressure = num_ghosts_pressure[0];
        const int num_ghosts_1_pressure = num_ghosts_pressure[1];
        const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
        
        const int num_ghosts_0_thermo_properties = num_ghosts_thermo_properties[0];
        const int num_ghosts_1_thermo_properties = num_ghosts_thermo_properties[1];
        const int ghostcell_dim_0_thermo_properties = ghostcell_dims_thermo_properties[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_enthalpy = (i + num_ghosts_0_enthalpy) +
                    (j + num_ghosts_1_enthalpy)*ghostcell_dim_0_enthalpy;
                
                const int idx_density = (i + num_ghosts_0_density) +
                    (j + num_ghosts_1_density)*ghostcell_dim_0_density;
                
                const int idx_pressure = (i + num_ghosts_0_pressure) +
                    (j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure;
                
                const int idx_thermo_properties = (i + num_ghosts_0_thermo_properties) +
                    (j + num_ghosts_1_thermo_properties)*ghostcell_dim_0_thermo_properties;
                
                h[idx_enthalpy] = gamma[idx_thermo_properties]*p[idx_pressure]/
                    ((gamma[idx_thermo_properties] - double(1))*rho[idx_density]);
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
        
        const int num_ghosts_0_enthalpy = num_ghosts_enthalpy[0];
        const int num_ghosts_1_enthalpy = num_ghosts_enthalpy[1];
        const int num_ghosts_2_enthalpy = num_ghosts_enthalpy[2];
        const int ghostcell_dim_0_enthalpy = ghostcell_dims_enthalpy[0];
        const int ghostcell_dim_1_enthalpy = ghostcell_dims_enthalpy[1];
        
        const int num_ghosts_0_density = num_ghosts_density[0];
        const int num_ghosts_1_density = num_ghosts_density[1];
        const int num_ghosts_2_density = num_ghosts_density[2];
        const int ghostcell_dim_0_density = ghostcell_dims_density[0];
        const int ghostcell_dim_1_density = ghostcell_dims_density[1];
        
        const int num_ghosts_0_pressure = num_ghosts_pressure[0];
        const int num_ghosts_1_pressure = num_ghosts_pressure[1];
        const int num_ghosts_2_pressure = num_ghosts_pressure[2];
        const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
        const int ghostcell_dim_1_pressure = ghostcell_dims_pressure[1];
        
        const int num_ghosts_0_thermo_properties = num_ghosts_thermo_properties[0];
        const int num_ghosts_1_thermo_properties = num_ghosts_thermo_properties[1];
        const int num_ghosts_2_thermo_properties = num_ghosts_thermo_properties[2];
        const int ghostcell_dim_0_thermo_properties = ghostcell_dims_thermo_properties[0];
        const int ghostcell_dim_1_thermo_properties = ghostcell_dims_thermo_properties[1];
        
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
                    const int idx_enthalpy = (i + num_ghosts_0_enthalpy) +
                        (j + num_ghosts_1_enthalpy)*ghostcell_dim_0_enthalpy +
                        (k + num_ghosts_2_enthalpy)*ghostcell_dim_0_enthalpy*
                            ghostcell_dim_1_enthalpy;
                    
                    const int idx_density = (i + num_ghosts_0_density) +
                        (j + num_ghosts_1_density)*ghostcell_dim_0_density +
                        (k + num_ghosts_2_density)*ghostcell_dim_0_density*
                            ghostcell_dim_1_density;
                    
                    const int idx_pressure = (i + num_ghosts_0_pressure) +
                        (j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure +
                        (k + num_ghosts_2_pressure)*ghostcell_dim_0_pressure*
                            ghostcell_dim_1_pressure;
                    
                    const int idx_thermo_properties = (i + num_ghosts_0_thermo_properties) +
                        (j + num_ghosts_1_thermo_properties)*ghostcell_dim_0_thermo_properties +
                        (k + num_ghosts_2_thermo_properties)*ghostcell_dim_0_thermo_properties*
                            ghostcell_dim_1_thermo_properties;
                    
                    h[idx_enthalpy] = gamma[idx_thermo_properties]*p[idx_pressure]/
                        ((gamma[idx_thermo_properties] - double(1))*rho[idx_density]);
                }
            }
        }
    }
}


/*
 * Compute the temperature.
 */
void
EquationOfStateIdealGas::computeTemperature(
    double* const T,
    const double* const rho,
    const double* const p,
    const double& gamma,
    const double& c_v,
    const hier::IntVector& num_ghosts_temperature,
    const hier::IntVector& num_ghosts_density,
    const hier::IntVector& num_ghosts_pressure,
    const hier::IntVector& ghostcell_dims_temperature,
    const hier::IntVector& ghostcell_dims_density,
    const hier::IntVector& ghostcell_dims_pressure,
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
        
        const int num_ghosts_0_temperature = num_ghosts_temperature[0];
        const int num_ghosts_0_density = num_ghosts_density[0];
        const int num_ghosts_0_pressure = num_ghosts_pressure[0];
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_temperature = i + num_ghosts_0_temperature;
            const int idx_density = i + num_ghosts_0_density;
            const int idx_pressure = i + num_ghosts_0_pressure;
            
            T[idx_temperature] = p[idx_pressure]/((gamma - double(1))*c_v*rho[idx_density]);
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
        
        const int num_ghosts_0_temperature = num_ghosts_temperature[0];
        const int num_ghosts_1_temperature = num_ghosts_temperature[1];
        const int ghostcell_dim_0_temperature = ghostcell_dims_temperature[0];
        
        const int num_ghosts_0_density = num_ghosts_density[0];
        const int num_ghosts_1_density = num_ghosts_density[1];
        const int ghostcell_dim_0_density = ghostcell_dims_density[0];
        
        const int num_ghosts_0_pressure = num_ghosts_pressure[0];
        const int num_ghosts_1_pressure = num_ghosts_pressure[1];
        const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_temperature = (i + num_ghosts_0_temperature) +
                    (j + num_ghosts_1_temperature)*ghostcell_dim_0_temperature;
                
                const int idx_density = (i + num_ghosts_0_density) +
                    (j + num_ghosts_1_density)*ghostcell_dim_0_density;
                
                const int idx_pressure = (i + num_ghosts_0_pressure) +
                    (j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure;
                
                T[idx_temperature] = p[idx_pressure]/((gamma - double(1))*c_v*rho[idx_density]);
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
        
        const int num_ghosts_0_temperature = num_ghosts_temperature[0];
        const int num_ghosts_1_temperature = num_ghosts_temperature[1];
        const int num_ghosts_2_temperature = num_ghosts_temperature[2];
        const int ghostcell_dim_0_temperature = ghostcell_dims_temperature[0];
        const int ghostcell_dim_1_temperature = ghostcell_dims_temperature[1];
        
        const int num_ghosts_0_density = num_ghosts_density[0];
        const int num_ghosts_1_density = num_ghosts_density[1];
        const int num_ghosts_2_density = num_ghosts_density[2];
        const int ghostcell_dim_0_density = ghostcell_dims_density[0];
        const int ghostcell_dim_1_density = ghostcell_dims_density[1];
        
        const int num_ghosts_0_pressure = num_ghosts_pressure[0];
        const int num_ghosts_1_pressure = num_ghosts_pressure[1];
        const int num_ghosts_2_pressure = num_ghosts_pressure[2];
        const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
        const int ghostcell_dim_1_pressure = ghostcell_dims_pressure[1];
        
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
                    const int idx_temperature = (i + num_ghosts_0_temperature) +
                        (j + num_ghosts_1_temperature)*ghostcell_dim_0_temperature +
                        (k + num_ghosts_2_temperature)*ghostcell_dim_0_temperature*
                            ghostcell_dim_1_temperature;
                    
                    const int idx_density = (i + num_ghosts_0_density) +
                        (j + num_ghosts_1_density)*ghostcell_dim_0_density +
                        (k + num_ghosts_2_density)*ghostcell_dim_0_density*
                            ghostcell_dim_1_density;
                    
                    const int idx_pressure = (i + num_ghosts_0_pressure) +
                        (j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure +
                        (k + num_ghosts_2_pressure)*ghostcell_dim_0_pressure*
                            ghostcell_dim_1_pressure;
                    
                    T[idx_temperature] = p[idx_pressure]/((gamma - double(1))*c_v*rho[idx_density]);
                }
            }
        }
    }
}


/*
 * Compute the temperature.
 */
void
EquationOfStateIdealGas::computeTemperature(
    double* const T,
    const double* const rho,
    const double* const p,
    const double* const gamma,
    const double* const c_v,
    const hier::IntVector& num_ghosts_temperature,
    const hier::IntVector& num_ghosts_density,
    const hier::IntVector& num_ghosts_pressure,
    const hier::IntVector& num_ghosts_thermo_properties,
    const hier::IntVector& ghostcell_dims_temperature,
    const hier::IntVector& ghostcell_dims_density,
    const hier::IntVector& ghostcell_dims_pressure,
    const hier::IntVector& ghostcell_dims_thermo_properties,
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
        
        const int num_ghosts_0_temperature = num_ghosts_temperature[0];
        const int num_ghosts_0_density = num_ghosts_density[0];
        const int num_ghosts_0_pressure = num_ghosts_pressure[0];
        const int num_ghosts_0_thermo_properties = num_ghosts_thermo_properties[0];
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_temperature = i + num_ghosts_0_temperature;
            const int idx_density = i + num_ghosts_0_density;
            const int idx_pressure = i + num_ghosts_0_pressure;
            const int idx_thermo_properties = i + num_ghosts_0_thermo_properties;
            
            T[idx_temperature] = p[idx_pressure]/((gamma[idx_thermo_properties] - double(1))*
                c_v[idx_thermo_properties]*rho[idx_density]);
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
        
        const int num_ghosts_0_temperature = num_ghosts_temperature[0];
        const int num_ghosts_1_temperature = num_ghosts_temperature[1];
        const int ghostcell_dim_0_temperature = ghostcell_dims_temperature[0];
        
        const int num_ghosts_0_density = num_ghosts_density[0];
        const int num_ghosts_1_density = num_ghosts_density[1];
        const int ghostcell_dim_0_density = ghostcell_dims_density[0];
        
        const int num_ghosts_0_pressure = num_ghosts_pressure[0];
        const int num_ghosts_1_pressure = num_ghosts_pressure[1];
        const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
        
        const int num_ghosts_0_thermo_properties = num_ghosts_thermo_properties[0];
        const int num_ghosts_1_thermo_properties = num_ghosts_thermo_properties[1];
        const int ghostcell_dim_0_thermo_properties = ghostcell_dims_thermo_properties[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_temperature = (i + num_ghosts_0_temperature) +
                    (j + num_ghosts_1_temperature)*ghostcell_dim_0_temperature;
                
                const int idx_density = (i + num_ghosts_0_density) +
                    (j + num_ghosts_1_density)*ghostcell_dim_0_density;
                
                const int idx_pressure = (i + num_ghosts_0_pressure) +
                    (j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure;
                
                const int idx_thermo_properties = (i + num_ghosts_0_thermo_properties) +
                    (j + num_ghosts_1_thermo_properties)*ghostcell_dim_0_thermo_properties;
                
                T[idx_temperature] = p[idx_pressure]/((gamma[idx_thermo_properties] - double(1))*
                    c_v[idx_thermo_properties]*rho[idx_density]);
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
        
        const int num_ghosts_0_temperature = num_ghosts_temperature[0];
        const int num_ghosts_1_temperature = num_ghosts_temperature[1];
        const int num_ghosts_2_temperature = num_ghosts_temperature[2];
        const int ghostcell_dim_0_temperature = ghostcell_dims_temperature[0];
        const int ghostcell_dim_1_temperature = ghostcell_dims_temperature[1];
        
        const int num_ghosts_0_density = num_ghosts_density[0];
        const int num_ghosts_1_density = num_ghosts_density[1];
        const int num_ghosts_2_density = num_ghosts_density[2];
        const int ghostcell_dim_0_density = ghostcell_dims_density[0];
        const int ghostcell_dim_1_density = ghostcell_dims_density[1];
        
        const int num_ghosts_0_pressure = num_ghosts_pressure[0];
        const int num_ghosts_1_pressure = num_ghosts_pressure[1];
        const int num_ghosts_2_pressure = num_ghosts_pressure[2];
        const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
        const int ghostcell_dim_1_pressure = ghostcell_dims_pressure[1];
        
        const int num_ghosts_0_thermo_properties = num_ghosts_thermo_properties[0];
        const int num_ghosts_1_thermo_properties = num_ghosts_thermo_properties[1];
        const int num_ghosts_2_thermo_properties = num_ghosts_thermo_properties[2];
        const int ghostcell_dim_0_thermo_properties = ghostcell_dims_thermo_properties[0];
        const int ghostcell_dim_1_thermo_properties = ghostcell_dims_thermo_properties[1];
        
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
                    const int idx_temperature = (i + num_ghosts_0_temperature) +
                        (j + num_ghosts_1_temperature)*ghostcell_dim_0_temperature +
                        (k + num_ghosts_2_temperature)*ghostcell_dim_0_temperature*
                            ghostcell_dim_1_temperature;
                    
                    const int idx_density = (i + num_ghosts_0_density) +
                        (j + num_ghosts_1_density)*ghostcell_dim_0_density +
                        (k + num_ghosts_2_density)*ghostcell_dim_0_density*
                            ghostcell_dim_1_density;
                    
                    const int idx_pressure = (i + num_ghosts_0_pressure) +
                        (j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure +
                        (k + num_ghosts_2_pressure)*ghostcell_dim_0_pressure*
                            ghostcell_dim_1_pressure;
                    
                    const int idx_thermo_properties = (i + num_ghosts_0_thermo_properties) +
                        (j + num_ghosts_1_thermo_properties)*ghostcell_dim_0_thermo_properties +
                        (k + num_ghosts_2_thermo_properties)*ghostcell_dim_0_thermo_properties*
                            ghostcell_dim_1_thermo_properties;
                    
                    T[idx_temperature] = p[idx_pressure]/((gamma[idx_thermo_properties] - double(1))*
                        c_v[idx_thermo_properties]*rho[idx_density]);
                }
            }
        }
    }
}


/*
 * Compute the specific internal energy from temperature.
 */
void
EquationOfStateIdealGas::computeInternalEnergyFromTemperature(
    double* const epsilon,
    const double* const T,
    const double& c_v,
    const hier::IntVector& num_ghosts_internal_energy,
    const hier::IntVector& num_ghosts_temperature,
    const hier::IntVector& ghostcell_dims_internal_energy,
    const hier::IntVector& ghostcell_dims_temperature,
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
        
        const int num_ghosts_0_internal_energy = num_ghosts_internal_energy[0];
        const int num_ghosts_0_temperature = num_ghosts_temperature[0];
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_internal_energy = i + num_ghosts_0_internal_energy;
            const int idx_temperature = i + num_ghosts_0_temperature;
            
            epsilon[idx_internal_energy] = c_v*T[idx_temperature];
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
        
        const int num_ghosts_0_internal_energy = num_ghosts_internal_energy[0];
        const int num_ghosts_1_internal_energy = num_ghosts_internal_energy[1];
        const int ghostcell_dim_0_internal_energy = ghostcell_dims_internal_energy[0];
        
        const int num_ghosts_0_temperature = num_ghosts_temperature[0];
        const int num_ghosts_1_temperature = num_ghosts_temperature[1];
        const int ghostcell_dim_0_temperature = ghostcell_dims_temperature[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_internal_energy = (i + num_ghosts_0_internal_energy) +
                    (j + num_ghosts_1_internal_energy)*ghostcell_dim_0_internal_energy;
                
                const int idx_temperature = (i + num_ghosts_0_temperature) +
                    (j + num_ghosts_1_temperature)*ghostcell_dim_0_temperature;
                
                epsilon[idx_internal_energy] = c_v*T[idx_temperature];
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
        
        const int num_ghosts_0_internal_energy = num_ghosts_internal_energy[0];
        const int num_ghosts_1_internal_energy = num_ghosts_internal_energy[1];
        const int num_ghosts_2_internal_energy = num_ghosts_internal_energy[2];
        const int ghostcell_dim_0_internal_energy = ghostcell_dims_internal_energy[0];
        const int ghostcell_dim_1_internal_energy = ghostcell_dims_internal_energy[1];
        
        const int num_ghosts_0_temperature = num_ghosts_temperature[0];
        const int num_ghosts_1_temperature = num_ghosts_temperature[1];
        const int num_ghosts_2_temperature = num_ghosts_temperature[2];
        const int ghostcell_dim_0_temperature = ghostcell_dims_temperature[0];
        const int ghostcell_dim_1_temperature = ghostcell_dims_temperature[1];
        
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
                    const int idx_internal_energy = (i + num_ghosts_0_internal_energy) +
                        (j + num_ghosts_1_internal_energy)*ghostcell_dim_0_internal_energy +
                        (k + num_ghosts_2_internal_energy)*ghostcell_dim_0_internal_energy*
                            ghostcell_dim_1_internal_energy;
                    
                    const int idx_temperature = (i + num_ghosts_0_temperature) +
                        (j + num_ghosts_1_temperature)*ghostcell_dim_0_temperature +
                        (k + num_ghosts_2_temperature)*ghostcell_dim_0_temperature*
                            ghostcell_dim_1_temperature;
                    
                    epsilon[idx_internal_energy] = c_v*T[idx_temperature];
                }
            }
        }
    }
}


/*
 * Compute the specific internal energy from temperature.
 */
void
EquationOfStateIdealGas::computeInternalEnergyFromTemperature(
    double* const epsilon,
    const double* const T,
    const double* const c_v,
    const hier::IntVector& num_ghosts_internal_energy,
    const hier::IntVector& num_ghosts_temperature,
    const hier::IntVector& num_ghosts_thermo_properties,
    const hier::IntVector& ghostcell_dims_internal_energy,
    const hier::IntVector& ghostcell_dims_temperature,
    const hier::IntVector& ghostcell_dims_thermo_properties,
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
        
        const int num_ghosts_0_internal_energy = num_ghosts_internal_energy[0];
        const int num_ghosts_0_temperature = num_ghosts_temperature[0];
        const int num_ghosts_0_thermo_properties = num_ghosts_thermo_properties[0];
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_internal_energy = i + num_ghosts_0_internal_energy;
            const int idx_temperature = i + num_ghosts_0_temperature;
            const int idx_thermo_properties = i + num_ghosts_0_thermo_properties;
            
            epsilon[idx_internal_energy] = c_v[idx_thermo_properties]*T[idx_temperature];
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
        
        const int num_ghosts_0_internal_energy = num_ghosts_internal_energy[0];
        const int num_ghosts_1_internal_energy = num_ghosts_internal_energy[1];
        const int ghostcell_dim_0_internal_energy = ghostcell_dims_internal_energy[0];
        
        const int num_ghosts_0_temperature = num_ghosts_temperature[0];
        const int num_ghosts_1_temperature = num_ghosts_temperature[1];
        const int ghostcell_dim_0_temperature = ghostcell_dims_temperature[0];
        
        const int num_ghosts_0_thermo_properties = num_ghosts_thermo_properties[0];
        const int num_ghosts_1_thermo_properties = num_ghosts_thermo_properties[1];
        const int ghostcell_dim_0_thermo_properties = ghostcell_dims_thermo_properties[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_internal_energy = (i + num_ghosts_0_internal_energy) +
                    (j + num_ghosts_1_internal_energy)*ghostcell_dim_0_internal_energy;
                
                const int idx_temperature = (i + num_ghosts_0_temperature) +
                    (j + num_ghosts_1_temperature)*ghostcell_dim_0_temperature;
                
                const int idx_thermo_properties = (i + num_ghosts_0_thermo_properties) +
                    (j + num_ghosts_1_thermo_properties)*ghostcell_dim_0_thermo_properties;
                
                epsilon[idx_internal_energy] = c_v[idx_thermo_properties]*T[idx_temperature];
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
        
        const int num_ghosts_0_internal_energy = num_ghosts_internal_energy[0];
        const int num_ghosts_1_internal_energy = num_ghosts_internal_energy[1];
        const int num_ghosts_2_internal_energy = num_ghosts_internal_energy[2];
        const int ghostcell_dim_0_internal_energy = ghostcell_dims_internal_energy[0];
        const int ghostcell_dim_1_internal_energy = ghostcell_dims_internal_energy[1];
        
        const int num_ghosts_0_temperature = num_ghosts_temperature[0];
        const int num_ghosts_1_temperature = num_ghosts_temperature[1];
        const int num_ghosts_2_temperature = num_ghosts_temperature[2];
        const int ghostcell_dim_0_temperature = ghostcell_dims_temperature[0];
        const int ghostcell_dim_1_temperature = ghostcell_dims_temperature[1];
        
        const int num_ghosts_0_thermo_properties = num_ghosts_thermo_properties[0];
        const int num_ghosts_1_thermo_properties = num_ghosts_thermo_properties[1];
        const int num_ghosts_2_thermo_properties = num_ghosts_thermo_properties[2];
        const int ghostcell_dim_0_thermo_properties = ghostcell_dims_thermo_properties[0];
        const int ghostcell_dim_1_thermo_properties = ghostcell_dims_thermo_properties[1];
        
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
                    const int idx_internal_energy = (i + num_ghosts_0_internal_energy) +
                        (j + num_ghosts_1_internal_energy)*ghostcell_dim_0_internal_energy +
                        (k + num_ghosts_2_internal_energy)*ghostcell_dim_0_internal_energy*
                            ghostcell_dim_1_internal_energy;
                    
                    const int idx_temperature = (i + num_ghosts_0_temperature) +
                        (j + num_ghosts_1_temperature)*ghostcell_dim_0_temperature +
                        (k + num_ghosts_2_temperature)*ghostcell_dim_0_temperature*
                            ghostcell_dim_1_temperature;
                    
                    const int idx_thermo_properties = (i + num_ghosts_0_thermo_properties) +
                        (j + num_ghosts_1_thermo_properties)*ghostcell_dim_0_thermo_properties +
                        (k + num_ghosts_2_thermo_properties)*ghostcell_dim_0_thermo_properties*
                            ghostcell_dim_1_thermo_properties;
                    
                    epsilon[idx_internal_energy] = c_v[idx_thermo_properties]*T[idx_temperature];
                }
            }
        }
    }
}


/*
 * Compute the isochoric specific heat capacity.
 */
void
EquationOfStateIdealGas::computeIsochoricSpecificHeatCapacity(
    double* const c_v,
    const double& c_v_src,
    const hier::IntVector& num_ghosts_isochoric_specific_heat_capacity,
    const hier::IntVector& ghostcell_dims_isochoric_specific_heat_capacity,
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
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_isochoric_specific_heat_capacity =
                i + num_ghosts_0_isochoric_specific_heat_capacity;
            
            c_v[idx_isochoric_specific_heat_capacity] = c_v_src;
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
        
        const int num_ghosts_0_isochoric_specific_heat_capacity =
            num_ghosts_isochoric_specific_heat_capacity[0];
        const int num_ghosts_1_isochoric_specific_heat_capacity =
            num_ghosts_isochoric_specific_heat_capacity[1];
        const int ghostcell_dim_0_isochoric_specific_heat_capacity =
            ghostcell_dims_isochoric_specific_heat_capacity[0];
        
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
                
                c_v[idx_isochoric_specific_heat_capacity] = c_v_src;
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
        
        const int num_ghosts_0_isochoric_specific_heat_capacity =
            num_ghosts_isochoric_specific_heat_capacity[0];
        const int num_ghosts_1_isochoric_specific_heat_capacity =
            num_ghosts_isochoric_specific_heat_capacity[1];
        const int num_ghosts_2_isochoric_specific_heat_capacity =
            num_ghosts_isochoric_specific_heat_capacity[2];
        const int ghostcell_dim_0_isochoric_specific_heat_capacity =
            ghostcell_dims_isochoric_specific_heat_capacity[0];
        const int ghostcell_dim_1_isochoric_specific_heat_capacity =
            ghostcell_dims_isochoric_specific_heat_capacity[1];
        
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
                    
                    c_v[idx_isochoric_specific_heat_capacity] = c_v_src;
                }
            }
        }
    }
}


/*
 * Compute the isochoric specific heat capacity.
 */
void
EquationOfStateIdealGas::computeIsochoricSpecificHeatCapacity(
    double* const c_v,
    const double* const c_v_src,
    const hier::IntVector& num_ghosts_isochoric_specific_heat_capacity,
    const hier::IntVector& num_ghosts_thermo_properties,
    const hier::IntVector& ghostcell_dims_isochoric_specific_heat_capacity,
    const hier::IntVector& ghostcell_dims_thermo_properties,
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
        
        const int num_ghosts_0_thermo_properties = num_ghosts_thermo_properties[0];
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_isochoric_specific_heat_capacity =
                i + num_ghosts_0_isochoric_specific_heat_capacity;
            
            const int idx_thermo_properties = i + num_ghosts_0_thermo_properties;
            
            c_v[idx_isochoric_specific_heat_capacity] = c_v_src[idx_thermo_properties];
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
        
        const int num_ghosts_0_isochoric_specific_heat_capacity =
            num_ghosts_isochoric_specific_heat_capacity[0];
        const int num_ghosts_1_isochoric_specific_heat_capacity =
            num_ghosts_isochoric_specific_heat_capacity[1];
        const int ghostcell_dim_0_isochoric_specific_heat_capacity =
            ghostcell_dims_isochoric_specific_heat_capacity[0];
        
        const int num_ghosts_0_thermo_properties = num_ghosts_thermo_properties[0];
        const int num_ghosts_1_thermo_properties = num_ghosts_thermo_properties[1];
        const int ghostcell_dim_0_thermo_properties = ghostcell_dims_thermo_properties[0];
        
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
                
                const int idx_thermo_properties = (i + num_ghosts_0_thermo_properties) +
                    (j + num_ghosts_1_thermo_properties)*ghostcell_dim_0_thermo_properties;
                
                c_v[idx_isochoric_specific_heat_capacity] = c_v_src[idx_thermo_properties];
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
        
        const int num_ghosts_0_isochoric_specific_heat_capacity =
            num_ghosts_isochoric_specific_heat_capacity[0];
        const int num_ghosts_1_isochoric_specific_heat_capacity =
            num_ghosts_isochoric_specific_heat_capacity[1];
        const int num_ghosts_2_isochoric_specific_heat_capacity =
            num_ghosts_isochoric_specific_heat_capacity[2];
        const int ghostcell_dim_0_isochoric_specific_heat_capacity =
            ghostcell_dims_isochoric_specific_heat_capacity[0];
        const int ghostcell_dim_1_isochoric_specific_heat_capacity =
            ghostcell_dims_isochoric_specific_heat_capacity[1];
        
        const int num_ghosts_0_thermo_properties = num_ghosts_thermo_properties[0];
        const int num_ghosts_1_thermo_properties = num_ghosts_thermo_properties[1];
        const int num_ghosts_2_thermo_properties = num_ghosts_thermo_properties[2];
        const int ghostcell_dim_0_thermo_properties = ghostcell_dims_thermo_properties[0];
        const int ghostcell_dim_1_thermo_properties = ghostcell_dims_thermo_properties[1];
        
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
                    
                    const int idx_thermo_properties = (i + num_ghosts_0_thermo_properties) +
                        (j + num_ghosts_1_thermo_properties)*ghostcell_dim_0_thermo_properties +
                        (k + num_ghosts_2_thermo_properties)*ghostcell_dim_0_thermo_properties*
                            ghostcell_dim_1_thermo_properties;
                    
                    c_v[idx_isochoric_specific_heat_capacity] = c_v_src[idx_thermo_properties];
                }
            }
        }
    }
}


/*
 * Compute the isobaric specific heat capacity.
 */
void
EquationOfStateIdealGas::computeIsobaricSpecificHeatCapacity(
    double* const c_p,
    const double& c_p_src,
    const hier::IntVector& num_ghosts_isobaric_specific_heat_capacity,
    const hier::IntVector& ghostcell_dims_isobaric_specific_heat_capacity,
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
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_isobaric_specific_heat_capacity =
                i + num_ghosts_0_isobaric_specific_heat_capacity;
            
            c_p[idx_isobaric_specific_heat_capacity] = c_p_src;
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
        
        const int num_ghosts_0_isobaric_specific_heat_capacity =
            num_ghosts_isobaric_specific_heat_capacity[0];
        const int num_ghosts_1_isobaric_specific_heat_capacity =
            num_ghosts_isobaric_specific_heat_capacity[1];
        const int ghostcell_dim_0_isobaric_specific_heat_capacity =
            ghostcell_dims_isobaric_specific_heat_capacity[0];
        
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
                
                c_p[idx_isobaric_specific_heat_capacity] = c_p_src;
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
        
        const int num_ghosts_0_isobaric_specific_heat_capacity =
            num_ghosts_isobaric_specific_heat_capacity[0];
        const int num_ghosts_1_isobaric_specific_heat_capacity =
            num_ghosts_isobaric_specific_heat_capacity[1];
        const int num_ghosts_2_isobaric_specific_heat_capacity =
            num_ghosts_isobaric_specific_heat_capacity[2];
        const int ghostcell_dim_0_isobaric_specific_heat_capacity =
            ghostcell_dims_isobaric_specific_heat_capacity[0];
        const int ghostcell_dim_1_isobaric_specific_heat_capacity =
            ghostcell_dims_isobaric_specific_heat_capacity[1];
        
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
                    
                    c_p[idx_isobaric_specific_heat_capacity] = c_p_src;
                }
            }
        }
    }
}


/*
 * Compute the isobaric specific heat capacity.
 */
void
EquationOfStateIdealGas::computeIsobaricSpecificHeatCapacity(
    double* const c_p,
    const double* const c_p_src,
    const hier::IntVector& num_ghosts_isobaric_specific_heat_capacity,
    const hier::IntVector& num_ghosts_thermo_properties,
    const hier::IntVector& ghostcell_dims_isobaric_specific_heat_capacity,
    const hier::IntVector& ghostcell_dims_thermo_properties,
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
        
        const int num_ghosts_0_thermo_properties = num_ghosts_thermo_properties[0];
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_isobaric_specific_heat_capacity =
                i + num_ghosts_0_isobaric_specific_heat_capacity;
            
            const int idx_thermo_properties = i + num_ghosts_0_thermo_properties;
            
            c_p[idx_isobaric_specific_heat_capacity] = c_p_src[idx_thermo_properties];
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
        
        const int num_ghosts_0_isobaric_specific_heat_capacity =
            num_ghosts_isobaric_specific_heat_capacity[0];
        const int num_ghosts_1_isobaric_specific_heat_capacity =
            num_ghosts_isobaric_specific_heat_capacity[1];
        const int ghostcell_dim_0_isobaric_specific_heat_capacity =
            ghostcell_dims_isobaric_specific_heat_capacity[0];
        
        const int num_ghosts_0_thermo_properties = num_ghosts_thermo_properties[0];
        const int num_ghosts_1_thermo_properties = num_ghosts_thermo_properties[1];
        const int ghostcell_dim_0_thermo_properties = ghostcell_dims_thermo_properties[0];
        
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
                
                const int idx_thermo_properties = (i + num_ghosts_0_thermo_properties) +
                    (j + num_ghosts_1_thermo_properties)*ghostcell_dim_0_thermo_properties;
                
                c_p[idx_isobaric_specific_heat_capacity] = c_p_src[idx_thermo_properties];
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
        
        const int num_ghosts_0_isobaric_specific_heat_capacity =
            num_ghosts_isobaric_specific_heat_capacity[0];
        const int num_ghosts_1_isobaric_specific_heat_capacity =
            num_ghosts_isobaric_specific_heat_capacity[1];
        const int num_ghosts_2_isobaric_specific_heat_capacity =
            num_ghosts_isobaric_specific_heat_capacity[2];
        const int ghostcell_dim_0_isobaric_specific_heat_capacity =
            ghostcell_dims_isobaric_specific_heat_capacity[0];
        const int ghostcell_dim_1_isobaric_specific_heat_capacity =
            ghostcell_dims_isobaric_specific_heat_capacity[1];
        
        const int num_ghosts_0_thermo_properties = num_ghosts_thermo_properties[0];
        const int num_ghosts_1_thermo_properties = num_ghosts_thermo_properties[1];
        const int num_ghosts_2_thermo_properties = num_ghosts_thermo_properties[2];
        const int ghostcell_dim_0_thermo_properties = ghostcell_dims_thermo_properties[0];
        const int ghostcell_dim_1_thermo_properties = ghostcell_dims_thermo_properties[1];
        
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
                    
                    const int idx_thermo_properties = (i + num_ghosts_0_thermo_properties) +
                        (j + num_ghosts_1_thermo_properties)*ghostcell_dim_0_thermo_properties +
                        (k + num_ghosts_2_thermo_properties)*ghostcell_dim_0_thermo_properties*
                            ghostcell_dim_1_thermo_properties;
                    
                    c_p[idx_isobaric_specific_heat_capacity] = c_p_src[idx_thermo_properties];
                }
            }
        }
    }
}


/*
 * Compute the partial derivative of pressure w.r.t. density under constant specific internal energy.
 */
void
EquationOfStateIdealGas::computePartialPressurePartialDensity(
    double* const Psi,
    const double* const rho,
    const double* const p,
    const hier::IntVector& num_ghosts_partial_pressure_partial_density,
    const hier::IntVector& num_ghosts_density,
    const hier::IntVector& num_ghosts_pressure,
    const hier::IntVector& ghostcell_dims_partial_pressure_partial_density,
    const hier::IntVector& ghostcell_dims_density,
    const hier::IntVector& ghostcell_dims_pressure,
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
        
        const int num_ghosts_0_partial_pressure_partial_density = num_ghosts_partial_pressure_partial_density[0];
        const int num_ghosts_0_density = num_ghosts_density[0];
        const int num_ghosts_0_pressure = num_ghosts_pressure[0];
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_partial_pressure_partial_density = i + num_ghosts_0_partial_pressure_partial_density;
            const int idx_density = i + num_ghosts_0_density;
            const int idx_pressure = i + num_ghosts_0_pressure;
            
            Psi[idx_partial_pressure_partial_density] = p[idx_pressure]/rho[idx_density];
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
        
        const int num_ghosts_0_partial_pressure_partial_density = num_ghosts_partial_pressure_partial_density[0];
        const int num_ghosts_1_partial_pressure_partial_density = num_ghosts_partial_pressure_partial_density[1];
        const int ghostcell_dim_0_partial_pressure_partial_density = ghostcell_dims_partial_pressure_partial_density[0];
        
        const int num_ghosts_0_density = num_ghosts_density[0];
        const int num_ghosts_1_density = num_ghosts_density[1];
        const int ghostcell_dim_0_density = ghostcell_dims_density[0];
        
        const int num_ghosts_0_pressure = num_ghosts_pressure[0];
        const int num_ghosts_1_pressure = num_ghosts_pressure[1];
        const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_partial_pressure_partial_density = (i + num_ghosts_0_partial_pressure_partial_density) +
                    (j + num_ghosts_1_partial_pressure_partial_density)*ghostcell_dim_0_partial_pressure_partial_density;
                
                const int idx_density = (i + num_ghosts_0_density) +
                    (j + num_ghosts_1_density)*ghostcell_dim_0_density;
                
                const int idx_pressure = (i + num_ghosts_0_pressure) +
                    (j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure;
                
                Psi[idx_partial_pressure_partial_density] = p[idx_pressure]/rho[idx_density];
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
        
        const int num_ghosts_0_partial_pressure_partial_density = num_ghosts_partial_pressure_partial_density[0];
        const int num_ghosts_1_partial_pressure_partial_density = num_ghosts_partial_pressure_partial_density[1];
        const int num_ghosts_2_partial_pressure_partial_density = num_ghosts_partial_pressure_partial_density[2];
        const int ghostcell_dim_0_partial_pressure_partial_density = ghostcell_dims_partial_pressure_partial_density[0];
        const int ghostcell_dim_1_partial_pressure_partial_density = ghostcell_dims_partial_pressure_partial_density[1];
        
        const int num_ghosts_0_density = num_ghosts_density[0];
        const int num_ghosts_1_density = num_ghosts_density[1];
        const int num_ghosts_2_density = num_ghosts_density[2];
        const int ghostcell_dim_0_density = ghostcell_dims_density[0];
        const int ghostcell_dim_1_density = ghostcell_dims_density[1];
        
        const int num_ghosts_0_pressure = num_ghosts_pressure[0];
        const int num_ghosts_1_pressure = num_ghosts_pressure[1];
        const int num_ghosts_2_pressure = num_ghosts_pressure[2];
        const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
        const int ghostcell_dim_1_pressure = ghostcell_dims_pressure[1];
        
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
                    const int idx_partial_pressure_partial_density = (i + num_ghosts_0_partial_pressure_partial_density) +
                        (j + num_ghosts_1_partial_pressure_partial_density)*ghostcell_dim_0_partial_pressure_partial_density +
                        (k + num_ghosts_2_partial_pressure_partial_density)*ghostcell_dim_0_partial_pressure_partial_density*
                            ghostcell_dim_1_partial_pressure_partial_density;
                    
                    const int idx_density = (i + num_ghosts_0_density) +
                        (j + num_ghosts_1_density)*ghostcell_dim_0_density +
                        (k + num_ghosts_2_density)*ghostcell_dim_0_density*
                            ghostcell_dim_1_density;
                    
                    const int idx_pressure = (i + num_ghosts_0_pressure) +
                        (j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure +
                        (k + num_ghosts_2_pressure)*ghostcell_dim_0_pressure*
                            ghostcell_dim_1_pressure;
                    
                    Psi[idx_partial_pressure_partial_density] = p[idx_pressure]/rho[idx_density];
                }
            }
        }
    }
}


/*
 * Compute the Gruneisen parameter (partial derivative of pressure w.r.t. specific internal energy under
 * constant density divided by density).
 */
void
EquationOfStateIdealGas::computeGruneisenParameter(
    double* const Gamma,
    const double& gamma,
    const hier::IntVector& num_ghosts_gruneisen_parameter,
    const hier::IntVector& ghostcell_dims_gruneisen_parameter,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims) const
{
    const double Gamma_src = gamma - double(1);
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and numbers of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int num_ghosts_0_gruneisen_parameter =
            num_ghosts_gruneisen_parameter[0];
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_gruneisen_parameter =
                i + num_ghosts_0_gruneisen_parameter;
            
            Gamma[idx_gruneisen_parameter] = Gamma_src;
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
        
        const int num_ghosts_0_gruneisen_parameter =
            num_ghosts_gruneisen_parameter[0];
        const int num_ghosts_1_gruneisen_parameter =
            num_ghosts_gruneisen_parameter[1];
        const int ghostcell_dim_0_gruneisen_parameter =
            ghostcell_dims_gruneisen_parameter[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_gruneisen_parameter =
                    (i + num_ghosts_0_gruneisen_parameter) +
                    (j + num_ghosts_1_gruneisen_parameter)*
                        ghostcell_dim_0_gruneisen_parameter;
                
                Gamma[idx_gruneisen_parameter] = Gamma_src;
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
        
        const int num_ghosts_0_gruneisen_parameter =
            num_ghosts_gruneisen_parameter[0];
        const int num_ghosts_1_gruneisen_parameter =
            num_ghosts_gruneisen_parameter[1];
        const int num_ghosts_2_gruneisen_parameter =
            num_ghosts_gruneisen_parameter[2];
        const int ghostcell_dim_0_gruneisen_parameter =
            ghostcell_dims_gruneisen_parameter[0];
        const int ghostcell_dim_1_gruneisen_parameter =
            ghostcell_dims_gruneisen_parameter[1];
        
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
                    const int idx_gruneisen_parameter =
                        (i + num_ghosts_0_gruneisen_parameter) +
                        (j + num_ghosts_1_gruneisen_parameter)*
                            ghostcell_dim_0_gruneisen_parameter +
                        (k + num_ghosts_2_gruneisen_parameter)*
                            ghostcell_dim_0_gruneisen_parameter*
                            ghostcell_dim_1_gruneisen_parameter;
                    
                    Gamma[idx_gruneisen_parameter] = Gamma_src;
                }
            }
        }
    }
}


/*
 * Compute the Gruneisen parameter (partial derivative of pressure w.r.t. specific internal energy under
 * constant density divided by density).
 */
void
EquationOfStateIdealGas::computeGruneisenParameter(
    double* const Gamma,
    const double* const gamma,
    const hier::IntVector& num_ghosts_gruneisen_parameter,
    const hier::IntVector& num_ghosts_thermo_properties,
    const hier::IntVector& ghostcell_dims_gruneisen_parameter,
    const hier::IntVector& ghostcell_dims_thermo_properties,
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
        
        const int num_ghosts_0_gruneisen_parameter =
            num_ghosts_gruneisen_parameter[0];
        
        const int num_ghosts_0_thermo_properties = num_ghosts_thermo_properties[0];
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_gruneisen_parameter =
                i + num_ghosts_0_gruneisen_parameter;
            
            const int idx_thermo_properties = i + num_ghosts_0_thermo_properties;
            
            Gamma[idx_gruneisen_parameter] = gamma[idx_thermo_properties] - double(1);
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
        
        const int num_ghosts_0_gruneisen_parameter =
            num_ghosts_gruneisen_parameter[0];
        const int num_ghosts_1_gruneisen_parameter =
            num_ghosts_gruneisen_parameter[1];
        const int ghostcell_dim_0_gruneisen_parameter =
            ghostcell_dims_gruneisen_parameter[0];
        
        const int num_ghosts_0_thermo_properties = num_ghosts_thermo_properties[0];
        const int num_ghosts_1_thermo_properties = num_ghosts_thermo_properties[1];
        const int ghostcell_dim_0_thermo_properties = ghostcell_dims_thermo_properties[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_gruneisen_parameter =
                    (i + num_ghosts_0_gruneisen_parameter) +
                    (j + num_ghosts_1_gruneisen_parameter)*
                        ghostcell_dim_0_gruneisen_parameter;
                
                const int idx_thermo_properties = (i + num_ghosts_0_thermo_properties) +
                    (j + num_ghosts_1_thermo_properties)*ghostcell_dim_0_thermo_properties;
                
                Gamma[idx_gruneisen_parameter] = gamma[idx_thermo_properties] - double(1);
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
        
        const int num_ghosts_0_gruneisen_parameter =
            num_ghosts_gruneisen_parameter[0];
        const int num_ghosts_1_gruneisen_parameter =
            num_ghosts_gruneisen_parameter[1];
        const int num_ghosts_2_gruneisen_parameter =
            num_ghosts_gruneisen_parameter[2];
        const int ghostcell_dim_0_gruneisen_parameter =
            ghostcell_dims_gruneisen_parameter[0];
        const int ghostcell_dim_1_gruneisen_parameter =
            ghostcell_dims_gruneisen_parameter[1];
        
        const int num_ghosts_0_thermo_properties = num_ghosts_thermo_properties[0];
        const int num_ghosts_1_thermo_properties = num_ghosts_thermo_properties[1];
        const int num_ghosts_2_thermo_properties = num_ghosts_thermo_properties[2];
        const int ghostcell_dim_0_thermo_properties = ghostcell_dims_thermo_properties[0];
        const int ghostcell_dim_1_thermo_properties = ghostcell_dims_thermo_properties[1];
        
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
                    const int idx_gruneisen_parameter =
                        (i + num_ghosts_0_gruneisen_parameter) +
                        (j + num_ghosts_1_gruneisen_parameter)*
                            ghostcell_dim_0_gruneisen_parameter +
                        (k + num_ghosts_2_gruneisen_parameter)*
                            ghostcell_dim_0_gruneisen_parameter*
                            ghostcell_dim_1_gruneisen_parameter;
                    
                    const int idx_thermo_properties = (i + num_ghosts_0_thermo_properties) +
                        (j + num_ghosts_1_thermo_properties)*ghostcell_dim_0_thermo_properties +
                        (k + num_ghosts_2_thermo_properties)*ghostcell_dim_0_thermo_properties*
                            ghostcell_dim_1_thermo_properties;
                    
                    Gamma[idx_gruneisen_parameter] = gamma[idx_thermo_properties] - double(1);
                }
            }
        }
    }
}


/*
 * Compute the partial derivative of internal energy w.r.t. pressure under constant density.
 */
void
EquationOfStateIdealGas::computeIsochoricPartialInternalEnergyPartialPressure(
    double* const xi,
    const double& gamma,
    const hier::IntVector& num_ghosts_partial_internal_energy_partial_pressure,
    const hier::IntVector& ghostcell_dims_partial_internal_energy_partial_pressure,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims) const
{
    const double xi_src = double(1)/(gamma - double(1));
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and numbers of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int num_ghosts_0_partial_internal_energy_partial_pressure =
            num_ghosts_partial_internal_energy_partial_pressure[0];
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_partial_internal_energy_partial_pressure =
                i + num_ghosts_0_partial_internal_energy_partial_pressure;
            
            xi[idx_partial_internal_energy_partial_pressure] = xi_src;
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
        
        const int num_ghosts_0_partial_internal_energy_partial_pressure =
            num_ghosts_partial_internal_energy_partial_pressure[0];
        const int num_ghosts_1_partial_internal_energy_partial_pressure =
            num_ghosts_partial_internal_energy_partial_pressure[1];
        const int ghostcell_dim_0_partial_internal_energy_partial_pressure =
            ghostcell_dims_partial_internal_energy_partial_pressure[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_partial_internal_energy_partial_pressure =
                    (i + num_ghosts_0_partial_internal_energy_partial_pressure) +
                    (j + num_ghosts_1_partial_internal_energy_partial_pressure)*
                        ghostcell_dim_0_partial_internal_energy_partial_pressure;
                
                xi[idx_partial_internal_energy_partial_pressure] = xi_src;
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
        
        const int num_ghosts_0_partial_internal_energy_partial_pressure =
            num_ghosts_partial_internal_energy_partial_pressure[0];
        const int num_ghosts_1_partial_internal_energy_partial_pressure =
            num_ghosts_partial_internal_energy_partial_pressure[1];
        const int num_ghosts_2_partial_internal_energy_partial_pressure =
            num_ghosts_partial_internal_energy_partial_pressure[2];
        const int ghostcell_dim_0_partial_internal_energy_partial_pressure =
            ghostcell_dims_partial_internal_energy_partial_pressure[0];
        const int ghostcell_dim_1_partial_internal_energy_partial_pressure =
            ghostcell_dims_partial_internal_energy_partial_pressure[1];
        
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
                    const int idx_partial_internal_energy_partial_pressure =
                        (i + num_ghosts_0_partial_internal_energy_partial_pressure) +
                        (j + num_ghosts_1_partial_internal_energy_partial_pressure)*
                            ghostcell_dim_0_partial_internal_energy_partial_pressure +
                        (k + num_ghosts_2_partial_internal_energy_partial_pressure)*
                            ghostcell_dim_0_partial_internal_energy_partial_pressure*
                            ghostcell_dim_1_partial_internal_energy_partial_pressure;
                    
                    xi[idx_partial_internal_energy_partial_pressure] = xi_src;
                }
            }
        }
    }
}


/*
 * Compute the partial derivative of internal energy w.r.t. pressure under constant density.
 */
void
EquationOfStateIdealGas::computeIsochoricPartialInternalEnergyPartialPressure(
    double* const xi,
    const double* const gamma,
    const hier::IntVector& num_ghosts_partial_internal_energy_partial_pressure,
    const hier::IntVector& num_ghosts_thermo_properties,
    const hier::IntVector& ghostcell_dims_partial_internal_energy_partial_pressure,
    const hier::IntVector& ghostcell_dims_thermo_properties,
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
        
        const int num_ghosts_0_partial_internal_energy_partial_pressure =
            num_ghosts_partial_internal_energy_partial_pressure[0];
        
        const int num_ghosts_0_thermo_properties = num_ghosts_thermo_properties[0];
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_partial_internal_energy_partial_pressure =
                i + num_ghosts_0_partial_internal_energy_partial_pressure;
            
            const int idx_thermo_properties = i + num_ghosts_0_thermo_properties;
            
            xi[idx_partial_internal_energy_partial_pressure] =
                double(1)/(gamma[idx_thermo_properties] - double(1));
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
        
        const int num_ghosts_0_partial_internal_energy_partial_pressure =
            num_ghosts_partial_internal_energy_partial_pressure[0];
        const int num_ghosts_1_partial_internal_energy_partial_pressure =
            num_ghosts_partial_internal_energy_partial_pressure[1];
        const int ghostcell_dim_0_partial_internal_energy_partial_pressure =
            ghostcell_dims_partial_internal_energy_partial_pressure[0];
        
        const int num_ghosts_0_thermo_properties = num_ghosts_thermo_properties[0];
        const int num_ghosts_1_thermo_properties = num_ghosts_thermo_properties[1];
        const int ghostcell_dim_0_thermo_properties = ghostcell_dims_thermo_properties[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_partial_internal_energy_partial_pressure =
                    (i + num_ghosts_0_partial_internal_energy_partial_pressure) +
                    (j + num_ghosts_1_partial_internal_energy_partial_pressure)*
                        ghostcell_dim_0_partial_internal_energy_partial_pressure;
                
                const int idx_thermo_properties = (i + num_ghosts_0_thermo_properties) +
                    (j + num_ghosts_1_thermo_properties)*ghostcell_dim_0_thermo_properties;
                
                xi[idx_partial_internal_energy_partial_pressure] =
                    double(1)/(gamma[idx_thermo_properties] - double(1));
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
        
        const int num_ghosts_0_partial_internal_energy_partial_pressure =
            num_ghosts_partial_internal_energy_partial_pressure[0];
        const int num_ghosts_1_partial_internal_energy_partial_pressure =
            num_ghosts_partial_internal_energy_partial_pressure[1];
        const int num_ghosts_2_partial_internal_energy_partial_pressure =
            num_ghosts_partial_internal_energy_partial_pressure[2];
        const int ghostcell_dim_0_partial_internal_energy_partial_pressure =
            ghostcell_dims_partial_internal_energy_partial_pressure[0];
        const int ghostcell_dim_1_partial_internal_energy_partial_pressure =
            ghostcell_dims_partial_internal_energy_partial_pressure[1];
        
        const int num_ghosts_0_thermo_properties = num_ghosts_thermo_properties[0];
        const int num_ghosts_1_thermo_properties = num_ghosts_thermo_properties[1];
        const int num_ghosts_2_thermo_properties = num_ghosts_thermo_properties[2];
        const int ghostcell_dim_0_thermo_properties = ghostcell_dims_thermo_properties[0];
        const int ghostcell_dim_1_thermo_properties = ghostcell_dims_thermo_properties[1];
        
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
                    const int idx_partial_internal_energy_partial_pressure =
                        (i + num_ghosts_0_partial_internal_energy_partial_pressure) +
                        (j + num_ghosts_1_partial_internal_energy_partial_pressure)*
                            ghostcell_dim_0_partial_internal_energy_partial_pressure +
                        (k + num_ghosts_2_partial_internal_energy_partial_pressure)*
                            ghostcell_dim_0_partial_internal_energy_partial_pressure*
                            ghostcell_dim_1_partial_internal_energy_partial_pressure;
                    
                    const int idx_thermo_properties = (i + num_ghosts_0_thermo_properties) +
                        (j + num_ghosts_1_thermo_properties)*ghostcell_dim_0_thermo_properties +
                        (k + num_ghosts_2_thermo_properties)*ghostcell_dim_0_thermo_properties*
                            ghostcell_dim_1_thermo_properties;
                    
                    xi[idx_partial_internal_energy_partial_pressure] =
                        double(1)/(gamma[idx_thermo_properties] - double(1));
                }
            }
        }
    }
}


/*
 * Compute the partial derivative of internal energy w.r.t. density under constant pressure.
 */
void
EquationOfStateIdealGas::computeIsobaricPartialInternalEnergyPartialDensity(
    double* const delta,
    const hier::IntVector& num_ghosts_partial_internal_energy_partial_density,
    const hier::IntVector& ghostcell_dims_partial_internal_energy_partial_density,
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
        
        const int num_ghosts_0_partial_internal_energy_partial_density =
            num_ghosts_partial_internal_energy_partial_density[0];
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_partial_internal_energy_partial_density =
                i + num_ghosts_0_partial_internal_energy_partial_density;
            
            delta[idx_partial_internal_energy_partial_density] = double(0);
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
        
        const int num_ghosts_0_partial_internal_energy_partial_density =
            num_ghosts_partial_internal_energy_partial_density[0];
        const int num_ghosts_1_partial_internal_energy_partial_density =
            num_ghosts_partial_internal_energy_partial_density[1];
        const int ghostcell_dim_0_partial_internal_energy_partial_density =
            ghostcell_dims_partial_internal_energy_partial_density[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_partial_internal_energy_partial_density =
                    (i + num_ghosts_0_partial_internal_energy_partial_density) +
                    (j + num_ghosts_1_partial_internal_energy_partial_density)*
                        ghostcell_dim_0_partial_internal_energy_partial_density;
                
                delta[idx_partial_internal_energy_partial_density] = double(0);
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
        
        const int num_ghosts_0_partial_internal_energy_partial_density =
            num_ghosts_partial_internal_energy_partial_density[0];
        const int num_ghosts_1_partial_internal_energy_partial_density =
            num_ghosts_partial_internal_energy_partial_density[1];
        const int num_ghosts_2_partial_internal_energy_partial_density =
            num_ghosts_partial_internal_energy_partial_density[2];
        const int ghostcell_dim_0_partial_internal_energy_partial_density =
            ghostcell_dims_partial_internal_energy_partial_density[0];
        const int ghostcell_dim_1_partial_internal_energy_partial_density =
            ghostcell_dims_partial_internal_energy_partial_density[1];
        
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
                    const int idx_partial_internal_energy_partial_density =
                        (i + num_ghosts_0_partial_internal_energy_partial_density) +
                        (j + num_ghosts_1_partial_internal_energy_partial_density)*
                            ghostcell_dim_0_partial_internal_energy_partial_density +
                        (k + num_ghosts_2_partial_internal_energy_partial_density)*
                            ghostcell_dim_0_partial_internal_energy_partial_density*
                            ghostcell_dim_1_partial_internal_energy_partial_density;
                    
                    delta[idx_partial_internal_energy_partial_density] = double(0);
                }
            }
        }
    }
}


/*
 * Compute the density.
 */
void
EquationOfStateIdealGas::computeDensity(
    double* const rho,
    const double* const p,
    const double* const T,
    const double& R,
    const hier::IntVector& num_ghosts_density,
    const hier::IntVector& num_ghosts_pressure,
    const hier::IntVector& num_ghosts_temperature,
    const hier::IntVector& ghostcell_dims_density,
    const hier::IntVector& ghostcell_dims_pressure,
    const hier::IntVector& ghostcell_dims_temperature,
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
        
        const int num_ghosts_0_density = num_ghosts_density[0];
        const int num_ghosts_0_pressure = num_ghosts_pressure[0];
        const int num_ghosts_0_temperature = num_ghosts_temperature[0];
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_density = i + num_ghosts_0_density;
            const int idx_pressure = i + num_ghosts_0_pressure;
            const int idx_temperature = i + num_ghosts_0_temperature;
            
            rho[idx_density] = p[idx_pressure]/(R*T[idx_temperature]);
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
        
        const int num_ghosts_0_density = num_ghosts_density[0];
        const int num_ghosts_1_density = num_ghosts_density[1];
        const int ghostcell_dim_0_density = ghostcell_dims_density[0];
        
        const int num_ghosts_0_pressure = num_ghosts_pressure[0];
        const int num_ghosts_1_pressure = num_ghosts_pressure[1];
        const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
        
        const int num_ghosts_0_temperature = num_ghosts_temperature[0];
        const int num_ghosts_1_temperature = num_ghosts_temperature[1];
        const int ghostcell_dim_0_temperature = ghostcell_dims_temperature[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_density = (i + num_ghosts_0_density) +
                    (j + num_ghosts_1_density)*ghostcell_dim_0_density;
                
                const int idx_pressure = (i + num_ghosts_0_pressure) +
                    (j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure;
                
                const int idx_temperature = (i + num_ghosts_0_temperature) +
                    (j + num_ghosts_1_temperature)*ghostcell_dim_0_temperature;
                
                rho[idx_density] = p[idx_pressure]/(R*T[idx_temperature]);
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
        
        const int num_ghosts_0_density = num_ghosts_density[0];
        const int num_ghosts_1_density = num_ghosts_density[1];
        const int num_ghosts_2_density = num_ghosts_density[2];
        const int ghostcell_dim_0_density = ghostcell_dims_density[0];
        const int ghostcell_dim_1_density = ghostcell_dims_density[1];
        
        const int num_ghosts_0_pressure = num_ghosts_pressure[0];
        const int num_ghosts_1_pressure = num_ghosts_pressure[1];
        const int num_ghosts_2_pressure = num_ghosts_pressure[2];
        const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
        const int ghostcell_dim_1_pressure = ghostcell_dims_pressure[1];
        
        const int num_ghosts_0_temperature = num_ghosts_temperature[0];
        const int num_ghosts_1_temperature = num_ghosts_temperature[1];
        const int num_ghosts_2_temperature = num_ghosts_temperature[2];
        const int ghostcell_dim_0_temperature = ghostcell_dims_temperature[0];
        const int ghostcell_dim_1_temperature = ghostcell_dims_temperature[1];
        
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
                    const int idx_density = (i + num_ghosts_0_density) +
                        (j + num_ghosts_1_density)*ghostcell_dim_0_density +
                        (k + num_ghosts_2_density)*ghostcell_dim_0_density*
                            ghostcell_dim_1_density;
                    
                    const int idx_pressure = (i + num_ghosts_0_pressure) +
                        (j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure +
                        (k + num_ghosts_2_pressure)*ghostcell_dim_0_pressure*
                            ghostcell_dim_1_pressure;
                    
                    const int idx_temperature = (i + num_ghosts_0_temperature) +
                        (j + num_ghosts_1_temperature)*ghostcell_dim_0_temperature +
                        (k + num_ghosts_2_temperature)*ghostcell_dim_0_temperature*
                            ghostcell_dim_1_temperature;
                    
                    rho[idx_density] = p[idx_pressure]/(R*T[idx_temperature]);
                }
            }
        }
    }
}


/*
 * Compute the density.
 */
void
EquationOfStateIdealGas::computeDensity(
    double* const rho,
    const double* const p,
    const double* const T,
    const double* const R,
    const hier::IntVector& num_ghosts_density,
    const hier::IntVector& num_ghosts_pressure,
    const hier::IntVector& num_ghosts_temperature,
    const hier::IntVector& num_ghosts_thermo_properties,
    const hier::IntVector& ghostcell_dims_density,
    const hier::IntVector& ghostcell_dims_pressure,
    const hier::IntVector& ghostcell_dims_temperature,
    const hier::IntVector& ghostcell_dims_thermo_properties,
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
        
        const int num_ghosts_0_density = num_ghosts_density[0];
        const int num_ghosts_0_pressure = num_ghosts_pressure[0];
        const int num_ghosts_0_temperature = num_ghosts_temperature[0];
        const int num_ghosts_0_thermo_properties = num_ghosts_thermo_properties[0];
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_density = i + num_ghosts_0_density;
            const int idx_pressure = i + num_ghosts_0_pressure;
            const int idx_temperature = i + num_ghosts_0_temperature;
            const int idx_thermo_properties = i + num_ghosts_0_thermo_properties;
            
            rho[idx_density] = p[idx_pressure]/(R[idx_thermo_properties]*T[idx_temperature]);
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
        
        const int num_ghosts_0_density = num_ghosts_density[0];
        const int num_ghosts_1_density = num_ghosts_density[1];
        const int ghostcell_dim_0_density = ghostcell_dims_density[0];
        
        const int num_ghosts_0_pressure = num_ghosts_pressure[0];
        const int num_ghosts_1_pressure = num_ghosts_pressure[1];
        const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
        
        const int num_ghosts_0_temperature = num_ghosts_temperature[0];
        const int num_ghosts_1_temperature = num_ghosts_temperature[1];
        const int ghostcell_dim_0_temperature = ghostcell_dims_temperature[0];
        
        const int num_ghosts_0_thermo_properties = num_ghosts_thermo_properties[0];
        const int num_ghosts_1_thermo_properties = num_ghosts_thermo_properties[1];
        const int ghostcell_dim_0_thermo_properties = ghostcell_dims_thermo_properties[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_density = (i + num_ghosts_0_density) +
                    (j + num_ghosts_1_density)*ghostcell_dim_0_density;
                
                const int idx_pressure = (i + num_ghosts_0_pressure) +
                    (j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure;
                
                const int idx_temperature = (i + num_ghosts_0_temperature) +
                    (j + num_ghosts_1_temperature)*ghostcell_dim_0_temperature;
                
                const int idx_thermo_properties = (i + num_ghosts_0_thermo_properties) +
                    (j + num_ghosts_1_thermo_properties)*ghostcell_dim_0_thermo_properties;
                
                rho[idx_density] = p[idx_pressure]/(R[idx_thermo_properties]*T[idx_temperature]);
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
        
        const int num_ghosts_0_density = num_ghosts_density[0];
        const int num_ghosts_1_density = num_ghosts_density[1];
        const int num_ghosts_2_density = num_ghosts_density[2];
        const int ghostcell_dim_0_density = ghostcell_dims_density[0];
        const int ghostcell_dim_1_density = ghostcell_dims_density[1];
        
        const int num_ghosts_0_pressure = num_ghosts_pressure[0];
        const int num_ghosts_1_pressure = num_ghosts_pressure[1];
        const int num_ghosts_2_pressure = num_ghosts_pressure[2];
        const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
        const int ghostcell_dim_1_pressure = ghostcell_dims_pressure[1];
        
        const int num_ghosts_0_temperature = num_ghosts_temperature[0];
        const int num_ghosts_1_temperature = num_ghosts_temperature[1];
        const int num_ghosts_2_temperature = num_ghosts_temperature[2];
        const int ghostcell_dim_0_temperature = ghostcell_dims_temperature[0];
        const int ghostcell_dim_1_temperature = ghostcell_dims_temperature[1];
        
        const int num_ghosts_0_thermo_properties = num_ghosts_thermo_properties[0];
        const int num_ghosts_1_thermo_properties = num_ghosts_thermo_properties[1];
        const int num_ghosts_2_thermo_properties = num_ghosts_thermo_properties[2];
        const int ghostcell_dim_0_thermo_properties = ghostcell_dims_thermo_properties[0];
        const int ghostcell_dim_1_thermo_properties = ghostcell_dims_thermo_properties[1];
        
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
                    const int idx_density = (i + num_ghosts_0_density) +
                        (j + num_ghosts_1_density)*ghostcell_dim_0_density +
                        (k + num_ghosts_2_density)*ghostcell_dim_0_density*
                            ghostcell_dim_1_density;
                    
                    const int idx_pressure = (i + num_ghosts_0_pressure) +
                        (j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure +
                        (k + num_ghosts_2_pressure)*ghostcell_dim_0_pressure*
                            ghostcell_dim_1_pressure;
                    
                    const int idx_temperature = (i + num_ghosts_0_temperature) +
                        (j + num_ghosts_1_temperature)*ghostcell_dim_0_temperature +
                        (k + num_ghosts_2_temperature)*ghostcell_dim_0_temperature*
                            ghostcell_dim_1_temperature;
                    
                    const int idx_thermo_properties = (i + num_ghosts_0_thermo_properties) +
                        (j + num_ghosts_1_thermo_properties)*ghostcell_dim_0_thermo_properties +
                        (k + num_ghosts_2_thermo_properties)*ghostcell_dim_0_thermo_properties*
                            ghostcell_dim_1_thermo_properties;
                    
                    rho[idx_density] = p[idx_pressure]/(R[idx_thermo_properties]*T[idx_temperature]);
                }
            }
        }
    }
}
