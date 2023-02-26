#include "util/mixing_rules/equations_of_state/stiffened_gas/EquationOfStateStiffenedGas.hpp"

/*
 * Print all characteristics of the equation of state class.
 */
void
EquationOfStateStiffenedGas::printClassData(
    std::ostream& os) const
{
    os << "\nPrint EquationOfStateStiffenedGas object..."
       << std::endl;
       
    os << std::endl;
    
    os << "EquationOfStateStiffenedGas: this = "
       << (EquationOfStateStiffenedGas *)this
       << std::endl;
    
    os << "d_object_name = "
       << d_object_name
       << std::endl;
}


/*
 * Compute the pressure.
 */
Real
EquationOfStateStiffenedGas::getPressure(
    const Real* const density,
    const Real* const internal_energy,
    const std::vector<const Real*>& thermo_properties) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 2);
#endif
    
    const Real& gamma = *(thermo_properties[0]);
    const Real& p_inf = *(thermo_properties[1]);
    
    const Real& rho = *density;
    const Real& epsilon = *internal_energy;
    
    return ((gamma - Real(1))*rho*epsilon - gamma*p_inf); // Return p.
}


/*
 * Compute the pressure.
 */
void
EquationOfStateStiffenedGas::computePressure(
    HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_internal_energy,
    const std::vector<const Real*>& thermo_properties,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_internal_energy);
    
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 2);
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_pressure = data_pressure->getGhostBox();
    const hier::IntVector ghostcell_dims_pressure = ghost_box_pressure.numberCells();
    
    const hier::Box ghost_box_density = data_density->getGhostBox();
    const hier::IntVector ghostcell_dims_density = ghost_box_density.numberCells();
    
    const hier::Box ghost_box_internal_energy = data_internal_energy->getGhostBox();
    const hier::IntVector ghostcell_dims_internal_energy = ghost_box_internal_energy.numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_pressure(d_dim);
    hier::IntVector offset_density(d_dim);
    hier::IntVector offset_internal_energy(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_internal_energy = data_internal_energy->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_pressure->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_internal_energy->getBox().isSpatiallyEqual(interior_box));
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_pressure;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_internal_energy, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_pressure = num_ghosts_pressure;
        offset_density = num_ghosts_density;
        offset_internal_energy = num_ghosts_internal_energy;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_internal_energy->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_pressure = domain.lower() - ghost_box_pressure.lower();
        offset_density = domain.lower() - ghost_box_density.lower();
        offset_internal_energy = domain.lower() - ghost_box_internal_energy.lower();
    }
    
    /*
     * Get the pointers to the cell data.
     */
    
    Real* const p = data_pressure->getPointer(0);
    const Real* const rho = data_density->getPointer(0);
    const Real* const epsilon = data_internal_energy->getPointer(0);
    
    const Real& gamma = *(thermo_properties[0]);
    const Real& p_inf = *(thermo_properties[1]);
    
    computePressure(
        p,
        rho,
        epsilon,
        gamma,
        p_inf,
        offset_pressure,
        offset_density,
        offset_internal_energy,
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
EquationOfStateStiffenedGas::computePressure(
    HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_internal_energy,
    const std::vector<const Real*>& thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_internal_energy);
    
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 2);
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_pressure = data_pressure->getGhostBox();
    hier::IntVector ghostcell_dims_pressure = ghost_box_pressure.numberCells();
    
    const hier::Box ghost_box_density = data_density->getGhostBox();
    hier::IntVector ghostcell_dims_density = ghost_box_density.numberCells();
    
    const hier::Box ghost_box_internal_energy = data_internal_energy->getGhostBox();
    hier::IntVector ghostcell_dims_internal_energy = ghost_box_internal_energy.numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_pressure(d_dim);
    hier::IntVector offset_density(d_dim);
    hier::IntVector offset_internal_energy(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_internal_energy = data_internal_energy->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_pressure->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_internal_energy->getBox().isSpatiallyEqual(interior_box));
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_pressure;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_internal_energy, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_pressure = num_ghosts_pressure;
        offset_density = num_ghosts_density;
        offset_internal_energy = num_ghosts_internal_energy;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_internal_energy->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_pressure = domain.lower() - ghost_box_pressure.lower();
        offset_density = domain.lower() - ghost_box_density.lower();
        offset_internal_energy = domain.lower() - ghost_box_internal_energy.lower();
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
    
    Real* const p = data_pressure->getPointer(side_normal, 0);
    const Real* const rho = data_density->getPointer(side_normal, 0);
    const Real* const epsilon = data_internal_energy->getPointer(side_normal, 0);
    
    const Real& gamma = *(thermo_properties[0]);
    const Real& p_inf = *(thermo_properties[1]);
    
    computePressure(
        p,
        rho,
        epsilon,
        gamma,
        p_inf,
        offset_pressure,
        offset_density,
        offset_internal_energy,
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
EquationOfStateStiffenedGas::computePressure(
    HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_internal_energy,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_thermo_properties,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_internal_energy);
    TBOX_ASSERT(data_thermo_properties);
    
    TBOX_ASSERT(data_thermo_properties->getDepth() >= 2);
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_pressure = data_pressure->getGhostBox();
    const hier::IntVector ghostcell_dims_pressure = ghost_box_pressure.numberCells();
    
    const hier::Box ghost_box_density = data_density->getGhostBox();
    const hier::IntVector ghostcell_dims_density = ghost_box_density.numberCells();
    
    const hier::Box ghost_box_internal_energy = data_internal_energy->getGhostBox();
    const hier::IntVector ghostcell_dims_internal_energy = ghost_box_internal_energy.numberCells();
    
    const hier::Box ghost_box_thermo_properties = data_thermo_properties->getGhostBox();
    const hier::IntVector ghostcell_dims_thermo_properties = ghost_box_thermo_properties.numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_pressure(d_dim);
    hier::IntVector offset_density(d_dim);
    hier::IntVector offset_internal_energy(d_dim);
    hier::IntVector offset_thermo_properties(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_internal_energy = data_internal_energy->getGhostCellWidth();
        const hier::IntVector num_ghosts_thermo_properties = data_thermo_properties->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_pressure->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_internal_energy->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_thermo_properties->getBox().isSpatiallyEqual(interior_box));
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_pressure;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_internal_energy, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_thermo_properties, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_pressure = num_ghosts_pressure;
        offset_density = num_ghosts_density;
        offset_internal_energy = num_ghosts_internal_energy;
        offset_thermo_properties = num_ghosts_thermo_properties;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_internal_energy->getGhostBox().contains(domain));
        TBOX_ASSERT(data_thermo_properties->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_pressure = domain.lower() - ghost_box_pressure.lower();
        offset_density = domain.lower() - ghost_box_density.lower();
        offset_internal_energy = domain.lower() - ghost_box_internal_energy.lower();
        offset_thermo_properties = domain.lower() - ghost_box_thermo_properties.lower();
    }
    
    /*
     * Get the pointers to the cell data.
     */
    
    Real* p = data_pressure->getPointer(0);
    const Real* const rho = data_density->getPointer(0);
    const Real* const epsilon = data_internal_energy->getPointer(0);
    const Real* const gamma = data_thermo_properties->getPointer(0);
    const Real* const p_inf = data_thermo_properties->getPointer(1);
    
    computePressure(
        p,
        rho,
        epsilon,
        gamma,
        p_inf,
        offset_pressure,
        offset_density,
        offset_internal_energy,
        offset_thermo_properties,
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
EquationOfStateStiffenedGas::computePressure(
    HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_internal_energy,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_internal_energy);
    TBOX_ASSERT(data_thermo_properties);
    
    TBOX_ASSERT(data_thermo_properties->getDepth() >= 2);
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_pressure = data_pressure->getGhostBox();
    hier::IntVector ghostcell_dims_pressure = ghost_box_pressure.numberCells();
    
    const hier::Box ghost_box_density = data_density->getGhostBox();
    hier::IntVector ghostcell_dims_density = ghost_box_density.numberCells();
    
    const hier::Box ghost_box_internal_energy = data_internal_energy->getGhostBox();
    hier::IntVector ghostcell_dims_internal_energy = ghost_box_internal_energy.numberCells();
    
    const hier::Box ghost_box_thermo_properties = data_thermo_properties->getGhostBox();
    hier::IntVector ghostcell_dims_thermo_properties = ghost_box_thermo_properties.numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_pressure(d_dim);
    hier::IntVector offset_density(d_dim);
    hier::IntVector offset_internal_energy(d_dim);
    hier::IntVector offset_thermo_properties(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_internal_energy = data_internal_energy->getGhostCellWidth();
        const hier::IntVector num_ghosts_thermo_properties = data_thermo_properties->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_pressure->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_internal_energy->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_thermo_properties->getBox().isSpatiallyEqual(interior_box));
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_pressure;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_internal_energy, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_thermo_properties, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_pressure = num_ghosts_pressure;
        offset_density = num_ghosts_density;
        offset_internal_energy = num_ghosts_internal_energy;
        offset_thermo_properties = num_ghosts_thermo_properties;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_internal_energy->getGhostBox().contains(domain));
        TBOX_ASSERT(data_thermo_properties->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_pressure = domain.lower() - ghost_box_pressure.lower();
        offset_density = domain.lower() - ghost_box_density.lower();
        offset_internal_energy = domain.lower() - ghost_box_internal_energy.lower();
        offset_thermo_properties = domain.lower() - ghost_box_thermo_properties.lower();
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
    
    Real* p = data_pressure->getPointer(side_normal, 0);
    const Real* const rho = data_density->getPointer(side_normal, 0);
    const Real* const epsilon = data_internal_energy->getPointer(side_normal, 0);
    const Real* const gamma = data_thermo_properties->getPointer(side_normal, 0);
    const Real* const p_inf = data_thermo_properties->getPointer(side_normal, 1);
    
    computePressure(
        p,
        rho,
        epsilon,
        gamma,
        p_inf,
        offset_pressure,
        offset_density,
        offset_internal_energy,
        offset_thermo_properties,
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
Real
EquationOfStateStiffenedGas::getSoundSpeed(
    const Real* const density,
    const Real* const pressure,
    const std::vector<const Real*>& thermo_properties) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 2);
#endif
    
    const Real& gamma = *(thermo_properties[0]);
    const Real& p_inf = *(thermo_properties[1]);
    
    const Real& rho = *density;
    const Real& p = *pressure;
    
    return std::sqrt(gamma*(p + p_inf)/rho); // Return c.
}


/*
 * Compute the sound speed.
 */
void
EquationOfStateStiffenedGas::computeSoundSpeed(
    HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_sound_speed,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
    const std::vector<const Real*>& thermo_properties,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_sound_speed);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 2);
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_sound_speed = data_sound_speed->getGhostBox();
    const hier::IntVector ghostcell_dims_sound_speed = ghost_box_sound_speed.numberCells();
    
    const hier::Box ghost_box_density = data_density->getGhostBox();
    const hier::IntVector ghostcell_dims_density = ghost_box_density.numberCells();
    
    const hier::Box ghost_box_pressure = data_pressure->getGhostBox();
    const hier::IntVector ghostcell_dims_pressure = ghost_box_pressure.numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_sound_speed(d_dim);
    hier::IntVector offset_density(d_dim);
    hier::IntVector offset_pressure(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_sound_speed = data_sound_speed->getGhostCellWidth();
        const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_sound_speed->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_pressure->getBox().isSpatiallyEqual(interior_box));
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_sound_speed;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_sound_speed = num_ghosts_sound_speed;
        offset_density = num_ghosts_density;
        offset_pressure = num_ghosts_pressure;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_sound_speed->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_sound_speed = domain.lower() - ghost_box_sound_speed.lower();
        offset_density = domain.lower() - ghost_box_density.lower();
        offset_pressure = domain.lower() - ghost_box_pressure.lower();
    }
    
    /*
     * Get the pointers to the cell data.
     */
    
    Real* const c = data_sound_speed->getPointer(0);
    const Real* const rho = data_density->getPointer(0);
    const Real* const p = data_pressure->getPointer(0);
    
    const Real& gamma = *(thermo_properties[0]);
    const Real& p_inf = *(thermo_properties[1]);
    
    computeSoundSpeed(
        c,
        rho,
        p,
        gamma,
        p_inf,
        offset_sound_speed,
        offset_density,
        offset_pressure,
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
EquationOfStateStiffenedGas::computeSoundSpeed(
    HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_sound_speed,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
    const std::vector<const Real*>& thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_sound_speed);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 2);
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_sound_speed = data_sound_speed->getGhostBox();
    hier::IntVector ghostcell_dims_sound_speed = ghost_box_sound_speed.numberCells();
    
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
    
    hier::IntVector offset_sound_speed(d_dim);
    hier::IntVector offset_density(d_dim);
    hier::IntVector offset_pressure(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_sound_speed = data_sound_speed->getGhostCellWidth();
        const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_sound_speed->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_pressure->getBox().isSpatiallyEqual(interior_box));
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_sound_speed;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_sound_speed = num_ghosts_sound_speed;
        offset_density = num_ghosts_density;
        offset_pressure = num_ghosts_pressure;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_sound_speed->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_sound_speed = domain.lower() - ghost_box_sound_speed.lower();
        offset_density = domain.lower() - ghost_box_density.lower();
        offset_pressure = domain.lower() - ghost_box_pressure.lower();
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
    
    Real* const c = data_sound_speed->getPointer(side_normal, 0);
    const Real* const rho = data_density->getPointer(side_normal, 0);
    const Real* const p = data_pressure->getPointer(side_normal, 0);
    
    const Real& gamma = *(thermo_properties[0]);
    const Real& p_inf = *(thermo_properties[1]);
    
    computeSoundSpeed(
        c,
        rho,
        p,
        gamma,
        p_inf,
        offset_sound_speed,
        offset_density,
        offset_pressure,
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
EquationOfStateStiffenedGas::computeSoundSpeed(
    HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_sound_speed,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_thermo_properties,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_sound_speed);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_thermo_properties);
    
    TBOX_ASSERT(data_thermo_properties->getDepth() >= 2);
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_sound_speed = data_sound_speed->getGhostBox();
    const hier::IntVector ghostcell_dims_sound_speed = ghost_box_sound_speed.numberCells();
    
    const hier::Box ghost_box_density = data_density->getGhostBox();
    const hier::IntVector ghostcell_dims_density = ghost_box_density.numberCells();
    
    const hier::Box ghost_box_pressure = data_pressure->getGhostBox();
    const hier::IntVector ghostcell_dims_pressure = ghost_box_pressure.numberCells();
    
    const hier::Box ghost_box_thermo_properties = data_thermo_properties->getGhostBox();
    const hier::IntVector ghostcell_dims_thermo_properties = ghost_box_thermo_properties.numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_sound_speed(d_dim);
    hier::IntVector offset_density(d_dim);
    hier::IntVector offset_pressure(d_dim);
    hier::IntVector offset_thermo_properties(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_sound_speed = data_sound_speed->getGhostCellWidth();
        const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        const hier::IntVector num_ghosts_thermo_properties = data_thermo_properties->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_sound_speed->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_pressure->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_thermo_properties->getBox().isSpatiallyEqual(interior_box));
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_sound_speed;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_thermo_properties, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_sound_speed = num_ghosts_sound_speed;
        offset_density = num_ghosts_density;
        offset_pressure = num_ghosts_pressure;
        offset_thermo_properties = num_ghosts_thermo_properties;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_sound_speed->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_thermo_properties->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_sound_speed = domain.lower() - ghost_box_sound_speed.lower();
        offset_density = domain.lower() - ghost_box_density.lower();
        offset_pressure = domain.lower() - ghost_box_pressure.lower();
        offset_thermo_properties = domain.lower() - ghost_box_thermo_properties.lower();
    }
    
    /*
     * Get the pointers to the cell data.
     */
    
    Real* const c = data_sound_speed->getPointer(0);
    const Real* const rho = data_density->getPointer(0);
    const Real* const p = data_pressure->getPointer(0);
    const Real* const gamma = data_thermo_properties->getPointer(0);
    const Real* const p_inf = data_thermo_properties->getPointer(1);
    
    computeSoundSpeed(
        c,
        rho,
        p,
        gamma,
        p_inf,
        offset_sound_speed,
        offset_density,
        offset_pressure,
        offset_thermo_properties,
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
EquationOfStateStiffenedGas::computeSoundSpeed(
    HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_sound_speed,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_sound_speed);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_thermo_properties);
    
    TBOX_ASSERT(data_thermo_properties->getDepth() >= 2);
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_sound_speed = data_sound_speed->getGhostBox();
    hier::IntVector ghostcell_dims_sound_speed = ghost_box_sound_speed.numberCells();
    
    const hier::Box ghost_box_density = data_density->getGhostBox();
    hier::IntVector ghostcell_dims_density = ghost_box_density.numberCells();
    
    const hier::Box ghost_box_pressure = data_pressure->getGhostBox();
    hier::IntVector ghostcell_dims_pressure = ghost_box_pressure.numberCells();
    
    const hier::Box ghost_box_thermo_properties = data_thermo_properties->getGhostBox();
    hier::IntVector ghostcell_dims_thermo_properties = ghost_box_thermo_properties.numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_sound_speed(d_dim);
    hier::IntVector offset_density(d_dim);
    hier::IntVector offset_pressure(d_dim);
    hier::IntVector offset_thermo_properties(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_sound_speed = data_sound_speed->getGhostCellWidth();
        const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        const hier::IntVector num_ghosts_thermo_properties = data_thermo_properties->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_sound_speed->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_pressure->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_thermo_properties->getBox().isSpatiallyEqual(interior_box));
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_sound_speed;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_thermo_properties, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_sound_speed = num_ghosts_sound_speed;
        offset_density = num_ghosts_density;
        offset_pressure = num_ghosts_pressure;
        offset_thermo_properties = num_ghosts_thermo_properties;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_sound_speed->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_thermo_properties->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_sound_speed = domain.lower() - ghost_box_sound_speed.lower();
        offset_density = domain.lower() - ghost_box_density.lower();
        offset_pressure = domain.lower() - ghost_box_pressure.lower();
        offset_thermo_properties = domain.lower() - ghost_box_thermo_properties.lower();
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
    
    Real* const c = data_sound_speed->getPointer(side_normal, 0);
    const Real* const rho = data_density->getPointer(side_normal, 0);
    const Real* const p = data_pressure->getPointer(side_normal, 0);
    const Real* const gamma = data_thermo_properties->getPointer(side_normal, 0);
    const Real* const p_inf = data_thermo_properties->getPointer(side_normal, 1);
    
    computeSoundSpeed(
        c,
        rho,
        p,
        gamma,
        p_inf,
        offset_sound_speed,
        offset_density,
        offset_pressure,
        offset_thermo_properties,
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
Real
EquationOfStateStiffenedGas::getInternalEnergy(
    const Real* const density,
    const Real* const pressure,
    const std::vector<const Real*>& thermo_properties) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 2);
#endif
    
    const Real& gamma = *(thermo_properties[0]);
    const Real& p_inf = *(thermo_properties[1]);
    
    const Real& rho = *density;
    const Real& p = *pressure;
    
    return (p + gamma*p_inf)/((gamma - Real(1))*rho); // Return epsilon.
}


/*
 * Compute the specific internal energy.
 */
void
EquationOfStateStiffenedGas::computeInternalEnergy(
    HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_internal_energy,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
    const std::vector<const Real*>& thermo_properties,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_internal_energy);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 2);
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_internal_energy = data_internal_energy->getGhostBox();
    const hier::IntVector ghostcell_dims_internal_energy = ghost_box_internal_energy.numberCells();
    
    const hier::Box ghost_box_density = data_density->getGhostBox();
    const hier::IntVector ghostcell_dims_density = ghost_box_density.numberCells();
    
    const hier::Box ghost_box_pressure = data_pressure->getGhostBox();
    const hier::IntVector ghostcell_dims_pressure = ghost_box_pressure.numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_internal_energy(d_dim);
    hier::IntVector offset_density(d_dim);
    hier::IntVector offset_pressure(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_internal_energy = data_internal_energy->getGhostCellWidth();
        const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_internal_energy->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_pressure->getBox().isSpatiallyEqual(interior_box));
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_internal_energy;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_internal_energy = num_ghosts_internal_energy;
        offset_density = num_ghosts_density;
        offset_pressure = num_ghosts_pressure;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_internal_energy->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_internal_energy = domain.lower() - ghost_box_internal_energy.lower();
        offset_density = domain.lower() - ghost_box_density.lower();
        offset_pressure = domain.lower() - ghost_box_pressure.lower();
    }
    
    /*
     * Get the pointers to the cell data.
     */
    
    Real* const epsilon = data_internal_energy->getPointer(0);
    const Real* const rho = data_density->getPointer(0);
    const Real* const p = data_pressure->getPointer(0);
    
    const Real& gamma = *(thermo_properties[0]);
    const Real& p_inf = *(thermo_properties[1]);
    
    computeInternalEnergy(
        epsilon,
        rho,
        p,
        gamma,
        p_inf,
        offset_internal_energy,
        offset_density,
        offset_pressure,
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
EquationOfStateStiffenedGas::computeInternalEnergy(
    HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_internal_energy,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
    const std::vector<const Real*>& thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_internal_energy);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 2);
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_internal_energy = data_internal_energy->getGhostBox();
    hier::IntVector ghostcell_dims_internal_energy = ghost_box_internal_energy.numberCells();
    
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
    
    hier::IntVector offset_internal_energy(d_dim);
    hier::IntVector offset_density(d_dim);
    hier::IntVector offset_pressure(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_internal_energy = data_internal_energy->getGhostCellWidth();
        const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_internal_energy->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_pressure->getBox().isSpatiallyEqual(interior_box));
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_internal_energy;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_internal_energy = num_ghosts_internal_energy;
        offset_density = num_ghosts_density;
        offset_pressure = num_ghosts_pressure;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_internal_energy->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_internal_energy = domain.lower() - ghost_box_internal_energy.lower();
        offset_density = domain.lower() - ghost_box_density.lower();
        offset_pressure = domain.lower() - ghost_box_pressure.lower();
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
    
    Real* const epsilon = data_internal_energy->getPointer(side_normal, 0);
    const Real* const rho = data_density->getPointer(side_normal, 0);
    const Real* const p = data_pressure->getPointer(side_normal, 0);
    
    const Real& gamma = *(thermo_properties[0]);
    const Real& p_inf = *(thermo_properties[1]);
    
    computeInternalEnergy(
        epsilon,
        rho,
        p,
        gamma,
        p_inf,
        offset_internal_energy,
        offset_density,
        offset_pressure,
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
EquationOfStateStiffenedGas::computeInternalEnergy(
    HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_internal_energy,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_thermo_properties,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_internal_energy);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_thermo_properties);
    
    TBOX_ASSERT(data_thermo_properties->getDepth() >= 2);
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_internal_energy = data_internal_energy->getGhostBox();
    const hier::IntVector ghostcell_dims_internal_energy = ghost_box_internal_energy.numberCells();
    
    const hier::Box ghost_box_density = data_density->getGhostBox();
    const hier::IntVector ghostcell_dims_density = ghost_box_density.numberCells();
    
    const hier::Box ghost_box_pressure = data_pressure->getGhostBox();
    const hier::IntVector ghostcell_dims_pressure = ghost_box_pressure.numberCells();
    
    const hier::Box ghost_box_thermo_properties = data_thermo_properties->getGhostBox();
    const hier::IntVector ghostcell_dims_thermo_properties = ghost_box_thermo_properties.numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_internal_energy(d_dim);
    hier::IntVector offset_density(d_dim);
    hier::IntVector offset_pressure(d_dim);
    hier::IntVector offset_thermo_properties(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_internal_energy = data_internal_energy->getGhostCellWidth();
        const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        const hier::IntVector num_ghosts_thermo_properties = data_thermo_properties->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_internal_energy->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_pressure->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_thermo_properties->getBox().isSpatiallyEqual(interior_box));
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_internal_energy;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_thermo_properties, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_internal_energy = num_ghosts_internal_energy;
        offset_density = num_ghosts_density;
        offset_pressure = num_ghosts_pressure;
        offset_thermo_properties = num_ghosts_thermo_properties;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_internal_energy->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_thermo_properties->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_internal_energy = domain.lower() - ghost_box_internal_energy.lower();
        offset_density = domain.lower() - ghost_box_density.lower();
        offset_pressure = domain.lower() - ghost_box_pressure.lower();
        offset_thermo_properties = domain.lower() - ghost_box_thermo_properties.lower();
    }
    
    /*
     * Get the pointers to the cell data.
     */
    
    Real* const epsilon = data_internal_energy->getPointer(0);
    const Real* const rho = data_density->getPointer(0);
    const Real* const p = data_pressure->getPointer(0);
    const Real* const gamma = data_thermo_properties->getPointer(0);
    const Real* const p_inf = data_thermo_properties->getPointer(1);
    
    computeInternalEnergy(
        epsilon,
        rho,
        p,
        gamma,
        p_inf,
        offset_internal_energy,
        offset_density,
        offset_pressure,
        offset_thermo_properties,
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
EquationOfStateStiffenedGas::computeInternalEnergy(
    HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_internal_energy,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_internal_energy);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_thermo_properties);
    
    TBOX_ASSERT(data_thermo_properties->getDepth() >= 2);
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_internal_energy = data_internal_energy->getGhostBox();
    hier::IntVector ghostcell_dims_internal_energy = ghost_box_internal_energy.numberCells();
    
    const hier::Box ghost_box_density = data_density->getGhostBox();
    hier::IntVector ghostcell_dims_density = ghost_box_density.numberCells();
    
    const hier::Box ghost_box_pressure = data_pressure->getGhostBox();
    hier::IntVector ghostcell_dims_pressure = ghost_box_pressure.numberCells();
    
    const hier::Box ghost_box_thermo_properties = data_thermo_properties->getGhostBox();
    hier::IntVector ghostcell_dims_thermo_properties = ghost_box_thermo_properties.numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_internal_energy(d_dim);
    hier::IntVector offset_density(d_dim);
    hier::IntVector offset_pressure(d_dim);
    hier::IntVector offset_thermo_properties(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_internal_energy = data_internal_energy->getGhostCellWidth();
        const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        const hier::IntVector num_ghosts_thermo_properties = data_thermo_properties->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_internal_energy->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_pressure->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_thermo_properties->getBox().isSpatiallyEqual(interior_box));
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_internal_energy;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_thermo_properties, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_internal_energy = num_ghosts_internal_energy;
        offset_density = num_ghosts_density;
        offset_pressure = num_ghosts_pressure;
        offset_thermo_properties = num_ghosts_thermo_properties;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_internal_energy->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_thermo_properties->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_internal_energy = domain.lower() - ghost_box_internal_energy.lower();
        offset_density = domain.lower() - ghost_box_density.lower();
        offset_pressure = domain.lower() - ghost_box_pressure.lower();
        offset_thermo_properties = domain.lower() - ghost_box_thermo_properties.lower();
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
    
    Real* const epsilon = data_internal_energy->getPointer(side_normal, 0);
    const Real* const rho = data_density->getPointer(side_normal, 0);
    const Real* const p = data_pressure->getPointer(side_normal, 0);
    const Real* const gamma = data_thermo_properties->getPointer(side_normal, 0);
    const Real* const p_inf = data_thermo_properties->getPointer(side_normal, 1);
    
    computeInternalEnergy(
        epsilon,
        rho,
        p,
        gamma,
        p_inf,
        offset_internal_energy,
        offset_density,
        offset_pressure,
        offset_thermo_properties,
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
Real
EquationOfStateStiffenedGas::getEnthalpy(
    const Real* const density,
    const Real* const pressure,
    const std::vector<const Real*>& thermo_properties) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 2);
#endif
    
    const Real& gamma = *(thermo_properties[0]);
    const Real& p_inf = *(thermo_properties[1]);
    
    const Real& rho = *density;
    const Real& p = *pressure;
    
    return gamma*(p + p_inf)/((gamma - Real(1))*rho); // Return h.
}


/*
 * Compute the specific enthalpy.
 */
void
EquationOfStateStiffenedGas::computeEnthalpy(
    HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_enthalpy,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
    const std::vector<const Real*>& thermo_properties,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_enthalpy);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 2);
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_enthalpy = data_enthalpy->getGhostBox();
    const hier::IntVector ghostcell_dims_enthalpy = ghost_box_enthalpy.numberCells();
    
    const hier::Box ghost_box_density = data_density->getGhostBox();
    const hier::IntVector ghostcell_dims_density = ghost_box_density.numberCells();
    
    const hier::Box ghost_box_pressure = data_pressure->getGhostBox();
    const hier::IntVector ghostcell_dims_pressure = ghost_box_pressure.numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_enthalpy(d_dim);
    hier::IntVector offset_density(d_dim);
    hier::IntVector offset_pressure(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_enthalpy = data_enthalpy->getGhostCellWidth();
        const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_enthalpy->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_pressure->getBox().isSpatiallyEqual(interior_box));
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_enthalpy;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_enthalpy = num_ghosts_enthalpy;
        offset_density = num_ghosts_density;
        offset_pressure = num_ghosts_pressure;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_enthalpy->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_enthalpy = domain.lower() - ghost_box_enthalpy.lower();
        offset_density = domain.lower() - ghost_box_density.lower();
        offset_pressure = domain.lower() - ghost_box_pressure.lower();
    }
    
    /*
     * Get the pointers to the cell data.
     */
    
    Real* const h = data_enthalpy->getPointer(0);
    const Real* const rho = data_density->getPointer(0);
    const Real* const p = data_pressure->getPointer(0);
    
    const Real& gamma = *(thermo_properties[0]);
    const Real& p_inf = *(thermo_properties[1]);
    
    computeEnthalpy(
        h,
        rho,
        p,
        gamma,
        p_inf,
        offset_enthalpy,
        offset_density,
        offset_pressure,
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
EquationOfStateStiffenedGas::computeEnthalpy(
    HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_enthalpy,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
    const std::vector<const Real*>& thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_enthalpy);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 2);
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_enthalpy = data_enthalpy->getGhostBox();
    hier::IntVector ghostcell_dims_enthalpy = ghost_box_enthalpy.numberCells();
    
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
    
    hier::IntVector offset_enthalpy(d_dim);
    hier::IntVector offset_density(d_dim);
    hier::IntVector offset_pressure(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_enthalpy = data_enthalpy->getGhostCellWidth();
        const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_enthalpy->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_pressure->getBox().isSpatiallyEqual(interior_box));
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_enthalpy;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_enthalpy = num_ghosts_enthalpy;
        offset_density = num_ghosts_density;
        offset_pressure = num_ghosts_pressure;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_enthalpy->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_enthalpy = domain.lower() - ghost_box_enthalpy.lower();
        offset_density = domain.lower() - ghost_box_density.lower();
        offset_pressure = domain.lower() - ghost_box_pressure.lower();
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
    
    Real* const h = data_enthalpy->getPointer(side_normal, 0);
    const Real* const rho = data_density->getPointer(side_normal, 0);
    const Real* const p = data_pressure->getPointer(side_normal, 0);
    
    const Real& gamma = *(thermo_properties[0]);
    const Real& p_inf = *(thermo_properties[1]);
    
    computeEnthalpy(
        h,
        rho,
        p,
        gamma,
        p_inf,
        offset_enthalpy,
        offset_density,
        offset_pressure,
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
EquationOfStateStiffenedGas::computeEnthalpy(
    HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_enthalpy,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_thermo_properties,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_enthalpy);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_thermo_properties);
    
    TBOX_ASSERT(data_thermo_properties->getDepth() >= 2);
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_enthalpy = data_enthalpy->getGhostBox();
    const hier::IntVector ghostcell_dims_enthalpy = ghost_box_enthalpy.numberCells();
    
    const hier::Box ghost_box_density = data_density->getGhostBox();
    const hier::IntVector ghostcell_dims_density = ghost_box_density.numberCells();
    
    const hier::Box ghost_box_pressure = data_pressure->getGhostBox();
    const hier::IntVector ghostcell_dims_pressure = ghost_box_pressure.numberCells();
    
    const hier::Box ghost_box_thermo_properties = data_thermo_properties->getGhostBox();
    const hier::IntVector ghostcell_dims_thermo_properties = ghost_box_thermo_properties.numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_enthalpy(d_dim);
    hier::IntVector offset_density(d_dim);
    hier::IntVector offset_pressure(d_dim);
    hier::IntVector offset_thermo_properties(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_enthalpy = data_enthalpy->getGhostCellWidth();
        const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        const hier::IntVector num_ghosts_thermo_properties = data_thermo_properties->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_enthalpy->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_pressure->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_thermo_properties->getBox().isSpatiallyEqual(interior_box));
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_enthalpy;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_thermo_properties, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_enthalpy = num_ghosts_enthalpy;
        offset_density = num_ghosts_density;
        offset_pressure = num_ghosts_pressure;
        offset_thermo_properties = num_ghosts_thermo_properties;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_enthalpy->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_thermo_properties->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_enthalpy = domain.lower() - ghost_box_enthalpy.lower();
        offset_density = domain.lower() - ghost_box_density.lower();
        offset_pressure = domain.lower() - ghost_box_pressure.lower();
        offset_thermo_properties = domain.lower() - ghost_box_thermo_properties.lower();
    }
    
    /*
     * Get the pointers to the cell data.
     */
    
    Real* const h = data_enthalpy->getPointer(0);
    const Real* const rho = data_density->getPointer(0);
    const Real* const p = data_pressure->getPointer(0);
    const Real* const gamma = data_thermo_properties->getPointer(0);
    const Real* const p_inf = data_thermo_properties->getPointer(1);
    
    computeEnthalpy(
        h,
        rho,
        p,
        gamma,
        p_inf,
        offset_enthalpy,
        offset_density,
        offset_pressure,
        offset_thermo_properties,
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
EquationOfStateStiffenedGas::computeEnthalpy(
    HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_enthalpy,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_enthalpy);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_thermo_properties);
    
    TBOX_ASSERT(data_thermo_properties->getDepth() >= 2);
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_enthalpy = data_enthalpy->getGhostBox();
    hier::IntVector ghostcell_dims_enthalpy = ghost_box_enthalpy.numberCells();
    
    const hier::Box ghost_box_density = data_density->getGhostBox();
    hier::IntVector ghostcell_dims_density = ghost_box_density.numberCells();
    
    const hier::Box ghost_box_pressure = data_pressure->getGhostBox();
    hier::IntVector ghostcell_dims_pressure = ghost_box_pressure.numberCells();
    
    const hier::Box ghost_box_thermo_properties = data_thermo_properties->getGhostBox();
    hier::IntVector ghostcell_dims_thermo_properties = ghost_box_thermo_properties.numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_enthalpy(d_dim);
    hier::IntVector offset_density(d_dim);
    hier::IntVector offset_pressure(d_dim);
    hier::IntVector offset_thermo_properties(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_enthalpy = data_enthalpy->getGhostCellWidth();
        const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        const hier::IntVector num_ghosts_thermo_properties = data_thermo_properties->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_enthalpy->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_pressure->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_thermo_properties->getBox().isSpatiallyEqual(interior_box));
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_enthalpy;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_thermo_properties, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_enthalpy = num_ghosts_enthalpy;
        offset_density = num_ghosts_density;
        offset_pressure = num_ghosts_pressure;
        offset_thermo_properties = num_ghosts_thermo_properties;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_enthalpy->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_thermo_properties->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_enthalpy = domain.lower() - ghost_box_enthalpy.lower();
        offset_density = domain.lower() - ghost_box_density.lower();
        offset_pressure = domain.lower() - ghost_box_pressure.lower();
        offset_thermo_properties = domain.lower() - ghost_box_thermo_properties.lower();
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
    
    Real* const h = data_enthalpy->getPointer(side_normal, 0);
    const Real* const rho = data_density->getPointer(side_normal, 0);
    const Real* const p = data_pressure->getPointer(side_normal, 0);
    const Real* const gamma = data_thermo_properties->getPointer(side_normal, 0);
    const Real* const p_inf = data_thermo_properties->getPointer(side_normal, 1);
    
    computeEnthalpy(
        h,
        rho,
        p,
        gamma,
        p_inf,
        offset_enthalpy,
        offset_density,
        offset_pressure,
        offset_thermo_properties,
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
Real
EquationOfStateStiffenedGas::getTemperature(
    const Real* const density,
    const Real* const pressure,
    const std::vector<const Real*>& thermo_properties) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) == 4);
#endif
    
    const Real& gamma = *(thermo_properties[0]);
    const Real& p_inf = *(thermo_properties[1]);
    const Real& c_v = *(thermo_properties[3]);
    
    const Real& rho = *density;
    const Real& p = *pressure;
    
    return (p + p_inf)/((gamma - Real(1))*c_v*rho); // Return T.
}


/*
 * Compute the temperature.
 */
void
EquationOfStateStiffenedGas::computeTemperature(
    HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_temperature,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
    const std::vector<const Real*>& thermo_properties,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_temperature);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) == 4);
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_temperature = data_temperature->getGhostBox();
    const hier::IntVector ghostcell_dims_temperature = ghost_box_temperature.numberCells();
    
    const hier::Box ghost_box_density = data_density->getGhostBox();
    const hier::IntVector ghostcell_dims_density = ghost_box_density.numberCells();
    
    const hier::Box ghost_box_pressure = data_pressure->getGhostBox();
    const hier::IntVector ghostcell_dims_pressure = ghost_box_pressure.numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_temperature(d_dim);
    hier::IntVector offset_density(d_dim);
    hier::IntVector offset_pressure(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
        const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_temperature->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_pressure->getBox().isSpatiallyEqual(interior_box));
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_temperature;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_temperature = num_ghosts_temperature;
        offset_density = num_ghosts_density;
        offset_pressure = num_ghosts_pressure;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_temperature->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_temperature = domain.lower() - ghost_box_temperature.lower();
        offset_density = domain.lower() - ghost_box_density.lower();
        offset_pressure = domain.lower() - ghost_box_pressure.lower();
    }
    
    /*
     * Get the pointers to the cell data.
     */
    
    Real* const T = data_temperature->getPointer(0);
    const Real* const rho = data_density->getPointer(0);
    const Real* const p = data_pressure->getPointer(0);
    
    const Real& gamma = *(thermo_properties[0]);
    const Real& p_inf = *(thermo_properties[1]);
    const Real& c_v = *(thermo_properties[3]);
    
    computeTemperature(
        T,
        rho,
        p,
        gamma,
        p_inf,
        c_v,
        offset_temperature,
        offset_density,
        offset_pressure,
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
EquationOfStateStiffenedGas::computeTemperature(
    HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_temperature,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
    const std::vector<const Real*>& thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_temperature);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) == 4);
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_temperature = data_temperature->getGhostBox();
    hier::IntVector ghostcell_dims_temperature = ghost_box_temperature.numberCells();
    
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
    
    hier::IntVector offset_temperature(d_dim);
    hier::IntVector offset_density(d_dim);
    hier::IntVector offset_pressure(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
        const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_temperature->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_pressure->getBox().isSpatiallyEqual(interior_box));
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_temperature;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_temperature = num_ghosts_temperature;
        offset_density = num_ghosts_density;
        offset_pressure = num_ghosts_pressure;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_temperature->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_temperature = domain.lower() - ghost_box_temperature.lower();
        offset_density = domain.lower() - ghost_box_density.lower();
        offset_pressure = domain.lower() - ghost_box_pressure.lower();
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
    
    Real* const T = data_temperature->getPointer(side_normal, 0);
    const Real* const rho = data_density->getPointer(side_normal, 0);
    const Real* const p = data_pressure->getPointer(side_normal, 0);
    
    const Real& gamma = *(thermo_properties[0]);
    const Real& p_inf = *(thermo_properties[1]);
    const Real& c_v = *(thermo_properties[3]);
    
    computeTemperature(
        T,
        rho,
        p,
        gamma,
        p_inf,
        c_v,
        offset_temperature,
        offset_density,
        offset_pressure,
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
EquationOfStateStiffenedGas::computeTemperature(
    HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_temperature,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_thermo_properties,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_temperature);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_thermo_properties);
    
    TBOX_ASSERT(data_thermo_properties->getDepth() == 4);
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_temperature = data_temperature->getGhostBox();
    const hier::IntVector ghostcell_dims_temperature = ghost_box_temperature.numberCells();
    
    const hier::Box ghost_box_density = data_density->getGhostBox();
    const hier::IntVector ghostcell_dims_density = ghost_box_density.numberCells();
    
    const hier::Box ghost_box_pressure = data_pressure->getGhostBox();
    const hier::IntVector ghostcell_dims_pressure = ghost_box_pressure.numberCells();
    
    const hier::Box ghost_box_thermo_properties = data_thermo_properties->getGhostBox();
    const hier::IntVector ghostcell_dims_thermo_properties = ghost_box_thermo_properties.numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_temperature(d_dim);
    hier::IntVector offset_density(d_dim);
    hier::IntVector offset_pressure(d_dim);
    hier::IntVector offset_thermo_properties(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
        const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        const hier::IntVector num_ghosts_thermo_properties = data_thermo_properties->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_temperature->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_pressure->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_thermo_properties->getBox().isSpatiallyEqual(interior_box));
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_temperature;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_thermo_properties, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_temperature = num_ghosts_temperature;
        offset_density = num_ghosts_density;
        offset_pressure = num_ghosts_pressure;
        offset_thermo_properties = num_ghosts_thermo_properties;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_temperature->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_thermo_properties->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_temperature = domain.lower() - ghost_box_temperature.lower();
        offset_density = domain.lower() - ghost_box_density.lower();
        offset_pressure = domain.lower() - ghost_box_pressure.lower();
        offset_thermo_properties = domain.lower() - ghost_box_thermo_properties.lower();
    }
    
    /*
     * Get the pointers to the cell data.
     */
    
    Real* const T = data_temperature->getPointer(0);
    const Real* const rho = data_density->getPointer(0);
    const Real* const p = data_pressure->getPointer(0);
    const Real* const gamma = data_thermo_properties->getPointer(0);
    const Real* const p_inf = data_thermo_properties->getPointer(1);
    const Real* const c_v = data_thermo_properties->getPointer(3);
    
    computeTemperature(
        T,
        rho,
        p,
        gamma,
        p_inf,
        c_v,
        offset_temperature,
        offset_density,
        offset_pressure,
        offset_thermo_properties,
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
EquationOfStateStiffenedGas::computeTemperature(
    HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_temperature,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_thermo_properties,
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
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_temperature = data_temperature->getGhostBox();
    hier::IntVector ghostcell_dims_temperature = ghost_box_temperature.numberCells();
    
    const hier::Box ghost_box_density = data_density->getGhostBox();
    hier::IntVector ghostcell_dims_density = ghost_box_density.numberCells();
    
    const hier::Box ghost_box_pressure = data_pressure->getGhostBox();
    hier::IntVector ghostcell_dims_pressure = ghost_box_pressure.numberCells();
    
    const hier::Box ghost_box_thermo_properties = data_thermo_properties->getGhostBox();
    hier::IntVector ghostcell_dims_thermo_properties = ghost_box_thermo_properties.numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_temperature(d_dim);
    hier::IntVector offset_density(d_dim);
    hier::IntVector offset_pressure(d_dim);
    hier::IntVector offset_thermo_properties(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
        const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        const hier::IntVector num_ghosts_thermo_properties = data_thermo_properties->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_temperature->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_pressure->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_thermo_properties->getBox().isSpatiallyEqual(interior_box));
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_temperature;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_thermo_properties, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_temperature = num_ghosts_temperature;
        offset_density = num_ghosts_density;
        offset_pressure = num_ghosts_pressure;
        offset_thermo_properties = num_ghosts_thermo_properties;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_temperature->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_thermo_properties->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_temperature = domain.lower() - ghost_box_temperature.lower();
        offset_density = domain.lower() - ghost_box_density.lower();
        offset_pressure = domain.lower() - ghost_box_pressure.lower();
        offset_thermo_properties = domain.lower() - ghost_box_thermo_properties.lower();
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
    
    Real* const T = data_temperature->getPointer(side_normal, 0);
    const Real* const rho = data_density->getPointer(side_normal, 0);
    const Real* const p = data_pressure->getPointer(side_normal, 0);
    const Real* const gamma = data_thermo_properties->getPointer(side_normal, 0);
    const Real* const p_inf = data_thermo_properties->getPointer(side_normal, 1);
    const Real* const c_v = data_thermo_properties->getPointer(side_normal, 3);
    
    computeTemperature(
        T,
        rho,
        p,
        gamma,
        p_inf,
        c_v,
        offset_temperature,
        offset_density,
        offset_pressure,
        offset_thermo_properties,
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
Real
EquationOfStateStiffenedGas::getInternalEnergyFromTemperature(
    const Real* const density,
    const Real* const temperature,
    const std::vector<const Real*>& thermo_properties) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) == 4);
#endif
    
    const Real& p_inf = *(thermo_properties[1]);
    const Real& c_v = *(thermo_properties[3]);
    
    const Real& rho = *density;
    const Real& T = *temperature;
    
    return (c_v*T + p_inf/rho); // Return epsilon.
}


/*
 * Compute the specific internal energy from temperature.
 */
void
EquationOfStateStiffenedGas::computeInternalEnergyFromTemperature(
    HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_internal_energy,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_temperature,
    const std::vector<const Real*>& thermo_properties,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_internal_energy);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_temperature);
    
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) == 4);
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_internal_energy = data_internal_energy->getGhostBox();
    const hier::IntVector ghostcell_dims_internal_energy = ghost_box_internal_energy.numberCells();
    
    const hier::Box ghost_box_density = data_density->getGhostBox();
    const hier::IntVector ghostcell_dims_density = ghost_box_density.numberCells();
    
    const hier::Box ghost_box_temperature = data_temperature->getGhostBox();
    const hier::IntVector ghostcell_dims_temperature = ghost_box_temperature.numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_internal_energy(d_dim);
    hier::IntVector offset_density(d_dim);
    hier::IntVector offset_temperature(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_internal_energy = data_internal_energy->getGhostCellWidth();
        const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_internal_energy->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_temperature->getBox().isSpatiallyEqual(interior_box));
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_internal_energy;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_temperature, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_internal_energy = num_ghosts_internal_energy;
        offset_density = num_ghosts_density;
        offset_temperature = num_ghosts_temperature;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_internal_energy->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_temperature->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_internal_energy = domain.lower() - ghost_box_internal_energy.lower();
        offset_density = domain.lower() - ghost_box_density.lower();
        offset_temperature = domain.lower() - ghost_box_temperature.lower();
    }
    
    /*
     * Get the pointers to the cell data.
     */
    
    Real* const epsilon = data_internal_energy->getPointer(0);
    const Real* const rho = data_density->getPointer(0);
    const Real* const T = data_temperature->getPointer(0);
    
    const Real& p_inf = *(thermo_properties[1]);
    const Real& c_v = *(thermo_properties[3]);
    
    computeInternalEnergyFromTemperature(
        epsilon,
        rho,
        T,
        p_inf,
        c_v,
        offset_internal_energy,
        offset_density,
        offset_temperature,
        ghostcell_dims_internal_energy,
        ghostcell_dims_density,
        ghostcell_dims_temperature,
        domain_lo,
        domain_dims);
}


/*
 * Compute the specific internal energy from temperature.
 */
void
EquationOfStateStiffenedGas::computeInternalEnergyFromTemperature(
    HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_internal_energy,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_temperature,
    const std::vector<const Real*>& thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_internal_energy);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_temperature);
    
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) == 4);
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_internal_energy = data_internal_energy->getGhostBox();
    hier::IntVector ghostcell_dims_internal_energy = ghost_box_internal_energy.numberCells();
    
    const hier::Box ghost_box_density = data_density->getGhostBox();
    hier::IntVector ghostcell_dims_density = ghost_box_density.numberCells();
    
    const hier::Box ghost_box_temperature = data_temperature->getGhostBox();
    hier::IntVector ghostcell_dims_temperature = ghost_box_temperature.numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_internal_energy(d_dim);
    hier::IntVector offset_density(d_dim);
    hier::IntVector offset_temperature(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_internal_energy = data_internal_energy->getGhostCellWidth();
        const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_internal_energy->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_temperature->getBox().isSpatiallyEqual(interior_box));
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_internal_energy;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_temperature, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_internal_energy = num_ghosts_internal_energy;
        offset_density = num_ghosts_density;
        offset_temperature = num_ghosts_temperature;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_internal_energy->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_temperature->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_internal_energy = domain.lower() - ghost_box_internal_energy.lower();
        offset_density = domain.lower() - ghost_box_density.lower();
        offset_temperature = domain.lower() - ghost_box_temperature.lower();
    }
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(side_normal < d_dim.getValue());
    
    TBOX_ASSERT(data_internal_energy->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_density->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_temperature->getDirectionVector()[side_normal] > 0);
#endif
    
    ghostcell_dims_internal_energy[side_normal]++;
    ghostcell_dims_density[side_normal]++;
    ghostcell_dims_temperature[side_normal]++;
    domain_dims[side_normal]++;
    
    /*
     * Get the pointers to the cell data.
     */
    
    Real* const epsilon = data_internal_energy->getPointer(side_normal, 0);
    const Real* const rho = data_density->getPointer(side_normal, 0);
    const Real* const T = data_temperature->getPointer(side_normal, 0);
    
    const Real& p_inf = *(thermo_properties[1]);
    const Real& c_v = *(thermo_properties[3]);
    
    computeInternalEnergyFromTemperature(
        epsilon,
        rho,
        T,
        p_inf,
        c_v,
        offset_internal_energy,
        offset_density,
        offset_temperature,
        ghostcell_dims_internal_energy,
        ghostcell_dims_density,
        ghostcell_dims_temperature,
        domain_lo,
        domain_dims);
}


/*
 * Compute the specific internal energy from temperature.
 */
void
EquationOfStateStiffenedGas::computeInternalEnergyFromTemperature(
    HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_internal_energy,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_temperature,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_thermo_properties,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_internal_energy);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_temperature);
    TBOX_ASSERT(data_thermo_properties);
    
    TBOX_ASSERT(data_thermo_properties->getDepth() == 4);
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_internal_energy = data_internal_energy->getGhostBox();
    const hier::IntVector ghostcell_dims_internal_energy = ghost_box_internal_energy.numberCells();
    
    const hier::Box ghost_box_density = data_density->getGhostBox();
    const hier::IntVector ghostcell_dims_density = ghost_box_density.numberCells();
    
    const hier::Box ghost_box_temperature = data_temperature->getGhostBox();
    const hier::IntVector ghostcell_dims_temperature = ghost_box_temperature.numberCells();
    
    const hier::Box ghost_box_thermo_properties = data_thermo_properties->getGhostBox();
    const hier::IntVector ghostcell_dims_thermo_properties = ghost_box_thermo_properties.numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_internal_energy(d_dim);
    hier::IntVector offset_density(d_dim);
    hier::IntVector offset_temperature(d_dim);
    hier::IntVector offset_thermo_properties(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_internal_energy = data_internal_energy->getGhostCellWidth();
        const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
        const hier::IntVector num_ghosts_thermo_properties = data_thermo_properties->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_internal_energy->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_temperature->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_thermo_properties->getBox().isSpatiallyEqual(interior_box));
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_internal_energy;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_temperature, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_thermo_properties, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_internal_energy = num_ghosts_internal_energy;
        offset_density = num_ghosts_density;
        offset_temperature = num_ghosts_temperature;
        offset_thermo_properties = num_ghosts_thermo_properties;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_internal_energy->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_temperature->getGhostBox().contains(domain));
        TBOX_ASSERT(data_thermo_properties->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_internal_energy = domain.lower() - ghost_box_internal_energy.lower();
        offset_density = domain.lower() - ghost_box_density.lower();
        offset_temperature = domain.lower() - ghost_box_temperature.lower();
        offset_thermo_properties = domain.lower() - ghost_box_thermo_properties.lower();
    }
    
    /*
     * Get the pointers to the cell data.
     */
    
    Real* const epsilon = data_internal_energy->getPointer(0);
    const Real* const rho = data_density->getPointer(0);
    const Real* const T = data_temperature->getPointer(0);
    const Real* const p_inf = data_thermo_properties->getPointer(1);
    const Real* const c_v = data_thermo_properties->getPointer(3);
    
    computeInternalEnergyFromTemperature(
        epsilon,
        rho,
        T,
        p_inf,
        c_v,
        offset_internal_energy,
        offset_density,
        offset_temperature,
        offset_thermo_properties,
        ghostcell_dims_internal_energy,
        ghostcell_dims_density,
        ghostcell_dims_temperature,
        ghostcell_dims_thermo_properties,
        domain_lo,
        domain_dims);
}


/*
 * Compute the specific internal energy from temperature.
 */
void
EquationOfStateStiffenedGas::computeInternalEnergyFromTemperature(
    HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_internal_energy,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_temperature,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_internal_energy);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_temperature);
    TBOX_ASSERT(data_thermo_properties);
    
    TBOX_ASSERT(data_thermo_properties->getDepth() == 4);
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_internal_energy = data_internal_energy->getGhostBox();
    hier::IntVector ghostcell_dims_internal_energy = ghost_box_internal_energy.numberCells();
    
    const hier::Box ghost_box_density = data_density->getGhostBox();
    hier::IntVector ghostcell_dims_density = ghost_box_density.numberCells();
    
    const hier::Box ghost_box_temperature = data_temperature->getGhostBox();
    hier::IntVector ghostcell_dims_temperature = ghost_box_temperature.numberCells();
    
    const hier::Box ghost_box_thermo_properties = data_thermo_properties->getGhostBox();
    hier::IntVector ghostcell_dims_thermo_properties = ghost_box_thermo_properties.numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_internal_energy(d_dim);
    hier::IntVector offset_density(d_dim);
    hier::IntVector offset_temperature(d_dim);
    hier::IntVector offset_thermo_properties(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_internal_energy = data_internal_energy->getGhostCellWidth();
        const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
        const hier::IntVector num_ghosts_thermo_properties = data_thermo_properties->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_internal_energy->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_temperature->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_thermo_properties->getBox().isSpatiallyEqual(interior_box));
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_internal_energy;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_temperature, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_thermo_properties, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_internal_energy = num_ghosts_internal_energy;
        offset_density = num_ghosts_density;
        offset_temperature = num_ghosts_temperature;
        offset_thermo_properties = num_ghosts_thermo_properties;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_internal_energy->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_temperature->getGhostBox().contains(domain));
        TBOX_ASSERT(data_thermo_properties->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_internal_energy = domain.lower() - ghost_box_internal_energy.lower();
        offset_density = domain.lower() - ghost_box_density.lower();
        offset_temperature = domain.lower() - ghost_box_temperature.lower();
        offset_thermo_properties = domain.lower() - ghost_box_thermo_properties.lower();
    }
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(side_normal < d_dim.getValue());
    
    TBOX_ASSERT(data_internal_energy->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_density->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_temperature->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_thermo_properties->getDirectionVector()[side_normal] > 0);
#endif
    
    ghostcell_dims_internal_energy[side_normal]++;
    ghostcell_dims_density[side_normal]++;
    ghostcell_dims_temperature[side_normal]++;
    ghostcell_dims_thermo_properties[side_normal]++;
    domain_dims[side_normal]++;
    
    /*
     * Get the pointers to the cell data.
     */
    
    Real* const epsilon = data_internal_energy->getPointer(side_normal, 0);
    const Real* const rho = data_density->getPointer(side_normal, 0);
    const Real* const T = data_temperature->getPointer(side_normal, 0);
    const Real* const p_inf = data_thermo_properties->getPointer(side_normal, 1);
    const Real* const c_v = data_thermo_properties->getPointer(side_normal, 3);
    
    computeInternalEnergyFromTemperature(
        epsilon,
        rho,
        T,
        p_inf,
        c_v,
        offset_internal_energy,
        offset_density,
        offset_temperature,
        offset_thermo_properties,
        ghostcell_dims_internal_energy,
        ghostcell_dims_density,
        ghostcell_dims_temperature,
        ghostcell_dims_thermo_properties,
        domain_lo,
        domain_dims);
}


/*
 * Compute the isochoric specific heat capacity.
 */
Real
EquationOfStateStiffenedGas::getIsochoricSpecificHeatCapacity(
    const Real* const density,
    const Real* const pressure,
    const std::vector<const Real*>& thermo_properties) const
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
EquationOfStateStiffenedGas::computeIsochoricSpecificHeatCapacity(
    HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_isochoric_specific_heat_capacity,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
    const std::vector<const Real*>& thermo_properties,
    const hier::Box& domain) const
{
    NULL_USE(data_density);
    NULL_USE(data_pressure);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_isochoric_specific_heat_capacity);
    
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) == 4);
#endif
    
    // Get the dimensions of the ghost cell box.
    const hier::Box ghost_box_isochoric_specific_heat_capacity =
        data_isochoric_specific_heat_capacity->getGhostBox();
    const hier::IntVector ghostcell_dims_isochoric_specific_heat_capacity =
        ghost_box_isochoric_specific_heat_capacity.numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offset.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_isochoric_specific_heat_capacity(d_dim);
    
    if (domain.empty())
    {
        // Get the number of ghost cells.
        const hier::IntVector num_ghosts_isochoric_specific_heat_capacity =
            data_isochoric_specific_heat_capacity->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_isochoric_specific_heat_capacity->getBox();
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_isochoric_specific_heat_capacity);
        
        domain_lo = -num_ghosts_isochoric_specific_heat_capacity;
        domain_dims = ghost_box.numberCells();
        
        offset_isochoric_specific_heat_capacity = num_ghosts_isochoric_specific_heat_capacity;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_isochoric_specific_heat_capacity->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_isochoric_specific_heat_capacity = domain.lower() - ghost_box_isochoric_specific_heat_capacity.lower();
    }
    
    /*
     * Get the pointer to the cell data.
     */
    
    Real* const c_v = data_isochoric_specific_heat_capacity->getPointer(0);
    
    const Real& c_v_src = *(thermo_properties[3]);
    
    computeIsochoricSpecificHeatCapacity(
        c_v,
        c_v_src,
        offset_isochoric_specific_heat_capacity,
        ghostcell_dims_isochoric_specific_heat_capacity,
        domain_lo,
        domain_dims);
}


/*
 * Compute the isochoric specific heat capacity.
 */
void
EquationOfStateStiffenedGas::computeIsochoricSpecificHeatCapacity(
    HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_isochoric_specific_heat_capacity,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
    const std::vector<const Real*>& thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
    NULL_USE(data_density);
    NULL_USE(data_pressure);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_isochoric_specific_heat_capacity);
    
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) == 4);
#endif
    
    // Get the dimensions of the ghost cell box.
    const hier::Box ghost_box_isochoric_specific_heat_capacity =
        data_isochoric_specific_heat_capacity->getGhostBox();
    hier::IntVector ghostcell_dims_isochoric_specific_heat_capacity =
        ghost_box_isochoric_specific_heat_capacity.numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offset.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_isochoric_specific_heat_capacity(d_dim);
    
    if (domain.empty())
    {
        // Get the number of ghost cells.
        const hier::IntVector num_ghosts_isochoric_specific_heat_capacity =
            data_isochoric_specific_heat_capacity->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_isochoric_specific_heat_capacity->getBox();
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_isochoric_specific_heat_capacity);
        
        domain_lo = -num_ghosts_isochoric_specific_heat_capacity;
        domain_dims = ghost_box.numberCells();
        
        offset_isochoric_specific_heat_capacity = num_ghosts_isochoric_specific_heat_capacity;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_isochoric_specific_heat_capacity->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_isochoric_specific_heat_capacity = domain.lower() - ghost_box_isochoric_specific_heat_capacity.lower();
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
    
    Real* const c_v = data_isochoric_specific_heat_capacity->getPointer(side_normal, 0);
    
    const Real& c_v_src = *(thermo_properties[3]);
    
    computeIsochoricSpecificHeatCapacity(
        c_v,
        c_v_src,
        offset_isochoric_specific_heat_capacity,
        ghostcell_dims_isochoric_specific_heat_capacity,
        domain_lo,
        domain_dims);
}


/*
 * Compute the isochoric specific heat capacity.
 */
void
EquationOfStateStiffenedGas::computeIsochoricSpecificHeatCapacity(
    HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_isochoric_specific_heat_capacity,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_thermo_properties,
    const hier::Box& domain) const
{
    NULL_USE(data_density);
    NULL_USE(data_pressure);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_isochoric_specific_heat_capacity);
    TBOX_ASSERT(data_thermo_properties);
    
    TBOX_ASSERT(data_thermo_properties->getDepth() == 4);
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_isochoric_specific_heat_capacity =
        data_isochoric_specific_heat_capacity->getGhostBox();
    const hier::IntVector ghostcell_dims_isochoric_specific_heat_capacity =
        ghost_box_isochoric_specific_heat_capacity.numberCells();
    
    const hier::Box ghost_box_thermo_properties = data_thermo_properties->getGhostBox();
    const hier::IntVector ghostcell_dims_thermo_properties = ghost_box_thermo_properties.numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_isochoric_specific_heat_capacity(d_dim);
    hier::IntVector offset_thermo_properties(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_isochoric_specific_heat_capacity =
            data_isochoric_specific_heat_capacity->getGhostCellWidth();
        const hier::IntVector num_ghosts_thermo_properties = data_thermo_properties->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_isochoric_specific_heat_capacity->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_thermo_properties->getBox().isSpatiallyEqual(interior_box));
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_isochoric_specific_heat_capacity;
        num_ghosts_min = hier::IntVector::min(num_ghosts_thermo_properties, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_isochoric_specific_heat_capacity = num_ghosts_isochoric_specific_heat_capacity;
        offset_thermo_properties = num_ghosts_thermo_properties;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_isochoric_specific_heat_capacity->getGhostBox().contains(domain));
        TBOX_ASSERT(data_thermo_properties->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_isochoric_specific_heat_capacity = domain.lower() - ghost_box_isochoric_specific_heat_capacity.lower();
        offset_thermo_properties = domain.lower() - ghost_box_thermo_properties.lower();
    }
    
    /*
     * Get the pointers to the cell data.
     */
    
    Real* const c_v = data_isochoric_specific_heat_capacity->getPointer(0);
    const Real* const c_v_src = data_thermo_properties->getPointer(3);
    
    computeIsochoricSpecificHeatCapacity(
        c_v,
        c_v_src,
        offset_isochoric_specific_heat_capacity,
        offset_thermo_properties,
        ghostcell_dims_isochoric_specific_heat_capacity,
        ghostcell_dims_thermo_properties,
        domain_lo,
        domain_dims);
}


/*
 * Compute the isochoric specific heat capacity.
 */
void
EquationOfStateStiffenedGas::computeIsochoricSpecificHeatCapacity(
    HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_isochoric_specific_heat_capacity,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_thermo_properties,
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
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_isochoric_specific_heat_capacity =
        data_isochoric_specific_heat_capacity->getGhostBox();
    hier::IntVector ghostcell_dims_isochoric_specific_heat_capacity =
        ghost_box_isochoric_specific_heat_capacity.numberCells();
    
    const hier::Box ghost_box_thermo_properties = data_thermo_properties->getGhostBox();
    hier::IntVector ghostcell_dims_thermo_properties = ghost_box_thermo_properties.numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_isochoric_specific_heat_capacity(d_dim);
    hier::IntVector offset_thermo_properties(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_isochoric_specific_heat_capacity =
            data_isochoric_specific_heat_capacity->getGhostCellWidth();
        const hier::IntVector num_ghosts_thermo_properties = data_thermo_properties->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_isochoric_specific_heat_capacity->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_thermo_properties->getBox().isSpatiallyEqual(interior_box));
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_isochoric_specific_heat_capacity;
        num_ghosts_min = hier::IntVector::min(num_ghosts_thermo_properties, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_isochoric_specific_heat_capacity = num_ghosts_isochoric_specific_heat_capacity;
        offset_thermo_properties = num_ghosts_thermo_properties;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_isochoric_specific_heat_capacity->getGhostBox().contains(domain));
        TBOX_ASSERT(data_thermo_properties->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_isochoric_specific_heat_capacity = domain.lower() - ghost_box_isochoric_specific_heat_capacity.lower();
        offset_thermo_properties = domain.lower() - ghost_box_thermo_properties.lower();
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
    
    Real* const c_v = data_isochoric_specific_heat_capacity->getPointer(side_normal, 0);
    const Real* const c_v_src = data_thermo_properties->getPointer(side_normal, 3);
    
    computeIsochoricSpecificHeatCapacity(
        c_v,
        c_v_src,
        offset_isochoric_specific_heat_capacity,
        offset_thermo_properties,
        ghostcell_dims_isochoric_specific_heat_capacity,
        ghostcell_dims_thermo_properties,
        domain_lo,
        domain_dims);
}


/*
 * Compute the isobaric specific heat capacity.
 */
Real
EquationOfStateStiffenedGas::getIsobaricSpecificHeatCapacity(
    const Real* const density,
    const Real* const pressure,
    const std::vector<const Real*>& thermo_properties) const
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
EquationOfStateStiffenedGas::computeIsobaricSpecificHeatCapacity(
    HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_isobaric_specific_heat_capacity,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
    const std::vector<const Real*>& thermo_properties,
    const hier::Box& domain) const
{
    NULL_USE(data_density);
    NULL_USE(data_pressure);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_isobaric_specific_heat_capacity);
    
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 3);
#endif
    
    // Get the dimensions of the ghost cell box.
    const hier::Box ghost_box_isobaric_specific_heat_capacity =
        data_isobaric_specific_heat_capacity->getGhostBox();
    const hier::IntVector ghostcell_dims_isobaric_specific_heat_capacity =
        ghost_box_isobaric_specific_heat_capacity.numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offset.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_isobaric_specific_heat_capacity(d_dim);
    
    if (domain.empty())
    {
        // Get the number of ghost cells.
        const hier::IntVector num_ghosts_isobaric_specific_heat_capacity =
            data_isobaric_specific_heat_capacity->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_isobaric_specific_heat_capacity->getBox();
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_isobaric_specific_heat_capacity);
        
        domain_lo = -num_ghosts_isobaric_specific_heat_capacity;
        domain_dims = ghost_box.numberCells();
        
        offset_isobaric_specific_heat_capacity = num_ghosts_isobaric_specific_heat_capacity;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_isobaric_specific_heat_capacity->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_isobaric_specific_heat_capacity = domain.lower() - ghost_box_isobaric_specific_heat_capacity.lower();
    }
    
    /*
     * Get the pointer to the cell data.
     */
    
    Real* const c_p = data_isobaric_specific_heat_capacity->getPointer(0);
    
    const Real& c_p_src = *(thermo_properties[2]);
    
    computeIsobaricSpecificHeatCapacity(
        c_p,
        c_p_src,
        offset_isobaric_specific_heat_capacity,
        ghostcell_dims_isobaric_specific_heat_capacity,
        domain_lo,
        domain_dims);
}


/*
 * Compute the isobaric specific heat capacity.
 */
void
EquationOfStateStiffenedGas::computeIsobaricSpecificHeatCapacity(
    HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_isobaric_specific_heat_capacity,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
    const std::vector<const Real*>& thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
    NULL_USE(data_density);
    NULL_USE(data_pressure);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_isobaric_specific_heat_capacity);
    
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 3);
#endif
    
    // Get the dimensions of the ghost cell box.
    const hier::Box ghost_box_isobaric_specific_heat_capacity =
        data_isobaric_specific_heat_capacity->getGhostBox();
    hier::IntVector ghostcell_dims_isobaric_specific_heat_capacity =
        ghost_box_isobaric_specific_heat_capacity.numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offset.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_isobaric_specific_heat_capacity(d_dim);
    
    if (domain.empty())
    {
        // Get the number of ghost cells.
        const hier::IntVector num_ghosts_isobaric_specific_heat_capacity =
            data_isobaric_specific_heat_capacity->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_isobaric_specific_heat_capacity->getBox();
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_isobaric_specific_heat_capacity);
        
        domain_lo = -num_ghosts_isobaric_specific_heat_capacity;
        domain_dims = ghost_box.numberCells();
        
        offset_isobaric_specific_heat_capacity = num_ghosts_isobaric_specific_heat_capacity;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_isobaric_specific_heat_capacity->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_isobaric_specific_heat_capacity = domain.lower() - ghost_box_isobaric_specific_heat_capacity.lower();
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
    
    Real* const c_p = data_isobaric_specific_heat_capacity->getPointer(side_normal, 0);
    
    const Real& c_p_src = *(thermo_properties[2]);
    
    computeIsobaricSpecificHeatCapacity(
        c_p,
        c_p_src,
        offset_isobaric_specific_heat_capacity,
        ghostcell_dims_isobaric_specific_heat_capacity,
        domain_lo,
        domain_dims);
}


/*
 * Compute the isobaric specific heat capacity.
 */
void
EquationOfStateStiffenedGas::computeIsobaricSpecificHeatCapacity(
    HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_isobaric_specific_heat_capacity,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_thermo_properties,
    const hier::Box& domain) const
{
    NULL_USE(data_density);
    NULL_USE(data_pressure);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_isobaric_specific_heat_capacity);
    TBOX_ASSERT(data_thermo_properties);
    
    TBOX_ASSERT(data_thermo_properties->getDepth() >= 3);
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_isobaric_specific_heat_capacity =
        data_isobaric_specific_heat_capacity->getGhostBox();
    const hier::IntVector ghostcell_dims_isobaric_specific_heat_capacity =
        ghost_box_isobaric_specific_heat_capacity.numberCells();
    
    const hier::Box ghost_box_thermo_properties = data_thermo_properties->getGhostBox();
    const hier::IntVector ghostcell_dims_thermo_properties = ghost_box_thermo_properties.numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_isobaric_specific_heat_capacity(d_dim);
    hier::IntVector offset_thermo_properties(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_isobaric_specific_heat_capacity =
            data_isobaric_specific_heat_capacity->getGhostCellWidth();
        const hier::IntVector num_ghosts_thermo_properties = data_thermo_properties->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_isobaric_specific_heat_capacity->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_thermo_properties->getBox().isSpatiallyEqual(interior_box));
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_isobaric_specific_heat_capacity;
        num_ghosts_min = hier::IntVector::min(num_ghosts_thermo_properties, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_isobaric_specific_heat_capacity = num_ghosts_isobaric_specific_heat_capacity;
        offset_thermo_properties = num_ghosts_thermo_properties;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_isobaric_specific_heat_capacity->getGhostBox().contains(domain));
        TBOX_ASSERT(data_thermo_properties->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_isobaric_specific_heat_capacity = domain.lower() - ghost_box_isobaric_specific_heat_capacity.lower();
        offset_thermo_properties = domain.lower() - ghost_box_thermo_properties.lower();
    }
    
    /*
     * Get the pointers to the cell data.
     */
    
    Real* const c_p = data_isobaric_specific_heat_capacity->getPointer(0);
    const Real* const c_p_src = data_thermo_properties->getPointer(2);
    
    computeIsobaricSpecificHeatCapacity(
        c_p,
        c_p_src,
        offset_isobaric_specific_heat_capacity,
        offset_thermo_properties,
        ghostcell_dims_isobaric_specific_heat_capacity,
        ghostcell_dims_thermo_properties,
        domain_lo,
        domain_dims);
}


/*
 * Compute the isobaric specific heat capacity.
 */
void
EquationOfStateStiffenedGas::computeIsobaricSpecificHeatCapacity(
    HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_isobaric_specific_heat_capacity,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_thermo_properties,
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
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_isobaric_specific_heat_capacity =
        data_isobaric_specific_heat_capacity->getGhostBox();
    hier::IntVector ghostcell_dims_isobaric_specific_heat_capacity =
        ghost_box_isobaric_specific_heat_capacity.numberCells();
    
    const hier::Box ghost_box_thermo_properties = data_thermo_properties->getGhostBox();
    hier::IntVector ghostcell_dims_thermo_properties = ghost_box_thermo_properties.numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_isobaric_specific_heat_capacity(d_dim);
    hier::IntVector offset_thermo_properties(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_isobaric_specific_heat_capacity =
            data_isobaric_specific_heat_capacity->getGhostCellWidth();
        const hier::IntVector num_ghosts_thermo_properties = data_thermo_properties->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_isobaric_specific_heat_capacity->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_thermo_properties->getBox().isSpatiallyEqual(interior_box));
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_isobaric_specific_heat_capacity;
        num_ghosts_min = hier::IntVector::min(num_ghosts_thermo_properties, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_isobaric_specific_heat_capacity = num_ghosts_isobaric_specific_heat_capacity;
        offset_thermo_properties = num_ghosts_thermo_properties;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_isobaric_specific_heat_capacity->getGhostBox().contains(domain));
        TBOX_ASSERT(data_thermo_properties->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_isobaric_specific_heat_capacity = domain.lower() - ghost_box_isobaric_specific_heat_capacity.lower();
        offset_thermo_properties = domain.lower() - ghost_box_thermo_properties.lower();
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
    
    Real* const c_p = data_isobaric_specific_heat_capacity->getPointer(side_normal, 0);
    const Real* const c_p_src = data_thermo_properties->getPointer(side_normal, 2);
    
    computeIsobaricSpecificHeatCapacity(
        c_p,
        c_p_src,
        offset_isobaric_specific_heat_capacity,
        offset_thermo_properties,
        ghostcell_dims_isobaric_specific_heat_capacity,
        ghostcell_dims_thermo_properties,
        domain_lo,
        domain_dims);
}


/*
 * Compute the Gruneisen parameter (partial derivative of pressure w.r.t. specific internal energy under
 * constant density divided by density).
 */
Real
EquationOfStateStiffenedGas::getGruneisenParameter(
    const Real* const density,
    const Real* const pressure,
    const std::vector<const Real*>& thermo_properties) const
{
    NULL_USE(density);
    NULL_USE(pressure);
    
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 1);
#endif
    
    const Real& gamma = *(thermo_properties[0]);
    
    return (gamma - Real(1));
}


/*
 * Compute the Gruneisen parameter (partial derivative of pressure w.r.t. specific internal energy under
 * constant density divided by density).
 */
void
EquationOfStateStiffenedGas::computeGruneisenParameter(
    HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_gruneisen_parameter,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
    const std::vector<const Real*>& thermo_properties,
    const hier::Box& domain) const
{
    NULL_USE(data_density);
    NULL_USE(data_pressure);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_gruneisen_parameter);
    
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 1);
#endif
    
    // Get the dimensions of the ghost cell box.
    const hier::Box ghost_box_gruneisen_parameter = data_gruneisen_parameter->getGhostBox();
    const hier::IntVector ghostcell_dims_gruneisen_parameter = ghost_box_gruneisen_parameter.numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offset.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_gruneisen_parameter(d_dim);
    
    if (domain.empty())
    {
        // Get the number of ghost cells.
        const hier::IntVector num_ghosts_gruneisen_parameter = data_gruneisen_parameter->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_gruneisen_parameter->getBox();
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_gruneisen_parameter);
        
        domain_lo = -num_ghosts_gruneisen_parameter;
        domain_dims = ghost_box.numberCells();
        
        offset_gruneisen_parameter = num_ghosts_gruneisen_parameter;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_gruneisen_parameter->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_gruneisen_parameter = domain.lower() - ghost_box_gruneisen_parameter.lower();
    }
    
    /*
     * Get the pointer to the cell data.
     */
    
    Real* const Gamma = data_gruneisen_parameter->getPointer(0);
    
    const Real& gamma = *(thermo_properties[0]);
    
    computeGruneisenParameter(
        Gamma,
        gamma,
        offset_gruneisen_parameter,
        ghostcell_dims_gruneisen_parameter,
        domain_lo,
        domain_dims);
}


/*
 * Compute the Gruneisen parameter (partial derivative of pressure w.r.t. specific internal energy under
 * constant density divided by density).
 */
void
EquationOfStateStiffenedGas::computeGruneisenParameter(
    HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_gruneisen_parameter,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
    const std::vector<const Real*>& thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
    NULL_USE(data_density);
    NULL_USE(data_pressure);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_gruneisen_parameter);
    
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 1);
#endif
    
    // Get the dimensions of the ghost cell box.
    const hier::Box ghost_box_gruneisen_parameter = data_gruneisen_parameter->getGhostBox();
    hier::IntVector ghostcell_dims_gruneisen_parameter = ghost_box_gruneisen_parameter.numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offset.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_gruneisen_parameter(d_dim);
    
    if (domain.empty())
    {
        // Get the number of ghost cells.
        const hier::IntVector num_ghosts_gruneisen_parameter = data_gruneisen_parameter->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_gruneisen_parameter->getBox();
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_gruneisen_parameter);
        
        domain_lo = -num_ghosts_gruneisen_parameter;
        domain_dims = ghost_box.numberCells();
        
        offset_gruneisen_parameter = num_ghosts_gruneisen_parameter;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_gruneisen_parameter->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_gruneisen_parameter = domain.lower() - ghost_box_gruneisen_parameter.lower();
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
    
    Real* const Gamma = data_gruneisen_parameter->getPointer(side_normal, 0);
    
    const Real& gamma = *(thermo_properties[0]);
    
    computeGruneisenParameter(
        Gamma,
        gamma,
        offset_gruneisen_parameter,
        ghostcell_dims_gruneisen_parameter,
        domain_lo,
        domain_dims);
}


/*
 * Compute the Gruneisen parameter (partial derivative of pressure w.r.t. specific internal energy under
 * constant density divided by density).
 */
void
EquationOfStateStiffenedGas::computeGruneisenParameter(
    HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_gruneisen_parameter,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_thermo_properties,
    const hier::Box& domain) const
{
    NULL_USE(data_density);
    NULL_USE(data_pressure);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_gruneisen_parameter);
    TBOX_ASSERT(data_thermo_properties);
    
    TBOX_ASSERT(data_thermo_properties->getDepth() >= 1);
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_gruneisen_parameter = data_gruneisen_parameter->getGhostBox();
    const hier::IntVector ghostcell_dims_gruneisen_parameter = ghost_box_gruneisen_parameter.numberCells();
    
    const hier::Box ghost_box_thermo_properties = data_thermo_properties->getGhostBox();
    const hier::IntVector ghostcell_dims_thermo_properties = ghost_box_thermo_properties.numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_gruneisen_parameter(d_dim);
    hier::IntVector offset_thermo_properties(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_gruneisen_parameter = data_gruneisen_parameter->getGhostCellWidth();
        const hier::IntVector num_ghosts_thermo_properties = data_thermo_properties->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_gruneisen_parameter->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_thermo_properties->getBox().isSpatiallyEqual(interior_box));
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_gruneisen_parameter;
        num_ghosts_min = hier::IntVector::min(num_ghosts_thermo_properties, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_gruneisen_parameter = num_ghosts_gruneisen_parameter;
        offset_thermo_properties = num_ghosts_thermo_properties;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_gruneisen_parameter->getGhostBox().contains(domain));
        TBOX_ASSERT(data_thermo_properties->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_gruneisen_parameter = domain.lower() - ghost_box_gruneisen_parameter.lower();
        offset_thermo_properties = domain.lower() - ghost_box_thermo_properties.lower();
    }
    
    /*
     * Get the pointers to the cell data.
     */
    
    Real* const Gamma = data_gruneisen_parameter->getPointer(0);
    const Real* const gamma = data_thermo_properties->getPointer(0);
    
    computeGruneisenParameter(
        Gamma,
        gamma,
        offset_gruneisen_parameter,
        offset_thermo_properties,
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
EquationOfStateStiffenedGas::computeGruneisenParameter(
    HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_gruneisen_parameter,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_thermo_properties,
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
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_gruneisen_parameter = data_gruneisen_parameter->getGhostBox();
    hier::IntVector ghostcell_dims_gruneisen_parameter = ghost_box_gruneisen_parameter.numberCells();
    
    const hier::Box ghost_box_thermo_properties = data_thermo_properties->getGhostBox();
    hier::IntVector ghostcell_dims_thermo_properties = ghost_box_thermo_properties.numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_gruneisen_parameter(d_dim);
    hier::IntVector offset_thermo_properties(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_gruneisen_parameter = data_gruneisen_parameter->getGhostCellWidth();
        const hier::IntVector num_ghosts_thermo_properties = data_thermo_properties->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_gruneisen_parameter->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_thermo_properties->getBox().isSpatiallyEqual(interior_box));
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_gruneisen_parameter;
        num_ghosts_min = hier::IntVector::min(num_ghosts_thermo_properties, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_gruneisen_parameter = num_ghosts_gruneisen_parameter;
        offset_thermo_properties = num_ghosts_thermo_properties;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_gruneisen_parameter->getGhostBox().contains(domain));
        TBOX_ASSERT(data_thermo_properties->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_gruneisen_parameter = domain.lower() - ghost_box_gruneisen_parameter.lower();
        offset_thermo_properties = domain.lower() - ghost_box_thermo_properties.lower();
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
    
    Real* const Gamma = data_gruneisen_parameter->getPointer(side_normal, 0);
    const Real* const gamma = data_thermo_properties->getPointer(side_normal, 0);
    
    computeGruneisenParameter(
        Gamma,
        gamma,
        offset_gruneisen_parameter,
        offset_thermo_properties,
        ghostcell_dims_gruneisen_parameter,
        ghostcell_dims_thermo_properties,
        domain_lo,
        domain_dims);
}


/*
 * Compute the partial derivative of pressure w.r.t. density under constant specific internal energy.
 */
Real
EquationOfStateStiffenedGas::getPressureDerivativeWithDensity(
    const Real* const density,
    const Real* const pressure,
    const std::vector<const Real*>& thermo_properties) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) == 2);
#endif
    
    const Real& gamma = *(thermo_properties[0]);
    const Real& p_inf = *(thermo_properties[1]);
    
    const Real& rho = *density;
    const Real& p = *pressure;
    
    return (p + gamma*p_inf)/rho;
}


/*
 * Compute the partial derivative of pressure w.r.t. density under constant specific internal energy.
 */
void
EquationOfStateStiffenedGas::computePressureDerivativeWithDensity(
    HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_partial_pressure_partial_density,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
    const std::vector<const Real*>& thermo_properties,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_partial_pressure_partial_density);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 2);
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_partial_pressure_partial_density =
        data_partial_pressure_partial_density->getGhostBox();
    const hier::IntVector ghostcell_dims_partial_pressure_partial_density =
        ghost_box_partial_pressure_partial_density.numberCells();
    
    const hier::Box ghost_box_density = data_density->getGhostBox();
    const hier::IntVector ghostcell_dims_density = ghost_box_density.numberCells();
    
    const hier::Box ghost_box_pressure = data_pressure->getGhostBox();
    const hier::IntVector ghostcell_dims_pressure = ghost_box_pressure.numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_partial_pressure_partial_density(d_dim);
    hier::IntVector offset_density(d_dim);
    hier::IntVector offset_pressure(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_partial_pressure_partial_density =
            data_partial_pressure_partial_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_partial_pressure_partial_density->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_pressure->getBox().isSpatiallyEqual(interior_box));
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_partial_pressure_partial_density;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_partial_pressure_partial_density = num_ghosts_partial_pressure_partial_density;
        offset_density = num_ghosts_density;
        offset_pressure = num_ghosts_pressure;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_partial_pressure_partial_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_partial_pressure_partial_density = domain.lower() - ghost_box_partial_pressure_partial_density.lower();
        offset_density = domain.lower() - ghost_box_density.lower();
        offset_pressure = domain.lower() - ghost_box_pressure.lower();
    }
    
    /*
     * Get the pointers to the cell data.
     */
    
    Real* const Psi = data_partial_pressure_partial_density->getPointer(0);
    const Real* const rho = data_density->getPointer(0);
    const Real* const p = data_pressure->getPointer(0);
    
    const Real& gamma = *(thermo_properties[0]);
    const Real& p_inf = *(thermo_properties[1]);
    
    computePressureDerivativeWithDensity(
        Psi,
        rho,
        p,
        gamma,
        p_inf,
        offset_partial_pressure_partial_density,
        offset_density,
        offset_pressure,
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
EquationOfStateStiffenedGas::computePressureDerivativeWithDensity(
    HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_partial_pressure_partial_density,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
    const std::vector<const Real*>& thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_partial_pressure_partial_density);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) >= 2);
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_partial_pressure_partial_density =
        data_partial_pressure_partial_density->getGhostBox();
    hier::IntVector ghostcell_dims_partial_pressure_partial_density =
        ghost_box_partial_pressure_partial_density.numberCells();
    
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
    
    hier::IntVector offset_partial_pressure_partial_density(d_dim);
    hier::IntVector offset_density(d_dim);
    hier::IntVector offset_pressure(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_partial_pressure_partial_density =
            data_partial_pressure_partial_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_partial_pressure_partial_density->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_pressure->getBox().isSpatiallyEqual(interior_box));
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_partial_pressure_partial_density;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_partial_pressure_partial_density = num_ghosts_partial_pressure_partial_density;
        offset_density = num_ghosts_density;
        offset_pressure = num_ghosts_pressure;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_partial_pressure_partial_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_partial_pressure_partial_density = domain.lower() - ghost_box_partial_pressure_partial_density.lower();
        offset_density = domain.lower() - ghost_box_density.lower();
        offset_pressure = domain.lower() - ghost_box_pressure.lower();
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
    
    Real* const Psi = data_partial_pressure_partial_density->getPointer(side_normal, 0);
    const Real* const rho = data_density->getPointer(side_normal, 0);
    const Real* const p = data_pressure->getPointer(side_normal, 0);
    
    const Real& gamma = *(thermo_properties[0]);
    const Real& p_inf = *(thermo_properties[1]);
    
    computePressureDerivativeWithDensity(
        Psi,
        rho,
        p,
        gamma,
        p_inf,
        offset_partial_pressure_partial_density,
        offset_density,
        offset_pressure,
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
EquationOfStateStiffenedGas::computePressureDerivativeWithDensity(
    HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_partial_pressure_partial_density,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_thermo_properties,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_partial_pressure_partial_density);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_thermo_properties);
    
    TBOX_ASSERT(data_thermo_properties->getDepth() >= 2);
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_partial_pressure_partial_density =
        data_partial_pressure_partial_density->getGhostBox();
    const hier::IntVector ghostcell_dims_partial_pressure_partial_density =
        ghost_box_partial_pressure_partial_density.numberCells();
    
    const hier::Box ghost_box_density = data_density->getGhostBox();
    const hier::IntVector ghostcell_dims_density = ghost_box_density.numberCells();
    
    const hier::Box ghost_box_pressure = data_pressure->getGhostBox();
    const hier::IntVector ghostcell_dims_pressure = ghost_box_pressure.numberCells();
    
    const hier::Box ghost_box_thermo_properties = data_thermo_properties->getGhostBox();
    const hier::IntVector ghostcell_dims_thermo_properties = ghost_box_thermo_properties.numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_partial_pressure_partial_density(d_dim);
    hier::IntVector offset_density(d_dim);
    hier::IntVector offset_pressure(d_dim);
    hier::IntVector offset_thermo_properties(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_partial_pressure_partial_density =
            data_partial_pressure_partial_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        const hier::IntVector num_ghosts_thermo_properties = data_thermo_properties->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_partial_pressure_partial_density->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_pressure->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_thermo_properties->getBox().isSpatiallyEqual(interior_box));
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_partial_pressure_partial_density;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_thermo_properties, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_partial_pressure_partial_density = num_ghosts_partial_pressure_partial_density;
        offset_density = num_ghosts_density;
        offset_pressure = num_ghosts_pressure;
        offset_thermo_properties = num_ghosts_thermo_properties;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_partial_pressure_partial_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_thermo_properties->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_partial_pressure_partial_density = domain.lower() - ghost_box_partial_pressure_partial_density.lower();
        offset_density = domain.lower() - ghost_box_density.lower();
        offset_pressure = domain.lower() - ghost_box_pressure.lower();
        offset_thermo_properties = domain.lower() - ghost_box_thermo_properties.lower();
    }
    
    /*
     * Get the pointers to the cell data.
     */
    
    Real* const Psi = data_partial_pressure_partial_density->getPointer(0);
    const Real* const rho = data_density->getPointer(0);
    const Real* const p = data_pressure->getPointer(0);
    const Real* const gamma = data_thermo_properties->getPointer(0);
    const Real* const p_inf = data_thermo_properties->getPointer(1);
    
    computePressureDerivativeWithDensity(
        Psi,
        rho,
        p,
        gamma,
        p_inf,
        offset_partial_pressure_partial_density,
        offset_density,
        offset_pressure,
        offset_thermo_properties,
        ghostcell_dims_partial_pressure_partial_density,
        ghostcell_dims_density,
        ghostcell_dims_pressure,
        ghostcell_dims_thermo_properties,
        domain_lo,
        domain_dims);
}


/*
 * Compute the partial derivative of pressure w.r.t. density under constant specific internal energy.
 */
void
EquationOfStateStiffenedGas::computePressureDerivativeWithDensity(
    HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_partial_pressure_partial_density,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_partial_pressure_partial_density);
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_thermo_properties);
    
    TBOX_ASSERT(data_thermo_properties->getDepth() >= 2);
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_partial_pressure_partial_density =
        data_partial_pressure_partial_density->getGhostBox();
    hier::IntVector ghostcell_dims_partial_pressure_partial_density =
        ghost_box_partial_pressure_partial_density.numberCells();
    
    const hier::Box ghost_box_density = data_density->getGhostBox();
    hier::IntVector ghostcell_dims_density = ghost_box_density.numberCells();
    
    const hier::Box ghost_box_pressure = data_pressure->getGhostBox();
    hier::IntVector ghostcell_dims_pressure = ghost_box_pressure.numberCells();
    
    const hier::Box ghost_box_thermo_properties = data_thermo_properties->getGhostBox();
    hier::IntVector ghostcell_dims_thermo_properties = ghost_box_thermo_properties.numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_partial_pressure_partial_density(d_dim);
    hier::IntVector offset_density(d_dim);
    hier::IntVector offset_pressure(d_dim);
    hier::IntVector offset_thermo_properties(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_partial_pressure_partial_density =
            data_partial_pressure_partial_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        const hier::IntVector num_ghosts_thermo_properties = data_thermo_properties->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_partial_pressure_partial_density->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_pressure->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_thermo_properties->getBox().isSpatiallyEqual(interior_box));
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_partial_pressure_partial_density;
        num_ghosts_min = hier::IntVector::min(num_ghosts_density, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_thermo_properties, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_partial_pressure_partial_density = num_ghosts_partial_pressure_partial_density;
        offset_density = num_ghosts_density;
        offset_pressure = num_ghosts_pressure;
        offset_thermo_properties = num_ghosts_thermo_properties;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_partial_pressure_partial_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_thermo_properties->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_partial_pressure_partial_density = domain.lower() - ghost_box_partial_pressure_partial_density.lower();
        offset_density = domain.lower() - ghost_box_density.lower();
        offset_pressure = domain.lower() - ghost_box_pressure.lower();
        offset_thermo_properties = domain.lower() - ghost_box_thermo_properties.lower();
    }
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(side_normal < d_dim.getValue());
    
    TBOX_ASSERT(data_partial_pressure_partial_density->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_density->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_pressure->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_thermo_properties->getDirectionVector()[side_normal] > 0);
#endif
    
    ghostcell_dims_partial_pressure_partial_density[side_normal]++;
    ghostcell_dims_density[side_normal]++;
    ghostcell_dims_pressure[side_normal]++;
    ghostcell_dims_thermo_properties[side_normal]++;
    domain_dims[side_normal]++;
    
    /*
     * Get the pointers to the cell data.
     */
    
    Real* const Psi = data_partial_pressure_partial_density->getPointer(side_normal, 0);
    const Real* const rho = data_density->getPointer(side_normal, 0);
    const Real* const p = data_pressure->getPointer(side_normal, 0);
    const Real* const gamma = data_thermo_properties->getPointer(0);
    const Real* const p_inf = data_thermo_properties->getPointer(1);
    
    computePressureDerivativeWithDensity(
        Psi,
        rho,
        p,
        gamma,
        p_inf,
        offset_partial_pressure_partial_density,
        offset_density,
        offset_pressure,
        offset_thermo_properties,
        ghostcell_dims_partial_pressure_partial_density,
        ghostcell_dims_density,
        ghostcell_dims_pressure,
        ghostcell_dims_thermo_properties,
        domain_lo,
        domain_dims);
}


/*
 * Compute the density.
 */
Real
EquationOfStateStiffenedGas::getDensity(
    const Real* const pressure,
    const Real* const temperature,
    const std::vector<const Real*>& thermo_properties) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) == 4);
#endif
    
    const Real& gamma = *(thermo_properties[0]);
    const Real& p_inf = *(thermo_properties[1]);
    const Real& c_v = *(thermo_properties[3]);
    
    const Real& p = *pressure;
    const Real& T = *temperature;
    
    return (p + p_inf)/((gamma - Real(1))*c_v*T); // Return rho.
}


/*
 * Compute the density.
 */
void
EquationOfStateStiffenedGas::computeDensity(
    HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_temperature,
    const std::vector<const Real*>& thermo_properties,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_temperature);
    
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) == 4);
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_density = data_density->getGhostBox();
    const hier::IntVector ghostcell_dims_density = ghost_box_density.numberCells();
    
    const hier::Box ghost_box_pressure = data_pressure->getGhostBox();
    const hier::IntVector ghostcell_dims_pressure = ghost_box_pressure.numberCells();
    
    const hier::Box ghost_box_temperature = data_temperature->getGhostBox();
    const hier::IntVector ghostcell_dims_temperature = ghost_box_temperature.numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_density(d_dim);
    hier::IntVector offset_pressure(d_dim);
    hier::IntVector offset_temperature(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_density->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_pressure->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_temperature->getBox().isSpatiallyEqual(interior_box));
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_density;
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_temperature, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_density = num_ghosts_density;
        offset_pressure = num_ghosts_pressure;
        offset_temperature = num_ghosts_temperature;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_temperature->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_density = domain.lower() - ghost_box_density.lower();
        offset_pressure = domain.lower() - ghost_box_pressure.lower();
        offset_temperature = domain.lower() - ghost_box_temperature.lower();
    }
    
    /*
     * Get the pointers to the cell data.
     */
    
    Real* const rho = data_density->getPointer(0);
    const Real* const p = data_pressure->getPointer(0);
    const Real* const T = data_temperature->getPointer(0);
    
    const Real& gamma = *(thermo_properties[0]);
    const Real& p_inf = *(thermo_properties[1]);
    const Real& c_v = *(thermo_properties[3]);
    
    computeDensity(
        rho,
        p,
        T,
        gamma,
        p_inf,
        c_v,
        offset_density,
        offset_pressure,
        offset_temperature,
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
EquationOfStateStiffenedGas::computeDensity(
    HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_temperature,
    const std::vector<const Real*>& thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_temperature);
    
    TBOX_ASSERT(static_cast<int>(thermo_properties.size()) == 4);
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_density = data_density->getGhostBox();
    hier::IntVector ghostcell_dims_density = ghost_box_density.numberCells();
    
    const hier::Box ghost_box_pressure = data_pressure->getGhostBox();
    hier::IntVector ghostcell_dims_pressure = ghost_box_pressure.numberCells();
    
    const hier::Box ghost_box_temperature = data_temperature->getGhostBox();
    hier::IntVector ghostcell_dims_temperature = ghost_box_temperature.numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_density(d_dim);
    hier::IntVector offset_pressure(d_dim);
    hier::IntVector offset_temperature(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_density->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_pressure->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_temperature->getBox().isSpatiallyEqual(interior_box));
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_density;
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_temperature, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_density = num_ghosts_density;
        offset_pressure = num_ghosts_pressure;
        offset_temperature = num_ghosts_temperature;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_temperature->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_density = domain.lower() - ghost_box_density.lower();
        offset_pressure = domain.lower() - ghost_box_pressure.lower();
        offset_temperature = domain.lower() - ghost_box_temperature.lower();
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
    
    Real* const rho = data_density->getPointer(side_normal, 0);
    const Real* const p = data_pressure->getPointer(side_normal, 0);
    const Real* const T = data_temperature->getPointer(side_normal, 0);
    
    const Real& gamma = *(thermo_properties[0]);
    const Real& p_inf = *(thermo_properties[1]);
    const Real& c_v = *(thermo_properties[3]);
    
    computeDensity(
        rho,
        p,
        T,
        gamma,
        p_inf,
        c_v,
        offset_density,
        offset_pressure,
        offset_temperature,
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
EquationOfStateStiffenedGas::computeDensity(
    HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_temperature,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_thermo_properties,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_temperature);
    TBOX_ASSERT(data_thermo_properties);
    
    TBOX_ASSERT(data_thermo_properties->getDepth() == 4);
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_density = data_density->getGhostBox();
    const hier::IntVector ghostcell_dims_density = ghost_box_density.numberCells();
    
    const hier::Box ghost_box_pressure = data_pressure->getGhostBox();
    const hier::IntVector ghostcell_dims_pressure = ghost_box_pressure.numberCells();
    
    const hier::Box ghost_box_temperature = data_temperature->getGhostBox();
    const hier::IntVector ghostcell_dims_temperature = ghost_box_temperature.numberCells();
    
    const hier::Box ghost_box_thermo_properties = data_thermo_properties->getGhostBox();
    const hier::IntVector ghostcell_dims_thermo_properties = ghost_box_thermo_properties.numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_density(d_dim);
    hier::IntVector offset_pressure(d_dim);
    hier::IntVector offset_temperature(d_dim);
    hier::IntVector offset_thermo_properties(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
        const hier::IntVector num_ghosts_thermo_properties = data_thermo_properties->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_density->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_pressure->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_temperature->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_thermo_properties->getBox().isSpatiallyEqual(interior_box));
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_density;
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_temperature, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_thermo_properties, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_density = num_ghosts_density;
        offset_pressure = num_ghosts_pressure;
        offset_temperature = num_ghosts_temperature;
        offset_thermo_properties = num_ghosts_thermo_properties;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_temperature->getGhostBox().contains(domain));
        TBOX_ASSERT(data_thermo_properties->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_density = domain.lower() - ghost_box_density.lower();
        offset_pressure = domain.lower() - ghost_box_pressure.lower();
        offset_temperature = domain.lower() - ghost_box_temperature.lower();
        offset_thermo_properties = domain.lower() - ghost_box_thermo_properties.lower();
    }
    
    /*
     * Get the pointers to the cell data.
     */
    
    Real* const rho = data_density->getPointer(0);
    const Real* const p = data_pressure->getPointer(0);
    const Real* const T = data_temperature->getPointer(0);
    const Real* const gamma = data_thermo_properties->getPointer(0);
    const Real* const p_inf = data_thermo_properties->getPointer(1);
    const Real* const c_v = data_thermo_properties->getPointer(3);
    
    computeDensity(
        rho,
        p,
        T,
        gamma,
        p_inf,
        c_v,
        offset_density,
        offset_pressure,
        offset_temperature,
        offset_thermo_properties,
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
EquationOfStateStiffenedGas::computeDensity(
    HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_density,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_temperature,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& data_thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_density);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_temperature);
    TBOX_ASSERT(data_thermo_properties);
    
    TBOX_ASSERT(data_thermo_properties->getDepth() == 4);
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_density = data_density->getGhostBox();
    hier::IntVector ghostcell_dims_density = ghost_box_density.numberCells();
    
    const hier::Box ghost_box_pressure = data_pressure->getGhostBox();
    hier::IntVector ghostcell_dims_pressure = ghost_box_pressure.numberCells();
    
    const hier::Box ghost_box_temperature = data_temperature->getGhostBox();
    hier::IntVector ghostcell_dims_temperature = ghost_box_temperature.numberCells();
    
    const hier::Box ghost_box_thermo_properties = data_thermo_properties->getGhostBox();
    hier::IntVector ghostcell_dims_thermo_properties = ghost_box_thermo_properties.numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_density(d_dim);
    hier::IntVector offset_pressure(d_dim);
    hier::IntVector offset_temperature(d_dim);
    hier::IntVector offset_thermo_properties(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
        const hier::IntVector num_ghosts_thermo_properties = data_thermo_properties->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_density->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_pressure->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_temperature->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_thermo_properties->getBox().isSpatiallyEqual(interior_box));
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_density;
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_temperature, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_thermo_properties, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_density = num_ghosts_density;
        offset_pressure = num_ghosts_pressure;
        offset_temperature = num_ghosts_temperature;
        offset_thermo_properties = num_ghosts_thermo_properties;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_temperature->getGhostBox().contains(domain));
        TBOX_ASSERT(data_thermo_properties->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_density = domain.lower() - ghost_box_density.lower();
        offset_pressure = domain.lower() - ghost_box_pressure.lower();
        offset_temperature = domain.lower() - ghost_box_temperature.lower();
        offset_thermo_properties = domain.lower() - ghost_box_thermo_properties.lower();
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
    
    Real* const rho = data_density->getPointer(side_normal, 0);
    const Real* const p = data_pressure->getPointer(side_normal, 0);
    const Real* const T = data_temperature->getPointer(side_normal, 0);
    const Real* const gamma = data_thermo_properties->getPointer(side_normal, 0);
    const Real* const p_inf = data_thermo_properties->getPointer(side_normal, 1);
    const Real* const c_v = data_thermo_properties->getPointer(side_normal, 3);
    
    computeDensity(
        rho,
        p,
        T,
        gamma,
        p_inf,
        c_v,
        offset_density,
        offset_pressure,
        offset_temperature,
        offset_thermo_properties,
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
EquationOfStateStiffenedGas::computePressure(
    Real* const p,
    const Real* const rho,
    const Real* const epsilon,
    const Real& gamma,
    const Real& p_inf,
    const hier::IntVector& offset_pressure,
    const hier::IntVector& offset_density,
    const hier::IntVector& offset_internal_energy,
    const hier::IntVector& ghostcell_dims_pressure,
    const hier::IntVector& ghostcell_dims_density,
    const hier::IntVector& ghostcell_dims_internal_energy,
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
        
        const int offset_0_pressure = offset_pressure[0];
        const int offset_0_density = offset_density[0];
        const int offset_0_internal_energy = offset_internal_energy[0];
        
        HAMERS_PRAGMA_SIMD
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_pressure = i + offset_0_pressure;
            const int idx_density = i + offset_0_density;
            const int idx_internal_energy = i + offset_0_internal_energy;
            
            p[idx_pressure] = (gamma - Real(1))*rho[idx_density]*epsilon[idx_internal_energy] -
                gamma*p_inf;
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
        
        const int offset_0_pressure = offset_pressure[0];
        const int offset_1_pressure = offset_pressure[1];
        const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
        
        const int offset_0_density = offset_density[0];
        const int offset_1_density = offset_density[1];
        const int ghostcell_dim_0_density = ghostcell_dims_density[0];
        
        const int offset_0_internal_energy = offset_internal_energy[0];
        const int offset_1_internal_energy = offset_internal_energy[1];
        const int ghostcell_dim_0_internal_energy = ghostcell_dims_internal_energy[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_pressure = (i + offset_0_pressure) +
                    (j + offset_1_pressure)*ghostcell_dim_0_pressure;
                
                const int idx_density = (i + offset_0_density) +
                    (j + offset_1_density)*ghostcell_dim_0_density;
                
                const int idx_internal_energy = (i + offset_0_internal_energy) +
                    (j + offset_1_internal_energy)*ghostcell_dim_0_internal_energy;
                
                p[idx_pressure] = (gamma - Real(1))*rho[idx_density]*epsilon[idx_internal_energy] -
                    gamma*p_inf;
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
        
        const int offset_0_pressure = offset_pressure[0];
        const int offset_1_pressure = offset_pressure[1];
        const int offset_2_pressure = offset_pressure[2];
        const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
        const int ghostcell_dim_1_pressure = ghostcell_dims_pressure[1];
        
        const int offset_0_density = offset_density[0];
        const int offset_1_density = offset_density[1];
        const int offset_2_density = offset_density[2];
        const int ghostcell_dim_0_density = ghostcell_dims_density[0];
        const int ghostcell_dim_1_density = ghostcell_dims_density[1];
        
        const int offset_0_internal_energy = offset_internal_energy[0];
        const int offset_1_internal_energy = offset_internal_energy[1];
        const int offset_2_internal_energy = offset_internal_energy[2];
        const int ghostcell_dim_0_internal_energy = ghostcell_dims_internal_energy[0];
        const int ghostcell_dim_1_internal_energy = ghostcell_dims_internal_energy[1];
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_pressure = (i + offset_0_pressure) +
                        (j + offset_1_pressure)*ghostcell_dim_0_pressure +
                        (k + offset_2_pressure)*ghostcell_dim_0_pressure*
                            ghostcell_dim_1_pressure;
                    
                    const int idx_density = (i + offset_0_density) +
                        (j + offset_1_density)*ghostcell_dim_0_density +
                        (k + offset_2_density)*ghostcell_dim_0_density*
                            ghostcell_dim_1_density;
                    
                    const int idx_internal_energy = (i + offset_0_internal_energy) +
                        (j + offset_1_internal_energy)*ghostcell_dim_0_internal_energy +
                        (k + offset_2_internal_energy)*ghostcell_dim_0_internal_energy*
                            ghostcell_dim_1_internal_energy;
                    
                    p[idx_pressure] = (gamma - Real(1))*rho[idx_density]*epsilon[idx_internal_energy] -
                        gamma*p_inf;
                }
            }
        }
    }
}


/*
 * Compute the pressure.
 */
void
EquationOfStateStiffenedGas::computePressure(
    Real* const p,
    const Real* const rho,
    const Real* const epsilon,
    const Real* const gamma,
    const Real* const p_inf,
    const hier::IntVector& offset_pressure,
    const hier::IntVector& offset_density,
    const hier::IntVector& offset_internal_energy,
    const hier::IntVector& offset_thermo_properties,
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
         * Get the local lower index, numbers of cells in each dimension and offsets.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int offset_0_pressure = offset_pressure[0];
        const int offset_0_density = offset_density[0];
        const int offset_0_internal_energy = offset_internal_energy[0];
        const int offset_0_thermo_properties = offset_thermo_properties[0];
        
        HAMERS_PRAGMA_SIMD
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_pressure = i + offset_0_pressure;
            const int idx_density = i + offset_0_density;
            const int idx_internal_energy = i + offset_0_internal_energy;
            const int idx_thermo_properties = i + offset_0_thermo_properties;
            
            p[idx_pressure] = (gamma[idx_thermo_properties] - Real(1))*rho[idx_density]*
                epsilon[idx_internal_energy] - gamma[idx_thermo_properties]*p_inf[idx_thermo_properties];
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
        
        const int offset_0_pressure = offset_pressure[0];
        const int offset_1_pressure = offset_pressure[1];
        const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
        
        const int offset_0_density = offset_density[0];
        const int offset_1_density = offset_density[1];
        const int ghostcell_dim_0_density = ghostcell_dims_density[0];
        
        const int offset_0_internal_energy = offset_internal_energy[0];
        const int offset_1_internal_energy = offset_internal_energy[1];
        const int ghostcell_dim_0_internal_energy = ghostcell_dims_internal_energy[0];
        
        const int offset_0_thermo_properties = offset_thermo_properties[0];
        const int offset_1_thermo_properties = offset_thermo_properties[1];
        const int ghostcell_dim_0_thermo_properties = ghostcell_dims_thermo_properties[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_pressure = (i + offset_0_pressure) +
                    (j + offset_1_pressure)*ghostcell_dim_0_pressure;
                
                const int idx_density = (i + offset_0_density) +
                    (j + offset_1_density)*ghostcell_dim_0_density;
                
                const int idx_internal_energy = (i + offset_0_internal_energy) +
                    (j + offset_1_internal_energy)*ghostcell_dim_0_internal_energy;
                
                const int idx_thermo_properties = (i + offset_0_thermo_properties) +
                    (j + offset_1_thermo_properties)*ghostcell_dim_0_thermo_properties;
                
                p[idx_pressure] = (gamma[idx_thermo_properties] - Real(1))*rho[idx_density]*
                    epsilon[idx_internal_energy] - gamma[idx_thermo_properties]*p_inf[idx_thermo_properties];
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
        
        const int offset_0_pressure = offset_pressure[0];
        const int offset_1_pressure = offset_pressure[1];
        const int offset_2_pressure = offset_pressure[2];
        const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
        const int ghostcell_dim_1_pressure = ghostcell_dims_pressure[1];
        
        const int offset_0_density = offset_density[0];
        const int offset_1_density = offset_density[1];
        const int offset_2_density = offset_density[2];
        const int ghostcell_dim_0_density = ghostcell_dims_density[0];
        const int ghostcell_dim_1_density = ghostcell_dims_density[1];
        
        const int offset_0_internal_energy = offset_internal_energy[0];
        const int offset_1_internal_energy = offset_internal_energy[1];
        const int offset_2_internal_energy = offset_internal_energy[2];
        const int ghostcell_dim_0_internal_energy = ghostcell_dims_internal_energy[0];
        const int ghostcell_dim_1_internal_energy = ghostcell_dims_internal_energy[1];
        
        const int offset_0_thermo_properties = offset_thermo_properties[0];
        const int offset_1_thermo_properties = offset_thermo_properties[1];
        const int offset_2_thermo_properties = offset_thermo_properties[2];
        const int ghostcell_dim_0_thermo_properties = ghostcell_dims_thermo_properties[0];
        const int ghostcell_dim_1_thermo_properties = ghostcell_dims_thermo_properties[1];
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_pressure = (i + offset_0_pressure) +
                        (j + offset_1_pressure)*ghostcell_dim_0_pressure +
                        (k + offset_2_pressure)*ghostcell_dim_0_pressure*
                            ghostcell_dim_1_pressure;
                    
                    const int idx_density = (i + offset_0_density) +
                        (j + offset_1_density)*ghostcell_dim_0_density +
                        (k + offset_2_density)*ghostcell_dim_0_density*
                            ghostcell_dim_1_density;
                    
                    const int idx_internal_energy = (i + offset_0_internal_energy) +
                        (j + offset_1_internal_energy)*ghostcell_dim_0_internal_energy +
                        (k + offset_2_internal_energy)*ghostcell_dim_0_internal_energy*
                            ghostcell_dim_1_internal_energy;
                    
                    const int idx_thermo_properties = (i + offset_0_thermo_properties) +
                        (j + offset_1_thermo_properties)*ghostcell_dim_0_thermo_properties +
                        (k + offset_2_thermo_properties)*ghostcell_dim_0_thermo_properties*
                            ghostcell_dim_1_thermo_properties;
                    
                    p[idx_pressure] = (gamma[idx_thermo_properties] - Real(1))*rho[idx_density]*
                        epsilon[idx_internal_energy] - gamma[idx_thermo_properties]*p_inf[idx_thermo_properties];
                }
            }
        }
    }
}


/*
 * Compute the sound speed.
 */
void
EquationOfStateStiffenedGas::computeSoundSpeed(
    Real* const c,
    const Real* const rho,
    const Real* const p,
    const Real& gamma,
    const Real& p_inf,
    const hier::IntVector& offset_sound_speed,
    const hier::IntVector& offset_density,
    const hier::IntVector& offset_pressure,
    const hier::IntVector& ghostcell_dims_sound_speed,
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
        
        const int offset_0_sound_speed = offset_sound_speed[0];
        const int offset_0_density = offset_density[0];
        const int offset_0_pressure = offset_pressure[0];
        
        HAMERS_PRAGMA_SIMD
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_sound_speed = i + offset_0_sound_speed;
            const int idx_density = i + offset_0_density;
            const int idx_pressure = i + offset_0_pressure;
            
            c[idx_sound_speed] = std::sqrt(gamma*(p[idx_pressure] + p_inf)/rho[idx_density]);
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
        
        const int offset_0_sound_speed = offset_sound_speed[0];
        const int offset_1_sound_speed = offset_sound_speed[1];
        const int ghostcell_dim_0_sound_speed = ghostcell_dims_sound_speed[0];
        
        const int offset_0_density = offset_density[0];
        const int offset_1_density = offset_density[1];
        const int ghostcell_dim_0_density = ghostcell_dims_density[0];
        
        const int offset_0_pressure = offset_pressure[0];
        const int offset_1_pressure = offset_pressure[1];
        const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_sound_speed = (i + offset_0_sound_speed) +
                    (j + offset_1_sound_speed)*ghostcell_dim_0_sound_speed;
                
                const int idx_density = (i + offset_0_density) +
                    (j + offset_1_density)*ghostcell_dim_0_density;
                
                const int idx_pressure = (i + offset_0_pressure) +
                    (j + offset_1_pressure)*ghostcell_dim_0_pressure;
                
                c[idx_sound_speed] = std::sqrt(gamma*(p[idx_pressure] + p_inf)/rho[idx_density]);
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
        
        const int offset_0_sound_speed = offset_sound_speed[0];
        const int offset_1_sound_speed = offset_sound_speed[1];
        const int offset_2_sound_speed = offset_sound_speed[2];
        const int ghostcell_dim_0_sound_speed = ghostcell_dims_sound_speed[0];
        const int ghostcell_dim_1_sound_speed = ghostcell_dims_sound_speed[1];
        
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
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_sound_speed = (i + offset_0_sound_speed) +
                        (j + offset_1_sound_speed)*ghostcell_dim_0_sound_speed +
                        (k + offset_2_sound_speed)*ghostcell_dim_0_sound_speed*
                            ghostcell_dim_1_sound_speed;
                    
                    const int idx_density = (i + offset_0_density) +
                        (j + offset_1_density)*ghostcell_dim_0_density +
                        (k + offset_2_density)*ghostcell_dim_0_density*
                            ghostcell_dim_1_density;
                    
                    const int idx_pressure = (i + offset_0_pressure) +
                        (j + offset_1_pressure)*ghostcell_dim_0_pressure +
                        (k + offset_2_pressure)*ghostcell_dim_0_pressure*
                            ghostcell_dim_1_pressure;
                    
                    c[idx_sound_speed] = std::sqrt(gamma*(p[idx_pressure] + p_inf)/rho[idx_density]);
                }
            }
        }
    }
}


/*
 * Compute the sound speed.
 */
void
EquationOfStateStiffenedGas::computeSoundSpeed(
    Real* const c,
    const Real* const rho,
    const Real* const p,
    const Real* const gamma,
    const Real* const p_inf,
    const hier::IntVector& offset_sound_speed,
    const hier::IntVector& offset_density,
    const hier::IntVector& offset_pressure,
    const hier::IntVector& offset_thermo_properties,
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
         * Get the local lower index, numbers of cells in each dimension and offsets.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int offset_0_sound_speed = offset_sound_speed[0];
        const int offset_0_density = offset_density[0];
        const int offset_0_pressure = offset_pressure[0];
        const int offset_0_thermo_properties = offset_thermo_properties[0];
        
        HAMERS_PRAGMA_SIMD
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_sound_speed = i + offset_0_sound_speed;
            const int idx_density = i + offset_0_density;
            const int idx_pressure = i + offset_0_pressure;
            const int idx_thermo_properties = i + offset_0_thermo_properties;
            
            c[idx_sound_speed] = std::sqrt(gamma[idx_thermo_properties]*(p[idx_pressure] +
                p_inf[idx_thermo_properties])/rho[idx_density]);
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
        
        const int offset_0_sound_speed = offset_sound_speed[0];
        const int offset_1_sound_speed = offset_sound_speed[1];
        const int ghostcell_dim_0_sound_speed = ghostcell_dims_sound_speed[0];
        
        const int offset_0_density = offset_density[0];
        const int offset_1_density = offset_density[1];
        const int ghostcell_dim_0_density = ghostcell_dims_density[0];
        
        const int offset_0_pressure = offset_pressure[0];
        const int offset_1_pressure = offset_pressure[1];
        const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
        
        const int offset_0_thermo_properties = offset_thermo_properties[0];
        const int offset_1_thermo_properties = offset_thermo_properties[1];
        const int ghostcell_dim_0_thermo_properties = ghostcell_dims_thermo_properties[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_sound_speed = (i + offset_0_sound_speed) +
                    (j + offset_1_sound_speed)*ghostcell_dim_0_sound_speed;
                
                const int idx_density = (i + offset_0_density) +
                    (j + offset_1_density)*ghostcell_dim_0_density;
                
                const int idx_pressure = (i + offset_0_pressure) +
                    (j + offset_1_pressure)*ghostcell_dim_0_pressure;
                
                const int idx_thermo_properties = (i + offset_0_thermo_properties) +
                    (j + offset_1_thermo_properties)*ghostcell_dim_0_thermo_properties;
                
                c[idx_sound_speed] = std::sqrt(gamma[idx_thermo_properties]*(p[idx_pressure] +
                    p_inf[idx_thermo_properties])/rho[idx_density]);
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
        
        const int offset_0_sound_speed = offset_sound_speed[0];
        const int offset_1_sound_speed = offset_sound_speed[1];
        const int offset_2_sound_speed = offset_sound_speed[2];
        const int ghostcell_dim_0_sound_speed = ghostcell_dims_sound_speed[0];
        const int ghostcell_dim_1_sound_speed = ghostcell_dims_sound_speed[1];
        
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
        
        const int offset_0_thermo_properties = offset_thermo_properties[0];
        const int offset_1_thermo_properties = offset_thermo_properties[1];
        const int offset_2_thermo_properties = offset_thermo_properties[2];
        const int ghostcell_dim_0_thermo_properties = ghostcell_dims_thermo_properties[0];
        const int ghostcell_dim_1_thermo_properties = ghostcell_dims_thermo_properties[1];
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_sound_speed = (i + offset_0_sound_speed) +
                        (j + offset_1_sound_speed)*ghostcell_dim_0_sound_speed +
                        (k + offset_2_sound_speed)*ghostcell_dim_0_sound_speed*
                            ghostcell_dim_1_sound_speed;
                    
                    const int idx_density = (i + offset_0_density) +
                        (j + offset_1_density)*ghostcell_dim_0_density +
                        (k + offset_2_density)*ghostcell_dim_0_density*
                            ghostcell_dim_1_density;
                    
                    const int idx_pressure = (i + offset_0_pressure) +
                        (j + offset_1_pressure)*ghostcell_dim_0_pressure +
                        (k + offset_2_pressure)*ghostcell_dim_0_pressure*
                            ghostcell_dim_1_pressure;
                    
                    const int idx_thermo_properties = (i + offset_0_thermo_properties) +
                        (j + offset_1_thermo_properties)*ghostcell_dim_0_thermo_properties +
                        (k + offset_2_thermo_properties)*ghostcell_dim_0_thermo_properties*
                            ghostcell_dim_1_thermo_properties;
                    
                    c[idx_sound_speed] = std::sqrt(gamma[idx_thermo_properties]*(p[idx_pressure] +
                        p_inf[idx_thermo_properties])/rho[idx_density]);
                }
            }
        }
    }
}


/*
 * Compute the specific internal energy.
 */
void
EquationOfStateStiffenedGas::computeInternalEnergy(
    Real* const epsilon,
    const Real* const rho,
    const Real* const p,
    const Real& gamma,
    const Real& p_inf,
    const hier::IntVector& offset_internal_energy,
    const hier::IntVector& offset_density,
    const hier::IntVector& offset_pressure,
    const hier::IntVector& ghostcell_dims_internal_energy,
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
        
        const int offset_0_internal_energy = offset_internal_energy[0];
        const int offset_0_density = offset_density[0];
        const int offset_0_pressure = offset_pressure[0];
        
        HAMERS_PRAGMA_SIMD
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_internal_energy = i + offset_0_internal_energy;
            const int idx_density = i + offset_0_density;
            const int idx_pressure = i + offset_0_pressure;
            
            epsilon[idx_internal_energy] = (p[idx_pressure] + gamma*p_inf)/((gamma - Real(1))*
                rho[idx_density]);
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
        
        const int offset_0_internal_energy = offset_internal_energy[0];
        const int offset_1_internal_energy = offset_internal_energy[1];
        const int ghostcell_dim_0_internal_energy = ghostcell_dims_internal_energy[0];
        
        const int offset_0_density = offset_density[0];
        const int offset_1_density = offset_density[1];
        const int ghostcell_dim_0_density = ghostcell_dims_density[0];
        
        const int offset_0_pressure = offset_pressure[0];
        const int offset_1_pressure = offset_pressure[1];
        const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_internal_energy = (i + offset_0_internal_energy) +
                    (j + offset_1_internal_energy)*ghostcell_dim_0_internal_energy;
                
                const int idx_density = (i + offset_0_density) +
                    (j + offset_1_density)*ghostcell_dim_0_density;
                
                const int idx_pressure = (i + offset_0_pressure) +
                    (j + offset_1_pressure)*ghostcell_dim_0_pressure;
                
                epsilon[idx_internal_energy] = (p[idx_pressure] + gamma*p_inf)/((gamma - Real(1))*
                    rho[idx_density]);
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
        
        const int offset_0_internal_energy = offset_internal_energy[0];
        const int offset_1_internal_energy = offset_internal_energy[1];
        const int offset_2_internal_energy = offset_internal_energy[2];
        const int ghostcell_dim_0_internal_energy = ghostcell_dims_internal_energy[0];
        const int ghostcell_dim_1_internal_energy = ghostcell_dims_internal_energy[1];
        
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
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_internal_energy = (i + offset_0_internal_energy) +
                        (j + offset_1_internal_energy)*ghostcell_dim_0_internal_energy +
                        (k + offset_2_internal_energy)*ghostcell_dim_0_internal_energy*
                            ghostcell_dim_1_internal_energy;
                    
                    const int idx_density = (i + offset_0_density) +
                        (j + offset_1_density)*ghostcell_dim_0_density +
                        (k + offset_2_density)*ghostcell_dim_0_density*
                            ghostcell_dim_1_density;
                    
                    const int idx_pressure = (i + offset_0_pressure) +
                        (j + offset_1_pressure)*ghostcell_dim_0_pressure +
                        (k + offset_2_pressure)*ghostcell_dim_0_pressure*
                            ghostcell_dim_1_pressure;
                    
                    epsilon[idx_internal_energy] = (p[idx_pressure] + gamma*p_inf)/((gamma - Real(1))*
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
EquationOfStateStiffenedGas::computeInternalEnergy(
    Real* const epsilon,
    const Real* const rho,
    const Real* const p,
    const Real* const gamma,
    const Real* const p_inf,
    const hier::IntVector& offset_internal_energy,
    const hier::IntVector& offset_density,
    const hier::IntVector& offset_pressure,
    const hier::IntVector& offset_thermo_properties,
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
         * Get the local lower index, numbers of cells in each dimension and offsets.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int offset_0_internal_energy = offset_internal_energy[0];
        const int offset_0_density = offset_density[0];
        const int offset_0_pressure = offset_pressure[0];
        const int offset_0_thermo_properties = offset_thermo_properties[0];
        
        HAMERS_PRAGMA_SIMD
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_internal_energy = i + offset_0_internal_energy;
            const int idx_density = i + offset_0_density;
            const int idx_pressure = i + offset_0_pressure;
            const int idx_thermo_properties = i + offset_0_thermo_properties;
            
            epsilon[idx_internal_energy] = (p[idx_pressure] + gamma[idx_thermo_properties]*
                p_inf[idx_thermo_properties])/((gamma[idx_thermo_properties] - Real(1))*rho[idx_density]);
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
        
        const int offset_0_internal_energy = offset_internal_energy[0];
        const int offset_1_internal_energy = offset_internal_energy[1];
        const int ghostcell_dim_0_internal_energy = ghostcell_dims_internal_energy[0];
        
        const int offset_0_density = offset_density[0];
        const int offset_1_density = offset_density[1];
        const int ghostcell_dim_0_density = ghostcell_dims_density[0];
        
        const int offset_0_pressure = offset_pressure[0];
        const int offset_1_pressure = offset_pressure[1];
        const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
        
        const int offset_0_thermo_properties = offset_thermo_properties[0];
        const int offset_1_thermo_properties = offset_thermo_properties[1];
        const int ghostcell_dim_0_thermo_properties = ghostcell_dims_thermo_properties[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_internal_energy = (i + offset_0_internal_energy) +
                    (j + offset_1_internal_energy)*ghostcell_dim_0_internal_energy;
                
                const int idx_density = (i + offset_0_density) +
                    (j + offset_1_density)*ghostcell_dim_0_density;
                
                const int idx_pressure = (i + offset_0_pressure) +
                    (j + offset_1_pressure)*ghostcell_dim_0_pressure;
                
                const int idx_thermo_properties = (i + offset_0_thermo_properties) +
                    (j + offset_1_thermo_properties)*ghostcell_dim_0_thermo_properties;
                
                epsilon[idx_internal_energy] = (p[idx_pressure] + gamma[idx_thermo_properties]*
                    p_inf[idx_thermo_properties])/((gamma[idx_thermo_properties] - Real(1))*rho[idx_density]);
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
        
        const int offset_0_internal_energy = offset_internal_energy[0];
        const int offset_1_internal_energy = offset_internal_energy[1];
        const int offset_2_internal_energy = offset_internal_energy[2];
        const int ghostcell_dim_0_internal_energy = ghostcell_dims_internal_energy[0];
        const int ghostcell_dim_1_internal_energy = ghostcell_dims_internal_energy[1];
        
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
        
        const int offset_0_thermo_properties = offset_thermo_properties[0];
        const int offset_1_thermo_properties = offset_thermo_properties[1];
        const int offset_2_thermo_properties = offset_thermo_properties[2];
        const int ghostcell_dim_0_thermo_properties = ghostcell_dims_thermo_properties[0];
        const int ghostcell_dim_1_thermo_properties = ghostcell_dims_thermo_properties[1];
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_internal_energy = (i + offset_0_internal_energy) +
                        (j + offset_1_internal_energy)*ghostcell_dim_0_internal_energy +
                        (k + offset_2_internal_energy)*ghostcell_dim_0_internal_energy*
                            ghostcell_dim_1_internal_energy;
                    
                    const int idx_density = (i + offset_0_density) +
                        (j + offset_1_density)*ghostcell_dim_0_density +
                        (k + offset_2_density)*ghostcell_dim_0_density*
                            ghostcell_dim_1_density;
                    
                    const int idx_pressure = (i + offset_0_pressure) +
                        (j + offset_1_pressure)*ghostcell_dim_0_pressure +
                        (k + offset_2_pressure)*ghostcell_dim_0_pressure*
                            ghostcell_dim_1_pressure;
                    
                    const int idx_thermo_properties = (i + offset_0_thermo_properties) +
                        (j + offset_1_thermo_properties)*ghostcell_dim_0_thermo_properties +
                        (k + offset_2_thermo_properties)*ghostcell_dim_0_thermo_properties*
                            ghostcell_dim_1_thermo_properties;
                    
                    epsilon[idx_internal_energy] = (p[idx_pressure] + gamma[idx_thermo_properties]*
                        p_inf[idx_thermo_properties])/((gamma[idx_thermo_properties] - Real(1))*rho[idx_density]);
                }
            }
        }
    }
}


/*
 * Compute the specific enthalpy.
 */
void
EquationOfStateStiffenedGas::computeEnthalpy(
    Real* const h,
    const Real* const rho,
    const Real* const p,
    const Real& gamma,
    const Real& p_inf,
    const hier::IntVector& offset_enthalpy,
    const hier::IntVector& offset_density,
    const hier::IntVector& offset_pressure,
    const hier::IntVector& ghostcell_dims_enthalpy,
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
        
        const int offset_0_enthalpy = offset_enthalpy[0];
        const int offset_0_density = offset_density[0];
        const int offset_0_pressure = offset_pressure[0];
        
        HAMERS_PRAGMA_SIMD
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_enthalpy = i + offset_0_enthalpy;
            const int idx_density = i + offset_0_density;
            const int idx_pressure = i + offset_0_pressure;
            
            h[idx_enthalpy] = gamma*(p[idx_pressure] + p_inf)/((gamma - Real(1))*rho[idx_density]);
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
        
        const int offset_0_enthalpy = offset_enthalpy[0];
        const int offset_1_enthalpy = offset_enthalpy[1];
        const int ghostcell_dim_0_enthalpy = ghostcell_dims_enthalpy[0];
        
        const int offset_0_density = offset_density[0];
        const int offset_1_density = offset_density[1];
        const int ghostcell_dim_0_density = ghostcell_dims_density[0];
        
        const int offset_0_pressure = offset_pressure[0];
        const int offset_1_pressure = offset_pressure[1];
        const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_enthalpy = (i + offset_0_enthalpy) +
                    (j + offset_1_enthalpy)*ghostcell_dim_0_enthalpy;
                
                const int idx_density = (i + offset_0_density) +
                    (j + offset_1_density)*ghostcell_dim_0_density;
                
                const int idx_pressure = (i + offset_0_pressure) +
                    (j + offset_1_pressure)*ghostcell_dim_0_pressure;
                
                h[idx_enthalpy] = gamma*(p[idx_pressure] + p_inf)/((gamma - Real(1))*rho[idx_density]);
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
        
        const int offset_0_enthalpy = offset_enthalpy[0];
        const int offset_1_enthalpy = offset_enthalpy[1];
        const int offset_2_enthalpy = offset_enthalpy[2];
        const int ghostcell_dim_0_enthalpy = ghostcell_dims_enthalpy[0];
        const int ghostcell_dim_1_enthalpy = ghostcell_dims_enthalpy[1];
        
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
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_enthalpy = (i + offset_0_enthalpy) +
                        (j + offset_1_enthalpy)*ghostcell_dim_0_enthalpy +
                        (k + offset_2_enthalpy)*ghostcell_dim_0_enthalpy*
                            ghostcell_dim_1_enthalpy;
                    
                    const int idx_density = (i + offset_0_density) +
                        (j + offset_1_density)*ghostcell_dim_0_density +
                        (k + offset_2_density)*ghostcell_dim_0_density*
                            ghostcell_dim_1_density;
                    
                    const int idx_pressure = (i + offset_0_pressure) +
                        (j + offset_1_pressure)*ghostcell_dim_0_pressure +
                        (k + offset_2_pressure)*ghostcell_dim_0_pressure*
                            ghostcell_dim_1_pressure;
                    
                    h[idx_enthalpy] = gamma*(p[idx_pressure] + p_inf)/((gamma - Real(1))*rho[idx_density]);
                }
            }
        }
    }
}


/*
 * Compute the specific enthalpy.
 */
void
EquationOfStateStiffenedGas::computeEnthalpy(
    Real* const h,
    const Real* const rho,
    const Real* const p,
    const Real* const gamma,
    const Real* const p_inf,
    const hier::IntVector& offset_enthalpy,
    const hier::IntVector& offset_density,
    const hier::IntVector& offset_pressure,
    const hier::IntVector& offset_thermo_properties,
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
         * Get the local lower index, numbers of cells in each dimension and offsets.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int offset_0_enthalpy = offset_enthalpy[0];
        const int offset_0_density = offset_density[0];
        const int offset_0_pressure = offset_pressure[0];
        const int offset_0_thermo_properties = offset_thermo_properties[0];
        
        HAMERS_PRAGMA_SIMD
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_enthalpy = i + offset_0_enthalpy;
            const int idx_density = i + offset_0_density;
            const int idx_pressure = i + offset_0_pressure;
            const int idx_thermo_properties = i + offset_0_thermo_properties;
            
            h[idx_enthalpy] = gamma[idx_thermo_properties]*(p[idx_pressure] + p_inf[idx_thermo_properties])/
                ((gamma[idx_thermo_properties] - Real(1))*rho[idx_density]);
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
        
        const int offset_0_enthalpy = offset_enthalpy[0];
        const int offset_1_enthalpy = offset_enthalpy[1];
        const int ghostcell_dim_0_enthalpy = ghostcell_dims_enthalpy[0];
        
        const int offset_0_density = offset_density[0];
        const int offset_1_density = offset_density[1];
        const int ghostcell_dim_0_density = ghostcell_dims_density[0];
        
        const int offset_0_pressure = offset_pressure[0];
        const int offset_1_pressure = offset_pressure[1];
        const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
        
        const int offset_0_thermo_properties = offset_thermo_properties[0];
        const int offset_1_thermo_properties = offset_thermo_properties[1];
        const int ghostcell_dim_0_thermo_properties = ghostcell_dims_thermo_properties[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_enthalpy = (i + offset_0_enthalpy) +
                    (j + offset_1_enthalpy)*ghostcell_dim_0_enthalpy;
                
                const int idx_density = (i + offset_0_density) +
                    (j + offset_1_density)*ghostcell_dim_0_density;
                
                const int idx_pressure = (i + offset_0_pressure) +
                    (j + offset_1_pressure)*ghostcell_dim_0_pressure;
                
                const int idx_thermo_properties = (i + offset_0_thermo_properties) +
                    (j + offset_1_thermo_properties)*ghostcell_dim_0_thermo_properties;
                
                h[idx_enthalpy] = gamma[idx_thermo_properties]*(p[idx_pressure] + p_inf[idx_thermo_properties])/
                    ((gamma[idx_thermo_properties] - Real(1))*rho[idx_density]);
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
        
        const int offset_0_enthalpy = offset_enthalpy[0];
        const int offset_1_enthalpy = offset_enthalpy[1];
        const int offset_2_enthalpy = offset_enthalpy[2];
        const int ghostcell_dim_0_enthalpy = ghostcell_dims_enthalpy[0];
        const int ghostcell_dim_1_enthalpy = ghostcell_dims_enthalpy[1];
        
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
        
        const int offset_0_thermo_properties = offset_thermo_properties[0];
        const int offset_1_thermo_properties = offset_thermo_properties[1];
        const int offset_2_thermo_properties = offset_thermo_properties[2];
        const int ghostcell_dim_0_thermo_properties = ghostcell_dims_thermo_properties[0];
        const int ghostcell_dim_1_thermo_properties = ghostcell_dims_thermo_properties[1];
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_enthalpy = (i + offset_0_enthalpy) +
                        (j + offset_1_enthalpy)*ghostcell_dim_0_enthalpy +
                        (k + offset_2_enthalpy)*ghostcell_dim_0_enthalpy*
                            ghostcell_dim_1_enthalpy;
                    
                    const int idx_density = (i + offset_0_density) +
                        (j + offset_1_density)*ghostcell_dim_0_density +
                        (k + offset_2_density)*ghostcell_dim_0_density*
                            ghostcell_dim_1_density;
                    
                    const int idx_pressure = (i + offset_0_pressure) +
                        (j + offset_1_pressure)*ghostcell_dim_0_pressure +
                        (k + offset_2_pressure)*ghostcell_dim_0_pressure*
                            ghostcell_dim_1_pressure;
                    
                    const int idx_thermo_properties = (i + offset_0_thermo_properties) +
                        (j + offset_1_thermo_properties)*ghostcell_dim_0_thermo_properties +
                        (k + offset_2_thermo_properties)*ghostcell_dim_0_thermo_properties*
                            ghostcell_dim_1_thermo_properties;
                    
                    h[idx_enthalpy] = gamma[idx_thermo_properties]*(p[idx_pressure] + p_inf[idx_thermo_properties])/
                        ((gamma[idx_thermo_properties] - Real(1))*rho[idx_density]);
                }
            }
        }
    }
}


/*
 * Compute the temperature.
 */
void
EquationOfStateStiffenedGas::computeTemperature(
    Real* const T,
    const Real* const rho,
    const Real* const p,
    const Real& gamma,
    const Real& p_inf,
    const Real& c_v,
    const hier::IntVector& offset_temperature,
    const hier::IntVector& offset_density,
    const hier::IntVector& offset_pressure,
    const hier::IntVector& ghostcell_dims_temperature,
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
        
        const int offset_0_temperature = offset_temperature[0];
        const int offset_0_density = offset_density[0];
        const int offset_0_pressure = offset_pressure[0];
        
        HAMERS_PRAGMA_SIMD
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_temperature = i + offset_0_temperature;
            const int idx_density = i + offset_0_density;
            const int idx_pressure = i + offset_0_pressure;
            
            T[idx_temperature] = (p[idx_pressure] + p_inf)/((gamma - Real(1))*c_v*rho[idx_density]);
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
        
        const int offset_0_temperature = offset_temperature[0];
        const int offset_1_temperature = offset_temperature[1];
        const int ghostcell_dim_0_temperature = ghostcell_dims_temperature[0];
        
        const int offset_0_density = offset_density[0];
        const int offset_1_density = offset_density[1];
        const int ghostcell_dim_0_density = ghostcell_dims_density[0];
        
        const int offset_0_pressure = offset_pressure[0];
        const int offset_1_pressure = offset_pressure[1];
        const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_temperature = (i + offset_0_temperature) +
                    (j + offset_1_temperature)*ghostcell_dim_0_temperature;
                
                const int idx_density = (i + offset_0_density) +
                    (j + offset_1_density)*ghostcell_dim_0_density;
                
                const int idx_pressure = (i + offset_0_pressure) +
                    (j + offset_1_pressure)*ghostcell_dim_0_pressure;
                
                T[idx_temperature] = (p[idx_pressure] + p_inf)/((gamma - Real(1))*c_v*rho[idx_density]);
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
        
        const int offset_0_temperature = offset_temperature[0];
        const int offset_1_temperature = offset_temperature[1];
        const int offset_2_temperature = offset_temperature[2];
        const int ghostcell_dim_0_temperature = ghostcell_dims_temperature[0];
        const int ghostcell_dim_1_temperature = ghostcell_dims_temperature[1];
        
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
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_temperature = (i + offset_0_temperature) +
                        (j + offset_1_temperature)*ghostcell_dim_0_temperature +
                        (k + offset_2_temperature)*ghostcell_dim_0_temperature*
                            ghostcell_dim_1_temperature;
                    
                    const int idx_density = (i + offset_0_density) +
                        (j + offset_1_density)*ghostcell_dim_0_density +
                        (k + offset_2_density)*ghostcell_dim_0_density*
                            ghostcell_dim_1_density;
                    
                    const int idx_pressure = (i + offset_0_pressure) +
                        (j + offset_1_pressure)*ghostcell_dim_0_pressure +
                        (k + offset_2_pressure)*ghostcell_dim_0_pressure*
                            ghostcell_dim_1_pressure;
                    
                    T[idx_temperature] = (p[idx_pressure] + p_inf)/((gamma - Real(1))*c_v*rho[idx_density]);
                }
            }
        }
    }
}


/*
 * Compute the temperature.
 */
void
EquationOfStateStiffenedGas::computeTemperature(
    Real* const T,
    const Real* const rho,
    const Real* const p,
    const Real* const gamma,
    const Real* const p_inf,
    const Real* const c_v,
    const hier::IntVector& offset_temperature,
    const hier::IntVector& offset_density,
    const hier::IntVector& offset_pressure,
    const hier::IntVector& offset_thermo_properties,
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
         * Get the local lower index, numbers of cells in each dimension and offsets.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int offset_0_temperature = offset_temperature[0];
        const int offset_0_density = offset_density[0];
        const int offset_0_pressure = offset_pressure[0];
        const int offset_0_thermo_properties = offset_thermo_properties[0];
        
        HAMERS_PRAGMA_SIMD
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_temperature = i + offset_0_temperature;
            const int idx_density = i + offset_0_density;
            const int idx_pressure = i + offset_0_pressure;
            const int idx_thermo_properties = i + offset_0_thermo_properties;
            
            T[idx_temperature] = (p[idx_pressure] + p_inf[idx_thermo_properties])/
                ((gamma[idx_thermo_properties] - Real(1))*c_v[idx_thermo_properties]*rho[idx_density]);
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
        
        const int offset_0_temperature = offset_temperature[0];
        const int offset_1_temperature = offset_temperature[1];
        const int ghostcell_dim_0_temperature = ghostcell_dims_temperature[0];
        
        const int offset_0_density = offset_density[0];
        const int offset_1_density = offset_density[1];
        const int ghostcell_dim_0_density = ghostcell_dims_density[0];
        
        const int offset_0_pressure = offset_pressure[0];
        const int offset_1_pressure = offset_pressure[1];
        const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
        
        const int offset_0_thermo_properties = offset_thermo_properties[0];
        const int offset_1_thermo_properties = offset_thermo_properties[1];
        const int ghostcell_dim_0_thermo_properties = ghostcell_dims_thermo_properties[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_temperature = (i + offset_0_temperature) +
                    (j + offset_1_temperature)*ghostcell_dim_0_temperature;
                
                const int idx_density = (i + offset_0_density) +
                    (j + offset_1_density)*ghostcell_dim_0_density;
                
                const int idx_pressure = (i + offset_0_pressure) +
                    (j + offset_1_pressure)*ghostcell_dim_0_pressure;
                
                const int idx_thermo_properties = (i + offset_0_thermo_properties) +
                    (j + offset_1_thermo_properties)*ghostcell_dim_0_thermo_properties;
                
                T[idx_temperature] = (p[idx_pressure] + p_inf[idx_thermo_properties])/
                    ((gamma[idx_thermo_properties] - Real(1))*c_v[idx_thermo_properties]*rho[idx_density]);
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
        
        const int offset_0_temperature = offset_temperature[0];
        const int offset_1_temperature = offset_temperature[1];
        const int offset_2_temperature = offset_temperature[2];
        const int ghostcell_dim_0_temperature = ghostcell_dims_temperature[0];
        const int ghostcell_dim_1_temperature = ghostcell_dims_temperature[1];
        
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
        
        const int offset_0_thermo_properties = offset_thermo_properties[0];
        const int offset_1_thermo_properties = offset_thermo_properties[1];
        const int offset_2_thermo_properties = offset_thermo_properties[2];
        const int ghostcell_dim_0_thermo_properties = ghostcell_dims_thermo_properties[0];
        const int ghostcell_dim_1_thermo_properties = ghostcell_dims_thermo_properties[1];
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_temperature = (i + offset_0_temperature) +
                        (j + offset_1_temperature)*ghostcell_dim_0_temperature +
                        (k + offset_2_temperature)*ghostcell_dim_0_temperature*
                            ghostcell_dim_1_temperature;
                    
                    const int idx_density = (i + offset_0_density) +
                        (j + offset_1_density)*ghostcell_dim_0_density +
                        (k + offset_2_density)*ghostcell_dim_0_density*
                            ghostcell_dim_1_density;
                    
                    const int idx_pressure = (i + offset_0_pressure) +
                        (j + offset_1_pressure)*ghostcell_dim_0_pressure +
                        (k + offset_2_pressure)*ghostcell_dim_0_pressure*
                            ghostcell_dim_1_pressure;
                    
                    const int idx_thermo_properties = (i + offset_0_thermo_properties) +
                        (j + offset_1_thermo_properties)*ghostcell_dim_0_thermo_properties +
                        (k + offset_2_thermo_properties)*ghostcell_dim_0_thermo_properties*
                            ghostcell_dim_1_thermo_properties;
                    
                    T[idx_temperature] = (p[idx_pressure] + p_inf[idx_thermo_properties])/
                        ((gamma[idx_thermo_properties] - Real(1))*c_v[idx_thermo_properties]*rho[idx_density]);
                }
            }
        }
    }
}


/*
 * Compute the specific internal energy from temperature.
 */
void
EquationOfStateStiffenedGas::computeInternalEnergyFromTemperature(
    Real* const epsilon,
    const Real* const rho,
    const Real* const T,
    const Real& p_inf,
    const Real& c_v,
    const hier::IntVector& offset_internal_energy,
    const hier::IntVector& offset_density,
    const hier::IntVector& offset_temperature,
    const hier::IntVector& ghostcell_dims_internal_energy,
    const hier::IntVector& ghostcell_dims_density,
    const hier::IntVector& ghostcell_dims_temperature,
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
        
        const int offset_0_internal_energy = offset_internal_energy[0];
        const int offset_0_density = offset_density[0];
        const int offset_0_temperature = offset_temperature[0];
        
        HAMERS_PRAGMA_SIMD
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_internal_energy = i + offset_0_internal_energy;
            const int idx_density = i + offset_0_density;
            const int idx_temperature = i + offset_0_temperature;
            
            epsilon[idx_internal_energy] = c_v*T[idx_temperature] + p_inf/rho[idx_density];
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
        
        const int offset_0_internal_energy = offset_internal_energy[0];
        const int offset_1_internal_energy = offset_internal_energy[1];
        const int ghostcell_dim_0_internal_energy = ghostcell_dims_internal_energy[0];
        
        const int offset_0_density = offset_density[0];
        const int offset_1_density = offset_density[1];
        const int ghostcell_dim_0_density = ghostcell_dims_density[0];
        
        const int offset_0_temperature = offset_temperature[0];
        const int offset_1_temperature = offset_temperature[1];
        const int ghostcell_dim_0_temperature = ghostcell_dims_temperature[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_internal_energy = (i + offset_0_internal_energy) +
                    (j + offset_1_internal_energy)*ghostcell_dim_0_internal_energy;
                
                const int idx_density = (i + offset_0_density) +
                    (j + offset_1_density)*ghostcell_dim_0_density;
                
                const int idx_temperature = (i + offset_0_temperature) +
                    (j + offset_1_temperature)*ghostcell_dim_0_temperature;
                
                epsilon[idx_internal_energy] = c_v*T[idx_temperature] + p_inf/rho[idx_density];
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
        
        const int offset_0_internal_energy = offset_internal_energy[0];
        const int offset_1_internal_energy = offset_internal_energy[1];
        const int offset_2_internal_energy = offset_internal_energy[2];
        const int ghostcell_dim_0_internal_energy = ghostcell_dims_internal_energy[0];
        const int ghostcell_dim_1_internal_energy = ghostcell_dims_internal_energy[1];
        
        const int offset_0_density = offset_density[0];
        const int offset_1_density = offset_density[1];
        const int offset_2_density = offset_density[2];
        const int ghostcell_dim_0_density = ghostcell_dims_density[0];
        const int ghostcell_dim_1_density = ghostcell_dims_density[1];
        
        const int offset_0_temperature = offset_temperature[0];
        const int offset_1_temperature = offset_temperature[1];
        const int offset_2_temperature = offset_temperature[2];
        const int ghostcell_dim_0_temperature = ghostcell_dims_temperature[0];
        const int ghostcell_dim_1_temperature = ghostcell_dims_temperature[1];
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_internal_energy = (i + offset_0_internal_energy) +
                        (j + offset_1_internal_energy)*ghostcell_dim_0_internal_energy +
                        (k + offset_2_internal_energy)*ghostcell_dim_0_internal_energy*
                            ghostcell_dim_1_internal_energy;
                    
                    const int idx_density = (i + offset_0_density) +
                        (j + offset_1_density)*ghostcell_dim_0_density +
                        (k + offset_2_density)*ghostcell_dim_0_density*
                            ghostcell_dim_1_density;
                    
                    const int idx_temperature = (i + offset_0_temperature) +
                        (j + offset_1_temperature)*ghostcell_dim_0_temperature +
                        (k + offset_2_temperature)*ghostcell_dim_0_temperature*
                            ghostcell_dim_1_temperature;
                    
                    epsilon[idx_internal_energy] = c_v*T[idx_temperature] + p_inf/rho[idx_density];
                }
            }
        }
    }
}


/*
 * Compute the specific internal energy from temperature.
 */
void
EquationOfStateStiffenedGas::computeInternalEnergyFromTemperature(
    Real* const epsilon,
    const Real* const rho,
    const Real* const T,
    const Real* const p_inf,
    const Real* const c_v,
    const hier::IntVector& offset_internal_energy,
    const hier::IntVector& offset_density,
    const hier::IntVector& offset_temperature,
    const hier::IntVector& offset_thermo_properties,
    const hier::IntVector& ghostcell_dims_internal_energy,
    const hier::IntVector& ghostcell_dims_density,
    const hier::IntVector& ghostcell_dims_temperature,
    const hier::IntVector& ghostcell_dims_thermo_properties,
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
        
        const int offset_0_internal_energy = offset_internal_energy[0];
        const int offset_0_density = offset_density[0];
        const int offset_0_temperature = offset_temperature[0];
        const int offset_0_thermo_properties = offset_thermo_properties[0];
        
        HAMERS_PRAGMA_SIMD
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_internal_energy = i + offset_0_internal_energy;
            const int idx_density = i + offset_0_density;
            const int idx_temperature = i + offset_0_temperature;
            const int idx_thermo_properties = i + offset_0_thermo_properties;
            
            epsilon[idx_internal_energy] = c_v[idx_thermo_properties]*T[idx_temperature] +
                p_inf[idx_thermo_properties]/rho[idx_density];
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
        
        const int offset_0_internal_energy = offset_internal_energy[0];
        const int offset_1_internal_energy = offset_internal_energy[1];
        const int ghostcell_dim_0_internal_energy = ghostcell_dims_internal_energy[0];
        
        const int offset_0_density = offset_density[0];
        const int offset_1_density = offset_density[1];
        const int ghostcell_dim_0_density = ghostcell_dims_density[0];
        
        const int offset_0_temperature = offset_temperature[0];
        const int offset_1_temperature = offset_temperature[1];
        const int ghostcell_dim_0_temperature = ghostcell_dims_temperature[0];
        
        const int offset_0_thermo_properties = offset_thermo_properties[0];
        const int offset_1_thermo_properties = offset_thermo_properties[1];
        const int ghostcell_dim_0_thermo_properties = ghostcell_dims_thermo_properties[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_internal_energy = (i + offset_0_internal_energy) +
                    (j + offset_1_internal_energy)*ghostcell_dim_0_internal_energy;
                
                const int idx_density = (i + offset_0_density) +
                    (j + offset_1_density)*ghostcell_dim_0_density;
                
                const int idx_temperature = (i + offset_0_temperature) +
                    (j + offset_1_temperature)*ghostcell_dim_0_temperature;
                
                const int idx_thermo_properties = (i + offset_0_thermo_properties) +
                    (j + offset_1_thermo_properties)*ghostcell_dim_0_thermo_properties;
                
                epsilon[idx_internal_energy] = c_v[idx_thermo_properties]*T[idx_temperature] +
                    p_inf[idx_thermo_properties]/rho[idx_density];
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
        
        const int offset_0_internal_energy = offset_internal_energy[0];
        const int offset_1_internal_energy = offset_internal_energy[1];
        const int offset_2_internal_energy = offset_internal_energy[2];
        const int ghostcell_dim_0_internal_energy = ghostcell_dims_internal_energy[0];
        const int ghostcell_dim_1_internal_energy = ghostcell_dims_internal_energy[1];
        
        const int offset_0_density = offset_density[0];
        const int offset_1_density = offset_density[1];
        const int offset_2_density = offset_density[2];
        const int ghostcell_dim_0_density = ghostcell_dims_density[0];
        const int ghostcell_dim_1_density = ghostcell_dims_density[1];
        
        const int offset_0_temperature = offset_temperature[0];
        const int offset_1_temperature = offset_temperature[1];
        const int offset_2_temperature = offset_temperature[2];
        const int ghostcell_dim_0_temperature = ghostcell_dims_temperature[0];
        const int ghostcell_dim_1_temperature = ghostcell_dims_temperature[1];
        
        const int offset_0_thermo_properties = offset_thermo_properties[0];
        const int offset_1_thermo_properties = offset_thermo_properties[1];
        const int offset_2_thermo_properties = offset_thermo_properties[2];
        const int ghostcell_dim_0_thermo_properties = ghostcell_dims_thermo_properties[0];
        const int ghostcell_dim_1_thermo_properties = ghostcell_dims_thermo_properties[1];
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_internal_energy = (i + offset_0_internal_energy) +
                        (j + offset_1_internal_energy)*ghostcell_dim_0_internal_energy +
                        (k + offset_2_internal_energy)*ghostcell_dim_0_internal_energy*
                            ghostcell_dim_1_internal_energy;
                    
                    const int idx_density = (i + offset_0_density) +
                        (j + offset_1_density)*ghostcell_dim_0_density +
                        (k + offset_2_density)*ghostcell_dim_0_density*
                            ghostcell_dim_1_density;
                    
                    const int idx_temperature = (i + offset_0_temperature) +
                        (j + offset_1_temperature)*ghostcell_dim_0_temperature +
                        (k + offset_2_temperature)*ghostcell_dim_0_temperature*
                            ghostcell_dim_1_temperature;
                    
                    const int idx_thermo_properties = (i + offset_0_thermo_properties) +
                        (j + offset_1_thermo_properties)*ghostcell_dim_0_thermo_properties +
                        (k + offset_2_thermo_properties)*ghostcell_dim_0_thermo_properties*
                            ghostcell_dim_1_thermo_properties;
                    
                    epsilon[idx_internal_energy] = c_v[idx_thermo_properties]*T[idx_temperature] +
                        p_inf[idx_thermo_properties]/rho[idx_density];
                }
            }
        }
    }
}


/*
 * Compute the isochoric specific heat capacity.
 */
void
EquationOfStateStiffenedGas::computeIsochoricSpecificHeatCapacity(
    Real* const c_v,
    const Real& c_v_src,
    const hier::IntVector& offset_isochoric_specific_heat_capacity,
    const hier::IntVector& ghostcell_dims_isochoric_specific_heat_capacity,
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
        
        HAMERS_PRAGMA_SIMD
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_isochoric_specific_heat_capacity =
                i + offset_0_isochoric_specific_heat_capacity;
            
            c_v[idx_isochoric_specific_heat_capacity] = c_v_src;
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
        
        const int offset_0_isochoric_specific_heat_capacity =
            offset_isochoric_specific_heat_capacity[0];
        const int offset_1_isochoric_specific_heat_capacity =
            offset_isochoric_specific_heat_capacity[1];
        const int ghostcell_dim_0_isochoric_specific_heat_capacity =
            ghostcell_dims_isochoric_specific_heat_capacity[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_isochoric_specific_heat_capacity =
                    (i + offset_0_isochoric_specific_heat_capacity) +
                    (j + offset_1_isochoric_specific_heat_capacity)*
                        ghostcell_dim_0_isochoric_specific_heat_capacity;
                
                c_v[idx_isochoric_specific_heat_capacity] = c_v_src;
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
        
        const int offset_0_isochoric_specific_heat_capacity =
            offset_isochoric_specific_heat_capacity[0];
        const int offset_1_isochoric_specific_heat_capacity =
            offset_isochoric_specific_heat_capacity[1];
        const int offset_2_isochoric_specific_heat_capacity =
            offset_isochoric_specific_heat_capacity[2];
        const int ghostcell_dim_0_isochoric_specific_heat_capacity =
            ghostcell_dims_isochoric_specific_heat_capacity[0];
        const int ghostcell_dim_1_isochoric_specific_heat_capacity =
            ghostcell_dims_isochoric_specific_heat_capacity[1];
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
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
EquationOfStateStiffenedGas::computeIsochoricSpecificHeatCapacity(
    Real* const c_v,
    const Real* const c_v_src,
    const hier::IntVector& offset_isochoric_specific_heat_capacity,
    const hier::IntVector& offset_thermo_properties,
    const hier::IntVector& ghostcell_dims_isochoric_specific_heat_capacity,
    const hier::IntVector& ghostcell_dims_thermo_properties,
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
        
        const int offset_0_thermo_properties = offset_thermo_properties[0];
        
        HAMERS_PRAGMA_SIMD
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_isochoric_specific_heat_capacity =
                i + offset_0_isochoric_specific_heat_capacity;
            
            const int idx_thermo_properties = i + offset_0_thermo_properties;
            
            c_v[idx_isochoric_specific_heat_capacity] = c_v_src[idx_thermo_properties];
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
        
        const int offset_0_isochoric_specific_heat_capacity =
            offset_isochoric_specific_heat_capacity[0];
        const int offset_1_isochoric_specific_heat_capacity =
            offset_isochoric_specific_heat_capacity[1];
        const int ghostcell_dim_0_isochoric_specific_heat_capacity =
            ghostcell_dims_isochoric_specific_heat_capacity[0];
        
        const int offset_0_thermo_properties = offset_thermo_properties[0];
        const int offset_1_thermo_properties = offset_thermo_properties[1];
        const int ghostcell_dim_0_thermo_properties = ghostcell_dims_thermo_properties[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_isochoric_specific_heat_capacity =
                    (i + offset_0_isochoric_specific_heat_capacity) +
                    (j + offset_1_isochoric_specific_heat_capacity)*
                        ghostcell_dim_0_isochoric_specific_heat_capacity;
                
                const int idx_thermo_properties = (i + offset_0_thermo_properties) +
                    (j + offset_1_thermo_properties)*ghostcell_dim_0_thermo_properties;
                
                c_v[idx_isochoric_specific_heat_capacity] = c_v_src[idx_thermo_properties];
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
        
        const int offset_0_isochoric_specific_heat_capacity =
            offset_isochoric_specific_heat_capacity[0];
        const int offset_1_isochoric_specific_heat_capacity =
            offset_isochoric_specific_heat_capacity[1];
        const int offset_2_isochoric_specific_heat_capacity =
            offset_isochoric_specific_heat_capacity[2];
        const int ghostcell_dim_0_isochoric_specific_heat_capacity =
            ghostcell_dims_isochoric_specific_heat_capacity[0];
        const int ghostcell_dim_1_isochoric_specific_heat_capacity =
            ghostcell_dims_isochoric_specific_heat_capacity[1];
        
        const int offset_0_thermo_properties = offset_thermo_properties[0];
        const int offset_1_thermo_properties = offset_thermo_properties[1];
        const int offset_2_thermo_properties = offset_thermo_properties[2];
        const int ghostcell_dim_0_thermo_properties = ghostcell_dims_thermo_properties[0];
        const int ghostcell_dim_1_thermo_properties = ghostcell_dims_thermo_properties[1];
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
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
                    
                    const int idx_thermo_properties = (i + offset_0_thermo_properties) +
                        (j + offset_1_thermo_properties)*ghostcell_dim_0_thermo_properties +
                        (k + offset_2_thermo_properties)*ghostcell_dim_0_thermo_properties*
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
EquationOfStateStiffenedGas::computeIsobaricSpecificHeatCapacity(
    Real* const c_p,
    const Real& c_p_src,
    const hier::IntVector& offset_isobaric_specific_heat_capacity,
    const hier::IntVector& ghostcell_dims_isobaric_specific_heat_capacity,
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
        
        HAMERS_PRAGMA_SIMD
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_isobaric_specific_heat_capacity =
                i + offset_0_isobaric_specific_heat_capacity;
            
            c_p[idx_isobaric_specific_heat_capacity] = c_p_src;
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
        
        const int offset_0_isobaric_specific_heat_capacity =
            offset_isobaric_specific_heat_capacity[0];
        const int offset_1_isobaric_specific_heat_capacity =
            offset_isobaric_specific_heat_capacity[1];
        const int ghostcell_dim_0_isobaric_specific_heat_capacity =
            ghostcell_dims_isobaric_specific_heat_capacity[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_isobaric_specific_heat_capacity =
                    (i + offset_0_isobaric_specific_heat_capacity) +
                    (j + offset_1_isobaric_specific_heat_capacity)*
                        ghostcell_dim_0_isobaric_specific_heat_capacity;
                
                c_p[idx_isobaric_specific_heat_capacity] = c_p_src;
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
        
        const int offset_0_isobaric_specific_heat_capacity =
            offset_isobaric_specific_heat_capacity[0];
        const int offset_1_isobaric_specific_heat_capacity =
            offset_isobaric_specific_heat_capacity[1];
        const int offset_2_isobaric_specific_heat_capacity =
            offset_isobaric_specific_heat_capacity[2];
        const int ghostcell_dim_0_isobaric_specific_heat_capacity =
            ghostcell_dims_isobaric_specific_heat_capacity[0];
        const int ghostcell_dim_1_isobaric_specific_heat_capacity =
            ghostcell_dims_isobaric_specific_heat_capacity[1];
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
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
EquationOfStateStiffenedGas::computeIsobaricSpecificHeatCapacity(
    Real* const c_p,
    const Real* const c_p_src,
    const hier::IntVector& offset_isobaric_specific_heat_capacity,
    const hier::IntVector& offset_thermo_properties,
    const hier::IntVector& ghostcell_dims_isobaric_specific_heat_capacity,
    const hier::IntVector& ghostcell_dims_thermo_properties,
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
        
        const int offset_0_thermo_properties = offset_thermo_properties[0];
        
        HAMERS_PRAGMA_SIMD
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_isobaric_specific_heat_capacity =
                i + offset_0_isobaric_specific_heat_capacity;
            
            const int idx_thermo_properties = i + offset_0_thermo_properties;
            
            c_p[idx_isobaric_specific_heat_capacity] = c_p_src[idx_thermo_properties];
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
        
        const int offset_0_isobaric_specific_heat_capacity =
            offset_isobaric_specific_heat_capacity[0];
        const int offset_1_isobaric_specific_heat_capacity =
            offset_isobaric_specific_heat_capacity[1];
        const int ghostcell_dim_0_isobaric_specific_heat_capacity =
            ghostcell_dims_isobaric_specific_heat_capacity[0];
        
        const int offset_0_thermo_properties = offset_thermo_properties[0];
        const int offset_1_thermo_properties = offset_thermo_properties[1];
        const int ghostcell_dim_0_thermo_properties = ghostcell_dims_thermo_properties[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_isobaric_specific_heat_capacity =
                    (i + offset_0_isobaric_specific_heat_capacity) +
                    (j + offset_1_isobaric_specific_heat_capacity)*
                        ghostcell_dim_0_isobaric_specific_heat_capacity;
                
                const int idx_thermo_properties = (i + offset_0_thermo_properties) +
                    (j + offset_1_thermo_properties)*ghostcell_dim_0_thermo_properties;
                
                c_p[idx_isobaric_specific_heat_capacity] = c_p_src[idx_thermo_properties];
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
        
        const int offset_0_isobaric_specific_heat_capacity =
            offset_isobaric_specific_heat_capacity[0];
        const int offset_1_isobaric_specific_heat_capacity =
            offset_isobaric_specific_heat_capacity[1];
        const int offset_2_isobaric_specific_heat_capacity =
            offset_isobaric_specific_heat_capacity[2];
        const int ghostcell_dim_0_isobaric_specific_heat_capacity =
            ghostcell_dims_isobaric_specific_heat_capacity[0];
        const int ghostcell_dim_1_isobaric_specific_heat_capacity =
            ghostcell_dims_isobaric_specific_heat_capacity[1];
        
        const int offset_0_thermo_properties = offset_thermo_properties[0];
        const int offset_1_thermo_properties = offset_thermo_properties[1];
        const int offset_2_thermo_properties = offset_thermo_properties[2];
        const int ghostcell_dim_0_thermo_properties = ghostcell_dims_thermo_properties[0];
        const int ghostcell_dim_1_thermo_properties = ghostcell_dims_thermo_properties[1];
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
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
                    
                    const int idx_thermo_properties = (i + offset_0_thermo_properties) +
                        (j + offset_1_thermo_properties)*ghostcell_dim_0_thermo_properties +
                        (k + offset_2_thermo_properties)*ghostcell_dim_0_thermo_properties*
                            ghostcell_dim_1_thermo_properties;
                    
                    c_p[idx_isobaric_specific_heat_capacity] = c_p_src[idx_thermo_properties];
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
EquationOfStateStiffenedGas::computeGruneisenParameter(
    Real* const Gamma,
    const Real& gamma,
    const hier::IntVector& offset_gruneisen_parameter,
    const hier::IntVector& ghostcell_dims_gruneisen_parameter,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims) const
{
    const Real Gamma_src = gamma - Real(1);
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and offsets.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int offset_0_gruneisen_parameter =
            offset_gruneisen_parameter[0];
        
        HAMERS_PRAGMA_SIMD
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_gruneisen_parameter =
                i + offset_0_gruneisen_parameter;
            
            Gamma[idx_gruneisen_parameter] = Gamma_src;
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
        
        const int offset_0_gruneisen_parameter =
            offset_gruneisen_parameter[0];
        const int offset_1_gruneisen_parameter =
            offset_gruneisen_parameter[1];
        const int ghostcell_dim_0_gruneisen_parameter =
            ghostcell_dims_gruneisen_parameter[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_gruneisen_parameter =
                    (i + offset_0_gruneisen_parameter) +
                    (j + offset_1_gruneisen_parameter)*
                        ghostcell_dim_0_gruneisen_parameter;
                
                Gamma[idx_gruneisen_parameter] = Gamma_src;
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
        
        const int offset_0_gruneisen_parameter =
            offset_gruneisen_parameter[0];
        const int offset_1_gruneisen_parameter =
            offset_gruneisen_parameter[1];
        const int offset_2_gruneisen_parameter =
            offset_gruneisen_parameter[2];
        const int ghostcell_dim_0_gruneisen_parameter =
            ghostcell_dims_gruneisen_parameter[0];
        const int ghostcell_dim_1_gruneisen_parameter =
            ghostcell_dims_gruneisen_parameter[1];
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_gruneisen_parameter =
                        (i + offset_0_gruneisen_parameter) +
                        (j + offset_1_gruneisen_parameter)*
                            ghostcell_dim_0_gruneisen_parameter +
                        (k + offset_2_gruneisen_parameter)*
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
EquationOfStateStiffenedGas::computeGruneisenParameter(
    Real* const Gamma,
    const Real* const gamma,
    const hier::IntVector& offset_gruneisen_parameter,
    const hier::IntVector& offset_thermo_properties,
    const hier::IntVector& ghostcell_dims_gruneisen_parameter,
    const hier::IntVector& ghostcell_dims_thermo_properties,
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
        
        const int offset_0_gruneisen_parameter =
            offset_gruneisen_parameter[0];
        
        const int offset_0_thermo_properties = offset_thermo_properties[0];
        
        HAMERS_PRAGMA_SIMD
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_gruneisen_parameter =
                i + offset_0_gruneisen_parameter;
            
            const int idx_thermo_properties = i + offset_0_thermo_properties;
            
            Gamma[idx_gruneisen_parameter] = gamma[idx_thermo_properties] - Real(1);
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
        
        const int offset_0_gruneisen_parameter =
            offset_gruneisen_parameter[0];
        const int offset_1_gruneisen_parameter =
            offset_gruneisen_parameter[1];
        const int ghostcell_dim_0_gruneisen_parameter =
            ghostcell_dims_gruneisen_parameter[0];
        
        const int offset_0_thermo_properties = offset_thermo_properties[0];
        const int offset_1_thermo_properties = offset_thermo_properties[1];
        const int ghostcell_dim_0_thermo_properties = ghostcell_dims_thermo_properties[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_gruneisen_parameter =
                    (i + offset_0_gruneisen_parameter) +
                    (j + offset_1_gruneisen_parameter)*
                        ghostcell_dim_0_gruneisen_parameter;
                
                const int idx_thermo_properties = (i + offset_0_thermo_properties) +
                    (j + offset_1_thermo_properties)*ghostcell_dim_0_thermo_properties;
                
                Gamma[idx_gruneisen_parameter] = gamma[idx_thermo_properties] - Real(1);
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
        
        const int offset_0_gruneisen_parameter =
            offset_gruneisen_parameter[0];
        const int offset_1_gruneisen_parameter =
            offset_gruneisen_parameter[1];
        const int offset_2_gruneisen_parameter =
            offset_gruneisen_parameter[2];
        const int ghostcell_dim_0_gruneisen_parameter =
            ghostcell_dims_gruneisen_parameter[0];
        const int ghostcell_dim_1_gruneisen_parameter =
            ghostcell_dims_gruneisen_parameter[1];
        
        const int offset_0_thermo_properties = offset_thermo_properties[0];
        const int offset_1_thermo_properties = offset_thermo_properties[1];
        const int offset_2_thermo_properties = offset_thermo_properties[2];
        const int ghostcell_dim_0_thermo_properties = ghostcell_dims_thermo_properties[0];
        const int ghostcell_dim_1_thermo_properties = ghostcell_dims_thermo_properties[1];
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_gruneisen_parameter =
                        (i + offset_0_gruneisen_parameter) +
                        (j + offset_1_gruneisen_parameter)*
                            ghostcell_dim_0_gruneisen_parameter +
                        (k + offset_2_gruneisen_parameter)*
                            ghostcell_dim_0_gruneisen_parameter*
                            ghostcell_dim_1_gruneisen_parameter;
                    
                    const int idx_thermo_properties = (i + offset_0_thermo_properties) +
                        (j + offset_1_thermo_properties)*ghostcell_dim_0_thermo_properties +
                        (k + offset_2_thermo_properties)*ghostcell_dim_0_thermo_properties*
                            ghostcell_dim_1_thermo_properties;
                    
                    Gamma[idx_gruneisen_parameter] = gamma[idx_thermo_properties] - Real(1);
                }
            }
        }
    }
}


/*
 * Compute the partial derivative of pressure w.r.t. density under constant specific internal energy.
 */
void
EquationOfStateStiffenedGas::computePressureDerivativeWithDensity(
    Real* const Psi,
    const Real* const rho,
    const Real* const p,
    const Real& gamma,
    const Real& p_inf,
    const hier::IntVector& offset_partial_pressure_partial_density,
    const hier::IntVector& offset_density,
    const hier::IntVector& offset_pressure,
    const hier::IntVector& ghostcell_dims_partial_pressure_partial_density,
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
        
        const int offset_0_partial_pressure_partial_density = offset_partial_pressure_partial_density[0];
        const int offset_0_density = offset_density[0];
        const int offset_0_pressure = offset_pressure[0];
        
        HAMERS_PRAGMA_SIMD
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_partial_pressure_partial_density = i + offset_0_partial_pressure_partial_density;
            const int idx_density = i + offset_0_density;
            const int idx_pressure = i + offset_0_pressure;
            
            Psi[idx_partial_pressure_partial_density] = (p[idx_pressure] + gamma*p_inf)/rho[idx_density];
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
        
        const int offset_0_partial_pressure_partial_density = offset_partial_pressure_partial_density[0];
        const int offset_1_partial_pressure_partial_density = offset_partial_pressure_partial_density[1];
        const int ghostcell_dim_0_partial_pressure_partial_density = ghostcell_dims_partial_pressure_partial_density[0];
        
        const int offset_0_density = offset_density[0];
        const int offset_1_density = offset_density[1];
        const int ghostcell_dim_0_density = ghostcell_dims_density[0];
        
        const int offset_0_pressure = offset_pressure[0];
        const int offset_1_pressure = offset_pressure[1];
        const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_partial_pressure_partial_density = (i + offset_0_partial_pressure_partial_density) +
                    (j + offset_1_partial_pressure_partial_density)*ghostcell_dim_0_partial_pressure_partial_density;
                
                const int idx_density = (i + offset_0_density) +
                    (j + offset_1_density)*ghostcell_dim_0_density;
                
                const int idx_pressure = (i + offset_0_pressure) +
                    (j + offset_1_pressure)*ghostcell_dim_0_pressure;
                
                Psi[idx_partial_pressure_partial_density] = (p[idx_pressure] + gamma*p_inf)/rho[idx_density];
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
        
        const int offset_0_partial_pressure_partial_density = offset_partial_pressure_partial_density[0];
        const int offset_1_partial_pressure_partial_density = offset_partial_pressure_partial_density[1];
        const int offset_2_partial_pressure_partial_density = offset_partial_pressure_partial_density[2];
        const int ghostcell_dim_0_partial_pressure_partial_density = ghostcell_dims_partial_pressure_partial_density[0];
        const int ghostcell_dim_1_partial_pressure_partial_density = ghostcell_dims_partial_pressure_partial_density[1];
        
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
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_partial_pressure_partial_density = (i + offset_0_partial_pressure_partial_density) +
                        (j + offset_1_partial_pressure_partial_density)*ghostcell_dim_0_partial_pressure_partial_density +
                        (k + offset_2_partial_pressure_partial_density)*ghostcell_dim_0_partial_pressure_partial_density*
                            ghostcell_dim_1_partial_pressure_partial_density;
                    
                    const int idx_density = (i + offset_0_density) +
                        (j + offset_1_density)*ghostcell_dim_0_density +
                        (k + offset_2_density)*ghostcell_dim_0_density*
                            ghostcell_dim_1_density;
                    
                    const int idx_pressure = (i + offset_0_pressure) +
                        (j + offset_1_pressure)*ghostcell_dim_0_pressure +
                        (k + offset_2_pressure)*ghostcell_dim_0_pressure*
                            ghostcell_dim_1_pressure;
                    
                    Psi[idx_partial_pressure_partial_density] = (p[idx_pressure] + gamma*p_inf)/rho[idx_density];
                }
            }
        }
    }
}


/*
 * Compute the partial derivative of pressure w.r.t. density under constant specific internal energy.
 */
void
EquationOfStateStiffenedGas::computePressureDerivativeWithDensity(
    Real* const Psi,
    const Real* const rho,
    const Real* const p,
    const Real* const gamma,
    const Real* const p_inf,
    const hier::IntVector& offset_partial_pressure_partial_density,
    const hier::IntVector& offset_density,
    const hier::IntVector& offset_pressure,
    const hier::IntVector& offset_thermo_properties,
    const hier::IntVector& ghostcell_dims_partial_pressure_partial_density,
    const hier::IntVector& ghostcell_dims_density,
    const hier::IntVector& ghostcell_dims_pressure,
    const hier::IntVector& ghostcell_dims_thermo_properties,
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
        
        const int offset_0_partial_pressure_partial_density = offset_partial_pressure_partial_density[0];
        const int offset_0_density = offset_density[0];
        const int offset_0_pressure = offset_pressure[0];
        const int offset_0_thermo_properties = offset_thermo_properties[0];
        
        HAMERS_PRAGMA_SIMD
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_partial_pressure_partial_density = i + offset_0_partial_pressure_partial_density;
            const int idx_density = i + offset_0_density;
            const int idx_pressure = i + offset_0_pressure;
            const int idx_thermo_properties = i + offset_0_thermo_properties;
            
            Psi[idx_partial_pressure_partial_density] = (p[idx_pressure] + gamma[idx_thermo_properties]*
                p_inf[idx_thermo_properties])/rho[idx_density];
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
        
        const int offset_0_partial_pressure_partial_density = offset_partial_pressure_partial_density[0];
        const int offset_1_partial_pressure_partial_density = offset_partial_pressure_partial_density[1];
        const int ghostcell_dim_0_partial_pressure_partial_density = ghostcell_dims_partial_pressure_partial_density[0];
        
        const int offset_0_density = offset_density[0];
        const int offset_1_density = offset_density[1];
        const int ghostcell_dim_0_density = ghostcell_dims_density[0];
        
        const int offset_0_pressure = offset_pressure[0];
        const int offset_1_pressure = offset_pressure[1];
        const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
        
        const int offset_0_thermo_properties = offset_thermo_properties[0];
        const int offset_1_thermo_properties = offset_thermo_properties[1];
        const int ghostcell_dim_0_thermo_properties = ghostcell_dims_thermo_properties[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_partial_pressure_partial_density = (i + offset_0_partial_pressure_partial_density) +
                    (j + offset_1_partial_pressure_partial_density)*ghostcell_dim_0_partial_pressure_partial_density;
                
                const int idx_density = (i + offset_0_density) +
                    (j + offset_1_density)*ghostcell_dim_0_density;
                
                const int idx_pressure = (i + offset_0_pressure) +
                    (j + offset_1_pressure)*ghostcell_dim_0_pressure;
                
                const int idx_thermo_properties = (i + offset_0_thermo_properties) +
                    (j + offset_1_thermo_properties)*ghostcell_dim_0_thermo_properties;
                
                Psi[idx_partial_pressure_partial_density] = (p[idx_pressure] + gamma[idx_thermo_properties]*
                    p_inf[idx_thermo_properties])/rho[idx_density];
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
        
        const int offset_0_partial_pressure_partial_density = offset_partial_pressure_partial_density[0];
        const int offset_1_partial_pressure_partial_density = offset_partial_pressure_partial_density[1];
        const int offset_2_partial_pressure_partial_density = offset_partial_pressure_partial_density[2];
        const int ghostcell_dim_0_partial_pressure_partial_density = ghostcell_dims_partial_pressure_partial_density[0];
        const int ghostcell_dim_1_partial_pressure_partial_density = ghostcell_dims_partial_pressure_partial_density[1];
        
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
        
        const int offset_0_thermo_properties = offset_thermo_properties[0];
        const int offset_1_thermo_properties = offset_thermo_properties[1];
        const int offset_2_thermo_properties = offset_thermo_properties[2];
        const int ghostcell_dim_0_thermo_properties = ghostcell_dims_thermo_properties[0];
        const int ghostcell_dim_1_thermo_properties = ghostcell_dims_thermo_properties[1];
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_partial_pressure_partial_density = (i + offset_0_partial_pressure_partial_density) +
                        (j + offset_1_partial_pressure_partial_density)*ghostcell_dim_0_partial_pressure_partial_density +
                        (k + offset_2_partial_pressure_partial_density)*ghostcell_dim_0_partial_pressure_partial_density*
                            ghostcell_dim_1_partial_pressure_partial_density;
                    
                    const int idx_density = (i + offset_0_density) +
                        (j + offset_1_density)*ghostcell_dim_0_density +
                        (k + offset_2_density)*ghostcell_dim_0_density*
                            ghostcell_dim_1_density;
                    
                    const int idx_pressure = (i + offset_0_pressure) +
                        (j + offset_1_pressure)*ghostcell_dim_0_pressure +
                        (k + offset_2_pressure)*ghostcell_dim_0_pressure*
                            ghostcell_dim_1_pressure;
                    
                    const int idx_thermo_properties = (i + offset_0_thermo_properties) +
                        (j + offset_1_thermo_properties)*ghostcell_dim_0_thermo_properties +
                        (k + offset_2_thermo_properties)*ghostcell_dim_0_thermo_properties*
                            ghostcell_dim_1_thermo_properties;
                    
                    Psi[idx_partial_pressure_partial_density] = (p[idx_pressure] + gamma[idx_thermo_properties]*
                        p_inf[idx_thermo_properties])/rho[idx_density];
                }
            }
        }
    }
}


/*
 * Compute the density.
 */
void
EquationOfStateStiffenedGas::computeDensity(
    Real* const rho,
    const Real* const p,
    const Real* const T,
    const Real& gamma,
    const Real& p_inf,
    const Real& c_v,
    const hier::IntVector& offset_density,
    const hier::IntVector& offset_pressure,
    const hier::IntVector& offset_temperature,
    const hier::IntVector& ghostcell_dims_density,
    const hier::IntVector& ghostcell_dims_pressure,
    const hier::IntVector& ghostcell_dims_temperature,
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
        
        const int offset_0_density = offset_density[0];
        const int offset_0_pressure = offset_pressure[0];
        const int offset_0_temperature = offset_temperature[0];
        
        HAMERS_PRAGMA_SIMD
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_density = i + offset_0_density;
            const int idx_pressure = i + offset_0_pressure;
            const int idx_temperature = i + offset_0_temperature;
            
            rho[idx_density] = (p[idx_pressure] + p_inf)/((gamma - Real(1))*c_v*T[idx_temperature]);
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
        
        const int offset_0_density = offset_density[0];
        const int offset_1_density = offset_density[1];
        const int ghostcell_dim_0_density = ghostcell_dims_density[0];
        
        const int offset_0_pressure = offset_pressure[0];
        const int offset_1_pressure = offset_pressure[1];
        const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
        
        const int offset_0_temperature = offset_temperature[0];
        const int offset_1_temperature = offset_temperature[1];
        const int ghostcell_dim_0_temperature = ghostcell_dims_temperature[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_density = (i + offset_0_density) +
                    (j + offset_1_density)*ghostcell_dim_0_density;
                
                const int idx_pressure = (i + offset_0_pressure) +
                    (j + offset_1_pressure)*ghostcell_dim_0_pressure;
                
                const int idx_temperature = (i + offset_0_temperature) +
                    (j + offset_1_temperature)*ghostcell_dim_0_temperature;
                
                rho[idx_density] = (p[idx_pressure] + p_inf)/((gamma - Real(1))*c_v*T[idx_temperature]);
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
        
        const int offset_0_temperature = offset_temperature[0];
        const int offset_1_temperature = offset_temperature[1];
        const int offset_2_temperature = offset_temperature[2];
        const int ghostcell_dim_0_temperature = ghostcell_dims_temperature[0];
        const int ghostcell_dim_1_temperature = ghostcell_dims_temperature[1];
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_density = (i + offset_0_density) +
                        (j + offset_1_density)*ghostcell_dim_0_density +
                        (k + offset_2_density)*ghostcell_dim_0_density*
                            ghostcell_dim_1_density;
                    
                    const int idx_pressure = (i + offset_0_pressure) +
                        (j + offset_1_pressure)*ghostcell_dim_0_pressure +
                        (k + offset_2_pressure)*ghostcell_dim_0_pressure*
                            ghostcell_dim_1_pressure;
                    
                    const int idx_temperature = (i + offset_0_temperature) +
                        (j + offset_1_temperature)*ghostcell_dim_0_temperature +
                        (k + offset_2_temperature)*ghostcell_dim_0_temperature*
                            ghostcell_dim_1_temperature;
                    
                    rho[idx_density] = (p[idx_pressure] + p_inf)/((gamma - Real(1))*c_v*T[idx_temperature]);
                }
            }
        }
    }
}


/*
 * Compute the density.
 */
void
EquationOfStateStiffenedGas::computeDensity(
    Real* const rho,
    const Real* const p,
    const Real* const T,
    const Real* const gamma,
    const Real* const p_inf,
    const Real* const c_v,
    const hier::IntVector& offset_density,
    const hier::IntVector& offset_pressure,
    const hier::IntVector& offset_temperature,
    const hier::IntVector& offset_thermo_properties,
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
         * Get the local lower index, numbers of cells in each dimension and offsets.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int offset_0_density = offset_density[0];
        const int offset_0_pressure = offset_pressure[0];
        const int offset_0_temperature = offset_temperature[0];
        const int offset_0_thermo_properties = offset_thermo_properties[0];
        
        HAMERS_PRAGMA_SIMD
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_density = i + offset_0_density;
            const int idx_pressure = i + offset_0_pressure;
            const int idx_temperature = i + offset_0_temperature;
            const int idx_thermo_properties = i + offset_0_thermo_properties;
            
            rho[idx_density] = (p[idx_pressure] + p_inf[idx_thermo_properties])/
                ((gamma[idx_thermo_properties] - Real(1))*c_v[idx_thermo_properties]*T[idx_temperature]);
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
        
        const int offset_0_density = offset_density[0];
        const int offset_1_density = offset_density[1];
        const int ghostcell_dim_0_density = ghostcell_dims_density[0];
        
        const int offset_0_pressure = offset_pressure[0];
        const int offset_1_pressure = offset_pressure[1];
        const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
        
        const int offset_0_temperature = offset_temperature[0];
        const int offset_1_temperature = offset_temperature[1];
        const int ghostcell_dim_0_temperature = ghostcell_dims_temperature[0];
        
        const int offset_0_thermo_properties = offset_thermo_properties[0];
        const int offset_1_thermo_properties = offset_thermo_properties[1];
        const int ghostcell_dim_0_thermo_properties = ghostcell_dims_thermo_properties[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_density = (i + offset_0_density) +
                    (j + offset_1_density)*ghostcell_dim_0_density;
                
                const int idx_pressure = (i + offset_0_pressure) +
                    (j + offset_1_pressure)*ghostcell_dim_0_pressure;
                
                const int idx_temperature = (i + offset_0_temperature) +
                    (j + offset_1_temperature)*ghostcell_dim_0_temperature;
                
                const int idx_thermo_properties = (i + offset_0_thermo_properties) +
                    (j + offset_1_thermo_properties)*ghostcell_dim_0_thermo_properties;
                
                rho[idx_density] = (p[idx_pressure] + p_inf[idx_thermo_properties])/
                    ((gamma[idx_thermo_properties] - Real(1))*c_v[idx_thermo_properties]*T[idx_temperature]);
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
        
        const int offset_0_temperature = offset_temperature[0];
        const int offset_1_temperature = offset_temperature[1];
        const int offset_2_temperature = offset_temperature[2];
        const int ghostcell_dim_0_temperature = ghostcell_dims_temperature[0];
        const int ghostcell_dim_1_temperature = ghostcell_dims_temperature[1];
        
        const int offset_0_thermo_properties = offset_thermo_properties[0];
        const int offset_1_thermo_properties = offset_thermo_properties[1];
        const int offset_2_thermo_properties = offset_thermo_properties[2];
        const int ghostcell_dim_0_thermo_properties = ghostcell_dims_thermo_properties[0];
        const int ghostcell_dim_1_thermo_properties = ghostcell_dims_thermo_properties[1];
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_density = (i + offset_0_density) +
                        (j + offset_1_density)*ghostcell_dim_0_density +
                        (k + offset_2_density)*ghostcell_dim_0_density*
                            ghostcell_dim_1_density;
                    
                    const int idx_pressure = (i + offset_0_pressure) +
                        (j + offset_1_pressure)*ghostcell_dim_0_pressure +
                        (k + offset_2_pressure)*ghostcell_dim_0_pressure*
                            ghostcell_dim_1_pressure;
                    
                    const int idx_temperature = (i + offset_0_temperature) +
                        (j + offset_1_temperature)*ghostcell_dim_0_temperature +
                        (k + offset_2_temperature)*ghostcell_dim_0_temperature*
                            ghostcell_dim_1_temperature;
                    
                    const int idx_thermo_properties = (i + offset_0_thermo_properties) +
                        (j + offset_1_thermo_properties)*ghostcell_dim_0_thermo_properties +
                        (k + offset_2_thermo_properties)*ghostcell_dim_0_thermo_properties*
                            ghostcell_dim_1_thermo_properties;
                    
                    rho[idx_density] = (p[idx_pressure] + p_inf[idx_thermo_properties])/
                        ((gamma[idx_thermo_properties] - Real(1))*c_v[idx_thermo_properties]*T[idx_temperature]);
                }
            }
        }
    }
}
